module cemodel

    implicit none

    real*8 :: minstellarmass = 0.1d0, maxstellarmass = 100.0d0
    real*8 :: maxtime = 1.3d10, initialgasmass = 1.4d5
    integer :: maxmodelnum = 10**6
    real*8 :: stellarmasslow, stellarmasshigh !used for integrals
    real*8 :: sfrefficiency = 0.0d0, maxchangemgas = 0.10, maxchangeser = 0.10
    real*8 :: sfrstartmgasabove = 1.0d2, sfrendmgasbelow = 1.0d-6
    real*8 :: initialcompscalez = 1.0d0
    real*8 :: timestepgrowthrate = 0.03d0
    character(len=50), parameter :: filespath = "./"
    real*8 :: initialtimestep = 1.0d2, mintimestep = 1.0d2, maxtimestep = 6.0d6
    real*8 :: maxerrorser = 5.d-6
    integer :: serintegralsteps = 80000
    real*8 :: imfnormfactor = 1.0 ! to be replaced during initialisation

    integer:: stepnum
    real*8 :: timestep

    type stellarmodeldata
        real*8 :: mass
        real*8 :: remnantmass
        real*8 :: stellarlifetime
        real*8 :: lifetimemassexponent ! lifetime = C * mass ^ exponent to match upper model set automatically on initialisation
        real*8, allocatable :: yield(:)
    end type

    type(stellarmodeldata), allocatable :: stellarmodel(:)

    type speciestype
        character (len=6) :: name ! h, d, he4, etc
        character (len=3) :: symbol ! h, he, c, fe, etc
        integer :: z ! number of protons (atomic number)
        integer :: a ! number of nucleons (mass number)
        logical :: yieldisrelative
        logical :: foundyields = .false.
        logical :: foundinitialvalue = .false.
    end type

    type(speciestype), allocatable :: species(:)

    type modelstate
        real*8 :: time
        real*8 :: gasmass
        real*8 :: starsmass
        real*8 :: sfr         ! star formation rate (Msun / year)
        real*8 :: ser         ! stellar ejection rate (Msun / year)
        real*8, allocatable :: speciesmassfrac(:)
    end type

    type(modelstate), allocatable :: model(:)

interface
   pure real*8 function fi(x, realargin, intargin)
    real*8, intent(in) :: x
    real*8, intent(in), optional :: realargin
    integer, intent(in), optional :: intargin
   end function fi
end interface

contains

pure real*8 function pastspeciesmassfrac(time,species)
    real*8, intent(in) :: time
    integer, intent(in) :: species
    integer :: i
    real*8 :: weight2

    pastspeciesmassfrac = model(1)%speciesmassfrac(species)
    do i = 2,stepnum
        if (model(i)%time > time) then
            weight2 = (time -  model(i-1)%time) / (model(i)%time - model(i-1)%time)
            pastspeciesmassfrac = (1-weight2) * model(i-1)%speciesmassfrac(species) +&
                weight2 * model(i)%speciesmassfrac(species)
            exit
        end if
    end do
end function

pure integer function findmodelnum(time)
    real*8, intent(in) :: time
    integer :: low, middle, high

    low = 1
    high = stepnum

    do while (high - low > 1)
        middle = (high + low) / 2
        if (model(middle)%time > time) then
            high = middle
        else
            low = middle
        end if
    end do
    findmodelnum = low
end function findmodelnum

! calculate log epsilon abundance of an element in the current timestep
pure real*8 function logepsilon(symbol)
    character (len=*), intent(in) :: symbol
    real*8 :: eldensitysum
    integer :: s, hydrogenindex

    hydrogenindex = 1
    eldensitysum = 0.0d0
    do s = 1,size(species)
        if (species(s)%z == 1 .and. species(s)%a == 1) then
            hydrogenindex = s
        end if

        if (trim(adjustl(species(s)%symbol)) == trim(adjustl(symbol))) then
            eldensitysum = eldensitysum + (model(stepnum)%speciesmassfrac(s) / species(s)%a)
        end if
    end do

    logepsilon = log10(eldensitysum / model(stepnum)%speciesmassfrac(hydrogenindex)) + 12.0d0
end function logepsilon

! calculate log epsilon abundance of C+N+O in the current timestep
pure real*8 function logepsiloncno()
    real*8 :: eldensitysum
    integer :: s, hydrogenindex

    hydrogenindex = 1
    eldensitysum = 0.0d0
    do s = 1,size(species)
        if (species(s)%z == 1 .and. species(s)%a == 1) then
            hydrogenindex = s
        end if

        if (trim(adjustl(species(s)%symbol)) == trim(adjustl('c')) .or. trim(adjustl(species(s)%symbol)) == trim(adjustl('n')) &
            .or. trim(adjustl(species(s)%symbol)) == trim(adjustl('o'))) then
            eldensitysum = eldensitysum + (model(stepnum)%speciesmassfrac(s) / species(s)%a)
        end if
    end do

    logepsiloncno = log10(eldensitysum / model(stepnum)%speciesmassfrac(hydrogenindex)) + 12.0d0
end function logepsiloncno


! calculate mass fraction of an element in the current timestep
pure real*8 function elmassfrac(symbol)
    character (len=*), intent(in) :: symbol
    real*8 :: elmassfracsum
    integer :: s

    elmassfracsum = 0.0d0
    do s = 1,size(species)
        if (trim(adjustl(species(s)%symbol)) == trim(adjustl(symbol))) then
            elmassfracsum = elmassfracsum + model(stepnum)%speciesmassfrac(s)
        end if
    end do

    elmassfrac = elmassfracsum
end function elmassfrac

! returns main sequence lifetime in years of a star
! with initial mass in solar masses
pure real*8 function stellarlifetime(initmass)
    implicit none
    real*8, intent(in) :: initmass
    integer :: i
    integer :: modelindex ! reference model to interpolate from
!    real*8 :: c

!    stellarlifetime = 1.d+10 * (mass ** (-3.1d0))

    ! interpolate lifetime from ev. model data
    if (initmass < stellarmodel(1)%mass) then
        modelindex = 1
    else if (initmass > stellarmodel(size(stellarmodel))%mass) then
!        stellarlifetime = stellarmodel(size(stellarmodel))%stellarlifetime * &
!            (initmass / stellarmodel(size(stellarmodel))%mass) ** (-3.5)
        modelindex = size(stellarmodel)
    else
        do i = 2, size(stellarmodel)
            if (stellarmodel(i)%mass >= initmass) then
                ! linear interpolation
!                c = (initmass - stellarmodel(i-1)%mass) / (stellarmodel(i)%mass - stellarmodel(i-1)%mass)
!                stellarlifetime = c * stellarmodel(i)%stellarlifetime + (1-c) * stellarmodel(i-1)%stellarlifetime
                modelindex = i - 1
                exit
            end if
        end do
    end if
    stellarlifetime = stellarmodel(modelindex)%stellarlifetime * &
            (initmass/stellarmodel(modelindex)%mass) ** stellarmodel(modelindex)%lifetimemassexponent
end function stellarlifetime

! remnant mass as a function of initial mass in solar masses
pure real*8 function remnantmass(initmass)
    implicit none
    real*8, intent(in) :: initmass
    integer :: i
    real*8 :: c

!   Iben & Tutukov 1984, in Pagel 2009 after eqn 7.10
!    if (mass <= 0.506) then
!        remnantmass = mass
!    else if (mass <= 9.5) then
!        remnantmass = 0.45d0 + 0.11d0 * mass
!    else
!        remnantmass = 1.5d0
!    end if

    if (initmass < stellarmodel(1)%mass) then
        c = (initmass - stellarmodel(1)%mass) / (stellarmodel(2)%mass - stellarmodel(1)%mass)
        remnantmass = (1-c) * stellarmodel(1)%remnantmass + c * stellarmodel(2)%remnantmass
        ! scale lowest mass end with initial mass
!        remnantmass = stellarmodel(1)%remnantmass * initmass / stellarmodel(1)%mass
    else if (initmass > stellarmodel(size(stellarmodel))%mass) then
        c = (stellarmodel(size(stellarmodel))%remnantmass - stellarmodel(size(stellarmodel)-1)%remnantmass) / &
            (stellarmodel(size(stellarmodel))%mass - stellarmodel(size(stellarmodel)-1)%mass)
        remnantmass = stellarmodel(size(stellarmodel))%remnantmass + c * (initmass - stellarmodel(size(stellarmodel))%mass)

!        remnantmass = stellarmodel(size(stellarmodel))%remnantmass
    else
    ! interpolate remnant mass from ev model data
        do i = 2, size(stellarmodel)
            if (stellarmodel(i)%mass >= initmass) then
                c = (initmass - stellarmodel(i-1)%mass) / (stellarmodel(i)%mass - stellarmodel(i-1)%mass)
                remnantmass = c * stellarmodel(i)%remnantmass + (1-c) * stellarmodel(i-1)%remnantmass
                exit
            end if
        end do
    end if
end function remnantmass

! relative or absolute yield of species in solar masses from a star with
! initial mass in solar masses
pure real*8 function yield(species, initmass)
    implicit none
    real*8, intent(in) :: initmass
    integer, intent(in) :: species
    integer :: i
    real*8 :: c

    yield = 0.0

    if ((initmass >= stellarmodel(1)%mass) .and. (initmass <= stellarmodel(size(stellarmodel))%mass)) then
        do i = 2, size(stellarmodel)
            if (stellarmodel(i)%mass >= initmass) then
                c = (initmass - stellarmodel(i-1)%mass) / (stellarmodel(i)%mass - stellarmodel(i-1)%mass)
                yield = (1-c) * stellarmodel(i-1)%yield(species) + c * stellarmodel(i)%yield(species)
                exit
            end if
        end do
    end if
end function yield

! stellar ejection rate in solar masses per year
pure real*8 function stellarejectionrate(time)
    implicit none
    real*8, intent(in) :: time

    stellarejectionrate = integral(dserbydmass, stellarmasshigh**(-3), stellarmasslow**(-3), serintegralsteps,&
        maxerrorser, realarg=time, label='ser mass integral')
end function stellarejectionrate

! D stellar ejection rate / D mass
pure real*8 function dserbydmass(mdashexp, time, dummyint)
    implicit none
    real*8, intent(in) :: mdashexp
    real*8, intent(in), optional :: time
    integer, intent(in), optional :: dummyint
    real*8 :: timeatbirth, mdash

    mdash = mdashexp ** (-1.d0/3.d0)

    timeatbirth = time - stellarlifetime(mdash)

    if (timeatbirth >= 0.0d0) then
        dserbydmass = (mdash - remnantmass(mdash)) * &
                    starformationrate(timeatbirth) * imf(mdash)
    else
        dserbydmass = 0.0d0
    end if
    dserbydmass = dserbydmass * (1.d0/3.d0) * mdashexp ** (-4.d0/3.d0) !required for integral change of variables
end function dserbydmass

! stellar ejection rate of species speciesnum in solar masses per year
pure real*8 function stellarejectionrateofspecies(time, speciesnum)
    implicit none
    real*8, intent(in) :: time
    integer, intent(in) :: speciesnum

    stellarejectionrateofspecies = integral(dserbydmassofspecies, stellarmasshigh**(-3), stellarmasslow**(-3), serintegralsteps,&
        maxerrorser, intarg=speciesnum, realarg=time, label='ser species mass integral')
end function

! D stellar ejection rate / D mass of species
pure real*8 function dserbydmassofspecies(mdashexp, time, speciesnum)
    implicit none
    real*8, intent(in) :: mdashexp
    real*8, intent(in), optional :: time
    integer, intent(in), optional :: speciesnum
    real*8 :: timeatbirth, mdash

    mdash = mdashexp ** (-1.d0/3.d0)

    timeatbirth = time - stellarlifetime(mdash)
    if (timeatbirth >= 0.0d0) then
        if (species(speciesnum)%yieldisrelative .eqv. .true.) then
            ! use relative yields
            dserbydmassofspecies = ((mdash - remnantmass(mdash)) * &
                pastspeciesmassfrac(timeatbirth,speciesnum) + &
                yield(speciesnum, mdash)) * &
                starformationrate(timeatbirth) * imf(mdash)
!            if (dserbydmassofspecies < 0.0d0) then
                !write(14,*) 'WARNING: Negative ejection rate, species:',speciesnum,', mass:',mdash
                !dserbydmassofspecies = 0.0d0
!            end if
        else
            ! use absolute yields
            dserbydmassofspecies = yield(speciesnum, mdash) * &
                starformationrate(timeatbirth) * imf(mdash)
        end if
    else
        dserbydmassofspecies = 0.0d0
    end if
    dserbydmassofspecies = dserbydmassofspecies * (1.d0/3.d0) * mdashexp ** (-4.d0/3.d0) !required for integral change of variables
end function

!integrate with Simpson's rule (and mid-point to measure error)
pure recursive function integral(f,x1,x2,numsteps,maxerror,intarg,realarg,label,recurdepth) result(integralresult)
    implicit none
!    real*8, external :: f
    procedure(fi) :: f
    integer, intent(in), optional :: intarg
    character(*), intent(in), optional :: label
    real*8, intent(in), optional :: realarg
    real*8, intent(in) :: x1,x2
    real*8, intent(in) :: maxerror
    integer, intent(in) :: numsteps
    integer, intent(in), optional:: recurdepth
    integer :: precurdepth ! passed to the next function
    integer :: i
    real*8 :: stepsize,sum,sumloworder,loworderintegral,xdash,error
    real*8 :: integralresult

    if (present(recurdepth)) then
        precurdepth = recurdepth
    else
        precurdepth = 0
    end if

    stepsize = (x2 - x1) / numsteps

    ! add up the boundary terms and the first midpoint
    if (present(intarg) .and. present(realarg)) then
        sumloworder = f(x1,realarg,intarg) + f(x2,realarg,intarg)
        sum = sumloworder + 4 * f(x1 + stepsize * 0.5d0,realarg,intarg)
    else if (present(realarg)) then
        sumloworder = f(x1,realarg) + f(x2,realarg)
        sum = sumloworder + 4 * f(x1 + stepsize * 0.5d0,realarg)
    else
        sumloworder = f(x1) + f(x2)
        sum = sumloworder + 4 * f(x1 + stepsize * 0.5d0)
    end if

    !!$omp parallel do private (i,xdash) reduction (+: sum,sumloworder)
    do i = 1,numsteps-1
        xdash = x1 + stepsize * i
        if (present(intarg) .and. present(realarg)) then
            sum = sum + 2 * f(xdash,realarg,intarg) + 4 * f(xdash + stepsize * 0.5d0,realarg,intarg)
            sumloworder = sumloworder + 2 * f(xdash,realarg,intarg)
        else if (present(realarg)) then
            sum = sum + 2 * f(xdash,realarg) + 4 * f(xdash + stepsize * 0.5d0,realarg)
            sumloworder = sumloworder + 2 * f(xdash,realarg)
        else
            sum = sum + 2 * f(xdash) + 4 * f(xdash + stepsize * 0.5d0)
            sumloworder = sumloworder + 2 * f(xdash)
        end if
    end do
    !!$omp end parallel do

    integralresult = sum * stepsize * (1.d0/6.d0)
    loworderintegral = sumloworder * stepsize * 0.5d0
    error = abs(loworderintegral/integralresult - 1.0d0)
    if (error > maxerror .and. maxerror > 1e-14) then
        if (precurdepth >= 19) then
!            if (present(label)) then
!                write(14,'(A,A,I8,A,ES10.3,A,ES10.3,A,ES10.3)') label,' reached maximum depth of ',precurdepth,', ',&
!                    integralresult,' error ',error,' >',maxerror
!            else
!                write(14,'(A,I8,A,ES10.3,A,ES10.3,A,ES10.3)') 'Integrator recursing to ',precurdepth,', ',&
!                    integralresult,' error',error,' >',maxerror
!            end if
        else
        !loworderintegral = integralresult
!        if (present(label)) then
!            write(14,'(A,A,I8,A,ES10.3,A,ES10.3,A,ES10.3)') label,' recursing to ',numsteps*2**(precurdepth+1),' steps ',&
!                integralresult,' error ',error,' >',maxerror
!        else
!            write(14,'(A,I8,A,ES10.3,A,ES10.3,A,ES10.3)') 'Integrator recursing to ',numsteps*2**(precurdepth+1),' steps ',&
!                integralresult,' error',error,' >',maxerror
!        end if
            if (present(intarg) .and. present(realarg)) then
                integralresult = integral(f,x1,x1+(x2-x1)*0.5,numsteps,maxerror,intarg=intarg,&
                    realarg=realarg,label=label,recurdepth=precurdepth+1) +&
                    integral(f,(x2+x1)*0.5,x2,numsteps,maxerror,intarg=intarg,&
                    realarg=realarg,label=label,recurdepth=precurdepth+1)
            else if (present(intarg)) then
                integralresult = integral(f,x1,x1+(x2-x1)*0.5,numsteps,maxerror,intarg=intarg,&
                    label=label,recurdepth=precurdepth+1) +&
                    integral(f,(x2+x1)*0.5,x2,numsteps,maxerror,intarg=intarg,&
                    label=label,recurdepth=precurdepth+1)
            else if (present(realarg)) then
                integralresult = integral(f,x1,x1+(x2-x1)*0.5,numsteps,maxerror,realarg=realarg,&
                    label=label,recurdepth=precurdepth+1) +&
                    integral(f,(x2+x1)*0.5,x2,numsteps,maxerror,realarg=realarg,&
                    label=label,recurdepth=precurdepth+1)
            else
                integralresult = integral(f,x1,x1+(x2-x1)*0.5,numsteps,maxerror,&
                        label=label,recurdepth=precurdepth+1) +&
                    integral(f,(x2+x1)*0.5,x2,numsteps,maxerror,&
                    label=label,recurdepth=precurdepth+1)
            end if
        end if
    end if
end function integral

! sort stellarmodel array entries by initial mass, required for interpolation to work
subroutine sortstellarmodels()
    integer :: i,j,minkeyposition
    real*8 :: minkeyvalue
    type(stellarmodeldata) :: temp

    do i = 1,size(stellarmodel)-1
        minkeyposition = i
        minkeyvalue = stellarmodel(i)%mass
        do j = i+1,size(stellarmodel)
            if (stellarmodel(i)%mass == stellarmodel(j)%mass) then
                write(14,*),"ERROR: multiple ev. models with same mass",stellarmodel(i)%mass
                stop
            end if
            if (stellarmodel(j)%mass < minkeyvalue) then
                minkeyposition = j
                minkeyvalue = stellarmodel(j)%mass
            end if
        end do
        if (minkeyvalue < stellarmodel(i)%mass) then
            temp = stellarmodel(i)
            stellarmodel(i) = stellarmodel(minkeyposition)
            stellarmodel(minkeyposition) = temp
        end if
    end do
end subroutine sortstellarmodels

subroutine initconfig()
    integer :: ios
    character(len=45) :: configline

    write(14,'(A)') 'reading config file...'
    open(unit=114, file=trim(filespath) // "config.txt", action="read", status="old")

    do while (.true.)
        read(114,'(A45)',iostat=ios) configline
        if (IS_IOSTAT_END(ios)) exit

        if (configline(1:25) == 'initialgasmass           ') then
            read(configline,'(25X,D20.7)') initialgasmass
            write(14,'(A,ES13.6)') 'initialgasmass =     ',initialgasmass

        else if (configline(1:25) == 'minstellarmass           ') then
            read(configline,'(25X,D20.7)') minstellarmass
            write(14,'(A,F13.2)') 'minstellarmass =       ',minstellarmass

        else if (configline(1:25) == 'maxstellarmass           ') then
            read(configline,'(25X,D20.7)') maxstellarmass
            write(14,'(A,F13.2)') 'minstellarmass =       ',maxstellarmass

        else if (configline(1:25) == 'sfrefficiency            ') then
            read(configline,'(25X,D20.7)') sfrefficiency
            write(14,'(A,ES13.6)') 'sfrefficiency =      ',sfrefficiency

        else if (configline(1:25) == 'sfrstartmgasabove        ') then
            read(configline,'(25X,D20.7)') sfrstartmgasabove
            write(14,'(A,ES13.6)') 'sfrstartmgasabove =  ',sfrstartmgasabove

        else if (configline(1:25) == 'sfrendmgasbelow          ') then
            read(configline,'(25X,D20.7)') sfrendmgasbelow
            write(14,'(A,ES13.6)') 'sfrendmgasbelow =    ',sfrendmgasbelow

        else if (configline(1:25) == 'maxtime                  ') then
            read(configline,'(25X,D20.7)') maxtime
            write(14,'(A,ES13.6)') 'maxtime =            ',maxtime

        else if (configline(1:25) == 'initialtimestep          ') then
            read(configline,'(25X,D20.7)') initialtimestep
            write(14,'(A,ES13.6)') 'initialtimestep =    ',initialtimestep

        else if (configline(1:25) == 'mintimestep              ') then
            read(configline,'(25X,D20.7)') mintimestep
            write(14,'(A,ES13.6)') 'mintimestep =        ',mintimestep

        else if (configline(1:25) == 'maxtimestep              ') then
            read(configline,'(25X,D20.7)') maxtimestep
            write(14,'(A,ES13.6)') 'maxtimestep =        ',maxtimestep

        else if (configline(1:25) == 'timestepgrowthrate       ') then
            read(configline,'(25X,D20.7)') timestepgrowthrate
            write(14,'(A,ES13.6)') 'timestepgrowthrate = ',timestepgrowthrate

        else if (configline(1:25) == 'maxmodelnum              ') then
            read(configline,'(25X,I20)') maxmodelnum
            write(14,'(A,I13)') 'maxmodelnum =            ',maxmodelnum

        else if (configline(1:25) == 'maxchangemgas            ') then
            read(configline,'(25X,D20.7)') maxchangemgas
            write(14,'(A,ES13.6)') 'maxchangemgas =      ',maxchangemgas

        else if (configline(1:25) == 'maxchangeser             ') then
            read(configline,'(25X,D20.7)') maxchangeser
            write(14,'(A,ES13.6)') 'maxchangeser  =      ',maxchangeser

        else if (configline(1:25) == 'serintegralsteps         ') then
            read(configline,'(25X,I20)') serintegralsteps
            write(14,'(A,I13)') 'serintegralsteps =   ',serintegralsteps

        else if (configline(1:25) == 'initialcompscalez        ') then
            read(configline,'(25X,D20.7)') initialcompscalez
            write(14,'(A,ES13.6)') 'initialcompscalez =  ',initialcompscalez

        end if
    end do

    if (mintimestep > maxtimestep) then
        write(14,'(A)') 'STOPPING. mintimestep > maxtimestep'
        stop
    end if
    if (sfrendmgasbelow > sfrstartmgasabove) then
        write(14,'(A)') 'STOPPING. sfrendmgasbelow > sfrstartmgasabove'
        stop
    end if

    close(114)
    write(14,'(A)') 'finished reading config file.'
    write(14,*)
end subroutine initconfig

subroutine initspecies()
    integer :: s, speciescount, neutrons

    write(14,'(A)') 'reading species.dat...'
    open(unit=114, file=trim(filespath) // "species.dat", action="read", status="old")
    read(114,'(I16)') speciescount
    allocate(species(1:speciescount))
    write(14,'(I4,A)') speciescount, ' species'
    write(14,'(A6,A4,A4,A4)') 'name','el','Z','A'
    do s=1,size(species)
        read(114,'(I6,A6,A6,I6)') species(s)%a, species(s)%symbol, species(s)%name, neutrons
        species(s)%z = species(s)%a - neutrons
        write(14,'(A6,A4,I4,I4)') species(s)%name,species(s)%symbol,species(s)%z,species(s)%a
    end do
    close(114)
    write(14,'(A)') 'finished reading species.dat.'
    write(14,*)
end subroutine initspecies

subroutine inityields()
    integer :: stellarmodelcount
    integer :: ios, i, s
    character(len=25) :: stellarmodelname
    character(len=6) :: spname !species name
    character(len=8) :: startyieldlist !line in yield file
    character(len=15) :: startmodellist !line in yield file
    character(len=10) :: absrel !absolute or relative yield
    real*8,allocatable :: yieldrow(:)

    write(14,'(A)') 'reading yields.dat...'
    open(unit=7, file=trim(filespath) // 'yields.txt', action="read", status="old", access="sequential", form="formatted")
    read(7,*,iostat=ios)
    startmodellist = ""
    do while (startmodellist /= "[stellarmodels]")
        read(7,*) startmodellist
    end do
    read(7,*) stellarmodelcount
    allocate(stellarmodel(stellarmodelcount))
    write(14,'(I4,A)') size(stellarmodel), ' stellar models'

    write(14,'(A25,1X,A6,A11,A14)') 'ModelName','Mini','Mremnant','lifetime'
    do i = 1, size(stellarmodel)
        read(7,'(A25,E14.2,14X,E14.2,E14.2)') stellarmodelname,stellarmodel(i)%mass,stellarmodel(i)%remnantmass,&
            stellarmodel(i)%stellarlifetime
        write(14,'(A25,1X,F6.2,F11.3,ES14.3)') stellarmodelname,stellarmodel(i)%mass,stellarmodel(i)%remnantmass,&
            stellarmodel(i)%stellarlifetime
        allocate(stellarmodel(i)%yield(size(species)))
        stellarmodel(i)%yield(:) = 0.0d0
    end do

    startyieldlist = ""
    do while (startyieldlist /= "[yields]")
        read(7,*) startyieldlist
    end do

    allocate(yieldrow(stellarmodelcount))
    do while (.true.)
        read(7,'(A6,2X,A10,*(E14.6))',iostat=ios) spname, absrel, yieldrow
        if (IS_IOSTAT_END(ios)) exit
!        write(14,*) name,absrel,yieldrow

        do s = 1, size(species)
            if (trim(adjustl(species(s)%name))==trim(adjustl(spname))) then
                write(14,'(A,A)') 'loading yields of ',species(s)%name
                do i = 1, size(stellarmodel)
                    if (absrel=='  relative') then
                        species(s)%yieldisrelative = .true.
                    else
                        species(s)%yieldisrelative = .false.
                    end if
                        species(s)%foundyields = .true.
                    stellarmodel(i)%yield(s) = yieldrow(i)
                    !write(14,'(A,F5.2,A,E14.7)') 'M=',stellarmodel(i)%mass,', yield=',stellarmodel(i)%yield(s)
                end do
            end if
        end do
    end do
    deallocate(yieldrow)
    close(7)

    do s = 1, size(species)
        if (species(s)%foundyields .eqv. .false.) then
            write(14,'(A,A,A)') 'no yields found for ',species(s)%name,', setting to 0.0'
            ! have already been initialised to zero earlier on
        end if
    end do

    call sortstellarmodels()

    ! get the exponents for lifetime = A * mass ^ B, based on each model and the following (higher mass) one
    do i = 1, size(stellarmodel)-1
        stellarmodel(i)%lifetimemassexponent = log(stellarmodel(i+1)%stellarlifetime / stellarmodel(i)%stellarlifetime) / &
            log(stellarmodel(i+1)%mass / stellarmodel(i)%mass)
    end do
    stellarmodel(size(stellarmodel))%lifetimemassexponent = stellarmodel(size(stellarmodel)-1)%lifetimemassexponent

    ! debugging
!    do i = 1, 100
!        m0 = 1.0d0 + i*0.5d0
!        write(14,*) i,stellarmodel(i)%mass,stellarmodel(i)%stellarlifetime,&
!            stellarmodel(i)%remnantmass,stellarmodel(i)%yield(23)
!        write(14,'(F6.2,F7.3,*(ES12.3))') m0, remnantmass(m0), stellarlifetime(m0)!, yield(23,stellarmodel(i)%mass)
!    end do
!    stop

    write(14,'(A)') 'finished reading yields.dat.'
    write(14,*)

end subroutine inityields

subroutine loadinitialcomposition ()
    integer :: ios, s, inmassnum
    character(len=6) :: inspname
    real*8 :: inabund

    write(14,'(A)') 'reading initial_comp.dat...'
    open(unit=115, file=trim(filespath) // "initial_comp.dat", action="read", status="old")
    do while (.true.)
        read(115,'(1X,A6,E14.4,I7)',iostat=ios) inspname, inabund, inmassnum
        if (IS_IOSTAT_END(ios)) exit
        !write(14,'(A6,E14.4,I7)') inspname, inabund, inmassnum
        do s = 1, size(species)
            if (trim(adjustl(species(s)%name))==trim(adjustl(inspname))) then
                species(s)%foundinitialvalue = .true.
                if (species(s)%z > 2) then
                    model(1)%speciesmassfrac(s) = inabund * float(inmassnum) * initialcompscalez
                else
                    model(1)%speciesmassfrac(s) = inabund * float(inmassnum)
                end if
                write(14,'(A,A6,A,ES14.7)') 'setting initial mass fraction of ',species(s)%name,' to ',model(1)%speciesmassfrac(s)
            end if
        end do
    end do
    close(115)

    do s = 1, size(species)
        if (species(s)%foundinitialvalue .eqv. .false.) then
            write(14,'(A,A,A)') 'no initial mass frac for ',species(s)%name,', setting to 0.0'
        end if
    end do

    write(14,'(A)') 'finished reading initial_comp.dat.'
    write(14,*)
end subroutine loadinitialcomposition

! returns star formation rate in solar masses at time in years
pure real*8 function starformationrate(time)
    real*8, intent(in) :: time
    integer :: modelnumber
    real*8 :: c

    if (time < 0.d0) then
        starformationrate = 0.d0
    else
        modelnumber = findmodelnum(time)
        if (modelnumber >= stepnum - 1) then
            starformationrate = model(modelnumber-1)%sfr !current model doesn't have an SFR yet
        elseif (modelnumber >= 2) then
            starformationrate = model(modelnumber)%sfr
            c = (time - model(modelnumber)%time) / &
                (model(modelnumber+1)%time - model(modelnumber)%time)
            ! this is necessary to conserve mass. The starformation rate was euler integrated
            ! to get the change in stellar/gass mass and the area must be kept constant
            ! under the interpolation
            starformationrate = 0.5d0 * (c * (model(modelnumber)%sfr + model(modelnumber+1)%sfr) +&
                (1-c) * (model(modelnumber-1)%sfr + model(modelnumber)%sfr))
        else
            starformationrate = 0.0d0
        end if
    end if

    !    starformationrate = 0.4d1  * exp(-time/3.2d4)
end function starformationrate

! initial mass function by number (dN/dM)
! normalised s.t. integral of m*imf(m) from m=minstellarmass to m=maxstellarmass Msun is 1.0
! i.e. the normalised function gives dN/dM per Msun of star formation
pure real*8 function imf(initmass)
    implicit none
    real*8, intent(in) :: initmass

    ! Modified Kroupa, Tout, Gilmore 1993 IMF
    !if (mass >= minstellarmass .AND. mass <= maxstellarmass) then
    if (initmass > 1.0d0) then
        imf = imfnormfactor * initmass ** (-2.7d0)
    elseif (initmass > 0.5d0) then
        imf = imfnormfactor * initmass ** (-2.2d0)
    elseif (initmass > 0.08d0) then
        imf = imfnormfactor * initmass ** (-1.50d0) / (0.5 ** (0.7d0))
    else
        imf = 0.0d0
    end if

    ! Kroupa 2001 IMF
    !if (mass >= minstellarmass .AND. mass <= maxstellarmass) then
    !    if (mass < 0.08d0) then
    !        imf = imfnormfactor * mass ** (-0.3d0)
    !    elseif (mass < 0.5d0) then
    !        imf = imfnormfactor * 0.08d0 * mass ** (-1.3d0)
    !    else
    !        imf = imfnormfactor * 0.08d0 * 0.5d0 * mass ** (-2.3d0)
    !    end if
    !else
    !    imf = 0.0d0
    !end if
end function imf

! dM/dM, the mass of stars from M to M + dM,
! normalised to a total stellar mass of 1 Msun
pure real*8 function mimf(initmass, dummyreal, dummyint)
    implicit none
    real*8, intent(in) :: initmass
    real*8, intent(in), optional :: dummyreal
    integer, intent(in), optional :: dummyint

    mimf = initmass * imf(initmass)
end function mimf

subroutine closefiles()
    flush(12)
    close(12)
    flush(13)
    close(13)
    flush(14)
    close(14)
end subroutine

end module
