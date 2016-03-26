module cemodel
    use cemodelsetup
    use stringutility

    implicit none

    integer(kind=4):: stepnum
    real(kind=8) :: timestep

    type evmodeldata
        real(kind=8) :: mass
        real(kind=8) :: remnantmass
        real(kind=8), allocatable :: yield(:)
    end type

    type(evmodeldata), allocatable :: evmodel(:)

    type speciestype
        character (len=3) :: symbol ! He, C, Fe, etc
        integer :: z ! number of protons
        integer :: a ! number of nucleons
        logical :: yieldisrelative
    end type

    type(speciestype), allocatable :: species(:)

    type modelstate
        real(kind=8) ::  time, gasmass, starsmass
        real(kind=8), allocatable :: speciesmassfrac(:)
    end type

    type(modelstate), allocatable :: model(:)

contains

real(kind=8) function pastspeciesmassfrac(time,species)
    real(kind=8), intent(in) :: time
    integer, intent(in) :: species
    integer :: i
    real(kind=8) :: weight2

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

!here
real(kind=8) function logepsilon(symbol)
    character (len=*), intent(in) :: symbol
    real(kind=8) :: eldensitysum
    integer :: s, hydrogenindex

    hydrogenindex = 1
    eldensitysum = 0.0d0
    do s = 1,size(species)
        if (species(s)%z == 1 .and. species(s)%a == 1) then
            hydrogenindex = s
        end if

        if (trim(adjustl(species(s)%symbol)) == symbol) then
            eldensitysum = eldensitysum + (model(stepnum)%speciesmassfrac(s) / species(s)%a)
        end if
    end do
    logepsilon = log10(eldensitysum / model(stepnum)%speciesmassfrac(hydrogenindex)) + 12.0d0
end function

! remnant mass as a function of initial mass in solar masses
real(kind=8) function remnantmass(initmass)
    implicit none
    real(kind=8), intent(in) :: initmass
    integer :: i
    real(kind=8) :: c

!   Iben & Tutukov 1984, in Pagel 2009 after eqn 7.10
!    if (mass <= 0.506) then
!        remnantmass = mass
!    else if (mass <= 9.5) then
!        remnantmass = 0.45d0 + 0.11d0 * mass
!    else
!        remnantmass = 1.5d0
!    end if

    if (initmass < evmodel(1)%mass) then
        remnantmass = evmodel(1)%remnantmass * initmass / evmodel(1)%mass ! scale lowest mass end with initial mass
    else if (initmass > evmodel(size(evmodel))%mass) then
        remnantmass = evmodel(size(evmodel))%remnantmass
    else
    ! interpolate remnant mass from ev. model data
        do i = 2, size(evmodel)
            if (evmodel(i)%mass >= initmass) then
                c = (initmass - evmodel(i-1)%mass) / (evmodel(i)%mass - evmodel(i-1)%mass)
                remnantmass = c * evmodel(i)%remnantmass + (1-c) * evmodel(i-1)%remnantmass
                exit
            end if
        end do
    end if
end function

! relative or absolute yield of species in solar masses from a star with
! initial mass in solar masses
real(kind=8) function yield(species, initmass)
    real(kind=8), intent(in) :: initmass
    integer, intent(in) :: species
    integer :: i
    real(kind=8) :: c

    yield = 0.0

    if (initmass < evmodel(1)%mass) then
        ! scale the lowest mass yield with the ejecta mass
        yield = evmodel(1)%yield(species) * (initmass - remnantmass(initmass)) / (evmodel(1)%mass - evmodel(1)%remnantmass)
    else if (initmass > evmodel(size(evmodel))%mass) then
        ! scale the highest mass yield with the ejecta mass
        yield = evmodel(size(evmodel))%yield(species) * (initmass - remnantmass(initmass)) / &
                (evmodel(size(evmodel))%mass - evmodel(size(evmodel))%remnantmass)
    else
        do i = 2, size(evmodel)
            if (evmodel(i)%mass >= initmass) then
                c = (initmass - evmodel(i-1)%mass) / (evmodel(i)%mass - evmodel(i-1)%mass)
                yield = (1-c) * evmodel(i-1)%yield(species) + c * evmodel(i)%yield(species)
                exit
            end if
        end do
    end if
end function

! stellar ejection rate in solar masses per year
real(kind=8) function stellarejectionrate(time)
    real(kind=8), intent(in) :: time

    stellarejectionrate = integral(dserbydmass, stellarmasslow, stellarmasshigh, 10 ** 4, sermaxerror, &
        realarg=time, label='ser mass integral')
end function

! D stellar ejection rate / D mass
real(kind=8) function dserbydmass(mdash, time)
    real(kind=8), intent(in) :: mdash, time
    real(kind=8) :: timeatbirth

    timeatbirth = time - stellarlifetime(mdash)

    if (timeatbirth >= 0.0d0) then
        dserbydmass = (mdash - remnantmass(mdash)) * &
                    (starformationrate(timeatbirth) / mdash) * imf(mdash)
    else
        dserbydmass = 0.0d0
    end if
end function

! stellar ejection rate of species speciesnum in solar masses per year
real(kind=8) function stellarejectionrateofspecies(time, speciesnum)
    real(kind=8), intent(in) :: time
    integer, intent(in) :: speciesnum

    stellarejectionrateofspecies = integral(dserbydmassofspecies, stellarmasslow, stellarmasshigh, 10 ** 4,&
        sermaxerror, intarg=speciesnum, realarg=time, label='ser species mass integral')
end function

! D stellar ejection rate / D mass of species 
real(kind=8) function dserbydmassofspecies(mdash, speciesnum, time)
    real(kind=8), intent(in) :: mdash, time
    integer, intent(in) :: speciesnum
    real(kind=8) :: timeatbirth

    timeatbirth = time - stellarlifetime(mdash)
    if (timeatbirth >= 0.0d0) then
        if (species(speciesnum)%yieldisrelative .eqv. .TRUE.) then
        ! use relative yields
            dserbydmassofspecies = ((mdash - remnantmass(mdash)) * &
                pastspeciesmassfrac(timeatbirth,speciesnum) + &
                yield(speciesnum, mdash)) * &
                (starformationrate(timeatbirth) / mdash) * imf(mdash)
            if (dserbydmassofspecies < 0.0d0) then
                !write(*,*) 'WARNING: Negative ejection rate, species:',speciesnum,', mass:',mdash
                !dserbydmassofspecies = 0.0d0
            end if
        else
            !use absolute yields
            dserbydmassofspecies = yield(speciesnum, mdash) * &
                (starformationrate(timeatbirth) / mdash) * imf(mdash)
        end if
    else
        dserbydmassofspecies = 0.0d0
    end if
end function

!integrate with Simpson's rule (and mid-point to measure error)
real(kind=8) function integral(f,x1,x2,numstepsguess,maxerror,errorout,intarg,realarg,label)
    real(kind=8), external :: f
    integer, intent(in), optional :: intarg
    character(*), intent(in), optional :: label
    real(kind=8), intent(in), optional :: realarg
    real(kind=8), intent(in) :: x1,x2
    real(kind=8), intent(in) :: maxerror
    real(kind=8), intent(out), optional :: errorout
    integer(kind=4), intent(in):: numstepsguess
    integer(kind=4) :: i, numsteps
    real(kind=8) :: stepsize,sum,sumloworder,loworderintegral,xdash,error

    numsteps = min(numstepsguess,10 ** 8)
    error = maxerror

    do while ((error > maxerror .or. numsteps == numstepsguess) .or. numsteps > 10 ** 8)
        stepsize = (x2-x1) / numsteps

        if (present(intarg) .and. present(realarg)) then
            sumloworder = f(x1,intarg,realarg) + f(x2,intarg,realarg) ! just the boundaries
            sum = sumloworder + 4 * f(x1 + stepsize * 0.5d0,intarg,realarg)
        else if (present(intarg)) then
            sumloworder = f(x1,intarg) + f(x2,intarg) ! just the boundaries
            sum = sumloworder + 4 * f(x1 + stepsize * 0.5d0,intarg) ! first midpoint
        else if (present(realarg)) then
            sumloworder = f(x1,realarg) + f(x2,realarg) ! just the boundaries
            sum = sumloworder + 4 * f(x1 + stepsize * 0.5d0,realarg) ! first midpoint
        else
            sumloworder = f(x1) + f(x2) ! just the boundaries
            sum = sumloworder + 4 * f(x1 + stepsize * 0.5d0) ! first midpoint
        end if

    !$omp parallel do private (i,xdash) reduction (+: sum,sumloworder)
        do i = 1,numsteps-1
            xdash = x1 + stepsize * i
            if (present(intarg) .and. present(realarg)) then
                sum = sum + 2 * f(xdash,intarg,realarg) + 4 * f(xdash + stepsize * 0.5d0,intarg,realarg)
                sumloworder = sumloworder + 2 * f(xdash,intarg,realarg)
            else if (present(intarg)) then
                sum = sum + 2 * f(xdash,intarg) + 4 * f(xdash + stepsize * 0.5d0,intarg)
                sumloworder = sumloworder + 2 * f(xdash,intarg)
            else if (present(realarg)) then
                sum = sum + 2 * f(xdash,realarg) + 4 * f(xdash + stepsize * 0.5d0,realarg)
                sumloworder = sumloworder + 2 * f(xdash,realarg)
            else
                sum = sum + 2 * f(xdash) + 4 * f(xdash + stepsize * 0.5d0)
                sumloworder = sumloworder + 2 * f(xdash)
            end if
        end do
    !$omp end parallel do
    
        integral = sum * stepsize * (1.d0/6.d0)
        loworderintegral = sumloworder * stepsize * 0.5d0
        error = abs(loworderintegral - integral)
        numsteps = numsteps * 2
        if (integral == 0 .and. loworderintegral==0) then
            exit
        end if
    end do
    if (numsteps /= numstepsguess * 2 .and. .false.) then
        if (present(label)) then
            write(*,*),label,' used',numsteps,'steps instead of',numstepsguess,integral,'+/-',error
        else
            write(*,*) 'Integrator used',numsteps,'steps instead of',numstepsguess,integral,'+/-',error
        end if
        if (present(realarg)) then
            write(*,*),realarg
        end if
    end if
    if (present(errorout)) then
        errorout = error
    end if
end function

! sort evmodel array entries by initial mass, required for interpolation to work
subroutine sortevmodels()
    integer :: i,j,minkeyposition
    real(kind=8) :: minkeyvalue
    type(evmodeldata) :: temp

    do i = 1,size(evmodel)-1
        minkeyvalue = evmodel(i)%mass
        do j = i+1,size(evmodel)
            if (evmodel(i)%mass == evmodel(j)%mass) then
                write(*,*),"ERROR: multiple ev. models with same mass",evmodel(i)%mass
                stop
            end if
            if (evmodel(j)%mass < minkeyvalue) then
                minkeyposition = j
                minkeyvalue = evmodel(j)%mass
            end if
        end do
        if (minkeyvalue < evmodel(i)%mass) then
            temp = evmodel(i)
            evmodel(i) = evmodel(minkeyposition)
            evmodel(minkeyposition) = temp
        end if
    end do
end subroutine

subroutine initspecies()
    integer :: s, speciescount, neutrons
    character(len=3) :: stra

    open(unit=14, file=trim(filespath) // "species.dat", action="read", status="old")
    read(14,'(I4)') speciescount
    allocate(species(1:speciescount))
    do s=1,size(species)
        read(14,'(I6,A6,A6,I6)') species(s)%a, species(s)%symbol, stra, neutrons
        species(s)%z = species(s)%a - neutrons
    end do
    flush(14)
end subroutine

subroutine initevmodels()
    integer :: evmodelcount
    integer :: ios, i, s
    character(len=14) :: evmodelname
    real(kind=8) :: m0, finalmass
    character(len=8) :: name,nameb !species names for comparison
    character(len=3) :: stra
    character(len=8) :: startyieldlist
    character(len=10) :: startmodellist
    character(len=10) :: absrel !absolute or relative yield
    real(kind=8),allocatable :: yieldrow(:)

    open(unit=7, file=trim(filespath) // 'yields.txt', action="read", status="old", access="sequential", form="formatted")
    read(7,*,iostat=ios)
    startmodellist = ""
    do while (startmodellist /= "[evmodels]")
        read(7,*) startmodellist
    end do
    read(7,*) evmodelcount
    allocate(evmodel(evmodelcount))

    do i = 1, size(evmodel)
        read(7,'(A25,E14.2,14X,E14.2)') evmodelname,m0,finalmass
!        write(*,*) evmodelname,m0,finalmass
        evmodel(i)%mass = m0
        evmodel(i)%remnantmass = finalmass

        allocate(evmodel(i)%yield(size(species)))
        evmodel(i)%yield = 0.0d0
    end do

    startyieldlist = ""
    do while (startyieldlist /= "[yields]")
        read(7,*) startyieldlist
    end do

    allocate(yieldrow(evmodelcount))
    do while (.true.)
        read(7,'(A8,A10,*(E14.6))',iostat=ios) name, absrel, yieldrow
        if (IS_IOSTAT_END(ios)) exit
!        write(*,*) name,absrel,yieldrow
!        write(*,*) absrel=='  relative'

        do s = 1, size(species)
            write(stra,'(I3)') species(s)%a
            if (species(s)%z == 1) then
                nameb = species(s)%symbol
            else
                nameb = trim(adjustl(StrLowCase(species(s)%symbol)))//trim(adjustl(stra))
            end if
            if (name==nameb) then
!                write(*,*)name,nameb,StrLowCase(name)==nameb
                do i = 1, size(evmodel)
                    if (absrel=='  relative') then
                        species(s)%yieldisrelative = .true.
                    else
                        species(s)%yieldisrelative = .false.
                    end if
                    evmodel(i)%yield(s) = yieldrow(i)
                end do
            end if
        end do
    end do

    close(7)

    call sortevmodels()

    !do i = 1, 4
    !    m0 = 10.0d0 + i*0.5
    !    write(*,*),evmodel(i)%mass,species(23)%z,evmodel(i)%yield(23),yield(23,evmodel(i)%mass)
    !end do
    !stop

end subroutine

end module
