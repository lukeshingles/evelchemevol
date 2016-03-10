module cemodel
    use cemodelsetup
    use stringutility

    implicit none

    integer(kind=4):: stepnum
    real(kind=8) :: timestep

    type evmodelatonemassonemetallicitydata
        real(kind=8) :: fetoh
        real(kind=8) :: remnantmass
        real(kind=8), allocatable :: yield(:)
        real(kind=8), allocatable :: speciesmasslost(:)
    end type

    type evmodeldataonemass
        type(evmodelatonemassonemetallicitydata), allocatable :: Z(:)
        real(kind=8) :: mass
    end type

    type(evmodeldataonemass), allocatable :: evmodel(:)

    type speciestype
        character (len=3) :: symbol ! He, C, Fe, etc
        integer :: z ! number of protons
        integer :: a ! number of nucleons
    end type

    type(speciestype), allocatable :: species(:)

    type modelstate
        real(kind=8) :: time, gasmass, starsmass, fetoh
        real(kind=8), allocatable :: speciesmassfrac(:)
    end type

    type(modelstate), allocatable :: model(:)

contains

real(kind=8) function yieldatfetoh(species, modelnum, fetoh)
    integer :: modelnum, i, Zlownum, Zhighnum
    integer, intent(in) :: species
    real(kind=8) :: fetoh, weight2

    Zlownum = 1
    Zhighnum = 1
    do i = 1, size(evmodel(modelnum)%Z)
        if (evmodel(modelnum)%Z(i)%fetoh > evmodel(modelnum)%Z(Zlownum)%fetoh .and. evmodel(modelnum)%Z(i)%fetoh < fetoh) then
            Zlownum = i
        end if
        if (evmodel(modelnum)%Z(i)%fetoh < evmodel(modelnum)%Z(Zhighnum)%fetoh .and. evmodel(modelnum)%Z(i)%fetoh > fetoh) then
            Zhighnum = i
        end if
    end do

    if (Zlownum /= Zhighnum) then
        weight2 = (fetoh - evmodel(modelnum)%Z(Zlownum)%fetoh) /&
                (evmodel(modelnum)%Z(Zhighnum)%fetoh - evmodel(modelnum)%Z(Zlownum)%fetoh)
    else
        weight2 = 0.0d0
    end if

    yieldatfetoh = (1-weight2) * evmodel(modelnum)%Z(Zlownum)%yield(species) +&
            weight2 * evmodel(modelnum)%Z(Zhighnum)%yield(species)
end function

real(kind=8) function speciesmasslostatfetoh(species, modelnum, fetoh)
    integer :: modelnum, i, Zlownum, Zhighnum
    integer, intent(in) :: species
    real(kind=8) :: fetoh, weight2

    Zlownum = 1
    Zhighnum = 1
    do i = 1, size(evmodel(modelnum)%Z)
        if (evmodel(modelnum)%Z(i)%fetoh > evmodel(modelnum)%Z(Zlownum)%fetoh .and. evmodel(modelnum)%Z(i)%fetoh < fetoh) then
            Zlownum = i
        end if
        if (evmodel(modelnum)%Z(i)%fetoh < evmodel(modelnum)%Z(Zhighnum)%fetoh .and. evmodel(modelnum)%Z(i)%fetoh > fetoh) then
            Zhighnum = i
        end if
    end do

    if (Zlownum /= Zhighnum) then
        weight2 = (fetoh - evmodel(modelnum)%Z(Zlownum)%fetoh) /&
                (evmodel(modelnum)%Z(Zhighnum)%fetoh - evmodel(modelnum)%Z(Zlownum)%fetoh)
    else
        weight2 = 0.0d0
    end if

    speciesmasslostatfetoh = (1-weight2) * evmodel(modelnum)%Z(Zlownum)%speciesmasslost(species) +&
            weight2 * evmodel(modelnum)%Z(Zhighnum)%speciesmasslost(species)
end function

real(kind=8) function pastfetoh(time)
    real(kind=8), intent(in) :: time
    integer :: i, j
    real(kind=8) :: weight2

    ! replace with binary search for better performance
    do i = 2,stepnum
        if (model(i)%time > time) then
            exit
        end if
    end do
    pastfetoh = model(i-1)%fetoh
end function

real(kind=8) function pastspeciesmassfrac(time,species)
    real(kind=8), intent(in) :: time
    integer, intent(in) :: species
    integer :: i
    real(kind=8) :: weight2
    ! replace with binary search for better performance
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
        remnantmass = evmodel(1)%Z(1)%remnantmass * initmass / evmodel(1)%mass ! scale lowest mass end with initial mass
    else if (initmass > evmodel(size(evmodel))%mass) then
        remnantmass = evmodel(size(evmodel))%Z(1)%remnantmass
    else
    ! interpolate remnant mass from ev. model data
        do i = 2, size(evmodel)
            if (evmodel(i)%mass >= initmass) then
                c = (initmass - evmodel(i-1)%mass) / (evmodel(i)%mass - evmodel(i-1)%mass)
                remnantmass = c * evmodel(i)%Z(1)%remnantmass + (1-c) * evmodel(i-1)%Z(1)%remnantmass
                exit
            end if
        end do
    end if
end function

! absolute yield of species in solar masses from a star with
! initial mass in solar masses
real(kind=8) function speciesmasslost(species, initmass, fetoh)
    real(kind=8), intent(in) :: initmass, fetoh
    integer, intent(in) :: species
    integer :: i
    real(kind=8) :: c

    speciesmasslost = 0.0

    if (initmass < evmodel(1)%mass) then
        ! scale the lowest mass yield with the ejecta mass
        speciesmasslost = speciesmasslostatfetoh(species,1,fetoh) * &
                (initmass - remnantmass(initmass)) / (evmodel(1)%mass - evmodel(1)%Z(1)%remnantmass)
    else if (initmass > evmodel(size(evmodel))%mass) then
        ! scale the highest mass yield with the ejecta mass
        speciesmasslost = speciesmasslostatfetoh(species,size(evmodel),fetoh) * (initmass - remnantmass(initmass)) / &
                (evmodel(size(evmodel))%mass - evmodel(size(evmodel))%Z(1)%remnantmass)
    else
        do i = 2, size(evmodel)
            if (evmodel(i)%mass >= initmass) then
                c = (initmass - evmodel(i-1)%mass) / (evmodel(i)%mass - evmodel(i-1)%mass)
                speciesmasslost = (1-c) * speciesmasslostatfetoh(species,i-1,fetoh) +&
                        c * speciesmasslostatfetoh(species,i,fetoh)
                exit
            end if
        end do
    end if
end function

! relative yield of species in solar masses from a star with
! initial mass in solar masses
real(kind=8) function yield(species, initmass, fetoh)
    real(kind=8), intent(in) :: initmass, fetoh
    integer, intent(in) :: species
    integer :: i
    real(kind=8) :: c

    yield = 0.0

    if (initmass < evmodel(1)%mass) then
        ! scale the lowest mass yield with the ejecta mass
        yield = yieldatfetoh(species,1,fetoh) * (initmass - remnantmass(initmass)) / (evmodel(1)%mass - evmodel(1)%Z(1)%remnantmass)
    else if (initmass > evmodel(size(evmodel))%mass) then
        ! scale the highest mass yield with the ejecta mass
        yield = yieldatfetoh(species,size(evmodel),fetoh) * (initmass - remnantmass(initmass)) / &
                (evmodel(size(evmodel))%mass - evmodel(size(evmodel))%Z(1)%remnantmass)
    else
        do i = 2, size(evmodel)
            if (evmodel(i)%mass >= initmass) then
                c = (initmass - evmodel(i-1)%mass) / (evmodel(i)%mass - evmodel(i-1)%mass)
                yield = (1-c) * yieldatfetoh(species,i-1,fetoh) + c * yieldatfetoh(species,i,fetoh)
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
    real(kind=8) :: timeatbirth, fetohatbirth

    timeatbirth = time - stellarlifetime(mdash)
    fetohatbirth = pastfetoh(timeatbirth)
    if (timeatbirth >= 0.0d0) then
        if (species(speciesnum)%z > 26) then
        ! use relative yields
            dserbydmassofspecies = ((mdash - remnantmass(mdash)) * &
                pastspeciesmassfrac(timeatbirth,speciesnum) + &
                yield(speciesnum, mdash, fetohatbirth)) * &
                (starformationrate(timeatbirth) / mdash) * imf(mdash)
            if (dserbydmassofspecies < 0.0d0) then
                !write(*,*) 'WARNING: Negative ejection rate, species:',speciesnum,', mass:',mdash
                !dserbydmassofspecies = 0.0d0
            end if
        else
            !use absolute yields
            dserbydmassofspecies = speciesmasslost(speciesnum, mdash, fetohatbirth) * &
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
    type(evmodeldataonemass) :: temp

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
end module
