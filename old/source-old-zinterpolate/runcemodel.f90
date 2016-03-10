program runcemodel
    use cemodelsetup
    use cemodel

    implicit none

    integer(kind=4) :: steplastoutput = 0
    real(kind=8) :: serdt, sfrdt, error
    real(kind=8) :: fetoh, rbtofe,srtofe,ytofe,zrtofe,motofe,batofe,latofe,cetofe,&
                    prtofe,ndtofe,smtofe,eutofe,hftofe,pbtofe
    real(kind=8), allocatable :: serspeciesdt(:)
    integer :: speciesnum
    logical :: acceptableerrors

    call initcemodel

    allocate(serspeciesdt(size(species)))

    open(unit=12, file="dataout/cemodel.txt", action="write", status="replace")
    open(unit=13, file="dataout/cemodel-abund.txt", action="write", status="replace")
23  format (E13.4','E10.3','E10.3','E10.3','E10.3','F6.2','E7.1','E10.3) !console output

    write(12,*) '#time  Mgas    Mstars  SFR SER'
    write(13,*) '#time    [Fe/H]  [Rb/Fe]   [Sr/Fe] [Y/Fe]  [Zr/Fe] [Mo/Fe] [Ba/Fe] [La/Fe]   [Ce/Fe] &
            [Pr/Fe] [Nd/Fe]   [Sm/Fe]   [Hf/Fe] [Pb/Fe]'

    timestep = 1.d2
    do stepnum = 2, size(model)
        acceptableerrors = .false.
        do while (acceptableerrors .eqv. .false.)
            acceptableerrors = .true.
            model(stepnum)%time = model(stepnum-1)%time + timestep

            ! all masses
            stellarmasslow = minstellarmass
            stellarmasshigh = maxstellarmass

            sfrdt = integral(starformationrate,model(stepnum-1)%time,model(stepnum)%time,1,1.d0,&
                    error,label='sfr time integral')
            if (error/timestep > sermaxerror .and. timestep > mintimestep) then
                write(*,*),'star formation rate error too high:',error/timestep
                goto 9000
            end if

            serdt = integral(stellarejectionrate,model(stepnum-1)%time,model(stepnum)%time,1,1.d0,error,&
                    label='ser-ret time integral')
            if (error/timestep > sermaxerror .and. timestep > mintimestep) then
                write(*,*),'stellar ejection rate retained error too high:',error/timestep
                goto 9000
            end if

            do speciesnum = 1, size(species)
                serspeciesdt(speciesnum) = integral(stellarejectionrateofspecies,model(stepnum-1)%time,&
                        model(stepnum)%time,1,1.d0,error,intarg=speciesnum,label='ser-ret species time integral')
                if (error/timestep > sermaxerror .and. timestep > mintimestep) then
                    write(*,*),'species-stellar ejection rate error too high:',speciesnum,error/timestep
                    goto 9000
                end if
            end do

            model(stepnum)%starsmass = model(stepnum-1)%starsmass + sfrdt - serdt
            model(stepnum)%gasmass = model(stepnum-1)%gasmass - sfrdt + serdt

            if (model(stepnum)%gasmass < 0.0d0) then
                write(*,*) 'Quitting: negative gas mass!'
                stop
            end if

            do speciesnum = 1, size(species)
                model(stepnum)%speciesmassfrac(speciesnum) = (1.0d0 / model(stepnum)%gasmass) * &
                    (model(stepnum-1)%speciesmassfrac(speciesnum) * model(stepnum-1)%gasmass + &
                    serspeciesdt(speciesnum) - model(stepnum-1)%speciesmassfrac(speciesnum) * sfrdt)
            end do

            timestep = min(timestep * 1.03,maxtimestep)
            exit
    9000    continue
                acceptableerrors = .false.
                timestep = max(timestep * 0.5d0, mintimestep)
                cycle
        end do

        if (stepnum - steplastoutput >= 1) then
            steplastoutput = stepnum

            ! Asplund et al. (ARAA, 2009)
            fetoh = logepsilon('fe') - 7.50d0
            rbtofe = logepsilon('rb') - logepsilon('fe') - (2.52d0 - 7.50d0)
            srtofe = logepsilon('sr') - logepsilon('fe') - (2.87d0 - 7.50d0)
            ytofe  = logepsilon('y')  - logepsilon('fe') - (2.21d0 - 7.50d0)
            zrtofe = logepsilon('zr') - logepsilon('fe') - (2.58d0 - 7.50d0)
            motofe = logepsilon('mo') - logepsilon('fe') - (1.88d0 - 7.50d0)
            batofe = logepsilon('ba') - logepsilon('fe') - (2.18d0 - 7.50d0)
            latofe = logepsilon('la') - logepsilon('fe') - (1.10d0 - 7.50d0)
            cetofe = logepsilon('ce') - logepsilon('fe') - (1.58d0 - 7.50d0)
            prtofe = logepsilon('pr') - logepsilon('fe') - (0.72d0 - 7.50d0)
            ndtofe = logepsilon('nd') - logepsilon('fe') - (1.42d0 - 7.50d0)
            smtofe = logepsilon('sm') - logepsilon('fe') - (0.96d0 - 7.50d0)
            hftofe = logepsilon('hf') - logepsilon('fe') - (0.85d0 - 7.50d0)
            pbtofe = logepsilon('pb') - logepsilon('fe') - (1.75d0 - 7.50d0)


            model(stepnum)%fetoh = fetoh

            if (fetoh > -1.2) then
                stop
            end if

            write(12,*) model(stepnum)%time,model(stepnum)%gasmass,model(stepnum)%starsmass,&
            sfrdt/timestep,serdt/timestep

            flush(12)

            write(13,*) model(stepnum)%time,fetoh,rbtofe,srtofe,ytofe,zrtofe,motofe,batofe,latofe,cetofe,&
                    prtofe,ndtofe,smtofe,hftofe,pbtofe

            flush(13)

            write(*,23) model(stepnum)%time,model(stepnum)%gasmass,model(stepnum)%starsmass,sfrdt/timestep,&
                    serdt/timestep,fetoh,serspeciesdt(3),timestep
        end if
    end do

contains
subroutine initcemodel()
    real(kind=8) :: error
    integer :: s

    imfnormfactor = 1.0d0
    imfnormfactor = 1.0d0 / integral(imf,minstellarmass,maxstellarmass,100,1.d-7,error,label='IMF normalisation')

    call initspecies
    call initevmodels

    allocate(model(1:10**6))
    do s=1,size(model)
        allocate(model(s)%speciesmassfrac(1:size(species)))
        model(s)%speciesmassfrac = 0.0d0
    end do

    model(1)%fetoh = -50
    model(1)%time = 0.0d0
    model(1)%gasmass = 1.4d5
    model(1)%starsmass = 0.0d0

    model(1)%speciesmassfrac(1) = 0.73d0
    !do s = 2, size(species)
    !    model(1)%speciesmassfrac(1) = 0.0d0
    !end do
end subroutine

subroutine initspecies()
    integer :: s, speciescount, neutrons
    character(len=3) :: stra

    open(unit=14, file="species.dat", action="read", status="old")
    read(14,'(I4)') speciescount
    allocate(species(1:speciescount))
    do s=1,size(species)
        read(14,'(I6,A6,A6,I6)') species(s)%a, species(s)%symbol, stra, neutrons
        species(s)%z = species(s)%a - neutrons
    end do
    flush(14)
end subroutine

subroutine initevmodels()
    integer :: status, i, j, s, a
    real(kind=8) :: m0, m1, miyield, milost, finalmass
    character(len=8) :: name,nameb
    character(len=3) :: stra
    character(len=1) :: mix
    character(len=3) :: el
    real(kind=8),dimension(29) :: yieldrow
    real(kind=8), dimension(8) :: eightreals
    real(kind=8) :: r !rotation parameter
    real(kind=8) :: cf88on10factor

    allocate(evmodel(4))

    do i=1,size(evmodel)
        allocate(evmodel(i)%Z(2))
        do j=1,size(evmodel(i)%Z)
            allocate(evmodel(i)%Z(j)%yield(size(species)))
            allocate(evmodel(i)%Z(j)%speciesmasslost(size(species)))
        end do
        evmodel(i)%Z(1)%yield(s) = 0.0d0
        evmodel(i)%Z(1)%speciesmasslost(s) = 0.0d0
    end do

    evmodel(1)%mass = 15.0d0
    evmodel(1)%Z(1)%fetoh = -3.8d0
    evmodel(1)%Z(1)%remnantmass = 1.74d0
    evmodel(1)%Z(2)%fetoh = -1.8d0
    evmodel(1)%Z(2)%remnantmass = 1.77d0

    evmodel(2)%mass = 20.0d0
    evmodel(2)%Z(1)%remnantmass = 2.21d0
    evmodel(2)%Z(1)%fetoh = -3.8d0
    evmodel(2)%Z(2)%fetoh = -1.8d0
    evmodel(2)%Z(2)%remnantmass = 2.22d0

    evmodel(3)%mass = 25.0d0
    evmodel(3)%Z(1)%fetoh = -3.8d0
    evmodel(3)%Z(1)%remnantmass = 2.69d0
    evmodel(3)%Z(2)%fetoh = -1.8d0
    evmodel(3)%Z(2)%remnantmass = 2.81d0

    evmodel(4)%mass = 40.0d0
    evmodel(4)%Z(1)%fetoh = -3.8d0
    evmodel(4)%Z(1)%remnantmass = 4.36d0
    evmodel(4)%Z(2)%fetoh = -1.8d0
    evmodel(4)%Z(2)%remnantmass = 4.57d0

15  format (A7,E10.2,E10.2,E10.2,E10.2,E10.2,E10.2,E10.2,E10.2)
    open(unit=7, file= "data/limongi2012.txt")
    status = 0
    do while (status==0 .or. status==5010)
        read(7,15,iostat = status) name, eightreals
        if (status==0) then
            !write(*,*) name,eightreals
            do s = 1, size(species)
                if (species(s)%z > 26) then
                    cycle
                end if
                write(stra,'(I3)') species(s)%a
                if (species(s)%z == 1) then
                    nameb = species(s)%symbol
                else
                    nameb = '^'//trim(adjustl(stra))//'^'//adjustl(species(s)%symbol)
                end if
                !write(*,*)name,nameb,StrLowCase(name)==nameb
                if (StrLowCase(name)==nameb .or. (name=='^56^Ni' .and. nameb=='^56^fe') .or. &
                    (name=='^1^H' .and. trim(adjustl(nameb))=='p')) then
                    !at Z=0, relative yields are absolute yields (except H,He)
                    evmodel(1)%Z(1)%speciesmasslost(s) = evmodel(1)%Z(1)%speciesmasslost(s) + eightreals(2)
                    evmodel(2)%Z(1)%speciesmasslost(s) = evmodel(2)%Z(1)%speciesmasslost(s) + eightreals(3)
                    evmodel(3)%Z(1)%speciesmasslost(s) = evmodel(3)%Z(1)%speciesmasslost(s) + eightreals(4)
                    evmodel(4)%Z(1)%speciesmasslost(s) = evmodel(4)%Z(1)%speciesmasslost(s) + &
                            eightreals(6)/3.0d0 + eightreals(7)*2.0d0/3.0d0

                    evmodel(1)%Z(2)%speciesmasslost(s) = evmodel(1)%Z(1)%speciesmasslost(s)
                    evmodel(2)%Z(2)%speciesmasslost(s) = evmodel(2)%Z(1)%speciesmasslost(s)
                    evmodel(3)%Z(2)%speciesmasslost(s) = evmodel(3)%Z(1)%speciesmasslost(s)
                    evmodel(4)%Z(2)%speciesmasslost(s) = evmodel(4)%Z(1)%speciesmasslost(s)
                end if
            end do
        end if
    end do

    ! Hirschi et al. rotating massive pre-SN yields
    open(unit=7, file= "data/hirschi/yields.tab")
    status = 0
    do i = 1, 62
        read(7,*)
    end do
    do while (status==0 .or. status==5010)
        read(7,*,iostat = status) name,yieldrow
        !write(*,*)name,yieldrow
        if (status==0) then
            do s = 1, size(species)
                if (species(s)%z <= 26) then
                    ! don't get light element yields from rotating models
                    cycle
                end if
                write(stra,'(I3)') species(s)%a
                nameb = trim(adjustl(StrLowCase(species(s)%symbol)))//trim(adjustl(stra))
                !write(*,*)name,nameb,trim(adjustl(name))==nameb
                if (trim(adjustl(name))==nameb) then
                    !r=0 -> no rotation, r=1.0 -> max rotation
                    r = 1.0d0
                    cf88on10factor = yieldrow(19) / yieldrow(18) !use only with rotation
                    !cf88on10factor = 1.0d0

                    ! G015zm5S003 and G015zm5S413
                    evmodel(1)%Z(1)%yield(s) = (1-r) * yieldrow(5)  + r * yieldrow(6) * cf88on10factor
                    ! G020zm5S003 and G020zm5S413
                    evmodel(2)%Z(1)%yield(s) = (1-r) * yieldrow(11) + r * yieldrow(12) * cf88on10factor
                    ! G025zm5S003 and G025zm5S413
                    evmodel(3)%Z(1)%yield(s) = (1-r) * yieldrow(17) + r * yieldrow(18) * cf88on10factor
                    ! G040zm5S418
                    evmodel(4)%Z(1)%yield(s) = yieldrow(29) * cf88on10factor

                    ! G015z01S013 and G015z01S413
                    evmodel(1)%Z(2)%yield(s) = (1-r) * yieldrow(3)  + r * yieldrow(4) * cf88on10factor
                    ! G020z01S013 and G020z01S413
                    evmodel(2)%Z(2)%yield(s) = (1-r) * yieldrow(9) + r * yieldrow(10) * cf88on10factor
                    ! G025z01S013 and G025z01S413
                    evmodel(3)%Z(2)%yield(s) = (1-r) * yieldrow(15) + r * yieldrow(16) * cf88on10factor
                    ! G040z01S413
                    evmodel(4)%Z(2)%yield(s) = yieldrow(28) * cf88on10factor

                    !write(*,*)'enhancement',species(s)%symbol,species(s)%a, yieldrow(18), yieldrow(19), cf88on10factor
                end if
            end do
        end if
    end do

    call sortevmodels()
   ! do i = 1, 20
    !    m0 = 15.0d0 + i*2.0d0
        !fetoh = 0.0d0
        !write(*,*) evmodel(i)%mass,speciesmasslostatfetoh(2,i,-1.0d0), evmodel(i)%Z(1)%speciesmasslost(2)
        !write(*,*) speciesmasslost(2, evmodel(i)%mass, 0.0d0)
     !   write(*,*) m0, speciesmasslost(3, m0, -1.0d0), remnantmass(m0)
        !write(*,*) m0, remnantmass(m0)
   ! end do
  !  stop

end subroutine

end program


