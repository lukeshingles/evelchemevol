program runcemodel
    use cemodelsetup
    use cemodel

    implicit none

    integer(kind=4) :: steplastoutput = 0
    real(kind=8) :: serdt, sfrdt, error
    real(kind=8) :: fetoh
    real(kind=8), allocatable :: serspeciesdt(:)
    integer :: speciesnum
    logical :: acceptableerrors

    call initcemodel

    allocate(serspeciesdt(size(species)))

    open(unit=12, file=trim(filespath) // "out-cemodel.txt", action="write", status="replace")
    open(unit=13, file=trim(filespath) // "out-abund.txt", action="write", status="replace")
23  format (E13.4','E10.3','E10.3','E10.3','E10.3','F6.2','E7.1','E10.3) !console output

    write(12,*) '#time  Mgas    Mstars  SFR SER'
    write(13,*) '#time    [Fe/H]  Fe  Rb  Sr  Y  Zr  Mo  Ba  La  Ce &
            Pr  Nd  Sm  Eu  Gd  Hf  Pb'

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
                    label='ser time integral')
            if (error/timestep > sermaxerror .and. timestep > mintimestep) then
                write(*,*),'stellar ejection rate retained error too high:',error/timestep
                goto 9000
            end if

            do speciesnum = 1, size(species)
                serspeciesdt(speciesnum) = integral(stellarejectionrateofspecies,model(stepnum-1)%time,&
                        model(stepnum)%time,1,1.d0,error,intarg=speciesnum,label='ser species time integral')
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

            if (fetoh > -1.2) then
                stop
            end if

            write(12,*) model(stepnum)%time,model(stepnum)%gasmass,model(stepnum)%starsmass,&
            sfrdt/timestep,serdt/timestep

            flush(12)

            write(13,*) model(stepnum)%time,fetoh,logepsilon('fe'),logepsilon('rb'),logepsilon('sr'),&
                logepsilon('y'),logepsilon('zr'),logepsilon('mo'),logepsilon('ba'),logepsilon('la'),&
                logepsilon('ce'),logepsilon('pr'),logepsilon('nd'),logepsilon('sm'),&
                logepsilon('eu'),logepsilon('gd'),logepsilon('hf'),logepsilon('pb')

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

    model(1)%time = 0.0d0
    model(1)%gasmass = 1.4d5
    model(1)%starsmass = 0.0d0

    model(1)%speciesmassfrac(1) = 0.73d0
    !do s = 2, size(species)
    !    model(1)%speciesmassfrac(1) = 0.0d0
    !end do
end subroutine


end program


