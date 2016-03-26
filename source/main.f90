program echemevol
    use cemodel

    implicit none

    integer(kind=4) :: steplastoutput = 1
    real*8 :: serdt, sfrdt
    real*8 :: fetoh
    real*8 :: systemtimeelapsed
    integer*8 :: sysclockstart = 0, sysclocknow, sysclockcountrate, sysclockcountmax
    real*8, allocatable :: serspeciesdt(:)
    integer :: speciesnum
    integer :: iterations
    logical :: acceptableerrors

    open(unit=14, file=trim(filespath) // "out-log.txt", action="write", status="replace")

    call initcemodel

    allocate(serspeciesdt(size(species)))

    open(unit=12, file=trim(filespath) // "out-cemodel.txt", action="write", status="replace")
    open(unit=13, file=trim(filespath) // "out-abundances.txt", action="write", status="replace")

    write(12,'(A)') '#stepnum         time         Mgas       Mstars         SFR          SER        [Fe/H]'
    write(13,'(A)') '#stepnum         time [Fe/H]      O     Na     Fe     Rb     Sr      Y     Zr     Ba     La     Ce&
           &     Pr     Nd     Sm     Eu     Pb'

    timestep = initialtimestep
    do stepnum = 1, size(model)
        iterations = 0
        acceptableerrors = .false.
        do while ((acceptableerrors .eqv. .false.) .and. (stepnum > 1)) !skip first model
            iterations = iterations + 1
            acceptableerrors = .true.
            model(stepnum)%time = model(stepnum-1)%time + timestep

            if (model(stepnum)%time > maxtime) then
                write(14,'(A,ES11.4,A)') 't > ',maxtime,', finished!'
                goto 950
            end if

            ! all masses
            stellarmasslow = minstellarmass
            stellarmasshigh = maxstellarmass

            !more accurate if SFR varies between model steps
            !sfrdt = integral(starformationrate,model(stepnum-1)%time,model(stepnum)%time,1,1.d0,&
            !        error,label='sfr time integral')
            !if (error/timestep > maxerrorser .and. timestep > mintimestep) then
            !    write(14,*),'Restepping: star formation rate error too high:',error/timestep
            !    goto 900
            !end if

            sfrdt = model(stepnum-1)%sfr * timestep

            !serdt = integral(stellarejectionrate,model(stepnum-1)%time,model(stepnum)%time,10,&
            !        max(maxchangeser*model(stepnum-1)%ser,1d-8),error,label='ser time integral')
            !model(stepnum)%ser = serdt/timestep

            ! this variable is read by stellarejectionrate and stellarejectionrateofspecies
            ! make sure it is accurate enough to determined whether dSER/SER > maxchangeser
            maxerrorser = maxchangeser / 2.d0

            model(stepnum)%ser = stellarejectionrate(model(stepnum)%time)
            serdt = model(stepnum)%ser * timestep

            if (model(stepnum-1)%ser > 1.0d-8 .and. model(stepnum)%ser > 1.0e-8) then
                if (abs(model(stepnum)%ser/model(stepnum-1)%ser - 1) > maxchangeser &
                        .and. timestep > mintimestep) then
                    write(14,'(A,ES9.2,A,ES9.2)') 'Restepping: dSER/SER = ',&
                            abs(model(stepnum)%ser/model(stepnum-1)%ser - 1),' > ',maxchangeser
                    goto 900
                end if
            end if

            !$omp parallel do simd
            do speciesnum = 1, size(species)
                !serspeciesdt(speciesnum) = integral(stellarejectionrateofspecies,model(stepnum-1)%time,&
                !        model(stepnum)%time,10,maxchangeser,intarg=speciesnum,&
                !        label='ser species time integral')
                serspeciesdt(speciesnum) = stellarejectionrateofspecies(model(stepnum)%time, speciesnum) * timestep
            end do
            !$omp end parallel do simd

            model(stepnum)%starsmass = model(stepnum-1)%starsmass + sfrdt - serdt
            model(stepnum)%gasmass = model(stepnum-1)%gasmass - sfrdt + serdt

            if (model(stepnum)%gasmass < 0.0d0 .and. timestep > mintimestep) then
                write(14,'(A,ES11.4,A,ES9.2)') 'Restepping: Negative gas mass. Mgas=',model(stepnum)%gasmass,', SFR*dt=',sfrdt
                goto 900
            end if

            if (abs(model(stepnum)%gasmass/model(stepnum-1)%gasmass - 1) > maxchangemgas .and. timestep > mintimestep) then
                write(14,'(A,ES9.2,A,ES9.2)') 'Restepping: dMgas/Mgas = ',&
                        abs(model(stepnum)%gasmass/model(stepnum-1)%gasmass - 1),' > ',maxchangemgas
                goto 900
            end if

            do concurrent (speciesnum = 1:size(species))
                model(stepnum)%speciesmassfrac(speciesnum) = (1.0d0 / model(stepnum)%gasmass) * &
                    (model(stepnum-1)%speciesmassfrac(speciesnum) * model(stepnum-1)%gasmass + &
                    serspeciesdt(speciesnum) - model(stepnum-1)%speciesmassfrac(speciesnum) * sfrdt)
            end do

            timestep = min(timestep * (1.00d0 + timestepgrowthrate),maxtimestep)
            exit
    900     continue
            acceptableerrors = .false.
            timestep = max(timestep * 0.5d0, mintimestep)
            cycle
        end do

        ! will be overwritten with zero afterwards if SF is switched off
        model(stepnum)%sfr = model(stepnum)%gasmass * sfrefficiency

        if (stepnum == 1) then
            if (model(stepnum)%gasmass > sfrstartmgasabove) then
                write(14,'(A,ES9.2)') 'Starting star formation, Mgas >',sfrstartmgasabove
            else
                model(stepnum)%sfr = 0.0d0
            end if
        elseif (model(stepnum-1)%sfr > 0.0d0) then !SF was active in the previous step
            if (model(stepnum)%gasmass < sfrendmgasbelow) then
                write(14,'(A,ES9.2)') 'Ending star formation, Mgas <',sfrendmgasbelow
                model(stepnum)%sfr = 0.0d0
            end if
        elseif (model(stepnum-1)%sfr <= 0.0d0) then !no SF in previous step
            if (model(stepnum)%gasmass > sfrstartmgasabove) then
                write(14,'(A,ES9.2)') 'Starting star formation, Mgas >',sfrstartmgasabove
            else
                model(stepnum)%sfr = 0.0d0
            end if
        end if

        if (stepnum - steplastoutput >= 1 .or. stepnum == 1) then

            ! Asplund et al. (ARAA, 2009)
            fetoh = logepsilon('fe') - 7.50d0

            call SYSTEM_CLOCK(sysclocknow, sysclockcountrate, sysclockcountmax)
            if (sysclockstart /= 0) then
                systemtimeelapsed = float(sysclocknow - sysclockstart) / float(sysclockcountrate)
            else
                systemtimeelapsed = 0.0d0
            end if
            sysclockstart = sysclocknow

            ! structure output
            write(12,'(I7,1X,*(ES13.4))') stepnum,model(stepnum)%time,model(stepnum)%gasmass,model(stepnum)%starsmass,&
            model(stepnum)%sfr,serdt/timestep,fetoh
            flush(12)

            ! abundance output
            write(13,'(I7,1X,ES13.4,*(F7.3))') stepnum,model(stepnum)%time,fetoh,logepsilon('o'),logepsilon('na'),&
                logepsilon('fe'),logepsilon('rb'),&
                logepsilon('sr'),logepsilon('y'),logepsilon('zr'),logepsilon('ba'),&
                logepsilon('la'),logepsilon('ce'),logepsilon('pr'),logepsilon('nd'),logepsilon('sm'),&
                logepsilon('eu'),logepsilon('pb')
            flush(13)

            ! log output
            write(14,'(A)')'-------------------------'
            write(14,'(A,ES11.4,A,I6,A,ES11.4)') 't=',model(stepnum)%time,', StepNum=',stepnum,', dt=',timestep
            write(14,'(A,ES11.4,A,ES11.4)') 'Mgas=',model(stepnum)%gasmass,', Mstars=',model(stepnum)%starsmass
            write(14,'(A,ES11.4,A,ES11.4)') 'SFR=',model(stepnum)%sfr,', SER=',model(stepnum)%ser
            write(14,'(A,F6.2,A,ES11.4,A,ES11.4)') '[Fe/H]=',fetoh,', X(Sr)=',elmassfrac('sr'),', X(Ba)=',elmassfrac('ba')
            write(14,'(A,F6.2)') 'e(CNO)=',logepsiloncno()
            if (stepnum > 1) then
                write(14,'(A,F9.4)') 'performance (models/sec)=',(stepnum - steplastoutput)/systemtimeelapsed
                write(14,'(A,ES11.2)') 'performance (years/sec)=',&
                    (model(stepnum)%time - model(steplastoutput)%time)/systemtimeelapsed
                !if (iterations > 1) then
                !    write(14,'(A,I4)') 'iterations=',iterations
                !end if
            end if

            flush(14)
            steplastoutput = stepnum
        end if
    end do

950 call closefiles

contains
subroutine initcemodel()
    integer :: s

    write(14,'(A)',advance="no") 'Normalising IMF...'
    flush(14)
    imfnormfactor = 1.0d0 ! a value is required to calculate the following integral
    imfnormfactor = 1.0d0 / integral(mimf,0.1d0,100.0d0,100,1.d-6,label='IMF normalisation')
    write(14,'(A)') 'done.'
    flush(14)

    call initconfig
    call initspecies
    call inityields

    allocate(model(1:maxmodelnum))
    do s=1,size(model)
        allocate(model(s)%speciesmassfrac(1:size(species)))
        model(s)%speciesmassfrac = 0.0d0
    end do

    call loadinitialcomposition

    model(1)%time = 0.0d0
    model(1)%gasmass = initialgasmass
    model(1)%starsmass = 0.0d0
    model(1)%ser = 0.0d0

end subroutine

end program


