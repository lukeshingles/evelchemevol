module cemodelsetup
    real(kind=8) :: minstellarmass = 15.0d0, maxstellarmass = 40.0d0
    real(kind=8) :: stellarmasslow, stellarmasshigh
    character(len=50), parameter :: filespath = "./"
    real(kind=8), parameter :: mintimestep = 1.0d1, maxtimestep = 6.0d3
    real(kind=8), parameter :: sermaxerror=1.d-6
    real(kind=8) :: imfnormfactor = 1.0 ! to be replaced during initialisation

contains

! returns star formation rate in solar masses at time in years
real(kind=8) function starformationrate(time)
    real(kind=8), intent(in) :: time

    if (time < 0.d0 .or. time > 1.d8) then
        starformationrate = 0.d0
    else
        starformationrate = 0.10d0  * exp(-time/3.0d4)
    end if
end function

! returns main sequence lifetime in years of a star
! with an initial mass in solar masses
real(kind=8) function stellarlifetime(mass)
    implicit none
    real(kind=8), intent(in) :: mass
    stellarlifetime = 1.d+10 * (mass ** (-2.9d0))
end function

! initial mass function (dN/dM) normalised s.t. integral of imf(m) from m=minstellarmass to m=maxstellarmass Msun is 1.0
real(kind=8) function imf(mass)
    real(kind=8), intent(in) :: mass

    ! Modified Kroupa, Tout, Gilmore 1993 IMF
    if (mass >= minstellarmass .AND. mass <= maxstellarmass) then
        if (mass > 1.0d0) then
            imf = imfnormfactor * mass ** (-2.7d0)
        elseif (mass > 0.5d0) then
            imf = imfnormfactor * mass ** (-2.2d0)
        elseif (mass > 0.08d0) then
            imf = imfnormfactor * mass ** (-1.50d0) / (0.5 ** (0.7d0))
        else
            imf = 0.0d0
        end if
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
end function

end module