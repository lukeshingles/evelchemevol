rm runcemodel > /dev/null
gfortran -fopenmp -Ofast -o runcemodel stringutility.f90 cemodelsetup.f90 cemodel.f90 runcemodel.f90
./runcemodel