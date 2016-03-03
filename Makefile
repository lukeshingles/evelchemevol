FC=gfortran
FFLAGS=-fopenmp -Ofast -Wunused -floop-parallelize-all -ftree-parallelize-loops=8 -fcheck=all
FDEBUGFLAGS=-fcheck=all

all:	clean echemevol

clean:
	rm -f *.mod echemevol

echemevol:
	$(FC) $(FFLAGS) -o echemevol cemodel.f90 main.f90

echemevolifort:
	ifort -openmp -O3 -parallel -o echemevol cemodel.f90 main.f90

echemevol-debug:
	$(FC) $(FFLAGS) $(FDEBUGFLAGS) -o echemevol cemodel.f90 main.f90
