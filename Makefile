# This makefile is used to produce quickly the code 
# We successfully implemented the random number generator.
#
#

FC=gfortran-mp-4.9
LIBPATH=-L$(HOME)/.lib/
LIBS=$(LIBPATH) -llapack -lblas

all: matrixlapack.o GFEOM.o   
	$(FC)  matrixlapack.o GFEOM.o $(LIBS) -o EXE_GFEOM

matrixlapack.o:
	$(FC) -c matrixlapack.f90
 
GFEOM.o: matrixlapack.o
	$(FC) -c GFEOM.f90

clean:
	rm -f *\.mod *\.o *~
	
