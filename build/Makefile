VPATH=../src:../test:../ParseText

F90FLAGS=
F90=gfortran
ifeq ($(F90),gfortran)
F90FLAGS+= -ffree-line-length-none
else ifeq ($(F90),ifort)
F90FLAGS+= -O3 -axT
endif
DEBUG=F
ifeq ($(DEBUG),T)
F90FLAGS+= -g
endif

%.o : %.f90
	$(F90) $(F90FLAGS) -c $<

test_group.o : ParseText.o group.o

MPCDMD_test : ParseText.o sys.o group.o MPCD.o LJ.o MD.o MPCDMD_test.o
	$(F90) $(F90FLAGS) ParseText.o sys.o group.o MPCD.o LJ.o MD.o MPCDMD_test.o -o MPCDMD_test

clean:
	rm *.o *.mod
