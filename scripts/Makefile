VPATH=../src:../test:../ParseText:../mtprng

F90OPT=
F90FLAGS=$(F90OPT)
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

MPCDMD_test : stdtypes.o mtprng.o ParseText.o sys.o group.o MPCD.o LJ.o MD.o MPCDMD_test.o
	$(F90) $(F90FLAGS) stdtypes.o mtprng.o ParseText.o sys.o group.o MPCD.o LJ.o MD.o MPCDMD_test.o -o MPCDMD_test

MPCD_test : stdtypes.o mtprng.o ParseText.o sys.o group.o MPCD.o MPCD_test.o
	$(F90) $(F90FLAGS) stdtypes.o mtprng.o ParseText.o sys.o group.o MPCD.o MPCD_test.o -o MPCD_test 

clean:
	rm *.o *.mod
