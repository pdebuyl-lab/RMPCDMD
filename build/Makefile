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

#F90LIBS=-L. -ldcmt

libdcmt.a :
	@cd ../dcmt/lib ; $(MAKE) lib
	@cp ../dcmt/lib/libdcmt.a ./

%.o : %.f90
	$(F90) $(F90FLAGS) -c $< $(F90LIBS)

MPCDMD_test : stdtypes.o mtprng.o ParseText.o sys.o group.o MPCD.o LJ.o MD.o MPCDMD_test.o
	$(F90) $(F90FLAGS) stdtypes.o mtprng.o ParseText.o sys.o group.o MPCD.o LJ.o MD.o MPCDMD_test.o -o MPCDMD_test $(F90LIBS)

MPCD_test : stdtypes.o mtprng.o ParseText.o sys.o group.o MPCD.o MPCD_test.o
	$(F90) $(F90FLAGS) stdtypes.o mtprng.o ParseText.o sys.o group.o MPCD.o MPCD_test.o -o MPCD_test $(F90LIBS)

clean:
	rm *.o *.mod
