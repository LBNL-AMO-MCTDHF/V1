

include ../Makefile.header

HTWO= init_sinc_new.o    sinc_extraproject.o  sinc_params.o sincproject.o SineDVR.o coreproject.o readit.o tooth.o toeplitz_subs.o ftcore.o  cooleytukey_shared.o cooleytukey_mpi.o

DEPEND= ../Definitions.INC  ../Definitions.ALL ../Makefile.header

default: sincprojectar.a

clean:
	rm *.o *.a *genmod* *.mod

sincprojectar.a:  $(HTWO) 
	$(AR) -r sincprojectar.a    $(HTWO) 


sinc_params.o: $(DEPEND) sinc_params.F90
	$(F90) $(FFLAGS)  -c sinc_params.F90

toeplitz_subs.o: $(DEPEND) toeplitz_subs.F90
	$(F90) $(FFLAGS) $(MPIFLAG)  -c toeplitz_subs.F90 -I$(MYINCLUDE)

ftcore.o: $(DEPEND) ftcore.F90
	$(F90) $(FFLAGS) $(FFTFLAG) $(MPIFLAG)  -c ftcore.F90 -I$(MYINCLUDE)

cooleytukey_shared.o: $(DEPEND) cooleytukey_shared.F90
	$(F90) $(FFLAGS) $(MPIFLAG)  -c cooleytukey_shared.F90 -I$(MYINCLUDE)

cooleytukey_mpi.o: $(DEPEND) cooleytukey_mpi.F90 cooleytukey_shared.o
	$(F90) $(FFLAGS)  -c cooleytukey_mpi.F90 -I$(MYINCLUDE)

readit.o: $(DEPEND) sinc_params.o sincproject.o tooth.o readit.F90
	$(F90) $(FFLAGS)  -c readit.F90 

#runtooth.o: $(DEPEND) sinc_params.o sincproject.o tooth.o runtooth.F90
#	$(F90) $(FFLAGS)  -c runtooth.F90 

tooth.o: $(DEPEND) tooth.f90
	$(F90) $(FFLAGS)  -c tooth.f90 

sincproject.o: $(DEPEND) sincproject.F90 sinc_params.o
	$(F90) $(FFLAGS)  -c sincproject.F90

init_sinc_new.o: $(DEPEND)  sinc_params.o sincproject.o init_sinc_new.F90
	$(F90) $(FFLAGS)  -c init_sinc_new.F90

coreproject.o: $(DEPEND)  sinc_params.o sincproject.o coreproject.F90
	$(F90) $(FFLAGS)  $(MPIFLAG) -c coreproject.F90

sinc_extraproject.o: $(DEPEND) sinc_params.o sincproject.o sinc_extraproject.F90
	$(F90) $(FFLAGS) $(MPIFLAG) $(PGFFLAG) -c sinc_extraproject.F90 -o sinc_extraproject.o


SineDVR.o: SineDVR.F90 
	$(F90) $(FFLAGS) -c SineDVR.F90 

