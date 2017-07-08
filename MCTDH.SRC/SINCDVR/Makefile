

include ../Makefile.header

HTWO= mpimod.o sinc_params.o sincproject.o init_sinc_new.o    sinc_extraproject.o   SineDVR.o coreproject.o \
	tooth.o toeplitz_subs.o ftcore.o cooleytukey_mpi.o

DEPEND= ../Definitions.INC  ../Definitions.ALL ../Makefile.header sinc_params.o mpimod.o Makefile

default: sincprojectar.a

clean:
	rm *.o *.a *genmod* *.mod

../MPI.o: ../MPI.F90
	cd ..; make MPI.o

../proputils.o: ../proputils.F90
	cd ..; make proputils.o

../parameters.o: ../parameters.F90
	cd ..; make parameters.o

../blocklanczos.o: ../blocklanczos.F90
	cd ..; make blocklanczos.o

sincprojectar.a: $(HTWO) 
	$(AR) -r sincprojectar.a    $(HTWO) 

sinc_params.o: sinc_params.F90
	$(F90) $(FFLAGS) $(MPIFLAG) -c sinc_params.F90

mpimod.o: mpimod.f90 ../Definitions.INC  ../Definitions.ALL ../Makefile.header sinc_params.o Makefile
	$(F90) $(FFLAGS) $(MPIFLAG) -c mpimod.f90

toeplitz_subs.o: $(DEPEND) toeplitz_subs.F90 ftcore.o cooleytukey_mpi.o ../MPI.o ../parameters.o
	$(F90) $(FFLAGS) $(MPIFLAG) -I.. -c toeplitz_subs.F90

ftcore.o: $(DEPEND) ftcore.F90
	$(F90) $(FFLAGS) $(MPIFLAG) $(FFTFLAG) -c ftcore.F90 -I$(MYINCLUDE)

cooleytukey_mpi.o: $(DEPEND) cooleytukey_mpi.F90 ftcore.o ../MPI.o
	$(F90) $(FFLAGS) $(MPIFLAG) -I..  -c cooleytukey_mpi.F90

#readit.o: $(DEPEND) sinc_params.o sincproject.o tooth.o readit.F90 ../MPI.o
#	$(F90) $(FFLAGS) $(MPIFLAG) -I.. -c readit.F90 

#runtooth.o: $(DEPEND) sinc_params.o sincproject.o tooth.o runtooth.F90
#	$(F90) $(FFLAGS)  -c runtooth.F90 

tooth.o: $(DEPEND) tooth.F90
	$(F90) $(FFLAGS) $(MPIFLAG) -c tooth.F90 

sincproject.o: $(DEPEND) sincproject.F90 sinc_params.o tooth.o ../MPI.o
	$(F90) $(FFLAGS) $(MPIFLAG) -I.. -c sincproject.F90

init_sinc_new.o: $(DEPEND) sinc_params.o sincproject.o init_sinc_new.F90 coreproject.o ../MPI.o ../proputils.o ../utils.o ../blocklanczos.o
	$(F90) $(FFLAGS) $(MPIFLAG) -I.. -c init_sinc_new.F90

coreproject.o: $(DEPEND) sinc_params.o sincproject.o coreproject.F90 toeplitz_subs.o ../MPI.o ../proputils.o ../parameters.o
	$(F90) $(FFLAGS) $(MPIFLAG) -I.. -c coreproject.F90

sinc_extraproject.o: $(DEPEND) sinc_params.o sincproject.o sinc_extraproject.F90
	$(F90) $(FFLAGS) $(MPIFLAG) $(PGFFLAG) -c sinc_extraproject.F90 -o sinc_extraproject.o

SineDVR.o: $(DEPEND) SineDVR.F90 
	$(F90) $(FFLAGS) $(MPIFLAG) -c SineDVR.F90 

