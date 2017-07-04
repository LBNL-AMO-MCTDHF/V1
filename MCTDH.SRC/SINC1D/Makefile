

include ../Makefile.header

HTWO= mpimod.o sinc1d_params.o sincproject.o init_sinc_new.o    sinc_extraproject.o   SineDVR.o coreproject.o \
	toeplitz_subs.o ftcore.o cooleytukey_mpi.o

DEPEND= ../Definitions.INC  ../Definitions.ALL ../Makefile.header sinc1d_params.o mpimod.o Makefile

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

../eigen.o: ../eigen.F90
	cd ..; make eigen.o

sincprojectar.a:  $(HTWO) 
	$(AR) -r sincprojectar.a    $(HTWO) 

sinc1d_params.o: sinc1d_params.F90
	$(F90) $(FFLAGS) $(MPIFLAG) -c sinc1d_params.F90

mpimod.o: mpimod.f90 ../Definitions.INC  ../Definitions.ALL ../Makefile.header sinc1d_params.o Makefile
	$(F90) $(FFLAGS) $(MPIFLAG) -c mpimod.f90

toeplitz_subs.o: $(DEPEND) toeplitz_subs.F90 cooleytukey_mpi.o ftcore.o ../parameters.o
	$(F90) $(FFLAGS) $(MPIFLAG) -I.. -c toeplitz_subs.F90

cooleytukey_mpi.o: $(DEPEND) cooleytukey_mpi.F90 ftcore.o ../MPI.o
	$(F90) $(FFLAGS) $(MPIFLAG) -I.. -c cooleytukey_mpi.F90

ftcore.o: $(DEPEND) ftcore.F90
	$(F90) $(FFLAGS) $(MPIFLAG) $(FFTFLAG) -c ftcore.F90 -I$(MYINCLUDE)

sincproject.o: $(DEPEND) sincproject.F90 sinc1d_params.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c sincproject.F90

init_sinc_new.o: $(DEPEND) sinc1d_params.o sincproject.o coreproject.o init_sinc_new.F90 ../MPI.o ../proputils.o ../parameters.o ../blocklanczos.o ../eigen.o
	$(F90) $(FFLAGS) $(MPIFLAG) -I.. -c init_sinc_new.F90

coreproject.o: $(DEPEND)  sinc1d_params.o sincproject.o coreproject.F90 toeplitz_subs.o ../MPI.o ../proputils.o ../parameters.o
	$(F90) $(FFLAGS) $(MPIFLAG) -I.. -c coreproject.F90

sinc_extraproject.o: $(DEPEND) sinc1d_params.o sincproject.o sinc_extraproject.F90
	$(F90) $(FFLAGS) $(MPIFLAG) $(PGFFLAG) -c sinc_extraproject.F90 -o sinc_extraproject.o

SineDVR.o: $(DEPEND) SineDVR.F90 
	$(F90) $(FFLAGS) $(MPIFLAG) -c SineDVR.F90 

