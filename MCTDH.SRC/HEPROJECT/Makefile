

include ../Makefile.header

HTWO= init_HE_new.o    he_extraproject.o  He_params.o heproject.o GetJacobiKE.o GetLobatto.o coreproject.o 

DEPEND= ../Definitions.INC  ../Definitions.ALL  ../Makefile.header

default: heprojectar.a

clean:
	rm *.o *.a *.mod *genmod*

../proputils.o: ../proputils.F90
	cd ..; make proputils.o

../eigen.o: ../eigen.F90
	cd ..; make eigen.o

../parameters.o: ../parameters.F90
	cd ..; make parameters.o

heprojectar.a:  $(HTWO) 
	$(AR) -r heprojectar.a    $(HTWO) 


He_params.o: $(DEPEND) He_params.F90
	$(F90) $(FFLAGS)  -c He_params.F90


heproject.o: $(DEPEND) heproject.F90 He_params.o
	$(F90) $(FFLAGS)  -c heproject.F90

init_HE_new.o: $(DEPEND)  He_params.o heproject.o coreproject.o ../eigen.o init_HE_new.F90
	$(F90) $(FFLAGS) -I.. -c init_HE_new.F90

coreproject.o: $(DEPEND)  He_params.o heproject.o coreproject.F90 ../proputils.o ../parameters.o
	$(F90) $(FFLAGS) -I.. -c coreproject.F90

he_extraproject.o: $(DEPEND) He_params.o heproject.o coreproject.o he_extraproject.F90
	$(F90) $(FFLAGS) $(PGFFLAG) -c he_extraproject.F90 -o he_extraproject.o


GetJacobiKE.o: GetJacobiKE.F90 
	$(F90) $(FFLAGS) -c GetJacobiKE.F90 

GetLobatto.o: $(DEPEND) GetLobatto.F90
	$(F90) $(FFLAGS) -c GetLobatto.F90


jacobi.o: jacobi.f 
	$(FC) -O0 -c jacobi.f
#	$(FC) $(F77FLAGS) -c jacobi.f

gaussq.o: gaussq.f 
	$(FC) $(F77FLAGS) -c gaussq.f

mclpmn.o: mclpmn.f
	$(FC) -O0 -c mclpmn.f

mclqmn.o: mclqmn.f
	$(FC) -O0 -c mclqmn.f


newfluxcall.o: newfluxcall.F90 
	$(F90) $(FFLAGS)  -c newfluxcall.F90 

newflux.o: $(DEPEND) newflux.F90 
	$(F90) $(FFLAGS)  -c newflux.F90 
