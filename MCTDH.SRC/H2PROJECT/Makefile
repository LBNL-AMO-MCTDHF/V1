

include ../Makefile.header

HTWO= init_H2_new.o  mclpmn.o mclqmn.o  h2_extraproject.o PSC-newest.o H2_params.o h2project.o lpm.o coreproject.o 

DEPEND= ../Definitions.INC  ../Definitions.ALL ../Makefile.header

default: h2projectar.a

clean:
	rm *.o *.a *.mod *genmod*

../proputils.o: ../proputils.F90
	cd ..; make proputils.o

../utils.o: ../utils.F90
	cd ..; make utils.o

../eigen.o: ../eigen.F90
	cd ..; make eigen.o

../parameters.o: ../parameters.F90
	cd ..; make parameters.o

h2projectar.a:  $(HTWO) 
	$(AR) -r h2projectar.a    $(HTWO) 

H2_params.o: $(DEPEND) H2_params.F90
	$(F90) $(FFLAGS)  -c H2_params.F90

h2project.o: $(DEPEND) h2project.F90 H2_params.o ../utils.o
	$(F90) $(FFLAGS)  -I.. -c h2project.F90

init_H2_new.o: $(DEPEND)  H2_params.o h2project.o coreproject.o PSC-newest.o ../eigen.o init_H2_new.F90
	$(F90) $(FFLAGS) -I.. -c init_H2_new.F90

coreproject.o: $(DEPEND)  H2_params.o h2project.o coreproject.F90 ../proputils.o ../parameters.o
	$(F90) $(FFLAGS) -I.. -c coreproject.F90

coretwoe.o: $(DEPEND)  HE_params.o heproject.o coretwoe.F90
	$(F90) $(FFLAGS)  -c coretwoe.F90

PSC-newest.o: $(DEPEND)  H2_params.o h2project.o PSC-newest.F90 lpm.o
	$(F90) $(FFLAGS)  -c PSC-newest.F90

lpm.o: $(DEPEND)  H2_params.o h2project.o lpm.F90
	$(F90) $(FFLAGS)  -c lpm.F90


h2_extraproject.o: $(DEPEND) H2_params.o h2project.o coreproject.o h2_extraproject.F90
	$(F90) $(FFLAGS) $(PGFFLAG) -c h2_extraproject.F90 -o h2_extraproject.o


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
