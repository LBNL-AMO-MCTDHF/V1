
# for checkmodules to work, all .o instructions must have object.o:  etc.  not object.o : etc.

include Makefile.header

include Name.Txt

HEPROJECT=HEPROJECT/heprojectar.a
H2PROJECT=H2PROJECT/h2projectar.a
SINCPROJECT=SINCDVR/sincprojectar.a

include Name.Txt

include DFFTPACK/Dfftfiles.Inc

include DGMRES/Dgmfiles.Inc


bessrcso = bessel/d1mach.o bessel/d9lgmc.o bessel/dasyjy.o bessel/dbesj.o bessel/dcsevl.o bessel/dgamlm.o bessel/dgamln.o bessel/dgamma.o bessel/djairy.o bessel/dlngam.o bessel/fdump.o bessel/i1mach.o bessel/initds.o bessel/j4save.o bessel/xercnt.o bessel/xerhlt.o bessel/xermsg.o bessel/xerprn.o bessel/xersve.o bessel/xgetua.o bessel/drc3jj.o


SRCS= parameters.o  gaussq.o  eigen.o jacobi.o  spin.o rkf45.o   spfs.o  main_modules.o spinwalks.o ovlsub.o autocall.o utils.o odex.o   getparams.o MPI.o   actions.o   readactions.o povactions.o saveactions.o natprojaction.o orbvectoractions.o   expokit.o  expo_driver.o quad.o   dipolecall.o    dipolesub.o  proputils.o autosub.o  derivs.o prop.o matel.o  projeflux.o loadstuff.o configload.o blocklanczos.o psistats.o second_derivs.o keprojector.o opfunctions.o electronflux.o walks.o walkmult.o  newconfig.o configstuff.o biortho.o    dfconstrain.o denmat.o mean.o driving.o $(ARGO) 


DEPEND=  main_modules.o parameters.o 
DEPEND= Definitions.INC main_modules.o parameters.o Definitions.ALL Makefile.header 

default: $(MYDEFAULT)

clean:
	rm *.o *.mod *~  $(NAME)_atom $(NAME)_diatom $(NAME)_sinc *genmod*

##double colon
#linkme::
#	./linkme $(NAME)

STUFFS= mctdhf.o $(SRCS) $(bessrcso) $(DFFTFILES) $(dgmfiles) $(dgmparfiles)

$(NAME)_diatom:   $(STUFFS)  $(H2PROJECT) 
	$(F90) $(LOADFLAGS) -o $(NAME)_diatom $(STUFFS) $(H2PROJECT)  $(LAPACK)

$(NAME)_atom:   $(STUFFS)  $(HEPROJECT)
	$(F90) $(LOADFLAGS) -o $(NAME)_atom  $(STUFFS) $(HEPROJECT)  $(LAPACK)

$(NAME)_sinc:   $(STUFFS)  $(SINCPROJECT)
	$(F90) $(LOADFLAGS) -o $(NAME)_sinc  $(STUFFS) $(SINCPROJECT)  $(LAPACK)

parameters.o:   Definitions.INC Definitions.ALL parameters.F90 
	$(F90) $(FFLAGS) $(MPIFLAG) -c parameters.F90

main_modules.o: Definitions.INC Definitions.ALL parameters.o main_modules.F90 
	$(F90) $(FFLAGS) $(PGFFLAG) -c main_modules.F90

biortho.o:      Definitions.INC  Definitions.ALL parameters.o main_modules.o biortho.F90 
	$(F90) $(FFLAGS)  -c biortho.F90 

biortho2.o:      Definitions.INC  Definitions.ALL parameters.o main_modules.o biortho2.F90 
	$(F90) $(FFLAGS)  -c biortho2.F90 



second_derivs.o: $(DEPEND) second_derivs.F90
	$(F90) $(FFLAGS)  -c second_derivs.F90 

keprojector.o: $(DEPEND) keprojector.F90
	$(F90) $(FFLAGS)  -c keprojector.F90 

loadstuff.o: $(DEPEND) loadstuff.F90
	$(F90) $(FFLAGS)  -c loadstuff.F90 

configload.o: $(DEPEND) configload.F90
	$(F90) $(FFLAGS)  -c configload.F90 

psistats.o: $(DEPEND) psistats.F90
	$(F90) $(FFLAGS)  -c psistats.F90 

dfconstrain.o: $(DEPEND) dfconstrain.F90
	$(F90) $(FFLAGS)  -c dfconstrain.F90 

expo_driver.o: $(DEPEND) expo_driver.F90
	$(F90) $(FFLAGS)  -c expo_driver.F90 

quad.o: $(DEPEND) quad.F90
	$(F90) $(FFLAGS)  -c quad.F90 

dgmres_driver.o: $(DEPEND) implicit.o dgmres_driver.F90
	$(F90) $(FFLAGS)  -c dgmres_driver.F90 

BCG_driver.o: $(DEPEND) implicit.o BCG_driver.F90
	$(F90) $(FFLAGS)  -c BCG_driver.F90 




mctdhf.o: $(DEPEND) mctdhf.F90 
	$(F90) $(FFLAGS)  -c mctdhf.F90 

MCSCF-assemble.o: $(DEPEND) MCSCF-assemble.F90 
	$(F90) $(FFLAGS)  -c MCSCF-assemble.F90 

MCSCF-matel.o: $(DEPEND) MCSCF-matel.F90 
	$(F90) $(FFLAGS) $(MPIFLAG) -c MCSCF-matel.F90 


actions.o: $(DEPEND) actions.F90 
	$(F90) $(FFLAGS) $(MATHEMAT) $(POVRAY) $(NOCFLAG) -c actions.F90 

readactions.o: $(DEPEND) readactions.F90 
	$(F90) $(FFLAGS) $(MATHEMAT) $(POVRAY) $(NOCFLAG) -c readactions.F90 

saveactions.o: $(DEPEND) saveactions.F90 
	$(F90) $(FFLAGS) $(MATHEMAT) $(POVRAY) $(NOCFLAG) -c saveactions.F90 

natprojaction.o: $(DEPEND) natprojaction.F90 
	$(F90) $(FFLAGS) $(MATHEMAT) $(POVRAY) $(NOCFLAG) -c natprojaction.F90 

orbvectoractions.o: $(DEPEND) orbvectoractions.F90 
	$(F90) $(FFLAGS) $(MATHEMAT) $(POVRAY) $(NOCFLAG) -c orbvectoractions.F90 

povactions.o: $(DEPEND) povactions.F90 
	$(F90) $(FFLAGS) $(MATHEMAT) $(POVRAY) $(NOCFLAG) -c povactions.F90 


prop.o: $(DEPEND) prop.F90  
	$(F90) $(FFLAGS) -c prop.F90 


xarr.o: $(DEPEND) xarr.F90  
	$(F90) $(FFLAGS) -c xarr.F90 


proputils.o: $(DEPEND) proputils.F90  
	$(F90) $(FFLAGS) -c proputils.F90 


derivs.o: $(DEPEND) derivs.F90  
	$(F90) $(FFLAGS) -c derivs.F90 

mcscf.o: $(DEPEND) mcscf.F90  
	$(F90) $(FFLAGS) -c mcscf.F90 


MPI.o:  $(DEPEND) MPI.F90
	$(F90) $(MPIFLAG) $(FFLAGS) -c MPI.F90


mean.o: $(DEPEND) mean.F90  
	$(F90) $(FFLAGS)  -c mean.F90 

matel.o: $(DEPEND)  matel.F90
	$(F90) $(FFLAGS) -c matel.F90 

transposematel.o: $(DEPEND)  transposematel.F90
	$(F90) $(FFLAGS) -c transposematel.F90 


spin.o: $(DEPEND)  spin.F90
	$(F90) $(FFLAGS) -c spin.F90 

spfs.o: $(DEPEND) spfs.F90
	$(F90) $(FFLAGS)  -c spfs.F90

denmat.o: $(DEPEND) denmat.F90 biortho.o
	$(F90) $(FFLAGS)  -c denmat.F90

driving.o: $(DEPEND) driving.F90 biortho.o
	$(F90) $(FFLAGS)  -c driving.F90

newconfig.o: $(DEPEND) newconfig.F90
	$(F90) $(FFLAGS)  -c newconfig.F90

configstuff.o: $(DEPEND) configstuff.F90
	$(F90) $(FFLAGS)  -c configstuff.F90

walks.o: $(DEPEND) walks.F90 
	$(F90) $(FFLAGS)  -c walks.F90

walkmult.o: $(DEPEND) walkmult.F90 
	$(F90) $(FFLAGS)  -c walkmult.F90

spinwalks.o: $(DEPEND) spinwalks.F90
	$(F90) $(FFLAGS)  -c spinwalks.F90

ovlsub.o: $(DEPEND) ovlsub.F90
	$(F90) $(FFLAGS)  -c ovlsub.F90


autosub.o: $(DEPEND) autosub.F90 biortho.o
	$(F90) $(FFLAGS)  -c autosub.F90

autocall.o: $(DEPEND)  autocall.F90
	$(F90) $(FFLAGS)  -c autocall.F90

dipolesub.o: $(DEPEND) dipolesub.F90
	$(F90) $(FFLAGS)  -c dipolesub.F90

dipolecall.o: $(DEPEND)  dipolecall.F90
	$(F90) $(FFLAGS)  -c dipolecall.F90


blocklanczos.o: blocklanczos.F90
	$(F90) $(FFLAGS) -c blocklanczos.F90

eigen.o: eigen.F90  Definitions.ALL
	$(F90) $(FFLAGS) -c eigen.F90 

utils.o: $(DEPEND) utils.F90
	$(F90) $(FFLAGS) -c utils.F90 

getparams.o: Definitions.ALL Definitions.INC parameters.o main_modules.o getparams.F90 
	$(F90) $(FFLAGS) $(SAFEFLAG) $(PGFFLAG) -c getparams.F90 

electronflux.o: $(DEPEND) electronflux.F90  biortho.o
	$(F90) $(FFLAGS) -c electronflux.F90

projeflux.o: $(DEPEND) projeflux.F90 biortho.o
	$(F90) $(FFLAGS) -c projeflux.F90




odex.o: odex.f 
	$(FC) -O2 -c odex.f
#	$(FC) $(F77FLAGS) -c odex.f

expokit.o: expokit.f 
	$(FC) $(F77FLAGS) -c expokit.f
#2014 try again	$(FC) -O0 -c expokit.f


rkf45.o: rkf45.f 
	$(FC) $(F77FLAGS) -c rkf45.f

opfunctions.o: opfunctions.f 
	$(FC) $(F77FLAGS) -c opfunctions.f

jacobi.o: jacobi.f 
	$(FC) -O0 -c jacobi.f
#	$(FC) $(F77FLAGS) -c jacobi.f

gaussq.o: gaussq.f 
	$(FC) $(F77FLAGS) -c gaussq.f

arg.o: arg.c 
	$(CC) $(PGFFLAG) -c arg.c



.F90.o: 
	$(F90) $(FFLAGS) -c $<




