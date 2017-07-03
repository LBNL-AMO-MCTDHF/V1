
######                                 ######
######     LBNL-AMO-MCTDHF MAKEFILE    ######
######                                 ######

include Makefile.header

include Name.Txt

HEPROJECT=HEPROJECT/heprojectar.a
H2PROJECT=H2PROJECT/h2projectar.a
SINCPROJECT=SINCDVR/sincprojectar.a
SINC1DPROJECT=SINC1D/sincprojectar.a

include Name.Txt

include DFFTPACK/Dfftfiles.Inc

bessrcso = bessel/d1mach.o bessel/d9lgmc.o bessel/dasyjy.o bessel/dbesj.o bessel/dcsevl.o bessel/dgamlm.o bessel/dgamln.o bessel/dgamma.o bessel/djairy.o bessel/dlngam.o bessel/fdump.o bessel/i1mach.o bessel/initds.o bessel/j4save.o bessel/xercnt.o bessel/xerhlt.o bessel/xermsg.o bessel/xerprn.o bessel/xersve.o bessel/xgetua.o bessel/drc3jj.o

SRCS= parameters.o  gaussq.o  eigen.o jacobi.o  spin.o rkf45.o   spfs.o  main_modules.o spinwalks.o ovlsub.o utils.o odex.o   getparams.o MPI.o   actions.o   readactions.o povactions.o saveactions.o natprojaction.o orbvectoractions.o   expokit.o  expo_driver.o quad.o dipolesub.o  proputils.o pulse.o autosub.o  derivs.o prop.o matel.o  projeflux.o loadstuff.o configload.o blocklanczos.o psistats.o second_derivs.o keprojector.o electronflux.o walks.o walkmult.o  newconfig.o configstuff.o biortho.o    dfconstrain.o denmat.o mean.o driving.o fock.o all_dgmres.o $(ARGO)  

DEPEND= Definitions.INC Definitions.ALL Makefile.header Makefile MPI.o main_modules.o parameters.o

default: $(MYDEFAULT)

clean:
	rm *.o *.mod *~  $(NAME)_atom $(NAME)_diatom $(NAME)_sinc $(NAME)_sinc1d *genmod*

##double colon
#linkme::
#	./linkme $(NAME)

STUFFS= mctdhf.o $(SRCS) $(bessrcso) $(DFFTFILES)

$(NAME)_diatom:   $(STUFFS)  $(H2PROJECT) 
	$(F90) $(LOADFLAGS) -o $(NAME)_diatom $(STUFFS) $(H2PROJECT)  $(LAPACK)

$(NAME)_atom:   $(STUFFS)  $(HEPROJECT)
	$(F90) $(LOADFLAGS) -o $(NAME)_atom  $(STUFFS) $(HEPROJECT)  $(LAPACK)

$(NAME)_sinc:   $(STUFFS)  $(SINCPROJECT)
	$(F90) $(LOADFLAGS) -o $(NAME)_sinc  $(STUFFS) $(SINCPROJECT)  $(LAPACK)

$(NAME)_sinc1d:   $(STUFFS)  $(SINC1DPROJECT)
	$(F90) $(LOADFLAGS) -o $(NAME)_sinc1d  $(STUFFS) $(SINC1DPROJECT)  $(LAPACK)

all_dgmres.o: all_dgmres.F90
	$(F90) $(FFLAGS) $(MPIFLAG) -c all_dgmres.F90 

parameters.o: Definitions.INC Definitions.ALL Makefile.header Makefile parameters.F90 
	$(F90) $(FFLAGS) $(MPIFLAG) -c parameters.F90

MPI.o:  Definitions.INC Definitions.ALL Makefile.header Makefile MPI.F90 parameters.o 
	$(F90) $(FFLAGS) $(MPIFLAG) -c MPI.F90

main_modules.o: Definitions.INC Definitions.ALL Makefile.header Makefile MPI.o parameters.o main_modules.F90
	$(F90) $(FFLAGS) $(MPIFLAG) $(PGFFLAG) -c main_modules.F90

biortho.o: $(DEPEND) biortho.F90 newconfig.o spin.o utils.o expokit.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c biortho.F90 

second_derivs.o: $(DEPEND) second_derivs.F90 derivs.o proputils.o spin.o expo_driver.o spfs.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c second_derivs.F90 

keprojector.o: $(DEPEND) keprojector.F90 walkmult.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c keprojector.F90 

loadstuff.o: $(DEPEND) loadstuff.F90 configload.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c loadstuff.F90 

configload.o: $(DEPEND) configload.F90 newconfig.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c configload.F90 

psistats.o: $(DEPEND) psistats.F90 walkmult.o biortho.o autosub.o denmat.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c psistats.F90 

dfconstrain.o: $(DEPEND) dfconstrain.F90 walkmult.o spin.o utils.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c dfconstrain.F90 

expo_driver.o: $(DEPEND) expo_driver.F90 derivs.o proputils.o spin.o pulse.o expokit.o spfs.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c expo_driver.F90 

quad.o: $(DEPEND) quad.F90 blocklanczos.o expo_driver.o derivs.o walkmult.o spin.o matel.o utils.o spfs.o all_dgmres.o 
	$(F90) $(FFLAGS) $(MPIFLAG) -c quad.F90 

BCG_driver.o: $(DEPEND) implicit.o BCG_driver.F90
	$(F90) $(FFLAGS) $(MPIFLAG) -c BCG_driver.F90 

mctdhf.o: $(DEPEND) mctdhf.F90 newconfig.o configload.o quad.o configstuff.o utils.o denmat.o proputils.o prop.o fock.o loadstuff.o mean.o psistats.o spfs.o spin.o spinwalks.o walks.o getparams.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c mctdhf.F90 

actions.o: $(DEPEND) actions.F90 autosub.o denmat.o dfconstrain.o proputils.o dipolesub.o electronflux.o projeflux.o keprojector.o ovlsub.o psistats.o
	$(F90) $(FFLAGS) $(MPIFLAG) $(MATHEMAT) $(POVRAY) $(NOCFLAG) -c actions.F90 

readactions.o: $(DEPEND) readactions.F90 
	$(F90) $(FFLAGS) $(MPIFLAG) $(MATHEMAT) $(POVRAY) $(NOCFLAG) -c readactions.F90 

saveactions.o: $(DEPEND) saveactions.F90 
	$(F90) $(FFLAGS) $(MPIFLAG) $(MATHEMAT) $(POVRAY) $(NOCFLAG) -c saveactions.F90 

natprojaction.o: $(DEPEND) natprojaction.F90 walkmult.o denmat.o
	$(F90) $(FFLAGS) $(MPIFLAG) $(MATHEMAT) $(POVRAY) $(NOCFLAG) -c natprojaction.F90 

orbvectoractions.o: $(DEPEND) orbvectoractions.F90
	$(F90) $(FFLAGS) $(MPIFLAG) $(MATHEMAT) $(POVRAY) $(NOCFLAG) -c orbvectoractions.F90 

povactions.o: $(DEPEND) povactions.F90 
	$(F90) $(FFLAGS) $(MPIFLAG) $(MATHEMAT) $(POVRAY) $(NOCFLAG) -c povactions.F90 

prop.o: $(DEPEND) prop.F90 derivs.o walkmult.o spin.o pulse.o odex.o quad.o configstuff.o second_derivs.o denmat.o dfconstrain.o proputils.o driving.o fock.o loadstuff.o mean.o spfs.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c prop.F90 

proputils.o: $(DEPEND) proputils.F90 pulse.o dfconstrain.o denmat.o mean.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c proputils.F90 

pulse.o: $(DEPEND) pulse.F90  
	$(F90) $(FFLAGS) $(MPIFLAG) -c pulse.F90 

derivs.o: $(DEPEND) derivs.F90 proputils.o walkmult.o utils.o pulse.o spin.o spfs.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c derivs.F90 

mean.o: $(DEPEND) mean.F90 MPI.o denmat.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c mean.F90 

matel.o: $(DEPEND)  matel.F90 spin.o proputils.o pulse.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c matel.F90 

spin.o: $(DEPEND)  spin.F90 newconfig.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c spin.F90 

spfs.o: $(DEPEND) spfs.F90 utils.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c spfs.F90

denmat.o: $(DEPEND) denmat.F90 biortho.o utils.o eigen.o spfs.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c denmat.F90

fock.o: $(DEPEND) fock.F90 biortho.o walkmult.o matel.o utils.o newconfig.o eigen.o spfs.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c fock.F90

driving.o: $(DEPEND) driving.F90 biortho.o walkmult.o denmat.o proputils.o matel.o pulse.o dipolesub.o mean.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c driving.F90

newconfig.o: $(DEPEND) newconfig.F90 
	$(F90) $(FFLAGS) $(MPIFLAG) -c newconfig.F90

configstuff.o: $(DEPEND) configstuff.F90 spin.o utils.o matel.o expokit.o expo_driver.o newconfig.o blocklanczos.o eigen.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c configstuff.F90

walks.o: $(DEPEND) walks.F90 newconfig.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c walks.F90

walkmult.o: $(DEPEND) walkmult.F90 pulse.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c walkmult.F90

spinwalks.o: $(DEPEND) spinwalks.F90 newconfig.o spin.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c spinwalks.F90

ovlsub.o: $(DEPEND) ovlsub.F90 configload.o biortho.o autosub.o utils.o loadstuff.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c ovlsub.F90

autosub.o: $(DEPEND) autosub.F90 biortho.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c autosub.F90

dipolesub.o: $(DEPEND) dipolesub.F90 biortho.o walkmult.o pulse.o utils.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c dipolesub.F90

blocklanczos.o: $(DEPEND) blocklanczos.F90 walkmult.o spin.o eigen.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c blocklanczos.F90

eigen.o: eigen.F90  Definitions.ALL
	$(F90) $(FFLAGS) $(MPIFLAG) -c eigen.F90 

utils.o: $(DEPEND) utils.F90
	$(F90) $(FFLAGS) $(MPIFLAG) -c utils.F90 

getparams.o: Definitions.ALL Definitions.INC parameters.o main_modules.o getparams.F90 pulse.o actions.o utils.o
	$(F90) $(FFLAGS) $(MPIFLAG) $(SAFEFLAG) $(PGFFLAG) -c getparams.F90 

electronflux.o: $(DEPEND) electronflux.F90 biortho.o walkmult.o proputils.o matel.o pulse.o utils.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c electronflux.F90

projeflux.o: $(DEPEND) projeflux.F90 biortho.o newconfig.o pulse.o utils.o configload.o electronflux.o loadstuff.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c projeflux.F90

expokit.o: expokit.F90
	$(F90) $(FFLAGS) -c expokit.F90

#odex.o: odex.f
#	$(FC) $(F77FLAGS) -c odex.f

odex.o: odex.F90
	$(F90) $(FFLAGS) -c odex.F90

####  F77 #####

rkf45.o: rkf45.f 
	$(FC) $(F77FLAGS) -c rkf45.f

jacobi.o: jacobi.f 
	$(FC) -O0 -c jacobi.f
#	$(FC) $(F77FLAGS) -c jacobi.f

gaussq.o: gaussq.f 
	$(FC) $(F77FLAGS) -c gaussq.f

arg.o: arg.c 
	$(CC) $(PGFFLAG) -c arg.c

