
!! ALL ONE MODULE

!!! MAIN SUBROUTINE FOR READING NAMELIST INPUT AND COMMAND LINE OPTIONS
  
#include "Definitions.INC"

module getparammod
contains 

subroutine getinpfile()
  use parameters
  use mpimod
  implicit none

  integer :: nargs, getlen, i, len
#ifdef PGFFLAG
  integer :: myiargc
#endif
  character (len=SLN) :: buffer

#ifdef PGFFLAG
  nargs=myiargc()
#else
  nargs=iargc()
#endif
  do i=1,nargs
     buffer=nullbuff;     call getarg(i,buffer);     len=getlen(buffer)
     if (buffer(1:4) .eq. 'Inp=') then
        inpfile=nullbuff;        inpfile(1:len-4)=buffer(5:len)
        OFLWR "Inpfile is ", inpfile(1:len-4+1); CFL
     endif
  enddo
end subroutine getinpfile

subroutine getparams()
  use parameters
  use denreg_parameters
  use bio_parameters
  use constraint_parameters
  use lanparameters
  use output_parameters
  use mpimod
  use orblabelmod
  use pulse_parameters   !! numpulses
  use actionlistmod      !! ftwindowpower,fttriwindow
  implicit none
  integer :: nargs, getlen, i, len,  ishell, ispf,j, myiostat, iiflag,needpulse,ipulse
#ifdef PGFFLAG
  integer :: myiargc
#endif
  character (len=SLN) :: buffer
  integer :: shelltop(100)=-1            !! greater than zero: not full CI.  Number of orbitals in core shell.  Must be even.
  integer :: avectorhole(1000)=-1001
  integer :: avectorexcitefrom(1000)=-1001
  integer :: avectorexciteto(1000)=-1001
  real*8 :: tempreal,mymax

!! DUMMIES
  integer :: restrictms=0,  dfrestrictflag=0, allspinproject=1

  NAMELIST/parinp/  noftflag, biodim,biotol, rdenflag,cdenflag, notiming, littlesteps, &
       expotol, eground,  ceground, maxexpodim, numloadfrozen, numholecombo, numholes, excitations, &
       excitecombos, jacsymflag, jacprojorth,jacgmatthird,  fluxoptype, timefac, threshflag, &
       timefacforce, avectoroutfile, spfoutfile,  autopermthresh, messamount, numshells,   &
       lanthresh, lanczosorder,  lioreg, &   !! rcond
       autonormthresh,  saveflag, save_every, &

!!$ IMPLEMENT ME (DEPRECATE fluxinterval as namelist input)      fluxtimestep,   

       mrestrictflag, mrestrictval, autotimestep, &
       nosparseforce,  allspinproject,  numfluxfiles,  verletnum, &
       constraintflag, dfrestrictflag, improvedrelaxflag, mcscfnum,  improvednatflag, avectorfile, spffile, &
       quadprecon,  improvedquadflag, quadtol, improvedfockflag, &
       plotxyrange,plotrange,plotskip,pm3d, plotterm, plotmodulus,plotpause, numfluxcurves,  & 
       plotres, nspf ,  spfrestrictflag,  spfmvals, spfugvals, spfugrestrict, ugrestrictflag, ugrestrictval, &
       restrictflag,  restrictms,   loadspfflag,  loadavectorflag,  avectorhole, &
       par_timestep,  stopthresh ,  cmf_flag,  intopt, timedepexpect,  avector_flag, &
       numelec,  relerr,  myrelerr,  spf_flag,  denreg,  timingout, tdflag, finaltime, actions, numactions, &
       messflag,  sparseconfigflag,  aorder, maxaorder, aerror, shelltop, numexcite, povres, povrange,&
       numpovranges, povsparse,  povmult, plotnum,lancheckstep,  plotview1, plotview2, &
       computeFlux, FluxInterval, FluxSkipMult, &
       numfrozen, nucfluxopt, natplotbin, spfplotbin, denplotbin, denprojplotbin, &
       natprojplotbin, rnatplotbin, dendatfile, denrotfile, rdendatfile,   avectorexcitefrom, avectorexciteto,&
       numovlfiles, ovlavectorfiles, ovlspffiles, outovl,fluxmofile,fluxafile, spifile, astoptol,  &
       xdipfile,ydipfile,zdipfile,   xtworkfile,ytworkfile,ztworkfile,  xytworkfile,yztworkfile,zxtworkfile,&
       xdftfile,ydftfile,zdftfile,   xydftfile,xzdftfile,yxdftfile,yzdftfile,zxdftfile,zydftfile, &
       xoworkfile,yoworkfile,zoworkfile,&
       xyoworkfile,xzoworkfile,yxoworkfile,yzoworkfile,zxoworkfile,zyoworkfile,&
       xophotonfile, yophotonfile,zophotonfile,&
       xyophotonfile,xzophotonfile,yxophotonfile,yzophotonfile,zxophotonfile, zyophotonfile,&
       act21circ,       ftwindowpower, ftdiff,  &
       diptime, fluxafile2,fluxmofile2, minocc,maxocc, corrdatfile,corrftfile,numavectorfiles,projfluxfile, &
       expodim,timingdir, hanningflag, numspffiles, condamp,conway, &
       mrestrictmin,mrestrictmax,lntol,invtol,psistatsfile, psistatfreq, &
       parorbsplit,maxbiodim, nkeproj,keprojminenergy,keprojenergystep,keprojminrad,keprojmaxrad, &
       debugflag, drivingflag,drivingproportion, drivingmethod, eigprintflag, &
       avecloadskip,nonsparsepropmode,sparseopt,lanprintflag,dipmodtime,conprop,&
       orbcompact,spin_restrictval,mshift,numskiporbs,orbskip,debugfac,denmatfciflag,&
       walkwriteflag,iprintconfiglist,timestepfac,max_timestep,expostepfac, maxquadnorm,quadstarttime,&
       reinterp_orbflag,spf_gridshift,load_avector_product,projspifile,readfullvector,walksinturn,&
       turnbatchsize,energyshift, pulseft_estep, finalstatsfile, fluxtsumfile, projfluxtsumfile,&
       sparsedfflag,sparseprime,sparsesummaflag, par_consplit, fttriwindow,&
       pulsewindowtoo,redobra,dipolesumstart,dipolesumend,outmatel,numcatfiles,&
       catspffiles,catavectorfiles,aquadstarttime,quadorthflag,normboflag,logbranch,nzflag,&
       shuffle_dfwalktype,maxdgdim, messavec, messaamount,holeflag, angularflag, angprojspifile,&
       prepropflag, step_flag, postpropflag, scalarflag, angprojfluxtsumfile, &
       catfacs, flux_subtract, jacsymquad, exact_exchange, jacquaddir, tentmode, followflag, &
       exchange_mode, gaugefluxflag, strongcatflag, strongcatspffiles, strongcatavectorfiles, &
       catenergies, nonuc_checkflag, autostart, auto_subtract

  OFL
  write(mpifileptr, *)
  write(mpifileptr, *) " *************************  COMMAND LINE OPTIONS  ***************************"
  write(mpifileptr, *) 
  CFL
#ifdef PGFFLAG
  nargs=myiargc()
#else
  nargs=iargc()
#endif

!! Defaults not set in declarations

  numfluxcurves=0;  numfluxcurves(1)=1

#ifdef REALGO
  conway=0
#else
  conway=3
#endif

!! F.T. windowing function, defaults for each action
!! use linear on [0:t], cosine-squared on [-t:t]

  fttriwindow(:)=1
  ftwindowpower(:)=1

  fttriwindow(1)=0
  ftwindowpower(1)=2
  fttriwindow(16)=0
  ftwindowpower(16)=2
  fttriwindow(17)=0
  ftwindowpower(17)=2

  open(971,file=inpfile, status="old", iostat=myiostat)
  if (myiostat/=0) then
     OFLWR "No Input.Inp found for parinp, iostat=",myiostat; CFL
  else
     OFLWR "Reading ",inpfile; CFL
     read(971,nml=parinp,iostat=myiostat)
     call checkiostat(myiostat," reading namelist PARINP")
     restrict_ms=restrictms
     all_spinproject=allspinproject
     df_restrictflag=dfrestrictflag

!! input-dependent defaults

!!$ prepropflag.eq.1 default
!!$     if (constraintflag.ne.0) then
!!$        prepropflag=1
!!$     endif

     if (numfrozen.gt.0) then 
        jacsymquad=0
     endif

     jacquaddir = -1
     if ((jacsymflag.eq.0.and.jacsymquad.eq.0).or.(numfrozen.gt.0.and.exact_exchange.eq.0)) then
        jacquaddir = 1
     endif

     if (improvedrelaxflag.ne.0) then
        maxexpodim=max(300,maxexpodim)
        expodim=max(40,expodim)
     endif

     if (intopt.eq.4) then
        expotol=min(1d-9,expotol)
     endif
     if (intopt.eq.1) then
        relerr=1d-8
        myrelerr=1d-8
     endif

     if (improvedrelaxflag.ne.0) then
        expotol=min(expotol,1d-9)
     endif

     if (spin_restrictval.lt.abs(restrictms)) then
        spin_restrictval=abs(restrictms)
     endif

!!$ denreg now 1d-10 default...
!!$     if (improvedrelaxflag.ne.0) then    !! not good.  reprogram this whole thing later.  Defaults here should
!!$        denreg=1.d-9                     !!   go after a first parinp AND command line argument parse; then repeat
!!$     endif

     if (improvedrelaxflag.ne.0) then
        aorder=max(300,aorder)
        maxaorder=max(aorder,maxaorder)
     endif

     if (spfrestrictflag.eq.0) then
        orbcompact=0
     endif

     close(971)     
     open(971,file=inpfile, status="old")
     read(971,nml=parinp,iostat=myiostat)
     call checkiostat(myiostat," reading namelist PARINP")
  endif
  close(971)

  
  !!   ************************************************************************************
  !!
  !!    Coord-dependent namelist input and command line options may be overridden by
  !!    options of the main program
  !!
  !!   ************************************************************************************


  !!   NOW TAKE MCTDHF NAMELIST AND COMMAND LINE INPUT

  call openfile()

  !! NOTE THAT LBIG IS NOT A COMMAND LINE OPTION; IS IN LOOP IN GETH2OPTS.  So as of now it is set.  change later.

  do i=1,nargs
     buffer=nullbuff;     call getarg(i,buffer);     len=getlen(buffer)
     myiostat=0
     if (buffer(1:5) .eq. 'Nspf=') then
        read(buffer(6:len),*,iostat=myiostat) nspf;   
        write(mpifileptr,*) "Nspf set to  ", nspf, " by command line option."
     endif

     if (buffer(1:5) .eq. 'Holes') then
        holeflag=1
        write(mpifileptr, *) "indexing holes not electrons by command line option"
     endif

     if (buffer(1:8) .eq. 'NucFlux=') then
        read(buffer(9:len),*,iostat=myiostat) nucfluxopt
        write(mpifileptr, *) "nucfluxopt set to ", nucfluxopt, " by command line option."
     endif

     if (buffer(1:9) .eq. 'NoTiming=') then
        read(buffer(10:len),*,iostat=myiostat) notiming
        write(mpifileptr, *) "notiming variable set to ",notiming," by command line input."
     endif
     if (buffer(1:10) .eq. 'TimingDir=') then
        read(buffer(11:len),*,iostat=myiostat) timingdir
        write(mpifileptr, *) "timing directory set to ",timingdir," by command line input."
     endif
     if (buffer(1:7) .eq. 'Timing=') then
        read(buffer(8:len),*,iostat=myiostat) j
        notiming=2-j
        write(mpifileptr, *) "notiming variable set to ",notiming," by command line input."
     endif
     if (buffer(1:4) .eq. 'Rel=') then
        read(buffer(5:len),*,iostat=myiostat) relerr
        write(mpifileptr, *) "Relative error for spf prop set to ", relerr, " by command line option."
     endif
     if (buffer(1:9) .eq. 'FluxSkip=') then
        read(buffer(10:len),*,iostat=myiostat) fluxskipmult
        write(mpifileptr, *) "Fluxskipmult set to ", fluxskipmult, " by command line option."
     endif
     if (buffer(1:6) .eq. 'Myrel=') then
        read(buffer(7:len),*,iostat=myiostat) myrelerr
        write(mpifileptr, *) "Absolute error (myrelerr) for spf prop set to ", myrelerr, " by command line option."
     endif
     if (buffer(1:9) .eq. 'PovRange=') then
        read(buffer(10:len),*,iostat=myiostat) tempreal
        write(mpifileptr, *) "Povrange set to ", tempreal
        povrange=tempreal  !! all of them
     endif
     if (buffer(1:2) .eq. 'M=') then
        mrestrictflag=1
        read(buffer(3:len),*,iostat=myiostat) mrestrictval
        write(mpifileptr, *) "Restricting configs to ", mrestrictval, " by command line option."
     endif
     if (buffer(1:5) .eq. 'Debug') then
        if (.not.buffer(1:6) .eq. 'Debug=') then
           WRFL "Please specify debug flag option N with command line Debug=N not just Debug"; CFLST
        endif
        read(buffer(7:len),*,iostat=myiostat) debugflag
        write(mpifileptr, *) "Debugflag set to ",debugflag," by command line option"
     endif
     if (buffer(1:5) .eq. 'Walks') then
        walkwriteflag=1
        write(mpifileptr, *) "Walks will be written by command line option"
     endif
     if (buffer(1:7) .eq. 'NoWalks') then
        walkwriteflag=0
        write(mpifileptr, *) "Walks will not be written by command line option"
     endif
     if (buffer(1:2) .eq. 'A=') then
        avectorfile(:)=nullbuff
        avectorfile(:)(1:len-2)=buffer(3:len)
        write(mpifileptr, *) "Avector file is ", avectorfile(1)(1:len-2+1)
        write(mpifileptr, *) "    Loadavectorflag turned on."
        loadavectorflag=1
     endif
     if (buffer(1:4) .eq. 'Spf=') then
        spffile=nullbuff
        spffile(:)(1:len-4)=buffer(5:len)
        write(mpifileptr, *) "Spf file is ", spffile(1)(1:len-4+1)
        write(mpifileptr, *) "    Loadspfflag turned on."
        loadspfflag=1
     endif
     if (buffer(1:5) .eq. 'Pulse') then
        write(mpifileptr, *) "Turning pulse on by command line option."
        tdflag=1
     endif
     if (buffer(1:6) .eq. 'SaveOn') then
        write(mpifileptr, *) "Saving wave function by command line option."
        saveflag=1
     endif
     if (buffer(1:7) .eq. 'SaveOff') then
        write(mpifileptr, *) "NOT Saving wave function by command line option."
        saveflag=0
     endif
#ifndef REALGO
     if (buffer(1:4) .eq. 'Prop') then
        write(mpifileptr, *) "  Forcing propagation in real time by command line option."
        improvedrelaxflag=0
     endif
#endif
     if (buffer(1:6) .eq. 'Relax=') then
        read(buffer(7:len),*,iostat=myiostat) improvedrelaxflag
        write(mpifileptr, *) "  Forcing improved relaxation by command line option.  Relaxing to state ", improvedrelaxflag
     else
        if (buffer(1:5) .eq. 'Relax') then
           write(mpifileptr, *) "  Forcing improved relaxation to ground state by command line option."
           improvedrelaxflag=max(1,improvedrelaxflag)
        endif
     endif

! so you can input numactions=0 but specify the actions you would want in the input file, and then turn them on with:

     if (buffer(1:5) .eq. 'Noact') then
        numactions=0        
     endif
     if (buffer(1:6) .eq. 'Allact') then
        do while (actions(numactions+1).ne.0)
           numactions=numactions+1
        enddo
        write(mpifileptr,*) " Setting all specified actions on: they are ", actions(1:numactions)
     endif
     if (buffer(1:4) .eq. 'Act=') then
        numactions=numactions+1;        read(buffer(5:len),*,iostat=myiostat) actions(numactions)
     endif
     if (buffer(1:4) .eq. 'Mess') then
        messflag=1;        messamount=1.d-3
        if (buffer(1:5) .eq. 'Mess=') then
           read(buffer(6:len),*,iostat=myiostat) messamount
        endif
     endif
     if (buffer(1:9) .eq. 'Autoperm=') then
        read(buffer(10:len),*,iostat=myiostat) autopermthresh
        write(mpifileptr,*) "Permutation threshold for autocorr set to ", autopermthresh
     endif

     if (buffer(1:3) .eq. 'VMF') then
        cmf_flag=0;        write(mpifileptr, *) "VMF set by command line option"
     endif

     if (buffer(1:5) .eq. 'Step=') then
        read(buffer(6:len),*,iostat=myiostat) par_timestep
        write(mpifileptr, *) "Timestep set to  ", par_timestep, " by command line option."
     endif

     if (buffer(1:10) .eq. 'Numfrozen=') then
        numshells=2
        read(buffer(11:len),*,iostat=myiostat) shelltop(1)
        shelltop(2)=nspf
        write(mpifileptr, *) "Numshells set to 2 with  ", shelltop(1), " in the first shell by command line option."
     endif
     if (buffer(1:10) .eq. 'Numexcite=') then
        read(buffer(11:len),*,iostat=myiostat) numexcite(1)
        write(mpifileptr, *) "Numexcite for first shell set to  ", numexcite(1), " by command line option."
     endif

     if (buffer(1:9) .eq. 'PlotSkip=') then
        read(buffer(10:len),*,iostat=myiostat) plotskip
        write(mpifileptr, *) "Plotskip set to ", plotskip, " by command line option."
     endif
     if (buffer(1:7) .eq. 'PlotXY=') then
        read(buffer(8:len),*,iostat=myiostat) plotxyrange
        write(mpifileptr, *) "Plot xy-range set to ", plotxyrange, " bohr by command line option."
     endif

     if (buffer(1:4) .eq. 'PM3D') then
        pm3d=1;        write(mpifileptr, *) "PM3D set to on."
     endif
     if (buffer(1:6) .eq. 'PlotZ=') then
        read(buffer(7:len),*,iostat=myiostat) plotrange
        write(mpifileptr, *) "Plot z-range set to ", plotrange, " bohr by command line option."
     endif
     if (buffer(1:8) .eq. 'PlotNum=') then
        read(buffer(9:len),*,iostat=myiostat) plotnum
        write(mpifileptr, *) "Plotnum set to ", plotnum, " bohr by command line option."
     endif
     if (buffer(1:8) .eq. 'PlotRes=') then
        read(buffer(9:len),*,iostat=myiostat) plotres
        write(mpifileptr, *) "Plot resolution set to ", plotres, " bohr by command line option."
     endif
     if (buffer(1:10) .eq. 'PlotPause=') then
        read(buffer(11:len),*,iostat=myiostat) plotpause
        write(mpifileptr, *) "Plotpause set to ", plotpause, " seconds by command line option."
     endif
     if (buffer(1:11) .eq. 'Stopthresh=') then
        read(buffer(12:len),*,iostat=myiostat) stopthresh
        write(mpifileptr, *) "Stopthresh set to  ", stopthresh, " by command line option."
     endif
     if (buffer(1:6) .eq. 'Sparse') then
        write(mpifileptr, *) "Sparseconfigflag turned on by command line option."
        sparseconfigflag=1
     endif
     if (buffer(1:8) .eq. 'NoSparse') then
        write(mpifileptr, *) "Sparseconfigflag turned OFF by command line option."
        sparseconfigflag=0
     endif
     if (buffer(1:7) .eq. 'Denreg=') then
        read(buffer(8:len),*,iostat=myiostat) denreg
        write(mpifileptr, *) "Denreg set by command line option to ", denreg
     endif
     if (buffer(1:3) .eq. 'GBS') then
        intopt=1
        write(mpifileptr, *) "GBS integration set by command line option."
     endif
     if (buffer(1:2) .eq. 'RK') then
        intopt=0
        write(mpifileptr, *) "RK integration set by command line option."
     endif
     if (buffer(1:8) .eq. 'Eground=') then
        read(buffer(9:len),*,iostat=myiostat) eground
        write(mpifileptr, *) "Eground for autoft set to  ", eground, " by command line option."
     endif

     call checkiostat(myiostat," command line argument "//buffer)
     
  enddo
  write(mpifileptr, *) " ****************************************************************************"     
  write(mpifileptr,*);  call closefile()
  
  if (constraintflag > 2) then
     OFLWR "Constraintflag not supported: ", constraintflag;     CFLST
  endif
  
  if (improvedrelaxflag.eq.0) then
     improvednatflag=0
     improvedfockflag=0
  endif

  if (improvedrelaxflag.ne.0.and.constraintflag.eq.1) then
     OFLWR "FOR DEN CONSTRAINT, USE IMPROVEDNATFLAG FOR RELAX, NO CONSTRAINTFLAG."; CFLST
  endif

  if ((sparseconfigflag.ne.0).and.(stopthresh.lt.lanthresh).and.(improvedquadflag.eq.0.or.improvedquadflag.eq.2)) then
     OFLWR "Enforcing lanthresh.le.stopthresh"
     lanthresh=stopthresh
     write(mpifileptr,*) "    --> lanthresh now  ", lanthresh; CFL
  endif
  if (intopt.eq.4) then
     if ((constraintflag.ne.0)) then
        OFLWR "Verlet not available with pulse nor constraint yet."; CFLST
     endif
  endif
  if ((intopt.eq.3).or.(intopt.eq.4)) then  
     OFLWR "Enforcing CMF defaults for Verlet or EXPO."; CFL;     cmf_flag=1
  endif

  if (timefacforce.eq.0) then
     timefac=(0.0d0, -1.0d0)
  endif
  if (improvedrelaxflag.ne.0) then
     threshflag=1
  endif
  if (threshflag.ne.0) then
     if (timefacforce.eq.0) then
        timefac=(-1.0d0, 0.0d0)
     endif
  endif

  if (numshells.lt.1) then
     OFLWR "Shell error ", numshells; CFLST
  endif

  allshelltop(0)=0;  allshelltop(numshells)=nspf
  allshelltop(1:numshells-1)=shelltop(1:numshells-1)
  do ishell=1,numshells
     if (allshelltop(ishell).gt.nspf) then
        allshelltop(ishell)=nspf
     endif
  enddo
  do ishell=numshells,2,-1
     if (allshelltop(ishell).le.allshelltop(ishell-1)) then
        numshells=numshells-1
        allshelltop(ishell:numshells)=allshelltop(ishell+1:numshells+1)
        numexcite(ishell:numshells)=numexcite(ishell+1:numshells+1)
     endif
  enddo

  do i=1,numshells
     if (allshelltop(i).le.allshelltop(i-1)) then
        OFLWR "allShell error ", allshelltop(i), allshelltop(i-1); CFLST
     endif
  enddo

  do i=2,numshells
     if (numexcite(i).lt.numexcite(i-1)) then
        OFLWR "numexcite error ", numexcite(i), numexcite(i-1); CFLST
     endif
  enddo

  !! define shells

  ishell=1
  do ispf=1,nspf
     shells(ispf)=ishell
     if (ispf.eq.allshelltop(ishell)) then
        ishell=ishell+1
     endif
  enddo

  call get_skipflag_from_actions(skipflag)

  needpulse=0
  if (tdflag.ne.0) then
     needpulse=1
  else 
     if (noftflag.eq.0) then
        do j=1,numactions
           if ( (actions(j).eq.14).or.&
                (actions(j).eq.16).or.&
                (actions(j).eq.17)) then
              needpulse=1
           endif
        enddo
     endif
  endif

  if (needpulse.ne.0) then
     call getpulse(0)
     if (dipolesumend.le.0d0) then
        do ipulse=1,numpulses
           if (pulsetype(ipulse).eq.1) then
              mymax=omega(ipulse)*3d0
           else
              mymax=omega2(ipulse)*3d0
           endif
           if (dipolesumend.lt.mymax) then
              dipolesumend=mymax
           endif
        enddo
     endif
  endif

  do i=1,nargs
     buffer=nullbuff
     myiostat=0
     !!     call mygetarg(i,buffer)
     call getarg(i,buffer)
     len=getlen(buffer)
     if (buffer(1:2) .eq. 'T=') then
        read(buffer(3:len),*,iostat=myiostat) finaltime 
        OFLWR "Finaltime set by command line option to ", finaltime
        numpropsteps=floor((finaltime+0.0000000001d0)/par_timestep) +1
        finaltime=numpropsteps*par_timestep
        write(mpifileptr, *) "     numpropsteps now   ", numpropsteps
        write(mpifileptr, *) "     finaltime    now   ", finaltime;        call closefile()
     endif
     call checkiostat(myiostat,"Reading finaltime command line input "//buffer)
  enddo

  numpropsteps=floor((finaltime+0.0000000001d0)/par_timestep)

  autosteps=floor(max(1.d0,autotimestep/par_timestep)); 

  autotimestep = autosteps * par_timestep   !! OOPS forgot 6-16

  autosize=numpropsteps/autosteps+1

  if (autostart.lt.0d0) then
     autostart = (ceiling(lastfinish/par_timestep/autosteps))*par_timestep*autosteps
     OFLWR "   autostart set to ",autostart
  endif

!!$ IMPLEMENT ME (DEPRECATE fluxinterval as namelist input)   
!!$ fluxsteps=floor(max(1.d0,fluxtimestep/par_timestep));  fluxtimestep=par_timestep*fluxsteps

  do i=1,nargs
     buffer=nullbuff;     call getarg(i,buffer)
     myiostat=0
     len=getlen(buffer)
     if (buffer(1:12) .eq. 'PlotModulus=') then
        read(buffer(13:len),*,iostat=myiostat) plotmodulus
        write(mpifileptr, *) "Plotmodulus set to ", plotmodulus, " steps by command line option."
     endif
     call checkiostat(myiostat," command line argument plotmodulus "//buffer)
  enddo

  if (holeflag.eq.0) then
     numpart=numelec
  else
     numpart=2*nspf-numelec
  endif
  num2part=2*numpart

  if (stopthresh.lt.1.d-12) then
     OFLWR "Error, stopthresh cannot be less than 1d-12"; CFLST  !! then would send hgram 1d-14
  endif

  if (constraintflag.eq.1.and.improvedrelaxflag.ne.0) then 
     OFLWR " Removing denmat constraint for relaxation. Not allowed."; CFL
     constraintflag=0
  endif
  if (spfrestrictflag.eq.0) then
     mrestrictflag=0
  endif
  if (ceground.eq.(0d0,0d0)) then
     ceground=eground
  endif

!! 121912
!! if numholes or numexcite is not set, define from avectorhole etc. input for backwards compatibility

  if (numholes.eq.0) then
     if (numholecombo.ne.1) then
        OFLWR "If setting numholecombo, set numholes."; CFLST
     endif
     do while (avectorhole(numholes+1).ne.-1001)
        numholes=numholes+1
     enddo
     if (numholes.ne.0) then
        OFLWR "You did not set numholes; I am setting it to ", numholes; CFL
     endif
  endif
  if (excitations.eq.0) then
     if (excitecombos.ne.1) then
        OFLWR "If setting excitecombos, set excitations."; CFLST
     endif
     do while (avectorexcitefrom(excitations+1).ne.-1001)
        print *, "blah ",excitations,avectorexcitefrom(excitations+1)
        excitations=excitations+1
     enddo
     if (excitations.ne.0) then
        OFLWR "You did not set excitations; I am setting it to ", excitations; CFL
     endif
  endif
  if (numholes.gt.0) then
     excitations=0
     allocate(myavectorhole(numholes,numholecombo,mcscfnum))
     myavectorhole=RESHAPE(avectorhole(1:numholecombo*numholes*mcscfnum),(/ numholes,numholecombo,mcscfnum/))
  else
     allocate(myavectorhole(1,1,1))
  endif
  if (excitations.gt.0) then
     numholes=0
     allocate(myavectorexcitefrom(excitations,excitecombos,mcscfnum))
     myavectorexcitefrom=RESHAPE(avectorexcitefrom(1:excitations*excitecombos*mcscfnum),(/excitations,excitecombos,mcscfnum/))
     allocate(myavectorexciteto(excitations,excitecombos,mcscfnum))
     myavectorexciteto=RESHAPE(avectorexciteto(1:excitations*excitecombos*mcscfnum),(/excitations,excitecombos,mcscfnum/))
  else
     allocate(myavectorexcitefrom(1,1,1), myavectorexciteto(1,1,1))
  endif

  if (nprocs.eq.1) then
     par_consplit=0
     parorbsplit=0
  endif

  if (sparseconfigflag.eq.0) then
     sparseopt=0
     par_consplit=0
  endif

#ifndef MPIFLAG
  if (par_consplit.ne.0.or.parorbsplit.eq.3) then
     OFLWR "Error, this is not MPI chmctdhf; cannot use par_consplit.ne.0 or parorbsplit.eq.3"; CFLST
  endif
#endif

  if (numavectorfiles.gt.MXF.or.numspffiles.gt.MXF) then
     OFLWR "PROGRAMMER REDIM M X F",numavectorfiles,numspffiles,MXF; CFLST
  endif

  call openfile()
  write(mpifileptr, *)
  write(mpifileptr, *) " ****************************************************************************"     
  write(mpifileptr, *) "*****************************************************************************"
  write(mpifileptr, *) 
  write(mpifileptr, *) "Number of states in propagation=", mcscfnum
  write(mpifileptr, *)
  write(mpifileptr, *) " Parameters: electronic"
  write(mpifileptr, *)
  write(mpifileptr,'(A40,I5,2F10.5)') "Number of electrons ", numelec
  write(mpifileptr,'(A40,3I5)') "Nuclear KE flag (nonuc_checkflag):     ",  nonuc_checkflag
  write(mpifileptr, *)
  
  call printmyopts()
     
  write(mpifileptr, *)
  write(mpifileptr, *) "************************   Parameters: config/spf    ************************   "
  write(mpifileptr, *)
  
  if (numshells.eq.1) then
     write(mpifileptr, *) "Doing full CI: numshells=1.  Shells:"
  else
     write(mpifileptr, *) "Shells:"
  endif
  do ishell=1,numshells
     if (ishell.eq.numshells) then
        write(mpifileptr,*) "Shell ", ishell
     else
        write(mpifileptr,*) "Shell ", ishell,   "Excitation level: ", numexcite(ishell)
     endif
     write(mpifileptr,'(12A6)')  (orblabel(i),i=allshelltop(ishell-1)*2+1,allshelltop(ishell)*2)
  enddo
  if (df_restrictflag.gt.0) then
     write(mpifileptr,*) " DF restrictflag is on for constraintflag=2 or other purposes = ",df_restrictflag
  endif
  write(mpifileptr, *) 
  write(mpifileptr,'(A30,I5)')    "   Number of unfrozen spfs:  ", nspf
  write(mpifileptr,'(A30,I5)')    "   Number of frozen spfs:    ", numfrozen
  write(mpifileptr,'(A30,100I4)') "   Spfs start in m=  ", spfmvals(1:nspf)
  write(mpifileptr, *) 
  if (spfrestrictflag.eq.1) then
     write(mpifileptr, *) "Spfs will be restricted to their original m-values."
  endif
  if (spfugrestrict.eq.1) then
     write(mpifileptr, *) "Spfs will be restricted to their original parity values."
  endif
  if (mrestrictflag==1) then
     write(mpifileptr, *) "Configurations will be restricted to total M= ", mrestrictval
  endif
  if (ugrestrictflag==1) then
     write(mpifileptr, *) "Configurations will be restricted to total parity= ", ugrestrictval
  endif
  if (restrictflag.eq.1) then
     write(mpifileptr, *) "Configurations will be restricted to spin projection", restrictms, "/2"
  endif
  if (all_spinproject.ne.0) then
     write(mpifileptr,*) " Configurations will be restricted to spin ",spin_restrictval 
  endif
  write(mpifileptr, *) 
  write(mpifileptr, *) "***********************    Initial state      ***********************   "
  write(mpifileptr, *) 
  if (loadspfflag.eq.1) then
     write(mpifileptr, *) "Spfs will be loaded from files "
     if (numskiporbs.gt.0) then
        write(mpifileptr,*) "Skipping orbitals on file.  Orbitals ", (orbskip(i),i=1,numskiporbs)
     endif
  else
     write(mpifileptr, *) "Spfs will be one-electron eigfuncts."
  endif
  write(mpifileptr, *) 
  if (threshflag.ne.0) then

     if (improvedquadflag.gt.1) then
        write(mpifileptr,*) "Spf Quad flag is ON (quadflag>1).  Start time ", quadstarttime
     else
        write(mpifileptr,*) "Sfp Quad flag is OFF."
     endif
     if (improvednatflag.ne.0) then
        write(mpifileptr,*) "Improvednatflag is ON."
     elseif (improvedfockflag.ne.0) then
        write(mpifileptr,*) "Improvedfockflag is ON."
     else
        write(mpifileptr,*) "Improvednatflag and improvedfockflag are OFF."
     endif
     write(mpifileptr, *) 
     if (improvedquadflag.eq.1.or.improvedquadflag.eq.3) then
        write(mpifileptr,*) "Avector Quad flag is ON (quadflag=1,3)."
     else
        write(mpifileptr,*) "Avector Quad flag is OFF."
     endif
  end if
  if (loadavectorflag.eq.1) then
     write(mpifileptr, *) "Avector will be loaded from files.  Number of files= ",numavectorfiles
     if (numholes.gt.0) then
        OFLWR; WRFL "We have holes:", numholes, " concurrent holes with ", numholecombo," wfns combined "; WRFL; CFL
        excitations=0
     endif
     if (excitations.gt.0) then
        OFLWR; WRFL "We have exitations: ",excitations, " concurrent excitations with ", excitecombos," wfns combined "; WRFL; CFL
     endif
  else
     write(mpifileptr, *) "Avector will be obtained from diagonalization."
  endif
  write(mpifileptr, *) 
  write(mpifileptr, *) "***********************    Parameters: propagation    ***********************   "
  write(mpifileptr, *)
  write(mpifileptr,*)  " PAR_TIMESTEP IS ", par_timestep, " LITTLESTEPS IS ", littlesteps
  write(mpifileptr,*)  "  prepropflag is ", prepropflag, " postpropflag is ", postpropflag
  write(mpifileptr,*)  "  step_flag is ", step_flag
  write(mpifileptr,*)

  if (messflag.ne.0) then
     write(mpifileptr,*) "MESSFLAG is on -- messing with spfs.  Messamount=", messamount; write(mpifileptr,*)
  endif
  if (spf_flag /= 1) then
     write(mpifileptr, *) "Spfs will be held CONSTANT. (except for constraint)";    write(mpifileptr,*)
  endif

  if (threshflag.eq.1) then
     write(mpifileptr, *) "Calculation will be stopped with threshold ", stopthresh, "; timestep is ", par_timestep
  else
     write(mpifileptr,'(A40,I10,2F24.8)') "# of steps, final time:",   numpropsteps,  finaltime
  endif

  write(mpifileptr,*);  write(mpifileptr,*)
  if (cmf_flag.ne.0) then
     write(mpifileptr, *) "*******  USING POLYNOMIAL MEAN FIELDS/MAGNUS A-VECTOR PREDICTOR/CORRECTOR   *******  "
     write(mpifileptr,*)
  else
     write(mpifileptr, *) "*************   Variable mean fields    *************   "
     write(mpifileptr,*)
  endif
  if (sparseconfigflag==1) then
     write(mpifileptr, *)  "    Will use sparse configuration routines."
     if (sparseopt.eq.0) then
        write(mpifileptr,*) "   Direct CI, sparseopt=0"
     else
        write(mpifileptr,*) "   Explicit sparse matrix, sparseopt=1"
     endif
     if (par_consplit.ne.0) then
        write(mpifileptr,*) "   Vectors are distributed, par_consplit.ne.0" 
     else
        write(mpifileptr,*) "   Vectors are not distributed, par_consplit=0" 
     endif
     write(mpifileptr,*)    "          Lanczosorder is ", lanczosorder
     write(mpifileptr,*)    "          Lanthresh is    ", lanthresh
     write(mpifileptr,*)    "          Aorder is       ", aorder
     write(mpifileptr,*)    "          Maxaorder is    ", maxaorder
     write(mpifileptr,*)    "          Aerror is       ", aerror
  else
     write(mpifileptr,*)    "    Using nonsparse configuration routines."
  endif
  iiflag=0

  write(mpifileptr,*) " Jacobian options:"
  write(mpifileptr,*) "    Jacprojorth=", jacprojorth
  write(mpifileptr,*) "    Jacsymflag=", jacsymflag
  if (constraintflag.ne.0) then
     write(mpifileptr,*) "    Jacgmatthird=", jacgmatthird
  endif

  select case (intopt)
  case(4)
     write(mpifileptr,*) " Using VERLET integration, expo first step."
     write(mpifileptr,*) "     Verletnum= ", verletnum
     iiflag=1
  case(3)
     write(mpifileptr,*) " Using EXPONENTIAL integration."
     write(mpifileptr,*) "    Expotol    =", expotol
     write(mpifileptr,*) "    Maxexpodim=", maxexpodim
  case(0)
     write(mpifileptr, *) "RK integration.  Recommend errors 1.d-8"
     write(mpifileptr, *) "Relerr, Myrelerr=", relerr, myrelerr;      iiflag=1
  case(1)
     write(mpifileptr, *) "GBS integration.  Recommend errors 1.d-8";     iiflag=1
     write(mpifileptr, *) "Relerr, Myrelerr=", relerr, myrelerr;      iiflag=1
  case(2)
     write(mpifileptr, *) "DLSODPK integration disabled (intopt=2) !!!"; CFLST
  case default
     write(mpifileptr,*) "Intopt not recognized: ", intopt; CFLST
  end select
  write(mpifileptr,*) 
  write(mpifileptr, '(A40, E10.3)') " Density matrix regularized with denreg= ", denreg
  write(mpifileptr,*) 

  if (constraintflag.ne.0) then
     select case (constraintflag)
     case(1)
        write(mpifileptr, *) "Using constraintflag=1, density matrices with full lioville solve (assume full, constant off-block diagonal)."
     case(2)
        OFLWR "Using true Dirac-Frenkel equation for constraint."
        WRFL "     dfrestrictflag = ", df_restrictflag;      
        select case(conway)
        case(0)
           WRFL "McLachlan constraint"
        case(1)
           WRFL "50/50 SVD"
        case(2)
           WRFL "Lagrangian constraint"
        case(3)
           WRFL "Lagrangian with epsilon times McLachlan.  epsilon=conprop=",conprop
        case default
           WRFL "Conway not supported ", conway; CFLST
        end select
     case(0)
        write(mpifileptr, *) "Using zero constraint, constraintflag=0"
     case default
        write(mpifileptr,*) "Constraintflag error ", constraintflag;     call mpistop()
     end select
     if (lioreg.le.0.d0) then
        write(mpifileptr,*) "No regularization of lioville solve."
     else
        write(mpifileptr,*) "Lioville solve regularized with lioreg=", lioreg
     endif
  endif
  

  write(mpifileptr, *)
  write(mpifileptr, *) "****************************************************************************"
  write(mpifileptr, *)
  write(mpifileptr,*) " Autotimestep= ", autotimestep
  write(mpifileptr,*) " Numpropsteps= ", numpropsteps

!!$ IMPLEMENT ME (DEPRECATE fluxinterval as namelist input)  
!!$ write(mpifileptr,*) "Fluxsteps is ", fluxsteps," Fluxtimestep is ", fluxtimestep

  write(mpifileptr, *)
  write(mpifileptr, *) "*****************************************************************************"
  write(mpifileptr, *)
  if (skipflag.ne.0) then
     write(mpifileptr,*) "   ****************************************"
     write(mpifileptr,*) "     SKIPPING CALCULATION!  Doing analysis."
     write(mpifileptr,*) "   ****************************************"
  endif
  write(mpifileptr,*) ;  call closefile()
  call write_actions()
  
end subroutine getparams


subroutine getpulse(no_error_exit_flag)   !! if flag is 0, will exit if &pulse is not read
  use parameters
  use pulse_parameters
  use mpimod
  use pulsesubmod
  use utilmod
  implicit none
  NAMELIST /pulse/ omega,pulsestart,pulsestrength, velflag, omega2,phaseshift,intensity,pulsetype, &
       pulsetheta,pulsephi, longstep, numpulses, reference_pulses, minpulsetime, maxpulsetime, chirp, ramp
  real*8 ::  time, fac, pulse_end, estep
  DATATYPE :: pots1(3),pots2(3),pots3(3), pots4(3), pots5(3), csumx,csumy,csumz
  integer :: i, myiostat, ipulse,no_error_exit_flag
  character (len=12) :: line
  real*8, parameter :: epsilon=1d-4
  integer, parameter :: neflux=10000
  complex*16 :: lenpot(0:neflux,3),velpot(0:neflux,3)
  real*8 :: pulseftsq(0:neflux), vpulseftsq(0:neflux)

  open(971,file=inpfile, status="old", iostat=myiostat)
  if (myiostat/=0) then
     OFLWR "No Input.Inp found, not reading pulse. iostat=",myiostat; CFL
  else
     read(971,nml=pulse,iostat=myiostat)
     if (myiostat.ne.0.and.no_error_exit_flag.eq.0) then
        OFLWR "Need &pulse namelist input!!"; CFLST
     endif
  endif
  close(971)
  call openfile()
  if (tdflag.ne.0) then
     line="PULSE IS ON:"
  else
     line="READ PULSE: "
  endif
  select case (velflag)
  case (0)
     write(mpifileptr, *) line,"   length."
  case(1)
     write(mpifileptr, *) line,"   velocity, usual way."
  case(2)
     write(mpifileptr, *) line,"   velocity, DJH way."
  end select
  if (no_error_exit_flag.ne.0) then    !! mcscf.  just need velflag.
     return
  endif
  write(mpifileptr, *) ;  write(mpifileptr, *) "NUMBER OF PULSES:  ", numpulses;  write(mpifileptr, *) 

  lastfinish=0.d0
  do ipulse=1,numpulses
     if (pulsephi(ipulse).ne.0d0) then
        offaxispulseflag=1
     endif

     write(mpifileptr, *) "    -----> Pulse ", ipulse," : "

     if (pulsetype(ipulse).eq.1.or.pulsetype(ipulse).eq.5) then
        fac=omega(ipulse)
     else
        fac=omega2(ipulse)
     endif
     if (intensity(ipulse).ne.-1.d0) then !! overrides pulsestrength
        pulsestrength(ipulse) = sqrt(intensity(ipulse)/3.5094478d0)/fac
     else
        intensity(ipulse) = (fac*pulsestrength(ipulse))**2 * 3.5094478d0  !! just output
     endif
     select case (pulsetype(ipulse))
     case (1)
        write(mpifileptr,*) "Pulse type is 1: single sine squared envelope"
        write(mpifileptr, *) "Omega, pulsestart, pulsefinish, pulsestrength:"
        write(mpifileptr, '(8F18.12)') omega(ipulse), pulsestart(ipulse), pulsestart(ipulse) + pi/omega(ipulse), pulsestrength(ipulse)
     case (2,3)
        write(mpifileptr,*) "Pulse type is 2 or 3: envelope with carrier"
        write(mpifileptr, *) "   chirp:           ", chirp(ipulse)
        write(mpifileptr, *) "   ramp:           ", ramp(ipulse), " Hartree"
        write(mpifileptr, *) "   Envelope omega:  ", omega(ipulse)
        write(mpifileptr, *) "   Pulse omega:     ", omega2(ipulse)
        write(mpifileptr, *) "   Pulsestart:      ",pulsestart(ipulse)
        write(mpifileptr, *) "   Pulsestrength:   ",pulsestrength(ipulse)
        write(mpifileptr, *) "   Intensity:       ",intensity(ipulse), " x 10^16 W cm^-2"
        write(mpifileptr, *) "   Pulsetheta:      ",pulsetheta(ipulse)
        write(mpifileptr, *) "   Pulsephi:      ",pulsephi(ipulse)
        write(mpifileptr, *) "   Pulsefinish:     ",pulsestart(ipulse) + pi/omega(ipulse)
        if (pulsetype(ipulse).eq.3) then
           write(mpifileptr,*) "---> Pulsetype 3; longstep = ", longstep(ipulse)
        else if (pulsetype(ipulse).eq.2) then
           write(mpifileptr,*) "---> Pulsetype 2."
        endif
     case (4)
        WRFL "Pulse type is 4, cw"
        write(mpifileptr, *) "   Duration omega:  ", omega(ipulse)
        write(mpifileptr, *) "   Pulse omega:     ", omega2(ipulse)
        write(mpifileptr, *) "   Pulsestrength:   ",pulsestrength(ipulse)
        write(mpifileptr, *) "   Intensity:       ",intensity(ipulse), " x 10^16 W cm^-2"
        write(mpifileptr, *) "   Pulsetheta:      ",pulsetheta(ipulse)
        write(mpifileptr, *) "   Pulsephi:      ",pulsephi(ipulse)

     end select
     if (lastfinish.lt.pulsestart(ipulse)+pi/omega(ipulse)) then
        lastfinish=pulsestart(ipulse)+pi/omega(ipulse)
     endif
  end do
  
#ifdef REALGO
  write(mpifileptr,*) "Pulse not available with real prop.";  call closefile();  call mpistop
#endif

  if (tdflag.ne.0) then
     finaltime=lastfinish+par_timestep*4

!! I think this is right

     pulse_end = 2*pi/(pulseft_estep*(neflux+1)/neflux)

     if (pulse_end.lt.finaltime) then
        pulse_end=finaltime
     endif

     estep = 2*pi/(pulse_end*(neflux+1)/neflux)

     if (finaltime.lt.minpulsetime) then
        write(mpifileptr,*) " Enforcing minpulsetime =        ", minpulsetime
        finaltime=minpulsetime
     endif
     if (finaltime.gt.maxpulsetime) then
        write(mpifileptr,*) " Enforcing maxpulsetime =        ", maxpulsetime
        finaltime=maxpulsetime
     endif
     numpropsteps=floor((finaltime+0.0000000001d0)/par_timestep) +1
     finaltime=numpropsteps*par_timestep
     write(mpifileptr, *) "    ---> Resetting finaltime to ", finaltime
     write(mpifileptr, *) "                numpropsteps to ", numpropsteps
     write(mpifileptr, *) "   ... Writing pulse to file"
     CFL

     if (myrank.eq.1) then
        open(886, file="Dat/Pulse.Datx", status="unknown",iostat=myiostat)
        call checkiostat(myiostat," opening Dat/Pulse.Datx")
        open(887, file="Dat/Pulse.Daty", status="unknown")
        open(888, file="Dat/Pulse.Datz", status="unknown")

        csumx=0; csumy=0; csumz=0
        do i=0,neflux
           time=i*pulse_end/neflux

!! checking that E(t) = d/dt A(t)   and  A(t) = integral E(t)

           call vectdpot(time,                  0,pots1,1)
           call vectdpot(time,                  1,pots2,1)
           call vectdpot(time-epsilon,          1,pots3,1)
           call vectdpot(time+epsilon,          1,pots4,1)
           call vectdpot(time+pulse_end/neflux ,0,pots5,1)

           lenpot(i,:)=pots1(:)
           velpot(i,:)=pots2(:)

           write(886,'(100F15.10)') time, pots1(1), pots2(1), (pots4(1)-pots3(1))/2/epsilon, csumx
           write(887,'(100F15.10)') time, pots1(2), pots2(2), (pots4(2)-pots3(2))/2/epsilon, csumy
           write(888,'(100F15.10)') time, pots1(3), pots2(3), (pots4(3)-pots3(3))/2/epsilon, csumz

           csumx=csumx+ pulse_end/neflux * pots5(1)
           csumy=csumy+ pulse_end/neflux * pots5(2)
           csumz=csumz+ pulse_end/neflux * pots5(3)
        enddo
        close(886);    close(887); close(888)

        do i=1,3
           call zfftf_wrap(neflux+1,lenpot(:,i))
           call zfftf_wrap(neflux+1,velpot(:,i))
        enddo

        lenpot(:,:)=lenpot(:,:)*pulse_end/neflux

        velpot(:,:)=velpot(:,:)*pulse_end/neflux

        pulseftsq(:)=abs(lenpot(:,1)**2)+abs(lenpot(:,2)**2)+abs(lenpot(:,3)**2)

        vpulseftsq(:)=abs(velpot(:,1)**2)+abs(velpot(:,2)**2)+abs(velpot(:,3)**2)

        open(885, file="Dat/Pulseftsq.Dat", status="unknown")
        open(886, file="Dat/Pulseft.Datx", status="unknown")
        open(887, file="Dat/Pulseft.Daty", status="unknown")
        open(888, file="Dat/Pulseft.Datz", status="unknown")

!! PREVIOUS OUTPUT IN Pulseft.Dat was vpulseftsq / 4

        do i=0,neflux
           write(885,'(100F30.12)',iostat=myiostat) i*estep, pulseftsq(i), vpulseftsq(i)
           call checkiostat(myiostat," wrting Dat/Pulse.Dat file")
           write(886,'(100F30.12)') i*estep, lenpot(i,1),velpot(i,1)
           write(887,'(100F30.12)') i*estep, lenpot(i,2),velpot(i,2)
           write(888,'(100F30.12)') i*estep, lenpot(i,3),velpot(i,3)
        enddo

        close(885);   close(886);   close(887);   close(888)

     endif

     call mpibarrier()

     OFLWR "       ... done Writing pulse to file"; CFL

  endif

end subroutine getpulse

end module getparammod
