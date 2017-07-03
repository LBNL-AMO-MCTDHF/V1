
#include "Definitions.INC"

!! ALL MODULES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! This file contains exponential propagation routines... They are separated nicely by the modules they
!!!   contain, mutually exclusively.  In order, in this file:
!!!   1) expomod files, for orbital propagation.  These are for expoprop, called from prop.f90, driver
!!!     for expokit subroutines.  These expokit subroutines for the orbital equation call 
!!!   2) jacmod files, for operation of jacobian, called by expokit, for orbital equation, and improvedquadflag=2,3
!!!   3) configexpomod files, for A-vector propagation
!!!   4) utilities that do not use any modules (besides mpimod and parameters)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   JACMOD SUBROUTINES, ETC 
!!   JACOPERATE IS THE SUBROUTINE PASSED TO EXPOKIT
!!   FOR ORBITAL PROPAGATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module jacmod
  implicit none
  DATATYPE, allocatable :: jacvect(:,:)     !! orbitals used to construct jacobian
  DATATYPE, allocatable :: jacvectout(:,:)  !! rho^-1 W jacvect = Q jacvect
  real*8 :: jactime
  integer :: allocated=0

end module


module jacopmod
  use mpisubmod

contains

  subroutine jacoperate00(lowspf,highspf,dentimeflag,conflag,inspfs,outspfs)
    use parameters
    use jacmod
    use mpimod
    use xxxmod    !! drivingorbs.... hmmm could just make wrapper but whatever
    use linearmod
    use jactimingmod
    use orbprojectmod
    use derivativemod
    use orbgathersubmod
    use pulsesubmod
    implicit none
    integer,intent(in) :: lowspf,highspf,dentimeflag,conflag
    DATATYPE,intent(in) ::  inspfs(spfsize,nspf)
    DATATYPE,intent(out) :: outspfs(spfsize,lowspf:highspf)
    integer :: ii,itop,getlen,numspf,myiostat,itime,jtime,jjj
    DATATYPE :: csum, nulldouble(2),pots(3)
    real*8 :: facs(0:1),rsum
    DATATYPE :: bigwork(spfsize,nspf),   workspfs(spfsize,lowspf:highspf),& 
         tempspfs(spfsize,lowspf:highspf), tempmat(nspf,lowspf:highspf)       !! AUTOMATIC

    numcalledhere=numcalledhere+1

    if (highspf.lt.lowspf) then
       call waitawhile()
       print *, "ERRROROR JACOP00",highspf,lowspf;
       call waitawhile()
       stop
    endif

!!$ 06-16 something is wonky with jacsymflag.. debugging frozenspfs
!!$ still allowing it though
!!$    if (jacsymflag.ne.0) then
!!$       OFLWR "something is wonky, set jacsymflag.eq.0"; CFLST
!!$    endif

    bigwork=0; workspfs=0; tempspfs=0; tempmat=0

    numspf=highspf-lowspf+1

    if (effective_cmf_linearflag.eq.1) then
       itop=1;
       facs(0)=(jactime-firsttime)/(lasttime-firsttime);     facs(1)=1d0-facs(0)
    else
       itop=0;     facs(0)=1d0;     facs(1)=0d0
    endif

    outspfs(:,lowspf:highspf)=0.d0

!! ** term with inspfs on far right ** !!

    do ii=0,itop
     
       if (jacsymflag.ne.0) then

          call myclock(itime)
          call project00(lowspf,highspf,inspfs(:,lowspf:highspf),&
               bigwork(:,lowspf:highspf),jacvect) 
          call myclock(jtime); times(2)=times(2)+jtime-itime;   itime=jtime
          if (parorbsplit.eq.1) then
             call mpiorbgather_nz(bigwork,spfsize)
          endif
          call myclock(jtime); times(4)=times(4)+jtime-itime;   itime=jtime

          call actreduced00(lowspf,highspf,dentimeflag,jactime,bigwork,nulldouble,&
               workspfs,ii,0,0)

          outspfs(:,:)=outspfs(:,:)+workspfs(:,lowspf:highspf)*facs(ii)

          call actreduced00(lowspf,highspf,dentimeflag,jactime,inspfs,&
               nulldouble,tempspfs,ii,0,0)
          call myclock(jtime); times(1)=times(1)+jtime-itime; 

       else
           
          call myclock(itime)
          call actreduced00(lowspf,highspf,dentimeflag,jactime, inspfs,&
               nulldouble, workspfs,ii,0,0)
!!tempspfs used later
          tempspfs(:,lowspf:highspf)=workspfs(:,lowspf:highspf)
          outspfs(:,:)=outspfs(:,:)+workspfs(:,lowspf:highspf)*facs(ii)

          call myclock(jtime); times(1)=times(1)+jtime-itime;      
       endif

       call myclock(itime)
       call project00(lowspf,highspf,tempspfs(:,lowspf:highspf),&
            workspfs(:,lowspf:highspf),jacvect)
       outspfs(:,:)=outspfs(:,:)-workspfs(:,lowspf:highspf)*facs(ii)

       call myclock(jtime); times(2)=times(2)+jtime-itime;

!! terms from projector
!! timing in derproject
       call derproject00(lowspf,highspf,jacvectout(:,lowspf:highspf),&
            workspfs,jacvect,inspfs)
       outspfs(:,:)=outspfs(:,:)-workspfs(:,lowspf:highspf)*facs(ii)

       if (jacsymflag.ne.0) then
          call derproject00(lowspf,highspf,jacvect(:,lowspf:highspf),&
               bigwork(:,min(lowspf,nspf):highspf),jacvect,inspfs)

          call myclock(itime)
          if (parorbsplit.eq.1) then
             call mpiorbgather_nz(bigwork,spfsize)
          endif
          call myclock(jtime); times(4)=times(4)+jtime-itime;   itime=jtime

          call actreduced00(lowspf,highspf,dentimeflag,jactime,bigwork,nulldouble,&
               workspfs,ii,0,0)
          outspfs(:,:)=outspfs(:,:)+workspfs(:,lowspf:highspf)*facs(ii)

          call myclock(jtime); times(1)=times(1)+jtime-itime;  
       endif
        
!!  CONSTRAINT!!  FORGOTTEN I GUESS !! APR 2014

       if (constraintflag.ne.0.and.conflag.ne.0) then
          call myclock(itime)
          call op_gmat00(lowspf,highspf,inspfs,&
               workspfs(:,lowspf:highspf),ii,jactime,jacvect)
          outspfs(:,:)=outspfs(:,:)+workspfs(:,lowspf:highspf)*facs(ii)
          if (jacgmatthird.ne.0) then
             call der_gmat00(lowspf,highspf,jacvect,&
                  workspfs(:,lowspf:highspf),ii,jactime,jacvect,inspfs)
             outspfs(:,:)=outspfs(:,:)+workspfs(:,lowspf:highspf)*facs(ii)
          endif
          call myclock(jtime); times(5)=times(5)+jtime-itime;
       endif


       if (numfrozen.gt.0.and.exact_exchange.eq.0) then

!! EXCHANGE
          select case(exchange_mode)
          case(0)   !! original (1-P(phi)) frozenexchinvr

             call derproject00(1,nspf,yyy%frozenexchinvr(:,:,ii),bigwork,jacvect,inspfs)

          case(1)   !! with frozenexchmat

             call myclock(itime)
             call MYGEMM('N','N',spfsize,nspf,nspf,DATAONE,&
                  inspfs,spfsize,  yyy%frozenexchmat(:,:,ii),nspf,&
                  DATAZERO, bigwork(:,:), spfsize)
             call myclock(jtime); times(6)=times(6)+jtime-itime;

          case default
             OFLWR "exchange_mode not supported",exchange_mode; CFLST
          end select

          call myclock(itime)
          if (dentimeflag.ne.0) then
!! TIMEFAC and facs HERE
             csum=timefac*facs(ii)
             call MYGEMM('N','N', spfsize,numspf,nspf,csum, &
                  bigwork(:,:),spfsize, &
                  yyy%invdenmat(:,lowspf:highspf,ii), nspf, DATAZERO, &
                  tempspfs(:,lowspf:highspf), spfsize)
          else
             tempspfs(:,lowspf:highspf)= & !! no more factor (-1) * &
                  bigwork(:,lowspf:highspf)*facs(ii)
          endif
          call myclock(jtime); times(6)=times(6)+jtime-itime;

          outspfs(:,lowspf:highspf)=outspfs(:,lowspf:highspf) &
               - tempspfs(:,lowspf:highspf)
       endif

!! DRIVING (PSI-PRIME)

       if (drivingflag.ne.0) then
          if (dentimeflag.eq.0) then
             OFLWR "error, no driving for quad!!"; CFLST !! invdenmat already in drivingorbs
          endif
          rsum=0
          call vectdpot(jactime,velflag,pots,-1)
          do jjj=1,3
             rsum=rsum+abs(pots(jjj))**2
          enddo
          if (rsum.ne.0d0) then
             workspfs(:,lowspf:highspf)=&
                  pots(1)*yyy%drivingorbsxx(:,lowspf:highspf,ii)+&
                  pots(2)*yyy%drivingorbsyy(:,lowspf:highspf,ii)+&
                  pots(3)*yyy%drivingorbszz(:,lowspf:highspf,ii)
             call derproject00(lowspf,highspf,workspfs,tempspfs,jacvect,inspfs)
             outspfs(:,:)=outspfs(:,:)-tempspfs(:,lowspf:highspf)*facs(ii)*timefac
          endif
       endif

    enddo  !! do ii=0,itop

    if ((myrank.eq.1).and.(notiming.eq.0)) then
       if (numcalledhere==1) then
          open(8577, file=timingdir(1:getlen(timingdir))//"/jacoperate.time.dat", &
               status="unknown",iostat=myiostat)
          call checkiostat(myiostat,"opening jacoperate timing file")
          write(8577,'(T16,100A9)',iostat=myiostat) &
               "actreduced ", &
               "project", &
               "derproject", &
               "MPI", &
               "constraint", &
               "frozen", &
               "total"

          close(8577)
          call checkiostat(myiostat,"writing jacoperate timing file")
       endif
       if (mod(numcalledhere,timingout).eq.0) then
          open(8577, file=timingdir(1:getlen(timingdir))//"/jacoperate.time.dat", &
               status="unknown", position="append",iostat=myiostat)
          call checkiostat(myiostat,"opening jacoperate timing file")
          write(8577,'(A3,F12.4,15I9)',iostat=myiostat) "T= ", jactime,  times(1:7)/1000
          call checkiostat(myiostat,"writing jacoperate timing file")
          close(8577)
       endif
    endif

  end subroutine jacoperate00

  subroutine jacoperate0(dentimeflag,conflag,inspfs,outspfs)
    use parameters
    use jactimingmod
    use orbgathersubmod
    implicit none
    integer,intent(in) :: dentimeflag,conflag
    DATATYPE,intent(in) ::  inspfs(spfsize,nspf)
    DATATYPE,intent(out) :: outspfs(spfsize,nspf)
    integer :: lowspf,highspf,itime,jtime,atime,btime

    call myclock(atime)

    lowspf=1; highspf=nspf
    if (parorbsplit.eq.1) then
       call getOrbSetRange(lowspf,highspf)
    endif

    if (highspf.ge.lowspf) then
       call jacoperate00(lowspf,highspf,dentimeflag,conflag,inspfs,&
            outspfs(:,lowspf:highspf))
    endif
    if (parorbsplit.eq.1) then
       call myclock(itime)
       call mpiorbgather(outspfs,spfsize)
       call myclock(jtime);      times(4)=times(4)+jtime-itime
    endif

    call myclock(btime); times(7)=times(7)+btime-atime;
  
  end subroutine jacoperate0


!! SUBROUTINE PASSED TO EXPOKIT FOR ORBITAL PROPAGATION

  subroutine jacopcompact(com_inspfs,com_outspfs)
    use parameters
    use spfsubmod
    implicit none
    DATATYPE,intent(in) :: com_inspfs(spfsmallsize,nspf)
    DATATYPE,intent(out) :: com_outspfs(spfsmallsize,nspf)
    DATATYPE ::  inspfs(spfsize,nspf), outspfs(spfsize,nspf)  !! AUTOMATIC

    inspfs=0; outspfs=0

    call spfs_expand(com_inspfs,inspfs)
    call jacoperate0(1,1,inspfs,outspfs)
    call spfs_compact(outspfs,com_outspfs)
  end subroutine jacopcompact

!! SUBROUTINE PASSED TO EXPOKIT FOR ORBITAL PROPAGATION

  subroutine jacoperate(inspfs,outspfs)
    use parameters
    implicit none
    DATATYPE,intent(in) :: inspfs(spfsize,nspf)
    DATATYPE,intent(out) :: outspfs(spfsize,nspf)
    call jacoperate0(1,1,inspfs,outspfs)
  end subroutine jacoperate


  subroutine parjacoperate0(dentimeflag,conflag,inspfs,outspfs)
    use parameters
    use jactimingmod
    use mpi_orbsetmod
    use orbgathersubmod
    implicit none
    integer,intent(in) :: dentimeflag,conflag
    DATATYPE,intent(in) ::  inspfs(spfsize,firstmpiorb:firstmpiorb+orbsperproc-1)
    DATATYPE,intent(out) :: outspfs(spfsize,firstmpiorb:firstmpiorb+orbsperproc-1)
    DATATYPE :: workspfs(spfsize,nspf)   !! AUTOMATIC
    integer :: lowspf,highspf,itime,jtime,atime,btime

    if (parorbsplit.ne.1) then
       OFLWR "ERROR, parjacoperate called but parorbsplit.ne.1",parorbsplit
       CFLST
    endif

    call myclock(atime)

    lowspf=1; highspf=nspf
    if (parorbsplit.eq.1) then
       call getOrbSetRange(lowspf,highspf)
    endif

    if (highspf.lt.lowspf) then
       call waitawhile()
       print *, "ERRROROR parjacop",highspf,lowspf;
       call waitawhile()
       stop
    endif

    workspfs=0
    workspfs(:,lowspf:highspf)=inspfs(:,lowspf:highspf)

    call myclock(itime)
    call mpiorbgather_nz(workspfs,spfsize)
    call myclock(jtime);      times(4)=times(4)+jtime-itime

    outspfs=0d0    !! padded

    call jacoperate00(lowspf,highspf,dentimeflag,conflag,&
         workspfs(:,:),outspfs(:,lowspf:highspf))

    call myclock(btime); times(7)=times(7)+btime-atime;
  
  end subroutine parjacoperate0


!! SUBROUTINE PASSED TO EXPOKIT FOR ORBITAL PROPAGATION

  subroutine parjacopcompact(com_inspfs,com_outspfs)
    use parameters
    use mpi_orbsetmod
    use spfsubmod
    implicit none
    DATATYPE,intent(in) :: com_inspfs(spfsmallsize,firstmpiorb:firstmpiorb+orbsperproc-1)
    DATATYPE,intent(out) :: com_outspfs(spfsmallsize,firstmpiorb:firstmpiorb+orbsperproc-1)
    DATATYPE ::  inspfs(spfsize,firstmpiorb:firstmpiorb+orbsperproc-1),  &
         outspfs(spfsize,firstmpiorb:firstmpiorb+orbsperproc-1)  !! AUTOMATIC

    inspfs=0; outspfs=0

    call spfs_expand_local(com_inspfs,inspfs)
    call parjacoperate0(1,1,inspfs,outspfs)
    call spfs_compact_local(outspfs,com_outspfs)
  end subroutine parjacopcompact

!! SUBROUTINE PASSED TO EXPOKIT FOR ORBITAL PROPAGATION

  subroutine parjacoperate(inspfs,outspfs)
    use parameters
    use mpi_orbsetmod
    implicit none
    DATATYPE,intent(in) :: inspfs(spfsize,firstmpiorb:firstmpiorb+orbsperproc-1)
    DATATYPE,intent(out) :: outspfs(spfsize,firstmpiorb:firstmpiorb+orbsperproc-1)
    call parjacoperate0(1,1,inspfs,outspfs)
  end subroutine parjacoperate

!! attempt for quad

  subroutine jacorth(inspfs,outspfs)
    use parameters
    use jacmod
    use orbprojectmod
    implicit none
    DATATYPE,intent(in) ::  inspfs(spfsize,nspf)
    DATATYPE,intent(out) :: outspfs(spfsize,nspf)
    call project(inspfs,outspfs,jacvect)
    outspfs(:,:)=inspfs(:,:)-outspfs(:,:)
  end subroutine jacorth
  
end module jacopmod


module jacinitmod
contains

!! SETS UP JACOBIAN FOR ORBITAL EXPO PROP

subroutine jacinit(inspfs, thistime)
  use parameters
  implicit none
  DATATYPE,intent(in) :: inspfs(spfsize,nspf) 
  real*8,intent(in) :: thistime
  call jacinit0(1,inspfs,thistime)
end subroutine jacinit


subroutine jacinit0(dentimeflag,inspfs, thistime)
  use parameters
  use linearmod
  use jacmod
  use derivativemod
  implicit none
  integer,intent(in) :: dentimeflag
  DATATYPE,intent(in) :: inspfs(spfsize,nspf) 
  real*8,intent(in) :: thistime
  DATATYPE :: nulldouble(2)
  DATATYPE,allocatable :: jactemp(:,:)
  real*8 :: gridtime

  if (allocated.eq.0) then
     allocate(jacvect(spfsize,nspf), jacvectout(spfsize,nspf))
     jacvect=0;jacvectout=0
  endif
  allocated=1

  jactime=thistime;  jacvect=inspfs; 

!! ONLY GOOD FOR CMF, LMF !!

  call actreduced0(dentimeflag,jactime,jacvect,nulldouble,jacvectout,0,0,0)

  if (effective_cmf_linearflag.eq.1) then
     allocate( jactemp(spfsize,nspf) );   jactemp=0
     gridtime=(jactime-firsttime)/(lasttime-firsttime) 
     if ((gridtime.lt.0.d0).or.(gridtime.gt.1)) then
        print *, "GGGRIDTIME ERR ", gridtime, jactime, firsttime, lasttime
        stop
     endif
     jacvectout=jacvectout*gridtime
     call actreduced0(dentimeflag,jactime,jacvect,jacvect,jactemp,1,0,0)
     jacvectout=jacvectout+(1.d0-gridtime)*jactemp
     deallocate(jactemp)
  endif

end subroutine jacinit0

end module jacinitmod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   EXPOPROP SUBROUTINE FOR ORBITAL PROPAGATION
!!   Driver for expokit call (z/dgphiv) for orbital equation
!!   z/dghpiv uses jacoperate and jacmod subroutines
!!   (at bottom of this file; expomod are at top)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module expospfpropmod
contains

subroutine expospfprop(time1,time2,in_inspfs, numiters)
  use parameters
  use mpimod
  use jacopmod
  use jacinitmod
  use orbdermod
  use mpi_orbsetmod
  use orbgathersubmod
  use mpisubmod
  use expokitmod, only: dgexpthird, dgexpthirdxxx2, dgphiv, dgphivxxx2
  use spfsubmod
  implicit none
  real*8,intent(in) :: time1,time2
  DATATYPE,intent(inout) :: in_inspfs(spfsize,nspf)
  integer,intent(out) :: numiters
  real*8, save :: tempstepsize = -1d0
  real*8 :: midtime, tdiff, error, norm
  integer :: itrace, iflag,getlen,minflag
  integer :: expofileptr=805

!! lowers thisexpodim until number of internal expokit steps is 2 or less, 
!!  then goes back and sets
!!  exposet to 1.  once exposet is 1, does not decrease further but may
!!  increase.

  integer, save :: thisexpodim=0, icalled=0
  integer, save :: exposet=0   
!! 0 = not set; reduce to minimum      !! 1 = minimum found; don't decrease

  DATATYPE,allocatable :: aspfs(:,:), proppspfs(:,:),   outspfs(:,:), &
       com_aspfs(:,:), com_proppspfs(:,:),   com_outspfs(:,:),  tempspfs(:,:),&
       inspfs(:,:)
  integer :: idim, liwsp=0, lwsp=0,ttott,myiostat,maxnorbs,lastmpiorb
  real*8, allocatable :: wsp(:)
  integer, allocatable :: iwsp(:)

  lastmpiorb=firstmpiorb+orbsperproc-1
  if (parorbsplit.eq.1) then
!!TEMP? was    maxnorbs=nzprocsperset*orbsperproc
     maxnorbs=maxprocsperset*orbsperproc
  else
     maxnorbs=nspf
  endif

  allocate(aspfs(spfsize,maxnorbs), proppspfs(spfsize,maxnorbs),&
       outspfs(spfsize,maxnorbs),tempspfs(spfsize,maxnorbs),&
       inspfs(spfsize,maxnorbs))
  aspfs=0; proppspfs=0; outspfs=0; tempspfs=0; inspfs=0

  inspfs(:,1:nspf)=in_inspfs(:,:)

  if (orbcompact.ne.0) then
     allocate(com_aspfs(spfsmallsize,maxnorbs), com_proppspfs(spfsmallsize,maxnorbs),   &
          com_outspfs(spfsmallsize,maxnorbs))
  else
     allocate(com_aspfs(1,maxnorbs), com_proppspfs(1,maxnorbs), com_outspfs(1,maxnorbs))
  endif
  com_aspfs=0; com_proppspfs=0; com_outspfs=0

  tdiff=(time2-time1)
  midtime=(time2+time1)/2d0

  call jacinit(inspfs,midtime)

!!$ PREPARE FOR CALL TO DGEXPVxxx2

  if (parorbsplit.eq.1) then
     if (orbcompact.ne.0) then
        ttott=spfsmallsize*orbsperproc
     else
        ttott=spfsize*orbsperproc
     endif
  else
     if (orbcompact.ne.0) then
        ttott=spfsmallsize*nspf
     else
        ttott=totspfdim
     endif
  endif

#ifndef REALGO
  ttott=ttott*2
#endif

  idim=ttott
  if (parorbsplit.eq.1.or.parorbsplit.eq.3) then
     ttott=ttott*nprocs
  endif
  if (maxexpodim.gt.ttott-1) then
     OFLWR "Maxexpodim too big, resetting to", ttott-1; CFL
     maxexpodim=ttott-1
  endif

  icalled=icalled+1
  if (icalled.eq.1) then
     tempstepsize=par_timestep/littlesteps     
     if (expodim.gt.maxexpodim) then
        expodim=maxexpodim
     endif
     thisexpodim=expodim
     if ((myrank.eq.1).and.(notiming==0)) then
        open(expofileptr,file=timingdir(1:getlen(timingdir))//"/expo.dat",&
             status="unknown",iostat=myiostat)
        call checkiostat(myiostat,"opening expo.dat timing file")
        write(expofileptr,*,iostat=myiostat) " Exponential propagator output.  expotol=",expotol
        call checkiostat(myiostat,"writing expo.dat timing file")
        write(expofileptr,*);        close(expofileptr)
     endif
  endif

  if (mod(icalled,12).eq.0) then
     exposet=0
  endif
  numiters=0;  

  if (thisexpodim.lt.maxexpodim) then
     tempstepsize=tempstepsize*4
  else
!!     tempstepsize=tempstepsize*1.1
     tempstepsize=tempstepsize*expostepfac
  endif

  if ((myrank.eq.1).and.(notiming.eq.0)) then
     open(expofileptr,file=timingdir(1:getlen(timingdir))//"/expo.dat",&
          status="old", position="append",iostat=myiostat)
     call checkiostat(myiostat,"opening expo.dat timing file")
     write(expofileptr,*,iostat=myiostat) "Go Orbital Expoprop.  Tinit=", time1, &
          " thisexpodim=",thisexpodim, " step ", min(par_timestep/littlesteps,tempstepsize)
     call checkiostat(myiostat,"writing expo.dat timing file")
     close(expofileptr)
  endif

  if (myrank.eq.1.and.notiming.eq.0) then
     open(expofileptr,file=timingdir(1:getlen(timingdir))//"/expo.dat",&
          status="old", position="append",iostat=myiostat)
     call checkiostat(myiostat,"opening expo.dat timing file")
  else
!!$ opening /dev/null multiple times not allowed :<       
!!$ open(expofileptr,file="/dev/null",status="unknown")

     expofileptr=nullfileptr
  endif

  itrace=0   !! 1=print info
  norm=1.d0
  
  lwsp = idim*(thisexpodim+4)  + 6*(thisexpodim+3)**2 + 100
  if (lwsp.lt.0) then
     OFLWR "Oops, integer overflow, reduce maxexpodim", thisexpodim,lwsp
     WRFL "      it's ",maxexpodim," now."; CFLST
  endif
  allocate(wsp(lwsp)); wsp=0
  liwsp=thisexpodim + 100
  allocate(iwsp(liwsp)); iwsp=0
  
!! nodgexpthirdflag=1 HARDWIRE 10-2015  NOT SURE ABOUT DGEXPTHIRD SUBROUTINES

  iflag=0
  if  ((jacsymflag.ne.0).and.(jacprojorth.eq.0).and.&
       (constraintflag.eq.0.or.jacgmatthird.ne.0).and.(nodgexpthirdflag.eq.0)) then

!! homogeneous third order.  Simpler expression for propagator.

     aspfs=inspfs

     if (parorbsplit.eq.3) then
        call dgexpthirdxxx2(idim, thisexpodim, tdiff, aspfs, inspfs, expotol, norm, &
             wsp, lwsp, iwsp, liwsp, jacoperate, itrace, iflag,expofileptr,&
             tempstepsize,realpardotsub_par3,ttott)
     elseif (parorbsplit.eq.1) then
        if (firstmpiorb.le.nspf) then
           call dgexpthirdxxx2(idim, thisexpodim, tdiff, &
                aspfs(:,firstmpiorb:lastmpiorb), inspfs(:,firstmpiorb:lastmpiorb),&
                expotol, norm, wsp,lwsp, iwsp, liwsp, parjacoperate, itrace, iflag, &
                expofileptr,tempstepsize,realpardotsub_par1,ttott)
        endif
        call mpiorbgather(inspfs,spfsize)
     else
        call dgexpthird(idim, thisexpodim, tdiff, aspfs, inspfs, expotol, norm, &
             wsp, lwsp, iwsp, liwsp, jacoperate, itrace, iflag,expofileptr,tempstepsize)
     endif

  else
     !! this does phi_1 = phi_0 + (exp(J t)-1)/J f(phi_0)

     call spf_linear_derivs(midtime,inspfs,aspfs) 

     call apply_spf_constraints(aspfs)

     proppspfs=0.d0

     if (orbcompact.ne.0) then
        call spfs_compact(aspfs,com_aspfs); 
        call spfs_compact(proppspfs,com_proppspfs);
        if (parorbsplit.eq.3) then
           call dgphivxxx2(idim, thisexpodim, tdiff, com_aspfs, &
                com_proppspfs, com_outspfs, expotol, norm, &
                wsp, lwsp, iwsp, liwsp, jacopcompact, itrace, iflag,&
                expofileptr,tempstepsize,realpardotsub_par3,ttott)
        elseif (parorbsplit.eq.1) then
           if (firstmpiorb.le.nspf) then
              call dgphivxxx2(idim, thisexpodim, tdiff, &
                   com_aspfs(:,firstmpiorb:lastmpiorb), &
                   com_proppspfs(:,firstmpiorb:lastmpiorb),&
                   com_outspfs(:,firstmpiorb:lastmpiorb), expotol, norm, wsp, &
                   lwsp, iwsp, liwsp, parjacopcompact, itrace, iflag, expofileptr,&
                   tempstepsize,realpardotsub_par1,ttott)
           endif
           call mpiorbgather(com_outspfs,spfsmallsize)
        else
           call dgphiv(idim, thisexpodim, tdiff, com_aspfs, &
                com_proppspfs, com_outspfs, expotol, norm, &
                wsp, lwsp, iwsp, liwsp, jacopcompact, itrace, iflag,&
                expofileptr,tempstepsize)
        endif
        call spfs_expand(com_outspfs,outspfs);
     else
        if (parorbsplit.eq.3) then
           call dgphivxxx2(idim, thisexpodim, tdiff, aspfs, proppspfs, &
                outspfs, expotol, norm, wsp, lwsp, iwsp, liwsp, jacoperate, &
                itrace, iflag,expofileptr,tempstepsize,realpardotsub_par3,ttott)
        elseif (parorbsplit.eq.1) then
           if (firstmpiorb.le.nspf) then
              call dgphivxxx2(idim, thisexpodim, tdiff, &
                   aspfs(:,firstmpiorb:lastmpiorb),proppspfs(:,firstmpiorb:lastmpiorb), &
                   outspfs(:,firstmpiorb:lastmpiorb), expotol, norm, wsp, lwsp, iwsp, &
                   liwsp, parjacoperate,itrace, iflag,expofileptr,&
                   tempstepsize,realpardotsub_par1,ttott)
           endif
           call mpiorbgather(outspfs,spfsize)
        else
           call dgphiv(idim, thisexpodim, tdiff, aspfs, proppspfs, &
                outspfs, expotol, norm, wsp, lwsp, iwsp, liwsp, jacoperate, &
                itrace, iflag,expofileptr,tempstepsize)
        endif
     endif

     inspfs=inspfs+outspfs
     
  endif

  if (expofileptr.ne.nullfileptr) then
     close(expofileptr)
  endif

  minflag=iflag
  call mympiimax(iflag)
  call mympiimin(minflag)
  if (iflag/=0.or.minflag/=0) then
     OFLWR "Expo error spf ", iflag,minflag; CFLST
  endif

!! should check to see if any numbers are unequal among processors.

  call mympiimax(iwsp(4))

  if (iwsp(4).gt.1) then
     thisexpodim=min(maxexpodim,max(thisexpodim+5,ceiling(thisexpodim*1.2)))
     exposet=1
  endif

  numiters=numiters+iwsp(1)

  call apply_spf_constraints(inspfs)

  call spf_orthogit(inspfs, error)
  
  in_inspfs(:,:)=inspfs(:,1:nspf)

  if ((myrank.eq.1).and.(notiming.eq.0)) then
     open(expofileptr,file=timingdir(1:getlen(timingdir))//"/expo.dat",&
          status="old", position="append",iostat=myiostat)
     call checkiostat(myiostat,"opening expo.dat timing file")
     write(expofileptr,*,iostat=myiostat) &
          "   End expo. Orthog error, steps, stepsize, iterations",&
          error,iwsp(4),tempstepsize,iwsp(1); 
     call checkiostat(myiostat,"writing expo.dat timing file")
     if (abs(error).gt.1d-5) then
        write(expofileptr,*) "    **** THATSA BIG ERROR I THINK! ***"
     endif
     write(expofileptr,*);     close(expofileptr)
  endif

  if ((exposet==0)) then
     thisexpodim=max(5,floor(thisexpodim*0.95))
  endif


  deallocate(wsp,iwsp)
  deallocate(aspfs, proppspfs, outspfs, tempspfs, inspfs)
  deallocate(com_aspfs, com_proppspfs, com_outspfs)
  return
  return  !! END SUBROUTINE
  return

contains
  subroutine realpardotsub_par3(one,two,n,out)
    implicit none
    integer,intent(in) :: n
    real*8,intent(in) :: one(n),two(n)
    real*8,intent(out) :: out
    real*8 :: sum
    sum=DOT_PRODUCT(one,two)
    call mympirealreduceone(sum)
    out=sum
  end subroutine realpardotsub_par3

  subroutine realpardotsub_par1(one,two,n,out)
    implicit none
    integer,intent(in) :: n
    real*8,intent(in) :: one(n),two(n)
    real*8,intent(out) :: out
    real*8 :: sum
    sum=DOT_PRODUCT(one,two)
    call mympirealreduceone_local(sum,NZ_COMM_ORB(myorbset))
    out=sum
  end subroutine realpardotsub_par1

end subroutine expospfprop

end module expospfpropmod

!!$!! KEEPME
!!$function checknan2(input,size)
!!$  implicit none
!!$  integer :: size,i
!!$  logical :: checknan2
!!$  DATATYPE :: input(size)
!!$  do i=1,size
!!$  if (input(i).eq.input(i)+1) then
!!$    checknan2=.true.
!!$    return
!!$  endif
!!$  enddo
!!$  checknan2=.false.
!!$end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  CONFIGEXPOMOD SUBROUTINES: A-vector propagation using expokit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module expoavecpropmod
contains

subroutine expoavecprop(inavector,outavector,time,imc,numiters)
  use parameters
  use configpropmod
  use configmod
  use mpimod
  use sparse_parameters   !! nzflag
  use basissubmod
  use mpisubmod
  use expokitmod, only: dgexpvxxx2, dgphivxxx2
  implicit none
  integer :: times(100)=0, iitime= -1 , jjtime = -1
  integer,intent(in) :: imc
  integer,intent(out) :: numiters
  DATATYPE,intent(in) :: inavector(numr,www%firstconfig:www%lastconfig)
  DATATYPE,intent(out) :: outavector(numr,www%firstconfig:www%lastconfig)
                                                                !!  AUTOMATIC
  DATATYPE :: localvector(numr,www%maxdfbasisperproc),&         !!       DFWW
       drivingavecdf(numr,www%botdfbasis:www%topdfbasis+1),&    !!       DFWW
       smallvector(numr,dwwptr%maxdfbasisperproc),&             !!       FDWW
       smallvectorout(numr,dwwptr%maxdfbasisperproc), &         !!       FDWW
       zerovector(numr,dwwptr%maxdfbasisperproc),&              !!       FDWW
       smallvectortemp(numr,dwwptr%maxdfbasisperproc), &        !!       FDWW
       workdrivingavecdf(numr,dwwptr%botdfbasis:dwwptr%topdfbasis+1) !!  FDWW

  real*8 :: one,time
  real*8, save :: tempstepsize=-1d0
  integer :: itrace, iflag, numsteps,expofileptr=61142, liwsp=0, lwsp=0,getlen,&
       myiostat,ixx,minflag
#ifdef REALGO
  integer, parameter :: zzz=1
#else
  integer, parameter ::   zzz=2
#endif

  real*8, allocatable :: wsp(:)
  integer, allocatable :: iwsp(:)
  integer, save :: thisexpodim=-1, icalled=0,my_maxaorder=0
!! 0 = not set; reduce to minimum 
!! 1 = minimum found; don't decrease
  integer, save :: exposet=0

  if (sparseconfigflag.eq.0) then
     OFLWR "must use sparseconfigflag.ne.0 for expoavecprop"; CFLST
  endif

  call avectortimeset()

  smallvector=0;   smallvectorout=0;   zerovector=0;  smallvectortemp=0
  workdrivingavecdf=0; localvector=0; drivingavecdf=0;

  icalled=icalled+1

  ixx = zzz * dwwptr%maxdfbasisperproc * numr

  if (icalled.eq.1) then
     tempstepsize=par_timestep/littlesteps

     my_maxaorder=max(1,min(zzz*www%numdfbasis*numr-1,maxaorder))

     thisexpodim=min(my_maxaorder,aorder)
  endif
  
  if (mod(icalled,20).eq.0) then
     exposet=0
  endif
  
  call configexpoinit(time,imc)

  if (thisexpodim.lt.my_maxaorder) then
     tempstepsize=tempstepsize*4
  else
     tempstepsize=tempstepsize*1.1
  endif

  if ((myrank.eq.1).and.(notiming==0)) then
     if (icalled.eq.1) then
        open(expofileptr,file=timingdir(1:getlen(timingdir))//"/avecexpo.dat",&
             status="unknown",iostat=myiostat)
        call checkiostat(myiostat,"opening avecexpo timing file")
        write(expofileptr,*,iostat=myiostat) " Avector lanczos propagator.  Order =",&
             thisexpodim," Aerror= ",aerror
        call checkiostat(myiostat,"writing avecexpo timing file")
        write(expofileptr,*);        close(expofileptr)
     endif

     open(expofileptr,file=timingdir(1:getlen(timingdir))//"/avecexpo.dat",&
          status="old", position="append",iostat=myiostat)
     call checkiostat(myiostat,"opening avecexpo timing file")
     write(expofileptr,*,iostat=myiostat) "Go Avector Lanczos.  time=", time, &
          " Order =",thisexpodim, "step ",min(par_timestep/littlesteps,tempstepsize)
     call checkiostat(myiostat,"writing avecexpo timing file")
     close(expofileptr)
  endif

  if (myrank.eq.1.and.notiming.eq.0) then
     open(expofileptr,file=timingdir(1:getlen(timingdir))//"/avecexpo.dat",&
          status="old", position="append",iostat=myiostat)
     call checkiostat(myiostat,"opening avecexpo timing file")
  else
!!$ opening /dev/null multiple times not allowed :<       
!!$ open(expofileptr,file="/dev/null",status="unknown")

     expofileptr=nullfileptr
  endif
  
  one=1.d0; itrace=0 ! running mode
  
  lwsp =  ixx*(thisexpodim+4) + 6*(thisexpodim+3)**2 + 100

  liwsp=thisexpodim+100
  allocate(wsp(lwsp),iwsp(liwsp));    
  wsp=0.d0; iwsp=0

  localvector(:,:)=0
  if (www%topconfig.ge.www%botconfig) then
     call basis_transformto_local(www,numr,inavector(:,www%botconfig),localvector(:,:))
  endif

  call avectortime(1)

  if (drivingflag.ne.0.and.www%topdfbasis.ge.www%botdfbasis) then
     call basis_transformto_local(www,numr,&
          workdrivingavec(:,www%botconfig:www%topconfig),&
          drivingavecdf(:,www%botdfbasis:www%topdfbasis))
  endif

  call mpibarrier()
  if (use_dfwalktype.and.shuffle_dfwalktype) then
     call basis_shuffle(numr,www,localvector,dwwptr,smallvector)
     if (drivingflag.ne.0) then
        call basis_shuffle(numr,www,drivingavecdf,dwwptr,workdrivingavecdf)
     endif
  else
     smallvector(:,:)=localvector(:,:)
     if (drivingflag.ne.0) then
        workdrivingavecdf(:,:)=drivingavecdf(:,:)
     endif
  endif

!! par_timestep is a-norm estimate, ok, whatever

  iflag=0
  if (nzflag.eq.0.or.dwwptr%nzrank.gt.0) then
     if (drivingflag.ne.0) then
        smallvectortemp(:,:)=0d0
        call parconfigexpomult_padded(smallvector,smallvectortemp)
        if (dwwptr%topdfbasis.ge.dwwptr%botdfbasis) then
           smallvectortemp(:,1:dwwptr%topdfbasis-dwwptr%botdfbasis+1)=  &
                smallvectortemp(:,1:dwwptr%topdfbasis-dwwptr%botdfbasis+1) + &
                workdrivingavecdf(:,dwwptr%botdfbasis:dwwptr%topdfbasis) * timefac 
        endif
        zerovector(:,:)=0d0
        call DGPHIVxxx2( ixx, thisexpodim, one, smallvectortemp,&
             zerovector, smallvectorout, aerror, par_timestep/littlesteps, wsp, lwsp, &
             iwsp, liwsp, parconfigexpomult_padded, itrace, iflag,expofileptr,&
             tempstepsize,realpardotsub,my_maxaorder+1)
        
        smallvectorout(:,:)=smallvectorout(:,:)+smallvector(:,:)
     else
        call DGEXPVxxx2( ixx, thisexpodim, one, smallvector, &
             smallvectorout, aerror, par_timestep/littlesteps, wsp, lwsp, &
             iwsp, liwsp, parconfigexpomult_padded, itrace, iflag,expofileptr,&
             tempstepsize,realpardotsub,my_maxaorder+1)
     endif
  else
     smallvectorout(:,:)=0d0 !! should be zero and stay zero but doing this anyway
  endif
  call mpibarrier()

  minflag=iflag
  call mympiimax(iflag)
  call mympiimin(minflag)
  if (iflag .ne. 0.or.minflag.ne.0) then
     OFLWR "Error expoavecprop: ", iflag,minflag; CFLST
  endif

  if (use_dfwalktype.and.shuffle_dfwalktype) then
     call basis_shuffle(numr,dwwptr,smallvectorout,www,localvector)
  else
     localvector(:,:)=smallvectorout(:,:)
  endif

  call avectortime(3)

  if (www%lastconfig.ge.www%firstconfig) then
     outavector(:,:)=0d0;   
     if (www%topconfig.ge.www%botconfig) then
        call basis_transformfrom_local(www,numr,localvector,outavector(:,www%botconfig))
     endif
  endif

  if (www%parconsplit.eq.0) then
     call mpiallgather(outavector,www%numconfig*numr,&
          www%configsperproc*numr,www%maxconfigsperproc*numr)
  endif

  if (expofileptr.ne.nullfileptr) then
     close(expofileptr)
  endif
   
  numsteps=iwsp(4);  numiters=iwsp(1)
   
  if ((exposet==0).and.(numsteps.eq.1).and.(thisexpodim.gt.2)) then
     thisexpodim=max(2,thisexpodim-1)
  endif

  if (myrank.eq.1.and.notiming.eq.0) then
     open(expofileptr,file=timingdir(1:getlen(timingdir))//"/avecexpo.dat",&
          status="old", position="append",iostat=myiostat)
     call checkiostat(myiostat,"opening avecexpo timing file")
     write(expofileptr,*,iostat=myiostat) &
          "   End avector prop.  steps, iterations, stepsize", &
          iwsp(4),iwsp(1),tempstepsize; 
     call checkiostat(myiostat,"writing avecexpo timing file")
     write(expofileptr,*);    
     call avectortimewrite(expofileptr)
     close(expofileptr)
  endif

  if (numsteps.gt.1) then
     exposet=1
     thisexpodim=min(floor(thisexpodim*1.2)+5, my_maxaorder)
  endif

  deallocate(wsp,iwsp)
 
  call avectortime(1)

contains

  subroutine avectortimeset()
    implicit none
    call myclock(iitime)
  end subroutine avectortimeset

  subroutine avectortimewrite(fileptr)
    implicit none
    integer,intent(in) :: fileptr
    integer :: myiostat
    write(fileptr,'(10(A13,I15))',iostat=myiostat) "Start/end",times(1)/1000,&
         "mult",times(2)/1000,"expokit",times(3)/1000
    call checkiostat(myiostat," writing in subroutine avectortimewrite")
  end subroutine avectortimewrite

  subroutine avectortime(which)
    implicit none
    integer,intent(in) :: which
    
!! times(1) = miscellaneous.  
!! times(2) = parconfigexpomult
!! times(3) = in-between (expokit)

    call myclock(jjtime); times(which)=times(which)+jjtime-iitime; iitime=jjtime
  end subroutine avectortime

  subroutine realpardotsub(one,two,n,out)
    implicit none
    integer,intent(in) :: n
    real*8,intent(in) :: one(n),two(n)
    real*8,intent(out) :: out
    real*8 :: sum
    sum=DOT_PRODUCT(one,two)
    call mympirealreduceone_local(sum,dwwptr%NZ_COMM)
    out=sum
  end subroutine realpardotsub

!! MULTIPLY BY A-VECTOR HAMILTONIAN MATRIX

!! NOTE BOUNDS !!  PADDED

  subroutine parconfigexpomult_padded0_gather(wwin,inconfigpointer,&
       insparsepointer,inavector,outavector)
    use fileptrmod
    use r_parameters
    use sparse_parameters
    use ham_parameters   !! timefac
    use mpimod
    use configexpotimemod
    use walkmod
    use configptrmod
    use sparseptrmod
    use sparsemultmod
    use basissubmod
    implicit none
    type(walktype),intent(in) :: wwin
    type(CONFIGPTR),intent(in) :: inconfigpointer
    type(SPARSEPTR),intent(in) :: insparsepointer
    DATATYPE,intent(in) :: inavector(numr,wwin%botdfbasis:wwin%botdfbasis+wwin%maxdfbasisperproc-1)
    DATATYPE,intent(out) :: outavector(numr,wwin%botdfbasis:wwin%botdfbasis+wwin%maxdfbasisperproc-1)
    DATATYPE,allocatable :: intemp(:,:)
    DATATYPE :: outtemp(numr,wwin%botconfig:wwin%topconfig+1)  !! AUTOMATIC

    call avectortime(3)

    if (sparseconfigflag.eq.0) then
       OFLWR "error, must use sparse for parconfigexpomult"; CFLST
    endif

    allocate(intemp(numr,wwin%numconfig));    intemp(:,:)=0d0;   outtemp=0

!! transform second to reduce communication?
!!   no, spin transformations done locally now.

    if (wwin%topdfbasis.ge.wwin%botdfbasis) then
       call basis_transformfrom_local(wwin,numr,inavector,&
            intemp(:,wwin%botconfig:wwin%topconfig))
    endif

    call mpiallgather_local(intemp,wwin%numconfig*numr,&
         wwin%nzconfsperproc(:)*numr, wwin%maxconfigsperproc*numr,&
         wwin%NZ_COMM,wwin%nzprocs,wwin%nzrank)

    outavector(:,:)=0d0    !!     inavector(:,:)   !! PADDED

    if (wwin%topconfig.ge.wwin%botconfig) then
       call sparseconfigmult_byproc(1,nprocs,wwin,intemp,outtemp, &
            inconfigpointer, insparsepointer, 1,1,1,1,&
            configexpotime,0,1,numr,0,imc)
    endif

    if (wwin%topdfbasis.ge.wwin%botdfbasis) then
       call basis_transformto_local(wwin,numr,outtemp,outavector)
    endif

    outavector=outavector*timefac

    deallocate(intemp)

    call avectortime(2)

  end subroutine parconfigexpomult_padded0_gather

#ifdef MPIFLAG

  subroutine parconfigexpomult_padded0_summa(wwin,inconfigpointer,&
       insparsepointer,inavector,outavector)
    use fileptrmod
    use r_parameters
    use sparse_parameters
    use ham_parameters   !! timefac
    use mpimod
    use configexpotimemod
    use walkmod
    use configptrmod
    use sparseptrmod
    use sparsemultmod
    use basissubmod
    implicit none
    type(walktype),intent(in) :: wwin
    type(CONFIGPTR),intent(in) :: inconfigpointer
    type(SPARSEPTR),intent(in) :: insparsepointer
    DATATYPE,intent(in) :: inavector(numr,wwin%botdfbasis:wwin%botdfbasis+wwin%maxdfbasisperproc-1)
    DATATYPE,intent(out) :: outavector(numr,wwin%botdfbasis:wwin%botdfbasis+wwin%maxdfbasisperproc-1)
    DATATYPE :: intemp(numr,wwin%maxconfigsperproc), outwork(numr,wwin%botconfig:wwin%topconfig+1),&
         outtemp(numr,wwin%botconfig:wwin%topconfig+1)   !! AUTOMATIC
    integer :: iproc,iiproc

    call avectortime(3)

    if (sparseconfigflag.eq.0) then
       OFLWR "error, must use sparse for parconfigexpomult"; CFLST
    endif
  
    outwork=0d0; outtemp=0; intemp=0

    do iiproc=1,wwin%nzprocs

       iproc=wwin%nzproclist(iiproc)

!! transform second to reduce communication?
!!   no, spin transformations done locally now.

       if (wwin%alltopconfigs(iproc).ge.wwin%allbotconfigs(iproc)) then
          if (myrank.eq.iproc) then
             intemp=0
             if (wwin%topdfbasis.ge.wwin%botdfbasis) then
                call basis_transformfrom_local(wwin,numr,inavector,intemp)
             endif
          endif

          call mympibcast_local(intemp,iiproc,wwin%nzconfsperproc(iiproc)*numr,&
               wwin%NZ_COMM)

          call sparseconfigmult_byproc(iproc,iproc,wwin,intemp,outtemp, &
               inconfigpointer, insparsepointer, 1,1,1,1,&
               configexpotime,0,1,numr,0,imc)
     
          outwork(:,:)=outwork(:,:)+outtemp(:,:)

       endif
    enddo

    outavector(:,:)=0d0   !!     inavector(:,:)   !! PADDED

    if (wwin%topdfbasis.ge.wwin%botdfbasis) then
       call basis_transformto_local(wwin,numr,outwork,outavector)
    endif

    outavector=outavector*timefac
  
    call avectortime(2)

  end subroutine parconfigexpomult_padded0_summa


  subroutine parconfigexpomult_padded0_circ(wwin,inconfigpointer,&
       insparsepointer,inavector,outavector)
    use fileptrmod
    use r_parameters
    use sparse_parameters
    use ham_parameters   !! timefac
    use mpimod
    use configexpotimemod
    use walkmod
    use configptrmod
    use sparseptrmod
    use sparsemultmod
    use basissubmod
    implicit none
    type(walktype),intent(in) :: wwin
    type(CONFIGPTR),intent(in) :: inconfigpointer
    type(SPARSEPTR),intent(in) :: insparsepointer
    DATATYPE,intent(in) :: inavector(numr,wwin%botdfbasis:wwin%botdfbasis+wwin%maxdfbasisperproc-1)
    DATATYPE,intent(out) :: outavector(numr,wwin%botdfbasis:wwin%botdfbasis+wwin%maxdfbasisperproc-1)
    DATATYPE :: workvector(numr,wwin%maxconfigsperproc), workvector2(numr,wwin%maxconfigsperproc),&
         outwork(numr,wwin%botconfig:wwin%topconfig+1),  &
         outtemp(numr,wwin%botconfig:wwin%topconfig+1)             !! AUTOMATIC
    integer :: iproc,prevproc,nextproc,deltaproc,iiproc

    call avectortime(3)

    if (sparseconfigflag.eq.0) then
       OFLWR "error, must use sparse for parconfigexpomult"; CFLST
    endif

!! doing circ mult slightly different than e.g. SINCDVR/coreproject.f90 
!!     and ftcore.f90, holding hands in a circle, prevproc and nextproc, 
!!     each chunk gets passed around the circle

    prevproc=mod(wwin%nzprocs+wwin%nzrank-2,wwin%nzprocs)+1
    nextproc=mod(wwin%nzrank,wwin%nzprocs)+1

    outwork=0d0; outtemp=0; workvector=0; workvector2=0

    if (wwin%topdfbasis.ge.wwin%botdfbasis) then
       call basis_transformfrom_local(wwin,numr,inavector,workvector)
    endif

    do deltaproc=0,wwin%nzprocs-1

!! PASSING BACKWARD (plus deltaproc)
       iiproc=mod(wwin%nzrank-1+deltaproc,wwin%nzprocs)+1
       iproc=wwin%nzproclist(iiproc)

       if (wwin%alltopconfigs(iproc).ge.wwin%allbotconfigs(iproc)) then
          call sparseconfigmult_byproc(iproc,iproc,wwin,workvector,outtemp, &
               inconfigpointer, insparsepointer, &
               1,1,1,1,configexpotime,0,1,numr,0,imc)
          outwork(:,:)=outwork(:,:)+outtemp(:,:)
       endif

!! PASSING BACKWARD
!! mympisendrecv(sendbuf,recvbuf,dest,source,...)

       call mympisendrecv_local(workvector,workvector2,prevproc,&
            nextproc,deltaproc,numr*wwin%maxconfigsperproc,wwin%NZ_COMM)
       workvector(:,:)=workvector2(:,:)
    enddo

    outavector(:,:)=0d0    !!     inavector(:,:)   !! PADDED

    if (wwin%topdfbasis.ge.wwin%botdfbasis) then
       call basis_transformto_local(wwin,numr,outwork,outavector)
    endif

    outavector=outavector*timefac
  
    call avectortime(2)

  end subroutine parconfigexpomult_padded0_circ

#endif

  subroutine parconfigexpomult_padded(inavector,outavector)
    use fileptrmod
    use r_parameters
    use sparse_parameters
    use ham_parameters   !! timefac
    use mpimod
    use configexpotimemod
    use configpropmod
    use configmod
    implicit none
    DATATYPE,intent(in) :: inavector(numr,dwwptr%maxdfbasisperproc)
    DATATYPE,intent(out) :: outavector(numr,dwwptr%maxdfbasisperproc)

#ifdef MPIFLAG
    select case (sparsesummaflag)
    case(0)
#endif
       call parconfigexpomult_padded0_gather(dwwptr,workconfigpointer,&
            worksparsepointerptr,inavector,outavector)
#ifdef MPIFLAG
    case(1)
       call parconfigexpomult_padded0_summa(dwwptr,workconfigpointer,&
            worksparsepointerptr,inavector,outavector)
    case(2)
       call parconfigexpomult_padded0_circ(dwwptr,workconfigpointer,&
            worksparsepointerptr,inavector,outavector)
    case default
       OFLWR "Error sparsesummaflag ",sparsesummaflag; CFLST
    end select
#endif

  end subroutine parconfigexpomult_padded
 
end subroutine expoavecprop

end module expoavecpropmod
