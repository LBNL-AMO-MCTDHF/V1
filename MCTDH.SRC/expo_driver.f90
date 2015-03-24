
#include "Definitions.INC"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! This file contains exponential propagation routines... They are separated nicely by the modules they
!!!   contain, mutually exclusively.  In order, in this file:
!!!   1) expomod files, for orbital propagation.  These are for expoprop, called from prop.f90, driver
!!!     for expokit subroutines.  These expokit subroutines for the orbital equation call 
!!!   2) jacmod files, for operation of jacobian, called by expokit, for orbital equation.
!!!   3) configexpomod files, for A-vector propagation
!!!   4) utilities that do not use any modules (besides mpimod and parameters)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   EXPOMOD SUBROUTINES, ETC !!
!!   Driver for expokit call (z/dgphiv) for orbital equation
!!   z/dghpiv uses jacoperate and jacmod subroutines
!!   (at bottom of this file; expomod are at top)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!! EXPONENTIAL PROPAGATION OF ORBITALS !!

subroutine expoprop(time1,time2,inspfs, numiters)
  use parameters
  use mpimod
  implicit none

  real*8 :: time1,time2,midtime, tdiff, error, norm
  DATATYPE :: inspfs(spfsize,nspf),tempspfs(spfsize,nspf)
  real*8, save :: tempstepsize = -1d0
  integer :: itrace, iflag,  numiters,getlen
  integer, parameter :: expofileptr=805
  external :: jacoperate,jacopcompact ,realpardotsub

!! lowers thisexpodim until number of internal expokit steps is 2 or less, 
!!  then goes back and sets
!!  exposet to 1.  once exposet is 1, does not decrease further but may
!!  increase.

  integer, save :: thisexpodim=0, icalled=0
  integer, save :: exposet=0   
!! 0 = not set; reduce to minimum      !! 1 = minimum found; don't decrease

  DATATYPE :: aspfs(spfsize,nspf), propspfs(spfsize,nspf),   outspfs(spfsize,nspf), &
       com_aspfs(spfsmallsize,nspf), com_propspfs(spfsmallsize,nspf),   com_outspfs(spfsmallsize,nspf)
  integer :: idim, liwsp=0, lwsp=0,ttott
  real*8, allocatable :: wsp(:)
  integer, allocatable :: iwsp(:)

  if (orbcompact.ne.0) then
     ttott=spfsmallsize*nspf
  else
     ttott=totspfdim
  endif

#ifdef REALGO
  idim=ttott  !!totspfdim
#else
  idim=2*ttott  !!2*totspfdim
#endif

  if (maxexpodim.gt.ttott-1) then
     OFLWR "Maxexpodim too big, resetting to", ttott-1, " totspfdim is ", totspfdim; CFL
     maxexpodim=ttott-1
  endif
  lwsp = ( idim*(maxexpodim+3) + idim + (maxexpodim+3)**2 + 5*(maxexpodim+3)**2+6+1 );  
  allocate(wsp(lwsp)); wsp=0
  liwsp=maxexpodim*2+3;   
  allocate(iwsp(liwsp)); iwsp=0


  icalled=icalled+1
  if (icalled.eq.1) then
     tempstepsize=par_timestep/littlesteps     
     if (expodim.gt.maxexpodim) then
        expodim=maxexpodim
     endif
     thisexpodim=expodim
     if ((myrank.eq.1).and.(notiming==0)) then
        open(expofileptr,file=timingdir(1:getlen(timingdir)-1)//"/expo.dat",status="unknown")
        write(expofileptr,*) " Exponential propagator output.  expotol=",expotol
        write(expofileptr,*);        close(expofileptr)
     endif
  endif

  if (mod(icalled,12).eq.0) then
     exposet=0
  endif
  numiters=0;  

  tdiff=(time2-time1)
  midtime=(time2+time1)/2d0

  if (thisexpodim.lt.maxexpodim) then
     tempstepsize=tempstepsize*4
  else
!!     tempstepsize=tempstepsize*1.1
     tempstepsize=tempstepsize*expostepfac
  endif

  if ((myrank.eq.1).and.(notiming.eq.0)) then
     open(expofileptr,file=timingdir(1:getlen(timingdir)-1)//"/expo.dat",status="old", position="append")
     write(expofileptr,*) "Go Orbital Expoprop.  Tinit=", time1, " thisexpodim=",thisexpodim, " step ", min(par_timestep/littlesteps,tempstepsize)
     close(expofileptr)
  endif

  call jacinit(inspfs,midtime)
  

  itrace=0   !! 1=print info
  norm=1.d0

  if  ((jacsymflag.ne.0).and.(jacprojorth.eq.0)) then  !! homogeneous third order.  Simpler expression for propagator.

     call noparorbsupport("dgexpthird")

     aspfs=inspfs
     call dgexpthird(idim, thisexpodim, tdiff, aspfs, inspfs, expotol, norm, &
          wsp, lwsp, iwsp, liwsp, jacoperate, itrace, iflag)
  else
     !! this does phi_1 = phi_0 + (exp(J t)-1)/J f(phi_0)


     call spf_linear_derivs(midtime,inspfs,aspfs) 

     if (drivingflag.eq.1.or.drivingflag.eq.2) then
        call driving_linear_derivs(midtime,inspfs,tempspfs)
        aspfs(:,:)=aspfs(:,:)+tempspfs(:,:)
     endif

     call apply_spf_constraints(aspfs)

     propspfs=0.d0

     if (myrank.eq.1.and.notiming.eq.0) then
        open(expofileptr,file=timingdir(1:getlen(timingdir)-1)//"/expo.dat",status="old", position="append")
     else
        open(expofileptr,file="/dev/null",status="unknown")
     endif
     if (orbcompact.ne.0) then
        call spfs_compact(aspfs,com_aspfs); call spfs_compact(propspfs,com_propspfs);
        if (parorbsplit.eq.3) then
           call dgphivxxx2(idim, thisexpodim, tdiff, com_aspfs, com_propspfs, com_outspfs, expotol, norm, &
                wsp, lwsp, iwsp, liwsp, jacopcompact, itrace, iflag,expofileptr,tempstepsize,realpardotsub,idim*nprocs)
        else
           call dgphiv(idim, thisexpodim, tdiff, com_aspfs, com_propspfs, com_outspfs, expotol, norm, &
                wsp, lwsp, iwsp, liwsp, jacopcompact, itrace, iflag,expofileptr,tempstepsize)
        endif
        call spfs_expand(com_outspfs,outspfs);
     else
        if (parorbsplit.eq.3) then
           call dgphivxxx2(idim, thisexpodim, tdiff, aspfs, propspfs, outspfs, expotol, norm, &
                wsp, lwsp, iwsp, liwsp, jacoperate, itrace, iflag,expofileptr,tempstepsize,realpardotsub,idim*nprocs)
        else
           call dgphiv(idim, thisexpodim, tdiff, aspfs, propspfs, outspfs, expotol, norm, &
                wsp, lwsp, iwsp, liwsp, jacoperate, itrace, iflag,expofileptr,tempstepsize)
        endif
     endif

     close(expofileptr)
     inspfs=inspfs+outspfs

     
  endif

  if (iflag/=0) then
     OFLWR "Expo error spf ", iflag; CFLST
  endif

  if (iwsp(4).gt.1) then
        thisexpodim=min(maxexpodim,max(thisexpodim+5,ceiling(thisexpodim*1.2)))
     exposet=1
  endif

  numiters=numiters+iwsp(1)

  call apply_spf_constraints(inspfs)

  call spf_orthogit(inspfs, error)
  
  if ((myrank.eq.1).and.(notiming.eq.0)) then
     open(expofileptr,file=timingdir(1:getlen(timingdir)-1)//"/expo.dat",status="old", position="append")
     write(expofileptr,*) "   End expo. Orthog error, steps, stepsize, iterations", error,iwsp(4),tempstepsize,iwsp(1); 
     if (abs(error).gt.1d-5) then
        write(expofileptr,*) "    **** THATSA BIG ERROR I THINK! ***"
     endif
     write(expofileptr,*);     close(expofileptr)
  endif

  if ((exposet==0)) then
     thisexpodim=max(5,floor(thisexpodim*0.95))
  endif


  deallocate(wsp,iwsp)

end subroutine expoprop


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


!! SUBROUTINE PASSED TO EXPOKIT FOR ORBITAL PROPAGATION

subroutine jacopcompact(com_inspfs,com_outspfs)
  use parameters
  implicit none

  DATATYPE ::  inspfs(spfsize,nspf), outspfs(spfsize,nspf), com_inspfs(spfsmallsize,nspf), com_outspfs(spfsmallsize,nspf)

  call spfs_expand(com_inspfs,inspfs)
  call jacoperate(inspfs,outspfs)
  call spfs_compact(outspfs,com_outspfs)

end subroutine jacopcompact

     

subroutine jacoperate(inspfs,outspfs)
  use parameters
  use jacmod
  use mpimod
  use xxxmod    !! drivingorbs.... hmmm could just make wrapper but whatever
  use linearmod
  implicit none

  integer :: ii,ibot,getlen
  integer, save :: times(20), numcalledhere=0,itime,jtime
  DATATYPE ::  inspfs(spfsize,nspf), outspfs(spfsize,nspf), nulldouble(2),pots(3)
  real*8 :: facs(0:1)
  DATATYPE,save,allocatable :: jactemp3(:,:),   jactemp2(:,:),  tempspfs(:,:),temporbs(:,:)
  integer,save :: allochere=0

  if (allochere.eq.0) then
     allocate(jactemp3(spfsize,nspf),   jactemp2(spfsize,nspf),  tempspfs(spfsize,nspf),temporbs(spfsize,nspf))
  endif
  allochere=1


  numcalledhere=numcalledhere+1

  call system_clock(itime)

  if (effective_cmf_linearflag.eq.1) then
     ibot=0
     facs(0)=(jactime-firsttime)/(lasttime-firsttime);     facs(1)=1d0-facs(0)
  else
     ibot=1;     facs(0)=0d0;     facs(1)=1d0
  endif

  outspfs=0.d0

  call system_clock(jtime); times(1)=times(1)+jtime-itime; 

!! ** term with inspfs on far right ** !!

  do ii=1,ibot,-1

     if (effective_cmf_spfflag.ne.0) then

     if (jacsymflag.ne.0) then
        call system_clock(itime)

        call project(inspfs,jactemp2,jacvect) 

        call system_clock(jtime); times(1)=times(1)+jtime-itime;   call system_clock(itime)
        
        call actreduced0(jactime,jactemp2,nulldouble,jactemp3,ii,0,0)
        
        call system_clock(jtime); times(2)=times(2)+jtime-itime;   call system_clock(itime)
        
        outspfs=outspfs+jactemp3*facs(ii)

        call system_clock(jtime); times(1)=times(1)+jtime-itime
     else

        call system_clock(itime)

        call actreduced0(jactime, inspfs,nulldouble, jactemp3,ii,0,0)

        call system_clock(jtime); times(2)=times(2)+jtime-itime;   call system_clock(itime)

        outspfs=outspfs+jactemp3*facs(ii)

        call system_clock(jtime); times(1)=times(1)+jtime-itime;      
     endif

     call system_clock(itime)

     call actreduced0(jactime,inspfs,nulldouble,jactemp3,ii,0,0)

     call system_clock(jtime); times(2)=times(2)+jtime-itime;   call system_clock(itime)

     call project(jactemp3,jactemp2,jacvect)
     outspfs=outspfs-jactemp2*facs(ii)

!! terms from projector

     call derproject(jacvectout,jactemp2,jacvect,inspfs)
     outspfs=outspfs-jactemp2*facs(ii)

     call system_clock(jtime); times(1)=times(1)+jtime-itime

     if (jacsymflag.ne.0) then
        call system_clock(itime)
        call derproject(jacvect,jactemp3,jacvect,inspfs) 
        call system_clock(jtime); times(1)=times(1)+jtime-itime;   call system_clock(itime)
        call actreduced0(jactime,jactemp3,nulldouble,jactemp2,ii,0,0)
        call system_clock(jtime); times(2)=times(2)+jtime-itime;   call system_clock(itime)
        outspfs=outspfs+jactemp2*facs(ii)
        call system_clock(jtime); times(1)=times(1)+jtime-itime;  
     endif

     endif !! spfflag



!!  CONSTRAINT!!  FORGOTTEN I GUESS !! APR 2014

     if (constraintflag.ne.0) then
        call system_clock(itime)
        call tauop(inspfs,jactemp2,ii,jactime)
        outspfs=outspfs+jactemp2*facs(ii)
        call system_clock(jtime); times(1)=times(1)+jtime-itime;
     endif

     if (effective_cmf_spfflag.ne.0) then

     if (drivingflag.eq.1.or.drivingflag.eq.3) then
        call system_clock(itime)

        call vectdpot(jactime,pots)

        temporbs(:,:)=pots(1)*yyy%drivingorbsxx(:,:,ii)+pots(2)*yyy%drivingorbsyy(:,:,ii)+pots(3)*yyy%drivingorbszz(:,:,ii)
        call derproject(temporbs,tempspfs,jacvect,inspfs)

        outspfs(:,:)=outspfs(:,:)-tempspfs(:,:)*facs(ii)*timefac
        call system_clock(jtime); times(3)=times(3)+jtime-itime;
     endif

     endif

  enddo





  if ((myrank.eq.1).and.(notiming.eq.0)) then
     if (numcalledhere==1) then
        open(8577, file=timingdir(1:getlen(timingdir)-1)//"/jacoperate.time.dat", status="unknown")
        write(8577,'(T16,100A9)') "other ", "actreduced ", "driving"; close(8577)
     endif

     if (mod(numcalledhere,timingout).eq.0) then
        open(8577, file=timingdir(1:getlen(timingdir)-1)//"/jacoperate.time.dat", status="unknown", position="append")
        write(8577,'(A3,F12.4,15I9)') "T= ", jactime,  times(1:3)/1000
        close(8577);        close(8577)
     endif
  endif


end subroutine jacoperate

!! keepme
function checknan2(input,size)
  implicit none
  integer :: size
  logical :: checknan2
  DATATYPE :: input(size)
!  do i=1,size
!  if (input(i).eq.input(i)+1) then
!    checknan2=.true.
!    return
!  endif
!  enddo
  checknan2=.false.
end function


!! keepme
function checknan2real(input,size)
  implicit none
  integer :: size
  logical :: checknan2real
  DATATYPE :: input(size)
!  do i=1,size
!  if (input(i).eq.input(i)+1) then
!    checknan2=.true.
!    return
!  endif
!  enddo
  checknan2real=.false.
end function



subroutine jacmodalloc()
  use jacmod
  use parameters
  implicit none
end subroutine jacmodalloc


!! SETS UP JACOBIAN FOR ORBITAL EXPO PROP


subroutine jacinit(inspfs, thistime) !!, timestep)
  use parameters
  use linearmod
  use jacmod
  implicit none

  DATATYPE :: inspfs(spfsize,nspf) 
  DATATYPE :: nulldouble(2), jactemp(spfsize,nspf)
  real*8 :: thistime,gridtime

  if (allocated.eq.0) then
     allocate(jacvect(spfsize,nspf), jacvectout(spfsize,nspf))
  endif
  allocated=1

  jactime=thistime;  jacvect=inspfs; 

  gridtime=(jactime-firsttime)/(lasttime-firsttime) 

  if ((gridtime.lt.0.d0).or.(gridtime.gt.1)) then
     print *, "GGGRIDTIME ERR ", gridtime, jactime, firsttime, lasttime
     stop
  endif

!! ONLY GOOD FOR CMF, LMF !!

  call actreduced0(jactime,jacvect,nulldouble,jacvectout,1,0,0)

  if (effective_cmf_linearflag.ne.0) then
     jacvectout=jacvectout*(1.d0-gridtime)
     call actreduced0(jactime,jacvect,jacvect,jactemp,0,0,0)
     jacvectout=jacvectout+gridtime*jactemp
  endif


end subroutine jacinit



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  CONFIGEXPOMOD SUBROUTINES: A-vector propagation using expokit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module avectortimemod
  implicit none
  integer :: times(100)=0, iitime= -1 , jjtime = -1
end module avectortimemod


subroutine avectortimeset()
  use avectortimemod
  implicit none
  call system_clock(iitime)
end subroutine


subroutine avectortimewrite(fileptr)
  use avectortimemod
  implicit none
  integer :: fileptr
  write(fileptr,'(10(A13,I15))') "Start/end",times(1)/1000,"mult",times(2)/1000,"expokit",times(3)/1000
end subroutine avectortimewrite


subroutine avectortime(which)
  use avectortimemod
  implicit none
  integer :: which

!! times(1) = miscellaneous.  
!! times(2) = parconfigexpomult_transpose or parconfigexpomult_transpose_spin.
!! times(3) = in-between (expokit)

   call system_clock(jjtime); times(which)=times(which)+jjtime-iitime; iitime=jjtime

end subroutine





subroutine expoconfigprop(inavector0,outavector,time)
  use parameters
  use configpropmod
  use mpimod
  implicit none
  DATATYPE :: inavector(numconfig,numr), outavector(numconfig,numr), &
       inavector0(numconfig,numr), inavectorspin(spintotrank,numr)
  external :: parconfigexpomult_transpose_spin,parconfigexpomult_transpose,realdotsub,realpardotsub
  real*8 :: one,time
  real*8, save :: tempstepsize=-1d0
  integer :: itrace, iflag, numsteps, numiters,expofileptr=61142, liwsp=0, lwsp=0,qqq,getlen,mytop,mybot,ixx,ii
#ifdef REALGO
  integer, parameter :: zzz=1
#else
  integer, parameter ::   zzz=2
#endif

  DATATYPE, allocatable :: wsp(:), smallvectortr(:,:),smallvectortrout(:,:),smallvectortrtemp(:,:),zerovectortr(:,:), &
       outavectortr(:,:),outavectorspin(:,:)
  integer, allocatable :: iwsp(:)
  integer, save :: thisexpodim=-1, icalled=0
!! 0 = not set; reduce to minimum 
!! 1 = minimum found; don't decrease
  integer, save :: exposet=0

  call avectortimeset()

  inavector(:,:)=inavector0(:,:)

  if (allspinproject.ne.0) then
     call configspin_projectall(inavector,0)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict(inavector,numr)
  endif

  icalled=icalled+1

  if (icalled.eq.1) then
     tempstepsize=par_timestep/littlesteps
     if (allspinproject.ne.0) then
!!        ii=spintotrank*numr
        ii=maxspinrank*nprocs*numr
     else
!!        ii=totadim
        ii=maxconfigsperproc*nprocs*numr
     endif
!why was this? 061814     if (maxaorder.gt.ii-1) then
!        maxaorder=ii-1
     if (maxaorder.gt.ii) then
        maxaorder=ii
     endif
     if (aorder.gt.maxaorder) then
        aorder=maxaorder
     endif
     if (aorder.lt.3) then
        OFLWR "Error, you must have like no configs, use nonsparse fool"; CFLST
     endif
     thisexpodim=aorder
  endif
  
  if (mod(icalled,20).eq.0) then
     exposet=0
  endif
  
  call configexpotimeinit(time)

  if (thisexpodim.lt.maxaorder) then
     tempstepsize=tempstepsize*4
  else
     tempstepsize=tempstepsize*1.1
  endif

  if ((myrank.eq.1).and.(notiming==0)) then
     if (icalled.eq.1) then
        open(expofileptr,file=timingdir(1:getlen(timingdir)-1)//"/avecexpo.dat",status="unknown")
        write(expofileptr,*) " Avector lanczos propagator.  Aorder =",aorder," Aerror= ",aerror
        write(expofileptr,*);        close(expofileptr)
     endif

     open(expofileptr,file=timingdir(1:getlen(timingdir)-1)//"/avecexpo.dat",status="old", position="append")
     write(expofileptr,*) "Go Avector Lanczos.  time=", time, " thisexpodim=",thisexpodim, "step ",min(par_timestep/littlesteps,tempstepsize)
     close(expofileptr)
  endif

  if (myrank.eq.1.and.notiming.eq.0) then
     open(expofileptr,file=timingdir(1:getlen(timingdir)-1)//"/avecexpo.dat",status="old", position="append")
  else
     open(expofileptr,file="/dev/null",status="unknown")
  endif


  one=1.d0; itrace=0 ! running mode

!!$  #ifdef REALGO     
!!$    OFLWR "REPROGRAM REAL VALUED PROP SORRY"; CFLST
!!$  !  if (nprocs.gt.1) then
!!$  !     OFLWR "parconfigsplit not implemented for real valued (mctdhf_diatom, _atom)"; CFLST
!!$  !  endif
!!$  !  lwsp =  ( totadim*(thisexpodim+3) + totadim + (thisexpodim+3)**2 + 5*(thisexpodim+3)**2+6+1 )
!!$  !  liwsp=thisexpodim+3
!!$  !  allocate(wsp(lwsp),iwsp(liwsp));    wsp=0.d0; iwsp=0
!!$  !
!!$  !!! par_timestep is a-norm estimate, ok, whatever
!!$  !  call EXPSPARSE( totadim, thisexpodim, one, inavector, outavector, aerror, par_timestep/littlesteps, wsp, lwsp, iwsp, liwsp, configexpomult, itrace, iflag,expofileptr)
!!$  #else
!!$  




  if (allspinproject.ne.0) then
     ixx=maxspinrank; mybot=spinstart; mytop=spinend ; qqq=spintotrank
  else
     ixx=maxconfigsperproc; mybot=botwalk; mytop=topwalk ; qqq=numconfig
  endif


  allocate(smallvectortr(numr,ixx),smallvectortrout(numr,ixx), &
       zerovectortr(numr,ixx),smallvectortrtemp(numr,ixx), outavectortr(numr,qqq), outavectorspin(qqq,numr))

   lwsp =  ( ixx*numr*(thisexpodim+3) + ixx*numr + (thisexpodim+3)**2 + 5*(thisexpodim+3)**2+6+1 )
   liwsp=max(10,thisexpodim+3)
   allocate(wsp(lwsp),iwsp(liwsp));    wsp=0.d0; iwsp=0

   smallvectortr(:,:)=0; smallvectortrout(:,:)=0d0

   if (allspinproject.ne.0) then

      call configspin_transformto(numr,inavector,inavectorspin)

      smallvectortr(:,1:spinend-spinstart+1)=TRANSPOSE(inavectorspin(spinstart:spinend,:))

   else
      smallvectortr(:,1:topwalk-botwalk+1)=TRANSPOSE(inavector(botwalk:topwalk,:))
   endif

   lwsp=zzz*lwsp

   call avectortime(1)

   
   if (drivingflag.ne.0) then
      zerovectortr(:,:)=0d0
      if (allspinproject.ne.0) then
         call parconfigexpomult_transpose_spin(smallvectortr,smallvectortrtemp)

         smallvectortrtemp(:,1:spinend-spinstart+1)=  smallvectortrtemp(:,1:spinend-spinstart+1) + TRANSPOSE(workdrivingavecspin(spinstart:spinend,:)) * timefac 

!! par_timestep is a-norm estimate, ok, whatever

         call DGPHIVxxx2( zzz*ixx*numr, thisexpodim, one, smallvectortrtemp, zerovectortr, smallvectortrout, aerror, par_timestep/littlesteps, wsp, lwsp, &
              iwsp, liwsp, parconfigexpomult_transpose_spin, itrace, iflag,expofileptr,tempstepsize,realpardotsub,zzz*qqq*numr)
         
      else
         call parconfigexpomult_transpose(smallvectortr,smallvectortrtemp)

         smallvectortrtemp(:,1:topwalk-botwalk+1)=  smallvectortrtemp(:,1:topwalk-botwalk+1) + TRANSPOSE(workdrivingavec(botwalk:topwalk,:)) * timefac 

         call DGPHIVxxx2( zzz*ixx*numr, thisexpodim, one, smallvectortrtemp, zerovectortr, smallvectortrout, aerror, par_timestep/littlesteps, wsp, lwsp, &
              iwsp, liwsp, parconfigexpomult_transpose, itrace, iflag,expofileptr,tempstepsize,realpardotsub,zzz*qqq*numr)
         
      endif

      smallvectortrout(:,:)=smallvectortrout(:,:)+smallvectortr(:,:)
      
   else

      if (allspinproject.ne.0) then

         call DGEXPVxxx2( zzz*ixx*numr, thisexpodim, one, smallvectortr, smallvectortrout, aerror, par_timestep/littlesteps, wsp, lwsp, &
              iwsp, liwsp, parconfigexpomult_transpose_spin, itrace, iflag,expofileptr,tempstepsize,realpardotsub,zzz*qqq*numr)

      else
         
         call DGEXPVxxx2( zzz*ixx*numr, thisexpodim, one, smallvectortr, smallvectortrout, aerror, par_timestep/littlesteps, wsp, lwsp, &
              iwsp, liwsp, parconfigexpomult_transpose, itrace, iflag,expofileptr,tempstepsize,realpardotsub,zzz*qqq*numr)

      endif

   endif

   call avectortime(3)


   lwsp=lwsp/zzz
   
   outavectortr(:,:)=0d0;   
   outavectortr(:,mybot:mytop)=smallvectortrout(:,1:mytop-mybot+1)

   if (allspinproject.ne.0) then
      call mpiallgather(outavectortr,spintotrank*numr,allspinranks*numr,maxspinrank*numr)

      outavectorspin(:,:)=TRANSPOSE(outavectortr(:,:))

      call configspin_transformfrom(numr,outavectorspin,outavector)

   else
      call mpiallgather(outavectortr,numconfig*numr,configsperproc*numr,maxconfigsperproc*numr)
      outavector(:,:)=TRANSPOSE(outavectortr(:,:))
   endif

   deallocate(smallvectortr,smallvectortrout,  zerovectortr,smallvectortrtemp, outavectortr, outavectorspin)

!!$   #endif

  if (allspinproject.ne.0) then
     call configspin_projectall(outavector,0)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict(outavector,numr)
  endif


  close(expofileptr)

  if (iflag .ne. 0) then
     OFLWR "Error expoconfigprop: ", iflag; CFLST
  endif
  numsteps=iwsp(4);  numiters=iwsp(1)

  if ((exposet==0).and.(numsteps.eq.1).and.(thisexpodim.gt.2)) then
     thisexpodim=max(2,thisexpodim-1)
  endif

  if (myrank.eq.1.and.notiming.eq.0) then
     open(expofileptr,file=timingdir(1:getlen(timingdir)-1)//"/avecexpo.dat",status="old", position="append")
     write(expofileptr,*) "   End avector prop.  steps, iterations, stepsize", iwsp(4),iwsp(1),tempstepsize; 
     write(expofileptr,*);    
     call avectortimewrite(expofileptr)
 close(expofileptr)
  endif

  if (numsteps.gt.1) then
     exposet=1
     thisexpodim=min(floor(thisexpodim*1.2)+5, maxaorder)
  endif

  deallocate(wsp,iwsp)

   call avectortime(1)

end subroutine expoconfigprop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   ALL SUBROUTINES BELOW HERE US NO MODULES BUT MPIMOD AND PARAMETERS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! for derivative of PROJECTOR using derivative of spfs.     
!!      on call inspfs is for example jacvectout

subroutine derproject(inspfs, outspfs, prospfs, prospfderivs)
  use parameters
  implicit none
  DATATYPE :: inspfs(spfsize, nspf), outspfs(spfsize, nspf),  prospfs(spfsize, nspf), prospfderivs(spfsize, nspf)
  call derproject0(inspfs, outspfs, prospfs, prospfderivs,0)  ! no conjg
end subroutine derproject


subroutine derproject0(inspfs, outspfs, prospfs, prospfderivs,inconjgflag)
  use parameters
  implicit none

  DATATYPE :: inspfs(spfsize, nspf), outspfs(spfsize, nspf), prospfs(spfsize, nspf),  prospfderivs(spfsize, nspf)
  integer :: i,j, inconjgflag
  DATATYPE :: dot
#ifdef CNORMFLAG
  DATATYPE :: hermdot
#endif
  integer,save :: allochere=0
  DATATYPE,save,allocatable ::  mydot(:,:), prodot(:,:), multdot(:,:), derdot(:,:),mydot2(:,:),prodot2(:,:),multdot2(:,:)

  if (allochere.eq.0) then
     allocate(       mydot(nspf,nspf), prodot(nspf,nspf), multdot(nspf,nspf), derdot(nspf,nspf),mydot2(nspf,nspf),prodot2(nspf,nspf),multdot2(nspf,nspf))
  endif
  allochere=1


  outspfs(:,:)=0.d0
  do i=1,nspf
     do j=1,nspf
        mydot(i,j) = dot(prospfs(:,i),inspfs(:,j),spfsize)

#ifdef CNORMFLAG
        if (inconjgflag==1) then
           derdot(i,j) = conjg(hermdot(prospfderivs(:,i),inspfs(:,j),spfsize))
        else
#endif
           derdot(i,j) = dot(prospfderivs(:,i),inspfs(:,j),spfsize)
#ifdef CNORMFLAG
        endif
#endif

     enddo
  enddo

  if (parorbsplit.eq.3) then
     call mympireduce(mydot,nspf**2)
     call mympireduce(derdot,nspf**2)
  endif


  if (inconjgflag==1) then
     mydot=ANTICON(mydot)    !! in derprojectconjg call, this is herm conjg of that in derproject call
  endif
  do i=1,nspf
     if (inconjgflag==1) then
        do j=1,nspf
           outspfs(:,i) = outspfs(:,i) + ANTICON(prospfs(:,j)) *                derdot(j,i)
        enddo
     else
        do j=1,nspf
           outspfs(:,i) = outspfs(:,i) + prospfs(:,j) *           derdot(j,i)
        enddo
     endif
     do j=1,nspf
        outspfs(:,i) = outspfs(:,i) + (prospfderivs(:,j)) *             mydot(j,i)
     enddo
  enddo

  if (jacprojorth.ne.0) then

     call noparorbsupport("jacprojorth dude!")

!! Proj in always-orthogonal-derivative form,
!!
!!  P = sum_ij | prospf_i > (S^-1)_ij < prospf_j |
!!
!!  where S=delta_ij = <jacvect_i | jacvect_j>
!!
!!  (dS)_ij = <dphi_i | jacvect_j> + <jacvect_i | dphi_j>
!!
!!  (dS^-1)_ij = - <dphi_i | jacvect_j> - <jacvect_i | dphi_j>  at S=1
!!

!        prodot is     (pro/proder,pro/proder)
!        mydot2 is     (pro,in)
!        conjg: want  (in/proder, in/proder)
!        and mydot2 is  (in,pro)

     prodot2=0.d0
     do i=1,nspf
        do j=1,nspf
           if (inconjgflag.eq.1) then
              prodot(i,j) = dot(ANTICON(inspfs(:,i)),prospfderivs(:,j),spfsize) 
              prodot2(i,j) = (dot(ALLCON(inspfs(:,i)),CONJUGATE(prospfderivs(:,j)),spfsize))
              mydot2(i,j) = dot(ALLCON(prospfs(:,i)),ALLCON(inspfs(:,j)),spfsize)
           else
              prodot(i,j) = dot(prospfs(:,i),prospfderivs(:,j),spfsize) 
              prodot(i,j) = prodot(i,j) + dot(prospfderivs(:,i),prospfs(:,j),spfsize)
           endif
        enddo
     enddo

!prodot  (pro,der) x mydot (pro,in)   : pro is in;    in is out
!prodot2  (der,pro) x mydot (pro,in)   : der is in;    in is out
!
! conjg : given in as pro=in,in=pro. mydot already cc.
!
!  mydot (ANTI(in),ANTI(pro))  x prodot (ANTI(in),der)  
!  prodot2(der,ALC(in))    x  mydot (ANTI(in),ANTI(pro))
! 
! multiply prodot -> multdot(pro,in)  by gemm(n,n,prom,my,mult)
!
! conjg:  mydot(in,pro)* prodot (pro,der)   der is out p
!         
     if (inconjgflag.eq.1) then
!! first term like the second, but having ANTICONNED all but prospfderivs
        call MYGEMM('N', 'N', nspf, nspf, nspf, DATANEGONE, prodot, nspf, mydot, nspf, DATAZERO, multdot, nspf)

        call MYGEMM('N', 'N', nspf, nspf, nspf, DATANEGONE, prodot2, nspf, mydot2, nspf, DATAZERO, multdot2, nspf)
     else
        call MYGEMM('N', 'N', nspf, nspf, nspf, DATANEGONE, prodot, nspf, mydot, nspf, DATAZERO, multdot, nspf)
     endif
     do i=1,nspf
        do j=1,nspf
           if (inconjgflag.eq.1) then
              outspfs(:,i) = outspfs(:,i) + ANTICON(inspfs(:,j)) * multdot(j,i)
              outspfs(:,i) = outspfs(:,i) + ANTICON(inspfs(:,j)) * (multdot2(i,j))
           else
              outspfs(:,i) = outspfs(:,i) + prospfs(:,j) * multdot(j,i)        !!multdot 
           endif
        enddo
     enddo
  endif

end subroutine derproject0

