
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
  real*8,intent(in) :: time1,time2
  DATATYPE,intent(inout) :: inspfs(spfsize,nspf)
  integer,intent(out) :: numiters
  real*8, save :: tempstepsize = -1d0
  real*8 :: midtime, tdiff, error, norm
  integer :: itrace, iflag,getlen
  integer :: expofileptr=805
  external :: jacoperate,jacopcompact ,realpardotsub

!! lowers thisexpodim until number of internal expokit steps is 2 or less, 
!!  then goes back and sets
!!  exposet to 1.  once exposet is 1, does not decrease further but may
!!  increase.

  integer, save :: thisexpodim=0, icalled=0
  integer, save :: exposet=0   
!! 0 = not set; reduce to minimum      !! 1 = minimum found; don't decrease

  DATATYPE,allocatable :: aspfs(:,:), propspfs(:,:),   outspfs(:,:), &
       com_aspfs(:,:), com_propspfs(:,:),   com_outspfs(:,:),  tempspfs(:,:)
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
  if (lwsp.lt.0) then
     OFLWR "Oops, integer overflow, reduce maxexpodim", lwsp
     WRFL "      it's ",maxexpodim," now."; CFLST
  endif
  allocate(wsp(lwsp)); wsp=0
  liwsp=maxexpodim*2+3;   
  allocate(iwsp(liwsp)); iwsp=0

  allocate(aspfs(spfsize,nspf), propspfs(spfsize,nspf),   outspfs(spfsize,nspf), tempspfs(spfsize,nspf))
  aspfs=0; propspfs=0; outspfs=0; tempspfs=0

  if (orbcompact.ne.0) then
     allocate(com_aspfs(spfsmallsize,nspf), com_propspfs(spfsmallsize,nspf),   com_outspfs(spfsmallsize,nspf))
     com_aspfs=0; com_propspfs=0; com_outspfs=0
  endif

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


  if (myrank.eq.1.and.notiming.eq.0) then
     open(expofileptr,file=timingdir(1:getlen(timingdir)-1)//"/expo.dat",status="old", position="append")
  else
!!$ opening /dev/null multiple times not allowed :<       open(expofileptr,file="/dev/null",status="unknown")
     expofileptr=nullfileptr
  endif
  
!! nodgexpthirdflag=1 HARDWIRE 10-2015  NOT SURE ABOUT DGEXPTHIRD SUBROUTINES

  if  ((jacsymflag.ne.0).and.(jacprojorth.eq.0).and.(constraintflag.eq.0.or.jacgmatthird.ne.0).and.(nodgexpthirdflag.eq.0)) then  !! homogeneous third order.  Simpler expression for propagator.
     aspfs=inspfs

     if (parorbsplit.eq.3) then
        call dgexpthirdxxx2(idim, thisexpodim, tdiff, aspfs, inspfs, expotol, norm, &
             wsp, lwsp, iwsp, liwsp, jacoperate, itrace, iflag,expofileptr,tempstepsize,realpardotsub,idim*nprocs)
     else
        call dgexpthird(idim, thisexpodim, tdiff, aspfs, inspfs, expotol, norm, &
             wsp, lwsp, iwsp, liwsp, jacoperate, itrace, iflag,expofileptr,tempstepsize)
     endif

  else
     !! this does phi_1 = phi_0 + (exp(J t)-1)/J f(phi_0)


     call spf_linear_derivs(midtime,inspfs,aspfs) 

     if (drivingflag.ne.0) then

        call driving_linear_derivs(midtime,inspfs,tempspfs)
        aspfs(:,:)=aspfs(:,:)+tempspfs(:,:)

     endif

     call apply_spf_constraints(aspfs)

     propspfs=0.d0

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

     inspfs=inspfs+outspfs

     
  endif

  if (expofileptr.ne.nullfileptr) then
     close(expofileptr)
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
  deallocate(aspfs, propspfs, outspfs, tempspfs)
  if (orbcompact.ne.0) then
     deallocate(com_aspfs, com_propspfs, com_outspfs)
  endif

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


!! attempt for quad

subroutine jacorth(inspfs,outspfs)
  use parameters
  use jacmod
  implicit none
  DATATYPE,intent(in) ::  inspfs(spfsize,nspf)
  DATATYPE,intent(out) :: outspfs(spfsize,nspf)
  call oneminusproject(inspfs,outspfs,jacvect)
end subroutine jacorth


!! SUBROUTINE PASSED TO EXPOKIT FOR ORBITAL PROPAGATION

subroutine jacopcompact(com_inspfs,com_outspfs)
  use parameters
  implicit none
  DATATYPE,intent(in) :: com_inspfs(spfsmallsize,nspf)
  DATATYPE,intent(out) :: com_outspfs(spfsmallsize,nspf)
  DATATYPE ::  inspfs(spfsize,nspf), outspfs(spfsize,nspf)  !! AUTOMATIC

  call spfs_expand(com_inspfs,inspfs)
  call jacoperate0(1,inspfs,outspfs)
  call spfs_compact(outspfs,com_outspfs)

end subroutine jacopcompact


subroutine jacoperate(inspfs,outspfs)
  use parameters
  implicit none
  DATATYPE,intent(in) :: inspfs(spfsize,nspf)
  DATATYPE,intent(out) :: outspfs(spfsize,nspf)

  call jacoperate0(1,inspfs,outspfs)

end subroutine jacoperate

     

subroutine jacoperate0(dentimeflag,inspfs,outspfs)
  use parameters
  use jacmod
  use mpimod
  use xxxmod    !! drivingorbs.... hmmm could just make wrapper but whatever
  use linearmod
  implicit none
  integer,intent(in) :: dentimeflag
  DATATYPE,intent(in) ::  inspfs(spfsize,nspf)
  DATATYPE,intent(out) :: outspfs(spfsize,nspf)
  integer :: ii,ibot,getlen
  integer, save :: times(20), numcalledhere=0,itime,jtime
  DATATYPE ::  nulldouble(2),pots(3)
  real*8 :: facs(0:1)
  DATATYPE :: jactemp3(spfsize,nspf),   jactemp2(spfsize,nspf),&   !! AUTOMATIC
       tempspfs(spfsize,nspf),temporbs(spfsize,nspf)

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
     
     if (jacsymflag.ne.0) then
        call system_clock(itime)
           
        call project(inspfs,jactemp2,jacvect) 
           
        call system_clock(jtime); times(1)=times(1)+jtime-itime;   call system_clock(itime)
           
        call actreduced0(dentimeflag,jactime,jactemp2,nulldouble,jactemp3,ii,0,0)
           
        call system_clock(jtime); times(2)=times(2)+jtime-itime;   call system_clock(itime)
           
        outspfs=outspfs+jactemp3*facs(ii)
           
        call system_clock(jtime); times(1)=times(1)+jtime-itime
     else
           
        call system_clock(itime)
           
        call actreduced0(dentimeflag,jactime, inspfs,nulldouble, jactemp3,ii,0,0)
           
        call system_clock(jtime); times(2)=times(2)+jtime-itime;   call system_clock(itime)
           
        outspfs=outspfs+jactemp3*facs(ii)
           
        call system_clock(jtime); times(1)=times(1)+jtime-itime;      
     endif
        
     call system_clock(itime)
        
     call actreduced0(dentimeflag,jactime,inspfs,nulldouble,jactemp3,ii,0,0)
        
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
        call actreduced0(dentimeflag,jactime,jactemp3,nulldouble,jactemp2,ii,0,0)
        call system_clock(jtime); times(2)=times(2)+jtime-itime;   call system_clock(itime)
        outspfs=outspfs+jactemp2*facs(ii)
        call system_clock(jtime); times(1)=times(1)+jtime-itime;  
     endif
        
!!  CONSTRAINT!!  FORGOTTEN I GUESS !! APR 2014

     if (constraintflag.ne.0) then
        call system_clock(itime)
        call op_gmat(inspfs,jactemp2,ii,jactime,jacvect)
        outspfs=outspfs+jactemp2*facs(ii)
        if (jacgmatthird.ne.0) then
           call der_gmat(jacvect,jactemp2,ii,jactime,jacvect,inspfs)
           outspfs=outspfs+jactemp2*facs(ii)
        endif
        call system_clock(jtime); times(1)=times(1)+jtime-itime;
     endif

     if (drivingflag.ne.0) then
        call system_clock(itime)
           
        call vectdpot(jactime,velflag,pots,-1)
        temporbs(:,:)=&
             pots(1)*yyy%drivingorbsxx(:,:,ii)+&
             pots(2)*yyy%drivingorbsyy(:,:,ii)+&
             pots(3)*yyy%drivingorbszz(:,:,ii)
        call derproject(temporbs,tempspfs,jacvect,inspfs)
           
        outspfs(:,:)=outspfs(:,:)-tempspfs(:,:)*facs(ii)*timefac
        call system_clock(jtime); times(3)=times(3)+jtime-itime;
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

end subroutine jacoperate0


!! KEEPME
function checknan2(input,size)
  implicit none
  integer :: size,i
  logical :: checknan2
  DATATYPE :: input(size)
  do i=1,size
  if (input(i).eq.input(i)+1) then
    checknan2=.true.
    return
  endif
  enddo
  checknan2=.false.
end function


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
  implicit none
  integer,intent(in) :: dentimeflag
  DATATYPE,intent(in) :: inspfs(spfsize,nspf) 
  real*8,intent(in) :: thistime
  DATATYPE :: nulldouble(2)
  DATATYPE,allocatable :: jactemp(:,:)
  real*8 :: gridtime

  if (allocated.eq.0) then
     allocate(jacvect(spfsize,nspf), jacvectout(spfsize,nspf))
  endif
  allocated=1

  jactime=thistime;  jacvect=inspfs; 

!! ONLY GOOD FOR CMF, LMF !!

  call actreduced0(dentimeflag,jactime,jacvect,nulldouble,jacvectout,1,0,0)

  if (effective_cmf_linearflag.ne.0) then
     allocate( jactemp(spfsize,nspf) )
     gridtime=(jactime-firsttime)/(lasttime-firsttime) 
     if ((gridtime.lt.0.d0).or.(gridtime.gt.1)) then
        print *, "GGGRIDTIME ERR ", gridtime, jactime, firsttime, lasttime
        stop
     endif
     jacvectout=jacvectout*(1.d0-gridtime)
     call actreduced0(dentimeflag,jactime,jacvect,jacvect,jactemp,0,0,0)
     jacvectout=jacvectout+gridtime*jactemp
     deallocate(jactemp)
  endif

end subroutine jacinit0


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
  integer,intent(in) :: fileptr
  write(fileptr,'(10(A13,I15))') "Start/end",times(1)/1000,"mult",times(2)/1000,"expokit",times(3)/1000
end subroutine avectortimewrite

subroutine avectortime(which)
  use avectortimemod
  implicit none
  integer,intent(in) :: which

!! times(1) = miscellaneous.  
!! times(2) = parconfigexpomult
!! times(3) = in-between (expokit)

   call system_clock(jjtime); times(which)=times(which)+jjtime-iitime; iitime=jjtime

end subroutine


subroutine exposparseprop(www,inavector,outavector,time,imc)
  use parameters
  use configpropmod
  use walkmod
  use mpimod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: imc
  DATATYPE,intent(in) :: inavector(numr,www%firstconfig:www%lastconfig)
  DATATYPE,intent(out) :: outavector(numr,www%firstconfig:www%lastconfig)
  DATATYPE :: smallvector(numr,www%maxdfbasisperproc),smallvectorout(numr,www%maxdfbasisperproc), &
       zerovector(numr,www%maxdfbasisperproc),smallvectortemp(numr,www%maxdfbasisperproc), &
       workdrivingavecdfbasis(numr,www%botdfbasis:www%topdfbasis)    !! AUTOMATIC
  external :: parconfigexpomult_padded,realpardotsub
  real*8 :: one,time
  real*8, save :: tempstepsize=-1d0
  integer :: itrace, iflag, numsteps, numiters,expofileptr=61142, liwsp=0, lwsp=0,getlen
#ifdef REALGO
  integer, parameter :: zzz=1
#else
  integer, parameter ::   zzz=2
#endif

  DATATYPE, allocatable :: wsp(:)
  integer, allocatable :: iwsp(:)
  integer, save :: thisexpodim=-1, icalled=0,my_maxaorder=0
!! 0 = not set; reduce to minimum 
!! 1 = minimum found; don't decrease
  integer, save :: exposet=0

  if (sparseconfigflag.eq.0) then
     OFLWR "must use sparseconfigflag.ne.0 for exposparseprop"; CFLST
  endif

  call avectortimeset()

  icalled=icalled+1

  if (icalled.eq.1) then
     tempstepsize=par_timestep/littlesteps
     my_maxaorder=min(www%maxdfbasisperproc*nprocs*numr,maxaorder)
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
        open(expofileptr,file=timingdir(1:getlen(timingdir)-1)//"/avecexpo.dat",status="unknown")
        write(expofileptr,*) " Avector lanczos propagator.  Order =",thisexpodim," Aerror= ",aerror
        write(expofileptr,*);        close(expofileptr)
     endif

     open(expofileptr,file=timingdir(1:getlen(timingdir)-1)//"/avecexpo.dat",status="old", position="append")
     write(expofileptr,*) "Go Avector Lanczos.  time=", time, " Order =",thisexpodim, "step ",min(par_timestep/littlesteps,tempstepsize)
     close(expofileptr)
  endif

  if (myrank.eq.1.and.notiming.eq.0) then
     open(expofileptr,file=timingdir(1:getlen(timingdir)-1)//"/avecexpo.dat",status="old", position="append")
  else
!!$ opening /dev/null multiple times not allowed :<       open(expofileptr,file="/dev/null",status="unknown")
     expofileptr=nullfileptr
  endif
  
  one=1.d0; itrace=0 ! running mode
  
  lwsp =  ( www%maxdfbasisperproc*numr*(thisexpodim+3) + www%maxdfbasisperproc*numr + (thisexpodim+3)**2 + 5*(thisexpodim+3)**2+6+1 )
  liwsp=max(10,thisexpodim+3)
  allocate(wsp(lwsp),iwsp(liwsp));    wsp=0.d0; iwsp=0
  
  smallvector(:,:)=0; smallvectorout(:,:)=0d0

  call basis_transformto_local(www,numr,inavector(:,www%botconfig),smallvector(:,:))

  lwsp=zzz*lwsp

  call avectortime(1)

!! par_timestep is a-norm estimate, ok, whatever

  if (drivingflag.ne.0) then
     smallvectortemp(:,:)=0d0
     call parconfigexpomult_padded(smallvector,smallvectortemp)

     call basis_transformto_local(www,numr,workdrivingavec(:,www%botconfig),workdrivingavecdfbasis(:,:))
     smallvectortemp(:,1:www%topdfbasis-www%botdfbasis+1)=  smallvectortemp(:,1:www%topdfbasis-www%botdfbasis+1) + workdrivingavecdfbasis(:,:) * timefac 

     zerovector(:,:)=0d0
     call DGPHIVxxx2( zzz*www%maxdfbasisperproc*numr, thisexpodim, one, smallvectortemp, zerovector, smallvectorout, aerror, par_timestep/littlesteps, wsp, lwsp, &
           iwsp, liwsp, parconfigexpomult_padded, itrace, iflag,expofileptr,tempstepsize,realpardotsub,zzz*www%numdfbasis*numr)

      smallvectorout(:,:)=smallvectorout(:,:)+smallvector(:,:)
   else
      call DGEXPVxxx2( zzz*www%maxdfbasisperproc*numr, thisexpodim, one, smallvector, smallvectorout, aerror, par_timestep/littlesteps, wsp, lwsp, &
           iwsp, liwsp, parconfigexpomult_padded, itrace, iflag,expofileptr,tempstepsize,realpardotsub,zzz*www%numdfbasis*numr)
   endif

   call avectortime(3)

   lwsp=lwsp/zzz
   
   outavector(:,:)=0d0;   

   call basis_transformfrom_local(www,numr,smallvectorout,outavector(:,www%botconfig))

   if (www%parconsplit.eq.0) then
      call mpiallgather(outavector,www%numconfig*numr,www%configsperproc*numr,www%maxconfigsperproc*numr)
   endif

   if (expofileptr.ne.nullfileptr) then
      close(expofileptr)
   endif
   
   if (iflag .ne. 0) then
      OFLWR "Error exposparseprop: ", iflag; CFLST
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
      thisexpodim=min(floor(thisexpodim*1.2)+5, my_maxaorder)
   endif
   
   deallocate(wsp,iwsp)
   
   call avectortime(1)
   
end subroutine exposparseprop



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   ALL SUBROUTINES BELOW HERE US NO MODULES BUT MPIMOD AND PARAMETERS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! for derivative of PROJECTOR using derivative of spfs.     
!!      on call inspfs is for example jacvectout

subroutine derproject(inspfs, outspfs, prospfs, prospfderivs)
  use parameters
  implicit none
  DATATYPE, intent(in) :: inspfs(spfsize, nspf), prospfs(spfsize, nspf),  prospfderivs(spfsize, nspf)
  DATATYPE, intent(out) :: outspfs(spfsize, nspf)
  integer :: i,j
  DATATYPE :: dot
  DATATYPE ::    mydot(nspf,nspf), prodot(nspf,nspf), multdot(nspf,nspf), derdot(nspf,nspf) !! AUTOMATIC

  outspfs(:,:)=0.d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
  do i=1,nspf
     do j=1,nspf
        mydot(i,j) = dot(prospfs(:,i),inspfs(:,j),spfsize)
        derdot(i,j) = dot(prospfderivs(:,i),inspfs(:,j),spfsize)
     enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

  if (parorbsplit.eq.3) then
     call mympireduce(mydot,nspf**2)
     call mympireduce(derdot,nspf**2)
  endif

  call MYGEMM('N','N',spfsize,nspf,nspf,DATAONE,prospfs,     spfsize,derdot,nspf,DATAONE,outspfs,spfsize)
  call MYGEMM('N','N',spfsize,nspf,nspf,DATAONE,prospfderivs,spfsize,mydot, nspf,DATAONE,outspfs,spfsize)

  if (jacprojorth.ne.0) then

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

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
     do i=1,nspf
        do j=1,nspf
           prodot(i,j) = dot(prospfs(:,i),prospfderivs(:,j),spfsize) 
           prodot(i,j) = prodot(i,j) + dot(prospfderivs(:,i),prospfs(:,j),spfsize)
        enddo
     enddo
!$OMP END DO
!$OMP END PARALLEL

     if (parorbsplit.eq.3) then
        call mympireduce(prodot,nspf**2)
     endif

     call MYGEMM('N', 'N', nspf, nspf, nspf, DATAONE, prodot, nspf, mydot, nspf, DATAZERO, multdot, nspf)
     call MYGEMM('N','N',spfsize,nspf,nspf,DATANEGONE,prospfs,spfsize,multdot,nspf,DATAONE,outspfs,spfsize)

  endif

end subroutine derproject



subroutine der_gmat(inspfs, outspfs, ireduced,thistime,prospfs, prospfderivs)
  use parameters
  implicit none
  integer, intent(in) :: ireduced
  real*8, intent(in) :: thistime
  DATATYPE, intent(in) :: inspfs(spfsize, nspf), prospfs(spfsize, nspf),  prospfderivs(spfsize, nspf)
  DATATYPE, intent(out) :: outspfs(spfsize, nspf)
  integer :: i,j
  DATATYPE :: dot
  DATATYPE ::    mydot(nspf,nspf), derdot(nspf,nspf), mydot0(nspf,nspf), &  !!  AUTOMATIC
       derdot0(nspf,nspf), conmat(nspf,nspf)

  outspfs(:,:)=0.d0

  if (constraintflag.eq.0) then
     return
  endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
  do i=1,nspf
     do j=1,nspf
        mydot0(i,j) = dot(prospfs(:,i),inspfs(:,j),spfsize)
        derdot0(i,j) = dot(prospfderivs(:,i),inspfs(:,j),spfsize)
     enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

  if (parorbsplit.eq.3) then
     call mympireduce(mydot0,nspf**2)
     call mympireduce(derdot0,nspf**2)
  endif

  call getconmat(thistime,ireduced,conmat)

  call MYGEMM('N','N',nspf,nspf,nspf,DATAONE,conmat,nspf,mydot0,nspf,DATAZERO,mydot,nspf)
  call MYGEMM('N','N',nspf,nspf,nspf,DATAONE,conmat,nspf,derdot0,nspf,DATAZERO,derdot,nspf)

  call MYGEMM('N','N',spfsize,nspf,nspf,DATAONE,prospfs,     spfsize,derdot,nspf,DATAONE,outspfs,spfsize)
  call MYGEMM('N','N',spfsize,nspf,nspf,DATAONE,prospfderivs,spfsize,mydot, nspf,DATAONE,outspfs,spfsize)

end subroutine der_gmat



