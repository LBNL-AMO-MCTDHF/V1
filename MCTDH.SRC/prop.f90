
!! ALL MODULES  

!! SUBROUTINES FOR WAVE FUNCTION PROPAGATION including core prop_loop subroutine
!! prop_loop calls cmf_prop_wfn which calls cmf_prop_avector and propspfs
!! for regular forward time propagation (improvedrelaxflag=0)

#include "Definitions.INC"


module prop_loop_sub_mod
contains

!! ***************************************************************************************** !!
!!                        EXACT (VMF) PROPAGATION ROUTINE 
!! ***************************************************************************************** !!

  subroutine prop_wfn(tin, tout)
    use parameters
    use xxxmod
    use mpimod
    use orbdermod
    use configmod
    use basissubmod
    use mpisubmod
    use odexmod
    use getstuffmod
    use spfsubmod
    implicit none
    real*8,intent(in) :: tout, tin
    real*8 :: mytime,nullreal, nulldouble,gbsstepsize
    integer ::  iflag,zzz, nullint, idid, itime, jtime, time=0,  time2=0 , numiters=0,rkworkdim,&
         rkiworkdim,psilength
    real*8, allocatable :: rkwork(:)
    integer, allocatable :: rkiwork(:)
    DATATYPE,allocatable :: psivec(:)

    call myclock(itime)

    if (numfrozen.gt.0) then
       print *, "program frozen in all_derivs for VMF.";     stop
    endif
    mytime=tin

    zzz=2
#ifdef REALGO
    zzz=1
#endif

    psilength=totspfdim+tot_adim*mcscfnum
    allocate(psivec(psilength))
    psivec(tot_adim*mcscfnum+1:tot_adim*mcscfnum+totspfdim)=yyy%cmfspfs(:,0)
    if (tot_adim.gt.0) then
       psivec(1:tot_adim*mcscfnum)=RESHAPE(yyy%cmfavec(:,:,0),(/tot_adim*mcscfnum/))
    endif

    if (nprocs.gt.1.and.(parorbsplit.ne.0.or.par_consplit.ne.0)) then
       OFLWR "error, for parallel options can't do VMF."
       WRFL  "set parorbsplit=0 and par_consplit=0 in namelist &parinp"; CFLST
    endif

    if (intopt==0) then
       rkworkdim=(3+6*psilength)*zzz; rkiworkdim=100
       allocate(rkwork(rkworkdim),rkiwork(rkiworkdim))
       rkiwork(:)=0; rkwork=0
       rkwork(psilength*zzz+1)=0.00001d0
       iflag=1
       call rkf45(all_derivs, psilength*zzz, psivec(:), mytime, tout, &
            relerr, abserr, iflag, rkwork, rkiwork)
       if (iflag/=2) then
          OFLWR; WRFL "RK45F ERR   Iflag = ", iflag, abserr, rkiwork(1); CFLST
       endif
       numiters=rkiwork(1)
       deallocate(rkwork,rkiwork)
    else if (intopt.eq.1) then
       rkworkdim=20*psilength*zzz; rkiworkdim=20*psilength*zzz
       allocate(rkwork(rkworkdim),rkiwork(rkiworkdim))
       rkwork=0; rkiwork=0
       call odexq(psilength*zzz,gbs_derivs,mytime,psivec(:),tout,gbsstepsize,relerr,abserr,&
            0,dummysub,0,rkwork,rkworkdim,rkiwork,rkiworkdim,idid)
       if (.not.(idid.eq.1)) then
          OFLWR; WRFL "ERR ODEX", idid; CFLST
       endif
       numiters=rkiwork(17)
       deallocate(rkwork,rkiwork)
    else
       OFLWR "INTOPT NOT REC", intopt; CFLST
    endif
!!    OFLWR "       Number of evaluations: ", numiters;  CFL

    yyy%cmfspfs(:,0)=psivec(tot_adim*mcscfnum+1:tot_adim*mcscfnum+totspfdim)
    if (tot_adim.gt.0) then
       yyy%cmfavec(:,:,0)=RESHAPE(psivec(1:tot_adim*mcscfnum),(/tot_adim,mcscfnum/))
    endif
    deallocate(psivec)

    call myclock(jtime);  time=time+jtime-itime;  itime=jtime

    call basis_project(www,numr,yyy%cmfavec(:,:,0))

    call apply_spf_constraints(yyy%cmfspfs(:,0))

    call spf_orthogit(yyy%cmfspfs(:,0), nulldouble)

!! prevent drift.. need timings
    if (parorbsplit.ne.3) then
       call mympibcast(yyy%cmfspfs(:,0),1,totspfdim)
    endif
    if (par_consplit.eq.0) then
       call mympibcast(yyy%cmfavec(:,:,0),1,tot_adim*mcscfnum)
    endif

    call get_stuff(tout) 

    call myclock(jtime);  time2=time2+jtime-itime

! rkiwork!
!  if ((myrank.eq.1).and.(notiming.eq.0)) then
!     open(853, file=timingdir(1:getlen(timingdir))//"/vmf_prop.time.dat", &
!             status="unknown", position="append")
!     write(853,'(A3,F12.3,100I15)') "T=",  tout, time/1000, time2/1000,rkiwork(1);     
!     close(853)
!  endif

  contains
    subroutine dummysub()
    end subroutine dummysub

  end subroutine prop_wfn


  subroutine propspfs(inspfs,outspfs,tin, tout,inlinearflag,numiters)
    use parameters
    implicit none
    DATATYPE,intent(in) :: inspfs(spfsize,nspf)
    DATATYPE,intent(out) :: outspfs(spfsize,nspf)
    DATATYPE,allocatable :: tempspfs2(:,:)
    integer :: inlinearflag,numiters,jj,k
    real*8 :: tout, tin,timea,timeb

    outspfs(:,:)=inspfs(:,:)
    numiters=0

    allocate(tempspfs2(spfsize,nspf));   tempspfs2(:,:)=0

    do k=1,littlesteps
       timea=tin+(tout-tin)*(k-1)/littlesteps;     timeb=tin+(tout-tin)*k/littlesteps

       call propspfs0(outspfs,tempspfs2,timea,timeb,inlinearflag,jj)
       numiters=numiters+jj
       outspfs(:,:)=tempspfs2(:,:)
    enddo

    deallocate(tempspfs2)

  end subroutine propspfs

  subroutine propspfs0(inspfs,outspfs,tin, tout,inlinearflag,numiters)
    use parameters
    use linearmod
    use orbdermod
    use odexmod
    use expospfpropmod
    use verletmod
    use spfsubmod
    implicit none
    DATATYPE,intent(in) :: inspfs(spfsize,nspf)
    DATATYPE,intent(out) :: outspfs(spfsize,nspf)
    integer ::  numiters,iflag,zzz, nullint, idid, inlinearflag,rkworkdim,rkiworkdim
    real*8 :: tout, tin,nullreal, time1, time2, nulldouble,gbsstepsize
    real*8, allocatable :: rkwork(:)
    integer, allocatable :: rkiwork(:)

    time1=tin;     time2=tout; numiters=0

    zzz=2
#ifdef REALGO
    zzz=1
#endif

    effective_cmf_linearflag=inlinearflag     !!LINEARMOD

    outspfs(:,:)=inspfs(:,:)

    if (intopt==0) then
       rkworkdim=(3+6*totspfdim)*zzz;      rkiworkdim=100
       allocate(rkwork(rkworkdim),rkiwork(rkiworkdim))
       iflag=1    ; rkwork=0;  rkiwork=0       
       call rkf45(spf_linear_derivs, totspfdim*zzz, outspfs, time1, time2, &
            relerr, abserr, iflag, rkwork, rkiwork)
       if (iflag/=2) then
          OFLWR "RK45F ERR   Iflag = ", iflag, abserr, relerr; CFLST
       endif
       numiters=rkiwork(1)
       deallocate(rkwork,rkiwork)
    else if (intopt==1) then 
       gbsstepsize=1.d-6
       rkworkdim=20*totspfdim*zzz; rkiworkdim=20*totspfdim*zzz
       allocate(rkwork(rkworkdim),rkiwork(rkiworkdim))
       rkwork=0; rkiwork=0
       call odexq(totspfdim*zzz,gbs_linear_derivs,time1,outspfs,time2,gbsstepsize,relerr,abserr,&
            0,dummysub,0,rkwork,rkworkdim,rkiwork,rkiworkdim,idid)
       numiters=rkiwork(17)
       if (.not.(idid.eq.1)) then
          OFLWR "ERR ODEX", idid; CFLST
       endif
       deallocate(rkwork,rkiwork)
    else if (intopt.eq.3) then
       call expospfprop(time1,time2,outspfs,numiters)
    else if (intopt.eq.4) then
       call verlet(outspfs,time1,time2,numiters)
    else
       OFLWR "Intopt not recognized: ", intopt; CFLST
    endif  !! intopt

    call apply_spf_constraints(outspfs)

    call spf_orthogit(outspfs, nulldouble)

  contains
    subroutine dummysub()
    end subroutine dummysub

  end subroutine propspfs0

!! ***************************************************************************************** !!
!!                        MEAN FIELD (CMF) PROPAGATION ROUTINE 
!! ***************************************************************************************** !!

  subroutine cmf_prop_wfn(tin, tout)
    use parameters
    use linearmod
    use xxxmod
    use configmod
    use mpimod
    use mpisubmod
    use derivativemod !! conpropspfs
    use quadavecmod
    use quadspfsmod
    use configstuffmod
    use repnatmod
    use denutilmod
    use getstuffmod
    use dfconsubmod
    use fockrepsubmod
    use meansubmod
    implicit none
    real*8,intent(in) :: tout, tin
    integer,save :: times(0:20)=0
    integer, save :: xxcount=0 
    integer :: numaiters,numiters,linearflag,imc,printflag=1,getlen,qq,&
         myiostat,ii,itime,jtime
    real*8 :: time1, time2
    DATAECS :: myvalues(mcscfnum)

    firsttime=tin;   lasttime=tout

    numiters=0;   numaiters=0

    xxcount=xxcount+1

    call myclock(itime)

!!     MOVE SAVED DATA (MEAN FIELDS ETC.) ONE TIME SLOT BACK. 

!! should clean this up into derived type

    yyy%invdenmat(:,:,1) = yyy%invdenmat(:,:,0)
    yyy%denmat(:,:,1) = yyy%denmat(:,:,0)
    yyy%reducedinvrsq(:,:,1)=yyy%reducedinvrsq(:,:,0)
    yyy%reducedinvr(:,:,1)=yyy%reducedinvr(:,:,0)
    yyy%reducedr(:,:,1)=yyy%reducedr(:,:,0)
    if (nonuc_checkflag/=1) then
       yyy%reducedproderiv(:,:,1)=yyy%reducedproderiv(:,:,0)
    endif
    yyy%reducedpot(:,:,:,1) = yyy%reducedpot(:,:,:,0)
    yyy%reducedpottally(:,:,:,:,1) = yyy%reducedpottally(:,:,:,:,0)
    if (drivingflag.ne.0.and.tot_adim.gt.0) then
       yyy%drivingavectorsxx(:,:,:,1)=yyy%drivingavectorsxx(:,:,:,0)
       yyy%drivingavectorsyy(:,:,:,1)=yyy%drivingavectorsyy(:,:,:,0)
       yyy%drivingavectorszz(:,:,:,1)=yyy%drivingavectorszz(:,:,:,0)
       yyy%drivingorbsxx(:,:,1)=yyy%drivingorbsxx(:,:,0)
       yyy%drivingorbsyy(:,:,1)=yyy%drivingorbsyy(:,:,0)
       yyy%drivingorbszz(:,:,1)=yyy%drivingorbszz(:,:,0)
    endif

    if (numfrozen.ne.0) then
       yyy%frozenexchinvr(:,:,1) = yyy%frozenexchinvr(:,:,0)
       yyy%frozenexchmat(:,:,1) = yyy%frozenexchmat(:,:,0)
    endif

    if (use_fockmatrix) then
       yyy%fockmatrix(:,:,1) = yyy%fockmatrix(:,:,0)
    endif

    call assign_cptr(yyy%cptr(1),yyy%cptr(0),DATAONE)

    if (sparseopt.ne.0) then
       call assign_sptr(yyysptr(1),yyysptr(0),DATAONE,www)
    endif
    if (use_dfwalktype) then
       if (sparseopt.ne.0) then
          if (shuffle_dfwalktype) then
             call assign_sptr(yyysfdptr(1),yyysfdptr(0),DATAONE,fdww)
          else
             call assign_sptr(yyysdfptr(1),yyysdfptr(0),DATAONE,dfww)
          endif
       endif
    endif
    yyy%cmfspfs(:,1)= yyy%cmfspfs(:,0)
    if (tot_adim.gt.0) then
       yyy%cmfavec(:,:,1)= yyy%cmfavec(:,:,0)
    endif

    call myclock(jtime);  times(6)=times(6)+jtime-itime

    if ( (threshflag.ne.0.and.improvedrelaxflag.le.0).and.&
         ( ( (improvedquadflag.eq.1.or.improvedquadflag.eq.3).and.tin.ge.aquadstarttime ).or.&
         ( (improvedquadflag.ge.2).and.tin.ge.quadstarttime ) ) ) then
       OFLWR " ** SETTING IMPROVEDRELAXFLAG=1 ** "; CFL
       improvedrelaxflag=max(1,(-1)*improvedrelaxflag)
    endif

    if (improvedrelaxflag.gt.0) then

       call myclock(itime)

       if (spf_flag.ne.0) then
          if (improvedquadflag.gt.1.and.tin.ge.quadstarttime) then
             call quadspfs(yyy%cmfspfs(:,0), qq)
             numiters=numiters+qq
          else
             time1=tin;        time2=tout
             call propspfs(yyy%cmfspfs(:,1), yyy%cmfspfs(:,0), time1,time2,0,qq)
             numiters=numiters+qq
          endif
          call myclock(jtime);     times(4)=times(4)+jtime-itime;     itime=jtime
!! prevent drift
          if (parorbsplit.ne.3) then
             call mympibcast(yyy%cmfspfs(:,0),1,totspfdim)
          endif
          call myclock(jtime);     times(9)=times(9)+jtime-itime;   itime=jtime
       endif

       if (improvednatflag.ne.0) then
          call replace_withnat(printflag)
       elseif (improvedfockflag.ne.0) then
          call replace_withfock(printflag)
       endif
       call myclock(jtime);     times(7)=times(7)+jtime-itime;   itime=jtime

       call all_matel()
       call myclock(jtime);     times(1)=times(1)+jtime-itime;   itime=jtime

       if (avector_flag.ne.0) then
          if ((improvedquadflag.eq.1.or.improvedquadflag.eq.3).and.tin.ge.aquadstarttime) then
             call quadavector(yyy%cmfavec(:,:,0),qq)
             numaiters=numaiters+qq
          else
             if (followflag.ne.0) then
                call myconfigeig(yyy%cptr(0),yyy%cmfavec(:,:,0),&
                     myvalues,mcscfnum,eigprintflag,2,0d0,max(0,improvedrelaxflag-1))
             else
                call myconfigeig(yyy%cptr(0),yyy%cmfavec(:,:,0),&
                     myvalues,mcscfnum,eigprintflag,1,0d0,max(0,improvedrelaxflag-1))
             endif
          endif
          call myclock(jtime);     times(5)=times(5)+jtime-itime;    itime=jtime
!! prevent drift
          if (par_consplit.eq.0) then
             call mympibcast(yyy%cmfavec(:,:,0),1,tot_adim*mcscfnum)
          endif
          call myclock(jtime);     times(9)=times(9)+jtime-itime;    itime=jtime
       endif

       call get_allden()
       call myclock(jtime);        times(2)=times(2)+jtime-itime;     itime=jtime

       if (use_fockmatrix) then
          call get_fockmatrix()
       endif     
       call get_reducedpot()
       if (numfrozen.gt.0) then
          call get_frexchange()
       endif
       call myclock(jtime);        times(3)=times(3)+jtime-itime;       itime=jtime
     
       if (constraintflag.ne.0) then    !! probably need denmat right
          call get_constraint(tout);      
          call myclock(jtime);        times(7)=times(7)+jtime-itime
       endif

    else  !! IMPROVEDRELAX

!! ******* CMF PREDICTOR ***** !!
!!            and              !!
!! ******* LMF CORRECTOR ***** !!

       linearflag=0

       call myclock(itime)

!! PRE-PROP

       if (prepropflag.ne.0) then
          if (constraintflag.ne.0) then
             time1=tin;        time2=tout
             call conpropspfs(yyy%cmfspfs(:,1),yyy%cmfspfs(:,0), time1,time2)
             call myclock(jtime);     times(4)=times(4)+jtime-itime;   itime=jtime
!! prevent drift
             if (parorbsplit.ne.3) then
                call mympibcast(yyy%cmfspfs(:,0),1,totspfdim)
                call myclock(jtime);     times(9)=times(9)+jtime-itime;   itime=jtime
             endif
             call all_matel()
             call myclock(jtime);  times(1)=times(1)+jtime-itime;   itime=jtime
          endif

          if(avector_flag.ne.0) then
             call cmf_prop_avector_and_stuff()

             call myclock(itime)
             call get_allden()
             call myclock(jtime);  times(2)=times(2)+jtime-itime;   itime=jtime
          endif

          call end_stuff()

          linearflag=1

       endif  !! prepropflag

!! MAIN STEPS

       if (step_flag.lt.0) then
          OFLWR "what? step_flag", step_flag; CFLST
       endif

       do ii=1,step_flag
          call propspfs_and_stuff()
          call cmf_prop_avector_and_stuff()
          call get_stuff0(tout,times)
          linearflag=1
       enddo

!! POST-PROP

       if (postpropflag.eq.1.and.spf_flag.ne.0) then

          call propspfs_and_stuff()

          call myclock(itime)
          call all_matel()
          call myclock(jtime);  times(1)=times(1)+jtime-itime

          call end_stuff()

       endif

       if (postpropflag.eq.2.and.avector_flag.ne.0) then

          call cmf_prop_avector_and_stuff()

          call myclock(itime)
          call get_allden()
          call myclock(jtime);  times(2)=times(2)+jtime-itime

          call end_stuff()

       endif  !! postpropflag

    endif  !!improvedrelax


!!!!! ******** !!!!!!!
!!!!! ******** !!!!!!!
!!!!! ******** !!!!!!!


    if (myrank.eq.1.and.notiming.le.1) then
       call output_denmat(1,tout)
    endif

    if ((myrank.eq.1).and.(notiming.le.1)) then
       if (xxcount==1) then
          open(853, file=timingdir(1:getlen(timingdir))//"/cmf_prop.time.dat", &
               status="unknown",iostat=myiostat)
          call checkiostat(myiostat," opening cmf_prop_time.dat")
          write(853,'(A15,100A11)',iostat=myiostat) &
               "Time", &      !! 
               "matel", &     !! (1)
               "denmat", &    !! (2)
               "reduced",&    !! (3)
               "spfprop",&    !! (4)
               "aprop", &     !! (5)
               "advance", &   !! (6)
               "constrain", & !! (7)
               "driving",&    !! (8)
               "MPI",&        !! (9)
               "#SDERIVS",&   !!
               "#ADERIVS"     !!
          call checkiostat(myiostat," writing cmf_prop_time.dat")
          close(853)
       endif
       open(853, file=timingdir(1:getlen(timingdir))//"/cmf_prop.time.dat", &
            status="unknown", position="append",iostat=myiostat)
       call checkiostat(myiostat," opening cmf_prop_time.dat")
       write(853,'(A3,F12.3,T16, 100I11)',iostat=myiostat)  "T=", tout, &
            times(1:9)/1000, numiters,numaiters
       call checkiostat(myiostat," writing cmf_prop_time.dat")
       close(853)
    endif

  contains
    subroutine propspfs_and_stuff()
      implicit none
      if (spf_flag.ne.0) then
         call myclock(itime)
         time1=tin;        time2=tout
         call propspfs(yyy%cmfspfs(:,1),yyy%cmfspfs(:,0), time1,time2,linearflag,qq)
         numiters=numiters+qq
         call myclock(jtime);     times(4)=times(4)+jtime-itime;   itime=jtime
!! prevent drift
         if (parorbsplit.ne.3) then
            call mympibcast(yyy%cmfspfs(:,0),1,totspfdim)
            call myclock(jtime);     times(9)=times(9)+jtime-itime
         endif
      endif
    end subroutine propspfs_and_stuff

    subroutine cmf_prop_avector_and_stuff()
      implicit none
      if(avector_flag.ne.0) then
         call myclock(itime)
         do imc=1,mcscfnum
            call cmf_prop_avector(yyy%cmfavec(:,imc,1),  &
                 yyy%cmfavec(:,imc,0), linearflag,tin,tout,imc,qq)
            numaiters=numaiters+qq
         enddo
         call myclock(jtime);     times(5)=times(5)+jtime-itime;  itime=jtime
!! prevent drift
         if (par_consplit.eq.0) then
            call mympibcast(yyy%cmfavec(:,:,0),1,tot_adim*mcscfnum)
            call myclock(jtime);     times(9)=times(9)+jtime-itime
         endif
      endif
    end subroutine cmf_prop_avector_and_stuff

    subroutine end_stuff()
      implicit none
      if (constraintflag.ne.0) then
         call myclock(itime)
         call get_constraint(tout)
         call myclock(jtime); times(7)=times(7)+jtime-itime
      endif

      if (drivingflag.ne.0) then
         call myclock(itime)
         call drivingtrans(tout)
         call myclock(jtime); times(8)=times(8)+jtime-itime
      endif

      call myclock(itime)
      if (use_fockmatrix) then
         call get_fockmatrix()
      endif
      call get_reducedpot()
      if (numfrozen.gt.0) then
         call get_frexchange()
      endif
      call myclock(jtime);     times(3)=times(3)+jtime-itime
    end subroutine end_stuff

  end subroutine cmf_prop_wfn


!!!   A-VECTOR PROPAGATION ROUTINE FOR CMF

!! now just cmf and lmf.

  subroutine cmf_prop_avector(avectorin,avectorout,linearflag,time1,time2,imc,numiters)
    use parameters
    use mpisubmod
    implicit none
    DATATYPE,intent(in) :: avectorin(tot_adim)
    DATATYPE,intent(out) :: avectorout(tot_adim)
    integer,intent(in) :: linearflag,imc
    integer,intent(out) :: numiters
    real*8,intent(in) :: time1,time2
    DATATYPE :: csum
    DATATYPE,allocatable :: tempvector(:)
    integer :: k,qq
    real*8 :: timea,timeb

    if (avector_flag.eq.0) then
       if (tot_adim.gt.0) then
          avectorout(:)=avectorin(:)
       endif
       return
    endif

    allocate(tempvector(tot_adim))
    if (tot_adim.gt.0) then
       tempvector(:)=avectorin(:)
    endif
    numiters=0
    do k=1,littlesteps
       timea=time1+(time2-time1)*(k-1)/littlesteps
       timeb=time1+(time2-time1)*k/littlesteps

       call cmf_prop_avector0(tempvector,avectorout,linearflag,timea,timeb,imc,qq)
       if (tot_adim.gt.0) then
          tempvector(:)=avectorout(:)
       endif
       numiters=numiters+qq
    enddo

    if (threshflag.ne.0) then
       csum=0d0
       if (tot_adim.gt.0) then
          csum=dot(avectorout,avectorout,tot_adim)
       endif
       if (par_consplit.ne.0) then
          call mympireduceone(csum)
       endif
       if (tot_adim.gt.0) then
          avectorout(:)=avectorout(:)/sqrt(csum)
       endif
    endif

    deallocate(tempvector)

  contains
    subroutine cmf_prop_avector0(avectorin,avectorout,linearflag,time1,time2,imc,numiters)
      use linearmod
      use parameters
      use mpimod
      use xxxmod
      use configmod
      use configpropmod
      use pulsesubmod
      use configstuffmod
      implicit none
      DATATYPE,intent(in) :: avectorin(tot_adim)
      DATATYPE,intent(out) :: avectorout(tot_adim)
      integer,intent(in) :: imc,linearflag
      integer,intent(out) :: numiters
      real*8, intent(in) :: time1,time2
      DATATYPE :: sum1,sum0,pots(3)=0d0
      DATATYPE :: nullvector1(numr),nullvector2(numr)
      integer :: itime,jtime,getlen,ii,iflag,myiostat
      real*8 :: thisstep,midtime,rsum
      integer, save :: times(2)=0, icalled=0

      nullvector1(:)=0;nullvector2(:)=0

      call myclock(itime)

      if (sparseopt.eq.0) then
         call zero_cptr(workconfigpointer)
      else
         call zero_sptr(worksparsepointer,www)
         if (use_dfwalktype) then
            if (shuffle_dfwalktype) then
               call zero_sptr(workfdsparsepointer,fdww)
            else
               call zero_sptr(workdfsparsepointer,dfww)
            endif
         endif
      endif

      thisstep=time2-time1; 

!! 062714
      midtime=(time2+time1)/2d0

      iflag=0
      if (drivingflag.ne.0) then
         call vectdpot(midtime,velflag,pots,imc)
         rsum=0d0
         do ii=1,3
            rsum=rsum+abs(pots(ii))**2
         enddo
         if (rsum.eq.0d0) then
            iflag=0
         else
            iflag=1
         endif
      endif

      if (linearflag.eq.0) then

         if (sparseopt.eq.0) then
            call assign_cptr(workconfigpointer,yyy%cptr(0),thisstep*DATAONE)
         else
            call assign_sptr(worksparsepointer,yyysptr(0),thisstep*DATAONE,www)
            if (use_dfwalktype) then
               if (shuffle_dfwalktype) then
                  call assign_sptr(workfdsparsepointer,yyysfdptr(0),thisstep*DATAONE,fdww)
               else
                  call assign_sptr(workdfsparsepointer,yyysdfptr(0),thisstep*DATAONE,dfww)
               endif
            endif
         endif

         if (drivingflag.ne.0.and.tot_adim.gt.0) then
            workdrivingavec(:,:)=0d0

            if (iflag.eq.1) then
               workdrivingavec(:,:)=( yyy%drivingavectorsxx(:,:,imc,0)*pots(1) + &
                    yyy%drivingavectorsyy(:,:,imc,0)*pots(2) + &
                    yyy%drivingavectorszz(:,:,imc,0)*pots(3) )*thisstep
            endif
         endif
      else

!! MAJOR BUG 062714     midtime=(time2+time1)/2d0

         sum0=(midtime-firsttime)/(lasttime-firsttime) 
         sum1=(lasttime-midtime)/(lasttime-firsttime) 

         if (drivingflag.ne.0.and.tot_adim.gt.0) then
            workdrivingavec(:,:)=0d0
            if (iflag.eq.1) then
               workdrivingavec(:,:) = thisstep*( &
   (sum1*yyy%drivingavectorsxx(:,:,imc,1)+ sum0*yyy%drivingavectorsxx(:,:,imc,0))*pots(1) + &
   (sum1*yyy%drivingavectorsyy(:,:,imc,1)+ sum0*yyy%drivingavectorsyy(:,:,imc,0))*pots(2) + &
   (sum1*yyy%drivingavectorszz(:,:,imc,1)+ sum0*yyy%drivingavectorszz(:,:,imc,0))*pots(3) )
            endif
         endif

         if (sparseopt.eq.0) then
            call add_cptr(yyy%cptr(1),yyy%cptr(0),workconfigpointer,sum1*thisstep,sum0*thisstep)
         else
            call add_sptr(yyysptr(1),yyysptr(0),worksparsepointer,sum1*thisstep,sum0*thisstep,www)
            if (use_dfwalktype) then
               if (shuffle_dfwalktype) then
                  call add_sptr(yyysfdptr(1),yyysfdptr(0),workfdsparsepointer,sum1*thisstep,&
                       sum0*thisstep,fdww)
               else
                  call add_sptr(yyysdfptr(1),yyysdfptr(0),workdfsparsepointer,sum1*thisstep,&
                       sum0*thisstep,dfww)
               endif
            endif
         endif
      endif
      call myclock(jtime); times(1)=times(1)+jtime-itime; itime=jtime

      if (tot_adim.gt.0) then
         call myconfigprop(avectorin(:),avectorout(:),midtime,imc,numiters)
      else
         call myconfigprop(nullvector1(:),nullvector2(:),midtime,imc,numiters)
      endif
      call myclock(jtime); times(2)=times(2)+jtime-itime;

      if (notiming.eq.0.and.myrank.eq.1) then

         if (icalled.eq.0) then
            open(11766, file=timingdir(1:getlen(timingdir))//"/aprop.time.dat", &
                 status="unknown",iostat=myiostat)
            call checkiostat(myiostat," writing aprop.time.dat")
            write(11766,'(100A11)',iostat=myiostat)  "init", "prop"
            call checkiostat(myiostat," writing aprop.time.dat")
            close(11766)
         endif
         open(11766, file=timingdir(1:getlen(timingdir))//"/aprop.time.dat", &
              status="unknown", position="append",iostat=myiostat)
         call checkiostat(myiostat," writing aprop.time.dat")
         write(11766,'(100I11)',iostat=myiostat)  times(1:2)
         call checkiostat(myiostat," writing aprop.time.dat")
         close(11766)
      endif

      icalled=1

    end subroutine cmf_prop_avector0

  end subroutine cmf_prop_avector

end module prop_loop_sub_mod

module proploopmod
contains

subroutine prop_loop( starttime)
  use parameters
  use mpimod
  use xxxmod
  use configmod
  use derivativemod
  use sparsemultmod
  use basissubmod
  use savenormmod
  use mpisubmod
  use prop_loop_sub_mod
  use orbdermod
  use getstuffmod
  use loadstuffmod
  use spfsubmod
  implicit none
  integer ::  jj,flag, itime, jtime, times(20)=0, qq,imc,getlen,myiostat
  DATAECS :: thisenergy(mcscfnum), lastenergy(mcscfnum) ,thisenergyavg,&
       lastenergyavg,startenergy(mcscfnum)
  CNORMTYPE :: norms(mcscfnum)
  real*8 :: thistime, starttime, thattime,error=1d10,avecerror=1d10
  DATATYPE :: sum2,sum,drivingoverlap(mcscfnum)
  DATATYPE, allocatable :: avectorp(:),outspfs(:)

  thistime=starttime;  flag=0;    call zero_mpi_times()

  allocate(avectorp(tot_adim),outspfs(totspfdim))
  outspfs=0
  if (tot_adim.gt.0) then
     avectorp=0; 
  endif

  lastenergy(:)=1.d+3;  lastenergyavg=1.d+3

  call myclock(itime)

  do imc=1,mcscfnum
     call basis_project(www,numr,yyy%cmfavec(:,imc,0))
  enddo

  if (normboflag.ne.0) then
     OFLWR "    Enforcing BO norms due to normboflag"; CFL
     call enforce_bonorms(mcscfnum,  yyy%cmfavec(:,:,0),savenorms(:,:))
  endif

  call get_stuff(0.0d0)

  do imc=1,mcscfnum

     call sparseconfigmult(www,yyy%cmfavec(:,imc,0),avectorp,&
          yyy%cptr(0),yyysptr(0),1,1,1,0,0d0,imc)

     sum=0;     sum2=0
     if (tot_adim.gt.0) then
        sum = dot(yyy%cmfavec(:,imc,0),yyy%cmfavec(:,imc,0),tot_adim)
        sum2=dot(yyy%cmfavec(:,imc,0),avectorp,tot_adim)
     endif
     if (par_consplit.ne.0) then
        call mympireduceone(sum);        call mympireduceone(sum2)
     endif
     OFLWR "IN PROP: VECTOR NORM ",sqrt(sum);CFL
     startenergy(imc)=sum2/sum
     OFLWR "         ENERGY ", startenergy(imc); CFL

!!$  moving to main
!!$     if (normboflag.ne.0) then
!!$        OFLWR "    ... will enforce BO norms due to normboflag"; CFL
!!$        call get_bonorms(1,yyy%cmfavec(:,imc,0),savenorms(:,imc))
!!$     endif

  enddo

  if (drivingflag.ne.0) then
     OFLWR "call drivinginit"; CFL
     call drivinginit(startenergy)
     call get_stuff(0.0d0)
     OFLWR "called drivinginit"; CFL
  endif

  if (debugflag.eq.956) then
     OFLWR "Stopping due to debugflag=956"; CFLST
  endif

  if ((myrank.eq.1).and.(notiming.le.1)) then
     call system("echo -n > "//timingdir(1:getlen(timingdir))//"/abstiming.dat")
     open(853, file=timingdir(1:getlen(timingdir))//"/Main.time.dat", &
          status="unknown",iostat=myiostat)
     call checkiostat(myiostat," opening Main.time.dat")
     write(853,'(T16,100A15)')  &
          "Prop ", &     !! (1)
          "Act ", &      !! (2)
          "After ", &    !! (3)
          "Init", &      !! (4)
          "Save",&       !! (5)
          "MPI",&        !! (6)
          "Non MPI"      !! (7)
     close(853)
  endif

  call myclock(jtime)  ;  times(4)=times(4)+jtime-itime

  jj=0
  do while (flag==0)      !! BEGIN MAIN LOOP
     jj=jj+1

     call myclock(itime)

     if (save_every.ne.0.and.mod(jj,save_every).eq.0) then
        call save_vector(yyy%cmfavec(:,:,0),yyy%cmfspfs(:,0),avectoroutfile,spfoutfile)  
     endif

     if ((jj==numpropsteps) .and. (threshflag/=1)) then
        flag=1
     endif

     call checkstopfile()

     thattime=thistime+par_timestep

     call myclock(jtime)  ;  times(5)=times(5)+jtime-itime;    itime=jtime

                                 !! ************************************************************* !!
     call actionsub( thistime)   !! ACTIONS.  IF ACTION CHANGES PSI THEN CALL GET_STUFF AFTER ACTION.
                                 !! ************************************************************* !!

     call myclock(jtime)  ;     times(2)=times(2)+jtime-itime;     itime=jtime

!! norms not gotten yet, tempfix 08-2016 for odex
!     rsum=0d0
!     do imc=1,mcscfnum
!        rsum=rsum+norms(imc)**2
!     enddo
!     rsum=sqrt(rsum);     abserr=sqrt(rsum)*myrelerr
!! tempfix:
     abserr=myrelerr

     if ((cmf_flag==1)) then
        call cmf_prop_wfn(thistime, thattime)
     else
        call prop_wfn(thistime, thattime)
     endif

     do imc=1,mcscfnum
        call basis_project(www,numr,yyy%cmfavec(:,imc,0))
     enddo

     if (normboflag.ne.0) then
!!        OFLWR "    ... enforcing BO norms due to normboflag"; CFL
        call enforce_bonorms(mcscfnum,  yyy%cmfavec(:,:,0),savenorms(:,:))
     endif

     call myclock(jtime)  ;  times(1)=times(1)+jtime-itime;    itime=jtime

     do imc=1,mcscfnum

        call sparseconfigmult(www,yyy%cmfavec(:,imc,0),avectorp,&
             yyy%cptr(0),yyysptr(0),1,1,timedepexpect,0,thattime,imc)


        call basis_project(www,numr,avectorp)

        sum=0; sum2=0
        if (tot_adim.gt.0) then
           sum=dot(yyy%cmfavec(:,imc,0),avectorp(:),tot_adim)
           sum2=dot(yyy%cmfavec(:,imc,0),yyy%cmfavec(:,imc,0),tot_adim)
        endif
        if (par_consplit.ne.0) then
           call mympireduceone(sum); call mympireduceone(sum2)
        endif

        thisenergy(imc) = sum/sum2   !! ok conversion
        norms(imc)=sqrt(sum2)   !! ok conversion
        
        OFL
#ifdef ECSFLAG     
#ifdef CNORMFLAG
        write(mpifileptr,'(A3,F16.5, 2(A10, 2E18.10),E18.10)') "T= ",thattime,&
             "Energy: ", thisenergy(imc), "Norm: ", abs(norms(imc)),norms(imc)
#else
        write(mpifileptr,'(A3,F16.5, 2(A10, 2E18.10))') "T= ",thattime, &
             "Energy: ", thisenergy(imc), "Norm: ", norms(imc)
#endif
#else 
        write(mpifileptr,'(A3,F16.5, 2(A10, E18.10))') "T= ", thattime, &
             "Energy: ", thisenergy(imc), "Norm: ", norms(imc)
#endif     
        CFL

        if (drivingflag.ne.0) then
          call getdrivingoverlap(drivingoverlap,mcscfnum)

          OFL
#ifdef ECSFLAG     
           write(mpifileptr,'(A3,F16.5, 2(A10, 2E18.10))') "t= ",thattime, "DEnergy ", &
#else 
           write(mpifileptr,'(A3,F16.5, 2(A10, E18.10))') "t= ", thattime, "DEnergy ", &
#endif

    (    (drivingenergies(imc)*drivingproportion**2 + CONJUGATE(drivingenergies(imc) * drivingoverlap(imc)) + &
    drivingenergies(imc)*drivingoverlap(imc) + thisenergy(imc)*norms(imc)**2)     )/  &
    (drivingproportion**2 + drivingoverlap(imc) + CONJUGATE(drivingoverlap(imc)) + norms(imc)**2), &
    "DNorm ", (sqrt &
    (drivingproportion**2 + drivingoverlap(imc) + CONJUGATE(drivingoverlap(imc)) + norms(imc)**2))

           CFL

!! ok so < Psi(t) | H | Psi(t) > =

!!   < Psi_0 | H | Psi_0 > +          drivingenergy * drivingproportion**2
!!   < Psi_0 | H | Psi' > +          = CONJUGATE(drivingenergy) * < Psi_0 | Psi' >  !! factor in psi0 in code * drivingproportion
!!   < Psi' | H | Psi_0 > +          = drivingenergy * < Psi_0 | Psi' >  !! factor in psi0 in code * drivingproportion
!!   < Psi' | H | Psi' >               in code

!! ok so < Psi(t) | Psi(t) > =

!!   < Psi_0 | Psi_0 > +          drivingproportion**2
!!   < Psi_0 |  Psi' > +          = < Psi_0 | Psi' >    !! factor in psi_0 in code * drivingproportion
!!   < Psi' | Psi_0 > +          = < Psi' | Psi_0 >     !!  ditto   * drivingproportion
!!   < Psi' | Psi' >               in code

        endif

     enddo  !! imc

     thisenergyavg=0
     do imc=1,mcscfnum
        thisenergyavg=thisenergyavg+thisenergy(imc)/mcscfnum
     enddo
     
     if ((myrank.eq.1).and.(notiming.le.1)) then
        call system("date >> "//timingdir(1:getlen(timingdir))//"/abstiming.dat")
        open(853, file=timingdir(1:getlen(timingdir))//"/Main.time.dat", &
             status="old", position="append",iostat=myiostat)
        call checkiostat(myiostat," opening Main.time.dat")
        write(853,'(F13.3,T16,100I15)',iostat=myiostat)  thistime, &
             times(1:5)/1000, mpitime/1000, nonmpitime/1000
        call checkiostat(myiostat," writing Main.time.dat")
        close(853)
     endif

     if (debugflag.eq.42) then
        OFLWR "Stopping due to debugflag=42"; CFLST
     endif

     thistime=thattime

     if (threshflag.eq.1) then

!! 06-16 why change to actreduced
!! not spf_linear_derivs?  need exchange
!!        call actreduced0(1,0d0,yyy%cmfspfs(:,0),yyy%cmfspfs(:,0),outspfs,0,1,1)

        call spf_linear_derivs0(0,1,0d0,yyy%cmfspfs(:,0),outspfs,1,1)

        call apply_spf_constraints(outspfs)

        avecerror=abs(thisenergyavg-lastenergyavg)/par_timestep

!!$ OOPS 04-24-15        error=sqrt(abs(hermdot(outspfs,outspfs,totspfdim)))
!                        if (parorbsplit.eq.3) then
!                           call mympirealreduceone(error)
!                        endif

        error=abs(hermdot(outspfs,outspfs,totspfdim))
        if (parorbsplit.eq.3) then
           call mympirealreduceone(error)
        endif
        error=sqrt(error)

        OFL; write(mpifileptr,'(A24,2E10.2,A10,2E10.2)') &
             " STOPTEST : ORBITALS ", error,stopthresh, " AVECTOR ",avecerror,astoptol;CFL
        
        if ( error.lt.stopthresh .or.spf_flag.eq.0) then
           if (avecerror.gt.astoptol) then
              OFL; write(mpifileptr,'(A67,2E12.5)') &
                   "   Orbitals Converged, Avector not necessarily converged "; CFL
           else
              flag=1
              OFLWR; WRFL "   ***  CONVERGED *** "; WRFL
              do qq=1,mcscfnum
                 write(mpifileptr,'(T5,10F20.12)') thisenergy(qq), thisenergy(qq)-lastenergy(qq)
              enddo
              write(mpifileptr, *) "   ***   ";              write(mpifileptr, *);CFL
           endif
        endif
        do qq=1,mcscfnum
           lastenergy(qq)=thisenergy(qq)
        enddo
        lastenergyavg=thisenergyavg
        
        par_timestep=par_timestep*timestepfac
        if (par_timestep.gt.max_timestep) then
           par_timestep=max_timestep
        endif

     endif

     call myclock(jtime)  ;        times(3)=times(3)+jtime-itime

  enddo    !! END MAIN LOOP

  norms=0
  if (tot_adim.gt.0) then
     do imc=1,mcscfnum
        norms(imc)=dot(yyy%cmfavec(:,imc,0),& !! ok conversion
             yyy%cmfavec(:,imc,0),tot_adim)   !! ok conversion
     enddo
  endif
  if (par_consplit.ne.0) then
#ifndef REALGO
#ifndef CNORMFLAG     
     call mympirealreduce(norms,mcscfnum)
#else
     call mympireduce(norms,mcscfnum)
#endif
#else
     call mympireduce(norms,mcscfnum)
#endif
  endif
  norms(:)=sqrt(norms(:))

  if (threshflag.eq.1) then
     if (tot_adim.gt.0) then
        do imc=1,mcscfnum
           yyy%cmfavec(:,imc,0)=yyy%cmfavec(:,imc,0)/norms(imc)
        enddo
     endif
     norms(:)=1d0 
  endif
  
  call actions_final()

  OFLWR "   ...done prop..."; CFL
  call mpibarrier()

  if (saveflag.ne.0) then
     OFLWR "  ...saving vector..."; CFL
     call mpibarrier()
     call save_vector(yyy%cmfavec(:,:,0),yyy%cmfspfs(:,0),avectoroutfile,spfoutfile)  
  endif
  OFLWR "   ...end prop..."; CFL
  call mpibarrier()

  deallocate(avectorp,outspfs)

end subroutine prop_loop

end module proploopmod

