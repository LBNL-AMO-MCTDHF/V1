  
#include "Definitions.INC"

!! SUBROUTINES FOR WAVE FUNCTION PROPAGATION
!! THESE ROUTINES  do not have any dummy psi(:)'s .  They use yyy%cmfpsivec(:,0) and for cmf, the others yyy%cmfpsivec(:,0:).
    
!! at the start of cmf_prop_wfn and prop_wfn the convention is that everything must be initialized to the proper values at that time.  
  
subroutine prop_loop( starttime)
  use parameters
  use mpimod
  use xxxmod
  use configmod
  implicit none

  integer ::  jj,flag,  iii, itime, jtime, times(20)=0, qq,imc,getlen
  DATAECS :: thisenergy(mcscfnum), lastenergy(mcscfnum) ,thisenergyavg,lastenergyavg,startenergy(mcscfnum)
  CNORMTYPE :: norms(mcscfnum)
  real*8 :: thistime, starttime, thattime,error=1d10,rsum,avecerror=1d10
  DATATYPE :: dot,  sum2,sum, psip(psilength),outspfs(totspfdim),hermdot,drivingoverlap(mcscfnum)

  thistime=starttime;  flag=0

!!GOING TO FULL ORDER.
!  if (improvedrelaxflag.ne.0.and.spf_fl ag.ne.0.and.sparseconfigflag.ne.0) then
!     lanagain=1   !! max number of builds of krylov space for eigen (#restarts + 1)
!  else
     lanagain=-1
!  endif

  lastenergy(:)=1.d+3;  lastenergyavg=1.d+3

  if (skipflag.eq.0) then
     do imc=1,mcscfnum
        call basis_project(www,numr,yyy%cmfpsivec(astart(imc),0))
     enddo

!!$ 06-2015
!!$     call all_matel()      ! initial
!!$     if (drivingflag.ne.0.and.1==1) then
!!$        OFLWR "call drivinginit"; CFL
!!$        call drivinginit(startenergy)
!!$        OFLWR "called drivinginit"; CFL
!!$     endif

     call get_stuff(0.0d0)  !initial

!!$ 06-2015
     if (drivingflag.ne.0.and.1==1) then
        OFLWR "call drivinginit"; CFL
        call drivinginit(startenergy)
        OFLWR "called drivinginit"; CFL
        call get_stuff(0.0d0)  ! for actions_initial
     endif

     if ((myrank.eq.1).and.(notiming.le.1)) then
        call system("echo -n > "//timingdir(1:getlen(timingdir)-1)//"/abstiming.dat")
        !! times(1:5)
        open(853, file=timingdir(1:getlen(timingdir)-1)//"/Main.time.dat", status="unknown")
        write(853,'(T16,100A15)')  "Spfs ", "Prop ", "Act ", "Final ", "MPI", "Non MPI";        close(853)
     endif
  endif

  call system_clock(itime)
  call actions_initial()                 !!   ACTIONS_INITIAL: PUT ANALYSIS ROUTINES HERE; THEN TERMINATES.
                                         !! ***************************************************************
  call system_clock(jtime)  ;  times(3)=times(3)+jtime-itime

  if (skipflag.ne.0) then
     OFLWR "Stopping"; CFLST
  endif

!! 06-2015 ADDING THIS HERE
  if (improvedrelaxflag.ne.0.and.improvednatflag.ne.0) then
     call system_clock(itime)
     call replace_withnat(1)
     call system_clock(jtime);        times(7)=times(7)+jtime-itime;
  endif
!! allowing that actions_initial could change psi.  OR replace_withnat!  drivinginit ok.

  call get_stuff(0.d0)

!! 06-2015 moved this down here
  do imc=1,mcscfnum

     call sparseconfigmult(www,yyy%cmfpsivec(astart(imc),0),psip(astart(imc)),yyy%cptr(0),yyy%sptr(0),1,1,1,0,0d0)

     sum = dot(     yyy%cmfpsivec(astart(imc),0),yyy%cmfpsivec(astart(imc),0),tot_adim)
     if (par_consplit.ne.0) then
        call mympireduceone(sum)
     endif
     OFLWR "IN PROP: VECTOR NORM ",sqrt(sum);CFL

     sum2=dot(     yyy%cmfpsivec(astart(imc),0),psip(astart(imc)),tot_adim)
     if (par_consplit.ne.0) then
        call mympireduceone(sum2)
     endif

     startenergy(imc)=sum2/sum
     OFLWR "         ENERGY ", startenergy(imc); CFL
     
  enddo

  if (debugflag.eq.956) then
     OFLWR "Stopping due to debugflag=956"; CFLST
  endif

  jj=0
  do while (flag==0)
     jj=jj+1

     if (save_every.ne.0.and.mod(jj,save_every).eq.0) then
        call save_vector(yyy%cmfpsivec(:,0),avectoroutfile,spfoutfile)  
     endif

     if ((jj==numpropsteps) .and. (threshflag/=1)) then
        flag=1
     endif
     open(222,file="stop",status="old", iostat=iii)
     if (iii==0) then
        close(222)
     endif
     if (iii==0) then
        flag=1
        call openfile();        write(mpifileptr, *) "Stopping due to stopfile!";CFLST
     endif

     thattime=thistime+par_timestep
     call system_clock(itime)  

                                  !! *****************************************************************************
                                  !!    ACTIONS.
     call actionsub( thistime)    !! actions may not change yyy%cmfpsivec   !! what? screw it, trying action_replacenat
                                  !! *****************************************************************************
     call system_clock(jtime)  
     times(3)=times(3)+jtime-itime

!!! (Not used for exponential propagation default - abserr thus myrelerr not used then)

     rsum=0d0
     do imc=1,mcscfnum
        rsum=rsum+norms(imc)**2
     enddo
     rsum=sqrt(rsum);     abserr=sqrt(rsum)*myrelerr

     call system_clock(itime)  

     if ((cmf_flag==1)) then
        call cmf_prop_wfn(thistime, thattime)
     else
        call prop_wfn(thistime, thattime)
     endif

     do imc=1,mcscfnum
        call basis_project(www,numr,yyy%cmfpsivec(astart(imc),0))
     enddo

     call system_clock(jtime)  ;     times(2)=times(2)+jtime-itime;     call system_clock(itime)  

     do imc=1,mcscfnum

        call sparseconfigmult(www,yyy%cmfpsivec(astart(imc),0),psip(astart(imc)),yyy%cptr(0),yyy%sptr(0),1,1,timedepexpect,0,thattime)


        call basis_project(www,numr,psip(astart(imc)))

        sum=dot(yyy%cmfpsivec(astart(imc),0),psip(astart(imc)),tot_adim)
        sum2=dot(yyy%cmfpsivec(astart(imc),0),yyy%cmfpsivec(astart(imc),0),tot_adim)
        if (par_consplit.ne.0) then
           call mympireduceone(sum); call mympireduceone(sum2)
        endif

        thisenergy(imc) = sum/sum2   !! ok implicit pmctdh
        norms(imc)=sqrt(sum2)   !! ok implicit p/chmctdh
        
!        if (threshflag.eq.1) then
!           yyy%cmfpsivec(astart(imc):aend(imc),0)=yyy%cmfpsivec(astart(imc):aend(imc),0)/norms(imc)
!           norms(imc)=1.d0
!        endif


        OFL
#ifdef ECSFLAG     
        write(mpifileptr,'(A3,F16.5, 2(A10, 2E18.10))') "T= ",thattime, "Energy: ", thisenergy(imc), "Norm: ", norms(imc)
#else 
        write(mpifileptr,'(A3,F16.5, 2(A10, E18.10))') "T= ", thattime, "Energy: ", thisenergy(imc), "Norm: ", norms(imc)
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

           CFL
        endif

        call system_clock(jtime)  ;        times(4)=times(4)+jtime-itime
     enddo  !! imc

     thisenergyavg=0
     do imc=1,mcscfnum
        thisenergyavg=thisenergyavg+thisenergy(imc)/mcscfnum
     enddo
     
     if ((myrank.eq.1).and.(notiming.le.1)) then
        call system("date >> "//timingdir(1:getlen(timingdir)-1)//"/abstiming.dat")
        open(853, file=timingdir(1:getlen(timingdir)-1)//"/Main.time.dat", status="old", position="append")
        write(853,'(F13.3,T16,100I15)')  thistime, times(1:4)/1000, mpitime/1000, nonmpitime/1000;     close(853)
     endif

     if (debugflag.eq.42) then
        OFLWR "Stopping due to debugflag=42"; CFLST
     endif

     thistime=thattime

     if (threshflag.eq.1) then

        call actreduced0(0d0,yyy%cmfpsivec(spfstart,0),yyy%cmfpsivec(spfstart,0),outspfs,0,1,1)

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
        
        if ( real(thisenergyavg).gt.real(lastenergyavg).and.lanagain.ne.-1 ) then
           OFLWR "  !! Incrementing lanagain !! ", lanagain, lanagain*2; CFL
           lanagain = lanagain * 2
           if (lanagain.gt.128) then
              lanagain=-1
           endif
        endif
        if (( error.lt.stopthresh.or.spf_flag.eq.0) .and. lanagain.ne.-1) then
           OFLWR "  !! --> Orbitals near converged, going to full A-vector order"; CFL
           lanagain=-1
        else if ( error.lt.stopthresh .or.spf_flag.eq.0) then
           if (avecerror.gt.astoptol) then
              OFL; write(mpifileptr,'(A67,2E12.5)') "   Orbitals Converged, Avector not necessarily converged "; CFL
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
  enddo
  do imc=1,mcscfnum
     norms(imc)=dot(yyy%cmfpsivec(astart(imc),0),& !! ok implicit
          yyy%cmfpsivec(astart(imc),0),tot_adim)   !! ok implicit
  enddo
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
     do imc=1,mcscfnum
        yyy%cmfpsivec(astart(imc):aend(imc),0)=yyy%cmfpsivec(astart(imc):aend(imc),0)/norms(imc)
        norms(imc)=1.d0
     enddo
  endif
  
  call actions_final()

  OFLWR "   ...done prop..."; CFL
  call mpibarrier()

  if (saveflag.ne.0) then
     OFLWR "  ...saving vector..."; CFL
     call mpibarrier()
     call save_vector(yyy%cmfpsivec(:,0),avectoroutfile,spfoutfile)  
  endif
  OFLWR "   ...end prop..."; CFL
  call mpibarrier()

  if (improvedrelaxflag.ne.0) then
     lanagain=-1 !! reset to original value -- should already be -1 in all cases, redundant
  endif

end subroutine prop_loop


!! ***************************************************************************************** !!
!!                        EXACT (VMF) PROPAGATION ROUTINE 
!! ***************************************************************************************** !!

subroutine prop_wfn(tin, tout)
  use parameters
  use xxxmod
  use mpimod
  implicit none

  real*8, external :: all_derivs, gbs_derivs, dummysub 
  real*8 :: tout, tin, mytime,nullreal, nulldouble,gbsstepsize
  integer ::  iflag,zzz, nullint, idid, itime, jtime, time=0,  time2=0 , numiters=0,rkworkdim,rkiworkdim
  real*8, allocatable :: rkwork(:)
  integer, allocatable :: rkiwork(:)

  if (numfrozen.gt.0) then
     print *, "program frozen in all_derivs for VMF.";     stop
  endif
  mytime=tin

  zzz=2
#ifdef REALGO
  zzz=1
#endif

  if (intopt==0) then
     rkworkdim=(3+6*psilength)*zzz; rkiworkdim=100
     allocate(rkwork(rkworkdim),rkiwork(rkiworkdim))
     rkiwork(:)=0
     rkwork(psilength*zzz+1)=0.00001d0
     iflag=1
     call rkf45(all_derivs, psilength*zzz, yyy%cmfpsivec(:,0), mytime, tout, relerr, abserr, iflag, rkwork, rkiwork)
     if (iflag/=2) then
        OFLWR; WRFL "RK45F ERR   Iflag = ", iflag, abserr, rkiwork(1); CFLST
     endif
     numiters=rkiwork(1)
     deallocate(rkwork,rkiwork)
  else if (intopt.eq.1) then
     rkworkdim=20*psilength*zzz; rkiworkdim=20*psilength*zzz
     allocate(rkwork(rkworkdim),rkiwork(rkiworkdim))
     call odex(psilength*zzz,gbs_derivs,mytime,yyy%cmfpsivec(:,0),tout,gbsstepsize,relerr,abserr,0,dummysub,0,rkwork,rkworkdim,rkiwork,rkiworkdim,nullreal,nullint,idid)
     if (.not.(idid.eq.1)) then
       OFLWR; WRFL "ERR ODEX", idid; CFLST
     endif
     numiters=rkiwork(17)
     deallocate(rkwork,rkiwork)
  else
     OFLWR "INTOPT NOT REC", intopt; CFLST
  endif
  OFLWR "       Number of evaluations: ", numiters;  CFL

  call system_clock(jtime);  time=time+jtime-itime;  call system_clock(itime)

  call apply_spf_constraints(yyy%cmfpsivec(spfstart,0))

!! when would you have onlyspfflag with more than orbital?  with one this just norms; don't want to do that
  call spf_orthogit(yyy%cmfpsivec(spfstart,0), nulldouble)
  call get_stuff(tout) 

  call system_clock(jtime);  time2=time2+jtime-itime

! rkiwork!
!  if ((myrank.eq.1).and.(notiming.eq.0)) then
!     open(853, file=timingdir(1:getlen(timingdir)-1)//"/vmf_prop.time.dat", status="unknown", position="append")
!     write(853,'(A3,F12.3,100I15)') "T=",  tout, time/1000, time2/1000,rkiwork(1);     close(853)
!  endif

end subroutine prop_wfn


subroutine propspfs(inspfs,outspfs,tin, tout,inlinearflag,inspfflag,numiters)
  use parameters
  implicit none
  DATATYPE :: inspfs(spfsize,nspf),outspfs(spfsize,nspf), tempspfs2(spfsize,nspf), tempspfs(spfsize,nspf)
  integer :: inlinearflag,numiters,jj,k,inspfflag
  real*8 :: tout, tin,timea,timeb

  if (spf_flag.eq.0.and.constraintflag.eq.0) then     !! no tau
     outspfs(:,:)=inspfs(:,:)
     return
  endif

  numiters=0
  tempspfs(:,:)=inspfs(:,:)

  do k=1,littlesteps
     timea=tin+(tout-tin)*(k-1)/littlesteps;     timeb=tin+(tout-tin)*k/littlesteps

     call propspfs0(tempspfs,tempspfs2,timea,timeb,inlinearflag,inspfflag,jj)
     numiters=numiters+jj
     tempspfs(:,:)=tempspfs2(:,:)
  enddo
  outspfs(:,:)=tempspfs2(:,:)

end subroutine propspfs



subroutine propspfs0(inspfs,outspfs,tin, tout,inlinearflag,inspfflag,numiters)
  use parameters
  use linearmod
  implicit none

  DATATYPE :: inspfs(spfsize,nspf),outspfs(spfsize,nspf)
  integer ::  numiters,iflag,zzz, nullint, idid, inlinearflag,inspfflag,rkworkdim,rkiworkdim
  real*8, external ::   spf_linear_derivs , gbs_linear_derivs , dummysub
  real*8 :: tout, tin,nullreal, time1, time2, nulldouble,gbsstepsize
  real*8, allocatable :: rkwork(:)
  integer, allocatable :: rkiwork(:)

  time1=tin;     time2=tout; numiters=0

  zzz=2
#ifdef REALGO
  zzz=1
#endif

  effective_cmf_linearflag=inlinearflag     !!LINEARMOD

  effective_cmf_spfflag=inspfflag

  outspfs(:,:)=inspfs(:,:)

  if (intopt==0) then
     rkworkdim=(3+6*totspfdim)*zzz;      rkiworkdim=20
     allocate(rkwork(rkworkdim),rkiwork(rkiworkdim))
     iflag=1           
     call rkf45(spf_linear_derivs, totspfdim*zzz, outspfs, time1, time2, relerr, abserr, iflag, rkwork, rkiwork)
     if (iflag/=2) then
        OFLWR "RK45F ERR   Iflag = ", iflag, abserr, relerr; CFLST
     endif
     numiters=rkiwork(1)
     deallocate(rkwork,rkiwork)
  else if (intopt==1) then 
     gbsstepsize=1.d-6
     rkworkdim=20*psilength*zzz; rkiworkdim=20*psilength*zzz
     allocate(rkwork(rkworkdim),rkiwork(rkiworkdim))
     call odex(totspfdim*zzz,gbs_linear_derivs,time1,outspfs,time2,gbsstepsize,relerr,abserr,0,dummysub,0,rkwork,rkworkdim,rkiwork,rkiworkdim,nullreal,nullint,idid)
     numiters=rkiwork(17)
     if (.not.(idid.eq.1)) then
        OFLWR "ERR ODEX", idid; CFLST
     endif
     deallocate(rkwork,rkiwork)
  else if (intopt.eq.3) then
     call expoprop(time1,time2,outspfs,numiters)
  else if (intopt.eq.4) then
     call verlet(outspfs,time1,time2,numiters)
  else
     OFLWR "Intopt not recognized: ", intopt; CFLST
  endif  !! intopt

  call apply_spf_constraints(outspfs)

  call spf_orthogit(outspfs, nulldouble)

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
  implicit none

  integer ::  itime,jtime,times(0:20)=0,numiters=0,linearflag,imc,printflag=1,getlen,qq
  real*8 :: tout, tin, time1, time2
  integer, save :: xxcount=0 
  DATATYPE :: myvalues(mcscfnum)

  firsttime=tin;   lasttime=tout

  numiters=0

  xxcount=xxcount+1

  call system_clock(itime)

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
  if (drivingflag.ne.0) then
     yyy%drivingavectorsxx(:,:,:,1)=yyy%drivingavectorsxx(:,:,:,0)
     yyy%drivingavectorsyy(:,:,:,1)=yyy%drivingavectorsyy(:,:,:,0)
     yyy%drivingavectorszz(:,:,:,1)=yyy%drivingavectorszz(:,:,:,0)
     yyy%drivingorbsxx(:,:,1)=yyy%drivingorbsxx(:,:,0)
     yyy%drivingorbsyy(:,:,1)=yyy%drivingorbsyy(:,:,0)
     yyy%drivingorbszz(:,:,1)=yyy%drivingorbszz(:,:,0)
  endif

  if (numfrozen.ne.0) then
     yyy%frozenexchange(:,:,1) = yyy%frozenexchange(:,:,0)
  endif

  call assign_cptr(yyy%cptr(1),yyy%cptr(0),DATAONE)

  if (sparseopt.ne.0) then
     call assign_sptr(yyy%sptr(1),yyy%sptr(0),DATAONE)
  endif
  if (df_restrictflag.ne.0.and.sparsedfflag.ne.0) then
     if (sparseopt.ne.0) then
        call assign_sptr(yyy%sdfptr(1),yyy%sdfptr(0),DATAONE)
     endif
  endif
  yyy%cmfpsivec(:,1)= yyy%cmfpsivec(:,0)

  call system_clock(jtime);  times(6)=times(6)+jtime-itime

  if (improvedrelaxflag.ne.0) then

     call system_clock(itime)
     if (improvedquadflag.gt.1.and.tin.gt.quadstarttime.and.spf_flag.ne.0) then
        call quadspfs(yyy%cmfpsivec(spfstart,0), qq)
        numiters=numiters+qq
     else
        time1=tin;        time2=tout
        call propspfs(yyy%cmfpsivec(spfstart,1), yyy%cmfpsivec(spfstart,0), time1,time2,0,spf_flag,qq)
        numiters=numiters+qq
     endif
     call system_clock(jtime);     times(4)=times(4)+jtime-itime;     

     if (improvednatflag.ne.0) then
        call system_clock(itime)
!! IT WORKS
        call replace_withnat(printflag)
!! 06-2015        call system_clock(jtime);        times(2)=times(2)+jtime-itime;
        call system_clock(jtime);        times(7)=times(7)+jtime-itime;
     endif

     call system_clock(itime)
     call all_matel()
     call system_clock(jtime);     times(1)=times(1)+jtime-itime;    

     call system_clock(itime)
     if (avector_flag.ne.0) then
        if (improvedquadflag.eq.1.or.improvedquadflag.eq.3) then
           call quadavector(yyy%cmfpsivec(astart(1),0),qq)
        else
           call myconfigeig(www,yyy%cptr(0),yyy%cmfpsivec(astart(1),0),myvalues,mcscfnum,eigprintflag,1,0d0,max(0,improvedrelaxflag-1))
        endif
     endif
     call system_clock(jtime);     times(5)=times(5)+jtime-itime;     call system_clock(itime)

     call get_allden()
     call system_clock(jtime);        times(2)=times(2)+jtime-itime;             call system_clock(itime)
     
     call get_reducedpot()
     if (numfrozen.gt.0) then
        call get_frexchange()
     endif
     call system_clock(jtime);        times(3)=times(3)+jtime-itime;           
     
     if (constraintflag.ne.0) then    !! probably need denmat right
        call system_clock(itime)
        call get_constraint(tout);      
        call system_clock(jtime);        times(7)=times(7)+jtime-itime
     endif

  else  !! IMPROVEDRELAX

        
!! ******* CMF PREDICTOR ***** !!

     linearflag=0
     
     if (threshflag.eq.0.and.(constraintflag.ne.0.or.drivingflag.ne.0)) then   !!.and.tempnoconguess.eq.0) then
        call system_clock(itime)
     
        time1=tin;        time2=tout
        call propspfs(yyy%cmfpsivec(spfstart,1),yyy%cmfpsivec(spfstart,0), time1,time2,linearflag,0,qq)
        numiters=numiters+qq
        call system_clock(jtime);     times(4)=times(4)+jtime-itime;     call system_clock(itime)

        do imc=1,mcscfnum
           call cmf_prop_avector(yyy%cmfpsivec(astart(imc),1), yyy%cmfpsivec(astart(imc),0), linearflag,tin,tout,imc)
        enddo
        
        call system_clock(jtime);        times(5)=times(5)+jtime-itime

        call system_clock(itime)
        call all_matel()
        call system_clock(jtime);     times(1)=times(1)+jtime-itime;   
        
        call system_clock(itime)
        call get_allden()
        call system_clock(jtime);     times(2)=times(2)+jtime-itime
        
        if (constraintflag.ne.0) then
           call system_clock(itime)
           call get_constraint(tout);          
           call system_clock(jtime);        times(7)=times(7)+jtime-itime
        endif
        if (drivingflag.ne.0) then
           call system_clock(itime)
           call drivingtrans(tout);         
           call system_clock(jtime);        times(8)=times(8)+jtime-itime
        endif
     
        call system_clock(itime)
        call get_reducedpot()
        if (numfrozen.gt.0) then
           call get_frexchange()
        endif
        call system_clock(jtime);     times(3)=times(3)+jtime-itime;     call system_clock(itime)  

        linearflag=1

     endif
     

     call system_clock(itime)     
     
     time1=tin;        time2=tout
     call propspfs(yyy%cmfpsivec(spfstart,1),yyy%cmfpsivec(spfstart,0), time1,time2,linearflag,spf_flag,qq)
     numiters=numiters+qq
     call system_clock(jtime);     times(4)=times(4)+jtime-itime;     

     call system_clock(itime)
     call all_matel()
     call system_clock(jtime);     times(1)=times(1)+jtime-itime

     if (threshflag.ne.0) then

        call system_clock(itime)
        do imc=1,mcscfnum
           call cmf_prop_avector(yyy%cmfpsivec(astart(imc),1), yyy%cmfpsivec(astart(imc),0), linearflag,tin,tout,imc)
        enddo
        call system_clock(jtime);     times(5)=times(5)+jtime-itime;     call system_clock(itime)
        call get_allden()
        call system_clock(jtime);     times(2)=times(2)+jtime-itime
        
     endif

     if (constraintflag.ne.0) then
        call system_clock(itime)
        call get_constraint(tout);         
        call system_clock(jtime);        times(7)=times(7)+jtime-itime
     endif
     
     if (drivingflag.ne.0) then
        call system_clock(itime)
        call drivingtrans(tout);         
        call system_clock(jtime);        times(8)=times(8)+jtime-itime
     endif

     
!! ******* LMF CORRECTOR ***** !!

     if (threshflag.ne.0) then

        call get_reducedpot()
        if (numfrozen.gt.0) then
           call get_frexchange()
        endif
        call system_clock(jtime);     times(3)=times(3)+jtime-itime;     call system_clock(itime)  
        
     else

        linearflag=1
        
        call system_clock(itime)
        do imc=1,mcscfnum
           call cmf_prop_avector(yyy%cmfpsivec(astart(imc),1), yyy%cmfpsivec(astart(imc),0), linearflag,tin,tout,imc)
        enddo
        call system_clock(jtime);     times(5)=times(5)+jtime-itime;     call system_clock(itime)
        
        call get_allden()
        call system_clock(jtime);     times(2)=times(2)+jtime-itime
        
        if (constraintflag.ne.0) then
           call system_clock(itime)
           call get_constraint(tout);          
           call system_clock(jtime);        times(7)=times(7)+jtime-itime
        endif
        if (drivingflag.ne.0) then
           call system_clock(itime)
           call drivingtrans(tout);         
           call system_clock(jtime);        times(8)=times(8)+jtime-itime
        endif
        
        call system_clock(itime)

     
        call get_reducedpot()
        if (numfrozen.gt.0) then
           call get_frexchange()
        endif
        call system_clock(jtime);     times(3)=times(3)+jtime-itime;     call system_clock(itime)  
        

        time1=tin;        time2=tout
        call propspfs(yyy%cmfpsivec(spfstart,1),yyy%cmfpsivec(spfstart,0), time1,time2,linearflag,spf_flag,qq)
        numiters=numiters+qq
        
        call system_clock(jtime);     times(4)=times(4)+jtime-itime;     call system_clock(itime)
        
        call all_matel()
        
        call system_clock(jtime);     times(1)=times(1)+jtime-itime;   
        
        if (constraintflag.ne.0) then
           call system_clock(itime)
           call get_constraint(tout);           
           call system_clock(jtime);     times(7)=times(7)+jtime-itime
        endif
        if (drivingflag.ne.0) then
           call system_clock(itime)
           call drivingtrans(tout);         
           call system_clock(jtime);        times(8)=times(8)+jtime-itime
        endif
        
        call system_clock(itime)
        
        call get_reducedpot()
        if (numfrozen.gt.0) then
           call get_frexchange()
        endif
        call system_clock(jtime);     times(3)=times(3)+jtime-itime
     
     endif  !!threshflag
  endif  !!improvedrelax


!!!!! ******** !!!!!!!
!!!!! ******** !!!!!!!
!!!!! ******** !!!!!!!


  if (myrank.eq.1.and.notiming.le.1) then
     call output_denmat(1,tout)
  endif

  if ((myrank.eq.1).and.(notiming.le.1)) then
     if (xxcount==1) then
        open(853, file=timingdir(1:getlen(timingdir)-1)//"/cmf_prop.time.dat", status="unknown")
        write(853,'(A15,100A11)')  "Time",  "matel", "denmat", "reduced", "spfprop", "aprop",  "advance", "constrain", "driving", "#DERIVS"
        close(853)
     endif
        open(853, file=timingdir(1:getlen(timingdir)-1)//"/cmf_prop.time.dat", status="unknown", position="append")
        write(853,'(A3,F12.3,T16, 100I11)')  "T=", tout, times(1:8)/1000, numiters;        close(853)

  endif

end subroutine cmf_prop_wfn


!!!   A-VECTOR PROPAGATION ROUTINE FOR CMF

!! now just cmf and lmf.

subroutine cmf_prop_avector(avectorin,avectorout,linearflag,time1,time2,imc)
  use parameters
  implicit none
  DATATYPE :: avectorin(tot_adim), avectorout(tot_adim),tempvector(tot_adim),dot,csum
  integer :: k, linearflag,imc
  real*8 :: time1,time2,timea,timeb

  if (avector_flag.eq.0) then
     avectorout(:)=avectorin(:)
     return
  endif

  tempvector(:)=avectorin(:)
  do k=1,littlesteps
     timea=time1+(time2-time1)*(k-1)/littlesteps;     timeb=time1+(time2-time1)*k/littlesteps
     call cmf_prop_avector0(tempvector,avectorout,linearflag,timea,timeb,imc)
     tempvector(:)=avectorout(:)
  enddo

  if (threshflag.ne.0) then
     csum=dot(avectorout,avectorout,tot_adim)
     if (par_consplit.ne.0) then
        call mympireduceone(csum)
     endif
     avectorout(:)=avectorout(:)/sqrt(csum)
  endif

end subroutine cmf_prop_avector

subroutine cmf_prop_avector0(avectorin,avectorout,linearflag,time1,time2,imc)
  use linearmod
  use parameters
  use mpimod
  use xxxmod
  use configmod
  use configpropmod
  implicit none

  DATATYPE :: avectorin(tot_adim), avectorout(tot_adim),sum1,sum0,pots(3)=0d0
  integer :: linearflag,imc,itime,jtime,getlen,ii,iflag
  real*8 :: time1,time2,thisstep,midtime,rsum
  integer, save :: times(2)=0, icalled=0

  call system_clock(itime)

  if (sparseopt.eq.0) then
     call zero_cptr(workconfigpointer)
  else
     call zero_sptr(worksparsepointer)
     if (df_restrictflag.ne.0.and.sparsedfflag.ne.0) then
        call zero_sptr(workdfsparsepointer)
     endif
  endif

  thisstep=time2-time1; 

!! 062714
  midtime=(time2+time1)/2d0

  iflag=0
  if (drivingflag.ne.0) then
     call vectdpot(midtime,velflag,pots)
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
        call assign_cptr(workconfigpointer,yyy%cptr(1),thisstep*DATAONE)
     else
        call assign_sptr(worksparsepointer,yyy%sptr(1),thisstep*DATAONE)
        if (df_restrictflag.ne.0.and.sparsedfflag.ne.0) then
           call assign_sptr(workdfsparsepointer,yyy%sdfptr(1),thisstep*DATAONE)
        endif
     endif

     if (drivingflag.ne.0) then
        workdrivingavec(:,:)=0d0

        if (iflag.eq.1) then
           workdrivingavec(:,:)=( yyy%drivingavectorsxx(:,:,imc,1)*pots(1) + &
                yyy%drivingavectorsyy(:,:,imc,1)*pots(2) + yyy%drivingavectorszz(:,:,imc,1)*pots(3) )*thisstep

        endif
     endif
  else

!! MAJOR BUG 062714     midtime=(time2+time1)/2d0

     sum0=(midtime-firsttime)/(lasttime-firsttime) 
     sum1=(lasttime-midtime)/(lasttime-firsttime) 

     if (drivingflag.ne.0) then
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
        call add_sptr(yyy%sptr(1),yyy%sptr(0),worksparsepointer,sum1*thisstep,sum0*thisstep)
        if (df_restrictflag.ne.0.and.sparsedfflag.ne.0) then
           call add_sptr(yyy%sdfptr(1),yyy%sdfptr(0),workdfsparsepointer,sum1*thisstep,sum0*thisstep)
        endif
     endif
  endif

  call system_clock(jtime); times(1)=times(1)+jtime-itime; itime=jtime

  call myconfigprop(www,avectorin,avectorout,midtime)

  call system_clock(jtime); times(2)=times(2)+jtime-itime;

  if (notiming.eq.0.and.myrank.eq.1) then

     if (icalled.eq.0) then
        open(11766, file=timingdir(1:getlen(timingdir)-1)//"/aprop.time.dat", status="unknown")
        write(11766,'(100A11)')  "init", "prop"
        close(11766)
     endif
     open(11766, file=timingdir(1:getlen(timingdir)-1)//"/aprop.time.dat", status="unknown", position="append")
     write(11766,'(100I11)')  times(1:2)
     close(11766)
  endif

  icalled=1

end subroutine cmf_prop_avector0


