
!! ALL MODULES

!! SUBROUTINES FOR COMPUTING MEAN FIELD !!

!! REDUCEDPOT is usual notation bra, ket   other reduced mats opposite notation for ease of actreducedconjg

#include "Definitions.INC"


module tworeducedxmod
contains

subroutine get_tworeducedx(www,reducedpottally,avector1,in_avector2,numvects)
  use r_parameters
  use walkmod
  use dotmod
  use densubmod
  use mpisubmod
  implicit none
  integer,intent(in) :: numvects
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: avector1(numr,www%firstconfig:www%lastconfig,numvects),  &
       in_avector2(numr,www%firstconfig:www%lastconfig,numvects)
  DATATYPE,intent(out) :: reducedpottally(www%nspf,www%nspf,www%nspf,www%nspf)
  DATATYPE,allocatable :: avector2(:,:,:), mytally(:,:,:,:), &
       holeden(:,:), tempvec(:,:,:)
  DATATYPE ::  a1(numr,numvects), a2(numr,numvects), csum
  DATAECS :: rvalues(numr)
  integer ::   ispf, jspf, iispf, jjspf ,  config2, config1,dirphase, iwalk,ii,ihop

!! DO SUMMA (parconsplit.ne.0 and sparsesummaflag.eq.2, "circ")

  allocate(avector2(numr,www%numconfig,numvects))

  if (www%lastconfig.ge.www%firstconfig) then
     avector2(:,www%firstconfig:www%lastconfig,:)=in_avector2(:,:,:)
  endif

  if (www%parconsplit.ne.0) then
     do ii=1,numvects
        call mpiallgather(avector2(:,:,ii),www%numconfig*numr,&
             www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
     enddo
  endif

  reducedpottally(:,:,:,:)=0.d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rvalues,config1,a1,a2,iwalk,dirphase,config2,ii,ispf,iispf,jspf,jjspf,ihop,csum,mytally)

!! REDUCTION(+:reducedpottally) caused crashes!  OMP CRITICAL instead below.

  allocate(mytally(www%nspf,www%nspf,www%nspf,www%nspf))
  mytally(:,:,:,:)=0d0

  rvalues(:)=1d0/bondpoints(:)

        !! doubly off diagonal walks

!$OMP DO SCHEDULE(DYNAMIC)
  do config1=www%botconfig,www%topconfig

     a1(:,:)=avector1(:,config1,:)
        
     do ihop=1,www%numdoublehops(config1)
        config2=www%doublehop(ihop,config1)

        do ii=1,numvects
           a2(:,ii)=avector2(:,config2,ii) * rvalues(:)
        enddo

        csum=dot(a1,a2,numvects*numr)           !! 1/R factor above

        do iwalk=www%doublehopwalkstart(ihop,config1),www%doublehopwalkend(ihop,config1)
           
           dirphase=www%doublewalkdirphase(iwalk+www%dcol(config1))

!! switched 2-2016, keep this logic the same
           ispf=www%doublewalkdirspf(1,iwalk+www%dcol(config1))   !KET2 
           jspf=www%doublewalkdirspf(2,iwalk+www%dcol(config1))   !BRA2 (walk)
           iispf=www%doublewalkdirspf(3,iwalk+www%dcol(config1))  !KET1
           jjspf=www%doublewalkdirspf(4,iwalk+www%dcol(config1))  !BRA1 (walk)

           mytally(ispf,jspf,iispf,jjspf) =  &    
                mytally(ispf,jspf,iispf,jjspf) +  &
                dirphase*csum           !! 1/R factor above
           
           mytally(iispf,jjspf,ispf,jspf) =  &
                mytally(iispf,jjspf,ispf,jspf) + &
                dirphase*csum           !! 1/R factor above

        enddo

     enddo
  enddo   ! config1
!$OMP END DO
!$OMP CRITICAL
  reducedpottally(:,:,:,:)=reducedpottally(:,:,:,:)+mytally(:,:,:,:)
!$OMP END CRITICAL
  deallocate(mytally)
!$OMP END PARALLEL

  deallocate(avector2)

  call mympireduce(reducedpottally(:,:,:,:), www%nspf**4)

!!  OFLWR "TALLY ", reducedpottally(1,1,1,1); CFL

  if (www%holeflag.ne.0) then
     allocate(tempvec(numr,www%firstconfig:www%lastconfig,numvects),&
          holeden(www%nspf,www%nspf))
     csum=0
     rvalues(:)=1d0/bondpoints(:)

     if (www%lastconfig.ge.www%firstconfig) then
        do ii=1,numvects
           do config1=www%firstconfig,www%lastconfig
              tempvec(:,config1,ii) = in_avector2(:,config1,ii) * rvalues
           enddo
        enddo
        csum=dot(avector1,tempvec,numr*www%localnconfig*numvects)
     endif
     if (www%parconsplit.ne.0) then
        call mympireduceone(csum)
     endif

!! second avector(avector1 here) is conjugated in getdenmat00
     call getdenmat00(0,www,in_avector2,avector1,rvalues,holeden,numr,numvects)

     do ispf=1,www%nspf
        do jspf=1,www%nspf
           reducedpottally(ispf,ispf,jspf,jspf) = reducedpottally(ispf,ispf,jspf,jspf) + 4d0*csum
           reducedpottally(ispf,jspf,jspf,ispf) = reducedpottally(ispf,jspf,jspf,ispf) - 2d0*csum
        enddo
        reducedpottally(:,:,ispf,ispf) = reducedpottally(:,:,ispf,ispf) - holeden(:,:) * 2 
        reducedpottally(ispf,ispf,:,:) = reducedpottally(ispf,ispf,:,:) - holeden(:,:) * 2 
        reducedpottally(:,ispf,ispf,:) = reducedpottally(:,ispf,ispf,:) + holeden(:,:) 
        reducedpottally(ispf,:,:,ispf) = reducedpottally(ispf,:,:,ispf) + TRANSPOSE(holeden(:,:)) 
     enddo
     deallocate(holeden,tempvec)
  endif


!!  OFLWR "FIRST"
!!#ifdef REALGO
!!  write(mpifileptr,'(F20.10)') reducedpottally(:,1,2,3)
!!#else
!!  write(mpifileptr,'(2F20.10)') reducedpottally(:,1,2,3)
!!#endif
!!  CFLST

end subroutine get_tworeducedx

end module tworeducedxmod


module reducedhamsubmod
contains

subroutine get_reducedham()
  use parameters
  use configmod
  use tworeducedxmod
  use xxxmod
  implicit none
  integer :: itime,jtime,times(100)

  call myclock(itime)
  call get_tworeducedx(www,yyy%reducedpottally(:,:,:,:,0),&
       yyy%cmfavec(:,:,0),yyy%cmfavec(:,:,0),mcscfnum)
  call myclock(jtime);  times(4)=times(4)+jtime-itime

  call myclock(itime)
  call get_reducedproderiv(www,yyy%reducedproderiv(:,:,0),&
       yyy%cmfavec(:,:,0),yyy%cmfavec(:,:,0),mcscfnum)
  call myclock(jtime);  times(6)=times(6)+jtime-itime

  call myclock(itime)

  if (numr.eq.1) then
     yyy%reducedinvr(:,:,0) = yyy%denmat(:,:,0) / bondpoints(1)
     yyy%reducedinvrsq(:,:,0) = yyy%denmat(:,:,0) / bondpoints(1)**2
     yyy%reducedr(:,:,0) = yyy%denmat(:,:,0) * bondpoints(1)
  else
     call get_reducedr(www,yyy%reducedinvr,yyy%reducedinvrsq,yyy%reducedr,&
          yyy%cmfavec(:,:,0),yyy%cmfavec(:,:,0),mcscfnum)
  endif

  call myclock(itime);  times(7)=times(7)+itime-jtime

contains

subroutine get_reducedproderiv(www,reducedproderiv,avector1,in_avector2,numvects)
  use fileptrmod
  use r_parameters
  use opmod   !! rkemod, proderivmod
  use walkmod
  use dotmod
  use mpisubmod
  implicit none
  integer,intent(in) :: numvects
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: avector1(numr,www%firstconfig:www%lastconfig,numvects), &
       in_avector2(numr,www%firstconfig:www%lastconfig,numvects)
  DATATYPE,intent(out) :: reducedproderiv(www%nspf,www%nspf)
  DATATYPE,allocatable :: avector2(:,:,:), tempvec(:,:,:)
  DATATYPE :: a1(numr,numvects), a2(numr,numvects), a2mult(numr,numvects), mypro(numr,numr), &
       myredpro(www%nspf,www%nspf),csum
  integer ::  config1,config2,  ispf,jspf,  dirphase,     iwalk,ii,ihop

!! DO SUMMA (parconsplit.ne.0 and sparsesummaflag.eq.2, "circ")

  allocate(avector2(numr,www%numconfig,numvects))
  avector2(:,:,:)=0d0
  if (www%lastconfig.ge.www%firstconfig) then
     avector2(:,www%firstconfig:www%lastconfig,:)=in_avector2(:,:,:)
  endif

  if (www%parconsplit.ne.0) then
     do ii=1,numvects
        call mpiallgather(avector2(:,:,ii),www%numconfig*numr,&
             www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
     enddo
  endif

  reducedproderiv(:,:)=0.d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(a1,a2,mypro,myredpro,config1,iwalk,config2,dirphase,ispf,jspf,a2mult,ihop,csum) 

!! removing REDUCTION(+:reducedproderiv) due to get_twoereducedx problems OMP CRITICAL instead

  mypro(:,:)=proderivmod(:,:)
  myredpro(:,:)=0d0

!$OMP DO SCHEDULE(DYNAMIC)
  do config1=www%botconfig,www%topconfig

     do ihop=1,www%numsinglehops(config1)
        config2=www%singlehop(ihop,config1)

        a1(:,:)=avector1(:,config1,:)
        a2(:,:)=avector2(:,config2,:)

        call MYGEMM('N','N',numr,numvects,numr,DATAONE,&
             mypro(:,:),numr,a2(:,:),numr,DATAZERO,a2mult(:,:),numr)

        csum=dot(a1,a2mult,numvects*numr)

        do iwalk=www%singlehopwalkstart(ihop,config1),www%singlehopwalkend(ihop,config1)

           dirphase=www%singlewalkdirphase(iwalk+www%scol(config1))
           ispf=www%singlewalkopspf(1,iwalk+www%scol(config1))
           jspf=www%singlewalkopspf(2,iwalk+www%scol(config1))

           myredpro(jspf,ispf)=myredpro(jspf,ispf)+ &
                dirphase*csum

        enddo

! PREVIOUS (jj,ii private)
!        do jj=1,numr
!           do ii=1,numr
!              a1(:)=avector1(ii,config1,:)
!              a2(:)=avector2(jj,config2,:)
!              
!              reducedproderiv(jspf,ispf)=reducedproderiv(jspf,ispf)+ &
!                   dirphase*dot(a1,a2,numvects)*proderivmod(ii,jj)
!           enddo
!        enddo

     enddo
  enddo
!$OMP END DO
!$OMP CRITICAL
  reducedproderiv(:,:)=reducedproderiv(:,:)+myredpro(:,:)
!$OMP END CRITICAL
!$OMP END PARALLEL

  deallocate(avector2)

  call mympireduce(reducedproderiv(:,:), www%nspf**2)

  if (www%holeflag.ne.0) then

     allocate(tempvec(numr,www%firstconfig:www%lastconfig,numvects))

     do config1=www%botconfig,www%topconfig

        a2(:,:)=in_avector2(:,config1,:)

        call MYGEMM('N','N',numr,numvects,numr,DATAONE,&
             proderivmod(:,:),numr,a2(:,:),numr,DATAZERO,a2mult(:,:),numr)
        tempvec(:,config1,:)=a2mult(:,:)
     enddo

     if (www%parconsplit.eq.0) then
        do ii=1,numvects
           call mpiallgather(tempvec(:,:,ii),www%numconfig*numr,&
                www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
        enddo
     endif

     csum=0
     if (www%lastconfig.ge.www%firstconfig) then
        csum=dot(avector1,tempvec,numr*www%localnconfig*numvects)
     endif
     if (www%parconsplit.ne.0) then
        call mympireduceone(csum)
     endif
     reducedproderiv(:,:) = reducedproderiv(:,:) * (-1)
     do ispf=1,www%nspf
        reducedproderiv(ispf,ispf) = reducedproderiv(ispf,ispf) + csum * 2d0
     enddo
     deallocate(tempvec)
  endif

end subroutine get_reducedproderiv


subroutine get_reducedr(www,reducedinvr,reducedinvrsq,reducedr,avector1,in_avector2,numvects)
  use fileptrmod
  use r_parameters
  use walkmod
  use dotmod
  use mpisubmod
  implicit none
  integer,intent(in) :: numvects
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: avector1(numr,www%firstconfig:www%lastconfig,numvects), &
       in_avector2(numr,www%firstconfig:www%lastconfig,numvects)
  DATATYPE,intent(out) :: reducedinvr(www%nspf,www%nspf),reducedr(www%nspf,www%nspf), &
       reducedinvrsq(www%nspf,www%nspf)
  DATATYPE,allocatable :: avector2(:,:,:),tempr(:,:,:),tempinvr(:,:,:),tempinvrsq(:,:,:)
  DATATYPE ::  a1(numr,numvects), a2(numr,numvects), a2r(numr,numvects), &
       a2inv(numr,numvects), a2invsq(numr,numvects), rdot,invdot,invsqdot
  integer ::  config1,config2,   ispf,jspf,  dirphase,    iwalk,ii, ihop
  DATAECS ::  invrvalues(numr),invrsqvalues(numr),rvalues(numr)
  DATATYPE :: myinvr(www%nspf,www%nspf),myr(www%nspf,www%nspf),  myinvrsq(www%nspf,www%nspf),&
       csumr,csuminvr,csuminvrsq

  if (numr.eq.1) then
     OFLWR "programmer fail, don't call get_reducedr if numr.eq.1"; CFLST
  endif

  reducedinvr(:,:)=0.d0;  reducedr(:,:)=0.d0;  reducedinvrsq(:,:)=0.d0

!! DO SUMMA (parconsplit.ne.0 and sparsesummaflag.eq.2, "circ")

  allocate(avector2(numr,www%numconfig,numvects))
  avector2(:,:,:)=0d0

  if (www%lastconfig.ge.www%firstconfig) then
     avector2(:,www%firstconfig:www%lastconfig,:) = in_avector2(:,:,:)
  endif
  if (www%parconsplit.ne.0) then
     do ii=1,numvects
        call mpiallgather(avector2(:,:,ii),www%numconfig*numr,&
             www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
     enddo
  endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rvalues,invrvalues,invrsqvalues,config1,a1,a2,a2r,a2inv,a2invsq,iwalk,config2,dirphase,ispf,jspf,myinvr,myr,myinvrsq,ihop,invdot,rdot,invsqdot) 

!! removing REDUCTION(+:reducedinvr,reducedr,reducedinvrsq) due to get_twoereducedx problems CRITICAL instead

  rvalues(:)=bondpoints(:)
  invrvalues(:)=(1.d0/bondpoints(:))
  invrsqvalues(:)=(1.d0/bondpoints(:)**2)

  myinvr(:,:)=0.d0;  myr(:,:)=0.d0;  myinvrsq(:,:)=0.d0

!$OMP DO SCHEDULE(DYNAMIC)
  do config1=www%botconfig,www%topconfig
     a1(:,:)=avector1(:,config1,:)

     do ihop=1,www%numsinglehops(config1)
        config2=www%singlehop(ihop,config1)

        a2(:,:)=avector2(:,config2,:)
        do ii=1,numvects
           a2r(:,ii)=a2(:,ii)*rvalues(:)
           a2inv(:,ii)=a2(:,ii)*invrvalues(:)
           a2invsq(:,ii)=a2(:,ii)*invrsqvalues(:)
        enddo

        invdot=dot(a1,a2inv,numvects*numr)
        rdot=dot(a1,a2r,numvects*numr)
        invsqdot=dot(a1,a2invsq,numvects*numr)

        do iwalk=www%singlehopwalkstart(ihop,config1),www%singlehopwalkend(ihop,config1)

           dirphase=www%singlewalkdirphase(iwalk+www%scol(config1))
           ispf=www%singlewalkopspf(1,iwalk+www%scol(config1))
           jspf=www%singlewalkopspf(2,iwalk+www%scol(config1))

           myinvr(jspf,ispf)=myinvr(jspf,ispf)+         dirphase*invdot
           myr(jspf,ispf)=myr(jspf,ispf)+         dirphase*rdot
           myinvrsq(jspf,ispf)=myinvrsq(jspf,ispf)+     dirphase*invsqdot

        enddo

     enddo
  enddo
!$OMP END DO
!$OMP CRITICAL
  reducedinvr(:,:)=reducedinvr(:,:)+myinvr(:,:)
  reducedinvrsq(:,:)=reducedinvrsq(:,:)+myinvrsq(:,:)
  reducedr(:,:)=reducedr(:,:)+myr(:,:)
!$OMP END CRITICAL
!$OMP END PARALLEL

  deallocate(avector2)
  
  call mympireduce(reducedr(:,:), www%nspf**2)
  call mympireduce(reducedinvr(:,:), www%nspf**2)
  call mympireduce(reducedinvrsq(:,:), www%nspf**2)

  if (www%holeflag.ne.0) then
     csumr=0; csuminvr=0; csuminvrsq=0

     allocate(tempr(numr,www%firstconfig:www%lastconfig,numvects),&
          tempinvr(numr,www%firstconfig:www%lastconfig,numvects),&
          tempinvrsq(numr,www%firstconfig:www%lastconfig,numvects))

     if (www%lastconfig.ge.www%firstconfig) then

        tempr=0; tempinvr=0; tempinvrsq=0
        do ii=1,numvects
           do config1=www%firstconfig,www%lastconfig
              tempr(:,config1,ii) = in_avector2(:,config1,ii) * bondpoints(:)
              tempinvr(:,config1,ii) = in_avector2(:,config1,ii) / bondpoints(:)
              tempinvrsq(:,config1,ii) = in_avector2(:,config1,ii) / bondpoints(:)**2
           enddo
        enddo
        csumr     =dot(avector1,tempr,numr*www%localnconfig*numvects)
        csuminvr  =dot(avector1,tempinvr,numr*www%localnconfig*numvects)
        csuminvrsq=dot(avector1,tempinvrsq,numr*www%localnconfig*numvects)
     endif
     if (www%parconsplit.ne.0) then
        call mympireduceone(csumr)
        call mympireduceone(csuminvr)
        call mympireduceone(csuminvrsq)
     endif
     reducedr(:,:) = reducedr(:,:) * (-1)
     reducedinvr(:,:) = reducedinvr(:,:) * (-1)
     reducedinvrsq(:,:) = reducedinvrsq(:,:) * (-1)
     do ispf=1,www%nspf
        reducedr(ispf,ispf) = reducedr(ispf,ispf) + csumr * 2d0
        reducedinvr(ispf,ispf) = reducedinvr(ispf,ispf) + csuminvr * 2d0
        reducedinvrsq(ispf,ispf) = reducedinvrsq(ispf,ispf) + csuminvrsq * 2d0
     enddo

     deallocate(tempr,tempinvr,tempinvrsq)

  endif

end subroutine get_reducedr

end subroutine get_reducedham

end module reducedhamsubmod


module meansubmod
contains

subroutine get_allden()
  use denmatxmod
  use reducedhamsubmod
  implicit none

!!  call get_reducedham();  call getdenmatx();

  call getdenmatx();
  call get_reducedham();  

end subroutine get_allden


!! IN THE END reducedpottally is bra2,ket2,bra1,ket1 ; 
!! I sum it up by conjugating ket2 bra2 ket1 bra1.  Notice a1 is conjg not a2.
!! reducedpottally is the TWO ELECTRON REDUCED DENSITY MATRIX.

subroutine get_reducedpot()
  use xxxmod
  use configmod
  use opmod
  use parameters
  use mpi_orbsetmod
  implicit none

  if (parorbsplit.eq.1) then
     call get_reducedpot0(www,yyy%reducedpottally(:,:,:,:,0), yyy%reducedpot(:,:,:,0),&
          twoereduced(:,:,:),firstmpiorb,firstmpiorb+orbsperproc-1)
  else
     call get_reducedpot0(www,yyy%reducedpottally(:,:,:,:,0), yyy%reducedpot(:,:,:,0),&
          twoereduced(:,:,:),1,nspf)
  endif
  call mpibarrier()

contains

subroutine get_reducedpot0(www,intwoden,outpot,twoereduced,firstorb,lastorb)
  use walkmod
  use spfsize_parameters
  use mpi_orbsetmod
  use mpisubmod
  implicit none
  integer,intent(in) :: firstorb,lastorb
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: twoereduced(reducedpotsize, www%nspf,firstorb:lastorb),&
       intwoden(www%nspf,www%nspf,www%nspf,www%nspf)
  DATATYPE,intent(out) :: outpot(reducedpotsize,www%nspf,firstorb:lastorb)
  integer :: lowspf,highspf
#ifdef MPIFLAG
  DATATYPE,allocatable :: workpot(:,:,:),workpot2(:,:,:),temptwoden(:,:,:,:)
  integer :: numspf,norbs,iproc,ispf,jspf
#endif

  lowspf=1; highspf=www%nspf

#ifdef MPIFLAG
  if (parorbsplit.eq.1) then
     call getorbsetrange(lowspf,highspf)
  endif
  numspf=highspf-lowspf+1
  norbs=lastorb-firstorb+1
  if (numspf.le.0) then
     return
  endif
#endif

! I have bra2,ket2,bra1,ket1. (chemists' notation for contraction)

! e.g. M= -2, -1, 2, 1 all different classes. 
!   no zeroed entries even if constraining, in general, I believe.

!  if(deb ugflag.ne.0) then
!     temptwoden(:,:,:,:)=intwoden(:,:,:,:)
!     do ispf=1,www%nspf
!        do jspf=1,www%nspf
!           if (orbclass(ispf).ne.orbclass(jspf)) then
!              temptwoden(:,ispf,:,jspf)=0d0
!              temptwoden(ispf,:,jspf,:)=0d0
!           endif
!        enddo
!     enddo
!     call MYGEMM('N','N', reducedpotsize,www%nspf**2,www%nspf**2,(1.0d0,0.d0), &
!         twoereduced,reducedpotsize,temptwoden,www%nspf**2, (0.d0,0.d0),outpot,&
!         reducedpotsize)
!  else

#ifdef MPIFLAG
  if (parorbsplit.ne.1) then
#endif
     call MYGEMM('N','N', reducedpotsize,www%nspf**2,www%nspf**2,DATAONE, &
          twoereduced,reducedpotsize,intwoden,www%nspf**2, DATAZERO,&
          outpot,reducedpotsize)
#ifdef MPIFLAG
  else
     allocate(temptwoden(www%nspf,firstorb:lastorb,www%nspf,norbs),&
          workpot(reducedpotsize,www%nspf,norbs),workpot2(reducedpotsize,www%nspf,norbs))
     outpot=0
     ispf=0
     do iproc=1,nzprocsperset
        jspf=min(ispf+orbsperproc,www%nspf)
        if (jspf.ge.ispf+1) then

           temptwoden=0d0;  workpot=0
           temptwoden(:,lowspf:highspf,:,1:jspf-ispf)=intwoden(:,lowspf:highspf,:,ispf+1:jspf)

           call MYGEMM('N','N', reducedpotsize,www%nspf*(jspf-ispf),www%nspf*numspf,DATAONE,&
                twoereduced,reducedpotsize,temptwoden,&
                www%nspf*norbs, DATAZERO,workpot,reducedpotsize)
           call mympireduceto_local(workpot,workpot2,reducedpotsize*www%nspf*(jspf-ispf),iproc,&
                NZ_COMM_ORB(myorbset))

           if (orbrank.eq.iproc) then
              outpot(:,:,:)=outpot(:,:,:)+workpot2(:,:,1:numspf)
           endif
        endif
        ispf=ispf+orbsperproc
     enddo
     deallocate(workpot,temptwoden,workpot2)
  endif
#endif

end subroutine get_reducedpot0

end subroutine get_reducedpot

end module meansubmod

