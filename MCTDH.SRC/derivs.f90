
!! OPERATION OF MEAN FIELD ON ORBITALS : "DERIVS" ARE DERIVATIVES
!!  OF ORBITALS WRT TIME.  SPF_LINEAR_DERIVS IS MAIN GUY
!!  (DOES LMF -- QMF ETC PROBABLY BROKEN)
!! 
!! IF YOU WANT TO OPERATE WITH REDUCED HAMILTONIAN YOURSELF,
!!  USE ACTREDUCED.

!! ALL MODULES

#include "Definitions.INC"

module derivtimingmod
  implicit none
  integer :: times(20)=0,numcalledhere=0
end module derivtimingmod


module derivativemod
contains

!! WITH TIMEFAC

subroutine getconmat(thistime,ireduced,conmat)
  use parameters
  use xxxmod
  use pulsesubmod
  implicit none
  integer, intent(in) :: ireduced
  real*8, intent(in) ::  thistime 
  DATATYPE, intent(out) :: conmat(nspf,nspf)
  DATATYPE :: pots(3)

  if (constraintflag.eq.0) then
     conmat(:,:)=0d0
     return
  endif

  conmat(:,:) =   yyy%cptr(ireduced)%xconmatel(:,:) * timefac
  if (tdflag.ne.0) then
     call vectdpot(thistime,velflag,pots,-1)
     conmat(:,:) =   conmat(:,:) + &
          yyy%cptr(ireduced)%xconmatelxx(:,:) *pots(1) * timefac + &
          yyy%cptr(ireduced)%xconmatelyy(:,:) *pots(2) * timefac + &
          yyy%cptr(ireduced)%xconmatelzz(:,:) *pots(3) * timefac
  endif

end subroutine getconmat

!! for derivative of PROJECTOR using derivative of spfs.     
!!      on call inspfs is for example jacvectout

  subroutine derproject00(lowspf,highspf,inspfs, outspfs, prospfs, prospfderivs)
    use parameters
    use jactimingmod
    use orbgathersubmod
    use mpisubmod
    implicit none
    integer,intent(in) :: lowspf,highspf
    DATATYPE, intent(in) :: inspfs(spfsize, lowspf:highspf), &
         prospfs(spfsize, nspf),  prospfderivs(spfsize, nspf)
    DATATYPE, intent(out) :: outspfs(spfsize, lowspf:highspf)
    DATATYPE :: csum
    DATATYPE :: mydot(nspf,lowspf:highspf), prodot(nspf,nspf), &
         derdot(nspf,lowspf:highspf) !! AUTOMATIC
    integer :: i,j,numspf,itime,jtime

    numspf=highspf-lowspf+1

    if (numspf.le.0) then
       call waitawhile()
       print *, "OOOGA EdddERRRRR"
       call waitawhile()
       stop
    endif

    call myclock(itime)

    mydot(:,:)=0d0; derdot(:,:)=0d0; prodot(:,:)=0d0
    outspfs(:,:)=0.d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i=lowspf,highspf
       do j=1,nspf
          mydot(j,i) = dot(prospfs(:,j),inspfs(:,i),spfsize)
          derdot(j,i) = dot(prospfderivs(:,j),inspfs(:,i),spfsize)
       enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

    call myclock(jtime); times(3)=times(3)+jtime-itime;     itime=jtime

    if (parorbsplit.eq.3) then
       call mympireduce(mydot,nspf**2)
       call mympireduce(derdot,nspf**2)
    endif

    call myclock(jtime); times(4)=times(4)+jtime-itime;     itime=jtime

    call MYGEMM('N','N',spfsize,numspf,nspf,DATAONE,prospfs,spfsize,&
         derdot(:,lowspf:highspf),nspf,DATAONE,outspfs,spfsize)
    call MYGEMM('N','N',spfsize,numspf,nspf,DATAONE,prospfderivs,spfsize,&
         mydot(:,lowspf:highspf), nspf,DATAONE,outspfs,spfsize)

    call myclock(jtime); times(3)=times(3)+jtime-itime;

    if (jacprojorth.ne.0) then

       call myclock(itime)

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

! need all nspf^2 even if parorbsplit.eq.1

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,csum)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
       do i=lowspf,highspf
          do j=1,nspf
             csum = dot(prospfs(:,j),prospfderivs(:,i),spfsize) + &
                  dot(prospfderivs(:,j),prospfs(:,i),spfsize)
             prodot(j,i) = csum
          enddo
       enddo
!$OMP END DO
!$OMP END PARALLEL

       call myclock(jtime); times(3)=times(3)+jtime-itime;     itime=jtime

       if (parorbsplit.eq.1) then
          call mpiorbgather_nz(prodot,nspf)
       endif
       if (parorbsplit.eq.3) then
          call mympireduce(prodot,nspf**2)
       endif

       call myclock(jtime); times(4)=times(4)+jtime-itime;     itime=jtime

       call MYGEMM('N', 'N', nspf, numspf, nspf, DATAONE, prodot, nspf, &
            mydot(:,lowspf:highspf), nspf, DATAZERO, derdot(:,lowspf:highspf), nspf)
          
       call MYGEMM('N', 'N', spfsize, numspf, nspf, DATANEGONE, prospfs, spfsize, &
            derdot(:,lowspf:highspf), nspf, DATAONE, outspfs, spfsize)

       call myclock(jtime); times(3)=times(3)+jtime-itime

    endif

  end subroutine derproject00

!! subroutine derproject not used currently only derproject00

  subroutine derproject(inspfs, outspfs, prospfs, prospfderivs)
    use parameters
    use orbgathersubmod
    implicit none
    DATATYPE, intent(in) :: inspfs(spfsize, nspf), prospfs(spfsize, nspf),&
         prospfderivs(spfsize, nspf)
    DATATYPE, intent(out) :: outspfs(spfsize, nspf)
    integer :: lowspf,highspf

    lowspf=1; highspf=nspf
    if (parorbsplit.eq.1) then
       call getOrbSetRange(lowspf,highspf)
    endif
    if (highspf.ge.lowspf) then
       call derproject00(lowspf,highspf,inspfs(:,lowspf:highspf),&
            outspfs(:,lowspf:highspf),prospfs,prospfderivs)
    endif
    if (parorbsplit.eq.1) then
       call mpiorbgather(outspfs,spfsize)
    endif

  end subroutine derproject

  subroutine der_gmat00(lowspf,highspf,inspfs, outspfs, &
       ireduced,thistime,prospfs, prospfderivs)
    use parameters
    use mpisubmod
    implicit none
    integer, intent(in) :: ireduced,lowspf,highspf
    real*8, intent(in) :: thistime
    DATATYPE, intent(in) :: inspfs(spfsize, nspf), prospfs(spfsize, nspf),  &
         prospfderivs(spfsize, nspf)
    DATATYPE, intent(out) :: outspfs(spfsize, lowspf:highspf)
    integer :: i,j,numspf
    DATATYPE :: mydot(nspf,lowspf:highspf), &
         derdot(nspf,lowspf:highspf), &            !!  AUTOMATIC
         mydot0(nspf,lowspf:highspf), &
         derdot0(nspf,lowspf:highspf), conmat(nspf,nspf)

    numspf=highspf-lowspf+1

    if (numspf.le.0) then
       call waitawhile()
       print *, "OOOGA EERRRRR"
       call waitawhile()
       stop
    endif

    outspfs(:,:)=0.d0

    if (constraintflag.eq.0) then
       return
    endif

    mydot0(:,:)=0d0; derdot0(:,:)=0d0; mydot(:,:)=0d0; derdot(:,:)=0d0; conmat(:,:)=0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i=lowspf,highspf
       do j=1,nspf
          mydot0(j,i) = dot(prospfs(:,j),inspfs(:,i),spfsize)
          derdot0(j,i) = dot(prospfderivs(:,j),inspfs(:,i),spfsize)
       enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

    if (parorbsplit.eq.3) then
       call mympireduce(mydot0,nspf**2)
       call mympireduce(derdot0,nspf**2)
    endif

    call getconmat(thistime,ireduced,conmat)
    
    call MYGEMM('N','N',nspf,numspf,nspf,DATAONE,conmat,nspf,&
         mydot0(:,lowspf:highspf),nspf,DATAZERO,mydot(:,lowspf:highspf),nspf)
    call MYGEMM('N','N',nspf,numspf,nspf,DATAONE,conmat,nspf,&
         derdot0(:,lowspf:highspf),nspf,DATAZERO,derdot(:,lowspf:highspf),nspf)

    call MYGEMM('N','N',spfsize,numspf,nspf,DATAONE,prospfs,     spfsize,&
         derdot(:,lowspf:highspf),nspf,DATAZERO,outspfs,spfsize)
    call MYGEMM('N','N',spfsize,numspf,nspf,DATAONE,prospfderivs,spfsize,&
         mydot(:,lowspf:highspf), nspf,DATAONE,outspfs,spfsize)

  end subroutine der_gmat00

  subroutine der_gmat(inspfs, outspfs, ireduced,thistime,prospfs, prospfderivs)
    use parameters
    use orbgathersubmod
    implicit none
    integer, intent(in) :: ireduced
    real*8, intent(in) :: thistime
    DATATYPE, intent(in) :: inspfs(spfsize, nspf), prospfs(spfsize, nspf),  &
         prospfderivs(spfsize, nspf)
    DATATYPE, intent(out) :: outspfs(spfsize, nspf)
    integer :: lowspf,highspf,numspf
    
    lowspf=1; highspf=nspf
    if (parorbsplit.eq.1) then
       call getOrbSetRange(lowspf,highspf)
    endif

    numspf=highspf-lowspf+1

    if (numspf.gt.0) then
       call der_gmat00(lowspf,highspf,inspfs,outspfs(:,lowspf:highspf),&
            ireduced,thistime,prospfs,prospfderivs)
    endif

    if (parorbsplit.eq.1) then
       call mpiorbgather(outspfs,spfsize)
    endif

  end subroutine der_gmat

!! WITH TIMEFAC

  subroutine op_gmat_firstorder00(lowspf,highspf,inspfs, outspfs, ireduced,thistime)
    use parameters
    implicit none
    integer,intent(in) :: lowspf,highspf,ireduced
    real*8, intent(in) ::  thistime 
    DATATYPE, intent(in) :: inspfs(spfsize,nspf)
    DATATYPE, intent(out) :: outspfs(spfsize,lowspf:highspf)
    DATATYPE :: conmat(nspf,nspf)        !! AUTOMATIC
    integer :: numspf

    numspf=highspf-lowspf+1

    if (numspf.le.0) then
       call waitawhile()
       print *, "OOOGA EERRRRR"
       call waitawhile()
       stop
    endif
    if (constraintflag.eq.0) then
       outspfs(:,:)=0d0
       return
    endif

!! with timefac
    conmat=0d0
    call getconmat(thistime,ireduced,conmat)

    call MYGEMM('N','N',spfsize,numspf,nspf,DATAONE,inspfs,&
         spfsize,conmat(:,lowspf:highspf),&
         nspf,DATAZERO,outspfs(:,lowspf:highspf),spfsize)

  end subroutine op_gmat_firstorder00

  subroutine op_gmat_thirdorder00(lowspf,highspf,inspfs, outspfs, &
       ireduced,thistime,projspfs)
    use parameters
    use mpisubmod
    implicit none
    integer,intent(in) :: lowspf,highspf,ireduced
    real*8, intent(in) ::  thistime 
    DATATYPE, intent(in) :: inspfs(spfsize,nspf), projspfs(spfsize,nspf)
    DATATYPE, intent(out) :: outspfs(spfsize,lowspf:highspf)
    DATATYPE :: conmat(nspf,nspf), mydot(nspf,nspf), mymat(nspf,nspf)    !! AUTOMATIC
    integer :: i,j,numspf

    numspf=highspf-lowspf+1

    if (numspf.le.0) then
       call waitawhile()
       print *, "OOOGA EERRRRR"
       call waitawhile()
       stop
    endif
    if (constraintflag.eq.0) then
       outspfs(:,:)=0d0
       return
    endif

    mydot(:,:)=0d0;  conmat(:,:)=0d0;  mymat(:,:)=0d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i=1,lowspf,highspf
       do j=1,nspf
          mydot(j,i)=dot(projspfs(:,j),inspfs(:,i),spfsize)
       enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

    if (parorbsplit.eq.3) then
       call mympireduce(mydot,nspf**2)
    endif

!! with timefac
    call getconmat(thistime,ireduced,conmat)

    call MYGEMM('N','N',nspf,numspf,nspf,DATAONE,conmat,nspf,&
         mydot(:,lowspf:highspf),nspf,DATAZERO,mymat(:,lowspf:highspf),nspf)

    call MYGEMM('N','N',spfsize,numspf,nspf,DATAONE,projspfs,spfsize,&
         mymat(:,lowspf:highspf),nspf,DATAZERO,outspfs(:,lowspf:highspf),&
         spfsize)

  end subroutine op_gmat_thirdorder00

  subroutine op_gmat00(lowspf,highspf,inspfs, outspfs, ireduced,thistime,projspfs)
    use parameters
    implicit none
    integer, intent(in) :: lowspf,highspf,ireduced
    real*8, intent(in) ::  thistime 
    DATATYPE, intent(in) :: inspfs(spfsize,nspf), projspfs(spfsize,nspf)
    DATATYPE, intent(out) :: outspfs(spfsize,lowspf:highspf)

    if (jacgmatthird.eq.0) then
       call op_gmat_firstorder00(lowspf,highspf,inspfs, outspfs, ireduced, thistime)
    else
       call op_gmat_thirdorder00(lowspf,highspf,inspfs, outspfs, &
            ireduced, thistime,projspfs)
    endif

  end subroutine op_gmat00

!! WITH TIMEFAC

  subroutine op_gmat(inspfs, outspfs, ireduced,thistime,projspfs)
    use parameters
    use orbgathersubmod
    implicit none
    integer, intent(in) :: ireduced
    real*8, intent(in) ::  thistime 
    DATATYPE, intent(in) :: inspfs(spfsize,nspf), projspfs(spfsize,nspf)
    DATATYPE, intent(out) :: outspfs(spfsize,nspf)
    integer :: lowspf,highspf,numspf

    lowspf=1; highspf=nspf
    if (parorbsplit.eq.1) then
       call getOrbSetRange(lowspf,highspf)
    endif
    
    numspf=highspf-lowspf+1
    if (numspf.gt.0) then
       call op_gmat00(lowspf,highspf,inspfs, outspfs(:,lowspf:highspf), &
            ireduced,thistime,projspfs)
    endif

    if (parorbsplit.eq.1) then
       call mpiorbgather(outspfs,spfsize)
    endif

  end subroutine op_gmat

!! conpropspfs FOR prepropflag.ne.0 !!

  subroutine conpropspfs00(lowspf,highspf,inspfs, outspfs, time1,time2)
    use parameters
    use expsubmod
    implicit none
    integer,intent(in) :: lowspf,highspf
    real*8, intent(in) ::  time1,time2
    DATATYPE, intent(in) :: inspfs(spfsize,nspf)
    DATATYPE, intent(out) :: outspfs(spfsize,lowspf:highspf)
    DATATYPE :: conmat(nspf,nspf)        !! AUTOMATIC
    real*8 :: midtime
    integer :: numspf

    numspf=highspf-lowspf+1

    if (numspf.le.0) then
       call waitawhile()
       print *, "OOOGA EErrrrr 77 42 42"
       call waitawhile()
       stop
    endif
    if (constraintflag.eq.0) then
       OFLWR "what? why call conpropspfs?"; CFLST
       outspfs(:,:)=inspfs(:,:)
       return
    endif

!! with timefac
    midtime=(time1+time2)/2d0

    conmat=0d0
    call getconmat(midtime,0,conmat)

    conmat(:,:)=conmat(:,:) * (time2-time1)

    call expmat(conmat,nspf)

    call MYGEMM('N','N',spfsize,numspf,nspf,DATAONE,inspfs,&
         spfsize,conmat(:,lowspf:highspf),&
         nspf,DATAZERO,outspfs(:,lowspf:highspf),spfsize)

  end subroutine conpropspfs00

  subroutine conpropspfs(inspfs, outspfs, time1,time2)
    use parameters
    use orbgathersubmod
    use spfsubmod
    implicit none
    real*8, intent(in) ::  time1,time2
    DATATYPE, intent(in) :: inspfs(spfsize,nspf)
    DATATYPE, intent(out) :: outspfs(spfsize,nspf)
    integer :: lowspf,highspf,numspf
    real*8 :: orthogerror

    lowspf=1; highspf=nspf
    if (parorbsplit.eq.1) then
       call getOrbSetRange(lowspf,highspf)
    endif
    
    numspf=highspf-lowspf+1
    if (numspf.gt.0) then
       call conpropspfs00(lowspf,highspf,inspfs, outspfs(:,lowspf:highspf),time1,time2)
    endif

    if (parorbsplit.eq.1) then
       call mpiorbgather(outspfs,spfsize)
    endif

    call spf_orthogit(outspfs,orthogerror)
    if (orthogerror.gt.1d-10) then
       OFLWR "WARNING - orthog error in conpropspfs",orthogerror; CFL
    endif

  end subroutine conpropspfs

!! Wmult subroutines for multiplication by reduced operator W, not used at the moment.
!! When scalar terms in W are implemented (see LBNL-AMO-MCTDHF version 1 to do list
!! on the github wiki) then this can be used for various purposes like enforcing
!! Fock eigenfunctions e.g. Hartree-Fock orbitals, like we currently do with
!! enforcing natural orbitals with improvenatflag=1.

  subroutine wmult00(lowspf,highspf,inspfs, outspfs, ireduced)
    use parameters
    use xxxmod
    use orbmultsubmod
    implicit none
    integer,intent(in) :: lowspf,highspf,ireduced
    DATATYPE, intent(in) :: inspfs(spfsize,nspf)
    DATATYPE, intent(out) :: outspfs(spfsize,lowspf:highspf)
    DATATYPE :: workmult(spfsize,lowspf:highspf), &
         spfinvr( spfsize,lowspf:highspf ),&  !! AUTOMATIC
         spfinvrsq(  spfsize,lowspf:highspf), &
         spfproderiv(  spfsize,lowspf:highspf )
    integer :: numspf

    numspf=highspf-lowspf+1

    if (numspf.le.0) then
       call waitawhile()
       print *, "OOOGA EERReeeeRRR"
       call waitawhile()
       stop
    endif

    workmult=0; spfinvr=0; spfinvrsq=0; spfproderiv=0
    outspfs(:,:)=0d0


!! sum over fast index reduced matrices, because doing spfinvrsq= reducedinvrsq * inspfs 
!! BUT 1) store in transposed order and 2) have to reverse the call in BLAS

    call MYGEMM('N', 'N', spfsize,numspf,nspf,DATAONE, inspfs, spfsize, &
         yyy%reducedinvrsq(:,lowspf:highspf,ireduced),&
         nspf,DATAZERO, spfinvrsq(:,lowspf:highspf), spfsize)
    call MYGEMM('N', 'N', spfsize,numspf,nspf,DATAONE, inspfs, spfsize, &
         yyy%reducedinvr(:,lowspf:highspf,ireduced),&
         nspf, DATAZERO, spfinvr(:,lowspf:highspf), spfsize)
    if ((nonuc_checkflag.eq.0)) then
       call MYGEMM('N', 'N', spfsize,numspf,nspf,DATAONE, inspfs, spfsize, &
            yyy%reducedproderiv(:,lowspf:highspf,ireduced),&
            nspf, DATAZERO, spfproderiv(:,lowspf:highspf), spfsize)
    endif
  
    call mult_ke(spfinvrsq(:,lowspf:highspf),outspfs(:,lowspf:highspf),&
         numspf,timingdir,notiming)

    call mult_pot(numspf,spfinvr(:,lowspf:highspf),workmult(:,lowspf:highspf))
    outspfs(:,lowspf:highspf)=outspfs(:,lowspf:highspf)+workmult(:,lowspf:highspf)

    if ((nonuc_checkflag.eq.0)) then
       call op_yderiv(numspf,spfproderiv(:,lowspf:highspf),workmult(:,lowspf:highspf))
       outspfs(:,lowspf:highspf)=outspfs(:,lowspf:highspf) + workmult(:,lowspf:highspf)
    endif

    call mult_reducedpot(lowspf,highspf,inspfs,&
         workmult(:,lowspf:highspf),yyy%reducedpot(:,:,lowspf:highspf,ireduced))
    outspfs(:,lowspf:highspf)=outspfs(:,lowspf:highspf)+workmult(:,lowspf:highspf)
    
  end subroutine wmult00

  subroutine wmult(inspfs, outspfs, ireduced)
    use parameters
    use orbgathersubmod
    implicit none
    DATATYPE, intent(in) :: inspfs(spfsize,nspf)
    DATATYPE, intent(out) :: outspfs(spfsize,nspf)
    integer,intent(in) :: ireduced
    integer :: lowspf,highspf,numspf

    lowspf=1; highspf=nspf
    if (parorbsplit.eq.1) then
       call getOrbSetRange(lowspf,highspf)
    endif

    numspf=highspf-lowspf+1
    if (numspf.gt.0) then
       call wmult00(lowspf,highspf,inspfs, outspfs(:,lowspf:highspf), ireduced)
    endif
    if (parorbsplit.eq.1) then
       call mpiorbgather(outspfs,spfsize)
    endif

  end subroutine wmult

  subroutine denmult00(lowspf,highspf,inspfs, outspfs, ireduced)
    use parameters
    use xxxmod
    implicit none
    integer,intent(in) :: lowspf,highspf,ireduced
    DATATYPE,intent(in) :: inspfs(spfsize,nspf)
    DATATYPE,intent(out) :: outspfs(spfsize,lowspf:highspf)
    integer :: numspf

    numspf=highspf-lowspf+1

    if (numspf.le.0) then
       call waitawhile()
       print *, "DENDEN EERReeeeRRR",lowspf,highspf
       call waitawhile()
       stop
    endif

    call MYGEMM('N','N', spfsize,numspf,nspf,DATAONE, inspfs,spfsize, &
         yyy%denmat(:,lowspf:highspf,ireduced), nspf, DATAZERO, outspfs, &
         spfsize) 
    outspfs=outspfs/numelec

  end subroutine denmult00

  subroutine denmult(inspfs, outspfs, ireduced)
    use parameters
    use orbgathersubmod
    implicit none
    DATATYPE, intent(in) :: inspfs(spfsize,nspf)
    DATATYPE, intent(out) :: outspfs(spfsize,nspf)
    integer,intent(in) :: ireduced
    integer :: lowspf,highspf,numspf

    lowspf=1; highspf=nspf
    if (parorbsplit.eq.1) then
       call getOrbSetRange(lowspf,highspf)
    endif

    numspf=highspf-lowspf+1
    if (numspf.gt.0) then
       call denmult00(lowspf,highspf,inspfs, outspfs(:,lowspf:highspf), ireduced)
    endif

    if (parorbsplit.eq.1) then
       call mpiorbgather(outspfs,spfsize)
    endif

  end subroutine denmult

  subroutine actreduced00(lowspf,highspf,dentimeflag,thistime,inspfs, projspfs, &
       outspfs, ireduced, projflag,conflag)
    use parameters
    use mpimod
    use xxxmod
    use derivtimingmod
    use orbprojectmod
    use orbgathersubmod
    use orbmultsubmod
    use pulsesubmod
    use opmod   !! frozenreduced, hatomreduced
    implicit none
    integer, intent(in) :: lowspf,highspf,dentimeflag,ireduced,projflag,conflag
    real*8, intent(in) :: thistime
    DATATYPE, intent(in) :: inspfs(spfsize, nspf), projspfs(spfsize,nspf)
    DATATYPE,intent(out) :: outspfs(spfsize,lowspf:highspf)
    integer :: itime, jtime, getlen,numspf,myiostat,ispf,jspf
    DATATYPE :: myxtdpot=0,  myytdpot=0, myztdpot=0, pots(3)=0d0
    DATATYPE :: spfmult(spfsize,nspf),workmult(spfsize,lowspf:highspf), &  !! AUTOMATIC
         spfinvr( spfsize,lowspf:highspf), spfr( spfsize,lowspf:highspf ),  &
         spfinvrsq(  spfsize,lowspf:highspf),spfproderiv(  spfsize,lowspf:highspf )
  
    numspf=highspf-lowspf+1

    if (numspf.le.0) then
       call waitawhile()
       print *, "ACT00 EERReeeeRRR",lowspf,highspf
       call waitawhile()
       stop
    endif

    if (tdflag.ne.0) then
       call vectdpot(thistime,velflag,pots,-1)
       myxtdpot=pots(1);  myytdpot=pots(2);  myztdpot=pots(3);
    endif

    numcalledhere=numcalledhere+1

    spfmult=0; workmult=0; spfinvr=0; spfr=0; spfinvrsq=0; spfproderiv=0

!! sum over fast index reduced matrices, because doing spfinvrsq= reducedinvrsq * inspfs 
!!   BUT 1) store in transposed order and 2) have to reverse the call in BLAS

    call myclock(itime)
    if (numr.eq.1) then
       call MYGEMM('N', 'N', spfsize,numspf,nspf,DATAONE, inspfs, spfsize, &
            yyy%denmat(:,lowspf:highspf,ireduced),nspf, DATAZERO, &
            spfinvrsq(:,lowspf:highspf), spfsize)

       spfinvr(:,lowspf:highspf)=spfinvrsq(:,lowspf:highspf)/bondpoints(1)

       if (tdflag.ne.0) then
          spfr(:,lowspf:highspf)=spfinvrsq(:,lowspf:highspf)*bondpoints(1)
       endif

       if (nonuc_checkflag.eq.0) then
          spfproderiv(:,lowspf:highspf)=spfinvrsq(:,lowspf:highspf)*bondpoints(1)
       endif

       spfinvrsq(:,lowspf:highspf)=spfinvrsq(:,lowspf:highspf)/bondpoints(1)**2

    else
       call MYGEMM('N', 'N', spfsize,numspf,nspf,DATAONE, inspfs, spfsize, &
            yyy%reducedinvrsq(:,lowspf:highspf,ireduced),&
            nspf, DATAZERO, spfinvrsq(:,lowspf:highspf), spfsize)
       call MYGEMM('N', 'N', spfsize,numspf,nspf,DATAONE, inspfs, spfsize, &
            yyy%reducedinvr(:,lowspf:highspf,ireduced),&
            nspf, DATAZERO, spfinvr(:,lowspf:highspf), spfsize)
        
       if (tdflag.ne.0) then
          call MYGEMM('N', 'N', spfsize,numspf,nspf,DATAONE, inspfs, spfsize, &
               yyy%reducedr(:,lowspf:highspf,ireduced),&
               nspf, DATAZERO, spfr(:,lowspf:highspf), spfsize)
       endif
     
       if (nonuc_checkflag.eq.0) then
          call MYGEMM('N', 'N', spfsize,numspf,nspf,DATAONE, inspfs, spfsize, &
               yyy%reducedproderiv(:,lowspf:highspf,ireduced),nspf, DATAZERO, &
               spfproderiv(:,lowspf:highspf), spfsize)
       endif
    endif
    call myclock(jtime);  times(1)=times(1)+jtime-itime; itime=jtime

    call mult_ke(spfinvrsq(:,lowspf:highspf),workmult(:,lowspf:highspf),&
         numspf,timingdir,notiming)
    spfmult(:,lowspf:highspf) = spfmult(:,lowspf:highspf) + workmult(:,lowspf:highspf)
    call myclock(jtime);  times(2)=times(2)+jtime-itime;      itime=jtime

!!    OFLWR "CHECKMULT1  ",spfmult(1,1); CFL

    call mult_pot(numspf,spfinvr(:,lowspf:highspf),workmult(:,lowspf:highspf))
    spfmult(:,lowspf:highspf) = spfmult(:,lowspf:highspf) + workmult(:,lowspf:highspf)
    call hatom_op(numspf,spfinvr(:,lowspf:highspf),workmult(:,lowspf:highspf),hatomreduced(:))
    spfmult(:,lowspf:highspf)=spfmult(:,lowspf:highspf)+workmult(:,lowspf:highspf)

!!    OFLWR "CHECKMULT2  ",spfmult(1,1); CFL

    if (numfrozen.gt.0) then
!! DIRECT ONLY in linear operator actreduced.  Exchange treated like driving term.
       call op_frozenreduced(numspf,spfinvr(:,lowspf:highspf),workmult(:,lowspf:highspf),frozenreduced(:))
       spfmult(:,lowspf:highspf)=spfmult(:,lowspf:highspf)+workmult(:,lowspf:highspf)

!!$ exact exchange (too slow), use frozenexchinvr instead

       if (exact_exchange.ne.0) then
          call op_frozen_exchange(lowspf,highspf,spfinvr(:,lowspf:highspf),workmult(:,lowspf:highspf))
          spfmult(:,lowspf:highspf)=spfmult(:,lowspf:highspf)+workmult(:,lowspf:highspf)
       endif

    endif
    call myclock(jtime);     times(3)=times(3)+jtime-itime;      itime=jtime

    if (tdflag.ne.0) then
       select case (velflag)
       case (0)
          call lenmultiply(numspf,spfr(:,lowspf:highspf),workmult(:,lowspf:highspf), &
               myxtdpot,myytdpot,myztdpot)
          spfmult(:,lowspf:highspf)=spfmult(:,lowspf:highspf)+workmult(:,lowspf:highspf)
       case default
          call velmultiply(numspf,spfinvr(:,lowspf:highspf),workmult(:,lowspf:highspf), &
               myxtdpot,myytdpot,myztdpot)
          spfmult(:,lowspf:highspf)=spfmult(:,lowspf:highspf)+workmult(:,lowspf:highspf)
       end select
    endif  !! tdpot
    call myclock(jtime);        times(4)=times(4)+jtime-itime;        itime=jtime
  
    if (nonuc_checkflag.eq.0) then
       call op_yderiv(numspf,spfproderiv(:,lowspf:highspf),workmult(:,lowspf:highspf))
       spfmult(:,lowspf:highspf)=spfmult(:,lowspf:highspf)+workmult(:,lowspf:highspf)
    endif
    call myclock(jtime);     times(5)=times(5)+jtime-itime;         itime=jtime

!!    OFLWR "CHECKMULT2BBB  ",spfmult(1,1); CFL

    call mult_reducedpot(lowspf,highspf,inspfs,workmult(:,lowspf:highspf),&
         yyy%reducedpot(:,:,lowspf:highspf,ireduced))
    spfmult(:,lowspf:highspf)=spfmult(:,lowspf:highspf) + workmult(:,lowspf:highspf)
    call myclock(jtime);  times(6)=times(6)+jtime-itime;  

    if (scalarflag.ne.0) then
       if (.not.use_fockmatrix) then
          OFLWR "programmer error use_fockmatrix"; CFLST
       endif
       do jspf=lowspf,highspf
          do ispf=1,nspf
             spfmult(:,jspf)=spfmult(:,jspf) + yyy%fockmatrix(jspf,ispf,ireduced) * inspfs(:,ispf)
          enddo
       enddo
    endif

!!    OFLWR "CHECKMULT3  ",spfmult(1,1); CFL

!! WITH TIMEFAC
    if (dentimeflag.ne.0) then
       call myclock(itime)
       if (parorbsplit.eq.1) then
          call mpiorbgather_nz(spfmult,spfsize)
       endif
       call myclock(jtime);        times(10)=times(10)+jtime-itime;      itime=jtime
       call MYGEMM('N','N', spfsize,numspf,nspf,timefac, spfmult,spfsize, &
            yyy%invdenmat(:,lowspf:highspf,ireduced), nspf, DATAZERO, workmult, spfsize)
       spfmult(:,lowspf:highspf)=workmult(:,lowspf:highspf)

       call myclock(jtime);  times(7)=times(7)+jtime-itime
!! no more factor -1
!!    else
!!       spfmult(:,lowspf:highspf) = spfmult(:,lowspf:highspf) * (-1)
    endif

!!    OFLWR "CHECKMULT4  ",spfmult(1,1); CFL

    if (projflag==1) then
       call myclock(itime)
       if (parorbsplit.eq.1) then
          call mpiorbgather_nz(spfmult,spfsize)
       endif
       call myclock(jtime);  times(10)=times(10)+jtime-itime;    itime=jtime
       call project00(lowspf,highspf,spfmult(:,lowspf:highspf), workmult, projspfs)
       outspfs(:,:) = spfmult(:,lowspf:highspf) - workmult(:,lowspf:highspf)
       call myclock(jtime);     times(8)=times(8)+jtime-itime;      
    else
       outspfs(:,:)=spfmult(:,lowspf:highspf)
    endif

!!    OFLWR "CHECKMULT5  ",spfmult(1,1); CFLST
     
    if (constraintflag/=0.and.conflag.ne.0) then
       call myclock(itime)
       call op_gmat00(lowspf,highspf,inspfs,workmult,ireduced,thistime,projspfs)
       outspfs(:,:)=outspfs(:,:)+workmult(:,lowspf:highspf)
       call myclock(jtime);        times(9)=times(9)+jtime-itime
    endif

    if ((myrank.eq.1).and.(notiming.eq.0)) then
       if (numcalledhere==1) then
          open(853, file=timingdir(1:getlen(timingdir))//"/actreduced.time.dat", &
               status="unknown",iostat=myiostat)
          call checkiostat(myiostat,"opening actreduced timing file")
          write(853,'(T16,100A9)',iostat=myiostat) &
               " rmult ",&     !   (1)
               " ke ",&        !   (2)
               " pot ",&       !   (3)
               " pulse",&       !   (4)
               " nuc ",&       !   (5)
               " twoe ",&      !   (6)
               " invdenmat ",& !   (7)
               " project ", &  !   (8)
               " constrain ", &!   (9)
               " MPI "         !  (10)
          close(853)
          call checkiostat(myiostat,"writing actreduced timing file")
       endif

       if (mod(numcalledhere,timingout).eq.0) then
          open(853, file=timingdir(1:getlen(timingdir))//"/actreduced.time.dat", &
               status="unknown", position="append",iostat=myiostat)
          call checkiostat(myiostat,"opening actreduced timing file")
          write(853,'(A3,F12.4,15I9)',iostat=myiostat) "T= ", thistime,  times(1:10)/1000
          call checkiostat(myiostat,"writing actreduced timing file")
          close(853)
       endif
    endif

  end subroutine actreduced00

!! MAIN ROUTINE TO OPERATE WITH REDUCED HAMILTONIAN (TIMES TIMEFAC!)
!!
!!  outspfs = (1-P) H inspfs
!! where P is projector onto opspfs
!! H is inverse denmat times reducedham

  subroutine actreduced0(dentimeflag,thistime,inspfs, projspfs, &
       outspfs, ireduced, projflag,conflag)
    use parameters
    use derivtimingmod
    use orbgathersubmod
    implicit none
    integer, intent(in) :: dentimeflag,ireduced,projflag,conflag
    real*8, intent(in) :: thistime
    DATATYPE, intent(in) :: inspfs(spfsize, nspf), projspfs(spfsize,nspf)
    DATATYPE,intent(out) :: outspfs(spfsize,nspf)
    integer :: lowspf,highspf,itime,jtime

    lowspf=1; highspf=nspf
    if (parorbsplit.eq.1) then
       call getOrbSetRange(lowspf,highspf)
    endif
    if (highspf.ge.lowspf) then
       call actreduced00(lowspf,highspf,dentimeflag,thistime,inspfs, projspfs, &
            outspfs(:,lowspf:highspf), ireduced, projflag,conflag)
    endif
    if (parorbsplit.eq.1) then
       call myclock(itime)
       call mpiorbgather(outspfs,spfsize)
       call myclock(jtime);        times(10)=times(10)+jtime-itime
    endif

  end subroutine actreduced0

!! only used for relax ; time not needed    
!!     this is actually the conjugate-transpose of the operator, even if cmctdh
!!     (not used)

  subroutine actreducedconjg0(thistime,inspfs, projspfs, outspfs, &
       ireduced, projflag,conflag)
    use parameters
    implicit none
    DATATYPE, intent(in) :: inspfs(spfsize,nspf),projspfs(spfsize,nspf)
    DATATYPE, intent(out) ::  outspfs(spfsize,nspf)
    integer,intent(in) :: ireduced,projflag,conflag
    real*8,intent(in) :: thistime
    DATATYPE :: ttempspfs(spfsize,nspf) !! AUTOMATIC

    ttempspfs(:,:)=ALLCON(inspfs(:,:))
    call actreduced0(1,thistime,ttempspfs, projspfs, outspfs, &
         ireduced, projflag,conflag)
    outspfs=ALLCON(outspfs)

  end subroutine actreducedconjg0

end module derivativemod



module orbdermod
  use derivativemod

contains

!! FOR ODEX PROPAGATION (NOT DEFAULT)

  subroutine gbs_derivs(notusedint,thistime,psi,psip,notuseddbl,notusedint2)
    use parameters
    implicit none
    integer,intent(in) :: notusedint,notusedint2
    real*8,intent(in) :: thistime,notuseddbl
    DATATYPE,intent(in) :: psi(tot_adim*mcscfnum+totspfdim)
    DATATYPE,intent(out) :: psip(tot_adim*mcscfnum+totspfdim)
    call all_derivs(thistime,psi,psip)
  end subroutine gbs_derivs


!! VMF SUBROUTINE (possibly broken, not used for a while)
  subroutine all_derivs(thistime,in_xpsi, out_xpsip)
    use parameters
    use mpimod
    use configmod
    use xxxmod
    use sparsemultmod
    use basissubmod
    use getstuffmod
    implicit none
    DATATYPE,intent(in) :: in_xpsi(tot_adim*mcscfnum+totspfdim)
    DATATYPE,intent(out) :: out_xpsip(tot_adim*mcscfnum+totspfdim)
    DATATYPE :: xspfs(totspfdim),xavec(tot_adim,mcscfnum),&
         xspfsp(totspfdim),xavecp(tot_adim,mcscfnum)
    DATATYPE :: avector(tot_adim)   !! AUTOMATIC
    real*8 :: thistime
    integer :: itime,jtime,getlen,myiostat,spfstart,spfend
    integer, save :: times(20)=0, numcalledhere=0,imc

    numcalledhere=numcalledhere+1

    xspfsp(:)=0d0
    if (tot_adim.gt.0) then
       xavec(:,:)=RESHAPE(in_xpsi(1:tot_adim*mcscfnum),(/tot_adim,mcscfnum/))
       yyy%cmfavec(:,:,0)=xavec(:,:)
       avector=0
       xavecp(:,:)=0d0
    endif

    spfstart=tot_adim*mcscfnum+1;     spfend=tot_adim*mcscfnum+totspfdim

    xspfs(:)=in_xpsi(spfstart:spfend)
    yyy%cmfspfs(:,0)=xspfs(:)

    call get_stuff0(thistime,times)

!! ireduced should be zero right   07-2015

    if (spf_flag.ne.0) then
       call myclock(itime)
       call actreduced0(1,thistime,xspfs,xspfs,xspfsp,0,1,1)
       call myclock(jtime);     times(5)=times(5)+jtime-itime
    endif

  !! AVECTOR PART.

    if (avector_flag.ne.0) then
       call myclock(itime)
       do imc=1,mcscfnum
          if (tot_adim.gt.0) then
             avector(:)=xavec(:,imc)
          endif
          call basis_project(www,numr,avector)
          call sparseconfigmult(www,avector,xavecp(:,imc),&
               yyy%cptr(0),yyysptr(0),1,1,1,1,thistime,imc)
          call basis_project(www,numr,xavecp(:,imc))
       enddo
    endif

    if (tot_adim.gt.0) then
       xavecp(:,:)=xavecp(:,:)*timefac
    endif

    call myclock(jtime);  times(6)=times(6)+jtime-itime
  
    if ((myrank.eq.1).and.(notiming.eq.0)) then
       if (numcalledhere==1) then
          open(853, file=timingdir(1:getlen(timingdir))//"/all_deriv.time.dat", &
               status="unknown",iostat=myiostat)
          call checkiostat(myiostat,"opening all_deriv timing file")
          write(853,'( T16, 100A15)') " init ", " matel ", " denmat ", " getredham ", &
               " actreduced ", " amult "
          close(853)
       endif
       if (mod(numcalledhere,timingout).eq.0) then
          open(853, file=timingdir(1:getlen(timingdir))//"/all_deriv.time.dat", &
               status="unknown", position="append",iostat=myiostat)
          call checkiostat(myiostat,"writing all_deriv timing file")
          write(853,'(A3,F12.3,100I15)') "T= ", thistime,  times(1:6)/1000
          close(853)
       endif
    endif


    out_xpsip(spfstart:spfend)=xspfsp(:)
    out_xpsip(1:tot_adim*mcscfnum)=RESHAPE(xavecp(:,:),(/tot_adim*mcscfnum/))

  end subroutine all_derivs


!! MAIN DERIVATIVE ROUTINE FOR ODEX

  subroutine gbs_linear_derivs(notusedint,thistime,spfsin,spfsout,&
       notuseddbl,notusedint2)
    use parameters
    implicit none
    DATATYPE,intent(in) :: spfsin(spfsize, nspf)
    DATATYPE,intent(out) :: spfsout(spfsize, nspf)
    integer,intent(in) :: notusedint,notusedint2
    real*8,intent(in) :: thistime,notuseddbl
    call spf_linear_derivs(thistime,spfsin,spfsout)
  end subroutine gbs_linear_derivs


!! MAIN DERIVATIVE ROUTINE FOR RK, EXPO

  subroutine spf_linear_derivs(thistime,spfsin,spfsout)
    use parameters
    use linearmod   !! effective_cmf_linearflag
    implicit none
    DATATYPE,intent(in) :: spfsin(spfsize, nspf)
    DATATYPE,intent(out) :: spfsout(spfsize, nspf)
    real*8,intent(in) :: thistime
    call spf_linear_derivs0(effective_cmf_linearflag,1,thistime,spfsin,spfsout,1,1)
  end subroutine spf_linear_derivs


  subroutine spf_linear_derivs0(inlinearflag,dentimeflag,thistime,spfsin,&
       spfsout, projflag,conflag)
    use derivtimingmod
    use parameters
    use linearmod    !! firsttime,lasttime
    use xxxmod  !! frozenexchange and driving orbs
    use orbprojectmod
    use orbgathersubmod
    use pulsesubmod
    implicit none
    integer,intent(in) :: inlinearflag,dentimeflag,projflag,conflag
    real*8,intent(in) :: thistime
    DATATYPE,intent(in) :: spfsin(spfsize, nspf)
    DATATYPE,intent(out) :: spfsout(spfsize, nspf)
    DATATYPE :: facs(0:1),csum,pots(3)
    DATATYPE,allocatable :: tempspfs(:,:),workspfs(:,:)
    real*8 :: rsum
    integer ::  jjj, itop,lowspf,highspf,numspf,itime,jtime

    if (inlinearflag.eq.1) then
       itop=1;     
       facs(0)=(thistime-firsttime)/(lasttime-firsttime); 
       facs(1)=1d0-facs(0)
    else
       itop=0;     facs(0)=1d0;     facs(1)=0d0
    endif

    lowspf=1; highspf=nspf
    if (parorbsplit.eq.1) then
       call getOrbSetRange(lowspf,highspf)
    endif
    numspf=highspf-lowspf+1

    spfsout(:,:)=0d0; 

    if (numspf.gt.0) then

       allocate(tempspfs(spfsize,nspf), workspfs(spfsize,nspf))
       tempspfs(:,:)=0d0; workspfs(:,:)=0d0

       do jjj=0,itop
          call actreduced00(lowspf,highspf,dentimeflag,thistime,spfsin,spfsin,&
               workspfs(:,lowspf:highspf),jjj, projflag,conflag)
          spfsout(:,lowspf:highspf)=spfsout(:,lowspf:highspf) + &
               workspfs(:,lowspf:highspf)*facs(jjj)
       enddo

!! EXCHANGE IS TREATED AS A DRIVING TERM.  DOES NOT GO IN LINEAR OPERATOR ACTREDUCED.
!!  (it is prohibitive to call op_frozen_exchange repeatedly)
!! EXCHANGE AND DRIVING CONTRIBUTE TO JACOBIAN (jacoperate) via projector.

       if (numfrozen.gt.0.and.exact_exchange.eq.0) then

          do jjj=0,itop

             call myclock(itime)

             if (projflag.ne.0) then
                select case (exchange_mode)
                case(0)
                   call project00(1,nspf,yyy%frozenexchinvr(:,1:nspf,jjj),&
                        workspfs(:,1:nspf),spfsin)
                case(1)
                   call MYGEMM('N','N',spfsize,nspf,nspf,DATAONE,&
                        spfsin(:,:), spfsize, yyy%frozenexchmat(:,1:nspf,jjj), nspf,&
                        DATAZERO, workspfs(:,1:nspf),spfsize)
                case default
                   OFLWR "error exchange_mode derivs"; CFLST
                end select
                tempspfs(:,1:nspf) = yyy%frozenexchinvr(:,1:nspf,jjj) - &
                     workspfs(:,1:nspf)
             else
                tempspfs(:,1:nspf) = yyy%frozenexchinvr(:,1:nspf,jjj)
             endif
             call myclock(jtime);        times(8)=times(8)+jtime-itime

             if (dentimeflag.ne.0) then
!! TIMEFAC and facs HERE
                csum=timefac*facs(jjj)
                call MYGEMM('N','N', spfsize,numspf,nspf,csum, &
                     tempspfs(:,1:nspf),spfsize, &
                     yyy%invdenmat(:,lowspf:highspf,jjj), nspf, DATAZERO, &
                     workspfs(:,lowspf:highspf), spfsize)
             else
                workspfs(:,lowspf:highspf)= &   !! no more factor -1
                     tempspfs(:,lowspf:highspf)*facs(jjj)
             endif

             spfsout(:,lowspf:highspf)=spfsout(:,lowspf:highspf) + workspfs(:,lowspf:highspf)

             call myclock(jtime);    times(7)=times(7)+jtime-itime;
          enddo
       endif

!! DRIVING (PSI-PRIME)

       if (drivingflag.ne.0) then
          if (dentimeflag.eq.0) then
             OFLWR "Error, no drivingflag for quad"; CFLST !! invdenmat already in drivingorbs
          endif
          call vectdpot(thistime,velflag,pots,-1)
          rsum=0
          do jjj=1,3
             rsum=rsum+abs(pots(jjj))**2
          enddo
          if (rsum.ne.0d0) then
             tempspfs(:,:)=0d0
             do jjj=0,itop
                tempspfs(:,lowspf:highspf)=tempspfs(:,lowspf:highspf)+ ( &
                     yyy%drivingorbsxx(:,lowspf:highspf,jjj) * pots(1) + &
                     yyy%drivingorbsyy(:,lowspf:highspf,jjj) * pots(2) + &
                     yyy%drivingorbszz(:,lowspf:highspf,jjj) * pots(3) ) &
                     *facs(jjj) * timefac                             !! WITH TIMEFAC
             enddo
             call myclock(itime)
             if (projflag.ne.0) then
                call project00(lowspf,highspf,tempspfs(:,lowspf:highspf),&
                     workspfs(:,lowspf:highspf),spfsin)
                spfsout(:,lowspf:highspf)=spfsout(:,lowspf:highspf)+&
                     tempspfs(:,lowspf:highspf)-workspfs(:,lowspf:highspf)
             else
                spfsout(:,lowspf:highspf)=spfsout(:,lowspf:highspf) + &
                     tempspfs(:,lowspf:highspf)
             endif
             call myclock(jtime);        times(8)=times(8)+jtime-itime
          endif
       endif

       deallocate(tempspfs,workspfs)

    endif

    if (parorbsplit.eq.1) then
       call myclock(itime)
       call mpiorbgather(spfsout,spfsize)
       call myclock(jtime);        times(10)=times(10)+jtime-itime
    endif

  end subroutine spf_linear_derivs0

end module orbdermod


