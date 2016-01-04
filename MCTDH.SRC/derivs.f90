
!! OPERATION OF MEAN FIELD ON ORBITALS : "DERIVS" ARE DERIVATIVES
!!  OF ORBITALS WRT TIME.  SPF_LINEAR_DERIVS IS MAIN GUY
!!  (DOES LMF -- QMF ETC PROBABLY BROKEN)
!! 
!! IF YOU WANT TO OPERATE WITH REDUCED HAMILTONIAN YOURSELF,
!!  USE ACTREDUCED.

#include "Definitions.INC"

!! NEVERMIND FOR CONFIGURATIONS -- NEWCOORDS 011514
!!!!   yyy%cmfpsivec  should NOT appear in these routines -- called with psi arguments
!!! nor should iglobalprop (iglobalmod)

!! FOR ODEX PROPAGATION (NOT DEFAULT)

subroutine gbs_derivs(notusedint,thistime,psi,psip,notuseddbl,notusedint2)
  use parameters
  implicit none
  integer,intent(in) :: notusedint,notusedint2
  real*8,intent(in) :: thistime,notuseddbl
  DATATYPE,intent(in) :: psi(psilength)
  DATATYPE,intent(out) :: psip(psilength)
  call all_derivs(thistime,psi,psip)
end subroutine gbs_derivs


subroutine all_derivs(thistime,xpsi, xpsip)
  use parameters
  use mpimod
  use configmod
  use xxxmod
  implicit none
  DATATYPE,intent(in) :: xpsi(psilength)
  DATATYPE,intent(out) :: xpsip(psilength)
  DATATYPE :: avector(tot_adim)   !! AUTOMATIC
  real*8 :: thistime
  integer :: itime,jtime,getlen
  integer, save :: times(20)=0, numcalledhere=0,imc

  numcalledhere=numcalledhere+1

  yyy%cmfpsivec(:,0)=xpsi(:)

  call system_clock(itime)
  call get_allden()
  call system_clock(jtime);  times(3)=times(3)+jtime-itime
  
  call system_clock(itime)
  call all_matel()
  call system_clock(jtime);  times(2)=times(2)+jtime-itime
  
     call system_clock(itime)     
     call get_reducedpot()
     call system_clock(jtime);     times(4)=times(4)+jtime-itime

     call system_clock(itime)
     if (constraintflag.eq.1) then
        OFLWR "CHECK/FIX GETDENCONSTRAINT"; CFLST
!!        call get_den const raint(xpsi(astart(1)))
     endif
     if (constraintflag.eq.2) then
        OFLWR "CHECK/FIX GETDFCONSTRAINT"; CFLST
!!        call get_dfc onstra int(xpsi(astart(1)))
     endif
!! ireduced should be zero right   07-2015
     call actreduced0(thistime,xpsi(spfstart),xpsi(spfstart),xpsip(spfstart),0,1,1)
     call system_clock(jtime);     times(5)=times(5)+jtime-itime

  call system_clock(itime)

  !! AVECTOR PART.

  do imc=1,mcscfnum
     avector(:)=xpsi(astart(imc):aend(imc))

     call basis_project(www,numr,avector)

     call sparseconfigmult(www,avector,xpsip(astart(imc)),yyy%cptr(0),yyy%sptr(0),1,1,1,1,thistime,imc)
     call basis_project(www,numr,xpsip(astart(imc)))

  enddo
  
  xpsip(astart(1):aend(mcscfnum))=xpsip(astart(1):aend(mcscfnum))*timefac

  call system_clock(jtime);  times(6)=times(6)+jtime-itime
  
  if ((myrank.eq.1).and.(notiming.eq.0)) then
     if (numcalledhere==1) then
        open(853, file=timingdir(1:getlen(timingdir)-1)//"/all_deriv.time.dat", status="unknown")

        write(853,'( T16, 100A15)') " init ", " matel ", " denmat ", " getredham ", " actreduced ", " amult "
        close(853)
     endif
     if (mod(numcalledhere,timingout).eq.0) then
        open(853, file=timingdir(1:getlen(timingdir)-1)//"/all_deriv.time.dat", status="unknown", position="append")
        write(853,'(A3,F12.3,100I15)') "T= ", thistime,  times(1:6)/1000;        close(853)
     endif
  endif
end subroutine all_derivs


subroutine gbs_linear_derivs(notusedint,thistime,spfsin,spfsout,notuseddbl,notusedint2)
  use parameters
  implicit none
  DATATYPE,intent(in) :: spfsin(spfsize, nspf)
  DATATYPE,intent(out) :: spfsout(spfsize, nspf)
  integer,intent(in) :: notusedint,notusedint2
  real*8,intent(in) :: thistime,notuseddbl
  call spf_linear_derivs(thistime,spfsin,spfsout)
end subroutine gbs_linear_derivs


subroutine spf_linear_derivs(thistime,spfsin,spfsout)
  use parameters
  implicit none
  DATATYPE,intent(in) :: spfsin(spfsize, nspf)
  DATATYPE,intent(out) :: spfsout(spfsize, nspf)
  real*8,intent(in) :: thistime
  call spf_linear_derivs0(thistime,spfsin,spfsout,1)
end subroutine spf_linear_derivs


!! MAIN DERIVATIVE ROUTINE FOR RK, EXPO

subroutine spf_linear_derivs0(thistime,spfsin,spfsout, allflag)
  use linearmod
  use parameters
  use xxxmod  !! frozenexchange
  implicit none
  DATATYPE,intent(in) :: spfsin(spfsize, nspf)
  DATATYPE,intent(out) :: spfsout(spfsize, nspf)
  integer,intent(in) :: allflag
  real*8,intent(in) :: thistime
  DATATYPE :: facs(0:1)
  integer ::  jjj, ibot
  DATATYPE ::   spfmult(spfsize,nspf), spfmult2(spfsize,nspf)  !!AUTOMATIC

  if (effective_cmf_linearflag.eq.1) then
     ibot=0;     facs(0)=(thistime-firsttime)/(lasttime-firsttime);     facs(1)=1d0-facs(0)
  else
     ibot=1;     facs(0)=0d0;     facs(1)=1d0
  endif
  spfsout(:,:) = 0.d0

  if (effective_cmf_spfflag.eq.0) then
     if (constraintflag==0) then
        return
     endif
     do jjj=ibot,1
        call op_gmat(spfsin(:,:),spfmult(:,:),jjj,thistime,spfsin(:,:))
        spfsout(:,:)=spfsout(:,:)+spfmult(:,:)*facs(jjj)
     enddo
     return
  endif

  do jjj=ibot,1
     call actreduced0(thistime,spfsin,spfsin,spfmult(:,:),jjj, allflag,allflag)
     spfsout(:,:)=spfsout(:,:)+spfmult(:,:)*facs(jjj)
  enddo

!! EXCHANGE IS TREATED AS A DRIVING TERM.  DOES NOT GO IN LINEAR OPERATOR ACTREDUCED.
!!  (it is prohibitive to call op_frozen_exchange repeatedly)
!!  DOES NOT CONTRIBUTE TO JACOBIAN (jacoperate)

  if (numfrozen.gt.0) then

     spfmult2(:,:)=0d0
     do jjj=ibot,1

! previous
!        spfmult(:,:)=0d0
!        do ispf=1,nspf
!           do jspf=1,nspf
!              spfmult(:,ispf) = spfmult(:,ispf) + &
!                   yyy%frozenexchange(:,jspf,jjj) * yyy%reducedinvr(jspf,ispf,jjj) 
!           enddo
!        enddo

!! facs(jjj) here
        call MYGEMM('N','N',spfsize,nspf,nspf,facs(jjj), yyy%frozenexchange(:,:,jjj),spfsize,&
             yyy%reducedinvr(:,:,jjj),nspf, DATAZERO, spfmult,spfsize)
!! TIMEFAC HERE
        call MYGEMM('N','N', spfsize,nspf,nspf,timefac, spfmult(:,:),spfsize, &
             yyy%invdenmat(:,:,jjj), nspf, DATAONE, spfmult2(:,:), spfsize)
     enddo
     call oneminusproject(spfmult2(:,:),spfmult(:,:),spfsin)
     spfsout(:,:)=spfsout(:,:)+spfmult(:,:)
  endif

end subroutine spf_linear_derivs0


subroutine driving_linear_derivs(thistime,spfsin,spfsout)
  use linearmod
  use parameters
  use xxxmod  !! driving orbs
  implicit none
  DATATYPE,intent(in) :: spfsin(spfsize, nspf)
  DATATYPE,intent(out) :: spfsout(spfsize, nspf)
  real*8,intent(in) :: thistime
  real*8 :: facs(0:1),rsum
  integer ::  jjj, ibot
  DATATYPE :: spfmult(spfsize,nspf),pots(3)  !! AUTOMATIC

  spfsout(:,:)=0d0

  if (effective_cmf_spfflag.eq.0) then
     return
  endif

  call vectdpot(thistime,velflag,pots,-1)
  rsum=0
  do jjj=1,3
     rsum=rsum+abs(pots(jjj))**2
  enddo
  if (rsum.eq.0d0) then
     return
  endif

  if (effective_cmf_linearflag.eq.1) then
     ibot=0;     facs(0)=(thistime-firsttime)/(lasttime-firsttime);     facs(1)=1d0-facs(0)
  else
     ibot=1;     facs(0)=0d0;     facs(1)=1d0
  endif

  spfmult(:,:)=0d0
  do jjj=ibot,1
     spfmult(:,:)=spfmult(:,:)+ ( &
          yyy%drivingorbsxx(:,:,jjj) * pots(1) + &
          yyy%drivingorbsyy(:,:,jjj) * pots(2) + &
          yyy%drivingorbszz(:,:,jjj) * pots(3) ) &
          *facs(jjj) * timefac
  enddo

  call oneminusproject(spfmult(:,:),spfsout(:,:),spfsin(:,:))

end subroutine driving_linear_derivs


!! MAIN ROUTINE TO OPERATE WITH REDUCED HAMILTONIAN (TIMES TIMEFAC!)
!!
!!  outspfs = (1-P) H inspfs
!! where P is projector onto opspfs
!! H is inverse denmat times reducedham

subroutine actreduced0(thistime,inspfs0, projspfs, outspfs, ireduced, projflag,conflag)
  use parameters
  use mpimod
  use configptrmod
  use xxxmod
  implicit none
  integer, intent(in) :: ireduced,projflag,conflag
  real*8, intent(in) :: thistime
  DATATYPE, intent(in) :: inspfs0(spfsize, nspf), projspfs(spfsize,nspf)
  DATATYPE,intent(out) :: outspfs(spfsize,nspf)
  integer :: itime, jtime, getlen, lowspf,highspf,numspf
  integer, save :: times(0:20)=0,numcalledhere=0
  DATATYPE :: myxtdpot=0,  myytdpot=0, myztdpot=0, pots(3)=0d0
  DATATYPE ::   inspfs(spfsize,nspf), workmult(spfsize,nspf),spfmult(spfsize,nspf),&  !! AUTOMATIC
       spfinvr( spfsize,nspf ), spfr( spfsize,nspf ),  &
       spfinvrsq(  spfsize,nspf),spfproderiv(  spfsize,nspf )

  if (tdflag.eq.1) then
     call vectdpot(thistime,velflag,pots,-1)
     myxtdpot=pots(1);  myytdpot=pots(2);  myztdpot=pots(3);
  endif

  inspfs(:,:)=inspfs0(:,:)

  numcalledhere=numcalledhere+1
  call system_clock(itime)

 !! sum over fast index reduced matrices, because doing spfinvrsq= reducedinvrsq * inspfs BUT 1) store in transposed order and 2) have to reverse the call in BLAS

  if (numr.eq.1) then
     call MYGEMM('N', 'N', spfsize,nspf,nspf,DATAONE, inspfs, spfsize, yyy%denmat(:,:,ireduced),nspf, DATAZERO, spfinvrsq(:,:), spfsize)

     spfinvr(:,:)=spfinvrsq(:,:)/bondpoints(1)

     if (tdflag.ne.0) then
        spfr(:,:)=spfinvrsq(:,:)*bondpoints(1)
     endif

     spfinvrsq(:,:)=spfinvrsq(:,:)/bondpoints(1)**2
  else
     call MYGEMM('N', 'N', spfsize,nspf,nspf,DATAONE, inspfs, spfsize, yyy%reducedinvrsq(:,:,ireduced),nspf, DATAZERO, spfinvrsq(:,:), spfsize)
     call MYGEMM('N', 'N', spfsize,nspf,nspf,DATAONE, inspfs, spfsize, yyy%reducedinvr(:,:,ireduced),nspf, DATAZERO, spfinvr(:,:), spfsize)

     if (tdflag.eq.1) then
        call MYGEMM('N', 'N', spfsize,nspf,nspf,DATAONE, inspfs, spfsize, yyy%reducedr(:,:,ireduced),nspf, DATAZERO, spfr(:,:), spfsize)
     endif
     
     if (nonuc_checkflag.eq.0) then
        call MYGEMM('N', 'N', spfsize,nspf,nspf,DATAONE, inspfs, spfsize, yyy%reducedproderiv(:,:,ireduced),nspf, DATAZERO, spfproderiv(:,:), spfsize)
     endif
  endif

  call system_clock(jtime);  times(1)=times(1)+jtime-itime

  lowspf=1; highspf=nspf
  if (parorbsplit.eq.1.and.nprocs.gt.1) then
     lowspf=1; highspf=nspf
     call getOrbSetRange(lowspf,highspf)
  endif
  numspf=highspf-lowspf+1

  spfmult=0.d0

     call system_clock(itime)

     call mult_ke(spfinvrsq(:,lowspf:highspf),workmult(:,lowspf:highspf),numspf,timingdir,notiming)
     spfmult(:,lowspf:highspf) = spfmult(:,lowspf:highspf) + workmult(:,lowspf:highspf)
     call system_clock(jtime);  times(2)=times(2)+jtime-itime;  call system_clock(itime)

     call mult_pot(numspf,spfinvr(:,lowspf:highspf),workmult(:,lowspf:highspf))
     spfmult(:,lowspf:highspf) = spfmult(:,lowspf:highspf) + workmult(:,lowspf:highspf)
     call hatom_op(numspf,spfinvr(:,lowspf:highspf),workmult(:,lowspf:highspf))
     spfmult(:,lowspf:highspf)=spfmult(:,lowspf:highspf)+workmult(:,lowspf:highspf)
     if (numfrozen.gt.0) then
!! DIRECT ONLY.
        call op_frozenreduced(numspf,spfinvr(:,lowspf:highspf),workmult(:,lowspf:highspf))
        spfmult(:,lowspf:highspf)=spfmult(:,lowspf:highspf)+workmult(:,lowspf:highspf)
     endif
     call system_clock(jtime);     times(3)=times(3)+jtime-itime;         call system_clock(itime)

     if (tdflag.eq.1) then
        select case (velflag)
        case (0)
           call lenmultiply(numspf,spfr(:,lowspf:highspf),workmult(:,lowspf:highspf), myxtdpot,myytdpot,myztdpot)
           spfmult(:,lowspf:highspf)=spfmult(:,lowspf:highspf)+workmult(:,lowspf:highspf)
        case default
           call velmultiply(numspf,spfinvr(:,lowspf:highspf),workmult(:,lowspf:highspf), myxtdpot,myytdpot,myztdpot)
           spfmult(:,lowspf:highspf)=spfmult(:,lowspf:highspf)+workmult(:,lowspf:highspf)
        end select
     endif  !! tdpot

     call system_clock(jtime);        times(4)=times(4)+jtime-itime;        call system_clock(itime)
  
     if (nonuc_checkflag.eq.0) then
        call noparorbsupport("another call mult_yderiv!!")
        call op_yderiv(numspf,spfproderiv(:,lowspf:highspf),workmult(:,lowspf:highspf))
        spfmult(:,lowspf:highspf)=spfmult(:,lowspf:highspf)+workmult(:,lowspf:highspf)
     endif
     call system_clock(jtime);     times(5)=times(5)+jtime-itime; call system_clock(itime)

!! NOW ONLY OUTPUTS ONE, TAKES ALL
     call mult_reducedpot(lowspf,highspf,inspfs(:,:),workmult(:,lowspf:highspf),yyy%reducedpot(:,:,:,ireduced))
     spfmult(:,lowspf:highspf)=spfmult(:,lowspf:highspf) + workmult(:,lowspf:highspf)
     call system_clock(jtime);  times(6)=times(6)+jtime-itime;  



  if (parorbsplit.eq.1.and.nprocs.gt.1) then
     call system_clock(itime)
     call mpiorbgather(spfmult,spfsize)
     call system_clock(jtime);  times(10)=times(10)+jtime-itime
  endif

!! WITH TIMEFAC

  call system_clock(itime)
  call MYGEMM('N','N', spfsize,nspf,nspf,timefac, spfmult,spfsize, yyy%invdenmat(:,:,ireduced), nspf, DATAZERO, outspfs, spfsize)
  call system_clock(jtime);  times(7)=times(7)+jtime-itime


  if (projflag==1) then
     call system_clock(itime)
     call project(outspfs, spfmult, projspfs)
     outspfs(:,:) = outspfs(:,:) - spfmult(:,:)
     call system_clock(jtime);     times(8)=times(8)+jtime-itime;      
  endif  
     
  if (constraintflag/=0.and.conflag.ne.0) then
     call system_clock(itime)
     call op_gmat(inspfs,spfmult,ireduced,thistime,projspfs)
     outspfs(:,:)=outspfs(:,:)+spfmult(:,:)
     call system_clock(jtime);        times(9)=times(9)+jtime-itime
  endif

  if ((myrank.eq.1).and.(notiming.eq.0)) then
     if (numcalledhere==1) then
        open(853, file=timingdir(1:getlen(timingdir)-1)//"/actreduced.time.dat", status="unknown")
        write(853,'(T16,100A9)') " rmult "," ke "," pot ","pulse", " nuc "," twoe "," invdenmat ",&
             " project ", " constrain ", " MPI ";        close(853)
     endif

     if (mod(numcalledhere,timingout).eq.0) then
        open(853, file=timingdir(1:getlen(timingdir)-1)//"/actreduced.time.dat", status="unknown", position="append")
        write(853,'(A3,F12.4,15I9)') "T= ", thistime,  times(1:10)/1000
        close(853);        close(853)
     endif
  endif

end subroutine actreduced0


!! WITH TIMEFAC

subroutine op_gmat(inspfs, outspfs, ireduced,thistime,projspfs)
  use parameters
  implicit none
  DATATYPE, intent(in) :: inspfs(spfsize,nspf), projspfs(spfsize,nspf)
  DATATYPE, intent(out) :: outspfs(spfsize,nspf)
  integer, intent(in) :: ireduced
  real*8, intent(in) ::  thistime 
  if (jacgmatthird.eq.0) then
     call op_gmat_firstorder(inspfs, outspfs, ireduced, thistime)
  else
     call op_gmat_thirdorder(inspfs, outspfs, ireduced, thistime,projspfs)
  endif
end subroutine op_gmat


!! WITH TIMEFAC

subroutine getconmat(thistime,ireduced,conmat)
  use parameters
  use xxxmod
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
  if (tdflag.eq.1) then
     call vectdpot(thistime,velflag,pots,-1)
     conmat(:,:) =   conmat(:,:) + &
          yyy%cptr(ireduced)%xconmatelxx(:,:) *pots(1) * timefac + &
          yyy%cptr(ireduced)%xconmatelyy(:,:) *pots(2) * timefac + &
          yyy%cptr(ireduced)%xconmatelzz(:,:) *pots(3) * timefac
  endif

end subroutine getconmat


subroutine op_gmat_firstorder(inspfs, outspfs, ireduced,thistime)
  use parameters
  implicit none
  DATATYPE, intent(in) :: inspfs(spfsize,nspf)
  DATATYPE, intent(out) :: outspfs(spfsize,nspf)
  integer, intent(in) :: ireduced
  real*8, intent(in) ::  thistime 
  DATATYPE :: conmat(nspf,nspf)        !! AUTOMATIC

  if (constraintflag.eq.0) then
     outspfs(:,:)=0d0
     return
  endif

!! with timefac
  call getconmat(thistime,ireduced,conmat)

  call MYGEMM('N','N',spfsize,nspf,nspf,DATAONE,inspfs,spfsize,conmat,nspf,DATAZERO,outspfs,spfsize)

!!  do ispf=1,nspf
!!     do jspf=1,nspf
!!        outspfs(:,ispf) = outspfs(:,ispf) + inspfs(:,jspf) * conmat(jspf,ispf)
!!     enddo
!!  enddo

end subroutine op_gmat_firstorder



subroutine op_gmat_thirdorder(inspfs, outspfs, ireduced,thistime,projspfs)
  use parameters
  implicit none
  DATATYPE, intent(in) :: inspfs(spfsize,nspf), projspfs(spfsize,nspf)
  DATATYPE, intent(out) :: outspfs(spfsize,nspf)
  integer, intent(in) :: ireduced
  real*8, intent(in) ::  thistime 
  DATATYPE :: conmat(nspf,nspf), mydot(nspf,nspf), mymat(nspf,nspf)    !! AUTOMATIC
  DATATYPE :: dot
  integer :: i,j

  if (constraintflag.eq.0) then
     outspfs(:,:)=0d0
     return
  endif

  do i=1,nspf
     do j=1,nspf
        mydot(i,j)=dot(projspfs(:,i),inspfs(:,j),spfsize)
     enddo
  enddo

  if (parorbsplit.eq.3) then
     call mympireduce(mydot,nspf**2)
  endif

!! with timefac
  call getconmat(thistime,ireduced,conmat)

  call MYGEMM('N','N',nspf,nspf,nspf,DATAONE,conmat,nspf,mydot,nspf,DATAZERO,mymat,nspf)

  call MYGEMM('N','N',spfsize,nspf,nspf,DATAONE,projspfs,spfsize,mymat,nspf,DATAZERO,outspfs,spfsize)

end subroutine op_gmat_thirdorder





!! only used for relax ; time not needed    
!!     this is actually the conjugate-transpose of the operator, even if cmctdh
!!     (not used)

subroutine actreducedconjg0(thistime,inspfs, projspfs, outspfs, ireduced, projflag,conflag)
  use parameters
  implicit none
  DATATYPE, intent(in) :: inspfs(spfsize,nspf),projspfs(spfsize,nspf)
  DATATYPE, intent(out) ::  outspfs(spfsize,nspf)
  integer,intent(in) :: ireduced,projflag,conflag
  real*8,intent(in) :: thistime
  DATATYPE :: ttempspfs(spfsize,nspf) !! AUTOMATIC

  ttempspfs(:,:)=ALLCON(inspfs(:,:))
  call actreduced0(thistime,ttempspfs, projspfs, outspfs, ireduced, projflag,conflag)
  outspfs=ALLCON(outspfs)

end subroutine actreducedconjg0


subroutine wmult(inspfs, outspfs, ireduced)
  use parameters
  use xxxmod
  implicit none
  DATATYPE, intent(in) :: inspfs(spfsize,nspf)
  DATATYPE, intent(out) :: outspfs(spfsize,nspf)
  integer,intent(in) :: ireduced
  DATATYPE :: spfmult(spfsize,nspf), workmult(spfsize,nspf),& !!AUTOMATIC
       spfinvr( spfsize,nspf ), spfinvrsq(  spfsize,nspf),spfproderiv(  spfsize,nspf )

!! sum over fast index reduced matrices, because doing spfinvrsq= reducedinvrsq * inspfs BUT 1) store in transposed order and 2) have to reverse the call in BLAS

  OFLWR "HMM CHECK WMULT.  CONMATEL ETC. also allocate arrays."; CFLST
  
  call MYGEMM('N', 'N', spfsize,nspf,nspf,DATAONE, inspfs, spfsize, yyy%reducedinvrsq(:,:,ireduced),nspf, DATAZERO, spfinvrsq(:,:), spfsize)
  call MYGEMM('N', 'N', spfsize,nspf,nspf,DATAONE, inspfs, spfsize, yyy%reducedinvr(:,:,ireduced),nspf, DATAZERO, spfinvr(:,:), spfsize)

  if ((nonuc_checkflag/=1)) then
     call MYGEMM('N', 'N', spfsize,nspf,nspf,DATAONE, inspfs, spfsize, yyy%reducedproderiv(:,:,ireduced),nspf, DATAZERO, spfproderiv(:,:), spfsize)
  endif
  
  spfmult=0.d0

  call mult_ke(spfinvrsq(:,:),spfmult(:,:),nspf,timingdir,notiming)

  call mult_pot(nspf,spfinvr(:,:),workmult(:,:))
  spfmult(:,:)=spfmult(:,:)+workmult(:,:)

  if ((nonuc_checkflag/=1)) then
     call noparorbsupport("another call op_yderiv")

     call op_yderiv(nspf,spfproderiv(:,:),workmult(:,:))
     spfmult(:,:)=spfmult(:,:) + workmult(:,:)

  endif

  call mult_reducedpot(1,nspf,inspfs(:,:),workmult(:,:),yyy%reducedpot(:,:,:,ireduced))

  spfmult(:,:)=spfmult(:,:)+workmult(:,:)

  call MYGEMM('N','N', spfsize,nspf,nspf,timefac, spfmult,spfsize, yyy%invdenmat(:,:,ireduced), nspf, DATAZERO, outspfs, spfsize)

end subroutine wmult

subroutine denmult(inspfs, outspfs, ireduced)
  use parameters
  use xxxmod
  implicit none
  integer,intent(in) :: ireduced
  DATATYPE,intent(in) :: inspfs(spfsize,nspf)
  DATATYPE,intent(out) :: outspfs(spfsize,nspf)

  call MYGEMM('N','N', spfsize,nspf,nspf,DATAONE, inspfs,spfsize, yyy%denmat(:,:,ireduced), nspf, DATAZERO, outspfs, spfsize)  !! ADDS TO OUTSPFS

  outspfs=outspfs/numelec

end subroutine denmult

