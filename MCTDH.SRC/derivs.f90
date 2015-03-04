
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

subroutine gbs_derivs(nullint,thistime,psi,psip,nulldbl,nullint2)
  use parameters
  implicit none
  integer :: nullint,nullint2
  real*8 :: thistime,nulldbl
  DATATYPE :: psi(psilength), psip(psilength)
  call all_derivs(thistime,psi,psip)
end subroutine gbs_derivs


subroutine all_derivs(thistime,xpsi, xpsip)
  use parameters
  use mpimod
  use xxxmod
  implicit none

  DATATYPE :: xpsi(psilength), xpsip(psilength)
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
     call actreduced(thistime,xpsi(spfstart),xpsi(spfstart),xpsip(spfstart),1)
     call system_clock(jtime);     times(5)=times(5)+jtime-itime

  call system_clock(itime)

  !! AVECTOR PART.

  do imc=1,mcscfnum
     call sparseconfigmult(xpsi(astart(imc)),xpsip(astart(imc)),yyy%cptr(0),yyy%sptr(0),1,1,1,1,thistime)
     if (allspinproject.ne.0) then
        call configspin_projectall(xpsip(astart(imc)),0)
     endif
     if (dfrestrictflag.ne.0) then
        call dfrestrict(xpsip(astart(imc)),numr)
     endif
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



subroutine gbs_linear_derivs(nullint,thistime,spfsin,spfsout,nulldbl,nullint2)
  use parameters
  implicit none
  DATATYPE :: spfsin(spfsize, nspf), spfsout(spfsize, nspf)
  integer :: nullint,nullint2
  real*8 :: thistime,nulldbl
  call spf_linear_derivs(thistime,spfsin,spfsout)
end subroutine gbs_linear_derivs

subroutine spf_linear_derivs(thistime,spfsin,spfsout)
  use parameters
  implicit none

  DATATYPE :: spfsin(spfsize, nspf), spfsout(spfsize, nspf)
  real*8 :: thistime
  call spf_linear_derivs0(thistime,spfsin,spfsout,1)

end subroutine spf_linear_derivs


!! MAIN DERIVATIVE ROUTINE FOR RK, EXPO

subroutine spf_linear_derivs0(thistime,spfsin,spfsout, allflag)
  use linearmod
  use parameters
  use xxxmod  !! frozenexchange
  implicit none

  DATATYPE :: spfsin(spfsize, nspf), spfsout(spfsize, nspf)
  real*8 :: thistime,facs(0:1)
  integer ::  jjj, allflag, ispf,jspf,ibot
  integer, save :: allochere=0
!  DATATYPE ::   spfmult(spfsize,nspf), spfmult2(spfsize,nspf), spfmult3(spfsize,nspf)
  DATATYPE,save,allocatable ::   spfmult(:,:),spfmult2(:,:),spfmult3(:,:)

  if (allochere.eq.0) then
     allocate(spfmult(spfsize,nspf), spfmult2(spfsize,nspf), spfmult3(spfsize,nspf))
  endif
  allochere=1


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
        call tauop(spfsin(:,:),spfmult(:,:), jjj,thistime)
        spfsout(:,:)=spfsout(:,:)+spfmult(:,:)*facs(jjj)
     enddo
     return
  endif


  do jjj=ibot,1
     call actreduced0(thistime,spfsin,spfsin,spfmult(:,:),jjj, allflag,allflag)
     spfsout(:,:)=spfsout(:,:)+spfmult(:,:)*facs(jjj)
  enddo

!! NOW INTERPOLATING EXCHANGE -- UNTESTED -- BEFORE HAD CMF EXCHANGE FOR STABILITY

  if (numfrozen.gt.0) then
     spfmult3(:,:)=0d0

     do jjj=ibot,1
        spfmult(:,:)=0d0
        do ispf=1,nspf
           do jspf=1,nspf
              spfmult(:,ispf) = spfmult(:,ispf) + &
                   yyy%frozenexchange(:,jspf,jjj) * yyy%reducedinvr(jspf,ispf,jjj) 
           enddo
        enddo

!! TIMEFAC HERE
        call MYGEMM('N','N', spfsize,nspf,nspf,timefac, spfmult(:,:),spfsize, yyy%invdenmat(:,:,jjj), nspf, DATAONE, spfmult2(:,:), spfsize)

        spfmult3(:,:)=spfmult3(:,:)+spfmult2(:,:)*facs(jjj)
     enddo

     call oneminusproject(spfmult3(:,:),spfmult(:,:),spfsin)
     spfsout(:,:)=spfsout(:,:)+spfmult(:,:)
  endif

end subroutine spf_linear_derivs0








subroutine driving_linear_derivs(thistime,spfsin,spfsout)
  use linearmod
  use parameters
  use xxxmod  !! driving orbs
  implicit none

  DATATYPE :: spfsin(spfsize, nspf),spfsout(spfsize, nspf),pots(3)
  real*8 :: thistime,facs(0:1),rsum
  integer ::  jjj, ibot
!  DATATYPE :: spfmult(spfsize,nspf)
  DATATYPE,save,allocatable :: spfmult(:,:)
  integer, save :: allochere=0

  if (allochere.eq.0) then
     allocate(spfmult(spfsize,nspf))
  endif
  allochere=1

  spfsout(:,:)=0d0

  if (effective_cmf_spfflag.eq.0) then
     return
  endif

  call vectdpot(thistime,pots)
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



!!  outspfs = (1-P) H inspfs

!! where P is projector onto opspfs
!! H is inverse denmat times reducedham

subroutine actreduced(thistime,inspfs,projspfs, outspfs, ireduced)
  use parameters
  implicit none
  real*8 :: thistime
  DATATYPE :: inspfs(spfsize,nspf),  outspfs(spfsize,nspf),  projspfs(spfsize,nspf)
  integer :: ireduced
  call actreduced0(thistime,inspfs, projspfs, outspfs, ireduced, 1,1)
end subroutine actreduced


!! MAIN ROUTINE TO OPERATE WITH REDUCED HAMILTONIAN (YOURSELF, NOT
!!  IN DIFF EQ SOLVER)


subroutine actreduced0(thistime,inspfs0, projspfs, outspfs, ireduced, projflag,conflag)
  use parameters
  use xxxmod
  implicit none

  DATATYPE :: outspfs(spfsize,nspf),inspfs0(spfsize, nspf), projspfs(spfsize,nspf)
  integer :: ireduced    ,    projflag,conflag
  real*8 :: thistime

  call actreduced00(thistime,inspfs0,projspfs,outspfs,projflag,conflag,1,0, yyy%cptr(ireduced),yyy%reducedinvr(:,:,ireduced), &
       yyy%reducedinvrsq(:,:,ireduced),yyy%reducedr(:,:,ireduced),yyy%reducedproderiv(:,:,ireduced), &
       yyy%reducedpot(:,:,:,ireduced),yyy%denmat(:,:,ireduced),yyy%invdenmat(:,:,ireduced))

end subroutine actreduced0

subroutine actreduced00(thistime,inspfs0, projspfs, outspfs, projflag,conflag, denandtimeflag, onlypulseflag, matrix_ptr, &
     reducedinvr,reducedinvrsq,reducedr,reducedproderiv, reducedpot,denmat,invdenmat)
  use parameters
  use mpimod
  use configptrmod
  implicit none

  DATATYPE, intent(in) :: inspfs0(spfsize, nspf), projspfs(spfsize,nspf)
  DATATYPE,intent(in) :: reducedinvr(nspf,nspf),reducedinvrsq(nspf,nspf),reducedr(nspf,nspf),reducedpot(reducedpotsize,nspf,nspf), &
       denmat(nspf,nspf), invdenmat(nspf,nspf), reducedproderiv(nspf,nspf)
  DATATYPE :: outspfs(spfsize,nspf)
  integer :: ispf, jspf,   itime, jtime, &
        projflag,getlen, conflag,lowspf,highspf, denandtimeflag,onlypulseflag
  integer, save :: times(0:20)=0,numcalledhere=0
  DATATYPE :: myxtdpot=0,  myytdpot=0, myztdpot=0, pots(3)=0d0
  real*8 :: thistime
  Type(CONFIGPTR) :: matrix_ptr

  integer, save :: allochere=0
  DATATYPE,save,allocatable ::  inspfs(:,:), tempmult(:),spfmult(:,:),myspf(:), &
       spfinvr( :,: ), spfr( :,: ),  spfinvrsq(  :,:),spfproderiv(  :,: )


  if (allochere.eq.0) then
     allocate(  inspfs(spfsize,nspf), tempmult(spfsize),spfmult(spfsize,nspf*2),myspf(spfsize), &
       spfinvr( spfsize,nspf ), spfr( spfsize,nspf ),  spfinvrsq(  spfsize,nspf),spfproderiv(  spfsize,nspf ))
  endif
  allochere=1


  if (tdflag.eq.1) then
     call vectdpot(thistime,pots)
     myxtdpot=pots(1);  myytdpot=pots(2);  myztdpot=pots(3);

  endif

  inspfs(:,:)=inspfs0(:,:)

  numcalledhere=numcalledhere+1
  call system_clock(itime)

 !! sum over fast index reduced matrices, because doing spfinvrsq= reducedinvrsq * inspfs BUT 1) store in transposed order and 2) have to reverse the call in BLAS

  if (numr.eq.1) then

     call MYGEMM('N', 'N', spfsize,nspf,nspf,DATAONE, inspfs, spfsize, denmat(:,:),nspf, DATAZERO, spfinvrsq(:,:), spfsize)

     spfinvr(:,:)=spfinvrsq(:,:)/bondpoints(1)

     if (tdflag.ne.0) then
        spfr(:,:)=spfinvrsq(:,:)*bondpoints(1)
     endif

     spfinvrsq(:,:)=spfinvrsq(:,:)/bondpoints(1)**2

  else


     call MYGEMM('N', 'N', spfsize,nspf,nspf,DATAONE, inspfs, spfsize, reducedinvrsq(:,:),nspf, DATAZERO, spfinvrsq(:,:), spfsize)
     call MYGEMM('N', 'N', spfsize,nspf,nspf,DATAONE, inspfs, spfsize, reducedinvr(:,:),nspf, DATAZERO, spfinvr(:,:), spfsize)

     
     if (tdflag.eq.1) then
        call MYGEMM('N', 'N', spfsize,nspf,nspf,DATAONE, inspfs, spfsize, reducedr(:,:),nspf, DATAZERO, spfr(:,:), spfsize)
     endif
     
     if ((nonuc_checkflag/=1).and.onlypulseflag.eq.0) then
        call MYGEMM('N', 'N', spfsize,nspf,nspf,DATAONE, inspfs, spfsize, reducedproderiv(:,:),nspf, DATAZERO, spfproderiv(:,:), spfsize)
     endif

  endif

  call system_clock(jtime);  times(1)=times(1)+jtime-itime

  lowspf=1; highspf=nspf
  if (parorbsplit.eq.1) then
     lowspf=firstmpiorb;     highspf=min(nspf,firstmpiorb+orbsperproc-1)
     if (multmanyflag.ne.0) then
        OFLWR "PARORBSPLIT=1/MULTMANY not compatible"; CFLST
     endif
  endif

  spfmult=0.d0

  if (multmanyflag.ne.0) then
     if (onlypulseflag.eq.0) then
        call system_clock(itime)
        call mult_ke(spfinvrsq(:,:),spfmult,nspf,timingdir,notiming)
        call system_clock(jtime);  times(2)=times(2)+jtime-itime
     endif
  endif

  do ispf=lowspf,highspf
     call system_clock(itime)
     if (onlypulseflag.eq.0) then
        if (multmanyflag.eq.0) then
           call mult_ke(spfinvrsq(:,ispf),tempmult,1,timingdir,notiming)
           spfmult(:,ispf) = spfmult(:,ispf) + tempmult(:)
           call system_clock(jtime);  times(2)=times(2)+jtime-itime;  call system_clock(itime)
        endif
        call mult_pot(spfinvr(:,ispf),tempmult)
        spfmult(:,ispf) = spfmult(:,ispf) + tempmult(:)
        call hatom_op(spfinvr(:,ispf),tempmult)
        spfmult(:,ispf)=spfmult(:,ispf)+tempmult(:)
        if (numfrozen.gt.0) then
           call op_frozenreduced(spfinvr(:,ispf),spfmult(:,ispf))  !! adds to spfmult DIRECT ONLY
        endif
     endif
     call system_clock(jtime);     times(3)=times(3)+jtime-itime;         call system_clock(itime)

     if (tdflag.eq.1) then

        select case (velflag)
        case (0)
           call lenmultiply(spfr(:,ispf),myspf(:), myxtdpot,myytdpot,myztdpot)
           spfmult(:,ispf)=spfmult(:,ispf)+myspf(:)
        case default
           call velmultiply(spfinvr(:,ispf),myspf(:), myxtdpot,myytdpot,myztdpot)
           spfmult(:,ispf)=spfmult(:,ispf)+myspf(:)
        end select
        
     endif  !! tdpot

     call system_clock(jtime);        times(4)=times(4)+jtime-itime;        call system_clock(itime)
  
     if ((nonuc_checkflag/=1).and.onlypulseflag.eq.0) then

        call noparorbsupport("another call mult_yderiv!!")
        call op_yderiv(spfproderiv(:,ispf),tempmult(:))
        spfmult(:,ispf)=spfmult(:,ispf)+tempmult(:)
        
     endif
     call system_clock(jtime);     times(5)=times(5)+jtime-itime; call system_clock(itime)

     if (onlypulseflag.eq.0) then

!! NOW ONLY OUTPUTS ONE, TAKES ALL
        call mult_reducedpot(inspfs(:,:),tempmult(:),ispf,reducedpot(:,:,:))
        spfmult(:,ispf)=spfmult(:,ispf) + tempmult(:)
     endif
     call system_clock(jtime);  times(6)=times(6)+jtime-itime;  
  enddo

  if (parorbsplit.eq.1) then
     call system_clock(itime)
     call mpiorbgather(spfmult,spfsize)
     call system_clock(jtime);  times(10)=times(10)+jtime-itime
  endif

!! WITH TIMEFAC
  if (denandtimeflag.ne.0) then
     call system_clock(itime)
     call MYGEMM('N','N', spfsize,nspf,nspf,timefac, spfmult,spfsize, invdenmat(:,:), nspf, DATAZERO, outspfs, spfsize)
     call system_clock(jtime);  times(7)=times(7)+jtime-itime
  endif

  if (projflag==1) then
     call system_clock(itime)
     call project(outspfs, spfmult, projspfs)
     outspfs(:,:) = outspfs(:,:) - spfmult(:,1:nspf)
     call system_clock(jtime);     times(8)=times(8)+jtime-itime;      
  endif  
     
  if (constraintflag/=0.and.conflag.ne.0) then
     call system_clock(itime)
     do ispf=1,nspf
        do jspf=1,nspf
           outspfs(:,ispf) = outspfs(:,ispf) + inspfs(:,jspf) * matrix_ptr%xconmatel(jspf,ispf) * timefac 
           if (tdflag.eq.1) then
              outspfs(:,ispf) = outspfs(:,ispf) + inspfs(:,jspf) * timefac * ( &
                   matrix_ptr%xconmatelxx(jspf,ispf) *myxtdpot + &
                   matrix_ptr%xconmatelyy(jspf,ispf) *myytdpot + &
                   matrix_ptr%xconmatelzz(jspf,ispf) *myztdpot )
           endif
        enddo
     enddo
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

end subroutine actreduced00




subroutine tauop(inspfs, outspfs, ireduced,thistime)
  use parameters
  use mpimod
  use xxxmod
  implicit none

  DATATYPE :: inspfs(spfsize,nspf), outspfs(spfsize,nspf)
  integer :: ispf, jspf, ireduced
  real*8 ::  thistime 
  DATATYPE :: myxtdpot=0,  myytdpot=0, myztdpot=0,pots(3)=0d0

  if (tdflag.eq.1) then
     call vectdpot(thistime,pots)
     myxtdpot=pots(1);     myytdpot=pots(2);     myztdpot=pots(3);
  endif

  outspfs(:,:)=0d0

  if (constraintflag.eq.0) then
     return
  endif

  do ispf=1,nspf
     do jspf=1,nspf
        outspfs(:,ispf) = outspfs(:,ispf) + inspfs(:,jspf) * yyy%cptr(ireduced)%xconmatel(jspf,ispf) * timefac 
        if (tdflag.eq.1) then
           outspfs(:,ispf) = outspfs(:,ispf) + inspfs(:,jspf) * timefac * ( &
                yyy%cptr(ireduced)%xconmatelxx(jspf,ispf) *myxtdpot + &
                yyy%cptr(ireduced)%xconmatelyy(jspf,ispf) *myytdpot + &
                yyy%cptr(ireduced)%xconmatelzz(jspf,ispf) *myztdpot )
        endif
     enddo
  enddo

end subroutine tauop


!! only used for relax ; time not needed    
!!     this is actually the conjugate-transpose of the operator, even if cmctdh
!!     (not used)

subroutine actreducedconjg0(thistime,inspfs, projspfs, outspfs, ireduced, projflag,conflag)
  use parameters
  implicit none

  DATATYPE, intent(in) :: inspfs(spfsize,nspf),projspfs(spfsize,nspf)
  DATATYPE ::  outspfs(spfsize,nspf)
  integer :: ireduced,projflag,conflag
  real*8 :: thistime
  integer, save :: allochere=0
  DATATYPE,save,allocatable ::  ttempspfs(:,:)

  if (allochere.eq.0) then
     allocate(ttempspfs(spfsize,nspf))
  endif
  allochere=1

!!!  OFLWR "You should check me!! actreducedconjg.  Maybe ok."; CFLST  why would it not be ok, i forget

  ttempspfs(:,:)=ALLCON(inspfs(:,:))
  call actreduced0(thistime,ttempspfs, projspfs, outspfs, ireduced, projflag,conflag)
  outspfs=ALLCON(outspfs)

end subroutine actreducedconjg0

subroutine wmult(inspfs, outspfs, ireduced)
  use parameters
  use xxxmod
  implicit none

  DATATYPE, intent(in) :: inspfs(spfsize,nspf)
  DATATYPE ::  outspfs(spfsize,nspf),  spfmult(spfsize,nspf), tempbigmult(spfsize,nspf), tempmult(spfsize)
  integer :: ispf,ireduced
  DATATYPE :: spfinvr( spfsize,nspf ), spfinvrsq(  spfsize,nspf),spfproderiv(  spfsize,nspf )

!! sum over fast index reduced matrices, because doing spfinvrsq= reducedinvrsq * inspfs BUT 1) store in transposed order and 2) have to reverse the call in BLAS

OFLWR "HMM CHECK WMULT.  CONMATEL ETC. also allocate arrays."; CFLST
  
  call MYGEMM('N', 'N', spfsize,nspf,nspf,DATAONE, inspfs, spfsize, yyy%reducedinvrsq(:,:,ireduced),nspf, DATAZERO, spfinvrsq(:,:), spfsize)
  call MYGEMM('N', 'N', spfsize,nspf,nspf,DATAONE, inspfs, spfsize, yyy%reducedinvr(:,:,ireduced),nspf, DATAZERO, spfinvr(:,:), spfsize)

  if ((nonuc_checkflag/=1)) then
     call MYGEMM('N', 'N', spfsize,nspf,nspf,DATAONE, inspfs, spfsize, yyy%reducedproderiv(:,:,ireduced),nspf, DATAZERO, spfproderiv(:,:), spfsize)
  endif
  
  spfmult=0.d0

  if (multmanyflag.ne.0) then
     call mult_ke(spfinvrsq(:,:),spfmult(:,:),nspf,timingdir,notiming); 
  else
     do ispf=1,nspf
        call mult_ke(spfinvrsq(:,ispf),spfmult(:,ispf),1,timingdir,notiming); 
     enddo
  endif

  do ispf=1,nspf
     call mult_pot(spfinvr(:,ispf),tempmult(:));     spfmult(:,ispf)=spfmult(:,ispf)+tempmult(:)
  enddo

  if ((nonuc_checkflag/=1)) then

     call noparorbsupport("another call op_yderiv")

     do ispf=1,nspf
        call op_yderiv(spfproderiv(:,ispf),tempmult(:));        spfmult(:,ispf)=spfmult(:,ispf) + tempmult(:)
     enddo
  endif

  do ispf=1,nspf
!! NOW OUTPUTS ONLY ONE (TAKES ALL)
     call mult_reducedpot(inspfs,tempbigmult(:,ispf),ispf,yyy%reducedpot(:,:,:,ireduced))
  enddo
  spfmult(:,:)=spfmult(:,:)+tempbigmult(:,:)

  if (whichquad.ne.1) then
     outspfs=spfmult
  else
     call MYGEMM('N','N', spfsize,nspf,nspf,timefac, spfmult,spfsize, yyy%invdenmat(:,:,ireduced), nspf, DATAZERO, outspfs, spfsize)
  endif

end subroutine wmult


subroutine denmult(inspfs, outspfs, ireduced)
  use parameters
  use xxxmod
  implicit none

  integer :: ireduced
  DATATYPE  :: inspfs(spfsize,nspf),outspfs(spfsize,nspf)

  call MYGEMM('N','N', spfsize,nspf,nspf,DATAONE, inspfs,spfsize, yyy%denmat(:,:,ireduced), nspf, DATAZERO, outspfs, spfsize)  !! ADDS TO OUTSPFS

  outspfs=outspfs/numelec

end subroutine denmult

