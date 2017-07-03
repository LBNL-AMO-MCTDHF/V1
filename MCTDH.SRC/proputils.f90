
!! ALL MODULES

#include "Definitions.INC"


module orbprojectmod
contains

  subroutine project00(lowspf,highspf,inspfs, outspfs, prospfs)
    use parameters
    use opmod !! frozenspfs
    use mpisubmod
    implicit none
    integer,intent(in) :: lowspf,highspf
    DATATYPE,intent(in) :: inspfs(spfsize,lowspf:highspf),  prospfs(spfsize,nspf)
    DATATYPE,intent(out) :: outspfs(spfsize,lowspf:highspf)
    integer :: i,j,numspf
    DATATYPE :: mydot(nspf+numfrozen,lowspf:highspf+1),&
         tempprospfs(spfsize,nspf+numfrozen)            !! AUTOMATIC

    numspf=highspf-lowspf+1
    if (numspf.le.0) then
       call waitawhile()
       print *, "ACK ERRR PROJ00"
       call waitawhile()
       stop
    endif

    tempprospfs(:,1:nspf)=prospfs(:,:)

    if (numfrozen.ne.0) then
       tempprospfs(:,nspf+1:nspf+numfrozen)=frozenspfs(:,:)
    endif

!  if (noorthogflag.eq.0) then   !hardwire.  eliminated realproject.
!     call spf_orthogi t(tempprospfs,nspf+numfrozen, nulldouble)
!  endif

    mydot(:,:)=0d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i=lowspf,highspf
       do j=1,nspf+numfrozen
          mydot(j,i)=   dot(tempprospfs(:,j),inspfs(:,i),spfsize)
       enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

    if (parorbsplit.eq.3) then
       call mympireduce(mydot,nspf*(nspf+numfrozen))
    endif

    call MYGEMM('N','N',spfsize,numspf,nspf+numfrozen,DATAONE,&
         tempprospfs(:,:),spfsize,mydot(:,lowspf:highspf),&
         nspf+numfrozen,DATAZERO,outspfs(:,lowspf:highspf),spfsize)
     
  end subroutine project00


  subroutine project(inspfs, outspfs, prospfs)
    use parameters
    use orbgathersubmod
    implicit none
    DATATYPE,intent(in) :: inspfs(spfsize,nspf),  prospfs(spfsize,nspf)
    DATATYPE,intent(out) :: outspfs(spfsize,nspf)
    integer :: lowspf,highspf

    lowspf=1; highspf=nspf
    if (parorbsplit.eq.1) then
       call getOrbSetRange(lowspf,highspf)
    endif
    if (highspf-lowspf+1.gt.0) then
       call  project00(lowspf,highspf,inspfs(:,lowspf:highspf), &
            outspfs(:,lowspf:highspf), prospfs)
    endif
    if (parorbsplit.eq.1) then
       call mpiorbgather(outspfs,spfsize)
    endif

  end subroutine project


  subroutine project_onfrozen(inspf, outspf)
    use parameters
    use opmod !! frozenspfs
    use mpisubmod
    implicit none
    DATATYPE,intent(in) :: inspf(spfsize)
    DATATYPE,intent(out) :: outspf(spfsize)
    integer :: j
    DATATYPE :: mydot(numfrozen)

    outspf(:)=0;   mydot(:)=0

    if (numfrozen.gt.0) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j)
!$OMP DO SCHEDULE(STATIC)
       do j=1,numfrozen
          mydot(j)=   dot(frozenspfs(:,j),inspf(:),spfsize)
       enddo
!$OMP END DO
!$OMP END PARALLEL

       if (parorbsplit.eq.3) then
          call mympireduce(mydot,numfrozen)
       endif

       do j=1,numfrozen
          outspf(:) = outspf(:) + frozenspfs(:,j) * mydot(j)
       enddo
    endif

  end subroutine project_onfrozen

end module orbprojectmod


module orbmultsubmod
contains

!! needs factor of 1/r  for hamiltonian

  subroutine mult_pot(howmany,in, out)
    use parameters
    use opmod 
    implicit none
    integer,intent(in) :: howmany
    DATATYPE,intent(in) :: in(spfsize,howmany)
    DATATYPE,intent(out) :: out(spfsize,howmany)
    integer :: ii

    call op_contact(howmany,in,out)
    do ii=1,howmany
       out(:,ii)=out(:,ii)+in(:,ii)*pot(:)   !! NO INTERNUCLEAR REPULSION !!
    enddo

  end subroutine mult_pot

!! needs factor of 1/r  for hamiltonian

  subroutine mult_halfniumpot(howmany,in, out)
    use parameters
    use opmod  
    implicit none
    integer,intent(in) :: howmany
    DATATYPE,intent(in) :: in(spfsize,howmany)
    DATATYPE,intent(out) :: out(spfsize,howmany)
    integer :: ii
    do ii=1,howmany
       out(:,ii)=in(:,ii)*halfniumpot(:)
    enddo
  end subroutine mult_halfniumpot

  subroutine lenmultiply(howmany,spfin,spfout, myxtdpot,myytdpot,myztdpot)
    use mpimod
    use parameters
    implicit none
    integer,intent(in) :: howmany
    DATATYPE,intent(in) :: spfin(spfsize,howmany)
    DATATYPE,intent(out) :: spfout(spfsize,howmany)
    DATATYPE :: ttempspf(spfsize,howmany)              !! AUTOMATIC
    DATATYPE :: myxtdpot,myztdpot,myytdpot

    spfout(:,:)=0d0

    if (abs(myztdpot).ne.0d0) then
       call mult_zdipole(howmany,spfin(:,:),ttempspf(:,:),0)
       spfout(:,:)=spfout(:,:)+ttempspf(:,:)*myztdpot
    endif
    if (abs(myytdpot).ne.0d0) then
       call mult_ydipole(howmany,spfin(:,:),ttempspf(:,:),0)
       spfout(:,:)=spfout(:,:)+ttempspf(:,:)*myytdpot
    endif
    if (abs(myxtdpot).ne.0d0) then
       call mult_xdipole(howmany,spfin(:,:),ttempspf(:,:),0)
       spfout(:,:)=spfout(:,:)+ttempspf(:,:)*myxtdpot
    endif

  end subroutine lenmultiply

  subroutine op_frozen_exchange(lowspf,highspf,inspfs,outspfs)
    use opmod
    use parameters
    implicit none
    integer,intent(in) :: lowspf,highspf
    DATATYPE,intent(in) :: inspfs(spfsize,lowspf:highspf)
    DATATYPE,intent(out) :: outspfs(spfsize,lowspf:highspf)
    integer :: numspf
    numspf=highspf-lowspf+1
    if (numspf.gt.0) then
       call op_frozen_exchange0(numspf,inspfs,outspfs,frozenspfs,&
            numfrozen,spfmvals(lowspf:highspf))
    endif
  end subroutine op_frozen_exchange

!! transform to the other gauge (not exactly, taking real part of
!! dipole operator so that transformation is unitary)

  subroutine gauge_transform(invelflag,intime,numspf,inspfs,outspfs)
    use spfsize_parameters
    use pulsesubmod
    implicit none
    real*8,intent(in) :: intime
    integer,intent(in) :: numspf,invelflag
    DATATYPE,intent(in) :: inspfs(spfsize,numspf)
    DATATYPE,intent(out) :: outspfs(spfsize,numspf)
    DATATYPE,allocatable :: tempspf(:),tempspf2(:)
    DATATYPE :: pots(3)
    complex*16 :: phase
    integer :: ispf

    allocate(tempspf(spfsize),tempspf2(spfsize))
    tempspf=0; tempspf2=0
    outspfs(:,:)=inspfs(:,:)

    call vectdpot(intime,1,pots,-1)  !! A-vector velocity gauge, real part for unitary

    if (invelflag.ne.0) then
       phase=(0d0,-1d0)    !! velocity
    else
       phase=(0d0, 1d0)    !! length
    endif

    if (abs(pots(1)).gt.0d0) then
       tempspf=1d0
       call mult_xdipole(1,tempspf(:),tempspf2(:),1)
       tempspf(:)=exp(phase*pots(1)*tempspf2(:))
       do ispf=1,numspf
          outspfs(:,ispf)=outspfs(:,ispf) * tempspf(:)
       enddo
    endif

    if (abs(pots(2)).gt.0d0) then
       tempspf=1d0
       call mult_ydipole(1,tempspf(:),tempspf2(:),1)
       tempspf(:)=exp(phase*pots(2)*tempspf2(:))
       do ispf=1,numspf
          outspfs(:,ispf)=outspfs(:,ispf) * tempspf(:)
       enddo
    endif

    if (abs(pots(3)).gt.0d0) then
       tempspf=1d0
       call mult_zdipole(1,tempspf(:),tempspf2(:),1)
       tempspf(:)=exp(phase*pots(3)*tempspf2(:))
       do ispf=1,numspf
          outspfs(:,ispf)=outspfs(:,ispf) * tempspf(:)
       enddo
    endif

    deallocate(tempspf,tempspf2)

  end subroutine gauge_transform

end module orbmultsubmod


module getstuffmod
contains

subroutine get_stuff0(thistime,times)
  use parameters
  use dfconsubmod
  use meansubmod
  implicit none
  real*8,intent(in) :: thistime
  integer,intent(inout) :: times(20)
  integer :: itime,jtime

  call myclock(itime)

  call get_allden()
  call myclock(jtime);  times(2)=times(2)+jtime-itime;    itime=jtime

  call all_matel()
  call myclock(jtime);  times(1)=times(1)+jtime-itime;    itime=jtime

  if (constraintflag.ne.0) then
     call get_constraint(thistime)
  endif
  call myclock(jtime); times(7)=times(7)+jtime-itime;     itime=jtime

  if (drivingflag.ne.0) then
     call drivingtrans(thistime)
  endif
  call myclock(jtime); times(8)=times(8)+jtime-itime;     itime=jtime

  if (use_fockmatrix) then
     call get_fockmatrix()
  endif
  call get_reducedpot()
  if (numfrozen.gt.0) then
     call get_frexchange()
  endif
  call myclock(jtime);     times(3)=times(3)+jtime-itime

end subroutine get_stuff0


subroutine get_stuff(thistime)
  implicit none
  real*8,intent(in) :: thistime
  integer :: times(20)=0
  call get_stuff0(thistime,times)
end subroutine get_stuff

subroutine get_frexchange()
  use parameters
  use xxxmod !! frozenexchinvr
  use orbgathersubmod
  use orbmultsubmod
  use mpisubmod
  use opmod
  implicit none
  integer :: lowspf,highspf,numspf,ispf,jspf
  DATATYPE,allocatable :: frozenexchange(:,:)
  DATATYPE :: frodot

  if (numfrozen.eq.0) then
     return
  endif

  yyy%frozenexchinvr(:,:,0)=0

  if (exact_exchange.ne.0) then
     return
  endif

  allocate(frozenexchange(spfsize,nspf));     frozenexchange=0
  lowspf=1; highspf=nspf
  if (parorbsplit.eq.1) then
     call getorbsetrange(lowspf,highspf)
  endif
  numspf=highspf-lowspf+1
  frozenexchange(:,:)=0

  if (numspf.gt.0) then
     call MYGEMM('N','N',spfsize,numspf,nspf,DATAONE, yyy%cmfspfs(:,0), &
          spfsize,yyy%reducedinvr(:,lowspf:highspf,0),nspf, DATAZERO, &
          frozenexchange(:,lowspf:highspf),spfsize)
  endif
!     if (parorbsplit.eq.1) then
!        call mpiorbgather(frozenexchange(:,:),spfsize)
!     endif

  if (numspf.gt.0) then
     call op_frozen_exchange(lowspf,highspf,&
          frozenexchange(:,lowspf:highspf),&
          yyy%frozenexchinvr(:,lowspf:highspf,0))
  endif

  deallocate(frozenexchange)

! orthog to frozen !
  do ispf=lowspf,highspf
     do jspf=1,numfrozen
        frodot = dot(frozenspfs(:,jspf), yyy%frozenexchinvr(:,ispf,0), spfsize)
        if (parorbsplit.eq.3) then
           call mympireduceone(frodot)
        endif
        yyy%frozenexchinvr(:,ispf,0) = yyy%frozenexchinvr(:,ispf,0) - &
             frodot * frozenspfs(:,jspf) 
     enddo
  enddo

  do jspf=1,nspf
     do ispf=lowspf,highspf
        yyy%frozenexchmat(jspf,ispf,0) = &
             dot(yyy%cmfspfs((jspf-1)*spfsize+1:jspf*spfsize,0), yyy%frozenexchinvr(:,ispf,0), spfsize)
     enddo
  enddo
  if (parorbsplit.eq.3) then
     call mympireduce(yyy%frozenexchmat(:,:,0),nspf**2)
  endif

  if (parorbsplit.eq.1) then
     call mpiorbgather(yyy%frozenexchinvr(:,:,0),spfsize)
     call mpiorbgather(yyy%frozenexchmat(:,:,0),nspf)
  endif

end subroutine get_frexchange

end module getstuffmod




