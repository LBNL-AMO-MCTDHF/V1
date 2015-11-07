

#include "Definitions.INC"


subroutine oneminusproject(inspfs, outspfs, prospfs)
  use parameters
  implicit none  
  DATATYPE :: inspfs(spfsize,nspf), &
       outspfs(spfsize,nspf), &
       prospfs(spfsize,nspf)
  call project(inspfs, outspfs, prospfs)
  outspfs=inspfs-outspfs
end subroutine


subroutine project(inspfs, outspfs, prospfs)
  use parameters
  use opmod !! frozenspfs
  implicit none

  DATATYPE :: inspfs(spfsize, nspf), &
       outspfs(spfsize, nspf), &
       prospfs(spfsize, nspf)
  integer :: i,j
  DATATYPE :: dot,mydot(nspf+numfrozen,nspf)
  DATATYPE :: tempprospfs(spfsize,nspf+numfrozen)

  tempprospfs(:,1:nspf)=prospfs(:,:)

  if (numfrozen.ne.0) then
     tempprospfs(:,nspf+1:nspf+numfrozen)=frozenspfs(:,:)
  endif

!  if (noorthogflag.eq.0) then   !hardwire.  eliminated realproject.
!     call spf_orthogi t(tempprospfs,nspf+numfrozen, nulldouble)
!  endif

!!$  do i=1,nspf
!!$     outspfs(:,i)=0.d0
!!$     do j=1,nspf+numfrozen
!!$        outspfs(:,i) = outspfs(:,i) + tempprospfs(:,j) * &
!!$             dot(tempprospfs(:,j),inspfs(:,i),spfsize)
!!$     enddo
!!$  enddo

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
  do i=1,nspf
     do j=1,nspf+numfrozen
        mydot(j,i)=   dot(tempprospfs(:,j),inspfs(:,i),spfsize)
     enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

  if (parorbsplit.eq.3) then
     call mympireduce(mydot,nspf*(nspf+numfrozen))
  endif

  outspfs(:,:)=0
  do i=1,nspf
     do j=1,nspf+numfrozen
        outspfs(:,i) = outspfs(:,i) + tempprospfs(:,j) * &
             mydot(j,i)
     enddo
  enddo

end subroutine project



subroutine project_onfrozen(inspf, outspf)
  use parameters
  use opmod !! frozenspfs
  implicit none
  DATATYPE :: inspf(spfsize),outspf(spfsize)
  integer :: j
  DATATYPE :: dot,mydot(numfrozen)

  outspf(:)=0

  if (numfrozen.gt.1) then
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



subroutine get_frexchange()
  use parameters
  use xxxmod !! frozenexchange
  implicit none
  if (numfrozen.gt.0) then
     yyy%frozenexchange(:,:,0)=0
     call op_frozen_exchange(yyy%cmfpsivec(spfstart,0),yyy%frozenexchange(:,:,0))
  endif
end subroutine get_frexchange

subroutine get_stuff(thistime)
  implicit none
  real*8 :: thistime
  integer :: times(20)=0
  call get_stuff0(thistime,times)
end subroutine get_stuff


subroutine get_stuff0(thistime,times)
  use parameters
  use xxxmod !! frozenexchange
  implicit none

  real*8 :: thistime
  integer :: times(20),itime,jtime

  call system_clock(itime)

  call get_allden()
  call system_clock(jtime);  times(2)=times(2)+jtime-itime;  call system_clock(itime)

  call all_matel()

  call system_clock(jtime);  times(1)=times(1)+jtime-itime;

  if (constraintflag.ne.0) then
     call system_clock(itime)
     call get_constraint(thistime)
     call system_clock(jtime); times(7)=times(7)+jtime-itime;     
  endif
  if (drivingflag.ne.0) then
     call system_clock(itime)
     call drivingtrans(thistime)
     call system_clock(jtime); times(8)=times(8)+jtime-itime;     
  endif

  call system_clock(itime)

  call get_reducedpot()
  if (numfrozen.gt.0) then
     call get_frexchange()
  endif
  call system_clock(jtime);     times(3)=times(3)+jtime-itime


end subroutine get_stuff0




!! needs factor of 1/r  for hamiltonian

subroutine mult_impot(in, out)
  use parameters
  use opmod 
  implicit none
  DATATYPE :: in(spfsize), out(spfsize)
  out(:)=in(:)*imag((0d0,0d0)+pot(:))   !! NO INTERNUCLEAR REPULSION !!
end subroutine mult_impot



!! needs factor of 1/r  for hamiltonian

subroutine mult_repot(in, out)
  use parameters
  use opmod 
  implicit none
  DATATYPE :: out(spfsize), in(spfsize)
  out(:)=in(:)*real(pot(:),8)   !! NO INTERNUCLEAR REPULSION !!
end subroutine mult_repot


subroutine mult_pot(in, out)
  use parameters
  use opmod 
  implicit none
  DATATYPE :: in(spfsize),out(spfsize)
  out(:)=in(:)*pot(:)   !! NO INTERNUCLEAR REPULSION !!
end subroutine mult_pot


!! needs factor of 1/r  for hamiltonian

subroutine mult_imhalfniumpot(in, out)
  use parameters
  use opmod  
  implicit none
  DATATYPE :: in(spfsize),out(spfsize)
  out(:)=in(:)*imag((0d0,0d0)+halfniumpot(:))
end subroutine mult_imhalfniumpot


subroutine mult_rehalfniumpot(in, out)
  use parameters
  use opmod  
  implicit none
  DATATYPE :: in(spfsize),out(spfsize)
  out(:)=in(:)*real(halfniumpot(:),8)
end subroutine mult_rehalfniumpot


subroutine lenmultiply(spfin,spfout, myxtdpot,myytdpot,myztdpot)
  use mpimod
  use parameters
  implicit none
  DATATYPE :: ttempspf(spfsize), spfin(spfsize), spfout(spfsize)
  DATATYPE :: myxtdpot,myztdpot,myytdpot

  spfout(:)=0d0

  if (abs(myztdpot).ne.0d0) then
     call mult_zdipole(spfin(:),ttempspf(:))
     spfout(:)=spfout(:)+ttempspf(:)*myztdpot
  endif
  if (abs(myytdpot).ne.0d0) then
     call mult_ydipole(spfin(:),ttempspf(:))
     spfout(:)=spfout(:)+ttempspf(:)*myytdpot
  endif
  if (abs(myxtdpot).ne.0d0) then
     call mult_xdipole(spfin(:),ttempspf(:))
     spfout(:)=spfout(:)+ttempspf(:)*myxtdpot
  endif

end subroutine lenmultiply


subroutine call_frozen_matels()
  use opmod
  use parameters
  implicit none
  call call_frozen_matels0(frozenspfs(:,:),numfrozen,frozenkediag,frozenpotdiag)  !! returns diags; has matels in twoemod
end subroutine call_frozen_matels


subroutine call_frozen_exchange(inspfs,outspfs)
  use opmod
  use parameters
  implicit none
  DATATYPE :: inspfs(totspfdim), outspfs(totspfdim)
  call call_frozen_exchange0(inspfs,outspfs,frozenspfs,numfrozen)
end subroutine call_frozen_exchange





