

#include "Definitions.INC"


subroutine oneminusproject(inspfs, outspfs, prospfs)
  use parameters
  implicit none  
  DATATYPE,intent(in) :: inspfs(spfsize,nspf),  prospfs(spfsize,nspf)
  DATATYPE,intent(out) :: outspfs(spfsize,nspf)

  call project(inspfs, outspfs, prospfs)
  outspfs=inspfs-outspfs
end subroutine


subroutine project(inspfs, outspfs, prospfs)
  use parameters
  use opmod !! frozenspfs
  implicit none

  DATATYPE,intent(in) :: inspfs(spfsize,nspf),  prospfs(spfsize,nspf)
  DATATYPE,intent(out) :: outspfs(spfsize,nspf)
  integer :: i,j,lowspf,highspf
  DATATYPE :: dot,mydot(nspf+numfrozen,nspf)
  DATATYPE :: tempprospfs(spfsize,nspf+numfrozen)  !! AUTOMATIC

  tempprospfs(:,1:nspf)=prospfs(:,:)

  if (numfrozen.ne.0) then
     tempprospfs(:,nspf+1:nspf+numfrozen)=frozenspfs(:,:)
  endif

!  if (noorthogflag.eq.0) then   !hardwire.  eliminated realproject.
!     call spf_orthogi t(tempprospfs,nspf+numfrozen, nulldouble)
!  endif

  lowspf=1; highspf=nspf
  if (parorbsplit.eq.1) then
     call getOrbSetRange(lowspf,highspf)
  endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
  do i=lowspf,highspf
     do j=1,nspf+numfrozen
        mydot(j,i)=   dot(tempprospfs(:,j),inspfs(:,i),spfsize)
     enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

  if (parorbsplit.eq.1) then
     call mpiorbgather(mydot,nspf+numfrozen)
  endif

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
  DATATYPE,intent(in) :: inspf(spfsize)
  DATATYPE,intent(out) :: outspf(spfsize)
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
  integer :: lowspf,highspf
  if (numfrozen.gt.0) then
     lowspf=1; highspf=nspf
     if (parorbsplit.eq.1) then
        call getorbsetrange(lowspf,highspf)
     endif
     yyy%frozenexchange(:,:,0)=0
     call op_frozen_exchange(lowspf,highspf,&
          yyy%cmfpsivec(spfstart+(lowspf-1)*spfsize:spfstart+highspf*spfsize-1,0),&
          yyy%frozenexchange(:,lowspf:highspf,0))
     if (parorbsplit.eq.1) then
        call mpiorbgather(yyy%frozenexchange(:,:,0),spfsize)
     endif
  endif
  
end subroutine get_frexchange

subroutine get_stuff(thistime)
  implicit none
  real*8,intent(in) :: thistime
  integer :: times(20)=0
  call get_stuff0(thistime,times)
end subroutine get_stuff


subroutine get_stuff0(thistime,times)
  use parameters
  use xxxmod !! frozenexchange
  implicit none
  real*8,intent(in) :: thistime
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


subroutine mult_reke(howmany,in,out)
  use parameters
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  DATATYPE :: work(spfsize,howmany)
  work(:,:)=ALLCON(in(:,:))
  call mult_ke(work,out,howmany,"booga",2)
  work=ALLCON(out)
  call mult_ke(in,out,howmany,"booga",2)  
  work=work+out
  out=work/2d0
end subroutine mult_reke


subroutine mult_imke(howmany,in,out)
  use parameters
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  DATATYPE :: work(spfsize,howmany)
  work(:,:)=ALLCON(in(:,:))
  call mult_ke(work,out,howmany,"booga",2)
  work=ALLCON(out)
  call mult_ke(in,out,howmany,"booga",2)  
  work=out-work
  out=work/(0d0,2d0)
end subroutine mult_imke

subroutine op_reyderiv(howmany,in,out)
  use parameters
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  DATATYPE :: work(spfsize,howmany)
  work(:,:)=ALLCON(in(:,:))
  call op_yderiv(howmany,work,out)
  work=ALLCON(out)
  call op_yderiv(howmany,in,out)
  work=work+out
  out=work/2d0
end subroutine op_reyderiv


subroutine op_imyderiv(howmany,in,out)
  use parameters
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  DATATYPE :: work(spfsize,howmany)
  work(:,:)=ALLCON(in(:,:))
  call op_yderiv(howmany,work,out)
  work=ALLCON(out)
  call op_yderiv(howmany,in,out)
  work=out-work
  out=work/(0d0,2d0)
end subroutine op_imyderiv


!! needs factor of 1/r  for hamiltonian

subroutine mult_impot(howmany,in, out)
  use parameters
  use opmod 
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  integer :: ii
  do ii=1,howmany
     out(:,ii)=in(:,ii)*imag((0d0,0d0)+pot(:))   !! NO INTERNUCLEAR REPULSION !!
  enddo
end subroutine mult_impot



!! needs factor of 1/r  for hamiltonian

subroutine mult_repot(howmany,in, out)
  use parameters
  use opmod 
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  integer :: ii
  do ii=1,howmany
     out(:,ii)=in(:,ii)*real(pot(:),8)   !! NO INTERNUCLEAR REPULSION !!
  enddo
end subroutine mult_repot


subroutine mult_pot(howmany,in, out)
  use parameters
  use opmod 
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  integer :: ii
  do ii=1,howmany
     out(:,ii)=in(:,ii)*pot(:)   !! NO INTERNUCLEAR REPULSION !!
  enddo
end subroutine mult_pot


!! needs factor of 1/r  for hamiltonian

subroutine mult_imhalfniumpot(howmany,in, out)
  use parameters
  use opmod  
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  integer :: ii
  do ii=1,howmany
     out(:,ii)=in(:,ii)*imag((0d0,0d0)+halfniumpot(:))
  enddo
end subroutine mult_imhalfniumpot


subroutine mult_rehalfniumpot(howmany,in, out)
  use parameters
  use opmod  
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  integer :: ii
  do ii=1,howmany
     out(:,ii)=in(:,ii)*real(halfniumpot(:),8)
  enddo
end subroutine mult_rehalfniumpot


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
  call op_frozen_exchange0(numspf,inspfs,outspfs,frozenspfs,numfrozen,spfmvals(lowspf:highspf))
end subroutine op_frozen_exchange





