

#include "Definitions.INC"


subroutine project(inspfs, outspfs, prospfs)
  use parameters
  use opmod !! frozenspfs
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
  

subroutine project00(lowspf,highspf,inspfs, outspfs, prospfs)
  use parameters
  use opmod !! frozenspfs
  implicit none
  integer,intent(in) :: lowspf,highspf
  DATATYPE,intent(in) :: inspfs(spfsize,lowspf:highspf),  prospfs(spfsize,nspf)
  DATATYPE,intent(out) :: outspfs(spfsize,lowspf:highspf)
  integer :: i,j,numspf
  DATATYPE :: dot
  DATATYPE :: mydot(nspf+numfrozen,lowspf:highspf+1),&
       tempprospfs(spfsize,nspf+numfrozen)            !! AUTOMATIC

  numspf=highspf-lowspf+1

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

  if (numspf.gt.0) then
     call MYGEMM('N','N',spfsize,numspf,nspf+numfrozen,DATAONE,&
          tempprospfs(:,:),spfsize,mydot(:,lowspf:highspf),&
          nspf+numfrozen,DATAZERO,outspfs(:,lowspf:highspf),spfsize)
  endif
     
end subroutine project00



subroutine project_onfrozen(inspf, outspf)
  use parameters
  use opmod !! frozenspfs
  implicit none
  DATATYPE,intent(in) :: inspf(spfsize)
  DATATYPE,intent(out) :: outspf(spfsize)
  integer :: j
  DATATYPE :: dot,mydot(numfrozen)

  outspf(:)=0;   mydot(:)=0

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
  use xxxmod !! frozenexchinvr
  implicit none
  integer :: lowspf,highspf,numspf
  DATATYPE,allocatable :: frozenexchange(:,:)

  if (numfrozen.gt.0) then

     allocate(frozenexchange(spfsize,nspf))

     lowspf=1; highspf=nspf
     if (parorbsplit.eq.1) then
        call getorbsetrange(lowspf,highspf)
     endif
     numspf=highspf-lowspf+1
     frozenexchange(:,:)=0
     if (numspf.gt.0) then
        call op_frozen_exchange(lowspf,highspf,&
             yyy%cmfspfs((lowspf-1)*spfsize+1:highspf*spfsize,0),&
             frozenexchange(:,lowspf:highspf))
     endif
     if (parorbsplit.eq.1) then
        call mpiorbgather(frozenexchange(:,:),spfsize)
     endif
     if (numspf.gt.0) then
        call MYGEMM('N','N',spfsize,numspf,nspf,DATAONE, frozenexchange(:,:),&
             spfsize,yyy%reducedinvr(:,lowspf:highspf,0),nspf, DATAZERO, &
             yyy%frozenexchinvr(:,lowspf:highspf,0),spfsize)
     endif
     if (parorbsplit.eq.1) then
        call mpiorbgather(yyy%frozenexchinvr(:,:,0),spfsize)
     endif

     deallocate(frozenexchange)

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
  implicit none
  real*8,intent(in) :: thistime
  integer,intent(inout) :: times(20)
  integer :: itime,jtime

  call system_clock(itime)

  call get_allden()
  call system_clock(jtime);  times(2)=times(2)+jtime-itime;    itime=jtime

  call all_matel()
  call system_clock(jtime);  times(1)=times(1)+jtime-itime;    itime=jtime

  if (constraintflag.ne.0) then
     call get_constraint(thistime)
  endif
  call system_clock(jtime); times(7)=times(7)+jtime-itime;     itime=jtime

  if (drivingflag.ne.0) then
     call drivingtrans(thistime)
  endif
  call system_clock(jtime); times(8)=times(8)+jtime-itime;     itime=jtime

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
  if (numspf.gt.0) then
     call op_frozen_exchange0(numspf,inspfs,outspfs,frozenspfs,&
          numfrozen,spfmvals(lowspf:highspf))
  endif
end subroutine op_frozen_exchange





