
!! ALL ONE MODULE

!! ROUTINES RELATED TO ORBITALS

#include "Definitions.INC"


module spfsubmod
contains

subroutine apply_spf_constraints(outspfs)
  use parameters
  implicit none
  DATATYPE,intent(inout) :: outspfs(spfsize,nspf)

  if (spfrestrictflag==1) then
     call restrict_spfs0(outspfs,nspf,spfmvals,0)   
     if (spfugrestrict.ne.0) then
        call ugrestrict_spfs0(outspfs,nspf,spfugvals,0)
     endif
  endif
end subroutine apply_spf_constraints


subroutine spfs_expand(inspfs,outspfs)
  use parameters
  implicit none
  DATATYPE,intent(in) :: inspfs(spfsmallsize,nspf)
  DATATYPE,intent(out) :: outspfs(spfsize,nspf)

  if (spfugrestrict.ne.0) then
     call bothexpand_spfs(inspfs,outspfs,nspf,spfmvals,spfugvals)
  else
     if (spfrestrictflag.ne.0) then
        call mexpand_spfs(inspfs,outspfs,nspf,spfmvals)
     else
        outspfs(:,:)=inspfs(:,:)
     endif
  end if
end subroutine spfs_expand

subroutine spfs_compact(inspfs,outspfs)
  use parameters
  implicit none
  DATATYPE,intent(in) :: inspfs(spfsize,nspf)
  DATATYPE,intent(out) :: outspfs(spfsmallsize,nspf)

  if (spfugrestrict.ne.0) then
     call bothcompact_spfs(inspfs,outspfs,nspf,spfmvals,spfugvals)
  else
     if (spfrestrictflag.ne.0) then
        call mcompact_spfs(inspfs,outspfs,nspf,spfmvals)
     else
        outspfs(:,:)=inspfs(:,:)
     endif
  end if
end subroutine spfs_compact


subroutine spfs_expand_local(inspfs,outspfs)
  use parameters
  use mpi_orbsetmod
  implicit none
  DATATYPE,intent(in) :: inspfs(spfsmallsize,firstmpiorb:firstmpiorb+orbsperproc-1)
  DATATYPE,intent(out) :: outspfs(spfsize,firstmpiorb:firstmpiorb+orbsperproc-1)
  integer :: lowspf,highspf,numspf

  outspfs=0

  if (firstmpiorb.le.nspf) then
     lowspf=firstmpiorb
     highspf=min(nspf,firstmpiorb+orbsperproc-1)
     numspf=highspf-lowspf+1
     if (spfugrestrict.ne.0) then
        call bothexpand_spfs(inspfs(:,lowspf:highspf),outspfs(:,lowspf:highspf),&
             numspf,spfmvals(lowspf:highspf),spfugvals(lowspf:highspf))
     else
        if (spfrestrictflag.ne.0) then
           call mexpand_spfs(inspfs(:,lowspf:highspf),outspfs(:,lowspf:highspf),&
                numspf,spfmvals(lowspf:highspf))
        else
           outspfs(:,lowspf:highspf)=inspfs(:,lowspf:highspf)
        endif
     end if
  endif

end subroutine spfs_expand_local

subroutine spfs_compact_local(inspfs,outspfs)
  use parameters
  use mpi_orbsetmod
  implicit none
  DATATYPE,intent(in) :: inspfs(spfsize,firstmpiorb:firstmpiorb+orbsperproc-1)
  DATATYPE,intent(out) :: outspfs(spfsmallsize,firstmpiorb:firstmpiorb+orbsperproc-1)
  integer :: lowspf,highspf,numspf

  outspfs=0

  if (firstmpiorb.le.nspf) then
     lowspf=firstmpiorb
     highspf=min(nspf,firstmpiorb+orbsperproc-1)
     numspf=highspf-lowspf+1
     if (spfugrestrict.ne.0) then
        call bothcompact_spfs(inspfs(:,lowspf:highspf),outspfs(:,lowspf:highspf),&
             numspf,spfmvals(lowspf:highspf),spfugvals(lowspf:highspf))
     else
        if (spfrestrictflag.ne.0) then
           call mcompact_spfs(inspfs(:,lowspf:highspf),outspfs(:,lowspf:highspf),&
                numspf,spfmvals(lowspf:highspf))
        else
           outspfs(:,lowspf:highspf)=inspfs(:,lowspf:highspf)
        endif
     end if
  endif

end subroutine spfs_compact_local


subroutine spf_orthogit(inspfs,error)
  use parameters
  use opmod !! frozenspfs
  use matsubmod
  use orbgathersubmod
  use mpisubmod
  use utilmod
  implicit none
  real*8,intent(out) :: error
  DATATYPE,intent(inout) :: inspfs(spfsize,nspf)
  DATATYPE :: tempspfs(spfsize,nspf), ovl(nspf,nspf)  !! AUTOMATIC
  DATATYPE :: sum
  logical :: orbparflag
  integer :: i,j,ii,lowspf,highspf,numspf

  orbparflag=.false.
  if (parorbsplit.eq.3) then
     orbparflag=.true.
  endif

!! this is the error in orthog on input, not the final error.

  error=0.d0

  if (numfrozen.gt.0) then
     do i=1,nspf
        call gramschmidt(spfsize, numfrozen, spfsize, frozenspfs(:,:), inspfs(:,i),orbparflag)
     enddo
  endif

  do ii=1,3

     do i=1,nspf
        do j=1,nspf
           ovl(i,j) = dot(inspfs(:,i),inspfs(:,j),spfsize)
        enddo
     enddo
     if (orbparflag) then
        call mympireduce(ovl,nspf**2)
     endif
     error=0d0
     do i=1,nspf
        do j=1,nspf
           sum=ovl(i,j)
           if (i==j) then 
              sum=sum-1.d0
           endif
           if (abs(sum).gt.error) then
              error=abs(sum)
           endif
        enddo
     enddo
     
     if (ii.lt.3) then
        call orthogmat(ovl,nspf)
        tempspfs=inspfs

        lowspf=1; highspf=nspf
        if (parorbsplit.eq.1) then
           call getOrbSetRange(lowspf,highspf)
        endif
        numspf=highspf-lowspf+1
        if (numspf.gt.0) then
           call MYGEMM('N', 'N', spfsize, numspf, nspf, DATAONE, tempspfs, spfsize, &
                ovl(:,lowspf:highspf), nspf, DATAZERO, inspfs(:,lowspf:highspf), spfsize)   
        endif
        if (parorbsplit.eq.1) then
           call mpiorbgather(inspfs,spfsize)
        endif
     endif
  enddo
  
!!$  if (abs(error) .gt. 0.000001d0) then
!!$     OFLWR "Err spf ovl!! ", error;  CFLST
!!$  endif
  
end subroutine spf_orthogit


subroutine spf_orthogit_gs(inspfs)
  use parameters
  use opmod !! frozenspfs
  use utilmod
  implicit none
  DATATYPE,intent(inout) :: inspfs(  spfsize, nspf   )
  integer :: i 
  logical :: orbparflag
  orbparflag=.false.
  if (parorbsplit.eq.3) then
     orbparflag=.true.
  endif

  if (numfrozen.gt.0) then
     do i=1,nspf
        call gramschmidt(spfsize, numfrozen, spfsize, frozenspfs(:,:), inspfs(:,i),orbparflag)
     enddo
  endif

  do i=1,nspf
     call gramschmidt(spfsize, i-1, spfsize, inspfs(:,:), inspfs(:,i),orbparflag)
  enddo


end subroutine spf_orthogit_gs


subroutine mess_with_spfs(inspfs)
  use parameters
  implicit none
  DATATYPE,intent(inout) :: inspfs(  spfsize, nspf   )
  integer :: ispf, ii
  real*8 :: nextran, nulldouble

!!  TO BE SAFE.  DON'T MESS.
  call rand_init(-0.0542899d0)

!! this would really screw things up for cnorm

!#ifndef CNORMFLAG  
!  inspfs(:,1) = inspfs(:,1) * (0.5d0, 0.88d0) !! ok for imp conv
!  if (nspf.gt.1) then
!     inspfs(:,2) = inspfs(:,2) * (-0.2d0, 0.88d0) !! ok for imp conv
!  endif
!#endif

  do ispf=1,nspf
     do ii=1,spfsize
        inspfs(ii,ispf)=inspfs(ii,ispf)+messamount*nextran()
     enddo
  enddo
  call spf_orthogit(inspfs, nulldouble)
  OFLWR "mess: orthog ",nulldouble; CFL

  call apply_spf_constraints(inspfs)

  call spf_orthogit(inspfs, nulldouble)
  OFLWR "mess: orthog2 ",nulldouble; CFL

end subroutine mess_with_spfs

end module spfsubmod


