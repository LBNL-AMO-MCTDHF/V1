
!! ROUTINES RELATED TO ORBITALS

#include "Definitions.INC"

subroutine noparorbsupport(intext)
  use parameters
  implicit none
  character*(*) :: intext
  if (parorbsplit.eq.3) then
     OFLWR "parorbsplit=3 not yet supported for ",intext; CFLST
  endif
end subroutine noparorbsupport



subroutine apply_spf_constraints(outspfs)
  use parameters
  implicit none
  DATATYPE :: outspfs(spfsize,nspf)

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
  DATATYPE :: outspfs(spfsize,nspf),inspfs(spfsmallsize,nspf)

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
  DATATYPE :: outspfs(spfsmallsize,nspf),inspfs(spfsize,nspf)

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


subroutine orthog_tofrozen(inspf)
  use parameters
  use opmod !! frozenspfs
  implicit none
  DATATYPE :: inspf(spfsize)
  logical :: orbparflag
  orbparflag=.false.
  if (parorbsplit.eq.3) then
     orbparflag=.true.
  endif
  if (numfrozen.gt.0) then
     call gramschmidt(spfsize, numfrozen, spfsize, frozenspfs(:,:), inspf(:),orbparflag)
  endif
end subroutine orthog_tofrozen


subroutine spf_orthogit(inspfs,error)
  use parameters
  use opmod !! frozenspfs
  implicit none

  integer :: i,j,ii
  real*8 :: error
  DATATYPE :: inspfs(spfsize,nspf),  dot, sum, &
       tempspfs(spfsize,nspf), ovl(nspf,nspf)
  logical :: orbparflag
  orbparflag=.false.
  if (parorbsplit.eq.3) then
     orbparflag=.true.
  endif

!! this is the error in orthog on input, not the final error.

  error=0.d0


  !! WTF.  if numfrozen is 0 this gives an error "Null pointer for frozenspfs"
  !!  on carver with pgf90.  Even though frozenspfs is never referenced.  Hmm.

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
        call MYGEMM('N', 'N', spfsize, nspf, nspf, DATAONE, tempspfs, spfsize, ovl, nspf, DATAZERO, inspfs, spfsize)   
     endif
  enddo
  

  if (abs(error) .gt. 0.000001d0) then
     OFLWR "Err spf ovl!! ", error;  CFLST
  endif
  
end subroutine spf_orthogit


subroutine spf_orthogit_gs(inspfs)
  use parameters
  use opmod !! frozenspfs
  implicit none
  DATATYPE :: inspfs(  spfsize, nspf   )
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
  DATATYPE :: inspfs(  spfsize, nspf   )
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

