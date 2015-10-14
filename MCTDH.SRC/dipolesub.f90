

#include "Definitions.INC"

module dipolemod
  implicit none
  DATATYPE, allocatable :: xdipoleexpect(:),ydipoleexpect(:),zdipoleexpect(:)
  integer :: calledflag=0,xcalledflag=0
end module dipolemod

subroutine dipolesub_initial()
  use dipolemod
  use parameters
  implicit none
  allocate(xdipoleexpect(0:autosize), ydipoleexpect(0:autosize), zdipoleexpect(0:autosize))
end subroutine

subroutine dipolesub()
  use dipolemod
  use parameters
  use xxxmod
  implicit none
  DATATYPE :: xx,yy,zz
  integer :: imc,sflag
  integer, save :: lastouttime=0
  real*8 :: thistime

  if (mod(xcalledflag,autosteps).eq.0) then

     xdipoleexpect(calledflag)=0d0
     ydipoleexpect(calledflag)=0d0
     zdipoleexpect(calledflag)=0d0

     do imc=1,mcscfnum
        call dipolesub_one(yyy%cmfpsivec(astart(imc),0),yyy%cmfpsivec(spfstart,0), xx,yy,zz)

!! 101414 REAL-VALUED FOR HERM.

#ifndef CNORMFLAG
        xdipoleexpect(calledflag) = xdipoleexpect(calledflag) + real(xx,8)
        ydipoleexpect(calledflag) = ydipoleexpect(calledflag) + real(yy,8)
        zdipoleexpect(calledflag) = zdipoleexpect(calledflag) + real(zz,8)
#else
        xdipoleexpect(calledflag) = xdipoleexpect(calledflag) + xx
        ydipoleexpect(calledflag) = ydipoleexpect(calledflag) + yy
        zdipoleexpect(calledflag) = zdipoleexpect(calledflag) + zz
#endif

     enddo

     if (mod(calledflag,dipmodtime).eq.0.and.calledflag.gt.0) then

        thistime=calledflag*par_timestep*autosteps
        sflag=0
        if (floor(thistime/diptime).gt.lastouttime) then
           lastouttime=floor(thistime/diptime)
           sflag=1
        endif

        call dipolecall(calledflag, xdipoleexpect(0:),xdipfile,xdftfile,1,sflag)
        call dipolecall(calledflag, ydipoleexpect(0:),ydipfile,ydftfile,2,sflag)
        call dipolecall(calledflag, zdipoleexpect(0:),zdipfile,zdftfile,3,sflag)

     endif
     calledflag=calledflag+1
  endif
  xcalledflag=xcalledflag+1

end subroutine dipolesub


subroutine dipolesub_one(avector,inspfs,xdipole_expect,ydipole_expect,zdipole_expect)
  use xxxmod
  use parameters
  implicit none

  DATATYPE ::  inspfs(  spfsize, nspf ),  zdipole_expect, ydipole_expect, xdipole_expect, &
       avector(numr,firstconfig:lastconfig),tempvector(numr,firstconfig:lastconfig)
  DATATYPE :: dot,nullcomplex(1),dipoles(3), &
       tempspfs(spfsize,nspf), dipolemat(nspf,nspf),csum
  DATAECS :: rvector(numr)
  integer :: i
 
!! independent of R for now.  multiply by R for prolate  (R set to 1 for atom)
  call nucdipvalue(nullcomplex,dipoles)

  do i=1,numr
     tempvector(i,:)=avector(i,:)*bondpoints(i)
  enddo
  csum=dot(avector,tempvector,totadim)
  if (parconsplit.ne.0) then
     call mympireduceone(csum)
  endif
  dipoles(:)=dipoles(:)*csum

!! Z DIPOLE

  do i=1,nspf
     call mult_zdipole(inspfs(:,i),tempspfs(:,i))
  enddo
  call MYGEMM(CNORMCHAR,'N',nspf,nspf,spfsize,DATAONE, inspfs, spfsize, tempspfs, spfsize, DATAZERO, dipolemat, nspf)
  if (parorbsplit.eq.3) then
     call mympireduce(dipolemat(:,:),nspf**2)
  endif

  rvector(:)=bondpoints(:)
  call arbitraryconfig_mult_singles(dipolemat,rvector,avector,tempvector,numr)
  zdipole_expect=dot(avector,tempvector,totadim)
  if (parconsplit.ne.0) then
     call mympireduceone(zdipole_expect)
  endif
  zdipole_expect=zdipole_expect + dipoles(3)


!! Y DIPOLE

  do i=1,nspf
     call mult_ydipole(inspfs(:,i),tempspfs(:,i))
  enddo
  call MYGEMM(CNORMCHAR,'N',nspf,nspf,spfsize,DATAONE, inspfs, spfsize, tempspfs, spfsize, DATAZERO, dipolemat, nspf)
  if (parorbsplit.eq.3) then
     call mympireduce(dipolemat(:,:),nspf**2)
  endif

  rvector(:)=bondpoints(:)
  call arbitraryconfig_mult_singles(dipolemat,rvector,avector,tempvector,numr)
  ydipole_expect=dot(avector,tempvector,totadim)
  if (parconsplit.ne.0) then
     call mympireduceone(ydipole_expect)
  endif
  ydipole_expect=ydipole_expect + dipoles(2)


!! X DIPOLE

  do i=1,nspf
     call mult_xdipole(inspfs(:,i),tempspfs(:,i))
  enddo
  call MYGEMM(CNORMCHAR,'N',nspf,nspf,spfsize,DATAONE, inspfs, spfsize, tempspfs, spfsize, DATAZERO, dipolemat, nspf)
  if (parorbsplit.eq.3) then
     call mympireduce(dipolemat(:,:),nspf**2)
  endif

  rvector(:)=bondpoints(:)
  call arbitraryconfig_mult_singles(dipolemat,rvector,avector,tempvector,numr)
  xdipole_expect=dot(avector,tempvector,totadim)
  if (parconsplit.ne.0) then
     call mympireduceone(xdipole_expect)
  endif
  xdipole_expect=xdipole_expect + dipoles(1)

end subroutine dipolesub_one


subroutine dipolesub_final()
  use parameters
  use dipolemod
  implicit none
  deallocate( zdipoleexpect, xdipoleexpect,  ydipoleexpect)

end subroutine dipolesub_final


