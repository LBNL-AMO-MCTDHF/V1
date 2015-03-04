

#include "Definitions.INC"

module dipolemod
  implicit none
  DATATYPE, allocatable :: xdipoleexpect(:),ydipoleexpect(:),zdipoleexpect(:)
  integer :: calledflags=0,xcalledflags=0
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

  if (mod(xcalledflags,autosteps).eq.0) then

     xdipoleexpect(calledflags)=0d0
     ydipoleexpect(calledflags)=0d0
     zdipoleexpect(calledflags)=0d0

     do imc=1,mcscfnum
        call dipolesub_one(yyy%cmfpsivec(astart(imc),0),yyy%cmfpsivec(spfstart,0), xx,yy,zz)

!! 101414 REAL-VALUED FOR HERM.

#ifndef CNORMFLAG
        xdipoleexpect(calledflags) = xdipoleexpect(calledflags) + real(xx,8)
        ydipoleexpect(calledflags) = ydipoleexpect(calledflags) + real(yy,8)
        zdipoleexpect(calledflags) = zdipoleexpect(calledflags) + real(zz,8)
#else
        xdipoleexpect(calledflags) = xdipoleexpect(calledflags) + xx
        ydipoleexpect(calledflags) = ydipoleexpect(calledflags) + yy
        zdipoleexpect(calledflags) = zdipoleexpect(calledflags) + zz
#endif

     enddo



     if (mod(calledflags,dipmodtime).eq.0) then

        thistime=calledflags*par_timestep*autosteps
        sflag=0
        if (floor(thistime/diptime).gt.lastouttime) then
           lastouttime=floor(thistime/diptime)
           sflag=1
        endif

        call dipolecall(calledflags, xdipoleexpect(0:),xdipfile,xdftfile,2,sflag)
        call dipolecall(calledflags, ydipoleexpect(0:),ydipfile,ydftfile,3,sflag)
        call dipolecall(calledflags, zdipoleexpect(0:),zdipfile,zdftfile,1,sflag)
     endif
     calledflags=calledflags+1
  endif
  xcalledflags=xcalledflags+1

end subroutine dipolesub


subroutine dipolesub_one(avector,inspfs,xdipole_expect,ydipole_expect,zdipole_expect)
  use xxxmod
  use parameters
  implicit none

  DATATYPE ::  inspfs(  spfsize, nspf ),  zdipole_expect, ydipole_expect, xdipole_expect, &
       avector(numconfig,numr),tempvector(numconfig,numr)
  DATATYPE :: dot,nullcomplex(1),dipoles(3), &
       tempspfs(spfsize,nspf), dipolemat(nspf,nspf)
  DATAECS :: rvector(numr)
  integer :: i
 
!! independent of R for now.  multiply by R for prolate  (R set to 1 for atom)
  call nucdipvalue(nullcomplex,dipoles)

  do i=1,numr
     tempvector(:,i)=avector(:,i)*bondpoints(i)
  enddo
  dipoles(:)=dipoles(:)*dot(avector,tempvector,numconfig*numr)

!! Z DIPOLE


  do i=1,nspf
     call mult_zdipole(inspfs(:,i),tempspfs(:,i))
  enddo
  call MYGEMM(CNORMCHAR,'N',nspf,nspf,spfsize,DATAONE, inspfs, spfsize, tempspfs, spfsize, DATAZERO, dipolemat, nspf)
  rvector(:)=bondpoints(:)
  call arbitraryconfig_mult(dipolemat,rvector,avector,tempvector,numr)
  zdipole_expect=dot(avector,tempvector,numconfig*numr)
  zdipole_expect=zdipole_expect + dipoles(3)


!! Y DIPOLE

  do i=1,nspf
     call mult_ydipole(inspfs(:,i),tempspfs(:,i))
  enddo
  call MYGEMM(CNORMCHAR,'N',nspf,nspf,spfsize,DATAONE, inspfs, spfsize, tempspfs, spfsize, DATAZERO, dipolemat, nspf)
  rvector(:)=bondpoints(:)
  call arbitraryconfig_mult(dipolemat,rvector,avector,tempvector,numr)
  ydipole_expect=dot(avector,tempvector,numconfig*numr)
  ydipole_expect=ydipole_expect + dipoles(2)


!! X DIPOLE

  do i=1,nspf
     call mult_xdipole(inspfs(:,i),tempspfs(:,i))
  enddo
  call MYGEMM(CNORMCHAR,'N',nspf,nspf,spfsize,DATAONE, inspfs, spfsize, tempspfs, spfsize, DATAZERO, dipolemat, nspf)
  rvector(:)=bondpoints(:)
  call arbitraryconfig_mult(dipolemat,rvector,avector,tempvector,numr)
  xdipole_expect=dot(avector,tempvector,numconfig*numr)
  xdipole_expect=xdipole_expect + dipoles(1)

end subroutine dipolesub_one


subroutine dipolesub_final()
  use parameters
  use dipolemod
  implicit none
  deallocate( zdipoleexpect, xdipoleexpect,  ydipoleexpect)

end subroutine dipolesub_final


