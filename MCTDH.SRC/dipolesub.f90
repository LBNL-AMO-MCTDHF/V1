

#include "Definitions.INC"

module dipolemod
  implicit none
  DATATYPE, allocatable :: xdipoleexpect(:),ydipoleexpect(:),zdipoleexpect(:),dipolenormsq(:)
  integer :: calledflag=0,xcalledflag=0
end module dipolemod

subroutine dipolesub_initial()
  use dipolemod
  use parameters
  implicit none
  allocate(xdipoleexpect(0:autosize), ydipoleexpect(0:autosize), zdipoleexpect(0:autosize),&
       dipolenormsq(0:autosize))
end subroutine

subroutine dipolesub()
  use dipolemod
  use parameters
  use configmod
  use xxxmod
  implicit none
  DATATYPE :: xx(mcscfnum),yy(mcscfnum),zz(mcscfnum),dd(mcscfnum),dot
  integer :: imc,sflag
  integer, save :: lastouttime=0
  real*8 :: thistime

  if (mod(xcalledflag,autosteps).eq.0) then

     if (conjgpropflag.ne.0) then
        if (mcscfnum.ne.2) then
           OFLWR "Whoot? conjgpropflag mcscfnum",mcscfnum; CFLST
        endif
        dipolenormsq(calledflag) = dot(yyy%cmfpsivec(astart(2),0),yyy%cmfpsivec(astart(1),0),tot_adim)
        if (par_consplit.ne.0) then
           call mympireduceone(dipolenormsq(calledflag))
        endif

        call dipolesub_one(www,yyy%cmfpsivec(astart(2),0),yyy%cmfpsivec(astart(1),0),yyy%cmfpsivec(spfstart,0), xdipoleexpect(calledflag),ydipoleexpect(calledflag),zdipoleexpect(calledflag))

     else

        do imc=1,mcscfnum
           dd(imc) = dot(yyy%cmfpsivec(astart(imc),0),yyy%cmfpsivec(astart(imc),0),tot_adim)
        if (par_consplit.ne.0) then
           call mympireduceone(dd(imc))
        endif

           call dipolesub_one(www,yyy%cmfpsivec(astart(imc),0),yyy%cmfpsivec(astart(imc),0),yyy%cmfpsivec(spfstart,0), xx(imc),yy(imc),zz(imc))
        enddo

        xdipoleexpect(calledflag)=0d0
        ydipoleexpect(calledflag)=0d0
        zdipoleexpect(calledflag)=0d0
        dipolenormsq(calledflag)=0d0

        do imc=1,mcscfnum

           dipolenormsq(calledflag) = dipolenormsq(calledflag) + dd(imc)

!! 101414 REAL-VALUED FOR HERM.
!! 1-2016 v1.17 should not be necessary with realflag in mult_zdipole(in,out,realflag) etc.
 
#ifndef CNORMFLAG
           xdipoleexpect(calledflag) = xdipoleexpect(calledflag) + real(xx(imc),8)
           ydipoleexpect(calledflag) = ydipoleexpect(calledflag) + real(yy(imc),8)
           zdipoleexpect(calledflag) = zdipoleexpect(calledflag) + real(zz(imc),8)
#else
           xdipoleexpect(calledflag) = xdipoleexpect(calledflag) + xx(imc)
           ydipoleexpect(calledflag) = ydipoleexpect(calledflag) + yy(imc)
           zdipoleexpect(calledflag) = zdipoleexpect(calledflag) + zz(imc)
#endif

        enddo

     endif  !! conjgpropflag complex Domcke

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

     if (conjgpropflag.ne.0) then
        OFLWR "   complex Domcke - off diagonal norm-squared ", dipolenormsq(calledflag)
     endif

     calledflag=calledflag+1

  endif
  xcalledflag=xcalledflag+1

end subroutine dipolesub


subroutine dipolesub_one(www,abra,aket,inspfs,xdipole_expect,ydipole_expect,zdipole_expect)
  use r_parameters
  use spfsize_parameters
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  DATATYPE, intent(in) :: inspfs(  spfsize, www%nspf ), abra(numr,www%firstconfig:www%lastconfig),&
       aket(numr,www%firstconfig:www%lastconfig)
  DATATYPE,intent(out) ::    zdipole_expect, ydipole_expect, xdipole_expect
  DATATYPE :: dot,nullcomplex(1),dipoles(3), dipolemat(www%nspf,www%nspf),csum
  DATATYPE,allocatable :: tempvector(:,:),tempspfs(:,:)
  DATAECS :: rvector(numr)
  integer :: i
 
  allocate(tempvector(numr,www%firstconfig:www%lastconfig), tempspfs(spfsize,www%nspf))

!! independent of R for now.  multiply by R for prolate  (R set to 1 for atom)
  call nucdipvalue(nullcomplex,dipoles)

  do i=1,numr
     tempvector(i,:)=aket(i,:)*bondpoints(i)
  enddo
  csum=dot(abra,tempvector,www%totadim)
  if (www%parconsplit.ne.0) then
     call mympireduceone(csum)
  endif
  dipoles(:)=dipoles(:)*csum

!! Z DIPOLE


  call mult_zdipole(www%nspf,inspfs(:,:),tempspfs(:,:),1)

  call MYGEMM(CNORMCHAR,'N',www%nspf,www%nspf,spfsize,DATAONE, inspfs, spfsize, tempspfs, spfsize, DATAZERO, dipolemat, www%nspf)
  if (parorbsplit.eq.3) then
     call mympireduce(dipolemat(:,:),www%nspf**2)
  endif

  rvector(:)=bondpoints(:)
  call arbitraryconfig_mult_singles(www,dipolemat,rvector,aket,tempvector,numr)
  zdipole_expect=dot(abra,tempvector,www%totadim)
  if (www%parconsplit.ne.0) then
     call mympireduceone(zdipole_expect)
  endif
  zdipole_expect=zdipole_expect + dipoles(3)


!! Y DIPOLE

  call mult_ydipole(www%nspf,inspfs(:,:),tempspfs(:,:),1)

  call MYGEMM(CNORMCHAR,'N',www%nspf,www%nspf,spfsize,DATAONE, inspfs, spfsize, tempspfs, spfsize, DATAZERO, dipolemat, www%nspf)
  if (parorbsplit.eq.3) then
     call mympireduce(dipolemat(:,:),www%nspf**2)
  endif

  rvector(:)=bondpoints(:)
  call arbitraryconfig_mult_singles(www,dipolemat,rvector,aket,tempvector,numr)
  ydipole_expect=dot(abra,tempvector,www%totadim)
  if (www%parconsplit.ne.0) then
     call mympireduceone(ydipole_expect)
  endif
  ydipole_expect=ydipole_expect + dipoles(2)


!! X DIPOLE

  call mult_xdipole(www%nspf,inspfs(:,:),tempspfs(:,:),1)

  call MYGEMM(CNORMCHAR,'N',www%nspf,www%nspf,spfsize,DATAONE, inspfs, spfsize, tempspfs, spfsize, DATAZERO, dipolemat, www%nspf)
  if (parorbsplit.eq.3) then
     call mympireduce(dipolemat(:,:),www%nspf**2)
  endif

  rvector(:)=bondpoints(:)
  call arbitraryconfig_mult_singles(www,dipolemat,rvector,aket,tempvector,numr)
  xdipole_expect=dot(abra,tempvector,www%totadim)
  if (www%parconsplit.ne.0) then
     call mympireduceone(xdipole_expect)
  endif
  xdipole_expect=xdipole_expect + dipoles(1)

  deallocate(tempvector,tempspfs)

end subroutine dipolesub_one


subroutine dipolesub_final()
  use dipolemod
  implicit none
  deallocate( zdipoleexpect, xdipoleexpect,  ydipoleexpect, dipolenormsq)

end subroutine dipolesub_final


