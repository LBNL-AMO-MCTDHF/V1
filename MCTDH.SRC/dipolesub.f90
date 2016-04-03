

!! ACTION 21 (emission/absorption) subroutine

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
  use mpisubmod
  use pulse_parameters !! conjgpropflag
  implicit none
  DATATYPE :: xx(mcscfnum),yy(mcscfnum),zz(mcscfnum),dd(mcscfnum),&
       axx(mcscfnum),ayy(mcscfnum),azz(mcscfnum),sxx(mcscfnum),syy(mcscfnum),&
       szz(mcscfnum),drivingoverlap(mcscfnum)
  integer :: imc,sflag
  integer, save :: lastouttime=0
  real*8 :: thistime
  xx=0;yy=0;zz=0;axx=0;ayy=0;azz=0;sxx=0;syy=0;szz=0;dd=0;drivingoverlap=0

  if (mod(xcalledflag,autosteps).eq.0) then

     if (conjgpropflag.ne.0) then
        if (mcscfnum.ne.2) then
           OFLWR "Whoot? conjgpropflag mcscfnum",mcscfnum; CFLST
        endif
        if (drivingflag.ne.0) then
           OFLWR "Driving not supported for conjprop yet"; CFLST
        endif
        dipolenormsq(calledflag)=0
        if (tot_adim.gt.0) then
           dipolenormsq(calledflag) = hermdot(yyy%cmfavec(:,2,0),yyy%cmfavec(:,1,0),tot_adim)
        endif
        if (par_consplit.ne.0) then
           call mympireduceone(dipolenormsq(calledflag))
        endif
        call dipolesub_one(www,bwwptr,yyy%cmfavec(:,2,0),yyy%cmfavec(:,1,0), yyy%cmfspfs(:,0),&
             xdipoleexpect(calledflag),ydipoleexpect(calledflag),zdipoleexpect(calledflag))

     else

        do imc=1,mcscfnum
           dd(imc)=0
           if (tot_adim.gt.0) then
              dd(imc) = hermdot(yyy%cmfavec(:,imc,0),yyy%cmfavec(:,imc,0),tot_adim)
           endif
           if (par_consplit.ne.0) then
              call mympireduceone(dd(imc))
           endif
           call dipolesub_one(www,bwwptr,yyy%cmfavec(:,imc,0),yyy%cmfavec(:,imc,0),&
                yyy%cmfspfs(:,0), xx(imc),yy(imc),zz(imc))

        enddo

        if (drivingflag.ne.0) then
#ifdef CNORMFLAG
           OFLWR "Error, driving with dipole not supported c-norm"; CFLST
#endif
           call dipolesub_driving(axx,ayy,azz,sxx,syy,szz,mcscfnum)
           xx(:)=xx(:)+axx(:)+CONJUGATE(axx(:))+sxx*drivingproportion**2
           yy(:)=yy(:)+ayy(:)+CONJUGATE(ayy(:))+syy*drivingproportion**2
           zz(:)=zz(:)+azz(:)+CONJUGATE(azz(:))+szz*drivingproportion**2

           call getdrivingoverlap(drivingoverlap,mcscfnum)   !! for time slot zero
           dd(:)=dd(:)+drivingproportion**2 + drivingoverlap(:) + &
                CONJUGATE(drivingoverlap(:))
        endif

        xdipoleexpect(calledflag)=0d0
        ydipoleexpect(calledflag)=0d0
        zdipoleexpect(calledflag)=0d0
        dipolenormsq(calledflag)=0d0

        do imc=1,mcscfnum

           dipolenormsq(calledflag) = dipolenormsq(calledflag) + dd(imc)

!! 101414 REAL-VALUED FOR HERM.
!! 1-2016 v1.17 should not be necessary with realflag in mult_zdipole(in,out,realflag) etc.
 
#ifndef CNORMFLAG
           xdipoleexpect(calledflag) = xdipoleexpect(calledflag) + real(xx(imc),8) / mcscfnum
           ydipoleexpect(calledflag) = ydipoleexpect(calledflag) + real(yy(imc),8) / mcscfnum
           zdipoleexpect(calledflag) = zdipoleexpect(calledflag) + real(zz(imc),8) / mcscfnum
#else
           xdipoleexpect(calledflag) = xdipoleexpect(calledflag) + xx(imc) / mcscfnum
           ydipoleexpect(calledflag) = ydipoleexpect(calledflag) + yy(imc) / mcscfnum
           zdipoleexpect(calledflag) = zdipoleexpect(calledflag) + zz(imc) / mcscfnum
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


subroutine checkorbsetrange(checknspf,flag)
  use parameters
  implicit none
  integer,intent(in) :: checknspf
  integer,intent(out) :: flag
  flag=0
  if (nspf.ne.checknspf) then
     flag=1
  endif
end subroutine checkorbsetrange


module dipbiomod
  use biorthotypemod
  implicit none
  type(biorthotype),target :: dipbiovar
end module dipbiomod


subroutine dipolesub_one(www,bioww,in_abra,&    !! ok unused bioww
     in_aket,inspfs,xdipole_expect,ydipole_expect,zdipole_expect)
  use r_parameters
  use spfsize_parameters
  use walkmod
  use fileptrmod
  use dotmod
  use dipbiomod
  use biorthomod
  use arbitrarymultmod
  use orbgathersubmod
  use mpisubmod
  implicit none
  type(walktype),intent(in) :: www,bioww
  DATATYPE, intent(in) :: inspfs(  spfsize, www%nspf ), &
       in_abra(numr,www%firstconfig:www%lastconfig),&
       in_aket(numr,www%firstconfig:www%lastconfig)
  DATATYPE,intent(out) ::    zdipole_expect, ydipole_expect, xdipole_expect
  DATATYPE,allocatable :: tempvector(:,:),tempspfs(:,:),abra(:,:),workspfs(:,:),&
       aket(:,:)
  DATATYPE :: nullcomplex(1),dipoles(3), dipolemat(www%nspf,www%nspf),csum
  DATAECS :: rvector(numr)
!!$  DATATYPE :: norm   !! datatype in case abra.ne.aket
  integer :: i,lowspf,highspf,numspf
#ifdef CNORMFLAG
  DATATYPE,target :: smo(www%nspf,www%nspf)
#endif

  lowspf=1; highspf=www%nspf
  if (parorbsplit.eq.1) then
     call checkorbsetrange(www%nspf,i)
     if (i.ne.0) then
        OFLWR "error exit, can't do dipolesub parorbsplit.eq.1 with",www%nspf,"orbitals"; CFLST
     endif
     call getOrbSetRange(lowspf,highspf)
  endif
  numspf=highspf-lowspf+1
 
  allocate(tempvector(numr,www%firstconfig:www%lastconfig+1), tempspfs(spfsize,lowspf:highspf+1),&
       abra(numr,www%firstconfig:www%lastconfig+1),workspfs(spfsize,www%nspf),&
       aket(numr,www%firstconfig:www%lastconfig+1))

  tempvector=0; tempspfs=0; workspfs=0; abra=0; aket=0

  if (www%lastconfig.ge.www%firstconfig) then
     abra(:,www%firstconfig:www%lastconfig)=in_abra(:,:)
     aket(:,www%firstconfig:www%lastconfig)=in_aket(:,:)
  endif

#ifndef CNORMFLAG
  workspfs(:,:)=inspfs(:,:)
#else
  call bioset(dipbiovar,smo,numr,bioww)
  dipbiovar%hermonly=.true.
  call biortho(inspfs,inspfs,workspfs,abra,dipbiovar)
  dipbiovar%hermonly=.false.
#endif


!!$  csum=dot(abra,aket,www%totadim)
!!$  if (www%parconsplit.ne.0) then
!!$     call mympireduceone(csum)
!!$  endif
!!$  norm=sqrt(csum)

!! independent of R for now.  multiply by R for prolate  (R set to 1 for atom)
  call nucdipvalue(nullcomplex,dipoles)

  do i=1,numr
     tempvector(i,:)=aket(i,:)*bondpoints(i)
  enddo
  csum=0d0
  if (www%totadim.gt.0) then
     csum=hermdot(abra,tempvector,www%totadim)
  endif
  if (www%parconsplit.ne.0) then
     call mympireduceone(csum)
  endif
  dipoles(:)=dipoles(:)*csum

!! Z DIPOLE

  dipolemat(:,:)=0d0
  if (numspf.gt.0) then
     call mult_zdipole(numspf,inspfs(:,lowspf:highspf),tempspfs(:,lowspf:highspf),1)
     call MYGEMM('C','N',www%nspf,numspf,spfsize,DATAONE, workspfs, spfsize, &
          tempspfs(:,lowspf:highspf), spfsize, DATAZERO, dipolemat(:,lowspf:highspf), www%nspf)
  endif
  if (parorbsplit.eq.1) then
     call mpiorbgather(dipolemat,www%nspf)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(dipolemat(:,:),www%nspf**2)
  endif

  rvector(:)=bondpoints(:)
  call arbitraryconfig_mult_singles(www,dipolemat,rvector,aket,tempvector,numr)
  zdipole_expect=0d0
  if (www%totadim.gt.0) then
     zdipole_expect=hermdot(abra,tempvector,www%totadim)
  endif
  if (www%parconsplit.ne.0) then
     call mympireduceone(zdipole_expect)
  endif
  zdipole_expect=zdipole_expect + dipoles(3)

!! Y DIPOLE

  dipolemat(:,:)=0d0
  if (numspf.gt.0) then
     call mult_ydipole(numspf,inspfs(:,lowspf:highspf),tempspfs(:,lowspf:highspf),1)

     call MYGEMM('C','N',www%nspf,numspf,spfsize,DATAONE, workspfs, spfsize, &
          tempspfs(:,lowspf:highspf), spfsize, DATAZERO, dipolemat(:,lowspf:highspf), www%nspf)
  endif
  if (parorbsplit.eq.1) then
     call mpiorbgather(dipolemat,www%nspf)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(dipolemat(:,:),www%nspf**2)
  endif

  rvector(:)=bondpoints(:)
  call arbitraryconfig_mult_singles(www,dipolemat,rvector,aket,tempvector,numr)
  ydipole_expect=0d0
  if (www%totadim.gt.0) then
     ydipole_expect=hermdot(abra,tempvector,www%totadim)
  endif
  if (www%parconsplit.ne.0) then
     call mympireduceone(ydipole_expect)
  endif
  ydipole_expect=ydipole_expect + dipoles(2)

!! X DIPOLE

  dipolemat(:,:)=0d0
  if (numspf.gt.0) then
     call mult_xdipole(numspf,inspfs(:,lowspf:highspf),tempspfs(:,lowspf:highspf),1)

     call MYGEMM('C','N',www%nspf,numspf,spfsize,DATAONE, workspfs, spfsize, &
          tempspfs(:,lowspf:highspf), spfsize, DATAZERO, dipolemat(:,lowspf:highspf), www%nspf)
  endif
  if (parorbsplit.eq.1) then
     call mpiorbgather(dipolemat,www%nspf)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(dipolemat(:,:),www%nspf**2)
  endif

  rvector(:)=bondpoints(:)
  call arbitraryconfig_mult_singles(www,dipolemat,rvector,aket,tempvector,numr)
  xdipole_expect=0d0
  if (www%totadim.gt.0) then
     xdipole_expect=hermdot(abra,tempvector,www%totadim)
  endif
  if (www%parconsplit.ne.0) then
     call mympireduceone(xdipole_expect)
  endif
  xdipole_expect=xdipole_expect + dipoles(1)

  deallocate(tempvector,tempspfs,abra,aket,workspfs)

!!$#ifdef CNORMFLAG
!!$  xdipole_expect=xdipole_expect*abs(norm)/norm
!!$  ydipole_expect=ydipole_expect*abs(norm)/norm
!!$  zdipole_expect=zdipole_expect*abs(norm)/norm
!!$#endif

end subroutine dipolesub_one


subroutine dipolesub_final()
  use dipolemod
  implicit none
  deallocate( zdipoleexpect, xdipoleexpect,  ydipoleexpect, dipolenormsq)

end subroutine dipolesub_final


