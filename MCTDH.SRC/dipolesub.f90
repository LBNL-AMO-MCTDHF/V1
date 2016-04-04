

!! ACTION 21 (emission/absorption) subroutine

#include "Definitions.INC"

module dipolemod
  implicit none
  DATATYPE, allocatable :: dipoleexpects(:,:,:),    dipolenormsq(:)
  integer :: calledflag=0,xcalledflag=0
end module dipolemod

subroutine dipolesub_initial()
  use dipolemod
  use parameters
  use pulse_parameters !! conjgpropflag
  implicit none

  if (conjgpropflag.eq.0) then
     allocate(dipoleexpects(0:autosize,3,1))   !! x,y,z
  else
     allocate(dipoleexpects(0:autosize,3,4))
  endif
  allocate(dipolenormsq(0:autosize))

end subroutine


subroutine dipolesub()
  use dipolemod
  use parameters
  use configmod
  use xxxmod
  use mpisubmod
  use pulse_parameters !! conjgpropflag
  use dipolecallsubmod
  implicit none

  DATATYPE :: myexpects(3), mcexpects(3,mcscfnum), dd(mcscfnum),&
       axx(mcscfnum),ayy(mcscfnum),azz(mcscfnum),sxx(mcscfnum),syy(mcscfnum),&
       szz(mcscfnum),drivingoverlap(mcscfnum)
  integer :: imc,sflag,getlen,ii
  integer, save :: lastouttime=0
  real*8 :: thistime
  character(len=2) :: tl(4) = (/ "BA", "AB", "AA", "AB" /)

  myexpects=0;mcexpects=0;axx=0;ayy=0;azz=0;sxx=0;syy=0;szz=0;dd=0;drivingoverlap=0

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

        call dipolesub_one(www,bwwptr,yyy%cmfavec(:,2,0),yyy%cmfavec(:,1,0), yyy%cmfspfs(:,0), myexpects(:))
        dipoleexpects(calledflag,:,1)=myexpects(:)
        call dipolesub_one(www,bwwptr,yyy%cmfavec(:,1,0),yyy%cmfavec(:,2,0), yyy%cmfspfs(:,0), myexpects(:))
        dipoleexpects(calledflag,:,2)=myexpects(:)
        call dipolesub_one(www,bwwptr,yyy%cmfavec(:,1,0),yyy%cmfavec(:,1,0), yyy%cmfspfs(:,0), myexpects(:))
        dipoleexpects(calledflag,:,3)=myexpects(:)
        call dipolesub_one(www,bwwptr,yyy%cmfavec(:,2,0),yyy%cmfavec(:,2,0), yyy%cmfspfs(:,0), myexpects(:))
        dipoleexpects(calledflag,:,4)=myexpects(:)

     else  !! conjgpropflag complex Domcke

        do imc=1,mcscfnum
           dd(imc)=0
           if (tot_adim.gt.0) then
              dd(imc) = hermdot(yyy%cmfavec(:,imc,0),yyy%cmfavec(:,imc,0),tot_adim)
           endif
           if (par_consplit.ne.0) then
              call mympireduceone(dd(imc))
           endif
           call dipolesub_one(www,bwwptr,yyy%cmfavec(:,imc,0),yyy%cmfavec(:,imc,0),yyy%cmfspfs(:,0),mcexpects(:,imc))

        enddo

        if (drivingflag.ne.0) then
#ifdef CNORMFLAG
           OFLWR "Error, driving with dipole not supported c-norm"; CFLST
#endif
           call dipolesub_driving(axx,ayy,azz,sxx,syy,szz,mcscfnum)
           mcexpects(1,:)=mcexpects(1,:)+axx(:)+CONJUGATE(axx(:))+sxx*drivingproportion**2
           mcexpects(2,:)=mcexpects(2,:)+ayy(:)+CONJUGATE(ayy(:))+syy*drivingproportion**2
           mcexpects(3,:)=mcexpects(3,:)+azz(:)+CONJUGATE(azz(:))+szz*drivingproportion**2

           call getdrivingoverlap(drivingoverlap,mcscfnum)   !! for time slot zero
           dd(:)=dd(:)+drivingproportion**2 + drivingoverlap(:) + &
                CONJUGATE(drivingoverlap(:))
        endif

        dipoleexpects(calledflag,:,1)=0d0
        dipolenormsq(calledflag)=0d0

        do imc=1,mcscfnum

           dipolenormsq(calledflag) = dipolenormsq(calledflag) + dd(imc)

!! 101414 REAL-VALUED FOR HERM.
!! 1-2016 v1.17 should not be necessary with realflag in mult_zdipole(in,out,realflag) etc.
 
#ifndef CNORMFLAG
           dipoleexpects(calledflag,:,1) = dipoleexpects(calledflag,:,1) + real(mcexpects(:,imc),8) / mcscfnum
#else
           dipoleexpects(calledflag,:,1) = dipoleexpects(calledflag,:,1) + mcexpects(:,imc) / mcscfnum
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

        if (conjgpropflag.eq.0) then
           call dipolecall(calledflag, dipoleexpects(:,:,1),   xdipfile(1:getlen(xdipfile)-1),&
                xdftfile(1:getlen(xdftfile)-1),            xoworkfile(1:getlen(xoworkfile)-1),&
                xtworkfile(1:getlen(xtworkfile)-1),    xophotonfile(1:getlen(xophotonfile)-1),&
                1,sflag)
           call dipolecall(calledflag, dipoleexpects(:,:,1),   ydipfile(1:getlen(ydipfile)-1),&
                ydftfile(1:getlen(ydftfile)-1),            yoworkfile(1:getlen(xoworkfile)-1),&
                ytworkfile(1:getlen(xtworkfile)-1),    yophotonfile(1:getlen(xophotonfile)-1),&
                2,sflag)
           call dipolecall(calledflag, dipoleexpects(:,:,1),   zdipfile(1:getlen(zdipfile)-1),&
                zdftfile(1:getlen(zdftfile)-1),            zoworkfile(1:getlen(xoworkfile)-1),&
                ztworkfile(1:getlen(xtworkfile)-1),    zophotonfile(1:getlen(xophotonfile)-1),&
                3,sflag)
           if (act21circ.ne.0) then
              call dipolecall(calledflag, dipoleexpects(:,:,1),  xydipfile(1:getlen(xydipfile)-1),&
                   xydftfile(1:getlen(xydftfile)-1),          xyoworkfile(1:getlen(xoworkfile)-1),&
                   xytworkfile(1:getlen(xtworkfile)-1),   xyophotonfile(1:getlen(xophotonfile)-1),&
                   4,sflag)
              call dipolecall(calledflag, dipoleexpects(:,:,1),  xzdipfile(1:getlen(xzdipfile)-1),&
                   xzdftfile(1:getlen(xzdftfile)-1),          xzoworkfile(1:getlen(xoworkfile)-1),&
                   xztworkfile(1:getlen(xtworkfile)-1),   xzophotonfile(1:getlen(xophotonfile)-1),&
                   5,sflag)
              call dipolecall(calledflag, dipoleexpects(:,:,1),  yxdipfile(1:getlen(yxdipfile)-1),&
                   yxdftfile(1:getlen(yxdftfile)-1),          yxoworkfile(1:getlen(xoworkfile)-1),&
                   yxtworkfile(1:getlen(xtworkfile)-1),   yxophotonfile(1:getlen(xophotonfile)-1),&
                   6,sflag)
              call dipolecall(calledflag, dipoleexpects(:,:,1),  yzdipfile(1:getlen(yzdipfile)-1),&
                   yzdftfile(1:getlen(yzdftfile)-1),          yzoworkfile(1:getlen(xoworkfile)-1),&
                   yztworkfile(1:getlen(xtworkfile)-1),   yzophotonfile(1:getlen(xophotonfile)-1),&
                   7,sflag)
              call dipolecall(calledflag, dipoleexpects(:,:,1),  zxdipfile(1:getlen(zxdipfile)-1),&
                   zxdftfile(1:getlen(zxdftfile)-1),          zxoworkfile(1:getlen(xoworkfile)-1),&
                   zxtworkfile(1:getlen(xtworkfile)-1),   zxophotonfile(1:getlen(xophotonfile)-1),&
                   8,sflag)
              call dipolecall(calledflag, dipoleexpects(:,:,1),  zydipfile(1:getlen(zydipfile)-1),&
                   zydftfile(1:getlen(zydftfile)-1),          zyoworkfile(1:getlen(xoworkfile)-1),&
                   zytworkfile(1:getlen(xtworkfile)-1),   zyophotonfile(1:getlen(xophotonfile)-1),&
                   9,sflag)
           endif
        else
           do ii=1,4
              call dipolecall(calledflag, dipoleexpects(:,:,1),      xdipfile(1:getlen(xdipfile)-1)//tl(ii),&
                   xdftfile(1:getlen(xdftfile)-1)//tl(ii),         xoworkfile(1:getlen(xoworkfile)-1)//tl(ii),&
                   xtworkfile(1:getlen(xtworkfile)-1)//tl(ii),   xophotonfile(1:getlen(xophotonfile)-1)//tl(ii),&
                   1,sflag)
              call dipolecall(calledflag, dipoleexpects(:,:,1),      ydipfile(1:getlen(ydipfile)-1)//tl(ii),&
                   ydftfile(1:getlen(ydftfile)-1)//tl(ii),         yoworkfile(1:getlen(xoworkfile)-1)//tl(ii),&
                   ytworkfile(1:getlen(xtworkfile)-1)//tl(ii),   yophotonfile(1:getlen(xophotonfile)-1)//tl(ii),&
                   2,sflag)
              call dipolecall(calledflag, dipoleexpects(:,:,1),     zdipfile(1:getlen(zdipfile)-1)//tl(ii),&
                   zdftfile(1:getlen(zdftfile)-1)//tl(ii),        zoworkfile(1:getlen(xoworkfile)-1)//tl(ii),&
                   ztworkfile(1:getlen(xtworkfile)-1)//tl(ii),  zophotonfile(1:getlen(xophotonfile)-1)//tl(ii),&
                   3,sflag)
              if (act21circ.ne.0) then
                 call dipolecall(calledflag, dipoleexpects(:,:,1),      xydipfile(1:getlen(xydipfile)-1)//tl(ii),&
                      xydftfile(1:getlen(xydftfile)-1)//tl(ii),       xyoworkfile(1:getlen(xoworkfile)-1)//tl(ii),&
                      xytworkfile(1:getlen(xtworkfile)-1)//tl(ii),  xyophotonfile(1:getlen(xophotonfile)-1)//tl(ii),&
                      4,sflag)
                 call dipolecall(calledflag, dipoleexpects(:,:,1),      xzdipfile(1:getlen(xzdipfile)-1)//tl(ii),&
                      xzdftfile(1:getlen(xzdftfile)-1)//tl(ii),       xzoworkfile(1:getlen(xoworkfile)-1)//tl(ii),&
                      xztworkfile(1:getlen(xtworkfile)-1)//tl(ii),  xzophotonfile(1:getlen(xophotonfile)-1)//tl(ii),&
                      5,sflag)
                 call dipolecall(calledflag, dipoleexpects(:,:,1),      yxdipfile(1:getlen(yxdipfile)-1)//tl(ii),&
                      yxdftfile(1:getlen(yxdftfile)-1)//tl(ii),       yxoworkfile(1:getlen(xoworkfile)-1)//tl(ii),&
                      yxtworkfile(1:getlen(xtworkfile)-1)//tl(ii),  yxophotonfile(1:getlen(xophotonfile)-1)//tl(ii),&
                      6,sflag)
                 call dipolecall(calledflag, dipoleexpects(:,:,1),      yzdipfile(1:getlen(yzdipfile)-1)//tl(ii),&
                      yzdftfile(1:getlen(yzdftfile)-1)//tl(ii),       yzoworkfile(1:getlen(xoworkfile)-1)//tl(ii),&
                      yztworkfile(1:getlen(xtworkfile)-1)//tl(ii),  yzophotonfile(1:getlen(xophotonfile)-1)//tl(ii),&
                      7,sflag)
                 call dipolecall(calledflag, dipoleexpects(:,:,1),      zxdipfile(1:getlen(zxdipfile)-1)//tl(ii),&
                      zxdftfile(1:getlen(zxdftfile)-1)//tl(ii),       zxoworkfile(1:getlen(xoworkfile)-1)//tl(ii),&
                      zxtworkfile(1:getlen(xtworkfile)-1)//tl(ii),  zxophotonfile(1:getlen(xophotonfile)-1)//tl(ii),&
                      8,sflag)
                 call dipolecall(calledflag, dipoleexpects(:,:,1),      zydipfile(1:getlen(zydipfile)-1)//tl(ii),&
                      zydftfile(1:getlen(zydftfile)-1)//tl(ii),       zyoworkfile(1:getlen(xoworkfile)-1)//tl(ii),&
                      zytworkfile(1:getlen(xtworkfile)-1)//tl(ii),  zyophotonfile(1:getlen(xophotonfile)-1)//tl(ii),&
                      9,sflag)
              endif
           enddo
        endif    !! conjgpropflag
     endif       !! calledflag (dipmodtime)

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
     in_aket,inspfs,dipole_expects)
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
  DATATYPE,intent(out) :: dipole_expects(3)
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
  dipole_expects(3)=0d0
  if (www%totadim.gt.0) then
     dipole_expects(3)=hermdot(abra,tempvector,www%totadim)
  endif
  if (www%parconsplit.ne.0) then
     call mympireduceone(dipole_expects(3))
  endif
  dipole_expects(3)=dipole_expects(3) + dipoles(3)

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
  dipole_expects(2)=0d0
  if (www%totadim.gt.0) then
     dipole_expects(2)=hermdot(abra,tempvector,www%totadim)
  endif
  if (www%parconsplit.ne.0) then
     call mympireduceone(dipole_expects(2))
  endif
  dipole_expects(2)=dipole_expects(2) + dipoles(2)

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
  dipole_expects(1)=0d0
  if (www%totadim.gt.0) then
     dipole_expects(1)=hermdot(abra,tempvector,www%totadim)
  endif
  if (www%parconsplit.ne.0) then
     call mympireduceone(dipole_expects(1))
  endif
  dipole_expects(1)=dipole_expects(1) + dipoles(1)

  deallocate(tempvector,tempspfs,abra,aket,workspfs)

!!$#ifdef CNORMFLAG
!!$  dipole_expects(1)=dipole_expects(1)*abs(norm)/norm
!!$  dipole_expects(2)=dipole_expects(2)*abs(norm)/norm
!!$  dipole_expects(3)=dipole_expects(3)*abs(norm)/norm
!!$#endif

end subroutine dipolesub_one


subroutine dipolesub_final()
  use dipolemod
  implicit none
  deallocate( dipoleexpects, dipolenormsq)

end subroutine dipolesub_final


