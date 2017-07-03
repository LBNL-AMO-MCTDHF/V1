
!! ALL MODULES

!! FOR finalstats.dat (after relaxation) or psistats.dat (action 25 for time-dependent calculation)

#include "Definitions.INC"

module orbmatsubmod
contains

subroutine get_orbmats( myspfs,  howmany,  ugmat,   &
     xdipmat,ydipmat,zdipmat,   xrefmat,yrefmat,zrefmat,conjgmat)
  use parameters
  use orbgathersubmod
  use mpisubmod
  implicit none
  integer, intent(in) :: howmany
  DATATYPE,intent(in) :: myspfs(spfsize,howmany)
  DATATYPE, intent(out) :: ugmat(howmany,howmany), xdipmat(howmany,howmany), &
       ydipmat(howmany,howmany), zdipmat(howmany,howmany), &
       xrefmat(howmany,howmany), yrefmat(howmany,howmany), zrefmat(howmany,howmany),&
       conjgmat(howmany,howmany)
  DATATYPE,allocatable ::  tempspfs(:,:),tempspfs2(:,:)
  integer :: i,lowspf,highspf,numspf

  allocate(tempspfs(spfsize,howmany),tempspfs2(spfsize,howmany))
  tempspfs(:,:)=0; tempspfs2(:,:)=0

  lowspf=1; highspf=howmany
  if (parorbsplit.eq.1) then
     if (howmany.ne.nspf) then
        OFLWR "In get_orbmats can't do parorbsplit.eq.1 with howmany.ne.nspf, error exit",howmany,nspf
        CFLST
     endif
     call getOrbSetRange(lowspf,highspf)
  endif
  numspf=highspf-lowspf+1

!! CONJUGATION
  do i=lowspf,highspf
     call op_conjg(myspfs(:,i),tempspfs(:,i))
  enddo

  if (numspf.gt.0) then
     call MYGEMM(CNORMCHAR,'N',howmany,numspf,spfsize,DATAONE, myspfs, spfsize, &
          tempspfs(:,lowspf:highspf), spfsize, DATAZERO, conjgmat(:,lowspf:highspf), howmany)
  endif
  if (parorbsplit.eq.1) then
     call mpiorbgather(conjgmat,nspf)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(conjgmat(:,:),howmany**2)
  endif

!! M & U/G

  do i=lowspf,highspf
     call op_reflectx(myspfs(:,i),tempspfs(:,i))
     call op_reflecty(tempspfs(:,i),tempspfs2(:,i))
     call op_reflectz(tempspfs2(:,i),tempspfs(:,i))
  enddo

  if (numspf.gt.0) then
     call MYGEMM(CNORMCHAR,'N',howmany,numspf,spfsize,DATAONE, myspfs, spfsize, &
          tempspfs(:,lowspf:highspf), spfsize, DATAZERO, ugmat(:,lowspf:highspf), howmany)
  endif
  if (parorbsplit.eq.1) then
     call mpiorbgather(ugmat,nspf)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(ugmat(:,:),howmany**2)
  endif

!! Z DIPOLE
  if (numspf.gt.0) then
     call mult_zdipole(numspf,myspfs(:,lowspf:highspf),tempspfs(:,lowspf:highspf),0)

     call MYGEMM(CNORMCHAR,'N',howmany,numspf,spfsize,DATAONE, myspfs, spfsize, &
          tempspfs(:,lowspf:highspf), spfsize, DATAZERO, zdipmat(:,lowspf:highspf), howmany)
  endif
  if (parorbsplit.eq.1) then
     call mpiorbgather(zdipmat,nspf)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(zdipmat(:,:),howmany**2)
  endif

!! Y DIPOLE
  if (numspf.gt.0) then
     call mult_ydipole(numspf,myspfs(:,lowspf:highspf),tempspfs(:,lowspf:highspf),0)

     call MYGEMM(CNORMCHAR,'N',howmany,numspf,spfsize,DATAONE, myspfs, spfsize, &
          tempspfs(:,lowspf:highspf), spfsize, DATAZERO, ydipmat(:,lowspf:highspf), howmany)
  endif
  if (parorbsplit.eq.1) then
     call mpiorbgather(ydipmat,nspf)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(ydipmat(:,:),howmany**2)
  endif

!! X DIPOLE
  if (numspf.gt.0) then
     call mult_xdipole(numspf,myspfs(:,lowspf:highspf),tempspfs(:,lowspf:highspf),0)

     call MYGEMM(CNORMCHAR,'N',howmany,numspf,spfsize,DATAONE, myspfs, spfsize, &
          tempspfs(:,lowspf:highspf), spfsize, DATAZERO, xdipmat(:,lowspf:highspf), howmany)
  endif
  if (parorbsplit.eq.1) then
     call mpiorbgather(xdipmat,nspf)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(xdipmat(:,:),howmany**2)
  endif



!! REFLECTIONS

  do i=lowspf,highspf
     call op_reflectz(myspfs(:,i),tempspfs(:,i))
  enddo
  if (numspf.gt.0) then
     call MYGEMM(CNORMCHAR,'N',howmany,numspf,spfsize,DATAONE, myspfs, spfsize, &
          tempspfs(:,lowspf:highspf), spfsize, DATAZERO, zrefmat(:,lowspf:highspf), howmany)
  endif
  if (parorbsplit.eq.1) then
     call mpiorbgather(zrefmat,nspf)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(zrefmat(:,:),howmany**2)
  endif

  do i=lowspf,highspf
     call op_reflecty(myspfs(:,i),tempspfs(:,i))
  enddo
  if (numspf.gt.0) then
     call MYGEMM(CNORMCHAR,'N',howmany,numspf,spfsize,DATAONE, myspfs, spfsize, &
          tempspfs(:,lowspf:highspf), spfsize, DATAZERO, yrefmat(:,lowspf:highspf), howmany)
  endif
  if (parorbsplit.eq.1) then
     call mpiorbgather(yrefmat,nspf)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(yrefmat(:,:),howmany**2)
  endif

  do i=lowspf,highspf
     call op_reflectx(myspfs(:,i),tempspfs(:,i))
  enddo
  if (numspf.gt.0) then
     call MYGEMM(CNORMCHAR,'N',howmany,numspf,spfsize,DATAONE, myspfs, spfsize, &
          tempspfs(:,lowspf:highspf), spfsize, DATAZERO, xrefmat(:,lowspf:highspf), howmany)
  endif
  if (parorbsplit.eq.1) then
     call mpiorbgather(xrefmat,nspf)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(xrefmat(:,:),howmany**2)
  endif

  deallocate(tempspfs,tempspfs2)

end subroutine get_orbmats

end module orbmatsubmod


module psistatsubmod
  use biorthotypemod
  implicit none
  integer,private :: psibioalloc=0, psibionumvec=(-1)
  type(biorthotype),target,allocatable,private :: ugbiovar(:),&
       xrefbiovar(:),yrefbiovar(:),zrefbiovar(:),conjgbiovar(:)

contains

subroutine psistats( thistime )
  use parameters
  use mpimod
  use configmod
  use xxxmod
  implicit none

  real*8,intent(in) :: thistime
  integer, save :: calledflag=0
  integer :: i,myiostat
  DATATYPE :: mexpect(mcscfnum),ugexpect(mcscfnum),m2expect(mcscfnum),&
       xreflect(mcscfnum),yreflect(mcscfnum),zreflect(mcscfnum),xdipole(mcscfnum),&
       ydipole(mcscfnum),zdipole(mcscfnum),conjugation(mcscfnum)

  if (calledflag==0) then
     if (myrank.eq.1) then
        open(662, file=psistatsfile, status="unknown",iostat=myiostat)
        call checkiostat(myiostat,"opening "//psistatsfile)
#ifdef REALGO        
        write(662,'(A16,200A14)',iostat=myiostat) &
             " Time        ", " M           ", " M2          ", " UG          ", &
             " XDipole     ", " YDipole     ", " ZDipole     ", " XReflect    ", &
             " YReflect    ", " ZReflect    ", " Conjugation ", " Abs(conjg)  "
#else
        write(662,'(A16,A14,200A28)',iostat=myiostat) &
             " Time        ", " M           ", " M2          ", " UG          ", &
             " XDipole     ", " YDipole     ", " ZDipole     ", " XReflect    ", &
             " YReflect    ", " ZReflect    ", " Conjugation ", " Abs(conjg)  "
#endif
        call checkiostat(myiostat,"writing "//psistatsfile)
        close(662)
     endif
  endif

  calledflag=1


  call get_psistats(www,bioww,yyy%cmfspfs(:,0),mcscfnum,yyy%cmfavec(:,:,0),&
       mexpect,m2expect,ugexpect,   xdipole,ydipole,zdipole,   &
       xreflect,yreflect,zreflect,conjugation)

  if (myrank.eq.1) then
     open(662, file=psistatsfile, status="old", position="append",iostat=myiostat)
     call checkiostat(myiostat,"opening "//psistatsfile)
     do i=1,mcscfnum
        write(662,'(F13.5,I3,2000F14.8)',iostat=myiostat) thistime, i, &
             mexpect(i),m2expect(i),ugexpect(i),  xdipole(i),ydipole(i),zdipole(i),  &
             xreflect(i),yreflect(i),zreflect(i), conjugation(i), abs(conjugation(i))+DATAZERO
     enddo
     call checkiostat(myiostat,"writing "//psistatsfile)
     close(662)
  endif

contains

subroutine get_psistats( wwin, bbin, myspfs, numvec, in_inavectors, mexpect,m2expect,ugexpect,&
     xdipole,ydipole,zdipole,   xreflect,yreflect,zreflect,conjgexpect)
  use parameters
  use walkmod
  use arbitrarymultmod
  use autocorrelate_one_mod
  use mpisubmod
  use orbmatsubmod
  implicit none
  type(walktype),intent(in) :: wwin, bbin
  integer, intent(in) :: numvec
  DATATYPE,intent(in) :: myspfs(spfsize,wwin%nspf), &
       in_inavectors(numr,wwin%firstconfig:wwin%lastconfig,numvec)       
  DATATYPE, intent(out) :: mexpect(numvec),ugexpect(numvec),m2expect(numvec),&
       xreflect(numvec),yreflect(numvec),zreflect(numvec),&
       xdipole(numvec),ydipole(numvec),zdipole(numvec),conjgexpect(numvec)
  DATATYPE,allocatable :: ugmat(:,:), xdipmat(:,:), &
       ydipmat(:,:), zdipmat(:,:), inavectors(:,:,:), &
       xrefmat(:,:), yrefmat(:,:), zrefmat(:,:),     &
       tempvector(:,:), tempspfs(:,:),tempspfs2(:,:),conjgmat(:,:)
  DATATYPE :: normsq(numvec),nullcomplex,dipoles(3),nucdipexpect(numvec,3),csum
  DATAECS :: rvector(numr)
  integer :: i,imc

  if (wwin%lastconfig.ge.wwin%firstconfig) then
     allocate(inavectors(numr,wwin%firstconfig:wwin%lastconfig,numvec))
     inavectors(:,:,:)=in_inavectors(:,:,:)
  else
     allocate(inavectors(numr,1,numvec))
     inavectors=0
  endif

  if (psibioalloc.ne.0.and.psibionumvec.ne.numvec) then
     deallocate(ugbiovar,xrefbiovar,yrefbiovar,zrefbiovar,conjgbiovar)
     psibioalloc=0
  endif
  if (psibioalloc.eq.0) then
     allocate(ugbiovar(numvec),xrefbiovar(numvec), yrefbiovar(numvec),zrefbiovar(numvec),&
          conjgbiovar(numvec))
  endif
  psibioalloc=1

  allocate(ugmat(wwin%nspf,wwin%nspf), xdipmat(wwin%nspf,wwin%nspf),&
       ydipmat(wwin%nspf,wwin%nspf), zdipmat(wwin%nspf,wwin%nspf), &
       xrefmat(wwin%nspf,wwin%nspf), yrefmat(wwin%nspf,wwin%nspf), &
       zrefmat(wwin%nspf,wwin%nspf),tempvector(numr,wwin%firstconfig:wwin%lastconfig), &
       tempspfs(spfsize,wwin%nspf),tempspfs2(spfsize,wwin%nspf),&
       conjgmat(wwin%nspf,wwin%nspf))
  
  ugmat=0; xdipmat=0; ydipmat=0; zdipmat=0; xrefmat=0; yrefmat=0; zrefmat=0
  call get_orbmats( myspfs,  wwin%nspf,  ugmat,   xdipmat,ydipmat,zdipmat,  &
       xrefmat,yrefmat,zrefmat, conjgmat)

  normsq(:)=0
  if (wwin%totadim.gt.0) then
     do imc=1,numvec
        normsq(imc)=dot(inavectors(:,:,imc),inavectors(:,:,imc),wwin%totadim)
     enddo
  endif
  if (wwin%parconsplit.ne.0) then
     call mympireduce(normsq,numvec)
  endif

!! CONJUGATION

  do i=1,nspf
     call op_conjg(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,numvec
     if (wwin%totadim.gt.0) then
        tempvector(:,:)=ALLCON(inavectors(:,:,imc))
     endif
     call autocorrelate_one(wwin,bbin,inavectors(:,:,imc),myspfs,&
          tempspfs,tempvector,conjgexpect(imc),numr,conjgbiovar(imc))
     conjgexpect(imc)=conjgexpect(imc)   /normsq(imc)
  enddo
  
!! M & U/G

  if (spfrestrictflag.eq.1) then
     mexpect(:)=0; ugexpect(:)=0; m2expect(:)=0
     if (wwin%totadim.gt.0) then
        do imc=1,numvec      
           do i=1,numr
              tempvector(i,:)=inavectors(i,:,imc) * &
                   wwin%configmvals(wwin%firstconfig:wwin%lastconfig)
           enddo
           mexpect(imc)=dot(inavectors(:,:,imc),tempvector(:,:),wwin%totadim)   /normsq(imc)
        
           do i=1,numr
              tempvector(i,:)=inavectors(i,:,imc) * &
                   wwin%configugvals(wwin%firstconfig:wwin%lastconfig)
           enddo
           ugexpect(imc)=dot(inavectors(:,:,imc),tempvector(:,:),wwin%totadim)   /normsq(imc)
        
           do i=1,numr
              tempvector(i,:)=inavectors(i,:,imc) * &
                   wwin%configmvals(wwin%firstconfig:wwin%lastconfig)**2
           enddo
           m2expect(imc)=dot(inavectors(:,:,imc),tempvector(:,:),wwin%totadim)   /normsq(imc)
        enddo
     endif
     if (wwin%parconsplit.ne.0) then
        call mympireduce(mexpect,numvec)
        call mympireduce(ugexpect,numvec)
        call mympireduce(m2expect,numvec)
     endif
  else
     mexpect=-99d0; m2expect=-99d0; !! could program these

     do i=1,nspf
        call op_reflectx(myspfs(:,i),tempspfs(:,i))
        call op_reflecty(tempspfs(:,i),tempspfs2(:,i))
        call op_reflectz(tempspfs2(:,i),tempspfs(:,i))
     enddo
     do imc=1,numvec
        call autocorrelate_one(wwin,bbin,inavectors(:,:,imc),myspfs,&
             tempspfs,inavectors(:,:,imc),ugexpect(imc),numr,ugbiovar(imc))
        ugexpect(imc)=ugexpect(imc)   /normsq(imc)
     enddo

  endif

!! like dipolesub

!! independent of R for now.  multiply by R for prolate  (R set to 1 for atom)
  call nucdipvalue(nullcomplex,dipoles)

  do imc=1,mcscfnum
     csum=0
     if (wwin%totadim.gt.0) then
        do i=1,numr
           tempvector(i,:)=inavectors(i,:,imc)*bondpoints(i)
        enddo
        csum=dot(inavectors(:,:,imc),tempvector,wwin%totadim)
     endif
     if (wwin%parconsplit.ne.0) then
        call mympireduceone(csum)
     endif
     nucdipexpect(imc,:)=dipoles(:)*csum
  enddo

  rvector(:)=bondpoints(:)

  zdipole(:)=0; ydipole(:)=0; xdipole(:)=0

!! Z DIPOLE

  do imc=1,numvec
     call arbitraryconfig_mult_singles(wwin,zdipmat,rvector,&
          inavectors(:,:,imc),tempvector,numr)
     if (wwin%totadim.gt.0) then
        zdipole(imc)=dot(inavectors(:,:,imc),tempvector,wwin%totadim)
     endif
  enddo
  if (wwin%parconsplit.ne.0) then
     call mympireduce(zdipole,numvec)
  endif
  zdipole(:)=(zdipole(:)+nucdipexpect(:,3))/normsq(:)

!! Y DIPOLE

  do imc=1,numvec
     call arbitraryconfig_mult_singles(wwin,ydipmat,rvector,&
          inavectors(:,:,imc),tempvector,numr)
     if (wwin%totadim.gt.0) then
        ydipole(imc)=dot(inavectors(:,:,imc),tempvector,wwin%totadim)
     endif
  enddo
  if (wwin%parconsplit.ne.0) then
     call mympireduce(ydipole,numvec)
  endif
  ydipole(:)=(ydipole(:)+nucdipexpect(:,2))/normsq(:)

!! X DIPOLE

  do imc=1,numvec
     call arbitraryconfig_mult_singles(wwin,xdipmat,rvector,&
          inavectors(:,:,imc),tempvector,numr)
     if (wwin%totadim.gt.0) then
        xdipole(imc)=dot(inavectors(:,:,imc),tempvector,wwin%totadim)
     endif
  enddo
  if (wwin%parconsplit.ne.0) then
     call mympireduce(xdipole,numvec)
  endif
  xdipole(:)=(xdipole(:)+nucdipexpect(:,1))/normsq(:)

!! REFLECTIONS

  do i=1,nspf
     call op_reflectz(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,numvec
     call autocorrelate_one(wwin,bbin,inavectors(:,:,imc),myspfs,tempspfs,&
          inavectors(:,:,imc),zreflect(imc),numr,zrefbiovar(imc))
     zreflect(imc)=zreflect(imc)   /normsq(imc)
  enddo

  do i=1,nspf
     call op_reflecty(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,numvec
     call autocorrelate_one(wwin,bbin,inavectors(:,:,imc),myspfs,tempspfs,&
          inavectors(:,:,imc),yreflect(imc),numr,yrefbiovar(imc))
     yreflect(imc)=yreflect(imc)   /normsq(imc)
  enddo

  do i=1,nspf
     call op_reflectx(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,numvec
     call autocorrelate_one(wwin,bbin,inavectors(:,:,imc),myspfs,tempspfs,&
          inavectors(:,:,imc),xreflect(imc),numr,xrefbiovar(imc))
     xreflect(imc)=xreflect(imc)   /normsq(imc)
  enddo

  deallocate(ugmat, xdipmat, ydipmat, zdipmat,     xrefmat, yrefmat, zrefmat, &
       tempvector,tempspfs,tempspfs2,inavectors, conjgmat)
  
end subroutine get_psistats

end subroutine psistats

end module psistatsubmod


module finalstatsubmod
contains

subroutine finalstats0(myspfs,in_inavectors,wwin,bbin )
  use parameters
  use mpimod
  use walkmod
  use arbitrarymultmod
  use biorthotypemod
  use autocorrelate_one_mod
  use mpisubmod
  use getoccmod
  use orbmatsubmod
  implicit none
  type(walktype),intent(in) :: wwin,bbin
  DATATYPE,intent(in) :: myspfs(spfsize,wwin%nspf), &
       in_inavectors(numr,wwin%firstconfig:wwin%lastconfig,mcscfnum)       
  type(biorthotype),target :: ugbiovar(mcscfnum,mcscfnum),xrefbiovar(mcscfnum,mcscfnum),&
       yrefbiovar(mcscfnum,mcscfnum),zrefbiovar(mcscfnum,mcscfnum),&
       conjgbiovar(mcscfnum,mcscfnum)
  DATATYPE :: mmatel(mcscfnum,mcscfnum),ugmatel(mcscfnum,mcscfnum),m2matel(mcscfnum,mcscfnum),&
       xrefmatel(mcscfnum,mcscfnum),yrefmatel(mcscfnum,mcscfnum),zrefmatel(mcscfnum,mcscfnum),&
       xdipmatel(mcscfnum,mcscfnum),ydipmatel(mcscfnum,mcscfnum),zdipmatel(mcscfnum,mcscfnum),&
       ovlmatel(mcscfnum,mcscfnum), nucdipmatel(mcscfnum,mcscfnum,3),&
       nullcomplex,dipoles(3), conjgmatel(mcscfnum,mcscfnum)
  DATATYPE,allocatable :: ugmat(:,:), xdipmat(:,:), ydipmat(:,:), zdipmat(:,:), inavectors(:,:,:), &
       xrefmat(:,:), yrefmat(:,:), zrefmat(:,:), tempvector(:,:), tempspfs(:,:),tempspfs2(:,:),&
       conjgmat(:,:)
  CNORMTYPE :: occupations(wwin%nspf,mcscfnum)
  DATAECS :: rvector(numr)
  integer :: i,j,imc,jmc,myiostat

  call mpibarrier()
  OFLWR "   ...GO finalstats."; CFL
  call mpibarrier()

  if (wwin%lastconfig.ge.wwin%firstconfig) then
     allocate(inavectors(numr,wwin%firstconfig:wwin%lastconfig,mcscfnum))
     inavectors(:,:,:)=in_inavectors(:,:,:)
  else
     allocate(inavectors(numr,1,mcscfnum))
     inavectors=0
  endif

  allocate(ugmat(wwin%nspf,wwin%nspf), xdipmat(wwin%nspf,wwin%nspf), &
       ydipmat(wwin%nspf,wwin%nspf), zdipmat(wwin%nspf,wwin%nspf), &
       xrefmat(wwin%nspf,wwin%nspf), yrefmat(wwin%nspf,wwin%nspf), &
       zrefmat(wwin%nspf,wwin%nspf),tempvector(numr,wwin%firstconfig:wwin%lastconfig), &
       tempspfs(spfsize,wwin%nspf),tempspfs2(spfsize,wwin%nspf),&
       conjgmat(wwin%nspf,wwin%nspf))

  if (wwin%totadim.gt.0) then
     tempvector(:,:)=0d0;
  endif
  tempspfs(:,:)=0d0; tempspfs2(:,:)=0d0

  ugmat=0; xdipmat=0; ydipmat=0; zdipmat=0; xrefmat=0; yrefmat=0; zrefmat=0

  call get_orbmats( myspfs,  wwin%nspf,  ugmat,   xdipmat,ydipmat,zdipmat,   &
       xrefmat,yrefmat,zrefmat,conjgmat)

  do imc=1,mcscfnum
     call getoccupations(wwin,inavectors(:,:,imc),numr,occupations(:,imc))
  enddo

!! CONJUGATION
  call mpibarrier()

  do i=1,nspf
     call op_conjg(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,mcscfnum
     if (wwin%totadim.gt.0) then
        tempvector(:,:)=ALLCON(inavectors(:,:,imc))
     endif
     do jmc=1,mcscfnum
        call autocorrelate_one(wwin,bbin,inavectors(:,:,jmc),myspfs,tempspfs,&
             tempvector,conjgmatel(jmc,imc),numr,conjgbiovar(jmc,imc))
     enddo
  enddo

!! M & U/G

  call mpibarrier()

  ovlmatel(:,:)=0
  if (wwin%totadim.gt.0) then
     do imc=1,mcscfnum
        do jmc=1,mcscfnum
           ovlmatel(jmc,imc)=dot(inavectors(:,:,jmc),inavectors(:,:,imc),wwin%totadim)
        enddo
     enddo
  endif
  if (wwin%parconsplit.ne.0) then
     call mympireduce(ovlmatel,mcscfnum**2)
  endif

  mmatel(:,:)=0; ugmatel(:,:)=0; m2matel(:,:)=0
  if (spfrestrictflag.eq.1) then
     if (wwin%totadim.gt.0) then
        do imc=1,mcscfnum
           do i=1,numr
              tempvector(i,:)=inavectors(i,:,imc) * wwin%configmvals(wwin%firstconfig:wwin%lastconfig)
           enddo
           do jmc=1,mcscfnum
              mmatel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector(:,:),wwin%totadim)
           enddo
           do i=1,numr
              tempvector(i,:)=inavectors(i,:,imc) * wwin%configugvals(wwin%firstconfig:wwin%lastconfig)
           enddo
           do jmc=1,mcscfnum
              ugmatel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector(:,:),wwin%totadim)
           enddo
           do i=1,numr
              tempvector(i,:)=inavectors(i,:,imc) * wwin%configmvals(wwin%firstconfig:wwin%lastconfig)**2
           enddo
           do jmc=1,mcscfnum
              m2matel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector(:,:),wwin%totadim)
           enddo
        enddo
     endif
     if (wwin%parconsplit.ne.0) then
        call mympireduce(mmatel,mcscfnum**2)
        call mympireduce(ugmatel,mcscfnum**2)
        call mympireduce(m2matel,mcscfnum**2)
     endif

  else

     mmatel=-99d0; m2matel=-99d0; !! could program these polyatomic

     do i=1,nspf
        call op_reflectx(myspfs(:,i),tempspfs(:,i))
        call op_reflecty(tempspfs(:,i),tempspfs2(:,i))
        call op_reflectz(tempspfs2(:,i),tempspfs(:,i))
     enddo
     do imc=1,mcscfnum
        do jmc=1,mcscfnum
           call autocorrelate_one(wwin,bbin,inavectors(:,:,jmc),myspfs,tempspfs,&
                inavectors(:,:,imc),ugmatel(jmc,imc),numr,ugbiovar(jmc,imc))
        enddo
     enddo
  endif

  call mpibarrier()

!! independent of R for now.  multiply by R for prolate  (R set to 1 for atom)
  call nucdipvalue(nullcomplex,dipoles)

  nucdipmatel(:,:,:)=0
  if (wwin%totadim.gt.0) then
     do imc=1,mcscfnum
        do i=1,numr
           tempvector(i,:)=inavectors(i,:,imc)*bondpoints(i)
        enddo
        do jmc=1,mcscfnum
           nucdipmatel(jmc,imc,:)=dipoles(:)*dot(inavectors(:,:,jmc),tempvector,wwin%totadim)
        enddo
     enddo
  endif
  if (wwin%parconsplit.ne.0) then
     call mympireduce(nucdipmatel,3*mcscfnum**2)
  endif

  rvector(:)=bondpoints(:)

  zdipmatel(:,:)=0; xdipmatel(:,:)=0; ydipmatel(:,:)=0

!! Z DIPOLE

  do imc=1,mcscfnum
     call arbitraryconfig_mult_singles(wwin,zdipmat,rvector,inavectors(:,:,imc),tempvector,numr)
     if (wwin%totadim.gt.0) then
        do jmc=1,mcscfnum
           zdipmatel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector,wwin%totadim)
        enddo
     endif
  enddo
  if (wwin%parconsplit.ne.0) then
     call mympireduce(zdipmatel(:,:),mcscfnum**2)
  endif
  zdipmatel(:,:)=zdipmatel(:,:)+nucdipmatel(:,:,3)

!! Y DIPOLE

  do imc=1,mcscfnum
     call arbitraryconfig_mult_singles(wwin,ydipmat,rvector,&
          inavectors(:,:,imc),tempvector,numr)
     if (wwin%totadim.gt.0) then
        do jmc=1,mcscfnum
           ydipmatel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector,wwin%totadim)
        enddo
     endif
  enddo
  if (wwin%parconsplit.ne.0) then
     call mympireduce(ydipmatel(:,:),mcscfnum**2)
  endif
  ydipmatel(:,:)=ydipmatel(:,:)+nucdipmatel(:,:,2)

!! X DIPOLE

  do imc=1,mcscfnum
     call arbitraryconfig_mult_singles(wwin,xdipmat,rvector,&
          inavectors(:,:,imc),tempvector,numr)
     if (wwin%totadim.gt.0) then
        do jmc=1,mcscfnum
           xdipmatel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector,wwin%totadim)
        enddo
     endif
  enddo
  if (wwin%parconsplit.ne.0) then
     call mympireduce(xdipmatel(:,:),mcscfnum**2)
  endif
  xdipmatel(:,:)=xdipmatel(:,:)+nucdipmatel(:,:,1)

!! REFLECTIONS

  call mpibarrier()

  do i=1,nspf
     call op_reflectz(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,mcscfnum
     do jmc=1,mcscfnum
        call autocorrelate_one(wwin,bbin,inavectors(:,:,jmc),myspfs,tempspfs,&
             inavectors(:,:,imc),zrefmatel(jmc,imc),numr,zrefbiovar(jmc,imc))
     enddo
  enddo

  do i=1,nspf
     call op_reflecty(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,mcscfnum
     do jmc=1,mcscfnum
        call autocorrelate_one(wwin,bbin,inavectors(:,:,jmc),myspfs,tempspfs,&
             inavectors(:,:,imc),yrefmatel(jmc,imc),numr,yrefbiovar(jmc,imc))
     enddo
  enddo

  do i=1,nspf
     call op_reflectx(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,mcscfnum
     do jmc=1,mcscfnum
        call autocorrelate_one(wwin,bbin,inavectors(:,:,jmc),myspfs,tempspfs,&
             inavectors(:,:,imc),xrefmatel(jmc,imc),numr,xrefbiovar(jmc,imc))
     enddo
  enddo

  if (myrank.eq.1) then
     open(662, file=finalstatsfile, status="unknown",iostat=myiostat)
     call checkiostat(myiostat,"opening "//finalstatsfile)
     write(662,*,iostat=myiostat);     
     call checkiostat(myiostat,"writing "//finalstatsfile)
     write(662,*)

     write(662,*) "--- EXPECTATION VALUES PSI VECTORS ---"

     write(662,*);     write(662,*)
     write(662,*) "  conjugation expectation values psi vectors, value and absolute value "
     do i=1,mcscfnum
        write(662,'(I10,2000F14.8)') i, conjgmatel(i,i),abs(conjgmatel(i,i))
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  inversion expectation values psi vectors "
     do i=1,mcscfnum
        write(662,'(I10,2000F14.8)') i, ugmatel(i,i)
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  M expectation values psi vectors "
     do i=1,mcscfnum
        write(662,'(I10,2000F14.8)') i, mmatel(i,i)
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  M^2 expectation values psi vectors "
     do i=1,mcscfnum
        write(662,'(I10,2000F14.8)') i, m2matel(i,i)
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  X, Y, Z dipole expectation values psi vectors"
     do i=1,mcscfnum
        write(662,'(I10,2000F14.8)') i, xdipmatel(i,i), ydipmatel(i,i), zdipmatel(i,i)
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  X, Y, Z reflection expectation values psi vectors"
     do i=1,mcscfnum
        write(662,'(I10,2000F14.8)') i, xrefmatel(i,i), yrefmatel(i,i), zrefmatel(i,i)
     enddo

     write(662,*);     write(662,*)
     write(662,*) "--- EXPECTATION VALUES ORBITALS ---"

     write(662,*);     write(662,*)
     write(662,*) "  conjugation expectation values orbitals, value and absolute value"
     do i=1,nspf
        write(662,'(I10,2000F14.8)') i, conjgmat(i,i), abs(conjgmat(i,i))
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  inversion expectation values orbitals"
     do i=1,nspf
        write(662,'(I10,2000F14.8)') i, ugmat(i,i)
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  X,Y,Z dipole expectation values orbitals"
     do i=1,nspf
        write(662,'(I10,2000F14.8)') i,xdipmat(i,i),ydipmat(i,i),zdipmat(i,i)
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  X,Y,Z reflection expectation values orbitals "
     do i=1,nspf
        write(662,'(I10, 2000F14.8)') i,xrefmat(i,i),yrefmat(i,i),zrefmat(i,i)
     enddo


     write(662,*);     write(662,*)
     write(662,*) "--- ORBITAL OCCUPATIONS BY STATE ---"
     write(662,'(A15,1000I15)') "Orbital State",(i,i=1,mcscfnum)
     do i=1,nspf
        write(662,'(I15,1000F15.10)') i, (occupations(i,j),j=1,mcscfnum)
     enddo


     write(662,*);     write(662,*)
     write(662,*) "---------------- MATRIX ELEMENTS ------------------"


     write(662,*);     write(662,*)
     write(662,*) "  conjugation matrix elements, absolute value "
     write(662,*) "--------------------------------"
     write(662,*) "     ORBITALS    "
     write(662,*) "--------------------------------"
     do i=1,nspf
        write(662,'(2000F14.8)') (abs(conjgmat(i,j)),j=1,nspf)
     enddo
     write(662,*) "--------------------------------"
     write(662,*) "     PSI VECTORS    "
     write(662,*) "--------------------------------"
     do i=1,mcscfnum
        write(662,'(2000F14.8)') (abs(conjgmatel(i,j)),j=1,mcscfnum)
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  inversion matrix elements "
     write(662,*) "--------------------------------"
     write(662,*) "     ORBITALS    "
     write(662,*) "--------------------------------"
     do i=1,nspf
        write(662,'(2000F14.8)') (ugmat(i,j),j=1,nspf)
     enddo
     write(662,*) "--------------------------------"
     write(662,*) "     PSI VECTORS    "
     write(662,*) "--------------------------------"
     do i=1,mcscfnum
        write(662,'(2000F14.8)') (ugmatel(i,j),j=1,mcscfnum)
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  M matrix elements "
     write(662,*) "--------------------------------"
     write(662,*) "     PSI VECTORS    "
     write(662,*) "--------------------------------"
     do i=1,mcscfnum
        write(662,'(2000F14.8)') (mmatel(i,j),j=1,mcscfnum)
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  M^2 matrix elements "
     write(662,*) "--------------------------------"
     write(662,*) "     PSI VECTORS    "
     write(662,*) "--------------------------------"
     do i=1,mcscfnum
        write(662,'(2000F14.8)') (m2matel(i,j),j=1,mcscfnum)
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  Z dipole matrix elements "
     write(662,*) "--------------------------------"
     write(662,*) "     ORBITALS    "
     write(662,*) "--------------------------------"
     do i=1,nspf
        write(662,'(2000F14.8)') (zdipmat(i,j),j=1,nspf)
     enddo
     write(662,*) "--------------------------------"
     write(662,*) "     PSI VECTORS    "
     write(662,*) "--------------------------------"
     do i=1,mcscfnum
        write(662,'(2000F14.8)') (zdipmatel(i,j),j=1,mcscfnum)
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  Y dipole matrix elements "
     write(662,*) "--------------------------------"
     write(662,*) "     ORBITALS    "
     write(662,*) "--------------------------------"
     do i=1,nspf
        write(662,'(2000F14.8)') (ydipmat(i,j),j=1,nspf)
     enddo
     write(662,*) "--------------------------------"
     write(662,*) "     PSI VECTORS    "
     write(662,*) "--------------------------------"
     do i=1,mcscfnum
        write(662,'(2000F14.8)') (ydipmatel(i,j),j=1,mcscfnum)
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  X dipole matrix elements "
     write(662,*) "--------------------------------"
     write(662,*) "     ORBITALS    "
     write(662,*) "--------------------------------"
     do i=1,nspf
        write(662,'(2000F14.8)') (xdipmat(i,j),j=1,nspf)
     enddo
     write(662,*) "--------------------------------"
     write(662,*) "     PSI VECTORS    "
     write(662,*) "--------------------------------"
     do i=1,mcscfnum
        write(662,'(2000F14.8)') (xdipmatel(i,j),j=1,mcscfnum)
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  Z reflection matrix elements "
     write(662,*) "--------------------------------"
     write(662,*) "     ORBITALS    "
     write(662,*) "--------------------------------"
     do i=1,nspf
        write(662,'(2000F14.8)') (zrefmat(i,j),j=1,nspf)
     enddo
     write(662,*) "--------------------------------"
     write(662,*) "     PSI VECTORS    "
     write(662,*) "--------------------------------"
     do i=1,mcscfnum
        write(662,'(2000F14.8)') (zrefmatel(i,j),j=1,mcscfnum)
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  Y reflection matrix elements "
     write(662,*) "--------------------------------"
     write(662,*) "     ORBITALS    "
     write(662,*) "--------------------------------"
     do i=1,nspf
        write(662,'(2000F14.8)') (yrefmat(i,j),j=1,nspf)
     enddo
     write(662,*) "--------------------------------"
     write(662,*) "     PSI VECTORS    "
     write(662,*) "--------------------------------"
     do i=1,mcscfnum
        write(662,'(2000F14.8)') (yrefmatel(i,j),j=1,mcscfnum)
     enddo

     write(662,*);     write(662,*)
     write(662,*) "  X reflection matrix elements "
     write(662,*) "--------------------------------"
     write(662,*) "     ORBITALS    "
     write(662,*) "--------------------------------"
     do i=1,nspf
        write(662,'(2000F14.8)') (xrefmat(i,j),j=1,nspf)
     enddo
     write(662,*) "--------------------------------"
     write(662,*) "     PSI VECTORS    "
     write(662,*) "--------------------------------"
     do i=1,mcscfnum
        write(662,'(2000F14.8)') (xrefmatel(i,j),j=1,mcscfnum)
     enddo

     close(662)

  endif

  deallocate(ugmat, xdipmat, ydipmat, zdipmat,     xrefmat, yrefmat, zrefmat, &
       tempvector,tempspfs,tempspfs2,inavectors, conjgmat)

  call mpibarrier()

end subroutine finalstats0


subroutine finalstats( )
  use parameters
  use configmod
  use xxxmod
  implicit none

  call finalstats0(yyy%cmfspfs(:,0),  yyy%cmfavec(:,:,0),www,bioww)

end subroutine finalstats


end module



