
#include "Definitions.INC"


subroutine get_orbmats( myspfs,  numspf,  ugmat,   xdipmat,ydipmat,zdipmat,   xrefmat,yrefmat,zrefmat)
  use parameters
  implicit none

  integer, intent(in) :: numspf
  DATATYPE,intent(in) :: myspfs(spfsize,numspf)
  DATATYPE, intent(out) :: ugmat(numspf,numspf), xdipmat(numspf,numspf), ydipmat(numspf,numspf), zdipmat(numspf,numspf), &
       xrefmat(numspf,numspf), yrefmat(numspf,numspf), zrefmat(numspf,numspf)
  DATATYPE ::  tempspfs(spfsize,numspf),tempspfs2(spfsize,numspf)
  integer :: i

!! M & U/G

  do i=1,numspf
     call op_reflectx(myspfs(:,i),tempspfs(:,i))
     call op_reflecty(tempspfs(:,i),tempspfs2(:,i))
     call op_reflectz(tempspfs2(:,i),tempspfs(:,i))
  enddo
  call MYGEMM(CNORMCHAR,'N',numspf,numspf,spfsize,DATAONE, myspfs, spfsize, tempspfs, spfsize, DATAZERO, ugmat, numspf)
  if (parorbsplit.eq.3) then
     call mympireduce(ugmat(:,:),numspf**2)
  endif


!! Z DIPOLE

  do i=1,numspf
     call mult_zdipole(myspfs(:,i),tempspfs(:,i))
  enddo
  call MYGEMM(CNORMCHAR,'N',numspf,numspf,spfsize,DATAONE, myspfs, spfsize, tempspfs, spfsize, DATAZERO, zdipmat, numspf)
  if (parorbsplit.eq.3) then
     call mympireduce(zdipmat(:,:),numspf**2)
  endif

!! Y DIPOLE

  do i=1,numspf
     call mult_ydipole(myspfs(:,i),tempspfs(:,i))
  enddo
  call MYGEMM(CNORMCHAR,'N',numspf,numspf,spfsize,DATAONE, myspfs, spfsize, tempspfs, spfsize, DATAZERO, ydipmat, numspf)
  if (parorbsplit.eq.3) then
     call mympireduce(ydipmat(:,:),numspf**2)
  endif

!! X DIPOLE

  do i=1,numspf
     call mult_xdipole(myspfs(:,i),tempspfs(:,i))
  enddo
  call MYGEMM(CNORMCHAR,'N',numspf,numspf,spfsize,DATAONE, myspfs, spfsize, tempspfs, spfsize, DATAZERO, xdipmat, numspf)
  if (parorbsplit.eq.3) then
     call mympireduce(xdipmat(:,:),numspf**2)
  endif

!! REFLECTIONS

  do i=1,numspf
     call op_reflectz(myspfs(:,i),tempspfs(:,i))
  enddo
  call MYGEMM(CNORMCHAR,'N',numspf,numspf,spfsize,DATAONE, myspfs, spfsize, tempspfs, spfsize, DATAZERO, zrefmat, numspf)
  if (parorbsplit.eq.3) then
     call mympireduce(zrefmat(:,:),numspf**2)
  endif

  do i=1,numspf
     call op_reflecty(myspfs(:,i),tempspfs(:,i))
  enddo
  call MYGEMM(CNORMCHAR,'N',numspf,numspf,spfsize,DATAONE, myspfs, spfsize, tempspfs, spfsize, DATAZERO, yrefmat, numspf)
  if (parorbsplit.eq.3) then
     call mympireduce(yrefmat(:,:),numspf**2)
  endif

  do i=1,numspf
     call op_reflectx(myspfs(:,i),tempspfs(:,i))
  enddo
  call MYGEMM(CNORMCHAR,'N',numspf,numspf,spfsize,DATAONE, myspfs, spfsize, tempspfs, spfsize, DATAZERO, xrefmat, numspf)
  if (parorbsplit.eq.3) then
     call mympireduce(xrefmat(:,:),numspf**2)
  endif

  
end subroutine get_orbmats




subroutine psistats( thistime )
  use parameters
  use mpimod
  use xxxmod
  implicit none

  real*8,intent(in) :: thistime
  integer, save :: calledflag=0
  integer :: i
  DATATYPE :: mexpect(mcscfnum),ugexpect(mcscfnum),m2expect(mcscfnum),&
       xreflect(mcscfnum),yreflect(mcscfnum),zreflect(mcscfnum),xdipole(mcscfnum),&
       ydipole(mcscfnum),zdipole(mcscfnum)

  if (calledflag==0) then
     if (myrank.eq.1) then
        open(662, file=psistatsfile, status="unknown")
#ifdef REALGO        
        write(662,'(A16,200A14)') &
             " Time     ", " M        ", " M2       ", " UG       ", " XDipole  ", " YDipole  ", " ZDipole  ", " XReflect ", " YReflect ", " ZReflect "
#else
        write(662,'(A16,A14,200A28)') &
             " Time     ", " M        ", " M2       ", " UG       ", " XDipole  ", " YDipole  ", " ZDipole  ", " XReflect ", " YReflect ", " ZReflect "
#endif
        close(662)
     endif
  endif

  calledflag=1

  call get_psistats(yyy%cmfpsivec(spfstart,0),mcscfnum,yyy%cmfpsivec(astart(1),0),&
       mexpect,m2expect,ugexpect,   xdipole,ydipole,zdipole,   xreflect,yreflect,zreflect)
       
  if (myrank.eq.1) then
     open(662, file=psistatsfile, status="old", position="append")
     do i=1,mcscfnum
        write(662,'(F13.5,I3,2000F14.8)') thistime, i, mexpect(i),m2expect(i),ugexpect(i),  xdipole(i),ydipole(i),zdipole(i),  xreflect(i),yreflect(i),zreflect(i)
     enddo
     close(662)
  endif

end subroutine psistats



subroutine get_psistats( myspfs, numvec, inavectors, mexpect,m2expect,ugexpect,   xdipole,ydipole,zdipole,   xreflect,yreflect,zreflect)
  use parameters
  implicit none

  integer, intent(in) :: numvec
  DATATYPE,intent(in) :: myspfs(spfsize,nspf), inavectors(numconfig,numr,numvec)       
  DATATYPE, intent(out) :: mexpect(numvec),ugexpect(numvec),m2expect(numvec),&
       xreflect(numvec),yreflect(numvec),zreflect(numvec),xdipole(numvec),ydipole(numvec),zdipole(numvec)
  DATATYPE :: ugmat(nspf,nspf), xdipmat(nspf,nspf), ydipmat(nspf,nspf), zdipmat(nspf,nspf), &
       xrefmat(nspf,nspf), yrefmat(nspf,nspf), zrefmat(nspf,nspf),&
       dot,normsq(numvec),nullcomplex,dipoles(3),tempvector(numconfig,numr),&
       tempspfs(spfsize,nspf),tempspfs2(spfsize,nspf), nucdipexpect(numvec,3)
  DATAECS :: rvector(numr)
  integer :: i,imc

  call get_orbmats( myspfs,  nspf,  ugmat,   xdipmat,ydipmat,zdipmat,   xrefmat,yrefmat,zrefmat)

!! M & U/G

  do imc=1,numvec
     normsq(imc)=dot(inavectors(:,:,imc),inavectors(:,:,imc),totadim)
  enddo

  if (spfrestrictflag.eq.1) then
     do imc=1,numvec      
        do i=1,numr
           tempvector(:,i)=inavectors(:,i,imc) * configmvals(:)
        enddo
        mexpect(imc)=dot(inavectors(:,:,imc),tempvector(:,:),totadim)   /normsq(imc)
        
        do i=1,numr
           tempvector(:,i)=inavectors(:,i,imc) * configugvals(:)
        enddo
        ugexpect(imc)=dot(inavectors(:,:,imc),tempvector(:,:),totadim)   /normsq(imc)
        
        do i=1,numr
           tempvector(:,i)=inavectors(:,i,imc) * configmvals(:)**2
        enddo
        m2expect(imc)=dot(inavectors(:,:,imc),tempvector(:,:),totadim)   /normsq(imc)
     enddo
  else
     mexpect=-99d0; m2expect=-99d0; !! could program these

     do i=1,nspf
        call op_reflectx(myspfs(:,i),tempspfs(:,i))
        call op_reflecty(tempspfs(:,i),tempspfs2(:,i))
        call op_reflectz(tempspfs2(:,i),tempspfs(:,i))
     enddo
     do imc=1,numvec
        call autocorrelate_one(inavectors(:,:,imc),myspfs,tempspfs,inavectors(:,:,imc),ugexpect(imc),numr)
        ugexpect(imc)=ugexpect(imc)   /normsq(imc)
     enddo

  endif

!! like dipolesub

!! independent of R for now.  multiply by R for prolate  (R set to 1 for atom)
  call nucdipvalue(nullcomplex,dipoles)

  do imc=1,mcscfnum
     do i=1,numr
        tempvector(:,i)=inavectors(:,i,imc)*bondpoints(i)
     enddo
     nucdipexpect(imc,:)=dipoles(:)*dot(inavectors(:,:,imc),tempvector,numconfig*numr)
  enddo

  rvector(:)=bondpoints(:)

!! Z DIPOLE

  do imc=1,numvec
     call arbitraryconfig_mult(zdipmat,rvector,inavectors(:,:,imc),tempvector,numr)
     zdipole(imc)=(dot(inavectors(:,:,imc),tempvector,numconfig*numr) + nucdipexpect(imc,3))   /normsq(imc)
  enddo

!! Y DIPOLE

  do imc=1,numvec
     call arbitraryconfig_mult(ydipmat,rvector,inavectors(:,:,imc),tempvector,numr)
     ydipole(imc)=(dot(inavectors(:,:,imc),tempvector,numconfig*numr) + nucdipexpect(imc,2))   /normsq(imc)
  enddo

!! X DIPOLE

  do imc=1,numvec
     call arbitraryconfig_mult(xdipmat,rvector,inavectors(:,:,imc),tempvector,numr)
     xdipole(imc)=(dot(inavectors(:,:,imc),tempvector,numconfig*numr) + nucdipexpect(imc,1))   /normsq(imc)
  enddo

!! REFLECTIONS

  do i=1,nspf
     call op_reflectz(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,numvec
     call autocorrelate_one(inavectors(:,:,imc),myspfs,tempspfs,inavectors(:,:,imc),zreflect(imc),numr)
     zreflect(imc)=zreflect(imc)   /normsq(imc)
  enddo

  do i=1,nspf
     call op_reflecty(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,numvec
     call autocorrelate_one(inavectors(:,:,imc),myspfs,tempspfs,inavectors(:,:,imc),yreflect(imc),numr)
     yreflect(imc)=yreflect(imc)   /normsq(imc)
  enddo

  do i=1,nspf
     call op_reflectx(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,numvec
     call autocorrelate_one(inavectors(:,:,imc),myspfs,tempspfs,inavectors(:,:,imc),xreflect(imc),numr)
     xreflect(imc)=xreflect(imc)   /normsq(imc)
  enddo
  
end subroutine get_psistats



subroutine finalstats( )
  use parameters
  use mpimod
  use xxxmod
  implicit none

  DATATYPE :: myspfs(spfsize,nspf), inavectors(numconfig,numr,mcscfnum)       
  DATATYPE :: mmatel(mcscfnum,mcscfnum),ugmatel(mcscfnum,mcscfnum),m2matel(mcscfnum,mcscfnum),&
       xrefmatel(mcscfnum,mcscfnum),yrefmatel(mcscfnum,mcscfnum),zrefmatel(mcscfnum,mcscfnum),&
       xdipmatel(mcscfnum,mcscfnum),ydipmatel(mcscfnum,mcscfnum),zdipmatel(mcscfnum,mcscfnum),&
       ovlmatel(mcscfnum,mcscfnum), nucdipmatel(mcscfnum,mcscfnum,3)
  DATATYPE :: ugmat(nspf,nspf), xdipmat(nspf,nspf), ydipmat(nspf,nspf), zdipmat(nspf,nspf), &
       xrefmat(nspf,nspf), yrefmat(nspf,nspf), zrefmat(nspf,nspf),&
       dot,nullcomplex,dipoles(3),tempvector(numconfig,numr),&
       tempspfs(spfsize,nspf),tempspfs2(spfsize,nspf)
  DATAECS :: rvector(numr)
  integer :: i,j,imc,jmc

  myspfs=RESHAPE(yyy%cmfpsivec(spfstart:spfend,0),(/spfsize,nspf/))

  inavectors(:,:,:)=RESHAPE(yyy%cmfpsivec(astart(1):aend(mcscfnum),0),(/numconfig,numr,mcscfnum/))

  call get_orbmats( myspfs,  nspf,  ugmat,   xdipmat,ydipmat,zdipmat,   xrefmat,yrefmat,zrefmat)

!! M & U/G

  do imc=1,mcscfnum
     do jmc=1,mcscfnum
        ovlmatel(jmc,imc)=dot(inavectors(:,:,jmc),inavectors(:,:,imc),totadim)
     enddo
  enddo

  if (spfrestrictflag.eq.1) then
     do imc=1,mcscfnum
        do i=1,numr
           tempvector(:,i)=inavectors(:,i,imc) * configmvals(:)
        enddo
        do jmc=1,mcscfnum
           mmatel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector(:,:),totadim)
        enddo
        do i=1,numr
           tempvector(:,i)=inavectors(:,i,imc) * configugvals(:)
        enddo
        do jmc=1,mcscfnum
           ugmatel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector(:,:),totadim)
        enddo
        do i=1,numr
           tempvector(:,i)=inavectors(:,i,imc) * configmvals(:)**2
        enddo
        do jmc=1,mcscfnum
           m2matel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector(:,:),totadim)
        enddo
     enddo
  else

     mmatel=-99d0; m2matel=-99d0; !! could program these polyatomic

     do i=1,nspf
        call op_reflectx(myspfs(:,i),tempspfs(:,i))
        call op_reflecty(tempspfs(:,i),tempspfs2(:,i))
        call op_reflectz(tempspfs2(:,i),tempspfs(:,i))
     enddo
     do imc=1,mcscfnum
        do jmc=1,mcscfnum
           call autocorrelate_one(inavectors(:,:,jmc),myspfs,tempspfs,inavectors(:,:,imc),ugmatel(jmc,imc),numr)
        enddo
     enddo
  endif


!! independent of R for now.  multiply by R for prolate  (R set to 1 for atom)
  call nucdipvalue(nullcomplex,dipoles)

  do imc=1,mcscfnum
     do i=1,numr
        tempvector(:,i)=inavectors(:,i,imc)*bondpoints(i)
     enddo
     do jmc=1,mcscfnum
        nucdipmatel(jmc,imc,:)=dipoles(:)*dot(inavectors(:,:,jmc),tempvector,numconfig*numr)
     enddo
  enddo

  rvector(:)=bondpoints(:)

!! Z DIPOLE

  do imc=1,mcscfnum
     call arbitraryconfig_mult(zdipmat,rvector,inavectors(:,:,imc),tempvector,numr)
     do jmc=1,mcscfnum
        zdipmatel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector,numconfig*numr) + nucdipmatel(jmc,imc,3)
     enddo
  enddo

!! Y DIPOLE

  do imc=1,mcscfnum
     call arbitraryconfig_mult(ydipmat,rvector,inavectors(:,:,imc),tempvector,numr)
     do jmc=1,mcscfnum
        ydipmatel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector,numconfig*numr) + nucdipmatel(jmc,imc,2)
     enddo
  enddo

!! X DIPOLE

  do imc=1,mcscfnum
     call arbitraryconfig_mult(xdipmat,rvector,inavectors(:,:,imc),tempvector,numr)
     do jmc=1,mcscfnum
        xdipmatel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector,numconfig*numr) + nucdipmatel(jmc,imc,1)
     enddo
  enddo

!! REFLECTIONS

  do i=1,nspf
     call op_reflectz(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,mcscfnum
     do jmc=1,mcscfnum
        call autocorrelate_one(inavectors(:,:,jmc),myspfs,tempspfs,inavectors(:,:,imc),zrefmatel(jmc,imc),numr)
     enddo
  enddo

  do i=1,nspf
     call op_reflecty(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,mcscfnum
     do jmc=1,mcscfnum
        call autocorrelate_one(inavectors(:,:,jmc),myspfs,tempspfs,inavectors(:,:,imc),yrefmatel(jmc,imc),numr)
     enddo
  enddo

  do i=1,nspf
     call op_reflectx(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,mcscfnum
     do jmc=1,mcscfnum
        call autocorrelate_one(inavectors(:,:,jmc),myspfs,tempspfs,inavectors(:,:,imc),xrefmatel(jmc,imc),numr)
     enddo
  enddo



  if (myrank.eq.1) then
     open(662, file=finalstatsfile, status="unknown")

     write(662,*)
     write(662,*)
     write(662,*) "  M matrix elements "

     write(662,*) "--------------------------------"
     write(662,*) "     PSI VECTORS    "
     write(662,*) "--------------------------------"

     do i=1,mcscfnum
        write(662,'(2000F14.8)') (mmatel(i,j),j=1,mcscfnum)
     enddo


     write(662,*)
     write(662,*)
     write(662,*) "  M^2 matrix elements "
     write(662,*) "--------------------------------"
     write(662,*) "     PSI VECTORS    "
     write(662,*) "--------------------------------"

     do i=1,mcscfnum
        write(662,'(2000F14.8)') (m2matel(i,j),j=1,mcscfnum)
     enddo


     write(662,*)
     write(662,*)
     write(662,*) "  UG matrix elements "
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


     write(662,*)
     write(662,*)
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


     write(662,*)
     write(662,*)
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


     write(662,*)
     write(662,*)
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


     write(662,*)
     write(662,*)
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
     write(662,*)
     write(662,*)
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
     write(662,*)
     write(662,*)
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
  
end subroutine finalstats





