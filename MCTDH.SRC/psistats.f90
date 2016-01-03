
#include "Definitions.INC"


subroutine get_orbmats( myspfs,  numspf,  ugmat,   xdipmat,ydipmat,zdipmat,   xrefmat,yrefmat,zrefmat)
  use parameters
  implicit none

  integer, intent(in) :: numspf
  DATATYPE,intent(in) :: myspfs(spfsize,numspf)
  DATATYPE, intent(out) :: ugmat(numspf,numspf), xdipmat(numspf,numspf), ydipmat(numspf,numspf), zdipmat(numspf,numspf), &
       xrefmat(numspf,numspf), yrefmat(numspf,numspf), zrefmat(numspf,numspf)
  DATATYPE,allocatable ::  tempspfs(:,:),tempspfs2(:,:)
  integer :: i

  allocate(tempspfs(spfsize,numspf),tempspfs2(spfsize,numspf))
  tempspfs(:,:)=0; tempspfs2(:,:)=0

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
     call mult_zdipole(myspfs(:,i),tempspfs(:,i),0)
  enddo
  call MYGEMM(CNORMCHAR,'N',numspf,numspf,spfsize,DATAONE, myspfs, spfsize, tempspfs, spfsize, DATAZERO, zdipmat, numspf)
  if (parorbsplit.eq.3) then
     call mympireduce(zdipmat(:,:),numspf**2)
  endif

!! Y DIPOLE

  do i=1,numspf
     call mult_ydipole(myspfs(:,i),tempspfs(:,i),0)
  enddo
  call MYGEMM(CNORMCHAR,'N',numspf,numspf,spfsize,DATAONE, myspfs, spfsize, tempspfs, spfsize, DATAZERO, ydipmat, numspf)
  if (parorbsplit.eq.3) then
     call mympireduce(ydipmat(:,:),numspf**2)
  endif

!! X DIPOLE

  do i=1,numspf
     call mult_xdipole(myspfs(:,i),tempspfs(:,i),0)
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

  deallocate(tempspfs,tempspfs2)

end subroutine get_orbmats




subroutine psistats( thistime )
  use parameters
  use mpimod
  use configmod
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


  call get_psistats(www,bwwptr,yyy%cmfpsivec(spfstart,0),mcscfnum,yyy%cmfpsivec(astart(1),0),&
       mexpect,m2expect,ugexpect,   xdipole,ydipole,zdipole,   xreflect,yreflect,zreflect)

  if (myrank.eq.1) then
     open(662, file=psistatsfile, status="old", position="append")
     do i=1,mcscfnum
        write(662,'(F13.5,I3,2000F14.8)') thistime, i, mexpect(i),m2expect(i),ugexpect(i),  xdipole(i),ydipole(i),zdipole(i),  xreflect(i),yreflect(i),zreflect(i)
     enddo
     close(662)
  endif

end subroutine psistats



subroutine get_psistats( www, bioww, myspfs, numvec, inavectors, mexpect,m2expect,ugexpect,   xdipole,ydipole,zdipole,   xreflect,yreflect,zreflect)
  use parameters
  use walkmod
  implicit none
  type(walktype),intent(in) :: www, bioww
  integer, intent(in) :: numvec
  DATATYPE,intent(in) :: myspfs(spfsize,www%nspf), inavectors(numr,www%firstconfig:www%lastconfig,numvec)       
  DATATYPE, intent(out) :: mexpect(numvec),ugexpect(numvec),m2expect(numvec),&
       xreflect(numvec),yreflect(numvec),zreflect(numvec),xdipole(numvec),ydipole(numvec),zdipole(numvec)
  DATATYPE,allocatable :: ugmat(:,:), xdipmat(:,:), ydipmat(:,:), zdipmat(:,:), &
       xrefmat(:,:), yrefmat(:,:), zrefmat(:,:),     tempvector(:,:), tempspfs(:,:),tempspfs2(:,:)
  DATATYPE ::  dot,normsq(numvec),nullcomplex,dipoles(3),nucdipexpect(numvec,3),csum
  DATAECS :: rvector(numr)
  integer :: i,imc

  allocate(ugmat(www%nspf,www%nspf), xdipmat(www%nspf,www%nspf), ydipmat(www%nspf,www%nspf), zdipmat(www%nspf,www%nspf), &
       xrefmat(www%nspf,www%nspf), yrefmat(www%nspf,www%nspf), zrefmat(www%nspf,www%nspf),&
       tempvector(numr,www%firstconfig:www%lastconfig), tempspfs(spfsize,www%nspf),tempspfs2(spfsize,www%nspf))
  
  ugmat=0; xdipmat=0; ydipmat=0; zdipmat=0; xrefmat=0; yrefmat=0; zrefmat=0
  call get_orbmats( myspfs,  www%nspf,  ugmat,   xdipmat,ydipmat,zdipmat,   xrefmat,yrefmat,zrefmat)

!! M & U/G

  do imc=1,numvec
     normsq(imc)=dot(inavectors(:,:,imc),inavectors(:,:,imc),www%totadim)
  enddo
  if (www%parconsplit.ne.0) then
     call mympireduce(normsq,numvec)
  endif

  if (spfrestrictflag.eq.1) then
     do imc=1,numvec      
        do i=1,numr
           tempvector(i,:)=inavectors(i,:,imc) * www%configmvals(www%firstconfig:www%lastconfig)
        enddo
        mexpect(imc)=dot(inavectors(:,:,imc),tempvector(:,:),www%totadim)   /normsq(imc)
        
        do i=1,numr
           tempvector(i,:)=inavectors(i,:,imc) * www%configugvals(www%firstconfig:www%lastconfig)
        enddo
        ugexpect(imc)=dot(inavectors(:,:,imc),tempvector(:,:),www%totadim)   /normsq(imc)
        
        do i=1,numr
           tempvector(i,:)=inavectors(i,:,imc) * www%configmvals(www%firstconfig:www%lastconfig)**2
        enddo
        m2expect(imc)=dot(inavectors(:,:,imc),tempvector(:,:),www%totadim)   /normsq(imc)
     enddo
     if (www%parconsplit.ne.0) then
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
        call autocorrelate_one(www,bioww,inavectors(:,:,imc),myspfs,tempspfs,inavectors(:,:,imc),ugexpect(imc),numr)
        ugexpect(imc)=ugexpect(imc)   /normsq(imc)
     enddo

  endif

!! like dipolesub

!! independent of R for now.  multiply by R for prolate  (R set to 1 for atom)
  call nucdipvalue(nullcomplex,dipoles)

  do imc=1,mcscfnum
     do i=1,numr
        tempvector(i,:)=inavectors(i,:,imc)*bondpoints(i)
     enddo
     csum=dot(inavectors(:,:,imc),tempvector,www%totadim)
     if (www%parconsplit.ne.0) then
        call mympireduceone(csum)
     endif
     nucdipexpect(imc,:)=dipoles(:)*csum
  enddo

  rvector(:)=bondpoints(:)

!! Z DIPOLE

  do imc=1,numvec
     call arbitraryconfig_mult_singles(www,zdipmat,rvector,inavectors(:,:,imc),tempvector,numr)
     zdipole(imc)=dot(inavectors(:,:,imc),tempvector,www%totadim)
  enddo
  if (www%parconsplit.ne.0) then
     call mympireduce(zdipole,numvec)
  endif
  zdipole(:)=(zdipole(:)+nucdipexpect(:,3))/normsq(:)

!! Y DIPOLE

  do imc=1,numvec
     call arbitraryconfig_mult_singles(www,ydipmat,rvector,inavectors(:,:,imc),tempvector,numr)
     ydipole(imc)=dot(inavectors(:,:,imc),tempvector,www%totadim)
  enddo
  if (www%parconsplit.ne.0) then
     call mympireduce(ydipole,numvec)
  endif
  ydipole(:)=(ydipole(:)+nucdipexpect(:,2))/normsq(:)

!! X DIPOLE

  do imc=1,numvec
     call arbitraryconfig_mult_singles(www,xdipmat,rvector,inavectors(:,:,imc),tempvector,numr)
     xdipole(imc)=dot(inavectors(:,:,imc),tempvector,www%totadim)
  enddo
  if (www%parconsplit.ne.0) then
     call mympireduce(xdipole,numvec)
  endif
  xdipole(:)=(xdipole(:)+nucdipexpect(:,1))/normsq(:)

!! REFLECTIONS

  do i=1,nspf
     call op_reflectz(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,numvec
     call autocorrelate_one(www,bioww,inavectors(:,:,imc),myspfs,tempspfs,inavectors(:,:,imc),zreflect(imc),numr)
     zreflect(imc)=zreflect(imc)   /normsq(imc)
  enddo

  do i=1,nspf
     call op_reflecty(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,numvec
     call autocorrelate_one(www,bioww,inavectors(:,:,imc),myspfs,tempspfs,inavectors(:,:,imc),yreflect(imc),numr)
     yreflect(imc)=yreflect(imc)   /normsq(imc)
  enddo

  do i=1,nspf
     call op_reflectx(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,numvec
     call autocorrelate_one(www,bioww,inavectors(:,:,imc),myspfs,tempspfs,inavectors(:,:,imc),xreflect(imc),numr)
     xreflect(imc)=xreflect(imc)   /normsq(imc)
  enddo

  deallocate(ugmat, xdipmat, ydipmat, zdipmat,     xrefmat, yrefmat, zrefmat,   tempvector,tempspfs,tempspfs2)
  
end subroutine get_psistats


subroutine finalstats( )
  use parameters
  use configmod
  use xxxmod
  implicit none

  call finalstats0(yyy%cmfpsivec(spfstart:spfend,0),yyy%cmfpsivec(astart(1):aend(mcscfnum),0),www,bwwptr)

end subroutine finalstats


subroutine finalstats0(myspfs,inavectors,www,bioww )
  use parameters
  use mpimod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www,bioww
  DATATYPE,intent(in) :: myspfs(spfsize,www%nspf), inavectors(numr,www%firstconfig:www%lastconfig,mcscfnum)       
  DATATYPE :: mmatel(mcscfnum,mcscfnum),ugmatel(mcscfnum,mcscfnum),m2matel(mcscfnum,mcscfnum),&
       xrefmatel(mcscfnum,mcscfnum),yrefmatel(mcscfnum,mcscfnum),zrefmatel(mcscfnum,mcscfnum),&
       xdipmatel(mcscfnum,mcscfnum),ydipmatel(mcscfnum,mcscfnum),zdipmatel(mcscfnum,mcscfnum),&
       ovlmatel(mcscfnum,mcscfnum), nucdipmatel(mcscfnum,mcscfnum,3),&
       dot,nullcomplex,dipoles(3)
  DATATYPE,allocatable :: ugmat(:,:), xdipmat(:,:), ydipmat(:,:), zdipmat(:,:), &
       xrefmat(:,:), yrefmat(:,:), zrefmat(:,:), tempvector(:,:), tempspfs(:,:),tempspfs2(:,:)
  CNORMTYPE :: occupations(www%nspf,mcscfnum)
  DATAECS :: rvector(numr)
  integer :: i,j,imc,jmc

  call mpibarrier()
  OFLWR "   ...GO finalstats."; CFL
  call mpibarrier()

  allocate(ugmat(www%nspf,www%nspf), xdipmat(www%nspf,www%nspf), ydipmat(www%nspf,www%nspf), zdipmat(www%nspf,www%nspf), &
       xrefmat(www%nspf,www%nspf), yrefmat(www%nspf,www%nspf), zrefmat(www%nspf,www%nspf),&
       tempvector(numr,www%firstconfig:www%lastconfig), tempspfs(spfsize,www%nspf),tempspfs2(spfsize,www%nspf))

  tempvector(:,:)=0d0; tempspfs(:,:)=0d0; tempspfs2(:,:)=0d0

  ugmat=0; xdipmat=0; ydipmat=0; zdipmat=0; xrefmat=0; yrefmat=0; zrefmat=0
  call get_orbmats( myspfs,  www%nspf,  ugmat,   xdipmat,ydipmat,zdipmat,   xrefmat,yrefmat,zrefmat)


  do imc=1,mcscfnum
     call getoccupations(www,inavectors(:,:,imc),numr,occupations(:,imc))
  enddo

!! M & U/G

  call mpibarrier()

  do imc=1,mcscfnum
     do jmc=1,mcscfnum
        ovlmatel(jmc,imc)=dot(inavectors(:,:,jmc),inavectors(:,:,imc),www%totadim)
     enddo
  enddo
  if (www%parconsplit.ne.0) then
     call mympireduce(ovlmatel,mcscfnum**2)
  endif

  if (spfrestrictflag.eq.1) then
     do imc=1,mcscfnum
        do i=1,numr
           tempvector(i,:)=inavectors(i,:,imc) * www%configmvals(:)
        enddo
        do jmc=1,mcscfnum
           mmatel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector(:,:),www%totadim)
        enddo
        do i=1,numr
           tempvector(i,:)=inavectors(i,:,imc) * www%configugvals(:)
        enddo
        do jmc=1,mcscfnum
           ugmatel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector(:,:),www%totadim)
        enddo
        do i=1,numr
           tempvector(i,:)=inavectors(i,:,imc) * www%configmvals(:)**2
        enddo
        do jmc=1,mcscfnum
           m2matel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector(:,:),www%totadim)
        enddo
     enddo
     if (www%parconsplit.ne.0) then
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
           call autocorrelate_one(www,bioww,inavectors(:,:,jmc),myspfs,tempspfs,inavectors(:,:,imc),ugmatel(jmc,imc),numr)
        enddo
     enddo
  endif

  call mpibarrier()

!! independent of R for now.  multiply by R for prolate  (R set to 1 for atom)
  call nucdipvalue(nullcomplex,dipoles)

  do imc=1,mcscfnum
     do i=1,numr
        tempvector(i,:)=inavectors(i,:,imc)*bondpoints(i)
     enddo
     do jmc=1,mcscfnum
        nucdipmatel(jmc,imc,:)=dipoles(:)*dot(inavectors(:,:,jmc),tempvector,www%totadim)
     enddo
  enddo
  if (www%parconsplit.ne.0) then
     call mympireduce(nucdipmatel,3*mcscfnum**2)
  endif

  rvector(:)=bondpoints(:)

!! Z DIPOLE

  do imc=1,mcscfnum
     call arbitraryconfig_mult_singles(www,zdipmat,rvector,inavectors(:,:,imc),tempvector,numr)
     do jmc=1,mcscfnum
        zdipmatel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector,www%totadim)
     enddo
  enddo
  if (www%parconsplit.ne.0) then
     call mympireduce(zdipmatel(:,:),mcscfnum**2)
  endif
  zdipmatel(:,:)=zdipmatel(:,:)+nucdipmatel(:,:,3)

!! Y DIPOLE

  do imc=1,mcscfnum
     call arbitraryconfig_mult_singles(www,ydipmat,rvector,inavectors(:,:,imc),tempvector,numr)
     do jmc=1,mcscfnum
        ydipmatel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector,www%totadim)
     enddo
  enddo
  if (www%parconsplit.ne.0) then
     call mympireduce(ydipmatel(:,:),mcscfnum**2)
  endif
  ydipmatel(:,:)=ydipmatel(:,:)+nucdipmatel(:,:,2)

!! X DIPOLE

  do imc=1,mcscfnum
     call arbitraryconfig_mult_singles(www,xdipmat,rvector,inavectors(:,:,imc),tempvector,numr)
     do jmc=1,mcscfnum
        xdipmatel(jmc,imc)=dot(inavectors(:,:,jmc),tempvector,www%totadim)
     enddo
  enddo
  if (www%parconsplit.ne.0) then
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
        call autocorrelate_one(www,bioww,inavectors(:,:,jmc),myspfs,tempspfs,inavectors(:,:,imc),zrefmatel(jmc,imc),numr)
     enddo
  enddo

  do i=1,nspf
     call op_reflecty(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,mcscfnum
     do jmc=1,mcscfnum
        call autocorrelate_one(www,bioww,inavectors(:,:,jmc),myspfs,tempspfs,inavectors(:,:,imc),yrefmatel(jmc,imc),numr)
     enddo
  enddo

  do i=1,nspf
     call op_reflectx(myspfs(:,i),tempspfs(:,i))
  enddo
  do imc=1,mcscfnum
     do jmc=1,mcscfnum
        call autocorrelate_one(www,bioww,inavectors(:,:,jmc),myspfs,tempspfs,inavectors(:,:,imc),xrefmatel(jmc,imc),numr)
     enddo
  enddo

  if (myrank.eq.1) then
     open(662, file=finalstatsfile, status="unknown")

     write(662,*);     write(662,*)
     write(662,*) "--- EXPECTATION VALUES PSI VECTORS ---"

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
     write(662,*) "  inversion expectation values psi vectors "

     do i=1,mcscfnum
        write(662,'(I10,2000F14.8)') i, ugmatel(i,i)
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

  deallocate(ugmat, xdipmat, ydipmat, zdipmat,     xrefmat, yrefmat, zrefmat,    tempvector,tempspfs,tempspfs2)
  call mpibarrier()

end subroutine finalstats0





