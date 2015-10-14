
!! SUBROUTINES FOR SPIN EIGENFUNCTION PROJECTOR
!!   SEE SPINWALKS.F90 ALSO

#include "Definitions.INC"


subroutine configspin_project(vector, iprint)
  use parameters
  implicit none
  integer :: iprint
  DATATYPE :: vector(numr,firstconfig:lastconfig)

  if (parconsplit.eq.0) then
     call configspin_project_all(vector(:,:),iprint)
  else
     call configspin_project_local(vector(:,:),iprint)
  endif

end subroutine configspin_project


subroutine configspin_project_all(vector, iprint)
  use spinwalkmod
  use parameters
  implicit none
  integer :: iprint
  DATATYPE :: vector(numr,numconfig)

  call configspin_project_local(vector(:,configstart),iprint)

  if (sparseconfigflag.ne.0) then
     call mpiallgather(vector,numconfig*numr, configsperproc*numr,maxconfigsperproc*numr)
  endif

end subroutine configspin_project_all



subroutine configspin_project_local(vector,iprint)
  use spinwalkmod
  use parameters
  implicit none
  integer :: iprint, iset, ii, isize
  real*8 :: normsq, normsq2
  DATATYPE :: hermdot, vector(numr,configstart:configend)
  DATATYPE :: smallvect(numr,maxspinsetsize), smalltemp(numr,maxspinsetsize), &
       outvector(numr,configstart:configend)

  isize=configend-configstart+1

  normsq=real(hermdot(vector,vector,numr*isize))  !! ok hermdot
  if (sparseconfigflag.ne.0) then
     call mympirealreduceone(normsq)
  endif

  outvector(:,:) = 0.d0
  do iset=1,numspinsets
     do ii=1,spinsetsize(iset)
        smallvect(:,ii)=vector(:,spinsets(ii,iset))
     enddo
     call MYGEMM('N', 'T', numr, spinsetsize(iset), spinsetsize(iset), DATAONE,  smallvect, numr, spinsetprojector(iset)%mat, spinsetsize(iset),DATAZERO,smalltemp, numr)
     do ii=1,spinsetsize(iset)
        outvector(:,spinsets(ii,iset)) = outvector(:,spinsets(ii,iset)) + smalltemp(:,ii)
     enddo
  enddo

  vector(:,:)=outvector(:,:)

  normsq2=real(hermdot(vector,vector,numr*isize))  !! ok hermdot
  if (sparseconfigflag.ne.0) then
     call mympirealreduceone(normsq2)
  endif

  if (iprint/=0) then
     if (abs(normsq/normsq2-1.d0).gt.1.d-7) then
        OFLWR "Warning, in configspin_project_local I lost norm: ", normsq, normsq2; CFL
     endif
  endif

end subroutine configspin_project_local


subroutine configspin_transformto(nblock,invector,outvector)
  use spinwalkmod
  use parameters
  implicit none
  integer,intent(in) :: nblock
  DATATYPE,intent(in) :: invector(nblock,localnconfig)
  DATATYPE,intent(out) :: outvector(nblock,localnspin)
  if (parconsplit.eq.0) then
     call configspin_transformto_all(nblock,invector,outvector)
  else
     call configspin_transformto_local(nblock,invector,outvector)
  endif
end subroutine configspin_transformto


!! ALL NUMCONFIG.

subroutine configspin_transformto_all(nblock,invector,outvector)
  use spinwalkmod
  use parameters
  implicit none
  integer :: nblock
  DATATYPE,intent(in) :: invector(nblock,numconfig)
  DATATYPE,intent(out) :: outvector(nblock,numspinconfig)

  outvector(:,:)=0d0
  call configspin_transformto_local(nblock,invector(:,configstart),outvector(:,spinstart))
  if (sparseconfigflag.ne.0) then
     call mpiallgather(outvector(:,:),numspinconfig*nblock, spinsperproc(:)*nblock,maxspinsperproc*nblock)
  endif

end subroutine configspin_transformto_all




subroutine configspin_transformto_local(nblock,invector,outvector)
  use spinwalkmod
  use parameters
  implicit none
  integer :: nblock, iset, iind,ii
  DATATYPE,intent(in) :: invector(nblock,configstart:configend)
  DATATYPE,intent(out) :: outvector(nblock,spinstart:spinend)

  DATATYPE :: smallvect(nblock,maxspinsetsize), smalltemp(nblock,maxspinsetsize)

  outvector(:,:)=0d0

  iind=spinstart
  do iset=1,numspinsets
     smallvect(:,:)=0d0
     do ii=1,spinsetsize(iset)
        smallvect(:,ii)=invector(:,spinsets(ii,iset))
     enddo
     call MYGEMM('N', 'N', nblock, spinsetrank(iset), spinsetsize(iset), DATAONE, smallvect,nblock, spinsetprojector(iset)%vects, spinsetsize(iset), DATAZERO,smalltemp, nblock)
     outvector(:,iind:iind+spinsetrank(iset)-1) = smalltemp(:,1:spinsetrank(iset))
     iind=iind+spinsetrank(iset)
  enddo

  if (iind.ne.spinend+1) then
     OFLWR "IIND ERRO", iind,spinstart,spinend; CFLST
  endif

end subroutine configspin_transformto_local



subroutine dfspin_transformto_local(nblock,invector,outvector)
  use spinwalkmod
  use parameters
  implicit none
  integer :: nblock, iset, iind,ii,jset
  DATATYPE,intent(in) :: invector(nblock,configdfstart:configdfend)
  DATATYPE,intent(out) :: outvector(nblock,spindfstart:spindfend)

  DATATYPE :: smallvect(nblock,maxspinsetsize), smalltemp(nblock,maxspinsetsize)

  outvector(:,:)=0d0

  iind=spindfstart
  do jset=1,numspindfsets
     iset=spindfsetindex(jset)

     smallvect(:,:)=0d0
     do ii=1,spinsetsize(iset)
        smallvect(:,ii)=invector(:,spindfsets(ii,jset))
     enddo
     call MYGEMM('N', 'N', nblock, spinsetrank(iset), spinsetsize(iset), DATAONE, smallvect,nblock, spinsetprojector(iset)%vects, spinsetsize(iset), DATAZERO,smalltemp, nblock)
     outvector(:,iind:iind+spinsetrank(iset)-1) = smalltemp(:,1:spinsetrank(iset))
     iind=iind+spinsetrank(iset)
  enddo

  if (iind.ne.spindfend+1) then
     OFLWR "IIND ERRO", iind,spindfstart,spindfend; CFLST
  endif

end subroutine dfspin_transformto_local



subroutine configspin_transformfrom(nblock,invector,outvector)
  use parameters
  implicit none
  integer :: nblock
  DATATYPE,intent(out) :: outvector(nblock,localnconfig)
  DATATYPE,intent(in) :: invector(nblock,localnspin)
  if (parconsplit.eq.0) then
     call configspin_transformfrom_all(nblock,invector,outvector)
  else
     call configspin_transformfrom_local(nblock,invector,outvector)
  endif
end subroutine configspin_transformfrom


subroutine configspin_transformfrom_all(nblock,invector,outvector)
  use spinwalkmod
  use parameters
  implicit none
  integer :: nblock
  DATATYPE,intent(out) :: outvector(nblock,numconfig)
  DATATYPE,intent(in) :: invector(nblock,numspinconfig)

  outvector(:,:)=0d0
  call configspin_transformfrom_local(nblock,invector(:,spinstart),outvector(:,configstart))
  if (sparseconfigflag.ne.0) then
        call mpiallgather(outvector(:,:),numconfig*nblock, configsperproc(:)*nblock,maxconfigsperproc*nblock)
  endif

end subroutine configspin_transformfrom_all



subroutine configspin_transformfrom_local(nblock,invector,outvector)
  use spinwalkmod
  use parameters
  implicit none
  integer :: nblock,iset, iind, ii
  DATATYPE,intent(out) :: outvector(nblock,configstart:configend)
  DATATYPE,intent(in) :: invector(nblock,spinstart:spinend)
  DATATYPE :: smallvect(nblock,maxspinsetsize), smalltemp(nblock,maxspinsetsize)

  outvector(:,:)=0d0

  iind=spinstart
  do iset=1,numspinsets
     smallvect(:,1:spinsetrank(iset))=invector(:,iind:iind+spinsetrank(iset)-1)
     call MYGEMM('N', 'T', nblock, spinsetsize(iset), spinsetrank(iset), DATAONE, smallvect, nblock, spinsetprojector(iset)%vects, spinsetsize(iset), DATAZERO,smalltemp, nblock)
     do ii=1,spinsetsize(iset)
        outvector(:,spinsets(ii,iset))=smalltemp(:,ii)
     enddo
     iind=iind+spinsetrank(iset)
  enddo

  if (iind.ne.spinend+1) then
     OFLWR "IIND ERROxx", iind, spinstart, spinend; CFLST
  endif

end subroutine configspin_transformfrom_local



subroutine dfspin_transformfrom_local(nblock,invector,outvector)
  use spinwalkmod
  use parameters
  implicit none
  integer :: nblock,iset, iind, ii, jset
  DATATYPE,intent(out) :: outvector(nblock,configdfstart:configdfend)
  DATATYPE,intent(in) :: invector(nblock,spindfstart:spindfend)
  DATATYPE :: smallvect(nblock,maxspinsetsize), smalltemp(nblock,maxspinsetsize)

  outvector(:,:)=0d0

  iind=spindfstart
  do jset=1,numspindfsets
     iset=spindfsetindex(jset)
     smallvect(:,1:spinsetrank(iset))=invector(:,iind:iind+spinsetrank(iset)-1)
     call MYGEMM('N', 'T', nblock, spinsetsize(iset), spinsetrank(iset), DATAONE, smallvect, nblock, spinsetprojector(iset)%vects, spinsetsize(iset), DATAZERO,smalltemp, nblock)
     do ii=1,spinsetsize(iset)
        outvector(:,spindfsets(ii,jset))=smalltemp(:,ii)
     enddo
     iind=iind+spinsetrank(iset)
  enddo

  if (iind.ne.spindfend+1) then
     OFLWR "IIND ERROxx", iind, spindfstart, spindfend; CFLST
  endif

end subroutine dfspin_transformfrom_local




function spinallowed(spinval)
  use parameters
  implicit none
  real*8 :: spinval
  logical :: spinallowed
  if (abs(spinval-(spinrestrictval/2.d0*(spinrestrictval/2.d0+1))).lt.1.d-3) then
     spinallowed=.true.
  else
     spinallowed=.false.
  endif
end function







