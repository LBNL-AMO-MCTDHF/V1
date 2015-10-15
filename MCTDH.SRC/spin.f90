
!! SUBROUTINES FOR SPIN EIGENFUNCTION PROJECTOR
!!   SEE SPINWALKS.F90 ALSO

#include "Definitions.INC"


subroutine configspin_project(nr, vector)
  use parameters
  implicit none
  integer,intent(in) :: nr
  DATATYPE :: vector(nr,firstconfig:lastconfig)

  if (parconsplit.eq.0) then
     call configspin_project_all(nr,vector(:,:))
  else
     call configspin_project_local(nr,vector(:,:))
  endif

end subroutine configspin_project


subroutine configspin_project_all(nr,vector)
  use parameters
  use mpimod   !! nprocs
  implicit none
  integer,intent(in) :: nr
  DATATYPE :: vector(nr,numconfig)

  call configspin_project_general(nr,vector,1,nprocs)

end subroutine configspin_project_all


subroutine configspin_project_local(nr,vector)
  use parameters
  use mpimod   !! myrank
  implicit none
  integer,intent(in) :: nr
  DATATYPE :: vector(nr,botconfig:topconfig)

  call configspin_project_general(nr,vector,myrank,myrank)

end subroutine configspin_project_local



subroutine configspin_project_general(nr,vector,iproc,jproc)
  use spinwalkmod
  use parameters
  implicit none
  integer :: iset, ii, iproc, jproc, nr, pp
  DATATYPE :: vector(nr,allbotconfigs(iproc):alltopconfigs(jproc))
  DATATYPE :: smallvect(nr,maxspinsetsize), smalltemp(nr,maxspinsetsize), &
       outvector(nr,allbotconfigs(iproc):alltopconfigs(jproc))

  outvector(:,:) = 0.d0

  do pp=iproc,jproc
     do iset=1,numspinsets(pp)
        do ii=1,spinsetsize(iset,pp)
           smallvect(:,ii)=vector(:,spinsets(ii,iset,pp))
        enddo
        call MYGEMM('N', 'T', nr, spinsetsize(iset,pp), spinsetsize(iset,pp), DATAONE,  smallvect, nr, spinsetprojector(iset,pp)%mat, spinsetsize(iset,pp),DATAZERO,smalltemp, nr)
        do ii=1,spinsetsize(iset,pp)
           outvector(:,spinsets(ii,iset,pp)) = outvector(:,spinsets(ii,iset,pp)) + smalltemp(:,ii)
        enddo
     enddo
  enddo

  vector(:,:)=outvector(:,:)


end subroutine configspin_project_general



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
  use parameters
  use mpimod !! nprocs
  implicit none
  integer :: nblock
  DATATYPE,intent(in) :: invector(nblock,numconfig)
  DATATYPE,intent(out) :: outvector(nblock,numspinconfig)

  call configspin_transformto_general(nblock,invector,outvector,1,nprocs)

end subroutine configspin_transformto_all



subroutine configspin_transformto_local(nblock,invector,outvector)
  use parameters
  use mpimod !! myrank
  implicit none
  integer :: nblock
  DATATYPE,intent(in) :: invector(nblock,botconfig:topconfig)
  DATATYPE,intent(out) :: outvector(nblock,botspin:topspin)

  call configspin_transformto_general(nblock,invector,outvector,myrank,myrank)

end subroutine configspin_transformto_local



subroutine configspin_transformto_general(nblock,invector,outvector,iproc,jproc)
  use spinwalkmod
  use parameters
  implicit none
  integer :: nblock, iset, iind,ii,iproc,jproc,pp
  DATATYPE,intent(in) :: invector(nblock,allbotconfigs(iproc):alltopconfigs(jproc))
  DATATYPE,intent(out) :: outvector(nblock,allbotspins(iproc):alltopspins(jproc))

  DATATYPE :: smallvect(nblock,maxspinsetsize), smalltemp(nblock,maxspinsetsize)

  outvector(:,:)=0d0

  do pp=iproc,jproc
     iind=allbotspins(pp)
     do iset=1,numspinsets(pp)
        smallvect(:,:)=0d0
        do ii=1,spinsetsize(iset,pp)
           smallvect(:,ii)=invector(:,spinsets(ii,iset,pp))
        enddo
        call MYGEMM('N', 'N', nblock, spinsetrank(iset,pp), spinsetsize(iset,pp), DATAONE, smallvect,nblock, spinsetprojector(iset,pp)%vects, spinsetsize(iset,pp), DATAZERO,smalltemp, nblock)
        outvector(:,iind:iind+spinsetrank(iset,pp)-1) = smalltemp(:,1:spinsetrank(iset,pp))
        iind=iind+spinsetrank(iset,pp)
     enddo
     
     if (iind.ne.alltopspins(pp)+1) then
        OFLWR "IIND ERROto", iind,pp,allbotspins(pp),alltopspins(pp); CFLST
     endif
  enddo

end subroutine configspin_transformto_general



subroutine dfspin_transformto_all(nblock,invector,outvector)
  use parameters
  use mpimod  !! nprocs
  implicit none
  integer :: nblock
  DATATYPE,intent(in) :: invector(nblock,numdfconfigs)
  DATATYPE,intent(out) :: outvector(nblock,numspindfconfig)

  call dfspin_transformto_general(nblock,invector,outvector,1,nprocs)

end subroutine dfspin_transformto_all



subroutine dfspin_transformto_local(nblock,invector,outvector)
  use parameters
  use mpimod  !! myrank
  implicit none
  integer :: nblock
  DATATYPE,intent(in) :: invector(nblock,botdfconfig:topdfconfig)
  DATATYPE,intent(out) :: outvector(nblock,botdfspin:topdfspin)

  call dfspin_transformto_general(nblock,invector,outvector,myrank,myrank)

end subroutine dfspin_transformto_local



subroutine dfspin_transformto_general(nblock,invector,outvector,iproc,jproc)
  use spinwalkmod
  use parameters
  implicit none
  integer :: nblock, iset, iind,ii,jset,iproc,jproc,pp
  DATATYPE,intent(in) :: invector(nblock,allbotdfconfigs(iproc):alltopdfconfigs(jproc))
  DATATYPE,intent(out) :: outvector(nblock,allbotspindfs(iproc):alltopspindfs(jproc))

  DATATYPE :: smallvect(nblock,maxspinsetsize), smalltemp(nblock,maxspinsetsize)

  outvector(:,:)=0d0

  do pp=iproc,jproc
     iind=allbotspindfs(pp)

     do jset=1,numspindfsets(pp)
        iset=spindfsetindex(jset,pp)

        smallvect(:,:)=0d0
        do ii=1,spinsetsize(iset,pp)
           smallvect(:,ii)=invector(:,spindfsets(ii,jset,pp))
        enddo
        call MYGEMM('N', 'N', nblock, spinsetrank(iset,pp), spinsetsize(iset,pp), DATAONE, smallvect,nblock, spinsetprojector(iset,pp)%vects, spinsetsize(iset,pp), DATAZERO,smalltemp, nblock)
        outvector(:,iind:iind+spinsetrank(iset,pp)-1) = smalltemp(:,1:spinsetrank(iset,pp))
        iind=iind+spinsetrank(iset,pp)
     enddo

     if (iind.ne.alltopspindfs(pp)+1) then
        OFLWR "IIND ERROdfto", iind,pp,allbotspindfs(pp),alltopspindfs(pp); CFLST
     endif
  end do

end subroutine dfspin_transformto_general



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
  use parameters
  use mpimod !! nprocs
  implicit none
  integer :: nblock
  DATATYPE,intent(out) :: outvector(nblock,numconfig)
  DATATYPE,intent(in) :: invector(nblock,numspinconfig)

  call configspin_transformfrom_general(nblock,invector,outvector,1,nprocs)

end subroutine configspin_transformfrom_all


subroutine configspin_transformfrom_local(nblock,invector,outvector)
  use parameters
  use mpimod !! nprocs
  implicit none
  integer :: nblock
  DATATYPE,intent(out) :: outvector(nblock,botconfig:topconfig)
  DATATYPE,intent(in) :: invector(nblock,botspin:topspin)

  call configspin_transformfrom_general(nblock,invector,outvector,myrank,myrank)

end subroutine configspin_transformfrom_local



subroutine configspin_transformfrom_general(nblock,invector,outvector,iproc,jproc)
  use spinwalkmod
  use parameters
  implicit none
  integer :: nblock,iset, iind, ii, iproc, jproc, pp
  DATATYPE,intent(out) :: outvector(nblock,allbotconfigs(iproc):alltopconfigs(jproc))
  DATATYPE,intent(in) :: invector(nblock,allbotspins(iproc):alltopspins(jproc))
  DATATYPE :: smallvect(nblock,maxspinsetsize), smalltemp(nblock,maxspinsetsize)

  outvector(:,:)=0d0

  do pp=iproc,jproc

     iind=allbotspins(pp)

     do iset=1,numspinsets(pp)
        smallvect(:,1:spinsetrank(iset,pp))=invector(:,iind:iind+spinsetrank(iset,pp)-1)
        call MYGEMM('N', 'T', nblock, spinsetsize(iset,pp), spinsetrank(iset,pp), DATAONE, smallvect, nblock, spinsetprojector(iset,pp)%vects, spinsetsize(iset,pp), DATAZERO,smalltemp, nblock)
        do ii=1,spinsetsize(iset,pp)
           outvector(:,spinsets(ii,iset,pp))=smalltemp(:,ii)
        enddo
        iind=iind+spinsetrank(iset,pp)
     enddo
     
     if (iind.ne.alltopspins(pp)+1) then
        OFLWR "IIND ERROfrom", iind, pp, allbotspins(pp), alltopspins(pp); CFLST
     endif

  enddo

end subroutine configspin_transformfrom_general



subroutine dfspin_transformfrom_all(nblock,invector,outvector)
  use parameters
  use mpimod  !! nprocs
  implicit none
  integer :: nblock
  DATATYPE,intent(in) :: invector(nblock,numspindfconfig)
  DATATYPE,intent(out) :: outvector(nblock,numdfconfigs)

  call dfspin_transformfrom_general(nblock,invector,outvector,1,nprocs)

end subroutine dfspin_transformfrom_all



subroutine dfspin_transformfrom_local(nblock,invector,outvector)
  use parameters
  use mpimod  !! myrank
  implicit none
  integer :: nblock
  DATATYPE,intent(in) :: invector(nblock,botdfspin:topdfspin)
  DATATYPE,intent(out) :: outvector(nblock,botdfconfig:topdfconfig)

  call dfspin_transformfrom_general(nblock,invector,outvector,myrank,myrank)

end subroutine dfspin_transformfrom_local


subroutine dfspin_transformfrom_general(nblock,invector,outvector,iproc,jproc)
  use spinwalkmod
  use parameters
  implicit none
  integer :: nblock,iset, iind, ii, jset, iproc, jproc, pp
  DATATYPE,intent(out) :: outvector(nblock,allbotdfconfigs(iproc):alltopdfconfigs(jproc))
  DATATYPE,intent(in) :: invector(nblock,allbotspindfs(iproc):alltopspindfs(jproc))
  DATATYPE :: smallvect(nblock,maxspinsetsize), smalltemp(nblock,maxspinsetsize)

  outvector(:,:)=0d0

  do pp=iproc,jproc

     iind=allbotspindfs(pp)
     do jset=1,numspindfsets(pp)
        iset=spindfsetindex(jset,pp)
        smallvect(:,1:spinsetrank(iset,pp))=invector(:,iind:iind+spinsetrank(iset,pp)-1)
        call MYGEMM('N', 'T', nblock, spinsetsize(iset,pp), spinsetrank(iset,pp), DATAONE, smallvect, nblock, spinsetprojector(iset,pp)%vects, spinsetsize(iset,pp), DATAZERO,smalltemp, nblock)
        do ii=1,spinsetsize(iset,pp)
           outvector(:,spindfsets(ii,jset,pp))=smalltemp(:,ii)
        enddo
        iind=iind+spinsetrank(iset,pp)
     enddo

     if (iind.ne.alltopspindfs(pp)+1) then
        OFLWR "IIND ERROdffrom", iind, pp, allbotspindfs(pp), alltopspindfs(pp)
     endif

  enddo

end subroutine dfspin_transformfrom_general









