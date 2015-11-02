
!! SUBROUTINES FOR SPIN EIGENFUNCTION PROJECTOR
!!   SEE SPINWALKS.F90 ALSO

#include "Definitions.INC"


subroutine configspin_project(www,nr, vector)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: nr
  DATATYPE :: vector(nr,www%firstconfig:www%lastconfig)

  if (www%parconsplit.eq.0) then
     call configspin_project_all(www,nr,vector(:,:))
  else
     call configspin_project_local(www,nr,vector(:,:))
  endif

end subroutine configspin_project


subroutine configspin_project_all(www,nr,vector)
  use walkmod
  use mpimod   !! nprocs
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: nr
  DATATYPE :: vector(nr,www%numconfig)

  call configspin_project_general(www,nr,vector,1,nprocs)

end subroutine configspin_project_all


subroutine configspin_project_local(www,nr,vector)
  use walkmod
  use mpimod   !! myrank
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: nr
  DATATYPE :: vector(nr,www%botconfig:www%topconfig)

  call configspin_project_general(www,nr,vector,myrank,myrank)

end subroutine configspin_project_local



subroutine configspin_project_general(www,nr,vector,iproc,jproc)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: iset, ii, iproc, jproc, nr, pp
  DATATYPE :: vector(nr,www%allbotconfigs(iproc):www%alltopconfigs(jproc))
  DATATYPE :: smallvect(nr,www%sss%maxspinsetsize), smalltemp(nr,www%sss%maxspinsetsize), &
       outvector(nr,www%allbotconfigs(iproc):www%alltopconfigs(jproc))

  outvector(:,:) = 0.d0

  do pp=iproc,jproc
     do iset=1,www%sss%numspinsets(pp)
        do ii=1,www%sss%spinsetsize(iset,pp)
           smallvect(:,ii)=vector(:,www%sss%spinsets(ii,iset,pp))
        enddo
        call MYGEMM('N', 'T', nr, www%sss%spinsetsize(iset,pp), www%sss%spinsetsize(iset,pp), &
             DATAONE,  smallvect, nr, www%sss%spinsetprojector(iset,pp)%mat, &
             www%sss%spinsetsize(iset,pp),DATAZERO,smalltemp, nr)
        do ii=1,www%sss%spinsetsize(iset,pp)
           outvector(:,www%sss%spinsets(ii,iset,pp)) = &
                outvector(:,www%sss%spinsets(ii,iset,pp)) + smalltemp(:,ii)
        enddo
     enddo
  enddo

  vector(:,:)=outvector(:,:)

end subroutine configspin_project_general



subroutine configspin_transformto(www,nblock,invector,outvector)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: nblock
  DATATYPE,intent(in) :: invector(nblock,www%localnconfig)
  DATATYPE,intent(out) :: outvector(nblock,www%sss%localnspin)
  if (www%parconsplit.eq.0) then
     call configspin_transformto_all(www,nblock,invector,outvector)
  else
     call configspin_transformto_local(www,nblock,invector,outvector)
  endif
end subroutine configspin_transformto


!! ALL NUMCONFIG.

subroutine configspin_transformto_all(www,nblock,invector,outvector)
  use walkmod
  use mpimod !! nprocs
  implicit none
  type(walktype),intent(in) :: www
  integer :: nblock
  DATATYPE,intent(in) :: invector(nblock,www%numconfig)
  DATATYPE,intent(out) :: outvector(nblock,www%sss%numspinconfig)

  call configspin_transformto_general(www,nblock,invector,outvector,1,nprocs)

end subroutine configspin_transformto_all



subroutine configspin_transformto_local(www,nblock,invector,outvector)
  use walkmod
  use mpimod !! myrank
  implicit none
  type(walktype),intent(in) :: www
  integer :: nblock
  DATATYPE,intent(in) :: invector(nblock,www%botconfig:www%topconfig)
  DATATYPE,intent(out) :: outvector(nblock,www%sss%botspin:www%sss%topspin)

  call configspin_transformto_general(www,nblock,invector,outvector,myrank,myrank)

end subroutine configspin_transformto_local



subroutine configspin_transformto_general(www,nblock,invector,outvector,iproc,jproc)
  use fileptrmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: nblock, iset, iind,ii,iproc,jproc,pp
  DATATYPE,intent(in) :: invector(nblock,www%allbotconfigs(iproc):www%alltopconfigs(jproc))
  DATATYPE,intent(out) :: outvector(nblock,www%sss%allbotspins(iproc):www%sss%alltopspins(jproc))

  DATATYPE :: smallvect(nblock,www%sss%maxspinsetsize), smalltemp(nblock,www%sss%maxspinsetsize)

  outvector(:,:)=0d0

  do pp=iproc,jproc
     iind=www%sss%allbotspins(pp)
     do iset=1,www%sss%numspinsets(pp)
        smallvect(:,:)=0d0
        do ii=1,www%sss%spinsetsize(iset,pp)
           smallvect(:,ii)=invector(:,www%sss%spinsets(ii,iset,pp))
        enddo
        call MYGEMM('N', 'N', nblock, www%sss%spinsetrank(iset,pp), www%sss%spinsetsize(iset,pp),&
             DATAONE, smallvect,nblock, www%sss%spinsetprojector(iset,pp)%vects, &
             www%sss%spinsetsize(iset,pp), DATAZERO,smalltemp, nblock)
        outvector(:,iind:iind+www%sss%spinsetrank(iset,pp)-1) = smalltemp(:,1:www%sss%spinsetrank(iset,pp))
        iind=iind+www%sss%spinsetrank(iset,pp)
     enddo
     
     if (iind.ne.www%sss%alltopspins(pp)+1) then
        OFLWR "IIND ERROto", iind,pp,www%sss%allbotspins(pp),www%sss%alltopspins(pp); CFLST
     endif
  enddo

end subroutine configspin_transformto_general


subroutine dfspin_transformto_all(www,nblock,invector,outvector)
  use walkmod
  use mpimod  !! nprocs
  implicit none
  type(walktype),intent(in) :: www
  integer :: nblock
  DATATYPE,intent(in) :: invector(nblock,www%numdfconfigs)
  DATATYPE,intent(out) :: outvector(nblock,www%sss%numspindfconfig)

  call dfspin_transformto_general(www,nblock,invector,outvector,1,nprocs)

end subroutine dfspin_transformto_all



subroutine dfspin_transformto_local(www,nblock,invector,outvector)
  use walkmod
  use mpimod  !! myrank
  implicit none
  type(walktype),intent(in) :: www
  integer :: nblock
  DATATYPE,intent(in) :: invector(nblock,www%botdfconfig:www%topdfconfig)
  DATATYPE,intent(out) :: outvector(nblock,www%sss%botdfspin:www%sss%topdfspin)

  call dfspin_transformto_general(www,nblock,invector,outvector,myrank,myrank)

end subroutine dfspin_transformto_local



subroutine dfspin_transformto_general(www,nblock,invector,outvector,iproc,jproc)
  use fileptrmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: nblock, iset, iind,ii,jset,iproc,jproc,pp
  DATATYPE,intent(in) :: invector(nblock,www%allbotdfconfigs(iproc):www%alltopdfconfigs(jproc))
  DATATYPE,intent(out) :: outvector(nblock,www%sss%allbotspindfs(iproc):www%sss%alltopspindfs(jproc))

  DATATYPE :: smallvect(nblock,www%sss%maxspinsetsize), smalltemp(nblock,www%sss%maxspinsetsize)

  outvector(:,:)=0d0

  do pp=iproc,jproc
     iind=www%sss%allbotspindfs(pp)

     do jset=1,www%sss%numspindfsets(pp)
        iset=www%sss%spindfsetindex(jset,pp)

        smallvect(:,:)=0d0
        do ii=1,www%sss%spinsetsize(iset,pp)
           smallvect(:,ii)=invector(:,www%sss%spindfsets(ii,jset,pp))
        enddo
        call MYGEMM('N', 'N', nblock, www%sss%spinsetrank(iset,pp), www%sss%spinsetsize(iset,pp),&
             DATAONE, smallvect,nblock, www%sss%spinsetprojector(iset,pp)%vects, &
             www%sss%spinsetsize(iset,pp), DATAZERO,smalltemp, nblock)
        outvector(:,iind:iind+www%sss%spinsetrank(iset,pp)-1) = smalltemp(:,1:www%sss%spinsetrank(iset,pp))
        iind=iind+www%sss%spinsetrank(iset,pp)
     enddo

     if (iind.ne.www%sss%alltopspindfs(pp)+1) then
        OFLWR "IIND ERROdfto", iind,pp,www%sss%allbotspindfs(pp),www%sss%alltopspindfs(pp); CFLST
     endif
  end do

end subroutine dfspin_transformto_general



subroutine configspin_transformfrom(www,nblock,invector,outvector)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: nblock
  DATATYPE,intent(out) :: outvector(nblock,www%localnconfig)
  DATATYPE,intent(in) :: invector(nblock,www%sss%localnspin)
  if (www%parconsplit.eq.0) then
     call configspin_transformfrom_all(www,nblock,invector,outvector)
  else
     call configspin_transformfrom_local(www,nblock,invector,outvector)
  endif
end subroutine configspin_transformfrom


subroutine configspin_transformfrom_all(www,nblock,invector,outvector)
  use walkmod
  use mpimod !! nprocs
  implicit none
  type(walktype),intent(in) :: www
  integer :: nblock
  DATATYPE,intent(out) :: outvector(nblock,www%numconfig)
  DATATYPE,intent(in) :: invector(nblock,www%sss%numspinconfig)

  call configspin_transformfrom_general(www,nblock,invector,outvector,1,nprocs)

end subroutine configspin_transformfrom_all


subroutine configspin_transformfrom_local(www,nblock,invector,outvector)
  use walkmod
  use mpimod !! nprocs
  implicit none
  type(walktype),intent(in) :: www
  integer :: nblock
  DATATYPE,intent(out) :: outvector(nblock,www%botconfig:www%topconfig)
  DATATYPE,intent(in) :: invector(nblock,www%sss%botspin:www%sss%topspin)

  call configspin_transformfrom_general(www,nblock,invector,outvector,myrank,myrank)

end subroutine configspin_transformfrom_local



subroutine configspin_transformfrom_general(www,nblock,invector,outvector,iproc,jproc)
  use fileptrmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: nblock,iset, iind, ii, iproc, jproc, pp
  DATATYPE,intent(out) :: outvector(nblock,www%allbotconfigs(iproc):www%alltopconfigs(jproc))
  DATATYPE,intent(in) :: invector(nblock,www%sss%allbotspins(iproc):www%sss%alltopspins(jproc))
  DATATYPE :: smallvect(nblock,www%sss%maxspinsetsize), smalltemp(nblock,www%sss%maxspinsetsize)

  outvector(:,:)=0d0

  do pp=iproc,jproc

     iind=www%sss%allbotspins(pp)

     do iset=1,www%sss%numspinsets(pp)
        smallvect(:,1:www%sss%spinsetrank(iset,pp))=invector(:,iind:iind+www%sss%spinsetrank(iset,pp)-1)
        call MYGEMM('N', 'T', nblock, www%sss%spinsetsize(iset,pp), www%sss%spinsetrank(iset,pp), &
             DATAONE, smallvect, nblock, www%sss%spinsetprojector(iset,pp)%vects, &
             www%sss%spinsetsize(iset,pp), DATAZERO,smalltemp, nblock)
        do ii=1,www%sss%spinsetsize(iset,pp)
           outvector(:,www%sss%spinsets(ii,iset,pp))=smalltemp(:,ii)
        enddo
        iind=iind+www%sss%spinsetrank(iset,pp)
     enddo
     
     if (iind.ne.www%sss%alltopspins(pp)+1) then
        OFLWR "IIND ERROfrom", iind, pp, www%sss%allbotspins(pp), www%sss%alltopspins(pp); CFLST
     endif

  enddo

end subroutine configspin_transformfrom_general



subroutine dfspin_transformfrom_all(www,nblock,invector,outvector)
  use walkmod
  use mpimod  !! nprocs
  implicit none
  type(walktype),intent(in) :: www
  integer :: nblock
  DATATYPE,intent(in) :: invector(nblock,www%sss%numspindfconfig)
  DATATYPE,intent(out) :: outvector(nblock,www%numdfconfigs)

  call dfspin_transformfrom_general(www,nblock,invector,outvector,1,nprocs)

end subroutine dfspin_transformfrom_all



subroutine dfspin_transformfrom_local(www,nblock,invector,outvector)
  use walkmod
  use mpimod  !! myrank
  implicit none
  type(walktype),intent(in) :: www
  integer :: nblock
  DATATYPE,intent(in) :: invector(nblock,www%sss%botdfspin:www%sss%topdfspin)
  DATATYPE,intent(out) :: outvector(nblock,www%botdfconfig:www%topdfconfig)

  call dfspin_transformfrom_general(www,nblock,invector,outvector,myrank,myrank)

end subroutine dfspin_transformfrom_local


subroutine dfspin_transformfrom_general(www,nblock,invector,outvector,iproc,jproc)
  use fileptrmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: nblock,iset, iind, ii, jset, iproc, jproc, pp
  DATATYPE,intent(out) :: outvector(nblock,www%allbotdfconfigs(iproc):www%alltopdfconfigs(jproc))
  DATATYPE,intent(in) :: invector(nblock,www%sss%allbotspindfs(iproc):www%sss%alltopspindfs(jproc))
  DATATYPE :: smallvect(nblock,www%sss%maxspinsetsize), smalltemp(nblock,www%sss%maxspinsetsize)

  outvector(:,:)=0d0

  do pp=iproc,jproc

     iind=www%sss%allbotspindfs(pp)
     do jset=1,www%sss%numspindfsets(pp)
        iset=www%sss%spindfsetindex(jset,pp)
        smallvect(:,1:www%sss%spinsetrank(iset,pp))=invector(:,iind:iind+www%sss%spinsetrank(iset,pp)-1)
        call MYGEMM('N', 'T', nblock, www%sss%spinsetsize(iset,pp), www%sss%spinsetrank(iset,pp),&
             DATAONE, smallvect, nblock, www%sss%spinsetprojector(iset,pp)%vects, &
             www%sss%spinsetsize(iset,pp), DATAZERO,smalltemp, nblock)
        do ii=1,www%sss%spinsetsize(iset,pp)
           outvector(:,www%sss%spindfsets(ii,jset,pp))=smalltemp(:,ii)
        enddo
        iind=iind+www%sss%spinsetrank(iset,pp)
     enddo

     if (iind.ne.www%sss%alltopspindfs(pp)+1) then
        OFLWR "IIND ERROdffrom", iind, pp, www%sss%allbotspindfs(pp), www%sss%alltopspindfs(pp)
     endif

  enddo

end subroutine dfspin_transformfrom_general









