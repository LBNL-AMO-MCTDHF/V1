
!! SUBROUTINES FOR SPIN EIGENFUNCTION PROJECTOR
!!   SEE SPINWALKS.F90 ALSO

#include "Definitions.INC"


subroutine configspin_project(www,nr, vector)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: nr
  DATATYPE,intent(inout) :: vector(nr,www%firstconfig:www%lastconfig)

  call configspin_project_local(www,nr,vector(:,www%botconfig))

  if (www%parconsplit.eq.0) then
     call mpiallgather(vector,www%numconfig*nr,www%configsperproc(:)*nr,www%maxconfigsperproc*nr)
  endif

end subroutine configspin_project



subroutine configspin_project_local(www,nr,vector)
  use walkmod
  use mpimod   !! myrank
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: nr
  DATATYPE,intent(inout) :: vector(nr,www%botconfig:www%topconfig)

  call configspin_project_general(www,nr,vector,myrank,myrank)

end subroutine configspin_project_local



subroutine configspin_project_general(www,nr,vector,iproc,jproc)
  use walkmod
  use fileptrmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: iset, ii, iproc, jproc, nr, pp
  DATATYPE,intent(inout) :: vector(nr,www%allbotconfigs(iproc):www%alltopconfigs(jproc))
  DATATYPE :: smallvect(nr,www%sss%maxspinsetsize), smalltemp(nr,www%sss%maxspinsetsize), &
       outvector(nr,www%allbotconfigs(iproc):www%alltopconfigs(jproc))      !! AUTOMATIC

  if (iproc.lt.www%startrank.or.jproc.gt.www%endrank) then
     OFLWR "STARTRANK ERR PROJGEN ",iproc,www%startrank,www%endrank; CFLST
  endif

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

  if (iproc.lt.www%startrank.or.jproc.gt.www%endrank) then
     OFLWR "STARTRANK ERR tRANSTO GEN ",iproc,www%startrank,www%endrank; CFLST
  endif

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

  if (iproc.lt.www%startrank.or.jproc.gt.www%endrank) then
     OFLWR "STARTRANK ERR DFTRANSTO GEN ",iproc,www%startrank,www%endrank; CFLST
  endif

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

  if (iproc.lt.www%startrank.or.jproc.gt.www%endrank) then
     OFLWR "STARTRANK ERR TRANSFROM GEN ",iproc,www%startrank,www%endrank; CFLST
  endif

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

  if (iproc.lt.www%startrank.or.jproc.gt.www%endrank) then
     OFLWR "STARTRANK ERR DFTRANSfrom GEN ",iproc,www%startrank,www%endrank; CFLST
  endif

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









