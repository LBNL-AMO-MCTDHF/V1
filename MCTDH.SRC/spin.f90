
!! ALL MODULES

!! SUBROUTINES FOR BASIS SET TRANSFORMATIONS AND PROJECTIONS (spin and restricted configuration space)

#include "Definitions.INC"


module spininitmod
contains

subroutine basis_set(www,innzflag)
  use r_parameters
  use walkmod
  use mpimod !! nprocs, MPI_COMM_WORLD
  implicit none
  type(walktype),intent(inout) :: www
  integer,intent(in) :: innzflag
  integer :: ii,jj

  allocate(www%basisperproc(nprocs),www%dfbasisperproc(nprocs))
  allocate(www%allbotbasis(nprocs),www%allbotdfbasis(nprocs))
  allocate(www%alltopbasis(nprocs),www%alltopdfbasis(nprocs))

  if (www%allspinproject.eq.0) then
     www%numbasis=www%numconfig
     www%basisperproc(:)=www%configsperproc(:)
     www%allbotbasis(:)=www%allbotconfigs(:)
     www%alltopbasis(:)=www%alltopconfigs(:)
     www%maxbasisperproc=www%maxconfigsperproc
     www%botbasis=www%botconfig;  www%topbasis=www%topconfig

     www%numdfbasis=www%numdfconfigs
     www%dfbasisperproc(:)=www%dfconfsperproc(:)
     www%allbotdfbasis(:)=www%allbotdfconfigs(:)
     www%alltopdfbasis(:)=www%alltopdfconfigs(:)
     www%maxdfbasisperproc=www%maxdfconfsperproc
     www%botdfbasis=www%botdfconfig;  www%topdfbasis=www%topdfconfig
        
  else
     www%numbasis=www%sss%numspinconfig
     www%basisperproc(:)=www%sss%spinsperproc(:)
     www%allbotbasis(:)=www%sss%allbotspins(:)
     www%alltopbasis(:)=www%sss%alltopspins(:)
     www%maxbasisperproc=www%sss%maxspinsperproc
     www%botbasis=www%sss%botspin;  www%topbasis=www%sss%topspin

     www%numdfbasis=www%sss%numspindfconfig
     www%dfbasisperproc(:)=www%sss%spindfsperproc(:)
     www%allbotdfbasis(:)=www%sss%allbotspindfs(:)
     www%alltopdfbasis(:)=www%sss%alltopspindfs(:)
     www%maxdfbasisperproc=www%sss%maxspindfsperproc
     www%botdfbasis=www%sss%botdfspin;  www%topdfbasis=www%sss%topdfspin
  endif

  www%totadim=www%localnconfig*numr

  allocate(www%nzconfsperproc(nprocs),www%nzproclist(nprocs))
  www%nzconfsperproc(:)=(-99); www%nzproclist(:)=(-99);  www%nzrank=(-99)

  jj=0
  do ii=1,nprocs
     if (www%dfbasisperproc(ii).gt.0.or.innzflag.eq.0) then
        jj=jj+1
        www%nzconfsperproc(jj)=www%configsperproc(ii)
        www%nzproclist(jj)=ii
        if (ii.eq.myrank) then
           www%nzrank=jj
        endif
     endif
  enddo
  www%nzprocs=jj
  if (innzflag.eq.0) then
#ifdef MPIFLAG
     www%NZ_COMM = MPI_COMM_WORLD
#else
     www%NZ_COMM = (-798)
#endif
     if (www%nzprocs.ne.nprocs) then
        print*, "FAILFAIL FAIL",myrank,www%nzprocs; stop
     endif
  else
     call make_mpi_comm(www%nzprocs,www%nzproclist(:),www%NZ_COMM)
  endif

end subroutine basis_set


! "included" configurations are those kept in the calculation.
! there are numdfconfigs INCLUDED configurations.

subroutine init_dfcon(www)
  use fileptrmod
  use walkmod
  use mpimod
  use configsubmod
  use mpisubmod
  implicit none
  type(walktype),intent(inout) :: www
  integer :: i,j,iconfig,nondfconfigs,dfrank,ii
  integer,allocatable :: dfnotconfigs(:)

  OFLWR "Allocating more arrays for Slater determinants..."; CFL
  call waitawhile()
  call mpibarrier()

  allocate(www%ddd%dfincludedmask(www%numconfig), www%ddd%dfincludedconfigs(www%numconfig),&
!!$SP       www%ddd%dfincludedindex(www%numconfig),&
       dfnotconfigs(www%numconfig),  www%dfconfsperproc(nprocs), &
       www%allbotdfconfigs(nprocs),www%alltopdfconfigs(nprocs))
  www%ddd%dfincludedmask=0; www%ddd%dfincludedconfigs=0; 
  dfnotconfigs=0; www%dfconfsperproc=0; www%allbotdfconfigs=0; www%alltopdfconfigs=0
!!$SP www%ddd%dfincludedindex=0;

  call waitawhile()
  call mpibarrier()
  OFLWR "     .. OK allocating."; CFL

  if (www%dfrestrictflag.lt.www%dflevel) then
     OFLWR "error, set dfrestrictflag .ge. dflevel",www%dfrestrictflag,www%dflevel; CFLST
  endif

  if (www%dfrestrictflag.eq.www%dflevel) then
     www%ddd%dfincludedmask(:)=1; dfnotconfigs(:)=(-1)
     do i=1,www%numconfig
        www%ddd%dfincludedconfigs(i)=i
!!$SP        www%ddd%dfincludedindex(i)=i
     enddo
     www%numdfconfigs=www%numconfig; nondfconfigs=0
     www%dfconfsperproc(:)=www%configsperproc(:)
     www%allbotdfconfigs(:)=www%allbotconfigs(:)
     www%alltopdfconfigs(:)=www%alltopconfigs(:)
     www%maxdfconfsperproc=www%maxconfigsperproc
     www%botdfconfig=www%botconfig;  www%topdfconfig=www%topconfig

  else

     www%ddd%dfincludedmask(:)=0;  www%ddd%dfincludedconfigs(:)=(-1);  dfnotconfigs(:)=(-1)
!!$SP     www%ddd%dfincludedindex(:)=(-999)

     www%numdfconfigs=0;  nondfconfigs=0

     do i=1,www%numconfig
        if (allowedconfig0(www,www%configlist(:,i),www%dfrestrictflag)) then 
           www%numdfconfigs=www%numdfconfigs+1
           www%ddd%dfincludedmask(i)=1;        www%ddd%dfincludedconfigs(www%numdfconfigs)=i
!!$SP           www%ddd%dfincludedindex(i)=www%numdfconfigs
        else
           nondfconfigs=nondfconfigs+1;        dfnotconfigs(nondfconfigs)=i
        endif
     enddo

     if (www%numdfconfigs.eq.0) then
        OFLWR "Error, no included DF configs!"; CFLST
     endif

     call waitawhile()
     call mpibarrier()
     OFLWR "Number of included configurations ", www%numdfconfigs
     WRFL  "Number of not included ", nondfconfigs;  WRFL; CFL

     dfrank=0
     do i=www%botconfig,www%topconfig
        if (allowedconfig0(www,www%configlist(:,i),www%dfrestrictflag)) then 
           dfrank=dfrank+1
        endif
     enddo

     www%dfconfsperproc(myrank)=dfrank
     www%maxdfconfsperproc=0
     ii=0
     do i=1,nprocs
        www%allbotdfconfigs(i)=ii+1
        call mympiibcastone(www%dfconfsperproc(i),i)
        ii=ii+www%dfconfsperproc(i)
        www%alltopdfconfigs(i)=ii
        if (www%dfconfsperproc(i).gt.www%maxdfconfsperproc) then
           www%maxdfconfsperproc=www%dfconfsperproc(i)
        endif
     enddo
     www%botdfconfig=www%allbotdfconfigs(myrank)
     www%topdfconfig=www%alltopdfconfigs(myrank)

  endif

  call mpibarrier()

  www%ddd%numdfwalks=0

  if (www%dfrestrictflag.eq.www%dflevel) then

     allocate(www%ddd%dfwalkfrom(1), www%ddd%dfwalkto(1), www%ddd%includedorb(1), &
          www%ddd%excludedorb(1), www%ddd%dfwalkphase(1))
     www%ddd%dfwalkfrom(:)=(-1); www%ddd%dfwalkto(:)=(-1); www%ddd%includedorb(:)=(-1)
     www%ddd%excludedorb(:)=(-1); www%ddd%dfwalkphase(:)=(-798)

  else

     do i=1,NONdfconfigs
        iconfig=dfNOTconfigs(i)
        if (iconfig.ge.www%configstart.and.iconfig.le.www%configend) then
           do j=1,www%numsinglewalks(iconfig)
              if (allowedconfig0(www,www%configlist(:,www%singlewalk(j+www%scol(iconfig))),www%dfrestrictflag)) then
                 www%ddd%numdfwalks=www%ddd%numdfwalks+1
              endif
           enddo
        endif
     enddo

     call waitawhile()
     call mpibarrier()
     OFLWR "numdfwalks:  ", www%ddd%numdfwalks," is total number DF single excitations on this processor"
     WRFL; CFL

     allocate(www%ddd%dfwalkfrom(www%ddd%numdfwalks), www%ddd%dfwalkto(www%ddd%numdfwalks),&
          www%ddd%includedorb(www%ddd%numdfwalks), www%ddd%excludedorb(www%ddd%numdfwalks), &
          www%ddd%dfwalkphase(www%ddd%numdfwalks))
     if (www%ddd%numdfwalks.gt.0) then
        www%ddd%dfwalkfrom=0; www%ddd%dfwalkto=0; www%ddd%includedorb=0; 
        www%ddd%excludedorb=0; www%ddd%dfwalkphase=0
     endif

     call waitawhile()
     call mpibarrier()
     OFLWR "       ...allocated more..."; CFL

     www%ddd%numdfwalks=0
     do i=1,NONdfconfigs
        iconfig=dfNOTconfigs(i)
        if (iconfig.ge.www%configstart.and.iconfig.le.www%configend) then
           do j=1,www%numsinglewalks(iconfig)
              if (allowedconfig0(www,www%configlist(:,www%singlewalk(j+www%scol(iconfig))),www%dfrestrictflag)) then
                 www%ddd%numdfwalks=www%ddd%numdfwalks+1
              
                 www%ddd%dfwalkTO(www%ddd%numdfwalks)=iconfig
                 www%ddd%dfwalkFROM(www%ddd%numdfwalks)=www%singlewalk(j+www%scol(iconfig))
                 www%ddd%EXcludedorb(www%ddd%numdfwalks)=www%singlewalkopspf(1,j+www%scol(iconfig))
                 www%ddd%INcludedorb(www%ddd%numdfwalks)=www%singlewalkopspf(2,j+www%scol(iconfig))
              
                 www%ddd%dfwalkphase(www%ddd%numdfwalks)=www%singlewalkdirphase(j+www%scol(iconfig))
              endif
           enddo
        endif
     enddo

     OFLWR "numdfwalks:  ", www%ddd%numdfwalks," is really number DF single excitations on this processor"
     WRFL; CFL
     
  endif

  deallocate(dfnotconfigs)

end subroutine init_dfcon

subroutine dfcondealloc(www)
  use walkmod
  implicit none
  type(walktype) :: www
  deallocate(www%ddd%dfincludedmask,www%ddd%dfincludedconfigs)
  deallocate(www%ddd%dfwalkfrom, www%ddd%dfwalkto, www%ddd%includedorb, &
       www%ddd%excludedorb,www%ddd%dfwalkphase)
  www%ddd%numdfwalks=-1;  www%numdfconfigs=-1
end subroutine dfcondealloc

end module spininitmod


module basissubmod
contains

!! Returns df configuration index (dfrestricted) up to numdfconfigs
!! given configuration index config1 (full list, numconfig)
!! if config1 is included in the dfrestricted list
!! or -1 if config1 is not in the dfrestricted config list

function getdfindex(www,config1)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: config1
  integer :: getdfindex, j,flag, dir,newdir, step

  getdfindex=-1

  if (www%ddd%dfincludedmask(config1).eq.0) then
     return
  endif
  
  dir=1;  j=1;  step=max(1,www%numdfconfigs/4);  flag=0

  do while (flag.eq.0)
     flag=1

     if (www%ddd%dfincludedconfigs(j) .ne. config1) then
        flag=0
        if (www%ddd%dfincludedconfigs(j).lt.config1) then
           newdir=1
        else
           newdir=-1
        endif
        if (newdir.ne.dir) then
           step=max(1,step/2)
        endif
        dir=newdir;           j=j+dir*step
        if (j.le.1) then
           j=1;              dir=1
        endif
        if (j.ge.www%numdfconfigs) then
           j=www%numdfconfigs;              dir=-1
        endif
     endif
  enddo

  getdfindex=j
  
end function getdfindex


subroutine configspin_project(www,nr, vector)
  use walkmod
  use mpisubmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: nr
  DATATYPE,intent(inout) :: vector(nr,www%firstconfig:www%lastconfig)

  if (www%configend.ge.www%configstart) then
     call configspin_project_general(www,nr,&
          vector(:,www%configstart:www%configend),www%startrank,www%endrank)
  endif
  if (www%parconsplit.eq.0.and.www%sparseconfigflag.ne.0) then
     call mpiallgather(vector,www%numconfig*nr,www%configsperproc(:)*nr,www%maxconfigsperproc*nr)
  endif

end subroutine configspin_project



subroutine configspin_project_general(www,nr,vector,iproc,jproc)
  use walkmod
  use fileptrmod
  implicit none
  integer,intent(in) :: nr,iproc,jproc
  type(walktype),intent(in) :: www
  DATATYPE,intent(inout) :: vector(nr,www%allbotconfigs(iproc):www%alltopconfigs(jproc))
  DATATYPE :: smallvect(nr,www%sss%maxspinsetsize), smalltemp(nr,www%sss%maxspinsetsize), &
       outvector(nr,www%allbotconfigs(iproc):www%alltopconfigs(jproc)+1)      !! AUTOMATIC
  integer :: iset, ii, pp

  if (iproc.lt.www%startrank.or.jproc.gt.www%endrank) then
     OFLWR "STARTRANK ERR PROJGEN ",iproc,www%startrank,www%endrank; CFLST
  endif

  outvector(:,:) = 0d0

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

  if (www%alltopconfigs(jproc).ge.www%allbotconfigs(iproc)) then
     vector(:,:)=outvector(:,www%allbotconfigs(iproc):www%alltopconfigs(jproc))
  endif

end subroutine configspin_project_general



subroutine configspin_transformto_general(www,nblock,invector,outvector,iproc,jproc)
  use fileptrmod
  use walkmod
  implicit none
  integer,intent(in) :: nblock,iproc,jproc
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: invector(nblock,www%allbotconfigs(iproc):www%alltopconfigs(jproc))
  DATATYPE,intent(out) :: outvector(nblock,www%sss%allbotspins(iproc):www%sss%alltopspins(jproc))
  DATATYPE :: smallvect(nblock,www%sss%maxspinsetsize), &
       smalltemp(nblock,www%sss%maxspinsetsize)  !! AUTOMATIC
  integer :: iset, iind,ii,pp

  if (iproc.lt.www%startrank.or.jproc.gt.www%endrank) then
     OFLWR "STARTRANK ERR tRANSTO GEN ",iproc,www%startrank,www%endrank; CFLST
  endif

  if (www%sss%alltopspins(jproc).ge.www%sss%allbotspins(iproc)) then
     outvector(:,:)=0d0
  endif

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
        outvector(:,iind:iind+www%sss%spinsetrank(iset,pp)-1) = &
             smalltemp(:,1:www%sss%spinsetrank(iset,pp))
        iind=iind+www%sss%spinsetrank(iset,pp)
     enddo
     
     if (iind.ne.www%sss%alltopspins(pp)+1) then
        OFLWR "IIND ERROto", iind,pp,www%sss%allbotspins(pp),www%sss%alltopspins(pp); CFLST
     endif
  enddo

end subroutine configspin_transformto_general



subroutine dfspin_transformto_general(www,nblock,invector,outvector,iproc,jproc)
  use fileptrmod
  use walkmod
  implicit none
  integer,intent(in) :: nblock,iproc,jproc
  type(walktype),intent(in) :: www
  integer :: iset, iind,ii,jset,pp
  DATATYPE,intent(in) :: invector(nblock,www%allbotdfconfigs(iproc):www%alltopdfconfigs(jproc))
  DATATYPE,intent(out) :: outvector(nblock,www%sss%allbotspindfs(iproc):www%sss%alltopspindfs(jproc))
  DATATYPE :: smallvect(nblock,www%sss%maxspinsetsize), &
       smalltemp(nblock,www%sss%maxspinsetsize)    !! AUTOMATIC

  if (iproc.lt.www%startrank.or.jproc.gt.www%endrank) then
     OFLWR "STARTRANK ERR DFTRANSTO GEN ",iproc,www%startrank,www%endrank; CFLST
  endif

  if (www%sss%alltopspindfs(jproc).ge.www%sss%allbotspindfs(iproc)) then
     outvector(:,:)=0d0
  endif

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
        outvector(:,iind:iind+www%sss%spinsetrank(iset,pp)-1) = &
             smalltemp(:,1:www%sss%spinsetrank(iset,pp))
        iind=iind+www%sss%spinsetrank(iset,pp)
     enddo

     if (iind.ne.www%sss%alltopspindfs(pp)+1) then
        OFLWR "IIND ERROdfto", iind,pp,www%sss%allbotspindfs(pp),www%sss%alltopspindfs(pp); CFLST
     endif
  end do

end subroutine dfspin_transformto_general



subroutine configspin_transformfrom_general(www,nblock,invector,outvector,iproc,jproc)
  use fileptrmod
  use walkmod
  implicit none
  integer,intent(in) :: nblock,iproc,jproc
  type(walktype),intent(in) :: www
  integer :: iset, iind, ii, pp
  DATATYPE,intent(out) :: outvector(nblock,www%allbotconfigs(iproc):www%alltopconfigs(jproc))
  DATATYPE,intent(in) :: invector(nblock,www%sss%allbotspins(iproc):www%sss%alltopspins(jproc))
  DATATYPE :: smallvect(nblock,www%sss%maxspinsetsize), &
       smalltemp(nblock,www%sss%maxspinsetsize)  !! AUTOMATIC

  if (iproc.lt.www%startrank.or.jproc.gt.www%endrank) then
     OFLWR "STARTRANK ERR TRANSFROM GEN ",iproc,www%startrank,www%endrank; CFLST
  endif

  if (www%alltopconfigs(jproc).ge.www%allbotconfigs(iproc)) then
     outvector(:,:)=0d0
  endif

  do pp=iproc,jproc

     iind=www%sss%allbotspins(pp)

     do iset=1,www%sss%numspinsets(pp)
        smallvect(:,1:www%sss%spinsetrank(iset,pp))=&
             invector(:,iind:iind+www%sss%spinsetrank(iset,pp)-1)
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



subroutine dfspin_transformfrom_general(www,nblock,invector,outvector,iproc,jproc)
  use fileptrmod
  use walkmod
  implicit none
  integer,intent(in) :: nblock,iproc,jproc
  type(walktype),intent(in) :: www
  integer :: iset, iind, ii, jset, pp
  DATATYPE,intent(out) :: outvector(nblock,www%allbotdfconfigs(iproc):www%alltopdfconfigs(jproc))
  DATATYPE,intent(in) :: invector(nblock,www%sss%allbotspindfs(iproc):www%sss%alltopspindfs(jproc))
  DATATYPE :: smallvect(nblock,www%sss%maxspinsetsize), &
       smalltemp(nblock,www%sss%maxspinsetsize)  !! AUTOMATIC

  if (iproc.lt.www%startrank.or.jproc.gt.www%endrank) then
     OFLWR "STARTRANK ERR DFTRANSfrom GEN ",iproc,www%startrank,www%endrank; CFLST
  endif

  if (www%alltopdfconfigs(jproc).ge.www%allbotdfconfigs(iproc)) then
     outvector(:,:)=0d0
  endif

  do pp=iproc,jproc

     iind=www%sss%allbotspindfs(pp)
     do jset=1,www%sss%numspindfsets(pp)
        iset=www%sss%spindfsetindex(jset,pp)
        smallvect(:,1:www%sss%spinsetrank(iset,pp))=&
             invector(:,iind:iind+www%sss%spinsetrank(iset,pp)-1)
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



subroutine df_project(www,howmany,avector)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: howmany
  DATATYPE,intent(inout) :: avector(howmany,www%firstconfig:www%lastconfig)

  if (www%parconsplit.eq.0) then
     call df_project_all(www,howmany,avector)
  else
     call df_project_local(www,howmany,avector)
  endif

end subroutine df_project


subroutine df_project_all(www,howmany,avector)
  use walkmod
  implicit none
  integer,intent(in) :: howmany
  type(walktype),intent(in) :: www
  DATATYPE,intent(inout) :: avector(howmany,www%numconfig)
  integer :: i
  do i=1,www%numconfig
     avector(:,i)=avector(:,i)*www%ddd%dfincludedmask(i)
  enddo
end subroutine df_project_all


subroutine df_project_local(www,howmany,avector)
  use walkmod
  implicit none
  integer,intent(in) :: howmany
  type(walktype),intent(in) :: www
  DATATYPE,intent(inout) :: avector(howmany,www%botconfig:www%topconfig)
  integer :: i
  do i=www%botconfig,www%topconfig
     avector(:,i)=avector(:,i)*www%ddd%dfincludedmask(i)
  enddo
end subroutine df_project_local


subroutine df_transformto_general(www,howmany,avectorin,avectorout,iproc,jproc)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: howmany,iproc,jproc
  DATATYPE,intent(in) :: avectorin(howmany,www%allbotconfigs(iproc):www%alltopconfigs(jproc))
  DATATYPE,intent(out) :: avectorout(howmany,www%allbotdfconfigs(iproc):www%alltopdfconfigs(jproc))
  integer :: i

  if (www%alltopdfconfigs(jproc).ge.www%allbotdfconfigs(iproc)) then
     avectorout(:,:)=0d0
  endif
  do i=www%allbotdfconfigs(iproc),www%alltopdfconfigs(jproc)
     avectorout(:,i)=avectorin(:,www%ddd%dfincludedconfigs(i))
  enddo

end subroutine df_transformto_general



subroutine df_transformfrom_general(www,howmany,avectorin,avectorout,iproc,jproc)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: howmany,iproc,jproc
  DATATYPE,intent(in) :: avectorin(howmany,www%allbotdfconfigs(iproc):www%alltopdfconfigs(jproc))
  DATATYPE,intent(out) :: avectorout(howmany,www%allbotconfigs(iproc):www%alltopconfigs(jproc))
  integer :: i

  if (www%alltopconfigs(jproc).ge.www%allbotconfigs(iproc)) then
     avectorout(:,:)=0d0
  endif
  do i=www%allbotdfconfigs(iproc),www%alltopdfconfigs(jproc)
     avectorout(:,www%ddd%dfincludedconfigs(i))=avectorin(:,i)
  enddo

end subroutine df_transformfrom_general


subroutine basis_project(www,howmany,avector)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: howmany
  DATATYPE,intent(inout) :: avector(howmany,www%firstconfig:www%lastconfig)

  if (www%allspinproject==1) then
     call configspin_project(www,howmany,avector)
  endif
  if (www%dfrestrictflag.ne.www%dflevel) then
     call df_project(www,howmany,avector)
  endif

end subroutine basis_project



subroutine basis_transformto_all(www,howmany,avectorin,avectorout)
  use walkmod
  use mpimod   !! myrank,nprocs ( not perfect )
  use mpisubmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,www%numconfig)
  DATATYPE,intent(out) :: avectorout(howmany,www%numdfbasis)

  avectorout(:,:)=0d0

  if (www%sparseconfigflag.ne.0) then
     if (www%topdfbasis.ge.www%botdfbasis) then
        call basis_transformto_local(www,howmany,avectorin(:,www%botconfig:www%topconfig),&
             avectorout(:,www%botdfbasis:www%topdfbasis))
     endif
     call mpiallgather(avectorout,www%numdfbasis*howmany,www%dfbasisperproc(:)*howmany,&
          www%maxdfbasisperproc*howmany)
  else
     call basis_transformto_general(www,howmany,avectorin,avectorout,1,nprocs)
  endif
end subroutine basis_transformto_all


subroutine basis_transformto_local(www,howmany,avectorin,avectorout)
  use walkmod
  use mpimod   !! myrank
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,www%botconfig:www%topconfig)
  DATATYPE,intent(out) :: avectorout(howmany,www%botdfbasis:www%topdfbasis)

  call basis_transformto_general(www,howmany,avectorin,avectorout,myrank,myrank)

end subroutine basis_transformto_local


subroutine basis_transformto_general(www,howmany,avectorin,avectorout,iproc,jproc)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: howmany,iproc,jproc
  DATATYPE,intent(in) :: avectorin(howmany,www%allbotconfigs(iproc):www%alltopconfigs(jproc))
  DATATYPE,intent(out) :: avectorout(howmany,www%allbotdfbasis(iproc):www%alltopdfbasis(jproc))
  DATATYPE :: workvec(howmany,www%allbotdfconfigs(iproc):www%alltopdfconfigs(jproc))   !! AUTOMATIC

  if (www%alltopdfbasis(jproc)-www%allbotdfbasis(iproc)+1.ne.0) then
     avectorout(:,:)=0d0
     if (www%dfrestrictflag.ne.www%dflevel) then
        if (www%allspinproject.ne.0) then
           call df_transformto_general(www,howmany,avectorin,workvec,iproc,jproc)
           call dfspin_transformto_general(www,howmany,workvec,avectorout,iproc,jproc)
        else
           call df_transformto_general(www,howmany,avectorin,avectorout,iproc,jproc)
        endif
     else
        if (www%allspinproject.ne.0) then
           call configspin_transformto_general(www,howmany,avectorin,avectorout,iproc,jproc)
        else
           if (www%alltopdfbasis(jproc).ge.www%allbotdfbasis(iproc)) then
              avectorout(:,:)=avectorin(:,:)
           endif
        endif
     endif
  endif

end subroutine basis_transformto_general



subroutine basis_transformfrom_all(www,howmany,avectorin,avectorout)
  use walkmod
  use mpimod    !! myrank,nprocs ( not perfect )
  use mpisubmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,www%numdfbasis)
  DATATYPE,intent(out) :: avectorout(howmany,www%numconfig)

  avectorout=0d0

  if (www%sparseconfigflag.ne.0) then
     if (www%topconfig.ge.www%botconfig) then
        call basis_transformfrom_local(www,howmany,avectorin(:,www%botdfbasis:www%topdfbasis),&
             avectorout(:,www%botconfig:www%topconfig))
     endif
     call mpiallgather(avectorout,www%numconfig*howmany,www%configsperproc(:)*howmany,&
          www%maxconfigsperproc*howmany)
  else
     call basis_transformfrom_general(www,howmany,avectorin,avectorout,1,nprocs)
  endif

end subroutine basis_transformfrom_all



subroutine basis_transformfrom_local(www,howmany,avectorin,avectorout)
  use walkmod
  use mpimod   !! myrank
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,www%botdfbasis:www%topdfbasis)
  DATATYPE,intent(out) :: avectorout(howmany,www%botconfig:www%topconfig)

  call basis_transformfrom_general(www,howmany,avectorin,avectorout,myrank,myrank)

end subroutine basis_transformfrom_local


subroutine basis_transformfrom_general(www,howmany,avectorin,avectorout,iproc,jproc)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: howmany,iproc,jproc
  DATATYPE,intent(in) :: avectorin(howmany,www%allbotdfbasis(iproc):www%alltopdfbasis(jproc))
  DATATYPE,intent(out) :: avectorout(howmany,www%allbotconfigs(iproc):www%alltopconfigs(jproc))
  DATATYPE :: workvec(howmany,www%allbotdfconfigs(iproc):www%alltopdfconfigs(jproc))   !! AUTOMATIC

  if (www%alltopconfigs(jproc)-www%allbotconfigs(iproc)+1.ne.0) then
     avectorout(:,:)=0d0
  endif

  if (www%alltopdfbasis(jproc)-www%allbotdfbasis(iproc)+1.ne.0) then
     if (www%dfrestrictflag.ne.www%dflevel) then
        if (www%allspinproject.ne.0) then
           call dfspin_transformfrom_general(www,howmany,avectorin,workvec,iproc,jproc)
           call df_transformfrom_general(www,howmany,workvec,avectorout,iproc,jproc)
        else
           call df_transformfrom_general(www,howmany,avectorin,avectorout,iproc,jproc)
        endif
     else
        if (www%allspinproject.ne.0) then
           call configspin_transformfrom_general(www,howmany,avectorin,avectorout,iproc,jproc)
        else
           if (www%alltopconfigs(jproc).ge.www%allbotconfigs(iproc)) then
              avectorout(:,:)=avectorin(:,:)
           endif
        endif
     endif
  endif

end subroutine basis_transformfrom_general


subroutine basis_shuffle(howmany,wwin,avectorin,wwout,avectorout)
  use timing_parameters
  use fileptrmod
  use walkmod
  use mpimod     !! nprocs
  use mpisubmod
  use clockmod
  implicit none
  type(walktype),intent(in) :: wwin,wwout
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,wwin%botdfbasis:wwin%topdfbasis)
  DATATYPE,intent(out) :: avectorout(howmany,wwout%botdfbasis:wwout%topdfbasis)
#ifndef MPIFLAG
  avectorout(:,:)=avectorin(:,:)
#else
  integer,allocatable :: pairs(:,:),pairlow(:),pairhigh(:),pairsize(:),&
       pairtag(:)
  integer :: iin, iout, numpairs, ipair, getlen, myiostat
  integer,save :: times(20)=0, atime = -99, btime = -99, icalled=0, &
       numcalled=0

  if (wwin%numdfbasis.ne.wwout%numdfbasis) then
     OFLWR "error, on shuffle basis sizes do not agree ",&
          wwin%numdfbasis,wwout%numdfbasis; CFLST
  endif

  if (nprocs.eq.1) then
     avectorout(:,:)=avectorin(:,:)
     return    !! RETURN
  endif

  if (icalled.eq.0) then
     call myclock(btime); atime=btime
     if (myrank.eq.1.and.notiming.le.1) then
        open(1777, file=timingdir(1:getlen(timingdir))//"/shuffletime.dat", &
             status="unknown",iostat=myiostat)
        call checkiostat(myiostat," opening shuffletime.dat")
        write(1777,'(T16,100A15)')  &
             "Init ", &     !! (1)
             "Alloc ", &    !! (2)
             "Pairs ", &    !! (3)
             "Sendrecv", &  !! (4)
             "Final", &     !! (5)
             "Away"         !! (6)
        close(1777)
     endif
  else
     call myclock(btime); times(6)=times(6)+btime-atime; atime=btime
  endif
  icalled=1; numcalled=numcalled+1

  if (wwout%topdfbasis.ge.wwout%botdfbasis) then
     avectorout(:,:)=0d0
  endif

  numpairs=0
  do iin=1,nprocs
     do iout=1,nprocs
        if (wwout%allbotdfbasis(iout).le.wwin%alltopdfbasis(iin).and.&
             wwout%alltopdfbasis(iout).ge.wwin%allbotdfbasis(iin).and.&
             wwin%alltopdfbasis(iin).ge.wwin%allbotdfbasis(iin).and.&
             wwout%alltopdfbasis(iout).ge.wwout%allbotdfbasis(iout)) then
           numpairs=numpairs+1
        endif
     enddo
  enddo

  call myclock(btime); times(1)=times(1)+ btime-atime; atime=btime

  allocate(pairs(2,numpairs),pairlow(numpairs),pairhigh(numpairs),&
       pairsize(numpairs),pairtag(numpairs))
  pairs=(-1); pairlow=(-1); pairhigh=(-1); pairsize=(-1); pairtag=(-1)

  call myclock(btime); times(2)=times(2)+ btime-atime; atime=btime

  numpairs=0
  do iin=1,nprocs
     do iout=1,nprocs
        if (wwout%allbotdfbasis(iout).le.wwin%alltopdfbasis(iin).and.&
             wwout%alltopdfbasis(iout).ge.wwin%allbotdfbasis(iin).and.&
             wwin%alltopdfbasis(iin).ge.wwin%allbotdfbasis(iin).and.&
             wwout%alltopdfbasis(iout).ge.wwout%allbotdfbasis(iout)) then
           numpairs=numpairs+1
           pairs(1,numpairs)=iin;     pairs(2,numpairs)=iout
           pairlow(numpairs)=max(wwout%allbotdfbasis(iout),wwin%allbotdfbasis(iin))
           pairhigh(numpairs)=min(wwout%alltopdfbasis(iout),wwin%alltopdfbasis(iin))
           pairsize(numpairs)=pairhigh(numpairs)-pairlow(numpairs)+1
           pairtag(numpairs)=iin+(iout-1)*nprocs
           if (pairsize(numpairs).lt.1) then
              print *, "progcheckkkk pairsize",pairsize(numpairs),&
                   wwin%allbotdfbasis(iin),wwin%alltopdfbasis(iin),&
                   wwout%allbotdfbasis(iout),wwout%alltopdfbasis(iout)
              stop
           endif
        endif
     enddo
  enddo

!!$  these things will happen and they are ok.
!!$  if (numpairs.eq.nprocs) then
!!$     OFLWR "What? it doesn't look like you need to shuffle."
!!$     write(mpifileptr,'(2I10)') pairs(:,:);     CFLST
!!$  elseif (numpairs.lt.nprocs) then
!!$     OFLWR "What? too few pairs",numpairs; CFLST
!!$  endif

  if (numpairs.lt.max(wwin%nzprocs,wwout%nzprocs)) then
     OFLWR "What? too few pairs",numpairs,wwin%nzprocs,wwout%nzprocs; CFLST
  endif
  call myclock(btime); times(3)=times(3)+ btime-atime; atime=btime

  do ipair=1,numpairs
     if (pairs(1,ipair).eq.myrank.and.pairs(2,ipair).eq.myrank) then
        avectorout(:,pairlow(ipair):pairhigh(ipair)) = &
             avectorin(:,pairlow(ipair):pairhigh(ipair))
     elseif (pairs(1,ipair).eq.myrank.and.pairs(2,ipair).ne.myrank) then
        call mympisend(avectorin(:,pairlow(ipair):pairhigh(ipair)),pairs(2,ipair),&
             pairtag(ipair), pairsize(ipair)*howmany)
     elseif (pairs(2,ipair).eq.myrank.and.pairs(1,ipair).ne.myrank) then
        call mympirecv(avectorout(:,pairlow(ipair):pairhigh(ipair)),pairs(1,ipair),&
             pairtag(ipair), pairsize(ipair)*howmany)
     endif
  enddo

  call myclock(btime); times(4)=times(4)+ btime-atime; atime=btime

  deallocate(pairs,pairlow,pairhigh,pairsize,pairtag)

  if (myrank.eq.1.and.notiming.le.1.and.mod(numcalled,timingout).eq.0) then
     open(1777, file=timingdir(1:getlen(timingdir))//"/shuffletime.dat", &
          status="old", position="append",iostat=myiostat)
     call checkiostat(myiostat," opening shuffletime.dat")
     write(1777,'(100I15)',iostat=myiostat) times(1:6)/1000
     call checkiostat(myiostat," writing shuffletime.dat")
     close(1777)
  endif

  call myclock(btime); times(5)=times(5)+ btime-atime; atime=btime

#endif

end subroutine basis_shuffle


subroutine basis_shuffle_several(howmany,wwin,avectorin,wwout,avectorout,numvec)
  use timing_parameters
  use fileptrmod
  use walkmod
  use mpimod     !! nprocs
  implicit none
  type(walktype),intent(in) :: wwin,wwout
  integer,intent(in) :: howmany,numvec
  DATATYPE,intent(in) :: avectorin(howmany,wwin%botdfbasis:wwin%topdfbasis,numvec)
  DATATYPE,intent(out) :: avectorout(howmany,wwout%botdfbasis:wwout%topdfbasis,numvec)
#ifndef MPIFLAG
  avectorout(:,:,:)=avectorin(:,:,:)
#else
  DATATYPE :: nullvector1(howmany),nullvector2(howmany)
  integer :: ii

  nullvector1=0; nullvector2=0
  do ii=1,numvec
     if (wwin%topdfbasis.ge.wwin%botdfbasis.and.wwout%topdfbasis.ge.wwout%botdfbasis) then
        call basis_shuffle(howmany,wwin,avectorin(:,:,ii),wwout,avectorout(:,:,ii))
     elseif (wwin%topdfbasis.ge.wwin%botdfbasis) then
        call basis_shuffle(howmany,wwin,avectorin(:,:,ii),wwout,nullvector2)
     elseif (wwout%topdfbasis.ge.wwout%botdfbasis) then
        call basis_shuffle(howmany,wwin,nullvector1,wwout,avectorout(:,:,ii))
     else
        call basis_shuffle(howmany,wwin,nullvector1,wwout,nullvector2)
     endif
  enddo

#endif

end subroutine basis_shuffle_several

end module basissubmod
