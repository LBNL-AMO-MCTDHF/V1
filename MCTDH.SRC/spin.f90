
!! SUBROUTINES FOR BASIS SET TRANSFORMATIONS AND PROJECTIONS (spin and restricted configuration space)

#include "Definitions.INC"


!! INDICATES WHETHER A GIVEN CONFIG is included in the restricted configuration list.

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


subroutine basis_set(www)
  use r_parameters
  use walkmod
  use mpimod !! nprocs
  implicit none
  type(walktype),intent(inout) :: www

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
end subroutine basis_set


! "included" configurations are those kept in the calculation.
! there are numdfconfigs INCLUDED configurations.

subroutine init_dfcon(www)
  use fileptrmod
  use walkmod
  use mpimod
  implicit none
  type(walktype),intent(inout) :: www
  integer :: i,j,iconfig,nondfconfigs,dfrank,ii
  logical :: allowedconfig0
  integer,allocatable :: dfnotconfigs(:)

  allocate(www%ddd%dfincludedmask(www%numconfig), www%ddd%dfincludedconfigs(www%numconfig), &
       www%ddd%dfincludedindex(www%numconfig),&
       dfnotconfigs(www%numconfig),  www%dfconfsperproc(nprocs), &
       www%allbotdfconfigs(nprocs),www%alltopdfconfigs(nprocs))

  if (www%dfrestrictflag.lt.www%dflevel) then
     OFLWR "error, set dfrestrictflag .ge. dflevel",www%dfrestrictflag,www%dflevel; CFLST
  endif

  if (www%dfrestrictflag.eq.www%dflevel) then
     www%ddd%dfincludedmask(:)=1; dfnotconfigs(:)=(-1)
     do i=1,www%numconfig
        www%ddd%dfincludedconfigs(i)=i
        www%ddd%dfincludedindex(i)=i
     enddo
     www%numdfconfigs=www%numconfig; nondfconfigs=0
     www%dfconfsperproc(:)=www%configsperproc(:)
     www%allbotdfconfigs(:)=www%allbotconfigs(:)
     www%alltopdfconfigs(:)=www%alltopconfigs(:)
     www%maxdfconfsperproc=www%maxconfigsperproc
     www%botdfconfig=www%botconfig;  www%topdfconfig=www%topconfig

  else

     www%ddd%dfincludedmask(:)=0;  www%ddd%dfincludedconfigs(:)=(-1);  dfnotconfigs(:)=(-1)
     www%ddd%dfincludedindex(:)=(-999)

     www%numdfconfigs=0;  nondfconfigs=0

     do i=1,www%numconfig
        if (allowedconfig0(www,www%configlist(:,i),www%dfrestrictflag)) then 
           www%numdfconfigs=www%numdfconfigs+1
           www%ddd%dfincludedmask(i)=1;        www%ddd%dfincludedconfigs(www%numdfconfigs)=i
           www%ddd%dfincludedindex(i)=www%numdfconfigs
        else
           nondfconfigs=nondfconfigs+1;        dfnotconfigs(nondfconfigs)=i
        endif
     enddo

     if (www%numdfconfigs.eq.0) then
        OFLWR "Error, no included DF configs!"; CFLST
     endif

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
              if (allowedconfig0(www,www%configlist(:,www%singlewalk(j,iconfig)),www%dfrestrictflag)) then
                 www%ddd%numdfwalks=www%ddd%numdfwalks+1
              endif
           enddo
        endif
     enddo

     OFLWR "numdfwalks:  ", www%ddd%numdfwalks," is total number DF single excitations on this processor"; WRFL; CFL

     allocate(www%ddd%dfwalkfrom(www%ddd%numdfwalks), www%ddd%dfwalkto(www%ddd%numdfwalks),&
          www%ddd%includedorb(www%ddd%numdfwalks), www%ddd%excludedorb(www%ddd%numdfwalks), &
          www%ddd%dfwalkphase(www%ddd%numdfwalks))

     www%ddd%numdfwalks=0
     do i=1,NONdfconfigs
        iconfig=dfNOTconfigs(i)
        if (iconfig.ge.www%configstart.and.iconfig.le.www%configend) then
           do j=1,www%numsinglewalks(iconfig)
              if (allowedconfig0(www,www%configlist(:,www%singlewalk(j,iconfig)),www%dfrestrictflag)) then
                 www%ddd%numdfwalks=www%ddd%numdfwalks+1
              
                 www%ddd%dfwalkTO(www%ddd%numdfwalks)=iconfig
                 www%ddd%dfwalkFROM(www%ddd%numdfwalks)=www%singlewalk(j,iconfig)
                 www%ddd%EXcludedorb(www%ddd%numdfwalks)=www%singlewalkopspf(1,j,iconfig)
                 www%ddd%INcludedorb(www%ddd%numdfwalks)=www%singlewalkopspf(2,j,iconfig)
              
                 www%ddd%dfwalkphase(www%ddd%numdfwalks)=www%singlewalkdirphase(j,iconfig)
              endif
           enddo
        endif
     enddo

     OFLWR "numdfwalks:  ", www%ddd%numdfwalks," is really number DF single excitations on this processor"; WRFL; CFL
     
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



subroutine configspin_project(www,nr, vector)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: nr
  DATATYPE,intent(inout) :: vector(nr,www%firstconfig:www%lastconfig)

  call configspin_project_general(www,nr,vector(:,www%configstart),www%startrank,www%endrank)

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
       outvector(nr,www%allbotconfigs(iproc):www%alltopconfigs(jproc))      !! AUTOMATIC
  integer :: iset, ii, pp

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



subroutine configspin_transformto_general(www,nblock,invector,outvector,iproc,jproc)
  use fileptrmod
  use walkmod
  implicit none
  integer,intent(in) :: nblock,iproc,jproc
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: invector(nblock,www%allbotconfigs(iproc):www%alltopconfigs(jproc))
  DATATYPE,intent(out) :: outvector(nblock,www%sss%allbotspins(iproc):www%sss%alltopspins(jproc))
  DATATYPE :: smallvect(nblock,www%sss%maxspinsetsize), smalltemp(nblock,www%sss%maxspinsetsize)  !! AUTOMATIC
  integer :: iset, iind,ii,pp

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



subroutine dfspin_transformto_general(www,nblock,invector,outvector,iproc,jproc)
  use fileptrmod
  use walkmod
  implicit none
  integer,intent(in) :: nblock,iproc,jproc
  type(walktype),intent(in) :: www
  integer :: iset, iind,ii,jset,pp
  DATATYPE,intent(in) :: invector(nblock,www%allbotdfconfigs(iproc):www%alltopdfconfigs(jproc))
  DATATYPE,intent(out) :: outvector(nblock,www%sss%allbotspindfs(iproc):www%sss%alltopspindfs(jproc))
  DATATYPE :: smallvect(nblock,www%sss%maxspinsetsize), smalltemp(nblock,www%sss%maxspinsetsize)    !! AUTOMATIC

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



subroutine configspin_transformfrom_general(www,nblock,invector,outvector,iproc,jproc)
  use fileptrmod
  use walkmod
  implicit none
  integer,intent(in) :: nblock,iproc,jproc
  type(walktype),intent(in) :: www
  integer :: iset, iind, ii, pp
  DATATYPE,intent(out) :: outvector(nblock,www%allbotconfigs(iproc):www%alltopconfigs(jproc))
  DATATYPE,intent(in) :: invector(nblock,www%sss%allbotspins(iproc):www%sss%alltopspins(jproc))
  DATATYPE :: smallvect(nblock,www%sss%maxspinsetsize), smalltemp(nblock,www%sss%maxspinsetsize)  !! AUTOMATIC

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



subroutine dfspin_transformfrom_general(www,nblock,invector,outvector,iproc,jproc)
  use fileptrmod
  use walkmod
  implicit none
  integer,intent(in) :: nblock,iproc,jproc
  type(walktype),intent(in) :: www
  integer :: iset, iind, ii, jset, pp
  DATATYPE,intent(out) :: outvector(nblock,www%allbotdfconfigs(iproc):www%alltopdfconfigs(jproc))
  DATATYPE,intent(in) :: invector(nblock,www%sss%allbotspindfs(iproc):www%sss%alltopspindfs(jproc))
  DATATYPE :: smallvect(nblock,www%sss%maxspinsetsize), smalltemp(nblock,www%sss%maxspinsetsize)  !! AUTOMATIC

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

  avectorout(:,:)=0d0

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
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,www%numconfig)
  DATATYPE,intent(out) :: avectorout(howmany,www%numdfbasis)

  if (www%sparseconfigflag.ne.0) then
     call basis_transformto_local(www,howmany,avectorin(:,www%botconfig),avectorout(:,www%botdfbasis))
     call mpiallgather(avectorout(:,:),www%numdfbasis*howmany,www%dfbasisperproc(:)*howmany,&
          www%maxdfbasisperproc*howmany)
  else
     call basis_transformto_general(www,howmany,avectorin(:,:),avectorout(:,:),1,nprocs)
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
     if (www%dfrestrictflag.ne.www%dflevel) then
        if (www%allspinproject.ne.0) then
           call df_transformto_general(www,howmany,avectorin(:,:),workvec(:,:),iproc,jproc)
           call dfspin_transformto_general(www,howmany,workvec(:,:),avectorout(:,:),iproc,jproc)
        else
           call df_transformto_general(www,howmany,avectorin(:,:),avectorout(:,:),iproc,jproc)
        endif
     else
        if (www%allspinproject.ne.0) then
           call configspin_transformto_general(www,howmany,avectorin(:,:),avectorout(:,:),iproc,jproc)
        else
           avectorout(:,:)=avectorin(:,:)
        endif
     endif
  endif

end subroutine basis_transformto_general



subroutine basis_transformfrom_all(www,howmany,avectorin,avectorout)
  use walkmod
  use mpimod    !! myrank,nprocs ( not perfect )
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,www%numdfbasis)
  DATATYPE,intent(out) :: avectorout(howmany,www%numconfig)

  if (www%sparseconfigflag.ne.0) then
     call basis_transformfrom_local(www,howmany,avectorin(:,www%botdfbasis),avectorout(:,www%botconfig))
     call mpiallgather(avectorout(:,:),www%numconfig*howmany,www%configsperproc(:)*howmany,&
          www%maxconfigsperproc*howmany)
  else
     call basis_transformfrom_general(www,howmany,avectorin(:,:),avectorout(:,:),1,nprocs)
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
           call dfspin_transformfrom_general(www,howmany,avectorin(:,:),workvec(:,:),iproc,jproc)
           call df_transformfrom_general(www,howmany,workvec(:,:),avectorout(:,:),iproc,jproc)
        else
           call df_transformfrom_general(www,howmany,avectorin(:,:),avectorout(:,:),iproc,jproc)
        endif
     else
        if (www%allspinproject.ne.0) then
           call configspin_transformfrom_general(www,howmany,avectorin(:,:),avectorout(:,:),iproc,jproc)
        else
           avectorout(:,:)=avectorin(:,:)
        endif
     endif
  endif

end subroutine basis_transformfrom_general




subroutine fullbasis_transformto_local(www,howmany,avectorin,avectorout)
  use walkmod
  use mpimod    !! myrank
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,www%botconfig:www%topconfig)
  DATATYPE,intent(out) :: avectorout(howmany,www%botbasis:www%topbasis)

  if (www%topbasis-www%botbasis+1.ne.0) then
     if (www%allspinproject.ne.0) then
        call configspin_transformto_general(www,howmany,avectorin(:,:),avectorout(:,:),myrank,myrank)
     else
        avectorout(:,:)=avectorin(:,:)
        endif
  endif

end subroutine fullbasis_transformto_local


subroutine fullbasis_transformfrom_local(www,howmany,avectorin,avectorout)
  use walkmod
  use mpimod    !! myrank
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,www%botbasis:www%topbasis)
  DATATYPE,intent(out) :: avectorout(howmany,www%botconfig:www%topconfig)

  if (www%topconfig-www%botconfig+1.ne.0) then
     avectorout(:,:)=0d0
  endif

  if (www%topbasis-www%botbasis+1.ne.0) then
     if (www%allspinproject.ne.0) then
        call configspin_transformfrom_general(www,howmany,avectorin(:,:),avectorout(:,:),myrank,myrank)
     else
        avectorout(:,:)=avectorin(:,:)
     endif
  endif

end subroutine fullbasis_transformfrom_local


