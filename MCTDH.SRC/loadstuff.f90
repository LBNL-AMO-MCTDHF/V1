
#include "Definitions.INC"

subroutine save_vector(psi,afile,sfile)
  use parameters
  use mpimod
  implicit none
  DATATYPE :: psi(psilength) 
  character :: afile*(*), sfile*(*)
  integer :: iprop,qqblocks(nprocs),qqend(nprocs),qqstart(nprocs),ii,ispf
  DATATYPE, allocatable :: parorbitals(:,:,:)

  if (myrank.ne.1.and.parorbsplit.ne.3) then
     OFLWR "   --> rank 1 will save. "; CFL;     return
  endif
  if (parorbsplit.eq.3) then
     if (myrank.eq.1) then
        OFLWR "allocating big orbitals for save!!"; CFL
        allocate(parorbitals(spfsize,nprocs,nspf)); parorbitals(:,:,:)=0d0
        OFLWR "   ...ok"; CFL
     else
        allocate(parorbitals(1,1,nspf))
     endif
     qqblocks(:)=spfsize
     do ii=1,nprocs
        qqend(ii)=ii*spfsize; qqstart(ii)=(ii-1)*spfsize+1
     enddo

     do ispf=1,nspf
#ifdef REALGO
        call mygatherv_real(psi(spfstart+(ispf-1)*spfsize),parorbitals(:,:,ispf),&
             spfsize*nprocs,qqstart(myrank),qqend(myrank),qqblocks(:),qqstart(:),.false.)
#else
        call mygatherv_complex(psi(spfstart+(ispf-1)*spfsize),parorbitals(:,:,ispf),&
             spfsize*nprocs,qqstart(myrank),qqend(myrank),qqblocks(:),qqstart(:),.false.)
#endif
     enddo
     if (myrank.ne.1) then
        OFLWR " ---> orbs sent to rank 1 for save!"; CFL
        deallocatE(parorbitals)
        return
     endif
  endif

  open(999,file=afile, status="unknown", form="unformatted")
  call avector_header_write(999,mcscfnum)

  open(998,file=sfile, status="unknown", form="unformatted")
  call spf_header_write(998,nspf+numfrozen)

  if (parorbsplit.eq.3) then
     call write_spf(998,parorbitals(:,:,:),spfsize*nprocs)
     deallocate(parorbitals)
  else
     call write_spf(998,psi(spfstart),spfsize)
  endif
  close(998)

  do iprop=1,mcscfnum
     call write_avector(999,psi(astart(iprop)))
  enddo
  close(999)


  OFLWR "rank 1 SAVED vectors!"; CFL 

end subroutine save_vector




subroutine write_spf(unit,spfin,inspfsize)
  use parameters
  use opmod !! frozenspfs
  implicit none
  integer :: unit,inspfsize
  DATATYPE :: allspf(inspfsize,nspf+numfrozen), spfin(inspfsize, nspf   )
  if (numfrozen.gt.0) then
     if (inspfsize.ne.spfsize) then
        OFLWR "AAACK SFFFRROZ"; CFLST
     endif
     allspf(:,1:numfrozen)=frozenspfs(:,:)
  endif
  allspf(:,numfrozen+1:nspf+numfrozen) = spfin(:,:)
  write (unit) allspf
end subroutine

subroutine spf_header_write(iunit,num)
  use mpimod
  use parameters
  implicit none
  integer ::  cflag,iunit,num,qqq(3)
#ifdef REALGO
  cflag=0
#else
  cflag=1
#endif
  qqq(:)=1
  if (parorbsplit.eq.3) then
     qqq(3)=nprocs
  endif
  write(iunit) spfdims(:)*qqq(:),num,cflag
end subroutine
subroutine spf_header_read(iunit,outdims,outnspf,cflag)
!! use fileptrmod  !! need nprocs and parorbsplit now
  use mpimod
  use parameters
  implicit none
  integer ::  cflag,outdims(3),outnspf,iunit,qdim
#ifdef REALGO
  cflag=0
#else
  cflag=1
#endif
  if (parorbsplit.eq.3) then
     qdim=nprocs
  else
     qdim=1
  endif

  read(iunit) outdims(:),outnspf,cflag

!!$  if ((outdims(3)/qdim)*qdim.ne.outdims(3)) then
!!$     OFLWR "parorb err header read",outdims(3),qdim; CFLST
!!$  endif
!!$  outdims(3)=outdims(3)/qdim
end subroutine

!subroutine loa d_spfs(inspfs, numloaded)    !! both output
!  use parameters
!  implicit none
!  integer :: numloaded
!  DATATYPE :: inspfs(spfsize,nspf)
!  call load_sp fs0(inspfs,spfdims,nspf,spfdimtype,spffile,numloaded)
!  OFLWR "Loaded ", numloaded, " spfs from file ",spffile; CFL
!end subroutine lo ad_spfs


subroutine load_spfs(inspfs, numloaded)    !! both output
  use parameters
  use mpimod
  implicit none
  integer :: numloaded,ii,jj,kk,ll,flag
  DATATYPE :: inspfs(spfsize,nspf+numfrozen), inspfs0(spfsize,nspf+numfrozen+numskiporbs)

  numloaded=0
  do ii=1,numspffiles

     call load_spfs0(inspfs0(:,numloaded+1),spfdims,nspf+numfrozen+numskiporbs-numloaded,spfdimtype,&
          spffile(ii),&
          jj)
     numloaded=numloaded+jj
     OFLWR "Loaded ", jj, " spfs from file ",spffile(ii); CFL
  enddo
  jj=0
  ll=numloaded
  do ii=1,ll
     flag=0
     do kk=1,numskiporbs
        if (orbskip(kk).eq.ii) then
           flag=1
           exit
        endif
     enddo
     if (flag.eq.0) then
        jj=jj+1
        if (jj.gt.nspf+numfrozen) then
           OFLWR "Dooga error.",ii,jj,nspf,numfrozen,numskiporbs,numloaded; CFLST
        endif
        inspfs(:,jj)=inspfs0(:,ii)
     else
        numloaded=numloaded-1
     endif
  enddo


end subroutine load_spfs

subroutine load_spfs0(inspfs, inspfdims, innspf, dimtypes, infile, numloaded)
  use fileptrmod
  implicit none
  character :: infile*(*)
  integer :: inspfdims(3), innspf,readdims(3),readnspf, dimtypes(3),numloaded,readcflag,myiostat
  DATATYPE :: inspfs(inspfdims(1),inspfdims(2),inspfdims(3),innspf)

  open(999,file=infile, status="unknown", form="unformatted",iostat=myiostat)

  if (myiostat.ne.0) then
     OFLWR "File ",infile," not found for spf read."; CFLST
  endif
  call spf_header_read(999,readdims,readnspf,readcflag)
  call spf_read0(999,innspf,inspfdims,readnspf,readdims,readcflag,dimtypes,inspfs)
  close(999)
  numloaded=min(readnspf,innspf)

end subroutine load_spfs0


! if parorbsplit=3:
!   outdims is local
!   readdims is global

subroutine spf_read0(iunit,outnspf,outdims,readnspf,bigreaddims,readcflag,dimtypes,outspfs)
!! use fileptrmod  !! need nprocs and parorbsplit now
  use mpimod
  use parameters

  implicit none
  integer :: iunit, dimtypes(3),outdims(3),readnspf,outnspf,readcflag,&
       bigreaddims(3),bigoutdims(3)  !! readdims(3),
  real*8,allocatable :: realspfs(:,:,:,:)
  complex*16,allocatable :: cspfs(:,:,:,:)
  real*8,allocatable :: bigoutrealspfs(:,:,:,:)
  complex*16,allocatable :: bigoutcspfs(:,:,:,:)
  real*8 :: outrealspfs(outdims(1),outdims(2),outdims(3),outnspf)
  complex*16 :: outcspfs(outdims(1),outdims(2),outdims(3),outnspf)

  DATATYPE :: outspfs(outdims(1),outdims(2),outdims(3),outnspf)
  integer :: numloaded,itop(3),ibot(3),otop(3),obot(3),idim
  integer :: qqblocks(nprocs),qqend(nprocs),qqstart(nprocs),ii,ispf,oosize

     numloaded = min(outnspf,readnspf)


!!$  bigreaddims(:)=readdims(:)
  bigoutdims(:)=outdims(:)
  if (parorbsplit.eq.3) then                !! see spf_header_read:
!!$     bigreaddims(3)=readdims(3)*nprocs      !! it is this variable bigreaddims not readdims that will not
     bigoutdims(3)=outdims(3)*nprocs        !!   change as nprocs is varied; likewise bigoutdims/outdims.
  endif 

  if (myrank.ne.1.and.parorbsplit.eq.3) then
     OFLWR "   --> rank 1 will read. "; CFL

     if (readcflag.eq.0) then
        allocate(realspfs(1,1,1,readnspf))
     else
        allocate(cspfs(1,1,1,readnspf))
     endif
  else
     OFLWR "Reading spfs."; CFL

     if (readcflag.eq.0) then
        allocate(realspfs(bigreaddims(1),bigreaddims(2),bigreaddims(3),readnspf))
     else
        allocate(cspfs(bigreaddims(1),bigreaddims(2),bigreaddims(3),readnspf))
     endif

  endif



  do idim=1,3
     select case(dimtypes(idim))
     case (0)
        ibot(idim)=1;        itop(idim)=min(bigreaddims(idim),bigoutdims(idim))
        obot(idim)=ibot(idim);        otop(idim)=itop(idim)
        
     case (1)
        if (mod(bigreaddims(idim)-bigoutdims(idim),2).eq.1) then
           OFLWR "On read must have both even or both odd points for dimension ",idim," they are ", &
                bigreaddims(idim),bigoutdims(idim); CFLST
        endif
        ibot(idim)=max(1,(bigreaddims(idim)-bigoutdims(idim))/2+1)
        obot(idim)=max(1,(bigoutdims(idim)-bigreaddims(idim))/2+1)
        itop(idim)=ibot(idim)+min(bigoutdims(idim),bigreaddims(idim))-1
        otop(idim)=obot(idim)+min(bigoutdims(idim),bigreaddims(idim))-1
     case (2)
        if (bigoutdims(idim).ne.bigreaddims(idim)) then
           OFLWR "For dimension ",idim," spfdimtype is 2 so dimensions must agree ", &
                bigreaddims(idim),bigoutdims(idim); CFLST
        endif
        itop(idim)=bigoutdims(idim);        otop(idim)=bigoutdims(idim)
        ibot(idim)=1;        obot(idim)=1
     case default
        OFLWR "spfdimtype not recognized for dim", idim," it is ", dimtypes(idim); CFLST
     end select
  end do
  

  if (myrank.eq.1.or.parorbsplit.ne.3) then
     if (readcflag.eq.0) then
        OFLWR "Read real spfs"; CFL
        read(iunit) realspfs(:,:,:, 1:numloaded )
     else
        OFLWR "Read complex spfs"; CFL
        read(iunit) cspfs(:,:,:, 1:numloaded)
     endif
     close(iunit)
     OFLWR "   ...okay"; CFL     


!!$      if (reinterp_orbflag.ne.0) then
!!$         if (oldspacing.le.0d0) then
!!$            OFLWR "must input oldspacing>0 in parinp for reinterp_orbflag"; CFLST
!!$         endif
!!$         if (readcflag.eq.0) then
!!$            call reinterpolate_orbs_real(realspfs(:,:,:,1:numloaded),bigreaddims,numloaded,oldspacing)
!!$         else
!!$            call reinterpolate_orbs_complex(cspfs(:,:,:,1:numloaded),bigreaddims,numloaded,oldspacing)
!!$         endif
!!$      endif

  endif

  outspfs(:,:,:,:)=0d0

  if (parorbsplit.ne.3) then
     if (readcflag.eq.0) then
        outspfs(obot(1):otop(1), obot(2):otop(2), obot(3):otop(3), 1:numloaded) =  &
             realspfs(ibot(1):itop(1), ibot(2):itop(2), ibot(3):itop(3), 1:numloaded )
     else
        outspfs(obot(1):otop(1), obot(2):otop(2), obot(3):otop(3), 1:numloaded) =  &
             cspfs(ibot(1):itop(1), ibot(2):itop(2), ibot(3):itop(3),  1:numloaded )
     endif
  else
     if (myrank.eq.1) then
        if (readcflag.eq.0) then
           allocate(bigoutrealspfs(bigoutdims(1),bigoutdims(2),bigoutdims(3),numloaded))
           bigoutrealspfs(:,:,:,:)=0d0
           bigoutrealspfs(obot(1):otop(1), obot(2):otop(2), obot(3):otop(3), 1:numloaded) =  &
                realspfs(ibot(1):itop(1), ibot(2):itop(2), ibot(3):itop(3), 1:numloaded )
        else
           allocate(bigoutcspfs(bigoutdims(1),bigoutdims(2),bigoutdims(3),numloaded))
           bigoutcspfs(:,:,:,:)=0d0
           bigoutcspfs(obot(1):otop(1), obot(2):otop(2), obot(3):otop(3), 1:numloaded) =  &
                cspfs(ibot(1):itop(1), ibot(2):itop(2), ibot(3):itop(3),  1:numloaded )
        endif
     else
        if (readcflag.eq.0) then
           allocate(bigoutrealspfs(1,1,1,numloaded))
        else
           allocate(bigoutcspfs(1,1,1,numloaded))
        endif
     endif
     oosize=outdims(1)*outdims(2)*outdims(3)
     qqblocks(:)=oosize
     do ii=1,nprocs
        qqend(ii)=ii*oosize; qqstart(ii)=(ii-1)*oosize+1
     enddo
     do ispf=1,numloaded
        if (readcflag.eq.0) then
           call myscatterv_real(bigoutrealspfs(:,:,:,ispf),outrealspfs(:,:,:,ispf),&
                oosize*nprocs,qqstart(myrank),qqend(myrank),qqblocks(:),qqstart(:))
           outspfs(:,:,:,ispf)= outrealspfs(:,:,:,ispf)
        else
           call myscatterv_complex(bigoutcspfs(:,:,:,ispf),outcspfs(:,:,:,ispf),&
                oosize*nprocs,qqstart(myrank),qqend(myrank),qqblocks(:),qqstart(:))
           outspfs(:,:,:,ispf)= outcspfs(:,:,:,ispf)
        endif
     enddo
     if (readcflag.eq.0) then
        deallocate(bigoutrealspfs)
     else
        deallocate(bigoutcspfs)
     endif
  endif
  if (readcflag.eq.0) then
     deallocate(realspfs)
  else
     deallocate(cspfs)
  endif
end subroutine spf_read0



!! initial version for parorbsplit=3 with multiple files
!!$subroutine save_vector(psi,afile,sfile)
!!$  use parameters
!!$  use mpimod
!!$  implicit none
!!$  DATATYPE :: psi(psilength) 
!!$  character :: afile*(*), sfile*(*)
!!$  integer :: iprop,getlen
!!$  character (len=3) :: iilab
!!$  character (len=4) :: iilab0
!!$
!!$  if (myrank.ne.1.and.parorbsplit.ne.3) then
!!$     OFLWR "   --> rank 1 will save. "; CFL;     return
!!$  endif
!!$
!!$  if (myrank.eq.1) then
!!$     open(999,file=afile, status="unknown", form="unformatted")
!!$     call avector_header_write(999,mcscfnum)
!!$  endif
!!$
!!$  if (parorbsplit.eq.3) then
!!$     write(iilab0,'(I4)') myrank+1000
!!$     iilab(:)=iilab0(2:4)
!!$     open(998,file=sfile(1:getlen(sfile)-1)//iilab, status="unknown", form="unformatted")
!!$  else
!!$     open(998,file=sfile, status="unknown", form="unformatted")
!!$  endif
!!$
!!$  call spf_head er_w rite(998,nspf+numfrozen)
!!$  call wr ite_spf(998,psi(spfstart));  close(998)
!!$
!!$  if (myrank.eq.1) then
!!$     do iprop=1,mcscfnum
!!$        call write_avector(999,psi(astart(iprop)))
!!$     enddo
!!$     close(999)
!!$  endif
!!$
!!$  OFLWR "SAVED vectors!"; CFL 
!!$end subroutine save_vector

!!$ old version for multiple files parorbsplit=3 (untested)
!!$subroutine load_spfs(inspfs, numloaded)    !! both output
!!$  use parameters
!!$  use mpimod
!!$  implicit none
!!$  integer :: numloaded,ii,jj,kk,ll,flag,getlen
!!$  DATATYPE :: inspfs(spfsize,nspf+numfrozen), inspfs0(spfsize,nspf+numfrozen+numskiporbs)
!!$  character (len=3) :: iilab
!!$  character (len=4) :: iilab0
!!$
!!$  write(iilab0,'(I4)') myrank+1000
!!$  iilab(:)=iilab0(2:4)
!!$
!!$  numloaded=0
!!$  do ii=1,numspffiles
!!$
!!$     if (parorbsplit.eq.3) then
!!$        call load_sp fs0(inspfs0(:,numloaded+1),spfdims,nspf+numfrozen+numskiporbs-numloaded,spfdimtype,&
!!$             spffile(ii)(1:getlen(spffile(ii))-1)//iilab,&
!!$             jj)
!!$     else
!!$        call loa d_spfs0(inspfs0(:,numloaded+1),spfdims,nspf+numfrozen+numskiporbs-numloaded,spfdimtype,&
!!$             spffile(ii),&
!!$             jj)
!!$     endif
!!$
!!$
!!$     numloaded=numloaded+jj
!!$     OFLWR "Loaded ", jj, " spfs from file ",spffile(ii); CFL
!!$  enddo
!!$  jj=0
!!$  ll=numloaded
!!$  do ii=1,ll
!!$     flag=0
!!$     do kk=1,numskiporbs
!!$        if (orbskip(kk).eq.ii) then
!!$           flag=1
!!$           exit
!!$        endif
!!$     enddo
!!$     if (flag.eq.0) then
!!$        jj=jj+1
!!$        if (jj.gt.nspf+numfrozen) then
!!$           OFLWR "Dooga error.",ii,jj,nspf,numfrozen,numskiporbs,numloaded; CFLST
!!$        endif
!!$        inspfs(:,jj)=inspfs0(:,ii)
!!$     else
!!$        numloaded=numloaded-1
!!$     endif
!!$  enddo
!!$
!!$
!!$end subroutine load_spfs


!!$ before new version for parorbsplit=3
!!$subroutine spf_r ead0(iunit,outnspf,outdims,readnspf,readdims,readcflag,dimtypes,outspfs)
!!$  use fileptrmod
!!$  implicit none
!!$  integer :: iunit, dimtypes(3),readdims(3),outdims(3),readnspf,outnspf,readcflag
!!$  real*8 :: realspfs(readdims(1),readdims(2),readdims(3),readnspf)
!!$  complex*16 :: cspfs(readdims(1),readdims(2),readdims(3),readnspf)
!!$  DATATYPE :: outspfs(outdims(1),outdims(2),outdims(3),outnspf)
!!$  integer :: numloaded,itop(3),ibot(3),otop(3),obot(3),idim
!!$
!!$  outspfs(:,:,:,:)=0d0
!!$  do idim=1,3
!!$     select case(dimtypes(idim))
!!$     case (0)
!!$        ibot(idim)=1;        itop(idim)=min(readdims(idim),outdims(idim))
!!$        obot(idim)=ibot(idim);        otop(idim)=itop(idim)
!!$
!!$     case (1)
!!$        if (mod(readdims(idim)-outdims(idim),2).eq.1) then
!!$           OFLWR "On read must have both even or both odd points for dimension ",idim," they are ", &
!!$                readdims(idim),outdims(idim); CFLST
!!$        endif
!!$        ibot(idim)=max(1,(readdims(idim)-outdims(idim))/2+1)
!!$        obot(idim)=max(1,(outdims(idim)-readdims(idim))/2+1)
!!$        itop(idim)=ibot(idim)+min(outdims(idim),readdims(idim))-1
!!$        otop(idim)=obot(idim)+min(outdims(idim),readdims(idim))-1
!!$     case (2)
!!$        if (outdims(idim).ne.readdims(idim)) then
!!$           OFLWR "For dimension ",idim," spfdimtype is 2 so dimensions must agree ", &
!!$                readdims(idim),outdims(idim); CFLST
!!$        endif
!!$        itop(idim)=outdims(idim);        otop(idim)=outdims(idim)
!!$        ibot(idim)=1;        obot(idim)=1
!!$     case default
!!$        OFLWR "spfdimtype not recognized for dim", idim," it is ", dimtypes(idim); CFLST
!!$     end select
!!$  end do
!!$ 
!!$  numloaded = min(outnspf,readnspf)
!!$  if (readcflag.eq.0) then
!!$     OFLWR "Read real spfs"; CFL
!!$     read(iunit) realspfs(:,:,:, 1:numloaded )
!!$     outspfs(obot(1):otop(1), obot(2):otop(2), obot(3):otop(3), 1:numloaded) =  &
!!$          realspfs(ibot(1):itop(1), ibot(2):itop(2), ibot(3):itop(3), 1:numloaded )
!!$  else
!!$     OFLWR "Read complex spfs"; CFL
!!$     read(iunit) cspfs(:,:,:, 1:numloaded)
!!$     outspfs(obot(1):otop(1), obot(2):otop(2), obot(3):otop(3), 1:numloaded) =  &
!!$          cspfs(ibot(1):itop(1), ibot(2):itop(2), ibot(3):itop(3), 1:numloaded )
!!$  endif
!!$  close(iunit)
!!$end subroutine spf_r ead0


