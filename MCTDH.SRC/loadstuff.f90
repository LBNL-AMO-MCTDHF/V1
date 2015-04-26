
#include "Definitions.INC"

subroutine save_vector(psi,afile,sfile)
  use parameters
  use mpimod
  implicit none
  DATATYPE :: psi(psilength) 
  character :: afile*(*), sfile*(*)
  integer :: iprop,qqblocks(nprocs),ispf
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

     do ispf=1,nspf
#ifdef REALGO
        call mygatherv_real(psi(spfstart+(ispf-1)*spfsize),parorbitals(:,:,ispf), qqblocks(:),.false.)
#else
        call mygatherv_complex(psi(spfstart+(ispf-1)*spfsize),parorbitals(:,:,ispf), qqblocks(:),.false.)
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
  use mpimod
  implicit none
  character :: infile*(*)
  integer :: inspfdims(3), innspf,readdims(3),readnspf, dimtypes(3),numloaded,readcflag,myiostat
  DATATYPE :: inspfs(inspfdims(1),inspfdims(2),inspfdims(3),innspf)

  if (myrank.eq.1) then
     open(999,file=infile, status="unknown", form="unformatted",iostat=myiostat)
  endif
  call mympiibcastone(myiostat,1)

  if (myiostat.ne.0) then
     OFLWR "File ",infile," not found for spf read."; CFLST
  endif
  if (myrank.eq.1) then
     call spf_header_read(999,readdims,readnspf,readcflag)
     call mympiibcast(readdims,1,3); call mympiibcastone(readnspf,1); call mympiibcastone(readcflag,1)
     call spf_read0(999,innspf,inspfdims,readnspf,readdims,readcflag,dimtypes,inspfs)
     close(999)
  else
     call mympiibcast(readdims,1,3); call mympiibcastone(readnspf,1); call mympiibcastone(readcflag,1)
     call spf_read0(-66666,innspf,inspfdims,readnspf,readdims,readcflag,dimtypes,inspfs)
  endif

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
  integer :: qqblocks(nprocs),ispf

     numloaded = min(outnspf,readnspf)


!!$  bigreaddims(:)=readdims(:)
  bigoutdims(:)=outdims(:)
  if (parorbsplit.eq.3) then                !! see spf_header_read:
!!$     bigreaddims(3)=readdims(3)*nprocs      !! it is this variable bigreaddims not readdims that will not
     bigoutdims(3)=outdims(3)*nprocs        !!   change as nprocs is varied; likewise bigoutdims/outdims.
  endif 

  if (myrank.ne.1.and.parorbsplit.eq.3) then
     OFLWR "   --> rank 1 will read. "; CFL
     allocate(realspfs(1,1,1,readnspf))
     allocate(cspfs(1,1,1,readnspf))
  else
     OFLWR "Reading spfs."; CFL
     if (readcflag.eq.0) then
        allocate(realspfs(bigreaddims(1),bigreaddims(2),bigreaddims(3),readnspf))
        allocate(cspfs(1,1,1,1))
     else
        allocate(cspfs(bigreaddims(1),bigreaddims(2),bigreaddims(3),readnspf))
        allocate(realspfs(1,1,1,1))
     endif
  endif
  realspfs=0; cspfs=0


  do idim=1,3
     select case(dimtypes(idim))
     case (0)
        ibot(idim)=1;        itop(idim)=min(bigreaddims(idim),bigoutdims(idim))
        obot(idim)=ibot(idim);        otop(idim)=itop(idim)
        
     case (1)

!XXARGH modulus function stupidly defined    if (mod(bigreaddims(idim)-bigoutdims(idim),2).eq.1) then

        if (reinterp_orbflag.eq.0) then
           if (mod(bigreaddims(idim)+bigoutdims(idim),2).eq.1) then
              OFLWR "On read must have both even or both odd points for dimension ",idim," they are ", &
                   bigreaddims(idim),bigoutdims(idim); CFLST
              
           endif

           ibot(idim)=max(1,(bigreaddims(idim)-bigoutdims(idim))/2+1)
           obot(idim)=max(1,(bigoutdims(idim)-bigreaddims(idim))/2+1)
           itop(idim)=ibot(idim)+min(bigoutdims(idim),bigreaddims(idim))-1
           otop(idim)=obot(idim)+min(bigoutdims(idim),bigreaddims(idim))-1

        else

           ibot(idim)=1    !! these are not used
           obot(idim)=1
           itop(idim)=bigoutdims(idim)
           otop(idim)=bigoutdims(idim)
        endif

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

!! sinc dvr only, half spacing reinterpolation only
     
  endif

  outspfs(:,:,:,:)=0d0

  if ((parorbsplit.eq.3.and.myrank.eq.1).or.(parorbsplit.ne.3.and.reinterp_orbflag.ne.0)) then
     if (readcflag.eq.0) then
        allocate(bigoutrealspfs(bigoutdims(1),bigoutdims(2),bigoutdims(3),numloaded))
        allocate(bigoutcspfs(1,1,1,numloaded))
     else
        allocate(bigoutcspfs(bigoutdims(1),bigoutdims(2),bigoutdims(3),numloaded))
        allocate(bigoutrealspfs(1,1,1,numloaded))
     endif
  else
     allocate(bigoutrealspfs(1,1,1,numloaded))
     allocate(bigoutcspfs(1,1,1,numloaded))
  endif
  bigoutrealspfs=0; bigoutcspfs=0

  
  if (reinterp_orbflag.ne.0.and.(myrank.eq.1.or.parorbsplit.ne.3)) then
     if (readcflag.eq.0) then
        call reinterpolate_orbs_real(realspfs(:,:,:,1:numloaded),bigreaddims,bigoutrealspfs,bigoutdims,numloaded)
     else
        call reinterpolate_orbs_complex(cspfs(:,:,:,1:numloaded),bigreaddims,bigoutcspfs,bigoutdims,numloaded)
     endif
  endif

  if (parorbsplit.ne.3) then
     if (reinterp_orbflag.ne.0) then
        if (readcflag.eq.0) then
           outspfs(:,:,:,1:numloaded)=bigoutrealspfs(:,:,:,:)
        else
           outspfs(:,:,:,1:numloaded)=bigoutcspfs(:,:,:,:)
        endif
     else
        if (readcflag.eq.0) then
           outspfs(obot(1):otop(1), obot(2):otop(2), obot(3):otop(3), 1:numloaded) =  &
                realspfs(ibot(1):itop(1), ibot(2):itop(2), ibot(3):itop(3), 1:numloaded )
        else
           outspfs(obot(1):otop(1), obot(2):otop(2), obot(3):otop(3), 1:numloaded) =  &
                cspfs(ibot(1):itop(1), ibot(2):itop(2), ibot(3):itop(3),  1:numloaded )
        endif
     endif
  else
     if (myrank.eq.1.and.reinterp_orbflag.eq.0) then
        if (readcflag.eq.0) then
           bigoutrealspfs(:,:,:,:)=0d0
           bigoutrealspfs(obot(1):otop(1), obot(2):otop(2), obot(3):otop(3), 1:numloaded) =  &
                realspfs(ibot(1):itop(1), ibot(2):itop(2), ibot(3):itop(3), 1:numloaded )
        else
           bigoutcspfs(:,:,:,:)=0d0
           bigoutcspfs(obot(1):otop(1), obot(2):otop(2), obot(3):otop(3), 1:numloaded) =  &
                cspfs(ibot(1):itop(1), ibot(2):itop(2), ibot(3):itop(3),  1:numloaded )
        endif
     endif

     qqblocks(:)=outdims(1)*outdims(2)*outdims(3)

     do ispf=1,numloaded
        if (readcflag.eq.0) then
           call myscatterv_real(bigoutrealspfs(:,:,:,ispf),outrealspfs(:,:,:,ispf),qqblocks(:))
           outspfs(:,:,:,ispf)= outrealspfs(:,:,:,ispf)
        else
           call myscatterv_complex(bigoutcspfs(:,:,:,ispf),outcspfs(:,:,:,ispf),qqblocks(:))
           outspfs(:,:,:,ispf)= outcspfs(:,:,:,ispf)
        endif
     enddo
  endif

  deallocate(realspfs,bigoutrealspfs,cspfs,bigoutcspfs)

end subroutine spf_read0



