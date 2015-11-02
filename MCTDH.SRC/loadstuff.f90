
#include "Definitions.INC"

subroutine save_vector(psi,afile,sfile)
  use parameters
  use opmod !! frozenspfs
  use mpimod
  implicit none
  DATATYPE :: psi(psilength) 
  character :: afile*(*), sfile*(*)
  integer :: iprop,ispf
  DATATYPE, allocatable :: parorbitals(:,:), parfrozen(:,:), paravec(:,:,:)

!! always allocate avoid warn bounds

  if (parorbsplit.eq.3.and.myrank.eq.1) then
     allocate(parorbitals(spfsize*nprocs,nspf)); parorbitals(:,:)=0d0
     allocate(parfrozen(spfsize*nprocs,max(numfrozen,1))); parfrozen(:,:)=0d0
  else
     allocate(parorbitals(1,nspf), parfrozen(1,max(numfrozen,1)))
  endif

  call mpibarrier()

  if (parorbsplit.eq.3) then     
     do ispf=1,nspf
        call splitgatherv(psi(spfstart+(ispf-1)*spfsize),parorbitals(:,ispf), .false.)
     enddo
     
     if (numfrozen.gt.0) then
        do ispf=1,numfrozen
           call splitgatherv(frozenspfs(:,ispf),parfrozen(:,ispf), .false.)
        enddo
     endif
  endif


  if (par_consplit.ne.0.and.myrank.eq.1) then
     allocate(paravec(numr,num_config,mcscfnum)); paravec(:,:,:)=0d0
  else
     allocate(paravec(1,1,mcscfnum))
  endif

  call mpibarrier()

  if (par_consplit.ne.0) then
     do iprop=1,mcscfnum
        call mygatherv(psi(astart(iprop)),paravec(:,:,iprop),configs_perproc(:)*numr,.false.)
     enddo
  endif

  if (myrank.eq.1) then
     open(998,file=sfile, status="unknown", form="unformatted")
     call spf_header_write(998,nspf+numfrozen)
     
     if (parorbsplit.eq.3) then
        call write_spf(998,parorbitals(:,:),parfrozen(:,:),spfsize*nprocs)
     else
        call write_spf(998,psi(spfstart),frozenspfs(:,:),spfsize)
     endif
     close(998)
     
     open(999,file=afile, status="unknown", form="unformatted")
     call avector_header_write(999,mcscfnum)
     if (par_consplit.ne.0) then
        do iprop=1,mcscfnum
           call write_avector(999,paravec(:,:,iprop))
        enddo
     else
        do iprop=1,mcscfnum
           call write_avector(999,psi(astart(iprop)))
        enddo
     endif
     close(999)
  endif
  
  call mpibarrier()

  deallocate(parorbitals,parfrozen)
  deallocate(paravec)

  if (myrank.eq.1) then
     OFLWR "   ...saved vectors!"; CFL
  else
     OFLWR "   ...rank 1 SAVED vectors!"; CFL 
  endif

end subroutine save_vector



subroutine write_spf(unit,spfin,frozenin,inspfsize)
  use parameters
  implicit none
  integer :: unit,inspfsize
  DATATYPE :: allspf(inspfsize,nspf+numfrozen), spfin(inspfsize, nspf   ),&
       frozenin(inspfsize, numfrozen)
  if (numfrozen.gt.0) then
     allspf(:,1:numfrozen)=frozenin(:,:)
  endif
  allspf(:,numfrozen+1:nspf+numfrozen) = spfin(:,:)
  write (unit) allspf
end subroutine


subroutine spf_header_write(iunit,num)
  use parameters
  implicit none
  integer ::  cflag,iunit,num,qqq(3),ppp(3)
#ifdef REALGO
  cflag=0
#else
  cflag=1
#endif
  qqq(:)=1
  if (parorbsplit.eq.3) then
     ppp(:)=1; call bigdimsub(ppp,qqq)
  endif
  write(iunit) spfdims(:)*qqq(:),num,cflag
end subroutine


subroutine spf_header_read(iunit,outdims,outnspf,cflag)
  implicit none
  integer ::  cflag,outdims(3),outnspf,iunit
#ifdef REALGO
  cflag=0
#else
  cflag=1
#endif
  read(iunit) outdims(:),outnspf,cflag
end subroutine



subroutine load_spfs(inspfs, numloaded)    !! both output
  use parameters
  implicit none
  integer :: numloaded,ii,jj,kk,ll,flag
  DATATYPE :: inspfs(spfsize,nspf+numfrozen), inspfs0(spfsize,nspf+numfrozen+numskiporbs)

  numloaded=0

  if (numspffiles.gt.100) then
     OFLWR "REDIM SPF_GRIDSHIFT?"; CFLST
  endif
  do ii=1,numspffiles
     call load_spfs0(inspfs0(:,numloaded+1),spfdims,nspf+numfrozen+numskiporbs-numloaded,spfdimtype,&
          spffile(ii),&
          jj,spf_gridshift(:,ii))
     eachloaded(ii)=jj
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


subroutine load_spfs0(inspfs, inspfdims, innspf, dimtypes, infile, numloaded,in_gridshift)
  use fileptrmod
  use mpimod
  implicit none
  character,intent(in) :: infile*(*)
  integer, intent(in) :: inspfdims(3), innspf, dimtypes(3),in_gridshift(3)
  DATATYPE,intent(out) :: inspfs(inspfdims(1),inspfdims(2),inspfdims(3),innspf)
  integer, intent(out) :: numloaded
  integer :: readcflag,myiostat, readdims(3),readnspf

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
     call spf_read0(999,innspf,inspfdims,readnspf,readdims,readcflag,dimtypes,inspfs,in_gridshift)
     close(999)
  else
     call mympiibcast(readdims,1,3); call mympiibcastone(readnspf,1); call mympiibcastone(readcflag,1)
     call spf_read0(-798,innspf,inspfdims,readnspf,readdims,readcflag,dimtypes,inspfs,in_gridshift)
  endif

  numloaded=min(readnspf,innspf)

end subroutine load_spfs0


! if parorbsplit=3:
!   outdims is local
!   readdims is global

subroutine spf_read0(iunit,outnspf,outdims,readnspf,bigreaddims,readcflag,dimtypes,outspfs,in_gridshift)
  use mpimod
  use parameters
  implicit none
  integer, intent(in) :: iunit, dimtypes(3),outdims(3),readnspf,outnspf,readcflag,&
       bigreaddims(3),in_gridshift(3)
  DATATYPE,intent(out) :: outspfs(outdims(1),outdims(2),outdims(3),outnspf)
  real*8,allocatable :: realspfs(:,:,:,:), bigoutrealspfs(:,:,:,:)
  complex*16,allocatable :: cspfs(:,:,:,:), bigoutcspfs(:,:,:,:)
  real*8 :: outrealspfs(outdims(1),outdims(2),outdims(3),outnspf)  !!AUTOMATIC
  complex*16 :: outcspfs(outdims(1),outdims(2),outdims(3),outnspf) !!AUTOMATIC
  integer :: numloaded,itop(3),ibot(3),otop(3),obot(3),idim,flag, amin(3),amax(3),bmin(3),bmax(3),&
       ispf, bigoutdims(3),ooshift,rrshift

  numloaded = min(outnspf,readnspf)

!!$  bigreaddims(:)=readdims(:)
  bigoutdims(:)=outdims(:)

  if (parorbsplit.eq.3) then                !! see spf_header_read:

!!!!$     bigreaddims(3)=readdims(3)*nprocs      !! it is this variable bigreaddims not readdims that will not
!!     bigoutdims(3)=outdims(3)*nprocs        !!   change as nprocs is varied; likewise bigoutdims/outdims.

     call bigdimsub(outdims,bigoutdims)

  endif 

  if (myrank.ne.1) then
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
           ooshift=0; rrshift=0
           if (mod(bigreaddims(idim)+bigoutdims(idim),2).eq.1) then

!!$              OFLWR "On read must have both even or both odd points for dimension ",idim," they are ", &
!!$                   bigreaddims(idim),bigoutdims(idim); CFLST

!!$ chop odd-valued grid at the positive end of xyz 
!!$ that means an even numbered grid with a nucleus at centershift=-1,-1,-1 (location 1/2,1/2,1/2 times spacing)
!!$ and an odd numbered grid with a nucleus at 0,0,0 
!!$ can be used back and forth, without invoking spf_gridshift

              rrshift=(-1)*mod(bigreaddims(idim),2)
              ooshift=(-1)*mod(bigoutdims(idim),2)
           endif

           ibot(idim)=max(1,((bigreaddims(idim)+rrshift)-(bigoutdims(idim)+ooshift))/2+1)
           obot(idim)=max(1,((bigoutdims(idim)+ooshift)-(bigreaddims(idim)+rrshift))/2+1)
           itop(idim)=ibot(idim)+min(bigoutdims(idim)+ooshift,bigreaddims(idim)+rrshift)-1
           otop(idim)=obot(idim)+min(bigoutdims(idim)+ooshift,bigreaddims(idim)+rrshift)-1

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

  if (myrank.eq.1) then
     if (readcflag.eq.0) then
        OFLWR "Read real spfs"; CFL
        read(iunit) realspfs(:,:,:, 1:numloaded )
     else
        OFLWR "Read complex spfs"; CFL
        read(iunit) cspfs(:,:,:, 1:numloaded)
     endif
     close(iunit)
     OFLWR "   ...rank 1 read from file"; CFL     

     do idim=1,3
        if (spfdimtype(idim).ne.1.and.in_gridshift(idim).ne.0) then
           OFLWR "spf_gridshift=",in_gridshift(idim)," not supported for dim",idim," because dimtype is ",spfdimtype(idim); CFLST
        endif
     enddo
     flag=0
     do idim=1,3
        if (in_gridshift(idim).ne.0) then
           flag=1
           exit
        endif
     enddo
     if (flag.ne.0) then
        OFLWR "SHIFTING ORBITALS, spf_gridshift=",in_gridshift(:); CFL

        amin(:)=max(1,1-in_gridshift(:)); amax(:)=min(bigreaddims(:),bigreaddims(:)-in_gridshift(:))
        bmin(:)=max(1,1+in_gridshift(:)); bmax(:)=min(bigreaddims(:),bigreaddims(:)+in_gridshift(:))
        
        if (readcflag.eq.0) then
           realspfs(amin(1):amax(1),amin(2):amax(2),amin(3):amax(3),1:numloaded)=&
                realspfs(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),1:numloaded)
           realspfs(1:amin(1)-1,:,:,1:numloaded)=0; realspfs(amax(1)+1:bigreaddims(1),:,:,1:numloaded)=0
           realspfs(:,1:amin(2)-1,:,1:numloaded)=0; realspfs(:,amax(2)+1:bigreaddims(2),:,1:numloaded)=0
           realspfs(:,:,1:amin(3)-1,1:numloaded)=0; realspfs(:,:,amax(3)+1:bigreaddims(3),1:numloaded)=0
        else
           cspfs(amin(1):amax(1),amin(2):amax(2),amin(3):amax(3),1:numloaded)=&
                cspfs(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),1:numloaded)
           cspfs(1:amin(1)-1,:,:,1:numloaded)=0; cspfs(amax(1)+1:bigreaddims(1),:,:,1:numloaded)=0
           cspfs(:,1:amin(2)-1,:,1:numloaded)=0; cspfs(:,amax(2)+1:bigreaddims(2),:,1:numloaded)=0
           cspfs(:,:,1:amin(3)-1,1:numloaded)=0; cspfs(:,:,amax(3)+1:bigreaddims(3),1:numloaded)=0
        endif
        OFLWR "    ..shifted orbitals."; CFL
     endif
  endif

  if (myrank.eq.1) then
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

!!$ (code above)
!!$     numloaded = min(outnspf,readnspf)
!!$

  if (myrank.eq.1) then
     if (reinterp_orbflag.ne.0) then

!! sinc dvr only, half spacing reinterpolation only

        if (readcflag.eq.0) then
           call reinterpolate_orbs_real(realspfs(:,:,:,1:numloaded),bigreaddims,bigoutrealspfs,bigoutdims,numloaded)
        else
           call reinterpolate_orbs_complex(cspfs(:,:,:,1:numloaded),bigreaddims,bigoutcspfs,bigoutdims,numloaded)
        endif
     else
        if (readcflag.eq.0) then
           bigoutrealspfs(obot(1):otop(1), obot(2):otop(2), obot(3):otop(3), 1:numloaded) =  &
                realspfs(ibot(1):itop(1), ibot(2):itop(2), ibot(3):itop(3), 1:numloaded )
        else
           bigoutcspfs(obot(1):otop(1), obot(2):otop(2), obot(3):otop(3), 1:numloaded) =  &
                cspfs(ibot(1):itop(1), ibot(2):itop(2), ibot(3):itop(3),  1:numloaded )
        endif
     endif
  endif


  outspfs(:,:,:,:)=0d0

  if (parorbsplit.ne.3) then
     do ispf=1,numloaded
        if (myrank.eq.1) then
           if (readcflag.eq.0) then
              outspfs(:,:,:,ispf)= bigoutrealspfs(:,:,:,ispf)
           else
              outspfs(:,:,:,ispf)= bigoutcspfs(:,:,:,ispf)
           endif
        endif
        call mympibcast(outspfs(:,:,:,ispf),1,outdims(1)*outdims(2)*outdims(3))
     enddo
  else
     do ispf=1,numloaded
        if (debugflag.ne.0) then
           OFLWR "call scatter orbital ", ispf; CFL
        endif
        if (readcflag.eq.0) then
           call splitscatterv_real(bigoutrealspfs(:,:,:,ispf),outrealspfs(:,:,:,ispf))
           outspfs(:,:,:,ispf)= outrealspfs(:,:,:,ispf)
        else
           call splitscatterv_complex(bigoutcspfs(:,:,:,ispf),outcspfs(:,:,:,ispf))
           outspfs(:,:,:,ispf)= outcspfs(:,:,:,ispf)
        endif

     enddo
  endif

  deallocate(realspfs,bigoutrealspfs,cspfs,bigoutcspfs)

end subroutine spf_read0



