
  
#include "Definitions.INC"
  
  
program mctdhf
  use mpimod
  use xxxmod
  use opmod !! frozenspfs
  use parameters
  implicit none

  integer :: i,spfsloaded,totread,ifile,readnum
  DATATYPE, allocatable ::  tempvals(:)
  DATATYPE, allocatable :: bigavector(:,:), bigspfs(:,:)
  logical :: logcheckpar

  spfsloaded=0
  pi=4.d0*atan(1.d0)

  call system("date")
  call system("mkdir Dat"); call system("mkdir Bin"); call system("mkdir Flux")
  call system("mkdir WALKS");  

  call MPIstart()
  
  call openfile()
  write(mpifileptr, *) " ****************************************************************************"     
  write(mpifileptr, *) " **********************  BEGIN LBNL AMO MCTDHF ******************************"
  write(mpifileptr, *) " ****************************************************************************"     
  write(mpifileptr, *) 
  write(mpifileptr, *) "   Atomic/Diatomic/Cartesian Polyatomic nonadiabatic MCTDHF code"
  write(mpifileptr, *) "     for ultrafast electronic and nuclear dynamics in strong laser fields"
  write(mpifileptr, *) "   AMO Theory Group, Lawrence Berkeley Laboratory, 2015"
  write(mpifileptr, *) "     DJ Haxton, C W McCurdy, T N Rescigno, K Lawler, B Abeln, X Li.."
  write(mpifileptr, *) "   VERSION 1.0 "
  write(mpifileptr, *) 
  write(mpifileptr, *) "   PROCESSOR ",myrank," OF ",nprocs
  write(mpifileptr, *)
#ifdef REALGO
  write(mpifileptr, *) "   This is MCTDH: Real-valued version."
#else
#ifdef CNORMFLAG
  write(mpifileptr, *) "   This is CMCTDH: C-NORMed ECS version."
#else
#ifdef ECSFLAG
  write(mpifileptr, *) "   This is CHMCTDH: HERM-normed ECS version."
#else
  write(mpifileptr, *) "   This is PMCTDH: HERM-normed non-ECS version."
#endif
#endif
#endif
  call closefile()

!! inpfile input; others output

  call getinpfile()
  call getmyparams(mpifileptr,inpfile,spfdims,spfdimtype,reducedpotsize,numr,nucrepulsion,nonuc_checkflag) 
  spfsize=spfdims(1)*spfdims(2)*spfdims(3)

  
  allocate(bondpoints(numr),bondweights(numr),elecweights(spfdims(1),spfdims(2),spfdims(3)),elecradii(spfsize))

  call getparams()
  call system("mkdir "//timingdir)  

  call transferparams(nspf,spfrestrictflag,spfmvals,spfugrestrict,spfugvals,spfsmallsize,logcheckpar);
  if (logcheckpar) then
     if (parorbsplit.ne.3) then
        OFLWR "Setting parorbsplit=3 - you have orbital parallelization in your project";CFL
        parorbsplit=3
     endif
  else
     if (parorbsplit.eq.3) then
        OFLWR "You don't have orbital parallelization in the coordinate system dependent routine."
        WRFL "   ..don't set parorbsplit=3"; CFLST
     endif
  endif

  if (spfsmallsize.gt.spfsize) then
     OFLWR "WTF."; CFLST
  endif
!  print *, "OK CHECKME TRANSFERPARMS (MOVED UP WAS TWOEALLOC)"; stop


  call getclasses()

  call mpiorbsets()

     call fast_newconfiglist(.false.);   


  if (numconfig.eq.0) then
     OFLWR " No configs!! "; CFLST
  endif

  totspfdim=nspf*spfsize; totadim=numconfig*numr;  
  psilength=totadim*mcscfnum+totspfdim

  do i=1,mcscfnum
     astart(i)=1+(i-1)*totadim;         ;     aend(i)=i*totadim;
  enddo
  spfstart=totadim*mcscfnum+1;           spfend=psilength

  if ((sparseconfigflag.eq.0).and.(numconfig.gt.1000).and.(nosparseforce.eq.0)) then
     OFLWR "You should really turn sparseconfigflag on, with ", &
          numconfig*numr, "configurations.";     CFLST
  endif

  call opalloc()

  call fast_newconfiglist(.true.);   

  if ((skipflag.lt.2)) then

  call walkalloc();             call walks()

     if (spinwalkflag==1) then
        call spinwalkalloc(); 
        call spinwalks()
        call spinsets_first()
        call configspin_matel()
        call configspinset_projector()
     endif

     if (walkwriteflag.ne.0) then
        OFLWR "Closing walks.BIN"; close(751); CFL
     endif

     call singlewalkwrite()

     call configalloc()
     if (dfrestrictflag.ne.0) then
        call getdfcon()
     endif

  endif

  call myprojectalloc()      !! Initialize coordinate-dependent arrays.

  call xalloc() !!   INITIALIZE XXX/YYY VECTORS!  

  allocate(bigspfs(spfsize,nspf+numfrozen))

  if (skipflag.eq.0) then
     if (loadspfflag.eq.1) then
        call load_spfs(bigspfs(:,:),spfsloaded)
        if (numfrozen.gt.0) then
           frozenspfs(:,:)=bigspfs(:,1:numfrozen)
           spfsloaded=spfsloaded-numfrozen
           bigspfs(:,1:nspf)=bigspfs(:,numfrozen+1:numfrozen+nspf)
        endif

!! hold it this is done below 041114           yyy%cmfpsivec(spfstart:spfend,0) = RESHAPE(bigspfs(:,:),(/totspfdim/))

     endif  ! loadspfflag
  endif  !! skipflag.eq.0


  OFLWR "Calling init_project",sizeof(pot); CFL
  call init_project(bigspfs,spfsloaded,pot,halfniumpot,rkemod,proderivmod,skipflag,&
       bondpoints,bondweights,elecweights,elecradii,numelec)
  OFLWR "   ....Called init_project."; CFL

!  print *, bigspfs
!  print *, "THAT WAS BIGSPSFS"; stop

!  NOW TRANSFERPARAMS
!  call twoealloc(nspf,spfrestrictflag,spfmvals,spfugrestrict,spfugvals,spfsmallsize);
!  if (spfsmallsize.gt.spfsize) then
!     OFLWR "WTF."; CFLST
!  endif


  if (skipflag.lt.2) then
   
     call spf_orthogit_gs(bigspfs(:,:))
     if (messflag.ne.0) then
        call mess_with_spfs(bigspfs)
        call spf_orthogit_gs(bigspfs(:,:))
     endif

     yyy%cmfpsivec(spfstart:spfend,0) = RESHAPE(bigspfs(:,1:nspf),(/totspfdim/))
     
     call apply_spf_constraints(yyy%cmfpsivec(spfstart,0))


     call spf_orthogit_gs(yyy%cmfpsivec(spfstart,0))
     
     if (numfrozen.ne.0) then
        call noparorbsupport("frozen_matels")
        call frozen_matels()
     endif
!! MAY 2014 now not doing load avector if skipflag=1 (flux)

  endif  !! skipflag.lt.2

  deallocate(bigspfs)

  yyy%cmfpsivec(astart(1):aend(mcscfnum),0) = 0d0

  if (skipflag.eq.0) then

     allocate(bigavector(totadim,mcscfnum));  bigavector(:,:)=0d0

     totread=0

     if (loadavectorflag.eq.1) then

        do ifile=1,numavectorfiles
           if (totread.lt.mcscfnum) then
              OFLWR "Reading avector..." ; CFL
              call load_avectors(avectorfile(ifile),bigavector(:,totread+1),mcscfnum-totread,readnum,avecloadskip(ifile))
           endif
           totread=totread+readnum
        enddo

        OFLWR "Read ",totread," A-vectors"; CFL

     endif

     if (totread.lt.mcscfnum) then
        if (improvedrelaxflag.eq.0.and.loadavectorflag.ne.0) then
           OFLWR "Loading avectors, but not enough on file for prop..."; CFLST
        endif

        OFLWR "Not enough avectors loaded! will diagonalize."; CFL
        do i=totread+1,mcscfnum
           call staticvector(bigavector(:,i),totadim)
           call gramschmidt(totadim,i-1,totadim,bigavector(:,:),bigavector(:,i),.false.)
        enddo
        
        call all_matel()
        allocate(tempvals(mcscfnum))
        call myconfigeig(bigavector,tempvals,mcscfnum,1,min(totread,1),0d0,max(0,improvedrelaxflag-1))
        deallocate(tempvals)

     endif  

     yyy%cmfpsivec(astart(1):aend(mcscfnum),0) = RESHAPE(bigavector(:,:),(/totadim*mcscfnum/))

     deallocate(bigavector)
  endif   !! skipflag lt 2

  if (cdenflag.ne.0) then
     call natprojalloc()
  endif

!! ****  CORE ROUTINE **** !!

  call prop_loop( 0.d0 )

!! *********************** !!

  call opdealloc()
!  call myprojectdealloc()  ;  
  call twoedealloc();     call xdealloc()

  deallocate(myavectorhole,myavectorexcitefrom,myavectorexciteto)

  if (dfrestrictflag.ne.0) then
     call dfcondealloc()
  endif
  
  call natprojdealloc(); 

  call system("date")
  call mpistop()

end program mctdhf



subroutine getclasses()
  use parameters
  implicit none
  integer :: iclass,ispf,jspf,flag

  allocate(orbclass(nspf))

  iclass=0
  do ispf=1,nspf
     flag=0
     do jspf=1,ispf-1
        
        if ( (spfrestrictflag.eq.0.or.spfmvals(jspf).eq.spfmvals(ispf)) .and. &
             (spfugrestrict.eq.0.or.spfugvals(jspf).eq.spfugvals(ispf)) ) then
           orbclass(ispf)=orbclass(jspf)
           flag=1
           exit
        endif
     enddo
     if (flag.eq.0) then
        iclass=iclass+1
        orbclass(ispf)=iclass
     endif
  enddo

!  OFLWR "ORBCLASS", orbclass(:); CFL

  numclasses=iclass
  allocate(classorb(nspf,numclasses),nperclass(numclasses))
  
  do iclass=1,numclasses
     jspf=0
     do ispf=1,nspf
        if (orbclass(ispf).eq.iclass) then
           jspf=jspf+1
           classorb(jspf,iclass)=ispf
        endif
     enddo
     nperclass(iclass)=jspf
  enddo
!  OFLWR "NPERCLASS", nperclass
!  WRFL
!  WRFL "CLASSORB"
!  do iclass=1,numclasses
!     WRFL classorb(1:nperclass(iclass),iclass)
!  enddo

end subroutine getclasses



