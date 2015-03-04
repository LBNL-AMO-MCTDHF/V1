
!! SEPARATE PROGRAM!  EXPERIMENTAL.

  
#include "Definitions.INC"
  
  
module aparams
  implicit none
  
  integer :: maxnumatomconfig, maxnumatomstates, maxnumatomspf, numatompoints,  numatomfiles=1, numalone, maxtransnumatomconfig, numpairstates

  real*8 :: atomstart, atomfinish  !! first and last radii on file
  character (len=200) :: atomfilenames(100)

  integer :: transnumatomconfig(2)

end module

  

module assemblemod
  implicit none

  DATATYPE, allocatable :: pairham(:,:), pairovl(:,:), pairvects(:,:), pairtemp(:,:)
  real*8, allocatable :: pairvals(:)

  integer, allocatable :: state_pairs(:,:)  !! 2, numpairstates

  integer, allocatable :: numatomstates(:), atom_spfmvals(:,:)

  integer, allocatable :: tempconfig(:), atomside(:)

  logical, allocatable :: isalone(:)
  integer, allocatable :: numatomconfig(:), numatomspf(:)
  !! remember ndof=2*numelec
  integer, allocatable :: atom_mval(:,:), atom_ms(:), atom_spin(:), atom_ndof(:)   

  integer, allocatable :: atomconfiglists(:,:,:)  !! ndof, maxnumatomconfig, numatomfiles

  integer, allocatable :: transatomconfiglists(:,:,:)  !! ndof, maxtransnumatomconfig, 2

  DATATYPE, allocatable :: denmat(:,:), denvects(:,:)
  CNORMTYPE, allocatable :: denvals(:)
  DATATYPE, allocatable :: atomvects(:,:,:) !! maxnumatomconfig, maxnumatomstates, numatomfiles
  DATATYPE, allocatable :: atomspfs(:,:,:) !! spfsize, maxnumatomspf, numatomfiles

  DATATYPE, allocatable :: finalspfs(:,:,:)  !! spfsize, nspf, numatompoints

  DATATYPE, allocatable :: finalvects(:,:,:)  !! numconfig, numpairstates, numatompoints

  DATATYPE, allocatable :: trans_atomvects(:,:,:)   !! maxtransnumatomconfig, maxnumatomstates, 2


end module

subroutine assemblealloc()
  use parameters
  use aparams
  use assemblemod
  implicit none

  allocate(numatomstates(numatomfiles), atomside(numatomfiles),  numatomconfig(numatomfiles), numatomspf(numatomfiles),  atom_ms(numatomfiles), atom_spin(numatomfiles), atom_ndof(numatomfiles))

  allocate(denmat(nspf,nspf), denvals(nspf), denvects(nspf,nspf))

  allocate(tempconfig(ndof), atom_spfmvals(nspf,numatomfiles))

end subroutine assemblealloc

subroutine assemblealloc2()
  use parameters
  use aparams
  use assemblemod
  implicit none

  allocate(state_pairs(2,maxnumatomstates**2))

  allocate(pairovl(maxnumatomstates**2, maxnumatomstates**2), pairham(maxnumatomstates**2, maxnumatomstates**2), pairvects(maxnumatomstates**2, maxnumatomstates**2), pairtemp(maxnumatomstates**2, maxnumatomstates**2), pairvals(maxnumatomstates**2))

  allocate( atomconfiglists(ndof,maxnumatomconfig,numatomfiles),  &
       atomvects(maxnumatomconfig,maxnumatomstates,numatomfiles), &
       atomspfs(spfsize,maxnumatomspf,numatomfiles) , &
       atom_mval(maxnumatomstates,numatomfiles))

  allocate(finalspfs(spfsize,nspf, numatompoints))

  allocate(transatomconfiglists(ndof,maxtransnumatomconfig,2))

  allocate(trans_atomvects(maxtransnumatomconfig, maxnumatomstates, 2))

  allocate(finalvects(numconfig, maxnumatomstates**2, numatompoints))

end subroutine assemblealloc2

subroutine atomread
  use assemblemod
  use aparams
  use mpimod
  use parameters
  implicit none

  integer :: myiostat, i
  character :: ext(9) = (/ "1", "2", "3", "4", "5", "6", "7", "8", "9" /) 

  namelist/assemble/ numatomfiles,atomfilenames

  do i=1,9
     atomfilenames="mcscf.bin."//ext(i)
  enddo

  open(971,file=inpfile, status="old", iostat=myiostat)
  if (myiostat/=0) then
     call openfile()
     write(mpifileptr,*) "No Input.Inp found for parinp, iostat=",myiostat
     call closefile()
     call mpistop()
  else
     read(971,nml=assemble)
  endif

end subroutine atomread


program mcscfassemble
  use mpimod
  use yyymod
  use parameters
  use iglobalmod
  use h2_params
  use opmod
  use h2projectmod
  use assemblemod
  use aparams
  use mcscfmod
  use configmod
  use mcloopmod
  use commonmod

  implicit none

  DATATYPE :: dot, a1, a2, nullspfs
  DATATYPE, allocatable :: ttvect(:)
  real*8 :: norm

  integer :: atomsum(-5:5,2), saveatomsum(-5:5,2)

  integer :: innumatompoints, innumerad, inlbig, inmbig, ifile, jfile, fileptr, &
       ipoint, istate, jstate, mymrestrictflag, myrestrictflag, myallspinproject, &
       myiostat, i, config1, ipair, nullint, info, lwork, &
       config2, config3, phase, reorder, getconfiguration, myspfrestrictflag, ii

  logical :: allowedconfig
  real*8 :: inatomstart, inatomfinish
  real*8 :: mynuccharge1, mynuccharge2, nullreal, nullperm, rsum, nulldouble
  real*8, allocatable :: rwork(:)
  DATATYPE, allocatable :: work(:)

  pi=4.d0*atan(1.d0)

  call MPIstart()

  call getparams()

  totspfdim=nspf*spfsize 

  call opalloc()


  if (numr.ne.1) then
     print *, "Error, set numr=1 (rgridpoints=3)"
     stop
  endif

  call atomread()

  if (numatomfiles.ne.2) then
     print *, "atomfiles ne 2 Not supported"
     stop
  endif

  call assemblealloc()
     

!!  allocate(psi(psilength))

     !! Initialize coordinate-dependent arrays.

  call myprojectalloc()

  if (atomflag==1) then
     print *, "Not supported"
     stop
  endif


  call init_h2(nullspfs,nspf) 
  
  call commonalloc()
  call commonset_h2()
  call twoealloc()
  call get_twoe_new()
  call myprojectdealloc1()


!  iglobalprop=1
!  call xxxptr_init()


!!!!!!!!!!!!!!!!!
!! GO ASSEMBLE !!
!!!!!!!!!!!!!!!!!


  maxnumatomconfig=0
  maxnumatomstates=0 
  maxnumatomspf=0 
  numatompoints=-1


  print *, "Reading header info from files."

  do ifile=1,numatomfiles

     print *, "File ", atomfilenames(ifile)

     fileptr=9954+ifile

     open(fileptr,file=atomfilenames(ifile), form="unformatted", status="old", iostat=myiostat)
     call checkiostat(myiostat)

!! dummy reads for now.
     read(fileptr) mynuccharge1, mynuccharge2, atom_ndof(ifile)

     if ( (mynuccharge1.eq.nuccharge1).and.(mynuccharge2.eq.0) ) then
        atomside(ifile)=1
     else if ( (mynuccharge2.eq.nuccharge2).and.(mynuccharge1.eq.0) ) then
        atomside(ifile)=2
     else
        call openfile()
        write(mpifileptr,*) "Error atom charges: file ", mynuccharge1, mynuccharge2
        write(mpifileptr,*) "                 Current ", nuccharge1, nuccharge2
        call closefile()
        call mpistop()
     endif

     print *, "Is on side ", atomside(ifile)

!! 052412: myallspinproject -> numr   (numr needed; allspinproject not relevant; different from mcnumr)
     read(fileptr) mymrestrictflag, myspfrestrictflag, myrestrictflag, myallspinproject

!!     if ((mymrestrictflag.ne.1).or.(myrestrictflag.ne.1).or.(myallspinproject.ne.1)) then

     if (myspfrestrictflag.ne.1) then
        
        call openfile()

!!        write(mpifileptr,*) "Error: all restrictflags must be on for atom input.  M, ms, S on file: ", mymrestrictflag, myrestrictflag, myallspinproject

        write(mpifileptr,*) "Error: spfrestrictflag must be on for atom input: ", myspfrestrictflag

        call closefile()
        call mpistop()
     endif

     read(fileptr) nullint, atom_ms(ifile), atom_spin(ifile)
!!     print *, "   Has mvalue ", nullint
     print *, "   Has ms     ", atom_ms(ifile)/2.d0
     print *, "     spin     ", atom_spin(ifile)/2.d0


     read(fileptr) numatomconfig(ifile), numatomspf(ifile)
     read(fileptr) atom_spfmvals(1:numatomspf(ifile),ifile)
     read(fileptr) numatomstates(ifile), innumatompoints, inatomstart, inatomfinish, innumerad, inlbig, inmbig

     print *, "    and has ", numatomstates(ifile), " states"
     print *, "    and     ", numatomspf(ifile), " spfs"
     if (ifile==1) then
        numatompoints=innumatompoints
        atomstart=inatomstart
        atomfinish=inatomfinish
     else
        if ((innumatompoints.ne.numatompoints).or.(abs(atomstart-inatomstart).gt.1.d-7).or.(abs(inatomfinish-atomfinish).gt.1.d-7))  then
           print *, "Error, files must be consistent "
           print *, "Atompoints", innumatompoints, numatompoints
           print *, "Start ", inatomstart, atomstart
           print *, "Finish ", inatomfinish, atomfinish
           stop
        endif
     endif

     if ((inlbig.ne.lbig) .or. (inmbig .ne. mbig) .or. (innumerad.ne.numerad) ) then
        call openfile()
        write(mpifileptr,*) "MCSCF vectors on file do not agree with current calculation.  Data on file, current:"
        write(mpifileptr,*) "numerad:         ", innumerad, numerad
        write(mpifileptr,*) "lbig:            ", inlbig, lbig
        write(mpifileptr,*) "mbig:            ", inmbig, mbig
        call closefile()
        call mpistop()
     endif

     if (maxnumatomspf.lt.numatomspf(ifile)) then
        maxnumatomspf=numatomspf(ifile)
     endif
     if (maxnumatomconfig.lt.numatomconfig(ifile)) then
        maxnumatomconfig=numatomconfig(ifile)
     endif
     if (maxnumatomstates.lt.numatomstates(ifile)) then
        maxnumatomstates=numatomstates(ifile)
     endif

  end do

  print *, "maxnumatomspf is ", maxnumatomspf
  print *, "maxnumatomconfig is ", maxnumatomconfig
  print *, "maxnumatomstates is ", maxnumatomstates


  ifile=1
  jfile=2

  spfmvals(1:numatomspf(ifile)) = atom_spfmvals(1:numatomspf(ifile),ifile)
  spfmvals(numatomspf(ifile)+1:numatomspf(ifile)+numatomspf(jfile)) = atom_spfmvals(1:numatomspf(jfile),jfile)




  call openfile()
  write(mpifileptr, *) "Getting configurations."
  call closefile()

  call get_numconfig()

  totadim=numconfig*numr
        
  psilength=totadim+totspfdim
        
  astart=1
  aend=totadim
  spfstart=totadim+1
  spfend=psilength
        
  call openfile()
  write(mpifileptr, *) "Num config: ", numconfig
        


  call get_configlist()


  call walkalloc()
  call walks()
  
  call spinwalkalloc()
  call spinwalks()
  
  call configalloc()

!!   INITIALIZE XXX/YYY VECTORS!  KEEP THIS AS LAST ALLOC STATEMENT.
  
  call xalloc()






  call get_numconfig0(transnumatomconfig(1), atom_ndof(ifile)/2, nspf)
  call get_numconfig0(transnumatomconfig(2), atom_ndof(jfile)/2, nspf)

  maxtransnumatomconfig=max(transnumatomconfig(1),transnumatomconfig(2))



  call assemblealloc2()





  print *
  print *, "Getting trans configlist"
  print *

  ifile=1
  jfile=2
  
  call get_configlist0(atom_ndof(ifile), ndof, nspf, transnumatomconfig(1), transatomconfiglists(:,:, 1))
  
  call get_configlist0(atom_ndof(jfile), ndof, nspf, transnumatomconfig(2), transatomconfiglists(:,:, 2))
  

  print *, "Done"







  call openfile()
  write(mpifileptr,*) "number of points on file:    ", numatompoints
  write(mpifileptr,*) "first point: ", atomstart
  write(mpifileptr,*) "last point:  ", atomfinish
  call closefile()


  print *, "Getting configs and spfs... "

  do ifile=1,numatomfiles

     print *, "File ", atomfilenames(ifile)

     fileptr=9954+ifile

     read(fileptr)  atomconfiglists(1:atom_ndof(ifile),1:numatomconfig(ifile),ifile)
  enddo


  print *, "mc alloc."

  mcnumr=numatompoints
  

  allocate(mcpoints(mcnumr))

  print *, "point loop"

!! NOW LOOP OVER POINTS.

  do ipoint=1,numatompoints

 

     do ifile=1,numatomfiles
        
!!        print *, "File ", atomfilenames(ifile)
        
        fileptr=9954+ifile
        
        if (ifile.eq.1) then
           read(fileptr)  mcpoints(ipoint)
        else
           read(fileptr) rsum
           if (abs(rsum-mcpoints(ipoint)).gt.1.d-7) then
              print *, "point error: point, files: ", ipoint, ifile, rsum,mcpoints(ipoint)
              stop
           endif
        endif

        if (ifile.eq.1) then
           print *, "Doing point ", ipoint, " of ", numatompoints," : ", mcpoints(ipoint)
        endif

        read(fileptr) atom_mval(1:numatomstates(ifile), ifile)

        atomsum=0
        do ii=1,numatomstates(ifile)
           atomsum(atom_mval(ii,ifile),ifile) = atomsum(atom_mval(ii,ifile),ifile) + 1
        enddo

        if (ipoint.eq.1) then
           saveatomsum(:,ifile)=atomsum(:,ifile)
        else
           do ii=-1,1  !!TEMP
              if (atomsum(ii,ifile).ne.saveatomsum(ii,ifile)) then
                 print *, "atomsum err ", ii, atomsum(ii,ifile), saveatomsum(ii,ifile)
                 stop
              endif
           enddo
        endif


!! THIS ONE - DONE!!xx
#ifdef REALGO
        call one_spf_real_read(fileptr,numerad,lbig+1,mbig,numatomspf(ifile),atomspfs(:,:,ifile),maxnumatomspf,myiostat)
#else
        call one_spf_complex_read(fileptr,numerad,lbig+1,mbig,numatomspf(ifile),atomspfs(:,:,ifile),maxnumatomspf,myiostat)
#endif
        call checkiostat(myiostat)

!!        read(fileptr)  atomspfs(:, 1:numatomspf(ifile), ifile) 


        do istate=1,numatomstates(ifile)

           read(fileptr) nullreal  !! eigenval
           read(fileptr) atomvects(1:numatomconfig(ifile) , istate, ifile)

        enddo

     enddo

     ii=0
     do istate=1,numatomstates(1)
        do jstate=1,numatomstates(2)
           if (atom_mval(istate,1)+atom_mval(jstate,2).eq.mrestrictval) then
              ii=ii+1
              state_pairs(1,ii) = istate
              state_pairs(2,ii) = jstate
           endif
        enddo
     enddo
     if (ipoint.eq.1) then
        numpairstates=ii
        print *, "Number of pair states is ", numpairstates, " out of possible ", numatomstates(1)*numatomstates(2)
     else
        if (ii.ne.numpairstates) then
           print *, "numpairstates error ", ii, numpairstates
           stop
        endif
     endif

!     do ii=1,numpairstates
!        print *, "   ", state_pairs(:,ii)
!     enddo
!     print *
     

     ifile=1
     jfile=2


     if (numatomspf(ifile)+numatomspf(jfile).ne.nspf) then
        print *, "Error, for mcscf vects must fill nspf slots."
        print *, "   atom 1, atom 2, nspf: ", numatomspf(ifile),numatomspf(jfile), nspf
        stop
     endif

     finalspfs(:,:,ipoint)=0
     finalspfs(:,1:numatomspf(ifile),ipoint)=atomspfs(:,1:numatomspf(ifile),ifile)
     finalspfs(:,numatomspf(ifile)+1:numatomspf(ifile)+numatomspf(jfile),ipoint)=atomspfs(:,1:numatomspf(jfile),jfile)

     finalvects(:,:,ipoint)=0.d0

!! lowdin orthogonalize joined orbital set
     call spf_orthogit(finalspfs(:,:,ipoint),nspf, nulldouble)
     
!! get new atomic wavefunctions 



     do istate=1,numatomstates(ifile)
        call permoverlaps(1,atom_ndof(ifile)/2, spfsize, finalspfs(:,:,ipoint),atomspfs(:,:,ifile), &
             trans_atomvects(:,istate,1), atomvects(:,istate,ifile), nullperm, 0, 0.0d0, 0.0d0, nspf,&
             numatomspf(ifile), transnumatomconfig(1), numatomconfig(ifile), &
             transatomconfiglists(:,:,1), ndof, atomconfiglists(:,:,ifile),ndof, 1)
        
!           print *, "Dot check: ", dot(atomvects(:,istate,ifile), atomvects(:,istate,ifile), numatomconfig(ifile)), dot(trans_atomvects(:,istate,1), trans_atomvects(:,istate,1), transnumatomconfig(1))

     enddo
     
     
     do istate=1,numatomstates(jfile)
        
        call permoverlaps(1,atom_ndof(jfile)/2, spfsize, finalspfs(:,:,ipoint),atomspfs(:,:,jfile), &
             trans_atomvects(:,istate,2), atomvects(:,istate,jfile), nullperm, 0, 0.0d0, 0.0d0, nspf,&
             numatomspf(jfile), transnumatomconfig(2), numatomconfig(jfile), &
             transatomconfiglists(:,:,2), ndof, atomconfiglists(:,:,jfile),ndof, 1)
        
!           print *, "Dot check: ", dot(atomvects(:,istate,jfile), atomvects(:,istate,jfile), numatomconfig(jfile)), dot(trans_atomvects(:,istate,2), trans_atomvects(:,istate,2), transnumatomconfig(2))


     enddo

!! FOR NOW: skipping raise/lower, clebsch gordan, etc.  total m_s must be less than calculation's restrictms.

!! COMBINE CONFIGS

        do config1=1,transnumatomconfig(1)
        do config2=1,transnumatomconfig(2)
           tempconfig(1:atom_ndof(ifile)) = transatomconfiglists(1:atom_ndof(ifile),config1,1)
           tempconfig(atom_ndof(ifile)+1:ndof) = transatomconfiglists(1:atom_ndof(jfile),config2,2)

!!           print *, "tempconfig is ", tempconfig(1:ndof)

           if (allowedconfig(tempconfig)) then
              phase=reorder(tempconfig)
              config3=getconfiguration(tempconfig)
!!              print *, "Allowed"

              do ipair=1,numpairstates
                 istate=state_pairs(1,ipair)
                 jstate=state_pairs(2,ipair)

                 a1=trans_atomvects(config1,istate,1)
                 a2=trans_atomvects(config2,jstate,2)
                 finalvects(config3,ipair,ipoint) = finalvects(config3,ipair,ipoint)+a1*a2*phase
              enddo
           endif
        enddo
        enddo

        print *, "DOTS BEFORE PROJECT (SHOULD BE 1 except if have same M and m_s)"
        print * 
        do ipair=1,numpairstates
           istate=state_pairs(1,ipair)
           jstate=state_pairs(2,ipair)

           print *, ipair, istate,jstate, atom_mval(istate,1), atom_mval(jstate,2), dot(finalvects(:,ipair,ipoint), finalvects(:,ipair,ipoint), numconfig)
        enddo

        print *, "DOTS AFTER PROJECT (SHOULD BE the same CG-coef squared except if have same M)"

        do ipair=1,numpairstates

           istate=state_pairs(1,ipair)
           jstate=state_pairs(2,ipair)

           call configspin_projectone(finalvects(:,ipair,ipoint))
           norm=real(dot(finalvects(:,ipair,ipoint), finalvects(:,ipair,ipoint), numconfig),8)
           print *, ipair,istate,jstate, atom_mval(istate,1), atom_mval(jstate,2),  norm
           finalvects(:,ipair,ipoint) = finalvects(:,ipair,ipoint) / sqrt(norm)
        enddo

        print *

!        print *, "tempstop"
!        stop




  end do  !! ipoint


!! NOW DIAGONALIZE TO GET MCSCF VECTORS for mctdhf run
!! requires mcscfnum to be set.  Needs atom spfs to fill
!! nspf slots.


  mcscfnum=numpairstates
  print *, "Getting MCSCF vects.  mcscfnum= ", mcscfnum

  call mcscfalloc(1)

  mcscfspfs(:,:,:) = finalspfs(:,:,:)

  allocate(ttvect(numconfig))

  ifile=1
  jfile=2
  


  
  print *, "doing mcscf"
  
  do ipoint=1,numatompoints
     
     mcwhichr=ipoint
     bondpoints=mcpoints(mcwhichr)
     
     call all_matel(finalspfs(:,:,ipoint), 1, -1.d0)

!! for full diag in space
     
!     do i=1,mcscfnum
!    need guessflag    call mcscf(ttvect, 1, i)
!     enddo

     print *, "get mats"

     pairham=0.d0
     pairovl=0.d0
     do jstate=1,numpairstates
        call configmultone(finalvects(:,jstate,ipoint),smallconfigvector,configpointer,1,0,0,1)
        do istate=1,numpairstates
           
           pairovl(istate,jstate)= dot(finalvects(:,istate,ipoint), finalvects(:,jstate,ipoint), numconfig)
           pairham(istate,jstate)= dot(finalvects(:,istate,ipoint), smallconfigvector, numconfig)
        enddo
     enddo

     print *,"call eigen"

!      SUBROUTINE ZHEGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
!     $                  LWORK, RWORK, INFO )

!      SUBROUTINE DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
!     $                  LWORK, INFO )

     lwork=numpairstates*3
     allocate(rwork(lwork), work(lwork))

     pairvects=pairham
     pairtemp=pairovl

#ifdef REALGO

     call DSYGV(1, 'V', 'U', numpairstates, pairvects, maxnumatomstates**2, pairtemp, maxnumatomstates**2, pairvals, work, lwork, info)

#else

#ifdef CNORMFLAG
     call ZGEGV(1, 'V', 'U', numpairstates, pairvects, maxnumatomstates**2, pairtemp, maxnumatomstates**2, pairvals, work, lwork, rwork, info)
#else

       !!!! well I should really do zgegv here too...  but will assume there is never ecs with chmctdh for the moment

     call ZHEGV(1, 'V', 'U', numpairstates, pairvects, maxnumatomstates**2, pairtemp, maxnumatomstates**2, pairvals, work, lwork, rwork, info)
#endif

#endif
     if (info/=0) then
        print *, "Info eigval", info
        stop
     endif

     deallocate(rwork,work)

     mcscfvals(:,ipoint)=pairvals(1:numpairstates)
     mcscfmvals(:,ipoint)=mrestrictval
     do istate=1,numpairstates
        mcscfvects(:,istate,ipoint)=0.d0
        do jstate=1,numpairstates
           mcscfvects(:,istate,ipoint) = mcscfvects(:,istate,ipoint) + pairvects(jstate,istate) * finalvects(:,jstate,ipoint)
        enddo
     enddo

!     do istate=1,numpairstates
!        call configmultone(mcscfvects(:,istate,ipoint),smallconfigvector,configpointer,1,0,0,1)
!        print *, "state 1, dot ", &
!             dot(mcscfvects(:,istate,ipoint), mcscfvects(:,istate,ipoint), numconfig), " eig ", &
!             dot(mcscfvects(:,istate,ipoint), smallconfigvector, numconfig), mcscfvals(istate,ipoint)
!     enddo
!     stop

  enddo
  
  print *, "done; output"
  call mcoutput()
  print *, "done"
  stop
  
  deallocate(ttvect)






     call commondealloc()

     call opdealloc()

     call myprojectdealloc2()
     call twoedealloc()
     call xdealloc()
     call walkdealloc()
     call configdealloc()
     call spinwalkdealloc()




end program mcscfassemble




