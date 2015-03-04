
!! MAIN PROGRAM IS HERE!   GOOD LUCK FOLLOWING THE LOGIC !!
  
  
#include "Definitions.INC"

  
program mcscf_matel
  use mpimod
  use parameters
  use configmod
  implicit none

  DATATYPE, allocatable :: projectorbs(:)  ! dummy

#ifdef MPIFLAG
  print *, "No, compile me without parallel option"
  stop
#endif

  pi=4.d0*atan(1.d0)
  call MPIstart()
  
  call getparams()

  skipflag=1
     call get_numconfig()
  
  totspfdim=nspf*spfsize; totadim=numconfig*numr;  psilength=totadim+totspfdim
  
  
!  if (1==1) then
!     astart=1;         aend=totadim;           spfstart=totadim+1;           spfend=psilength
!  else
!     spfstart=1;           spfend=totspfdim;           astart=totspfdim+1;           aend=psilength
!  endif
  
  OFLWR "Num config: ", numconfig
  
  if ((sparseconfigflag.eq.0).and.(numconfig.gt.1000).and.(nosparseforce.eq.0)) then
     write(mpifileptr, *) "You should really turn sparseconfigflag on, with ", &
          numconfig*numr, "configurations."
     CFLST
  endif
  CFL
  
     !! ALLOCATE.

     call opalloc()

     !! GET CONFIGURATONS.

  select case (psitype)
  case (1)
     call get_configlist()
     
        !!  GET WALKS (nonzero matrix elements and info)
     
     call walkalloc();              call walks()
     if (spinwalkflag==1) then
        call spinwalkalloc();              call spinwalks()
     endif
     call configalloc()
  case default
     OFLWR "PSITYPE ", psitype, "NOT SUPPORTED.";       CFLST
  end select

  OFLWR; CFL

  call curvealloc()

  call myprojectalloc()      !! Initialize coordinate-dependent arrays.

!!  call xalloc() !!   INITIALIZE XXX/YYY VECTORS!  KEEP THIS AS LAST ALLOC STATEMENT.


  allocate(projectorbs(totspfdim)); projectorbs=0d0


  
  if (atomflag==0) then
     call init_h2(projectorbs,nspf)        !! WITH POSSIBLE CORE EIGFUNCTS to init orbitals
     call commonalloc();        call commonset_h2()
  else
     call init_helium(projectorbs,nspf)      !! WITH POSSIBLE CORE EIGFUNCTS to init orbitals
     call commonalloc();        call commonset_atom()
  endif


  call mcmsub()

  deallocate(projectorbs)

  call mpistop()



end program mcscf_matel





!! SAVE INITIAL WAVE FUNCTION

module mcmmod
  implicit none
  DATATYPE, allocatable :: orig_spfs(:,:,:,:), orig_avectors(:,:,:),mcm_energies(:,:)
  integer :: mcmnumr(100)

  integer :: nummcm=0, nummcmfiles=0, mcmskip=0,mcmtot(40)=9999 !! MAX 40 FILES
  character(len=200) :: mcmfiles(100), mcmlabels(100)
  integer :: mcmtake(40,40)  !! MAX 40 FILES 40 STATES PER
!!  real*8 :: mcmzerothreshdip = -1, mcmzerothreshovl = -1d0
  real*8 :: mcmzerothreshdip = 1d-9, mcmzerothreshovl = 1d-9
  integer :: mcmoutmin=-1,mcminmax=1000
  integer :: mcmmultflag=0
  integer :: multny=2, multnx=2,multnz=2

end module mcmmod



subroutine mcmsub()
  use mcmmod
  use iglobalmod
  use parameters
  use mpimod
  use commonmod
  use yyymod
  implicit none

  integer :: tnumconfig, tndof, tnumerad, tlbig, tmbig, tnumr, dummyint(100), &
       tnspf, i, inummcm, jnummcm, ifile, index, iconfig, getconfiguration, ipoint, &
       myiostat,istate,jstate,jj,dirflag,jpoint,kk,ll,flag,numx,numy,numz,numovl,numt,j, &
       ii,k, it
  integer :: zdippairs(40,40) !! MAX 40 STATES TOTAL
  integer :: xdippairs(40,40) !! MAX 40 STATES TOTAL
  integer :: mdippairs(40,40) !! MAX 40 STATES TOTAL
  integer :: ydippairs(40,40) !! MAX 40 STATES TOTAL
  integer :: ovlpairs(40,40) !! MAX 40 STATES TOTAL


  real*8 :: rdummy,mcpoint,rdummy2, sum1,sum2,sum3, sum2b,eshift=0d0
  real*8, allocatable ::  ravec(:)
  complex*16 :: cdummy
  complex*16, allocatable ::  cavec(:)
  integer, allocatable :: dummyconfiglist(:,:)
  integer :: nocalc(1000)=0,nocalcstate(1000)=0,statefile(1000)=-1
  logical :: saveit
  DATATYPE, allocatable :: zdip(:,:,:), myovl(:,:,:), xdip(:,:,:), ydip(:,:,:), multdip(:,:,:,:,:,:),tempqspfs(:,:)
  character(len=200) :: mcmovlfile="MCMatel.Ovl.Dat"
  character(len=200) :: mcmxfile="MCMatel.XDipole.Dat"
  character(len=200) :: mcmyfile="MCMatel.YDipole.Dat"
  character(len=200) :: mcmzfile="MCMatel.ZDipole.Dat"
  character(len=200) :: mcmmfile="MCMatel.MDipole.Dat"
  character(len=200) :: mcmfile="MCMatel.Dat"

  EXTERNAL :: mult_zdipole,mult_ydipole,mult_xdipole, multdipoperate

  NAMELIST/mcminp/ nummcmfiles,mcmfiles,mcmskip,mcmfile,mcmzfile,mcmxfile,mcmyfile,mcmovlfile,mcmtake,mcmtot, mcmzerothreshovl, mcmzerothreshdip,mcminmax,mcmoutmin, mcmmultflag,multny,multnx,multnz,nocalc, eshift,mcmlabels

  do i=1,40
     mcmtake(i,:)=i
  enddo

  open(971,file=inpfile, status="old", iostat=myiostat)
  if (myiostat/=0) then
     OFLWR "No Input.Inp found, iostat=",myiostat; CFLST
  endif
  read(971,nml=mcminp,iostat=myiostat)
  if (myiostat/=0) then
     OFLWR "No mcminp namelist in Input.Inp, iostat=",myiostat; CFLST
  endif
  if (nummcmfiles.le.0) then
     OFLWR "Need to set nummcmfiles and mcmfiles in namelist mcminp"; CFLST
  endif

  mcmlabels(1:nummcmfiles)=mcmfiles(:nummcmfiles)

  close(971)
  open(971,file=inpfile, status="old", iostat=myiostat)
  read(971,nml=mcminp,iostat=myiostat)
  close(971)



  nummcm=0

  if (whichprojflux.eq.2) then
     call openfile(); write(mpifileptr,*) "Checking COMPLEX mcm files "; call closefile()
  else if (whichprojflux.eq.1) then
     call openfile(); write(mpifileptr,*) "Checking REAL mcm files "; call closefile()
  else
     call openfile(); write(mpifileptr,*) "Set whichprojflux to 1 or 2 depending on ";
     call openfile(); write(mpifileptr,*) " whether or not you have REAL or COMPLEX mcm.mcscf.bin file", ifile; call closefile()
  endif

  mcmnumr(:)=0



  do ifile=1,nummcmfiles
     if (nocalc(ifile).eq.0) then
        open(9954,file=mcmfiles(ifile),status="unknown",form="unformatted",position="rewind")
        read(9954) rdummy,rdummy,tndof  
        read(9954) i,i,i,tnumr
        read(9954) i,i,i
        read(9954) tnumconfig,tnspf
        read(9954) dummyint(1:tnspf) 
        read(9954) inummcm,mcmnumr(nummcm+1),rdummy,rdummy,tnumerad,tlbig,tmbig
        close(9954)
     else
        if (nummcm.eq.0) then
           OFLWR "can't do this kloodge with first guy"; CFLST
        endif
        mcmnumr(nummcm+1)=mcmnumr(nummcm)+1
     endif
     OFLWR "MCMNUMR ", nummcm+1,mcmnumr(nummcm+1); CFL

!! NOW 092512   mcmtake can be specified
!! if mcmtot is still 9999
!! then take all and check that mcmtake is not changed
!! otherwise check mcmtake(mcmtot)

     if (mcmtot(ifile).eq.9999) then
        mcmtot(ifile)=inummcm
        if (mcmtake(inummcm,ifile).ne.inummcm) then
           OFLWR "Error, mcmtot has not been specified and mcmtake has"; CFLST
        endif
     else
        if (mcmtake(mcmtot(ifile),ifile).gt.inummcm) then
           OFLWR "Error, mcmtake(mcmtot) is greater than # states on file", mcmtake(mcmtot(ifile),ifile),mcmtot(ifile),inummcm; CFLST
        endif
     endif

     nummcm=nummcm+mcmtot(ifile)
     mcmnumr(nummcm-mcmtot(ifile)+1:nummcm)=mcmnumr(nummcm-mcmtot(ifile)+1)-1


     call openfile()
     write(mpifileptr,*) "MCSCF Data on file, current:"
     write(mpifileptr,*) "lbig:          ", tlbig, lbig
     write(mpifileptr,*) "mbig:          ", tmbig, mbig
     write(mpifileptr,*) "nspf:          ", tnspf, nspf
     write(mpifileptr,*) "ndof:          ", tndof, ndof
     write(mpifileptr,*) "numconfig:     ", tnumconfig, numconfig
     write(mpifileptr,*) "numr:          ", mcmnumr(nummcm),numr
     if (dirflag.eq.(-1)) then
        write(mpifileptr,*) "  ... file is backwards."
     endif
     call closefile()
     
     if (tlbig.ne.lbig .or. tmbig.gt.mbig .or.tnspf.ne.nspf .or. tndof.ne.ndof .or. tnumconfig.gt.numconfig)  then
        call openfile()
        write(mpifileptr,*) " ...MCSCF vectors on file do not agree with current calculation. "
        call closefile();    call mpistop()
     endif
     
     if(tnumerad.gt.numerad) then
        call openfile()
        write(mpifileptr,*) "numerad on disk too big. on file, current:",tnumerad,numerad
        call closefile();    call mpistop()
     endif
     
     if(mcmnumr(nummcm).gt.numr) then
        call openfile()
        write(mpifileptr,*) "numr on disk too big. on file, current:",mcmnumr(nummcm),numr
        call closefile();    call mpistop()
     endif
     
     if(tnumerad.lt.numerad) then
        call openfile()
        write(mpifileptr,*) "numerad on disk smaller than current, adjusting. on file, current:",tnumerad,numerad
        call closefile()
     endif
  enddo


!!  print *, "MCMNUMR ", mcmnumr

  allocate(orig_spfs(spfsize,nspf,nummcm,0:numr),orig_avectors(numconfig,nummcm,0:numr),mcm_energies(nummcm,0:numr))

  orig_spfs=0d0
  orig_avectors=0d0
  mcm_energies=0d0

  call openfile
  if(whichprojflux.eq.1) then
    write(mpifileptr,*) "wavefunction on disk is REAL"
  else if(whichprojflux.eq.2) then
    write(mpifileptr,*) "wavefunction on disk is COMPLEX"
  endif
  call closefile


  jnummcm=0


  do ifile=1,nummcmfiles

     if (nocalc(ifile).eq.0) then
        open(9954,file=mcmfiles(ifile),status="unknown",form="unformatted",position="rewind")
        read(9954) rdummy,rdummy,tndof  
        read(9954) i,i,i,tnumr
        read(9954) i,i,i
        read(9954) tnumconfig,tnspf
        read(9954) dummyint(1:tnspf) 
        read(9954) inummcm,i,rdummy,rdummy2,tnumerad,tlbig,tmbig
     endif

     ! these are first and last points  - check for backwards
     if (rdummy.lt.rdummy2)  then   ! usual case
        dirflag=1
     else
        dirflag=(-1)
     endif

     print *, "points on file ", ifile, " are going in direction ", dirflag

     OFLWR "Reading mcm file ", ifile, "  Number of states on disk=",inummcm, "total used=",mcmtot(ifile); CFL

     if (tnumr.ne.1) then
        OFLWR "Error; it appears you have a vibrational mcscf file in ", mcmfiles(ifile)
        WRFL  "   because numr = ", tnumr, " on file."; CFLST
     endif


     allocate(dummyconfiglist(ndof,tnumconfig) )
     if (nocalc(ifile).eq.0) then
        read(9954) dummyconfiglist(:,:)
     endif

     do jpoint=0,mcmnumr(jnummcm+1)
        if (dirflag.eq.1) then
           ipoint=jpoint
        else
           !! will skip first point, jpoint=0
           ipoint=numr-jpoint+1
        endif
        if (jpoint.eq.0.and.dirflag.eq.(-1)) then
           saveit=.false.
        else
           saveit=.true.
        endif
        if (nocalc(ifile).eq.0) then
           read(9954) mcpoint
           read(9954) dummyint(1:inummcm) !! dummy mvalue
        endif

!! THIS ONE  DONE!! xx

        allocate(tempqspfs(spfsize,nspf))

        if (whichprojflux.eq.1) then
           call one_spf_real_read(9954,tnumerad,lbig+1,mbig,tnspf,tempqspfs,nspf,myiostat)
        else if (whichprojflux.eq.2) then
           call one_spf_complex_read(9954,tnumerad,lbig+1,mbig,tnspf,tempqspfs,nspf,myiostat)
        endif
        call checkiostat(myiostat)

        if (saveit) then
           orig_spfs(:,:,jnummcm+1:jnummcm+mcmtot(ifile),ipoint)=0d0
        endif

        kk=0
        do i=1,inummcm
           flag=0
           do ll=1,mcmtot(ifile)
              if (mcmtake(ll,ifile).eq.i) then
                 flag=1
                 exit
              endif
           enddo
           if (flag.eq.1) then
              kk=kk+1
              statefile(kk+jnummcm)=ifile
              nocalcstate(kk+jnummcm)=nocalc(ifile)
           endif
           if (saveit) then

              orig_spfs(:,:,jnummcm+kk,ipoint)=tempqspfs(:,:)

           endif

           if (nocalc(ifile).eq.0) then
              if (whichprojflux.eq.1) then
                 read(9954) rdummy ; if (flag.eq.1.and.saveit) mcm_energies(kk+jnummcm,ipoint)=rdummy
                 read(9954) ravec(:)
              else if (whichprojflux.eq.2) then
                 read(9954) cdummy ; if (flag.eq.1.and.saveit) mcm_energies(kk+jnummcm,ipoint)=cdummy
                 read(9954) cavec(:)
              endif

              if (saveit) then
                 if (flag.eq.1) then
                    do iconfig=1,tnumconfig
                       index=getconfiguration(dummyconfiglist(:,iconfig))
                       if (index.lt.1.or.index.gt.numconfig) then
                          print *, "whhoopsy error", index, numconfig
                          stop
                       endif
                       if (whichprojflux.eq.2) then
                          orig_avectors(index,kk+jnummcm,ipoint)=ravec(iconfig)
                       else if (whichprojflux.eq.2) then
                          orig_avectors(index,kk+jnummcm,ipoint)=cavec(iconfig)
                       endif
                    enddo
                 endif
              endif
           else
              mcm_energies(kk+jnummcm,ipoint)=0d0
              orig_avectors(:,kk+jnummcm,ipoint)=0d0
           endif
        enddo
        deallocate(tempqspfs)
     enddo     !jpoint

     deallocate(dummyconfiglist)
     deallocate(ravec, cavec)

     jnummcm=jnummcm+mcmtot(ifile)
     close(9954)

  enddo


  OFLWR "Have ", nummcm, " vectors"; CFL

  allocate(zdip(nummcm,nummcm,0:numr), myovl(nummcm,nummcm,0:numr), xdip(nummcm,nummcm,0:numr), ydip(nummcm,nummcm,0:numr))
  allocate(multdip(nummcm,nummcm,0:multnx,0:multny,0:multnz,0:numr))
  zdip=0d0; myovl=0d0; ydip=0d0; xdip=0d0

  numz=0;  numy=0; numx=0;  numovl=0; numt=0


  if (mcmskip.eq.0) then

     do istate=1,min(mcminmax,nummcm)
        jstate=istate
        
        print *, " istate ", istate

!           print *, "istate,jstate ", istate, jstate,nocalcstate(istate),nocalcstate(jstate), "UU"
!           print *, mcmnumr(istate),mcmnumr(jstate), "XsX"
        
        do ipoint=1,min(mcmnumr(istate),mcmnumr(jstate))

           if (nocalcstate(istate).ne.0 .or. nocalcstate(jstate).ne.0) then
              myovl(istate,jstate,ipoint)=-999.999d0
           else
              call autocorrelate_one(orig_avectors(:,istate,ipoint),orig_spfs(:,:,istate,ipoint),orig_spfs(:,:,jstate,ipoint),orig_avectors(:,jstate,ipoint), myovl(istate,jstate,ipoint),1)
           endif

           call mcmmultmbig(mbig)

           if (nocalcstate(istate).ne.0 .or. nocalcstate(jstate).ne.0) then
                 xdip(istate,jstate,ipoint)=-999.999d0
                 ydip(istate,jstate,ipoint)=-999.999d0
                 zdip(istate,jstate,ipoint)=-999.999d0
                 multdip(istate,jstate,:,:,:,ipoint)=-999.999d0
           else

              call onee_generalmatel(orig_spfs(:,:,istate,ipoint),orig_avectors(:,istate,ipoint), &
                   orig_spfs(:,:,jstate,ipoint),orig_avectors(:,jstate,ipoint), &
                   mult_zdipole,zdip(istate,jstate,ipoint))
              
              call onee_generalmatel(orig_spfs(:,:,istate,ipoint),orig_avectors(:,istate,ipoint), &
                   orig_spfs(:,:,jstate,ipoint),orig_avectors(:,jstate,ipoint), &
                   mult_xdipole,xdip(istate,jstate,ipoint))
              
              call onee_generalmatel(orig_spfs(:,:,istate,ipoint),orig_avectors(:,istate,ipoint), &
                   orig_spfs(:,:,jstate,ipoint),orig_avectors(:,jstate,ipoint), &
                   mult_ydipole,ydip(istate,jstate,ipoint))
           
              if (mcmmultflag.ne.0) then
                 print *, "go mcmmult"
                 do ii=0,multnx
                    do jj=0,multny
                       do kk=0,multnz

!! this is all wasting a lot of effort - 

                          call mcmmultinit(ii,jj,kk)

!! xmbig set in multdipoperate

           call onee_generalmatel(orig_spfs(:,:,istate,ipoint),orig_avectors(:,istate,ipoint), &
                orig_spfs(:,:,jstate,ipoint),orig_avectors(:,jstate,ipoint), &
                multdipoperate,multdip(istate,jstate,ii,jj,kk,ipoint))

!!           print *, "X--XXX  ", istate,jstate,ipoint

!!    print *, abs(multdip(istate,jstate,ii,jj,kk,ipoint)), "QQQ"
           

                 
                       enddo
                    enddo
                 enddo
              endif
           endif
        enddo
        
        do jstate=max(istate+1,mcmoutmin),nummcm
           
!           print *, "istate,jstate ", istate, jstate,nocalcstate(istate),nocalcstate(jstate), "UU"
!           print *, mcmnumr(istate),mcmnumr(jstate), "XsX"
           
           do ipoint=1,min(mcmnumr(istate),mcmnumr(jstate))
              
              if (nocalcstate(istate).ne.0 .or. nocalcstate(jstate).ne.0) then
                 myovl(istate,jstate,ipoint)=-999.999d0
              else
                 call autocorrelate_one(orig_avectors(:,istate,ipoint),orig_spfs(:,:,istate,ipoint),orig_spfs(:,:,jstate,ipoint),orig_avectors(:,jstate,ipoint), myovl(istate,jstate,ipoint),1)
              endif
              call mcmmultmbig(mbig)

              if (nocalcstate(istate).ne.0 .or. nocalcstate(jstate).ne.0) then

                 xdip(istate,jstate,ipoint)=-999.999d0
                 ydip(istate,jstate,ipoint)=-999.999d0
                 zdip(istate,jstate,ipoint)=-999.999d0
                 multdip(istate,jstate,:,:,:,ipoint)=-999.999d0

              else
                 
                 call onee_generalmatel(orig_spfs(:,:,istate,ipoint),orig_avectors(:,istate,ipoint), &
                      orig_spfs(:,:,jstate,ipoint),orig_avectors(:,jstate,ipoint), &
                      mult_zdipole,zdip(istate,jstate,ipoint))
                 
                 call onee_generalmatel(orig_spfs(:,:,istate,ipoint),orig_avectors(:,istate,ipoint), &
                      orig_spfs(:,:,jstate,ipoint),orig_avectors(:,jstate,ipoint), &
                      mult_xdipole,xdip(istate,jstate,ipoint))
                 
                 call onee_generalmatel(orig_spfs(:,:,istate,ipoint),orig_avectors(:,istate,ipoint), &
                      orig_spfs(:,:,jstate,ipoint),orig_avectors(:,jstate,ipoint), &
                      mult_ydipole,ydip(istate,jstate,ipoint))

                 if (mcmmultflag.ne.0) then
                    print *, "go mcmmult"
                    do ii=0,multnx
                       do jj=0,multny
                          do kk=0,multnz
                             
!! this is all wasting a lot of effort - 

                       call mcmmultinit(ii,jj,kk)

           call onee_generalmatel(orig_spfs(:,:,istate,ipoint),orig_avectors(:,istate,ipoint), &
                orig_spfs(:,:,jstate,ipoint),orig_avectors(:,jstate,ipoint), &
                multdipoperate,multdip(istate,jstate,ii,jj,kk,ipoint))
                 
                          enddo
                       enddo
                    enddo
                 endif
              endif

!           print *, "XXXX  ", istate,jstate,ipoint
!     print *, multdip(istate,jstate,:,:,:,ipoint)
!     stop
              
           enddo
           
        enddo
        
     enddo


     do istate=1,min(mcminmax,nummcm)
        do jstate=1,istate-1
           zdip(istate,jstate,:)=CONJG(zdip(jstate,istate,:))
           xdip(istate,jstate,:)=CONJG(xdip(jstate,istate,:))
           ydip(istate,jstate,:)=CONJG(ydip(jstate,istate,:))
           if (mcmmultflag.ne.0) then
              multdip(istate,jstate,:,:,:,:)=CONJG(multdip(jstate,istate,:,:,:,:))
           endif
           myovl(istate,jstate,:)=CONJG(myovl(jstate,istate,:))
        enddo
     enddo



!
     !! GET significant off diagonal elements
     
     do istate=1,min(nummcm,mcminmax)
        numz=numz+1;numy=numy+1;numx=numx+1;numovl=numovl+1
        zdippairs(1,numz)=istate; zdippairs(2,numz)=istate
        xdippairs(1,numx)=istate; xdippairs(2,numx)=istate
        mdippairs(1,numx)=istate; mdippairs(2,numx)=istate
        ydippairs(1,numy)=istate; ydippairs(2,numy)=istate
        ovlpairs(1,numovl)=istate; ovlpairs(2,numovl)=istate
        
        do jstate=max(mcmoutmin,istate+1),nummcm
           sum1=0d0; sum2=0d0;sum2b=0d0; sum3=0d0
           do ipoint=1,min(mcmnumr(istate),mcmnumr(jstate))
              sum1=sum1 + abs( zdip(istate,jstate,ipoint) )
              sum2=sum2 + abs( xdip(istate,jstate,ipoint) )
              sum2b=sum2b + abs( ydip(istate,jstate,ipoint) )
              if (istate.eq.jstate) then
                 sum3=sum3 + abs( myovl(istate,jstate,ipoint)-1d0 )
              else
                 sum3=sum3 + abs( myovl(istate,jstate,ipoint) )
              endif
           enddo
           if (sum1.gt.mcmzerothreshdip*min(mcmnumr(istate),mcmnumr(jstate))) then
              numz=numz+1
              zdippairs(1,numz)=istate; zdippairs(2,numz)=jstate
           endif
           if (sum2.gt.mcmzerothreshdip*min(mcmnumr(istate),mcmnumr(jstate))) then
              numx=numx+1
              xdippairs(1,numx)=istate; xdippairs(2,numx)=jstate
           endif
           if (sum2b.gt.mcmzerothreshdip*min(mcmnumr(istate),mcmnumr(jstate))) then
              numy=numy+1
              ydippairs(1,numy)=istate; ydippairs(2,numy)=jstate
           endif
           if (sum3.gt.mcmzerothreshovl*min(mcmnumr(istate),mcmnumr(jstate))) then
              numovl=numovl+1
              ovlpairs(1,numovl)=istate; ovlpairs(2,numovl)=jstate
           endif
           numt=numt+1
           mdippairs(1,numt)=istate; mdippairs(2,numt)=jstate
        enddo
     enddo
     
     open(6690,file=mcmovlfile, status="unknown")
     open(6691,file=mcmxfile, status="unknown")
     open(66911,file=mcmyfile, status="unknown")
     open(6692,file=mcmzfile, status="unknown")
     if (numr.ne.1) then
        write(6690,'(A10,1000(A8,I5,I5,A8))') " R ", (" ",ovlpairs(1,jj),ovlpairs(2,jj)," ", jj=1,numovl)
        write(6691,'(A10,1000(A8,I5,I5,A8))') " R ", (" ",xdippairs(1,jj),xdippairs(2,jj)," ", jj=1,numx)
        write(66911,'(A10,1000(A8,I5,I5,A8))') " R ", (" ",ydippairs(1,jj),ydippairs(2,jj)," ", jj=1,numy)
        write(6692,'(A10,1000(A8,I5,I5,A8))') " R ", (" ",zdippairs(1,jj),zdippairs(2,jj)," ", jj=1,numz)
     endif
  endif
  
  open(6693,file=mcmfile, status="unknown")
!  open(66930,file=mcmfile//".eV", status="unknown")
  
  if (numr.eq.1) then
     ipoint=1
     write(6693,*)
     write(6693,'(A10, F10.5)') "  R= ", real(bondpoints(ipoint),8)
     write(6693,*) "  ----------------  "
     write(6693,'(2000F13.6)')  mcm_energies(:,ipoint)

!     write(66930,*)
!     write(66930,'(A10, F10.5)') "  R= ", real(bondpoints(ipoint),8)
!     write(66930,*) "  ----------------  "
!     write(66930,'(2F13.6)')  27.2114*(mcm_energies(:,ipoint)-eshift)

     if (mcmskip.eq.0) then

        write(6690,*)
        write(6690,'(A10, F10.5)') "  R= ", real(bondpoints(ipoint),8)
        write(6690,*) "  ----------------  "
        
        write(6691,*)
        write(6691,'(A10, F10.5)') "  R= ", real(bondpoints(ipoint),8)
        write(6691,*) "  ----------------  "
        
        write(66911,*)
        write(66911,'(A10, F10.5)') "  R= ", real(bondpoints(ipoint),8)
        write(66911,*) "  ----------------  "
        
        write(6692,*)
        write(6692,'(A10, F10.5)') "  R= ", real(bondpoints(ipoint),8)
        write(6692,*) "  ----------------  "

        write(6690,'(2I5,2F13.6)') (ovlpairs(1,jj),ovlpairs(2,jj),myovl(ovlpairs(1,jj),ovlpairs(2,jj),ipoint),jj=1,numovl)

!! 051513  Ack!  had bug like this for all
!!        write(6691,'(2I5,2F13.6)') (xdippairs(1,jj),ovlpairs(2,jj),xdip(xdippairs(1,jj),xdippairs(2,jj),ipoint),jj=1,numx)

        write(6691,'(2I5,2F13.6,2F10.2,3F14.6,A5,2A60)') (xdippairs(1,jj),xdippairs(2,jj),&
             xdip(xdippairs(1,jj),xdippairs(2,jj),ipoint), &
             27.2114*(eshift-min(real(mcm_energies(xdippairs(1,jj),ipoint)),real(mcm_energies(xdippairs(2,jj),ipoint)))), &
             27.2114*abs(mcm_energies(xdippairs(1,jj),ipoint)-mcm_energies(xdippairs(2,jj),ipoint)), &
             eshift,real(mcm_energies(xdippairs(1,jj),ipoint)),real(mcm_energies(xdippairs(2,jj),ipoint)), &
             "   ",&
             mcmfiles(statefile(xdippairs(1,jj))),mcmfiles(statefile(xdippairs(2,jj))), &
             jj=1,numx)

        write(66911,'(2I5,2F13.6,2F10.2,3F14.6,A5,2A60)') (ydippairs(1,jj),ydippairs(2,jj),&
             ydip(ydippairs(1,jj),ydippairs(2,jj),ipoint), &
             27.2114*(eshift-min(real(mcm_energies(ydippairs(1,jj),ipoint)),real(mcm_energies(ydippairs(2,jj),ipoint)))), &
             27.2114*abs(mcm_energies(ydippairs(1,jj),ipoint)-mcm_energies(ydippairs(2,jj),ipoint)), &
             eshift,real(mcm_energies(ydippairs(1,jj),ipoint)),real(mcm_energies(ydippairs(2,jj),ipoint)), &
             "   ",&
             mcmfiles(statefile(ydippairs(1,jj))),mcmfiles(statefile(ydippairs(2,jj))), &
             jj=1,numy)

        write(6692,'(2I5,2F13.6,2F10.2,3F14.6,A5,2A60)') (zdippairs(1,jj),zdippairs(2,jj),&
             zdip(zdippairs(1,jj),zdippairs(2,jj),ipoint), &
             27.2114*(eshift-min(real(mcm_energies(zdippairs(1,jj),ipoint)),real(mcm_energies(zdippairs(2,jj),ipoint)))), &
             27.2114*abs(mcm_energies(zdippairs(1,jj),ipoint)-mcm_energies(zdippairs(2,jj),ipoint)), &
             eshift,real(mcm_energies(zdippairs(1,jj),ipoint)),real(mcm_energies(zdippairs(2,jj),ipoint)), &
             "   ",&
             mcmfiles(statefile(zdippairs(1,jj))),mcmfiles(statefile(zdippairs(2,jj))), &
             jj=1,numz)

        if (mcmmultflag.ne.0) then
           open(7692,file=mcmmfile, status="unknown")        

           do k=0,multnz
              do j=0,multny
                 do i=0,multnx
                    write(7692,*)
                    write(7692,*) "# R= ",real(bondpoints(ipoint),8),"  Multipole moment x^l y^m z^n , lmn=",i,j,k
                    write(7692,*) "#  istate, jstate, matrix element:"
                    write(7692,*)
!!                    write(7692,'(2I5,2F13.6)') (mdippairs(1,jj),mdippairs(2,jj),multdip(mdippairs(1,jj),mdippairs(2,jj),i,j,k,ipoint),jj=1,numt)

!                    write(7692,'(2I5,2E13.4)') (mdippairs(1,jj),mdippairs(2,jj),multdip(mdippairs(1,jj),mdippairs(2,jj),i,j,k,ipoint),jj=1,numt)

                    do jj=1,numt
                       if (abs(multdip(mdippairs(1,jj),mdippairs(2,jj),i,j,k,ipoint)).gt.mcmzerothreshdip) then
                          write(7692,'(2I5,2E13.4)') mdippairs(1,jj),mdippairs(2,jj),multdip(mdippairs(1,jj),mdippairs(2,jj),i,j,k,ipoint)
                       endif
                    enddo

                 enddo
              enddo
           enddo
           close(7692)
        endif
     endif

     if (1==0) then

        open(1534,file="MCMatrix.ZDipole.Dat",status="unknown")
        open(1535,file="MCMatrix.XDipole.Dat",status="unknown")
        open(15351,file="MCMatrix.YDipole.Dat",status="unknown")
        open(1536,file="MCMatrix.Ovl.Dat",status="unknown")
        write(1534,*)  "zdipole={"
        write(1535,*)  "xdipole={"
        write(15351,*)  "ydipole={"
        write(1536,*)  "ovl={"
        do i=1,nummcm
           write(1534,'(A2$)') "{"
           write(1535,'(A2$)') "{"
           write(15351,'(A2$)') "{"
           write(1536,'(A2$)') "{"
           do j=1,nummcm
              if (j.eq.nummcm) then
                 if (i.eq.nummcm) then
                    write(1534,'(A1,F15.10,A3,F15.10,A4)') &
                         "(",real(zdip(i,j,1)),"+I*",imag(zdip(i,j,1)),")}};"
                    write(1535,'(A1,F15.10,A3,F15.10,A4)') &
                         "(",real(xdip(i,j,1)),"+I*",imag(xdip(i,j,1)),")}};"
                    write(15351,'(A1,F15.10,A3,F15.10,A4)') &
                         "(",real(ydip(i,j,1)),"+I*",imag(ydip(i,j,1)),")}};"
                    write(1536,'(A1,F15.10,A3,F15.10,A4)') &
                         "(",real(myovl(i,j,1)),"+I*",imag(myovl(i,j,1)),")}};"
                 else
                    write(1534,'(A1,F15.10,A3,F15.10,A4)') &
                         "(",real(zdip(i,j,1)),"+I*",imag(zdip(i,j,1)),")},"
                    write(1535,'(A1,F15.10,A3,F15.10,A4)') &
                         "(",real(xdip(i,j,1)),"+I*",imag(xdip(i,j,1)),")},"
                    write(15351,'(A1,F15.10,A3,F15.10,A4)') &
                         "(",real(ydip(i,j,1)),"+I*",imag(ydip(i,j,1)),")},"
                    write(1536,'(A1,F15.10,A3,F15.10,A4)') &
                         "(",real(myovl(i,j,1)),"+I*",imag(myovl(i,j,1)),")},"

!!                 write(1534,'(2F15.10,A2)') zdip(i,j,1),"},"
                 endif
              else
                 write(1534,'(A1,F15.10,A3,F15.10,A2$)') &
                      "(",real(zdip(i,j,1)),"+I*",imag(zdip(i,j,1)),"), "
                 write(1535,'(A1,F15.10,A3,F15.10,A2$)') &
                      "(",real(xdip(i,j,1)),"+I*",imag(xdip(i,j,1)),"),"
                 write(15351,'(A1,F15.10,A3,F15.10,A2$)') &
                      "(",real(ydip(i,j,1)),"+I*",imag(ydip(i,j,1)),"),"
                 write(1536,'(A1,F15.10,A3,F15.10,A2$)') &
                      "(",real(myovl(i,j,1)),"+I*",imag(myovl(i,j,1)),"),"

!!              write(1534,'(2F15.10,A2$)') zdip(i,j,1),", "
                 
              endif
           enddo
        enddo
     endif
  else

     do ipoint=1,numr
        write(6693,'(F10.5, 1000F13.6)') real(bondpoints(ipoint),8), mcm_energies(:,ipoint)
!        write(66930,'(F10.5, 1000F13.6)') real(bondpoints(ipoint),8), 27.2114*(mcm_energies(:,ipoint)-eshift)
        if (mcmskip.eq.0) then
           write(6690,'(F10.5, 1000F13.6)') real(bondpoints(ipoint),8), (myovl(ovlpairs(1,jj),ovlpairs(2,jj),ipoint),jj=1,numovl)
           write(6691,'(F10.5, 1000F13.6)') real(bondpoints(ipoint),8), (xdip(xdippairs(1,jj),xdippairs(2,jj),ipoint),jj=1,numx)
           write(66911,'(F10.5, 1000F13.6)') real(bondpoints(ipoint),8), (ydip(ydippairs(1,jj),ydippairs(2,jj),ipoint),jj=1,numy)
           write(6692,'(F10.5, 1000F13.6)') real(bondpoints(ipoint),8), (zdip(zdippairs(1,jj),zdippairs(2,jj),ipoint),jj=1,numz)
        endif
     enddo
  endif


  close(6690)
  close(6691)
  close(66911)
  close(6692)
  close(6693)
!  close(66930)
  close(1534)
  close(1535)
  close(15351)
  close(1536)

  deallocate(zdip, myovl, xdip, ydip, multdip)
  deallocate(orig_spfs,orig_avectors,mcm_energies)

  print *, "DONE MCMSUB"

!!$
!!$  do istate=1,nummcm
!!$     write(*,'(100F10.5)') myovl(:,istate)
!!$  enddo
!!$
!!$  print *
!!$
!!$  do istate=1,nummcm
!!$     write(*,'(100F10.5)') zdip(:,istate)
!!$  enddo
!!$
!!$  print *
!!$
!!$  do istate=1,nummcm
!!$     write(*,'(100F10.5)') xydip(:,istate)
!!$  enddo

end subroutine mcmsub




module mcmmultmod
  integer :: ix=0,iy=0,iz=0,  xmbig=-1
end module mcmmultmod





!! wasting a LOT of effort by calling this repeatedly each time

subroutine multdipoperate(inspf,outspf)
  use mcmmultmod
  use parameters
  use opmod
  implicit none

  integer :: i,ii

  DATATYPE :: inspf(numerad,lbig+1,-mbig:mbig)
  DATATYPE :: outspf(numerad,lbig+1,-mbig:mbig)
  DATATYPE, allocatable :: tempspf(:,:,:),tempspf2(:,:,:)

!! THIS ONE  call mcmmultmbig(mbig+abs(ix)+abs(iy)+1)  !! sets xmbig

!! TEMP TEMP TEMP

 call mcmmultmbig(mbig+abs(ix)+abs(iy)-1)  !! sets xmbig.  -1 looks right (observation)

!  print *, "XMGBIG IS ", xmbig, mbig, "XX"; stop


  allocate(tempspf(numerad,lbig+1,-xmbig:xmbig), tempspf2(numerad,lbig+1,-xmbig:xmbig))

  tempspf=0d0; tempspf2=0d0
! debugging
!  tempspf(:,:,-mbig:mbig) = inspf(:,:,-mbig:mbig)

  tempspf(:,:,-int(min(xmbig,mbig)):int(min(xmbig,mbig))) = inspf(:,:,-int(min(xmbig,mbig)):int(min(xmbig,mbig)))

  do ii=1,iz
     call mult_zdipole0(tempspf,tempspf2,xmbig)
     tempspf=tempspf2
  enddo
  do ii=1,iy
     call mult_ydipole0(tempspf,tempspf2,xmbig)
     tempspf=tempspf2
  enddo
  do ii=1,ix
     call mult_xdipole0(tempspf,tempspf2,xmbig)
     tempspf=tempspf2
  enddo

!debug
!  outspf(:,:,-mbig:mbig)=tempspf(:,:,-mbig:mbig)
  outspf(:,:,-int(min(xmbig,mbig)):int(min(xmbig,mbig))) = tempspf(:,:,-int(min(xmbig,mbig)):int(min(xmbig,mbig)))
  deallocate(tempspf,tempspf2)

  call mcmmultmbig(mbig)

end subroutine multdipoperate



subroutine mcmmultmbig(inmbig)
  use mcmmultmod
  implicit none
  integer :: inmbig
  xmbig=inmbig
end subroutine mcmmultmbig

subroutine mcmmultinit(iix,iiy,iiz)
  use mcmmultmod
  implicit none
  integer :: iix,iiy,iiz
  ix=iix
  iy=iiy
  iz=iiz
end subroutine mcmmultinit


