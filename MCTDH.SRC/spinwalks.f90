

!! SUBROUTINE FOR "WALKS" (WHICH CONFIGURATIONS CONNECT TO WHICH) FOR SPARSE SPIN PROJECTOR.

#include "Definitions.INC"

module spinwalkinternal
  implicit none
  integer, allocatable :: numspinwalks(:),spinwalkdirphase(:,:),spinwalk(:,:),&
       unpaired(:,:),msvalue(:),numunpaired(:)
  real*8, allocatable :: configspinmatel(:,:)
  integer ::   maxspinwalks=0
end module


subroutine spinwalkinternal_dealloc()
  use parameters
  use spinwalkinternal
  implicit none
  deallocate(unpaired,numunpaired,msvalue,numspinwalks,spinwalk, spinwalkdirphase, configspinmatel)
end subroutine spinwalkinternal_dealloc



subroutine spinwalkinit()
  use parameters
  use spinwalkinternal
  use spinwalkmod
  use mpimod
  implicit none

  OFLWR "Go spinwalks. "; CFL
  
  allocate(numspinsets(nprocs),numspindfsets(nprocs))

  allocate(unpaired(numelec,numconfig), numunpaired(numconfig), msvalue(numconfig), numspinwalks(numconfig))

  call getnumspinwalks()

  allocate(spinwalk(maxspinwalks,numconfig),spinwalkdirphase(maxspinwalks,numconfig))
  allocate(configspinmatel(maxspinwalks+1,numconfig))

  allocate(spinsetsize(maxconfigsperproc,nprocs),spinsetrank(maxconfigsperproc,nprocs))

end subroutine spinwalkinit


subroutine spinwalkdealloc()
  use parameters
  use mpimod
  use spinwalkmod
  implicit none
  integer :: i,iproc

  deallocate(spinsets, spinsetsize, spinsetrank)

  do iproc=1,nprocs
     do i=1,numspinsets(iproc)
        deallocate(spinsetprojector(i,iproc)%mat,spinsetprojector(i,iproc)%vects)
     enddo
  enddo
  deallocate(spinsetprojector)

end subroutine spinwalkdealloc



subroutine spinwalks()
  use spinwalkmod
  use spinwalkinternal
  use configmod
  use parameters
  use mpimod !! nprocs
  use aarrmod
  implicit none

  integer ::     config1, config2, dirphase,  idof, jdof,iwalk , thisconfig(ndof),  &
       thatconfig(ndof), reorder,  getconfiguration, ii, jj, firstspin, secondspin, iproc
  logical :: allowedconfig !! extraconfig

  spinwalk=0;     spinwalkdirphase=0


  OFLWR "Calculating spin walks.";  CFL

  do iproc=1,nprocs

     do config1=allbotconfigs(iproc),alltopconfigs(iproc)

        iwalk=0
        do ii=1,numunpaired(config1)
           thisconfig=configlist(:,config1);        idof=unpaired(ii,config1)
           if (idof==0) then
              OFLWR "Unpaired error";CFLST
           endif
           firstspin=thisconfig(idof*2)
           thisconfig(idof*2)=mod(thisconfig(idof*2),2) + 1
           
           do jj=ii+1,numunpaired(config1)
              thatconfig=thisconfig;           jdof=unpaired(jj,config1)
              if (jdof==0) then
                 call openfile(); write(mpifileptr,*) "Unpaired error"
                 call closefile();              call mpistop()
              endif
              secondspin=thatconfig(jdof*2)

              if (secondspin.ne.firstspin) then
                 thatconfig(jdof*2)=mod(thatconfig(jdof*2),2) + 1
                 dirphase=reorder(thatconfig)
                 if (allowedconfig(thatconfig)) then
                    iwalk=iwalk+1
                    spinwalkdirphase(iwalk,config1)=dirphase
                    config2=getconfiguration(thatconfig)
                    spinwalk(iwalk,config1)=config2

  if ((config2.lt.allbotconfigs(iproc).or.config2.gt.alltopconfigs(iproc))) then
     OFLWR "BOT TOP NEWCONFIG ERR", config1,config2,iproc,allbotconfigs(iproc),alltopconfigs(iproc); CFLST
  endif
                 endif
              endif
           enddo   ! jj


        enddo  ! ii
        
        if (     numspinwalks(config1) /= iwalk ) then
           OFLWR "WALK ERROR SPIN.";CFLST
        endif

     enddo   ! config1
  enddo

  OFLWR "     ... done calculating spin walks.";  CFL

end subroutine spinwalks





subroutine spinsets_first()
  use dfconmod !! dfincludedmask
  use spinwalkmod
  use spinwalkinternal
  use configmod
  use parameters
  use aarrmod
  use mpimod
  implicit none

  integer ::  iwalk, jj,  iset, ilevel, currentnumwalks, prevnumwalks, flag, iflag, addwalks, &
       i, j, jwalk, jset, getdfindex, iproc
  integer, allocatable :: taken(:), tempwalks(:)

  OFLWR "   ... go spinsets ..."; CFL

  allocate(taken(numconfig), tempwalks(maxconfigsperproc))

  maxspinsetsize=0

  do jj=0,1

     taken(:)=0

     do iproc=1,nprocs

        iset=0;    jset=0

        do i=allbotconfigs(iproc),alltopconfigs(iproc)
           if (taken(i).ne.1) then
              if (dfincludedmask(i).ne.0) then
                 jset=jset+1
              endif
              taken(i)=1;           iset=iset+1
              ilevel=0;           flag=0
              tempwalks(1) = i;
              currentnumwalks = 1
              prevnumwalks = 0
              do while (flag==0)
                 ilevel=ilevel+1; flag=1;              addwalks=0
                 do j=prevnumwalks+1, currentnumwalks
                    do iwalk=1,numspinwalks(tempwalks(j))
                       iflag=0
                       do jwalk=1, currentnumwalks+addwalks
                          if (spinwalk(iwalk, tempwalks(j)) == tempwalks(jwalk)) then
                             iflag=1;                          exit
                          endif
                       enddo
                       if (iflag==0) then  ! walk is 
                          flag=0; addwalks=addwalks+1
                          tempwalks(currentnumwalks+addwalks)=spinwalk(iwalk,tempwalks(j))
                          if (taken(tempwalks(currentnumwalks+addwalks))==1) then
                             OFLWR "TAKEN ERROR!";CFLST
                          endif
                          taken(tempwalks(currentnumwalks+addwalks)) = 1
                       endif
                    enddo
                 enddo
                 prevnumwalks=currentnumwalks
                 currentnumwalks=currentnumwalks+addwalks
              enddo
              if (jj==0) then
                 spinsetsize(iset,iproc)=currentnumwalks
                 if (maxspinsetsize .lt. currentnumwalks) then
                    maxspinsetsize=currentnumwalks
                 endif
              else
                 if ((currentnumwalks.gt.maxspinsetsize).or.(spinsetsize(iset,iproc).ne.currentnumwalks)) then
                    OFLWR "WALK ERROR";CFLST
                 endif
                 spinsets(1:currentnumwalks,iset,iproc)=tempwalks(1:currentnumwalks)
                 if (dfincludedmask(i).ne.0) then
                    spindfsetindex(jset,iproc)=iset
                    do j=1,currentnumwalks
                       spindfsets(j,jset,iproc)=getdfindex(tempwalks(j))
                    enddo
                 endif
              endif
           endif
        enddo
        if (jj==0) then
           numspinsets(iproc)=iset
           numspindfsets(iproc)=jset
        else
           if (numspinsets(iproc).ne.iset) then
              OFLWR "NUMSPINSETS ERROR ", numspinsets(iproc), iset;CFLST
           endif
           if (numspindfsets(iproc).ne.jset) then
              OFLWR "NUMSPINdfSETS ERROR ", numspindfsets(iproc), jset;CFLST
           endif
        endif
     enddo  !! iproc

     do i=1,numconfig
        if (taken(i).ne.1) then
           OFLWR "TAKEN ERROR!!!!", i,taken(i);CFLST
        endif
     enddo

     j=0
     maxnumspinsets=0; maxnumspindfsets=0
     do iproc=1,nprocs
        if (maxnumspinsets.lt.numspinsets(iproc)) then
           maxnumspinsets=numspinsets(iproc)
        endif
        if (maxnumspindfsets.lt.numspindfsets(iproc)) then
           maxnumspindfsets=numspindfsets(iproc)
        endif
        do i=1,numspinsets(iproc)
           j=j+spinsetsize(i,iproc)
        enddo
     enddo
     if (j.ne.numconfig) then
        OFLWR "SPINSETSIZE ERROR!! ", j, numconfig; CFLST
     endif

     if (jj==0) then
        allocate(spinsets(maxspinsetsize,maxnumspinsets,nprocs),spindfsets(maxspinsetsize,maxnumspindfsets,nprocs),&
             spindfsetindex(maxnumspindfsets,nprocs))
     endif

  enddo

  deallocate(taken, tempwalks)

  OFLWR "Numspinsets this processor is ", numspinsets(myrank),"  maxspinset size is ", maxspinsetsize; CFL

end subroutine spinsets_first




subroutine getnumspinwalks()
  use spinwalkmod
  use spinwalkinternal
  use configmod
  use parameters
  implicit none

  integer ::   ispf,  config1, flag, idof, jdof,iwalk , thisconfig(ndof),  thatconfig(ndof), &
       ii, jj, dirphase, reorder, firstspin, secondspin
  real*8 :: avgspinwalks
  logical :: allowedconfig !! extraconfig


     OFLWR "Doing spin projector.";  CFL

     do config1=1,numconfig
        unpaired(:,config1)=0;    numunpaired(config1)=0;   msvalue(config1)=0
        thisconfig=configlist(:,config1)
        do idof=1,numelec
           msvalue(config1)=msvalue(config1) + thisconfig(idof*2)*2-3
           ispf=thisconfig(idof*2-1)
           flag=0
           do jdof=1,numelec
              if ((jdof.ne.idof).and.(thisconfig(jdof*2-1).eq.ispf)) then
                 flag=1
                 exit
              endif
           enddo
           if (flag==0) then
              numunpaired(config1)=numunpaired(config1)+1
              unpaired(numunpaired(config1),config1)=idof
           endif
        enddo
     enddo

     do config1=1,numconfig
        iwalk=0
        do ii=1,numunpaired(config1)
           thisconfig=configlist(:,config1)
           idof=unpaired(ii,config1)
           if (idof==0) then
              call openfile(); write(mpifileptr,*) "Unpaired error"; call closefile(); call mpistop()
           endif
           firstspin=thisconfig(idof*2)
           thisconfig(idof*2)=mod(thisconfig(idof*2),2) + 1
           do jj=ii+1,numunpaired(config1)
              thatconfig=thisconfig
              jdof=unpaired(jj,config1)
              if (jdof==0) then
                 call openfile(); write(mpifileptr,*) "Unpaired error"; call closefile(); call mpistop()
              endif
              secondspin=thatconfig(jdof*2)
              if (secondspin.ne.firstspin) then
                 thatconfig(jdof*2)=mod(thatconfig(jdof*2),2) + 1
                 dirphase=reorder(thatconfig)
                 if (allowedconfig(thatconfig)) then
                    iwalk=iwalk+1
                 endif   ! allowedconfig
              endif
           enddo   ! jj
        enddo  ! ii
        numspinwalks(config1) = iwalk 
     enddo   ! config1


  maxspinwalks=0
  avgspinwalks=0.d0
  do config1=1,numconfig
     avgspinwalks = avgspinwalks + numspinwalks(config1)
     if (maxspinwalks.lt.numspinwalks(config1)) then
        maxspinwalks=numspinwalks(config1)
     endif
  enddo

  avgspinwalks=avgspinwalks/numconfig

  OFLWR "Maximum number of spin walks= ",  maxspinwalks
  write(mpifileptr, *) "Avg number of spin walks= ",  avgspinwalks;CFL
  
end subroutine getnumspinwalks





subroutine configspin_matel()   
  use spinwalkmod
  use spinwalkinternal
  use parameters
  implicit none
  integer ::     config2, config1,   iwalk, myind

     configspinmatel(:,:)=0.d0
     do config1=1,numconfig
        myind=1
        
!! msvalue is 2x ms quantum number

        configspinmatel(myind,config1) = msvalue(config1)**2/4.d0 + numunpaired(config1)/2.d0
        
        do iwalk=1,numspinwalks(config1)
           config2=spinwalk(iwalk,config1)
           myind=1+iwalk
           configspinmatel(myind,config1) = configspinmatel(myind,config1) + &
                spinwalkdirphase(iwalk,config1)
        enddo
     enddo


end subroutine configspin_matel



function spinallowed(spinval)
  use parameters
  implicit none
  real*8 :: spinval
  logical :: spinallowed
  if (abs(spinval-(spinrestrictval/2.d0*(spinrestrictval/2.d0+1))).lt.1.d-3) then
     spinallowed=.true.
  else
     spinallowed=.false.
  endif
end function



subroutine configspinset_projector()   
  use spinwalkmod
  use spinwalkinternal
  use configmod   !! configlist for numspindfconfig
  use mpimod
  use parameters
  implicit none
  integer :: info, lwork,j,i,ii,iset,jj, elim, elimsets, flag, iwalk,&
       iproc
  real*8, allocatable :: spinvects(:,:), spinvals(:), work(:), realprojector(:,:)
  logical :: spinallowed,dfallowed
!  DATATYPE :: doublevects(maxspinsetsize**2)
!  real*8 :: doublevals(maxspinsetsize)
  
  OFLWR "Getting Spinset Projectors.  Numspinsets this process is ", numspinsets(myrank)
  WRFL "                                        maxspinsetsize is ", maxspinsetsize; CFL

  allocate(spinsetprojector(maxnumspinsets,nprocs))
  allocate(spinvects(maxspinsetsize,maxspinsetsize), spinvals(maxspinsetsize), &
       realprojector(maxspinsetsize,maxspinsetsize))

  lwork=10*maxspinsetsize;  allocate(work(lwork))
  
  allocate(spinsperproc(nprocs),spindfsperproc(nprocs), allbotspins(nprocs),alltopspins(nprocs),&
       allbotspindfs(nprocs), alltopspindfs(nprocs))

  spinsperproc(:)=0; spindfsperproc(:)=0

  elim=0;  elimsets=0;  

  do iproc=1,nprocs

     iset=1
  
     do while (iset.le.numspinsets(iproc))
        spinvects=0.d0
        do ii=1,spinsetsize(iset,iproc)
           do jj=1,spinsetsize(iset,iproc)
              
              if (ii.eq.jj) then
                 spinvects(ii,jj)=configspinmatel(1, spinsets(jj,iset,iproc))
              else
                 flag=0
                 do iwalk=1,numspinwalks(spinsets(jj,iset,iproc))
                    if (spinwalk(iwalk,spinsets(jj,iset,iproc)).eq.spinsets(ii,iset,iproc)) then
                       spinvects(ii,jj)=configspinmatel(iwalk+1, spinsets(jj,iset,iproc))
                       flag=1
                       exit
                    endif
                 enddo
              endif
           enddo
        enddo
     
        call dsyev('V','U', spinsetsize(iset,iproc), spinvects, maxspinsetsize, spinvals, work, lwork, info)
        if (info/=0) then
           OFLWR  "INFO SSYEV", info; CFLST
        endif
        j=0; 
        do i=1,spinsetsize(iset,iproc)
           if (spinallowed(spinvals(i))) then
              j=j+1;           spinvects(:,j)=spinvects(:,i)
           endif
        enddo
        spinsetrank(iset,iproc)=j

        spinsperproc(iproc)=spinsperproc(iproc)+j

        if (dfallowed(configlist(:,spinsets(1,iset,iproc)))) then
           spindfsperproc(iproc)=spindfsperproc(iproc)+j
        endif

        spinvects(:,j+1:maxspinsetsize)=0d0
        
        if (spinsetrank(iset,iproc)==0) then 
           elimsets=elimsets+1
           elim=elim+spinsetsize(iset,iproc)
           spinsetsize(iset:numspinsets(iproc)-1,iproc)=spinsetsize(iset+1:numspinsets(iproc),iproc)
           spinsets(:,iset:numspinsets(iproc)-1,iproc)=spinsets(:,iset+1:numspinsets(iproc),iproc)
           numspinsets(iproc)=numspinsets(iproc)-1
        else
           allocate(spinsetprojector(iset,iproc)%mat(spinsetsize(iset,iproc), spinsetsize(iset,iproc)))
           allocate(spinsetprojector(iset,iproc)%vects(spinsetsize(iset,iproc), spinsetrank(iset,iproc)))
           
           spinsetprojector(iset,iproc)%vects(:,:)=spinvects(1:spinsetsize(iset,iproc), 1:spinsetrank(iset,iproc))
           
           call dgemm('N', 'T', spinsetsize(iset,iproc), spinsetsize(iset,iproc),spinsetrank(iset,iproc),1.0d0, spinvects, maxspinsetsize, &
                spinvects, maxspinsetsize ,0.0d0,realprojector, maxspinsetsize)
           
           spinsetprojector(iset,iproc)%mat(:,:) = realprojector(1:spinsetsize(iset,iproc), 1:spinsetsize(iset,iproc))

!just checking right
!        call CONFIGHERM(spinsetprojector(iset,iproc)%mat,spinsetsize(iset,iproc),spinsetsize(iset,iproc), doublevects,doublevals)
!        do i=1,spinsetsize(iset,iproc)
!           if (abs(doublevals(i)*2-1.d0)-1.d0 .gt. 1.d-9) then
!              OFLWR "SPIN PROJECTOR ERROR", doublevals(i); CFLST
!           endif
!        enddo

           iset=iset+1

        endif
     enddo

  enddo

  deallocate(spinvects, spinvals, realprojector, work)
  
  OFLWR "...done.  Eliminated ", elimsets, " sets with total rank ", elim   ; CFL

  do iproc=1,nprocs
     i=0
     do ii=1,numspinsets(iproc)
        i=i+spinsetrank(ii,iproc)
     enddo
     if (i.ne.spinsperproc(iproc)) then
        OFLWR "CHECKMEERRR",iproc,i; CFLST
     endif
  enddo


!!!TEMP
!  print *, "TEMP OUTPUT",numspinsets(:)
!  do iproc=1,nprocs
!     do iset=1,numspinsets(iproc)
!        print *, iset, spinsetprojector(iset,iproc)%mat(:,:), spinsetprojector(iset,iproc)%vects(:,:)
!     enddo
!  enddo
!print *, "TEMPSTOP"
!call mpistop()




  maxspinsperproc=0
  ii=0
  do i=1,nprocs
     allbotspins(i)=ii+1
     ii=ii+spinsperproc(i)
     alltopspins(i)=ii
     if (spinsperproc(i).gt.maxspinsperproc) then
        maxspinsperproc=spinsperproc(i)
     endif
  enddo
  numspinconfig=ii
  botspin=allbotspins(myrank)
  topspin=alltopspins(myrank)

!! SPIN DF

  maxspindfsperproc=0
  ii=0
  do i=1,nprocs
     allbotspindfs(i)=ii+1
     ii=ii+spindfsperproc(i)
     alltopspindfs(i)=ii
     if (spindfsperproc(i).gt.maxspindfsperproc) then
        maxspindfsperproc=spindfsperproc(i)
     endif
  enddo
  numspindfconfig=ii
  botdfspin=allbotspindfs(myrank)
  topdfspin=alltopspindfs(myrank)

  if (parconsplit.eq.0) then
     firstspinconfig=1
     lastspinconfig=numspinconfig
     localnspin=numspinconfig
  else
     firstspinconfig=botspin
     lastspinconfig=topspin
     localnspin=spinsperproc(myrank)
  endif


  OFLWR
  WRFL  "This processor: "
  WRFL "      spin sets        ", numspinsets(myrank)
  WRFL "      spin rank            ", spinsperproc(myrank), " of ", configsperproc(myrank) 
  WRFL "      botspin,topspin:     ", botspin,topspin
  WRFL "      df spin rank         ", spindfsperproc(myrank), " of ", dfconfsperproc(myrank) 
  WRFL "      botdfspin,topdfspin: ", botdfspin,topdfspin
  WRFL "All processors:"
  WRFL "      spin rank, S^2 = ", (spinrestrictval/2.d0*(spinrestrictval/2.d0+1)), " is ", numspinconfig, " out of ", numconfig
  WRFL "      df spin rank         ", numspindfconfig, " of ", numdfconfigs
  WRFL; CFL

  
end subroutine configspinset_projector


