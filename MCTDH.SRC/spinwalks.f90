

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
  implicit none

  OFLWR "Go spinwalks. "; CFL
  
  allocate(unpaired(numelec,configstart:configend), numunpaired(configstart:configend), msvalue(configstart:configend), numspinwalks(configstart:configend))

  call getnumspinwalks()

  allocate(spinwalk(maxspinwalks,configstart:configend),spinwalkdirphase(maxspinwalks,configstart:configend))
  allocate(configspinmatel(maxspinwalks+1,configstart:configend))

  allocate(spinsetsize(configend-configstart+1),spinsetrank(configend-configstart+1))

end subroutine spinwalkinit


subroutine spinwalkdealloc()
  use parameters
  use spinwalkmod
  implicit none
  integer :: i

  deallocate(spinsets, spinsetsize, spinsetrank)

  do i=1,numspinsets
     deallocate(spinsetprojector(i)%mat,spinsetprojector(i)%vects)
  enddo
  deallocate(spinsetprojector)

end subroutine spinwalkdealloc



subroutine spinwalks()
  use spinwalkmod
  use spinwalkinternal
  use configmod
  use parameters
  use aarrmod
  implicit none

  integer ::     config1, config2, dirphase,  idof, jdof,iwalk , thisconfig(ndof),  &
       thatconfig(ndof), reorder,  getconfiguration, ii, jj, firstspin, secondspin, &
       myiostat
  logical :: allowedconfig !! extraconfig

  spinwalk=0;     spinwalkdirphase=0

  if (walksonfile.ne.0) then

     if (walksinturn) then
        call beforebarrier()
     endif

     OFLWR "   ...reading spinwalks..."; CFL
     read(751,iostat=myiostat) spinwalk(:,configstart:configend),spinwalkdirphase(:,configstart:configend) 
     OFLWR "   ...done reading spinwalks..."; CFL

     if (walksinturn) then
        call afterbarrier()
     endif


     call mympiimax(myiostat)
     if (myiostat.ne.0) then
        OFLWR "Read error for savewalks.BIN!  Delete it to recompute walks. 662", myiostat; CFLST
     endif

  else   !WALKSONFILE

     OFLWR "Calculating spin walks.";  call closefile()

     do config1=configstart,configend

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


                    if ((config2.lt.configstart.or.config2.gt.configend)) then
                       OFLWR "BOT TOP NEWCONFIG ERR", config1,config2,configstart,configend; CFLST
                    endif
                 endif
              endif
           enddo   ! jj


        enddo  ! ii
        
        if (     numspinwalks(config1) /= iwalk ) then
           OFLWR "WALK ERROR SPIN.";CFLST
        endif

     enddo   ! config1

     if (walkwriteflag.ne.0) then
        if (walksinturn) then
           call beforebarrier()
        endif
        OFLWR "   ...writing spinwalks..."; CFL
        write(751) spinwalk(:,configstart:configend),spinwalkdirphase(:,configstart:configend)  
        OFLWR "   ...ok, wrote spinwalks..."; CFL

        if (walksinturn) then
           call afterbarrier()
        endif
     endif

  endif  !! walksonfile

end subroutine spinwalks





subroutine spinsets_first()
  use dfconmod !! dfincludedmask
  use spinwalkmod
  use spinwalkinternal
  use configmod
  use parameters
  use aarrmod
  implicit none

  integer ::  iwalk, jj,  iset, ilevel, currentnumwalks, prevnumwalks, flag, iflag, addwalks,  i, j, jwalk, jset, getdfindex
  integer, allocatable :: taken(:), tempwalks(:)

  if (walksonfile.eq.0) then

     allocate(taken(configstart:configend), tempwalks(1:configend-configstart+1))
     maxspinsetsize=0

     do jj=0,1
        taken=0;        iset=0;    jset=0
        do i=configstart,configend
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
                 spinsetsize(iset)=currentnumwalks
                 if (maxspinsetsize .lt. currentnumwalks) then
                    maxspinsetsize=currentnumwalks
                 endif
              else
                 if ((currentnumwalks.gt.maxspinsetsize).or.(spinsetsize(iset).ne.currentnumwalks)) then
                    OFLWR "WALK ERROR";CFLST
                 endif
                 spinsets(1:currentnumwalks,iset)=tempwalks(1:currentnumwalks)
                 if (dfincludedmask(i).ne.0) then
                    spindfsetindex(jset)=iset
                    do j=1,currentnumwalks
                       spindfsets(j,jset)=getdfindex(tempwalks(j))
                    enddo
                 endif
              endif
           endif
        enddo
        if (jj==0) then
           numspinsets=iset
           numspindfsets=jset
        else
           if (numspinsets.ne.iset) then
              OFLWR "NUMSPINSETS ERROR ", numspinsets, iset;CFLST
           endif
           if (numspindfsets.ne.jset) then
              OFLWR "NUMSPINdfSETS ERROR ", numspindfsets, jset;CFLST
           endif
        endif
        do i=configstart,configend
           if (taken(i).ne.1) then
              OFLWR "TAKEN ERROR!!!!", i,taken(i);CFLST
           endif
        enddo
        j=0
        do i=1,numspinsets
           j=j+spinsetsize(i)
        enddo
        if (j.ne.configend-configstart+1) then
           OFLWR "SPINSETSIZE ERROR!! ", j, configend-configstart+1, numconfig;CFLST
        endif

        call mympiimax(maxspinsetsize)

        if (jj==0) then
           allocate(spinsets(maxspinsetsize,numspinsets),spindfsets(maxspinsetsize,numspindfsets),spindfsetindex(numspindfsets))
        endif
     enddo
     deallocate(taken, tempwalks)

     OFLWR "Numspinsets is ", numspinsets,"  maxspinset size is ", maxspinsetsize; CFL

  endif  !! walksonfile

end subroutine spinsets_first




subroutine getnumspinwalks()
  use spinwalkmod
  use spinwalkinternal
  use configmod
  use parameters
  implicit none

  integer ::   ispf,  config1, flag, idof, jdof,iwalk , thisconfig(ndof),  thatconfig(ndof), &
       ii, jj, dirphase, reorder, firstspin, secondspin, myiostat
  real*8 :: avgspinwalks
  logical :: allowedconfig !! extraconfig




  if (walksonfile.ne.0) then
     numunpaired=0; msvalue=0; numspinwalks=0; unpaired=0;

     if (walksinturn) then
        call beforebarrier()
     endif

     OFLWR "   ...reading spin projector....";  CFL
     read(751,iostat=myiostat)  numunpaired(configstart:configend), msvalue(configstart:configend), numspinwalks(configstart:configend), unpaired(:,configstart:configend)
     OFLWR "   ...done reading spin projector....";  CFL

     if (walksinturn) then
        call afterbarrier()
     endif

     call mympiimax(myiostat)
     if (myiostat.ne.0) then
        OFLWR "Read error for savewalks.BIN!  Delete it to recompute walks. 887", myiostat; CFLST
     endif
  else

     OFLWR "Doing spin projector.";  CFL

     do config1=configstart,configend
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

     do config1=configstart,configend
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

     if (walkwriteflag.ne.0) then
        if (walksinturn) then
           call beforebarrier()
        endif

        OFLWR "   ...writing spin info..."; CFL
        write(751)  numunpaired(configstart:configend), msvalue(configstart:configend), numspinwalks(configstart:configend), &
             unpaired(:,configstart:configend)
        OFLWR "   ...ok, wrote spin info..."; CFL

        if (walksinturn) then
           call afterbarrier()
        endif
     endif
  endif  !! walksonfile

  maxspinwalks=0
  avgspinwalks=0.d0
  do config1=configstart,configend
     avgspinwalks = avgspinwalks + numspinwalks(config1)
     if (maxspinwalks.lt.numspinwalks(config1)) then
        maxspinwalks=numspinwalks(config1)
     endif
  enddo

     avgspinwalks=avgspinwalks/(configend-configstart+1)

  OFLWR "Maximum number of spin walks= ",  maxspinwalks
  write(mpifileptr, *) "Avg number of spin walks= ",  avgspinwalks;CFL
  
end subroutine getnumspinwalks





subroutine configspin_matel()   
  use spinwalkmod
  use spinwalkinternal
  use parameters
  implicit none
  integer ::     config2, config1,   iwalk, myind,myiostat

  if (walksonfile.ne.0) then

     if (walksinturn) then 
        call beforebarrier()
     endif
     
     OFLWR "   ...reading configspinmatel..."; CFL
     read(751,iostat=myiostat) configspinmatel
     OFLWR "   ...ok, done reading configspinmatel..."; CFL

     if (walksinturn) then 
        call afterbarrier()
     endif

     call mympiimax(myiostat)
     if (myiostat.ne.0) then
        OFLWR "Read error for savewalks.BIN!  Delete it to recompute walks. 662", myiostat;CFLST
     endif

  else   !WALKSONFILE

     configspinmatel(:,:)=0.d0
     do config1=configstart,configend
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
     if (walkwriteflag.ne.0) then
        if (walksinturn) then
           call beforebarrier()
        endif

        OFLWR "   ...writing configspinmatel..."; CFL
        write(751) configspinmatel
        OFLWR "   ...ok, done writing configspinmatel..."; CFL

        if (walksinturn) then
           call afterbarrier()
        endif
     endif
  endif

end subroutine configspin_matel



subroutine configspinset_projector()   
  use spinwalkmod
  use spinwalkinternal
  use configmod   !! configlist for numspindfconfig
  use mpimod
  use parameters
  implicit none
  integer :: info, lwork,j,i,ii,iset,jj, elim, elimsets, flag, iwalk,&
       spindfrank,spinrank,spinlocalrank,spinlocaldfrank
  real*8, allocatable :: spinvects(:,:), spinvals(:), work(:), realprojector(:,:)
  logical :: spinallowed,dfallowed
  integer :: allbottemp(nprocs)
!  DATATYPE :: doublevects(maxspinsetsize**2)
!  real*8 :: doublevals(maxspinsetsize)
  
  OFLWR "Getting Spinset Projectors.  Numspinsets is ", numspinsets, " maxspinsetsize is ", maxspinsetsize; CFL


     allocate(spinsetprojector(numspinsets))

     allocate(spinvects(maxspinsetsize,maxspinsetsize), spinvals(maxspinsetsize), &
          realprojector(maxspinsetsize,maxspinsetsize))
     lwork=10*maxspinsetsize;  allocate(work(lwork))
     
     elim=0;  elimsets=0;  iset=1; spinrank=0; spindfrank=0
     spinlocalrank=0; spinlocaldfrank=0
  
     do while (iset.le.numspinsets)
        spinvects=0.d0
        do ii=1,spinsetsize(iset)
           do jj=1,spinsetsize(iset)
              
              if (ii.eq.jj) then
                 spinvects(ii,jj)=configspinmatel(1, spinsets(jj,iset))
              else
                 flag=0
                 do iwalk=1,numspinwalks(spinsets(jj,iset))
                    if (spinwalk(iwalk,spinsets(jj,iset)).eq.spinsets(ii,iset)) then
                       spinvects(ii,jj)=configspinmatel(iwalk+1, spinsets(jj,iset))
                       flag=1
                       exit
                    endif
                 enddo
              endif
           enddo
        enddo
     
        call dsyev('V','U', spinsetsize(iset), spinvects, maxspinsetsize, spinvals, work, lwork, info)
        if (info/=0) then
           OFLWR  "INFO SSYEV", info; CFLST
        endif
        j=0; 
        do i=1,spinsetsize(iset)
           if (spinallowed(spinvals(i))) then
              j=j+1;           spinvects(:,j)=spinvects(:,i)
           endif
        enddo
        spinsetrank(iset)=j

        spinrank=spinrank+j
        if (spinsets(1,iset).ge.botconfig.and.spinsets(1,iset).le.topconfig) then
           spinlocalrank=spinlocalrank+j
        endif
        if (dfallowed(configlist(:,spinsets(1,iset)))) then
           spindfrank=spindfrank+j
           if (spinsets(1,iset).ge.botconfig.and.spinsets(1,iset).le.topconfig) then
              spinlocaldfrank=spinlocaldfrank+j
           endif
        endif

        spinvects(:,j+1:maxspinsetsize)=0d0
        
        if (spinsetrank(iset)==0) then 
           elimsets=elimsets+1
           elim=elim+spinsetsize(iset)
           spinsetsize(iset:numspinsets-1)=spinsetsize(iset+1:numspinsets)
           spinsets(:,iset:numspinsets-1)=spinsets(:,iset+1:numspinsets)
           numspinsets=numspinsets-1
        else
           allocate(spinsetprojector(iset)%mat(spinsetsize(iset), spinsetsize(iset)))
           allocate(spinsetprojector(iset)%vects(spinsetsize(iset), spinsetrank(iset)))
           
           spinsetprojector(iset)%vects(:,:)=spinvects(1:spinsetsize(iset), 1:spinsetrank(iset))
           
           call dgemm('N', 'T', spinsetsize(iset), spinsetsize(iset),spinsetrank(iset),1.0d0, spinvects, maxspinsetsize, &
                spinvects, maxspinsetsize ,0.0d0,realprojector, maxspinsetsize)
           
           spinsetprojector(iset)%mat(:,:) = realprojector(1:spinsetsize(iset), 1:spinsetsize(iset))

!just checking right
!        call CONFIGHERM(spinsetprojector(iset)%mat,spinsetsize(iset),spinsetsize(iset), doublevects,doublevals)
!        do i=1,spinsetsize(iset)
!           if (abs(doublevals(i)*2-1.d0)-1.d0 .gt. 1.d-9) then
!              OFLWR "SPIN PROJECTOR ERROR", doublevals(i); CFLST
!           endif
!        enddo

           iset=iset+1

        endif
     enddo

     deallocate(spinvects, spinvals, realprojector, work)
  
     OFLWR "...done.  Eliminated ", elimsets, " sets with total rank ", elim   ; CFL

  i=0
  do ii=1,numspinsets
     i=i+spinsetrank(ii)
  enddo
  if (i.ne.spinrank) then
     OFLWR "CHECKMEERRR",i,spinrank; CFLST
  endif

  
  OFLWR  "Number of spin sets is now ", numspinsets
  WRFL "Number of spinvects with S^2 = ", (spinrestrictval/2.d0*(spinrestrictval/2.d0+1)), " is ", spinrank, " out of ", configend-configstart+1; CFL




  allocate(spinsperproc(nprocs)); spinsperproc=(-1)

  spinsperproc(myrank)=spinlocalrank
  maxspinsperproc=0
  ii=0
  do i=1,nprocs
     allbottemp(i)=ii+1
     call mympiibcastone(spinsperproc(i),i)
     ii=ii+spinsperproc(i)
     if (spinsperproc(i).gt.maxspinsperproc) then
        maxspinsperproc=spinsperproc(i)
     endif
  enddo
  numspinconfig=ii
  botspin=allbottemp(myrank)
  topspin=botspin+spinlocalrank-1


  if (sparseconfigflag.eq.0) then
     if (spinrank.ne.numspinconfig) then
        OFLWR "ACK CHECKME SPINNN"; CFLST
     endif
     spinstart=1
     spinend=numspinconfig
  else
     if (spinrank.ne.spinlocalrank) then
        OFLWR "ACK CHECKME SPINNNxxxx"; CFLST
     endif
     spinstart=botspin
     spinend=topspin
  endif

!! SPIN DF

  allocate(spindfsperproc(nprocs)); spindfsperproc=(-1)

  spindfsperproc(myrank)=spinlocaldfrank
  maxspindfsperproc=0
  ii=0
  do i=1,nprocs
     allbottemp(i)=ii+1
     call mympiibcastone(spindfsperproc(i),i)
     ii=ii+spindfsperproc(i)
     if (spindfsperproc(i).gt.maxspindfsperproc) then
        maxspindfsperproc=spindfsperproc(i)
     endif
  enddo
  numspindfconfig=ii
  botdfspin=allbottemp(myrank)
  topdfspin=botdfspin+spinlocaldfrank-1

  if (sparseconfigflag.eq.0) then
     if (numspindfconfig.ne.spindfrank) then
        OFLWR "CHECKME FAIL SPINDF ",numspindfconfig,spindfrank; CFLST
     endif
  endif

  if (sparseconfigflag.eq.0) then
     spindfstart=1
     spindfend=numspindfconfig
  else
     spindfstart=botdfspin
     spindfend=topdfspin
  endif

  if (parconsplit.eq.0) then
     firstspinconfig=1
     lastspinconfig=numspinconfig
     localnspin=numspinconfig
  else
     firstspinconfig=spinstart
     lastspinconfig=spinend
     localnspin=spinrank
  endif


  OFLWR "TOTAL (all processors) spin rank",numspinconfig," out of ", numconfig
  WRFL "     Spin DF rank (all processors) ", numspindfconfig
  WRFL "   This proc:  spinstart,spinend : ", spinstart,spinend
  WRFL; CFL

  
end subroutine configspinset_projector


