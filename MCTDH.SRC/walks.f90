

!! DETERMINES WHICH CONFIGURATIONS HAVE NONZERO MATRIX ELEMENTS WITH WHICH OTHERS, AND STORES INFORMATION
!!  ABOUT THE ORBITAL MATRIX ELEMENTS OF WHICH THEY ARE COMPRISED

#include "Definitions.INC"

function highspinorder(thisconfig)
  use parameters
  implicit none
  logical :: highspinorder
  integer :: thisconfig(ndof),ii,unpaired(numelec),flag,jj

  highspinorder=.true.

  unpaired(1:numelec)=1

  do ii=1,numelec
     do jj=1,numelec   !! WORKS
        if (jj.ne.ii) then   !!WORKS
! -xAVX error on lawrencium!  doesnt work this way.  compiler/instruction set bug.
!     do jj=ii+1,numelec   !!FAILS
           if (thisconfig(jj*2-1).eq.thisconfig(ii*2-1)) then
              unpaired(ii)=0
              unpaired(jj)=0
           endif
        endif     !!WORKS
     enddo
  enddo
  
  flag=0
  do ii=1,numelec
     if (unpaired(ii).eq.1) then
        if (thisconfig(ii*2).eq.1) then
           flag=1
        else
           if (flag==1) then
              highspinorder=.false.
              return
           endif
        endif
     endif
  enddo


end function highspinorder
   


function lowspinorder(thisconfig)
  use parameters
  implicit none
  logical :: lowspinorder
  integer :: thisconfig(ndof),ii,unpaired(numelec),flag,jj

  lowspinorder=.true.

  unpaired(:)=1

  do ii=1,numelec
!     do jj=ii+1,numelec    !!FAILS
     do jj=1,numelec        !!WORKS
        if (jj.ne.ii) then  !!WORKS
        if (thisconfig(jj*2-1).eq.thisconfig(ii*2-1)) then
           unpaired(ii)=0
           unpaired(jj)=0
        endif
        endif               !!WORKS
     enddo
  enddo

  flag=0
  do ii=1,numelec
     if (unpaired(ii).eq.1) then
        if (thisconfig(ii*2).eq.2) then
           flag=1
        else
           if (flag==1) then
              lowspinorder=.false.
           endif
        endif
     endif
  enddo


end function lowspinorder
        

subroutine walkalloc()
  use parameters
  use mpimod
  use configmod !! configlist for newconfigflag
  use walkmod
  implicit none
  integer :: ii,i
  logical :: highspinorder,lowspinorder

!! botwalk and topwalk in newconfig.f90.
!! now, with fast_newconfiglist, just check

  if (topwalk-botwalk.gt.0) then
     if (.not.highspinorder(configlist(:,topwalk))) then
        OFLWR "NOT HIGHSPIN",topwalk
        call printconfig(configlist(:,topwalk))
        CFLST
     endif
     if (.not.lowspinorder(configlist(:,botwalk))) then
        OFLWR "NOT LOWSPIN",botwalk
        call printconfig(configlist(:,topwalk))
        CFLST
     endif
  endif
     
  allocate(configsperproc(nprocs))

  if (sparseconfigflag.eq.0) then
     maxconfigsperproc=(-1); configsperproc(:)=(-1)       !! SHOULDN'T BE USED
  else
     configsperproc(myrank)=topwalk-botwalk+1
     ii=0
     maxconfigsperproc=0
     do i=1,nprocs
        call mympiibcast(configsperproc(i),i,1)
        ii=ii+configsperproc(i)
        if (configsperproc(i).gt.maxconfigsperproc) then
           maxconfigsperproc=configsperproc(i)
        endif
     enddo
     if (ii.ne.numconfig) then
        OFLWR "WTF EERRROROR", ii, numconfig; CFLST
     endif
     do i=1,nprocs
        if (configsperproc(i).lt.0) then
           OFLWR "Configs per proc lt 0 for proc ",i,"nprocs must be greater than numconfig???  Can't do."; CFLST
        endif
     enddo
  endif

  allocate( numsinglewalks(botwalk:topwalk) , numdoublewalks(botwalk:topwalk) )
  allocate( numsinglediagwalks(botwalk:topwalk) , numdoublediagwalks(botwalk:topwalk) )

  call getnumwalks()
  OFLWR "Allocating singlewalks"; CFL
  allocate( singlewalk(maxsinglewalks,botwalk:topwalk), singlediag(numelec,botwalk:topwalk) )
  singlewalk=-1
  allocate( singlewalkdirphase(maxsinglewalks,botwalk:topwalk) )
  singlewalkdirphase=0
  allocate( singlewalkopspf(1:2,maxsinglewalks,botwalk:topwalk) )
  singlewalkopspf=-1
  OFLWR "Allocating doublewalks"; CFL
  allocate( doublewalkdirspf(1:4,maxdoublewalks,botwalk:topwalk ) )
  doublewalkdirspf=-1
  allocate( doublewalkdirphase(maxdoublewalks,botwalk:topwalk) )
  doublewalkdirphase=0
  allocate( doublewalk(maxdoublewalks,botwalk:topwalk), doublediag(numelec*(numelec-1),botwalk:topwalk) )
  doublewalk=-1
  OFLWR "     ..done walkalloc."; CFL
end subroutine walkalloc


subroutine walkdealloc()
  use parameters
  use walkmod
  implicit none
  deallocate( numsinglewalks,numsinglediagwalks )
  deallocate( numdoublewalks,numdoublediagwalks )
  deallocate( singlewalk )
  deallocate( singlewalkdirphase )
  deallocate( singlewalkopspf )
  deallocate( doublewalkdirspf )
  deallocate( doublewalkdirphase )
  deallocate( doublewalk)
end subroutine walkdealloc


subroutine singlewalkwrite()
  use parameters
  use configmod
  use walkmod
  use mpimod
  implicit none

!! beforebarrier and afterbarrier in main

  if (myrank.eq.1) then
     open(1088,file=configlistfile,status="unknown",form="unformatted")
     write(1088) numconfig,ndof,maxsinglewalks
     write(1088) configlist(:,:)
     write(1088) numsinglewalks(:)
     write(1088) singlewalk(:,:)
     write(1088) singlewalkopspf(:,:,:)
     write(1088) singlewalkdirphase(:,:)
     close(1088)
  endif

end subroutine singlewalkwrite

subroutine singlewalkheaderread(iunit,readnumconfig,readndof,readmaxsinglewalks)
  implicit none
  integer :: iunit,readnumconfig,readndof,readmaxsinglewalks

  read(iunit) readnumconfig,readndof,readmaxsinglewalks

end subroutine singlewalkheaderread


subroutine singlewalkread(iunit,readnumconfig,readndof,readmaxsinglewalks,  readconfiglist,readnumsinglewalks, &
     readsinglewalks, readsinglewalkopspf, readsinglewalkdirphase)
  implicit none
  integer :: iunit,readnumconfig,readndof,readmaxsinglewalks, &
       readconfiglist(readndof,readnumconfig), readnumsinglewalks(readnumconfig), &
       readsinglewalks(readmaxsinglewalks,readnumconfig), &
       readsinglewalkopspf(2,readmaxsinglewalks,readnumconfig), &
       readsinglewalkdirphase(readmaxsinglewalks,readnumconfig)
  
  read(iunit) readconfiglist(:,:)
  read(iunit) readnumsinglewalks(:)
  read(iunit) readsinglewalks(:,:)
  read(iunit) readsinglewalkopspf(:,:,:)
  read(iunit) readsinglewalkdirphase(:,:)

end subroutine singlewalkread



subroutine walks()
  use walkmod
  use configmod
  use parameters
  use aarrmod
  implicit none

  integer :: iindex, iiindex, jindex, jjindex,  ispin, jspin, iispin, jjspin, ispf, jspf, iispf, jjspf, config2, config1, dirphase, &
       iind, flag, idof, iidof, jdof, iwalk, reorder, getconfiguration,myiostat,getmval,idiag
  logical :: allowedconfig   !! extraconfig
  integer :: thisconfig(ndof), thatconfig(ndof), temporb(2), temporb2(2), &
       newdoublediag(numelec*(numelec-1)), newsinglediag(numelec), listorder(maxdoublewalks+maxsinglewalks), i

  !!  ***********   SINGLES  **********

  if (walksonfile.ne.0) then

     if (walksinturn) then 
        call beforebarrier()
     endif

     OFLWR "Reading walks. Singles";  CFL

     read(751,iostat=myiostat) &
          singlewalkopspf, &
          singlewalkdirphase, &
          singlewalk, &
          numsinglediagwalks, &
          singlediag, &

          doublewalkdirspf, &
          doublewalkdirphase, &
          doublewalk, &
          numdoublediagwalks, &
          doublediag

     OFLWR "   ..read single walks this processor...";  CFL

     if (walksinturn) then 
        call afterbarrier()
     endif

     call mympiimax(myiostat)
     if (myiostat.ne.0) then
        OFLWR "Read error for savewalks.BIN!  Delete it to recompute walks!!", myiostat; CFLST
     else
        OFLWR "savewalks.BIN was found and read. SKIPPING WALK CALCULATION.";CFL
        return
!! RETURN
        return
!! RETURN
        return
!! RETURN
        return
!! RETURN
        return
     endif
  else
     OFLWR "No savewalks.BIN found.  Calculating walks.  Singles...";  CFL
  endif

  
  do config1=botwalk,topwalk
     if (mod(config1,1000).eq.0) then
        OFLWR config1, " out of ", topwalk;        call closefile()
     endif

     iwalk=0; idiag=0
     thisconfig=configlist(:,config1)

     do idof=1,numelec   !! position in thisconfig that we're walking 

        temporb=thisconfig((idof-1)*2+1 : idof*2)
        ispf=temporb(1)
        ispin=temporb(2)
        iindex=iind(temporb)
        
        do jindex=1,spftot   !! the walk

           temporb=aarr(jindex,nspf)
           jspf=temporb(1)
           jspin=temporb(2)

           if (ispin.ne.jspin) then
              cycle
           endif
           
           flag=0
           do jdof=1,numelec
              if (jdof.ne.idof) then !! INCLUDING DIAGONAL WALKS
                 if (iind(thisconfig((jdof-1)*2+1:jdof*2)) == jindex) then 
                    flag=1
                 endif
              endif
           enddo

           if (flag.ne.0) then    ! pauli dis allowed configuration.
              cycle
           endif

           thatconfig=thisconfig
           thatconfig((idof-1)*2+1  : idof*2)=temporb

           dirphase=reorder(thatconfig)
           
           if (.not.allowedconfig(thatconfig)) then
              cycle
           endif

           if (offaxispulseflag.eq.0.and.getmval(thatconfig).ne.getmval(thisconfig)) then
              cycle
           endif

           iwalk=iwalk+1
           singlewalkopspf(1:2,iwalk,config1)=[ ispf,jspf ]   !! ket, bra   bra is walk
           singlewalkdirphase(iwalk,config1)=dirphase
           
           config2=getconfiguration(thatconfig)
           
           singlewalk(iwalk,config1)=config2

           if (config2.eq.config1) then
              idiag=idiag+1
              singlediag(idiag,config1)=iwalk
           endif

        enddo   ! the walk
     enddo  ! position we're walking

     if (     numsinglewalks(config1) /= iwalk ) then
        OFLWR "WALK ERROR SINGLES.";        CFLST
     endif

     numsinglediagwalks(config1)=idiag

  enddo   ! config1


  OFLWR "Calculating walks.  Doubles...";  call closefile()

  !!   ***********  DOUBLES  ************

  do config1=botwalk,topwalk

     if (mod(config1,1000).eq.0) then
        OFLWR config1, " out of ", topwalk;        CFL
     endif

     iwalk=0; idiag=0
     thisconfig=configlist(:,config1)

     do idof=1,numelec         !! positions in thisconfig that we're walking 
        do iidof=idof+1,numelec   !! 

           temporb=thisconfig((idof-1)*2+1 : idof*2)
           ispf=temporb(1)
           ispin=temporb(2)
           iindex=iind(temporb)

           temporb=thisconfig((iidof-1)*2+1 : iidof*2)
           iispf=temporb(1)
           iispin=temporb(2)
           iiindex=iind(temporb)

           do jindex=1,spftot   !! the walk
              
              temporb=aarr(jindex,nspf)
              jspf=temporb(1) 
              jspin=temporb(2)
              
              if (.not.ispin.eq.jspin) then
                 cycle
              endif

!! no more exchange separately

              do jjindex=1,spftot
                 if (jjindex.eq.jindex) then
                    cycle
                 endif
                 
                 temporb2=aarr(jjindex,nspf)
                 jjspf=temporb2(1)
                 jjspin=temporb2(2)
                 
                 if (.not.iispin.eq.jjspin) then
                    cycle
                 endif

                 flag=0
                 do jdof=1,numelec
                    if (jdof.ne.idof.and.jdof.ne.iidof) then !! INCLUDING DIAGONAL AND SINGLE WALKS
                       if ((iind(thisconfig((jdof-1)*2+1:jdof*2)) == jindex).or. &
                            (iind(thisconfig((jdof-1)*2+1:jdof*2)) == jjindex)) then
                          flag=1
                          exit
                       endif
                    endif
                 enddo
                 
                 if (flag.ne.0) then    ! pauli dis allowed configuration.
                    cycle
                 endif

                 
                 thatconfig=thisconfig
                 thatconfig((idof-1)*2+1  : idof*2)=temporb
                 thatconfig((iidof-1)*2+1  : iidof*2)=temporb2

                 dirphase=reorder(thatconfig)

                 if (.not.allowedconfig(thatconfig)) then
                    cycle
                 endif

                 if (offaxispulseflag.eq.0.and.getmval(thatconfig).ne.getmval(thisconfig)) then
                    cycle
                 endif

                 
                 iwalk = iwalk+1
            
!!                                                      ket2   bra2   ket1   bra1
                 doublewalkdirspf(1:4,iwalk,config1)=[ iispf, jjspf, ispf, jspf ]
                 doublewalkdirphase(iwalk,config1)=dirphase
                 
                 config2=getconfiguration(thatconfig)
                 doublewalk(iwalk,config1)=config2

                 if (config2.eq.config1) then
                    idiag=idiag+1
                    doublediag(idiag,config1)=iwalk
                 endif

                 
              enddo   ! the walk
           enddo
           
        enddo  ! position we're walking
     enddo

     if (     numdoublewalks(config1) /= iwalk ) then
        OFLWR "WALK ERROR DOUBLES.",config1,numdoublewalks(config1),iwalk; CFLST
     endif

     numdoublediagwalks(config1)=idiag

  enddo   ! config1


!! SINGLEDIAG, DOUBLEDIAG NOT DEBUGGED

  if (sortwalks.ne.0) then

     OFLWR "Sorting walks..."; CFL
     do config1=botwalk,topwalk
        if (mod(config1,1000).eq.0) then
           OFLWR "   ...config ", config1," of ", topwalk; CFL
        endif
        
        call getlistorder(singlewalk(:,config1),listorder(:),numsinglewalks(config1))
        call listreorder(singlewalkdirphase(:,config1),listorder(:),numsinglewalks(config1),1)
        call listreorder(singlewalkopspf(:,:,config1),listorder(:),numsinglewalks(config1),2)
        call listreorder(singlewalk(:,config1),listorder(:),numsinglewalks(config1),1)
!!?? CORRECT ?? NOT DEBUGGED
        do i=1,numsinglediagwalks(config1)
           newsinglediag(i)=listorder(singlediag(i,config1))
        enddo
        singlediag(:,config1)=newsinglediag(:)

        call getlistorder(doublewalk(:,config1),listorder(:),numdoublewalks(config1))
        call listreorder(doublewalkdirphase(:,config1),listorder(:),numdoublewalks(config1),1)
        call listreorder(doublewalkdirspf(:,:,config1),listorder(:),numdoublewalks(config1),4)
        call listreorder(doublewalk(:,config1),listorder(:),numdoublewalks(config1),1)
!!?? CORRECT ?? NOT DEBUGGED
        do i=1,numdoublediagwalks(config1)
           newdoublediag(i)=listorder(doublediag(i,config1))
        enddo
        doublediag(:,config1)=newdoublediag(:)
     enddo
     OFLWR "    .... done sorting walks."; CFL
  endif
  
  if (walkwriteflag.ne.0) then
        if (walksinturn) then
           call beforebarrier()
        endif
        OFLWR "    ... writing double walks ..."; CFL
        write(751) &
             singlewalkopspf, &
             singlewalkdirphase, &
             singlewalk, &
             numsinglediagwalks, &
             singlediag, &
             
             doublewalkdirspf, &
             doublewalkdirphase, &
             doublewalk, &
             numdoublediagwalks, &
             doublediag
        OFLWR "    ... done writing double walks this processor..."; CFL
        if (walksinturn) then 
           call afterbarrier()
        endif
        OFLWR "    ... done writing double walks."; CFL
  endif

end subroutine walks



subroutine getnumwalks()
  use walkmod
  use configmod
  use parameters
  use mpimod
  use aarrmod
  implicit none

  integer :: iindex, iiindex, jindex, jjindex,  ispin, jspin, iispin, jjspin, ispf, iispf,  config1,innumconfig,  &
       dirphase, iind, flag, idof, iidof, jdof,iwalk , reorder, myiostat, inprocs , getmval

  logical :: allowedconfig !! extraconfig
  integer :: thisconfig(ndof), thatconfig(ndof), temporb(2), temporb2(2),totwalks
  character(len=3) :: iilab
  character(len=4) :: iilab0

  if (nprocs.gt.999) then
  print *, "redim getnumwalks";  call mpistop()
  endif

  write(iilab0,'(I4)') myrank+1000
  iilab(:)=iilab0(2:4)
  

  !!  ***********   SINGLES  **********

  call mpibarrier()

  flag=1
  
  open(751,file="WALKS/savewalks.BIN",status="old",iostat=myiostat, form="unformatted")
  if (myiostat==0) then
     if (myrank.ne.1) then
        close(751)
        open(751,file="WALKS/savewalks.BIN"//iilab,status="old",iostat=myiostat, form="unformatted")
        myiostat=myiostat*myrank
     endif
     call mympiimax(myiostat)
     if (myiostat.ne.0) then
        OFLWR "error savewalks ", myiostat; CFLST
     endif

     if (walksinturn) then 
        call beforebarrier()
     endif

     OFLWR "  ...reading walks...."; CFL
     read(751,iostat=myiostat) inprocs, innumconfig
     OFLWR "  ...ok reading walks...."; CFL

     if (walksinturn) then 
        call afterbarrier()
     endif

     call mympiimax(myiostat)

     if (myiostat==0) then
        if (inprocs.ne.nprocs) then
           OFLWR "savewalks.BIN exists but is for nprocs=", inprocs," not current nprocs=",nprocs      ;     CFLST
        endif
        if (innumconfig.ne.numconfig) then
           OFLWR "savewalks.BIN exists but has wrong number of configurations.  Current:", numconfig, " on file:", innumconfig;           CFLST
        endif
        OFLWR "savewalks.BIN exists; reading walks, not computing.";        CFL
        walksonfile=1;        flag=0
     endif

     if (walksinturn) then 
        call beforebarrier()
     endif

     read(751,iostat=myiostat) numsinglewalks,numdoublewalks

     if (walksinturn) then 
        call afterbarrier()
     endif

     call mympiimax(myiostat)
     if (myiostat.ne.0) then
        OFLWR "error savewalks xx", myiostat; CFLST
     endif
  endif

  if (flag==1) then

     OFLWR "Counting walks. Singles";  CFL

     do config1=botwalk,topwalk
        
        iwalk=0
        thisconfig=configlist(:,config1)
        
        do idof=1,numelec   !! position in thisconfig that we're walking 
           
           temporb=thisconfig((idof-1)*2+1 : idof*2)
           ispf=temporb(1)
           ispin=temporb(2)
           iindex=iind(temporb)
           
           do jindex=1,spftot   !! the walk
              
              temporb=aarr(jindex,nspf)
              jspin=temporb(2)
              if (ispin.ne.jspin) then  
                 cycle
              endif

              flag=0
              do jdof=1,numelec
                 if (jdof.ne.idof) then !! INCLUDING DIAGONAL WALKS
                    if (iind(thisconfig((jdof-1)*2+1:jdof*2)) == jindex) then 
                       flag=1
                    endif
                 endif
              enddo
              
              if (flag.ne.0) then    ! pauli dis allowed configuration.
                 cycle
              endif

              thatconfig=thisconfig
              thatconfig((idof-1)*2+1  : idof*2)=temporb

              dirphase=reorder(thatconfig)

              if (.not.allowedconfig(thatconfig)) then
                 cycle
              endif

              if (offaxispulseflag.eq.0.and.getmval(thatconfig).ne.getmval(thisconfig)) then
                 cycle
              endif
              
              iwalk=iwalk+1

           enddo   ! the walk
        enddo  ! position we're walking
        
        numsinglewalks(config1) = iwalk 
        
     enddo   ! config1

     OFLWR "Counting walks. Doubles"; CFL
     
  !!   ***********  DOUBLES  ************

     do config1=botwalk,topwalk
        if (mod(config1,1000).eq.0) then
           OFLWR config1, " out of ", topwalk;        call closefile()
        endif
        
        iwalk=0
        thisconfig=configlist(:,config1)
        
        do idof=1,numelec         !! positions in thisconfig that we're walking 
           do iidof=idof+1,numelec   !! 
              
              temporb=thisconfig((idof-1)*2+1 : idof*2)
              ispf=temporb(1)
              ispin=temporb(2)
              iindex=iind(temporb)
              
              temporb=thisconfig((iidof-1)*2+1 : iidof*2)
              iispf=temporb(1)
              iispin=temporb(2)
              iiindex=iind(temporb)
              
              do jindex=1,spftot   !! the walk

                 temporb=aarr(jindex,nspf)
                 jspin=temporb(2)

                 if (.not.ispin.eq.jspin) then
                    cycle
                 endif

!! no more exchange separately
                 do jjindex=1,spftot   !! the walk
                    if (jjindex.eq.jindex) then
                       cycle
                    endif

                    temporb2=aarr(jjindex,nspf)
                    jjspin=temporb2(2)

                    if (.not.iispin.eq.jjspin) then
                       cycle
                    endif

                    flag=0
                    do jdof=1,numelec
                       if (jdof.ne.idof.and.jdof.ne.iidof) then  !! INCLUDING DIAGONAL AND SINGLE WALKS
                          if ((iind(thisconfig((jdof-1)*2+1:jdof*2)) == jindex).or. &
                               (iind(thisconfig((jdof-1)*2+1:jdof*2)) == jjindex)) then
                             flag=1
                          endif
                       endif
                    enddo
                    
                    if (flag.ne.0) then    ! pauli dis allowed configuration.
                       cycle
                    endif
                    
                    thatconfig=thisconfig
                    thatconfig((idof-1)*2+1  : idof*2)=temporb
                    thatconfig((iidof-1)*2+1  : iidof*2)=temporb2
                    dirphase=reorder(thatconfig)
                    
                    if (allowedconfig(thatconfig)) then
                       if (offaxispulseflag.ne.0.or.getmval(thatconfig).eq.getmval(thisconfig)) then
                          iwalk = iwalk+1
                       endif
                    endif
                    
                 enddo   ! the walk
              enddo
           enddo  ! position we're walking
        enddo
        
        numdoublewalks(config1)=iwalk
        
     enddo   ! config1

     if (walkwriteflag.ne.0) then
        if (walksinturn) then
           OFLWR "OPENING WALKS and writing in turn...."; CFL
           call beforebarrier()
        else
           OFLWR "OPENING WALKS...."; CFL
        endif
        if (myrank.eq.1) then
           open(751,file="WALKS/walks.BIN",status="unknown", form="unformatted")
        else
           open(751,file="WALKS/walks.BIN"//iilab,status="unknown", form="unformatted")
        endif

        write(751) nprocs, numconfig;     write(751) numsinglewalks,numdoublewalks

        if (walksinturn) then
           call afterbarrier()
        endif
        OFLWR "    ...opened walks and wrote header info"; CFL
     endif
  endif

  maxsinglewalks=0;  maxdoublewalks=0

  totwalks=0
  do config1=botwalk,topwalk

     if (maxsinglewalks.lt.numsinglewalks(config1)) then
        maxsinglewalks=numsinglewalks(config1)
     endif
     if (maxdoublewalks.lt.numdoublewalks(config1)) then
        maxdoublewalks=numdoublewalks(config1)
     endif

     totwalks=totwalks+numsinglewalks(config1)+numdoublewalks(config1)

  enddo

  call mympiireduceone(totwalks)
  call mympiimax(maxsinglewalks);  call mympiimax(maxdoublewalks)


  OFLWR;  write(mpifileptr, *) "Maximum number of"
  write(mpifileptr, *) "           single walks= ",  maxsinglewalks
  write(mpifileptr, *) "           double walks= ",  maxdoublewalks;  
  WRFL "  TOTAL walks:", totwalks,"maxdoublewalks*numconfig",maxdoublewalks*numconfig
  WRFL; CFL

end subroutine getnumwalks




subroutine dfrestrictmatrix(inmatrix)
  use parameters
  implicit none
  DATATYPE :: inmatrix(numconfig,numconfig),tempmatrix(numconfig,numconfig)

  if (dfrestrictflag.eq.0) then
     return
  endif
  tempmatrix=TRANSPOSE(inmatrix)
  call dfrestrict(tempmatrix,numconfig)
  inmatrix=TRANSPOSE(tempmatrix)
  call dfrestrict(inmatrix,numconfig)
  
end subroutine dfrestrictmatrix
  




subroutine getlistorder(values, order,num)
  use fileptrmod
  implicit none
  integer :: num, values(num),taken(num), order(num)
  integer :: i,j,whichlowest, flag, lowval

  taken=0;  order=-1
  do j=1,num
     whichlowest=-1; flag=0;     lowval=10000000  !! is not used (see flag)
     do i=1,num
        if ( taken(i) .eq. 0 ) then
           if ((flag.eq.0) .or.(values(i) .le. lowval)) then
              flag=1;              lowval=values(i); whichlowest=i
           endif
        endif
     enddo
     if ((whichlowest.gt.num).or.(whichlowest.lt.1)) then
         OFLWR taken,"lowest ERROR, J=",j," WHICHLOWEST=", whichlowest;   CFLST
     endif
     if (taken(whichlowest).ne.0) then
        OFLWR "TAKENmm ERROR.";        CFLST
     endif
     taken(whichlowest)=1;     order(j)=whichlowest
  enddo

end subroutine getlistorder


subroutine listreorder(list, order,num,numper)
  implicit none
  integer :: num, numper, list(numper,num),order(num),newvals(numper,num),j

  do j=1,num
     newvals(:,j)=list(:,order(j))
  enddo

  list(:,:)=newvals(:,:)

end subroutine listreorder




