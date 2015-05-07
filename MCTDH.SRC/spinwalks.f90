

!! SUBROUTINE FOR "WALKS" (WHICH CONFIGURATIONS CONNECT TO WHICH) FOR SPARSE SPIN PROJECTOR.

#include "Definitions.INC"

module spinwalkinternal
  implicit none
  integer, allocatable :: unpaired(:,:)
end module

subroutine spinwalkalloc()
  use parameters
  use spinwalkinternal
  use spinwalkmod
  implicit none


  
  OFLWR "Go spinwalks. "; CFL
  
  allocate(unpaired(numelec,botwalk:topwalk), numunpaired(botwalk:topwalk), msvalue(botwalk:topwalk), numspinwalks(botwalk:topwalk))

  call getnumspinwalks()

  allocate(spinwalk(maxspinwalks,botwalk:topwalk),spinwalkdirphase(maxspinwalks,botwalk:topwalk))
  allocate(configspinmatel(maxspinwalks+1,botwalk:topwalk))
  allocate(spinsetsize(topwalk-botwalk+1),spinsetrank(topwalk-botwalk+1))

end subroutine spinwalkalloc


subroutine spinwalkdealloc()
  use parameters
  use spinwalkmod
  use spinwalkinternal
  implicit none
  integer :: i

  deallocate(unpaired,numunpaired,msvalue,numspinwalks,spinwalk, spinwalkdirphase, configspinmatel, &
       spinsets, spinsetsize, spinsetrank)

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
     read(751,iostat=myiostat) spinwalk(:,botwalk:topwalk),spinwalkdirphase(:,botwalk:topwalk) 
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

     do config1=botwalk,topwalk

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


                    if ((config2.lt.botwalk.or.config2.gt.topwalk)) then
                       OFLWR "BOT TOP NEWCONFIG ERR", config1,config2,botwalk,topwalk; CFLST
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
        write(751) spinwalk(:,botwalk:topwalk),spinwalkdirphase(:,botwalk:topwalk)  
        OFLWR "   ...ok, wrote spinwalks..."; CFL

        if (walksinturn) then
           call afterbarrier()
        endif
     endif

  endif  !! walksonfile

end subroutine spinwalks





subroutine spinsets_first()
  use spinwalkmod
  use spinwalkinternal
  use configmod
  use parameters
  use aarrmod
  implicit none

  integer ::  iwalk, jj,  iset, ilevel, currentnumwalks, prevnumwalks, flag, iflag, addwalks,  i, j, jwalk
  integer, allocatable :: taken(:), tempwalks(:)

  if (walksonfile.eq.0) then

     allocate(taken(botwalk:topwalk), tempwalks(1:topwalk-botwalk+1))
     maxspinsetsize=0

     do jj=0,1
        taken=0;        iset=0
        do i=botwalk,topwalk

           if (taken(i).ne.1) then
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
              endif
           endif
        enddo
        if (jj==0) then
           numspinsets=iset
        else
           if (numspinsets.ne.iset) then
              OFLWR "NUMSPINSETS ERROR ", numspinsets, iset;CFLST
           endif
        endif
        do i=botwalk,topwalk
           if (taken(i).ne.1) then
              OFLWR "TAKEN ERROR!!!!", i,taken(i);CFLST
           endif
        enddo
        j=0
        do i=1,numspinsets
           j=j+spinsetsize(i)
        enddo
        if (j.ne.topwalk-botwalk+1) then
           OFLWR "SPINSETSIZE ERROR!! ", j, topwalk-botwalk+1, numconfig;CFLST
        endif

        call mympiimax(maxspinsetsize)

        if (jj==0) then
           allocate(spinsets(maxspinsetsize,numspinsets))
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
     read(751,iostat=myiostat)  numunpaired(botwalk:topwalk), msvalue(botwalk:topwalk), numspinwalks(botwalk:topwalk), unpaired(:,botwalk:topwalk)
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

     do config1=botwalk,topwalk
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

     do config1=botwalk,topwalk
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
        write(751)  numunpaired(botwalk:topwalk), msvalue(botwalk:topwalk), numspinwalks(botwalk:topwalk), &
             unpaired(:,botwalk:topwalk)
        OFLWR "   ...ok, wrote spin info..."; CFL

        if (walksinturn) then
           call afterbarrier()
        endif
     endif
  endif  !! walksonfile

  maxspinwalks=0
  avgspinwalks=0.d0
  do config1=botwalk,topwalk
     avgspinwalks = avgspinwalks + numspinwalks(config1)
     if (maxspinwalks.lt.numspinwalks(config1)) then
        maxspinwalks=numspinwalks(config1)
     endif
  enddo

     avgspinwalks=avgspinwalks/(topwalk-botwalk+1)

  OFLWR "Maximum number of spin walks= ",  maxspinwalks
  write(mpifileptr, *) "Avg number of spin walks= ",  avgspinwalks;CFL
  
end subroutine getnumspinwalks

