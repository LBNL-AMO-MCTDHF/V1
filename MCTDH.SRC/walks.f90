


!! DETERMINES WHICH CONFIGURATIONS HAVE NONZERO MATRIX ELEMENTS WITH WHICH OTHERS, AND STORES INFORMATION
!!  ABOUT THE ORBITAL MATRIX ELEMENTS OF WHICH THEY ARE COMPRISED

#include "Definitions.INC"

function highspinorder(thisconfig,ndof,numelec)
  implicit none
  integer,intent(in) :: ndof,numelec
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
   


function lowspinorder(thisconfig,ndof,numelec)
  implicit none
  integer,intent(in) :: ndof,numelec
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
        

subroutine walkalloc(www)
  use fileptrmod
  use mpimod
  use walkmod
  implicit none
  type(walktype) :: www
  logical :: highspinorder,lowspinorder

!! training wheels

  if (www%topconfig-www%botconfig.gt.0) then
     if (.not.highspinorder(www%configlist(:,www%topconfig),www%ndof,www%numelec)) then
        OFLWR "NOT HIGHSPIN",www%topconfig
        call printconfig(www%configlist(:,www%topconfig),www)
        CFLST
     endif
     if (.not.lowspinorder(www%configlist(:,www%botconfig),www%ndof,www%numelec)) then
        OFLWR "NOT LOWSPIN",www%botconfig
        call printconfig(www%configlist(:,www%topconfig),www)
        CFLST
     endif
  endif

!! 06-2015 configpserproc also in newconfig.f90

  allocate( www%numsinglewalks(www%configstart:www%configend) , www%numdoublewalks(www%configstart:www%configend) )
  allocate( www%numsinglediagwalks(www%configstart:www%configend) , www%numdoublediagwalks(www%configstart:www%configend) )

  call getnumwalks(www)
  OFLWR "Allocating singlewalks"; CFL
  allocate( www%singlewalk(www%maxsinglewalks,www%configstart:www%configend))
  www%singlewalk=-1
  allocate(www%singlediag(www%numelec,www%configstart:www%configend) )
  www%singlediag=-1
  allocate( www%singlewalkdirphase(www%maxsinglewalks,www%configstart:www%configend) )
  www%singlewalkdirphase=0
  allocate( www%singlewalkopspf(1:2,www%maxsinglewalks,www%configstart:www%configend) )
  www%singlewalkopspf=-1
  OFLWR "Allocating doublewalks"; CFL
  allocate( www%doublewalkdirspf(1:4,www%maxdoublewalks,www%configstart:www%configend ) )
  www%doublewalkdirspf=-1
  allocate( www%doublewalkdirphase(www%maxdoublewalks,www%configstart:www%configend) )
  www%doublewalkdirphase=0
  allocate( www%doublewalk(www%maxdoublewalks,www%configstart:www%configend))
  www%doublewalk=-1
  allocate(www%doublediag(www%numelec*(www%numelec-1),www%configstart:www%configend))
  www%doublediag=-1
  OFLWR "     ..done walkalloc."; CFL
end subroutine walkalloc


subroutine walkdealloc(www)
  use walkmod
  implicit none
  type(walktype) :: www
  deallocate( www%numsinglewalks,www%numsinglediagwalks )
  deallocate( www%numdoublewalks,www%numdoublediagwalks )
  deallocate( www%singlewalk )
  deallocate( www%singlewalkdirphase )
  deallocate( www%singlewalkopspf )
  deallocate( www%doublewalkdirspf )
  deallocate( www%doublewalkdirphase )
  deallocate( www%doublewalk)
end subroutine walkdealloc


subroutine configlistwrite(www,inconfiglistfile)
  use walkmod
  use mpimod
  implicit none
  type(walktype),intent(in) :: www
  character :: inconfiglistfile*(*)
  if (myrank.eq.1) then
     open(1088,file=inconfiglistfile,status="unknown",form="unformatted")
     write(1088) www%numconfig,www%ndof
     write(1088) www%configlist(:,:)
     close(1088)
  endif

end subroutine configlistwrite

subroutine configlistheaderread(iunit,readnumconfig,readndof)
  implicit none
  integer :: iunit,readnumconfig,readndof

  read(iunit) readnumconfig,readndof

end subroutine configlistheaderread


subroutine configlistread(iunit,readnumconfig,readndof, readconfiglist)

  implicit none
  integer :: iunit,readnumconfig,readndof, readconfiglist(readndof,readnumconfig)
  
  read(iunit) readconfiglist(:,:)

end subroutine configlistread



subroutine walks(www)
  use fileptrmod
  use sparse_parameters
  use walkmod
  use ham_parameters !! offaxispulseflag
  use walks_parameters !! sortwalks
  use aarrmod
  implicit none
  type(walktype) :: www
  integer :: iindex, iiindex, jindex, jjindex,  ispin, jspin, iispin, jjspin, ispf, jspf, iispf, jjspf, config2, config1, dirphase, &
       iind, flag, idof, iidof, jdof, iwalk, reorder, getconfiguration,getmval,idiag
  logical :: allowedconfig0
  integer :: thisconfig(www%ndof), thatconfig(www%ndof), temporb(2), temporb2(2), isize, &
       listorder(www%maxdoublewalks+www%maxsinglewalks)

  !!  ***********   SINGLES  **********

  OFLWR "Calculating walks.  Singles..."; CFL
  
  do config1=www%botconfig,www%topconfig

     if (mod(config1,1000).eq.0) then
        OFLWR config1, " out of ", www%topconfig; CFL
     endif

     iwalk=0
     thisconfig=www%configlist(:,config1)

     do idof=1,www%numelec   !! position in thisconfig that we're walking 

        temporb=thisconfig((idof-1)*2+1 : idof*2)
        ispf=temporb(1)
        ispin=temporb(2)
        iindex=iind(temporb)

        do jindex=1,2*www%nspf   !! the walk

           temporb=aarr(jindex)
           jspf=temporb(1)
           jspin=temporb(2)

           if (ispin.ne.jspin) then
              cycle
           endif
           
           flag=0
           do jdof=1,www%numelec
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

           dirphase=reorder(thatconfig,www%numelec)

           if (.not.allowedconfig0(www,thatconfig,www%dfwalklevel)) then
              cycle
           endif

           if (offaxispulseflag.eq.0.and.getmval(www,thatconfig).ne.getmval(www,thisconfig)) then
              cycle
           endif

           iwalk=iwalk+1
           www%singlewalkopspf(1:2,iwalk,config1)=[ ispf,jspf ]   !! ket, bra   bra is walk
           www%singlewalkdirphase(iwalk,config1)=dirphase
           
           config2=getconfiguration(thatconfig,www)
           
           www%singlewalk(iwalk,config1)=config2

        enddo   ! the walk
     enddo  ! position we're walking

     if (     www%numsinglewalks(config1) /= iwalk ) then
        OFLWR "WALK ERROR SINGLES.";        CFLST
     endif

  enddo   ! config1


  OFLWR "Calculating walks.  Doubles...";  call closefile()

  !!   ***********  DOUBLES  ************

  do config1=www%botconfig,www%topconfig

     if (mod(config1,1000).eq.0) then
        OFLWR config1, " out of ", www%topconfig;        CFL
     endif

     iwalk=0

     thisconfig=www%configlist(:,config1)

     do idof=1,www%numelec         !! positions in thisconfig that we're walking 
        do iidof=idof+1,www%numelec   !! 

           temporb=thisconfig((idof-1)*2+1 : idof*2)
           ispf=temporb(1)
           ispin=temporb(2)
           iindex=iind(temporb)

           temporb=thisconfig((iidof-1)*2+1 : iidof*2)
           iispf=temporb(1)
           iispin=temporb(2)
           iiindex=iind(temporb)

           do jindex=1,2*www%nspf   !! the walk
              
              temporb=aarr(jindex)
              jspf=temporb(1) 
              jspin=temporb(2)
              
              if (.not.ispin.eq.jspin) then
                 cycle
              endif

!! no more exchange separately

              do jjindex=1,2*www%nspf
                 if (jjindex.eq.jindex) then
                    cycle
                 endif
                 
                 temporb2=aarr(jjindex)
                 jjspf=temporb2(1)
                 jjspin=temporb2(2)
                 
                 if (.not.iispin.eq.jjspin) then
                    cycle
                 endif

                 flag=0
                 do jdof=1,www%numelec
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

                 dirphase=reorder(thatconfig,www%numelec)

                 if (.not.allowedconfig0(www,thatconfig,www%dfwalklevel)) then
                    cycle
                 endif

                 if (offaxispulseflag.eq.0.and.getmval(www,thatconfig).ne.getmval(www,thisconfig)) then
                    cycle
                 endif

                 
                 iwalk = iwalk+1
            
!!                                                      ket2   bra2   ket1   bra1
                 www%doublewalkdirspf(1:4,iwalk,config1)=[ iispf, jjspf, ispf, jspf ]
                 www%doublewalkdirphase(iwalk,config1)=dirphase
                 
                 config2=getconfiguration(thatconfig,www)
                 www%doublewalk(iwalk,config1)=config2

              enddo   ! the walk
           enddo
           
        enddo  ! position we're walking
     enddo

     if (     www%numdoublewalks(config1) /= iwalk ) then
        OFLWR "WALK ERROR DOUBLES.",config1,www%numdoublewalks(config1),iwalk; CFLST
     endif

  enddo   ! config1

  call mpibarrier()

  if (sortwalks.ne.0) then

     OFLWR "Sorting walks..."; CFL
     do config1=www%botconfig,www%topconfig
        
        call getlistorder(www%singlewalk(:,config1),listorder(:),www%numsinglewalks(config1))
        call listreorder(www%singlewalkdirphase(:,config1),listorder(:),www%numsinglewalks(config1),1)
        call listreorder(www%singlewalkopspf(:,:,config1),listorder(:),www%numsinglewalks(config1),2)
        call listreorder(www%singlewalk(:,config1),listorder(:),www%numsinglewalks(config1),1)

        call getlistorder(www%doublewalk(:,config1),listorder(:),www%numdoublewalks(config1))
        call listreorder(www%doublewalkdirphase(:,config1),listorder(:),www%numdoublewalks(config1),1)
        call listreorder(www%doublewalkdirspf(:,:,config1),listorder(:),www%numdoublewalks(config1),4)
        call listreorder(www%doublewalk(:,config1),listorder(:),www%numdoublewalks(config1),1)
     enddo
     OFLWR "    .... done sorting walks."; CFL
  endif

  call mpibarrier()

  if (sparseconfigflag.eq.0.and.www%maxsinglewalks.ne.0) then
     isize=2*www%maxsinglewalks
     call mpiallgather_i(www%singlewalkopspf,   www%numconfig*isize,www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
     isize=www%maxsinglewalks
     call mpiallgather_i(www%singlewalkdirphase,www%numconfig*isize,www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
     call mpiallgather_i(www%singlewalk,        www%numconfig*isize,www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
  endif


  if (sparseconfigflag.eq.0.and.www%maxdoublewalks.ne.0) then
     isize=4*www%maxdoublewalks
     call mpiallgather_i(www%doublewalkdirspf,  www%numconfig*isize,www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
     isize=www%maxdoublewalks
     call mpiallgather_i(www%doublewalkdirphase,www%numconfig*isize,www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
     call mpiallgather_i(www%doublewalk,        www%numconfig*isize,www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
  endif




  do config1=www%configstart,www%configend
     idiag=0
     do iwalk=1,www%numsinglewalks(config1)
        if (www%singlewalk(iwalk,config1).eq.config1) then
           idiag=idiag+1
           www%singlediag(idiag,config1)=iwalk
        endif
     enddo
     www%numsinglediagwalks(config1)=idiag
     idiag=0
     do iwalk=1,www%numdoublewalks(config1)
        if (www%doublewalk(iwalk,config1).eq.config1) then
           idiag=idiag+1
           www%doublediag(idiag,config1)=iwalk
        endif
     enddo
     www%numdoublediagwalks(config1)=idiag
  enddo
     
  
end subroutine walks



subroutine getnumwalks(www)
  use fileptrmod
  use ham_parameters
  use sparse_parameters
  use walkmod
  use mpimod
  use aarrmod
  implicit none
  type(walktype) :: www
  integer :: iindex, iiindex, jindex, jjindex,  ispin, jspin, iispin, jjspin, ispf, iispf,  config1,  &
       dirphase, iind, flag, idof, iidof, jdof,iwalk , reorder, getmval
  logical :: allowedconfig0
  integer :: thisconfig(www%ndof), thatconfig(www%ndof), temporb(2), temporb2(2),totwalks
  character(len=3) :: iilab
  character(len=4) :: iilab0

  if (nprocs.gt.999) then
  print *, "redim getnumwalks";  call mpistop()
  endif

  write(iilab0,'(I4)') myrank+1000
  iilab(:)=iilab0(2:4)
  

  !!  ***********   SINGLES  **********

  call mpibarrier()

     OFLWR "Counting walks. Singles";  CFL

     do config1=www%botconfig,www%topconfig

        iwalk=0
        thisconfig=www%configlist(:,config1)

        do idof=1,www%numelec   !! position in thisconfig that we're walking 
           
           temporb=thisconfig((idof-1)*2+1 : idof*2)
           ispf=temporb(1)
           ispin=temporb(2)
           iindex=iind(temporb)
           
           do jindex=1,www%nspf * 2   !! the walk
              
              temporb=aarr(jindex)
              jspin=temporb(2)
              if (ispin.ne.jspin) then  
                 cycle
              endif

              flag=0
              do jdof=1,www%numelec
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

              dirphase=reorder(thatconfig,www%numelec)

              if (.not.allowedconfig0(www,thatconfig,www%dfwalklevel)) then
                 cycle
              endif

              if (offaxispulseflag.eq.0.and.getmval(www,thatconfig).ne.getmval(www,thisconfig)) then
                 cycle
              endif
              
              iwalk=iwalk+1

           enddo   ! the walk
        enddo  ! position we're walking
        
        www%numsinglewalks(config1) = iwalk 

     enddo   ! config1

     if (sparseconfigflag.eq.0) then
        call mpiallgather_i(www%numsinglewalks(:),www%numconfig,www%configsperproc(:),www%maxconfigsperproc)
     endif


     OFLWR "Counting walks. Doubles"; CFL
     
  !!   ***********  DOUBLES  ************

     do config1=www%botconfig,www%topconfig
        if (mod(config1,1000).eq.0) then
           OFLWR config1, " out of ", www%topconfig;  CFL
        endif
        
        iwalk=0
        thisconfig=www%configlist(:,config1)
        
        do idof=1,www%numelec            !! positions in thisconfig that we're walking 
           do iidof=idof+1,www%numelec   !! 
              
              temporb=thisconfig((idof-1)*2+1 : idof*2)
              ispf=temporb(1)
              ispin=temporb(2)
              iindex=iind(temporb)
              
              temporb=thisconfig((iidof-1)*2+1 : iidof*2)
              iispf=temporb(1)
              iispin=temporb(2)
              iiindex=iind(temporb)
              
              do jindex=1,2*www%nspf   !! the walk

                 temporb=aarr(jindex)
                 jspin=temporb(2)

                 if (.not.ispin.eq.jspin) then
                    cycle
                 endif

!! no more exchange separately
                 do jjindex=1,2*www%nspf   !! the walk
                    if (jjindex.eq.jindex) then
                       cycle
                    endif

                    temporb2=aarr(jjindex)
                    jjspin=temporb2(2)

                    if (.not.iispin.eq.jjspin) then
                       cycle
                    endif

                    flag=0
                    do jdof=1,www%numelec
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
                    dirphase=reorder(thatconfig,www%numelec)

                    if (allowedconfig0(www,thatconfig,www%dfwalklevel)) then
                       if (offaxispulseflag.ne.0.or.getmval(www,thatconfig).eq.getmval(www,thisconfig)) then
                          iwalk = iwalk+1
                       endif
                    endif
                    
                 enddo   ! the walk
              enddo
           enddo  ! position we're walking
        enddo
        
        www%numdoublewalks(config1)=iwalk

     enddo   ! config1

     if (sparseconfigflag.eq.0) then
        call mpiallgather_i(www%numdoublewalks(:),www%numconfig,www%configsperproc(:),www%maxconfigsperproc)
     endif


  www%maxsinglewalks=0;  www%maxdoublewalks=0

  totwalks=0
  do config1=www%configstart,www%configend

     if (www%maxsinglewalks.lt.www%numsinglewalks(config1)) then
        www%maxsinglewalks=www%numsinglewalks(config1)
     endif
     if (www%maxdoublewalks.lt.www%numdoublewalks(config1)) then
        www%maxdoublewalks=www%numdoublewalks(config1)
     endif

     totwalks=totwalks+www%numsinglewalks(config1)+www%numdoublewalks(config1)

  enddo

  if (sparseconfigflag.ne.0) then
     call mympiireduceone(totwalks)
     call mympiimax(www%maxsinglewalks);  call mympiimax(www%maxdoublewalks)
  endif

  OFLWR;  write(mpifileptr, *) "Maximum number of"
  write(mpifileptr, *) "           single walks= ",  www%maxsinglewalks
  write(mpifileptr, *) "           double walks= ",  www%maxdoublewalks;  
  WRFL "  TOTAL walks:", totwalks,"maxdoublewalks*numconfig",www%maxdoublewalks*www%numconfig
  WRFL; CFL

end subroutine getnumwalks





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




