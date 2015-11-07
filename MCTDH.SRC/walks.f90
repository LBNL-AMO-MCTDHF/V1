


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
  allocate( www%firstsinglewalkbyproc(nprocs,www%configstart:www%configend), www%lastsinglewalkbyproc(nprocs,www%configstart:www%configend) )
  allocate( www%firstdoublewalkbyproc(nprocs,www%configstart:www%configend), www%lastdoublewalkbyproc(nprocs,www%configstart:www%configend) )


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
  use mpimod !! nprocs
  use ham_parameters !! offaxispulseflag
  use aarrmod
  implicit none
  type(walktype) :: www
  integer :: iindex, iiindex, jindex, jjindex,  ispin, jspin, iispin, jjspin, ispf, jspf, iispf, jjspf, config2, config1, dirphase, &
       iind, flag, idof, iidof, jdof, iwalk, reorder, getconfiguration,getmval,idiag
  logical :: allowedconfig0
  integer :: thisconfig(www%ndof), thatconfig(www%ndof), temporb(2), temporb2(2), isize, iproc, &
       listorder(www%maxdoublewalks+www%maxsinglewalks)

  !!  ***********   SINGLES  **********

  OFLWR "Calculating walks.  Singles..."; CFL
  
  do config1=www%botconfig,www%topconfig

     if (mod(config1,1000).eq.0) then
        OFLWR config1, " out of ", www%topconfig; CFL
     endif

     iwalk=0

     if (www%singlewalkflag.ne.0) then

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
     endif  ! singlewalkflag

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

     if (www%doublewalkflag.ne.0) then

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
     endif  ! doublewalkflag

     if (     www%numdoublewalks(config1) /= iwalk ) then
        OFLWR "WALK ERROR DOUBLES.",config1,www%numdoublewalks(config1),iwalk; CFLST
     endif

  enddo   ! config1

  call mpibarrier()

  OFLWR "Sorting walks..."; CFL
  do config1=www%botconfig,www%topconfig

     if (www%numsinglewalks(config1).gt.1) then
        call getlistorder(www%singlewalk(:,config1),listorder(:),www%numsinglewalks(config1))
        call listreorder(www%singlewalkdirphase(:,config1),listorder(:),www%numsinglewalks(config1),1)
        call listreorder(www%singlewalkopspf(:,:,config1),listorder(:),www%numsinglewalks(config1),2)
        call listreorder(www%singlewalk(:,config1),listorder(:),www%numsinglewalks(config1),1)
     endif
     if (www%numdoublewalks(config1).gt.1) then
        call getlistorder(www%doublewalk(:,config1),listorder(:),www%numdoublewalks(config1))
        call listreorder(www%doublewalkdirphase(:,config1),listorder(:),www%numdoublewalks(config1),1)
        call listreorder(www%doublewalkdirspf(:,:,config1),listorder(:),www%numdoublewalks(config1),4)
        call listreorder(www%doublewalk(:,config1),listorder(:),www%numdoublewalks(config1),1)
     endif
  enddo
  OFLWR "    .... done sorting walks."; CFL


  call mpibarrier()
  
  do config1=www%botconfig,www%topconfig

     www%firstsinglewalkbyproc(1,config1)=1
     iproc=1
     do iwalk=1,www%numsinglewalks(config1)
        do while (www%singlewalk(iwalk,config1).gt.www%alltopconfigs(iproc))
           www%lastsinglewalkbyproc(iproc,config1)=iwalk-1
           iproc=iproc+1
           www%firstsinglewalkbyproc(iproc,config1)=iwalk
        enddo
     enddo
     www%lastsinglewalkbyproc(iproc,config1)=www%numsinglewalks(config1)
     www%firstsinglewalkbyproc(iproc+1:nprocs,config1)=www%numsinglewalks(config1)+1
     www%lastsinglewalkbyproc(iproc+1:nprocs,config1)=www%numsinglewalks(config1)

     www%firstdoublewalkbyproc(1,config1)=1
     iproc=1
     do iwalk=1,www%numdoublewalks(config1)
        do while (www%doublewalk(iwalk,config1).gt.www%alltopconfigs(iproc))
           www%lastdoublewalkbyproc(iproc,config1)=iwalk-1
           iproc=iproc+1
           www%firstdoublewalkbyproc(iproc,config1)=iwalk
        enddo
     enddo
     www%lastdoublewalkbyproc(iproc,config1)=www%numdoublewalks(config1)
     www%firstdoublewalkbyproc(iproc+1:nprocs,config1)=www%numdoublewalks(config1)+1
     www%lastdoublewalkbyproc(iproc+1:nprocs,config1)=www%numdoublewalks(config1)

  end do

  call mpibarrier()


  if (sparseconfigflag.eq.0) then
     isize=nprocs
     call mpiallgather_i(www%firstsinglewalkbyproc(:,:),   www%numconfig*isize,www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
     call mpiallgather_i(www%lastsinglewalkbyproc(:,:),   www%numconfig*isize,www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
     call mpiallgather_i(www%firstdoublewalkbyproc(:,:),   www%numconfig*isize,www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
     call mpiallgather_i(www%lastdoublewalkbyproc(:,:),   www%numconfig*isize,www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
  endif

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

  call mpibarrier()

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

        if (mod(config1,1000).eq.0) then
           OFLWR config1, " out of ", www%topconfig;  CFL
        endif

        iwalk=0

        if (www%singlewalkflag.ne.0) then

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
        endif  ! singlewalkflag

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

        if (www%doublewalkflag.ne.0) then

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
        endif   ! doublewalkflag

        www%numdoublewalks(config1)=iwalk

     enddo   ! config1

     if (sparseconfigflag.eq.0) then
        call mpiallgather_i(www%numdoublewalks(:),www%numconfig,www%configsperproc(:),www%maxconfigsperproc)
     endif


!!$  www%maxsinglewalks=0;  www%maxdoublewalks=0
  www%maxsinglewalks=1;  www%maxdoublewalks=1      !! ensure always allocate

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


subroutine hops(www)
  use fileptrmod
  use ham_parameters
  use sparse_parameters
  use walkmod
  use mpimod
  use aarrmod
  implicit none
  type(walktype) :: www
  integer :: ii,iwalk,iconfig,totsinglehops,totdoublehops,totsinglewalks,totdoublewalks,ihop,flag,iproc,isize

  integer :: numsinglehopsbyproc(nprocs), numdoublehopsbyproc(nprocs)


  allocate(www%numsinglehops(www%configstart:www%configend),&
       www%numdoublehops(www%configstart:www%configend))
  allocate( www%singlediaghop(www%configstart:www%configend),&
       www%doublediaghop(www%configstart:www%configend))
  allocate( www%singlehopdiagflag(www%configstart:www%configend),&
       www%doublehopdiagflag(www%configstart:www%configend))

  www%singlediaghop(:)=(-999)
  www%doublediaghop(:)=(-999)

  allocate( www%firstsinglehopbyproc(nprocs,www%configstart:www%configend), www%lastsinglehopbyproc(nprocs,www%configstart:www%configend) )
  allocate( www%firstdoublehopbyproc(nprocs,www%configstart:www%configend), www%lastdoublehopbyproc(nprocs,www%configstart:www%configend) )

  do ii=0,1

     if (ii.eq.0) then
!! avoid warn bounds
        allocate(www%singlehop(1,1),www%singlehopwalkstart(1,1),www%singlehopwalkend(1,1),&
             www%doublehop(1,1),www%doublehopwalkstart(1,1),www%doublehopwalkend(1,1))
     else
        deallocate(www%singlehop,www%singlehopwalkstart,www%singlehopwalkend,&
             www%doublehop,www%doublehopwalkstart,www%doublehopwalkend)
        allocate(www%singlehop(www%maxnumsinglehops,www%configstart:www%configend),&
             www%singlehopwalkstart(www%maxnumsinglehops,www%configstart:www%configend),&
             www%singlehopwalkend(www%maxnumsinglehops,www%configstart:www%configend),&
             www%doublehop(www%maxnumdoublehops,www%configstart:www%configend),&
             www%doublehopwalkstart(www%maxnumdoublehops,www%configstart:www%configend),&
             www%doublehopwalkend(www%maxnumdoublehops,www%configstart:www%configend))
     endif

     if (ii.eq.0) then
        OFLWR "Counting single hops..."; CFL
     else
        OFLWR "Getting single hops..."; CFL
     endif

     do iconfig=www%configstart,www%configend
        ihop=0
        if (www%numsinglewalks(iconfig).gt.0) then
           ihop=1
           if (ii.eq.1) then
              www%singlehop(1,iconfig)=www%singlewalk(1,iconfig)
              www%singlehopwalkstart(1,iconfig)=1
           endif
           do iwalk=2,www%numsinglewalks(iconfig)
              if (www%singlewalk(iwalk,iconfig).ne.www%singlewalk(iwalk-1,iconfig)) then
                 if (ii.eq.1) then
                    www%singlehopwalkend(ihop,iconfig)=iwalk-1
                 endif
                 ihop=ihop+1
                 if (ii.eq.1) then
                    www%singlehop(ihop,iconfig)=www%singlewalk(iwalk,iconfig)
                    www%singlehopwalkstart(ihop,iconfig)=iwalk
                 endif
              endif
           enddo
        endif ! if numsinglewalks.gt.0
        if (ii.eq.0) then
           www%numsinglehops(iconfig)=ihop
        else
           if (www%numsinglewalks(iconfig).gt.0) then
              www%singlehopwalkend(ihop,iconfig)=www%numsinglewalks(iconfig)
           endif
           if (www%numsinglehops(iconfig).ne.ihop) then
              OFLWR "CHECKME SINGLEHOPW",www%numsinglehops(iconfig),ihop,iconfig; CFLST
           endif
        endif
     enddo


     if (ii.eq.0) then
        OFLWR "Counting double hops..."; CFL
     else
        OFLWR "Getting double hops..."; CFL
     endif

     do iconfig=www%configstart,www%configend
        ihop=0
        if (www%numdoublewalks(iconfig).gt.0) then
           ihop=1
           if (ii.eq.1) then
              www%doublehop(1,iconfig)=www%doublewalk(1,iconfig)
              www%doublehopwalkstart(1,iconfig)=1
           endif
           do iwalk=2,www%numdoublewalks(iconfig)
              if (www%doublewalk(iwalk,iconfig).ne.www%doublewalk(iwalk-1,iconfig)) then
                 if (ii.eq.1) then
                    www%doublehopwalkend(ihop,iconfig)=iwalk-1
                 endif
                 ihop=ihop+1
                 if (ii.eq.1) then
                    www%doublehop(ihop,iconfig)=www%doublewalk(iwalk,iconfig)
                    www%doublehopwalkstart(ihop,iconfig)=iwalk
                 endif
              endif
           enddo
        endif  !! if numdoublewalks.gt.0
        if (ii.eq.0) then
           www%numdoublehops(iconfig)=ihop
        else
           if (www%numdoublewalks(iconfig).gt.0) then
              www%doublehopwalkend(ihop,iconfig)=www%numdoublewalks(iconfig)
           endif
           if (www%numdoublehops(iconfig).ne.ihop) then
              OFLWR "CHECKME DOUBLEHOPW",www%numdoublehops(iconfig),ihop,iconfig; CFLST
           endif
        endif
     enddo

     if (ii.eq.0) then
!!$        www%maxnumsinglehops=0
!!$        www%maxnumdoublehops=0
        www%maxnumsinglehops=1   !always allocate
        www%maxnumdoublehops=1
        totsinglehops=0; totsinglewalks=0
        totdoublehops=0; totdoublewalks=0
        do iconfig=www%configstart,www%configend
           totsinglehops=totsinglehops+www%numsinglehops(iconfig)
           totsinglewalks=totsinglewalks+www%numsinglewalks(iconfig)
           if (www%numsinglehops(iconfig).gt.www%maxnumsinglehops) then
              www%maxnumsinglehops=www%numsinglehops(iconfig)
           endif
           totdoublehops=totdoublehops+www%numdoublehops(iconfig)
           totdoublewalks=totdoublewalks+www%numdoublewalks(iconfig)
           if (www%numdoublehops(iconfig).gt.www%maxnumdoublehops) then
              www%maxnumdoublehops=www%numdoublehops(iconfig)
           endif
        enddo
     endif

  enddo

  do iconfig=www%configstart,www%configend

     flag=0
     if (www%numsinglehops(iconfig).gt.0) then
        do ihop=1,www%numsinglehops(iconfig)
           if (www%singlehop(ihop,iconfig).eq.iconfig) then
              if (flag.eq.1) then
                 OFLWR "EERRR HOPSSING"; CFLST
              else
                 flag=1
                 www%singlediaghop(iconfig)=ihop
              endif
           endif
        enddo
     endif
     www%singlehopdiagflag(iconfig)=flag

     flag=0
     if (www%numdoublehops(iconfig).gt.0) then
        do ihop=1,www%numdoublehops(iconfig)
           if (www%doublehop(ihop,iconfig).eq.iconfig) then
              if (flag.eq.1) then
                 OFLWR "EERRR HOPSDOUB"; CFLST
              else
                 flag=1
                 www%doublediaghop(iconfig)=ihop
              endif
           endif
        enddo
     endif
     www%doublehopdiagflag(iconfig)=flag
  enddo

  call mpibarrier()
  
  do iconfig=www%botconfig,www%topconfig

     www%firstsinglehopbyproc(1,iconfig)=1
     iproc=1
     do ihop=1,www%numsinglehops(iconfig)
        do while (www%singlehop(ihop,iconfig).gt.www%alltopconfigs(iproc))
           www%lastsinglehopbyproc(iproc,iconfig)=ihop-1
           iproc=iproc+1
           www%firstsinglehopbyproc(iproc,iconfig)=ihop
        enddo
     enddo
     www%lastsinglehopbyproc(iproc,iconfig)=www%numsinglehops(iconfig)
     www%firstsinglehopbyproc(iproc+1:nprocs,iconfig)=www%numsinglehops(iconfig)+1
     www%lastsinglehopbyproc(iproc+1:nprocs,iconfig)=www%numsinglehops(iconfig)

     www%firstdoublehopbyproc(1,iconfig)=1
     iproc=1
     do ihop=1,www%numdoublehops(iconfig)
        do while (www%doublehop(ihop,iconfig).gt.www%alltopconfigs(iproc))
           www%lastdoublehopbyproc(iproc,iconfig)=ihop-1
           iproc=iproc+1
           www%firstdoublehopbyproc(iproc,iconfig)=ihop
        enddo
     enddo
     www%lastdoublehopbyproc(iproc,iconfig)=www%numdoublehops(iconfig)
     www%firstdoublehopbyproc(iproc+1:nprocs,iconfig)=www%numdoublehops(iconfig)+1
     www%lastdoublehopbyproc(iproc+1:nprocs,iconfig)=www%numdoublehops(iconfig)

  end do

  call mpibarrier()

  if (sparseconfigflag.eq.0) then
     isize=nprocs
     call mpiallgather_i(www%firstsinglehopbyproc(:,:),   www%numconfig*isize,www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
     call mpiallgather_i(www%lastsinglehopbyproc(:,:),   www%numconfig*isize,www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
     call mpiallgather_i(www%firstdoublehopbyproc(:,:),   www%numconfig*isize,www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
     call mpiallgather_i(www%lastdoublehopbyproc(:,:),   www%numconfig*isize,www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
  endif


  if (sparseconfigflag.eq.0) then
     ii=1
     call mpiallgather_i(www%numdoublehops(:),www%numconfig*ii,www%configsperproc(:)*ii,www%maxconfigsperproc*ii)
     call mpiallgather_i(www%numsinglehops(:),www%numconfig*ii,www%configsperproc(:)*ii,www%maxconfigsperproc*ii)
     call mpiallgather_i(www%singlediaghop(:),www%numconfig*ii,www%configsperproc(:)*ii,www%maxconfigsperproc*ii)
     call mpiallgather_i(www%doublediaghop(:),www%numconfig*ii,www%configsperproc(:)*ii,www%maxconfigsperproc*ii)
     call mpiallgather_i(www%singlehopdiagflag(:),www%numconfig*ii,www%configsperproc(:)*ii,www%maxconfigsperproc*ii)
     call mpiallgather_i(www%doublehopdiagflag(:),www%numconfig*ii,www%configsperproc(:)*ii,www%maxconfigsperproc*ii)

     ii=www%maxnumsinglehops
     call mpiallgather_i(www%singlehop(:,:),www%numconfig*ii,www%configsperproc(:)*ii,www%maxconfigsperproc*ii)
     call mpiallgather_i(www%singlehopwalkstart(:,:),www%numconfig*ii,www%configsperproc(:)*ii,www%maxconfigsperproc*ii)
     call mpiallgather_i(www%singlehopwalkend(:,:),www%numconfig*ii,www%configsperproc(:)*ii,www%maxconfigsperproc*ii)
     ii=www%maxnumdoublehops
     call mpiallgather_i(www%doublehop(:,:),www%numconfig*ii,www%configsperproc(:)*ii,www%maxconfigsperproc*ii)
     call mpiallgather_i(www%doublehopwalkstart(:,:),www%numconfig*ii,www%configsperproc(:)*ii,www%maxconfigsperproc*ii)
     call mpiallgather_i(www%doublehopwalkend(:,:),www%numconfig*ii,www%configsperproc(:)*ii,www%maxconfigsperproc*ii)
  endif


  numsinglehopsbyproc(:)=0;   numdoublehopsbyproc(:)=0

  do iconfig=www%botconfig,www%topconfig
     numsinglehopsbyproc(:)=numsinglehopsbyproc(:) + &
          (www%lastsinglehopbyproc(:,iconfig)-www%firstsinglehopbyproc(:,iconfig)+1)
     numdoublehopsbyproc(:)=numdoublehopsbyproc(:) + &
          (www%lastdoublehopbyproc(:,iconfig)-www%firstdoublehopbyproc(:,iconfig)+1)
  enddo

  call mpibarrier()
  if (myrank.eq.1) then
     print *, "HOPS BY PROC ON PROCESSOR 1 :::::::::::::::::::::::"
     print *, "   singles:"
     write(*,'(I5,A2,1000I7)') myrank,": ",numsinglehopsbyproc(:)/1000
     print *, "   doubles:"
     write(*,'(I5,A2,1000I7)') myrank,": ",numdoublehopsbyproc(:)/1000
     print *
  endif
  call mpibarrier()


  OFLWR "GOT HOPS:  "
  WRFL " Single hops this processor ",totsinglehops, " of ", totsinglewalks
  WRFL " Double hops this processor ",totdoublehops, " of ", totdoublewalks; CFL
  if (sparseconfigflag.ne.0) then
     call mympiireduceone(totsinglehops);  call mympiireduceone(totdoublehops)
     call mympiireduceone(totsinglewalks);  call mympiireduceone(totdoublewalks)
     call mympiimax(www%maxnumsinglehops);  call mympiimax(www%maxnumdoublehops)
  endif
  OFLWR " Single hops total ",totsinglehops, " of ", totsinglewalks
  WRFL " Double hops total ",totdoublehops, " of ", totdoublewalks
  WRFL "    Max single hops ", www%maxnumsinglehops
  WRFL "    Max double hops ", www%maxnumdoublehops
  WRFL; CFL

end subroutine hops

subroutine set_matsize(www)
  use walkmod
  use sparse_parameters
  implicit none
  type(walktype),intent(inout) :: www

  if (sparseconfigflag.eq.0) then
     www%singlematsize=www%numconfig
     www%doublematsize=www%numconfig
  else
     www%singlematsize=www%maxnumsinglehops
     www%doublematsize=www%maxnumdoublehops
  endif
end subroutine set_matsize


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




