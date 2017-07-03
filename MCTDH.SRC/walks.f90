
!! ALL ONE MODULE

!!   Initialization subroutines for walks and hops which refer to the indexing of nonzero matrix 
!!   elements between pairs of configurations (matrix elements in the Hamiltonian matrix for the 
!!   Slater determinant representation) -- initialization of derived type walktype variables.


#include "Definitions.INC"

module walksubmod
contains        

subroutine walkalloc(www)
  use fileptrmod
  use mpimod
  use walkmod
  use configsubmod 
  implicit none
  type(walktype) :: www

!! training wheels

  if (www%topconfig-www%botconfig.gt.0) then
     if (.not.highspinorder(www%configlist(:,www%topconfig),www%numpart)) then
        OFLWR "NOT HIGHSPIN",www%topconfig
        call printconfig(www%configlist(:,www%topconfig),www)
        CFLST
     endif
     if (.not.lowspinorder(www%configlist(:,www%botconfig),www%numpart)) then
        OFLWR "NOT LOWSPIN",www%botconfig
        call printconfig(www%configlist(:,www%topconfig),www)
        CFLST
     endif
  endif

!! 06-2015 configpserproc also in newconfig.f90

  allocate( www%numsinglewalks(www%configstart:www%configend+1) , &
       www%numdoublewalks(www%configstart:www%configend+1) )
  www%numsinglewalks(:)=(-1);  www%numdoublewalks(:)=(-1)
  allocate( www%numsinglediagwalks(www%configstart:www%configend+1) , &
       www%numdoublediagwalks(www%configstart:www%configend+1) )
  www%numsinglediagwalks(:)=(-1);  www%numdoublediagwalks(:)=(-1)

  call getnumwalks(www)
  OFLWR "Allocating singlewalks"; CFL
  allocate( www%singlewalk(www%maxtotsinglewalks+1) )
  www%singlewalk=-1
  allocate(www%singlediag(www%numpart,www%configstart:www%configend+1) )
  www%singlediag=-1
  allocate( www%singlewalkdirphase(www%maxtotsinglewalks+1) )
  www%singlewalkdirphase=0
  allocate( www%singlewalkopspf(1:2,www%maxtotsinglewalks+1) )
  www%singlewalkopspf=-1
  OFLWR "Allocating doublewalks"; CFL
  allocate( www%doublewalkdirspf(1:4,www%maxtotdoublewalks+1) )
  www%doublewalkdirspf=-1
  allocate( www%doublewalkdirphase(www%maxtotdoublewalks+1) )
  www%doublewalkdirphase=0
  allocate( www%doublewalk(www%maxtotdoublewalks+1) )
  www%doublewalk=-1
  allocate(www%doublediag(www%numpart*(www%numpart-1),www%configstart:www%configend+1))
  www%doublediag=-1
  OFLWR "     ..done walkalloc."; CFL

contains

  function highspinorder(thisconfig,numpart)
    implicit none
    integer,intent(in) :: numpart,thisconfig(2*numpart)
    logical :: highspinorder
    integer :: ii,unpaired(numpart),flag,jj

    highspinorder=.true.

    unpaired(1:numpart)=1

    do ii=1,numpart
       do jj=1,numpart   !! WORKS
          if (jj.ne.ii) then   !!WORKS
! -xAVX error on lawrencium!  doesnt work this way.  compiler/instruction set bug.
!     do jj=ii+1,numpart   !!FAILS
             if (thisconfig(jj*2-1).eq.thisconfig(ii*2-1)) then
                unpaired(ii)=0
                unpaired(jj)=0
             endif
          endif     !!WORKS
       enddo
    enddo
  
    flag=0
    do ii=1,numpart
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

  function lowspinorder(thisconfig,numpart)
    implicit none
    integer,intent(in) :: numpart,thisconfig(2*numpart)
    logical :: lowspinorder
    integer :: ii,unpaired(numpart),flag,jj
    
    lowspinorder=.true.

    unpaired(:)=1

    do ii=1,numpart
!     do jj=ii+1,numpart    !!FAILS
       do jj=1,numpart        !!WORKS
          if (jj.ne.ii) then  !!WORKS
             if (thisconfig(jj*2-1).eq.thisconfig(ii*2-1)) then
                unpaired(ii)=0
                unpaired(jj)=0
             endif
          endif               !!WORKS
       enddo
    enddo

    flag=0
    do ii=1,numpart
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

end subroutine walkalloc



subroutine walks(www)
  use fileptrmod
  use walkmod
  use mpimod !! nprocs
  use aarrmod
  use configsubmod
  use mpisubmod
  implicit none
  type(walktype) :: www
  integer :: iindex, iiindex, jindex, jjindex,  ispin, jspin, iispin, jjspin, ispf, jspf, &
       iispf, jjspf, config2, config1, dirphase, flag, idof, iidof, jdof, iwalk, idiag
  integer :: thisconfig(www%num2part), thatconfig(www%num2part), temporb(2), temporb2(2),&
       qsize(nprocs)  !! AUTOMATIC
  integer, allocatable :: listorder(:)

  qsize=0

  !!  ***********   SINGLES  **********

  OFLWR "Calculating walks.  Singles..."; CFL
  
  do config1=www%botconfig,www%topconfig

     if (mod(config1,1000).eq.0) then
        OFLWR config1, " out of ", www%topconfig; CFL
     endif

     iwalk=0

     if (www%singlewalkflag.ne.0) then

        thisconfig=www%configlist(:,config1)

        do idof=1,www%numpart   !! position in thisconfig that we're walking 

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
              do jdof=1,www%numpart
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

              dirphase=reorder(thatconfig,www%numpart)

              if (.not.allowedconfig0(www,thatconfig,www%dfwalklevel)) then
                 cycle
              endif

              config2=getconfiguration(thatconfig,www)
           
              if (www%configtypes(config1).ne.www%configtypes(config2)) then
                 cycle
              endif

              iwalk=iwalk+1

!! ket, bra   bra is walk
if (www%holeflag.eq.0) then
              www%singlewalkopspf(1:2,iwalk+www%scol(config1))=[ ispf,jspf ] 
else
              www%singlewalkopspf(1:2,iwalk+www%scol(config1))=[ jspf,ispf ] 
endif

              www%singlewalkdirphase(iwalk+www%scol(config1))=dirphase
           
              www%singlewalk(iwalk+www%scol(config1))=config2

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

        do idof=1,www%numpart         !! positions in thisconfig that we're walking 
           do iidof=idof+1,www%numpart   !! 

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

!! INCLUDING DIAGONAL AND SINGLE WALKS
                    flag=0
                    do jdof=1,www%numpart
                       if (jdof.ne.idof.and.jdof.ne.iidof) then 
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

                    dirphase=reorder(thatconfig,www%numpart)

                    if (.not.allowedconfig0(www,thatconfig,www%dfwalklevel)) then
                       cycle
                    endif

                    config2=getconfiguration(thatconfig,www)

                    if (www%configtypes(config1).ne.www%configtypes(config2)) then
                       cycle
                    endif

                    iwalk = iwalk+1
            
!! switched 2-2016 was                                          ket2   bra2   ket1   bra1
!!                    www%doublewalkdirspf(1:4,iwalk,config1)=[ iispf, jjspf, ispf, jspf ]

if (www%holeflag.eq.0) then
!! now                                                         bra2   ket2   bra1   ket1
                    www%doublewalkdirspf(1:4,iwalk+www%dcol(config1))=[ jjspf, iispf, jspf, ispf ]
else
                    www%doublewalkdirspf(1:4,iwalk+www%dcol(config1))=[ iispf, jjspf, ispf, jspf ]
endif

                    www%doublewalkdirphase(iwalk+www%dcol(config1))=dirphase
                 
                    www%doublewalk(iwalk+www%dcol(config1))=config2

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

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(www,nprocs)

  allocate(listorder(www%singlemaxwalks+www%doublemaxwalks+1))

  listorder=0
!$OMP DO SCHEDULE(DYNAMIC)
  do config1=www%botconfig,www%topconfig

     if (www%numsinglewalks(config1).gt.1) then
        call getlistorder(www%singlewalk(www%scol(config1)+1:),listorder(:),www%numsinglewalks(config1))
        call listreorder(www%singlewalkdirphase(www%scol(config1)+1:),listorder(:),www%numsinglewalks(config1),1)
        call listreorder(www%singlewalkopspf(:,www%scol(config1)+1:),listorder(:),www%numsinglewalks(config1),2)
        call listreorder(www%singlewalk(www%scol(config1)+1:),listorder(:),www%numsinglewalks(config1),1)
     endif

     if (www%numdoublewalks(config1).gt.1) then
        call getlistorder(www%doublewalk(www%dcol(config1)+1:),listorder(:),www%numdoublewalks(config1))
        call listreorder(www%doublewalkdirphase(www%dcol(config1)+1:),listorder(:),www%numdoublewalks(config1),1)
        call listreorder(www%doublewalkdirspf(:,www%dcol(config1)+1:),listorder(:),www%numdoublewalks(config1),4)
        call listreorder(www%doublewalk(www%dcol(config1)+1:),listorder(:),www%numdoublewalks(config1),1)
     endif

  enddo
!$OMP END DO
  deallocate(listorder)
!$OMP END PARALLEL

  OFLWR "    .... done sorting walks."; CFL


#ifdef MPIFLAG

  call mpibarrier()

  if (www%sparseconfigflag.eq.0.and.www%maxtotsinglewalks.ne.0) then

     qsize(:) = www%scol(www%alltopconfigs(:)+1) - www%scol(www%allbotconfigs(:))
     
     call mpiallgather_i(www%singlewalkopspf,   2*www%maxtotsinglewalks,&
          2*qsize(:),-00420042)
     call mpiallgather_i(www%singlewalkdirphase,www%maxtotsinglewalks,&
          qsize(:),-79800798)
     call mpiallgather_i(www%singlewalk,        www%maxtotsinglewalks,&
          qsize(:),-798042)
  endif

  if (www%sparseconfigflag.eq.0.and.www%maxtotdoublewalks.ne.0) then

     qsize(:) = www%dcol(www%alltopconfigs(:)+1) - www%dcol(www%allbotconfigs(:))

     call mpiallgather_i(www%doublewalkdirspf,  4*www%maxtotdoublewalks,&
          4*qsize(:),-9994291)
     call mpiallgather_i(www%doublewalkdirphase,www%maxtotdoublewalks,&
          qsize(:),-9994291)
     call mpiallgather_i(www%doublewalk,        www%maxtotdoublewalks,&
          qsize(:),001234)
  endif

  call mpibarrier()

#endif

  do config1=www%configstart,www%configend
     idiag=0
     do iwalk=1,www%numsinglewalks(config1)
        if (www%singlewalk(iwalk+www%scol(config1)).eq.config1) then
           idiag=idiag+1
           www%singlediag(idiag,config1)=iwalk
        endif
     enddo
     www%numsinglediagwalks(config1)=idiag
     idiag=0
     do iwalk=1,www%numdoublewalks(config1)
        if (www%doublewalk(iwalk+www%dcol(config1)).eq.config1) then
           idiag=idiag+1
           www%doublediag(idiag,config1)=iwalk
        endif
     enddo
     www%numdoublediagwalks(config1)=idiag
  enddo

contains

  recursive subroutine getlistorder(values, order,num)
    implicit none
    integer,intent(in) :: num,values(num)
    integer,intent(out) :: order(num)
    integer :: i,j,whichlowest, flag, lowval
    integer, allocatable :: taken(:)

    allocate(taken(num))
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

    deallocate(taken)

  end subroutine getlistorder

  recursive subroutine listreorder(list, order,num,numper)
    implicit none
    integer,intent(in) :: num, numper, order(num)
    integer,intent(inout) :: list(numper,num)
    integer,allocatable :: newvals(:,:)
    integer :: j

    allocate(newvals(numper,num))
    newvals=0

    do j=1,num
       newvals(:,j)=list(:,order(j))
    enddo

    list(:,:)=newvals(:,:)

    deallocate(newvals)

  end subroutine listreorder
     
  
end subroutine walks



subroutine getnumwalks(www)
  use fileptrmod
  use walkmod
  use mpimod
  use aarrmod
  use configsubmod
  use mpisubmod
  implicit none
  type(walktype) :: www
  integer :: iindex, iiindex, jindex, jjindex,  ispin, jspin, iispin, jjspin, ispf, iispf,  config1,  &
       dirphase, flag, idof, iidof, jdof,iwalk ,config2, maxwalks
  integer :: thisconfig(www%num2part), thatconfig(www%num2part), temporb(2), temporb2(2)
  integer :: totwalks, totdoublewalks, totsinglewalks
  integer*8 :: allwalks
  character(len=3) :: iilab
  character(len=4) :: iilab0

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

           do idof=1,www%numpart   !! position in thisconfig that we're walking 
           
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
                 do jdof=1,www%numpart
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

                 dirphase=reorder(thatconfig,www%numpart)

                 if (.not.allowedconfig0(www,thatconfig,www%dfwalklevel)) then
                    cycle
                 endif

                 config2=getconfiguration(thatconfig,www)

                 if (www%configtypes(config1).ne.www%configtypes(config2)) then
                    cycle
                 endif

                 iwalk=iwalk+1

              enddo   ! the walk
           enddo  ! position we're walking
        endif  ! singlewalkflag

        www%numsinglewalks(config1) = iwalk 

     enddo   ! config1

#ifdef MPIFLAG
     if (www%sparseconfigflag.eq.0) then
        call mpiallgather_i(www%numsinglewalks(:),www%numconfig,&
             www%configsperproc(:),www%maxconfigsperproc)
     endif
#endif

     OFLWR "Counting walks. Doubles"; CFL
     
  !!   ***********  DOUBLES  ************

     do config1=www%botconfig,www%topconfig

        if (mod(config1,1000).eq.0) then
           OFLWR config1, " out of ", www%topconfig;  CFL
        endif
        
        iwalk=0

        if (www%doublewalkflag.ne.0) then

           thisconfig=www%configlist(:,config1)
        
           do idof=1,www%numpart            !! positions in thisconfig that we're walking 
              do iidof=idof+1,www%numpart   !! 
              
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

!! INCLUDING DIAGONAL AND SINGLE WALKS

                       flag=0
                       do jdof=1,www%numpart
                          if (jdof.ne.idof.and.jdof.ne.iidof) then  
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
                       dirphase=reorder(thatconfig,www%numpart)

                       if (.not.allowedconfig0(www,thatconfig,www%dfwalklevel)) then
                          cycle
                       endif

                       config2=getconfiguration(thatconfig,www)

                       if (www%configtypes(config1).ne.www%configtypes(config2)) then
                          cycle
                       endif

                       iwalk = iwalk+1
                    
                    enddo   ! the walk
                 enddo
              enddo  ! position we're walking
           enddo
        endif   ! doublewalkflag

        www%numdoublewalks(config1)=iwalk

     enddo   ! config1

#ifdef MPIFLAG
     if (www%sparseconfigflag.eq.0) then
        call mpiallgather_i(www%numdoublewalks(:),www%numconfig,www%configsperproc(:),&
             www%maxconfigsperproc)
     endif
#endif

  totwalks=0;  totsinglewalks=0;  totdoublewalks=0
  www%singlemaxwalks=0; www%doublemaxwalks=0

  allocate(www%scol(www%configstart:www%configend+1), www%dcol(www%configstart:www%configend+1))

  do config1=www%configstart,www%configend

     www%scol(config1) = totsinglewalks
     www%dcol(config1) = totdoublewalks

     totwalks=totwalks+www%numsinglewalks(config1)+www%numdoublewalks(config1)
     totsinglewalks=totsinglewalks + www%numsinglewalks(config1)
     totdoublewalks=totdoublewalks + www%numdoublewalks(config1)

     if (www%singlemaxwalks.lt.www%numsinglewalks(config1)) then
        www%singlemaxwalks=www%numsinglewalks(config1)
     endif
     if (www%doublemaxwalks.lt.www%numdoublewalks(config1)) then
        www%doublemaxwalks=www%numdoublewalks(config1)
     endif

  enddo

  www%scol(www%configend+1) = totsinglewalks
  www%dcol(www%configend+1) = totdoublewalks

  www%maxtotsinglewalks = totsinglewalks
  www%maxtotdoublewalks = totdoublewalks
  maxwalks = totwalks

  allwalks=totwalks
  if (www%sparseconfigflag.ne.0) then
     call mympii8reduceone(allwalks)
     call mympiimax(www%maxtotsinglewalks);  call mympiimax(www%maxtotdoublewalks);
     call mympiimax(maxwalks)
  endif

  OFLWR;  
  WRFL    "Maximum number of"
  WRFL    "           single walks=  ",  www%maxtotsinglewalks
  WRFL    "           double walks=  ",  www%maxtotdoublewalks;  
  WRFL    "            total walks=  ",  maxwalks;  
  WRFL    "TOTAL walks:    ", allwalks
  if (www%sparseconfigflag.ne.0) then
     WRFL "maxwalks*nprocs:",int(maxwalks,8)*nprocs
  endif
  WRFL; CFL

end subroutine getnumwalks


subroutine hops(www)
  use fileptrmod
  use walkmod
  use mpimod
  use aarrmod
  use configsubmod   !! allowedconfig0
  use mpisubmod
  implicit none
  type(walktype) :: www
  integer :: ii,iwalk,iconfig,ihop,flag,iproc,isize, &
       totsinglehops,totdoublehops, totsinglewalks,totdoublewalks
  integer*8 :: allsinglehops,alldoublehops, allsinglewalks,alldoublewalks

!!$  integer :: numsinglehopsbyproc(nprocs), numdoublehopsbyproc(nprocs)

  allocate(www%numsinglehops(www%configstart:www%configend+1),&
       www%numdoublehops(www%configstart:www%configend+1))
  allocate( www%singlediaghop(www%configstart:www%configend+1),&
       www%doublediaghop(www%configstart:www%configend+1))
  allocate( www%singlehopdiagflag(www%configstart:www%configend+1),&
       www%doublehopdiagflag(www%configstart:www%configend+1))

  www%numsinglehops(:)=(-99);  www%numdoublehops(:)=(-99)
  www%singlediaghop(:)=(-99);  www%doublediaghop(:)=(-99)
  www%singlehopdiagflag(:)=(-99);  www%doublehopdiagflag(:)=(-99)

  allocate( www%firstsinglehopbyproc(nprocs,www%configstart:www%configend+1), &
       www%lastsinglehopbyproc(nprocs,www%configstart:www%configend+1) )
  allocate( www%firstdoublehopbyproc(nprocs,www%configstart:www%configend+1), &
       www%lastdoublehopbyproc(nprocs,www%configstart:www%configend+1) )

  www%firstsinglehopbyproc(:,:)=(-99);  www%lastsinglehopbyproc(:,:)=(-99)
  www%firstdoublehopbyproc(:,:)=(-99);  www%lastdoublehopbyproc(:,:)=(-99)

  do ii=0,1

     if (ii.eq.0) then
!! avoid warn bounds
        allocate(www%singlehop(1,1),www%singlehopwalkstart(1,1),www%singlehopwalkend(1,1),&
             www%doublehop(1,1),www%doublehopwalkstart(1,1),www%doublehopwalkend(1,1))
     else
        deallocate(www%singlehop,www%singlehopwalkstart,www%singlehopwalkend,&
             www%doublehop,www%doublehopwalkstart,www%doublehopwalkend)
        allocate(www%singlehop(www%maxnumsinglehops,www%configstart:www%configend+1),&
             www%singlehopwalkstart(www%maxnumsinglehops,www%configstart:www%configend+1),&
             www%singlehopwalkend(www%maxnumsinglehops,www%configstart:www%configend+1),&
             www%doublehop(www%maxnumdoublehops,www%configstart:www%configend+1),&
             www%doublehopwalkstart(www%maxnumdoublehops,www%configstart:www%configend+1),&
             www%doublehopwalkend(www%maxnumdoublehops,www%configstart:www%configend+1))
        www%singlehop(:,:)=(-99)
        www%singlehopwalkstart(:,:)=1
        www%singlehopwalkend(:,:)=0
        www%doublehop(:,:)=(-99)
        www%doublehopwalkstart(:,:)=1
        www%doublehopwalkend(:,:)=0
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
              www%singlehop(1,iconfig)=www%singlewalk(1+www%scol(iconfig))
              www%singlehopwalkstart(1,iconfig)=1
           endif
           do iwalk=2,www%numsinglewalks(iconfig)
              if (www%singlewalk(iwalk+www%scol(iconfig)).ne.www%singlewalk(iwalk-1+www%scol(iconfig))) then
                 if (ii.eq.1) then
                    www%singlehopwalkend(ihop,iconfig)=iwalk-1
                 endif
                 ihop=ihop+1
                 if (ii.eq.1) then
                    www%singlehop(ihop,iconfig)=www%singlewalk(iwalk+www%scol(iconfig))
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
              www%doublehop(1,iconfig)=www%doublewalk(1+www%dcol(iconfig))
              www%doublehopwalkstart(1,iconfig)=1
           endif
           do iwalk=2,www%numdoublewalks(iconfig)
              if (www%doublewalk(iwalk+www%dcol(iconfig)).ne.www%doublewalk(iwalk-1+www%dcol(iconfig))) then
                 if (ii.eq.1) then
                    www%doublehopwalkend(ihop,iconfig)=iwalk-1
                 endif
                 ihop=ihop+1
                 if (ii.eq.1) then
                    www%doublehop(ihop,iconfig)=www%doublewalk(iwalk+www%dcol(iconfig))
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

  enddo   !! do ii=0,1


  do iconfig=www%configstart,www%configend

     if (allowedconfig0(www,www%configlist(:,iconfig),www%dfwalklevel) .and. www%holeflag.ne.0) then
        if (www%numsinglehops(iconfig).eq.0.and.www%singlewalkflag.ne.0) then
           www%numsinglehops(iconfig) = 1
           www%singlehop(1,iconfig)=iconfig
        endif
        if (www%numdoublehops(iconfig).eq.0.and.www%doublewalkflag.ne.0) then
           www%numdoublehops(iconfig)=1
           www%doublehop(1,iconfig)=iconfig
        endif
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

#ifdef MPIFLAG
  call mpibarrier()

  if (www%sparseconfigflag.eq.0) then
     isize=nprocs
     call mpiallgather_i(www%firstsinglehopbyproc(:,:),   www%numconfig*isize,&
          www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
     call mpiallgather_i(www%lastsinglehopbyproc(:,:),   www%numconfig*isize,&
          www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
     call mpiallgather_i(www%firstdoublehopbyproc(:,:),   www%numconfig*isize,&
          www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
     call mpiallgather_i(www%lastdoublehopbyproc(:,:),   www%numconfig*isize,&
          www%configsperproc(:)*isize,www%maxconfigsperproc*isize)
  endif

  if (www%sparseconfigflag.eq.0) then
     ii=1
     call mpiallgather_i(www%numdoublehops(:),www%numconfig*ii,www%configsperproc(:)*ii,&
          www%maxconfigsperproc*ii)
     call mpiallgather_i(www%numsinglehops(:),www%numconfig*ii,www%configsperproc(:)*ii,&
          www%maxconfigsperproc*ii)
     call mpiallgather_i(www%singlediaghop(:),www%numconfig*ii,www%configsperproc(:)*ii,&
          www%maxconfigsperproc*ii)
     call mpiallgather_i(www%doublediaghop(:),www%numconfig*ii,www%configsperproc(:)*ii,&
          www%maxconfigsperproc*ii)
     call mpiallgather_i(www%singlehopdiagflag(:),www%numconfig*ii,www%configsperproc(:)*ii,&
          www%maxconfigsperproc*ii)
     call mpiallgather_i(www%doublehopdiagflag(:),www%numconfig*ii,www%configsperproc(:)*ii,&
          www%maxconfigsperproc*ii)

     ii=www%maxnumsinglehops
     call mpiallgather_i(www%singlehop(:,:),www%numconfig*ii,www%configsperproc(:)*ii,&
          www%maxconfigsperproc*ii)
     call mpiallgather_i(www%singlehopwalkstart(:,:),www%numconfig*ii,www%configsperproc(:)*ii,&
          www%maxconfigsperproc*ii)
     call mpiallgather_i(www%singlehopwalkend(:,:),www%numconfig*ii,www%configsperproc(:)*ii,&
          www%maxconfigsperproc*ii)
     ii=www%maxnumdoublehops
     call mpiallgather_i(www%doublehop(:,:),www%numconfig*ii,www%configsperproc(:)*ii,&
          www%maxconfigsperproc*ii)
     call mpiallgather_i(www%doublehopwalkstart(:,:),www%numconfig*ii,www%configsperproc(:)*ii,&
          www%maxconfigsperproc*ii)
     call mpiallgather_i(www%doublehopwalkend(:,:),www%numconfig*ii,www%configsperproc(:)*ii,&
          www%maxconfigsperproc*ii)
  endif

!!$  numsinglehopsbyproc(:)=0;   numdoublehopsbyproc(:)=0
!!$
!!$  do iconfig=www%botconfig,www%topconfig
!!$     numsinglehopsbyproc(:)=numsinglehopsbyproc(:) + &
!!$          (www%lastsinglehopbyproc(:,iconfig)-www%firstsinglehopbyproc(:,iconfig)+1)
!!$     numdoublehopsbyproc(:)=numdoublehopsbyproc(:) + &
!!$          (www%lastdoublehopbyproc(:,iconfig)-www%firstdoublehopbyproc(:,iconfig)+1)
!!$  enddo
!!$
!!$  call mpibarrier()
!!$  if (myrank.eq.1) then
!!$     print *, "HOPS BY PROC ON PROCESSOR 1 :::::::::::::::::::::::"
!!$     print *, "   singles:"
!!$     write(*,'(I5,A2,1000I7)') myrank,": ",numsinglehopsbyproc(:)/1000
!!$     print *, "   doubles:"
!!$     write(*,'(I5,A2,1000I7)') myrank,": ",numdoublehopsbyproc(:)/1000
!!$     print *
!!$  endif

  call mpibarrier()

#endif

  OFLWR "GOT HOPS:  "
  WRFL " Single hops this processor ",totsinglehops, " of ", totsinglewalks
  WRFL " Double hops this processor ",totdoublehops, " of ", totdoublewalks; CFL

  allsinglehops=totsinglehops
  allsinglewalks=totsinglewalks
  alldoublehops=totdoublehops
  alldoublewalks=totdoublewalks

#ifdef MPIFLAG
  if (www%sparseconfigflag.ne.0) then
     call mympii8reduceone(allsinglehops);  call mympii8reduceone(alldoublehops)
     call mympii8reduceone(allsinglewalks);  call mympii8reduceone(alldoublewalks)
     call mympiimax(www%maxnumsinglehops);  call mympiimax(www%maxnumdoublehops)
  endif
#endif

  OFLWR " Single hops total ",allsinglehops, " of ", allsinglewalks
  WRFL " Double hops total ",alldoublehops, " of ", alldoublewalks
  WRFL "    Max single hops ", www%maxnumsinglehops
  WRFL "    Max double hops ", www%maxnumdoublehops
  WRFL; CFL

end subroutine hops

subroutine set_matsize(www)
  use walkmod
  implicit none
  type(walktype),intent(inout) :: www

  if (www%sparseconfigflag.eq.0) then
     www%singlematsize=www%numconfig
     www%doublematsize=www%numconfig
  else
     www%singlematsize=www%maxnumsinglehops
     www%doublematsize=www%maxnumdoublehops
  endif
end subroutine set_matsize


! deallocating doublewalks and singlewalks ahead of time
!subroutine walkdealloc(www)
!  use walkmod
!  implicit none
!  type(walktype) :: www
!  deallocate( www%numsinglewalks,www%numsinglediagwalks )
!  deallocate( www%numdoublewalks,www%numdoublediagwalks )
!  deallocate( www%singlewalk )
!  deallocate( www%singlewalkdirphase )
!  deallocate( www%singlewalkopspf )
!  deallocate( www%doublewalkdirspf )
!  deallocate( www%doublewalkdirphase )
!  deallocate( www%doublewalk)
!end subroutine walkdealloc


end module walksubmod


