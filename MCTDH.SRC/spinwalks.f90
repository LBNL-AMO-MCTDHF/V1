
!! ALL ONE MODULE

!! SUBROUTINE FOR "WALKS" (WHICH CONFIGURATIONS CONNECT TO WHICH) FOR SPARSE SPIN PROJECTOR.

#include "Definitions.INC"

module spinwalkinternal
  implicit none
  integer, allocatable,private :: numspinwalks(:),spinwalkdirphase(:,:),spinwalk(:,:),&
       unpaired(:,:),msvalue(:),numunpaired(:)
  real*8, allocatable,private :: configspinmatel(:,:)
  integer,private ::   maxspinwalks=0

contains

subroutine spinwalkinternal_dealloc()
  implicit none
  deallocate(unpaired,numunpaired,msvalue,numspinwalks,spinwalk, spinwalkdirphase, configspinmatel)
end subroutine spinwalkinternal_dealloc



subroutine spinwalkinit(www)
  use fileptrmod
  use walkmod
  implicit none
  type(walktype),intent(inout) :: www

  call mpibarrier()
  OFLWR "Go spinwalk init. "; CFL
  call mpibarrier()

  allocate(www%sss%numspinsets(www%startrank:www%endrank+1),&
       www%sss%numspindfsets(www%startrank:www%endrank+1))
  www%sss%numspinsets(:)=(-99);   www%sss%numspindfsets(:)=(-99)

  allocate(unpaired(www%numpart,www%configstart:www%configend+1), &
       numunpaired(www%configstart:www%configend+1), &
       msvalue(www%configstart:www%configend+1), &
       numspinwalks(www%configstart:www%configend+1))
  unpaired(:,:)=(-99);   numunpaired(:)=(-99)
  msvalue(:)=(-99);   numspinwalks(:)=(-99)

  call mpibarrier()
  OFLWR "Go get numspinwalks."; CFL
  
  call getnumspinwalks(www)

  allocate(spinwalk(maxspinwalks,www%configstart:www%configend+1), &
       spinwalkdirphase(maxspinwalks,www%configstart:www%configend+1),&
       configspinmatel(maxspinwalks+1,www%configstart:www%configend+1),&
       www%sss%spinsetsize(www%maxconfigsperproc,www%startrank:www%endrank+1), &
       www%sss%spinsetrank(www%maxconfigsperproc,www%startrank:www%endrank+1))

  spinwalk(:,:)=(-99);  spinwalkdirphase(:,:)=(-99)
  configspinmatel(:,:)=(-99)
  www%sss%spinsetsize(:,:)=(-99);   www%sss%spinsetrank(:,:)=(-99)

  call mpibarrier()

  OFLWR "Done spinwalk init."; CFL

end subroutine spinwalkinit


!subroutine spinwalkdealloc(sss)
!  use spinwalkmod
!  implicit none
!  type(spintype) :: sss
!  integer :: i,iproc
!
!  deallocate(sss%spinsets, sss%spinsetsize, sss%spinsetrank)
!
!  do iproc=www%startrank,www%endrank             !! oops
!     do i=1,sss%numspinsets(iproc)
!        deallocate(sss%spinsetprojector(i,iproc)%mat,sss%spinsetprojector(i,iproc)%vects)
!     enddo
!  enddo
!  deallocate(sss%spinsetprojector)
!
!end subroutine spinwalkdealloc



subroutine spinwalks(www)
  use fileptrmod
  use walkmod
  use aarrmod
  use configsubmod
  implicit none
  type(walktype),intent(in) :: www
  integer ::     config1, config2, dirphase,  idof, jdof,iwalk , thisconfig(www%num2part),  &
       thatconfig(www%num2part), ii, jj, firstspin, secondspin, iproc

  spinwalk=0;     spinwalkdirphase=0


  OFLWR "Calculating spin walks.";  CFL

  do iproc=www%startrank,www%endrank

     do config1=www%allbotconfigs(iproc),www%alltopconfigs(iproc)

        iwalk=0
        do ii=1,numunpaired(config1)
           thisconfig=www%configlist(:,config1);        idof=unpaired(ii,config1)
           if (idof==0) then
              OFLWR "Unpaired error";CFLST
           endif
           firstspin=thisconfig(idof*2)
           thisconfig(idof*2)=mod(thisconfig(idof*2),2) + 1
           
           do jj=ii+1,numunpaired(config1)
              thatconfig=thisconfig;           jdof=unpaired(jj,config1)
              if (jdof==0) then
                 OFLWR "UUnpaired error"; CFLST
              endif
              secondspin=thatconfig(jdof*2)

              if (secondspin.ne.firstspin) then
                 thatconfig(jdof*2)=mod(thatconfig(jdof*2),2) + 1
                 dirphase=reorder(thatconfig,www%numpart)
                 if (allowedconfig0(www,thatconfig,www%dflevel)) then
                    iwalk=iwalk+1
                    spinwalkdirphase(iwalk,config1)=dirphase
                    config2=getconfiguration(thatconfig,www)
                    spinwalk(iwalk,config1)=config2

  if ((config2.lt.www%allbotconfigs(iproc).or.config2.gt.www%alltopconfigs(iproc))) then
     OFLWR "BOT TOP NEWCONFIG ERR", config1,config2,iproc,www%allbotconfigs(iproc),&
          www%alltopconfigs(iproc);     CFLST
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





subroutine spinsets_first(www)
  use fileptrmod
  use walkmod
  use aarrmod
  use mpimod
  use basissubmod
  use mpisubmod
  implicit none
  type(walktype),intent(inout) :: www
  integer ::  iwalk, jj,  iset, ilevel, currentnumwalks, prevnumwalks, flag, iflag, &
       addwalks, i, j, jwalk, jset, iproc
  integer, allocatable :: taken(:), tempwalks(:)

  OFLWR "   ... go spinsets ..."; CFL

  allocate(taken(www%configstart:www%configend+1), tempwalks(www%maxconfigsperproc))
  
  www%sss%maxspinsetsize=0;    tempwalks=0

  do jj=0,1

     taken(:)=0

     do iproc=www%startrank,www%endrank

        iset=0;    jset=0

        do i=www%allbotconfigs(iproc),www%alltopconfigs(iproc)
           if (taken(i).ne.1) then
              if (www%ddd%dfincludedmask(i).ne.0) then
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
                 www%sss%spinsetsize(iset,iproc)=currentnumwalks
                 if (www%sss%maxspinsetsize .lt. currentnumwalks) then
                    www%sss%maxspinsetsize=currentnumwalks
                 endif
              else
                 if ((currentnumwalks.gt.www%sss%maxspinsetsize).or.&
                      (www%sss%spinsetsize(iset,iproc).ne.currentnumwalks)) then
                    OFLWR "WALK ERROR";CFLST
                 endif
                 www%sss%spinsets(1:currentnumwalks,iset,iproc)=tempwalks(1:currentnumwalks)
                 if (www%ddd%dfincludedmask(i).ne.0) then
                    www%sss%spindfsetindex(jset,iproc)=iset
                    do j=1,currentnumwalks
                       www%sss%spindfsets(j,jset,iproc)=getdfindex(www,tempwalks(j))
                    enddo
                 endif
              endif
           endif
        enddo
        if (jj==0) then
           www%sss%numspinsets(iproc)=iset
           www%sss%numspindfsets(iproc)=jset
        else
           if (www%sss%numspinsets(iproc).ne.iset) then
              OFLWR "NUMSPINSETS ERROR ", www%sss%numspinsets(iproc), iset;CFLST
           endif
           if (www%sss%numspindfsets(iproc).ne.jset) then
              OFLWR "NUMSPINdfSETS ERROR ", www%sss%numspindfsets(iproc), jset;CFLST
           endif
        endif
     enddo  !! iproc = startrank,endrank

     if (www%sparseconfigflag.ne.0) then
        call mympiimax(www%sss%maxspinsetsize)
     endif

     do i=www%configstart,www%configend
        if (taken(i).ne.1) then
           print *, "TAKEN ERROR!!!!", myrank,i,taken(i),www%startrank,www%endrank; stop
        endif
     enddo

     j=0
     www%sss%maxnumspinsets=0; www%sss%maxnumspindfsets=0
     do iproc=www%startrank,www%endrank
        if (www%sss%maxnumspinsets.lt.www%sss%numspinsets(iproc)) then
           www%sss%maxnumspinsets=www%sss%numspinsets(iproc)
        endif
        if (www%sss%maxnumspindfsets.lt.www%sss%numspindfsets(iproc)) then
           www%sss%maxnumspindfsets=www%sss%numspindfsets(iproc)
        endif
        do i=1,www%sss%numspinsets(iproc)
           j=j+www%sss%spinsetsize(i,iproc)
        enddo
     enddo
     if (www%sparseconfigflag.ne.0) then
        call mympiimax(www%sss%maxnumspinsets)
        call mympiimax(www%sss%maxnumspindfsets)
        call mympiireduceone(j)
     endif
     if (j.ne.www%numconfig) then
        OFLWR "SPINSETSIZE ERROR!! ", j, www%numconfig; CFLST
     endif

     if (jj==0) then
        allocate(www%sss%spinsets(www%sss%maxspinsetsize,www%sss%maxnumspinsets,www%startrank:www%endrank+1), &
             www%sss%spindfsets(www%sss%maxspinsetsize,www%sss%maxnumspindfsets,www%startrank:www%endrank+1),&
             www%sss%spindfsetindex(www%sss%maxnumspindfsets,www%startrank:www%endrank+1))
        www%sss%spinsets(:,:,:)=(-99)
        www%sss%spindfsets(:,:,:)=(-99)
        www%sss%spindfsetindex(:,:)=(-99)
     endif

  enddo  !! jj

  deallocate(taken, tempwalks)

  OFLWR "Numspinsets this processor is ", www%sss%numspinsets(myrank),"  maxspinset size is ", www%sss%maxspinsetsize
  CFL

end subroutine spinsets_first




subroutine getnumspinwalks(www)
  use fileptrmod
  use walkmod
  use configsubmod
  use mpisubmod
  implicit none
  type(walktype),intent(in) :: www
  integer ::   ispf,  config1, flag, idof, jdof,iwalk , thisconfig(www%num2part), &
       thatconfig(www%num2part), ii, jj, dirphase, firstspin, secondspin
  real*8 :: avgspinwalks

  OFLWR "Doing spin projector.";  CFL

  do config1=www%configstart,www%configend
     unpaired(:,config1)=0;    numunpaired(config1)=0;   msvalue(config1)=0
     thisconfig=www%configlist(:,config1)
     do idof=1,www%numpart
        msvalue(config1)=msvalue(config1) + thisconfig(idof*2)*2-3
        ispf=thisconfig(idof*2-1)
        flag=0
        do jdof=1,www%numpart
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

  do config1=www%configstart,www%configend
     iwalk=0
     do ii=1,numunpaired(config1)
        thisconfig=www%configlist(:,config1)
        idof=unpaired(ii,config1)
        if (idof==0) then
           OFLWR "Unpaired error"; CFLST
        endif
        firstspin=thisconfig(idof*2)
        thisconfig(idof*2)=mod(thisconfig(idof*2),2) + 1
        do jj=ii+1,numunpaired(config1)
           thatconfig=thisconfig
           jdof=unpaired(jj,config1)
           if (jdof==0) then
              OFLWR "Unpaired error"; CFLST
           endif
           secondspin=thatconfig(jdof*2)
           if (secondspin.ne.firstspin) then
              thatconfig(jdof*2)=mod(thatconfig(jdof*2),2) + 1
              dirphase=reorder(thatconfig,www%numpart)
              if (allowedconfig0(www,thatconfig,www%dflevel)) then
                 iwalk=iwalk+1
              endif   ! allowedconfig
           endif
        enddo   ! jj
     enddo  ! ii
     numspinwalks(config1) = iwalk 
  enddo   ! config1


  maxspinwalks=0
  avgspinwalks=0.d0
  do config1=www%configstart,www%configend
     avgspinwalks = avgspinwalks + numspinwalks(config1)
     if (maxspinwalks.lt.numspinwalks(config1)) then
        maxspinwalks=numspinwalks(config1)
     endif
  enddo

  if (www%sparseconfigflag.ne.0) then
     call mympirealreduceone(avgspinwalks)
     call mympiimax(maxspinwalks)
  endif

  avgspinwalks=avgspinwalks/www%numconfig

  OFLWR "Maximum number of spin walks= ",  maxspinwalks
  write(mpifileptr, *) "Avg number of spin walks= ",  avgspinwalks;CFL
  
end subroutine getnumspinwalks





subroutine configspin_matel(www)   
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer ::     config2, config1,   iwalk, myind

  configspinmatel(:,:)=0.d0

  do config1=www%configstart,www%configend
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



function spinallowed(spinval,sss)
  use spinwalkmod
  implicit none
  type(spintype),intent(in) :: sss
  real*8,intent(in) :: spinval
  logical :: spinallowed
  if (abs(spinval-(sss%spinrestrictval/2.d0*(sss%spinrestrictval/2.d0+1))).lt.1.d-3) then
     spinallowed=.true.
  else
     spinallowed=.false.
  endif
end function



subroutine configspinset_projector(www)   
  use fileptrmod
  use walkmod
  use mpimod
  use configsubmod
  use mpisubmod
  implicit none
  type(walktype),intent(inout) :: www
  integer :: info, lwork,j,i,ii,iset,jj, elim, elimsets, flag, iwalk,&
       iproc
  real*8, allocatable :: spinvects(:,:), spinvals(:), work(:), realprojector(:,:)
!  DATATYPE :: doublevects(maxspinsetsize**2)
!  real*8 :: doublevals(maxspinsetsize)

  call mpibarrier()
  OFLWR "Getting Spinset Projectors.  Numspinsets this process is ", www%sss%numspinsets(myrank);CFL
  call mpibarrier()
  allocate(www%sss%spinsetprojector(www%sss%maxnumspinsets,www%startrank:www%endrank+1))
  allocate(spinvects(www%sss%maxspinsetsize,www%sss%maxspinsetsize), spinvals(www%sss%maxspinsetsize), &
       realprojector(www%sss%maxspinsetsize,www%sss%maxspinsetsize))
  if (www%sss%maxspinsetsize.gt.0) then
     spinvects=0; spinvals=0; realprojector=0
  endif

  lwork=10*www%sss%maxspinsetsize;  allocate(work(lwork));    work=0
  
  allocate(www%sss%spinsperproc(nprocs),www%sss%spindfsperproc(nprocs), &
       www%sss%allbotspins(nprocs),www%sss%alltopspins(nprocs),&
       www%sss%allbotspindfs(nprocs), www%sss%alltopspindfs(nprocs))

  www%sss%spinsperproc(:)=0; www%sss%spindfsperproc(:)=0
  www%sss%allbotspins=0; www%sss%alltopspins=0
  www%sss%allbotspindfs=0; www%sss%alltopspindfs=0

  elim=0;  elimsets=0;  

  call mpibarrier()
  OFLWR "                                       maxspinsetsize is ", www%sss%maxspinsetsize; CFL
  call mpibarrier()

  do iproc=www%startrank,www%endrank

     iset=1
  
     do while (iset.le.www%sss%numspinsets(iproc))
        spinvects=0.d0
        do ii=1,www%sss%spinsetsize(iset,iproc)
           do jj=1,www%sss%spinsetsize(iset,iproc)
              
              if (ii.eq.jj) then
                 spinvects(ii,jj)=configspinmatel(1, www%sss%spinsets(jj,iset,iproc))
              else
                 flag=0
                 do iwalk=1,numspinwalks(www%sss%spinsets(jj,iset,iproc))
                    if (spinwalk(iwalk,www%sss%spinsets(jj,iset,iproc)).eq.www%sss%spinsets(ii,iset,iproc)) then
                       spinvects(ii,jj)=configspinmatel(iwalk+1, www%sss%spinsets(jj,iset,iproc))
                       flag=1
                       exit
                    endif
                 enddo
              endif
           enddo
        enddo
     
        call dsyev('V','U', www%sss%spinsetsize(iset,iproc), spinvects, &
             www%sss%maxspinsetsize, spinvals, work, lwork, info)
        if (info/=0) then
           OFLWR  "INFO SSYEV", info; CFLST
        endif
        j=0; 
        do i=1,www%sss%spinsetsize(iset,iproc)
           if (spinallowed(spinvals(i),www%sss)) then
              j=j+1;           spinvects(:,j)=spinvects(:,i)
           endif
        enddo
        www%sss%spinsetrank(iset,iproc)=j

        www%sss%spinsperproc(iproc)=www%sss%spinsperproc(iproc)+j

        if (allowedconfig0(www,www%configlist(:,www%sss%spinsets(1,iset,iproc)),www%dfrestrictflag)) then
           www%sss%spindfsperproc(iproc)=www%sss%spindfsperproc(iproc)+j
        endif

        spinvects(:,j+1:www%sss%maxspinsetsize)=0d0
        
        if (www%sss%spinsetrank(iset,iproc)==0) then 
           elimsets=elimsets+1
           elim=elim+www%sss%spinsetsize(iset,iproc)
           www%sss%spinsetsize(iset:www%sss%numspinsets(iproc)-1,iproc)=&
                www%sss%spinsetsize(iset+1:www%sss%numspinsets(iproc),iproc)
           www%sss%spinsets(:,iset:www%sss%numspinsets(iproc)-1,iproc)=&
                www%sss%spinsets(:,iset+1:www%sss%numspinsets(iproc),iproc)
           www%sss%numspinsets(iproc)=www%sss%numspinsets(iproc)-1
        else
           allocate(www%sss%spinsetprojector(iset,iproc)%mat(www%sss%spinsetsize(iset,iproc),&
                www%sss%spinsetsize(iset,iproc)))
           allocate(www%sss%spinsetprojector(iset,iproc)%vects(www%sss%spinsetsize(iset,iproc), &
                www%sss%spinsetrank(iset,iproc)))
           
           www%sss%spinsetprojector(iset,iproc)%vects(:,:)=spinvects(1:www%sss%spinsetsize(iset,iproc), &
                1:www%sss%spinsetrank(iset,iproc))
           
           call dgemm('N', 'T', www%sss%spinsetsize(iset,iproc), www%sss%spinsetsize(iset,iproc),&
                www%sss%spinsetrank(iset,iproc),1.0d0, spinvects, www%sss%maxspinsetsize, &
                spinvects, www%sss%maxspinsetsize ,0.0d0,realprojector, www%sss%maxspinsetsize)
           
           www%sss%spinsetprojector(iset,iproc)%mat(:,:) = realprojector(1:www%sss%spinsetsize(iset,iproc), &
                1:www%sss%spinsetsize(iset,iproc))

!just checking right
!        call CONFIGHERM(www%sss%spinsetprojector(iset,iproc)%mat,www%sss%spinsetsize(iset,iproc),&
!              www%sss%spinsetsize(iset,iproc), doublevects,doublevals)
!        do i=1,www%sss%spinsetsize(iset,iproc)
!           if (abs(doublevals(i)*2-1.d0)-1.d0 .gt. 1.d-9) then
!              OFLWR "SPIN PROJECTOR ERROR", doublevals(i); CFLST
!           endif
!        enddo

           iset=iset+1

        endif
     enddo

  enddo

  deallocate(spinvects, spinvals, realprojector, work)
  call mpibarrier()

  if (www%sparseconfigflag.ne.0) then
     call mympiireduceone(elim)
     call mympiireduceone(elimsets)
  endif

  call mpibarrier()
  OFLWR "...done.  Eliminated ", elimsets, " sets with total rank ", elim   ; CFL
  call mpibarrier()

  do iproc=www%startrank,www%endrank
     i=0
     do ii=1,www%sss%numspinsets(iproc)
        i=i+www%sss%spinsetrank(ii,iproc)
     enddo
     if (i.ne.www%sss%spinsperproc(iproc)) then
        print *, " SPIN CHECKMEERRR",myrank,iproc,i; call mpistop()
     endif
  enddo

  if (www%sparseconfigflag.ne.0) then
     do iproc=1,nprocs
        call mympiibcastone(www%sss%spinsperproc(iproc),iproc)
        call mympiibcastone(www%sss%spindfsperproc(iproc),iproc)
     enddo
  endif

  www%sss%maxspinsperproc=0
  ii=0
  do i=1,nprocs
     www%sss%allbotspins(i)=ii+1
     ii=ii+www%sss%spinsperproc(i)
     www%sss%alltopspins(i)=ii
     if (www%sss%spinsperproc(i).gt.www%sss%maxspinsperproc) then
        www%sss%maxspinsperproc=www%sss%spinsperproc(i)
     endif
  enddo
  www%sss%numspinconfig=ii
  www%sss%botspin=www%sss%allbotspins(myrank)
  www%sss%topspin=www%sss%alltopspins(myrank)

!! SPIN DF

  www%sss%maxspindfsperproc=0
  ii=0
  do i=1,nprocs
     www%sss%allbotspindfs(i)=ii+1
     ii=ii+www%sss%spindfsperproc(i)
     www%sss%alltopspindfs(i)=ii
     if (www%sss%spindfsperproc(i).gt.www%sss%maxspindfsperproc) then
        www%sss%maxspindfsperproc=www%sss%spindfsperproc(i)
     endif
  enddo
  www%sss%numspindfconfig=ii
  www%sss%botdfspin=www%sss%allbotspindfs(myrank)
  www%sss%topdfspin=www%sss%alltopspindfs(myrank)

  if (www%parconsplit.eq.0) then
     www%sss%firstspinconfig=1
     www%sss%lastspinconfig=www%sss%numspinconfig
     www%sss%localnspin=www%sss%numspinconfig
  else
     www%sss%firstspinconfig=www%sss%botspin
     www%sss%lastspinconfig=www%sss%topspin
     www%sss%localnspin=www%sss%spinsperproc(myrank)
  endif

  call mpibarrier()
  OFLWR
  WRFL  "This processor: "
  WRFL "      spin sets        ", www%sss%numspinsets(myrank)
  WRFL "      spin rank            ", www%sss%spinsperproc(myrank), " of ", www%configsperproc(myrank) 
  WRFL "      botspin,topspin:     ", www%sss%botspin,www%sss%topspin
  WRFL "      df spin rank         ", www%sss%spindfsperproc(myrank), " of ", www%dfconfsperproc(myrank) 
  WRFL "      botdfspin,topdfspin: ", www%sss%botdfspin,www%sss%topdfspin
  WRFL "All processors:"
  WRFL "      spin rank, S^2 = ", (www%sss%spinrestrictval/2.d0*(www%sss%spinrestrictval/2.d0+1)), &
       " is ", www%sss%numspinconfig, " out of ", www%numconfig
  WRFL "      df spin rank         ", www%sss%numspindfconfig, " of ", www%numdfconfigs
  WRFL; CFL
  call mpibarrier()
  
end subroutine configspinset_projector

end module spinwalkinternal

