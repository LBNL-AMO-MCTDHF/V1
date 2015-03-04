
#include "Definitions.INC"


subroutine print_excitations()
  use parameters
  implicit none
  integer :: i
  if (numholes.gt.0) then
     OFLWR
     do i=1,mcscfnum
        write(mpifileptr,'(A20,I3,A5,100I5)') "Holes for state ", i, "are ",myavectorhole(1:numholes,1:numholecombo,i)
     enddo
     WRFL; CFL
  endif
  if (excitations.gt.0) then
     OFLWR
     do i=1,mcscfnum
        write(mpifileptr,'(A20,I3,A5,I5,A5,I5,$)') "Excitations for state ", i, "are "
        WRFL myavectorexcitefrom(:,:,i);        WRFL " to "
        WRFL myavectorexciteto(:,:,i)
     enddo
     WRFL; CFL
  endif
end subroutine print_excitations



!!$ subroutine avector_header_read(iunit,numvects,outndof,outnumr,outnumconfig,outspinrestrictval,icomplex)
!!$   implicit none
!!$   integer :: outnumr,outnumconfig,icomplex,iunit,numvects,outndof,outspinrestrictval
!!$   read(iunit) numvects;  read(iunit) outndof,outnumr,outnumconfig,outspinrestrictval
!!$   read(iunit) icomplex
!!$ end subroutine
!!$ 
!!$ subroutine avector_header_write(iunit,numvects)
!!$   use parameters
!!$   implicit none
!!$   integer :: iunit,numvects
!!$ 
!!$   write(iunit) numvects;  write(iunit) ndof,numr,numconfig,spinrestrictval
!!$ #ifdef REALGO
!!$   write(iunit) 0
!!$ #else
!!$   write(iunit) 1
!!$ #endif
!!$ 
!!$ end subroutine


subroutine avector_header_read_simple(iunit,numvects,outndof,outnumr,outnumconfig,icomplex)
  implicit none
  integer :: outnumr,outnumconfig,icomplex,iunit,numvects,outndof,ierr,nullint

  call avector_header_read(iunit,numvects,outndof,outnumr,outnumconfig,nullint,nullint,nullint,nullint,icomplex,ierr)  
  if (ierr.ne.0) then
     call avector_header_read_old(iunit,numvects,outndof,outnumr,outnumconfig,nullint,icomplex)
  endif

end subroutine avector_header_read_simple


subroutine avector_header_read_old(iunit,numvects,outndof,outnumr,outnumconfig,outspinrestrictval,icomplex)
  implicit none
  integer :: outnumr,outnumconfig,icomplex,iunit,numvects,outndof,outspinrestrictval
  read(iunit) numvects;  read(iunit) outndof,outnumr,outnumconfig,outspinrestrictval
  read(iunit) icomplex
end subroutine


subroutine avector_header_read(iunit,numvects,outndof,outnumr,outnumconfig, &
     outrestrictflag, outrestrictms, outallspinproject, outspinrestrictval, &
     icomplex,ierr)
  implicit none
  integer :: outnumr,outnumconfig,icomplex,iunit,numvects,outndof,outspinrestrictval,ierr, &
       array(100),myiostat, outrestrictflag, outrestrictms,outallspinproject

  ierr=0

  read(iunit) numvects
  read(iunit,iostat=myiostat) array(:)
  if (myiostat.ne.0) then
     ierr=1; rewind(iunit); return
  endif
  outndof=array(1)
  outnumr=array(2)
  outnumconfig=array(3)
  outrestrictflag=array(4)
  outrestrictms=array(5)
  outallspinproject=array(6)
  outspinrestrictval=array(7)

  read(iunit) icomplex

end subroutine

subroutine avector_header_write(iunit,numvects)
  use parameters
  implicit none
  integer :: iunit,numvects, array(100)=0

  write(iunit) numvects

  array(1)=ndof
  array(2)=numr
  array(3)=numconfig
  array(4)=restrictflag
  array(5)=restrictms
  array(6)=allspinproject
  array(7)=spinrestrictval

  write(iunit) array(:)

#ifdef REALGO
  write(iunit) 0
#else
  write(iunit) 1
#endif

end subroutine







subroutine load_avectors(filename,myavectors,mynumvects,readnumvects,numskip)
  use parameters
  use mpimod
  implicit none
  character :: filename*(*)
  integer :: readnumvects,readndof,readnumr,readnumconfig,readcomplex,mynumvects,numskip
  DATATYPE :: myavectors(numconfig,numr,mynumvects)
  external :: readavectorsubroutine,readavectorsubsimple
  DATATYPE, allocatable :: readavectors(:,:,:)

  open(999,file=filename, status="unknown", form="unformatted")

!!  call avector_header_read(999,readnumvects,readndof,readnumr,readnumconfig,readspinrestrictval,readcomplex)  

 call avector_header_read_simple(999,readnumvects,readndof,readnumr,readnumconfig,readcomplex)  

  allocate(readavectors(numconfig,numr,readnumvects))

  readnumvects=min(mynumvects,readnumvects-numskip)


  close(999)

  OFL
  if (numr>readnumr) then
     write(mpifileptr, *) "WARNING, smaller numr.", numr, readnumr
  endif
  if (numconfig>readnumconfig) then 
    write(mpifileptr, *) "WARNING: smaller numconfigs.", numconfig,readnumconfig
  endif    
  if (numr<readnumr) then
     write(mpifileptr, *) "WARNING, bigger numr.", readnumr, numr
  endif
  if (numconfig<readnumconfig) then 
    write(mpifileptr, *) "WARNING: bigger numconfigs.", readnumconfig, numconfig
  endif    


  if (myrank.eq.1) then
     open(999,file=filename, status="unknown", form="unformatted")

!     call avector_header_read(999,readnumvects,readndof,readnumr,readnumconfig,readspinrestrictval,readcomplex)  

     call avector_header_read_simple(999,readnumvects,readndof,readnumr,readnumconfig,readcomplex)

     if (numholes.eq.0.and.excitations.eq.0) then

        call easy_load_avectors(999,readcomplex, readavectors(:,:,:), ndof, numr, readnumconfig, readnumvects)

!        call load_avectors0(999,readcomplex,myavectors(:,:,:),numr,numconfig,ndof           ,readnumr,readnumconfig, readavectorsubsimple,readnumvects)
     else
        call load_avectors0(999,readcomplex,readavectors(:,:,:),numr,numconfig,ndof+2*numholes,readnumr,readnumconfig, readavectorsubroutine,readnumvects)
     endif

     close(999)

     readnumvects=min(mynumvects,readnumvects-numskip)
     myavectors(:,:,1:readnumvects) = readavectors(:,:, 1+numskip : numskip+readnumvects)
     
!     myavectors(:,:,1:readnumvects)=readavectors(:,:,:)

  endif


  call mympibcast(myavectors(:,:,:),1,numconfig*numr*mynumvects)

!     OFLWR "DDDDDOT ", dot(myavectors(:,:,1),myavectors(:,:,1),totadim),dot(myavectors(:,:,2),myavectors(:,:,2),totadim); CFL
!  OFLWR "READ ", readnumvects, " VECTORS. ",numskip

end subroutine load_avectors


!! general subroutine, for loading vectors with some manipulating on read
!!     used in load_avectors (main load routine) and ovlsub.

subroutine load_avectors0(iunit, qq, myavectors, mynumr, mynumconfig, readndof, readnumr, readnumconfig, mysubroutine, mynumvects )
  use fileptrmod
  implicit none
 
  external :: mysubroutine
  integer :: mynumconfig, mynumr, mynumvects,iunit,i,  readndof, readnumr, readnumconfig
  integer :: qq, config1,  thatconfig(readndof),  myiostat, ivect
  DATATYPE :: myavectors(mynumconfig,mynumr,mynumvects),  mytempavector(mynumconfig), readvect(readnumr)
  real*8 :: rtempreadvect(readnumr)
  complex*16 :: ctempreadvect(readnumr)

  myavectors=0d0

  do ivect=1,mynumvects
     myavectors(:,:,ivect)=0d0
     do config1=1,readnumconfig
        if (qq==0) then
           read (iunit,iostat=myiostat) thatconfig(1:readndof), rtempreadvect(1:readnumr)
           readvect(:)=rtempreadvect(:)
        else
           read (iunit,iostat=myiostat) thatconfig(1:readndof), ctempreadvect(1:readnumr)
           readvect(:)=ctempreadvect(:)
        endif
        if (myiostat.ne.0) then
           OFLWR "err read config "; CFLST
        endif
!! mysubroutine takes input slater determinant with readndof electrons 
!!    mysubroutine may give a linear combo of slaters       so is programmed to return entire A-vector
!!  and checks validity and does its action and returns configuration index for myavectors (with possible +/- phase)

        call mysubroutine(thatconfig,mytempavector(:),ivect)
        do i=1,min(mynumr,readnumr)
           myavectors(:,i,ivect)= myavectors(:,i,ivect) +  readvect(i) * mytempavector(:)
        enddo
     enddo
  enddo
end subroutine load_avectors0
  

subroutine simple_load_avectors(iunit, qq, myavectors, myndof, mynumr, mynumconfig, mynumvects)
  use fileptrmod
  implicit none

  integer :: myndof, mynumconfig, mynumr,mynumvects,iunit,ivect,qq, config1, thatconfig(myndof), myiostat
  DATATYPE :: myavectors(mynumconfig,mynumr,mynumvects)
  real*8 :: rtempreadvect(mynumr)
  complex*16 :: ctempreadvect(mynumr)

  myavectors=0d0

  do ivect=1,mynumvects
     do config1=1,mynumconfig
        if (qq==0) then
           read (iunit,iostat=myiostat) thatconfig(1:myndof), rtempreadvect(1:mynumr)
           myavectors(config1,:,ivect)=rtempreadvect(:)
        else
           read (iunit,iostat=myiostat) thatconfig(1:myndof), ctempreadvect(1:mynumr)
           myavectors(config1,:,ivect)=ctempreadvect(:)
        endif
        if (myiostat.ne.0) then
           OFLWR "err read config "; CFLST
        endif
     enddo
  enddo
end subroutine simple_load_avectors


!! outavectors = numconfig

subroutine easy_load_avectors(iunit, qq, myavectors, myndof, mynumr, mynumconfig, mynumvects)
  use parameters
  implicit none

  integer :: myndof, mynumconfig, mynumr,mynumvects,iunit,ivect,qq, config1, thatconfig(myndof), myiostat,&
       phase,reorder,myconfig,getconfiguration
  DATATYPE :: myavectors(numconfig,mynumr,mynumvects)
  real*8 :: rtempreadvect(mynumr)
  logical :: allowedconfig
  complex*16 :: ctempreadvect(mynumr)

   myavectors=0d0

  do ivect=1,mynumvects
     do myconfig=1,mynumconfig
        if (qq==0) then
           read (iunit,iostat=myiostat) thatconfig(1:myndof), rtempreadvect(1:mynumr)
        else
           read (iunit,iostat=myiostat) thatconfig(1:myndof), ctempreadvect(1:mynumr)
        endif
        phase=reorder(thatconfig)
        if (allowedconfig(thatconfig)) then
           config1=getconfiguration(thatconfig)
           if (qq==0) then
              myavectors(config1,:,ivect)=rtempreadvect(:)*phase
           else
              myavectors(config1,:,ivect)=ctempreadvect(:)*phase
           endif
        else
           OFLWR "Config on file not allowed!!!"; CFLST
        endif
        if (myiostat.ne.0) then
           OFLWR "err read config "; CFLST
        endif
     enddo
  enddo
end subroutine easy_load_avectors
  


!! does conversion for initial a-vector read (excitations and holes).  N electron wfn output.

subroutine readavectorsubroutine(readconfig, outavector,ivect)
  use parameters
  use aarrmod
  implicit none
  integer :: i,iihole,iloop,jloop,ivect
  DATATYPE :: outavector(numconfig)
  integer ::  readconfig(ndof+2*numholes),xflag,iexcite,config2
  integer :: thisconfig(ndof+400)
  integer :: jj,iind,phase,reorder
  integer :: getconfiguration
  logical :: allowedconfig
  real*8 :: qfac,xphase

!!!!!!!!!!!!!!   HOLES   !!!!!!!!!!!!!!

  outavector(:)=0d0;  qfac=sqrt(1d0/numholecombo/excitecombos)

  do iloop=1,numholecombo   !! default numholecombo=1  with -1 entry for no avectorhole.
     thisconfig(1:ndof+2*numholes)=readconfig(1:ndof+2*numholes)
     xphase=1;     xflag=0

     do iihole=1,numholes
        xflag=1
        do jj=1,numelec+numholes-iihole+1
           if (iind(thisconfig(2*jj-1:2*jj)).eq.abs(myavectorhole(iihole,iloop,ivect))) then
              xflag=0
              xphase=xphase*myavectorhole(iihole,iloop,ivect)/abs(myavectorhole(iihole,iloop,ivect))
              thisconfig(2*jj-1:ndof+298)=thisconfig(2*jj+1:ndof+300)
              exit
           endif
        enddo
        if (xflag==1) then
           exit
        endif
     enddo
     if (xflag==1) then
        cycle                   !! cycle iloop
     endif
     do i=1,numelec+numholes
        thisconfig(i*2-1)   = thisconfig(i*2-1)-numloadfrozen   
     enddo

!!!!!!!!!!!!!!!!!!!      EXCITATIONS   !!!!!!!!!!!!!!!!!!!!!

     do jloop=1,excitecombos   
        xflag=0
        do iexcite=1,excitations  
           xflag=1
           do jj=1,numelec
              if (iind(thisconfig(2*jj-1:2*jj)).eq.abs(myavectorexcitefrom(iexcite,jloop,ivect))) then
                 xflag=0  
                 thisconfig(2*jj-1:2*jj)=aarr(abs(myavectorexciteto(iexcite,jloop,ivect)),nspf)
                 xphase=xphase*myavectorexciteto(iexcite,jloop,ivect)*myavectorexcitefrom(iexcite,jloop,ivect)/abs(myavectorexciteto(iexcite,jloop,ivect)*myavectorexcitefrom(iexcite,jloop,ivect))
                 exit
              endif
           enddo
           if (xflag.eq.1) then
              exit
           endif
        enddo
        if (xflag.eq.1) then
           cycle                     !! cycle jloop
        endif
        phase=reorder(thisconfig)
        if (allowedconfig(thisconfig)) then
           config2=getconfiguration(thisconfig)
           if (config2.eq.-1) then
              write(mpifileptr,*) 
              OFLWR "Hmm, looks like avector on file is ordered differently?"; CFLST
           endif
           outavector(config2) = outavector(config2) + phase * qfac * xphase
        endif
     enddo ! jloop (excitecombos)
  enddo ! iloop

end subroutine readavectorsubroutine

!! N electron wfn, no conversion

subroutine readavectorsubsimple(readconfig, outavector,nullint)
  use parameters
  use aarrmod
  implicit none
  integer :: i,nullint,phase,reorder,readconfig(ndof),config2,thisconfig(ndof),getconfiguration
  logical :: allowedconfig
  DATATYPE :: outavector(numconfig)

  outavector(:)=0d0;  thisconfig(:)=readconfig(:)
  do i=1,numelec
     thisconfig(i*2-1)   = thisconfig(i*2-1)-numloadfrozen   
  enddo
  phase=reorder(thisconfig)
  
  if (allowedconfig(thisconfig)) then
     config2=getconfiguration(thisconfig)
     if (config2.eq.-1) then
        write(mpifileptr,*) 
        OFLWR "Hmm, looks like avector on file is ordered differently?"; CFLST
     endif
     outavector(config2) = outavector(config2) + phase 
  endif
end subroutine readavectorsubsimple

subroutine write_avector(unit,avector)
  use configmod
  use parameters
  implicit none
  integer :: unit, config1
  DATATYPE :: avector(numconfig,numr)
  do config1=1,numconfig
     write (unit) configlist(:,config1), avector(config1,:)
  enddo
end subroutine write_avector

