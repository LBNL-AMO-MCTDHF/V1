
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
  array(3)=num_config
  array(4)=restrictflag
  array(5)=restrict_ms
  array(6)=all_spinproject
  array(7)=spin_restrictval

  write(iunit) array(:)

#ifdef REALGO
  write(iunit) 0
#else
  write(iunit) 1
#endif

end subroutine



!! PARCONSPLIT=0 ONLY RIGHT NOW

subroutine load_avector_productsub(myavector)
  use parameters
  use configmod
  use mpimod
  implicit none
  DATATYPE,intent(out) :: myavector(numr,first_config:last_config,mcscfnum)
  integer :: readnumvects(numavectorfiles),readndof(numavectorfiles),readnumr(numavectorfiles),&
       readnumconfig(numavectorfiles),readcomplex(numavectorfiles),readunit(numavectorfiles)
  integer :: configtable(num_config),configphase(num_config) !! AUTOMATIC
  type threemat
     DATATYPE, allocatable :: mat(:,:,:)
  end type threemat
  TYPE(threemat) :: readavectors(numavectorfiles)
  type inttwomat
     integer, allocatable :: mat(:,:)
  end type inttwomat
  type(inttwomat) :: readconfiglist(numavectorfiles)
  DATATYPE, allocatable :: productvector(:,:,:,:,:,:,:),productreshape(:,:)
  integer,allocatable :: newconfiglist(:,:)
  integer :: ifile,tot_ndof,tot_numconfig,iconfig,jj,kk,dirphase,reorder,num_allowed,jconfig,iitop(6),&
       dofsum,thisconfig(ndof),ir,getconfiguration,tot_wfns,iwfn
  logical :: allowedconfig0
  integer, target :: ii(6)
  integer, pointer :: ii1,ii2,ii3,ii4,ii5,ii6
!  DATATYPE :: dot

  if (numspffiles.ne.numavectorfiles) then
     OFLWR "numavectorfiles and numspffiles should be the same for load_avector_product"; CFLST
  endif
  do ifile=1,numavectorfiles
     if (eachloaded(ifile).le.0) then
        OFLWR "Error Eachloaded:", eachloaded(1:numavectorfiles); CFLST
     endif
  enddo

  if (myrank.eq.1) then
     do ifile=1,numavectorfiles
        readunit(ifile)=677+ifile
        open(readunit(ifile),file=avectorfile(ifile), status="unknown", form="unformatted")
        call avector_header_read_simple(readunit(ifile),readnumvects(ifile),readndof(ifile),readnumr(ifile),readnumconfig(ifile),readcomplex(ifile))  
        close(readunit(ifile))
     enddo
  endif
     
  call mympiibcast(readunit,1,numavectorfiles)
  call mympiibcast(readnumvects,1,numavectorfiles)
  call mympiibcast(readndof,1,numavectorfiles)
  call mympiibcast(readnumr,1,numavectorfiles)
  call mympiibcast(readnumconfig,1,numavectorfiles)
  call mympiibcast(readcomplex,1,numavectorfiles)

  do ifile=1,numavectorfiles
     if (readnumr(ifile).ne.numr) then
        OFLWR "error, for product all numr must agree", numr
        WRFL readnumr(1:numavectorfiles); CFLST
     endif
  enddo

  tot_ndof=0
  tot_numconfig=1
  tot_wfns=1
  do ifile=1,numavectorfiles
     tot_ndof=tot_ndof+readndof(ifile)
     tot_numconfig=tot_numconfig*readnumconfig(ifile)
     tot_wfns=tot_wfns*readnumvects(ifile)
  enddo

  if (tot_ndof.ne.ndof) then
     OFLWR "error, total number of electrons doesn't agree load_avector_product"
     WRFL tot_ndof,ndof,readndof(1:numavectorfiles); CFLST
  endif

  if (tot_wfns.ne.mcscfnum) then
     OFLWR "Load_avector_product: error: mcscfnum must equal total number of wave functions"; CFLST
  endif

  call mpibarrier()
  OFLWR "LOAD_AVECTOR_PRODUCT: tot_numconfig=",tot_numconfig; CFL
  OFLWR "     total number of wave functions = ",tot_wfns; CFL
  call mpibarrier()

  do ifile=1,numavectorfiles
     allocate(readavectors(ifile)%mat(numr,readnumconfig(ifile),readnumvects(ifile)))
     allocate(readconfiglist(ifile)%mat(readndof(ifile),readnumconfig(ifile)))
  enddo

  if (myrank.eq.1) then

     do ifile=1,numavectorfiles

        open(readunit(ifile),file=avectorfile(ifile), status="unknown", form="unformatted")
        call avector_header_read_simple(readunit(ifile),readnumvects(ifile),readndof(ifile),numr, &
             readnumconfig(ifile),readcomplex(ifile))  
        call simple_load_avectors(readunit(ifile),readcomplex(ifile), readavectors(ifile)%mat(:,:,:), &
             readndof(ifile), numr, readnumconfig(ifile), readnumvects(ifile))
        close(readunit(ifile))
        open(readunit(ifile),file=avectorfile(ifile), status="unknown", form="unformatted")

        call avector_header_read_simple(readunit(ifile),readnumvects(ifile),readndof(ifile),numr, &
             readnumconfig(ifile),readcomplex(ifile))  
        call get_avectorfile_configlist(readunit(ifile),readcomplex(ifile), readconfiglist(ifile)%mat(:,:), &
             readndof(ifile), numr, readnumconfig(ifile))
        close(readunit(ifile))
     enddo
  endif

  call mpibarrier()
  OFLWR "     ...A-Vectors read for product. broadcasting."; CFL
  call mpibarrier()

  do ifile=1,numavectorfiles
     call mympibcast(readavectors(ifile)%mat(:,:,:),1,readnumconfig(ifile)*numr*readnumvects(ifile))
     call mympiibcast(readconfiglist(ifile)%mat(:,:),1,readndof(ifile)*readnumconfig(ifile))
!     do iwfn=1,readnumvects(ifile)
!        OFLWR "    Vector norm on read",ifile,iwfn,dot(readavectors(ifile)%mat(:,:,iwfn), &
!             readavectors(ifile)%mat(:,:,iwfn), readnumconfig(ifile)*numr); CFL
!     enddo
  enddo

  if (numavectorfiles.gt.6) then
     OFLWR "REDIM LOAD_AVECTOR_PRODUCT",6,numavectorfiles; CFLST
  endif

  OFLWR "    ..go get newconfiglist load_avector_product"; CFL

  allocate(newconfiglist(ndof,tot_numconfig))

  iitop(:)=1
  iitop(1:numavectorfiles)=readnumconfig(1:numavectorfiles)

  ii1=>ii(1); ii2=>ii(2); ii3=>ii(3); ii4=>ii(4); ii5=>ii(5); ii6=>ii(6)

  iconfig=0

!! BUGFIX 12-2015 v1.16
  do ii6=1,iitop(6)
  do ii5=1,iitop(5)
  do ii4=1,iitop(4)
  do ii3=1,iitop(3)
  do ii2=1,iitop(2)
  do ii1=1,iitop(1)

     jj=0
     dofsum=0
     do ifile=1,numavectorfiles
        thisconfig(jj+1:jj+readndof(ifile))=readconfiglist(ifile)%mat(:,ii(ifile))
        do kk=1,readndof(ifile)-1,2
           thisconfig(jj+kk)=thisconfig(jj+kk)+    dofsum
        enddo
        dofsum=dofsum+eachloaded(ifile)
        jj=jj+readndof(ifile)
     enddo

     iconfig=iconfig+1
     if (iconfig.gt.tot_numconfig) then
        OFLWR "error iconfig (programmer fail)", iconfig,tot_numconfig,ii(:); CFLST
     endif
     newconfiglist(:,iconfig)=thisconfig(:)
!!$     iiarray(:,iconfig)=ii(:)

  enddo
  enddo
  enddo
  enddo
  enddo
  enddo

  if (iconfig.ne.tot_numconfig) then
     OFLWR "error iconfig  at end (programmer fail)", iconfig,tot_numconfig; CFLST
  endif

  OFLWR "    ..done get newconfiglist load_avector_product. go configtable"; CFL

  configtable(:)=(-99)

  num_allowed=0
  do iconfig=1,tot_numconfig
     dirphase=reorder(newconfiglist(:,iconfig),numelec)
     if (allowedconfig0(www,newconfiglist(:,iconfig),www%dflevel)) then
        num_allowed=num_allowed+1
        jconfig=getconfiguration(newconfiglist(:,iconfig),www)
        if (configtable(jconfig).ne.(-99)) then
           OFLWR "PROGFAIL CONFIGTABEEE"; CFLST
        endif
        configtable(jconfig)=iconfig
        configphase(jconfig)=dirphase
!     else
!        OFLWR "WARN NOT ALLOWED PRODUCT CONFIG"
!        call printconfig(newconfiglist(:,iconfig),www)
!        CFL
     endif
  enddo

  deallocate(newconfiglist)

  OFLWR "    ..done product table.  Go product vector"; CFL

  allocate(productvector(numr,iitop(1),iitop(2),iitop(3),iitop(4),iitop(5),iitop(6)))
  allocate(productreshape(numr,tot_numconfig))

  iwfn=0

  iitop(:)=1
  iitop(1:numavectorfiles)=readnumvects(1:numavectorfiles)

  do ii1=1,iitop(1)
  do ii2=1,iitop(2)
  do ii3=1,iitop(3)
  do ii4=1,iitop(4)
  do ii5=1,iitop(5)
  do ii6=1,iitop(6)

     iwfn=iwfn+1

     OFLWR "          ...construct wfn ",iwfn, " of ", mcscfnum; CFL

     if (iwfn.gt.mcscfnum) then
        OFLWR "PROGFAIL IWFN", iwfn,mcscfnum,ii(1:numavectorfiles);CFLST
     endif

     productvector=1d0
     do ir=1,numr
        do iconfig=1,readnumconfig(1)
           productvector(ir,iconfig,:,:,:,:,:)=productvector(ir,iconfig,:,:,:,:,:)*readavectors(1)%mat(ir,iconfig,ii(1))
        enddo
     enddo
     if (numavectorfiles.ge.2) then
        do ir=1,numr
           do iconfig=1,readnumconfig(2)
              productvector(ir,:,iconfig,:,:,:,:)=productvector(ir,:,iconfig,:,:,:,:)*readavectors(2)%mat(ir,iconfig,ii(2))
           enddo
        enddo
     endif
     if (numavectorfiles.ge.3) then
        do ir=1,numr
           do iconfig=1,readnumconfig(3)
              productvector(ir,:,:,iconfig,:,:,:)=productvector(ir,:,:,iconfig,:,:,:)*readavectors(3)%mat(ir,iconfig,ii(3))
           enddo
        enddo
     endif
     if (numavectorfiles.ge.4) then
        do ir=1,numr
           do iconfig=1,readnumconfig(4)
              productvector(ir,:,:,:,iconfig,:,:)=productvector(ir,:,:,:,iconfig,:,:)*readavectors(4)%mat(ir,iconfig,ii(4))
           enddo
        enddo
     endif
     if (numavectorfiles.ge.5) then
        do ir=1,numr
           do iconfig=1,readnumconfig(5)
              productvector(ir,:,:,:,:,iconfig,:)=productvector(ir,:,:,:,:,iconfig,:)*readavectors(5)%mat(ir,iconfig,ii(5))
           enddo
        enddo
     endif
     if (numavectorfiles.ge.6) then
        do ir=1,numr
           do iconfig=1,readnumconfig(6)
              productvector(ir,:,:,:,:,:,iconfig)=productvector(ir,:,:,:,:,:,iconfig)*readavectors(6)%mat(ir,iconfig,ii(6))
           enddo
        enddo
     endif
     
     productreshape(:,:)=RESHAPE(productvector,(/numr,tot_numconfig/))

!     OFLWR "check norm productreshape"
!     WRFL dot(productreshape,productreshape,numr*tot_numconfig)
!     CFL

     myavector(:,:,iwfn)=0d0

     do jconfig=first_config,last_config
        if (configtable(jconfig).gt.0) then
           myavector(:,jconfig,iwfn)=productreshape(:,configtable(jconfig))*configphase(jconfig)
        endif
     enddo

!     OFLWR "check norm myavector product"
!     WRFL dot(myavector,myavector,num_config*numr)
!     CFL

  enddo
  enddo
  enddo
  enddo
  enddo
  enddo

  do ifile=1,numavectorfiles
     deallocate(readavectors(ifile)%mat,readconfiglist(ifile)%mat)
  enddo

  deallocate(productvector,productreshape)

  OFLWR "   ... done load_avector_product."; CFL

end subroutine load_avector_productsub




subroutine get_avectorfile_configlist(iunit, qq, myconfiglist, myndof, mynumr, mynumconfig)
  use fileptrmod
  implicit none
  integer :: myndof, mynumconfig, mynumr,iunit,qq, config1, myiostat
  integer, intent(out) :: myconfiglist(myndof,mynumconfig)
  real*8 :: rtempreadvect(mynumr)
  complex*16 :: ctempreadvect(mynumr)

  do config1=1,mynumconfig
     if (qq==0) then
        read (iunit,iostat=myiostat) myconfiglist(:,config1), rtempreadvect(1:mynumr)
     else
        read (iunit,iostat=myiostat) myconfiglist(:,config1), ctempreadvect(1:mynumr)
     endif
     if (myiostat.ne.0) then
        OFLWR "err get_avectorfile_configlist "; CFLST
     endif
  enddo
end subroutine get_avectorfile_configlist





subroutine load_avectors(filename,myavectors,mynumvects,readnumvects,numskip)
  use parameters
  use mpimod
  implicit none
  character :: filename*(*)
  integer :: readnumvects,readndof,readnumr,readnumconfig,readcomplex,mynumvects,numskip,ii
  DATATYPE :: myavectors(numr,first_config:last_config,mynumvects)
  external :: readavectorsubroutine,readavectorsubsimple
  DATATYPE, allocatable :: readavectors(:,:,:)

  open(999,file=filename, status="unknown", form="unformatted")

 call avector_header_read_simple(999,readnumvects,readndof,readnumr,readnumconfig,readcomplex)  

 if (myrank.eq.1) then
    allocate(readavectors(numr,num_config,readnumvects))
 else
    allocate(readavectors(1,1,readnumvects))
 endif

  readnumvects=min(mynumvects,readnumvects-numskip)


  close(999)

  OFL
  if (numr>readnumr) then
     write(mpifileptr, *) "WARNING, smaller numr.", numr, readnumr
  endif
  if (num_config>readnumconfig) then 
    write(mpifileptr, *) "WARNING: smaller numconfigs.", num_config,readnumconfig
  endif    
  if (numr<readnumr) then
     write(mpifileptr, *) "WARNING, bigger numr.", readnumr, numr
  endif
  if (num_config<readnumconfig) then 
    write(mpifileptr, *) "WARNING: bigger numconfigs.", readnumconfig, num_config
  endif    


  if (myrank.eq.1) then
     open(999,file=filename, status="unknown", form="unformatted")

     call avector_header_read_simple(999,readnumvects,readndof,readnumr,readnumconfig,readcomplex)

     if (numholes.eq.0.and.excitations.eq.0) then

        call easy_load_avectors(999,readcomplex, readavectors(:,:,:), readnumr, readnumconfig, readnumvects)

     else
        call load_avectors0(999,readcomplex,readavectors(:,:,:),numr,num_config,ndof+2*numholes,readnumr,readnumconfig, readavectorsubroutine,readnumvects)
     endif

     close(999)

     readnumvects=min(mynumvects,readnumvects-numskip)

  endif

  if (par_consplit.eq.0) then
     if (myrank.eq.1) then
        myavectors(:,:,1:readnumvects) = readavectors(:,:, 1+numskip : numskip+readnumvects)
     endif
     call mympibcast(myavectors(:,:,:),1,num_config*numr*readnumvects)
  else
     do ii=1,readnumvects
        call myscatterv(readavectors(:,:,ii+numskip),myavectors(:,:,ii),configs_perproc(:)*numr)
     enddo
  endif

  deallocate(readavectors)

end subroutine load_avectors


!! general subroutine, for loading vectors with some manipulating on read
!!     used in load_avectors (main load routine) and ovlsub.

subroutine load_avectors0(iunit, qq, myavectors, mynumr, mynumconfig, readndof, readnumr, readnumconfig, mysubroutine, mynumvects )
  use fileptrmod
  implicit none
 
  external :: mysubroutine
  integer :: mynumconfig, mynumr, mynumvects,iunit,i,  readndof, readnumr, readnumconfig
  integer :: qq, config1,  thatconfig(readndof),  myiostat, ivect
  DATATYPE :: myavectors(mynumr,mynumconfig,mynumvects),  mytempavector(mynumconfig), readvect(readnumr)
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
           myavectors(i,:,ivect)= myavectors(i,:,ivect) +  readvect(i) * mytempavector(:)
        enddo
     enddo
  enddo
end subroutine load_avectors0
  

subroutine simple_load_avectors(iunit, qq, myavectors, myndof, mynumr, mynumconfig, mynumvects)
  use fileptrmod
  implicit none

  integer :: myndof, mynumconfig, mynumr,mynumvects,iunit,ivect,qq, config1, thatconfig(myndof), myiostat
  DATATYPE :: myavectors(mynumr,mynumconfig,mynumvects)
  real*8 :: rtempreadvect(mynumr)
  complex*16 :: ctempreadvect(mynumr)

  myavectors=0d0

  do ivect=1,mynumvects
     do config1=1,mynumconfig
        if (qq==0) then
           read (iunit,iostat=myiostat) thatconfig(1:myndof), rtempreadvect(1:mynumr)
           myavectors(:,config1,ivect)=rtempreadvect(:)
        else
           read (iunit,iostat=myiostat) thatconfig(1:myndof), ctempreadvect(1:mynumr)
           myavectors(:,config1,ivect)=ctempreadvect(:)
        endif
        if (myiostat.ne.0) then
           OFLWR "err read config "; CFLST
        endif
     enddo
  enddo
end subroutine simple_load_avectors


!! outavectors = numconfig   RANK 1 calls this subroutine, THEN DO PARCONSPLIT LOGIC

subroutine easy_load_avectors(iunit, qq, outavectors, mynumr, mynumconfig, mynumvects)
  use parameters
  use configmod
  use mpimod
  implicit none

  integer ::  mynumconfig, mynumr,mynumvects,iunit,ivect,qq, config1, thatconfig(ndof), myiostat,&
       phase,reorder,myconfig,getconfiguration
  DATATYPE :: outavectors(mynumr,num_config,mynumvects)
  real*8 :: rtempreadvect(mynumr)
  logical :: allowedconfig0
  complex*16 :: ctempreadvect(mynumr)

  if (mynumr.gt.numr) then
     OFLWR "error numr on file greater than calc",mynumr,numr; CFLST
  endif

  if (myrank.ne.1) then
     OFLWR "only call me rank 1 easy"; CFLST
  endif

  outavectors=0d0

  do ivect=1,mynumvects
     do myconfig=1,mynumconfig
        if (qq==0) then
           read (iunit,iostat=myiostat) thatconfig(1:ndof), rtempreadvect(1:mynumr)
        else
           read (iunit,iostat=myiostat) thatconfig(1:ndof), ctempreadvect(1:mynumr)
        endif
        phase=reorder(thatconfig,numelec)
        if (allowedconfig0(www,thatconfig,www%dflevel)) then
           config1=getconfiguration(thatconfig,www)
           if (qq==0) then
              outavectors(1:mynumr,config1,ivect)=rtempreadvect(:)*phase
           else
              outavectors(1:mynumr,config1,ivect)=ctempreadvect(:)*phase
           endif
        else
           if (readfullvector) then
              OFLWR "Config on file not allowed!!!"; CFLST
           endif
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
  use configmod
  use aarrmod
  implicit none
  integer :: i,iihole,iloop,jloop,ivect
  DATATYPE :: outavector(num_config)
  integer ::  readconfig(ndof+2*numholes),xflag,iexcite,config2
  integer :: thisconfig(ndof+400)
  integer :: jj,iind,phase,reorder
  integer :: getconfiguration
  logical :: allowedconfig0
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
                 thisconfig(2*jj-1:2*jj)=aarr(abs(myavectorexciteto(iexcite,jloop,ivect)))
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
        phase=reorder(thisconfig,numelec)
        if (allowedconfig0(www,thisconfig,www%dflevel)) then
           config2=getconfiguration(thisconfig,www)
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

subroutine readavectorsubsimple(readconfig, outavector,notusedint)
  use parameters
  use configmod
  use aarrmod
  implicit none
  integer :: i,notusedint,phase,reorder,readconfig(ndof),config2,thisconfig(ndof),getconfiguration
  logical :: allowedconfig0
  DATATYPE :: outavector(num_config)

  outavector(:)=0d0;  thisconfig(:)=readconfig(:)
  do i=1,numelec
     thisconfig(i*2-1)   = thisconfig(i*2-1)-numloadfrozen   
  enddo
  phase=reorder(thisconfig,numelec)
  
  if (allowedconfig0(www,thisconfig,www%dflevel)) then
     config2=getconfiguration(thisconfig,www)
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
  use mpimod
  implicit none
  integer :: unit, config1
  DATATYPE :: avector(numr,num_config)
  if (myrank.ne.1) then
     OFLWR "only call write_avector on process 1 right????"; CFLST
  endif
  do config1=1,num_config
     write (unit) www%configlist(:,config1), avector(:,config1)
  enddo
end subroutine write_avector

