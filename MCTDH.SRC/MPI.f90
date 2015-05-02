
!! PARALLEL SUBROUTINES.

#include "Definitions.INC"

subroutine getmyranknprocs(outrank,outnprocs)
  use mpimod
  implicit none
  integer :: outrank,outnprocs
  outrank=myrank;  outnprocs=nprocs
end subroutine getmyranknprocs


subroutine getWorldCommGroup(outcommunicator,outgroup)
  use mpimod
  implicit none
  integer, intent(out) :: outcommunicator,outgroup
#ifndef MPIFLAG
  outcommunicator= (-798)
  outgroup= 42
#else
  outcommunicator=MPI_COMM_WORLD
  outgroup=MPI_GROUP_WORLD
#endif
end subroutine getWorldCommGroup


subroutine mpiorbsets()
  use mpimod
  use mpi_orbsetmod
  use parameters
  implicit none
  integer :: iproc,iset,jproc,ierr
  integer, allocatable :: procsperset(:)
#ifdef MPIFLAG
  integer, allocatable :: process_ranks(:,:)
#endif
  if (nprocs.lt.nspf) then
     OFLWR " Get ORbsets.  Fewer procs than orbitals."; CFL
     norbsets=1
     orbsperproc=ceiling(real(nspf)/nprocs)
     allocate(procsperset(norbsets))
     procsperset(:)=nprocs
  else   
     norbsets=floor(real(nprocs)/nspf)
     orbsperproc=1
     allocate(procsperset(norbsets))     
     procsperset(:)=nspf
     procsperset(norbsets)= nprocs-nspf*(norbsets-1)
     iproc=0
     do iset=1,norbsets
        iproc=iproc+procsperset(iset)
     enddo
     if (iproc.ne.nprocs) then
        OFLWR "TTTYTUTU ERROR GUY."; CFLST
     endif
  endif
#ifndef MPIFLAG
  myorbset=1
  firstmpiorb=1
  jproc=0  ;   ierr=0     !! avoid warn unused
#else
  allocate(process_ranks(nspf*2,norbsets)); process_ranks(:,:)=(-1)
  allocate(MPI_GROUP_ORB(norbsets), MPI_COMM_ORB(norbsets)); MPI_GROUP_ORB(:)=0; MPI_COMM_ORB(:)=0
  iproc=0
  do iset=1,norbsets
     jproc=0
     do while (jproc.lt.procsperset(iset))
        jproc=jproc+1
        iproc=iproc+1
        process_ranks(jproc,iset)=iproc
        if (iproc.eq.myrank) then
           myorbset=iset
           firstmpiorb=(jproc-1)*orbsperproc+1
        endif
     enddo
     process_ranks(:,iset)=process_ranks(:,iset)-1
     call MPI_group_incl(MPI_GROUP_WORLD,procsperset(iset),process_ranks(:,iset),MPI_GROUP_ORB(iset),ierr)
     if (ierr.ne.0) then
        OFLWR "group_incl ERR ", ierr; CFLST
     endif
     call MPI_comm_create(MPI_COMM_WORLD, MPI_GROUP_ORB(iset),MPI_COMM_ORB(iset),ierr)
     if (ierr.ne.0) then
        OFLWR "group_create err ", ierr; CFLST
     endif
  enddo
  if (myorbset.le.0) then
     OFLWR "SSSETERRROR"; CFLST
  endif
  deallocate(process_ranks,procsperset)
#endif

  mpi_orbset_init=1
  mpi_orbset_init=1
  mpi_orbset_init=1

end subroutine mpiorbsets


#ifdef MPIFLAG

subroutine mpiorbgather(orbvector,insize)    !! insize=spfsize except debug
  use mpimod
  use mpi_orbsetmod
  use fileptrmod
  use parameters
  implicit none
  integer :: ierr,insize
  DATATYPE :: orbvector(insize,nspf*2)
  if (mpi_orbset_init.ne.1) then
     OFLWR "Programmer fail, mpiorbgather called but mpiorbsets() appears "
     WRFL "   not to have been called."; CFLST
  endif
  if (nprocs.eq.1) then
     return
  endif
  if (firstmpiorb+orbsperproc-1.gt.nspf*2 .or. orbsperproc*min(nspf,nprocs).gt.nspf*2) then
     OFLWR "YYY ERROR",firstmpiorb,orbsperproc; CFLST
  endif
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
#ifdef REALGO
  call mpi_allgather(orbvector(:,firstmpiorb:firstmpiorb+orbsperproc-1),&
       insize*orbsperproc,MPI_DOUBLE_PRECISION,orbvector(:,:),insize*orbsperproc,MPI_DOUBLE_PRECISION,MPI_COMM_ORB(myorbset),ierr)
#else
  call mpi_allgather(orbvector(:,firstmpiorb:firstmpiorb+orbsperproc-1),&
       insize*orbsperproc,MPI_DOUBLE_COMPLEX,orbvector(:,:),insize*orbsperproc,MPI_DOUBLE_COMPLEX,MPI_COMM_ORB(myorbset),ierr)
#endif
  if (ierr.ne.0) then
     OFLWR "ORBGATHER ERR ", ierr; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mpiorbgather


subroutine mpistart()
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr,  iilen !!,provided,required
  character(len=40) :: format
  integer :: nargs, i, flag
#ifdef PGFFLAG
  integer :: myiargc
#endif
  character (len=200) :: buffer
#ifdef PGFFLAG
  nargs=myiargc()
#else
  nargs=iargc()
#endif
  flag=0
  do i=1,nargs
     call getarg(i,buffer); 
     if (buffer(1:6) .eq. 'MPIOUT') then
        flag=1
     endif
  enddo
  call MPI_INIT(ierr)
  if (ierr/=0) then; print *,  "MPI ERR 1";  stop; 
  endif
  call MPI_Comm_group(MPI_COMM_WORLD,MPI_GROUP_WORLD,ierr)
  if (ierr.ne.0) then
     OFLWR "world comm group ERR ", ierr; CFLST
  endif
  call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr)   
  if (ierr/=0) then;   print *,  "MPI ERR 2";   stop; 
  endif

  myrank=myrank+1

  call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr)     
  if (ierr/=0) then;   print *,   "MPI ERR 3";   stop; 
  endif
  if (nprocs.eq.1.or.myrank.eq.1) then
     stdoutflag=1
  else
     stdoutflag=0
  endif
  if (stdoutflag==1) then
     mpifileptr=6
  else
     mpifileptr=987
     call system("mkdir MPIOUTS")
     if (mpioutfilelen.gt.9) then
        write(format,'(A2,I2,A2,I1,A1)') "(A",mpioutfilelen,",I",iilen(myrank),")"
     else
        write(format,'(A2,I1,A2,I1,A1)') "(A",mpioutfilelen,",I",iilen(myrank),")"
     endif
     write(mpioutfile,format) mpioutfilebase,myrank

!! nevermind writing out every one
     if (flag.eq.0) then
        if (myrank.ne.1) then
           mpioutfile="/dev/null"
        endif
     endif
     open(mpifileptr,file=mpioutfile,status="unknown"); WRFL; CFL
  endif
  call system_clock(mpibtime);  call system_clock(mpiatime)

end subroutine mpistart


subroutine mpibarrier()
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr
  call system_clock(mpiatime)
  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if (ierr/=0) then
     OFLWR "MPI ERR 3"; CFLST
  endif
  call system_clock(mpibtime)
  mpitime=mpitime+mpibtime-mpiatime
end subroutine mpibarrier


subroutine mpistop()
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr
  call mpi_finalize(ierr)
  if (ierr/=0) then
     OFLWR "MPI ERR 3"; CFL
  endif
  OFLWR "MCTDHF STOP!";CFL
  stop
end subroutine mpistop


subroutine mympireduce(input, isize)
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr, isize
  DATATYPE,intent(inout) :: input(isize)
  DATATYPE :: output(isize)           !! AUTOMATIC
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  ierr=0
  call MPI_allreduce( input, output, isize, MPIDATATYPE, MPI_SUM, MPI_COMM_WORLD , ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR mympireduce!";   CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympireduce


subroutine mympireduceto(input, output, isize, dest)
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr, isize ,dest,idest
  DATATYPE,intent(in) :: input(isize)
  DATATYPE,intent(out) :: output(isize) 
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  if (dest.lt.1.or.dest.gt.nprocs) then
     OFLWR "ERR DEST REDUCETO",dest,nprocs; CFLST
  endif
  ierr=0
  idest=dest-1
  call MPI_reduce( input, output, isize, MPIDATATYPE, MPI_SUM, idest, MPI_COMM_WORLD , ierr)
  if (ierr/=0) then
     OFLWR "ERR mpi_reduce!";   CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympireduceto


!! used in dGMRES parallel routines

subroutine mympirealreduceone(input)
  use mpimod
  use fileptrmod
  implicit none
  real*8 :: input,output
  integer :: ierr,isize=1
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  ierr=0
  call MPI_allreduce( input, output, isize, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD , ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR mympirealreduce!"; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympirealreduceone


subroutine mympiireduceone(input)
  use mpimod
  use fileptrmod
  implicit none
  integer :: input,output,ierr,isize=1
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  ierr=0
  call MPI_allreduce( input, output, isize, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD , ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR mympiireduce!"; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympiireduceone


subroutine mympiireduce(input, isize)
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr, isize
  integer :: input(isize), output(isize) 
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  ierr=0
  call MPI_allreduce( input, output, isize, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD , ierr)
  if (ierr/=0) then
     OFLWR "ERR mympiireduce!";   CFLST
  endif
  input=output
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympiireduce


subroutine mympireduceone(input)
  use mpimod
  use fileptrmod
  implicit none
  DATATYPE :: input,output
  integer :: ierr 
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  ierr=0
  call MPI_allreduce( input, output, 1, MPIDATATYPE, MPI_SUM, MPI_COMM_WORLD , ierr)
  if (ierr/=0) then
     OFLWR "ERR mympireduce!"; CFLST
  endif
  input=output
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympireduceone


subroutine mympibcast(input, source, isize)
  use mpimod
  implicit none
  integer,intent(in) :: isize, source
  DATATYPE,intent(inout) :: input(isize)
#ifdef REALGO
  call mympirealbcast_local(input, source, isize,MPI_COMM_WORLD)
#else
  call mympicomplexbcast_local(input, source, isize,MPI_COMM_WORLD)
#endif

end subroutine mympibcast


subroutine mympisendrecv(sendbuf, recvbuf, dest, source, tag, isize)
  use mpimod
  use fileptrmod
  implicit none
  integer, intent(in) :: dest,source,tag,isize
  integer :: ierr, idest,isource
  DATATYPE, intent(in) :: sendbuf(isize)
  DATATYPE, intent(out) :: recvbuf(isize)

  idest=dest-1
  isource=source-1
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_sendrecv(sendbuf,isize,MPIDATATYPE,idest,tag,&
       recvbuf,isize, MPIDATATYPE,isource,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympisendrecv",ierr; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympisendrecv


subroutine mympisendrecv_complex_local(sendbuf, recvbuf, dest, source, tag, isize,MPI_COMM_LOCAL)
  use mpimod
  use fileptrmod
  implicit none
  integer, intent(in) :: dest,source,tag,isize,MPI_COMM_LOCAL
  integer :: ierr, idest,isource
  complex*16, intent(in) :: sendbuf(isize)
  complex*16, intent(out) :: recvbuf(isize)
  idest=dest-1
  isource=source-1
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_sendrecv(sendbuf,isize,MPI_DOUBLE_COMPLEX,idest,tag,&
       recvbuf,isize, MPI_DOUBLE_COMPLEX,isource,tag,MPI_COMM_LOCAL,MPI_STATUS_IGNORE,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympisendrecv",ierr; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympisendrecv_complex_local



subroutine mympisend(input, dest, tag, isize)
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr, isize, dest, idest,tag
  DATATYPE :: input(isize)
  idest=dest-1
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_send(input,isize,MPIDATATYPE,idest,tag,MPI_COMM_WORLD,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympisend"; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympisend



subroutine mympirecv(input, source, tag, isize)
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr, isize, source, isource,tag
  DATATYPE :: input(isize)
  isource=source-1
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_recv(input,isize,MPIDATATYPE,isource,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympirecv"; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympirecv


subroutine mympirealbcast_local(input, source, isize,MPI_COMM_LOCAL)
  use mpimod
  use fileptrmod
  implicit none
  integer, intent(in) ::  isize, source, MPI_COMM_LOCAL
  real*8, intent(inout) :: input(isize)
  integer :: ierr, isource

  isource=source-1
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_bcast(input,isize,MPI_DOUBLE_PRECISION,isource,MPI_COMM_WORLD,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympirealbcast"; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympirealbcast_local


subroutine mympicomplexbcast_local(input, source, isize, MPI_COMM_LOCAL)
  use mpimod
  use fileptrmod
  implicit none
  integer, intent(in) ::  isize, source, MPI_COMM_LOCAL
  complex*16, intent(inout) :: input(isize)
  integer :: ierr, isource

  isource=source-1
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_bcast(input,isize,MPI_DOUBLE_COMPLEX,isource,MPI_COMM_LOCAL,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympicomplexbcast"; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympicomplexbcast_local


subroutine mympiibcast(input, source, isize)
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr, isize, source, isource, input(isize)
  isource=source-1
  call system_clock(mpiatime); nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_bcast(input,isize,MPI_INTEGER,isource,MPI_COMM_WORLD,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympibcast"; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympiibcast


subroutine mympirealbcastone(input, source)
  use mpimod
  use fileptrmod
  implicit none
  real*8 :: input
  integer :: ierr, source, isource
  isource=source-1
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_bcast(input,1,MPI_DOUBLE_PRECISION,isource,MPI_COMM_WORLD,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympibcast"; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympirealbcastone


subroutine mympiibcastone(input, source)
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr,  source, isource, input
  isource=source-1
  call system_clock(mpiatime); nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_bcast(input,1,MPI_INTEGER,isource,MPI_COMM_WORLD,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympibcast"; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympiibcastone


subroutine mympilogbcast(input, source, isize)
  use mpimod
  use fileptrmod
  implicit none
  logical :: input(isize)
  integer :: ierr, isize, source, isource
  isource=source-1
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_bcast(input,isize,MPI_LOGICAL,isource,MPI_COMM_WORLD,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympibcast"; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympilogbcast


subroutine mympiimax(input)
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr,input, output
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call MPI_reduce( input, output, 1, MPI_INTEGER, MPI_MAX,0, MPI_COMM_WORLD , ierr)
  if (ierr/=0) then
     OFLWR "ERR MPIMAX!"; CFLST
  endif
  call mpi_bcast(output,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR MPIIMAX!"; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympiimax


subroutine mympiimin(input)
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr,input, output
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call MPI_reduce( input, output, 1, MPI_INTEGER, MPI_MIN,0, MPI_COMM_WORLD , ierr)
  if (ierr/=0) then
     OFLWR "ERR MPIImin!"; CFLST
  endif
  call mpi_bcast(output,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR MPIImin!"; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympiimin


subroutine mympimax(input)
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr
  real*8 :: input, output
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call MPI_reduce( input, output, 1, MPI_DOUBLE_PRECISION, MPI_MAX,0, MPI_COMM_WORLD , ierr)
  if (ierr/=0) then
     OFLWR "ERR MPIMAX!"; CFLST
  endif
  call mpi_bcast(output,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR MPIIMAX!"; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympimax


subroutine mympimin(input)
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr
  real*8 :: input, output
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call MPI_reduce( input, output, 1, MPI_DOUBLE_PRECISION, MPI_MIN,0, MPI_COMM_WORLD , ierr)
  if (ierr/=0) then
     OFLWR "ERR MPImin!"; CFLST
  endif
  call mpi_bcast(output,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR MPIImin!"; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympimin




subroutine simpleallgather(vectorin,vectorsout,insize)
  use mpimod
  implicit none
  integer :: insize
  DATATYPE :: vectorin(insize),vectorsout(insize,nprocs)
#ifdef REALGO
  call simpleallgather_real(vectorin,vectorsout,insize)
#else
  call simpleallgather_complex(vectorin,vectorsout,insize)
#endif
end subroutine simpleallgather



subroutine simpleallgather_complex(vectorin,vectorsout,insize)
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr,insize
  complex*16 :: vectorin(insize),vectorsout(insize,nprocs)
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_allgather(vectorin(:),insize,MPI_DOUBLE_COMPLEX,vectorsout(:,:),insize,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "gather complex ERR ", ierr; call mpistop()
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine simpleallgather_complex


subroutine simpleallgather_real(vectorin,vectorsout,insize)
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr,insize
  real*8 :: vectorin(insize),vectorsout(insize,nprocs)
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_allgather(vectorin(:),insize,MPI_DOUBLE_PRECISION,vectorsout(:,:),insize,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "gather real ERR ", ierr; call mpistop()
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine simpleallgather_real


subroutine getgatherv_stuff(blocks,parsize,blockstart)
  use mpimod
  integer, intent(in) :: blocks(nprocs)
  integer, intent(out) :: parsize, blockstart(nprocs)
  integer :: ii, isum
  isum=blocks(1)
  blockstart(1)=1  
  do ii=2,nprocs
     blockstart(ii)=blockstart(ii-1)+blocks(ii-1)
     isum=isum+blocks(ii)
  enddo
  parsize=isum
end subroutine getgatherv_stuff


SUBROUTINE MYGATHERV_complex(V1,X1,BLOCKS,bcastflag)
  use mpimod
  IMPLICIT NONE
  logical, intent(in) :: bcastflag
  INTEGER, INTENT(IN) :: BLOCKS(nprocs)
  complex*16, INTENT(IN) :: V1(blocks(myrank))
  complex*16, INTENT(OUT) :: X1(*)
  INTEGER ::MPIERR=0, BLOCKSTART(nprocs),parsize

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  call getgatherv_stuff(blocks,parsize,blockstart)

  CALL MPI_GATHERV(V1, blocks(myrank), MPI_DOUBLE_COMPLEX, X1, BLOCKS, &
       BLOCKSTART(:)-1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, MPIERR)
  if (mpierr.ne.0) then
     print *, "GAtherv err",mpierr,myrank;          stop
  endif

  if (bcastflag) then
     CALL MPI_BCAST(X1,PARSIZE,MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, MPIERR)
     if (mpierr.ne.0) then
        print *, "bcast err",mpierr,myrank;          stop
     endif
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

END SUBROUTINE MYGATHERV_complex


SUBROUTINE MYGATHERV_real(V1,X1,BLOCKS,bcastflag)
  use mpimod
  IMPLICIT NONE
  logical, intent(in) :: bcastflag
  INTEGER, INTENT(IN) :: BLOCKS(nprocs)
  real*8, INTENT(IN) :: V1(blocks(myrank))
  real*8, INTENT(OUT) :: X1(*)
  INTEGER ::MPIERR=0, parsize, BLOCKSTART(nprocs)

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  call getgatherv_stuff(blocks,parsize,blockstart)

  CALL MPI_GATHERV(V1, blocks(myrank), MPI_DOUBLE_PRECISION, X1, BLOCKS, &
       BLOCKSTART(:)-1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPIERR)
  if (mpierr.ne.0) then
     print *, "GAtherv err",mpierr,myrank;          stop
  endif
  if (bcastflag) then
     CALL MPI_BCAST(X1,PARSIZE,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPIERR)
     if (mpierr.ne.0) then
        print *, "bcast err",mpierr,myrank;          stop
     endif
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

END SUBROUTINE MYGATHERV_real




SUBROUTINE MYSCATTERV_complex(X1,V1,BLOCKS)
  use mpimod
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: BLOCKS(nprocs)
  complex*16, INTENT(OUT) :: V1(blocks(myrank))
  complex*16, INTENT(IN) :: X1(*)
  INTEGER ::MPIERR=0, parsize, BLOCKSTART(nprocs)

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  call getgatherv_stuff(blocks,parsize,blockstart)

  CALL MPI_SCATTERV(X1, blocks, blockstart(:)-1, MPI_DOUBLE_COMPLEX, V1, BLOCKS(myrank), &
       MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, MPIERR)
  if (mpierr.ne.0) then
     print *, "Scatterv err",mpierr,myrank;          stop
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

END SUBROUTINE MYSCATTERV_complex


SUBROUTINE MYSCATTERV_real(X1,V1,BLOCKS)
  use mpimod
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: BLOCKS(nprocs)
  real*8, INTENT(OUT) :: V1(blocks(myrank))
  real*8, INTENT(IN) :: X1(*)
  INTEGER ::MPIERR=0, parsize, BLOCKSTART(nprocs)

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  call getgatherv_stuff(blocks,parsize,blockstart)

  CALL MPI_SCATTERV(X1, blocks, blockstart(:)-1, MPI_DOUBLE_PRECISION, V1, BLOCKS(myrank), &
       MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPIERR)
  if (mpierr.ne.0) then
     print *, "Scatterv err",mpierr,myrank;          stop
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

END SUBROUTINE MYSCATTERV_real


#else

subroutine mympiimax(input)
  implicit none
  integer :: input
  return
  input=input 
end subroutine mympiimax

subroutine mympiimin(input)
  implicit none
  integer :: input
  return
  input=input
end subroutine mympiimin

subroutine mympimax(input)
  implicit none
  real*8 :: input
  return
  input=input
end subroutine mympimax

subroutine mympimin(input)
  implicit none
  real*8 :: input
  return
  input=input
end subroutine mympimin



SUBROUTINE MYGATHERV_complex(V1,X1,BLOCKS,bcastflag)
  IMPLICIT NONE
  logical, intent(in) :: bcastflag
  INTEGER, INTENT(IN) :: BLOCKS(1)
  complex*16, INTENT(IN) :: V1(blocks(1))
  complex*16, INTENT(OUT) :: X1(blocks(1))
  logical :: ilog
  X1(:)=V1(:)
  return
  ilog=bcastflag   !!avoid warn unused
END SUBROUTINE MYGATHERV_complex


SUBROUTINE MYGATHERV_real(V1,X1,BLOCKS,bcastflag)
  use mpimod
  IMPLICIT NONE
  logical, intent(in) :: bcastflag
  INTEGER, INTENT(IN) :: BLOCKS(1)
  real*8, INTENT(IN) :: V1(blocks(1))
  real*8, INTENT(OUT) :: X1(blocks(1))
  logical :: ilog
  X1(:)=V1(:)
  return
  ilog=bcastflag   !!avoid warn unused
END SUBROUTINE MYGATHERV_real




SUBROUTINE MYSCATTERV_complex(X1,V1,BLOCKS)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: BLOCKS(1)
  complex*16, INTENT(OUT) :: V1(blocks(1))
  complex*16, INTENT(IN) :: X1(blocks(1))
  V1(:)=X1(:)
END SUBROUTINE MYSCATTERV_complex


SUBROUTINE MYSCATTERV_real(X1,V1,BLOCKS)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: BLOCKS(1)
  real*8, INTENT(OUT) :: V1(blocks(1))
  real*8, INTENT(IN) :: X1(blocks(1))
  V1(:)=X1(:)
END SUBROUTINE MYSCATTERV_real






subroutine mympirealreduceone(input)
  implicit none
  real*8 :: input
  return
  input=input
end subroutine

subroutine mympireduceone(input)
  implicit none
  DATATYPE :: input
  return
  input=input
end subroutine mympireduceone

subroutine mympireduce(input, isize)
  implicit none
  integer :: isize
  DATATYPE :: input(isize)
  return
  input(1)=input(1)
end subroutine mympireduce

subroutine mympiireduce(input, isize)
  implicit none
  integer :: isize, input(isize)
  return
  input(1)=input(1)
end subroutine mympiireduce

subroutine mympibcast(input, source, isize)
  implicit none
  integer :: source,isize
  DATATYPE :: input(isize)
  return
  source=source;  input(1)=input(1)
end subroutine mympibcast

subroutine mympirealbcast(input, source, isize)
  implicit none
  integer :: source,isize
  real*8 :: input(isize)
  return
  source=source;  input(1)=input(1)
end subroutine mympirealbcast

subroutine mympicomplexbcast(input, source, isize)
  implicit none
  integer :: source,isize
  complex*16 :: input(isize)
  return
  source=source;  input(1)=input(1)
end subroutine mympicomplexbcast

subroutine mympiibcast(input, source, isize)
  implicit none
  integer :: source,isize
  integer :: input(isize)
  return
  source=source;  input(1)=input(1)
end subroutine mympiibcast

subroutine mympirealbcastone(input, source)
  implicit none
  integer :: source
  real*8 :: input
  return
  source=source;  input=input
end subroutine mympirealbcastone

subroutine mympiibcastone(input, source)
  implicit none
  integer :: input,source
  return
  source=source; input=input
end subroutine mympiibcastone

subroutine mympilogbcast(input, source, isize)
  implicit none
  integer :: source,isize
  DATATYPE :: input(isize)
  return
  source=source; input(1)=input(1)
end subroutine mympilogbcast


subroutine mympiireduceone(input)
  implicit none
  integer :: input
  return
  input=input
end subroutine mympiireduceone


subroutine mpiorbgather(orbvector,insize)
  implicit none
  integer :: insize
  DATATYPE :: orbvector(insize)
  return
  orbvector(1)=orbvector(1)
end subroutine mpiorbgather


subroutine mympireduceto(input, output, isize, dest)
  use mpimod
  use fileptrmod
  implicit none
  integer :: isize ,dest
  DATATYPE :: input(isize), output(isize)
  if (dest.ne.1) then
     OFLWR "Error nonpar mympireduceto wrapper",dest; CFLST
  endif
  output(:)=input(:)
end subroutine mympireduceto


subroutine mympisendrecv(sendbuf, recvbuf, dest, source, tag, isize)
  use mpimod
  use fileptrmod
  implicit none
  integer :: isize, dest, source,tag
  DATATYPE :: sendbuf(isize),recvbuf(isize)
  if (dest.ne.1.or.source.ne.1) then
     OFLWR "Error nonpar mympireduceto wrapper",dest,source; CFLST
  endif
  recvbuf(:)=sendbuf(:)
  return
  tag=tag 
end subroutine mympisendrecv


subroutine mpistart
  use mpimod
  use fileptrmod
  implicit none
  integer ::  iilen
  character(len=40) :: format
  myrank=1
  nprocs=1
  stdoutflag=1
  if (stdoutflag==1) then
     mpifileptr=6
  else
     mpifileptr=987
     write(format,'(A2,I1,A2,I1,A1)') "(A",mpioutfilelen,",I",iilen(myrank),")"
     write(mpioutfile,format) mpioutfilebase,myrank
     open(mpifileptr,file=mpioutfile,status="unknown")
     write(mpifileptr,*);     close(mpifileptr)
  endif

end subroutine mpistart


subroutine mpiinit
end subroutine mpiinit


subroutine mpibarrier()
end subroutine mpibarrier


subroutine mpistop()
  stop
end subroutine mpistop


#endif


function iilen(which)
  use mpimod
  implicit none
  integer :: iilen, which
  iilen=0
  if (which < 10) then
     iilen=1
  else if (which < 100) then
     iilen=2
  else if (which < 1000) then
     iilen=3
  else if (which < 10000) then
     iilen=4
  else 
     write(*,*) "Whoa! reformat"
     call mpistop()
  endif
end function iilen


subroutine openfile()
  use mpimod
  use fileptrmod  
  implicit none
  if (stdoutflag==0) then
     open(mpifileptr,file=mpioutfile,status="old", position="append")
  endif
end subroutine openfile


subroutine closefile()
  use mpimod
  use fileptrmod
  implicit none
  if (stdoutflag==0)  then
     close(mpifileptr)
  endif
end subroutine closefile


subroutine mpiallgather(inout,totsize,blocksizes,notusedint)
  use mpimod
  use fileptrmod
  implicit none
  integer, intent(in) :: totsize,blocksizes(nprocs)
  DATATYPE :: inout(totsize)
  integer :: icount,notusedint
#ifndef MPIFLAG
  return
  inout(1)=inout(1); icount=blocksizes(1)
#else
  integer :: ierr,blockstart(nprocs)
  if (nprocs.eq.1) then
     return
  endif
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call getgatherv_stuff(blocksizes,icount,blockstart)
  if (icount.ne.totsize) then
     OFLWR "ALLgather ERR count", icount,totsize; 
     WRFL blocksizes;
     CFLST
  endif
  call mpi_allgatherv(inout(blockstart(myrank)),blocksizes(myrank),MPIDATATYPE,inout,blocksizes(:),blockstart(:)-1,MPIDATATYPE,MPI_COMM_WORLD,ierr)
  if (ierr.ne.0) then
     OFLWR "ALLGATHERv ERR ", ierr; CFLST
  endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
#endif
end subroutine mpiallgather



subroutine mympialltoall(input, output, count)
  use fileptrmod
  use mpimod
  implicit none
  integer :: ierr, count
  DATATYPE,intent(in) :: input(count,nprocs)
  DATATYPE, intent(out) :: output(count,nprocs)

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
#ifndef MPIFLAG
  output(:,:)=input(:,:)
  return
  ierr=ierr; count=count  !! avoid warn unused
#else
#ifdef REALGO
  call mpi_alltoall(input(:,:),count,MPI_DOUBLE_PRECISION,output(:,:),count,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
#else
  call mpi_alltoall(input(:,:),count,MPI_DOUBLE_COMPLEX,output(:,:),count,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
#endif
  if (ierr.ne.0) then
     OFLWR "ERROR ALLTOALL ", ierr; CFLST
  endif
#endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympialltoall



subroutine mympialltoall_complex(input, output, count)
  use fileptrmod
  use mpimod
  implicit none
  integer :: ierr, count
  complex*16, intent(in) :: input(count,nprocs)
  complex*16, intent(out) :: output(count,nprocs)

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
#ifndef MPIFLAG
  output(:,:)=input(:,:)
  return
  ierr=0;   !! avoid warn unused
#else
  call mpi_alltoall(input(:,:),count,MPI_DOUBLE_COMPLEX,output(:,:),count,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
  if (ierr.ne.0) then
     OFLWR "ERROR ALLTOALL ", ierr; CFLST
  endif
#endif
  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympialltoall_complex

