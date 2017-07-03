
!! PARALLEL SUBROUTINES.  NO DEPENDENCIES ON MAIN_MODULES !!  ONLY PARAMETERS.F90

#include "Definitions.INC"

subroutine zero_mpi_times()
  use mpimod
  use clockmod
  implicit none
  mpitime=0; nonmpitime=0
  call myclock(mpibtime);  mpiatime=mpibtime
end subroutine zero_mpi_times


subroutine getmyranknprocs(outrank,outnprocs)
  use mpimod
  implicit none
  integer,intent(out) :: outrank,outnprocs
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


subroutine make_mpi_comm(grpsize,inranks,OUTCOMM)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(in) :: grpsize, inranks(grpsize)
  integer,intent(out) :: OUTCOMM
#ifndef MPIFLAG
  outcomm=(-798)
  return
  outcomm=outcomm+inranks(1)   !! avoid warn unused
#else
  integer :: ranks(grpsize),ii,MYGROUP,ierr

  do ii=1,grpsize
     if (inranks(grpsize).gt.nprocs.or.inranks(grpsize).lt.1) then
        OFLWR "rank error",inranks(:); CFLST
     endif
  enddo

  ranks(:)=inranks(:)-1
  call MPI_group_incl(MPI_GROUP_WORLD,grpsize,ranks(:),MYGROUP,ierr)
  if (ierr.ne.0) then
     OFLWR "make_comm group_incl ERR ", ierr; CFLST
  endif
  call MPI_comm_create(MPI_COMM_WORLD,MYGROUP,OUTCOMM,ierr)
  if (ierr.ne.0) then
     OFLWR "make_comm group_create err ", ierr; CFLST
  endif
#endif

end subroutine make_mpi_comm


subroutine mpiorbsets()
  use mpimod
  use mpi_orbsetmod
  use parameters
  implicit none
  integer :: iproc,iset
  integer, allocatable :: procsperset(:)
#ifdef MPIFLAG
  integer :: jproc,ierr,MPI_GROUP_ORB
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
  myorbset=1;    firstmpiorb=1;  maxprocsperset=1; orbrank=1;
  nzprocsperset=1
#else
  allocate(process_ranks(nspf*2,norbsets)); process_ranks(:,:)=(-1)
  allocate(MPI_COMM_ORB(norbsets),NZ_COMM_ORB(norbsets))
  MPI_COMM_ORB(:)=0; NZ_COMM_ORB(:)=0

  nzprocsperset=(nspf+orbsperproc-1)/orbsperproc

  firstmpiorb=nspf+1
  iproc=0
  do iset=1,norbsets
     jproc=0
     do while (jproc.lt.procsperset(iset))
        jproc=jproc+1
        iproc=iproc+1
        process_ranks(jproc,iset)=iproc
        if (iproc.eq.myrank) then
           orbrank=jproc
           myorbset=iset
           firstmpiorb=(jproc-1)*orbsperproc+1
        endif
     enddo
     process_ranks(:,iset)=process_ranks(:,iset)-1
     call MPI_group_incl(MPI_GROUP_WORLD,procsperset(iset),process_ranks(:,iset),&
          MPI_GROUP_ORB,ierr)
     if (ierr.ne.0) then
        OFLWR "group_incl ERR ", ierr; CFLST
     endif
     call MPI_comm_create(MPI_COMM_WORLD, MPI_GROUP_ORB,MPI_COMM_ORB(iset),ierr)
     if (ierr.ne.0) then
        OFLWR "group_create err ", ierr; CFLST
     endif
     call MPI_group_incl(MPI_GROUP_WORLD,nzprocsperset,process_ranks(:,iset),&
          MPI_GROUP_ORB,ierr)
     if (ierr.ne.0) then
        OFLWR "group_incl ERR xx", ierr; CFLST
     endif
     call MPI_comm_create(MPI_COMM_WORLD, MPI_GROUP_ORB,NZ_COMM_ORB(iset),ierr)
     if (ierr.ne.0) then
        OFLWR "group_create err xx", ierr; CFLST
     endif
  enddo
  if (myorbset.le.0) then
     OFLWR "SSSETERRROR"; CFLST
  endif

  maxprocsperset=0
  do iset=1,norbsets
     if (procsperset(iset).ge.maxprocsperset) then
        maxprocsperset=procsperset(iset)
     endif
  enddo
  deallocate(process_ranks,procsperset)
#endif

  mpi_orbset_init=1
  mpi_orbset_init=1
  mpi_orbset_init=1

end subroutine mpiorbsets



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
  integer :: myiostat
  if (stdoutflag==0.and.mpifileptr.ne.nullfileptr) then
     open(mpifileptr,file=mpioutfile,status="old", position="append",iostat=myiostat)
     if (myiostat.ne.0) then
        print *, "IOSTAT OPENFILE",myiostat,myrank,mpifileptr,nullfileptr; stop
     endif
  endif
end subroutine openfile


subroutine closefile()
  use mpimod
  use fileptrmod
  implicit none
  if (stdoutflag==0.and.mpifileptr.ne.nullfileptr)  then
     close(mpifileptr)
  endif
end subroutine closefile


subroutine beforebarrier()
  use mpimod
  use parameters
  implicit none
  integer :: ii,qrank,qprocs
  qprocs=(nprocs-1)/turnbatchsize+1
  qrank=(myrank-1)/turnbatchsize+1
  do ii=1,qrank
     call mpibarrier()
  enddo
end subroutine beforebarrier


subroutine afterbarrier()
  use mpimod
  use parameters
  implicit none
  integer :: ii,qrank,qprocs
  qprocs=(nprocs-1)/turnbatchsize+1
  qrank=(myrank-1)/turnbatchsize+1
  do ii=qrank,qprocs
     call mpibarrier()
  enddo
end subroutine afterbarrier


#ifdef MPIFLAG


module orbgathersubmod
use clockmod
contains

subroutine mpiorbgather0(inoutvector,insize,onlynz)
  use mpimod
  use mpi_orbsetmod
  use fileptrmod
  use parameters
  implicit none
  integer,intent(in) :: insize,onlynz
  DATATYPE,intent(inout) :: inoutvector(insize,nspf)
  DATATYPE :: orbvector(insize,nspf*2),workvector(insize,orbsperproc)   !! AUTOMATIC
  integer :: ierr,lastmpiorb

  if (nprocs.eq.1.or.parorbsplit.ne.1) then
     return
  endif

  if (mpi_orbset_init.ne.1) then
     OFLWR "Programmer fail, mpiorbgather called but mpiorbsets() appears "
     WRFL "   not to have been called."; CFLST
  endif

  if (firstmpiorb+orbsperproc-1.gt.nspf*2 .or. orbsperproc*min(nspf,nprocs).gt.nspf*2) then
     OFLWR "YYY ERROR",firstmpiorb,orbsperproc; CFLST
  endif

  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  orbvector(:,:)=0d0;  workvector(:,:)=0d0

  if (firstmpiorb.le.nspf) then
     lastmpiorb=min(nspf,firstmpiorb+orbsperproc-1)
     workvector(:,1:lastmpiorb-firstmpiorb+1)=inoutvector(:,firstmpiorb:lastmpiorb)
  endif

  if (onlynz.eq.0) then
     call mpi_allgather(workvector(:,:),insize*orbsperproc,&
          MPIDATATYPE,orbvector(:,:),insize*orbsperproc,&
          MPIDATATYPE,MPI_COMM_ORB(myorbset),ierr)
  else
     call mpi_allgather(workvector(:,:),insize*orbsperproc,&
          MPIDATATYPE,orbvector(:,:),insize*orbsperproc,&
          MPIDATATYPE,NZ_COMM_ORB(myorbset),ierr)
  endif
  if (ierr.ne.0) then
     OFLWR "ORBGATHER ERR ", ierr; CFLST
  endif

  inoutvector(:,:)=orbvector(:,1:nspf)

  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mpiorbgather0


subroutine mpiorbgather(inoutvector,insize)
  use mpimod
  use mpi_orbsetmod
  use fileptrmod
  use parameters
  implicit none
  integer,intent(in) :: insize
  DATATYPE,intent(inout) :: inoutvector(insize,nspf)
  call mpiorbgather0(inoutvector,insize,0)
end subroutine mpiorbgather


subroutine mpiorbgather_nz(inoutvector,insize)
  use mpimod
  use mpi_orbsetmod
  use fileptrmod
  use parameters
  implicit none
  integer,intent(in) :: insize
  DATATYPE,intent(inout) :: inoutvector(insize,nspf)
  call mpiorbgather0(inoutvector,insize,1)
end subroutine mpiorbgather_nz


subroutine mpiorbreduce(input, isize)
  use mpimod
  use mpi_orbsetmod
  use fileptrmod
  implicit none
  integer,intent(in) :: isize
  DATATYPE,intent(inout) :: input(isize)
  DATATYPE,allocatable :: output(:)
  integer :: ierr

  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  allocate(output(isize))
  output=0;  ierr=0
  call MPI_allreduce( input, output, isize, MPIDATATYPE, MPI_SUM, &
       MPI_COMM_ORB(myorbset), ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR mympireduce!";   CFLST
  endif
  deallocate(output)
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

end subroutine mpiorbreduce

end module orbgathersubmod


subroutine mpistart()
  use mpimod
  use fileptrmod
  use clockmod
  implicit none
  integer :: ierr,  iilen !!,provided,required
  character(len=40) :: format
  integer :: nargs, i, mpioutflag
#ifdef PGFFLAG
  integer :: myiargc
#endif
  character (len=SLN) :: buffer
#ifdef PGFFLAG
  nargs=myiargc()
#else
  nargs=iargc()
#endif
  mpioutflag=0
  do i=1,nargs
     call getarg(i,buffer); 
     if (buffer(1:6) .eq. 'MPIOUT') then
        mpioutflag=1
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
     call system("mkdir -p MPIOUTS")
     if (mpioutfilelen.gt.9) then
        write(format,'(A2,I2,A2,I1,A1)') "(A",mpioutfilelen,",I",iilen(myrank),")"
     else
        write(format,'(A2,I1,A2,I1,A1)') "(A",mpioutfilelen,",I",iilen(myrank),")"
     endif
     write(mpioutfile,format) mpioutfilebase,myrank
!! nevermind writing out every one
     if (mpioutflag.eq.0.and.myrank.ne.1) then
        mpifileptr=nullfileptr
        mpioutfile="garmonbozia"
     else
        open(mpifileptr,file=mpioutfile,status="unknown"); WRFL; CFL
     endif
  endif
  OFLWR; CFL
  call myclock(mpibtime);  mpiatime=mpibtime

end subroutine mpistart



subroutine mpibarrier()
  use mpimod
  use fileptrmod
  use clockmod
  implicit none
  integer :: ierr
  call myclock(mpiatime)
  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if (ierr/=0) then
     OFLWR "MPI ERR 3"; CFLST
  endif
  call myclock(mpibtime)
  mpitime=mpitime+mpibtime-mpiatime
end subroutine mpibarrier


subroutine mpistop()
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr
  call waitawhile()
  call mpi_finalize(ierr)
  call waitawhile()
  if (ierr/=0) then
     OFLWR "MPI ERR 3"; CFL
  endif
  OFLWR "MCTDHF STOP!";CFL
  stop  !! STOP
end subroutine mpistop


subroutine mpiabort()
  use mpimod
  use fileptrmod
  implicit none
  integer :: ierr
  call mpi_abort(ierr)
  OFLWR "MCTDHF ABORT!";CFL
  stop
end subroutine mpiabort


module mpisubmod
  use clockmod

contains

subroutine mympireduce_local(input, isize, IN_COMM)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(in) :: isize,IN_COMM
  DATATYPE,intent(inout) :: input(isize)
  DATATYPE,allocatable :: output(:)
  integer :: ierr

  if (nprocs.eq.1) then
     return
  endif
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  allocate(output(isize))
  output=0; ierr=0
  call MPI_allreduce(input,output,isize,MPIDATATYPE,MPI_SUM,IN_COMM,ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR mympireduce!";   CFLST
  endif
  deallocate(output)
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

end subroutine mympireduce_local


subroutine mympireduce(input, isize)
  use mpimod
  implicit none
  integer,intent(in) :: isize
  DATATYPE,intent(inout) :: input(isize)
  call mympireduce_local(input, isize, MPI_COMM_WORLD)
end subroutine mympireduce


subroutine mympirealreduce(input, isize)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(in) :: isize
  real*8,intent(inout) :: input(isize)
  real*8,allocatable :: output(:)
  integer :: ierr
  if (nprocs.eq.1) then
     return
  endif
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  allocate(output(isize))
  ierr=0; output=0
  call MPI_allreduce( input, output, isize, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_WORLD , ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR mympireduce!";   CFLST
  endif
  deallocate(output)
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympirealreduce


subroutine mympireduceto_local(input, output, isize, dest, IN_COMM)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(in) :: isize,dest,IN_COMM
  DATATYPE,intent(in) :: input(isize)
  DATATYPE,intent(out) :: output(isize) 
  integer :: ierr, idest

  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  if (dest.lt.1.or.dest.gt.nprocs) then
     OFLWR "ERR DEST REDUCETO",dest,nprocs; CFLST
  endif
  ierr=0
  idest=dest-1
  call MPI_reduce( input, output, isize, MPIDATATYPE, MPI_SUM, idest, &
       IN_COMM, ierr)
  if (ierr/=0) then
     OFLWR "ERR mpi_reduce!";   CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympireduceto_local


subroutine mympireduceto(input, output, isize, dest)
  use mpimod
  implicit none
  integer,intent(in) :: isize,dest
  DATATYPE,intent(in) :: input(isize)
  DATATYPE,intent(out) :: output(isize) 
  call mympireduceto_local(input,output,isize,dest,MPI_COMM_WORLD)
end subroutine mympireduceto


!! used in dGMRES parallel routines

subroutine mympirealreduceone_local(input,IN_COMM)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(in) :: IN_COMM
  real*8,intent(inout) :: input
  real*8 :: output
  integer :: ierr
  if (nprocs.eq.1) then
     return
  endif
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  ierr=0
  call MPI_allreduce(input,output,1,MPI_DOUBLE_PRECISION,MPI_SUM,IN_COMM,ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR mympirealreduce!"; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympirealreduceone_local


subroutine mympiireduceone_local(input,IN_COMM)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(in) :: IN_COMM
  integer,intent(inout) :: input
  integer :: output,ierr
  if (nprocs.eq.1) then
     return
  endif
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  ierr=0
  call MPI_allreduce(input,output,1,MPI_INTEGER,MPI_SUM,IN_COMM,ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR mympiireduce!"; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympiireduceone_local


subroutine mympii8reduceone_local(input,IN_COMM)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(in) :: IN_COMM
  integer*8,intent(inout) :: input
  integer*8 :: output
  integer :: ierr
  if (nprocs.eq.1) then
     return
  endif
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  ierr=0
  call MPI_allreduce(input,output,1,MPI_INTEGER8,MPI_SUM,IN_COMM,ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR mympii8reduce!"; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympii8reduceone_local


subroutine mympirealreduceone(input)
  use mpimod
  implicit none
  real*8,intent(inout) :: input
  call mympirealreduceone_local(input,MPI_COMM_WORLD)
end subroutine mympirealreduceone


subroutine mympiireduceone(input)
  use mpimod
  implicit none
  integer,intent(inout) :: input
  call mympiireduceone_local(input,MPI_COMM_WORLD)
end subroutine mympiireduceone


subroutine mympii8reduceone(input)
  use mpimod
  implicit none
  integer*8,intent(inout) :: input
  call mympii8reduceone_local(input,MPI_COMM_WORLD)
end subroutine mympii8reduceone


subroutine mympiireduce(input, isize)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(in) :: isize
  integer,intent(inout) :: input(isize)
  integer :: output(isize) , ierr
  if (nprocs.eq.1) then
     return
  endif
  output=0
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  ierr=0
  call MPI_allreduce( input, output, isize, MPI_INTEGER, MPI_SUM, &
       MPI_COMM_WORLD , ierr)
  if (ierr/=0) then
     OFLWR "ERR mympiireduce!";   CFLST
  endif
  input=output
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympiireduce


subroutine mympireduceone_local(input,IN_COMM)
  use mpimod
  use fileptrmod
  implicit none
  integer, intent(in) :: IN_COMM
  DATATYPE, intent(inout) :: input
  DATATYPE :: output
  integer :: ierr 
  if (nprocs.eq.1) then
     return
  endif
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  ierr=0
  call MPI_allreduce( input, output, 1, MPIDATATYPE, MPI_SUM,IN_COMM, ierr)
  if (ierr/=0) then
     OFLWR "ERR mympireduce!"; CFLST
  endif
  input=output
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympireduceone_local


subroutine mympireduceone(input)
  use mpimod
  implicit none
  DATATYPE,intent(inout) :: input
  call mympireduceone_local(input,MPI_COMM_WORLD)
end subroutine mympireduceone


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
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_sendrecv(sendbuf,isize,MPI_DOUBLE_COMPLEX,idest,tag,&
       recvbuf,isize, MPI_DOUBLE_COMPLEX,isource,tag,MPI_COMM_LOCAL,MPI_STATUS_IGNORE,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympisendrecv",ierr; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympisendrecv_complex_local


subroutine mympisendrecv_real_local(sendbuf, recvbuf, dest, source, tag, isize,MPI_COMM_LOCAL)
  use mpimod
  use fileptrmod
  implicit none
  integer, intent(in) :: dest,source,tag,isize,MPI_COMM_LOCAL
  integer :: ierr, idest,isource
  real*8, intent(in) :: sendbuf(isize)
  real*8, intent(out) :: recvbuf(isize)
  idest=dest-1
  isource=source-1
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_sendrecv(sendbuf,isize,MPI_DOUBLE_PRECISION,idest,tag,&
       recvbuf,isize, MPI_DOUBLE_PRECISION,isource,tag,MPI_COMM_LOCAL,MPI_STATUS_IGNORE,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympisendrecv_real_local",ierr; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympisendrecv_real_local


subroutine mympisendrecv_complex(sendbuf, recvbuf, dest, source, tag, isize)
  use mpimod
  implicit none
  integer, intent(in) :: dest,source,tag,isize
  complex*16, intent(in) :: sendbuf(isize)
  complex*16, intent(out) :: recvbuf(isize)
  call mympisendrecv_complex_local(sendbuf, recvbuf, dest, source, tag, isize,MPI_COMM_WORLD)
end subroutine mympisendrecv_complex


subroutine mympisendrecv_local(sendbuf, recvbuf, dest, source, tag, isize,MPI_COMM_LOCAL)
  implicit none
  integer, intent(in) :: dest,source,tag,isize,MPI_COMM_LOCAL
  DATATYPE, intent(in) :: sendbuf(isize)
  DATATYPE, intent(out) :: recvbuf(isize)
#ifdef REALGO
  call mympisendrecv_real_local(sendbuf, recvbuf, dest, source, tag, isize,MPI_COMM_LOCAL)
#else
  call mympisendrecv_complex_local(sendbuf, recvbuf, dest, source, tag, isize,MPI_COMM_LOCAL)
#endif
end subroutine mympisendrecv_local


subroutine mympisendrecv(sendbuf, recvbuf, dest, source, tag, isize)
  use mpimod
  implicit none
  integer, intent(in) :: dest,source,tag,isize
  DATATYPE, intent(in) :: sendbuf(isize)
  DATATYPE, intent(out) :: recvbuf(isize)
  call mympisendrecv_local(sendbuf, recvbuf, dest, source, tag, isize,MPI_COMM_WORLD)
end subroutine mympisendrecv


subroutine mympisend(input, dest, tag, isize)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(in) :: dest,tag,isize
  DATATYPE,intent(in) :: input(isize)
  integer :: ierr, idest
  idest=dest-1
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_send(input,isize,MPIDATATYPE,idest,tag,MPI_COMM_WORLD,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympisend"; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympisend


subroutine mympirecv(output, source, tag, isize)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(in) :: source,tag,isize
  integer :: ierr,isource
  DATATYPE,intent(out) :: output(isize)
  isource=source-1
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_recv(output,isize,MPIDATATYPE,isource,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympirecv"; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympirecv


subroutine mympirealbcast_local(input, source, isize,MPI_COMM_LOCAL)
  use mpimod
  use fileptrmod
  implicit none
  integer, intent(in) ::  isize, source, MPI_COMM_LOCAL
  real*8, intent(inout) :: input(isize)
  integer :: ierr, isource

  isource=source-1
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_bcast(input,isize,MPI_DOUBLE_PRECISION,isource,MPI_COMM_LOCAL,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympirealbcast"; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympirealbcast_local


subroutine mympicomplexbcast_local(input, source, isize, MPI_COMM_LOCAL)
  use mpimod
  use fileptrmod
  implicit none
  integer, intent(in) ::  isize, source, MPI_COMM_LOCAL
  complex*16, intent(inout) :: input(isize)
  integer :: ierr, isource

  isource=source-1
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_bcast(input,isize,MPI_DOUBLE_COMPLEX,isource,MPI_COMM_LOCAL,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympicomplexbcast"; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympicomplexbcast_local


subroutine mympibcast_local(input, source, isize,INCOMM)
  use mpimod
  implicit none
  integer,intent(in) :: isize, source,INCOMM
  DATATYPE,intent(inout) :: input(isize)
#ifdef REALGO
  call mympirealbcast_local(input, source, isize,INCOMM)
#else
  call mympicomplexbcast_local(input, source, isize,INCOMM)
#endif

end subroutine mympibcast_local


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


subroutine mympirealbcast(input, source, isize)
  use mpimod
  implicit none
  integer,intent(in) :: isize, source
  real*8,intent(inout) :: input(isize)
  call mympirealbcast_local(input, source, isize,MPI_COMM_WORLD)
end subroutine mympirealbcast


subroutine mympicomplexbcast(input, source, isize)
  use mpimod
  implicit none
  integer,intent(in) :: isize, source
  complex*16,intent(inout) :: input(isize)
  call mympicomplexbcast_local(input, source, isize,MPI_COMM_WORLD)
end subroutine mympicomplexbcast


subroutine mympiibcast(input, source, isize)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(in) :: source,isize
  integer,intent(inout) :: input(isize) 
  integer :: ierr, isource
  isource=source-1
  call myclock(mpiatime); nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_bcast(input,isize,MPI_INTEGER,isource,MPI_COMM_WORLD,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympibcast"; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympiibcast


subroutine mympirealbcastone(input, source)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(in) :: source
  real*8,intent(inout) :: input
  integer :: ierr,isource
  isource=source-1
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_bcast(input,1,MPI_DOUBLE_PRECISION,isource,MPI_COMM_WORLD,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympibcast"; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympirealbcastone


subroutine mympiibcastone(input, source)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(in) :: source
  integer,intent(inout) :: input
  integer :: ierr, isource
  isource=source-1
  call myclock(mpiatime); nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_bcast(input,1,MPI_INTEGER,isource,MPI_COMM_WORLD,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympibcast"; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympiibcastone


subroutine mympilogbcast(input, source, isize)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(in) :: source,isize
  logical,intent(inout) :: input(isize)
  integer :: ierr, isource
  isource=source-1
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_bcast(input,isize,MPI_LOGICAL,isource,MPI_COMM_WORLD,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympibcast"; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympilogbcast


subroutine mympiimax(input)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(inout) :: input
  integer :: ierr,output
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call MPI_reduce( input, output, 1, MPI_INTEGER, MPI_MAX,0, MPI_COMM_WORLD , ierr)
  if (ierr/=0) then
     OFLWR "ERR MPIMAX!"; CFLST
  endif
  call mpi_bcast(output,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR MPIIMAX!"; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympiimax


subroutine mympiimin(input)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(inout) :: input
  integer :: ierr, output
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call MPI_reduce( input, output, 1, MPI_INTEGER, MPI_MIN,0, MPI_COMM_WORLD , ierr)
  if (ierr/=0) then
     OFLWR "ERR MPIImin!"; CFLST
  endif
  call mpi_bcast(output,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR MPIImin!"; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympiimin


subroutine mympimax_local(input,IN_COMM)
  use mpimod
  use fileptrmod
  implicit none
  integer,intent(in) :: IN_COMM
  real*8,intent(inout) :: input
  integer :: ierr
  real*8 :: output
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call MPI_reduce(input,output,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,IN_COMM,ierr)
  if (ierr/=0) then
     OFLWR "ERR MPIMAX!"; CFLST
  endif
  call mpi_bcast(output,1,MPI_DOUBLE_PRECISION,0,IN_COMM,ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR MPIMAX!"; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympimax_local


subroutine mympimax(input)
  use mpimod
  implicit none
  real*8,intent(inout) :: input
  call mympimax_local(input,MPI_COMM_WORLD)
end subroutine mympimax


subroutine mympimin(input)
  use mpimod
  use fileptrmod
  implicit none
  real*8,intent(inout) :: input
  integer :: ierr
  real*8 :: output
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call MPI_reduce( input, output, 1, MPI_DOUBLE_PRECISION, MPI_MIN,0, &
       MPI_COMM_WORLD , ierr)
  if (ierr/=0) then
     OFLWR "ERR MPImin!"; CFLST
  endif
  call mpi_bcast(output,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR MPIImin!"; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympimin




subroutine simpleallgather(vectorin,vectorsout,insize)
  use mpimod
  implicit none
  integer,intent(in) :: insize
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
  integer,intent(in) :: insize
  complex*16,intent(in) :: vectorin(insize)
  complex*16,intent(out) :: vectorsout(insize,nprocs)
  integer :: ierr
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_allgather(vectorin(:),insize,MPI_DOUBLE_COMPLEX,vectorsout(:,:),&
       insize, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "gather complex ERR ", ierr; call mpistop()
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine simpleallgather_complex


subroutine simpleallgather_real(vectorin,vectorsout,insize)
  use mpimod
  use fileptrmod
  implicit none

  integer,intent(in) :: insize
  real*8,intent(in) :: vectorin(insize)
  real*8,intent(out) :: vectorsout(insize,nprocs)
  integer :: ierr
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_allgather(vectorin(:),insize,MPI_DOUBLE_PRECISION,vectorsout(:,:),&
       insize,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "gather real ERR ", ierr; call mpistop()
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine simpleallgather_real


subroutine getgatherv_stuff(blocks,parsize,blockstart,numprocs)
  implicit none
  integer, intent(in) :: numprocs,blocks(numprocs)
  integer, intent(out) :: parsize, blockstart(numprocs)
  integer :: ii, isum
  isum=blocks(1)
  blockstart(1)=1  
  do ii=2,numprocs
     blockstart(ii)=blockstart(ii-1)+blocks(ii-1)
     isum=isum+blocks(ii)
  enddo
  parsize=isum
end subroutine getgatherv_stuff


SUBROUTINE MYGATHERV(V1,X1,BLOCKS,bcastflag)
  use mpimod
  IMPLICIT NONE
  logical, intent(in) :: bcastflag
  INTEGER, INTENT(IN) :: BLOCKS(nprocs)
  DATATYPE, INTENT(IN) :: V1(blocks(myrank))
  DATATYPE, INTENT(OUT) :: X1(*)
#ifdef REALGO
  call mygatherv_real(v1,x1,blocks,bcastflag)
#else
  call mygatherv_complex(v1,x1,blocks,bcastflag)
#endif

END SUBROUTINE MYGATHERV


SUBROUTINE MYGATHERV_complex(V1,X1,BLOCKS,bcastflag)
  use mpimod
  IMPLICIT NONE
  logical, intent(in) :: bcastflag
  INTEGER, INTENT(IN) :: BLOCKS(nprocs)
  complex*16, INTENT(IN) :: V1(blocks(myrank))
  complex*16, INTENT(OUT) :: X1(*)
  INTEGER ::MPIERR=0, BLOCKSTART(nprocs),parsize

  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  call getgatherv_stuff(blocks,parsize,blockstart,nprocs)

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
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

END SUBROUTINE MYGATHERV_complex


SUBROUTINE MYGATHERV_real(V1,X1,BLOCKS,bcastflag)
  use mpimod
  IMPLICIT NONE
  logical, intent(in) :: bcastflag
  INTEGER, INTENT(IN) :: BLOCKS(nprocs)
  real*8, INTENT(IN) :: V1(blocks(myrank))
  real*8, INTENT(OUT) :: X1(*)
  INTEGER ::MPIERR=0, parsize, BLOCKSTART(nprocs)

  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  call getgatherv_stuff(blocks,parsize,blockstart,nprocs)

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
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

END SUBROUTINE MYGATHERV_real



SUBROUTINE MYSCATTERV(X1,V1,BLOCKS)
  use mpimod
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: BLOCKS(nprocs)
  DATATYPE, INTENT(OUT) :: V1(blocks(myrank))
  DATATYPE, INTENT(IN) :: X1(*)
#ifdef REALGO
  call myscatterv_real(x1,v1,blocks)
#else
  call myscatterv_complex(x1,v1,blocks)
#endif
end SUBROUTINE MYSCATTERV


SUBROUTINE MYSCATTERV_complex(X1,V1,BLOCKS)
  use mpimod
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: BLOCKS(nprocs)
  complex*16, INTENT(OUT) :: V1(blocks(myrank))
  complex*16, INTENT(IN) :: X1(*)
  INTEGER ::MPIERR=0, parsize, BLOCKSTART(nprocs)

  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  call getgatherv_stuff(blocks,parsize,blockstart,nprocs)

  CALL MPI_SCATTERV(X1, blocks, blockstart(:)-1, MPI_DOUBLE_COMPLEX, V1, &
       BLOCKS(myrank), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, MPIERR)
  if (mpierr.ne.0) then
     print *, "Scatterv err",mpierr,myrank;          stop
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

END SUBROUTINE MYSCATTERV_complex


SUBROUTINE MYSCATTERV_real(X1,V1,BLOCKS)
  use mpimod
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: BLOCKS(nprocs)
  real*8, INTENT(OUT) :: V1(blocks(myrank))
  real*8, INTENT(IN) :: X1(*)
  INTEGER ::MPIERR=0, parsize, BLOCKSTART(nprocs)

  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  call getgatherv_stuff(blocks,parsize,blockstart,nprocs)

  CALL MPI_SCATTERV(X1, blocks, blockstart(:)-1, MPI_DOUBLE_PRECISION, V1, &
       BLOCKS(myrank), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPIERR)
  if (mpierr.ne.0) then
     print *, "Scatterv err",mpierr,myrank;          stop
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

END SUBROUTINE MYSCATTERV_real

!! MODULE mpisubmod continues

#else

module orbgathersubmod
contains

subroutine mpiorbgather(orbvector,insize)
  implicit none
  integer :: insize
  DATATYPE :: orbvector(insize)
  return
  orbvector(1)=orbvector(1)
end subroutine mpiorbgather


subroutine mpiorbgather_nz(orbvector,insize)
  implicit none
  integer :: insize
  DATATYPE :: orbvector(insize)
  return
  orbvector(1)=orbvector(1)
end subroutine mpiorbgather_nz


subroutine mpiorbreduce(orbvector,insize)
  implicit none
  integer :: insize
  DATATYPE :: orbvector(insize)
  return
  orbvector(1)=orbvector(1)
end subroutine mpiorbreduce

end module orbgathersubmod


subroutine mpistart
  use mpimod
  use fileptrmod
  implicit none
  myrank=1
  nprocs=1
  stdoutflag=1
  mpifileptr=6

end subroutine mpistart


subroutine mpibarrier()
end subroutine mpibarrier


subroutine mpistop()
  use fileptrmod
  implicit none
  call waitawhile()
  OFLWR "MCTDHF STOP!";CFL
  stop
end subroutine mpistop

subroutine mpiabort()
  use fileptrmod
  implicit none
  OFLWR "MCTDHF ABORT!";CFL
  stop
end subroutine mpiabort


module mpisubmod
contains

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

subroutine mympimax_local(input,incomm)
  implicit none
  real*8 :: input
  integer :: incomm
  return
  input=input; incomm=incomm
end subroutine mympimax_local

subroutine mympimin(input)
  implicit none
  real*8 :: input
  return
  input=input
end subroutine mympimin



SUBROUTINE MYGATHERV(V1,X1,BLOCKS,bcastflag)
  IMPLICIT NONE
  logical, intent(in) :: bcastflag
  INTEGER, INTENT(IN) :: BLOCKS(1)
  DATATYPE, INTENT(IN) :: V1(blocks(1))
  DATATYPE, INTENT(OUT) :: X1(blocks(1))
  logical :: ilog
  X1(:)=V1(:)
  return
  ilog=bcastflag   !!avoid warn unused
END SUBROUTINE MYGATHERV



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




SUBROUTINE MYSCATTERV(X1,V1,BLOCKS)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: BLOCKS(1)
  DATATYPE, INTENT(OUT) :: V1(blocks(1))
  DATATYPE, INTENT(IN) :: X1(blocks(1))
  V1(:)=X1(:)
END SUBROUTINE MYSCATTERV


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


subroutine mympirealreduceone_local(input,incomm)
  implicit none
  real*8 :: input
  integer :: incomm
  return
  input=input; incomm=incomm
end subroutine


subroutine mympireduceone(input)
  implicit none
  DATATYPE :: input
  return
  input=input
end subroutine mympireduceone


subroutine mympireduceone_local(input,incomm)
  implicit none
  DATATYPE :: input
  integer :: incomm
  return
  input=input; incomm=incomm
end subroutine mympireduceone_local


subroutine mympireduce(input, isize)
  implicit none
  integer :: isize
  DATATYPE :: input(isize)
  return
  input(1)=input(1)
end subroutine mympireduce


subroutine mympireduce_local(input, isize,incomm)
  implicit none
  integer :: isize,incomm
  DATATYPE :: input(isize)
  return
  input(1)=input(1); incomm=incomm
end subroutine mympireduce_local


subroutine mympirealreduce(input, isize)
  implicit none
  integer :: isize
  real*8 :: input(isize)
  return
  input(1)=input(1)
end subroutine mympirealreduce

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


subroutine mympibcast_local(input, source, isize, IN_COMM)
  implicit none
  integer :: source,isize, IN_COMM
  DATATYPE :: input(isize)
  return
  source=source;  input(1)=input(1); IN_COMM=IN_COMM
end subroutine mympibcast_local


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


subroutine mympiireduceone_local(input,incomm)
  implicit none
  integer :: input,incomm
  return
  input=input; incomm=incomm
end subroutine mympiireduceone_local


subroutine mympii8reduceone(input)
  implicit none
  integer*8 :: input
  return
  input=input
end subroutine mympii8reduceone


subroutine mympii8reduceone_local(input,incomm)
  implicit none
  integer*8 :: input
  integer :: incomm
  return
  input=input; incomm=incomm
end subroutine mympii8reduceone_local


subroutine mympisendrecv(sendbuf, recvbuf, dest, source, tag, isize)
  use mpimod
  use fileptrmod
  implicit none
  integer :: isize, dest, source,tag
  DATATYPE :: sendbuf(isize),recvbuf(isize)
  if (dest.ne.1.or.source.ne.1) then
     OFLWR "Error nonpar mympisendrecv wrapper",dest,source; CFLST
  endif
  recvbuf(:)=sendbuf(:)
  return
  tag=tag 
end subroutine mympisendrecv


subroutine mympisendrecv_local(sendbuf, recvbuf, dest, source, tag, isize,INCOMM)
  use mpimod
  use fileptrmod
  implicit none
  integer :: isize, dest, source,tag,incomm
  DATATYPE :: sendbuf(isize),recvbuf(isize)
  if (dest.ne.1.or.source.ne.1) then
     OFLWR "Error nonpar mympisendrecv_local wrapper",dest,source; CFLST
  endif
  recvbuf(:)=sendbuf(:)
  return
  tag=tag ; tag=incomm
end subroutine mympisendrecv_local

#endif

!! MODULE MPISUBMOD CONTINUES

subroutine mpiallgather_local(inout,totsize,blocksizes,notusedint,&
     IN_COMM,numprocs,locrank)
  use mpimod
  use fileptrmod
  implicit none
  integer, intent(in) :: numprocs,totsize,blocksizes(numprocs),&
       notusedint,IN_COMM,locrank
  DATATYPE,intent(inout) :: inout(totsize)
#ifdef MPIFLAG
  DATATYPE :: work(blocksizes(locrank))
  integer :: ierr,blockstart(numprocs),icount
  if (numprocs.eq.1) then
     return
  endif
  work=0; blockstart=0
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call getgatherv_stuff(blocksizes,icount,blockstart,numprocs)
  if (icount.ne.totsize) then
     OFLWR "ALLgather ERR count", icount,totsize; 
     WRFL blocksizes;
     CFLST
  endif
  work(:)=inout(blockstart(locrank):blockstart(locrank)+blocksizes(locrank)-1)
  call mpi_allgatherv(work,blocksizes(locrank),MPIDATATYPE,inout,blocksizes(:),&
       blockstart(:)-1,MPIDATATYPE,IN_COMM,ierr)
  if (ierr.ne.0) then
     OFLWR "ALLGATHERv ERR ", ierr; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
#endif
end subroutine mpiallgather_local


subroutine mpiallgather(inout,totsize,blocksizes,notusedint)
  use mpimod
  implicit none
  integer, intent(in) :: totsize,blocksizes(nprocs),notusedint
  DATATYPE,intent(inout) :: inout(totsize)
#ifdef MPIFLAG
  call mpiallgather_local(inout,totsize,blocksizes,notusedint,&
     MPI_COMM_WORLD,nprocs,myrank)
#endif
end subroutine mpiallgather


subroutine mpiallgather_i(inout,totsize,blocksizes,notusedint)
  use mpimod
  use fileptrmod
  implicit none
  integer, intent(in) :: totsize,blocksizes(nprocs),notusedint
  integer,intent(inout) :: inout(totsize)
#ifdef MPIFLAG
  integer :: icount, work(blocksizes(myrank))
  integer :: ierr,blockstart(nprocs)
  if (nprocs.eq.1) then
     return
  endif
  work=0; blockstart=0
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call getgatherv_stuff(blocksizes,icount,blockstart,nprocs)
  if (icount.ne.totsize) then
     OFLWR "ALLgather_i ERR count", icount,totsize; 
     WRFL blocksizes;
     CFLST
  endif
  work(:)=inout(blockstart(myrank):blockstart(myrank)+blocksizes(myrank)-1)
  call mpi_allgatherv(work,blocksizes(myrank),MPI_INTEGER,inout,blocksizes(:),&
       blockstart(:)-1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  if (ierr.ne.0) then
     OFLWR "ALLGATHERv ERR ", ierr; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
#endif
end subroutine mpiallgather_i


subroutine mympialltoall(input, output, count)
  use fileptrmod
  use mpimod
  implicit none
  integer,intent(in) :: count
  DATATYPE,intent(in) :: input(count,nprocs)
  DATATYPE, intent(out) :: output(count,nprocs)
#ifndef MPIFLAG
  output(:,:)=input(:,:)
#else
  integer :: ierr
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  call mpi_alltoall(input(:,:),count,MPIDATATYPE,output(:,:),count,&
       MPIDATATYPE,MPI_COMM_WORLD,ierr)

  if (ierr.ne.0) then
     OFLWR "ERROR ALLTOALL ", ierr; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
#endif
end subroutine mympialltoall



subroutine mympialltoall_complex(input, output, count)
  use fileptrmod
  use mpimod
  implicit none
  integer,intent(in) :: count
  complex*16, intent(in) :: input(count,nprocs)
  complex*16, intent(out) :: output(count,nprocs)
#ifndef MPIFLAG
  output(:,:)=input(:,:)
#else
  integer :: ierr
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_alltoall(input(:,:),count,MPI_DOUBLE_COMPLEX,output(:,:),count,&
       MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
  if (ierr.ne.0) then
     OFLWR "ERROR ALLTOALL ", ierr; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
#endif
end subroutine mympialltoall_complex



subroutine mympialltoall_complex_local(input, output, count,INCOMM)
  use fileptrmod
  use mpimod
  implicit none
  integer,intent(in) :: count,INCOMM
  complex*16, intent(in) :: input(count,nprocs)
  complex*16, intent(out) :: output(count,nprocs)
#ifndef MPIFLAG
  output(:,:)=input(:,:)
#else
  integer :: ierr
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_alltoall(input(:,:),count,MPI_DOUBLE_COMPLEX,output(:,:),count,&
       MPI_DOUBLE_COMPLEX,INCOMM,ierr)
  if (ierr.ne.0) then
     OFLWR "ERROR ALLTOALLlocal ", ierr; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
#endif
end subroutine mympialltoall_complex_local



subroutine mympialltoall_local(input, output, count,INCOMM)
  use fileptrmod
  use mpimod
  implicit none
  integer,intent(in) :: count,INCOMM
  DATATYPE, intent(in) :: input(count,nprocs)
  DATATYPE, intent(out) :: output(count,nprocs)
#ifndef MPIFLAG
  output(:,:)=input(:,:)
#else
  integer :: ierr
  call myclock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  call mpi_alltoall(input(:,:),count,MPIDATATYPE,output(:,:),count,&
       MPIDATATYPE,INCOMM,ierr)

  if (ierr.ne.0) then
     OFLWR "ERROR ALLTOALL   local ", ierr; CFLST
  endif
  call myclock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
#endif
end subroutine mympialltoall_local

end module mpisubmod


subroutine checkstopfile()
  use mpisubmod
  use fileptrmod
  implicit none
  integer :: iii
  open(222,file="stop",status="old", iostat=iii)
  if (iii==0) then
     close(222)
  endif
  call mympiimax(iii)
  if (iii==0) then
     OFLWR "Stopping due to stopfile!"; CFLST
  endif
end subroutine checkstopfile


