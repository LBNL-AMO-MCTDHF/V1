
!! PARALLEL SUBROUTINES.

#include "Definitions.INC"

subroutine getmyranknprocs(outrank,outnprocs)
  use mpimod
  implicit none
  integer :: outrank,outnprocs
  outrank=myrank
  outnprocs=nprocs
end subroutine getmyranknprocs

subroutine mpiorbsets()
  use mpimod
  use parameters
  implicit none
  integer :: iproc,iset,jproc,ierr
  integer, allocatable :: procsperset(:),process_ranks(:,:)

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
#else

  call MPI_Comm_group(MPI_COMM_WORLD,MPI_GROUP_WORLD,ierr)
  if (ierr.ne.0) then
     OFLWR "comm group ERR ", ierr; CFLST
  endif

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
!!???           firstmpiorb=real(jproc-1)*orbsperproc+1
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

end subroutine mpiorbsets




#ifdef MPIFLAG


subroutine mpiorbgather(orbvector,insize)    !! insize=spfsize except debug
  use mpimod
  use parameters
  implicit none

  integer :: ierr,insize
  DATATYPE :: orbvector(insize,nspf*2),tempvector(insize,firstmpiorb:firstmpiorb+orbsperproc-1)

  tempvector(:,:)=orbvector(:,firstmpiorb:firstmpiorb+orbsperproc-1)

#ifndef MPIFLAG
  return
#endif
  if (nprocs.eq.1) then
     return
  endif

  if (firstmpiorb+orbsperproc-1.gt.nspf*2 .or. orbsperproc*min(nspf,nprocs).gt.nspf*2) then
     OFLWR "YYY ERROR",firstmpiorb,orbsperproc; CFLST
  endif

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime


#ifdef REALGO
  call mpi_allgather(tempvector(:,:),insize*orbsperproc,MPI_DOUBLE_PRECISION,orbvector(:,:),insize*orbsperproc,MPI_DOUBLE_PRECISION,MPI_COMM_ORB(myorbset),ierr)
#else
  call mpi_allgather(tempvector(:,:),insize*orbsperproc,MPI_DOUBLE_COMPLEX,orbvector(:,:),insize*orbsperproc,MPI_DOUBLE_COMPLEX,MPI_COMM_ORB(myorbset),ierr)
#endif

  if (ierr.ne.0) then
     OFLWR "ORBGATHER ERR ", ierr; CFLST
  endif

  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime


end subroutine mpiorbgather






subroutine mpistart()
  use mpimod
  use parameters
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

!!$!  required=MPI_THREAD_FUNNELED
!!$  required=MPI_THREAD_SINGLE
!!$
!!$!$OMP PARALLEL
!!$!$OMP MASTER
!!$  call MPI_INIT_THREAD(required,provided,ierr)
!!$  if (provided.ne.required) then
!!$     print *, "MPI thread support not granted.  probably ok.  comment me out in code to continue.",&
!!$          provided,required; stop
!!$  endif
!!$  if (ierr/=0) then; print *,  "MPI ERR 1";  stop; 
!!$  endif
!!$!$OMP END MASTER
!!$!$OMP END PARALLEL

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
  use parameters
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
  use parameters
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
  use parameters
  implicit none
  DATATYPE :: input(isize)
  DATATYPE :: output(isize) 
  integer :: ierr, isize  !!, iroot

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  ierr=0

  call MPI_allreduce( input, output, isize, MPIDATATYPE, MPI_SUM, MPI_COMM_WORLD , ierr)

  input=output
  if (ierr/=0) then
     OFLWR "ERR mympireduce!";     call closefile(); CFLST
  endif

  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

end subroutine mympireduce


!! used in dGMRES parallel routines

subroutine mympirealreduceone(input)
  use mpimod
  use parameters
  implicit none
  real*8 :: input,output
  integer :: ierr,isize=1
  ierr=0
  call MPI_allreduce( input, output, isize, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD , ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR mympirealreduce!"; CFLST
  endif
end subroutine mympirealreduceone


subroutine mympirealreduce(input,isize)
  use mpimod
  use parameters
  implicit none
  integer :: ierr,isize
  real*8 :: input(isize),output(isize)

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  ierr=0
  call MPI_allreduce( input, output, isize, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD , ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR mympirealreduce!";  CFLST
  endif

  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympirealreduce



subroutine mympiireduceone(input)
  use mpimod
  use parameters
  implicit none
  integer :: input,output
  integer :: ierr,isize=1

  ierr=0
  call MPI_allreduce( input, output, isize, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD , ierr)
  input=output
  if (ierr/=0) then
     OFLWR "ERR mympiireduce!"; CFLST
  endif
end subroutine mympiireduceone


subroutine mympiireduce(input, isize)
  use mpimod
  use parameters
  implicit none
  integer :: ierr, isize  !!, iroot
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
  use parameters
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
  use parameters
  implicit none
  DATATYPE :: input(isize)
  integer :: ierr, isize, source, isource

  isource=source-1
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_bcast(input,isize,MPIDATATYPE,isource,MPI_COMM_WORLD,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympibcast"; CFLST
  endif

  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympibcast


subroutine mympisend(input, dest, tag, isize)
  use mpimod
  use parameters
  implicit none
  DATATYPE :: input(isize)
  integer :: ierr, isize, dest, idest,tag

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
  use parameters
  implicit none
  DATATYPE :: input(isize)
  integer :: ierr, isize, source, isource,tag

  isource=source-1
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_recv(input,isize,MPIDATATYPE,isource,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympirecv"; CFLST
  endif

  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympirecv




subroutine mympirealbcast(input, source, isize)
  use mpimod
  use parameters
  implicit none
  real*8 :: input(isize)
  integer :: ierr, isize, source, isource

  isource=source-1
  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
  call mpi_bcast(input,isize,MPI_DOUBLE_PRECISION,isource,MPI_COMM_WORLD,ierr)
  if (ierr/=0) then
     OFLWR "ERR mympibcast"; CFLST
  endif

  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
end subroutine mympirealbcast


subroutine mympiibcast(input, source, isize)
  use mpimod
  use parameters
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
  use parameters
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
  use parameters
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
  use parameters
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
  use parameters
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
  use parameters
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
  use parameters
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
  use parameters
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


SUBROUTINE MYGATHERV_complex(V1,X1,MINPUT,FIRST,LAST,BLOCKS,BLOCKSTART,bcastflag)
  use mpimod
  IMPLICIT NONE
  logical, intent(in) :: bcastflag
  INTEGER, INTENT(IN) :: MINPUT, FIRST, LAST
  INTEGER, INTENT(IN) :: BLOCKS(nprocs), BLOCKSTART(nprocs)
  complex*16, INTENT(IN) :: V1(first:last)
  complex*16, INTENT(OUT) :: X1(1:MINPUT)
  INTEGER ::MPIERR=0

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  if (last-first+1.ne.blocks(myrank)) then
     print *, "ACK BLOCKS",myrank,nprocs,last,first,blocks(myrank);     stop
  endif

  CALL MPI_GATHERV(V1(first:), blocks(myrank), MPI_DOUBLE_COMPLEX, X1(:), BLOCKS(:), &
       BLOCKSTART(:)-1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, MPIERR)
  if (mpierr.ne.0) then
     print *, "GAtherv err",mpierr,myrank;          stop
  endif

  if (bcastflag) then
     CALL MPI_BCAST(X1(:),MINPUT,MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, MPIERR)
     if (mpierr.ne.0) then
        print *, "bcast err",mpierr,myrank;          stop
     endif
  endif

  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

END SUBROUTINE MYGATHERV_complex


SUBROUTINE MYGATHERV_real(V1,X1,MINPUT,FIRST,LAST,BLOCKS,BLOCKSTART,bcastflag)
  use mpimod
  IMPLICIT NONE
  logical, intent(in) :: bcastflag
  INTEGER, INTENT(IN) :: MINPUT, FIRST, LAST
  INTEGER, INTENT(IN) :: BLOCKS(nprocs), BLOCKSTART(nprocs)
  real*8, INTENT(IN) :: V1(first:last)
  real*8, INTENT(OUT) :: X1(1:MINPUT)
  INTEGER ::MPIERR=0

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  if (last-first+1.ne.blocks(myrank)) then
     print *, "ACK BLOCKS",myrank,nprocs,last,first,blocks(myrank);     stop
  endif

  CALL MPI_GATHERV(V1(first:), blocks(myrank), MPI_DOUBLE_PRECISION, X1(:), BLOCKS(:), &
       BLOCKSTART(:)-1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPIERR)
  if (mpierr.ne.0) then
     print *, "GAtherv err",mpierr,myrank;          stop
  endif

  if (bcastflag) then
     CALL MPI_BCAST(X1(:),MINPUT,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPIERR)
     if (mpierr.ne.0) then
        print *, "bcast err",mpierr,myrank;          stop
     endif
  endif

  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

END SUBROUTINE MYGATHERV_real



SUBROUTINE MYGATHERV_integer(V1,X1,MINPUT,FIRST,LAST,BLOCKS,BLOCKSTART,bcastflag)
  use mpimod
  IMPLICIT NONE
  logical, intent(in) :: bcastflag
  INTEGER, INTENT(IN) :: MINPUT, FIRST, LAST
  INTEGER, INTENT(IN) :: BLOCKS(nprocs), BLOCKSTART(nprocs)
  integer, INTENT(IN) :: V1(first:last)
  integer, INTENT(OUT) :: X1(1:MINPUT)
  INTEGER ::MPIERR=0

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  if (last-first+1.ne.blocks(myrank)) then
     print *, "ACK BLOCKS",myrank,nprocs,last,first,blocks(myrank);     stop
  endif

  CALL MPI_GATHERV(V1(first:), blocks(myrank), MPI_INTEGER, X1(:), BLOCKS(:), &
       BLOCKSTART(:)-1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPIERR)
  if (mpierr.ne.0) then
     print *, "GAtherv err",mpierr,myrank;          stop
  endif

  if (bcastflag) then
     CALL MPI_BCAST(X1(:),MINPUT,MPI_INTEGER, 0, MPI_COMM_WORLD, MPIERR)
     if (mpierr.ne.0) then
        print *, "bcast err",mpierr,myrank;          stop
     endif
  endif

  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

END SUBROUTINE MYGATHERV_integer




SUBROUTINE MYSCATTERV_complex(X1,V1,MINPUT,FIRST,LAST,BLOCKS,BLOCKSTART)
  use mpimod
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: MINPUT, FIRST, LAST
  INTEGER, INTENT(IN) :: BLOCKS(nprocs), BLOCKSTART(nprocs)
  complex*16, INTENT(OUT) :: V1(first:last)
  complex*16, INTENT(IN) :: X1(1:MINPUT)
  INTEGER ::MPIERR=0

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  if (last-first+1.ne.blocks(myrank)) then
     print *, "ACK BLOCKS",myrank,nprocs,last,first,blocks(myrank);     stop
  endif

  CALL MPI_SCATTERV(X1(:), blocks(:), blockstart(:)-1, MPI_DOUBLE_COMPLEX, V1(first:), BLOCKS(myrank), &
       MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, MPIERR)
  if (mpierr.ne.0) then
     print *, "Scatterv err",mpierr,myrank;          stop
  endif

  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

END SUBROUTINE MYSCATTERV_complex


SUBROUTINE MYSCATTERV_real(X1,V1,MINPUT,FIRST,LAST,BLOCKS,BLOCKSTART)
  use mpimod
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: MINPUT, FIRST, LAST
  INTEGER, INTENT(IN) :: BLOCKS(nprocs), BLOCKSTART(nprocs)
  real*8, INTENT(OUT) :: V1(first:last)
  real*8, INTENT(IN) :: X1(1:MINPUT)
  INTEGER ::MPIERR=0

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  if (last-first+1.ne.blocks(myrank)) then
     print *, "ACK BLOCKS",myrank,nprocs,last,first,blocks(myrank);     stop
  endif

  CALL MPI_SCATTERV(X1(:), blocks(:), blockstart(:)-1, MPI_DOUBLE_PRECISION, V1(first:), BLOCKS(myrank), &
       MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPIERR)
  if (mpierr.ne.0) then
     print *, "Scatterv err",mpierr,myrank;          stop
  endif

  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime

END SUBROUTINE MYSCATTERV_real


#else

subroutine mympiimax(input)
end subroutine mympiimax
subroutine mympiimin(input)
end subroutine mympiimin
subroutine mympimax(input)
end subroutine mympimax
subroutine mympimin(input)
end subroutine mympimin



SUBROUTINE MYGATHERV_complex(V1,X1,MINPUT,FIRST,LAST,BLOCKS,BLOCKSTART,bcastflag)
  IMPLICIT NONE
  logical, intent(in) :: bcastflag
  INTEGER, INTENT(IN) :: MINPUT, FIRST, LAST, BLOCKS(1), BLOCKSTART(1)
  complex*16, INTENT(IN) :: V1(first:last)
  complex*16, INTENT(OUT) :: X1(1:MINPUT)
  if (first-last.ne.minput-1) then
     print *, "AAUAUAUAUAUAUA"; stop
  endif
  X1(:)=V1(:)
END SUBROUTINE MYGATHERV_complex


SUBROUTINE MYGATHERV_real(V1,X1,MINPUT,FIRST,LAST,BLOCKS,BLOCKSTART,bcastflag)
  use mpimod
  IMPLICIT NONE
  logical, intent(in) :: bcastflag
  INTEGER, INTENT(IN) :: MINPUT, FIRST, LAST, BLOCKS(1), BLOCKSTART(1)
  real*8, INTENT(IN) :: V1(first:last)
  real*8, INTENT(OUT) :: X1(1:MINPUT)
  if (first-last.ne.minput-1) then
     print *, "AAUAUAUAUAUAUA"; stop
  endif
  X1(:)=V1(:)
END SUBROUTINE MYGATHERV_real


SUBROUTINE MYGATHERV_integer(V1,X1,MINPUT,FIRST,LAST,BLOCKS,BLOCKSTART,bcastflag)
  IMPLICIT NONE
  logical, intent(in) :: bcastflag
  INTEGER, INTENT(IN) :: MINPUT, FIRST, LAST, BLOCKS(1), BLOCKSTART(1)
  integer, INTENT(IN) :: V1(first:last)
  integer, INTENT(OUT) :: X1(1:MINPUT)
  if (first-last.ne.minput-1) then
     print *, "AAUAUAUAUAUAUA"; stop
  endif
  X1(:)=V1(:)
END SUBROUTINE MYGATHERV_integer




SUBROUTINE MYSCATTERV_complex(X1,V1,MINPUT,FIRST,LAST,BLOCKS,BLOCKSTART)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: MINPUT, FIRST, LAST,BLOCKS(1), BLOCKSTART(1)
  complex*16, INTENT(OUT) :: V1(first:last)
  complex*16, INTENT(IN) :: X1(1:MINPUT)
  if (first-last.ne.minput-1) then
     print *, "AAUAUAUAUAUAUA"; stop
  endif
  V1(:)=X1(:)
END SUBROUTINE MYSCATTERV_complex

SUBROUTINE MYSCATTERV_real(X1,V1,MINPUT,FIRST,LAST,BLOCKS,BLOCKSTART)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: MINPUT, FIRST, LAST,BLOCKS(1), BLOCKSTART(1)
  real*8, INTENT(OUT) :: V1(first:last)
  real*8, INTENT(IN) :: X1(1:MINPUT)
  if (first-last.ne.minput-1) then
     print *, "AAUAUAUAUAUAUA"; stop
  endif
  V1(:)=X1(:)
END SUBROUTINE MYSCATTERV_real




subroutine mympirealreduce(input,isize)
end subroutine
subroutine mympirealreduceone(input)
end subroutine
subroutine mympireduceone(input)
end subroutine mympireduceone
subroutine mympireduce(input, isize)
end subroutine mympireduce
subroutine mympiireduce(input, isize)
end subroutine mympiireduce
subroutine mympibcast(input, source, isize)
end subroutine mympibcast
subroutine mympirealbcast(input, source, isize)
end subroutine mympirealbcast
subroutine mympiibcast(input, source, isize)
end subroutine mympiibcast

subroutine mympirealbcastone(input, source)
end subroutine mympirealbcastone
subroutine mympiibcastone(input, source)
end subroutine mympiibcastone

subroutine mympilogbcast(input, source, isize)
end subroutine mympilogbcast


subroutine mympiireduceone(input)
end subroutine mympiireduceone

subroutine mpiorbgather(orbvector,insize)
end subroutine mpiorbgather





subroutine mpistart
  use mpimod
  use parameters
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
  use parameters  
  implicit none
  if (stdoutflag==0) then
     open(mpifileptr,file=mpioutfile,status="old", position="append")
  endif
end subroutine openfile

subroutine closefile()
  use mpimod
  use parameters
  implicit none
  if (stdoutflag==0)  then
     close(mpifileptr)
  endif
end subroutine closefile






subroutine mpiallgather(inout,totsize,blocksizes,maxblocksize)
  use mpimod
  use parameters
  implicit none

  integer, intent(in) :: totsize,blocksizes(nprocs),maxblocksize
  DATATYPE :: inout(totsize)

#ifdef REALGO
 call mpiallgather_real(inout,totsize,blocksizes,maxblocksize)
#else
 call mpiallgather_complex(inout,totsize,blocksizes,maxblocksize)
#endif

end subroutine mpiallgather


subroutine mpiallgather_real(inout,totsize,blocksizes,maxblocksize)
  use mpimod
  use parameters
  implicit none

  integer, intent(in) :: totsize,blocksizes(nprocs),maxblocksize
  real*8 :: inout(totsize)
  real*8, allocatable :: intemp(:,:),outtemp(:,:)
  integer :: i,icount,ierr,blockstart(nprocs+1)

#ifndef MPIFLAG
  return
#else

if (nprocs.eq.1) then
   return
endif

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  icount=0
  blockstart(1)=1
  do i=1,nprocs
     blockstart(i+1)=blockstart(i)+blocksizes(i)
     if (blocksizes(i).gt.maxblocksize) then
        OFLWR "BBB ERR"; CFLST
     endif
     icount=icount + blocksizes(i)
  enddo
  if (icount.ne.totsize) then
     OFLWR "ALLgather ERR", icount,totsize; 
     WRFL blocksizes;
     CFLST
  endif

  allocate(intemp(maxblocksize,nprocs)); intemp(:,:) = 0d0
  allocate(outtemp(maxblocksize,nprocs)); outtemp(:,:) = 0d0

  intemp=0d0

!   intemp(1:blocksizes(myrank),myrank)=inout(blockstart(myrank):blockstart(myrank+1)-1)

  do i=1,  nprocs
     intemp(1:blocksizes(myrank),i)=inout(blockstart(myrank):blockstart(myrank+1)-1)
  enddo


  outtemp=0d0

  call mpi_allgather(intemp(:,:),maxblocksize,MPI_DOUBLE_PRECISION,outtemp(:,:),maxblocksize,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  if (ierr.ne.0) then
     OFLWR "ALLGATHER ERR ", ierr; CFLST
  endif

  inout(:)=0d0
  do i=1,nprocs
     inout(blockstart(i):blockstart(i+1)-1)=outtemp(1:blocksizes(i),i)
  enddo

  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
  deallocate(intemp,outtemp)

#endif

end subroutine mpiallgather_real


subroutine mpiallgather_complex(inout,totsize,blocksizes,maxblocksize)
  use mpimod
  use parameters
  implicit none

  integer, intent(in) :: totsize,blocksizes(nprocs),maxblocksize
  complex*16 :: inout(totsize)
  complex*16, allocatable :: intemp(:,:),outtemp(:,:)
  integer :: i,icount,ierr,blockstart(nprocs+1)

#ifndef MPIFLAG
  return
#else

if (nprocs.eq.1) then
   return
endif

  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime

  icount=0
  blockstart(1)=1
  do i=1,nprocs
     blockstart(i+1)=blockstart(i)+blocksizes(i)
     if (blocksizes(i).gt.maxblocksize) then
        OFLWR "BBB ERR"; CFLST
     endif
     icount=icount + blocksizes(i)
  enddo
  if (icount.ne.totsize) then
     OFLWR "ALLgather ERR", icount,totsize; 
     WRFL blocksizes;
     CFLST
  endif

  allocate(intemp(maxblocksize,nprocs)); intemp(:,:) = 0d0
  allocate(outtemp(maxblocksize,nprocs)); outtemp(:,:) = 0d0

  intemp=0d0

!   intemp(1:blocksizes(myrank),myrank)=inout(blockstart(myrank):blockstart(myrank+1)-1)

  do i=1,  nprocs
     intemp(1:blocksizes(myrank),i)=inout(blockstart(myrank):blockstart(myrank+1)-1)
  enddo


  outtemp=0d0

  call mpi_allgather(intemp(:,:),maxblocksize,MPI_DOUBLE_COMPLEX,outtemp(:,:),maxblocksize,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)

  if (ierr.ne.0) then
     OFLWR "C ALLGATHER ERR ", ierr; CFLST
  endif

  inout(:)=0d0
  do i=1,nprocs
     inout(blockstart(i):blockstart(i+1)-1)=outtemp(1:blocksizes(i),i)
  enddo

  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
  deallocate(intemp,outtemp)

#endif

end subroutine mpiallgather_complex


subroutine mympialltoall(input, output, count)
  use fileptrmod
  use mpimod
  implicit none
  integer :: ierr, count
  DATATYPE :: input(count,nprocs), output(count,nprocs)

#ifndef MPIFLAG
  output(:,:)=input(:,:)
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

end subroutine mympialltoall


!! probably just the same -- transferred above from boxes
!!subroutine mpiallgather(inout,totsize,blocksizes,maxblocksize)
!!  use mpimod
!!  use parameters
!!  implicit none
!!
!!  integer, intent(in) :: totsize,blocksizes(nprocs),maxblocksize
!!  DATATYPE :: inout(totsize)
!!  DATATYPE, allocatable :: intemp(:,:),outtemp(:,:)
!!  integer :: i,icount,ierr,blockstart(nprocs+1)
!!
!!#ifndef MPIFLAG
!!  return
!!#else
!!
!!if (nprocs.eq.1) then
!!   return
!!endif
!!
!!  call system_clock(mpiatime);  nonmpitime=nonmpitime+mpiatime-mpibtime
!!
!!  icount=0
!!  blockstart(1)=1
!!  do i=1,nprocs
!!     blockstart(i+1)=blockstart(i)+blocksizes(i)
!!     if (blocksizes(i).gt.maxblocksize) then
!!        OFLWR "BBB ERR"; CFLST
!!     endif
!!     icount=icount + blocksizes(i)
!!  enddo
!!  if (icount.ne.totsize) then
!!     OFLWR "ALLgather ERR", icount,totsize; 
!!     WRFL blocksizes;
!!     CFLST
!!  endif
!!
!!  allocate(intemp(maxblocksize,nprocs)); intemp(:,:) = 0d0
!!  allocate(outtemp(maxblocksize,nprocs)); outtemp(:,:) = 0d0
!!
!!  intemp=0d0
!!
!!!   intemp(1:blocksizes(myrank),myrank)=inout(blockstart(myrank):blockstart(myrank+1)-1)
!!
!!  do i=1,  nprocs
!!     intemp(1:blocksizes(myrank),i)=inout(blockstart(myrank):blockstart(myrank+1)-1)
!!  enddo
!!
!!
!!  outtemp=0d0
!!
!!#ifdef REALGO
!!  call mpi_allgather(intemp(:,:),maxblocksize,MPI_DOUBLE_PRECISION,outtemp(:,:),maxblocksize,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
!!#else
!!  call mpi_allgather(intemp(:,:),maxblocksize,MPI_DOUBLE_COMPLEX,outtemp(:,:),maxblocksize,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
!!#endif
!!
!!
!!  if (ierr.ne.0) then
!!     OFLWR "ALLGATHER ERR ", ierr; CFLST
!!  endif
!!
!!
!!  inout(:)=0d0
!!  do i=1,nprocs
!!     inout(blockstart(i):blockstart(i+1)-1)=outtemp(1:blocksizes(i),i)
!!  enddo
!!
!!
!!
!!  call system_clock(mpibtime);  mpitime=mpitime+mpibtime-mpiatime
!!  deallocate(intemp,outtemp)
!!
!!
!!#endif
!!
!!end subroutine mpiallgather
!!


