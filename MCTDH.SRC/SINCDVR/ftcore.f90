
!! myrank is indexed 1:nprocs

module bothblockmod
  integer :: nprocs=-1,dim=-1,myrank=-1
  integer, allocatable :: mpiblocks(:),mpiblockend(:),mpiblockstart(:)
end module bothblockmod


#ifdef FFTWFLAG

recursive subroutine myzfft1d(in,out,dim,howmany)
  implicit none
  integer, intent(in) :: dim,howmany
  complex*16, intent(in) :: in(dim,howmany)
  complex*16, intent(out) :: out(dim,howmany)
  call fftw1dfftsub(in,out,dim,howmany)
end subroutine myzfft1d


recursive subroutine myzfft3d(in,out,indim,howmany)
  use bothblockmod
  implicit none
  integer :: indim,howmany
  complex*16, intent(in) :: in(dim,dim,dim,howmany)
  complex*16, intent(out) :: out(dim,dim,dim,howmany)
  if (dim.ne.indim) then
     print *, "WRONG INIT FFTW3D",dim,indim;stop
  endif
  call fftw3dfftsub(in,out,howmany)
end subroutine myzfft3d


#else

recursive subroutine myzfft1d(in,out,dim,howmany)
  implicit none
  integer, intent(in) :: dim,howmany
  integer :: k
  complex*16, intent(in) :: in(dim,howmany)
  complex*16, intent(out) :: out(dim,howmany)
  complex*16 :: wsave(4*dim+15,howmany)   ! MAKE BIGGER IF SEGFAULT... iffy
  out(:,:)=in(:,:)
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(in,out,dim,howmany,wsave)
!$OMP DO SCHEDULE(STATIC)
  do k=1,howmany
     call zffti(dim,wsave(:,k))
     call zfftf(dim,out(:,k),wsave(:,k))
  enddo
!$OMP END DO
!$OMP END PARALLEL
end subroutine myzfft1d


subroutine myzfft3d(in,out,indim,howmany)
  use bothblockmod
  implicit none
  integer :: indim,howmany
  complex*16, intent(in) :: in(dim,dim,dim,howmany)
  complex*16, intent(out) :: out(dim,dim,dim,howmany)
  complex*16 :: work(dim,dim,dim,howmany)  !! AUTOMATIC
  integer :: ii,i
  if (dim.ne.indim) then
     print *, "WRONG INIT",dim,indim;stop
  endif
  out(:,:,:,:)=in(:,:,:,:)
  do ii=1,3
     call myzfft3d_oneblock(out,work,dim,howmany)
     do i=1,dim
        out(:,:,i,:)=work(i,:,:,:)
     enddo
  enddo
end subroutine myzfft3d

#endif


module littlestartmod
 integer :: mystart=(-1),mysize=(-1)
end module littlestartmod


subroutine setblock(innprocs,inmyrank,indims)
  use bothblockmod
  use littlestartmod
  implicit none
  integer :: innprocs,inmyrank,indims(3),i
  if (dim.ne.-1) then
     print *, "CALLME ONCE ONLY"; stop
  endif
  nprocs=innprocs
  myrank=inmyrank
  if (myrank.eq.0) then
     print *, "WRONG CONVENTION."; stop
  endif
  if (indims(2).ne.indims(3).or.indims(1).ne.indims(2)) then
     print *, "Only all dims equal for now", indims; stop
  endif
  dim=indims(1)
  allocate(mpiblocks(nprocs),mpiblockend(nprocs),mpiblockstart(nprocs))
  mpiblockstart(1)=1
  do i=1,nprocs
     mpiblockend(i)=(i*dim/nprocs)*dim**2
     if (i.lt.nprocs) then
        mpiblockstart(i+1)=mpiblockend(i)+1
     endif
  enddo
  do i=1,nprocs
     mpiblocks(i)=mpiblockend(i)-mpiblockstart(i)+1
  enddo
  mystart=mpiblockstart(myrank)
  mysize=mpiblocks(myrank)
end subroutine setblock


subroutine unsetblock()
  use bothblockmod
  implicit none
  if (dim.eq.(-1)) then
     if (myrank.eq.1) then 
        print *, "HMM UNSET ME BUT WHY?";
     endif
     return !!stop
  endif
!!  nprocs=(-1);  myrank=(-1);  
  dim=(-1)
  deallocate(mpiblocks,mpiblockend,mpiblockstart)
end subroutine unsetblock


recursive subroutine myzfft3d_oneblock(in,out,insize,howmany)
  use bothblockmod
  implicit none
  integer :: insize,howmany
  complex*16, intent(in) :: in(dim,dim,insize,howmany)
  complex*16, intent(out) :: out(dim,dim,insize,howmany)
  if (dim.eq.-1) then
     print *, "NEED TO INITIALIZE FFTBLOCK";stop
  endif
  call myzfft1d(in(:,:,:,:),out(:,:,:,:),dim,insize*dim*howmany)
end subroutine myzfft3d_oneblock


#ifdef FFTWFLAG


recursive subroutine fftw3dfftsub(in,out,howmany)
  use, intrinsic :: iso_c_binding
  use bothblockmod
  implicit none
  include "fftw3.f03"
  integer :: howmany,ii
  type(C_PTR),save :: plan
  integer,save :: icalled=0
  complex*16 :: in(dim,dim,dim,howmany),out(dim,dim,dim,howmany)
  if (icalled.eq.0) then
     plan=fftw_plan_dft_3d(dim,dim,dim,in(:,:,:,1),out(:,:,:,1),FFTW_FORWARD,FFTW_EXHAUSTIVE) 
  endif
  icalled=1
  do ii=1,howmany
     call fftw_execute_dft(plan, in(:,:,:,ii),out(:,:,:,ii))
  enddo

!!$  call fftw_destroy_plan(plan)

end subroutine fftw3dfftsub


#ifdef MPIFFTW

 not finished with howmany!  

subroutine fftw3dfftsub_mpi(in,out,indim,insize,howmany)
  use, intrinsic :: iso_c_binding
  use bothblockmod
  use littlestartmod
  include 'fftw3-mpi.f03'
  integer :: indim,insize,howmany
  integer, save :: icalled=0
  integer(C_INTPTR_T) :: alloc_local,LL,MM,SS
  complex*16,intent(in) :: in(dim,dim,mysize/dim**2,howmany)
  complex*16,intent(out) :: out(dim,dim,mysize/dim**2,howmany)
  type(C_PTR),save :: plan, cdata
  complex(C_DOUBLE_COMPLEX), pointer :: data(:,:,:)
  if (myrank.eq.1) then
     print *, "GO FFTW MPI SUB."
  endif

  LL=dim;MM=mysize/dim**2;SS=(mystart-1)/dim**2

  if (indim.ne.dim.or.insize.ne.mysize/dim**2) then
     print *, "AUGH FFTW MPI FFT ",indim,dim,insize,mysize; stop
  endif

!print *, "LMS ", LL,MM,SS
  alloc_local = fftw_mpi_local_size_3d(LL,LL,LL, MPI_COMM_WORLD, MM, SS)
!print *, "LMSnow ", LL,MM,SS

  cdata = fftw_alloc_complex(alloc_local)
  call c_f_pointer(cdata, data, [LL,LL,MM])
     
  if (icalled.eq.0) then

     !   create MPI plan for in-place forward DFT (note dimension reversal)

     plan = fftw_mpi_plan_dft_3d(LL,LL,LL, data, data, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_EXHAUSTIVE)
  endif
  icalled=1

     ! initialize data to some function my_function(i,j)

  data(:, :,:) = in(:,:,:)
  ! compute transform (as many times as desired)

  call fftw_mpi_execute_dft(plan, data, data)

!!  call fftw_execute_dft(plan, data, data)

  out(:,:,:)=data(:,:,:)

  if (myrank.eq.1) then
     print *, "ok done fftw3dfftsub_mpi"
  endif
end subroutine fftw3dfftsub_mpi

#endif

recursive subroutine fftw1dfftsub(in,out,dim,howmany)
  use, intrinsic :: iso_c_binding
  implicit none
  include "fftw3.f03"
  integer, intent(in) :: dim,howmany
  integer, parameter :: maxplans=3
  type(C_PTR),save :: plans(maxplans)
  integer, save :: plandims(maxplans)=-999, planhowmany(maxplans)=-999
  integer,save :: icalleds(maxplans)=0, numplans=0
  complex*16 :: in(dim,howmany),out(dim,howmany)
  integer :: ostride,istride,onembed(1),inembed(1),idist,odist, dims(1),iplan,thisplan

  inembed(1)=dim; onembed(1)=dim; idist=dim; odist=dim; istride=1; ostride=1; dims(1)=dim

  if (numplans.eq.0) then
     numplans=1
     thisplan=1
     plandims(thisplan)=dim; planhowmany(thisplan)=howmany
  else
     thisplan= -99
     do iplan=1,numplans
        if (plandims(iplan).eq.dim.and.planhowmany(iplan).eq.howmany) then
           if (icalleds(iplan).eq.0) then
              print *, "ERROR, plan not done ",iplan,dim,howmany; call mpistop()
           endif
           thisplan=iplan
           exit
        endif
     enddo
     if (thisplan.eq.-99) then
        if (numplans.eq.maxplans) then
           print *,  "all plans taken!", maxplans; call mpistop()
        endif
        numplans=numplans+1
        thisplan=numplans
        plandims(thisplan)=dim; planhowmany(thisplan)=howmany
     endif
  endif
  if (icalleds(thisplan).eq.0) then
     plans(thisplan) = fftw_plan_many_dft(1,dims,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,FFTW_FORWARD,FFTW_EXHAUSTIVE) 
  endif
  icalleds(thisplan)=1    

  call fftw_execute_dft(plans(thisplan), in,out)

!!$  call fftw_destroy_plan(plan)
end subroutine fftw1dfftsub

#endif



#ifdef MPIFLAG


!! times(1) = transpose   times(2) = mpi  times(3) = copy
!! This is 3D specific
!!   (123) -> (312)

module mytransposemod
contains
  recursive subroutine mytranspose(in,out,blocksize,howmany,times)
  use bothblockmod
  implicit none
  integer,intent(in) :: blocksize,howmany
  integer,intent(inout) :: times(3)
  complex*16,intent(in) :: in(nprocs*blocksize,nprocs*blocksize,blocksize,howmany)
  complex*16 :: inchop(nprocs*blocksize,blocksize,blocksize,howmany)  !!AUTOMATIC
  complex*16,intent(out) :: out(nprocs*blocksize,nprocs*blocksize,blocksize,howmany)
  integer :: atime,btime,i,count,ii,totsize,iproc
  complex*16 :: intranspose(blocksize,nprocs*blocksize,blocksize,howmany,nprocs)  !!AUTOMATIC
  complex*16 :: outtemp(blocksize,nprocs*blocksize,blocksize,howmany,nprocs)      !!AUTOMATIC

  totsize=blocksize*nprocs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!    (123)->(312)    !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call myclock(atime)
  
!!$  if (myrank.eq.1) print *, "TEMP NO OMP TRANSPOSE"

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,i,iproc)
  do iproc=1,nprocs
     inchop(:,:,:,:)=in(:,(iproc-1)*blocksize+1:iproc*blocksize,:,:)
!$OMP DO SCHEDULE(DYNAMIC)
     do ii=1,howmany
        do i=1,blocksize
           intranspose(:,:,i,ii,iproc)=TRANSPOSE(inchop(:,i,:,ii))
        enddo
     enddo
!$OMP END DO
  enddo
! (implied barrier at end)
!$OMP END PARALLEL

 call myclock(btime); times(1)=times(1)+btime-atime; atime=btime
  
  count=blocksize**2 * totsize * howmany
  
  call mympialltoall_complex(intranspose,outtemp,count)
  
 call myclock(btime); times(2)=times(2)+btime-atime; atime=btime
  
  do iproc=1,nprocs
     out((iproc-1)*blocksize+1:iproc*blocksize,:,:,:)=outtemp(:,:,:,:,iproc)
  enddo

  call myclock(btime); times(3)=times(3)+btime-atime;
  
end subroutine mytranspose
end module  
  

recursive subroutine myzfft3d_mpiwrap(in,out,indim,howmany)
  use littlestartmod
  use bothblockmod
  implicit none
  integer :: indim,nulltimes(10),howmany,ii
  complex*16, intent(in) :: in(dim**3,howmany)
  complex*16, intent(out) :: out(dim**3,howmany)

  if (dim.ne.indim) then
     print *, "WRONG INIT",dim,indim;stop
  endif

!!! ii=1,howmany loops are temporary(?) kloodge... 
!!!   this subroutine _mpiwrap (called when orbparflag=.false.) 
!!    should not really be used... hmm
!!! (3pm 04-13-2015 BUGFIX)

#ifdef MPIFFTW
  do ii=1,howmany
     call fftw3dfftsub_mpi(in(mystart,ii),out(mystart,ii),indim,mysize/dim**2,1)
  enddo
  return
#else
  do ii=1,howmany
     call myzfft3d_par(in(mystart,ii),out(mystart,ii),indim,mysize/dim**2,nulltimes,1)
  enddo

#define TENxxTEST
#ifdef TENTEST
  write(*,'(4E20.10,A10,I5)') DOT_PRODUCT(&
       in(mystart:mystart+mysize-1),&
       in(mystart:mystart+mysize-1)),&
       DOT_PRODUCT(&
       out(mystart:mystart+mysize-1),&
       out(mystart:mystart+mysize-1)), &
       "DDOTPAR",myrank
  call mybarrier()
  stop
#endif

  do ii=1,howmany
     call mygatherv_complex(out(mpiblockstart(myrank),ii),out(:,ii),dim**3,&
          mpiblockstart(myrank),&
          mpiblockend(myrank),&
          mpiblocks(:),&
          mpiblockstart(:),.true.)
  enddo
#endif
end subroutine myzfft3d_mpiwrap


!! adds to times

!!! times(1) = zero   times(2)=fourier
!!! from mytranspose times(3) = transpose   times(4) = mpi  times(5) = copy
  
recursive subroutine myzfft3d_par(in,out,indim,inblockdim,times,howmany)
  use bothblockmod
  use littlestartmod
  use mytransposemod
  implicit none
  integer, intent(in) :: indim,inblockdim,howmany
  complex*16, intent(in) :: in(dim,dim,mysize/dim**2,howmany)
  complex*16, intent(out) :: out(dim,dim,mysize/dim**2,howmany)
  integer, intent(inout) :: times(8)
  integer :: ii,atime,btime
  complex*16 :: mywork(dim,dim,mysize/dim**2,howmany),tempout(dim,dim,mysize/dim**2,howmany)
!!$  complex*16, allocatable, save :: mywork(:,:,:,:), tempout(:,:,:,:)
!!$  integer, save :: savehowmany = -999

  call myclock(atime)

!!$  if (dims(1)*dims(2)*(mysize/dims(1)/dims(2)).ne.mysize) then
  if (dim**2*(mysize/dim**2).ne.mysize) then
     print *, "WTF!!! 5656578",dim,mysize; stop
  endif
  if (dim.ne.indim) then
     print *, "WRONG INIT",dim,indim;stop
  endif
  if (dim**2*inblockdim.ne.mysize) then
     print *, "WRONG BLOCK",dim,indim," ",mysize,inblockdim,dim;stop
  endif
  if (mysize.ne.mpiblocks(myrank)) then
     print *, "MYSIZE/blocks disagree",mysize,mpiblocks(myrank),myrank,nprocs;stop
  endif

!!$  if (savehowmany.eq.-999) then
!!$     allocate(mywork(dim,dim,mysize/dim**2,howmany),tempout(dim,dim,mysize/dim**2,howmany))
!!$     mywork=0; tempout=0
!!$  else if (savehowmany.ne.howmany) then
!!$     if (myrank.eq.1) print *, "fDEALLOCATING",myrank,howmany,savehowmany 
!!$     deallocate(mywork,tempout)
!!$     allocate(mywork(dim,dim,mysize/dim**2,howmany),tempout(dim,dim,mysize/dim**2,howmany))
!!$     mywork=0; tempout=0
!!$  endif
!!$  savehowmany=howmany


#ifdef MPIFFTW
  call fftw3dfftsub_mpi(in,out,indim,mysize/dim**2,howmany)
  return
#else
  tempout(:,:,:,:)=in(:,:,:,:)
  call myclock(btime); times(1)=times(1)+btime-atime;

  do ii=1,3
     call myclock(atime)
#ifndef TENTEST
     call myzfft3d_oneblock( tempout, mywork, mysize/dim**2,howmany)
#else
     mywork(:,:,:,:)=tempout(:,:,:,:)*10
#endif
     call myclock(btime); times(2)=times(2)+btime-atime; atime=btime
     
!!! from mytranspose times(3) = transpose   times(4) = mpi  times(5) = copy
     call mytranspose(&
          mywork,  &
          tempout,  &
          mysize/dim**2, &
          howmany,times(3:))
  enddo
  out(:,:,:,:)=tempout(:,:,:,:)
#endif

end subroutine myzfft3d_par

#endif




