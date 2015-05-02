
!! myrank is 1:nprocs

#ifdef FFTWFLAG

recursive subroutine myzfft1d(in,out,dim,howmany)
  use, intrinsic :: iso_c_binding
  implicit none
  include "fftw3.f03"
  integer, intent(in) :: dim,howmany
  complex*16 :: in(dim,howmany)    !! cannot be declared intent(in)...hmmm...
  complex*16, intent(out) :: out(dim,howmany)
  integer, parameter :: maxplans=3
  type(C_PTR),save :: plans(maxplans)
  integer, save :: plandims(maxplans)=-999, planhowmany(maxplans)=-999
  integer,save :: icalleds(maxplans)=0, numplans=0
  integer :: ostride,istride,onembed(1),inembed(1),idist,odist, dims(1),iplan,thisplan
  integer :: myrank,nprocs
  call getmyranknprocs(myrank,nprocs)

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
     if (myrank.eq.1) then
        print *, "       Making a 1D FFT plan ", thisplan, dims, howmany
     endif
     plans(thisplan) = fftw_plan_many_dft(1,dims,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,FFTW_FORWARD,FFTW_EXHAUSTIVE) 
  endif
  icalleds(thisplan)=1    

  call fftw_execute_dft(plans(thisplan), in,out)

end subroutine myzfft1d


recursive subroutine myzfft3d(in,out,dim1,dim2,dim3,howmany)
  use, intrinsic :: iso_c_binding
  implicit none
  include "fftw3.f03"
  integer, intent(in) :: dim1,dim2,dim3,howmany
  complex*16 :: in(dim1,dim2,dim3,howmany)  !! cannot be declared intent(in)...hmmm...
  complex*16, intent(out) :: out(dim1,dim2,dim3,howmany)
  integer, parameter :: maxplans=3
  type(C_PTR),save :: plans(maxplans)
  integer, save :: plandims(3,maxplans)=-999, planhowmany(maxplans)=-999
  integer,save :: icalleds(maxplans)=0, numplans=0
  integer :: ostride,istride,onembed(3),inembed(3),idist,odist, dims(3),iplan,thisplan
  integer :: myrank,nprocs
  call getmyranknprocs(myrank,nprocs)

  dims(:)=(/dim3,dim2,dim1/)
  inembed(:)=dims(:); onembed(:)=dims(:); idist=dim1*dim2*dim3; odist=dim1*dim2*dim3; istride=1; ostride=1; 

  if (numplans.eq.0) then
     numplans=1
     thisplan=1
     plandims(:,thisplan)=dims(:); planhowmany(thisplan)=howmany
  else
     thisplan= -99
     do iplan=1,numplans
        if (plandims(1,iplan).eq.dims(1).and.&
             plandims(2,iplan).eq.dims(2).and.&
             plandims(3,iplan).eq.dims(3).and.&
             planhowmany(iplan).eq.howmany) then
           if (icalleds(iplan).eq.0) then
              print *, "ERROR, plan not done ",iplan,dims(:),howmany; call mpistop()
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
        plandims(:,thisplan)=dims(:); planhowmany(thisplan)=howmany
     endif
  endif
  if (icalleds(thisplan).eq.0) then
     if (myrank.eq.1) then
        print *, "       Making a 3D fft plan ", thisplan, dims, howmany
     endif
     plans(thisplan) = fftw_plan_many_dft(3,dims,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,FFTW_FORWARD,FFTW_EXHAUSTIVE) 
  endif
  icalleds(thisplan)=1    

  call fftw_execute_dft(plans(thisplan), in,out)

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


recursive subroutine myzfft3d(in,out,dim1,dim2,dim3,howmany)
  implicit none
  integer :: dim1,dim2,dim3,howmany
  complex*16, intent(in) :: in(dim1,dim2,dim3,howmany)
  complex*16, intent(out) :: out(dim1,dim2,dim3,howmany)
  out(:,:,:,:)=in(:,:,:,:)
  call fftblock_withtranspose(out,dim1,dim2,dim3,howmany)
  call fftblock_withtranspose(out,dim2,dim3,dim1,howmany)
  call fftblock_withtranspose(out,dim3,dim1,dim2,howmany)
end subroutine myzfft3d


recursive subroutine fftblock_withtranspose(inout,dim1,dim2,dim3,howmany)
  implicit none
  integer :: dim1,dim2,dim3,howmany
!!!!  is dimensioned (dim1,dim2,dim3) on input. !!!!
  complex*16,intent(inout) :: inout(dim2,dim3,dim1,howmany) 
  complex*16 :: work1(dim1,dim2,dim3,howmany)  !! AUTOMATIC
  integer :: i
  call myzfft1d(inout,work1,dim1,dim2*dim3*howmany)
  do i=1,dim1
     inout(:,:,i,:)=work1(i,:,:,:)
  enddo
end subroutine fftblock_withtranspose

#endif


#ifdef MPIFLAG

!! times(1) = transpose   times(2) = mpi  times(3) = copy
!! This is 3D specific
!!   (123) -> (312)

module mytransposemod
contains
  recursive subroutine mytranspose(in,out,blocksize,howmany,times,nprocs)
  implicit none
  integer,intent(in) :: blocksize,howmany,nprocs
  integer,intent(inout) :: times(3)
  complex*16,intent(in) :: in(nprocs*blocksize,nprocs*blocksize,blocksize,howmany)
  complex*16,intent(out) :: out(nprocs*blocksize,nprocs*blocksize,blocksize,howmany)
  integer :: atime,btime,i,count,ii,totsize,iproc,checkprocs,myrank,j
#define AUTOMATICBLOCKS
#ifdef AUTOMATICBLOCKS
  complex*16 :: intranspose(blocksize,nprocs*blocksize,blocksize,howmany,nprocs)  !!AUTOMATIC
  complex*16 :: outtemp(blocksize,nprocs*blocksize,blocksize,howmany,nprocs)      !!AUTOMATIC
#else
  complex*16, allocatable :: intranspose(:,:,:,:,:), outtemp(:,:,:,:,:)
#endif
#define AUTOMATICTRANS
#ifdef AUTOMATICTRANS
  complex*16 :: outone(blocksize,nprocs*blocksize,blocksize,nprocs)               !!AUTOMATIC
  complex*16 :: inchop(nprocs*blocksize,blocksize,blocksize,howmany)              !!AUTOMATIC
#else
!!
!! if outone and inchop are private to threads, then there is no need for barriers below.
!! 
  complex*16, allocatable, save :: outone(:,:,:,:),inchop(:,:,:,:)
!$OMP THREADPRIVATE(outone,inchop)
#endif

#ifndef AUTOMATICBLOCKS
  allocate( intranspose(blocksize,nprocs*blocksize,blocksize,howmany,nprocs), &
       outtemp(blocksize,nprocs*blocksize,blocksize,howmany,nprocs))
#endif

  call getmyranknprocs(myrank,checkprocs)
  if (nprocs.ne.checkprocs) then
     print *, "ACK CHECKPROCS",nprocs,checkprocs; call mpistop()
  endif

  totsize=blocksize*nprocs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!    (123)->(312)    !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call myclock(atime)

#ifdef AUTOMATICTRANS  
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,i)    !! IPROC IS SHARED (GOES WITH BARRIER, INCHOP SHARED)
#else
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,i,iproc)    !! IPROC IS private (no BARRIER, INCHOP threadprivate)
  allocate(inchop(nprocs*blocksize,blocksize,blocksize,howmany))
#endif

!$OMP MASTER
  intranspose(:,:,:,:,:)=0d0
!$OMP END MASTER
!! *** OMP BARRIER *** !!  make sure master touches first
!$OMP BARRIER

#ifndef AUTOMATICTRANS
!$OMP DO SCHEDULE(STATIC)
#endif
  do iproc=1,nprocs
     inchop(:,:,:,:)=in(:,(iproc-1)*blocksize+1:iproc*blocksize,:,:)
#ifdef AUTOMATICTRANS
!$OMP DO SCHEDULE(STATIC)
#endif
     do ii=1,howmany
        do i=1,blocksize
           intranspose(:,:,i,ii,iproc)=TRANSPOSE(inchop(:,i,:,ii))
        enddo
     enddo
#ifdef AUTOMATICTRANS
!$OMP END DO
!! *** OMP BARRIER *** !!   if inchop & iproc are shared
!$OMP BARRIER
  enddo
#else
  enddo
!$OMP END DO
  deallocate(inchop)
#endif
! (implied barrier at end)
!$OMP END PARALLEL

!! UMM does this do anything, not sure
!!
!$OcccMP PARALLEL DEFAULT(SHARED)
!$OcccMP MASTER
  call myclock(btime); times(1)=times(1)+btime-atime; atime=btime

  outtemp(:,:,:,:,:)=0d0
  
  count=blocksize**2 * totsize * howmany
  
  call mympialltoall_complex(intranspose,outtemp,count)
  call myclock(btime); times(2)=times(2)+btime-atime; atime=btime
!$OccccMP END MASTER
!! (barrier at end)
!$OccccMP END PARALLEL


!!  complex*16,intent(out) :: out(nprocs*blocksize,nprocs*blocksize,blocksize,howmany)
!!  complex*16 :: outtemp(blocksize,nprocs*blocksize,blocksize,howmany,nprocs)     
!!  complex*16 :: outone(blocksize,nprocs*blocksize,blocksize,nprocs)

!!$  OLD WAY - NEVER OMP'D WELL OBVIOUSLY
!!$  do iproc=1,nprocs
!!$     out((iproc-1)*blocksize+1:iproc*blocksize,:,:,:)=outtemp(:,:,:,:,iproc)
!!$  enddo

#ifdef AUTOMATICTRANS

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)  !! ii IS SHARED (OUTTEMP IS SHARED; BARRIER)
  do ii=1,howmany
     outone(:,:,:,:)=outtemp(:,:,:,ii,:)
!$OMP DO SCHEDULE(STATIC)
#else

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii)  !! ii IS PRIVATE (OUTTEMP IS THREADPRIVATE; NO BARRIER)
  allocate(outone(blocksize,nprocs*blocksize,blocksize,nprocs))
!$OMP DO SCHEDULE(STATIC)
  do ii=1,howmany
     outone(:,:,:,:)=outtemp(:,:,:,ii,:)

#endif
     do i=1,blocksize
        do j=1,nprocs*blocksize
           out(:,j,i,ii)=RESHAPE(outone(:,j,i,:),(/nprocs*blocksize/))
        enddo
     enddo

#ifdef AUTOMATICTRANS
!$OMP END DO
!! *** OMP BARRIER *** !!   if outone & ii are shared
!$OMP BARRIER
  enddo
#else
  enddo   !! no barrier because outone is threadprivate; ii is private
!$OMP END DO
  deallocate(outone)
#endif

! (implied barrier at end)
!$OMP END PARALLEL

  call myclock(btime); times(3)=times(3)+btime-atime;

#ifndef AUTOMATICBLOCKS
  deallocate(intranspose, outtemp)
#endif
  
end subroutine mytranspose
end module  
  

subroutine checkdivisible(number,divisor)
  implicit none
  integer :: number,divisor
  if ((number/divisor)*divisor.ne.number) then
     print *, "ACK NOT DIVISIBLE",number,divisor; call mpistop()
  endif
end subroutine checkdivisible


recursive subroutine myzfft3d_mpiwrap_forward(in,out,dim,howmany,placeopt)
  implicit none
  integer, intent(in) :: dim,howmany,placeopt
  complex*16, intent(in) :: in(*)
  complex*16, intent(out) :: out(*)
  call myzfft3d_mpiwrap0(in,out,dim,howmany,1,placeopt)
end subroutine myzfft3d_mpiwrap_forward

recursive subroutine myzfft3d_mpiwrap_backward(in,out,dim,howmany,placeopt)
  implicit none
  integer, intent(in) :: dim,howmany,placeopt
  complex*16, intent(in) :: in(*)
  complex*16, intent(out) :: out(*)
  call myzfft3d_mpiwrap0(in,out,dim,howmany,-1,placeopt)
end subroutine myzfft3d_mpiwrap_backward

recursive subroutine myzfft3d_mpiwrap0(in,out,dim,howmany,direction,placeopt)
  implicit none
  integer :: dim,nulltimes(10),howmany,ii,direction,placeopt
  complex*16, intent(in) :: in(dim**3,howmany)
  complex*16, intent(out) :: out(dim**3,howmany)
  complex*16,allocatable :: inlocal(:,:),outgather(:,:,:),outlocal(:,:)
  integer :: mystart, myend, mysize, myrank,nprocs

  call getmyranknprocs(myrank,nprocs)
  call checkdivisible(dim,nprocs)

  mystart=dim**3/nprocs*(myrank-1)+1
  myend=dim**3/nprocs*myrank
  mysize=dim**3/nprocs

  allocate(inlocal(mystart:myend,howmany), outlocal(mystart:myend,howmany),&
       outgather(mystart:myend,howmany,nprocs))

  inlocal(:,:)=in(mystart:myend,:)


  select case(direction)
  case(-1)
     if (placeopt.ne.1) then
        call ctdim(3)
        call cooleytukey_outofplace_backward_mpi(inlocal,outlocal,dim,dim,dim/nprocs,howmany)
     else
        call myzfft3d_par_backward(inlocal,outlocal,dim,nulltimes,howmany)
     endif
  case(1)
     if (placeopt.ne.1) then
        call ctdim(3)
        call cooleytukey_outofplace_forward_mpi(inlocal,outlocal,dim,dim,dim/nprocs,howmany)
     else
        call myzfft3d_par_forward(inlocal,outlocal,dim,nulltimes,howmany)
     endif
  case default
  print *, "ACK DIRECTION!!!!", direction; call mpistop()
end select

  call simpleallgather_complex(outlocal,outgather,dim**3/nprocs*howmany)
  do ii=1,howmany
     out(:,ii)=RESHAPE(outgather(:,ii,:),(/dim**3/))
  enddo
  deallocate(inlocal,outlocal,outgather)

end subroutine myzfft3d_mpiwrap0


recursive subroutine myzfft3d_par_forward(in,out,dim,times,howmany)
  implicit none
  integer, intent(in) :: dim,howmany
  complex*16, intent(in) :: in(*)
  complex*16, intent(out) :: out(*)
  integer, intent(inout) :: times(8)
  integer :: myrank,nprocs
  call getmyranknprocs(myrank,nprocs)
  call checkdivisible(dim,nprocs)
  call myzfft3d_par0(in,out,dim,times,howmany,nprocs,1)
end subroutine myzfft3d_par_forward

recursive subroutine myzfft3d_par_backward(in,out,dim,times,howmany)
  implicit none
  integer, intent(in) :: dim,howmany
  complex*16, intent(in) :: in(*)
  complex*16, intent(out) :: out(*)
  integer, intent(inout) :: times(8)
  integer :: myrank,nprocs
  call getmyranknprocs(myrank,nprocs)
  call checkdivisible(dim,nprocs)
  call myzfft3d_par0(in,out,dim,times,howmany,nprocs,-1)
end subroutine myzfft3d_par_backward


!!! adds to times

!!! times(1) = copy  times(2) = conjg  times(3) = ft
!!! from mytranspose times(4) = transpose   times(5) = mpi  times(6) = copy

recursive subroutine myzfft3d_par0(in,out,dim,times,howmany,nprocs,direction)
  use mytransposemod
  implicit none
  integer, intent(in) :: dim,howmany,nprocs,direction
  complex*16, intent(in) :: in(dim,dim,dim/nprocs,howmany)
  complex*16, intent(out) :: out(dim,dim,dim/nprocs,howmany)
  integer, intent(inout) :: times(8)
  integer :: ii,atime,btime
  complex*16 :: mywork(dim,dim,dim/nprocs,howmany) !! AUTOMATIC

  call checkdivisible(dim,nprocs)

  call myclock(atime)
  select case(direction)
  case(-1)
     out(:,:,:,:)=CONJG(in(:,:,:,:))
     call myclock(btime); times(2)=times(2)+btime-atime;
  case(1)
     out(:,:,:,:)=in(:,:,:,:)
     call myclock(btime); times(1)=times(1)+btime-atime;
  case default
     print *, "ACK PAR0 DIRECTION=",direction; call mpistop()
  end select

  do ii=1,3
     call myclock(atime)
     call myzfft1d( out, mywork, dim, dim**2/nprocs*howmany)
     call myclock(btime); times(3)=times(3)+btime-atime
     
!!! from mytranspose times(4) = transpose   times(5) = mpi  times(6) = copy
     call mytranspose(&
          mywork,  &
          out,  &
          dim/nprocs, &
          howmany,times(4:),nprocs)
  enddo

  call myclock(atime)
  select case(direction)
  case(-1)
     out(:,:,:,:)=CONJG(out(:,:,:,:))
     call myclock(btime); times(2)=times(2)+btime-atime
  case(1)
  case default
     print *, "ACK PAR0 DIRECTION=",direction; call mpistop()
  end select

end subroutine myzfft3d_par0

#endif
