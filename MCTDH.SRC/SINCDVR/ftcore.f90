
!! myrank is 1:nprocs

#ifdef FFTWFLAG

recursive subroutine myzfft1d(in,out,dim,howmany)
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
end subroutine myzfft1d


recursive subroutine myzfft3d(in,out,dim1,dim2,dim3,howmany)
  use, intrinsic :: iso_c_binding
  implicit none
  include "fftw3.f03"
  integer :: howmany,ii,dim1,dim2,dim3
  type(C_PTR),save :: plan
  integer,save :: icalled=0
  complex*16 :: in(dim1,dim2,dim3,howmany),out(dim1,dim2,dim3,howmany)
  integer, save :: savedim1,savedim2,savedim3

  if (icalled.eq.0) then
     plan=fftw_plan_dft_3d(dim1,dim2,dim3,in(:,:,:,1),out(:,:,:,1),FFTW_FORWARD,FFTW_EXHAUSTIVE) 
     savedim1=dim1; savedim2=dim2; savedim3=dim3
  else
     if (dim1.ne.savedim1.or.dim2.ne.savedim2.or.dim3.ne.savedim3) then
        print *, "ERROR SAVE DIMS MYZFFT3D FFTW ",dim1,dim2,dim3,savedim1,savedim2,savedim3; call mpistop()
     endif
  endif
  icalled=1
  do ii=1,howmany
     call fftw_execute_dft(plan, in(:,:,:,ii),out(:,:,:,ii))
  enddo

!!$  call fftw_destroy_plan(plan)

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


subroutine myzfft3d(in,out,dim1,dim2,dim3,howmany)
  implicit none
  integer :: dim1,dim2,dim3,howmany
  complex*16, intent(in) :: in(dim1,dim2,dim3,howmany)
  complex*16, intent(out) :: out(dim1,dim2,dim3,howmany)
  out(:,:,:,:)=in(:,:,:,:)
  call fftblock_withtranspose(out,dim1,dim2,dim3,howmany)
  call fftblock_withtranspose(out,dim2,dim3,dim1,howmany)
  call fftblock_withtranspose(out,dim3,dim1,dim2,howmany)
end subroutine myzfft3d


subroutine fftblock_withtranspose(inout,dim1,dim2,dim3,howmany)
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
  complex*16 :: inchop(nprocs*blocksize,blocksize,blocksize,howmany)  !!AUTOMATIC
  complex*16,intent(out) :: out(nprocs*blocksize,nprocs*blocksize,blocksize,howmany)
  integer :: atime,btime,i,count,ii,totsize,iproc,checkprocs,myrank,j
  complex*16 :: intranspose(blocksize,nprocs*blocksize,blocksize,howmany,nprocs)  !!AUTOMATIC
  complex*16 :: outtemp(blocksize,nprocs*blocksize,blocksize,howmany,nprocs)      !!AUTOMATIC
  complex*16 :: outone(blocksize,nprocs*blocksize,blocksize,nprocs)               !!AUTOMATIC
  
  call getmyranknprocs(myrank,checkprocs)
  if (nprocs.ne.checkprocs) then
     print *, "ACK CHECKPROCS",nprocs,checkprocs; call mpistop()
  endif

  totsize=blocksize*nprocs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!    (123)->(312)    !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call myclock(atime)
  
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

!! *** OMP BARRIER *** !!   if inchop is shared
!$OMP BARRIER

  enddo
! (implied barrier at end)
!$OMP END PARALLEL
  call myclock(btime); times(1)=times(1)+btime-atime; atime=btime
  
  count=blocksize**2 * totsize * howmany
  
  call mympialltoall_complex(intranspose,outtemp,count)
  call myclock(btime); times(2)=times(2)+btime-atime; atime=btime

!!  complex*16,intent(out) :: out(nprocs*blocksize,nprocs*blocksize,blocksize,howmany)
!!  complex*16 :: outtemp(blocksize,nprocs*blocksize,blocksize,howmany,nprocs)     
!!  complex*16 :: outone(blocksize,nprocs*blocksize,blocksize,nprocs)

  do iproc=1,nprocs
     out((iproc-1)*blocksize+1:iproc*blocksize,:,:,:)=outtemp(:,:,:,:,iproc)
  enddo

!  do ii=1,howmany
!     outone(:,:,:,:)=outtemp(:,:,:,ii,:)
!     do i=1,blocksize
!        do j=1,nprocs*blocksize
!           out(:,j,i,ii)=RESHAPE(outone(:,j,i,:),(/nprocs*blocksize/))
!        enddo
!     enddo
!  enddo


  call myclock(btime); times(3)=times(3)+btime-atime;
  
end subroutine mytranspose
end module  
  

subroutine checkdivisible(number,divisor)
  implicit none
  integer :: number,divisor
  if ((number/divisor)*divisor.ne.number) then
     print *, "ACK NOT DIVISIBLE",number,divisor; call mpistop()
  endif
end subroutine checkdivisible


recursive subroutine myzfft3d_mpiwrap(in,out,dim,howmany)
  implicit none
  integer :: dim,nulltimes(10),howmany,ii
  complex*16, intent(in) :: in(dim**3,howmany)
  complex*16, intent(out) :: out(dim**3,howmany)
  complex*16,allocatable :: inlocal(:,:),outgather(:,:,:),outlocal(:,:)
  integer :: mystart, myend, mysize, myrank, nprocs

  call getmyranknprocs(myrank,nprocs)
  call checkdivisible(dim,nprocs)

  mystart=dim**3/nprocs*(myrank-1)+1
  myend=dim**3/nprocs*myrank
  mysize=dim**3/nprocs

  allocate(inlocal(mystart:myend,howmany), outlocal(mystart:myend,howmany),&
       outgather(mystart:myend,howmany,nprocs))

  inlocal(:,:)=in(mystart:myend,:)

  call myzfft3d_par(inlocal,outlocal,dim,nulltimes,howmany)
  call simpleallgather_complex(outlocal,outgather,dim**3/nprocs*howmany)
  do ii=1,howmany
     out(:,ii)=RESHAPE(outgather(:,ii,:),(/dim**3/))
  enddo
  deallocate(inlocal,outlocal,outgather)

end subroutine myzfft3d_mpiwrap



!! adds to times

!!! times(1) = zero   times(2)=fourier
!!! from mytranspose times(3) = transpose   times(4) = mpi  times(5) = copy

recursive subroutine myzfft3d_par(in,out,dim,times,howmany)
  implicit none
  integer, intent(in) :: dim,howmany
  complex*16, intent(in) :: in(*)
  complex*16, intent(out) :: out(*)
  integer, intent(inout) :: times(8)
  integer :: myrank,nprocs
  call getmyranknprocs(myrank,nprocs)
  call checkdivisible(dim,nprocs)
  call myzfft3d_par0(in,out,dim,times,howmany,nprocs)
end subroutine myzfft3d_par


recursive subroutine myzfft3d_par0(in,out,dim,times,howmany,nprocs)
  use mytransposemod
  implicit none
  integer, intent(in) :: dim,howmany,nprocs
  complex*16, intent(in) :: in(dim,dim,dim/nprocs,howmany)
  complex*16, intent(out) :: out(dim,dim,dim/nprocs,howmany)
  integer, intent(inout) :: times(8)
  integer :: ii,atime,btime
  complex*16 :: mywork(dim,dim,dim/nprocs,howmany),tempout(dim,dim,dim/nprocs,howmany)

  call myclock(atime)

  call checkdivisible(dim,nprocs)

  tempout(:,:,:,:)=in(:,:,:,:)
  call myclock(btime); times(1)=times(1)+btime-atime;

  do ii=1,3
     call myclock(atime)
     call myzfft1d( tempout, mywork, dim, dim**2/nprocs*howmany)
     call myclock(btime); times(2)=times(2)+btime-atime; atime=btime
     
!!! from mytranspose times(3) = transpose   times(4) = mpi  times(5) = copy
     call mytranspose(&
          mywork,  &
          tempout,  &
          dim/nprocs, &
          howmany,times(3:),nprocs)
  enddo
  out(:,:,:,:)=tempout(:,:,:,:)

end subroutine myzfft3d_par0

#endif



