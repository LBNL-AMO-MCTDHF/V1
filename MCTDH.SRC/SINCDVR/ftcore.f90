
!! myrank is 1:nprocs

module ftoutmod
  implicit none
  integer :: ftoutflag=0
  integer :: ftfileptr=6
end module ftoutmod

subroutine ftset(inoutflag,infileptr)
  use ftoutmod
  integer, intent(in) :: inoutflag,infileptr
  ftoutflag=inoutflag; ftfileptr=infileptr
end subroutine ftset


#ifdef FFTWFLAG

!! Old version myzfft1d() for intel, should not be needed; see myzfft1d_not() below

subroutine myzfft1d(in,out,dim,howmany)
  use ftoutmod
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
     if (ftoutflag.ne.0) then
        print *, "       Making a 1D FFT plan ", thisplan, dims, howmany
     endif
     plans(thisplan) = fftw_plan_many_dft(1,dims,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,FFTW_FORWARD,FFTW_EXHAUSTIVE) 
     if (ftoutflag.ne.0) then
        print *, "       Done making a 1D FFT plan ", thisplan, dims, howmany
     endif
  endif
  icalleds(thisplan)=1    

!  if (ftoutflag.ne.0) then
!     print *, "       Doing a 1D FFT ", thisplan
!  endif

  call fftw_execute_dft(plans(thisplan), in,out)

!  if (ftoutflag.ne.0) then
!     print *, "          Done with a 1D FFT ", thisplan
!  endif

end subroutine myzfft1d


!! Not sure why this didn't work.  Old version myzfft1d() for intel, above, should not be needed


!!$  subroutine myzfft1d_not(in,out,dim,howmany)
!!$    use, intrinsic :: iso_c_binding
!!$    implicit none
!!$    include "fftw3.f03"
!!$    integer, intent(in) :: dim,howmany
!!$    complex*16,intent(in) :: in(dim,howmany)
!!$    complex*16, intent(out) :: out(dim,howmany)
!!$    call myzfft1d0(1,in,out,dim,howmany)
!!$  end subroutine myzfft1d_not
!!$  
!!$  
!!$  subroutine myzfft1d_slowindex_local(in,out,dim1,dim2,howmany)
!!$    implicit none
!!$    integer, intent(in) :: dim1,dim2,howmany
!!$    complex*16, intent(in) :: in(dim1,dim2,howmany)
!!$    complex*16, intent(out) :: out(dim1,dim2,howmany)
!!$    call myzfft1d0(dim1,in,out,dim2,howmany)
!!$  end subroutine myzfft1d_slowindex_local
!!$  
!!$  
!!$  subroutine myzfft1d0(blockdim,in,out,dim,howmany)
!!$    use ftoutmod
!!$    use, intrinsic :: iso_c_binding
!!$    implicit none
!!$    include "fftw3.f03"
!!$    integer, intent(in) :: dim,howmany,blockdim
!!$    complex*16 :: in(blockdim,dim,howmany)    !! cannot be declared intent(in)...hmmm...
!!$    complex*16, intent(out) :: out(blockdim,dim,howmany)
!!$    integer, parameter :: maxplans=3
!!$    type(C_PTR),save :: plans(maxplans)
!!$    integer, save :: plandims(maxplans)=-999, planhowmany(maxplans)=-999,&
!!$         planblockdim(maxplans)=-999
!!$    integer,save :: icalleds(maxplans)=0, numplans=0
!!$    integer :: ostride,istride,onembed(1),inembed(1),idist,odist, dims(1),iplan,thisplan
!!$  
!!$  !!$  KEEPME           EITHER WORK BLOCKDIM=1            KEEPME
!!$  !!$
!!$  !!$  inembed(1)=dim; onembed(1)=dim; idist=dim; odist=dim; istride=1; ostride=1; dims(1)=dim
!!$  !!$  inembed(1)=dim; onembed(1)=dim; idist=1; odist=1; istride=1; ostride=1; dims(1)=dim
!!$  !!$
!!$  
!!$    inembed(1)=dim; onembed(1)=dim; idist=1; odist=1; istride=blockdim; ostride=blockdim; dims(1)=dim
!!$  
!!$    if (numplans.eq.0) then
!!$       numplans=1
!!$       thisplan=1
!!$       plandims(thisplan)=dim; planhowmany(thisplan)=howmany;
!!$       planblockdim(thisplan)=blockdim
!!$    else
!!$       thisplan= -99
!!$       do iplan=1,numplans
!!$          if (plandims(iplan).eq.dim.and.planhowmany(iplan).eq.howmany&
!!$               .and.planblockdim(iplan).eq.blockdim) then
!!$             if (icalleds(iplan).eq.0) then
!!$                print *, "ERROR, plan not done ",iplan,dim,howmany; call mpistop()
!!$             endif
!!$             thisplan=iplan
!!$             exit
!!$          endif
!!$       enddo
!!$       if (thisplan.eq.-99) then
!!$          if (numplans.eq.maxplans) then
!!$             print *,  "all plans taken!", maxplans; call mpistop()
!!$          endif
!!$          numplans=numplans+1
!!$          thisplan=numplans
!!$          plandims(thisplan)=dim; planhowmany(thisplan)=howmany;
!!$          planblockdim(thisplan)=blockdim
!!$       endif
!!$    endif
!!$    if (icalleds(thisplan).eq.0) then
!!$       if (ftoutflag.ne.0) then
!!$          print *, "       Making a 1D FFT plan! ", thisplan,  howmany, blockdim
!!$          print *, "       ", dims
!!$       endif
!!$       plans(thisplan) = fftw_plan_many_dft(1,dims,howmany*blockdim,in,inembed,istride,idist,out,onembed,ostride,odist,FFTW_FORWARD,FFTW_EXHAUSTIVE) 
!!$       if (ftoutflag.ne.0) then
!!$          print *, "       Done making a 1D FFT plan! ", thisplan,  howmany, blockdim
!!$       endif
!!$    endif
!!$    icalleds(thisplan)=1    
!!$  
!!$  !  if (ftoutflag.ne.0) then
!!$  !     print *, "       Doing a 1D FFT! ", thisplan
!!$  !  endif
!!$  
!!$    call fftw_execute_dft(plans(thisplan), in,out)
!!$  
!!$  !  if (ftoutflag.ne.0) then
!!$  !     print *, "          Done with a 1D FFT! ", thisplan
!!$  !  endif
!!$  
!!$  end subroutine myzfft1d0

module fft3dsubmod
contains

!! called by recursive subroutines so making recursive

  recursive subroutine myzfft3d(in,out,dim1,dim2,dim3,howmany)
    use ftoutmod
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
       if (ftoutflag.ne.0) then
          print *, "       Making a 3D fft plan ", thisplan, dims, howmany
       endif
       plans(thisplan) = fftw_plan_many_dft(3,dims,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,FFTW_FORWARD,FFTW_EXHAUSTIVE) 
       if (ftoutflag.ne.0) then
          print *, "        ...ok, made a 3D fft plan ", thisplan, dims, howmany
       endif
    endif
    icalleds(thisplan)=1    

!  if (ftoutflag.ne.0) then
!     print *, "       Doing a 3D fft ", thisplan
!  endif
    call fftw_execute_dft(plans(thisplan), in,out)
!  if (ftoutflag.ne.0) then
!     print *, "          Done with a 3D fft ", thisplan
!  endif

  end subroutine myzfft3d

end module fft3dsubmod

#else


subroutine myzfft1d(in,out,dim,howmany)
  implicit none
  integer, intent(in) :: dim,howmany
  integer :: k
  complex*16, intent(in) :: in(dim,howmany)
  complex*16, intent(out) :: out(dim,howmany)
  complex*16 :: wsave(4*dim+15)    !! AUTOMATIC OpenMP
  out(:,:)=in(:,:)
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(in,out,dim,howmany)
  wsave(:)=0d0
!$OMP DO SCHEDULE(STATIC)
  do k=1,howmany
     call zffti(dim,wsave(:))
     call zfftf(dim,out(:,k),wsave(:))
  enddo
!$OMP END DO
!$OMP END PARALLEL
end subroutine myzfft1d



!!$  !! OBVIOUSLY UNSATISFACTORY WITH DFFTPACK ROUTINES
!!$  !! OBVIOUSLY UNSATISFACTORY WITH DFFTPACK ROUTINES
!!$  !! OBVIOUSLY UNSATISFACTORY WITH DFFTPACK ROUTINES USED CURRENTLY:
!!$  
!!$  subroutine myzfft1d_slowindex_local(in,out,dim1,dim2,howmany)
!!$    implicit none
!!$    integer, intent(in) :: dim1,dim2,howmany
!!$    complex*16, intent(in) :: in(dim1,dim2,howmany)
!!$    complex*16, intent(out) :: out(dim1,dim2,howmany)
!!$    complex*16 :: intrans(dim2,dim1,howmany),outtrans(dim2,dim1,howmany),&
!!$         loctrans21(dim2,dim1),loctrans12(dim1,dim2)    !! AUTOMATIC
!!$    integer :: ii
!!$  
!!$    intrans=0; outtrans=0
!!$  
!!$  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,loctrans21)
!!$    loctrans21(:,:)=0
!!$  !$OMP DO SCHEDULE(DYNAMIC)
!!$    do ii=1,howmany
!!$       loctrans21(:,:)=TRANSPOSE(in(:,:,ii))
!!$       intrans(:,:,ii)=loctrans21(:,:)
!!$    enddo
!!$  !$OMP END DO
!!$  !$OMP END PARALLEL
!!$  
!!$    call myzfft1d(intrans,outtrans,dim2,dim1*howmany)
!!$  
!!$  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,loctrans12)
!!$    loctrans12(:,:)=0
!!$  !$OMP DO SCHEDULE(DYNAMIC)
!!$    do ii=1,howmany
!!$       loctrans12(:,:)=TRANSPOSE(outtrans(:,:,ii))
!!$       out(:,:,ii)=loctrans12(:,:)
!!$    enddo
!!$  !$OMP END DO
!!$  !$OMP END PARALLEL
!!$  
!!$  end subroutine myzfft1d_slowindex_local


module fft3dsubmod
contains

!! called by recursive so making recursive
  recursive subroutine fftblock_withtranspose(inout,dim1,dim2,dim3,howmany)
    implicit none
    integer :: dim1,dim2,dim3,howmany
!!!!  is dimensioned (dim1,dim2,dim3) on input. !!!!
    complex*16,intent(inout) :: inout(dim2,dim3,dim1,howmany) 
    complex*16 :: work1(dim1,dim2,dim3,howmany)  !! AUTOMATIC
    integer :: i

    work1=0
    call myzfft1d(inout,work1,dim1,dim2*dim3*howmany)
    do i=1,dim1
       inout(:,:,i,:)=work1(i,:,:,:)
    enddo
  end subroutine fftblock_withtranspose

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

end module fft3dsubmod

#endif

#ifdef MPIFLAG

module fftparsubmod
contains

!! times(1) = transpose   times(2) = mpi  times(3) = copy
!!   (123) -> (231)

!! called by subroutine that's called by recursive subroutine so making recursive
!!  etc.

recursive subroutine mytranspose_complex(in,out,blocksize,howmany,times,nprocs1,nprocs2)
  use pmpimod  !! box_comm
  implicit none
  integer,intent(in) :: blocksize,howmany,nprocs1,nprocs2
  integer,intent(inout) :: times(3)
  complex*16,intent(in) :: in(nprocs1*blocksize,nprocs2*blocksize,blocksize,howmany)
  complex*16,intent(out) :: out(nprocs1*blocksize,nprocs2*blocksize,blocksize,howmany)
  integer :: atime,btime,i,count,ii,iproc,j
  complex*16 :: intranspose(nprocs2*blocksize,blocksize,blocksize,howmany,nprocs1)  !!AUTOMATIC
  complex*16 :: outtemp(nprocs2*blocksize,blocksize,blocksize,howmany,nprocs1)      !!AUTOMATIC
  complex*16 :: outone(nprocs2*blocksize,blocksize,blocksize,nprocs1)               !!AUTOMATIC
  complex*16 :: inchop(blocksize,nprocs2*blocksize,blocksize,howmany)              !!AUTOMATIC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!    (123)->(231)    !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call myclock(atime)

  intranspose=0d0; outtemp=0; outone=0; inchop=0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,i)    !! IPROC IS SHARED (GOES WITH BARRIER, INCHOP SHARED)

  do iproc=1,nprocs1

     inchop(:,:,:,:)=in((iproc-1)*blocksize+1:iproc*blocksize,:,:,:)

!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
     do ii=1,howmany
        do i=1,blocksize

!! (123)->(231)  (:,:,7) -> (:,7,:)
!! (123)->(231)  (:,5,7) -> (5,7,:)

           intranspose(:,:,i,ii,iproc)=inchop(i,:,:,ii)

        enddo
     enddo
!$OMP END DO
!! *** OMP BARRIER *** !!   if inchop & iproc are shared
!$OMP BARRIER
  enddo

!$OMP END PARALLEL

  call myclock(btime); times(1)=times(1)+btime-atime; atime=btime

  outtemp(:,:,:,:,:)=0d0
  
  count=blocksize**3 * nprocs2 * howmany

  if (nprocs1.eq.nprocs2) then  !! orbparlevel=3

!! 231  (:,7,:) -> (:,:,7)

     call mympialltoall_complex(intranspose,outtemp,count)

  else

!! 231  (5,7,:) -> (7,5,:)

     call mympisendrecv_complex(intranspose,outtemp,&
          rankbybox(boxrank(1),boxrank(3),boxrank(2)), rankbybox(boxrank(1),boxrank(3),boxrank(2)), &
          999,count*nprocs1)

     intranspose(:,:,:,:,:)=outtemp(:,:,:,:,:)

!! 231  (7,5,:) -> (:,5,7)

     call mympialltoall_complex_local(intranspose,outtemp,count,BOX_COMM(1,boxrank(2),3))

  endif


  call myclock(btime); times(2)=times(2)+btime-atime; atime=btime

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)  !! ii IS SHARED (OUTTEMP IS SHARED; BARRIER)
  do ii=1,howmany

     outone(:,:,:,:)=outtemp(:,:,:,ii,:)
     if (nprocs1.eq.nprocs2) then

!! (231) collecting middle index, 3.. have first index,2

!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
        do i=1,blocksize
           do j=1,nprocs1*blocksize
              out(j,:,i,ii)=RESHAPE(outone(j,:,i,:),(/nprocs1*blocksize/))
           enddo
        enddo
!$OMP END DO
!! *** OMP BARRIER *** !!   if outone & ii are shared
!$OMP BARRIER
     else

!! (231) collecting first index,2

!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
        do i=1,blocksize
           do j=1,blocksize
              out(:,j,i,ii)=RESHAPE(outone(:,j,i,:),(/nprocs1*blocksize/))
           enddo
        enddo
!$OMP END DO
!! *** OMP BARRIER *** !!   if outone & ii are shared
!$OMP BARRIER
        
     endif
  enddo
!$OMP END PARALLEL

  call myclock(btime); times(3)=times(3)+btime-atime;

end subroutine mytranspose_complex


recursive subroutine myzfft3d_par_forward(in,out,dim,times,howmany)
  use pmpimod
  implicit none
  integer, intent(in) :: dim,howmany
  complex*16, intent(in) :: in(*)
  complex*16, intent(out) :: out(*)
  integer, intent(inout) :: times(8)
  select case(orbparlevel)
  case(3)
     call myzfft3d_par0(in,out,dim,times,howmany,nprocs,1,1,orbparlevel)
  case(2)
     call myzfft3d_par0(in,out,dim,times,howmany,sqnprocs,sqnprocs,1,orbparlevel)
  case default
     print *, "ORBPARLEVEL NOT SUP", orbparlevel; call mpistop()
  end select
end subroutine myzfft3d_par_forward

recursive subroutine myzfft3d_par_backward(in,out,dim,times,howmany)
  use pmpimod
  implicit none
  integer, intent(in) :: dim,howmany
  complex*16, intent(in) :: in(*)
  complex*16, intent(out) :: out(*)
  integer, intent(inout) :: times(8)
  select case(orbparlevel)
  case(3)
     call myzfft3d_par0(in,out,dim,times,howmany,nprocs,1,-1,orbparlevel)
  case (2)
     call myzfft3d_par0(in,out,dim,times,howmany,sqnprocs,sqnprocs,-1,orbparlevel)
  case default
     print *,  "ORBPARLEVEL NOT SUP", orbparlevel; call mpistop()
  end select

end subroutine myzfft3d_par_backward


!!! adds to times

!!! times(1) = copy  times(2) = conjg  times(3) = ft
!!! from mytranspose times(4) = transpose   times(5) = mpi  times(6) = copy

recursive subroutine myzfft3d_par0(in,out,dim,times,howmany,nprocs1,nprocs2,direction,oplevel)
  implicit none
  integer, intent(in) :: dim,howmany,nprocs1,nprocs2,direction,oplevel
  complex*16, intent(in) :: in(dim**3/nprocs1/nprocs2,howmany)
  complex*16, intent(out) :: out(dim**3/nprocs1/nprocs2,howmany)
  integer, intent(inout) :: times(8)
  integer :: ii,atime,btime
  complex*16 :: mywork(dim**3/nprocs1/nprocs2,howmany) !! AUTOMATIC

  mywork=0

  if (oplevel.ne.2.and.oplevel.ne.3) then
     print *, "OPLEVEL NOT SUP",oplevel; call mpistop()
  endif

  call myclock(atime)
  select case(direction)
  case(-1)
     out(:,:)=CONJG(in(:,:))
     call myclock(btime); times(2)=times(2)+btime-atime;
  case(1)
     out(:,:)=in(:,:)
     call myclock(btime); times(1)=times(1)+btime-atime;
  case default
     print *, "ACK PAR0 DIRECTION=",direction; call mpistop()
  end select

  do ii=1,3
     call myclock(atime)
     call myzfft1d( out, mywork, dim, dim**2/nprocs1/nprocs2*howmany)
     call myclock(btime); times(3)=times(3)+btime-atime
     
!!! from mytranspose times(4) = transpose   times(5) = mpi  times(6) = copy

     select case(oplevel)
     case(3)
        call mytranspose_complex(&
             mywork,  &
             out,  &
             dim/nprocs1, &
             howmany,times(4:),nprocs1,nprocs1)
     case(2)
        call mytranspose_complex(&
             mywork,  &
             out,  &
             dim/nprocs1, &
             howmany,times(4:),nprocs1,1)
     case default
        print *, "NOOOT SUP", oplevel; call mpistop()
     end select
  enddo

  call myclock(atime)
  select case(direction)
  case(-1)
     out(:,:)=CONJG(out(:,:))
     call myclock(btime); times(2)=times(2)+btime-atime
  case(1)
  case default
     print *, "ACK PAR0 DIRECTION=",direction; call mpistop()
  end select

end subroutine myzfft3d_par0

end module fftparsubmod

#endif


!!$  subroutine checkdivisible(number,divisor)
!!$    implicit none
!!$    integer :: number,divisor
!!$    if ((number/divisor)*divisor.ne.number) then
!!$       print *, "ACK NOT DIVISIBLE",number,divisor; call mpistop()
!!$    endif
!!$  end subroutine checkdivisible
!!$  
!!$  
!!$  subroutine myzfft3d_mpiwrap_forward(in,out,dim,howmany,placeopt)
!!$    implicit none
!!$    integer, intent(in) :: dim,howmany,placeopt
!!$    complex*16, intent(in) :: in(*)
!!$    complex*16, intent(out) :: out(*)
!!$    call myzfft3d_mpiwrap0(in,out,dim,howmany,1,placeopt)
!!$  end subroutine myzfft3d_mpiwrap_forward
!!$  
!!$  subroutine myzfft3d_mpiwrap_backward(in,out,dim,howmany,placeopt)
!!$    implicit none
!!$    integer, intent(in) :: dim,howmany,placeopt
!!$    complex*16, intent(in) :: in(*)
!!$    complex*16, intent(out) :: out(*)
!!$    call myzfft3d_mpiwrap0(in,out,dim,howmany,-1,placeopt)
!!$  end subroutine myzfft3d_mpiwrap_backward
!!$  
!!$  subroutine myzfft3d_mpiwrap0(in,out,dim,howmany,direction,placeopt)
!!$    use pmpimod
!!$    implicit none
!!$    integer :: dim,nulltimes(10),howmany,ii,direction,placeopt
!!$    complex*16, intent(in) :: in(dim**3,howmany)
!!$    complex*16, intent(out) :: out(dim**3,howmany)
!!$    complex*16,allocatable :: inlocal(:,:),outgather(:,:,:),outlocal(:,:)
!!$    integer :: mystart, myend, mysize
!!$  
!!$    if (orbparlevel.ne.3) then
!!$       print *,  "orbparlevel .ne. 3 not supported mpiwrap"; call mpistop()
!!$    endif
!!$  
!!$    call checkdivisible(dim,nprocs)
!!$  
!!$    mystart=dim**3/nprocs*(myrank-1)+1
!!$    myend=dim**3/nprocs*myrank
!!$    mysize=dim**3/nprocs
!!$  
!!$    allocate(inlocal(mystart:myend,howmany), outlocal(mystart:myend,howmany),&
!!$         outgather(mystart:myend,howmany,nprocs))
!!$  
!!$    outlocal=0; outgather=0
!!$    inlocal(:,:)=in(mystart:myend,:)
!!$  
!!$  
!!$    select case(direction)
!!$    case(-1)
!!$       if (placeopt.ne.1) then
!!$  !!$        call ctdim(3)
!!$          call cooleytukey_outofplace_backward_mpi(inlocal,outlocal,dim/procsplit(1),dim/procsplit(2),dim/procsplit(3),howmany)
!!$       else
!!$          call myzfft3d_par_backward(inlocal,outlocal,dim,nulltimes,howmany)
!!$       endif
!!$    case(1)
!!$       if (placeopt.ne.1) then
!!$  !!$        call ctdim(3)
!!$          call cooleytukey_outofplace_forward_mpi(inlocal,outlocal,dim/procsplit(1),dim/procsplit(2),dim/procsplit(3),howmany)
!!$       else
!!$          call myzfft3d_par_forward(inlocal,outlocal,dim,nulltimes,howmany)
!!$       endif
!!$    case default
!!$       print *, "ACK DIRECTION!!!!", direction; call mpistop()
!!$    end select
!!$    
!!$    call simpleallgather_complex(outlocal,outgather,dim**3/nprocs*howmany)
!!$    do ii=1,howmany
!!$       out(:,ii)=RESHAPE(outgather(:,ii,:),(/dim**3/))
!!$    enddo
!!$    deallocate(inlocal,outlocal,outgather)
!!$  
!!$  end subroutine myzfft3d_mpiwrap0
