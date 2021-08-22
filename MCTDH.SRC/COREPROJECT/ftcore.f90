
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

module fft1dsubmod
contains


  subroutine myzfft1d(in,out,dim,howmany)
    use ftoutmod
    use, intrinsic :: iso_c_binding
    implicit none
    include "fftw3.f03"

!! CAN'T GET BLOCKDIM TO WORK. TOO BAD.
!!
!! Otherwise this subroutine would be myzfft1d_slowindex and myzfft1d would
!! call myzfft1d_slowindex instead of the other way around and we would not
!! need to transpose the data like myzfft1d_slowindex below now does.  I
!! simply want to 1D FT the middle column of a 3rd rank tensor and I can't
!! figure out how to do it with the onembed, idist, ostride, etc. variables.
!!
    integer,parameter :: blockdim=1
    integer, intent(in) :: dim,howmany
    complex*16 :: in(blockdim,dim,howmany)    !! cannot be declared intent(in)...hmmm...
    complex*16, intent(out) :: out(blockdim,dim,howmany)
    integer, parameter :: maxplans=127
    type(C_PTR),save :: plans(maxplans)
    integer, save :: plandims(maxplans)=-999, planhowmany(maxplans)=-999,&
         planblockdim(maxplans)=-999
    integer,save :: icalleds(maxplans)=0, numplans=0
    integer :: ostride,istride,onembed(1),inembed(1),idist,odist, dims(1),iplan,thisplan

    inembed(1)=dim; onembed(1)=dim; idist=dim; odist=dim; 
    istride=blockdim; ostride=blockdim; dims(1)=dim

    if (numplans.eq.0) then
       numplans=1
       thisplan=1
       plandims(thisplan)=dim
       planhowmany(thisplan)=howmany
       planblockdim(thisplan)=blockdim
    else
       thisplan= -99
       do iplan=1,numplans
          if (plandims(iplan).eq.dim.and.planhowmany(iplan).eq.howmany&
               .and.planblockdim(iplan).eq.blockdim) then
             if (icalleds(iplan).eq.0) then
                print *, "ERROR, plan not done ",iplan,blockdim,dim,howmany
                call mpistop()
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
          plandims(thisplan)=dim
          planhowmany(thisplan)=howmany
          planblockdim(thisplan)=blockdim
       endif
    endif
    if (icalleds(thisplan).eq.0) then
       if (ftoutflag.ne.0) then
          print *, "   Making 1D FFT plan ", thisplan, blockdim, dims, howmany
       endif
       plans(thisplan) = fftw_plan_many_dft(1,dims,howmany*blockdim,in,inembed,istride,idist,out,&
            onembed,ostride,odist,FFTW_FORWARD,FFTW_EXHAUSTIVE) 
!       if (ftoutflag.ne.0) then
!          print *, "       Done with 1D plan ", thisplan, blockdim, dims, howmany
!       endif
    endif
    icalleds(thisplan)=1    
    
    call fftw_execute_dft(plans(thisplan), in,out)

  end subroutine myzfft1d


  subroutine myzfft1d_slowindex(in,out,blocksize,size,howmany)
    implicit none
    integer,intent(in) :: blocksize,size,howmany
    complex*16,intent(in) :: in(blocksize,size,howmany)
    complex*16,intent(out) :: out(blocksize,size,howmany)
    complex*16,allocatable :: intrans(:,:,:), outtrans(:,:,:)
    integer :: ii

    if (blocksize.eq.1) then
       call myzfft1d(in,out,size,howmany)
    else
       allocate( intrans(size,blocksize,howmany), outtrans(size,blocksize,howmany) )
       intrans(:,:,:)=0; outtrans(:,:,:)=0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(DYNAMIC)
       do ii=1,howmany
          intrans(:,:,ii)=TRANSPOSE(in(:,:,ii))
       enddo
!$OMP END DO
!$OMP END PARALLEL

       call myzfft1d(intrans,outtrans,size,blocksize*howmany)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(DYNAMIC)
       do ii=1,howmany
          out(:,:,ii)=TRANSPOSE(outtrans(:,:,ii))
       enddo
!$OMP END DO
!$OMP END PARALLEL
       deallocate(intrans,outtrans)
    endif

  end subroutine myzfft1d_slowindex

end module fft1dsubmod

#else

module fft1dsubmod
contains

  subroutine myzfft1d(in,out,dim,howmany)
    implicit none
    integer, intent(in) :: dim,howmany
    integer :: k
    complex*16, intent(in) :: in(dim,howmany)
    complex*16, intent(out) :: out(dim,howmany)
    complex*16 :: csave(2*dim)    !! AUTOMATIC OpenMP
    real*8     :: rsave(2*dim)    !! AUTOMATIC OpenMP
    integer    :: ifac(15)        !! AUTOMATIC OpenMP
    out(:,:)=in(:,:)
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(in,out,dim,howmany)
    csave(:)=0d0; rsave(:)=0d0; ifac(:)=0
    call cffti1(dim,rsave,ifac)                    ! now here 08-2021
    !$OMP DO SCHEDULE(STATIC)
    do k=1,howmany
       ! note that initializer is inside the loop run every time
       !  there was a reason (probably a crash avoided).  08-2021 replacing
       !  zfftf with cfftf1 which may fix the crash.  bringing it outside
       call cffti1(dim,rsave,ifac)   ! 
       call cfftf1(dim,out(:,k),csave,rsave,ifac)  ! there was a reason.
    enddo
!$OMP END DO
!$OMP END PARALLEL
  end subroutine myzfft1d


  subroutine myzfft1d_slowindex(in,out,blocksize,size,howmany)
    implicit none
    integer,intent(in) :: blocksize,size,howmany
    complex*16,intent(in) :: in(blocksize,size,howmany)
    complex*16,intent(out) :: out(blocksize,size,howmany)
    complex*16,allocatable :: intrans(:,:,:), outtrans(:,:,:)
    integer :: ii

    if (blocksize.eq.1) then
       call myzfft1d(in,out,size,howmany)
    else
       allocate( intrans(size,blocksize,howmany), outtrans(size,blocksize,howmany) )
       intrans(:,:,:)=0; outtrans(:,:,:)=0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(DYNAMIC)
       do ii=1,howmany
          intrans(:,:,ii)=TRANSPOSE(in(:,:,ii))
       enddo
!$OMP END DO
!$OMP END PARALLEL

       call myzfft1d(intrans,outtrans,size,blocksize*howmany)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(DYNAMIC)
       do ii=1,howmany
          out(:,:,ii)=TRANSPOSE(outtrans(:,:,ii))
       enddo
!$OMP END DO
!$OMP END PARALLEL
       deallocate(intrans,outtrans)
    endif

  end subroutine myzfft1d_slowindex


end module fft1dsubmod

#endif



