
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
    integer, intent(in) :: dim,howmany
    complex*16 :: in(dim,howmany)    !! cannot be declared intent(in)...hmmm...
    complex*16, intent(out) :: out(dim,howmany)
    integer, parameter :: maxplans=10
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
       plans(thisplan) = fftw_plan_many_dft(1,dims,howmany,in,inembed,istride,idist,out,&
            onembed,ostride,odist,FFTW_FORWARD,FFTW_EXHAUSTIVE) 
       if (ftoutflag.ne.0) then
          print *, "       Done making a 1D FFT plan ", thisplan, dims, howmany
       endif
    endif
    icalleds(thisplan)=1    

    call fftw_execute_dft(plans(thisplan), in,out)

  end subroutine myzfft1d

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

end module fft1dsubmod

#endif
