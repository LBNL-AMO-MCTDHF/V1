
!! ACTION ROUTINES FOR CALLING DISSOCIATIVE FLUX ROUTINES

subroutine fluxcall(numdata, forwardovl, timestep, eground, fftrans)  !! returns estep and fftrans
  implicit none

  complex*16 :: forwardovl(1:numdata)
  integer :: i, numdata
  real*8 :: piover2,  timestep, eground
  complex*16 :: fftrans(numdata)
  real*8, allocatable :: wsave(:)

  allocate( wsave(20*numdata+15) )
  piover2=2.d0 * atan(1.d0)

  !! offset for fftrans!
  fftrans=forwardovl

  do i=1,numdata
     fftrans(i)=fftrans(i) * exp((0.d0,-1.d0) * (i-numdata/2-1)* timestep * eground * 2.d0)
  enddo

  call zffti(numdata,wsave)
  call zfftf(numdata,fftrans, wsave)

  fftrans=fftrans*4*timestep

!! NEW
!  do i=1,numdata
!0     fftrans(i)=fftrans(i)*exp((0.d0,1.d0)*(i-1)*(numdata-1)*pi/real(numdata))
! enddo

  deallocate(  wsave ) 

end subroutine fluxcall





