

#include "Definitions.INC"

function facfunct(myenergy)
  use parameters
  implicit none
  real*8 :: myenergy
  DATATYPE :: facfunct,ccsum
  ccsum=1d0
  if (diffdipoleflag.ne.0) then
     if (myenergy.ne.0d0) then
        ccsum=1d0/((0d0,1d0)*myenergy)
     else
        ccsum=0d0
     endif
  endif
  facfunct=ccsum
end function facfunct


!! actually have numdata+1 data points in indipolearray

subroutine dipolecall(numdata, indipolearray,outename,outftname,which ,sflag)   !! which=1,2,3  =  x,y,z
  use parameters
  use mpimod
  implicit none

  DATATYPE :: indipolearray(0:numdata), temparray(0:numdata),facfunct,pots(3)
  integer :: i, numdata, which,getlen,jj,sflag
  real*8 :: estep, thistime, myenergy,sum1,sum2
  character (len=7) :: number
  character :: outftname*(*), outename*(*)
  complex*16 ::  fftrans(0:autosize), eft(0:autosize)

  fftrans=0.d0; eft=0d0

#ifdef REALGO
  OFLWR "Cant use dipolesub for real valued code."; CFLST
#endif

  if (diffdipoleflag.ne.0) then
     temparray(:)=0d0
     do i=1,numdata-1
        jj=min(min(i,numdata-i),4)
        select case(jj)
        case(1)
           temparray(i)= &
                - 1d0/2d0 * indipolearray(i-1) &
                + 1d0/2d0 * indipolearray(i+1)
        case(2)
           temparray(i)= &
                1d0/12d0 * indipolearray(i-2) &
                - 2d0/3d0 * indipolearray(i-1) &
                + 2d0/3d0 * indipolearray(i+1) &
                - 1d0/12d0 * indipolearray(i+2)
        case(3)
           temparray(i)= &
                - 1d0/60d0 * indipolearray(i-3) &
                + 3d0/20d0 * indipolearray(i-2) &
                - 3d0/4d0 * indipolearray(i-1) &
                + 3d0/4d0 * indipolearray(i+1) &
                - 3d0/20d0 * indipolearray(i+2) &
                + 1d0/60d0 * indipolearray(i+3)
        case(4)
           temparray(i)= &
                1d0/280d0 * indipolearray(i-4) &
                - 4d0/105d0 * indipolearray(i-3) &
                + 1d0/5d0 * indipolearray(i-2) &
                - 4d0/5d0 * indipolearray(i-1) &
                + 4d0/5d0 * indipolearray(i+1) &
                - 1d0/5d0 * indipolearray(i+2) &
                + 4d0/105d0 * indipolearray(i+3) &
                - 1d0/280d0 * indipolearray(i+4)
        end select
     end do
     temparray(:)=temparray(:)/par_timestep/autosteps
     do i=0,numdata
        fftrans(i) = temparray(i)  *cos(pi/2d0 * real(i,8)/real(numdata,8))**dipolewindowpower
        call vectdpot(i*par_timestep*autosteps,0,pots)
        eft(i)=pots(which)
     enddo
  else
     do i=0,numdata
        fftrans(i) = (indipolearray(i)-indipolearray(0))  *cos(pi/2d0 * real(i,8)/real(numdata,8))**dipolewindowpower
        call vectdpot(i*par_timestep*autosteps,0,pots)
        eft(i)=pots(which)
     enddo
  endif

  if (myrank.eq.1) then
     open(171,file=outename,status="unknown");     write(171,*) "#   ", numdata
     do i=0,numdata
        write(171,'(F18.12, T22, 400E20.8)')  i*par_timestep*autosteps, fftrans(i),indipolearray(i),eft(i)
     enddo
     close(171)
  endif

  sum1=0d0
  do i=0,numdata
     sum1=sum1+abs(fftrans(i))**2
  enddo

  call zfftf_wrap(numdata+1,fftrans(0))
  call zfftf_wrap(numdata+1,eft(0))

  sum2=0d0
  do i=0,numdata
     sum2=sum2+abs(fftrans(i))**2
  enddo
  OFLWR "FFT : NORM BEFORE, AFTER ", sum1,sum2/(numdata+1); CFL

  Estep=2*pi/par_timestep/autosteps/(numdata+1)

  thistime=numdata*par_timestep*autosteps

  if (myrank.eq.1) then
     open(171,file=outftname,status="unknown")
     do i=0,numdata
        myenergy=i*Estep
        write(171,'(F18.12, T22, 400E20.8)')  myenergy, fftrans(i)*facfunct(myenergy), eft(i), fftrans(i)*facfunct(myenergy)*ALLCON(eft(i))
     enddo
     close(171)
     if (sflag.ne.0) then
        write(number,'(I7)') 1000000+floor(thistime)
        open(171,file=outftname(1:getlen(outftname)-1)//number(2:7),status="unknown")
        do i=0,numdata
           myenergy=i*Estep
           write(171,'(F18.12, T22, 400E20.8)')  myenergy, fftrans(i)*facfunct(myenergy), eft(i), fftrans(i)*facfunct(myenergy)*ALLCON(eft(i))
        enddo
        close(171)
     endif
  endif

end subroutine dipolecall

