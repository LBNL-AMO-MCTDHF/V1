

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


function windowfunct(i,numdata)
  use parameters
  implicit none
  real*8 :: windowfunct
  integer,intent(in) :: i,numdata

  if (i.lt.0.or.i.gt.numdata) then
     OFLWR "ERROR, windowfunct ",i,numdata; CFLST
  endif

  if (ftwindowlength.ge.0) then
     if (numdata-i.lt.ftwindowlength) then
        windowfunct = cos( pi/2d0 * (i-(numdata-ftwindowlength)) / real(ftwindowlength,8) )**ftwindowpower      !! **2
     else
        windowfunct=1d0
     endif
  else
     windowfunct = cos( pi/2d0 * i / real(numdata,8) )**ftwindowpower
  endif

end function windowfunct


!! actually have numdata+1 data points in indipolearray

subroutine dipolecall(numdata, indipolearray,outename,outftname,which ,sflag)   !! which=1,2,3  =  x,y,z
  use parameters
  use mpimod
  implicit none

  DATATYPE :: indipolearray(0:numdata), temparray(0:numdata),facfunct,pots(3)
  integer :: i, numdata, which,getlen,jj,sflag
  real*8 :: estep, thistime, myenergy,sum1,sum2,xsecunits, windowfunct
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


        fftrans(i) = temparray(i)  * windowfunct(i,numdata)
        call vectdpot(i*par_timestep*autosteps,0,pots)   !! LENGTH GAUGE
        eft(i)=pots(which)
     enddo
  else
     do i=0,numdata
        fftrans(i) = (indipolearray(i)-indipolearray(0))  * windowfunct(i,numdata)
        call vectdpot(i*par_timestep*autosteps,0,pots)   !! LENGTH GAUGE
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

        write(171,*) "## Photon energy (column 1); D(omega) (2,3); E(omega) (4,5); response (6,7); cross sect (9)" 
        write(171,*) "## UNITLESS RESPONSE FUNCTION FOR ABSORPTION/EMISSION ( 2 omega Im(D(omega)E(omega)*) ) IN COLUMN 7"
        write(171,*) "## QUANTUM MECHANICAL PHOTOABSORPTION/EMISSION CROSS SECTION IN MEGABARNS (no factor of 1/3) IN COLUMN 9"
        write(171,*)

     do i=0,numdata
        myenergy=i*Estep

!! LENGTH GAUGE (electric field) WAS FT'ed , OK with usual formula multiply by wfi
!! RESPONSE FUNCTION FOR ABSORPTION/EMISSION IN COLUMN 7
!! QUANTUM MECHANICAL PHOTOABSORPTION/EMISSION CROSS SECTION IN MEGABARNS (no factor of 1/3) IN COLUMN NINE

        xsecunits = 5.291772108d0**2 * 4d0 * PI / 1.37036d2 * myenergy

!! NOW FACTOR (2 omega) IN COLUMNS 6,7   v1.16 12-2015

        write(171,'(F18.12, T22, 400E20.8)')  myenergy, fftrans(i)*facfunct(myenergy), eft(i), fftrans(i)*facfunct(myenergy)*ALLCON(eft(i)) * 2 * myenergy, &
             fftrans(i)*facfunct(myenergy)*ALLCON(eft(i)) / abs(eft(i)**2) * xsecunits
     enddo
     close(171)
     if (sflag.ne.0) then
        write(number,'(I7)') 1000000+floor(thistime)
        open(171,file=outftname(1:getlen(outftname)-1)//number(2:7),status="unknown")

        write(171,*) "## Photon energy (column 1); D(omega) (2,3); E(omega) (4,5); response (6,7); cross sect (9)" 
        write(171,*) "## UNITLESS RESPONSE FUNCTION FOR ABSORPTION/EMISSION ( 2 omega Im(D(omega)E(omega)*) ) IN COLUMN 7"
        write(171,*) "## QUANTUM MECHANICAL PHOTOABSORPTION/EMISSION CROSS SECTION IN MEGABARNS (no factor of 1/3) IN COLUMN 9"
     write(171,*)

        do i=0,numdata
           myenergy=i*Estep

           xsecunits = 5.291772108d0**2 * 4d0 * PI / 1.37036d2 * myenergy

!! NOW FACTOR (2 omega) IN COLUMNS 6,7   v1.16 12-2015

           write(171,'(F18.12, T22, 400E20.8)')  myenergy, fftrans(i)*facfunct(myenergy), eft(i), fftrans(i)*facfunct(myenergy)*ALLCON(eft(i)) * 2 * myenergy, &
                fftrans(i)*facfunct(myenergy)*ALLCON(eft(i)) / abs(eft(i)**2) * xsecunits
        enddo
        close(171)
     endif
  endif

end subroutine dipolecall

