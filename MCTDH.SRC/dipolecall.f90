

!! Subroutines for action 21 emission/absorption.  Windowfunct for actions 1 16 17 21.
!!   See also dipolesub.f90

#include "Definitions.INC"

function windowfunct(i,numdata)
  use parameters
  implicit none
  real*8 :: windowfunct
  integer,intent(in) :: i,numdata

  if (i.lt.0.or.i.gt.numdata) then
     OFLWR "ERROR, windowfunct ",i,numdata; CFLST
  endif

  if (fttriwindow.ne.0) then

     windowfunct = ( real(numdata-i,8) / real(numdata,8) )**ftwindowpower

  else

     if (ftwindowlength.ge.0) then
        if (numdata-i.lt.ftwindowlength) then
           windowfunct = cos( pi/2d0 * (i-(numdata-ftwindowlength)) &
                / real(ftwindowlength,8) )**ftwindowpower
        else
           windowfunct=1d0
        endif
     else
        if (ftwindowpower.eq.0) then
           windowfunct = ( 1 - sin( pi/2d0 * i / real(numdata,8) ) )
        else
           windowfunct = cos( pi/2d0 * i / real(numdata,8) )**ftwindowpower
        endif
     endif
  endif

end function windowfunct


!! actually have numdata+1 data points in indipolearray
!! which=1,2,3  =  x,y,z   4,5,6,7,8,9 = x+iy, etc. see dipolesub

module dipolecallsubmod
contains

subroutine dipolecall(numdata, indipolearrays, outename, outftname, outoworkname, &
     outtworkname, outophotonname, which, sflag)
  use parameters
  use mpimod
  use pulse_parameters !! numpulses
  use pulsesubmod
  implicit none
  integer,intent(in) :: numdata, which, sflag
  DATATYPE,intent(in) :: indipolearrays(0:numdata,3)
  character,intent(in) :: outftname*(*), outename*(*), outoworkname*(*), outtworkname*(*),&
       outophotonname*(*)
  complex*16,allocatable ::  fftrans(:),eft(:), all_eft(:,:),dipole_diff(:)
  DATATYPE :: pots(3,numpulses)
  real*8 :: estep, thistime, myenergy,xsecunits, windowfunct
  real*8, allocatable :: worksums(:,:), worksum0(:,:), exsums(:,:), totworksums(:),&
       totworksum0(:), totexsums(:), xsums(:)
  character (len=7) :: number
  integer :: i,getlen,myiostat,ipulse

#ifdef REALGO
  OFLWR "Cant use dipolesub for real valued code."; CFLST
#endif

  pots=0
  allocate(fftrans(0:numdata), eft(0:numdata), all_eft(0:numdata,numpulses),dipole_diff(0:numdata))
  fftrans=0.d0; eft=0d0; all_eft=0d0; dipole_diff=0d0
  allocate(worksums(0:numdata,numpulses), worksum0(0:numdata,numpulses),exsums(0:numdata,numpulses),&
       totworksums(0:numdata), totworksum0(0:numdata),totexsums(0:numdata), xsums(0:numdata))
  worksums=0; worksum0=0; exsums=0; totworksums=0; totworksum0=0; totexsums=0; xsums=0

  do i=0,numdata
     do ipulse=1,numpulses
        call vectdpot0(i*par_timestep*autosteps,0,pots(:,ipulse),-1,ipulse,ipulse) !! LENGTH
     enddo
     select case(which)
     case(1,2,3)
        all_eft(i,:)=pots(which,:)
        fftrans(i) = (indipolearrays(i,which)-indipolearrays(0,which))
     case(4)
        all_eft(i,:)=pots(1,:) + (0d0,1d0) * pots(2,:)
        fftrans(i) = (indipolearrays(i,1)-indipolearrays(0,1)) + (0d0,1d0)*(indipolearrays(i,2)-indipolearrays(0,2))
     case(5)
        all_eft(i,:)=pots(1,:) + (0d0,1d0) * pots(3,:)
        fftrans(i) = (indipolearrays(i,1)-indipolearrays(0,1)) + (0d0,1d0)*(indipolearrays(i,3)-indipolearrays(0,3))
     case(6)
        all_eft(i,:)=pots(2,:) + (0d0,1d0) * pots(1,:)
        fftrans(i) = (indipolearrays(i,2)-indipolearrays(0,2)) + (0d0,1d0)*(indipolearrays(i,1)-indipolearrays(0,1))
     case(7)
        all_eft(i,:)=pots(2,:) + (0d0,1d0) * pots(3,:)
        fftrans(i) = (indipolearrays(i,2)-indipolearrays(0,2)) + (0d0,1d0)*(indipolearrays(i,3)-indipolearrays(0,3))
     case(8)
        all_eft(i,:)=pots(3,:) + (0d0,1d0) * pots(1,:)
        fftrans(i) = (indipolearrays(i,3)-indipolearrays(0,3)) + (0d0,1d0)*(indipolearrays(i,1)-indipolearrays(0,1))
     case(9)
        all_eft(i,:)=pots(3,:) + (0d0,1d0) * pots(2,:)
        fftrans(i) = (indipolearrays(i,3)-indipolearrays(0,3)) + (0d0,1d0)*(indipolearrays(i,2)-indipolearrays(0,2))
     case default
        OFLWR "ACK WHICH DIPOLECALL", which; CFLST
     end select
  enddo

  call mydiff(numdata+1,fftrans,dipole_diff,.false.)
  dipole_diff(:)=dipole_diff(:) / par_timestep / autosteps

!! dividing and multiplying for clarity not math
!! numbers are real-valued which=1,2,3 x,y,z
!! otherwise with which = 4 through 9 take real part with conjugate like below

  worksum0(0,:) = (-1) * real( dipole_diff(0) * conjg(all_eft(0,:)) , 8) * par_timestep * autosteps
  do i=1,numdata
     worksum0(i,:)=worksum0(i-1,:) - real( dipole_diff(i) * conjg(all_eft(i,:)) , 8) * par_timestep * autosteps
  enddo
  totworksum0=0d0
  do i=1,numpulses
     totworksum0(:)=totworksum0(:)+worksum0(:,i)
  enddo

  do i=0,numdata
     fftrans(i) = fftrans(i) * windowfunct(i,numdata)
  enddo

  if (pulsewindowtoo.ne.0) then
     do i=0,numdata
        all_eft(i,:)=all_eft(i,:) * windowfunct(i,numdata)
     enddo
  endif

  eft=0
  do ipulse=1,numpulses
     eft(:)=eft(:)+all_eft(:,ipulse)
  enddo

  if (myrank.eq.1) then
     open(171,file=outename,status="unknown",iostat=myiostat)
     call checkiostat(myiostat,"opening "//outename)
     write(171,*,iostat=myiostat) "#   ", numdata
     call checkiostat(myiostat,"writing "//outename)
     do i=0,numdata
        write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat)  i*par_timestep*autosteps, &
             fftrans(i),eft(i),all_eft(i,:)
     enddo
     call checkiostat(myiostat,"writing "//outename)
     close(171)

     open(171,file=outtworkname,status="unknown",iostat=myiostat)
     call checkiostat(myiostat,"opening "//outtworkname)
     do i=0,numdata
        write(171,'(A25,F10.5,400F15.10)') " EACH PULSE WORK T= ", i*par_timestep*autosteps,&
             totworksum0(i),worksum0(i,:)
     enddo
     close(171)
  endif

  call zfftf_wrap_diff(numdata+1,fftrans(:),ftdiff)
  call zfftf_wrap(numdata+1,eft(:))
  do ipulse=1,numpulses
     call zfftf_wrap(numdata+1,all_eft(:,ipulse))
  enddo

  fftrans(:)=fftrans(:)     * par_timestep * autosteps
  eft(:)=eft(:)             * par_timestep * autosteps
  all_eft(:,:)=all_eft(:,:) * par_timestep * autosteps

  Estep=2*pi/par_timestep/autosteps/(numdata+1)

  thistime=numdata*par_timestep*autosteps

  xsums(0)= Estep * imag(fftrans(0)*conjg(eft(0))) / abs(eft(0)**2) * myenergy * 2 / PI
  exsums(0,:) = Estep * imag(fftrans(0)*conjg(all_eft(0,:))) / PI
  worksums(0,:) = Estep * imag(fftrans(0)*conjg(all_eft(0,:))) / PI * myenergy

  do i=1,numdata
     myenergy=i*Estep

!! xsum sums to N for N electrons
     if (myenergy.ge.dipolesumstart.and.myenergy.le.dipolesumend) then
        xsums(i)=xsums(i-1) + Estep * imag(fftrans(i)*conjg(eft(i))) / abs(eft(i)**2) * myenergy * 2 / PI
        exsums(i,:)  =  exsums(i-1,:) + Estep * imag(fftrans(i)*conjg(all_eft(i,:))) / PI
        worksums(i,:)=worksums(i-1,:) + Estep * imag(fftrans(i)*conjg(all_eft(i,:))) / PI * myenergy
     endif
  enddo
  totexsums=0
  totworksums=0
  do i=1,numpulses
     totexsums(:)=totexsums(:)+exsums(:,i)
     totworksums(:)=totworksums(:)+worksums(:,i)
  enddo

  if (myrank.eq.1) then
     open(171,file=outftname,status="unknown",iostat=myiostat)
     call checkiostat(myiostat,"opening "//outftname)
     write(171,'(A120)',iostat=myiostat) &
          "## Photon energy (column 1); D(omega) (2,3); E(omega) (4,5); response (6,7); cross sect (9); integrated (10)" 
     call checkiostat(myiostat,"writing "//outftname)
     write(171,'(A120)') "## UNITLESS RESPONSE FUNCTION FOR ABSORPTION/EMISSION 2 omega im(D(omega)E(omega)^*) IN COLUMN 7"
     write(171,'(A120)') "## QUANTUM MECHANICAL PHOTOABSORPTION/EMISSION CROSS SECTION IN MEGABARNS (no factor of 1/3) IN COLUMN NINE"
     write(171,'(A120)') "## INTEGRATED DIFFERENTIAL OSCILLATOR STRENGTH (FOR SUM RULE) IN COLUMN 10"
     write(171,*)

     do i=0,numdata
        myenergy=i*Estep

!! LENGTH GAUGE (electric field) WAS FT'ed , OK with usual formula multiply by wfi
!! UNITLESS RESPONSE FUNCTION FOR ABSORPTION/EMISSION 2 omega im(D(omega)E(omega)^*) IN COLUMN 7
!! QUANTUM MECHANICAL PHOTOABSORPTION/EMISSION CROSS SECTION IN MEGABARNS (no factor of 1/3) IN COLUMN NINE
!! INTEGRATED DIFFERENTIAL OSCILLATOR STRENGTH (FOR SUM RULE) IN COLUMN 10

        xsecunits = 5.291772108d0**2 * 4d0 * PI / 1.37036d2 * myenergy

!! NOW FACTOR (2 omega) IN COLUMNS 6,7   v1.16 12-2015

        write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat)  myenergy, &
             fftrans(i), eft(i), fftrans(i)*conjg(eft(i)) * 2 * myenergy, &
             fftrans(i)*conjg(eft(i)) / abs(eft(i)**2) * xsecunits, xsums(i)
     enddo
     call checkiostat(myiostat,"writing "//outftname)
     close(171)

!!  NUMBER OF PHOTONS ABSORBED AND AND WORK DONE BY EACH PULSE
!!  worksum0 the time integral converges right after pulse is finished... others take longer

     open(171,file=outoworkname,status="unknown",iostat=myiostat)
     call checkiostat(myiostat,"opening "//outoworkname)
     do i=0,numdata
        write(171,'(A25,F10.5,400F15.10)') " WORK EACH PULSE E= ", i*Estep, totworksums(i),worksums(i,:)
     enddo
     close(171)

     open(171,file=outophotonname,status="unknown",iostat=myiostat)
     call checkiostat(myiostat,"opening "//outophotonname)
     do i=0,numdata
        write(171,'(A25,F10.5,400F15.10)') "PHOTONS EACH PULSE E= ", i*Estep, totexsums(i),exsums(i,:)
     enddo
     close(171)

     if (sflag.ne.0) then
        write(number,'(I7)') 1000000+floor(thistime)
        open(171,file=outftname(1:getlen(outftname)-1)//number(2:7),status="unknown",iostat=myiostat)
        call checkiostat(myiostat,"opening "//outftname)
        write(171,'(A120)',iostat=myiostat) &
             "## Photon energy (column 1); D(omega) (2,3); E(omega) (4,5); response (6,7); cross sect (9); integrated (10)" 
        call checkiostat(myiostat,"writing "//outftname)
        write(171,'(A120)') "## UNITLESS RESPONSE FUNCTION FOR ABSORPTION/EMISSION 2 omega im(D(omega)E(omega)^*) IN COLUMN 7"
        write(171,'(A120)') "## QUANTUM MECHANICAL PHOTOABSORPTION/EMISSION CROSS SECTION IN MEGABARNS (no factor of 1/3) IN COLUMN NINE"
        write(171,'(A120)') "## INTEGRATED DIFFERENTIAL OSCILLATOR STRENGTH (FOR SUM RULE) IN COLUMN 10"
        write(171,*)

        do i=0,numdata
           myenergy=i*Estep

           xsecunits = 5.291772108d0**2 * 4d0 * PI / 1.37036d2 * myenergy

!! NOW FACTOR (2 omega) IN COLUMNS 6,7   v1.16 12-2015

           write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat)  myenergy, &
                fftrans(i), eft(i), fftrans(i)*conjg(eft(i)) * 2 * myenergy, &
                fftrans(i)*conjg(eft(i)) / abs(eft(i)**2) * xsecunits, xsums(i)
        enddo
        call checkiostat(myiostat,"writing "//outftname)
        close(171)
     endif
  endif

  deallocate(fftrans,eft,all_eft,dipole_diff,worksums,worksum0,exsums,totworksums,&
       totworksum0,totexsums,xsums)

end subroutine dipolecall

end module dipolecallsubmod
