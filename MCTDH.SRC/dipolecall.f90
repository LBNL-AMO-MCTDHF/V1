

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
!! which=1,2,3  =  x,y,z

subroutine dipolecall(numdata, indipolearray,outename,outftname,which ,sflag)
  use parameters
  use mpimod
  use pulse_parameters !! numpulses
  use pulsesubmod
  implicit none

  DATATYPE :: indipolearray(0:numdata),pots(3)
  integer :: i, numdata, which,getlen,sflag,myiostat,ipulse
  real*8 :: estep, thistime, myenergy,xsecunits, windowfunct, xsum, &
       worksums(numpulses), worksum0(numpulses),exsums(numpulses)
  character (len=7) :: number
  character :: outftname*(*), outename*(*)
  complex*16,allocatable ::  fftrans(:),eft(:), all_eft(:,:),dipole_diff(:)
  character (len=3) :: xlabel(3) = (/ " X ", " Y ", " Z " /)

  allocate(fftrans(0:numdata), eft(0:numdata), all_eft(0:numdata,numpulses),dipole_diff(0:numdata))

  fftrans=0.d0; eft=0d0; all_eft=0d0; dipole_diff=0d0

#ifdef REALGO
  OFLWR "Cant use dipolesub for real valued code."; CFLST
#endif

  do i=0,numdata
     fftrans(i) = (indipolearray(i)-indipolearray(0))
     do ipulse=1,numpulses
        call vectdpot0(i*par_timestep*autosteps,0,pots,-1,ipulse,ipulse) !! LENGTH
        all_eft(i,ipulse)=pots(which)
     enddo
  enddo

  call mydiff(numdata+1,fftrans,dipole_diff,.false.)
  dipole_diff(:)=dipole_diff(:) / par_timestep / autosteps

  worksum0(:)=0
  do i=0,numdata
     worksum0(:)=worksum0(:) - real( dipole_diff(i) * all_eft(i,:) , 8) * & !! is real, ok
          par_timestep * autosteps        !! dividing and multiplying for clarity not math
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
             fftrans(i),indipolearray(i),eft(i),all_eft(i,:)
     enddo
     call checkiostat(myiostat,"writing "//outename)
     close(171)
  endif

  call zfftf_wrap_diff(numdata+1,fftrans(0:),ftdiff)
  call zfftf_wrap(numdata+1,eft(0:))
  do ipulse=1,numpulses
     call zfftf_wrap(numdata+1,all_eft(0:,ipulse))
  enddo

  fftrans(:)=fftrans(:)     * par_timestep * autosteps
  eft(:)=eft(:)             * par_timestep * autosteps
  all_eft(:,:)=all_eft(:,:) * par_timestep * autosteps

  Estep=2*pi/par_timestep/autosteps/(numdata+1)

  thistime=numdata*par_timestep*autosteps

  if (myrank.eq.1) then
     open(171,file=outftname,status="unknown",iostat=myiostat)
     call checkiostat(myiostat,"opening "//outftname)
     write(171,'(A120)',iostat=myiostat) &
          "## Photon energy (column 1); D(omega) (2,3); E(omega) (4,5); response (6,7); cross sect (9); integrated (10)" 
     call checkiostat(myiostat,"writing "//outftname)
     write(171,'(A120)') "## UNITLESS RESPONSE FUNCTION FOR ABSORPTION/EMISSION 2 omega im(D(omega)E(omega)^*) IN COLUMN 7"
     write(171,'(A120)') "## QUANTUM MECHANICAL PHOTOABSORPTION/EMISSION CROSS SECTION IN MEGABARNS (no factor of 1/3) IN COLUMN NINE"
     write(171,'(A120)') "## INTEGRATED DIFFERENTIAL OSCILLATOR STRENGTH (FOR SUM RULE) IN COLUMN 10"
     write(171,'(A120)') "## WORK DONE BY EACH PULSE (INTEGRATED) IN COLUMNS 11...."
     write(171,*)

     xsum=0d0
     worksums(:)=0d0
     exsums(:)=0d0

     do i=0,numdata
        myenergy=i*Estep

!! xsum sums to N for N electrons
        if (myenergy.ge.dipolesumstart.and.myenergy.le.dipolesumend) then
           xsum=xsum + Estep * imag(fftrans(i)*conjg(eft(i))) / abs(eft(i)**2) * myenergy * 2 / PI
           exsums(:)  =  exsums(:) + Estep * imag(fftrans(i)*conjg(all_eft(i,:))) / PI
           worksums(:)=worksums(:) + Estep * imag(fftrans(i)*conjg(all_eft(i,:))) / PI * myenergy
        endif

!! LENGTH GAUGE (electric field) WAS FT'ed , OK with usual formula multiply by wfi
!! UNITLESS RESPONSE FUNCTION FOR ABSORPTION/EMISSION 2 omega im(D(omega)E(omega)^*) IN COLUMN 7
!! QUANTUM MECHANICAL PHOTOABSORPTION/EMISSION CROSS SECTION IN MEGABARNS (no factor of 1/3) IN COLUMN NINE
!! INTEGRATED DIFFERENTIAL OSCILLATOR STRENGTH (FOR SUM RULE) IN COLUMN 10
!! WORK DONE BY EACH PULSE (INTEGRATED) IN COLUMNS 11....

        xsecunits = 5.291772108d0**2 * 4d0 * PI / 1.37036d2 * myenergy

!! NOW FACTOR (2 omega) IN COLUMNS 6,7   v1.16 12-2015

        write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat)  myenergy, &
             fftrans(i), eft(i), fftrans(i)*conjg(eft(i)) * 2 * myenergy, &
             fftrans(i)*conjg(eft(i)) / abs(eft(i)**2) * xsecunits, xsum, worksums(:)
     enddo
     call checkiostat(myiostat,"writing "//outftname)
     close(171)


!!  NUMBER OF PHOTONS ABSORBED AND AND WORK DONE BY EACH PULSE
!!  worksum0 the time integral converges right after pulse is finished... others take longer

     OFL;write(mpifileptr,'(A27,A3,400F15.10)') "     EACH PULSE WORK", xlabel(which), worksum0(:); CFL
     OFL;write(mpifileptr,'(A27,A3,400F15.10)') "     WORK EACH PULSE", xlabel(which), worksums(:); CFL
     OFL;write(mpifileptr,'(A27,A3,400F15.10)') "  PHOTONS EACH PULSE", xlabel(which), exsums(:); CFL

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
        write(171,'(A120)') "## WORK DONE BY EACH PULSE (INTEGRATED) IN COLUMNS 11...."
        write(171,*)

        xsum=0d0
        worksums(:)=0d0

        do i=0,numdata
           myenergy=i*Estep

!! xsum sums to N for N electrons
           if (myenergy.ge.dipolesumstart.and.myenergy.le.dipolesumend) then
              xsum=xsum + Estep * imag(fftrans(i)*conjg(eft(i))) / abs(eft(i)**2) * myenergy * 2 / PI
              worksums(:)=worksums(:) + Estep * imag(fftrans(i)*conjg(all_eft(i,:))) / PI * myenergy
           endif

           xsecunits = 5.291772108d0**2 * 4d0 * PI / 1.37036d2 * myenergy

!! NOW FACTOR (2 omega) IN COLUMNS 6,7   v1.16 12-2015

           write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat)  myenergy, &
                fftrans(i), eft(i), fftrans(i)*conjg(eft(i)) * 2 * myenergy, &
                fftrans(i)*conjg(eft(i)) / abs(eft(i)**2) * xsecunits, xsum, worksums(:)
        enddo
        call checkiostat(myiostat,"writing "//outftname)
        close(171)
     endif
  endif

  deallocate(fftrans,eft,all_eft,dipole_diff)

end subroutine dipolecall

