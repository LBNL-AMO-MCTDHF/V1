
#include "Definitions.INC"



subroutine vectdpot(myintime,invelflag,tdpotsout)
  use parameters
  implicit none
  real*8,intent(in) :: myintime
  integer, intent(in) :: invelflag
  DATATYPE,intent(out) :: tdpotsout(3)
  DATATYPE :: tdpotlen,tdpotvel

  if (invelflag.eq.0) then
     tdpotsout(1)=tdpotlen(myintime,1)  ;  tdpotsout(2)=tdpotlen(myintime,2);  tdpotsout(3)=tdpotlen(myintime,3)
  else
     tdpotsout(1)=tdpotvel(myintime,1)  ;  tdpotsout(2)=tdpotvel(myintime,2);  tdpotsout(3)=tdpotvel(myintime,3)
  endif

end subroutine vectdpot




function tdpotlen(myintime, which)
  use parameters
  implicit none
  integer :: which, ipulse
  real*8 :: myintime,fac
  DATATYPE :: tdpotlen, simplepulselen,pulselen, longpulselen,cwpulselen

  tdpotlen=0.d0

  do ipulse=1,numpulses

     if (which==3) then !! z component
        fac=cos(pulsetheta(ipulse))
     else if (which==1) then
        fac=sin(pulsetheta(ipulse))*cos(pulsephi(ipulse))
     else if (which==2) then
        fac=sin(pulsetheta(ipulse))*sin(pulsephi(ipulse))
     else 
        OFLWR "ACK which = ",which," not allowed tdpotlen"; CFLST
        fac=798d0   !! avoid warn unused
     endif

     if (fac.ne.0d0) then
        select case (pulsetype(ipulse))
        case (1)
           tdpotlen=tdpotlen+simplepulselen(myintime,ipulse) * fac
        case (2)
           tdpotlen=tdpotlen+pulselen(myintime,ipulse) * fac
        case (3)
           tdpotlen=tdpotlen+longpulselen(myintime,ipulse) * fac
        case (4)
           tdpotlen=tdpotlen+cwpulselen(myintime,ipulse) * fac
        case default
           OFLWR "Pulse type not supported: ", pulsetype(ipulse); CFLST
        end select
     endif
  enddo
  
end function tdpotlen


function tdpotvel(myintime,which)
  use parameters

  implicit none
  integer :: which, ipulse
  real*8 :: myintime,fac
  DATATYPE :: tdpotvel, simplepulsevel, pulsevel, longpulsevel,cwpulsevel

  tdpotvel=0.d0

  do ipulse=1,numpulses

     if (which==3) then !! z component
        fac=cos(pulsetheta(ipulse))
     else if (which==1) then
        fac=sin(pulsetheta(ipulse))*cos(pulsephi(ipulse))
     else if (which==2) then
        fac=sin(pulsetheta(ipulse))*sin(pulsephi(ipulse))
     else 
        OFLWR "ACK which = ",which," not allowed tdpotvel"; CFLST
        fac=798d0   !! avoid warn unused
     endif

     if (fac.ne.0d0) then
        select case (pulsetype(ipulse))
        case (1)
           tdpotvel=tdpotvel+simplepulsevel(myintime, ipulse) * fac
        case (2)
           tdpotvel=tdpotvel+pulsevel(myintime, ipulse) * fac
        case (3)
           tdpotvel=tdpotvel+longpulsevel(myintime, ipulse) * fac
        case (4)
           tdpotvel=tdpotvel+cwpulsevel(myintime, ipulse) * fac
        case default
           OFLWR "Pulse type not supported: ", pulsetype(ipulse); CFLST
        end select
     endif
  enddo

end function tdpotvel



function simplepulselen(myintime, ipulse)
  use parameters
  implicit none
  integer :: ipulse
  real*8 :: time, myintime
  DATATYPE :: simplepulselen

  simplepulselen=0.d0

  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)
     if (time.le.pi/omega(ipulse)) then

        simplepulselen = pulsestrength(ipulse) * omega(ipulse) * &
             ( sin(time*omega(ipulse))*cos(time*omega(ipulse) + phaseshift(ipulse)) + sin(time*omega(ipulse) + phaseshift(ipulse))*cos(time*omega(ipulse)) )

     endif
  endif
end function simplepulselen


function simplepulsevel(myintime, ipulse)
  use parameters
  implicit none
  integer :: ipulse
  real*8 :: time, myintime
  DATATYPE :: simplepulsevel

  simplepulsevel=0.d0

  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)
     if (time.le.pi/omega(ipulse)) then

        simplepulsevel = pulsestrength(ipulse) * sin(time*omega(ipulse)) * sin(time*omega(ipulse) + phaseshift(ipulse) )

     endif
  endif

end function simplepulsevel



function cwpulselen(myintime, ipulse)
  use parameters
  implicit none
  integer :: ipulse
  real*8 :: myintime
  DATATYPE :: cwpulselen

  cwpulselen=0d0
  if (myintime.lt.pi/omega(ipulse)) then
     cwpulselen = pulsestrength(ipulse) * omega2(ipulse)*cos(myintime*omega2(ipulse)+phaseshift(ipulse))
  endif

end function cwpulselen


function cwpulsevel(myintime, ipulse)
  use parameters
  implicit none
  integer :: ipulse
  real*8 :: myintime
  DATATYPE :: cwpulsevel

  cwpulsevel=0.d0
  if (myintime.lt.pi/omega(ipulse)) then
     cwpulsevel = pulsestrength(ipulse) * sin(myintime*omega2(ipulse) + phaseshift(ipulse))
  endif

end function cwpulsevel


function pulselen(myintime, ipulse)
  use parameters
  implicit none
  integer :: ipulse
  real*8 :: time, myintime
  DATATYPE :: pulselen

  pulselen=0.d0

  if (chirp(ipulse).ne.0d0.or.ramp(ipulse).ne.0d0) then
     OFLWR "Chirp / ramp not supported for length", chirp(ipulse), ipulse; CFLST
  endif
  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)
     if (time.le.pi/omega(ipulse)) then
        pulselen = pulsestrength(ipulse) * ( &
             2*omega(ipulse)*sin(time*omega(ipulse))*cos(time*omega(ipulse)) * sin((time-pi/omega(ipulse)/2)*omega2(ipulse) + phaseshift(ipulse)) &
             + sin(time*omega(ipulse))**2 * omega2(ipulse) * cos((time-pi/omega(ipulse)/2)*omega2(ipulse) + phaseshift(ipulse)) )
     endif
  endif

end function pulselen


function pulsevel(myintime, ipulse)
  use parameters
  implicit none
  integer :: ipulse
  real*8 :: time, myintime, thisomega2
  DATATYPE :: pulsevel, thisstrength

  pulsevel=0.d0

  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)

!!# notice denominator here , half the total pulse time: energy at start of pulse 
!!# is omega2-chirp; at 1/4 pulse omega2-chirp/2; at 3/4 omega2+chirp/2; end of pulse 
!!# omega2+chirp.  Therefore with given value of chirp, will span this range of 
!!# energies over FWHM.

     thisomega2=omega2(ipulse)+chirp(ipulse)*(time-pi/omega(ipulse)/2)/(pi/omega(ipulse)/2)   /2

     thisstrength=pulsestrength(ipulse)*(1d0+ramp(ipulse)*(time-pi/omega(ipulse)/2)/(pi/omega(ipulse)/2))

     if (time.le.pi/omega(ipulse)) then
        pulsevel = thisstrength * sin(time*omega(ipulse))**2 * sin((time-pi/omega(ipulse)/2)*thisomega2 + phaseshift(ipulse))
     endif
  endif

end function pulsevel

!! wtf with longstep...  why did I do it 2* longstep+1?  So for longstep 0, no constant part of
!!  envelope; for longstep 1, constant part is middle 2/3 of pulse.

function longpulselen(myintime, ipulse)
  use parameters
  implicit none
  integer :: ipulse
  real*8 :: time, myintime, fac, fac2
  DATATYPE :: longpulselen

  longpulselen=0.d0

  if (chirp(ipulse).ne.0d0.or.ramp(ipulse).ne.0d0) then
     OFLWR "Chirp not supported for length", chirp(ipulse), ipulse; CFLST
  endif
  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)

     if (time.le.pi/omega(ipulse)) then

        if ( (time.le.pi/2.d0/omega(ipulse)/(2*longstep(ipulse)+1)) .or. (time.ge.pi/omega(ipulse) - pi/2.d0/omega(ipulse)/(2*longstep(ipulse)+1)) ) then
           fac=2*omega(ipulse)*(2*longstep(ipulse)+1)*sin(time*omega(ipulse)*(2*longstep(ipulse)+1))*cos(time*omega(ipulse)*(2*longstep(ipulse)+1))
           fac2=sin(time*omega(ipulse)*(2*longstep(ipulse)+1))**2
        else
           fac=0.d0
           fac2=1.d0
        endif

        longpulselen = pulsestrength(ipulse) * ( &
             fac * sin(time*omega2(ipulse) + phaseshift(ipulse)) &
          + fac2 * omega2(ipulse) * cos(time*omega2(ipulse) + phaseshift(ipulse)) )
     endif
  endif

end function longpulselen


function longpulsevel(myintime, ipulse)
  use parameters
  implicit none
  integer :: ipulse
  real*8 :: time, myintime, thisomega2
  DATATYPE :: longpulsevel,thisstrength

  longpulsevel=0.d0

  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)

     if (time.le.pi/omega(ipulse)) then

!! this is right I think  goes to chirp/2  when 2*logstep+1 / 4*longstep+2 -way before half time

        thisomega2=omega2(ipulse)+chirp(ipulse)/2 *(time-pi/omega(ipulse)/2)/(pi/omega(ipulse)/2*(2*longstep(ipulse)+1)/(4*longstep(ipulse)+2))
        thisstrength=pulsestrength(ipulse)*(1d0 +ramp(ipulse) *(time-pi/omega(ipulse)/2)/(pi/omega(ipulse)/2))

!!(2*longstep(ipulse)+1)/(4*longstep(ipulse)+2)))

        if ( (time.le.pi/2.d0/omega(ipulse)/(2*longstep(ipulse)+1)) .or. (time.ge.pi/omega(ipulse) - pi/2.d0/omega(ipulse)/(2*longstep(ipulse)+1)) ) then
          longpulsevel = thisstrength  * sin((time-pi/omega(ipulse)/2)*thisomega2 + phaseshift(ipulse)) * sin(time*omega(ipulse)*(2*longstep(ipulse)+1))**2 
        else
           longpulsevel = thisstrength * sin((time-pi/omega(ipulse)/2)*thisomega2 + phaseshift(ipulse))
        endif
     endif
  endif
  
end function longpulsevel
