
#include "Definitions.INC"


module newpulsesubmod
contains

  ! function [aField,eField,dField] = GetField(calcTimes,fieldStrength,omega,duration,phaseShift,doplot,ENV_PWR)

  subroutine NewField(aField,eField,calcTime,omega,duration,phaseShift,envdernum,envpwr)
    use fileptrmod
    implicit none
    real*8, intent(in) :: calcTime, omega, duration, phaseShift
    integer, intent(in) :: envdernum, envpwr
    real*8, intent(out) :: aField,eField
    real*8 :: nothing1, nothing2, nothing3

    if (envdernum == 0) then
       call FieldFunction(calcTime,omega,duration,phaseShift,envpwr,nothing1,aField,eField,nothing2,nothing3)
       aField = aField / omega 
       eField = eField / omega 
    elseif (envdernum == 1) then
       call FieldFunction(calcTime,omega,duration,phaseShift,envpwr,nothing1,nothing2,aField,eField,nothing3)
       aField = aField / omega**2 
       eField = eField / omega**2
    elseif (envdernum == 2) then
       call FieldFunction(calcTime,omega,duration,phaseShift,envpwr,nothing1,nothing2,nothing3,aField,eField)
       aField = aField / omega**3 
       eField = eField / omega**3
    else
       OFLWR "envdernum Not supported ", envdernum; CFLST
    endif

    ! fac = sqrt(2*pi*fieldStrength) / omega / sqrt(duration);

    ! aField = aField * fac ;
    ! eField = eField * fac ;

  end subroutine NewField


  subroutine FieldFunction(intime,omega,duration,phaseShift,N,y,yp,ypp,yp3,yp4)
    use constant_parameters
    implicit none
    real*8, intent(in) :: intime,omega,duration,phaseShift
    integer, intent(in) :: N
    real*8, intent(out) :: y,yp,ypp,yp3,yp4
    real*8 :: O, funt, dfunt, ddfunt, d3funt, d4funt, env, denv, ddenv, d3env, d4env
    real*8 :: env0, denv0, ddenv0, d3env0, d4env0
    real*8 :: t

    t = intime;

    O = pi/duration;

    y = 0;
    yp = 0;
    ypp = 0;
    yp3 = 0;
    yp4 = 0;

    if (t<=0 .or. t>=duration) then
       return;
    endif

    ! NOW CENTERED
    t = t - duration / 2;

    !    fun   = @sin;
    !    dfun  = @(x) cos(x);
    !    ddfun = @(x) -sin(x);
    !    d3fun = @(x) -cos(x);
    !    d4fun = @(x) sin(x);
    !
    !    funt = fun(omega*t+phaseShift);
    !    dfunt = omega * dfun(omega*t+phaseShift);
    !    ddfunt = omega**2 * ddfun(omega*t+phaseShift);
    !    d3funt = omega**3 * d3fun(omega*t+phaseShift);
    !    d4funt = omega**4 * d4fun(omega*t+phaseShift);
    
    funt = sin(omega*t+phaseShift);
    dfunt = omega * cos(omega*t+phaseShift);
    ddfunt = -omega**2 * sin(omega*t+phaseShift);
    d3funt = -omega**3 * cos(omega*t+phaseShift);
    d4funt = omega**4 * sin(omega*t+phaseShift);

    if (1 == 0) then

       ! %$$ N = 3;   % N=3 is minimal value ensuring E is zero at t=0.
       ! %$$          % with N=3, d/dt E is not zero at t=0.

       ! N = 5;

       env = cos(O*t)**N;
       denv = -O * N * cos(O*t)**(N-1) * sin(O*t);
       ddenv =   &
            O**2 * N * (N-1) * cos(O*t)**(N-2) * sin(O*t)**2 &
            - O**2 * N * cos(O*t)**N ;

       ! - O**2 * N * cos(O*t)**(N-1) * cos(O*t) ;

       d3env =   &
            - O**3 * N * (N-1) * (N-2) * cos(O*t)**(N-3) * sin(O*t)**3 &
            + O**3 * N * (3*N-2)       * cos(O*t)**(N-1) * sin(O*t) ;

       ! + O**3 * N * (N-1) * 2     * cos(O*t)**(N-1) * sin(O*t) &
       ! + O**3 * N**2               * cos(O*t)**(N-1) * sin(O*t) ;

       d4env =   &
            + O**4 * N * (N-1) * (N-2) * (N-3) * cos(O*t)**(N-4) * sin(O*t)**4 &
            - O**4 * N * (N-1) * (6*N-8)   * cos(O*t)**(N-2) * sin(O*t)**2 & 
            + O**4 * N * (3*N-2)           * cos(O*t)**N  ;

       ! - O**4 * N * (N-1) * (N-2) * 3 * cos(O*t)**(N-2) * sin(O*t)**2 &
       ! - O**4 * N * (N-1) * (3*N-2)   * cos(O*t)**(N-2) * sin(O*t)**2 ;

    else

       ! N = 3;

       env0     = 1-(O*2/pi)**2*t**2;
       denv0    = -2*(O*2/pi)**2*t;
       ddenv0   = -2*(O*2/pi)**2*(t*0+1);
       d3env0   = t*0;
       d4env0   = t*0;


       env     = env0**N;
       denv    = N * env0**(N-1) * denv0;
       ddenv   = N * (N-1) * env0**(N-2) * denv0**2 + N * env0**(N-1) * ddenv0;
       
       !d3env   = N * (N-1) * (N-2) * env0**(N-3) * denv0**3 &
       !     + 2 * N * (N-1) * env0**(N-2) * denv0 * ddenv0 &
       !     +     N * (N-1) * env0**(N-2) * denv0 * ddenv0 &
       !     + N * env0**(N-1) * d3env0;

       d3env   = N * (N-1) * (N-2) * env0**(N-3) * denv0**3 &
            + 3 * N * (N-1) * env0**(N-2) * denv0 * ddenv0 &
            + N * env0**(N-1) * d3env0;

       !d4env   = N * (N-1) * (N-2) * (N-3) * env0**(N-4) * denv0**4 &
       !     + 3 * N * (N-1) * (N-2) * env0**(N-3) * ddenv0 * denv0**2 &
       !     + 3 * N * (N-1) * (N-2) * env0**(N-3) * denv0**2 * ddenv0 &
       !     + 3 * N * (N-1) * env0**(N-2) * ddenv0**2 &
       !     + 3 * N * (N-1) * env0**(N-2) * denv0 * d3env0 &
       !     + N * (N-1) * env0**(N-2) * denv0 * d3env0 &
       !     + N * env0**(N-1) * d4env0;

       d4env   = N * (N-1) * (N-2) * (N-3) * env0**(N-4) * denv0**4 &
            + 6 * N * (N-1) * (N-2) * env0**(N-3) * denv0**2 * ddenv0 &
            + N * (N-1) * env0**(N-2) * ( 3 * ddenv0**2 + 3 * denv0 * d3env0  + denv0 * d3env0 ) &
            + N * env0**(N-1) * d4env0;

    end if

    y = &
         env * funt ;

    yp = &
         denv * funt + &
         env * dfunt ;

    ypp = &
         ddenv    * funt + &
         2 * denv * dfunt + &
         env      * ddfunt ;

    yp3 = &
         d3env     * funt + &
         3 * ddenv * dfunt + &
         3 * denv  * ddfunt + &
         env       * d3funt ;

    yp4 = &
         d4env     * funt + &
         4 * d3env * dfunt + &
         6 * ddenv * ddfunt + &
         4 * denv  * d3funt + &
         env       * d4funt ;

  end subroutine FieldFunction

end module newpulsesubmod


module pulsesubmod
contains

subroutine vectdpot(myintime,invelflag,tdpotsout,imc)
  use pulse_parameters
  implicit none
  real*8,intent(in) :: myintime
  integer, intent(in) :: invelflag,imc
  DATATYPE,intent(out) :: tdpotsout(3)

  call vectdpot0(myintime,invelflag,tdpotsout,imc,1,numpulses)

end subroutine vectdpot


subroutine vectdpot0(myintime,invelflag,tdpotsout,imc,ilow,ihigh)
  use pulse_parameters
  use fileptrmod
  implicit none
  integer, intent(in) :: invelflag,imc,ilow,ihigh
  real*8,intent(in) :: myintime
  DATATYPE,intent(out) :: tdpotsout(3)

  if (invelflag.eq.0) then
     tdpotsout(1)=tdpotlen0(myintime,1,ilow,ihigh)
     tdpotsout(2)=tdpotlen0(myintime,2,ilow,ihigh)
     tdpotsout(3)=tdpotlen0(myintime,3,ilow,ihigh)
  else
     tdpotsout(1)=tdpotvel0(myintime,1,ilow,ihigh)
     tdpotsout(2)=tdpotvel0(myintime,2,ilow,ihigh)
     tdpotsout(3)=tdpotvel0(myintime,3,ilow,ihigh)
  endif

!  if (conjgpropflag.ne.0) then
!     select case (imc)
!     case(-1)
!        tdpotsout(:)=real(tdpotsout(:),8)
!     case(1)
!     case(2)
!        tdpotsout(:)=ALLCON(tdpotsout(:))
!     case default
!        OFLWR "WHOOPS?? conjgprop",imc; CFLST
!     end select
!  endif

end subroutine vectdpot0


function tdpotlen0(myintime, which,ilow,ihigh)
  use pulse_parameters
  use fileptrmod
  implicit none
  integer,intent(in) :: which,ilow,ihigh
  real*8,intent(in) :: myintime
  integer :: ipulse
  real*8 :: fac
  DATATYPE :: tdpotlen0

  tdpotlen0=0.d0

!!  do ipulse=1,numpulses

  do ipulse=ilow,ihigh

     if (which==3) then !! z component
        fac=cos(pulsetheta(ipulse))
     else if (which==1) then
        fac=sin(pulsetheta(ipulse))*cos(pulsephi(ipulse))
     else if (which==2) then
        fac=sin(pulsetheta(ipulse))*sin(pulsephi(ipulse))
     else 
        OFLWR "ACK which = ",which," not allowed tdpotlen0"; CFLST
        fac=798d0   !! avoid warn unused
     endif

     if (fac.ne.0d0) then
        select case (pulsetype(ipulse))
        case (1)
           tdpotlen0=tdpotlen0+simplepulselen(myintime,ipulse) * fac
        case (2)
           tdpotlen0=tdpotlen0+pulselen(myintime,ipulse) * fac
        !$$  case (3)
        !$$     tdpotlen0=tdpotlen0+longpulselen(myintime,ipulse) * fac
        !$$  case (4)
        !$$     tdpotlen0=tdpotlen0+cwpulselen(myintime,ipulse) * fac
        !$$  case (5)
        !$$    tdpotlen0=tdpotlen0+monopulselen(myintime,ipulse) * fac
        case (6)
           tdpotlen0=tdpotlen0+newpulselen(myintime,ipulse) * fac
        case (7)
           tdpotlen0=tdpotlen0+threepulselen(myintime,ipulse) * fac
        case default
           OFLWR "Pulse type not supported: ", pulsetype(ipulse); CFLST
        end select
     endif
  enddo
  
end function tdpotlen0


function tdpotvel0(myintime,which,ilow,ihigh)
  use pulse_parameters
  use fileptrmod
  implicit none
  integer,intent(in) :: which,ilow,ihigh
  real*8,intent(in) :: myintime
  integer :: ipulse
  real*8 :: fac
  DATATYPE :: tdpotvel0

  tdpotvel0=0.d0

!!  do ipulse=1,numpulses

  do ipulse=ilow,ihigh

     if (which==3) then !! z component
        fac=cos(pulsetheta(ipulse))
     else if (which==1) then
        fac=sin(pulsetheta(ipulse))*cos(pulsephi(ipulse))
     else if (which==2) then
        fac=sin(pulsetheta(ipulse))*sin(pulsephi(ipulse))
     else 
        OFLWR "ACK which = ",which," not allowed tdpotvel0"; CFLST
        fac=798d0   !! avoid warn unused
     endif

     if (fac.ne.0d0) then
        select case (pulsetype(ipulse))
        case (1)
           tdpotvel0=tdpotvel0+simplepulsevel(myintime, ipulse) * fac
        case (2)
           tdpotvel0=tdpotvel0+pulsevel(myintime, ipulse) * fac
        !$$  case (3)
        !$$     tdpotvel0=tdpotvel0+longpulsevel(myintime, ipulse) * fac
        !$$  case (4)
        !$$     tdpotvel0=tdpotvel0+cwpulsevel(myintime, ipulse) * fac
        !$$  case (5)
        !$$     tdpotvel0=tdpotvel0+monopulsevel(myintime, ipulse) * fac
        case (6)
           tdpotvel0=tdpotvel0+newpulsevel(myintime, ipulse) * fac
        case (7)
           tdpotvel0=tdpotvel0+threepulsevel(myintime, ipulse) * fac
        case default
           OFLWR "Pulse type not supported: ", pulsetype(ipulse); CFLST
        end select
     endif
  enddo

end function tdpotvel0





function newpulselen(myintime, ipulse)
  use pulse_parameters
  use constant_parameters
  use newpulsesubmod
  implicit none
  integer,intent(in) :: ipulse
  real*8,intent(in) :: myintime
  real*8 :: time, duration, aField, eField
  DATATYPE :: newpulselen

  newpulselen=0.d0

  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)
     duration = pi/omega(ipulse)
     if (time.le.duration) then
        call NewField(aField,eField,time,omega2(ipulse),duration,phaseshift(ipulse),envdernum,envpwr)
        newpulselen = pulsestrength(ipulse) * eField;
     endif
  endif
  
end function newpulselen



function newpulsevel(myintime, ipulse)
  use pulse_parameters
  use constant_parameters
  use newpulsesubmod
  implicit none
  integer,intent(in) :: ipulse
  real*8,intent(in) :: myintime
  real*8 :: time, duration, aField, eField
  DATATYPE :: newpulsevel

  newpulsevel=0.d0

  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)
     duration = pi/omega(ipulse)
     if (time.le.duration) then
        call NewField(aField,eField,time,omega2(ipulse),duration,phaseshift(ipulse),envdernum,envpwr)
        newpulsevel = pulsestrength(ipulse) * aField;
     endif
  endif
  
end function newpulsevel




function simplepulselen(myintime, ipulse)
  use pulse_parameters
  use constant_parameters
  implicit none
  integer,intent(in) :: ipulse
  real*8,intent(in) :: myintime
  real*8 :: time
  DATATYPE :: simplepulselen

  simplepulselen=0.d0

  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)
     if (time.le.pi/omega(ipulse)) then

        simplepulselen = pulsestrength(ipulse) * omega(ipulse) * &
             ( sin(time*omega(ipulse))*cos(time*omega(ipulse) + phaseshift(ipulse)) &
             + sin(time*omega(ipulse) + phaseshift(ipulse))*cos(time*omega(ipulse)) )

     endif
  endif
end function simplepulselen


function simplepulsevel(myintime, ipulse)
  use pulse_parameters
  use constant_parameters
  implicit none
  integer,intent(in) :: ipulse
  real*8,intent(in) :: myintime
  real*8 :: time
  DATATYPE :: simplepulsevel

  simplepulsevel=0.d0

  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)
     if (time.le.pi/omega(ipulse)) then

        simplepulsevel = pulsestrength(ipulse) * &
             sin(time*omega(ipulse)) * sin(time*omega(ipulse) + phaseshift(ipulse) )

     endif
  endif

end function simplepulsevel




function threepulsevel(intime, ipulse)
  use pulse_parameters
  use constant_parameters
  implicit none
  integer,intent(in) :: ipulse
  real*8,intent(in) :: intime
  real*8 :: time,dur,pstart,f,df,ddf
  DATATYPE :: pstren, threepulsevel
  integer :: opt1, negflag
  
  opt1    = pulse7opt1(ipulse) ;
  negflag = 0;
  select case(opt1)
  case(0,1)
  case(2,3)
     opt1 = opt1 - 2;
     negflag = 1;
  case default
  end select
  
  threepulsevel       = 0.d0
  dur                 = pulsedur(ipulse)
  pstart              = pulsestart(ipulse)
  pstren              = pulsestrength(ipulse)  
  if (negflag .ne. 0) then
     pstren           = pstren * (-1)  ! if flip time, mult deriv by (-1)
     pstren           = pstren * (-1)  ! flipping sign both length and velocity
  endif
  if (intime.ge.pstart) then
     time             = intime - pstart
     if (time.le.dur) then
        time          = time/dur
        if (negflag .ne. 0) then
           time = 1-time;
        endif
        call threepulsefunction(time,f,df,ddf,opt1)
        threepulsevel = pstren * df ;
     endif
  endif
end function threepulsevel




function threepulselen(intime, ipulse)
  use pulse_parameters
  use constant_parameters
  implicit none
  integer,intent(in) :: ipulse
  real*8,intent(in) :: intime
  real*8 :: time,dur,pstart,f,df,ddf
  DATATYPE :: pstren, threepulselen
  integer :: opt1, negflag
  
  opt1    = pulse7opt1(ipulse) ;
  negflag = 0;
  select case(opt1)
  case(0,1)
  case(2,3)
     opt1 = opt1 - 2;
     negflag = 1;
  case default
  end select
  
  threepulselen       = 0.d0
  dur                 = pulsedur(ipulse)
  pstart              = pulsestart(ipulse)
  pstren              = pulsestrength(ipulse)  
  if (negflag .ne. 0) then
     pstren           = pstren * (-1)  ! flipping sign both length and velocity
  endif
  if (intime.ge.pstart) then
     time             = intime - pstart
     if (time.le.dur) then
        time          = time/dur
        if (negflag .ne. 0) then
           time = 1-time;
        endif
        call threepulsefunction(time,f,df,ddf,opt1)
        threepulselen = pstren * ddf / dur;
     endif
  endif
end function threepulselen


subroutine threepulsefunction(time,f,df,ddf,OPT)
  use fileptrmod
  implicit none
  integer,intent(in) :: OPT
  real*8, intent(in) :: time
  real*8, intent(out) :: f,df,ddf
  real*8 :: aval=0, bval=0, nfac=0, &
       fac1=0,   fac2=0,   faca=0,   facb=0,  &
       dfac1=0,  dfac2=0,  dfaca=0,  dfacb=0,  &
       ddfac1=0, ddfac2=0, ddfaca=0, ddfacb=0
  
  if (time<0 .or. time>1) then
     OFLWR "TIME ERROR ", time; CFLST
  endif

  if (OPT == 0) then
     aval = 1/220d0 * ( 85 - sqrt(185d0) )
     bval = 1/115d0 * ( 45 + sqrt(185d0) )
     fac1 = time**8;
     dfac1 = 8*time**7;
     ddfac1 = 56*time**6;
     fac2 = (1-time)**13;
     dfac2 = -13*(1-time)**12;
     ddfac2 = 156*(1-time)**11;
  elseif (OPT == 1) then
     aval = 1/220d0 * ( 135 - sqrt(185d0) )
     bval = 1/115d0 * ( 70 + sqrt(185d0) )
     fac1 = time**13;
     dfac1 = 13*time**12;
     ddfac1 = 156*time**11;
     fac2 = (1-time)**8;
     dfac2 = -8*(1-time)**7;
     ddfac2 = 56*(1-time)**6;
  else
     OFLWR "programmer error threepulsefunction"; CFLST
  endif
 
  faca = (time-aval)**2;
  dfaca = 2*(time-aval);
  ddfaca = 2;
  facb = (time-bval);
  dfacb = 1;
  ddfacb = 0;

  f =  fac1 * fac2 * faca * facb;
  df = dfac1 * fac2 * faca * facb + &
       fac1 * dfac2 * faca * facb + &
       fac1 * fac2 * dfaca * facb + &
       fac1 * fac2 * faca * dfacb ;
  ddf= ddfac1 * fac2 * faca * facb + &
       dfac1 * dfac2 * faca * facb * 2 + &
       dfac1 * fac2 * dfaca * facb * 2 + &
       dfac1 * fac2 * faca * dfacb * 2 + &
       fac1 * ddfac2 * faca * facb + &
       fac1 * dfac2 * dfaca * facb * 2 + &
       fac1 * dfac2 * faca * dfacb * 2 + &
       fac1 * fac2 * ddfaca * facb + &
       fac1 * fac2 * dfaca * dfacb * 2 + &
       fac1 * fac2 * faca * ddfacb ;

  nfac = 3.069e-7;

  f   = f / nfac;
  df  = df / nfac;
  ddf = ddf / nfac;
  
end subroutine threepulsefunction




!!$function cwpulselen(myintime, ipulse)
!!$  use pulse_parameters
!!$  use constant_parameters
!!$  implicit none
!!$  integer,intent(in) :: ipulse
!!$  real*8,intent(in) :: myintime
!!$  DATATYPE :: cwpulselen
!!$
!!$  cwpulselen=0d0
!!$  if (myintime.lt.pi/omega(ipulse)) then
!!$     cwpulselen = pulsestrength(ipulse) * omega2(ipulse) * &
!!$          cos(myintime*omega2(ipulse)+phaseshift(ipulse))
!!$  endif
!!$
!!$end function cwpulselen
!!$
!!$
!!$function cwpulsevel(myintime, ipulse)
!!$  use pulse_parameters
!!$  use constant_parameters
!!$  implicit none
!!$  integer,intent(in) :: ipulse
!!$  real*8,intent(in) :: myintime
!!$  DATATYPE :: cwpulsevel
!!$
!!$  cwpulsevel=0.d0
!!$  if (myintime.lt.pi/omega(ipulse)) then
!!$     cwpulsevel = pulsestrength(ipulse) * sin(myintime*omega2(ipulse) + phaseshift(ipulse))
!!$  endif
!!$
!!$end function cwpulsevel


!!$function monopulselen(myintime, ipulse)
!!$  use pulse_parameters
!!$  use constant_parameters
!!$  implicit none
!!$  integer,intent(in) :: ipulse
!!$  real*8,intent(in) :: myintime
!!$  real*8 :: pptime
!!$  DATATYPE :: monopulselen
!!$
!!$  pptime=myintime*omega(ipulse) - pi/2
!!$
!!$  if (pptime.lt.pi/2) then
!!$     monopulselen = pulsestrength(ipulse) * (0.75d0*cos(pptime) + cos(3*pptime)/4d0)
!!$  else
!!$     monopulselen = 0
!!$  end if
!!$
!!$end function monopulselen
!!$
!!$
!!$function monopulsevel(myintime, ipulse)
!!$  use pulse_parameters
!!$  use constant_parameters
!!$  implicit none
!!$  integer,intent(in) :: ipulse
!!$  real*8,intent(in) :: myintime
!!$  real*8 :: pptime
!!$  DATATYPE :: monopulsevel
!!$
!!$  pptime=myintime*omega(ipulse) - pi/2
!!$
!!$  if (pptime.lt.pi/2) then
!!$     monopulsevel = pulsestrength(ipulse) * (0.75d0*sin(pptime) + sin(3*pptime)/12d0 + 2d0/3d0)
!!$  else
!!$     monopulsevel = pulsestrength(ipulse) * (4d0/3d0)
!!$  endif
!!$
!!$end function monopulsevel




function pulselen(myintime, ipulse)
  use pulse_parameters
  use constant_parameters
  use fileptrmod
  implicit none
  integer,intent(in) :: ipulse
  real*8,intent(in) :: myintime
  real*8 :: time
  DATATYPE :: pulselen

  pulselen=0.d0

  if (chirp(ipulse).ne.0d0.or.ramp(ipulse).ne.0d0) then
     OFLWR "Chirp / ramp not supported for length", chirp(ipulse), ipulse; CFLST
  endif
  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)
     if (time.le.pi/omega(ipulse)) then
        pulselen = pulsestrength(ipulse) * ( &
             2*omega(ipulse)*sin(time*omega(ipulse))*cos(time*omega(ipulse)) * &
             sin((time-pi/omega(ipulse)/2)*omega2(ipulse) + phaseshift(ipulse)) &
             + sin(time*omega(ipulse))**2 * omega2(ipulse) * &
             cos((time-pi/omega(ipulse)/2)*omega2(ipulse) + phaseshift(ipulse)) )
     endif
  endif

end function pulselen


function pulsevel(myintime, ipulse)
  use pulse_parameters
  use constant_parameters
  implicit none
  integer,intent(in) :: ipulse
  real*8,intent(in) :: myintime
  real*8 :: time,thisomega2
  DATATYPE :: pulsevel, thisstrength

  pulsevel=0.d0

  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)

!!# notice denominator here , half the total pulse time: energy at start of pulse 
!!# is omega2-chirp; at 1/4 pulse omega2-chirp/2; at 3/4 omega2+chirp/2; end of pulse 
!!# omega2+chirp.  Therefore with given value of chirp, will span this range of 
!!# energies over FWHM.

     thisomega2=omega2(ipulse)+chirp(ipulse)*(time-pi/omega(ipulse)/2)/(pi/omega(ipulse)/2)   /2

     thisstrength=pulsestrength(ipulse)*(1d0+ramp(ipulse) * &
          (time-pi/omega(ipulse)/2)/(pi/omega(ipulse)/2))

     if (time.le.pi/omega(ipulse)) then
        pulsevel = thisstrength * sin(time*omega(ipulse))**2 * &
             sin((time-pi/omega(ipulse)/2)*thisomega2 + phaseshift(ipulse))
     endif
  endif

end function pulsevel


!!$!! wtf with longstep...  why did I do it 2* longstep+1?  So for longstep 0, no constant part of
!!$!!  envelope; for longstep 1, constant part is middle 2/3 of pulse.
!!$
!!$function longpulselen(myintime, ipulse)
!!$  use pulse_parameters
!!$  use constant_parameters
!!$  use fileptrmod
!!$  implicit none
!!$  integer,intent(in) :: ipulse
!!$  real*8,intent(in) :: myintime
!!$  real*8 :: time, fac, fac2
!!$  DATATYPE :: longpulselen
!!$
!!$  longpulselen=0.d0
!!$
!!$  if (chirp(ipulse).ne.0d0.or.ramp(ipulse).ne.0d0) then
!!$     OFLWR "Chirp not supported for length", chirp(ipulse), ipulse; CFLST
!!$  endif
!!$  if (myintime.ge.pulsestart(ipulse)) then
!!$     time=myintime-pulsestart(ipulse)
!!$
!!$     if (time.le.pi/omega(ipulse)) then
!!$
!!$        if ( (time.le.pi/2.d0/omega(ipulse)/(2*longstep(ipulse)+1)) .or.&
!!$             (time.ge.pi/omega(ipulse) - pi/2.d0/omega(ipulse)/(2*longstep(ipulse)+1)) ) then
!!$           fac=2*omega(ipulse)*(2*longstep(ipulse)+1)*sin(time*omega(ipulse)*(2*longstep(ipulse)+1)) * &
!!$                cos(time*omega(ipulse)*(2*longstep(ipulse)+1))
!!$           fac2=sin(time*omega(ipulse)*(2*longstep(ipulse)+1))**2
!!$        else
!!$           fac=0.d0
!!$           fac2=1.d0
!!$        endif
!!$
!!$        longpulselen = pulsestrength(ipulse) * ( &
!!$             fac * sin(time*omega2(ipulse) + phaseshift(ipulse)) &
!!$          + fac2 * omega2(ipulse) * cos(time*omega2(ipulse) + phaseshift(ipulse)) )
!!$     endif
!!$  endif
!!$
!!$end function longpulselen
!!$
!!$
!!$function longpulsevel(myintime, ipulse)
!!$  use pulse_parameters
!!$  use constant_parameters
!!$  implicit none
!!$  integer,intent(in) :: ipulse
!!$  real*8,intent(in) :: myintime
!!$  real*8 :: time, thisomega2
!!$  DATATYPE :: longpulsevel,thisstrength
!!$
!!$  longpulsevel=0.d0
!!$
!!$  if (myintime.ge.pulsestart(ipulse)) then
!!$     time=myintime-pulsestart(ipulse)
!!$
!!$     if (time.le.pi/omega(ipulse)) then
!!$
!!$!! this is right I think  goes to chirp/2  when 2*logstep+1 / 4*longstep+2 -way before half time
!!$
!!$        thisomega2=omega2(ipulse)+chirp(ipulse)/2 *(time-pi/omega(ipulse)/2)/&
!!$             (pi/omega(ipulse)/2*(2*longstep(ipulse)+1)/(4*longstep(ipulse)+2))
!!$        thisstrength=pulsestrength(ipulse)*(1d0 +ramp(ipulse) * &
!!$             (time-pi/omega(ipulse)/2)/(pi/omega(ipulse)/2))
!!$
!!$!!(2*longstep(ipulse)+1)/(4*longstep(ipulse)+2)))
!!$
!!$        if ( (time.le.pi/2.d0/omega(ipulse)/(2*longstep(ipulse)+1)) .or. &
!!$             (time.ge.pi/omega(ipulse) - pi/2.d0/omega(ipulse)/(2*longstep(ipulse)+1)) ) then
!!$          longpulsevel = thisstrength  * &
!!$               sin((time-pi/omega(ipulse)/2)*thisomega2 + phaseshift(ipulse)) * &
!!$               sin(time*omega(ipulse)*(2*longstep(ipulse)+1))**2 
!!$        else
!!$           longpulsevel = thisstrength * &
!!$                sin((time-pi/omega(ipulse)/2)*thisomega2 + phaseshift(ipulse))
!!$        endif
!!$     endif
!!$  endif
!!$  
!!$end function longpulsevel


!! not useful because tentmode=1 isn't better than tentmode=0; tentmode=0 default;
!! keeping these subroutines anyway

function get_rtentsum(num,vec)
  use parameters   !! tentmode
  implicit none
  real*8 :: get_rtentsum
  integer,intent(in) :: num
  real*8, intent(in) :: vec(-num:num)
  integer :: itent,jtent

  if (num.lt.0) then
     OFLWR "rtentsumerror", num; CFLST
  elseif (num.eq.0) then
     get_rtentsum=vec(0)
     return
  endif

  if (tentmode.eq.0) then
     itent=0
     jtent=1
  else
     itent=(num+1)/2 * (-1)
     jtent=(num+1)/2
  endif

  get_rtentsum = SUM(vec(-num:itent)) + SUM(vec(jtent:num))

end function get_rtentsum

function get_ctentsum(num,vec)
  use parameters   !! tentmode
  implicit none
  complex*16 :: get_ctentsum
  integer,intent(in) :: num
  complex*16, intent(in) :: vec(-num:num)
  integer :: itent,jtent

  if (num.lt.0) then
     OFLWR "ctentsumerror", num; CFLST
  elseif (num.eq.0) then
     get_ctentsum=vec(0)
     return
  endif

  if (tentmode.eq.0) then
     itent=0
     jtent=1
  else
     itent=(num+1)/2 * (-1)
     jtent=(num+1)/2
  endif

  get_ctentsum = SUM(vec(-num:itent)) + SUM(vec(jtent:num))

end function get_ctentsum
  

end module pulsesubmod
