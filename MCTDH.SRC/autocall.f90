
!! AUTOCORRELATION SUBROUTINES.

#include "Definitions.INC"


!! actually have numdata+1 data points in forwardovl

subroutine autocall(numdata, forwardovl, sflag)
  use parameters
  use mpimod
  implicit none
  integer, intent(in) :: numdata,sflag
  DATATYPE,intent(in) :: forwardovl(0:autosize,mcscfnum)
  integer :: i, totdim, imc,getlen,ibot,myiostat
  real*8 ::   estep, thistime, myenergy, windowfunct
  character (len=7) :: number
  DATATYPE, allocatable :: fftrans(:,:),fftrans0(:,:)
  DATATYPE :: csums(mcscfnum), csums2(mcscfnum)   !! AUTOMATIC

#ifdef REALGO
  OFLWR "Error, autocall not supported for real-valued"; CFLST
#else

  if (numdata.gt.autosize+1) then
     OFLWR "Err numdata, autosize!  ", numdata, autosize; CFLST
  endif

  if (hanningflag.eq.4) then
     ibot=0
     totdim=numdata+1
  else
     ibot=(-numdata)
     totdim=2*numdata+1
  endif

  allocate( fftrans(ibot:numdata,mcscfnum),fftrans0(ibot:numdata,mcscfnum))

  fftrans=0.d0;  fftrans0=0.d0


  IF(hanningflag .NE. 4) THEN
     do i=0,numdata
        fftrans0(i,:) = forwardovl(i,:)  * exp((0.d0,-1.d0)*ceground*par_timestep*autosteps*i)

        fftrans(i,:) = fftrans0(i,:) * windowfunct(i,numdata,1)  !! action 1

     enddo
     do i=1,numdata
        fftrans(-i,:) = conjg(fftrans(i,:))
        fftrans0(-i,:) = conjg(fftrans0(i,:))
     enddo

  else
  
! Multiply the Hanning Window and the autocorrelation function
! Note the autocorrelation function is given as the complex conjugate of the literature definition.
! This is fixed.

     do i=0,numdata
        fftrans0(i,:) = conjg(forwardovl(i,:)) * &
             exp((0.d0,-1.d0)*ceground*par_timestep*autosteps*i)
        fftrans(i,:) = fftrans0(i,:) * (1-cos(pi*2*real(i,8)/real(numdata,8)))
     enddo

  ENDIF

!Specify the dimension of the transform

  if (myrank.eq.1) then
     open(171,file=corrdatfile,status="unknown",iostat=myiostat)
     call checkiostat(myiostat,"opening corrdatfile")
     write(171,*,iostat=myiostat) "#   ", totdim
     call checkiostat(myiostat,"writing corrdatfile")
     do i=ibot,numdata
        write(171,'(F18.12, T22, 400E20.8)')  i*par_timestep*autosteps, &
             (fftrans0(i,imc),fftrans(i,imc),imc=1,mcscfnum)
     enddo
     close(171)
  ENDIF

  do imc=1,mcscfnum

     if (hanningflag.ne.4) then
        call zfftf_wrap_diff(totdim,fftrans(:,imc),ftdiff)
        fftrans(:,imc)=fftrans(:,imc)*autosteps*par_timestep  /  2 / pi
        do i=-numdata,numdata
           fftrans(i,imc)=fftrans(i,imc) * &
                exp((0.d0,1.d0)*(numdata+i)*numdata*2*pi/real(2*numdata+1))
        enddo
!! for default ftwindowpower=2 for action 1
        if (ftdiff.eq.1) then
           fftrans(-numdata+1,:) = fftrans(-numdata+1,:) * 4d0
        elseif (ftdiff.eq.0) then
           fftrans(-numdata,:) = fftrans(-numdata,:) * sqrt(2d0)
        endif
     else
        CALL ZFFTB_wrap(totdim,fftrans(:,imc))
     endif

  enddo   !! IMC

  Estep=2*pi/par_timestep/autosteps/totdim

  if (myrank.eq.1) then

     if (sflag.ne.0) then
        thistime=numdata*par_timestep*autosteps
        write(number,'(I7)') 1000000+floor(thistime)
        open(1711,file=corrftfile(1:getlen(corrftfile))//number(2:7),&
             status="unknown",iostat=myiostat)
        call checkiostat(myiostat,"opening corrftfile")
        write(1711,*,iostat=myiostat) "#   ", totdim
        call checkiostat(myiostat,"opening corrftfile")
        csums(:)=0d0; csums2(:)=0d0
        do i=ibot,numdata
           myenergy=(i-ibot)*Estep
           if (myenergy.ge.dipolesumstart) then !! dipolesumstart only, not sumend
              csums(:)=csums(:) + fftrans(i,:) * Estep
              csums2(:)=csums2(:) + fftrans(i,:) * Estep * myenergy
           endif
           write(1711,'(F18.12, T22, 5000E20.8)')  myenergy, fftrans(i,:), csums, csums2
        enddo
        close(1711)
     endif

     open(171,file=corrftfile,status="unknown",iostat=myiostat)
     call checkiostat(myiostat,"opening corrftfile")
     write(171,*,iostat=myiostat) "#   ", totdim
     call checkiostat(myiostat,"opening corrftfile")
     csums(:)=0d0; csums2=0d0
     do i=ibot,numdata
        myenergy=(i-ibot)*Estep
        if (myenergy.ge.dipolesumstart) then !! dipolesumstart only, not sumend
           csums(:)=csums(:) + fftrans(i,:) * Estep
           csums2(:)=csums2(:) + fftrans(i,:) * Estep * myenergy
        endif
        write(171,'(F18.12, T22, 5000E20.8)')  myenergy, fftrans(i,:), csums(:), csums2(:)
     enddo
     close(171)

  endif

  deallocate( fftrans, fftrans0)

#endif

end subroutine autocall

