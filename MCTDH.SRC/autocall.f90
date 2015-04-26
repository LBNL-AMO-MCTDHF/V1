
!! AUTOCORRELATION SUBROUTINES.

#include "Definitions.INC"

module corrmod
  implicit none
  DATATYPE, allocatable ::  wsave(:), fftrans(:,:),fftrans0(:,:), fftranssave(:,:)
end module corrmod

subroutine ftalloc(insize)
  use corrmod
  use parameters
  implicit none
  integer :: insize

  autosize=insize
  allocate( wsave(9*autosize+15), fftrans(-autosize:autosize,mcscfnum), fftranssave(-autosize:autosize,mcscfnum),fftrans0(-autosize:autosize,mcscfnum))

  fftrans=0.d0;  fftrans0=0.d0;  fftranssave=0.d0
end subroutine ftalloc


subroutine ftdealloc()
  use corrmod
  implicit none
  deallocate( wsave, fftrans, fftranssave,fftrans0)
end subroutine ftdealloc

subroutine autocall(mynumdata, forwardovl, backwardovl,sflag)
  use parameters
  use mpimod
  use corrmod
  implicit none

  complex*16 :: tdpotft
  DATATYPE :: forwardovl(0:autosize,mcscfnum), backwardovl(0:autosize,mcscfnum), ffsum(mcscfnum)
  integer :: i, numdata, totdim,  mynumdata, status,imc,sflag,getlen
  real*8 :: piover2,  estep, thistime, myft, myenergy,  ftsum(mcscfnum), twopi
  character (len=200) :: syscommand
  character (len=7) :: number
  REAL*8, ALLOCATABLE :: hanningenergy(:)
  COMPLEX*16, ALLOCATABLE :: hanningft(:,:), hanningwork(:), hanningfftinp(:)

  piover2=2.d0 * atan(1.d0);  twopi=4*piover2

  !! offset for fftrans!

  numdata=mynumdata
  if (numdata.gt.autosize) then
     OFLWR "Err numdata, autosize!  ", numdata, autosize
     WRFL " Setting numdata=autosize and continuing."; CFLST
     numdata=autosize
  endif

  IF(hanningflag .NE. 3 .AND. hanningflag .NE. 4) THEN
     if (abs(conjg(backwardovl(0,1)+(0.d0,0.d0))-forwardovl(0,1)).gt.1.d-10) then
        OFLWR "Looks like forwards and backwards are from different calculations:"
        WRFL backwardovl(0,1), forwardovl(0,1); CFLST
     endif
  ENDIF

  totdim=2*numdata+1

  IF(hanningflag .EQ. 4) THEN
     totdim=numdata+1
  ENDIF

#ifdef REALGO
  OFLWR "Error, corrflag not supported for relax"; CFLST
#else

  IF(hanningflag .NE. 4) THEN
     do i=0,numdata
        fftrans(i,:) = forwardovl(i,:)           * cos(real(i,8)/real(numdata,8) * piover2) * exp((0.d0,-1.d0)*ceground*par_timestep*autosteps*i)
        fftrans0(i,:) = forwardovl(i,:)          * exp((0.d0,-1.d0)*ceground*par_timestep*autosteps*i)
     enddo
     
     do i=0,numdata
        fftrans(-i,:) = conjg(backwardovl(i,:))  * cos(real(i,8)/real(numdata,8) * piover2) * exp((0.d0,1.d0)*ceground*par_timestep*autosteps*i)
        fftrans0(-i,:) = conjg(backwardovl(i,:))  * exp((0.d0,1.d0)*ceground*par_timestep*autosteps*i)
     enddo
  ENDIF
  
! Multiply the Hanning Window and the autocorrelation function
! Note the autocorrelation function is given as the complex conjugate of the literature definition.
! This is fixed.

  IF(hanningflag .EQ. 4) THEN
     do i=0,numdata
        fftrans(i,:) = conjg(forwardovl(i,:))  *(1-cos(twopi*real(i,8)/real(numdata,8)))* exp((0.d0,-1.d0)*ceground*par_timestep*autosteps*i)
        fftrans0(i,:) = conjg(forwardovl(i,:)) * exp((0.d0,-1.d0)*ceground*par_timestep*autosteps*i)
     enddo
  ENDIF
#endif
  

!Specify the dimension of the transform

  if (myrank.eq.1 .AND. hanningflag .NE. 4) THEN
     open(171,file=corrdatfile,status="unknown")
     write(171,*) "#   ", 2*numdata+1
     
     do i=-numdata,numdata
        write(171,'(F18.12, T22, 400E20.8)')  i*par_timestep*autosteps, (fftrans0(i,imc),fftrans(i,imc),imc=1,mcscfnum)
     enddo
     close(171)
  ELSE IF(myrank.eq.1 .AND. hanningflag .EQ. 4) THEN
     open(171,file=corrdatfile,status="unknown")
     write(171,*) "#   ", numdata
     
     do i=0,numdata
        write(171,'(F18.12, T22, 400E20.8)')  i*par_timestep*autosteps, (fftrans0(i,imc),fftrans(i,imc),imc=1,mcscfnum)
     enddo
     close(171)
  ENDIF

!
!BRANT'S FFT using the Hanning Window.
  IF(hanningflag .EQ. 4)  THEN
     ALLOCATE(hanningenergy(1:numdata), hanningwork(1:numdata*4+19), hanningft(1:numdata,mcscfnum), hanningfftinp(1:numdata), STAT=status)
     IF(STATUS .NE. 0) THEN
        WRITE(222, *) "FAIL";        STOP
     ENDIF
  else
     ALLOCATE(hanningenergy(1), hanningwork(1), hanningft(1,1), hanningfftinp(1))
  endif

  do imc=1,mcscfnum
     IF(hanningflag .EQ. 4)  THEN
        DO i=1, numdata
           hanningfftinp(i) = fftrans(i,imc)
        ENDDO
!
!Perform FFT
!
        CALL ZFFTI(numdata, hanningwork)
        CALL ZFFTB(numdata, hanningfftinp, hanningwork)
        !
!Scale and translate the parts of the FFT correctly
!
        DO i=1, numdata
           IF(i .LE. (numdata)/2) THEN
              hanningenergy(i+(numdata)/2) = (i-1)*twopi/(par_timestep*autosteps*numdata)
              hanningft(i+(numdata)/2,imc) = hanningfftinp(i)*(par_timestep*autosteps)
           ELSE
              hanningenergy(i-(numdata)/2) = (i-1-numdata)*twopi/(par_timestep*autosteps*numdata)
              hanningft(i-(numdata)/2,imc) = hanningfftinp(i)*(par_timestep*autosteps)
           ENDIF
        ENDDO
     
!
! Write to hanningout.out
!
     ENDIF

     IF(hanningflag .NE. 4) THEN
        if (myrank.eq.1) then
#ifdef REALGO
           call dffti(totdim,wsave)
           call dfftf(totdim,fftrans(-numdata,imc), wsave)
#else
           call zffti(totdim,wsave)
           call zfftf(totdim,fftrans(-numdata,imc), wsave)
#endif
           fftrans(-numdata:numdata,imc)=fftrans(-numdata:numdata,imc)*autosteps*par_timestep
        endif
        call mympibcast(fftrans(-numdata,imc), 1, 2*numdata+1)
     ELSE
        if (myrank.eq.1) then
           call zffti(totdim,wsave)
           call zfftf(totdim,fftrans(0,imc), wsave)  
        endif
        call mympibcast(fftrans(-numdata,imc), 1, 2*numdata+1)
     ENDIF
  enddo   !! IMC


  Estep=4*piover2/par_timestep/autosteps/totdim

  thistime=numdata*par_timestep*autosteps
  syscommand=nullbuff

!! NEW 0603
! Edited to make sure the no errors for corrflag=4.  Note the output in the Correlation file is NOT CORRECT! 
! 2014 - still so???

  IF(hanningflag .NE. 4) THEN
     do i=-numdata,numdata
        fftrans(i,:)=fftrans(i,:)*exp((0.d0,1.d0)*(numdata+i)*numdata*2*pi/real(2*numdata+1))
     enddo
  ENDIF
  
  IF(hanningflag .NE. 4) THEN
     do i=1,numdata
        fftranssave(i-numdata-1,:) = fftrans(i,:)
     enddo
     do i=-numdata,0
        fftranssave(numdata+i,:)=fftrans(i,:)
     enddo
  ELSE
     do i=numdata/2+1,numdata
        fftranssave(i-((numdata)/2-1),:) = fftrans(i,:)
     enddo
     do i=0, numdata/2
        fftranssave((numdata)/2+i,:)=fftrans(i,:)
     enddo
  ENDIF
  
  ffsum(:)=0.d0
  do i=(-numdata),numdata
     ffsum(:)=ffsum(:)+fftranssave(i,:)
  enddo
  ffsum=ffsum*Estep
  
!!  OFL; write(mpifileptr,'(A12,100F15.10)') "AUTO SUM: ", ffsum(:)/2/pi; CFL
  
  if (myrank.eq.1) then
     IF(hanningflag .NE. 4) THEN

        if (sflag.ne.0) then
           write(number,'(I7)') 1000000+floor(thistime)
           open(1711,file=corrftfile(1:getlen(corrftfile)-1)//number(2:7),status="unknown")
           write(1711,*) "#   ", numdata
           do i=-numdata,0
              write(1711,'(F18.12, T22, 400E20.8)')  i*Estep, fftranssave(i,:)
           enddo
        endif
        open(171,file=corrftfile,status="unknown")
        write(171,*) "#   ", numdata
        do i=-numdata,0
           write(171,'(F18.12, T22, 400E20.8)')  i*Estep, fftranssave(i,:)
        enddo
        do i=1,numdata
           myenergy=i*Estep
           
           myft=0.25d0 * abs(tdpotft(myenergy)**2)*(myenergy)
           
           if ((myft.gt.0.d0).and.(myft+1.d0.ne.myft)) then
              ftsum(:)=abs(fftranssave(i,:))/myft
           else
              ftsum(:)=0.d0
           endif

!!           write(171,'(F18.12, T22, 400E20.8)')  myenergy, abs(fftranssave(i)), fftranssave(i), ftsum, myft
           write(171,'(F18.12, T22, 400E20.8)')  myenergy, myft, (abs(fftranssave(i,imc)), fftranssave(i,imc), ftsum(imc),imc=1,mcscfnum)
           if (sflag.ne.0) then
              write(1711,'(F18.12, T22, 400E20.8)')  myenergy, myft, (abs(fftranssave(i,imc)), fftranssave(i,imc), ftsum(imc),imc=1,mcscfnum)
           endif
        enddo
        close(171)
        if (sflag.ne.0) then
           close(1711)
        endif
     ELSE
        open(171,file=corrftfile,status="unknown")
        write(171,*) "#   ", numdata

        do i=0, numdata/2
           write(171,'(F18.12, T22, 400E20.8)')  i*Estep,( abs(fftranssave(i,imc)), fftranssave(i,imc),imc=1,mcscfnum)
        enddo
        
        do i=numdata/2+1, numdata
           myenergy=i*Estep
           myft=0.25d0 * abs(tdpotft(myenergy)**2)*(myenergy)
           
           if ((myft.gt.0.d0).and.(myft+1.d0.ne.myft)) then
              ftsum(:)=abs(fftranssave(i,:))/myft
           else
              ftsum(:)=0.d0
           endif
           
           write(171,'(F18.12, T22, 400E20.8)')  myenergy, ftsum, myft, (abs(fftranssave(i,imc)), fftranssave(i,imc),imc=1,mcscfnum)
        enddo
        close(171)

        OPEN(991, file='hanningout.out', status='unknown')
        DO i=1, numdata
           WRITE(991, '(F18.12, 3E20.8)') (hanningenergy(i), ABS(hanningft(i,imc)), DREAL(hanningft(i,imc)), DIMAG(hanningft(i,imc)),imc=1,mcscfnum)
        ENDDO
        CLOSE(991)
     ENDIF
  endif

  DEALLOCATE(hanningenergy, hanningft, hanningwork, hanningfftinp)

end subroutine autocall

