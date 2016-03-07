
!! AUTOCORRELATION SUBROUTINES.

#include "Definitions.INC"


subroutine autocall(mynumdata, forwardovl, backwardovl)
  use parameters
!  use automod
  use mpimod
  use corrmod
  implicit none

  complex*16 :: tdpotft
  DATATYPE :: forwardovl(0:autosize), backwardovl(0:autosize), ffsum
  integer :: i, numdata, totdim,  j, mynumdata
  real*8 :: piover2,  estep, thistime, myft, myenergy,  ftsum
  character (len=200) :: syscommand

  piover2=2.d0 * atan(1.d0)

  !! offset for fftrans!

  numdata=mynumdata
  if (numdata.gt.autosize) then
     call openfile()
     write(mpifileptr, *) "Err numdata, autosize!  ", numdata, autosize
     write(mpifileptr,*) " Setting numdata=autosize and continuing."
     call closefile();     call mpistop()
     numdata=autosize
  endif


  if (abs(conjg(backwardovl(0)+(0.d0,0.d0))-forwardovl(0)).gt.1.d-10) then
     call openfile()
     write(mpifileptr, *) "Looks like forwards and backwards are from different calculations:"
     write(mpifileptr, *) backwardovl(0), forwardovl(0)
     call closefile();     call mpistop()
  endif

  totdim=2*numdata+1

#ifdef REALGO
  call openfile()
  write(mpifileptr,*) "Error, corrflag not supported for relax"
  call closefile();  call mpistop()
#else

  do i=0,numdata
     fftrans(i) = forwardovl(i)           * cos(real(i,8)/real(numdata,8) * piover2) * exp((0.d0,-1.d0)*eground*par_timestep*autosteps*i)
  enddo

  do i=0,numdata
     fftrans(-i) = conjg(backwardovl(i))  * cos(real(i,8)/real(numdata,8) * piover2) * exp((0.d0,1.d0)*eground*par_timestep*autosteps*i)
  enddo

#endif

  if (myrank.eq.1) then
     open(171,file="Correlation.Dat",status="unknown")
     write(171,*) "#   ", numdata

     do i=-numdata,numdata
        write(171,'(F18.12, T22, 400E20.8)')  i*par_timestep*autosteps, fftrans(i)
     enddo
     close(171)
  endif

  if (myrank.eq.1) then
#ifdef REALGO
     call dffti(totdim,wsave)
     call dfftf(totdim,fftrans(-numdata), wsave)
#else
     call zffti(totdim,wsave)
     call zfftf(totdim,fftrans(-numdata), wsave)
#endif
     fftrans(-numdata:numdata)=fftrans(-numdata:numdata)*autosteps*par_timestep
  endif

  call mympibcast(fftrans(-numdata), 1, 2*numdata+1)

  Estep=4*piover2/par_timestep/autosteps/totdim

  thistime=numdata*par_timestep*autosteps
  syscommand=nullbuff

     
  numcorr=numcorr+1

  do j=10,2,-1
     fftranssave(:,j)=fftranssave(:,j-1)
  enddo

!! NEW 0603
  do i=-numdata,numdata
     fftrans(i)=fftrans(i)*exp((0.d0,1.d0)*(numdata+i)*numdata*2*pi/real(2*numdata+1))
  enddo

  do i=1,numdata
     fftranssave(i-numdata-1,1) = fftrans(i)
  enddo
  do i=-numdata,0
     fftranssave(numdata+i,1)=fftrans(i)
  enddo

  ffsum=0.d0
  do i=(-numdata),numdata
!!     ffsum=ffsum+abs(fftrans(i))
     ffsum=ffsum+fftranssave(i,1)
  enddo
  ffsum=ffsum*Estep

  call openfile();  write(mpifileptr, *) 
  write(mpifileptr, *) "AUTO SUM: ", ffsum/2/pi
  write(mpifileptr, *) "Timestep, delta-E:  ", autosteps, " x ",par_timestep, " = ", autosteps*par_timestep, Estep
  write(mpifileptr, *);  call closefile()


  if (myrank.eq.1) then

     open(171,file="Corrft.Dat",status="unknown")
     write(171,*) "#   ", numdata

     do i=-numdata,0
!!        write(171,'(F18.12, T22, 400E20.8)')  i*Estep, abs(fftranssave(i,1)), fftranssave(i,1)
        write(171,'(F18.12, T22, 400E20.8)')  i*Estep, fftranssave(i,1)
     enddo

     do i=1,numdata
        myenergy=i*Estep

        if (dipoleopflag==1) then
           if (velflag==0) then
              myft=1.d0/(myenergy)
           else
              myft=(myenergy)
           endif
        else
           myft=0.25d0 * abs(tdpotft(myenergy)**2)*(myenergy)
        endif

        if ((myft.gt.0.d0).and.(myft+1.d0.ne.myft)) then
           ftsum=abs(fftranssave(i,1))/myft
        else
           ftsum=0.d0
        endif

        write(171,'(F18.12, T22, 400E20.8)')  myenergy, abs(fftranssave(i,1)), fftranssave(i,1), ftsum, myft
     enddo
     close(171)
  endif

end subroutine autocall

