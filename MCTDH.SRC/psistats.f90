
#include "Definitions.INC"

subroutine psistats( thistime )
  use parameters
  use mpimod
  use xxxmod
  implicit none

  DATATYPE :: inspfs(spfsize,nspf), inavector(numconfig,numr), dot,sum
  DATATYPE :: avecttemp(numconfig,numr,10)
  DATATYPE :: mexpect(mcscfnum),ugexpect(mcscfnum),m2expect(mcscfnum)
  integer :: i,imc
  real*8 :: thistime
  integer, save :: calledflag=0

  inspfs=RESHAPE(yyy%cmfpsivec(spfstart:spfend,0),(/spfsize,nspf/))

  if (calledflag==0) then
     OFLWR "Writing psi stats."; CFL

     if (myrank.eq.1) then
        open(662, file=psistatsfile, status="unknown")
#ifdef REALGO        
        write(662,'(A6,20A14)') " Time ", " M ", " M2 ", " UG "
#else
        write(662,'(A6,20A28)') " Time ", " M ", " M2 ", " UG "
#endif
        close(662)
     endif
  endif

  calledflag=1

!! M & U/G

  avecttemp(:,:,:)=0d0

  
  if (spfrestrictflag.eq.1) then
     do imc=1,mcscfnum
        inavector(:,:)=RESHAPE(yyy%cmfpsivec(astart(imc):aend(imc),0),(/numconfig,numr/))

        do i=1,numr
           avecttemp(:,i,1)=inavector(:,i) * configmvals(:)
           avecttemp(:,i,2)=inavector(:,i) * configugvals(:)
           avecttemp(:,i,3)=inavector(:,i) * configmvals(:)**2
        enddo
        mexpect(imc)=dot(inavector,avecttemp(:,:,1),totadim)
        ugexpect(imc)=dot(inavector,avecttemp(:,:,2),totadim)
        m2expect(imc)=dot(inavector,avecttemp(:,:,3),totadim)
        sum=dot(inavector,inavector,totadim)

        mexpect(imc)=mexpect(imc)/sum;   m2expect(imc)=m2expect(imc)/sum;   ugexpect(imc)=ugexpect(imc)/sum
     enddo
  else
     mexpect=-99d0; m2expect=-99d0; ugexpect=-99d0
  endif


  if (myrank.eq.1) then
     open(662, file=psistatsfile, status="old", position="append")
     do i=1,mcscfnum
        write(662,'(F13.5,I3,2000F14.8)') thistime, i, mexpect(i),m2expect(i),ugexpect(i)
     enddo
     close(662)
  endif
end subroutine psistats





