
#include "Definitions.INC"


module drivingbiomod
  use biorthotypemod
  implicit none
  type(biorthotype),target :: drivingbiovar
  DATATYPE :: drivingoverlap(200)

end module

subroutine getdrivingoverlap(outdrivingoverlap)
  use drivingbiomod
  use parameters
  implicit none
  DATATYPE :: outdrivingoverlap(mcscfnum)
  outdrivingoverlap(:)= drivingoverlap(1:mcscfnum)
end subroutine getdrivingoverlap

subroutine drivingtrans(thistime)
  use parameters
  use drivingbiomod
  use configmod
  use xxxmod
  use walkmod
  use biorthomod
  use opmod
  implicit none

  DATATYPE :: dipmatxx(nspf,nspf),dipmatyy(nspf,nspf),dipmatzz(nspf,nspf),scalarpotyy(nspf,nspf),scalarpotzz(nspf,nspf)
  DATATYPE :: tempinvden(nspf,nspf),tempdenmat2(nspf,nspf), dreducedpottally(nspf,nspf,nspf,nspf)
  DATATYPE ::  tempdenmat(nspf,nspf), tempdrivingorbs(spfsize,nspf), &
       currentorbs(spfsize,nspf),tempdrivingavector(numr,firstconfig:lastconfig,mcscfnum), dot
  DATATYPE :: scalarpotxx(nspf,nspf), multorbsxx(spfsize,nspf), multorbsyy(spfsize,nspf), multorbszz(spfsize,nspf)
  DATATYPE :: dtwoereduced(reducedpotsize,nspf,nspf), pots(3)=0d0

  Type(CONFIGPTR) :: drivingptr     !!, nullptr
  integer ::  i, imc,j, nulltimes(10)
  real*8 :: thistime,rsum
  DATATYPE, target :: smo(nspf,nspf)
  DATAECS :: rvector(numr)

  yyy%drivingavectorsxx(:,:,:,0)=0;   yyy%drivingorbsxx(:,:,0) = 0.d0
  yyy%drivingavectorsyy(:,:,:,0)=0;   yyy%drivingorbsyy(:,:,0) = 0.d0
  yyy%drivingavectorszz(:,:,:,0)=0;   yyy%drivingorbszz(:,:,0) = 0.d0

  if (drivingflag==0) then
     OFLWR "Drivingflag is zero, bad"; CFLST
  endif

  if (drivingmethod.ne.0) then
     OFLWR "NOT SUPPORTED 558"; CFLST
  endif

  currentorbs(:,:)=RESHAPE(yyy%cmfpsivec(spfstart:spfend,0),(/spfsize,nspf/))
  
  do imc=1,mcscfnum

     tempdrivingavector(:,:,imc)=avector_driving(:,:,imc) * exp(timefac*drivingenergies(imc)*thistime)

     call bioset(drivingbiovar,smo,numr)
     call biortho(orbs_driving(:,:),currentorbs(:,:),tempdrivingorbs(:,:),tempdrivingavector(:,:,imc),drivingbiovar)

!! tempdrivingavector is vector in non-orthonormal tempdrivingorbs (which are bio to currentorbs)
!!    also is vector in currentorbs corresponding to insertion of identity
!!    so overlap is avector dot tempdrivingavector

     drivingoverlap(imc)=dot(yyy%cmfpsivec(astart(imc),0),tempdrivingavector(:,:,imc),totadim)

  enddo

  if (parconsplit.ne.0) then
     call mympireduce(drivingoverlap,mcscfnum)
  endif

  if (tdflag.ne.1) then
     return
  endif
  call vectdpot(thistime,velflag,pots)
  rsum=0d0
  do i=1,3
     rsum=rsum+abs(pots(i))**2
  enddo
  if (rsum.eq.0d0) then
     return
  endif


  call configptralloc(drivingptr)
  call zero_cptr(drivingptr)
  call all_matel0(drivingptr,currentorbs,tempdrivingorbs,dtwoereduced,nulltimes)

  if (velflag.eq.0) then
     rvector(:)=bondpoints(:)
  else
     rvector(:)=1/bondpoints(:)
  endif

  do imc=1,mcscfnum
     do i=1,nspf
        if (velflag.eq.0) then
           call lenmultiply(tempdrivingorbs(:,i),multorbsxx(:,i),DATAONE,DATAZERO,DATAZERO)
           call lenmultiply(tempdrivingorbs(:,i),multorbsyy(:,i),DATAZERO,DATAONE,DATAZERO)
           call lenmultiply(tempdrivingorbs(:,i),multorbszz(:,i),DATAZERO,DATAZERO,DATAONE)
        else
           call velmultiply(tempdrivingorbs(:,i),multorbsxx(:,i),DATAONE,DATAZERO,DATAZERO)
           call velmultiply(tempdrivingorbs(:,i),multorbsyy(:,i),DATAZERO,DATAONE,DATAZERO)
           call velmultiply(tempdrivingorbs(:,i),multorbszz(:,i),DATAZERO,DATAZERO,DATAONE)
        endif
     enddo
     do i=1,nspf
        do j=1,nspf
           dipmatxx(j,i)=dot(currentorbs(:,j),multorbsxx(:,i),spfsize)
           dipmatyy(j,i)=dot(currentorbs(:,j),multorbsyy(:,i),spfsize)
           dipmatzz(j,i)=dot(currentorbs(:,j),multorbszz(:,i),spfsize)
        enddo
     enddo

     if (parorbsplit.eq.3) then
        call mympireduce(dipmatxx,nspf**2)
        call mympireduce(dipmatyy,nspf**2)
        call mympireduce(dipmatzz,nspf**2)
     endif


     call arbitraryconfig_mult(dipmatxx,rvector,tempdrivingavector(:,:,imc),     yyy%drivingavectorsxx(:,:,imc,0),numr)
     call arbitraryconfig_mult(dipmatyy,rvector,tempdrivingavector(:,:,imc),     yyy%drivingavectorsyy(:,:,imc,0),numr)
     call arbitraryconfig_mult(dipmatzz,rvector,tempdrivingavector(:,:,imc),     yyy%drivingavectorszz(:,:,imc,0),numr)
  enddo

  call get_tworeducedx( dreducedpottally, yyy%cmfpsivec(astart(1),0), tempdrivingavector)

  call MYGEMM('N','N', 1,nspf**2,nspf**2,(1.0d0,0.d0), dipmatxx,1,dreducedpottally,nspf**2, (0.d0,0.d0),scalarpotxx,1)
  call MYGEMM('N','N', 1,nspf**2,nspf**2,(1.0d0,0.d0), dipmatyy,1,dreducedpottally,nspf**2, (0.d0,0.d0),scalarpotyy,1)
  call MYGEMM('N','N', 1,nspf**2,nspf**2,(1.0d0,0.d0), dipmatzz,1,dreducedpottally,nspf**2, (0.d0,0.d0),scalarpotzz,1)


!! NOT DENMAT, REDUCED OPERATOR FOR PULSE (rvector)

  call getdenmat00(tempdrivingavector(:,:,:),yyy%cmfpsivec(astart(1),0),rvector(:),tempdenmat(:,:),numr,mcscfnum)

!!$  tempdenmat(:,:)=0.d0
!!$
!!$  do imc=1,mcscfnum
!!$     bigavector(:,:)=0d0
!!$     bigavector(:,firstconfig:lastconfig)=RESHAPE(yyy%cmfpsivec(astart(imc):aend(imc),0),(/numr,localnconfig/))
!!$     if (parconsplit.ne.0) then
!!$        call mpiallgather(bigavector,numconfig*numr,configsperproc(:)*numr,maxconfigsperproc*numr)
!!$     endif
!!$
!!$     do config1=botconfig,topconfig
!!$
!!$        a1(:)=tempdrivingavector(:,config1,imc)   *rvector(:)
!!$
!!$        do iwalk=1,numsinglewalks(config1)
!!$           config2=singlewalk(iwalk,config1)
!!$           dirphase=singlewalkdirphase(iwalk,config1)
!!$           a2(:)=bigavector(:,config2)
!!$           ispf=singlewalkopspf(1,iwalk,config1)  !! goes with config1, driving
!!$           jspf=singlewalkopspf(2,iwalk,config1)  !! goes with config2, bra, current avector
!!$           
!!$           tempdenmat(ispf,jspf)=tempdenmat(ispf,jspf)+ &  
!!$                dirphase*dot(a2,a1,numr)
!!$        enddo
!!$     enddo
!!$
!!$  enddo
!!$
!!$  call mympireduce(tempdenmat,nspf**2)

!! (NO TIMEFAC)

  tempinvden(:,:)=yyy%invdenmat(:,:,0) 

  call MYGEMM('N','N', nspf,nspf,nspf,DATAONE, tempdenmat, nspf, tempinvden, nspf, DATAZERO, tempdenmat2, nspf)

  call MYGEMM('N','N', spfsize,nspf,nspf,DATAONE, multorbsxx,       spfsize, tempdenmat2, nspf, DATAZERO, yyy%drivingorbsxx(:,:,0),spfsize)
  call MYGEMM('N','N', spfsize,nspf,nspf,DATAONE, multorbsyy,       spfsize, tempdenmat2, nspf, DATAZERO, yyy%drivingorbsyy(:,:,0),spfsize)
  call MYGEMM('N','N', spfsize,nspf,nspf,DATAONE, multorbszz,       spfsize, tempdenmat2, nspf, DATAZERO, yyy%drivingorbszz(:,:,0),spfsize)


!! LIKE REDUCEDPOT, SCALARPOT SHOULD LEFT MULTIPLY INVECTORS-TRANSPOSE TO GET OUTVECTORS-TRANSPOSE

  call MYGEMM('T','N', nspf,nspf,nspf,DATAONE, scalarpotxx,  nspf, tempinvden, nspf, DATAZERO, tempdenmat2, nspf)
  call MYGEMM('N','N', spfsize,nspf,nspf,DATAONE, tempdrivingorbs,spfsize, tempdenmat2, nspf, DATAONE, yyy%drivingorbsxx(:,:,0),spfsize)

  call MYGEMM('T','N', nspf,nspf,nspf,DATAONE, scalarpotyy,  nspf, tempinvden, nspf, DATAZERO, tempdenmat2, nspf)
  call MYGEMM('N','N', spfsize,nspf,nspf,DATAONE, tempdrivingorbs,spfsize, tempdenmat2, nspf, DATAONE, yyy%drivingorbsyy(:,:,0),spfsize)

  call MYGEMM('T','N', nspf,nspf,nspf,DATAONE, scalarpotzz,  nspf, tempinvden, nspf, DATAZERO, tempdenmat2, nspf)
  call MYGEMM('N','N', spfsize,nspf,nspf,DATAONE, tempdrivingorbs,spfsize, tempdenmat2, nspf, DATAONE, yyy%drivingorbszz(:,:,0),spfsize)

  
end subroutine drivingtrans







subroutine drivinginit(inenergies)
  use parameters
  use xxxmod
  use opmod
  implicit none

  DATATYPE,intent(in) :: inenergies(mcscfnum)

  orbs_driving(:,:)=RESHAPE(yyy%cmfpsivec(spfstart:spfend,0),(/spfsize,nspf/))

  avector_driving(:,:,:)=RESHAPE(yyy%cmfpsivec(astart(1):aend(mcscfnum),0),(/numr,localnconfig,mcscfnum/))*drivingproportion

  yyy%cmfpsivec(astart(1):aend(mcscfnum),0) = yyy%cmfpsivec(astart(1):aend(mcscfnum),0) * (1d0-drivingproportion)

  drivingenergies(1:mcscfnum)=inenergies(:)

end subroutine drivinginit








