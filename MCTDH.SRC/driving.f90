

!! subroutines for PSI-PRIME method (still unstable) using drivingflag

#include "Definitions.INC"

module drivingbiomod
  use biorthotypemod
  implicit none
  type(biorthotype),target :: drivingbiovar
  DATATYPE :: drivingoverlap(200)=0,sxx(200)=0,syy(200)=0,szz(200)=0,&
       axx(200)=0,ayy(200)=0,azz(200)=0
end module

subroutine getdrivingoverlap(outdrivingoverlap,numvects)
  use drivingbiomod
  implicit none
  integer,intent(in) :: numvects
  DATATYPE,intent(out) :: outdrivingoverlap(numvects)
  outdrivingoverlap(:)= drivingoverlap(1:numvects)
end subroutine getdrivingoverlap


subroutine dipolesub_driving(oaxx,oayy,oazz,osxx,osyy,oszz,numvects)
  use drivingbiomod
  implicit none
  integer,intent(in) :: numvects
  DATATYPE,intent(out) :: oaxx(numvects),oayy(numvects),oazz(numvects),&
       osxx(numvects),osyy(numvects),oszz(numvects)
  oaxx(:)=axx(1:numvects)
  oayy(:)=ayy(1:numvects)
  oazz(:)=azz(1:numvects)
  osxx(:)=sxx(1:numvects)
  osyy(:)=syy(1:numvects)
  oszz(:)=szz(1:numvects)
end subroutine dipolesub_driving


subroutine drivingtrans(thistime)
  use parameters
  use drivingbiomod
  use configmod
  use xxxmod
  use biorthomod
  use opmod
  use arbitrarymultmod
  use mpi_orbsetmod
  use densubmod
  use mpisubmod
  use orbmultsubmod
  use allmat0mod
  use pulsesubmod
  use tworeducedxmod
  implicit none
  DATATYPE :: dipmatxx(nspf,nspf),dipmatyy(nspf,nspf),dipmatzz(nspf,nspf),&
       scalarpotyy(nspf,nspf),scalarpotzz(nspf,nspf),tempinvden(nspf,nspf),&
       tempdenmat2(nspf,nspf), scalarpotxx(nspf,nspf),  tempdenmat(nspf,nspf)  !! AUTOMATIC
  DATATYPE,allocatable :: currentorbs(:,:),tempdrivingorbs(:,:),&
       tempdrivingavector(:,:,:),dtwoereduced(:,:,:),&
       dreducedpottally(:,:,:,:),multorbsxx(:,:),multorbsyy(:,:),multorbszz(:,:)
  DATATYPE :: pots(3)=0d0
  Type(CONFIGPTR) :: drivingptr
  integer ::  i, imc,j, nulltimes(10),firstorb,lastorb
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

  allocate(currentorbs(spfsize,nspf), tempdrivingorbs(spfsize,nspf),&
       tempdrivingavector(numr,first_config:last_config,mcscfnum))

  currentorbs(:,:)=RESHAPE(yyy%cmfspfs(:,0),(/spfsize,nspf/))
  
  do imc=1,mcscfnum

     if (tot_adim.gt.0) then
        tempdrivingavector(:,:,imc)=avector_driving(:,:,imc) &
             * exp(timefac*drivingenergies(imc)*thistime)
     endif

     call bioset(drivingbiovar,smo,numr,bioww)
     call biortho(orbs_driving(:,:),currentorbs(:,:),&
          tempdrivingorbs(:,:),tempdrivingavector(:,:,imc),drivingbiovar)

!! tempdrivingavector is vector in non-orthonormal tempdrivingorbs (which are bio to currentorbs)
!!    also is vector in currentorbs corresponding to insertion of identity
!!    so overlap is avector dot tempdrivingavector
     drivingoverlap(imc)=0
     if (tot_adim.gt.0) then
        drivingoverlap(imc)=dot(yyy%cmfavec(:,imc,0),tempdrivingavector(:,:,imc),tot_adim)
     endif
  enddo

  if (par_consplit.ne.0) then
     call mympireduce(drivingoverlap,mcscfnum)
  endif

  if (tdflag.eq.0) then
     deallocate(currentorbs,tempdrivingorbs,tempdrivingavector)
     return
  endif

  rsum=0d0
  call vectdpot(thistime,velflag,pots,-1)
  do i=1,3
     rsum=rsum+abs(pots(i))**2
  enddo
  if (rsum.eq.0d0) then
     deallocate(currentorbs,tempdrivingorbs,tempdrivingavector)
     return
  endif

  if (parorbsplit.eq.1) then
     firstorb=firstmpiorb
     lastorb=firstmpiorb+orbsperproc-1
  else
     firstorb=1
     lastorb=nspf
  endif

  allocate(dtwoereduced(reducedpotsize,nspf,firstorb:lastorb),&
       dreducedpottally(nspf,nspf,nspf,nspf),&
       multorbsxx(spfsize,nspf), multorbsyy(spfsize,nspf), multorbszz(spfsize,nspf))
  dtwoereduced=0; dreducedpottally=0; multorbsxx=0; multorbsyy=0; multorbszz=0

  call configptralloc(drivingptr,www)
  call zero_cptr(drivingptr)
  if (holeflag.ne.0) then
     OFLWR "Not done holeflag driving"; CFLST
  endif
  call all_matel0(0,drivingptr,currentorbs,tempdrivingorbs,dtwoereduced,nulltimes,firstorb,lastorb)

  if (velflag.eq.0) then
     rvector(:)=bondpoints(:)
  else
     rvector(:)=1/bondpoints(:)
  endif

  do imc=1,mcscfnum
     if (velflag.eq.0) then
        call lenmultiply(nspf,tempdrivingorbs(:,:),multorbsxx(:,:),DATAONE,DATAZERO,DATAZERO)
        call lenmultiply(nspf,tempdrivingorbs(:,:),multorbsyy(:,:),DATAZERO,DATAONE,DATAZERO)
        call lenmultiply(nspf,tempdrivingorbs(:,:),multorbszz(:,:),DATAZERO,DATAZERO,DATAONE)
     else
        call velmultiply(nspf,tempdrivingorbs(:,:),multorbsxx(:,:),DATAONE,DATAZERO,DATAZERO)
        call velmultiply(nspf,tempdrivingorbs(:,:),multorbsyy(:,:),DATAZERO,DATAONE,DATAZERO)
        call velmultiply(nspf,tempdrivingorbs(:,:),multorbszz(:,:),DATAZERO,DATAZERO,DATAONE)
     endif

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


     call arbitraryconfig_mult_singles(www,dipmatxx,rvector,&
          tempdrivingavector(:,:,imc),     yyy%drivingavectorsxx(:,:,imc,0),numr)
     call arbitraryconfig_mult_singles(www,dipmatyy,rvector,&
          tempdrivingavector(:,:,imc),     yyy%drivingavectorsyy(:,:,imc,0),numr)
     call arbitraryconfig_mult_singles(www,dipmatzz,rvector,&
          tempdrivingavector(:,:,imc),     yyy%drivingavectorszz(:,:,imc,0),numr)

     axx(imc)=0; ayy(imc)=0; azz(imc)=0
     if (tot_adim.gt.0) then
        axx(imc)=dot(yyy%cmfavec(:,imc,0),yyy%drivingavectorsxx(:,:,imc,0),tot_adim)
        ayy(imc)=dot(yyy%cmfavec(:,imc,0),yyy%drivingavectorsyy(:,:,imc,0),tot_adim)
        azz(imc)=dot(yyy%cmfavec(:,imc,0),yyy%drivingavectorszz(:,:,imc,0),tot_adim)
     endif
  enddo

  if (par_consplit.ne.0) then
     call mympireduce(axx,mcscfnum)
     call mympireduce(ayy,mcscfnum)
     call mympireduce(azz,mcscfnum)
  endif

  call get_tworeducedx(www, dreducedpottally, yyy%cmfavec(:,:,0),  tempdrivingavector, mcscfnum)

  call MYGEMM('N','N', 1,nspf**2,nspf**2,(1.0d0,0.d0), dipmatxx,1,&
       dreducedpottally,nspf**2, (0.d0,0.d0),scalarpotxx,1)
  call MYGEMM('N','N', 1,nspf**2,nspf**2,(1.0d0,0.d0), dipmatyy,1,&
       dreducedpottally,nspf**2, (0.d0,0.d0),scalarpotyy,1)
  call MYGEMM('N','N', 1,nspf**2,nspf**2,(1.0d0,0.d0), dipmatzz,1,&
       dreducedpottally,nspf**2, (0.d0,0.d0),scalarpotzz,1)


!! NOT DENMAT, REDUCED OPERATOR FOR PULSE (rvector)

  call getdenmat00(1,www,tempdrivingavector(:,:,:),yyy%cmfavec(:,:,0),&
       rvector(:),tempdenmat(:,:),numr,mcscfnum)

!!$  tempdenmat(:,:)=0.d0
!!$
!!$  do imc=1,mcscfnum
!!$     bigavector(:,:)=0d0
!!$     bigavector(:,first_config:last_config)=RESHAPE(yyy%cmfavec(:,imc,0),(/numr,local_nconfig/))
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

  call MYGEMM('N','N', nspf,nspf,nspf,DATAONE, tempdenmat, nspf, &
       tempinvden, nspf, DATAZERO, tempdenmat2, nspf)

  call MYGEMM('N','N', spfsize,nspf,nspf,DATAONE, multorbsxx,     &
       spfsize, tempdenmat2, nspf, DATAZERO, yyy%drivingorbsxx(:,:,0),spfsize)
  call MYGEMM('N','N', spfsize,nspf,nspf,DATAONE, multorbsyy,     &
       spfsize, tempdenmat2, nspf, DATAZERO, yyy%drivingorbsyy(:,:,0),spfsize)
  call MYGEMM('N','N', spfsize,nspf,nspf,DATAONE, multorbszz,    &
       spfsize, tempdenmat2, nspf, DATAZERO, yyy%drivingorbszz(:,:,0),spfsize)


!! LIKE REDUCEDPOT, SCALARPOT SHOULD LEFT MULTIPLY INVECTORS-TRANSPOSE TO GET OUTVECTORS-TRANSPOSE

  call MYGEMM('T','N', nspf,nspf,nspf,DATAONE, scalarpotxx,  nspf, &
       tempinvden, nspf, DATAZERO, tempdenmat2, nspf)
  call MYGEMM('N','N', spfsize,nspf,nspf,DATAONE, tempdrivingorbs,&
       spfsize, tempdenmat2, nspf, DATAONE, yyy%drivingorbsxx(:,:,0),spfsize)

  call MYGEMM('T','N', nspf,nspf,nspf,DATAONE, scalarpotyy,  nspf, &
       tempinvden, nspf, DATAZERO, tempdenmat2, nspf)
  call MYGEMM('N','N', spfsize,nspf,nspf,DATAONE, tempdrivingorbs, &
       spfsize, tempdenmat2, nspf, DATAONE, yyy%drivingorbsyy(:,:,0),spfsize)

  call MYGEMM('T','N', nspf,nspf,nspf,DATAONE, scalarpotzz,  nspf, &
       tempinvden, nspf, DATAZERO, tempdenmat2, nspf)
  call MYGEMM('N','N', spfsize,nspf,nspf,DATAONE, tempdrivingorbs, &
       spfsize, tempdenmat2, nspf, DATAONE, yyy%drivingorbszz(:,:,0),spfsize)

  call configptrdealloc(drivingptr)
  deallocate(currentorbs,tempdrivingorbs,tempdrivingavector)
  deallocate(dtwoereduced, dreducedpottally,  multorbsxx,multorbsyy,multorbszz)
  
end subroutine drivingtrans



subroutine drivinginit(inenergies)
  use parameters
  use drivingbiomod
  use xxxmod
  use opmod
  use configmod
  use dipsubonemod
  implicit none
  DATAECS,intent(in) :: inenergies(mcscfnum)
  DATATYPE :: myexpects(3), blah
  integer :: imc

  orbs_driving(:,:)=RESHAPE(yyy%cmfspfs(:,0),(/spfsize,nspf/))

  do imc=1,mcscfnum
     call dipolesub_one(www,bioww,yyy%cmfavec(:,imc,0),yyy%cmfavec(:,imc,0),&
          yyy%cmfspfs(:,0),yyy%cmfspfs(:,0),.false.,myexpects(:),blah)
     sxx(imc)=myexpects(1)
     syy(imc)=myexpects(2)
     szz(imc)=myexpects(3)
  enddo

  if (tot_adim.gt.0) then
     avector_driving(:,:,:)=&
          RESHAPE(yyy%cmfavec(:,:,0),(/numr,local_nconfig,mcscfnum/))*drivingproportion
     yyy%cmfavec(:,:,0) = yyy%cmfavec(:,:,0) * (1d0-drivingproportion)
  endif

  drivingenergies(1:mcscfnum)=inenergies(:)


end subroutine drivinginit








