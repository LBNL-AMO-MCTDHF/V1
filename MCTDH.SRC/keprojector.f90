
!! ALL ONE MODULE

!! ACTION 24

#include "Definitions.INC"

module keprojsubmod
contains

subroutine keprojector(inavector,inspfs,infac,www)
  use mpimod
  use parameters
  use walkmod
  use arbitrarymultmod
  use orbgathersubmod
  implicit none
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: inavector(www%totadim,mcscfnum), inspfs(spfsize,nspf)
  DATATYPE :: rad, kedot(nspf),krval,csum,kedot2(nspf)
  DATATYPE, allocatable,save :: kevects(:,:),keproj(:,:,:),keproj2(:,:,:)
  real*8, allocatable, save ::  kesum(:),energy(:), kesum2(:)
  DATATYPE,allocatable :: tempelecweights(:),tempvector(:)
  real*8 :: infac
  DATAECS :: ones(numr)
  integer, save :: allocd=0
  integer :: ispf,ii,jj,jspf,lowspf,highspf,numspf,myiostat

  if(numr.gt.1) then
     OFLWR "KEPROJ NOT SUPPORTED NUMR>1"; CFLST
  endif

  allocate( tempelecweights(spfsize),tempvector(www%totadim) )

  if (allocd.eq.0) then
     allocd=1
     allocate(keproj(nspf,nspf,nkeproj),kevects(spfsize,nkeproj),&
          kesum2(nkeproj),kesum(nkeproj),energy(nkeproj), keproj2(nspf,nspf,nkeproj)); 
     kevects(:,:)=0d0; keproj(:,:,:)=0d0; kesum=0d0; kesum2=0d0; keproj2(:,:,:)=0d0

!! elec weights has spfdims for m value support
     tempelecweights=RESHAPE(elecweights(:,:,:,1)*elecweights(:,:,:,2)*elecweights(:,:,:,3),&
          (/spfsize/))
     
     do ii=1,nkeproj
        energy(ii)=keprojminenergy+ ii*keprojenergystep
        do jj=1,spfsize
           rad=elecradii(jj)
           if (real(rad).gt.keprojminrad.and.real(rad).lt.keprojmaxrad) then
              if (imag(rad+(0d0,0d0)).gt.1d-6) then
                 OFLWR "AAUGH keprojector is in scaled region i think", rad; CFLST
              endif
              krval=sqrt(2*energy(ii))*rad
              kevects(jj,ii)=sin((rad-keprojminrad)*pi/(keprojmaxrad-keprojminrad))**2 &
                   * exp((0d0,1d0)*krval)*sqrt(tempelecweights(jj))
           endif
        enddo
        kevects(:,ii)=kevects(:,ii)/sqrt(dot(kevects(:,ii),kevects(:,ii),spfsize))
!        OFLWR "VECT ",ii,dot(kevects(:,ii),kevects(:,ii),spfsize); CFL
     enddo
     OFLWR "KEVECTS GOTTEN."; CFL
  endif

  lowspf=1; highspf=nspf
  if (parorbsplit.eq.1) then
     call getOrbSetRange(lowspf,highspf)
  endif
  numspf=highspf-lowspf+1

  do ii=1,nkeproj
     do ispf=lowspf,highspf
        kedot(ispf)=dot(inspfs(:,ispf),kevects(:,ii),spfsize)
        kedot2(ispf)=dot(inspfs(:,ispf),ALLCON(kevects(:,ii)),spfsize)
!        OFLWR "KEDOT", ispf,kedot(ispf); CFL
     enddo
     if (parorbsplit.eq.1) then
        call mpiorbgather(kedot,1)
        call mpiorbgather(kedot2,1)
     endif
     if (parorbsplit.eq.3) then
        call mympireduce(kedot,nspf)
        call mympireduce(kedot2,nspf)
     endif
     do ispf=1,nspf
        do jspf=1,nspf
           keproj(ispf,jspf,ii)=kedot(ispf)*CONJUGATE(kedot(jspf)) 
           keproj2(ispf,jspf,ii)=kedot2(ispf)*CONJUGATE(kedot2(jspf)) 
        enddo
     enddo
  enddo

  OFLWR "KEPROJ GOTTEN"; CFL

  do ii=1,nkeproj
     ones(:)=1d0

     call arbitraryconfig_mult_singles(www,keproj(:,:,ii),ones,inavector,tempvector,numr)
     csum=dot(inavector(:,1),tempvector(:),www%totadim)
     if (www%parconsplit.ne.0) then
        call mympireduceone(csum)
     endif
!! should be positive real  (um but not for cmctdhf, so don't use it then)
     if (real(csum).lt.1d6*abs(imag((0d0,0d0)+csum))) then     
        OFLWR "SUMERRKEPRJOFS", csum; CFL
     endif
     kesum(ii)=kesum(ii)+real(csum)*infac

     call arbitraryconfig_mult_singles(www,keproj2(:,:,ii),ones,inavector,tempvector,numr)
     csum=dot(inavector(:,1),tempvector(:),www%totadim)
     if (www%parconsplit.ne.0) then
        call mympireduceone(csum)
     endif
!! should be positive real  (um but not for cmctdhf, so don't use it then)
     if (real(csum).lt.1d6*abs(imag((0d0,0d0)+csum))) then     
        OFLWR "SUMERRKEPRJOFS", csum; CFL
     endif
     kesum2(ii)=kesum2(ii)+real(csum)*infac

  enddo
  OFLWR "KEPROJ DONE.";CFL

  if (myrank.eq.1) then
     open(1986,file="Dat/Keprojector.Dat", status="unknown",iostat=myiostat)
     call checkiostat(myiostat," opening keprojector.dat")
     do ii=1,nkeproj
        write(1986,'(100F30.20)',iostat=myiostat) energy(ii), kesum(ii),kesum2(ii)
     enddo
     call checkiostat(myiostat," writing keprojector.dat")
     close(1986)
  endif

  deallocate(tempelecweights,tempvector)

end subroutine keprojector
  
end module keprojsubmod


