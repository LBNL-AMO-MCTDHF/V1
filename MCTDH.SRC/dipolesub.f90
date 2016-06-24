

!! ACTION 21 (emission/absorption) subroutine

#include "Definitions.INC"

module dipolemod
  implicit none
  DATATYPE, allocatable :: dipoleexpects(:,:,:),    dipolenormsq(:)
end module dipolemod


subroutine dipolesub_initial()
  use dipolemod
  use parameters
  use pulse_parameters !! conjgpropflag
  implicit none

  if (conjgpropflag.eq.0) then
     allocate(dipoleexpects(0:autosize,3,1))   !! x,y,z
  else
     allocate(dipoleexpects(0:autosize,3,4))
  endif
  allocate(dipolenormsq(0:autosize))

end subroutine


subroutine dipolesub_final()
  use dipolemod
  implicit none
  deallocate( dipoleexpects, dipolenormsq)

end subroutine dipolesub_final


subroutine checkorbsetrange(checknspf,flag)
  use parameters
  implicit none
  integer,intent(in) :: checknspf
  integer,intent(out) :: flag
  flag=0
  if (nspf.ne.checknspf) then
     flag=1
  endif
end subroutine checkorbsetrange


module dipbiomod
  use biorthotypemod
  implicit none
  type(biorthotype),target :: dipbiovar
end module dipbiomod


module dipsubonemod
contains

subroutine dipolesub_one(wwin,bbin,in_abra,&    !! ok unused bbin
     in_aket,inspfs,dipole_expects)
  use r_parameters
  use spfsize_parameters
  use walkmod
  use fileptrmod
  use dotmod
  use dipbiomod
  use biorthomod
  use arbitrarymultmod
  use orbgathersubmod
  use mpisubmod
  implicit none
  type(walktype),intent(in) :: wwin
  type(walktype),target,intent(in) :: bbin
  DATATYPE, intent(in) :: inspfs(  spfsize, wwin%nspf ), &
       in_abra(numr,wwin%firstconfig:wwin%lastconfig),&
       in_aket(numr,wwin%firstconfig:wwin%lastconfig)
  DATATYPE,intent(out) :: dipole_expects(3)
  DATATYPE,allocatable :: tempvector(:,:),tempspfs(:,:),abra(:,:),workspfs(:,:),&
       aket(:,:)
  DATATYPE :: nullcomplex(1),dipoles(3), dipolemat(wwin%nspf,wwin%nspf),csum
  DATAECS :: rvector(numr)
!!$  DATATYPE :: norm   !! datatype in case abra.ne.aket
  integer :: i,lowspf,highspf,numspf
#ifdef CNORMFLAG
  DATATYPE,target :: smo(wwin%nspf,wwin%nspf)
#endif

  lowspf=1; highspf=wwin%nspf
  if (parorbsplit.eq.1) then
     call checkorbsetrange(wwin%nspf,i)
     if (i.ne.0) then
        OFLWR "error exit, can't do dipolesub parorbsplit.eq.1 with",wwin%nspf,"orbitals"; CFLST
     endif
     call getOrbSetRange(lowspf,highspf)
  endif
  numspf=highspf-lowspf+1
 
  allocate(tempvector(numr,wwin%firstconfig:wwin%lastconfig+1), tempspfs(spfsize,lowspf:highspf+1),&
       abra(numr,wwin%firstconfig:wwin%lastconfig+1),workspfs(spfsize,wwin%nspf),&
       aket(numr,wwin%firstconfig:wwin%lastconfig+1))

  tempvector=0; tempspfs=0; workspfs=0; abra=0; aket=0

  if (wwin%lastconfig.ge.wwin%firstconfig) then
     abra(:,wwin%firstconfig:wwin%lastconfig)=in_abra(:,:)
     aket(:,wwin%firstconfig:wwin%lastconfig)=in_aket(:,:)
  endif

#ifndef CNORMFLAG
  workspfs(:,:)=inspfs(:,:)
#else
  call bioset(dipbiovar,smo,numr,bbin)
  dipbiovar%hermonly=.true.
  call biortho(inspfs,inspfs,workspfs,abra,dipbiovar)
  dipbiovar%hermonly=.false.
#endif


!!$  csum=dot(abra,aket,wwin%totadim)
!!$  if (wwin%parconsplit.ne.0) then
!!$     call mympireduceone(csum)
!!$  endif
!!$  norm=sqrt(csum)

!! independent of R for now.  multiply by R for prolate  (R set to 1 for atom)
  call nucdipvalue(nullcomplex,dipoles)

  do i=1,numr
     tempvector(i,:)=aket(i,:)*bondpoints(i)
  enddo
  csum=0d0
  if (wwin%totadim.gt.0) then
     csum=hermdot(abra,tempvector,wwin%totadim)
  endif
  if (wwin%parconsplit.ne.0) then
     call mympireduceone(csum)
  endif
  dipoles(:)=dipoles(:)*csum

!! Z DIPOLE

  dipolemat(:,:)=0d0
  if (numspf.gt.0) then
     call mult_zdipole(numspf,inspfs(:,lowspf:highspf),tempspfs(:,lowspf:highspf),1)
     call MYGEMM('C','N',wwin%nspf,numspf,spfsize,DATAONE, workspfs, spfsize, &
          tempspfs(:,lowspf:highspf), spfsize, DATAZERO, dipolemat(:,lowspf:highspf), wwin%nspf)
  endif
  if (parorbsplit.eq.1) then
     call mpiorbgather(dipolemat,wwin%nspf)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(dipolemat(:,:),wwin%nspf**2)
  endif

  rvector(:)=bondpoints(:)
  call arbitraryconfig_mult_singles(wwin,dipolemat,rvector,aket,tempvector,numr)
  dipole_expects(3)=0d0
  if (wwin%totadim.gt.0) then
     dipole_expects(3)=hermdot(abra,tempvector,wwin%totadim)
  endif
  if (wwin%parconsplit.ne.0) then
     call mympireduceone(dipole_expects(3))
  endif
  dipole_expects(3)=dipole_expects(3) + dipoles(3)

!! Y DIPOLE

  dipolemat(:,:)=0d0
  if (numspf.gt.0) then
     call mult_ydipole(numspf,inspfs(:,lowspf:highspf),tempspfs(:,lowspf:highspf),1)

     call MYGEMM('C','N',wwin%nspf,numspf,spfsize,DATAONE, workspfs, spfsize, &
          tempspfs(:,lowspf:highspf), spfsize, DATAZERO, dipolemat(:,lowspf:highspf), wwin%nspf)
  endif
  if (parorbsplit.eq.1) then
     call mpiorbgather(dipolemat,wwin%nspf)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(dipolemat(:,:),wwin%nspf**2)
  endif

  rvector(:)=bondpoints(:)
  call arbitraryconfig_mult_singles(wwin,dipolemat,rvector,aket,tempvector,numr)
  dipole_expects(2)=0d0
  if (wwin%totadim.gt.0) then
     dipole_expects(2)=hermdot(abra,tempvector,wwin%totadim)
  endif
  if (wwin%parconsplit.ne.0) then
     call mympireduceone(dipole_expects(2))
  endif
  dipole_expects(2)=dipole_expects(2) + dipoles(2)

!! X DIPOLE

  dipolemat(:,:)=0d0
  if (numspf.gt.0) then
     call mult_xdipole(numspf,inspfs(:,lowspf:highspf),tempspfs(:,lowspf:highspf),1)

     call MYGEMM('C','N',wwin%nspf,numspf,spfsize,DATAONE, workspfs, spfsize, &
          tempspfs(:,lowspf:highspf), spfsize, DATAZERO, dipolemat(:,lowspf:highspf), wwin%nspf)
  endif
  if (parorbsplit.eq.1) then
     call mpiorbgather(dipolemat,wwin%nspf)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(dipolemat(:,:),wwin%nspf**2)
  endif

  rvector(:)=bondpoints(:)
  call arbitraryconfig_mult_singles(wwin,dipolemat,rvector,aket,tempvector,numr)
  dipole_expects(1)=0d0
  if (wwin%totadim.gt.0) then
     dipole_expects(1)=hermdot(abra,tempvector,wwin%totadim)
  endif
  if (wwin%parconsplit.ne.0) then
     call mympireduceone(dipole_expects(1))
  endif
  dipole_expects(1)=dipole_expects(1) + dipoles(1)

  deallocate(tempvector,tempspfs,abra,aket,workspfs)

!!$#ifdef CNORMFLAG
!!$  dipole_expects(1)=dipole_expects(1)*abs(norm)/norm
!!$  dipole_expects(2)=dipole_expects(2)*abs(norm)/norm
!!$  dipole_expects(3)=dipole_expects(3)*abs(norm)/norm
!!$#endif

end subroutine dipolesub_one

end module dipsubonemod


module dipcallsubmod
contains

!! actually have numdata+1 data points in indipolearray

  subroutine dipolecall(numdata, indipolearrays, outenames, outtworknames, outangworknames, &
       outftnames, outoworknames, outophotonnames, sflag, referencepulses, npulses)
    use parameters
    use pulse_parameters
    use mpimod
    use pulsesubmod
    implicit none
    integer,intent(in) :: numdata, sflag, referencepulses, npulses
    DATATYPE,intent(in) :: indipolearrays(0:autosize,3)
    character(len=SLN),intent(in) :: outenames(3), outtworknames(3), outangworknames(3), &
         outftnames(9), outoworknames(9), outophotonnames(9)
    complex*16,allocatable :: fftrans(:,:),eft(:,:),each_eft(:,:,:),dipolearrays(:,:),&
         dipole_diff(:,:)
!! real valued variables -- field (taken real always); 
!!        integrals domega and circ polarization output left real valued
    real*8, allocatable :: &
         worksums(:,:,:), photsums(:,:,:), totworksums(:,:), totphotsums(:,:), &
         workpers(:,:,:), photpers(:,:,:), totworkpers(:,:), totphotpers(:,:), &
         sumrule(:,:), each_efield_ang(:,:,:), moment(:), dipole_ang(:,:), &
         angworksum0(:,:,:), totangworksum0(:,:),  efield(:,:), each_efield(:,:,:)
!! making integrals dt for work x y z complex valued for complex domcke wave mixing
    complex*16, allocatable :: worksum0(:,:,:), totworksum0(:,:)
    DATATYPE :: pots(3,npulses)
    real*8 :: estep, thistime, myenergy,xsecunits, windowfunct
    integer :: i,getlen,myiostat,ipulse,numft,ii
    character (len=7) :: number

#ifdef REALGO
    OFLWR "Cant use dipolesub for real valued code."; CFLST
#endif
    
    pots=0

    allocate(dipolearrays(0:numdata,3), efield(0:numdata,3), each_efield(0:numdata,3,npulses))
    dipolearrays=0; efield=0; each_efield=0

    if (sflag.ne.0) then
       thistime=numdata*par_timestep*autosteps
       write(number,'(I7)') 1000000+floor(thistime)
    endif

!! LENGTH GAUGE electric field vectdpot0

    do i=0,numdata
       if (referencepulses.eq.0) then
          do ipulse=1,npulses
             call vectdpot0(i*par_timestep*autosteps,0,pots(:,ipulse),-1,ipulse,ipulse)
          enddo
       else
          do ipulse=1,npulses
             call vectdpot0(i*par_timestep*autosteps,0,pots(:,ipulse),-1,ipulse+numpulses,ipulse+numpulses)
          enddo
       endif
       each_efield(i,:,:)=real(pots(:,:),8)  !! is real valued (imc = -1 above)
    enddo

    efield=0
    do ipulse=1,npulses
       efield(:,:)=efield(:,:)+each_efield(:,:,ipulse)
    enddo
    do i=0,numdata
       dipolearrays(i,:)=indipolearrays(i,:)-indipolearrays(0,:)
    enddo

    if (myrank.eq.1) then
       do ii=1,3
          open(171,file=outenames(ii),status="unknown",iostat=myiostat)
          call checkiostat(myiostat,"opening "//outenames(ii))
          write(171,*,iostat=myiostat) "#   ", numdata
          call checkiostat(myiostat,"writing "//outenames(ii))
          do i=0,numdata
             write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat)  i*par_timestep*autosteps, &
                  dipolearrays(i,ii),efield(i,ii),each_efield(i,ii,:)
          enddo
          call checkiostat(myiostat,"writing "//outenames(ii))
          close(171)
       enddo
       if (sflag.ne.0) then
          do ii=1,3
             open(171,file=outenames(ii)(1:getlen(outenames(ii)))//number(2:7),status="unknown",iostat=myiostat)
             call checkiostat(myiostat,"opening "//outenames(ii)(1:getlen(outenames(ii)))//number(2:7))
             write(171,*,iostat=myiostat) "#   ", numdata
             call checkiostat(myiostat,"writing "//outenames(ii)(1:getlen(outenames(ii)))//number(2:7))
             do i=0,numdata
                write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat)  i*par_timestep*autosteps, &
                     dipolearrays(i,ii),efield(i,ii),each_efield(i,ii,:)
             enddo
             call checkiostat(myiostat,"writing "//outenames(ii)(1:getlen(outenames(ii)))//number(2:7))
             close(171)
          enddo
       endif
    endif
    call mpibarrier()

!! work done by pulse integral dt

    allocate(dipole_diff(0:numdata,3), worksum0(0:numdata,3,npulses), totworksum0(0:numdata,3))
    dipole_diff=0d0;    worksum0=0;   totworksum0=0;  

    do i=1,3
       call complexdiff(numdata+1,dipolearrays(:,i),dipole_diff(:,i),.false.)
    enddo
    dipole_diff(:,:)= dipole_diff(:,:) / par_timestep / autosteps

    do ipulse=1,npulses
       worksum0(0,:,ipulse) = (-1) * dipole_diff(0,:) * each_efield(0,:,ipulse) * par_timestep * autosteps
       do i=1,numdata
          worksum0(i,:,ipulse)=worksum0(i-1,:,ipulse) - dipole_diff(i,:) * each_efield(i,:,ipulse) * par_timestep * autosteps
       enddo
    enddo
    totworksum0(:,:)=0d0
    do ipulse=1,npulses
       totworksum0(:,:)=totworksum0(:,:)+worksum0(:,:,ipulse)
    enddo

    if (myrank.eq.1) then
       do ii=1,3
          open(171,file=outtworknames(ii),status="unknown",iostat=myiostat)
          call checkiostat(myiostat,"opening "//outtworknames(ii))
          do i=0,numdata
             write(171,'(A25,F10.5,400F20.10)') " EACH PULSE WORK T= ", i*par_timestep*autosteps,&
                  totworksum0(i,ii),worksum0(i,ii,:)
          enddo
          close(171)
       enddo

       if (sflag.ne.0) then
          do ii=1,3
             open(171,file=outtworknames(ii)(1:getlen(outtworknames(ii)))//number(2:7),status="unknown",iostat=myiostat)
             call checkiostat(myiostat,"opening "//outtworknames(ii)(1:getlen(outtworknames(ii)))//number(2:7))
             do i=0,numdata
                write(171,'(A25,F10.5,400F20.10)') " EACH PULSE WORK T= ", i*par_timestep*autosteps,&
                     totworksum0(i,ii),worksum0(i,ii,:)
             enddo
             close(171)
          enddo
       endif
    endif
    call mpibarrier()
    deallocate(worksum0,totworksum0)

    if (act21circ.ne.0) then
       allocate(angworksum0(0:numdata,3,npulses), totangworksum0(0:numdata,3), dipole_ang(0:numdata,3), &
            each_efield_ang(0:numdata,3,npulses), moment(0:numdata))
       angworksum0=0;  totangworksum0=0; dipole_ang=0;   each_efield_ang=0;    moment=0

!! not worrying about complex values

       dipole_ang(:,1) = real(dipolearrays(:,2) * dipole_diff(:,3) - dipolearrays(:,3) * dipole_diff(:,2),8)
       dipole_ang(:,2) = real(dipolearrays(:,3) * dipole_diff(:,1) - dipolearrays(:,1) * dipole_diff(:,3),8)
       dipole_ang(:,3) = real(dipolearrays(:,1) * dipole_diff(:,2) - dipolearrays(:,2) * dipole_diff(:,1),8)

       moment(:) = real(dipolearrays(:,1))**2 + real(dipolearrays(:,2))**2 + real(dipolearrays(:,3))**2

       do i=1,npulses
          each_efield_ang(:,1,i) = real(dipolearrays(:,2) * each_efield(:,3,i) - dipolearrays(:,3) * each_efield(:,2,i),8)
          each_efield_ang(:,2,i) = real(dipolearrays(:,3) * each_efield(:,1,i) - dipolearrays(:,1) * each_efield(:,3,i),8)
          each_efield_ang(:,3,i) = real(dipolearrays(:,1) * each_efield(:,2,i) - dipolearrays(:,2) * each_efield(:,1,i),8)
       enddo

       angworksum0=0
       do ipulse=1,npulses

!! real-valued variables all
          if (moment(0).ne.0) then
             angworksum0(0,:,ipulse) = (-1) * dipole_ang(0,:) * each_efield_ang(0,:,ipulse) / moment(0) * par_timestep * autosteps
          endif
          do i=1,numdata
             if (moment(i).ne.0) then
                angworksum0(i,:,ipulse) = angworksum0(i-1,:,ipulse) - dipole_ang(i,:) * each_efield_ang(i,:,ipulse) / moment(i) &
                     * par_timestep * autosteps
             else
                angworksum0(i,:,ipulse) = angworksum0(i-1,:,ipulse)
             endif
          enddo
       enddo
       totangworksum0(:,:)=0d0
       do ipulse=1,npulses
          totangworksum0(:,:)=totangworksum0(:,:)+angworksum0(:,:,ipulse)
       enddo

       if (myrank.eq.1) then
          do ii=1,3
             open(171,file=outangworknames(ii),status="unknown",iostat=myiostat)
             call checkiostat(myiostat,"opening "//outangworknames(ii))
             do i=0,numdata

!! not worrying about complex values

                write(171,'(A25,F10.5,400F20.10)') " EACH PULSE WORK T= ", i*par_timestep*autosteps,&
                     totangworksum0(i,ii),angworksum0(i,ii,:)

             enddo
             close(171)
          enddo

          if (sflag.ne.0) then
             do ii=1,3
                open(171,file=outangworknames(ii)(1:getlen(outangworknames(ii)))//number(2:7),status="unknown",iostat=myiostat)
                call checkiostat(myiostat,"opening "//outangworknames(ii)(1:getlen(outangworknames(ii)))//number(2:7))
                do i=0,numdata

!! not worrying about complex values

                   write(171,'(A25,F10.5,400F20.10)') " EACH PULSE WORK T= ", i*par_timestep*autosteps,&
                        totangworksum0(i,ii),angworksum0(i,ii,:)

                enddo
                close(171)
             enddo
          endif
       endif
       call mpibarrier()
       deallocate(angworksum0,totangworksum0,dipole_ang, each_efield_ang,moment)
    endif
    deallocate(dipole_diff)

    if (act21circ.ne.0) then
       numft=9
    else
       numft=3
    endif

    allocate(fftrans(0:numdata,numft), eft(0:numdata,numft), each_eft(0:numdata,numft,npulses))
    fftrans=0; eft=0; each_eft=0

    fftrans(:,1:3)=dipolearrays(:,:)
    eft(:,1:3)=efield(:,:) 
    each_eft(:,1:3,:)=each_efield(:,:,:)

    deallocate(efield,each_efield,dipolearrays)

    do i=0,numdata
       fftrans(i,1:3) = fftrans(i,1:3) * windowfunct(i,numdata,21)  !! action 21
    enddo

    if (pulsewindowtoo.ne.0) then
       do i=0,numdata
          each_eft(i,1:3,:)=each_eft(i,1:3,:) * windowfunct(i,numdata,21)
          eft(i,1:3)=eft(i,1:3) * windowfunct(i,numdata,21)
       enddo
    endif

    do ii=1,3
       call zfftf_wrap_diff(numdata+1,fftrans(:,ii),ftdiff)
       call zfftf_wrap(numdata+1,eft(:,ii))
       do ipulse=1,npulses
          call zfftf_wrap(numdata+1,each_eft(:,ii,ipulse))
       enddo
    enddo

    fftrans(:,1:3)=fftrans(:,1:3)     * par_timestep * autosteps
    eft(:,1:3)=eft(:,1:3)             * par_timestep * autosteps
    each_eft(:,1:3,:)=each_eft(:,1:3,:) * par_timestep * autosteps

    if (act21circ.ne.0) then
!! XY, XZ, YX, YZ, ZX, ZY

       fftrans(:,4) = ( fftrans(:,1) + (0d0,1d0) * fftrans(:,2) )/sqrt(2d0)
       fftrans(:,5) = ( fftrans(:,1) + (0d0,1d0) * fftrans(:,3) )/sqrt(2d0)
       fftrans(:,6) = ( fftrans(:,2) + (0d0,1d0) * fftrans(:,1) )/sqrt(2d0)
       fftrans(:,7) = ( fftrans(:,2) + (0d0,1d0) * fftrans(:,3) )/sqrt(2d0)
       fftrans(:,8) = ( fftrans(:,3) + (0d0,1d0) * fftrans(:,1) )/sqrt(2d0)
       fftrans(:,9) = ( fftrans(:,3) + (0d0,1d0) * fftrans(:,2) )/sqrt(2d0)

       eft(:,4) = ( eft(:,1) + (0d0,1d0) * eft(:,2) )/sqrt(2d0)
       eft(:,5) = ( eft(:,1) + (0d0,1d0) * eft(:,3) )/sqrt(2d0)
       eft(:,6) = ( eft(:,2) + (0d0,1d0) * eft(:,1) )/sqrt(2d0)
       eft(:,7) = ( eft(:,2) + (0d0,1d0) * eft(:,3) )/sqrt(2d0)
       eft(:,8) = ( eft(:,3) + (0d0,1d0) * eft(:,1) )/sqrt(2d0)
       eft(:,9) = ( eft(:,3) + (0d0,1d0) * eft(:,2) )/sqrt(2d0)

       each_eft(:,4,:) = ( each_eft(:,1,:) + (0d0,1d0) * each_eft(:,2,:) )/sqrt(2d0)
       each_eft(:,5,:) = ( each_eft(:,1,:) + (0d0,1d0) * each_eft(:,3,:) )/sqrt(2d0)
       each_eft(:,6,:) = ( each_eft(:,2,:) + (0d0,1d0) * each_eft(:,1,:) )/sqrt(2d0)
       each_eft(:,7,:) = ( each_eft(:,2,:) + (0d0,1d0) * each_eft(:,3,:) )/sqrt(2d0)
       each_eft(:,8,:) = ( each_eft(:,3,:) + (0d0,1d0) * each_eft(:,1,:) )/sqrt(2d0)
       each_eft(:,9,:) = ( each_eft(:,3,:) + (0d0,1d0) * each_eft(:,2,:) )/sqrt(2d0)

    endif

    Estep=2*pi/par_timestep/autosteps/(numdata+1)

    allocate(&
         worksums(0:numdata,numft,npulses),photsums(0:numdata,numft,npulses),&
         totworksums(0:numdata,numft), totphotsums(0:numdata,numft), &
         workpers(0:numdata,numft,npulses),photpers(0:numdata,numft,npulses),&
         totworkpers(0:numdata,numft), totphotpers(0:numdata,numft), &
         sumrule(0:numdata,numft))
    worksums=0; photsums=0; totworksums=0; totphotsums=0; 
    workpers=0; photpers=0; totworkpers=0; totphotpers=0; 
    sumrule=0

    if (dipolesumstart.le.0d0) then
       do ipulse=1,npulses
          photpers(0,:,ipulse) = imag(fftrans(0,:)*conjg(each_eft(0,:,ipulse))) / PI / 2   !! /2 6/16
          photsums(0,:,ipulse) = Estep * photpers(0,:,ipulse)
       enddo
    endif
    do i=1,numdata
       myenergy=i*Estep

!! sumrule sums to N for N electrons
       if (myenergy.ge.dipolesumstart.and.myenergy.le.dipolesumend) then
          sumrule(i,:)=sumrule(i-1,:) + Estep * imag(fftrans(i,:)*conjg(eft(i,:))) / abs(eft(i,:)**2) * myenergy * 2 / PI
          do ipulse=1,npulses
             photpers(i,:,ipulse) = imag(fftrans(i,:)*conjg(each_eft(i,:,ipulse))) / PI / 2  !! /2 6/16
             workpers(i,:,ipulse) = imag(fftrans(i,:)*conjg(each_eft(i,:,ipulse))) / PI * myenergy
             photsums(i,:,ipulse) = photsums(i-1,:,ipulse) + Estep * photpers(i,:,ipulse)
             worksums(i,:,ipulse) = worksums(i-1,:,ipulse) + Estep * workpers(i,:,ipulse)
          enddo
       else
          sumrule(i,:)=sumrule(i-1,:)
          photsums(i,:,:)  =  photsums(i-1,:,:)
          worksums(i,:,:)=worksums(i-1,:,:)
       endif
    enddo
    totphotsums=0
    totworksums=0
    totphotpers=0
    totworkpers=0
    do ipulse=1,npulses
       totphotsums(:,:)=totphotsums(:,:)+photsums(:,:,ipulse)
       totworksums(:,:)=totworksums(:,:)+worksums(:,:,ipulse)
       totphotpers(:,:)=totphotpers(:,:)+photpers(:,:,ipulse)
       totworkpers(:,:)=totworkpers(:,:)+workpers(:,:,ipulse)
    enddo

    if (myrank.eq.1) then
       do ii=1,numft
          open(171,file=outftnames(ii),status="unknown",iostat=myiostat)
          call checkiostat(myiostat,"opening "//outftnames(ii))
          write(171,'(A120)',iostat=myiostat) &
               "## Photon energy (column 1); D(omega) (2,3); E(omega) (4,5); response (6,7); cross sect (9); integrated (10)" 
          call checkiostat(myiostat,"writing "//outftnames(ii))
          write(171,'(A120)') "## UNITLESS RESPONSE FUNCTION FOR ABSORPTION/EMISSION 2 omega im(D(omega)E(omega)^*) IN COLUMN 7"
          write(171,'(A120)') "## QUANTUM MECHANICAL PHOTOABSORPTION/EMISSION CROSS SECTION IN MEGABARNS (no factor of 1/3) IN COLUMN NINE"
          write(171,'(A120)') "## INTEGRATED DIFFERENTIAL OSCILLATOR STRENGTH (FOR SUM RULE) IN COLUMN 10"
          write(171,*)
       
          do i=0,numdata
             myenergy=i*Estep

!! LENGTH GAUGE (electric field) WAS FT'ed , OK with usual formula multiply by wfi
!! UNITLESS RESPONSE FUNCTION FOR ABSORPTION/EMISSION 2 omega im(D(omega)E(omega)^*) IN COLUMN 7
!! QUANTUM MECHANICAL PHOTOABSORPTION/EMISSION CROSS SECTION IN MEGABARNS (no factor of 1/3) IN COLUMN NINE
!! INTEGRATED DIFFERENTIAL OSCILLATOR STRENGTH (FOR SUM RULE) IN COLUMN 10

             xsecunits = 5.291772108d0**2 * 4d0 * PI / 1.37036d2 * myenergy

!! NOW FACTOR (2 omega) IN COLUMNS 6,7   v1.16 12-2015

             write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat)  myenergy, &
                  fftrans(i,ii), eft(i,ii), fftrans(i,ii)*conjg(eft(i,ii)) * 2 * myenergy, &
                  fftrans(i,ii)*conjg(eft(i,ii)) / abs(eft(i,ii)**2) * xsecunits, sumrule(i,ii)
          enddo
          call checkiostat(myiostat,"writing "//outftnames(ii))
          close(171)

!!  NUMBER OF PHOTONS ABSORBED AND AND WORK DONE BY EACH PULSE
!!  worksum0 the time integral converges right after pulse is finished... others take longer

          open(171,file=outoworknames(ii),status="unknown",iostat=myiostat)
          call checkiostat(myiostat,"opening "//outoworknames(ii))
          do i=0,numdata
             myenergy=i*Estep
             if (myenergy.ge.dipolesumstart.and.myenergy.le.dipolesumend) then
                write(171,'(A25,F10.5,400F20.10)') " WORK EACH PULSE E= ", myenergy, &
                     totworksums(i,ii), totworkpers(i,ii), worksums(i,ii,:), workpers(i,ii,:)
             endif
          enddo
          close(171)

          open(171,file=outophotonnames(ii),status="unknown",iostat=myiostat)
          call checkiostat(myiostat,"opening "//outophotonnames(ii))
          do i=0,numdata
             myenergy=i*Estep
             if (myenergy.ge.dipolesumstart.and.myenergy.le.dipolesumend) then
                write(171,'(A25,F10.5,400F20.10)') "PHOTONS EACH PULSE E= ", myenergy, &
                     totphotsums(i,ii), totphotpers(i,ii), photsums(i,ii,:), photpers(i,ii,:)
             endif
          enddo
          close(171)

          if (sflag.ne.0) then

             open(171,file=outoworknames(ii)(1:getlen(outoworknames(ii)))//number(2:7),status="unknown",iostat=myiostat)
             call checkiostat(myiostat,"opening "//outoworknames(ii)(1:getlen(outoworknames(ii)))//number(2:7))
             do i=0,numdata
                myenergy=i*Estep
                if (myenergy.ge.dipolesumstart.and.myenergy.le.dipolesumend) then
                   write(171,'(A25,F10.5,400F20.10)') " WORK EACH PULSE E= ", myenergy, &
                        totworksums(i,ii), totworkpers(i,ii), worksums(i,ii,:), workpers(i,ii,:)
                endif
             enddo
             close(171)

             open(171,file=outophotonnames(ii)(1:getlen(outophotonnames(ii)))//number(2:7),status="unknown",iostat=myiostat)
             call checkiostat(myiostat,"opening "//outophotonnames(ii)(1:getlen(outophotonnames(ii)))//number(2:7))
             do i=0,numdata
                myenergy=i*Estep
                if (myenergy.ge.dipolesumstart.and.myenergy.le.dipolesumend) then
                   write(171,'(A25,F10.5,400F20.10)') "PHOTONS EACH PULSE E= ", myenergy, &
                        totphotsums(i,ii), totphotpers(i,ii), photsums(i,ii,:), photpers(i,ii,:)
                endif
             enddo
             close(171)

             open(171,file=outftnames(ii)(1:getlen(outftnames(ii)))//number(2:7),status="unknown",iostat=myiostat)
             call checkiostat(myiostat,"opening "//outftnames(ii)(1:getlen(outftnames(ii)))//number(2:7))
             write(171,'(A120)',iostat=myiostat) &
                  "## Photon energy (column 1); D(omega) (2,3); E(omega) (4,5); response (6,7); cross sect (9); integrated (10)" 
             call checkiostat(myiostat,"writing "//outftnames(ii)(1:getlen(outftnames(ii)))//number(2:7))
             write(171,'(A120)') "## UNITLESS RESPONSE FUNCTION FOR ABSORPTION/EMISSION 2 omega im(D(omega)E(omega)^*) IN COLUMN 7"
             write(171,'(A120)') "## QUANTUM MECHANICAL PHOTOABSORPTION/EMISSION CROSS SECTION IN MEGABARNS (no factor of 1/3) IN COLUMN NINE"
             write(171,'(A120)') "## INTEGRATED DIFFERENTIAL OSCILLATOR STRENGTH (FOR SUM RULE) IN COLUMN 10"
             write(171,*)

             do i=0,numdata
                myenergy=i*Estep

                xsecunits = 5.291772108d0**2 * 4d0 * PI / 1.37036d2 * myenergy

!! NOW FACTOR (2 omega) IN COLUMNS 6,7   v1.16 12-2015

                write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat)  myenergy, &
                     fftrans(i,ii), eft(i,ii), fftrans(i,ii)*conjg(eft(i,ii)) * 2 * myenergy, &
                     fftrans(i,ii)*conjg(eft(i,ii)) / abs(eft(i,ii)**2) * xsecunits, sumrule(i,ii)
             enddo
             call checkiostat(myiostat,"writing "//outftnames(ii)(1:getlen(outftnames(ii)))//number(2:7))
             close(171)
          endif
       enddo    !! do ii=1,numft
    endif  !! myrank
    call mpibarrier()

    deallocate(fftrans,eft,each_eft)
    deallocate(worksums,photsums,totworksums,totphotsums,  workpers,photpers,totworkpers,totphotpers,&
         sumrule)

  end subroutine dipolecall

end module dipcallsubmod


   subroutine dipolesub0(inspfs,inavector)
     use dipolemod
     use dipsubonemod
     use dipcallsubmod
     use parameters
     use configmod
     use mpisubmod
     use pulse_parameters !! conjgpropflag
     implicit none
     DATATYPE,intent(in) :: inspfs(totspfdim), inavector(tot_adim,mcscfnum)
     DATATYPE :: myexpects(3), mcexpects(3,mcscfnum), dd(mcscfnum),&
          axx(mcscfnum),ayy(mcscfnum),azz(mcscfnum),sxx(mcscfnum),syy(mcscfnum),&
          szz(mcscfnum),drivingoverlap(mcscfnum)
     character(len=SLN) :: dipfiles(3), tworkfiles(3), angworkfiles(3), ftfiles(9), &
          oworkfiles(9), ophotonfiles(9)
     integer :: imc,sflag,getlen,ii
     integer, save :: lastouttime=0, calledflag=0
     real*8 :: thistime
     character(len=2) :: tl(4) = (/ "BA", "AB", "AA", "BB" /)

     myexpects=0;mcexpects=0;axx=0;ayy=0;azz=0;sxx=0;syy=0;szz=0;dd=0;drivingoverlap=0

     if (conjgpropflag.ne.0) then

        if (mcscfnum.ne.2) then
           OFLWR "Whoot? conjgpropflag mcscfnum",mcscfnum; CFLST
        endif
        if (drivingflag.ne.0) then
           OFLWR "Driving not supported for conjprop yet"; CFLST
        endif
        dipolenormsq(calledflag)=0
        if (tot_adim.gt.0) then
           dipolenormsq(calledflag) = hermdot(inavector(:,2),inavector(:,1),tot_adim)
        endif
        if (par_consplit.ne.0) then
           call mympireduceone(dipolenormsq(calledflag))
        endif

        OFLWR "   complex Domcke - off diagonal norm-squared ", dipolenormsq(calledflag)

        call dipolesub_one(www,bioww,inavector(:,2),inavector(:,1), inspfs(:), myexpects(:))
        dipoleexpects(calledflag,:,1)=myexpects(:)
        call dipolesub_one(www,bioww,inavector(:,1),inavector(:,2), inspfs(:), myexpects(:))
        dipoleexpects(calledflag,:,2)=myexpects(:)
        call dipolesub_one(www,bioww,inavector(:,1),inavector(:,1), inspfs(:), myexpects(:))
        dipoleexpects(calledflag,:,3)=myexpects(:)
        call dipolesub_one(www,bioww,inavector(:,2),inavector(:,2), inspfs(:), myexpects(:))
        dipoleexpects(calledflag,:,4)=myexpects(:)

     else  !! conjgpropflag complex Domcke

        do imc=1,mcscfnum
           dd(imc)=0
           if (tot_adim.gt.0) then
              dd(imc) = hermdot(inavector(:,imc),inavector(:,imc),tot_adim)
           endif
           if (par_consplit.ne.0) then
              call mympireduceone(dd(imc))
           endif
           call dipolesub_one(www,bioww,inavector(:,imc),inavector(:,imc),inspfs(:),mcexpects(:,imc))
        enddo

        if (drivingflag.ne.0) then
#ifdef CNORMFLAG
           OFLWR "Error, driving with dipole not supported c-norm"; CFLST
#endif
           call dipolesub_driving(axx,ayy,azz,sxx,syy,szz,mcscfnum)
           mcexpects(1,:)=mcexpects(1,:)+axx(:)+CONJUGATE(axx(:))+sxx*drivingproportion**2
           mcexpects(2,:)=mcexpects(2,:)+ayy(:)+CONJUGATE(ayy(:))+syy*drivingproportion**2
           mcexpects(3,:)=mcexpects(3,:)+azz(:)+CONJUGATE(azz(:))+szz*drivingproportion**2

!! for time slot zero
           call getdrivingoverlap(drivingoverlap,mcscfnum)
           dd(:)=dd(:)+drivingproportion**2 + drivingoverlap(:) + CONJUGATE(drivingoverlap(:))
        endif

        dipoleexpects(calledflag,:,1)=0d0
        dipolenormsq(calledflag)=0d0

        do imc=1,mcscfnum

           dipolenormsq(calledflag) = dipolenormsq(calledflag) + dd(imc)

!! 101414 REAL-VALUED FOR HERM.
!! 1-2016 v1.17 should not be necessary with realflag in mult_zdipole(in,out,realflag) etc.
 
#ifndef CNORMFLAG
           dipoleexpects(calledflag,:,1) = dipoleexpects(calledflag,:,1) + real(mcexpects(:,imc),8) / mcscfnum
#else
           dipoleexpects(calledflag,:,1) = dipoleexpects(calledflag,:,1) + mcexpects(:,imc) / mcscfnum
#endif

        enddo

     endif  !! conjgpropflag complex Domcke

     if (mod(calledflag,dipmodtime).eq.0.and.calledflag.gt.0) then

        OFLWR "Go emission/absorption action 21... "; CFL

        thistime=calledflag*par_timestep*autosteps
        sflag=0
        if (floor(thistime/diptime).gt.lastouttime) then
           lastouttime=floor(thistime/diptime)
           sflag=1
        endif

!!$  subroutine dipolecall(numdata, indipolearrays, outenames, outtworknames, outangworknames, &
!!$       outftnames, outoworknames, outophotonnames, sflag)

        if (conjgpropflag.eq.0) then
           dipfiles(:) = (/ xdipfile, ydipfile, zdipfile /)
           tworkfiles(:) = (/ xtworkfile, ytworkfile, ztworkfile /)
           angworkfiles(:) = (/ yztworkfile, zxtworkfile, xytworkfile /)
           ftfiles(:) = (/ xdftfile, ydftfile, zdftfile, &
                xydftfile, xzdftfile, yxdftfile, yzdftfile, zxdftfile, zydftfile /)
           oworkfiles(:) = (/ xoworkfile, yoworkfile, zoworkfile, &
                xyoworkfile, xzoworkfile, yxoworkfile, yzoworkfile, zxoworkfile, zyoworkfile /)
           ophotonfiles(:) = (/ xophotonfile, yophotonfile, zophotonfile, &
                xyophotonfile, xzophotonfile, yxophotonfile, yzophotonfile, zxophotonfile, zyophotonfile /)

           if (reference_pulses.eq.0) then
              call dipolecall(calledflag, dipoleexpects(:,:,1), &
                   dipfiles, tworkfiles, angworkfiles, ftfiles, oworkfiles, ophotonfiles, sflag, 0, numpulses)
           else
              call dipolecall(calledflag, dipoleexpects(:,:,1), &
                   dipfiles, tworkfiles, angworkfiles, ftfiles, oworkfiles, ophotonfiles, sflag, 1, reference_pulses)
           endif

        else   !! conjgpropflag complex Domcke:
           do ii=1,4
              dipfiles(:) = (/ &
                   xdipfile(1:getlen(xdipfile))//tl(ii),&
                   ydipfile(1:getlen(ydipfile))//tl(ii),&
                   zdipfile(1:getlen(zdipfile))//tl(ii) /)
              tworkfiles(:) = (/ &
                   xtworkfile(1:getlen(xtworkfile))//tl(ii),&
                   ytworkfile(1:getlen(ytworkfile))//tl(ii),&
                   ztworkfile(1:getlen(ztworkfile))//tl(ii) /)
              angworkfiles(:) = (/ &
                   yztworkfile(1:getlen(yztworkfile))//tl(ii),&
                   zxtworkfile(1:getlen(zxtworkfile))//tl(ii),&
                   xytworkfile(1:getlen(xytworkfile))//tl(ii) /)
              ftfiles(:) = (/ &
                   xdftfile(1:getlen(xdftfile))//tl(ii),&
                   ydftfile(1:getlen(ydftfile))//tl(ii),&
                   zdftfile(1:getlen(zdftfile))//tl(ii),&
                   xydftfile(1:getlen(xydftfile))//tl(ii),&
                   xzdftfile(1:getlen(xzdftfile))//tl(ii),&
                   yxdftfile(1:getlen(yxdftfile))//tl(ii),&
                   yzdftfile(1:getlen(yzdftfile))//tl(ii),&
                   zxdftfile(1:getlen(zxdftfile))//tl(ii),&
                   zydftfile(1:getlen(zydftfile))//tl(ii) /)
              oworkfiles(:) = (/ &
                   xoworkfile(1:getlen(xoworkfile))//tl(ii),&
                   yoworkfile(1:getlen(yoworkfile))//tl(ii),&
                   zoworkfile(1:getlen(zoworkfile))//tl(ii),&
                   xyoworkfile(1:getlen(xyoworkfile))//tl(ii),&
                   xzoworkfile(1:getlen(xzoworkfile))//tl(ii),&
                   yxoworkfile(1:getlen(yxoworkfile))//tl(ii),&
                   yzoworkfile(1:getlen(yzoworkfile))//tl(ii),&
                   zxoworkfile(1:getlen(zxoworkfile))//tl(ii),&
                   zyoworkfile(1:getlen(zyoworkfile))//tl(ii) /)
              ophotonfiles(:) = (/ &
                   xophotonfile(1:getlen(xophotonfile))//tl(ii),&
                   yophotonfile(1:getlen(yophotonfile))//tl(ii),&
                   zophotonfile(1:getlen(zophotonfile))//tl(ii),&
                   xyophotonfile(1:getlen(xyophotonfile))//tl(ii),&
                   xzophotonfile(1:getlen(xzophotonfile))//tl(ii),&
                   yxophotonfile(1:getlen(yxophotonfile))//tl(ii),&
                   yzophotonfile(1:getlen(yzophotonfile))//tl(ii),&
                   zxophotonfile(1:getlen(zxophotonfile))//tl(ii),&
                   zyophotonfile(1:getlen(zyophotonfile))//tl(ii) /)

              if (reference_pulses.eq.0) then
                 call dipolecall(calledflag, dipoleexpects(:,:,ii), &
                      dipfiles, tworkfiles, angworkfiles, ftfiles, oworkfiles, ophotonfiles, sflag, 0, numpulses)
              else
                 call dipolecall(calledflag, dipoleexpects(:,:,ii), &
                      dipfiles, tworkfiles, angworkfiles, ftfiles, oworkfiles, ophotonfiles, sflag, 1, reference_pulses)
              endif

           enddo !! do ii=1,4
        endif    !! conjgpropflag

        OFLWR "     ..done emission/absorption"; CFL

     endif       !! calledflag (dipmodtime)

     calledflag=calledflag+1

   end subroutine dipolesub0


subroutine dipolesub()   !! action 21
  use parameters
  use xxxmod
  implicit none
  integer,save :: xcalledflag=0

  if (mod(xcalledflag,autosteps).eq.0) then
     call dipolesub0(yyy%cmfspfs(:,0),yyy%cmfavec(:,:,0))
  endif
  xcalledflag=xcalledflag+1

end subroutine dipolesub


subroutine redo_dipolesub(alg)
  use dipolemod
  use dipsubonemod
  use dipcallsubmod
  use parameters
  use configmod
  use xxxmod
  use mpisubmod
  use pulse_parameters !! conjgpropflag
  use mpimod
  implicit none
  integer, intent(in) :: alg
  integer :: k,nt,i,molength,alength,  BatchSize,NBat,ketbat,ketreadsize, &
       kettime,imc, ispf, myiostat
  real*8 :: MemTot,MemVal, dt, MemNum
  DATATYPE :: nullvector(numr)
  DATATYPE, allocatable :: ketmo(:,:,:),ketavec(:,:,:,:),&
       read_ketmo(:,:),read_ketavec(:,:,:,:)

  nullvector=0

  dt=real(FluxInterval*FluxSkipMult,8)*par_timestep; 
  nt=floor(real(numpropsteps,8)/fluxinterval/fluxskipmult)

  if (abs(dt - par_timestep*autosteps)/dt .gt. 1d-3) then
     OFLWR "error, need fluxskipmult*fluxinterval*par_timestep.eq.autotimestep for action 29"
     WRFL "   fluxskipmult, fluxinterval, par_timestep, autosteps, autotimestep "
     WRFL fluxskipmult, fluxinterval, par_timestep, autosteps, autotimestep; CFLST
  endif

  call mpibarrier()
  OFLWR "Go recompute dipolemoment"; CFL

!! determine if we should do batching or not
!! 250,000 words/MB, real*8 2words/#, complex*16 4words/#

#ifdef REALGO
  MemVal = 1.25d5
#else
  MemVal = 6.25d4
#endif

!! memory based on rank 1 which right now stores the entire orbital vector, which is wasteful.
!! Unsatisfactory: need better MPI I/O.

  MemNum=0d0

!!$    if (parorbsplit.eq.3) then
!!$       MemNum = MemNum + spfsize*nspf*nprocs
!!$    else

  MemNum = MemNum + spfsize*nspf

!!$    endif

  if (par_consplit.ne.0) then
     MemNum = MemNum + www%numconfig*numr*mcscfnum
  else
     MemNum = MemNum + www%maxconfigsperproc*numr*mcscfnum
  endif

  call openfile()
  write(mpifileptr,'(A30,F9.3,A3)') " Guess at necessary memory is ",&
       (nt+1)*MemNum/MemVal," MB"
  if(alg.eq.0) then
     write(mpifileptr,*) "   ...will be computed with all of psi in core"
     BatchSize=nt+1
  else
     MemTot=real(alg,8)    
     write(mpifileptr,*) "    ... computed with all psi being read in batches"
     write(mpifileptr,'(A33,F9.3,A3)') "        Desired amount of memory ",MemTot," MB"
     BatchSize=floor(MemTot * MemVal / MemNum)
     if(BatchSize.lt.1) then
        write(mpifileptr,*) "Tiny amount of memory or huge wavefunction, Batchsize is 1" 
        BatchSize=1
     else if(BatchSize.ge.nt+1) then
        write(mpifileptr,*) "Good, there is enough memory for only one batch."
        BatchSize=nt+1
     else
        write(mpifileptr,*) "Batchsize is ",BatchSize,"/",(nt+1)
     endif
  endif
  call closefile()

  call mpibarrier()
  OFLWR "Allocating psi arrays for recomputing dipolemoment"; CFL

  allocate(ketmo(spfsize,nspf,BatchSize), &
       ketavec(numr,www%firstconfig:www%lastconfig,mcscfnum,BatchSize))
  ketmo=0
  if (www%lastconfig.ge.www%firstconfig) then
     ketavec=0
  endif
  if (myrank.eq.1) then
     if (parorbsplit.eq.3) then
        allocate(read_ketmo(spfsize*nprocs,nspf))
     else
        allocate(read_ketmo(spfsize,nspf))
     endif
  else
     allocate(read_ketmo(1,nspf))
  endif
  read_ketmo=0

  if (myrank.eq.1) then
     allocate(read_ketavec(numr,www%numconfig,mcscfnum,BatchSize))
  else
     allocate(read_ketavec(1,1,mcscfnum,BatchSize))
  endif
  read_ketavec=0

  NBat=ceiling(real(nt+1)/real(BatchSize))
  ketreadsize=0

  if (myrank.eq.1) then
     inquire (iolength=molength) read_ketmo(:,:)
     inquire (iolength=alength) read_ketavec(:,:,:,1)
  endif
  call mympiibcastone(molength,1); call mympiibcastone(alength,1)

  call openfile()
  write(mpifileptr,*) "MO record length is ",molength
  write(mpifileptr,*) "AVEC record length is ",alength
  call closefile()

!! begin the ket batch read loop
  do ketbat=1,NBat

     OFLWR "Reading batch ", ketbat, " of ", NBat; CFL
     ketreadsize=min(BatchSize,nt+1-(ketbat-1)*BatchSize)

!! ORBITALS
     if(myrank.eq.1) then
        open(1001,file=fluxmofile,status="old",form="unformatted",&
             access="direct",recl=molength,iostat=myiostat)
        call checkiostat(myiostat,"opening "//fluxmofile)
     endif
     do i=1,ketreadsize
        if(myrank.eq.1) then
           k=FluxSkipMult*((ketbat-1)*BatchSize+i-1)+1
           read(1001,rec=k,iostat=myiostat) read_ketmo(:,:) 
           call checkiostat(myiostat,"reading "//fluxmofile)
        endif
        if (parorbsplit.ne.3) then
           if (myrank.eq.1) then
              ketmo(:,:,i)=read_ketmo(:,:)
           endif
           call mympibcast(ketmo(:,:,i),1,totspfdim)
        else
           do ispf=1,nspf
              call splitscatterv(read_ketmo(:,ispf),ketmo(:,ispf,i))
           enddo
        endif
     enddo       !! do i=1,ketreadsize
     if (myrank.eq.1) then
        close(1001)
     endif
!! A-VECTOR
     if(myrank.eq.1) then
        open(1002,file=fluxafile,status="old",form="unformatted",&
             access="direct",recl=alength,iostat=myiostat)
        call checkiostat(myiostat,"opening "//fluxafile)
        do i=1,ketreadsize
           k=FluxSkipMult*((ketbat-1)*BatchSize+i-1)+1
           read(1002,rec=k,iostat=myiostat) read_ketavec(:,:,:,i) 
        enddo
        call checkiostat(myiostat,"reading "//fluxafile)
        close(1002)
     endif
     if (par_consplit.eq.0) then
        if (myrank.eq.1) then
           ketavec(:,:,:,1:ketreadsize)=read_ketavec(:,:,:,1:ketreadsize)
        endif
        call mympibcast(ketavec(:,:,:,1:ketreadsize),1,numr*www%numconfig*mcscfnum*ketreadsize)
     else
        do i=1,ketreadsize
           do imc=1,mcscfnum
              if (tot_adim.gt.0) then
                 call myscatterv(read_ketavec(:,:,imc,i),&
                      ketavec(:,:,imc,i),configs_perproc(:)*numr)
              else
                 call myscatterv(read_ketavec(:,:,imc,i),&
                      nullvector(:),configs_perproc(:)*numr)
              endif
           enddo
        enddo
     endif

!! loop over all time for the ket of the dipole integral
     do kettime=1,ketreadsize

        call dipolesub0(ketmo(:,:,kettime),ketavec(:,:,:,kettime))

     enddo

  enddo  !! ketbat

  deallocate(ketmo, ketavec, read_ketmo, read_ketavec)
  OFLWR "Done recomputing dipole, stopping"; CFLST

end subroutine redo_dipolesub
