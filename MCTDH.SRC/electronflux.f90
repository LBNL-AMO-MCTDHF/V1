
!! ALL MODULES

!! ELECTRON FLUX ROUTINES - ACTION 16 analysis total photoionization, 
!! ACTION 28 photoionization flux during calculation, ACTION 15 
!! save file, saved file may be used for actions 16 17 23)


#include "Definitions.INC"

module fluxutilmod
contains

subroutine mult_reke(howmany,in,out)
  use parameters
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  DATATYPE :: work(spfsize,howmany)
  work(:,:)=ALLCON(in(:,:))
  call mult_ke(work,out,howmany,"booga",2)
  work=ALLCON(out)
  call mult_ke(in,out,howmany,"booga",2)  
  work=work+out
  out=work/2d0
end subroutine mult_reke


subroutine mult_imke(howmany,in,out)
  use parameters
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  DATATYPE :: work(spfsize,howmany)
  work(:,:)=ALLCON(in(:,:))
  call mult_ke(work,out,howmany,"booga",2)
  work=ALLCON(out)
  call mult_ke(in,out,howmany,"booga",2)  
  work=out-work
  out=work/(0d0,2d0)
end subroutine mult_imke


subroutine op_reyderiv(howmany,in,out)
  use parameters
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  DATATYPE :: work(spfsize,howmany)
  work(:,:)=ALLCON(in(:,:))
  call op_yderiv(howmany,work,out)
  work=ALLCON(out)
  call op_yderiv(howmany,in,out)
  work=work+out
  out=work/2d0
end subroutine op_reyderiv


subroutine op_imyderiv(howmany,in,out)
  use parameters
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  DATATYPE :: work(spfsize,howmany)
  work(:,:)=ALLCON(in(:,:))
  call op_yderiv(howmany,work,out)
  work=ALLCON(out)
  call op_yderiv(howmany,in,out)
  work=out-work
  out=work/(0d0,2d0)
end subroutine op_imyderiv


!! needs factor of 1/r  for hamiltonian

subroutine mult_impot(howmany,in, out)
  use parameters
  use opmod 
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  integer :: ii
  do ii=1,howmany
     out(:,ii)=in(:,ii)*imag((0d0,0d0)+pot(:))   !! NO INTERNUCLEAR REPULSION !!
  enddo
end subroutine mult_impot



!! needs factor of 1/r  for hamiltonian

subroutine mult_repot(howmany,in, out)
  use parameters
  use opmod 
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  integer :: ii
  do ii=1,howmany
     out(:,ii)=in(:,ii)*real(pot(:),8)   !! NO INTERNUCLEAR REPULSION !!
  enddo
end subroutine mult_repot


subroutine mult_imhalfniumpot(howmany,in, out)
  use parameters
  use opmod  
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  integer :: ii
  do ii=1,howmany
     out(:,ii)=in(:,ii)*imag((0d0,0d0)+halfniumpot(:))
  enddo
end subroutine mult_imhalfniumpot


subroutine mult_rehalfniumpot(howmany,in, out)
  use parameters
  use opmod  
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(spfsize,howmany)
  DATATYPE,intent(out) :: out(spfsize,howmany)
  integer :: ii
  do ii=1,howmany
     out(:,ii)=in(:,ii)*real(halfniumpot(:),8)
  enddo
end subroutine mult_rehalfniumpot

end module fluxutilmod


module fluxgtaubiomod
  use biorthotypemod
  implicit none
  type(biorthotype),target :: fluxgtaubiovar
end module fluxgtaubiomod


module fluxduringmod
  implicit none
  integer :: allocated=0, curtime=0
  DATATYPE,allocatable :: gtausum(:), gtausave(:)
end module fluxduringmod


module fluxgtau0mod
contains

  subroutine fluxgtau_during0(wwin,bbin,in_ketmo,ketavec,dt)
    use fluxgtaubiomod
    use biorthomod
    use parameters
    use walkmod
    use mpimod
    use mpisubmod
    use pulsesubmod
    use fluxduringmod
    use orbmultsubmod   !! gauge_transform
    implicit none
!! FluxOpType:
!! 0       = use one-e potential and two-e contribution routines  (exact treatment)
!! 1       = replace one-e potential + two-e potential with halfnium one-e potential  (recommended)
!! 2       = use full one-e potential, no two-e 
!! other   = only KE

    type(walktype),target,intent(in) :: wwin,bbin
    type(walktype),pointer :: myww
    real*8,intent(in) :: dt
    DATATYPE,intent(in) :: ketavec(numr,wwin%firstconfig:wwin%lastconfig,mcscfnum), &
         in_ketmo(spfsize,nspf)
    integer :: imc, myiostat
    DATATYPE :: gtaunow(mcscfnum), nullvector(numr)
    DATATYPE, allocatable :: imke(:,:),impe(:,:),imV2(:,:,:,:), imyderiv(:,:),&
         imkeop(:,:),impeop(:,:),  imyop(:,:), reke(:,:),repe(:,:),reV2(:,:,:,:), reyderiv(:,:),&
         rekeop(:,:),repeop(:,:),  reyop(:,:), ketmo(:,:)

    gtaunow=0; nullvector=0

!! initial setup

    if (FluxOpType.eq.0) then    !! exact expression with two-electron
       myww=>wwin
    else
       myww=>bbin
    endif

    allocate(imke(nspf,nspf),impe(nspf,nspf),imV2(nspf,nspf,nspf,nspf),imyderiv(nspf,nspf),&
         reke(nspf,nspf),repe(nspf,nspf),reV2(nspf,nspf,nspf,nspf),reyderiv(nspf,nspf),&
         rekeop(spfsize,nspf),repeop(spfsize,nspf),reyop(spfsize,nspf),&
         imkeop(spfsize,nspf),impeop(spfsize,nspf),imyop(spfsize,nspf), ketmo(spfsize,nspf))

    imke=0; impe=0; imV2=0; imyderiv=0; imkeop=0; impeop=0; imyop=0;
    reke=0; repe=0; reV2=0; reyderiv=0; rekeop=0; repeop=0; reyop=0;
    ketmo=0

    if (allocated.eq.0) then
       allocated=1
       allocate(gtausum(mcscfnum),gtausave(mcscfnum))
       gtausum=0d0; gtausave=0
       if (myrank.eq.1) then
          open(454, file=fluxtsumfile, status="unknown",iostat=myiostat)
          call checkiostat(myiostat,"opening fluxtsumfile")
          write(454,*,iostat=myiostat) "#KVL flux sum: time, flux integral, flux"
          call checkiostat(myiostat,"writing fluxtsumfile")
          write(454,*)
          close(454)
       endif
    else
       curtime=curtime+1
    endif

!! gaugefluxflag available to attempt gauge-invariant calculation with photoionization
!! during the pulse.  Best performance seems to be calculate wfn in length gauge, transform
!! to velocity gauge for flux: velflag=0, gaugefluxflag=2 (Volkov phase is only a function of 
!! time in the velocity gauge, no preferred origin in the velocity gauge)

!! gaugefluxflag=1, transform velocity to length; gaugefluxflag=2, transform length to velocity
    if ((gaugefluxflag.eq.1.and.velflag.ne.0).or.(gaugefluxflag.eq.2.and.velflag.eq.0)) then
       call gauge_transform(velflag,curtime*dt,nspf,in_ketmo(:,:),ketmo(:,:))
    else
       ketmo(:,:)=in_ketmo(:,:)
    endif

!! get the one-e matrix elements for this ket time

    imkeop=0; impeop=0;  imyop=0
    if (nucfluxopt.ne.2) then
       call flux_op_onee(ketmo,imkeop,impeop,1)
       if (nonuc_checkflag.eq.0) then
          call flux_op_nuc(ketmo,imyop,1)
       endif
    endif

    rekeop=0; repeop=0;
    if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0) then
       call flux_op_onee(ketmo,rekeop,repeop,2)
       call flux_op_nuc(ketmo,reyop,2)
    endif

!! complete the one-e potential and kinetic energy matrix elements           

    imke=0;  impe=0;  imyderiv=0;  imke=0;  imke=0; 

    if (nucfluxopt.ne.2) then
       call MYGEMM('C','N',nspf,nspf,spfsize,DATAONE,ketmo,&
            spfsize,imkeop,spfsize,DATAZERO,imke,nspf)
       call MYGEMM('C','N',nspf,nspf,spfsize,DATAONE,ketmo,&
            spfsize,impeop,spfsize,DATAZERO,impe,nspf)
       if (nonuc_checkflag.eq.0) then
          call MYGEMM('C','N',nspf,nspf,spfsize,DATAONE,ketmo,&
               spfsize,imyop,spfsize,DATAZERO,imyderiv,nspf)
       endif
    endif

    reke=0; repe=0; reyderiv=0

    if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0) then
       call MYGEMM('C','N',nspf,nspf,spfsize,DATAONE,ketmo,&
            spfsize,rekeop,spfsize,DATAZERO,reke,nspf)
       call MYGEMM('C','N',nspf,nspf,spfsize,DATAONE,ketmo,&
            spfsize,repeop,spfsize,DATAZERO,repe,nspf)
       call MYGEMM('C','N',nspf,nspf,spfsize,DATAONE,ketmo,&
            spfsize,reyop,spfsize,DATAZERO,reyderiv,nspf)
    endif

    if (parorbsplit.eq.3) then
       if (nucfluxopt.ne.2) then
          call mympireduce(imke,nspf**2)
          call mympireduce(impe,nspf**2)
          if (nonuc_checkflag.eq.0) then
             call mympireduce(imyderiv,nspf**2)
          endif
       endif
       if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0)then
          call mympireduce(reke,nspf**2)
          call mympireduce(repe,nspf**2)
          call mympireduce(reyderiv,nspf**2)
       endif
    endif

!! get the two-e contribution for exact formula (fluxoptype=0)

    imV2=0d0; reV2=0
    if(FluxOpType.eq.0) then
       if (nucfluxopt.ne.2) then
          call flux_op_twoe(ketmo,ketmo,imV2,1)
       endif
       if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0)then
          call flux_op_twoe(ketmo,ketmo,reV2,2)
       endif
    endif

!! evaluate the actual g(tau) expression

    do imc=1,mcscfnum
       if (tot_adim.gt.0) then
          gtaunow(imc) = fluxeval(ketavec(:,:,imc),ketavec(:,:,imc),&
               imke,impe,imV2,imyderiv,reke,repe,reV2,reyderiv,myww)
       else
          gtaunow(imc) = fluxeval(nullvector(:),nullvector(:),&
               imke,impe,imV2,imyderiv,reke,repe,reV2,reyderiv,myww)
       endif
    enddo

    if (curtime.eq.0) then
       gtausave(:)=gtaunow(:)
    endif

    if (flux_subtract.ne.0) then
       gtaunow(:) = gtaunow(:) - gtausave(:)
    endif

    gtausum(:) = gtausum(:) + gtaunow(:) * dt / 2d0   !! 2 looks correct

    if (myrank.eq.1) then
       open(454, file=fluxtsumfile, status="old",  position="append",iostat=myiostat)
       call checkiostat(myiostat,"opening fluxtsumfile")
       write(454,'(F18.12, T22, 400E20.8)',iostat=myiostat) curtime*dt, &
            gtausum(:), gtaunow(:)
       call checkiostat(myiostat,"writing fluxtsumfile")
       close(454)
    endif

    call mpibarrier()

    deallocate(imke,impe,imV2,imyderiv,imkeop,impeop,imyop)
    deallocate(reke,repe,reV2,reyderiv,rekeop,repeop,reyop)
    deallocate(ketmo)

  end subroutine fluxgtau_during0


  subroutine fluxgtau0(alg,wwin,bbin)
!! actually compute the flux in a post processing kind of manner
!! input :
!! alg - determines how the memory management algorithm for loading up previous wavefunctions
    use fluxgtaubiomod
    use biorthomod
    use parameters
    use walkmod
    use mpimod
    use mpisubmod
    use pulsesubmod
    use orbmultsubmod   !! gauge_transform
    use utilmod
    implicit none
!! FluxOpType:
!! 0       = use one-e potential and two-e contribution routines  (exact treatment)
!! 1       = replace one-e potential + two-e potential with halfnium one-e potential  (recommended)
!! 2       = use full one-e potential, no two-e 
!! other   = only KE

    type(walktype),target,intent(in) :: wwin,bbin
    type(walktype),pointer :: myww
    integer,intent(in) :: alg
    integer :: curtime,oldtime,k,nt,i,molength,alength,  BatchSize,NBat,brabat,brareadsize, &
         bratime,ketbat,ketreadsize,kettime,bratop,itime,jtime,times(1:20)=0, &
         imc, tau, ispf, myiostat
    real*8 :: MemTot,MemVal,dt, myfac,wfi,estep,windowfunct,MemNum,tentsum
    complex*16, allocatable :: FTgtau(:,:), pulseft(:,:)
    real*8, allocatable :: pulseftsq(:), tentfunction(:)
    DATATYPE :: gtaunow(mcscfnum), ftgtausum(mcscfnum), gtausum(mcscfnum), &
         pots1(3)=0d0, nullvector(numr), gtausave(mcscfnum), csum
    DATATYPE, allocatable :: bramo(:,:,:),braavec(:,:,:,:), ketmo(:,:,:),ketavec(:,:,:,:),&
         gtau(:,:), mobio(:,:),abio(:,:,:), &
         read_bramo(:,:),read_braavec(:,:,:,:),read_ketmo(:,:),read_ketavec(:,:,:,:),&
         imke(:,:),impe(:,:),imV2(:,:,:,:), imyderiv(:,:), imkeop(:,:),impeop(:,:),  imyop(:,:),&
         reke(:,:),repe(:,:),reV2(:,:,:,:), reyderiv(:,:), rekeop(:,:),repeop(:,:),  reyop(:,:),&
         transxmo(:,:)
    DATATYPE,target :: smo(nspf,nspf)

    if (ceground.eq.(0d0,0d0)) then
       OFLWR "Eground is ZERO.  Are you sure?  If want zero just make it small."
       WRFL  "   Otherwise need eground: initial state energy."; CFLST
    endif

    gtaunow=0; ftgtausum=0; gtausum=0; gtaunow=0; nullvector=0; gtausave=0

!! initial setup

!!  dt=real(FluxInterval*FluxSkipMult,8)*par_timestep;  nt=floor(final time/dt)

    dt=real(FluxInterval*FluxSkipMult,8)*par_timestep;  nt=floor(real(numpropsteps,8)/fluxinterval/fluxskipmult)

    if (FluxOpType.eq.0) then    !! exact expression with two-electron
       myww=>wwin
    else
       myww=>bbin
    endif

    call mpibarrier()
    OFLWR "Initial allocation electronflux"; CFL
    
    allocate(gtau(0:nt,mcscfnum))
    gtau(:,:)=0d0

    allocate(imke(nspf,nspf),impe(nspf,nspf),imV2(nspf,nspf,nspf,nspf),imyderiv(nspf,nspf),&
         reke(nspf,nspf),repe(nspf,nspf),reV2(nspf,nspf,nspf,nspf),reyderiv(nspf,nspf),&
         rekeop(spfsize,nspf),repeop(spfsize,nspf),reyop(spfsize,nspf),&
         mobio(spfsize,nspf),abio(numr,myww%firstconfig:myww%lastconfig,mcscfnum), &
         imkeop(spfsize,nspf),impeop(spfsize,nspf),imyop(spfsize,nspf))

    imke=0; impe=0; imV2=0; imyderiv=0; imkeop=0; impeop=0; imyop=0;
    reke=0; repe=0; reV2=0; reyderiv=0; rekeop=0; repeop=0; reyop=0
    mobio=0;
    if (myww%lastconfig.ge.myww%firstconfig) then
       abio=0
    endif

!! determine if we should do batching or not
!! 250,000 words/MB, real*8 2words/#, complex*16 4words/#

#ifdef REALGO
    MemVal = 1.25d5
#else
    MemVal = 6.25d4
#endif

!! memory based on rank 1 which right now stores the entire a-vector, which is wasteful.
!! Unsatisfactory: need better MPI I/O.
!! For orbitals, read_bramo and read_ketmo are not dimensioned with BatchSize any more 
!! --> MemNum based on numconfig (whole vector) and spfsize (per processor, bramo and ketmo)

!! factors of 2 for bra and ket

    MemNum=0d0

!!$    if (parorbsplit.eq.3) then
!!$       MemNum = MemNum + spfsize*nspf*nprocs*2
!!$    else

    MemNum = MemNum + spfsize*nspf*2

!!$    endif

    if (par_consplit.ne.0) then
       MemNum = MemNum + myww%numconfig*numr*mcscfnum*2
    else
       MemNum = MemNum + myww%maxconfigsperproc*numr*mcscfnum*2
    endif

    call openfile()
    write(mpifileptr,'(A30,F9.3,A3)') " Guess at necessary memory is ",&
         (nt+1)*MemNum/MemVal," MB"
    if(alg.eq.0) then
       write(mpifileptr,*) "g(tau) will be computed with all of psi in core"
       BatchSize=nt+1
    else
       MemTot=real(alg,8)    
       write(mpifileptr,*) "g(tau) will be computed with all psi being read in batches"
       write(mpifileptr,'(A33,F9.3,A3)') " Desired amount of memory to use ",MemTot," MB"
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
    OFLWR "Allocating psi arrays for electronflux"; CFL

    allocate(ketmo(spfsize,nspf,BatchSize),ketavec(numr,myww%firstconfig:myww%lastconfig,mcscfnum,BatchSize),&
         bramo(spfsize,nspf,BatchSize),braavec(numr,myww%firstconfig:myww%lastconfig,mcscfnum,BatchSize))
    ketmo=0; bramo=0
    if (myww%lastconfig.ge.myww%firstconfig) then
       ketavec=0; braavec=0
    endif
    if (myrank.eq.1) then
       if (parorbsplit.eq.3) then
          allocate(read_bramo(spfsize*nprocs,nspf),read_ketmo(spfsize*nprocs,nspf))
       else
          allocate(read_bramo(spfsize,nspf),read_ketmo(spfsize,nspf))
       endif
    else
       allocate(read_bramo(1,nspf),read_ketmo(1,nspf))
    endif
    read_bramo=0; read_ketmo=0

    if (myrank.eq.1) then
       allocate(read_braavec(numr,myww%numconfig,mcscfnum,BatchSize),&
            read_ketavec(numr,myww%numconfig,mcscfnum,BatchSize))
    else
       allocate(read_braavec(1,1,mcscfnum,BatchSize),read_ketavec(1,1,mcscfnum,BatchSize))
    endif
    read_braavec=0; read_ketavec=0

    NBat=ceiling(real(nt+1)/real(BatchSize))
    ketreadsize=0;  brareadsize=0

    if (myrank.eq.1) then
       inquire (iolength=molength) read_ketmo(:,:)
       inquire (iolength=alength) read_ketavec(:,:,:,1)
    endif
    call mympiibcastone(molength,1); call mympiibcastone(alength,1)

    call openfile()
    write(mpifileptr,*) "MO record length is ",molength
    write(mpifileptr,*) "AVEC record length is ",alength
    call closefile()

    if (myrank.eq.1) then
       open(454, file=fluxtsumfile, status="unknown",iostat=myiostat)
       call checkiostat(myiostat,"opening fluxtsumfile")
       write(454,*,iostat=myiostat) "#KVL flux sum: time, flux integral, flux"
       call checkiostat(myiostat,"writing fluxtsumfile")
       write(454,*)
       close(454)
    endif

    gtausum(:)=0d0

!! begin the ket batch read loop
    do ketbat=1,NBat
       call myclock(itime)
       OFLWR "Reading ket batch ", ketbat, " of ", NBat; CFL
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

!! gaugefluxflag available to attempt gauge-invariant calculation with photoionization
!! during the pulse.  Best performance seems to be calculate wfn in length gauge, transform
!! to velocity gauge for flux: velflag=0, gaugefluxflag=2 (Volkov phase is only a function of 
!! time in the velocity gauge, no preferred origin in the velocity gauge)

!! gaugefluxflag=1, transform velocity to length; gaugefluxflag=2, transform length to velocity
       if ((gaugefluxflag.eq.1.and.velflag.ne.0).or.(gaugefluxflag.eq.2.and.velflag.eq.0)) then
          allocate(transxmo(spfsize,nspf))
          transxmo=0d0
          do i=1,ketreadsize
             curtime=(ketbat-1)*BatchSize+i-1 
             call gauge_transform(velflag,curtime*dt,nspf,ketmo(:,:,i),transxmo(:,:))
             ketmo(:,:,i) = transxmo(:,:)
          enddo
          deallocate(transxmo)
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
          call mympibcast(ketavec(:,:,:,1:ketreadsize),1,numr*myww%numconfig*mcscfnum*ketreadsize)
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

       call myclock(jtime)
       times(1)=times(1)+jtime-itime

!! begin the bra batch read loop
       do brabat=1,ketbat
          call myclock(itime)
          OFLWR "Reading bra batch ", brabat, " of ", ketbat; CFL

          brareadsize=min(BatchSize,nt+1-(brabat-1)*BatchSize)
          if(brabat.eq.ketbat) then
             bramo=ketmo;        
             if (tot_adim.gt.0) then
                braavec=ketavec
             endif
          else 
!! ORBITALS
             if(myrank.eq.1) then
                open(1001,file=fluxmofile,status="old",form="unformatted",&
                     access="direct",recl=molength,iostat=myiostat)
                call checkiostat(myiostat,"opening "//fluxmofile)
             endif
             do i=1,brareadsize
                if(myrank.eq.1) then
                   k=FluxSkipMult*((brabat-1)*BatchSize+i-1)+1
                   read(1001,rec=k,iostat=myiostat) read_bramo(:,:)
                   call checkiostat(myiostat,"reading "//fluxmofile)
                endif
                if (parorbsplit.ne.3) then
                   if (myrank.eq.1) then
                      bramo(:,:,i)=read_bramo(:,:)
                   endif
                   call mympibcast(bramo(:,:,i),1,totspfdim)
                else
                   do ispf=1,nspf
                      call splitscatterv(read_bramo(:,ispf),bramo(:,ispf,i))
                   enddo
                endif
             enddo               !! do i=1,brareadsize
             if (myrank.eq.1) then
                close(1001)
             endif

!! gaugefluxflag available to attempt gauge-invariant calculation with photoionization
!! during the pulse.  Best performance seems to be calculate wfn in length gauge, transform
!! to velocity gauge for flux: velflag=0, gaugefluxflag=2 (Volkov phase is only a function of 
!! time in the velocity gauge, no preferred origin in the velocity gauge)

!! gaugefluxflag=1, transform velocity to length; gaugefluxflag=2, transform length to velocity
             if ((gaugefluxflag.eq.1.and.velflag.ne.0).or.(gaugefluxflag.eq.2.and.velflag.eq.0)) then
                allocate(transxmo(spfsize,nspf))
                transxmo=0d0
                do i=1,brareadsize
                   oldtime=(brabat-1)*BatchSize+i-1
                   call gauge_transform(velflag,oldtime*dt,nspf,bramo(:,:,i),transxmo(:,:))
                   bramo(:,:,i) = transxmo(:,:)
                enddo
                deallocate(transxmo)
             endif

!! A-VECTOR
             if(myrank.eq.1) then
                open(1002,file=fluxafile,status="old",form="unformatted",&
                     access="direct",recl=alength,iostat=myiostat)
                call checkiostat(myiostat,"opening "//fluxafile)
                do i=1,brareadsize
                   k=FluxSkipMult*((brabat-1)*BatchSize+i-1)+1
                   read(1002,rec=k,iostat=myiostat) read_braavec(:,:,:,i) 
                enddo
                call checkiostat(myiostat,"reading "//fluxafile)
                close(1002)
             endif

             if (par_consplit.eq.0) then
                if (myrank.eq.1) then
                   braavec(:,:,:,1:brareadsize)=read_braavec(:,:,:,1:brareadsize)
                endif
                call mympibcast(braavec(:,:,:,1:brareadsize),1,&
                     numr*myww%numconfig*mcscfnum*brareadsize)
             else
                do i=1,brareadsize
                   do imc=1,mcscfnum
                      if (tot_adim.gt.0) then
                         call myscatterv(read_braavec(:,:,imc,i),&
                              braavec(:,:,imc,i),configs_perproc(:)*numr)
                      else
                         call myscatterv(read_braavec(:,:,imc,i),&
                              nullvector(:),configs_perproc(:)*numr)
                      endif
                   enddo
                enddo
             endif
          endif  !! if brabat.eq.ketbat
          call myclock(jtime)
          times(1)=times(1)+jtime-itime
        
!! loop over all time for the ket of the flux integral
          do kettime=1,ketreadsize

!! get the one-e half transformed matrix elements for this ket time
             call myclock(itime)
             curtime=(ketbat-1)*BatchSize+kettime-1 

             imkeop=0; impeop=0;  imyop=0
             if (nucfluxopt.ne.2) then
                call flux_op_onee(ketmo(:,:,kettime),imkeop,impeop,1)
                if (nonuc_checkflag.eq.0) then
                   call flux_op_nuc(ketmo(:,:,kettime),imyop,1)
                endif
             endif

             rekeop=0; repeop=0;
             if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0) then
                call flux_op_onee(ketmo(:,:,kettime),rekeop,repeop,2)
                call flux_op_nuc(ketmo(:,:,kettime),reyop,2)
             endif
             call myclock(jtime);        times(2)=times(2)+jtime-itime

!! determine bounds of loop over bras and setup doing in parallel with mpi!
             if(brabat.lt.ketbat) then
                bratop=brareadsize
             else if (brabat.gt.ketbat) then
                OFLWR "WHATBAT?",brabat,ketbat; CFLST
                bratop=-99
             else
                bratop=kettime
             endif
           
!! loop over all previous time for the bra of the flux integral
             do bratime=1,bratop

                oldtime=(brabat-1)*BatchSize+bratime-1
              
!! biortho this pair of times!        
                call myclock(itime)
                if (tot_adim.gt.0) then
                   abio(:,:,:)=braavec(:,:,:,bratime)
                endif
                call bioset(fluxgtaubiovar,smo,numr,bbin)
#ifdef CNORMFLAG
                fluxgtaubiovar%hermonly=.true.
#endif
                if (tot_adim.gt.0) then
                   call biortho(bramo(:,:,bratime),ketmo(:,:,kettime),mobio,abio(:,:,1),fluxgtaubiovar)
                   do imc=2,mcscfnum
                      call biotransform(bramo(:,:,bratime),ketmo(:,:,kettime),abio(:,:,imc),fluxgtaubiovar)
                   enddo
#ifdef CNORMFLAG
                   abio(:,:,:)=conjg(abio(:,:,:))
#endif
                else
                   call biortho(bramo(:,:,bratime),ketmo(:,:,kettime),mobio,nullvector(:),fluxgtaubiovar)
                   do imc=2,mcscfnum
                      call biotransform(bramo(:,:,bratime),ketmo(:,:,kettime),nullvector(:),fluxgtaubiovar)
                   enddo
                endif
                fluxgtaubiovar%hermonly=.false.   !! in case fluxgtaubiovar is reused

                call myclock(jtime);          times(3)=times(3)+jtime-itime;  itime=jtime

!! complete the one-e potential and kinetic energy matrix elements           

                imke=0;  impe=0;  imyderiv=0;  imke=0;  imke=0; 

                if (nucfluxopt.ne.2) then
                   call MYGEMM('C','N',nspf,nspf,spfsize,DATAONE,mobio,&
                        spfsize,imkeop,spfsize,DATAZERO,imke,nspf)
                   call MYGEMM('C','N',nspf,nspf,spfsize,DATAONE,mobio,&
                        spfsize,impeop,spfsize,DATAZERO,impe,nspf)
                   if (nonuc_checkflag.eq.0) then
                      call MYGEMM('C','N',nspf,nspf,spfsize,DATAONE,mobio,&
                           spfsize,imyop,spfsize,DATAZERO,imyderiv,nspf)
                   endif
                endif

                reke=0; repe=0; reyderiv=0

                if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0) then
                   call MYGEMM('C','N',nspf,nspf,spfsize,DATAONE,mobio,&
                        spfsize,rekeop,spfsize,DATAZERO,reke,nspf)
                   call MYGEMM('C','N',nspf,nspf,spfsize,DATAONE,mobio,&
                        spfsize,repeop,spfsize,DATAZERO,repe,nspf)
                   call MYGEMM('C','N',nspf,nspf,spfsize,DATAONE,mobio,&
                        spfsize,reyop,spfsize,DATAZERO,reyderiv,nspf)
                endif

                call myclock(jtime);          times(4)=times(4)+jtime-itime;  itime=jtime

                if (parorbsplit.eq.3) then
                   if (nucfluxopt.ne.2) then
                      call mympireduce(imke,nspf**2)
                      call mympireduce(impe,nspf**2)
                      if (nonuc_checkflag.eq.0) then
                         call mympireduce(imyderiv,nspf**2)
                      endif
                   endif
                   if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0)then
                      call mympireduce(reke,nspf**2)
                      call mympireduce(repe,nspf**2)
                      call mympireduce(reyderiv,nspf**2)
                   endif
                endif

                call myclock(jtime);          times(5)=times(5)+jtime-itime;  itime=jtime

!! get the two-e contribution for exact formula (fluxoptype=0)

                imV2=0d0; reV2=0
                if(FluxOpType.eq.0) then
                   if (nucfluxopt.ne.2) then
                      call flux_op_twoe(mobio,ketmo(:,:,kettime),imV2,1)
                   endif
                   if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0)then
                      call flux_op_twoe(mobio,ketmo(:,:,kettime),reV2,2)
                   endif
                endif
                call myclock(jtime);         times(6)=times(6)+jtime-itime;  itime=jtime

!! evaluate the actual g(tau) expression

                do imc=1,mcscfnum
                   if (tot_adim.gt.0) then
                      gtaunow(imc) = fluxeval(abio(:,:,imc),ketavec(:,:,imc,kettime),&
                           imke,impe,imV2,imyderiv,reke,repe,reV2,reyderiv,myww)
                   else
                      gtaunow(imc) = fluxeval(nullvector(:),nullvector(:),&
                           imke,impe,imV2,imyderiv,reke,repe,reV2,reyderiv,myww)
                   endif
                enddo

                oldtime=(brabat-1)*BatchSize+bratime-1
                tau=curtime-oldtime; 

                if (curtime.eq.0.and.oldtime.eq.0) then
                   gtausave(:) = gtaunow(:)
                endif

!! dt factor here
                gtau(tau,:) = gtau(tau,:) + gtaunow(:) * dt

                call myclock(jtime);          times(7)=times(7)+jtime-itime

                if (tau.eq.0) then

!! BRABAT = KETBAT : diagonal part, integral dt

!! nevermind KVLsum.dat and gtau.dat
!! now integrals dt and domega with fluxtsumfile in namelist &parinp
!! and additional column in spifile

                   call myclock(itime)

                   if (flux_subtract.ne.0) then
                      gtaunow(:) = gtaunow(:) - gtausave(:)
                   endif

                   gtausum(:) = gtausum(:) + gtaunow(:) * dt / 2d0   !! 2 looks correct

                   if (myrank.eq.1) then
                      call myclock(itime)
                      open(454, file=fluxtsumfile, status="old",  position="append",iostat=myiostat)
                      call checkiostat(myiostat,"opening fluxtsumfile")
                      write(454,'(F18.12, T22, 400E20.8)',iostat=myiostat) curtime*dt, &
                           gtausum(:), gtaunow(:)
                      call checkiostat(myiostat,"writing fluxtsumfile")
                      call myclock(jtime);        times(8)=times(8)+jtime-itime
                      close(454)
                   endif

                endif

             enddo !! do bratime

          enddo !! do kettime

          if (brabat.eq.ketbat) then
             if (notiming.ne.2) then
                write(mpifileptr,'(A28,F10.4)') " Timing statistics as of T= ",curtime*dt
                write(mpifileptr,'(100A10)') "Times: ", "Read", "One-e", "Biorth", &
                     "Matel", "moreMPI", "Two-e", "eval", "write"
                write(mpifileptr,'(A10,100I10)') " ", times(1:8)/1000; CFL
             endif
          endif
       enddo !! do brabat
    enddo !! do ketbat

    if (curtime.ne.nt) then
       OFLWR "DOOG CURTIME NT ERR",curtime,nt; CFLST
    endif

!! create the xsec as it should be now

    call mpibarrier()
    OFLWR " Taking the FT of g(tau) to get xsection at T= ",curtime*dt; CFL

    allocate(ftgtau(-curtime:curtime,mcscfnum), pulseft(-curtime:curtime,3), &
         pulseftsq(-curtime:curtime), tentfunction(-curtime:curtime))
    ftgtau(:,:)=0d0; pulseft(:,:)=0d0; pulseftsq(:)=0d0; tentfunction(:)=0d0

    do i=-curtime,curtime
       tentfunction(i) = (1+curtime-abs(i)) * windowfunct(abs(i),curtime,16)
    enddo

    tentsum = get_rtentsum(curtime,tentfunction(-curtime:curtime))

    do i=0,curtime
       ftgtau(i,:) = ALLCON(gtau(i,:))   * windowfunct(i,curtime,16) * &  !! action 16
            exp((0d0,-1d0)*ALLCON(ceground)*dt*i)
    enddo

    do i=1,curtime
       ftgtau(-i,:) = ALLCON(ftgtau(i,:))
    enddo

!! subtract tent function   (1+curtime-abs(i))/(curtime+1)^2   for better performance

    if (flux_subtract.ne.0) then
       do imc=1,mcscfnum
          csum = get_ctentsum(curtime,ftgtau(-curtime:curtime,imc))
          ftgtau(-curtime:curtime,imc) = ftgtau(-curtime:curtime,imc) - &
               csum/tentsum * tentfunction(-curtime:curtime)
       enddo
    endif

    do i=1,curtime
       call vectdpot(i*dt,0,pots1,-1)   !! LENGTH GAUGE.
       if (pulsewindowtoo == 0) then
          pulseft(i,:)=pots1(:)
       else
          pulseft(i,:)=pots1(:) * windowfunct(i,curtime,16)
       endif
    enddo

    OFLWR "   ....Go ft...."; CFL
    do imc=1,mcscfnum
       call zfftf_wrap_diff(2*curtime+1,ftgtau(-curtime:curtime,imc),ftdiff)
    enddo
    OFLWR "   ....Go ft pulse...."; CFL
    do i=1,3
       call zfftf_wrap(2*curtime+1,pulseft(-curtime:curtime,i))
    enddo
    OFLWR "   ....Done with ft...."; CFL

    ftgtau(:,:)=ftgtau(:,:)*dt
    pulseft(:,:)=pulseft(:,:)*dt

    do i=-curtime,curtime
       ftgtau(i,:)=ftgtau(i,:)*exp((0.d0,1.d0)*(curtime+i)*curtime*2*pi/real(2*curtime+1))
       pulseft(i,:)=pulseft(i,:)*exp((0.d0,1.d0)*(curtime+i)*curtime*2*pi/real(2*curtime+1))
    enddo

    pulseftsq(:) = abs(pulseft(:,1)**2) + abs(pulseft(:,2)**2) + abs(pulseft(:,3)**2)

    estep=2*pi/dt/(2*curtime+1)

    if (myrank.eq.1) then
       open(1004,file=spifile,status="replace",action="readwrite",&
            position="rewind",iostat=myiostat)
       call checkiostat(myiostat,"opening "//spifile)
       write(1004,*,iostat=myiostat)
       call checkiostat(myiostat,"writing "//spifile)
       write(1004,*) "# eight columns."

!! nevermind KVLsum.dat and gtau.dat
!! now integrals dt and domega with fluxtsumfile in namelist &parinp
!! and additional column in spifile

       write(1004,*) "# Omega (column 1); |pulse ft|^2 (2); cross section (Mb) (3); flux (5); flux integral domega (7)"
       write(1004,*)
       ftgtausum(:)=0d0
       do i=-curtime,curtime
          wfi=(i+curtime)*Estep

          ftgtausum(:)=ftgtausum(:)+ftgtau(i,:) * estep / PI / 4d0   !! 4 looks correct

!! LENGTH GAUGE WAS FT'ed multiply by wfi dont divide
!! NEVERMIND FACTOR OF 1/3
!!        myfac = 5.291772108d0**2 / 3d0 * 2d0 * PI / 1.37036d2 * wfi

!! WITH THIS FACTOR, NOW THE QUANTUM MECHANICAL CROSS SECTION IN 
!! MEGABARNS (10^-18 cm^2) IS IN COLUMN 3 REAL PART
          myfac = 5.291772108d0**2 * 2d0 * PI / 1.37036d2 * wfi

          write(1004,'(F8.4,100E18.6)',iostat=myiostat) wfi, pulseftsq(i), &
               FTgtau(i,:)/pulseftsq(i) * myfac, ftgtau(i,:) / 4 / PI, ftgtausum(:)
       enddo
       call checkiostat(myiostat,"writing "//spifile)
       close(1004)
    endif


    call mpibarrier()
    OFLWR "  finished electronflux, stopping."; CFL
    call waitawhile()
    call mpistop()

    deallocate(ftgtau,pulseft,pulseftsq,tentfunction)
    deallocate(bramo,ketmo,braavec,ketavec)
    deallocate(read_bramo,read_ketmo,read_braavec,read_ketavec)
    deallocate(gtau)
    deallocate(imke,impe,imV2,imyderiv,imkeop,impeop,imyop)
    deallocate(reke,repe,reV2,reyderiv,rekeop,repeop,reyop)
    deallocate(mobio,abio)

  end subroutine fluxgtau0

!! begin the flux matrix element and contraction routine section
!! flag=1, imag part; flag=2, real part; flag=0, all

  subroutine flux_op_onee(inspfs,keop,peop,flag)
    use parameters
    use orbgathersubmod
    use orbmultsubmod
    use fluxutilmod
    implicit none
    integer,intent(in) :: flag
    DATATYPE, intent(in) :: inspfs(spfsize,nspf)
    DATATYPE,intent(out) ::  keop(spfsize,nspf),peop(spfsize,nspf)
    integer :: lowspf,highspf,numspf

    lowspf=1; highspf=nspf
    if (parorbsplit.eq.1) then
       call getOrbSetRange(lowspf,highspf)
    endif
    numspf=highspf-lowspf+1

!! initialize
    keop=0d0; peop=0d0

!! the kinetic energy

    if (numspf.gt.0) then
       select case(flag)
       case(1)
          call mult_imke(numspf,inspfs(:,lowspf:highspf),keop(:,lowspf:highspf))
       case(2)
          call mult_reke(numspf,inspfs(:,lowspf:highspf),keop(:,lowspf:highspf))
       case default
          OFLWR "not supppported"; CFLST
          call mult_ke(inspfs(:,lowspf:highspf),keop(:,lowspf:highspf),numspf,"booga",2)
       end select

!! the one-e potential energy 

       select case(flag)
       case(1)
          if(FluxOpType.eq.0.or.FluxOpType.eq.2) then
             call mult_impot(numspf,inspfs(:,lowspf:highspf),peop(:,lowspf:highspf))
          else if(FluxOpType.eq.1) then
             call mult_imhalfniumpot(numspf,inspfs(:,lowspf:highspf),peop(:,lowspf:highspf))
          endif
       case(2)
          if(FluxOpType.eq.0.or.FluxOpType.eq.2) then
             call mult_repot(numspf,inspfs(:,lowspf:highspf),peop(:,lowspf:highspf))
          else if(FluxOpType.eq.1) then
             call mult_rehalfniumpot(numspf,inspfs(:,lowspf:highspf),peop(:,lowspf:highspf))
          endif
       case default
          OFLWR "Nottt supporteddd"; CFLST
          if(FluxOpType.eq.0.or.FluxOpType.eq.2) then
             call mult_pot(numspf,inspfs(:,lowspf:highspf),peop(:,lowspf:highspf))
          else if(FluxOpType.eq.1) then
             call mult_halfniumpot(numspf,inspfs(:,lowspf:highspf),peop(:,lowspf:highspf))
          endif
       end select

    endif

    if (parorbsplit.eq.1) then
       call mpiorbgather(peop,spfsize)
       call mpiorbgather(keop,spfsize)
    endif

  end subroutine flux_op_onee


!! operates with y derivative cross term.  
!! Assumes y^2 and l^2 are in the one electron hamiltonian.

!!  yop complex antisymmetric  hermitian part is imaginary   
!!  times antihermitian op in R (assuming no R ecs scaling, 
!!  as in all these routines)

!! fortran imag() returns real value so multiply by i

!! flag=1, flux (imag); flag=0, all    2 flux (imag)

  subroutine flux_op_nuc(inspfs,yop,flag)
    use parameters
    use orbgathersubmod
    use fluxutilmod
    implicit none
    DATATYPE, intent(in) :: inspfs(spfsize,nspf)
    DATATYPE,intent(out) ::  yop(spfsize,nspf)
    integer,intent(in) :: flag
    integer :: lowspf,highspf,numspf

    lowspf=1; highspf=nspf
    if (parorbsplit.eq.1) then
       call getOrbSetRange(lowspf,highspf)
    endif
    numspf=highspf-lowspf+1

!! the kinetic energy

    if (numspf.gt.0) then
       select case(flag)
       case(1)
          call op_imyderiv(numspf,inspfs(:,lowspf:highspf),yop(:,lowspf:highspf))
       case(2)
          call op_reyderiv(numspf,inspfs(:,lowspf:highspf),yop(:,lowspf:highspf))
       case default
          OFLWR "NNOT supported!"; CFLST
          call op_yderiv(numspf,inspfs(:,lowspf:highspf),yop(:,lowspf:highspf))
       end select

    endif

    if (parorbsplit.eq.1) then
       call mpiorbgather(yop,spfsize)
    endif

  end subroutine flux_op_nuc


  subroutine flux_op_twoe(mobra,moket,V2,flag) !! ok unused flux_op_twoe, unfinished
    use parameters
    use orbgathersubmod
    use mpisubmod
    implicit none
    DATATYPE,intent(in) :: mobra(spfsize,nspf),moket(spfsize,nspf)
    DATATYPE,intent(out) :: V2(nspf,nspf,nspf,nspf)
    integer,intent(in) :: flag
    DATATYPE, allocatable :: tempreduced(:,:,:)
    integer :: lowspf,highspf,numspf

    lowspf=1; highspf=nspf
    if (parorbsplit.eq.1) then
       call getOrbSetRange(lowspf,highspf)
    endif
    numspf=highspf-lowspf+1

    if (numspf.gt.0) then
       allocate(tempreduced(reducedpotsize,nspf,lowspf:highspf))
       tempreduced=0
       OFLWR "DOME TWOE_MATEL0000"; CFLST
!!$       call call_twoe_matel0000(lowspf,highspf,mobra,moket,V2(:,:,:,lowspf:highspf),tempreduced,"boogie",2,flag)
       deallocate(tempreduced)
    endif

    if (parorbsplit.eq.1) then
       call mpiorbgather(V2,nspf**3)
    endif
    if (parorbsplit.eq.3) then
       call mympireduce(V2,nspf**4)
    endif

  end subroutine flux_op_twoe


  function fluxeval(abra,aket,imke,impe,imV2,imyderiv,reke,repe,reV2,reyderiv,myww)
    use parameters 
    use walkmod
    use configptrmod
    use sparseptrmod
    use sparsemultmod
    use asssubmod
    use mpisubmod
    implicit none
    DATATYPE,intent(in) :: abra(numr,first_config:last_config),aket(numr,first_config:last_config),&
         imke(nspf,nspf),impe(nspf,nspf),imV2(nspf,nspf,nspf,nspf),imyderiv(nspf,nspf),&
         reke(nspf,nspf),repe(nspf,nspf),reV2(nspf,nspf,nspf,nspf),reyderiv(nspf,nspf)
    type(walktype),intent(in) :: myww
    DATATYPE, allocatable :: multket(:,:), ketwork(:,:), conjgket(:,:)
    DATATYPE :: fluxeval, outsum
    type(CONFIGPTR) :: matrix_ptr
    type(SPARSEPTR) :: sparse_ptr

    call configptralloc(matrix_ptr,myww)
    matrix_ptr%kefac = 0d0
    matrix_ptr%constfac = 0d0

    if (sparseopt.ne.0) then
       call sparseptralloc(sparse_ptr,myww)
    endif

    allocate(multket(numr,first_config:last_config), ketwork(numr,first_config:last_config),&
         conjgket(numr,first_config:last_config))
    if (last_config.ge.first_config) then
       multket=0; ketwork=0; conjgket=0
    endif

    outsum=0d0


!!!!   constant term - imaginary part of frozen, nuclear repulsion, energyshift  !!!!
!!!!        is NOT included since it does not define asymptotic region           !!!!
!!!!        constfac always set 0                                                !!!!

!!!!   imaginary part of electronic matrix elements with real part coefficients in r   !!!!

    if (nucfluxopt.ne.2) then

       matrix_ptr%kefac = 0d0
       matrix_ptr%constfac = 0d0

       matrix_ptr%xpotmatel(:,:) = impe(:,:)
       matrix_ptr%xopmatel(:,:)  = imke(:,:)
       matrix_ptr%xymatel(:,:)   = imyderiv(:,:)
       matrix_ptr%xtwoematel(:,:,:,:) = imV2(:,:,:,:)

! 06-16 umm shouldn't same logic be here as below?  nucfluxopt
!       if (sparseopt.ne.0) then
!          call assemble_sparsemats(myww,matrix_ptr,sparse_ptr,1,0,0,0)
!       endif
!
!       call sparseconfigmult(myww,aket,ketwork,matrix_ptr,sparse_ptr,1,0,0,0,0d0,-1)

       if (sparseopt.ne.0) then
          call assemble_sparsemats(myww,matrix_ptr,sparse_ptr,1,min(nucfluxopt,1),0,0)
       endif

       call sparseconfigmult(myww,aket,ketwork,matrix_ptr,sparse_ptr,1,min(nucfluxopt,1),0,0,0d0,-1)

       matrix_ptr%xpotmatel(:,:)      = ALLCON(matrix_ptr%xpotmatel(:,:))
       matrix_ptr%xopmatel(:,:)       = ALLCON(matrix_ptr%xopmatel(:,:))
       matrix_ptr%xymatel(:,:)        = ALLCON(matrix_ptr%xymatel(:,:))
       matrix_ptr%xtwoematel(:,:,:,:) = ALLCON(matrix_ptr%xtwoematel(:,:,:,:))

       if (sparseopt.ne.0) then
          call assemble_sparsemats(myww,matrix_ptr,sparse_ptr,1,min(nucfluxopt,1),0,0)
       endif

       if (tot_adim.gt.0) then
          conjgket(:,:) = ALLCON(aket(:,:))
       endif
       call sparseconfigmult(myww,conjgket,multket,matrix_ptr,sparse_ptr,1,min(nucfluxopt,1),0,0,0d0,-1)
       if (tot_adim.gt.0) then
          multket(:,:)=ALLCON(multket(:,:))
       endif
   
       if (tot_adim.gt.0.and.nucfluxopt.ne.2) then
          ketwork(:,:)=(ketwork(:,:)+multket(:,:)) / (-1d0)    !! -2x imag part overall
          outsum = outsum + hermdot(abra,ketwork,tot_adim)
       endif

    endif

!!!!   real part of electronic operators, imaginary part of nuclear operators   !!!!

    if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0) then

       matrix_ptr%kefac = 1d0
       matrix_ptr%constfac = 0d0

       matrix_ptr%xpotmatel(:,:) = repe(:,:)
       matrix_ptr%xopmatel(:,:)  = reke(:,:)
       matrix_ptr%xymatel(:,:)   = reyderiv(:,:)
       matrix_ptr%xtwoematel(:,:,:,:) = reV2(:,:,:,:)

       if (sparseopt.ne.0) then
          call assemble_sparsemats(myww,matrix_ptr,sparse_ptr,1,1,0,0)
       endif
       call sparseconfigmult(myww,aket,ketwork,matrix_ptr,sparse_ptr,1,1,0,0,0d0,-1)

       matrix_ptr%xpotmatel(:,:)      = ALLCON(matrix_ptr%xpotmatel(:,:))
       matrix_ptr%xopmatel(:,:)       = ALLCON(matrix_ptr%xopmatel(:,:))
       matrix_ptr%xymatel(:,:)        = ALLCON(matrix_ptr%xymatel(:,:))
       matrix_ptr%xtwoematel(:,:,:,:) = ALLCON(matrix_ptr%xtwoematel(:,:,:,:))

       if (sparseopt.ne.0) then
          call assemble_sparsemats(myww,matrix_ptr,sparse_ptr,1,1,0,0)
       endif

       if (tot_adim.gt.0) then
          conjgket(:,:) = ALLCON(aket(:,:))
       endif
       call sparseconfigmult(myww,conjgket,multket,matrix_ptr,sparse_ptr,1,1,0,0,0d0,-1)
       if (tot_adim.gt.0) then
          multket(:,:)=ALLCON(multket(:,:))
       endif

       if (tot_adim.gt.0) then
          ketwork(:,:)=(ketwork(:,:)-multket(:,:)) / (0d0,-1d0)   !! -2x imag part overall
          outsum = outsum + hermdot(abra,ketwork,tot_adim)
       endif

    endif

    if (par_consplit.ne.0) then
       call mympireduceone(outsum)
    endif
    fluxeval=outsum

    if (sparseopt.ne.0) then
       call sparseptrdealloc(sparse_ptr)
    endif
    call configptrdealloc(matrix_ptr)
    deallocate(multket,ketwork,conjgket)

  end function fluxeval

end module fluxgtau0mod


module fluxgtaumod
contains

subroutine fluxgtau(alg)
  use configmod
  use fluxgtau0mod
  implicit none
  integer,intent(in) :: alg

  call fluxgtau0(alg,www,bioww)

end subroutine fluxgtau


subroutine fluxgtau_during(ketmo,ketavec,dt)
  use parameters !! mcscfnum
  use configmod
  use fluxgtau0mod
  implicit none
  real*8,intent(in) :: dt
  DATATYPE,intent(in) :: ketavec(tot_adim,mcscfnum), ketmo(spfsize,nspf)
  call fluxgtau_during0(www,bioww,ketmo,ketavec,dt)
end subroutine fluxgtau_during
    

subroutine fluxwrite(curtime,in_xmo,in_xa)
!! Write out the wavefunction for the flux calculation
!! input :
!! curtime - the current time
!! xmo - this time's current orbitals
!! xa - this time's current A vector
  use parameters
  use mpimod
  use mpisubmod
  implicit none
  integer,intent(in) :: curtime
  DATATYPE,intent(in) :: in_xmo(spfsize,nspf),in_xa(numr,first_config:last_config,mcscfnum)
  DATATYPE,allocatable :: xmo(:,:), xa(:,:,:)
  DATATYPE :: nullvector(numr)
  integer :: molength,alength,ispf,ii,myiostat

  nullvector(:)=0

  if (myrank.eq.1) then
     if (parorbsplit.eq.3) then
        allocate(xmo(spfsize*nprocs,nspf))
     else
        allocate(xmo(spfsize,nspf))
     endif
  else
     allocate(xmo(1,nspf))
  endif
  xmo(:,:)=0d0

  if (parorbsplit.eq.3) then
     do ispf=1,nspf
        call splitgatherv(in_xmo(:,ispf),xmo(:,ispf),.false.)
     enddo
  elseif (myrank.eq.1) then
     xmo(:,:)=in_xmo(:,:)
  endif

  if (myrank.eq.1) then
     allocate(xa(numr,num_config,mcscfnum))
  else
     allocate(xa(1,1,mcscfnum))
  endif
  xa(:,:,:)=0d0

  if (par_consplit.ne.0) then
     do ii=1,mcscfnum
        if (tot_adim.gt.0) then
           call mygatherv(in_xa(:,:,ii),xa(:,:,ii), configs_perproc(:)*numr,.false.)
        else
           call mygatherv(nullvector(:),xa(:,:,ii), configs_perproc(:)*numr,.false.)
        endif
     enddo
  elseif(myrank.eq.1) then
     xa(:,:,:)=in_xa(:,:,:)
  endif

  if (myrank.eq.1) then
     inquire (iolength=molength) xmo
     inquire (iolength=alength) xa
     call openfile
     write(mpifileptr,'(A27,F11.4)') " Saving wavefunction at T= ",&
          curtime*FluxInterval*par_timestep
     call closefile


     open(1001,file=fluxmofile,status="unknown",form="unformatted",&
          access="direct",recl=molength,iostat=myiostat)
     call checkiostat(myiostat,"opening "//fluxmofile)
     write(1001,rec=curtime+1,iostat=myiostat) xmo
     call checkiostat(myiostat,"writing "//fluxmofile)
     close(1001)
     open(1002,file=fluxafile,status="unknown",form="unformatted",&
          access="direct",recl=alength,iostat=myiostat)
     call checkiostat(myiostat,"opening "//fluxafile)
     write(1002,rec=curtime+1,iostat=myiostat) xa
     call checkiostat(myiostat,"writing "//fluxafile)
     close(1002)
  endif

  deallocate(xmo,xa)
  
  call mpibarrier()

end subroutine fluxwrite

end module
