
!! KVL ELECTRON FLUX ROUTINES

#include "Definitions.INC"

subroutine fluxwrite(curtime,xmo,xa)
!! Write out the wavefunction for the flux calculation
!! input :
!! curtime - the current time
!! xmo - this time's current orbitals
!! xa - this time's current A vector
  use parameters
  implicit none
  integer :: curtime,molength,alength
  DATATYPE :: xmo(spfsize,nspf),xa(numconfig,numr,mcscfnum)
  inquire (iolength=molength) xmo
  inquire (iolength=alength) xa
  call openfile
  write(mpifileptr,'(A27,F11.4)') " Saving wavefunction at T= ",curtime*FluxInterval*par_timestep
  call closefile
  open(1001,file=fluxmofile,status="unknown",form="unformatted",access="direct",recl=molength)
  open(1002,file=fluxafile,status="unknown",form="unformatted",access="direct",recl=alength)

  write(1001,rec=curtime+1) xmo;  write(1002,rec=curtime+1) xa
  close(1001);  close(1002)
end subroutine fluxwrite


module fluxgtaubiomod
  use biorthotypemod
  implicit none
  type(biorthotype),target :: fluxgtaubiovar
end module fluxgtaubiomod


subroutine fluxgtau(alg)
!! actually compute the flux in a post processing kind of manner
!! input :
!! alg - determines how the memory management algorithm for loading up previous wavefunctions
  use fluxgtaubiomod
  use biorthomod
  use parameters
  use walkmod
  use mpimod
  implicit none
!! FluxOpType:
!! 0       = use one-e potential and two-e contribution routines  (exact treatment)
!! 1       = replace one-e potential + two-e potential with halfnium one-e potential  (recommended)
!! 2       = use full one-e potential, no two-e 
!! other   = only KE

  integer :: alg,curtime,oldtime,tau,k,nt,i,molength,alength,  BatchSize,NBat,brabat,brareadsize, &
       bratime,ketbat,ketreadsize,kettime,bratop, whichbin,atime,btime,itime,jtime,times(1:7)=0, clow, &
       chigh,jproc,cnum, flag, ipulse, imc
  integer, allocatable :: bintimes(:)
  real*8 :: MemTot,MemVal,dt,wfi,ffac, piover2, estep, myft,xyfac,zfac
!!092712 real*8 :: Energy
  complex*16 :: Energy,  tdpotft
  real*8, allocatable :: xsec(:)
  complex*16, allocatable :: FTgtau(:,:), FTgtausave(:,:)
  DATATYPE :: dot,fluxeval,fluxevalval(mcscfnum), pots1(3)=0d0, pots2(3)=0d0
  logical :: doFT
  DATATYPE, allocatable :: gtau(:,:,:),wsave(:), tempgtau(:,:), mobio(:,:),abio(:,:,:)
  DATATYPE, allocatable, target :: bramo(:,:,:),braavec(:,:,:,:), ketmo(:,:,:),ketavec(:,:,:,:)
  DATATYPE, pointer :: moket(:,:),mobra(:,:),aket(:,:,:)
  DATATYPE, allocatable :: ke(:,:),pe(:,:),V2(:,:,:,:), yderiv(:,:),zdip(:,:), xydip(:,:)
  DATATYPE, allocatable :: keop(:,:),peop(:,:),  yop(:,:), zdipop(:,:), xydipop(:,:)
  DATATYPE, allocatable :: reke(:,:),repe(:,:),reV2(:,:,:,:), reyderiv(:,:),rezdip(:,:), rexydip(:,:)
  DATATYPE, allocatable :: rekeop(:,:),repeop(:,:),  reyop(:,:), rezdipop(:,:), rexydipop(:,:)
  DATATYPE,target :: smo(nspf,nspf)

  if (ceground.eq.(0d0,0d0)) then
     OFLWR "Eground is ZERO.  Are you sure?  If want zero just make it small. \n     Otherwise need eground: initial state energy."; CFLST
  endif

  piover2=atan2(1d0,1d0)*2

!! initial setup

!!  dt=real(FluxInterval*FluxSkipMult,8)*par_timestep;  nt=floor(final time/dt)

  dt=real(FluxInterval*FluxSkipMult,8)*par_timestep;  nt=floor(real(numpropsteps,8)/fluxinterval/fluxskipmult)


  allocate(gtau(0:nt,FluxNBins,mcscfnum),bintimes(FluxNBins),xsec(nEFlux))
  allocate(xydip(nspf,nspf),zdip(nspf,nspf),ke(nspf,nspf),pe(nspf,nspf),V2(nspf,nspf,nspf,nspf),yderiv(nspf,nspf))
  allocate(rexydip(nspf,nspf),rezdip(nspf,nspf),reke(nspf,nspf),repe(nspf,nspf),reV2(nspf,nspf,nspf,nspf),reyderiv(nspf,nspf))
  gtau(:,:,:)=0d0
  call getFTpulse(xsec) 
  allocate(mobio(spfsize,nspf),abio(numconfig,numr,mcscfnum), &
       keop(spfsize,nspf),peop(spfsize,nspf),xydipop(spfsize,nspf),zdipop(spfsize,nspf), &
       rekeop(spfsize,nspf),repeop(spfsize,nspf),rexydipop(spfsize,nspf),rezdipop(spfsize,nspf))

  allocate(yop(spfsize,nspf));  allocate(reyop(spfsize,nspf))
  
  rexydip=0d0;rezdip=0d0;reke=0d0;repe=0d0;rev2=0d0;reyderiv=0d0;rekeop=0d0;repeop=0d0;rexydipop=0d0;rezdipop=0d0
  xydip=0d0;zdip=0d0;ke=0d0;pe=0d0;v2=0d0;yderiv=0d0;keop=0d0;peop=0d0;xydipop=0d0;zdipop=0d0


!! determine if we should do batching or not
!! 250,000 words/MB, real*8 2words/#, complex*16 4words/#
#ifdef REALGO
  MemVal = 1.25d5
#else
  MemVal = 6.25d4
#endif
  call openfile()
  write(mpifileptr,'(A30,F9.3,A3)') " Guess at necessary memory is ",2d0*real((nt+1)*(numconfig*numr*mcscfnum+spfsize*nspf),8)/MemVal," MB"
  if(alg.eq.0) then
    write(mpifileptr,*) "g(tau) will be computed with all of psi in core"
    BatchSize=nt+1
  else
    MemTot=real(alg,8)    
    write(mpifileptr,*) "g(tau) will be computed with all psi being read in batches"
    write(mpifileptr,'(A33,F9.3,A3)') " Desired amount of memory to use ",MemTot," MB"
    BatchSize=floor(MemTot * MemVal / (2d0*real(numconfig*numr*mcscfnum+spfsize*nspf,8)))
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


  allocate(ketmo(spfsize,nspf,BatchSize),ketavec(numconfig,numr,mcscfnum,BatchSize))
  allocate(bramo(spfsize,nspf,BatchSize),braavec(numconfig,numr,mcscfnum,BatchSize))


  NBat=ceiling(real(nt+1)/real(BatchSize))
  ketreadsize=0;  brareadsize=0

  inquire (iolength=molength) ketmo(:,:,1);  inquire (iolength=alength) ketavec(:,:,:,1)

  write(mpifileptr,*) "MO record length is ",molength
  write(mpifileptr,*) "AVEC record length is ",alength
  call closefile()

  open(454, file="Dat/KVLsum.dat", status="unknown")
  write(454,*) "#KVL flux sum: itime, time, flux sum"
  write(454,*);  close(454)
  
!! figure out where bins begin
  tau=floor(real(nt,8)/real(FluxNBins,8))
  k=0
  do i=FluxNBins,2,-1
     k=k+tau
     bintimes(i)=k
  enddo
  bintimes(1)=nt

!! begin the ket batch read loop
  do ketbat=1,NBat
     call system_clock(atime)
     OFLWR "Reading ket batch ", ketbat, " of ", NBat; CFL
     ketreadsize=min(BatchSize,nt+1-(ketbat-1)*BatchSize)
     if(myrank.eq.1) then
        open(1001,file=fluxmofile,status="old",form="unformatted",access="direct",recl=molength)
        open(1002,file=fluxafile,status="old",form="unformatted",access="direct",recl=alength)
        do i=1,ketreadsize
           k=FluxSkipMult*((ketbat-1)*BatchSize+i-1)+1
           read(1001,rec=k) ketmo(:,:,i) ;        read(1002,rec=k) ketavec(:,:,:,i) 
        enddo
        close(1001);      close(1002)
     endif

     call mympibcast(ketmo(:,:,1:ketreadsize),1,totspfdim*ketreadsize)
     call mympibcast(ketavec(:,:,:,1:ketreadsize),1,numr*numconfig*mcscfnum*ketreadsize)

     call system_clock(btime)
     times(1)=times(1)+btime-atime;    times(2)=times(2)+btime-atime
     
!! begin the bra batch read loop
     do brabat=1,ketbat
        call system_clock(atime)
        OFLWR "Reading bra batch ", brabat, " of ", ketbat; CFL
        brareadsize=min(BatchSize,nt+1-(brabat-1)*BatchSize)
        if(brabat.eq.ketbat) then
           bramo=ketmo;        braavec=ketavec
        else 
           if(myrank.eq.1) then
              open(1001,file=fluxmofile,status="old",form="unformatted",access="direct",recl=molength)
              open(1002,file=fluxafile,status="old",form="unformatted",access="direct",recl=alength)
              do i=1,brareadsize
                 k=FluxSkipMult*((brabat-1)*BatchSize+i-1)+1
                 read(1001,rec=k) bramo(:,:,i) ;            read(1002,rec=k) braavec(:,:,:,i) 
              enddo
              close(1001);          close(1002)
           endif
           call mympibcast(bramo(:,:,1:brareadsize),1,totspfdim*brareadsize)
           call mympibcast(braavec(:,:,:,1:brareadsize),1,numr*numconfig*mcscfnum*brareadsize)
        endif
        call system_clock(btime)
        times(1)=times(1)+btime-atime;      times(2)=times(2)+btime-atime
        
!! loop over all time for the ket of the flux integral
        do kettime=1,ketreadsize

!           OFLWR "KETTIME ", kettime,ketreadsize; CFL

!! get the one-e half transformed matrix elements for this ket time
           call system_clock(atime)
           curtime=(ketbat-1)*BatchSize+kettime-1 
           whichbin=1
           do k=2,FluxNBins
              if(curtime.le.bintimes(k)) whichbin=k
           enddo
           moket=>ketmo(:,:,kettime);        aket=>ketavec(:,:,:,kettime)
           
           call flux_op_onee(moket,keop,peop,zdipop,xydipop,1)  !! 1 means flux
           if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0)then
              call flux_op_onee(moket,rekeop,repeop,rezdipop,rexydipop,2)  !! 1 means flux
           endif
           
           if (nonuc_checkflag.eq.0) then
              call flux_op_nuc(moket,yop,1)  !! 1 means flux
              if (nucfluxopt.ne.0) then
                 
                 call flux_op_nuc(moket,reyop,2)  !! 1 means flux
              endif
           endif
!! determine bounds of loop over bras and setup doing in parallel with mpi!
           if(brabat.lt.ketbat) then
              bratop=brareadsize
           else
              bratop=kettime
           endif

!!MAY2014           clow = (myrank-1)*bratop/nprocs+1;        chigh = myrank*bratop/nprocs
  
           clow = 1;        chigh = bratop
         
           allocate(tempgtau(mcscfnum,1:bratop))
           tempgtau=0d0
           call system_clock(btime);        times(4)=times(4)+btime-atime
           
           !! loop over all previous time for the bra of the flux integral
           do bratime=clow,chigh
              oldtime=(brabat-1)*BatchSize+bratime-1
              
!! biortho this pair of times!        
              call system_clock(itime)
              mobra=>bramo(:,:,bratime)
              abio(:,:,:)=braavec(:,:,:,bratime)
              
              call bioset(fluxgtaubiovar,smo,numr)
              call biortho(mobra,moket,mobio,abio(:,:,1),fluxgtaubiovar)
              
              do imc=2,mcscfnum
                 call biotransform(mobra,mobio,abio(:,:,imc),fluxgtaubiovar)
              enddo
              
              call system_clock(jtime);          times(3)=times(3)+jtime-itime

!! complete the one-e potential and kinetic energy matrix elements           
              do i=1,nspf
                 do k=1,nspf
                    ke(k,i) = dot(mobio(:,k),keop(:,i),spfsize)
                    pe(k,i) = dot(mobio(:,k),peop(:,i),spfsize)
                    if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0)then
                       reke(k,i) = dot(mobio(:,k),rekeop(:,i),spfsize)
                       repe(k,i) = dot(mobio(:,k),repeop(:,i),spfsize)
                    endif
                    if (tdflag.ne.0) then
                       zdip(k,i) = dot(mobio(:,k),zdipop(:,i),spfsize)
                       xydip(k,i) = dot(mobio(:,k),xydipop(:,i),spfsize)
                       if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0)then
                          rezdip(k,i) = dot(mobio(:,k),rezdipop(:,i),spfsize)
                          rexydip(k,i) = dot(mobio(:,k),rexydipop(:,i),spfsize)
                       endif
                    endif
                 enddo
              enddo
              if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0) then
                 do i=1,nspf
                    do k=1,nspf
                       yderiv(k,i) = dot(mobio(:,k),yop(:,i),spfsize)
                    enddo
                 enddo
                 do i=1,nspf
                    do k=1,nspf
                       reyderiv(k,i) = dot(mobio(:,k),reyop(:,i),spfsize)
                    enddo
                 enddo
              endif

              if (parorbsplit.eq.3) then
                 call mympireduce(ke,nspf**2)
                 call mympireduce(pe,nspf**2)
                 if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0)then
                    call mympireduce(reke,nspf**2)
                    call mympireduce(repe,nspf**2)
                 endif
                 if (tdflag.ne.0) then
                    call mympireduce(zdip,nspf**2)
                    call mympireduce(xydip,nspf**2)
                    if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0)then
                       call mympireduce(rezdip,nspf**2)
                       call mympireduce(rexydip,nspf**2)
                    endif
                 endif
                 if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0) then
                    call mympireduce(yderiv,nspf**2)
                    call mympireduce(reyderiv,nspf**2)
                 endif
              endif

              call system_clock(itime);          times(4)=times(4)+itime-jtime

!! get the two-e contribution, boo this is slow and we don't like it!           
              
              V2=0d0
              if(FluxOpType.eq.0) then  !!.and.onee_checkflag/=1) then

                 call noparorbsupport("call flux_op_twoe")

                 call flux_op_twoe(mobio,moket,V2,1) !! one means flux
                 if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0)then
                    call flux_op_twoe(mobio,moket,reV2,2)  !! one means flux
                 endif
              endif
              call system_clock(jtime)
              times(5)=times(5)+jtime-itime
!! evaluate the actual g(tau) expression
              xyfac=0d0; zfac=0d0

!! for length at least, trying flux = [H(t),heaviside] = [H(t'),heaviside] = average    !! WAS BUG 01-2014 no dt
              if (tdflag.ne.1) then
!!! STILL NEED Y
!!! STILL NEED Y
!!! STILL NEED Y
!!! STILL NEED Y
                 call vectdpot(curtime*dt,pots1)
                 call vectdpot(oldtime*dt,pots2)
                 xyfac=real( pots1(1) + pots2(1) ,8) /2d0
                 zfac=real( pots1(3) + pots2(3) ,8) /2d0

!                 xyfac=( xtdpot(curtime*dt) + xtdpot(oldtime*dt) ) /2d0
!                 zfac= ( ztdpot(curtime*dt) + ztdpot(oldtime*dt) ) /2d0
              endif
              do imc=1,mcscfnum
                 fluxevalval(imc) = fluxeval(abio,aket,ke,pe,V2,zdip,zfac,xydip,xyfac,yderiv,1,imc) * dt  !! 1 means flux
              enddo
          
              if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0)then
                 do imc=1,mcscfnum
                    fluxevalval(imc) = fluxevalval(imc)+fluxeval(abio,aket,reke,repe,reV2,rezdip,zfac,rexydip,xyfac,reyderiv,2,imc) * dt  !! 2 means flux
                 enddo
              endif
              tempgtau(:,bratime) = tempgtau(:,bratime) + fluxevalval(:)
              nullify(mobra) !!,abra)
              call system_clock(itime);          times(6)=times(6)+itime-jtime
           enddo !! end loop over specific bras
           
!! broadcast the parallelized data
           call system_clock(itime)
           do jproc=1,nprocs

              clow=(jproc-1)*bratop/nprocs+1;          chigh=jproc*bratop/nprocs
              cnum=chigh-clow+1

              if (cnum.gt.0) then
                 call mympibcast(tempgtau(:,clow:),jproc,cnum*mcscfnum)
              endif
           enddo
!! sum up the parallelized data
           do bratime=1,bratop
              oldtime=(brabat-1)*BatchSize+bratime-1
              tau=curtime-oldtime; 
              do imc=1,mcscfnum
                 gtau(tau,1:whichbin,imc) = gtau(tau,1:whichbin,imc) + tempgtau(imc,bratime)
              enddo
           enddo
           deallocate(tempgtau);        nullify(moket,aket)
           call system_clock(jtime);        times(6)=times(6)+jtime-itime
           
!! create the xsec as it should be now
           doFT=.false.
           if(brabat.eq.ketbat) then
              do k=1,FluxNBins
                 if(curtime.eq.bintimes(k)) doFT=.true.
              enddo
           endif
           if(doFT) then
              call openfile
              write(mpifileptr,'(A47,F10.4)') " Taking the FT of g(tau) to get xsection at T= ",real(curtime,8)*dt; CFL
              if(myrank.eq.1) then
                 open(1004,file=spifile,status="replace",action="readwrite",position="rewind")
                 write(1004,*)
                 write(1004,*) "# Omega; pulse ft; flux at t= ... "
                 write(1004,'(A8, A8, 100F36.5)') "#  ", " ", real(bintimes(1:FluxNBins),8)*dt
                 write(1004,*)
                 
                 open(1005,file="Dat/xsec.alt.dat",status="replace",action="readwrite",position="rewind")
                 write(1005,*)
                 write(1005,*) "# Omega; pulse ft; flux at t= ... "
                 write(1005,'(A8, A8, 100F36.5)') "#  ", " ", real(bintimes(1:FluxNBins),8)*dt
                 write(1005,*)
              endif
!! get the FT of gtau
              allocate(FTgtau(FluxNBins,mcscfnum))
              do k=1,nEFlux
                 wfi=EFluxLo+(k-1)*dEFlux
                 FTgtau=0d0
!!092812            if(wfi.ge.0d0) then
                 Energy=wfi+ceground
                 FTgtau(:,:) = gtau(0,:,:) * dt 
                 do i=1,curtime
                    if(FluxSineOpt.eq.0) then
                       ffac = dt
                    else
                       ffac = (cos(real(i)/real(nt)*PI/2d0)**FluxSineOpt) * dt 
                    endif
                    FTgtau(:,:) = FTgtau(:,:) + exp((0d0,1d0)*Energy*i*dt) * gtau(i,:,:) * ffac
                    FTgtau(:,:) = FTgtau(:,:) + exp((0d0,1d0)*ALLCON(Energy)*(-i)*dt) * ALLCON(gtau(i,:,:)) * ffac
                 enddo
                 if(myrank.eq.1) then
                    write(1004,'(F8.4,100E18.6)') wfi, 1d0/xsec(k), xsec(k)*FTgtau(1:2,:)
                    write(1005,'(F8.4,100E18.6)') wfi, 1d0/xsec(k), xsec(k)*FTgtau(1,:), FTgtau(1,:)
                 endif
              enddo
              deallocate(FTgtau)
              if(myrank.eq.1) then
                 close(1004); close(1005)
              endif
           
              allocate(ftgtau(-curtime:curtime,mcscfnum),ftgtausave(-curtime:curtime,mcscfnum), wsave(20*(2*curtime+1)+30))
              ftgtau(:,:)=0d0; ftgtausave(:,:)=0d0
              
              do i=0,curtime
                 ftgtau(i,:) = ALLCON(gtau(i,1,:))   * cos(real(i,8)/real(curtime,8) * piover2) * exp((0.d0,-1.d0)*ALLCON(ceground)*par_timestep*FluxInterval*FluxSkipMult*i)
              enddo
              do i=0,curtime
                 ftgtau(-i,:) = (gtau(i,1,:))  * cos(real(i,8)/real(curtime,8) * piover2) * exp((0.d0,1.d0)*ceground*par_timestep*FluxInterval*FluxSkipMult*i)
              enddo

              if (myrank.eq.1) then
                 open(171,file="Dat/myGTau.Dat",status="unknown");          write(171,*) "#   ", curtime
                 do i=-curtime,curtime
                    write(171,'(F18.12, T22, 400E20.8)')  i*par_timestep*FluxInterval*FluxSkipMult, ftgtau(i,:)
                 enddo
                 close(171)
              endif
              do imc=1,mcscfnum

                 call zffti(2*curtime+1,wsave);          call zfftf(2*curtime+1,ftgtau(-curtime,imc),wsave)

              enddo

              ftgtau(-curtime:curtime,:)=ftgtau(-curtime:curtime,:)*par_timestep*FluxInterval*FluxSkipMult
              
              !! NEW 0603
              do i=-curtime,curtime
                 ftgtau(i,:)=ftgtau(i,:)*exp((0.d0,1.d0)*(curtime+i)*curtime*2*pi/real(2*curtime+1))
              enddo
              do i=1,curtime
                 ftgtausave(i-curtime-1,:) = ftgtau(i,:)
              enddo
              do i=-curtime,0
                 ftgtausave(curtime+i,:)=ftgtau(i,:)
              enddo
              
              estep=4*piover2/par_timestep/fluxinterval/fluxskipmult/(2*curtime+1)
              
              flag=0
              if (noftflag.eq.0) then
                 flag=1
                 do ipulse=1,numpulses
                    if ((pulsetype(ipulse).ne.1.and.pulsetype(ipulse).ne.2).or.chirp(ipulse).ne.0d0) then
                       flag=0
                    endif
                 enddo
              endif

              if (myrank.eq.1) then
                 open(171,file="Dat/FTGtau.Dat", status="unknown")
                 do i=-curtime,curtime
                    myft=1d0
                    if (flag.eq.1) then
                       myft=abs(tdpotft(i*Estep))**2
                    endif
!             write(171,'(F18.12, T22, 400E20.8)')  i*Estep, abs(ftgtausave(i)), ftgtausave(i), ftgtausave(i)/myft, myft
                    write(171,'(F18.12, T22, 400E20.8)')  i*Estep,  ftgtausave(i,:), ftgtausave(i,:)/myft, myft
                 enddo
                 close(171)
              endif
              deallocate(ftgtau,ftgtausave, wsave)
              call system_clock(btime);        times(7)=times(7)+btime-jtime;        times(1)=times(1)+btime-atime
!! write out times and totals
              if (brabat.eq.ketbat) then
                 call openfile
                 write(mpifileptr,'(A28,F10.4)') " Timing statistics as of T= ",real(curtime,8)*dt
                 write(mpifileptr,'(100A10)') "Times: ", "All", "Read","Biorth", "One-e", "Two-e", "Fluxeval", "FT gtau"
                 write(mpifileptr,'(A10,100I10)') " ", times(1:7)/100; CFL
              endif

           endif !! doFT
        
!! only do this after we are sure we've gone through every bra
           if(brabat.eq.ketbat) then
              if (myrank.eq.1) then
                 open(454, file="Dat/KVLsum.dat", status="old", position="append")
                 
                 write(454,'(I5,100F18.12)') curtime, curtime*dt, gtau(0,1,:);    
                 close(454)
              endif
           endif
        enddo !! end loop overspecific kets
     enddo !! end loop over bra batches
  enddo !! end loop over ket batches
  deallocate(bramo,braavec,ketmo,ketavec,mobio,abio,keop,peop,ke,pe,V2,gtau,xsec,bintimes,yderiv)
  deallocate(yop)

end subroutine fluxgtau


!! begin the flux matrix element and contraction routine section


subroutine getFTpulse(xsec)
  use parameters
  implicit none
  integer :: i,k,numstep
  real*8 :: time,wfi,xsec(nEFlux)
  DATATYPE :: alltdpot
  complex*16 :: FTpulse
!! Get the FT of the pulse and the coefficeints for the xsec 
  xsec=0d0

  if (noftflag.ne.0) then   !! for, for instance, propagating an ansatz core hole. no ft factor.
     xsec=1d0
     return
  endif

    do i=1,numpulses
!! DJH APR 2014
       numstep=ceiling(1d0/omega(i)/200d0/par_timestep)
    enddo

  do k=1,nEFlux
    wfi=EFluxLo+(k-1)*dEFlux
    FTpulse=0d0
!! get the FT of the pulse
    do i=1,numpulses
!! DJH APR 2014
       numstep=ceiling(1d0/omega(i)/200d0/par_timestep)
       time=0d0
       do while(time.lt.pi/omega(i)+pulsestart(i))
          FTpulse=FTpulse+par_timestep*numstep*alltdpot(time,4)*exp((0d0,1d0)*wfi*time)
          time=time+par_timestep*numstep
       enddo
    enddo
!! finish the cross section with the right constants
    if(velflag.eq.0) then !! length gauge
       xsec(k) = ( ((5.291772108d0**2) / 3d0) * (2d0 * PI / 1.37036d2) * wfi ) / real(conjg(FTpulse)*FTpulse,8)
    else !! velocity gauge
       xsec(k) = ( ((5.291772108d0**2) / 3d0) * (2d0 * PI / 1.37036d2) / wfi ) / real(conjg(FTpulse)*FTpulse,8)
    endif
  enddo
end subroutine getFTpulse


subroutine flux_op_onee(inspfs,keop,peop,zdipop,xydipop,flag) !! flag=1, flux (imag); flag=2, flux (real); flag=0, all  0 not used
  use parameters
  implicit none
  DATATYPE, intent(in) :: inspfs(spfsize,nspf)
  DATATYPE ::  keop(spfsize,nspf),peop(spfsize,nspf),zdipop(spfsize,nspf),xydipop(spfsize,nspf)
  integer :: ispf,flag

!! initialize
  keop=0d0; peop=0d0; zdipop=0d0; xydipop=0d0

!! the kinetic energy
  do ispf=1,nspf
     select case(flag)
     case(1)
        call mult_imke(inspfs(:,ispf),keop(:,ispf))
!        print *, "KEOP",keop(:,ispf); stop
     case(2)
        call mult_reke(inspfs(:,ispf),keop(:,ispf))
     case default
        call mult_ke(inspfs(:,ispf),keop(:,ispf),1,"booga",2)
     end select
  enddo

!! the one-e potential energy 
  do ispf=1,nspf
     select case(flag)
     case(1)
        if(FluxOpType.eq.0.or.FluxOpType.eq.2) then
           call mult_impot(inspfs(:,ispf),peop(:,ispf))
        else if(FluxOpType.eq.1) then
           call mult_imhalfniumpot(inspfs(:,ispf),peop(:,ispf))
        endif 
     case(2)
        if(FluxOpType.eq.0.or.FluxOpType.eq.2) then
           call mult_repot(inspfs(:,ispf),peop(:,ispf))
        else if(FluxOpType.eq.1) then
           call mult_rehalfniumpot(inspfs(:,ispf),peop(:,ispf))
        endif 
     case default
        call mult_pot(inspfs(:,ispf),peop(:,ispf))
     end select

     if (tdflag.ne.0) then
        if (flag.ne.1) then
           OFLWR "Programmmmm meeee"; CFLST
        endif

!! NO Y POLARIZATION YET
        select case (velflag)
        case (0)
           call mult_imzdipole(inspfs(:,ispf),zdipop(:,ispf))
           call mult_imxdipole(inspfs(:,ispf),xydipop(:,ispf))
        case (1)

           call noparorbsupport("call imvelmultiply")

           call imvelmultiply(inspfs(:,ispf),zdipop(:,ispf), 0d0,0d0,1d0)
           call imvelmultiply(inspfs(:,ispf),xydipop(:,ispf), 1d0,0d0,0d0)
        end select
     endif
  enddo

!! scale correctly and clear memory
  if(flag.ne.0) then
    zdipop=zdipop*(-2d0)
    xydipop=xydipop*(-2d0)
    peop=peop*(-2d0)
    keop=keop*(-2d0)
  endif
end subroutine flux_op_onee


!! operates with y derivative cross term.  Assumes y^2 and l^2 are in the one electron hamiltonian.

!!  yop complex antisymmetric  hermitian part is imaginary   times antihermitian op in R (assuming no R ecs scaling, 
!!    as in all these routines)

!! fortran imag() returns real value so multiply by i


subroutine flux_op_nuc(inspfs,yop,flag) !! flag=1, flux (imag); flag=0, all    2 flux (imag)
  use parameters
  implicit none
  DATATYPE, intent(in) :: inspfs(spfsize,nspf)
  DATATYPE ::  yop(spfsize,nspf)
  integer :: ispf,flag

  call noparorbsupport("in flux_op_nuc")

!! the kinetic energy
  do ispf=1,nspf
     select case(flag)
     case(1)
        call op_imyderiv(inspfs(:,ispf),yop(:,ispf))
     case(2)
        call op_reyderiv(inspfs(:,ispf),yop(:,ispf))
     case default
        call op_yderiv(inspfs(:,ispf),yop(:,ispf))
     end select
  end do

!! scale correctly and clear memory
  if(flag.ne.0) then
    yop=yop*(-2d0)
  endif
end subroutine flux_op_nuc


!! maybe double check if flag.ne.1.

subroutine flux_op_twoe(mobra,moket,V2,flag)  !! flag=1 means flux, otherwise whole op
  use parameters
  implicit none
  DATATYPE :: mobra(spfsize,nspf),moket(spfsize,nspf)
  DATATYPE :: V2(nspf,nspf,nspf,nspf)
  integer :: flag

  if (flag.ne.1) then
     OFLWR "Doublecheck flag.ne.1 ok in flux_op_twoe"; CFLST
  endif

  call noparorbsupport("call call_flux_op_twoe")

  call call_flux_op_twoe(mobra,moket,V2,flag)

end subroutine flux_op_twoe


!! with flag=1 is sent matrix elements of imaginary part of electronic operators (real valued; function "imaginary part" returns real)
!!      multiply these by the real part of the bond operators
!! with flag=2 vice versa
!! with flag=0 does full hamiltonian no real or imag part


function fluxeval(abra,aket,ke,pe,V2,zdip,zfac,xydip,xyfac,yderiv,flag,imc)   
  use parameters 
  implicit none
  integer :: flag,imc
  real*8 :: xyfac, zfac
  DATATYPE :: abra(numconfig,numr,mcscfnum),aket(numconfig,numr,mcscfnum),fluxeval,fluxeval00
  DATATYPE :: ke(nspf,nspf),pe(nspf,nspf),V2(nspf,nspf,nspf,nspf),yderiv(nspf,nspf),xydip(nspf,nspf),zdip(nspf,nspf)

  fluxeval=fluxeval00(abra(:,:,imc),aket(:,:,imc),ke,pe,V2,zdip,zfac,xydip,xyfac,yderiv,flag,nucfluxflag)

end function fluxeval


function fluxeval00(abra,aket,ke,pe,V2,zdip,zfac,xydip,xyfac,yderiv,flag,ipart)   
!! flag=1 means flux (take into account fluxoptype) otherwise whole thing

!! ipart=0 do all   1 elec only   (diagonal nuclear ke ; cap; herm(Y+3/2)*anti(d/dR term) 
!!   (herm=hermitian part (H + H^dag /2))  that's imag(Y+3/2)...   yes this is correct

!! determine the 2-electron matrix elements in the orbital basis for the flux operator i(H-H^{\dag}) 
!! input :
!! abra - the A vector for the bra space
!! aket - the A vector for the ket space
!! ke - the 1-electron matrix elements corresponding with kinetic energy (contract with 1/R**2) 
!! pe - the 1-electron matrix elements corresponding with potential energy (contract with 1/R) 
!! yderiv - the 1-electron matrix elements of Y+3/2 operator for KE cross term
!! V2 - the 2-electron matrix elements corresponding with potential energy (contract with 1/R) 
!! output : 
!! fluxeval - this is the actual value for <\Psi_{bra space}| i(H-H^{\dag}) |\Psi_{ket space}>

  use parameters 
  use opmod  !! rkemod & proderivmod
  use walkmod
  implicit none

  integer :: bra,ket,r,flag,rr, iiflag=0,ipart
  real*8 :: xyfac, zfac
  DATATYPE :: abra(numconfig,numr),aket(numconfig,numr),INVR,IIINVR,INVRSQ,fluxeval00,PRODERIV,BONDKE,PULSEFAC,csum
  DATATYPE :: ke(nspf,nspf),pe(nspf,nspf),V2(nspf,nspf,nspf,nspf),yderiv(nspf,nspf),xydip(nspf,nspf),zdip(nspf,nspf)

  if (flag.ne.1.or.fluxoptype.ne.1.or.ipart.ne.0) then
     OFLWR "maybe checkme debug"; CFLST
  endif

  iiflag=iiflag+1

  fluxeval00 = 0d0


  IF (1==0) THEN   !! 1==0 !! TEMP???



  do bra=botwalk,topwalk
!! A) do all diagonal terms

    IIINVR=0d0;    INVR=0d0;    INVRSQ=0d0;    PULSEFAC=0d0

    do r=1,numr
       select case(flag)
       case(1)
          INVR = INVR + CONJUGATE(abra(bra,r)) * aket(bra,r) * real(1.d0/bondpoints(r),8)
          IIINVR = IIINVR + CONJUGATE(abra(bra,r)) * aket(bra,r) * imag((0d0,0d0)+1.d0/bondpoints(r))
          INVRSQ = INVRSQ + CONJUGATE(abra(bra,r)) * aket(bra,r) * real(1d0/bondpoints(r)**2,8)
       case(2)
          INVR = INVR + CONJUGATE(abra(bra,r)) * aket(bra,r) * imag((0d0,0d0)+1d0/bondpoints(r)**2)
          INVRSQ = INVRSQ + CONJUGATE(abra(bra,r)) * aket(bra,r) * imag((0d0,0d0)+1d0/bondpoints(r)**2)
       case default
          INVR = INVR + CONJUGATE(abra(bra,r)) * aket(bra,r) / (bondpoints(r))
          IIINVR = IIINVR + CONJUGATE(abra(bra,r)) * aket(bra,r) / (bondpoints(r))
          INVRSQ = INVRSQ + CONJUGATE(abra(bra,r)) * aket(bra,r) / (bondpoints(r)**2)
       end select

       if (tdflag.ne.0) then

         if (velflag==0) then
            csum=bondpoints(r)
         else
            csum=1d0/bondpoints(r)
         endif

         select case(flag)
         case(1)
            PULSEFAC= PULSEFAC + CONJUGATE(abra(bra,r)) * aket(bra,r) * imag((0d0,0d0)+csum)
         case(2)
            PULSEFAC= PULSEFAC + CONJUGATE(abra(bra,r)) * aket(bra,r) * real(csum)
         case default
            PULSEFAC= PULSEFAC + CONJUGATE(abra(bra,r)) * aket(bra,r) * csum
         end select
       endif
    enddo

    BONDKE=0d0;    PRODERIV=0d0

    if (ipart.eq.0.or.ipart.eq.2) then

    if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0) then
       select case(flag)
       case(1)
          do rr=1,numr
             do r=1,numr
                BONDKE = BONDKE + CONJUGATE(abra(bra,r)) * aket(bra,rr) * imag((0d0,0d0)+rkemod(r,rr))
                PRODERIV = PRODERIV + CONJUGATE(abra(bra,r)) * aket(bra,rr) * real(proderivmod(r,rr),8)
             enddo
          enddo
       case(2)
          do rr=1,numr
             do r=1,numr
                PRODERIV = PRODERIV + CONJUGATE(abra(bra,r)) * aket(bra,rr) * imag((0d0,0d0)+proderivmod(r,rr))
             enddo
          enddo
       case default
          do rr=1,numr
             do r=1,numr
                BONDKE = BONDKE + CONJUGATE(abra(bra,r)) * aket(bra,rr) * rkemod(r,rr)
                PRODERIV = PRODERIV + CONJUGATE(abra(bra,r)) * aket(bra,rr) * proderivmod(r,rr) 
             enddo
          enddo
       end select

       if (flag.eq.1.or.flag.eq.2) then
          fluxeval00=fluxeval00+BONDKE * (-2d0)
       else
          fluxeval00=fluxeval00+BONDKE 
       endif

!! with flag=1 we should have imag(yderiv) thats hermitian part yderiv for nuclear part  (ipart=2)
!! yes this correct   likewise real(yderiv), imag(R derivative) flag=2   for elect
!! for flag=0 only do once, so do for ipart=1
!!
!! so then for capflag do imag part of yderiv times hermitian part of yderiv  (flag=1)
!!     but don't do flag.eq.2  (which part is zero; so no need for new logic )  SO WEIRD

! sparate diagonal REMOVED 09-2014 (WALKS SIMPLIFIED)
!       do i=1,numelec
!          !! A.1) 1 electron
!          ii=configlist((i-1)*2+1,bra);          is=configlist((i-1)*2+2,bra)
!          fluxeval00 = fluxeval00 + PRODERIV * yderiv(ii,ii) 
!       enddo

    endif
    endif  !! NONUC


    if (ipart.eq.0.or.ipart.eq.1) then
    select case(flag)
    case(1)
       fluxeval00=fluxeval00+IIINVR*nucrepulsion * (-2)  
    case(0)
       fluxeval00=fluxeval00+IIINVR*nucrepulsion
    end select



! separate diagonal REMOVED 09-2014 (walks simplified)
!    do i=1,numelec
!!! A.1) 1 electron
!      ii=configlist((i-1)*2+1,bra)
!      is=configlist((i-1)*2+2,bra)
!      fluxeval00 = fluxeval00 + INVRSQ * ke(ii,ii) + INVR * pe(ii,ii)
!      if (tdflag.ne.0) then
!         fluxeval00 = fluxeval00 + PULSEFAC * zdip(ii,ii) * zfac
!         fluxeval00 = fluxeval00 + PULSEFAC * xydip(ii,ii) * xyfac
!      endif
!!! A.2) 2 electron
!      if(((FluxOpType.eq.0).or.(flag.ne.1))) then  !!.and.onee_checkflag/=1) then
!        do j=i+1,numelec
!          jj=configlist((j-1)*2+1,bra)
!          js=configlist((j-1)*2+2,bra)
!          fluxeval00 = fluxeval00 + IIINVR * V2(ii,ii,jj,jj)
!          if(is.eq.js) fluxeval00 = fluxeval00 - IIINVR * V2(ii,jj,jj,ii)
!        enddo
!      endif
!    enddo

    endif   !! ipart

 enddo
endif   !! 1==0


!! B) do all the single walks

  do bra=botwalk,topwalk

    do ket=1,numsinglewalks(bra)

      INVR=0d0;      INVRSQ=0d0;      PULSEFAC=0d0
      do r=1,numr
         select case(flag)
         case(1)
            INVR = INVR + CONJUGATE(abra(bra,r)) * aket(singlewalk(ket,bra),r) * real(1d0/ (bondpoints(r)),8)
            INVRSQ = INVRSQ + CONJUGATE(abra(bra,r)) * aket(singlewalk(ket,bra),r) * real(1d0 / (bondpoints(r)**2),8)
         case(2)
            INVR = INVR + CONJUGATE(abra(bra,r)) * aket(singlewalk(ket,bra),r) *imag((0d0,0d0)+1d0/bondpoints(r))
            INVRSQ = INVRSQ + CONJUGATE(abra(bra,r)) * aket(singlewalk(ket,bra),r) * imag((0d0,0d0)+1d0/ (bondpoints(r)**2))
         case default
            INVR = INVR + CONJUGATE(abra(bra,r)) * aket(singlewalk(ket,bra),r) / (bondpoints(r))
            INVRSQ = INVRSQ + CONJUGATE(abra(bra,r)) * aket(singlewalk(ket,bra),r) / (bondpoints(r)**2)
         end select
        if (tdflag.ne.0) then
           if (velflag==0) then
              csum=bondpoints(r)
           else
              csum=1d0/bondpoints(r)
           endif
           select case(flag)
           case(1)
              PULSEFAC= PULSEFAC + CONJUGATE(abra(bra,r)) * aket(singlewalk(ket,bra),r) * real(csum,8)
           case(2)
              PULSEFAC= PULSEFAC + CONJUGATE(abra(bra,r)) * aket(singlewalk(ket,bra),r) * imag((0d0,0d0)+csum)
           case default
              PULSEFAC= PULSEFAC + CONJUGATE(abra(bra,r)) * aket(singlewalk(ket,bra),r) * csum
           end select
         endif
      enddo

      INVR = INVR * singlewalkdirphase(ket,bra)
      INVRSQ = INVRSQ * singlewalkdirphase(ket,bra)
      PULSEFAC = PULSEFAC * singlewalkdirphase(ket,bra)
      BONDKE=0d0;      PRODERIV=0d0

      if (ipart.eq.0.or.ipart.eq.2) then
      if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0) then
         select case(flag)
         case(1)
            do rr=1,numr
               do r=1,numr
                  BONDKE = BONDKE + CONJUGATE(abra(bra,r)) * aket(singlewalk(ket,bra),rr) * imag((0d0,0d0)+rkemod(r,rr))
                  PRODERIV = PRODERIV + CONJUGATE(abra(bra,r)) * aket(singlewalk(ket,bra),rr) * real(proderivmod(r,rr),8)
               enddo
            enddo
         case(2)
            do rr=1,numr
               do r=1,numr
                  PRODERIV = PRODERIV + CONJUGATE(abra(bra,r)) * aket(singlewalk(ket,bra),rr) * imag((0d0,0d0)+proderivmod(r,rr))
               enddo
            enddo
         case default
            do rr=1,numr
               do r=1,numr
                  BONDKE = BONDKE + CONJUGATE(abra(bra,r)) * aket(singlewalk(ket,bra),rr) * rkemod(r,rr)
                  PRODERIV = PRODERIV + CONJUGATE(abra(bra,r)) * aket(singlewalk(ket,bra),rr) * proderivmod(r,rr)
               enddo
            enddo
         end select

         BONDKE = BONDKE * singlewalkdirphase(ket,bra)
         PRODERIV = PRODERIV * singlewalkdirphase(ket,bra)

         if (flag.eq.1.or.flag.eq.2) then
            fluxeval00 = fluxeval00 + BONDKE   * (-2d0)
         else
            fluxeval00 = fluxeval00 + BONDKE  
         endif
         fluxeval00 = fluxeval00 + PRODERIV * yderiv(singlewalkopspf(1,ket,bra),singlewalkopspf(2,ket,bra))
      endif
      endif
      if (ipart.eq.1.or.ipart.eq.0) then

!! B.1) 1 electron
      fluxeval00 = fluxeval00 + INVRSQ * ke(singlewalkopspf(1,ket,bra),singlewalkopspf(2,ket,bra)) & !! conjugates OK
                          + INVR * pe(singlewalkopspf(1,ket,bra),singlewalkopspf(2,ket,bra))

      if (tdflag.ne.0) then
         fluxeval00 = fluxeval00 + PULSEFAC * zdip(singlewalkopspf(1,ket,bra),singlewalkopspf(2,ket,bra))  * zfac
         fluxeval00 = fluxeval00 + PULSEFAC * xydip(singlewalkopspf(1,ket,bra),singlewalkopspf(2,ket,bra))  * xyfac
      endif

!! B.2) 2 electron  
!! singlewalks removed REMOVED 09-2014 (walks simplified)

      endif  !! ipart
    enddo

    if(ipart.eq.0.or.ipart.eq.1) then

!! C) do all the double walks (only 2 electron)
    if(((FluxOpType.eq.0).or.(flag.ne.1))) then !!.and.onee_checkflag/=1) then

       OFLWR "MAYBE DOUBLE CHECKME"; CFLST

      do ket=1,numdoublewalks(bra)
        INVR=0d0
        select case(flag)
        case(1)
           do r=1,numr
              INVR = INVR + CONJUGATE(abra(bra,r)) * aket(doublewalk(ket,bra),r) *real(1d0/ (bondpoints(r)),8)
           enddo
        case(2)
           do r=1,numr
              INVR = INVR + CONJUGATE(abra(bra,r)) * aket(doublewalk(ket,bra),r) *imag((0d0,0d0)+1d0/bondpoints(r))
           enddo
        case default
           do r=1,numr
              INVR = INVR + CONJUGATE(abra(bra,r)) * aket(doublewalk(ket,bra),r) / (bondpoints(r))
           enddo
        end select
          fluxeval00 = fluxeval00 + INVR * V2(doublewalkdirspf(1,ket,bra), &
            doublewalkdirspf(2,ket,bra), &
            doublewalkdirspf(3,ket,bra), &
            doublewalkdirspf(4,ket,bra)) * &
            doublewalkdirphase(ket,bra)
      enddo
    endif
    endif  !! ipart
  enddo

  if (sparseconfigflag.ne.0) then
     call mympireduceone(fluxeval00)
  endif

end function fluxeval00
