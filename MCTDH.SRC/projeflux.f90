
!! PROJECTED FLUX BY KVL

#include "Definitions.INC"

module projefluxmod !! needed for cation walks and bi-orthonormalization

!! variables so we don't have to pass the single cation walks for the target states to the MCTDHF wavefunction 
!! target readable variables
  integer :: targetms, targetrestrictflag, targetspinproject, tnumconfig,tmaxsinglewalks,targetspinval,cgflag=0
  integer, allocatable :: tconfiglist(:,:),tnumsinglewalks(:),tsinglewalk(:,:),tsinglewalkopspf(:,:,:), &
       tsinglewalkdirphase(:,:), pphase1(:,:) ,pspf1(:,:,:) ,numpwalk1(:) ,pwalk1(:,:)  
  DATATYPE, allocatable :: tmo(:,:,:),ta(:,:,:)
!! cation walks variables
  integer :: maxpwalk1 !!,maxpwalk2

end module projefluxmod


!! here we are doing walks from our BO target state to our wavefunction \Psi(t)
subroutine projeflux_singlewalks()
  use parameters
  use aarrmod
  use projefluxmod
  use mpimod
  implicit none
  integer :: clow,chigh,jproc,cnum, iconfig,idof,iindex,iind,iwalk,flag,getconfiguration,reorder,dirphase
  integer :: tempconfig(ndof),temporb(2)
  logical :: allowedconfig

!! get the number of single walks from the target N-1 e- state to our regular N e- state
  OFLWR "Getting cation single walks"; CFL

  allocate(numpwalk1(tnumconfig))
  numpwalk1(:)=0
  clow=(myrank-1)*tnumconfig/nprocs+1;  chigh=myrank*tnumconfig/nprocs
  do iconfig=clow,chigh
    iwalk=0
    do iindex=1,spftot
      tempconfig(3:ndof)=tconfiglist(:,iconfig);      temporb=aarr(iindex,nspf)
      flag=0
      do idof=2,numelec 
        if(iind(tempconfig(idof*2-1:idof*2)).eq.iindex) flag=1
      enddo
      if(flag.eq.0) then
        tempconfig(1:2)=temporb(:);        dirphase=reorder(tempconfig)
        if(allowedconfig(tempconfig)) iwalk=iwalk+1
      endif
    enddo
    numpwalk1(iconfig)=iwalk
  enddo
!! mpi broadcasting of number of single walks
  do jproc=1,nprocs
    clow=(jproc-1)*tnumconfig/nprocs+1;    chigh=jproc*tnumconfig/nprocs;    cnum=chigh-clow+1
    call mympiibcast(numpwalk1(clow:),jproc,cnum)
  enddo
!! figure out the maximum number of single target walks
  maxpwalk1=0
  do iconfig=1,tnumconfig
    if(maxpwalk1.lt.numpwalk1(iconfig)) maxpwalk1=numpwalk1(iconfig)
  enddo
  OFLWR "Max # single walks from cation state is ",maxpwalk1;CFL

  allocate(pphase1(maxpwalk1,tnumconfig),pwalk1(maxpwalk1,tnumconfig),pspf1(2,maxpwalk1,tnumconfig))


  pwalk1=0;  pspf1=0;  pphase1=0
  clow=(myrank-1)*tnumconfig/nprocs+1;  chigh=myrank*tnumconfig/nprocs

  do iconfig=clow,chigh
    iwalk=0
    do iindex=1,spftot
      tempconfig(3:ndof)=tconfiglist(:,iconfig);      temporb=aarr(iindex,nspf);      flag=0
      do idof=2,numelec
        if(iind(tempconfig(idof*2-1:idof*2)).eq.iindex) flag=1
      enddo
      if(flag.eq.0) then
        tempconfig(1:2)=temporb(:);        dirphase=reorder(tempconfig)
        if(allowedconfig(tempconfig)) then
          iwalk=iwalk+1;          pwalk1(iwalk,iconfig)=getconfiguration(tempconfig)
          pspf1(:,iwalk,iconfig)=temporb;          pphase1(iwalk,iconfig)=dirphase
        endif
      endif
    enddo
!! check and make sure no bad walk error
    if(numpwalk1(iconfig).ne.iwalk) then
      OFLWR "TARGET SINGLE WALK ERROR";CFLST
    endif
  enddo

  OFLWR "DONE getting proj single walks, broadcasting"; CFL

!! mpi broadcasting of our actual single walks
  do jproc=1,nprocs
    clow=(jproc-1)*tnumconfig/nprocs+1;    chigh=jproc*tnumconfig/nprocs;    cnum=(chigh-clow+1)*maxpwalk1
    call mympiibcast(pwalk1(:,clow:),jproc,cnum);    call mympiibcast(pspf1(:,:,clow:),jproc,2*cnum)
    call mympiibcast(pphase1(:,clow:),jproc,cnum)
  enddo
  call mpibarrier()
  OFLWR "     ... done boroadcasting."; CFL

end subroutine projeflux_singlewalks


!! construct the one electron functions
subroutine projeflux_doproj(cata,neuta,mo,offset)
  use parameters
  use projefluxmod
  implicit none
  integer :: jconfig,iwalk,iconfig,ispf,ispin,iphase, offset
  DATATYPE :: cata(tnumconfig),neuta(numconfig),mo(spfsize,nspf), projwfn(spfsize,2)

!! make the single electron wfn
  projwfn=0d0
  do jconfig=1,tnumconfig
!! original spin loop
    do iwalk=1,numpwalk1(jconfig)
      iconfig=pwalk1(iwalk,jconfig);      ispf=pspf1(1,iwalk,jconfig)
      ispin=pspf1(2,iwalk,jconfig);      iphase=pphase1(iwalk,jconfig)

      projwfn(:,ispin) = projwfn(:,ispin) + CONJUGATE(cata(jconfig)) * neuta(iconfig) * mo(:,ispf) * iphase
    enddo
  enddo

!! write out when done for this time and r
  inquire (iolength=ispf) projwfn
!  open(1003,file="proj.flux.wfn.bin",status="unknown",form="unformatted",access="direct",recl=ispf)
  open(1003,file=projfluxfile,status="unknown",form="unformatted",access="direct",recl=ispf)
  write(1003,rec=offset) projwfn(:,:) ;  close(1003)

end subroutine projeflux_doproj


!! do the double time integral piece
subroutine projeflux_double_time_int(bintimes,xsec,mem,istate,nstate,imc,nt,dt)
  use parameters
  use projefluxmod  !! targetms, ..
  use mpimod
  implicit none  
  integer :: i,k,tlen,mem,istate,curtime,tau,whichbin,nt,ir ,imc,nstate, bintimes(FluxNBins),ipulse
  integer :: BatchSize,NBat,ketreadsize,brareadsize,ketbat,brabat,kettime,bratime,bratop,flag
 !! bintimes and xsec done outside so we don't have to redo over and over for each state
  real*8 :: doubleclebschsq,aa,bb,cc,xsec(nEFlux),MemTot,MemVal,Energy,dt,wfi,cgfac,ffac,myft,piover2,estep
  DATATYPE, allocatable,target :: bramo(:,:,:,:),ketmo(:,:,:,:),gtau(:,:)
  complex*16, allocatable :: ftgtau(:),ftgtausave(:),wsave(:)
  complex*16  :: myFTgtau(FluxNBins),tdpotft

  logical :: doFT
  DATATYPE :: dot
  character (len=4) :: xstate0,xmc0
  character (len=3) :: xstate1,xmc1

  if (ceground.eq.(0d0,0d0)) then
     OFLWR "Eground is ZERO.  Are you sure?  If want zero just make it small. \n     Otherwise need eground: initial state energy."; CFLST
  endif

  piover2=atan2(1d0,1d0)*2

              
!!!!  100912 NEED CG COEFFICIENT RATIO 
              
!! programmed this with restrictms, targetms.  Do have option to change spinrestrictval; should use this.
!!
!!    have high spin cation and high spin neutral.
!!
!!    | s s  >  =  sum_m_s'  < s s   | s' s' m_s'  >  | s' m_s' >  x  | 1 +/-1/2 >
!!
!!     neutral                                          cation       outgoing electron
!!     m_s = s
!!
!!    only have |s' s' > m_s'=s'  cation not |s' s'-1> m_s'=s'-1
!!
!!   
!!    if restrictms > targetms  then we couple high spin to high spin to get high spin, one CG coef = 1.
!!    
!!        otherwise we couple high spin target down to high spin neutral:
!!
!!         projection is spin down  x high spin target   e.g.   < s' s' |  s s  1/2 -1/2 > = A
!!         also have < s' s' | s s-1 1/2 1/2 > = B
!!  
!!         so mult cross section by ( A^2 + B^2 / A^2 )
!!    
!!   BUT if restrictflag is off, can't do this.  Must realize that you are only projecting one component of C.G. sum.  
!!     that would be a problem.  But if you were to calculate say the 3 lowest states for a triplet, you'd get all components.
!!     otherwise cross section is meaningless.
!! all arguments are x 2, integers for half spin.

!!$                 if (restrictms.lt.targetms.and.restrictflag.eq.1) then
!!$                    aa= doubleclebschsq(targetms,1,targetms,restrictms-targetms,restrictms)      !! what we're calculating  targetms-restrictms=+/-1
!!$                    bb= doubleclebschsq(targetms,1,targetms-2,restrictms-targetms+2,restrictms)      !! only one could be nonzero obviously
!!$                    cc= doubleclebschsq(targetms,1,targetms+2,restrictms-targetms-2,restrictms)      
!!$                    cgfac = (aa+bb+cc)/aa
!!$                 else
!!$                    cgfac=1
!!$                 endif

  if (cgflag.eq.0) then
     cgfac=1
  else
     aa= doubleclebschsq(targetspinval,1,targetms,restrictms-targetms,spinrestrictval)      !! what we're calculating  targetms-restrictms=+/-1
     bb= doubleclebschsq(targetspinval,1,targetms-2,restrictms-targetms+2,spinrestrictval)      !! only one could be nonzero obviously
     cc= doubleclebschsq(targetspinval,1,targetms+2,restrictms-targetms-2,spinrestrictval)      
     cgfac = (aa+bb+cc)/aa
  endif




!! determine if we should do batching or not
!! 250,000 words/MB, real*8 2words/#, complex*16 4words/#

  write(xstate0,'(I4)') istate+1000;  write(xmc0,'(I4)') imc+1000
  xstate1=xstate0(2:4);  xmc1=xmc0(2:4)

  OFLWR "Computing the CrossSection for state ",istate,"wfn",imc; CFL

#ifdef REALGO
  MemVal = 1.25d5
#else
  MemVal = 6.25d4
#endif
  call openfile()
  write(mpifileptr,'(A30,F9.3,A3)') " Guess at necessary memory is ",2d0*real((nt+1)*2*spfsize*numr,8)/MemVal," MB"
  if(mem.eq.0) then
     write(mpifileptr,*) "g(tau) will be computed with all of psi in core"
     BatchSize=nt+1
  else
     MemTot=real(mem,8)    
     write(mpifileptr,*) "g(tau) will be computed with all psi being read in batches"
     write(mpifileptr,'(A33,F9.3,A3)') " Desired amount of memory to use ",MemTot," MB"
     
     BatchSize=floor(MemTot * MemVal / (2d0*real(2*spfsize*numr,8))) !! for each time pair we need both spun orbitals in each r
     if(BatchSize.lt.1) then
        write(mpifileptr,*) "Tiny amount of memory or huge wavefunction, Batchsize is 1" 
        BatchSize=1
     else if(BatchSize.ge.nt+1) then
        write(mpifileptr,*) "Hooray, there is enough memory, doing it all in core" 
        BatchSize=nt+1
     else
        write(mpifileptr,*) "Batchsize is ",BatchSize,"/",(nt+1)
     endif
  endif

  allocate(gtau(0:nt,FluxNBins),ketmo(spfsize,numr,2,BatchSize),bramo(spfsize,numr,2,BatchSize))
  NBat=ceiling(real(nt+1)/real(BatchSize))
  ketreadsize=0;  brareadsize=0

  inquire (iolength=tlen) ketmo(:,1,:,1)
  write(mpifileptr,*) "Projected 1e- function record length is ",tlen;  call closefile()
!! set up the sum of gtau(0) for debugging monitoring purposes

  open(454,file="Dat/KVLsum."//xstate1//"_"//xmc1//".dat",status="unknown")
  write(454,*) "#KVL flux sum: itime, time, flux sum";  write(454,*)

!! lets go through this loop
!!  open(1003,file="proj.flux.wfn.bin",status="unknown",form="unformatted",access="direct",recl=tlen)

  open(1003,file=projfluxfile,status="unknown",form="unformatted",access="direct",recl=tlen)

  gtau=0d0

!! lets do this double time integral double batched monster loop here. 
!! as can be seen in total flux eveything but biortho is fast as hell

  do ketbat=1,NBat
     OFLWR "Reading ket batch ", ketbat, " of ", NBat," for state ",istate; CFL
     ketreadsize=min(BatchSize,nt+1-(ketbat-1)*BatchSize)
     
!! read the orbital |\psi(t)>
     if(myrank.eq.1) then
        do i=1,ketreadsize !! loop over times in this batch
!        do r=1,numr !! loop over all the r's for this time
!          k = (istate-1)*(nt+1)*numr + ((ketbat-1)*BatchSize+(i-1))*numr + r  

           do ir=1,numr !! loop over all the r's for this time           
              k = (imc-1)*nstate*(nt+1)*numr + (istate-1)*(nt+1)*numr + ((ketbat-1)*BatchSize+i-1)*numr + ir
              read(1003,rec=k) ketmo(:,ir,:,i)
           enddo
        enddo
     endif

     call mympibcast(ketmo(:,:,:,1:ketreadsize),1,spfsize*numr*2*ketreadsize)


!! change it to \hat{F}|\psi(t)>, the onee part
     do i=1,ketreadsize
        do k=1,2 !! loop over the spins...
           call projeflux_op_onee(ketmo(:,:,k,i))
        enddo
     enddo
!! begin the bra batch read loop
     do brabat=1,ketbat
        OFLWR "Reading bra batch ", brabat, " of ", ketbat," for state ",istate; CFL
        brareadsize=min(BatchSize,nt+1-(brabat-1)*BatchSize)
        if(myrank.eq.1) then
           do i=1,brareadsize
              do ir=1,numr 
!            k = (istate-1)*(nt+1)*numr + ((brabat-1)*BatchSize+(i-1))*numr + r 
             
                 k = (imc-1)*nstate*(nt+1)*numr + (istate-1)*(nt+1)*numr + ((brabat-1)*BatchSize+i-1)*numr + ir
                 read(1003,rec=k) bramo(:,ir,:,i)
              enddo
           enddo
        endif

     call mympibcast(bramo(:,:,:,1:brareadsize),1,spfsize*numr*2*brareadsize)

!! loop over the specific ket indices & over all previous times & evaluate the actual g(tau) expression           
        do kettime=1,ketreadsize
           curtime=(ketbat-1)*BatchSize+kettime-1 
           whichbin=1
           do k=2,FluxNBins
              if(curtime.le.bintimes(k)) whichbin=k
           enddo
           if(brabat.lt.ketbat) then
              bratop=brareadsize
           else
              bratop=kettime
           endif
           do bratime=1,bratop
              tau=curtime-((brabat-1)*BatchSize+bratime-1)
              gtau(tau,1:whichbin) = gtau(tau,1:whichbin) + dot(bramo(:,:,:,bratime),ketmo(:,:,:,kettime),2*spfsize*numr) * dt
           enddo
           
!! create the xsec as it should be now
           doFT=.false.
           if(brabat.eq.ketbat) then
              do k=1,FluxNBins
                 if(curtime.eq.bintimes(k)) doFT=.true.
              enddo
           endif
           if(doFT) then
              OFLWR "Taking the fourier transform of g(tau) to get cross section"; CFL
              
              if(myrank.eq.1) then
                 open(1004,file="Dat/xsec.proj."//xstate1//"_"//xmc1//".spi.dat",status="replace",action="readwrite",position="rewind")
                 write(1004,*);write(1004,*) "# Omega; pulse ft; projected flux at t= ... "
                 write(1004,'(A8, A8, 100F36.5)') "#  ", " ", bintimes(1:FluxNBins)*dt;            write(1004,*)
              endif

              do k=1,nEFlux
                 wfi=EFluxLo+(k-1)*dEFlux;                 Myftgtau=0d0
                 if(wfi.ge.0d0) then
                    Energy=wfi+eground;                    Myftgtau(:) = gtau(0,:) * dt
                    do i=1,curtime
                       if(FluxSineOpt.eq.0) then
                          ffac=1d0 * dt
                       else
                          ffac=(cos(real(i)/real(nt)*3.14159265d0/2d0)**FluxSineOpt) * dt
                       endif
                       Myftgtau(:) = Myftgtau(:) + exp((0d0,1d0)*Energy*i*dt) * gtau(i,:) * ffac
                       Myftgtau(:) = Myftgtau(:) + exp((0d0,1d0)*Energy*(-i)*dt) * ALLCON(gtau(i,:)) * ffac
                    enddo
                 endif

                 if(myrank.eq.1) then
                    write(1004,'(F8.4,100E18.6)') wfi, xsec(k), xsec(k)*Myftgtau(1:FluxNBins)* cgfac
                 endif
              enddo
              if(myrank.eq.1) then
                 close(1004)
              endif




              allocate(ftgtau(-curtime:curtime),ftgtausave(-curtime:curtime), wsave(20*(2*curtime+1)+30))
              ftgtau(:)=0d0; ftgtausave(:)=0d0
              
              do i=0,curtime
                 ftgtau(i) = ALLCON(gtau(i,1))   * cos(real(i,8)/real(curtime,8) * piover2) * exp((0.d0,-1.d0)*ALLCON(ceground)*par_timestep*FluxInterval*FluxSkipMult*i)
              enddo
              do i=0,curtime
                 ftgtau(-i) = (gtau(i,1))  * cos(real(i,8)/real(curtime,8) * piover2) * exp((0.d0,1.d0)*ceground*par_timestep*FluxInterval*FluxSkipMult*i)
              enddo

              if (myrank.eq.1) then
                 open(171,file="Dat/myGTau."//xstate1//"_"//xmc1//".Dat"  ,status="unknown");          write(171,*) "#   ", curtime
                 do i=-curtime,curtime
                    write(171,'(F18.12, T22, 400E20.8)')  i*par_timestep*FluxInterval*FluxSkipMult, ftgtau(i)
                 enddo
                 close(171)
              endif

              call zffti(2*curtime+1,wsave);          call zfftf(2*curtime+1,ftgtau(-curtime),wsave)


              ftgtau(-curtime:curtime)=ftgtau(-curtime:curtime)*par_timestep*FluxInterval*FluxSkipMult
              
              !! NEW 0603
              do i=-curtime,curtime
                 ftgtau(i)=ftgtau(i)*exp((0.d0,1.d0)*(curtime+i)*curtime*2*pi/real(2*curtime+1))
              enddo
              do i=1,curtime
                 ftgtausave(i-curtime-1) = ftgtau(i)
              enddo
              do i=-curtime,0
                 ftgtausave(curtime+i)=ftgtau(i)
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

                 open(171,file="Dat/FTGTau."//xstate1//"_"//xmc1//".Dat"  ,status="unknown");          write(171,*) "#   ", curtime
!!                 open(171,file="Dat/FTGtau.Dat", status="unknown")
                 do i=-curtime,curtime
                    myft=1d0
                    if (flag.eq.1) then
                       myft=abs(tdpotft(i*Estep))**2
                    endif
!             write(171,'(F18.12, T22, 400E20.8)')  i*Estep, abs(ftgtausave(i)), ftgtausave(i), ftgtausave(i)/myft, myft
                    write(171,'(F18.12, T22, 400E20.8)')  i*Estep,  ftgtausave(i), ftgtausave(i)/myft, myft
                 enddo
                 close(171)
              endif
              deallocate(ftgtau,ftgtausave, wsave)

           endif  !! doft

!! only do this after we are sure we've gone through every bra
           if(brabat.eq.ketbat) then
              if (myrank.eq.1) then
                 open(454, file="Dat/KVLsum."//xstate1//"_"//xmc1//".dat", status="old", position="append")
                 write(454,'(I5,100F18.12)') curtime, curtime*dt, gtau(0,1);    
                 close(454)
              endif
           endif


!! write out times
!        if(brabat.eq.ketbat.and.(mod(curtime,20).eq.0.or.curtime.eq.nt)) then
!          write(454,'(I5,100F18.12)') curtime, curtime*dt, gtau(0,1)
!          call openfile
!          write(mpifileptr,'(A50,F10.4)') " Timing Statistics: Computing gtau and xsec at T= ",curtime*dt
!          write(mpifileptr,'(100A10)') "Times: ", "All", "IO", "Walks", "Biorth", "1e Bld", "1e Op", "Fluxeval", "!FT gtau"
!          write(mpifileptr,'(A10,100I10)') " ", times(1:8)/100
!          call closefile
!        endif
        enddo
     enddo
  enddo
  close(1003);  close(454);  deallocate(gtau,ketmo,bramo)

end subroutine projeflux_double_time_int


!! get the contraction of the flux operator (iH-H^\dag) with our current set of orbitals F*inspfs=outspfs
subroutine projeflux_op_onee(inspfs)
  use parameters
  implicit none
  DATATYPE :: inspfs(spfsize,numr),outspfs(spfsize,numr),ttempspfs(spfsize)
  integer :: r

  outspfs=0d0
  do r=1,numr
     call mult_imke(inspfs(:,r),ttempspfs(:))
     outspfs(:,r)=outspfs(:,r)+ttempspfs(:);     outspfs(:,r) = outspfs(:,r) / bondpoints(r)
  enddo
!! scale correctly and clear memory
  inspfs=outspfs*(-2d0)

end subroutine projeflux_op_onee


module projbiomod
  use biorthotypemod
  implicit none
  type(biorthotype),target :: projbiovar
end module projbiomod

subroutine projeflux_single(mem)
  use projbiomod
  use biorthomod
  use parameters
  use walkmod
  use projefluxmod
  use mpimod
  implicit none
!! necessary working variables
  integer :: mem,tau, i,k,nt ,ir,tndof,tnspf,nstate,tnumr,istate,myiostat,ierr
  integer :: spfcomplex, acomplex, tdims(3),ttndof,ttnumconfig,imc,ii
  real*8 :: dt
  real*8, allocatable :: xsec(:)
  integer, allocatable :: bintimes(:)
  DATATYPE, allocatable :: readmo(:,:),readavec(:,:,:),mobio(:,:),abio(:,:)
!! mcscf specific read variables
  DATATYPE,target :: smo(nspf,nspf)

  OFLWR ;  WRFL   "   *** DOING PROJECTED FLUX. ***    ";  WRFL; CFL

!! read in the data from mcscf for our target cation state

  open(909,file="Bin/cation.spfs.bin",status="unknown",form="unformatted")
  open(910,file="Bin/cation.avector.bin",status="unknown",form="unformatted")

  call avector_header_read(910,nstate,tndof,tnumr,tnumconfig,targetrestrictflag,targetms,targetspinproject,targetspinval,acomplex,ierr)
  if (ierr.ne.0) then
     OFLWR "Avector header read error; redone 09/29/2014; recompute vector on disk."; CFLST
  endif

  call spf_header_read(909,tdims,tnspf,spfcomplex)

!! have to project on BO wfns.  Otherwise doesn't make sense in prolate.  Not supported anymore, will support
!!   again, mcscf mode w/many r's on file

  if (tnumr.ne.1) then
     OFLWR "projecting only one r value... input numr>1... numr>1 proj not supported yet"; CFLST
  endif
  if (tnspf.gt.nspf) then
     OFLWR "ERROR, for now can't do more orbs in projection than in calculation",tnspf,nspf; CFLST
  endif
  if (tndof.ne.ndof-2 ) then
     OFLWR "Vectors read are not n-1 electron functions"; CFLST
  endif

  cgflag=0
  if (restrictflag.eq.0) then
     OFLWR " WARNING: propagated state is not spin projection restricted.  No CG algebra."; CFL
  else
     if (targetrestrictflag.eq.0) then
        OFLWR "WARNING: N-1 electron state is not spin projection restricted (restrictflag=0).  No CG algebra.";CFL
     else
        if (abs(targetms-restrictms).ne.1) then
           OFLWR "Targetms should differ from restrictms by 1.  Cation and neutral have no dyson orbital"
           WRFL  "    Targetms=", targetms, " restrictms=",restrictms; CFLST
        endif
        if (allspinproject.eq.0) then
           OFLWR "WARNING: Wave function is not spin restricted (allspinproject=0).  No CG algebra.";CFL
        else
           if (targetspinproject.eq.0) then
              OFLWR "WARNING: N-1 electron state is not spin restricted (allspinproject=0).  No CG algebra.";CFL
           else
              if (abs(targetspinval-spinrestrictval).gt.1) then
                 OFLWR "Spin value of wave function and N-1 electron state differ by more than 1/2.",targetspinval,spinrestrictval; CFL !!ST
                 OFLWR "TEMP CONTINUE"; CFL
                 OFLWR "TEMP CONTINUE"; CFL
                 OFLWR "TEMP CONTINUE"; CFL
              else
                 OFLWR "Spin adapted wave functions for propagation and projection. doing CG algebra.";CFL
                 cgflag=1
              endif
           endif
        endif
     endif
  endif

  allocate(tmo(spfsize,nspf,numr),ta(tnumconfig,nstate,numr))

  OFLWR "Reading", nstate," Born-Oppenheimer states."; CFL

!! no should be ok  call noparorbsupport("in projeflux_single")

  tmo=0d0;  ta=0d0

  call simple_load_avectors(910,acomplex,ta(:,:,1),ndof-2,1,tnumconfig,nstate)

  call spf_read0(909,nspf,spfdims,tnspf,tdims,spfcomplex,spfdimtype,tmo(:,:,1),(/0,0,0/))


  do i=tnspf+1,nspf
     tmo(:,i,1)=0d0
     call staticvector(tmo(:,i,1),spfsize)
     call gramschmidt(spfsize,i-1,spfsize,tmo(:,:,1),tmo(:,i,1),.false.)
  enddo



!! only one r point on file supported..temporary
  do ir=1,numr
     tmo(:,:,ir)=tmo(:,:,1)
     ta(:,:,ir)=ta(:,:,1)
  enddo

  do ii=1,myrank
     call mpibarrier()
  enddo
  open(676,file="WALKS/cation.configlist.BIN",status="old",form="unformatted",iostat=myiostat)
  if (myiostat.ne.0) then
     OFLWR "iostat ",myiostat," in open of cation.configlist.BIN"; CFLST
  endif
  call singlewalkheaderread(676,ttnumconfig,ttndof,tmaxsinglewalks)

  if (ttnumconfig.ne.tnumconfig.or.tndof.ne.ttndof) then
     OFLWR "TTTYQ ERROR"; CFLST
  endif

  allocate(tconfiglist(tndof,tnumconfig),tnumsinglewalks(tnumconfig),tsinglewalk(tmaxsinglewalks,tnumconfig))
  allocate(tsinglewalkopspf(1:2,tmaxsinglewalks,tnumconfig),tsinglewalkdirphase(tmaxsinglewalks,tnumconfig))

  call singlewalkread(676,tnumconfig,tndof,tmaxsinglewalks,tconfiglist,tnumsinglewalks, &
       tsinglewalk,tsinglewalkopspf,tsinglewalkdirphase)
  close(676)

  do ii=myrank,nprocs
     call mpibarrier()
  enddo

!! do the walks from the target state into our final state

  call projeflux_singlewalks()

!! allocate all necessary extra memory and io params to do this looping business

!!  dt=real(FluxInterval*FluxSkipMult,8)*par_timestep;  nt=floor(final time/dt)

  dt=real(FluxInterval*FluxSkipMult,8)*par_timestep;  nt=floor(real(numpropsteps,8)/fluxinterval/fluxskipmult)

  allocate(mobio(spfsize,nspf),abio(numconfig,mcscfnum),readmo(spfsize,nspf),readavec(numconfig,numr,mcscfnum))
  inquire (iolength=i) readmo
  open(1001,file=fluxmofile,status="unknown",form="unformatted",access="direct",recl=i)
  inquire (iolength=i) readavec
  open(1002,file=fluxafile,status="unknown",form="unformatted",access="direct",recl=i)

!! do the loop over all time, ALL TIME and now in parallel-o-vision

  do tau=0,nt    
!! read in this time's wavefucntion

     read(1001,rec=FluxSkipMult*tau+1) readmo(:,:);     read(1002,rec=FluxSkipMult*tau+1) readavec(:,:,:)

!! do biortho and construct the single particle function

     do ir=1,numr

        abio(:,:)=readavec(:,ir,:)
        call bioset(projbiovar,smo,mcscfnum); 
        call biortho(readmo,tmo(:,:,ir),mobio,abio(:,:),projbiovar)

        do imc=1,mcscfnum
           do istate=1,nstate
!!      i=(istate-1)*(nt+1)*numr + tau*numr +r 

              i=(imc-1)*nstate*(nt+1)*numr + (istate-1)*(nt+1)*numr + tau*numr + ir
              call projeflux_doproj(ta(:,istate,ir),abio(:,imc),mobio(:,:),i)

           enddo
        enddo
     enddo

!! write out times cause we're bored
!    if(mod(tau,100).eq.0.or.tau.eq.nt) then
!      call openfile
!      write(mpifileptr,'(A52,F10.4)') " Timing Statistics: producing 1e- proj wfn as of T= ",tau*dt
!      write(mpifileptr,'(100A10)') "Times: ", "All", "IO", "Walks", "Biorth", "1e Bld", "1e Op", "Fluxeval", "FT g!tau"
!      write(mpifileptr,'(A10,100I10)') " ", times(1:8)/100
!    endif
  enddo
!! clean up
  close(1001);  close(1002)
  deallocate(mobio,abio,readmo,readavec,ta,tconfiglist,tmo)
  deallocate(tnumsinglewalks,tsinglewalk,tsinglewalkopspf,tsinglewalkdirphase)
  deallocate(numpwalk1,pwalk1,pspf1,pphase1) !!,num pw alk2,pwal k2,ps pf2,ppha se2)


!! do all the prep work and then do the double time integral
!! get FT of pulse and constants
  allocate(bintimes(FluxNBins),xsec(nEFlux))
  call getFTpulse(xsec)

!! figure out the bins
  tau=floor(real(nt,8)/real(FluxNBins,8))
  k=0
  do i=FluxNBins,2,-1
    k=k+tau
    bintimes(i)=k
  enddo
  bintimes(1)=nt

!! do the double time integral
  do imc=1,mcscfnum
  do istate=1,nstate
    call projeflux_double_time_int(bintimes,xsec,mem,istate,nstate,imc,nt,dt)
  enddo
  enddo
 
  deallocate(bintimes,xsec)
  OFLWR "Cross Section acquired, cleaning and closing";CFLST
end subroutine projeflux_single


