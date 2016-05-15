

!! PROJECTED FLUX (partial photoionization ACTION 17)


#include "Definitions.INC"

module projefluxmod !! needed for cation walks and bi-orthonormalization
!  implicit none
!  type PPTYPE
  integer :: tnumconfig
  integer, allocatable :: tconfiglist(:,:),pphase1(:,:) ,pspf1(:,:,:) ,numpwalk1(:) ,pwalk1(:,:)
  integer :: maxpwalk1 
!  end type PPTYPE
!  type(PPYPE),allocatable :: ppp(:)    !! numcatfiles
  integer :: nstate
  integer, allocatable :: eachstate(:) !! numcatfiles
end module projefluxmod


!! here we are doing walks from our BO target state to our wavefunction \Psi(t)
subroutine projeflux_singlewalks()
  use fileptrmod
  use basis_parameters !! numpart, num2part
  use aarrmod
  use projefluxmod
  use configmod
  use mpimod
  use configsubmod
  implicit none
  integer :: iconfig,jconfig,idof,iindex,iwalk,flag,dirphase
  integer :: tempconfig(num2part),temporb(2)

!! get the number of single walks from the target N-1 e- state to our regular N e- state
  OFLWR "Getting cation single walks"; CFL

  allocate(numpwalk1(tnumconfig));    numpwalk1(:)=0

  do iconfig=1,tnumconfig
     iwalk=0
     do iindex=1,2*www%nspf

       if (www%holeflag.eq.0) then
          tempconfig(3:num2part)=tconfiglist(:,iconfig)
          temporb=aarr(iindex)
          tempconfig(1:2)=temporb(:);        
          flag=0
          do idof=2,numpart 
             if(iind(tempconfig(idof*2-1:idof*2)).eq.iindex) flag=1
          enddo
       else   !! holeflag
          flag=1
          do idof=1,numpart+1
             if(iind(tconfiglist(idof*2-1:idof*2,iconfig)).eq.iindex) then
                flag=0
                tempconfig(1:(idof-1)*2) = tconfiglist(1:(idof-1)*2,iconfig)
                tempconfig(idof*2-1:num2part) = tconfiglist(idof*2+1:num2part+2,iconfig)
                temporb=aarr(iindex)
                exit
             endif
          enddo
       endif   !! holeflag
       if(flag.eq.0) then
          dirphase=reorder(tempconfig,www%numpart)
          if(allowedconfig0(www,tempconfig,www%dflevel)) then
             jconfig=getconfiguration(tempconfig,www)
             if (jconfig.ge.www%botconfig.and.jconfig.le.www%topconfig) then
                iwalk=iwalk+1
             endif
          endif
       endif
    enddo
    numpwalk1(iconfig)=iwalk
 enddo

!! figure out the maximum number of single target walks
  maxpwalk1=0
  do iconfig=1,tnumconfig
    if(maxpwalk1.lt.numpwalk1(iconfig)) maxpwalk1=numpwalk1(iconfig)
  enddo
  OFLWR "Max # single walks from cation state on this processor is ",maxpwalk1;CFL
  call mpibarrier()

  allocate(pphase1(maxpwalk1,tnumconfig),pwalk1(maxpwalk1,tnumconfig),&
       pspf1(2,maxpwalk1,tnumconfig))
  if (maxpwalk1.gt.0) then
     pwalk1=0;  pspf1=0;  pphase1=0
  endif

  do iconfig=1,tnumconfig
     iwalk=0
     do iindex=1,2*www%nspf
        dirphase=(-9999);        tempconfig(:)=(-9999);        temporb(:)=(-9999)
        if (www%holeflag.eq.0) then
           flag=0
           do idof=1,numpart-1
              if(iind(tconfiglist(idof*2-1:idof*2,iconfig)).eq.iindex) then
                 flag=1
                 exit
              endif
           enddo
           if (flag.eq.0) then
              temporb=aarr(iindex)
              tempconfig(1:2)=temporb(:)
              tempconfig(3:num2part)=tconfiglist(:,iconfig)
              dirphase=reorder(tempconfig,www%numpart)
           endif
        else    !! holeflag
           flag=1
           do idof=1,numpart+1
              if(iind(tconfiglist(idof*2-1:idof*2,iconfig)).eq.iindex) then
                 flag=0
                 temporb=aarr(iindex)
                 tempconfig(1:(idof-1)*2) = tconfiglist(1:(idof-1)*2,iconfig)
                 tempconfig(idof*2-1:num2part) = tconfiglist(idof*2+1:num2part+2,iconfig)
                 dirphase=(-1)**idof
                 exit
              endif
           enddo
        endif  !! holeflag
        if(flag.eq.0) then
           if(allowedconfig0(www,tempconfig,www%dflevel)) then
              jconfig=getconfiguration(tempconfig,www)
              if (jconfig.ge.www%botconfig.and.jconfig.le.www%topconfig) then
                 iwalk=iwalk+1;          pwalk1(iwalk,iconfig)=jconfig
                 pspf1(:,iwalk,iconfig)=temporb;          pphase1(iwalk,iconfig)=dirphase
              endif
           endif
        endif
     enddo
!! check and make sure no bad walk error
     if(numpwalk1(iconfig).ne.iwalk) then
        OFLWR "TARGET SINGLE WALK ERROR";CFLST
     endif
  enddo
  
  call mpibarrier()
  OFLWR "DONE getting proj single walks"; CFL

end subroutine projeflux_singlewalks


!! construct the one electron functions
subroutine projeflux_doproj(cata,neuta,mo,offset,cgfac)
  use parameters     !! nspf, projfluxfile etc.
  use projefluxmod
  use mpimod
  use mpisubmod
  implicit none
  integer,intent(in) :: offset 
  DATATYPE,intent(in) :: cata(tnumconfig),&
       neuta(first_config:last_config),mo(spfsize,nspf)
  real*8, intent(in) :: cgfac
  DATATYPE :: projwfn(spfsize,2),  projcoefs(nspf,2)
  integer :: jconfig,iwalk,iconfig,ispf,ispin,iphase,mylength,myiostat
  DATATYPE,allocatable:: bigprojwfn(:,:)

!! make the single electron wfn

  projcoefs(:,:)=0d0
  do jconfig=1,tnumconfig
    do iwalk=1,numpwalk1(jconfig)
      iconfig=pwalk1(iwalk,jconfig);      ispf=pspf1(1,iwalk,jconfig)
      ispin=pspf1(2,iwalk,jconfig);      iphase=pphase1(iwalk,jconfig)

      projcoefs(ispf,ispin)=projcoefs(ispf,ispin) + &
           CONJUGATE(cata(jconfig)) * neuta(iconfig) * iphase
    enddo
  enddo

  call mympireduce(projcoefs,nspf*2)

  projwfn(:,:)=0d0
  do ispin=1,2
     do ispf=1,nspf

!! sqrt(cgfac) to recover cgfac with bra and ket factor.  cgfac positive

        projwfn(:,ispin) = projwfn(:,ispin) + mo(:,ispf) * projcoefs(ispf,ispin) * sqrt(cgfac)
     enddo
  enddo

  if (myrank.eq.1) then
     if (parorbsplit.eq.3) then
        allocate(bigprojwfn(spfsize*nprocs,2));        bigprojwfn(:,:)=0d0
     else
        allocate(bigprojwfn(spfsize,2));        bigprojwfn(:,:)=projwfn(:,:)
     endif
  else
     allocate(bigprojwfn(1,2));   bigprojwfn(:,:)=0d0
  endif

  if (parorbsplit.eq.3) then
     do ispin=1,2
        call splitgatherv(projwfn(:,ispin),bigprojwfn(:,ispin),.false.)
     enddo
  endif

  if (myrank.eq.1) then
     inquire (iolength=mylength) bigprojwfn
     open(1003,file=projfluxfile,status="unknown",form="unformatted",&
          access="direct",recl=mylength,iostat=myiostat)
     call checkiostat(myiostat,"opening projfluxfile")
     write(1003,rec=offset,iostat=myiostat) bigprojwfn(:,:) 
     call checkiostat(myiostat,"writing projfluxfile")
     close(1003)
  endif

  deallocate(bigprojwfn)

end subroutine projeflux_doproj


!! do the double time integral piece    08-2015 now looping over istate and imc here

subroutine projeflux_double_time_int(mem,nstate,nt)
  use parameters    !! par_timestep and others
  use mpimod
  use mpisubmod
  use pulsesubmod
  implicit none
  integer,intent(in) :: mem,nstate,nt
  integer :: i,k,tlen,istate,curtime,tau,ir ,imc,ioffset, NUMANGLES, NUMERAD, il
  integer :: BatchSize,NBat,ketreadsize,brareadsize,ketbat,brabat,kettime,bratime,&
       bratop,getlen,myiostat
  real*8 :: MemTot,MemVal,wfi,estep,myfac,windowfunct,MemNum,dt
  DATATYPE, allocatable :: bramo(:,:,:,:),ketmo(:,:,:,:),gtau(:,:,:),gtau_ad(:,:,:,:),ketop(:,:,:,:),&
       read_bramo(:,:,:), read_ketmo(:,:,:), deweighted_bramo(:,:,:), ketmo_ad(:,:,:,:,:),&
       ketop_ad(:,:,:,:,:), gtaunow(:,:,:), gtaunow_ad(:,:,:,:), gtausum(:,:,:), gtausum_ad(:,:,:,:)
  complex*16, allocatable :: ftgtau(:),pulseft(:,:), total(:),ftgtau_ad(:,:),total_ad(:,:)
  complex*16 :: ftgtausum
  real*8, allocatable :: pulseftsq(:)
  DATATYPE :: pots1(3)
  character (len=4) :: xstate0,xmc0
  character (len=3) :: xstate1,xmc1

  dt=real(FluxInterval*FluxSkipMult,8)*par_timestep

  if (ceground.eq.(0d0,0d0)) then
     OFLWR "Eground is ZERO.  Are you sure?  If want zero just make it small."
     WRFL  "  Otherwise need eground: initial state energy."; CFLST
  endif

!! determine if we should do batching or not
!! 250,000 words/MB, real*8 2words/#, complex*16 4words/#

#ifdef REALGO
  MemVal = 1.25d5
#else
  MemVal = 6.25d4
#endif

!! rank 1 now only holds one entire orbital vector for parorbsplit.eq.3.
!! read_bramo and read_ketmo not dimensioned with BatchSize any more.
!! MemNum based on size of arrays bramo and ketmo instead.

!!$  if (parorbsplit.eq.3) then
!!$     MemNum=4*spfsize*numr*nprocs    !! 4 for bra and ket and both spun orbitals
!!$  else

  MemNum=4*spfsize*numr

!!$  endif

  call openfile()
  write(mpifileptr,'(A30,F9.3,A3)') " Guess at necessary memory is ",&
       (nt+1)*MemNum/MemVal," MB"
  if(mem.eq.0) then
     write(mpifileptr,*) "g(tau) will be computed with all of psi in core"
     BatchSize=nt+1
  else
     MemTot=real(mem,8)    
     write(mpifileptr,*) "g(tau) will be computed with all psi being read in batches"
     write(mpifileptr,'(A33,F9.3,A3)') " Desired amount of memory to use ",MemTot," MB"
     
     BatchSize=floor(MemTot * MemVal / MemNum) 
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
  call closefile()

  allocate(gtau(0:nt,nstate,mcscfnum),ketmo(spfsize,numr,2,BatchSize),ketop(spfsize,numr,2,BatchSize),&
       bramo(spfsize,numr,2,BatchSize),gtaunow(0:nt,nstate,mcscfnum),gtausum(0:nt,nstate,mcscfnum))
  gtau=0; ketmo=0; ketop=0; bramo=0; gtaunow=0; gtausum=0d0

  NUMANGLES=0; NUMERAD=0

  if (angularflag.ne.0) then    !! FULLY DIFFERENTIAL PHOTOIONIZATION (function of angle and energy)
     if (spfdimtype(2).ne.2) then
        OFLWR "For angularflag, fully differential photoionizatoin, must have spfdimtype(2).eq.2"
        WRFL "     (atom/diatom only)"; CFLST
     endif

     NUMANGLES=spfdims(2)
     NUMERAD=spfdims(1)

     allocate(gtau_ad(0:nt,nstate,mcscfnum,NUMANGLES),deweighted_bramo(spfsize,numr,2),&
          ketmo_ad(spfsize,numr,2,batchsize,NUMANGLES), ketop_ad(spfsize,numr,2,batchsize,NUMANGLES),&
          gtaunow_ad(0:nt,nstate,mcscfnum,NUMANGLES),gtausum_ad(0:nt,nstate,mcscfnum,NUMANGLES))
  else
     allocate(gtau_ad(0:0,1,1,1),deweighted_bramo(1,1,2),&
          ketmo_ad(1,1,2,1,1), ketop_ad(1,1,2,1,1),&
          gtaunow_ad(0:0,1,1,1),gtausum_ad(0:0,1,1,1))
  endif


  gtau_ad=0;   deweighted_bramo=0;  ketmo_ad=0; ketop_ad=0; gtaunow_ad=0; gtausum_ad=0

  if (myrank.eq.1) then
     if (parorbsplit.eq.3) then
        allocate(read_ketmo(spfsize*nprocs,numr,2),&
             read_bramo(spfsize*nprocs,numr,2))
     else
        allocate(read_ketmo(spfsize,numr,2),&
             read_bramo(spfsize,numr,2))
     endif
  else
     allocate(read_ketmo(1,numr,2),read_bramo(1,numr,2))
  endif
  read_ketmo=0; read_bramo=0

  NBat=ceiling(real(nt+1)/real(BatchSize))
  ketreadsize=0;  brareadsize=0

  if (myrank.eq.1) then
     inquire (iolength=tlen) read_ketmo(:,1,:)
  endif
  call mympiibcastone(tlen,1)

  OFLWR "Projected 1e- function record length is ",tlen;  CFL

  gtau(:,:,:)=0d0
  gtau_ad(:,:,:,:)=0d0

!! looping here now 08-2015

  do imc=1,mcscfnum

     do istate=1,nstate

        OFLWR "Computing the CrossSection for state ",istate,"wfn",imc; CFL

        if (myrank.eq.1) then
           open(1003,file=projfluxfile,status="unknown",form="unformatted",&
                access="direct",recl=tlen,iostat=myiostat)
           call checkiostat(myiostat,"opening "//projfluxfile)
        endif

!! lets do this double time integral double batched monster loop here. 
!! as can be seen in total flux eveything but biortho is fast as hell

        do ketbat=1,NBat
           OFLWR "Reading ket batch ", ketbat, " of ", NBat," for state ",istate; CFL
           ketreadsize=min(BatchSize,nt+1-(ketbat-1)*BatchSize)
     
!! read the orbital |\psi(t)>

           do i=1,ketreadsize !! loop over times in this batch
              if(myrank.eq.1) then
                 do ir=1,numr !! loop over all the r's for this time           

  ioffset=(istate-1)*mcscfnum*(nt+1)*numr + (imc-1)*(nt+1)*numr + ((ketbat-1)*BatchSize+i-1)*numr + ir

                    read(1003,rec=ioffset,iostat=myiostat) read_ketmo(:,ir,:)
                    call checkiostat(myiostat,"reading ketmo")
                 enddo
              endif
              if (parorbsplit.ne.3) then
                 if (myrank.eq.1) then
                    ketmo(:,:,:,i)=read_ketmo(:,:,:)
                 endif
                 call mympibcast(ketmo(:,:,:,i),1,spfsize*numr*2)
              else
                 do k=1,2
                    do ir=1,numr
                       call splitscatterv(read_ketmo(:,ir,k),ketmo(:,ir,k,i))
                    enddo
                 enddo
              endif
           enddo       !! do i=1,ketreadsize

           if (fluxoptype.eq.0) then
              OFLWR "PROGRAM ME FLUXOPTYPE 0 PROJEFLUX"; CFLST
           endif

           do i=1,ketreadsize
              do k=1,2 !! loop over the spins...
                 call projeflux_op_onee(ketmo(:,:,k,i), ketop(:,:,k,i))
              enddo
!! for integral dt
              curtime=(ketbat-1)*BatchSize+i-1
              gtaunow(curtime,istate,imc) = hermdot(ketmo(:,:,:,i),ketop(:,:,:,i),2*spfsize*numr)
           enddo
           if (parorbsplit.eq.3) then
              call mympireduce(gtaunow(:,istate,imc),ketreadsize)
           endif

           if (angularflag.ne.0) then
              do il=1,NUMANGLES
                 ketmo_ad(:,:,:,1:ketreadsize,il) = ketmo(:,:,:,1:ketreadsize)
                 call only_one_angle(2*numr*ketreadsize,il,ketmo_ad(:,:,:,:,il))
                 do i=1,ketreadsize
                    do k=1,2 !! loop over the spins...
                       call projeflux_op_onee(ketmo_ad(:,:,k,i,il), ketop_ad(:,:,k,i,il))
                    enddo
!! for integral dt
                    curtime=(ketbat-1)*BatchSize+i-1
                    gtaunow_ad(curtime,istate,imc,il) = hermdot(ketmo_ad(:,:,:,i,il),ketop_ad(:,:,:,i,il),2*spfsize*numr)
                 enddo
                 if (parorbsplit.eq.3) then
                    call mympireduce(gtaunow_ad(:,istate,imc,il),ketreadsize)
                 endif
              enddo
           endif

!! begin the bra batch read loop
           do brabat=1,ketbat
              OFLWR "Reading bra batch ", brabat, " of ", ketbat," for state ",istate; CFL
              brareadsize=min(BatchSize,nt+1-(brabat-1)*BatchSize)
              do i=1,brareadsize
                 if(myrank.eq.1) then
                    do ir=1,numr 

   ioffset = (istate-1)*mcscfnum*(nt+1)*numr + (imc-1)*(nt+1)*numr + ((brabat-1)*BatchSize+i-1)*numr + ir

                       read(1003,rec=ioffset,iostat=myiostat) read_bramo(:,ir,:)
                       call checkiostat(myiostat,"reading bramo")
                    enddo
                 endif
                 if (parorbsplit.ne.3) then
                    if (myrank.eq.1) then
                       bramo(:,:,:,i)=read_bramo(:,:,:)
                    endif
                    call mympibcast(bramo(:,:,:,i),1,spfsize*numr*2)
                 else
                    do k=1,2
                       do ir=1,numr
                          call splitscatterv(read_bramo(:,ir,k),bramo(:,ir,k,i))
                       enddo
                    enddo
                 endif
              enddo     !! do i=1,brareadsize

!! loop over the specific ket indices & over all previous 
!! times & evaluate the actual g(tau) expression           

              do kettime=1,ketreadsize

                 curtime=(ketbat-1)*BatchSize+kettime-1 

                 if(brabat.lt.ketbat) then
                    bratop=brareadsize
                 else
                    bratop=kettime
                 endif
                 do bratime=1,bratop
                    tau=curtime-((brabat-1)*BatchSize+bratime-1)
                    gtau(tau,istate,imc) = gtau(tau,istate,imc) + &
                         hermdot(bramo(:,:,:,bratime),ketop(:,:,:,kettime),2*spfsize*numr) * dt

!!      allocated   bramo(spfsize,numr,2,BatchSize)

                    if (angularflag.ne.0) then

!!$                       gtau_ad(tau,istate,imc,:) = gtau_ad(tau,istate,imc,:) + &
!!$                            myhermdots(bramo(:,:,:,bratime),ketop(:,:,:,kettime),2*spfsize*numr) * dt

                       do k=1,2
                          do ir=1,numr

!! with /2/pi, we have Mb per steradian column 4 angprojspifile... right?  column 5 is nonsense pmctdhf/chmctdhf
!! weights elecweights(:,:,:,3)=1d0.

                             deweighted_bramo(:,ir,k) = bramo(:,ir,k,bratime) / &
                                  RESHAPE(elecweights(:,:,:,2),(/spfsize/)) / 2 / pi
                          enddo
                       enddo
!!$                       gtau_ad(tau,istate,imc,:) = gtau_ad(tau,istate,imc,:) + &
!!$                            myhermdots(deweighted_bramo(:,:,:),ketop(:,:,:,kettime),2*spfsize*numr) * dt

!! better way, with symmetric flux operator.  total still equals integrated differential.
!! results (at threshold where there is oscillation and usually some negative values) actually 
!! look worse for H2 with this symmetrized operator, below, instead of what's above

                       gtau_ad(tau,istate,imc,:) = gtau_ad(tau,istate,imc,:) + 0.5d0 * &
                            myhermdots(deweighted_bramo(:,:,:),ketop(:,:,:,kettime),2*spfsize*numr) * dt

                       do il=1,NUMANGLES
                          gtau_ad(tau,istate,imc,il) = gtau_ad(tau,istate,imc,il) + 0.5d0 * &
                               hermdot(deweighted_bramo(:,:,:),ketop_ad(:,:,:,kettime,il),2*spfsize*numr) * dt
                       enddo

                    endif
                 enddo

              enddo !! do kettime
           enddo   !! do brabat
        enddo   !! do ketbat
        
        if (curtime.ne.nt) then
           OFLWR "DEEG CURTIME NT ERR", curtime, nt; CFLST
        endif

        if (myrank.eq.1) then
           close(1003)
        endif

     enddo
  enddo

  if (parorbsplit.eq.3) then
     call mympireduce(gtau(:,:,:), (nt+1)*nstate*mcscfnum)
     if (angularflag.ne.0) then
        call mympireduce(gtau_ad(:,:,:,:), (NUMANGLES)*(nt+1)*nstate*mcscfnum)
     endif
  endif

!! INTEGRAL DT

  gtausum(:,:,:)=0d0;  gtausum_ad(:,:,:,:)=0d0

  gtausum(0,:,:)=gtaunow(0,:,:)
  do i=1,nt
     gtausum(i,:,:)=gtausum(i-1,:,:) + gtaunow(i,:,:) * dt / 2d0   !! 2 looks correct
  enddo
  if (angularflag.ne.0) then
     gtausum_ad(0,:,:,:)=gtaunow_ad(0,:,:,:)
     do i=1,nt
        gtausum_ad(i,:,:,:)=gtausum_ad(i-1,:,:,:) + gtaunow_ad(i,:,:,:) * dt
     enddo
  endif

  if(myrank.eq.1) then

!! nevermind KVLsum.dat and gtau.dat
!! now integrals dt and domega with projfluxtsumfile, angprojfluxtsumfile in namelist &parinp
!! and additional column in projspifile

!! integral dt:

     do imc=1,mcscfnum
        do istate=1,nstate
           write(xstate0,'(I4)') istate+1000;  write(xmc0,'(I4)') imc+1000
           xstate1=xstate0(2:4);  xmc1=xmc0(2:4)

           open(171,file=projfluxtsumfile(1:getlen(projfluxtsumfile))//"_"//xstate1//"_"//xmc1//".dat", &
                status="unknown",iostat=myiostat)
           call checkiostat(myiostat,"opening projfluxtsumfile")
           write(171,*,iostat=myiostat) "#   ", curtime
           call checkiostat(myiostat,"writing projfluxtsumfile")
           do i=0,curtime
              write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat) &
                   i*dt,  gtausum(i,istate,imc), gtaunow(i,istate,imc)
           enddo
           call checkiostat(myiostat,"writing projfluxtsumfile")
           close(171)

           if (angularflag.ne.0) then
              open(171,file=angprojfluxtsumfile(1:getlen(angprojfluxtsumfile))//"_"//xstate1//"_"//xmc1//".dat", &
                   status="unknown",iostat=myiostat)
              call checkiostat(myiostat,"opening angprojfluxtsumfile")
              write(171,*,iostat=myiostat) "#   ", curtime
              call checkiostat(myiostat,"writing angprojfluxtsumfile")
              do i=0,curtime
                 write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat) &
                      i*dt,  gtausum_ad(i,istate,imc,:), gtaunow_ad(i,istate,imc,:)
              enddo
              call checkiostat(myiostat,"writing angprojfluxtsumfile")
              close(171)
           endif

        enddo !! istate
     enddo    !! imc
  endif       !! myrank

  call mpibarrier()

  OFLWR "Taking the fourier transform of g(tau) to get cross section at T= ",curtime*dt; CFL

  allocate(ftgtau(-curtime:curtime), pulseft(-curtime:curtime,3),&
       pulseftsq(-curtime:curtime),total(-curtime:curtime))
  ftgtau(:)=0d0; pulseft(:,:)=0d0; pulseftsq(:)=0d0; total(:)=0

  if (angularflag.ne.0) then
     allocate(ftgtau_ad(-curtime:curtime,NUMANGLES), total_ad(-curtime:curtime,NUMANGLES))
  else
     allocate(ftgtau_ad(0:0,NUMANGLES), total_ad(0:0,NUMANGLES))
  endif
  ftgtau_ad=0; total_ad=0

  do i=0,curtime
     call vectdpot(i*dt,0,pots1,-1)  !! LENGTH GAUGE
     if (pulsewindowtoo == 0) then
        pulseft(i,:)=pots1(:)
     else
        pulseft(i,:)=pots1(:) * windowfunct(i,curtime)
     endif
  enddo
  
  do i=1,3
     call zfftf_wrap(2*curtime+1,pulseft(-curtime:curtime,i))
  enddo

  pulseft(:,:)=pulseft(:,:)*dt

  do i=-curtime,curtime
     pulseft(i,:)=pulseft(i,:)*exp((0.d0,1.d0)*(curtime+i)*curtime*2*pi/real(2*curtime+1))
  enddo

  pulseftsq(:) = abs(pulseft(:,1)**2) + abs(pulseft(:,2)**2) + abs(pulseft(:,3)**2)

  estep=2*pi/dt/(2*curtime+1)
  
  do imc=1,mcscfnum

     total(:)=0d0
     total_ad(:,:)=0d0

     do istate=1,nstate

        ftgtau(:)=0d0;

        do i=0,curtime
           ftgtau(i) = ALLCON(gtau(i,istate,imc))   * windowfunct(i,curtime) * &
                exp((0.d0,-1.d0)*ALLCON(ceground)*dt*i)
        enddo
        do i=1,curtime
           ftgtau(-i) = ALLCON(ftgtau(i))
        enddo

        if (angularflag.ne.0) then
           do i=0,curtime
              ftgtau_ad(i,:) = ALLCON(gtau_ad(i,istate,imc,:))   * windowfunct(i,curtime) * &
                   exp((0.d0,-1.d0)*ALLCON(ceground)*dt*i)
           enddo
           do i=1,curtime
              ftgtau_ad(-i,:) = ALLCON(ftgtau_ad(i,:))
           enddo
        endif

        write(xstate0,'(I4)') istate+1000;  write(xmc0,'(I4)') imc+1000
        xstate1=xstate0(2:4);  xmc1=xmc0(2:4)

!! no more gtaufile

        call zfftf_wrap_diff(2*curtime+1,ftgtau(-curtime:curtime),ftdiff)
        
        ftgtau(:)=ftgtau(:)*dt

        do i=-curtime,curtime
           ftgtau(i)=ftgtau(i)*exp((0.d0,1.d0)*(curtime+i)*curtime*2*pi/real(2*curtime+1))
        enddo

        total(:)=total(:)+ftgtau(:)

        if (angularflag.ne.0) then
           do il=1,NUMANGLES
              call zfftf_wrap_diff(2*curtime+1,ftgtau_ad(-curtime:curtime,il),ftdiff)
           enddo

           ftgtau_ad(:,:)=ftgtau_ad(:,:)*dt

           do i=-curtime,curtime
              ftgtau_ad(i,:)=ftgtau_ad(i,:)*exp((0.d0,1.d0)*(curtime+i)*curtime*2*pi/real(2*curtime+1))
           enddo

           total_ad(:,:)=total_ad(:,:)+ftgtau_ad(:,:)
        endif

        if(myrank.eq.1) then

           open(1004,file=projspifile(1:getlen(projspifile))//"_"//xstate1//"_"//xmc1//".dat",&
                status="replace",action="readwrite",position="rewind",iostat=myiostat)
           call checkiostat(myiostat,"opening proj spi file")
           write(1004,*,iostat=myiostat)
           call checkiostat(myiostat,"writing proj spi file")
           write(1004,*) "# Omega; pulse ft; proj flux t= ",finaltime

!! nevermind KVLsum.dat and gtau.dat
!! now integrals dt and domega with projfluxtsumfile, angprojfluxtsumfile in namelist &parinp
!! and additional column in projspifile

!! integral domega in ftgtausum

           ftgtausum=0d0
           do i=-curtime,curtime
              wfi=(i+curtime)*estep

              ftgtausum=ftgtausum + ftgtau(i) * estep / PI / 4d0   !! 4 looks correct

!! LENGTH GAUGE WAS FT'ed multiply by wfi don't divide
!! NEVERMIND FACTOR OF 1/3
!!              myfac = 5.291772108d0**2 / 3d0 * 2d0 * PI / 1.37036d2 * wfi 

!! WITH THIS FACTOR, NOW THE QUANTUM MECHANICAL PHOTOIONIZATION CROSS SECTION IN 
!! MEGABARNS (10^-18 cm^2) IS IN COLUMN 3 REAL PART
              myfac = 5.291772108d0**2 * 2d0 * PI / 1.37036d2 * wfi 

              write(1004,'(F18.12, T22, 400E20.8)',iostat=myiostat)  wfi,  &
                   pulseftsq(i), ftgtau(i)/pulseftsq(i) * myfac, ftgtau(i), ftgtausum

           enddo
           call checkiostat(myiostat,"writing proj spi file")
           close(1004)

           if (angularflag.ne.0) then
              open(1004,file=angprojspifile(1:getlen(angprojspifile))//"_"//xstate1//"_"//xmc1//".dat",&
                   status="replace",action="readwrite",position="rewind",iostat=myiostat)
              call checkiostat(myiostat,"opening angproj spi file")
              write(1004,*,iostat=myiostat)
              call checkiostat(myiostat,"writing angproj spi file")
              write(1004,*) "# Omega; pulse ft; differential proj flux t= ",finaltime
              do i=-curtime,curtime
                 wfi=(i+curtime)*estep
                 myfac = 5.291772108d0**2 * 2d0 * PI / 1.37036d2 * wfi 

                 do il=1,NUMANGLES
                    write(1004,'(F18.12, I5, 1400E20.8)',iostat=myiostat)  wfi, il, &
                         pulseftsq(i), ftgtau_ad(i,il)/pulseftsq(i) * myfac, ftgtau_ad(i,il)
                 enddo
                 write(1004,*)
              enddo
              call checkiostat(myiostat,"writing angproj spi file")
              close(1004)
           endif
        endif  !! if myrank.eq.1
     enddo  !! do istate

     if(myrank.eq.1) then
        open(1004,file=projspifile(1:getlen(projspifile))//"_all_"//xmc1//".dat",&
             status="replace",action="readwrite",position="rewind",iostat=myiostat)
        call checkiostat(myiostat,"opening total proj spi file")
        write(1004,*,iostat=myiostat)
        call checkiostat(myiostat,"writing total proj spi file")
        write(1004,*) "# Omega; pulse ft; proj flux t= ",finaltime
        do i=-curtime,curtime
           wfi=(i+curtime)*estep
           myfac = 5.291772108d0**2 * 2d0 * PI / 1.37036d2 * wfi
           write(1004,'(F18.12, T22, 400E20.8)',iostat=myiostat)  wfi,  pulseftsq(i), &
                total(i)/pulseftsq(i) * myfac, total(i)
        enddo
        call checkiostat(myiostat,"writing total proj spi file")
        close(1004)

        if (angularflag.ne.0) then
           open(1004,file=angprojspifile(1:getlen(angprojspifile))//"_all_"//xmc1//".dat",&
                status="replace",action="readwrite",position="rewind",iostat=myiostat)
           call checkiostat(myiostat,"opening total ang proj spi file")
           write(1004,*,iostat=myiostat)
           call checkiostat(myiostat,"writing total ang proj spi file")
           write(1004,*) "# Omega; pulse ft; differential proj flux t= ",finaltime
           do i=-curtime,curtime
              wfi=(i+curtime)*estep
              myfac = 5.291772108d0**2 * 2d0 * PI / 1.37036d2 * wfi

              do il=1,NUMANGLES
                 write(1004,'(F18.12, I5, 400E20.8)',iostat=myiostat)  wfi, il, pulseftsq(i), &
                      total_ad(i,il)/pulseftsq(i) * myfac, total_ad(i,il)
              enddo
              write(1004,*)
           enddo
           call checkiostat(myiostat,"writing total ang proj spi file")
           close(1004)
        endif

     endif

  enddo  !! do imc

  deallocate(ftgtau,pulseft,pulseftsq,total)
  deallocate(read_bramo,read_ketmo)
  deallocate(gtau,ketmo,bramo,ketop,gtaunow,gtausum)
  deallocate(gtau_ad,ftgtau_ad,total_ad,deweighted_bramo,ketmo_ad,ketop_ad,gtaunow_ad,gtausum_ad)

contains

  subroutine only_one_angle(howmany,iangle,inout)
    integer,intent(in) :: howmany,iangle
    DATATYPE,intent(inout) :: inout(spfdims(1),spfdims(2),spfdims(3),howmany)
    if (NUMANGLES.ne.spfdims(2)) then
       OFLWR "WHAAAAAAATTTT?? onlyoneangle",spfdims(2),NUMANGLES; CFLST
    endif
    inout(:,1:iangle-1,:,:)=0
    inout(:,iangle+1:NUMANGLES,:,:)=0
  end subroutine only_one_angle

  function myhermdots(bra,ket,totsize)
    integer,intent(in) :: totsize
    DATATYPE,intent(in) :: bra(NUMERAD,NUMANGLES,totsize/NUMERAD/(NUMANGLES)),&
         ket(NUMERAD,NUMANGLES,totsize/NUMERAD/(NUMANGLES))
    DATATYPE :: myhermdots(NUMANGLES)
    DATATYPE :: tempbra(NUMERAD,totsize/NUMERAD/(NUMANGLES)), &
         tempket(NUMERAD,totsize/NUMERAD/(NUMANGLES)),&            !! AUTOMATIC
         temparray(NUMANGLES)
    integer :: ibb,il

    ibb=(totsize/NUMERAD/(NUMANGLES)) * NUMERAD * (NUMANGLES)
    if (ibb.ne.totsize) then
       OFLWR "MYHERMDOTS ERROR ",ibb,totsize,NUMERAD,NUMANGLES; CFLST
    endif

    temparray=0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(il,tempbra,tempket)
!$OMP DO SCHEDULE(DYNAMIC)
    do il=1,NUMANGLES
       tempbra(:,:)=bra(:,il,:);    tempket(:,:)=ket(:,il,:)
       temparray(il)=hermdot(tempbra,tempket,totsize/(NUMANGLES))
    enddo
!$OMP END DO
!$OMP END PARALLEL
    myhermdots(:) = temparray(:)

  end function myhermdots


!! get the contraction of the flux operator (iH-H^\dag) with our current set of orbitals 
!! F*inspfs=outspfs

  subroutine projeflux_op_onee(inspfs,outspfs)
    use r_parameters
    use spfsize_parameters
    use opmod
    implicit none
    DATATYPE,intent(in) :: inspfs(spfsize,numr)
    DATATYPE,intent(out) :: outspfs(spfsize,numr)
    DATATYPE :: workspfs(spfsize,numr), mymat(numr,numr)    !! AUTOMATIC
    integer :: r

    workspfs=0; outspfs=0; mymat=0

    if (nucfluxopt.ne.2) then

!! terms from antihermitian part of B-O hamiltonian and hermitian part of Y operator
!! which are nonzero in the external region for electrons and so characterize the
!! electronic flux

!! terms from antihermitian of B-O hamiltonian 

       call mult_imke(numr,inspfs(:,:),workspfs(:,:))
       do r=1,numr
          outspfs(:,r) = outspfs(:,r) + workspfs(:,r) * real(1d0 / bondpoints(r)**2,8)
       enddo
    
       if(FluxOpType.eq.0.or.FluxOpType.eq.2) then
          call mult_impot(numr,inspfs(:,:),workspfs(:,:))
       else if(FluxOpType.eq.1) then
          call mult_imhalfniumpot(numr,inspfs(:,:),workspfs(:,:))
       endif
       do r=1,numr
          outspfs(:,r) = outspfs(:,r) + workspfs(:,r) * real(1d0 / bondpoints(r),8)
       enddo

       if (nonuc_checkflag.eq.0) then

!! hermitian part of antisymmetric operator is imaginary valued
          call op_imyderiv(numr,inspfs,workspfs)

!! antihermitian part of antisymmetric operator is real valued
          mymat(:,:)=real(proderivmod(:,:),8)

!! transpose verified 2nd argument 'N' not 'T' (with 'T' results are way off)
          call MYGEMM('N','N',spfsize,numr,numr,DATAONE,workspfs,spfsize,mymat,numr,&
               DATAONE,outspfs,spfsize)

       endif  !! nonuc_checkflag

    endif   !! nucfluxopt.ne.2

    if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0) then

!! terms from imaginary part of inverse bond length, hermitian part of Y operator, and
!! antihermitian part of nuclear ke second derivative operator, all of which are 
!! nonzero only in external region for nuclei and so characterize the nuclear flux

       call mult_reke(numr,inspfs(:,:),workspfs(:,:))
       do r=1,numr
          outspfs(:,r) = outspfs(:,r) + workspfs(:,r) * imag((0d0,0d0)+ 1d0 / bondpoints(r)**2)
       enddo
    
       if(FluxOpType.eq.0.or.FluxOpType.eq.2) then
          call mult_repot(numr,inspfs(:,:),workspfs(:,:))
       else if(FluxOpType.eq.1) then
          call mult_rehalfniumpot(numr,inspfs(:,:),workspfs(:,:))
       endif
       do r=1,numr
          outspfs(:,r) = outspfs(:,r) + workspfs(:,r) * imag((0d0,0d0)+ 1d0 / bondpoints(r))
       enddo

!! real part of antisymmetric operator is antihermitian
       call op_reyderiv(numr,inspfs,workspfs)

!! imag part is hermitian
       mymat(:,:)=imag(proderivmod(:,:)+(0d0,0d0))

!! transpose verified 2nd argument 'N' not 'T' (with 'T' results are way off)
       call MYGEMM('N','N',spfsize,numr,numr,DATAONE,workspfs,spfsize,mymat,numr,&
            DATAONE,outspfs,spfsize)

       mymat(:,:)=imag(rkemod(:,:)+(0d0,0d0))
       call MYGEMM('N','N',spfsize,numr,numr,DATAONE,inspfs,spfsize,mymat,numr,&
            DATAONE,outspfs,spfsize)

    endif

!! scale correctly
    outspfs(:,:)=outspfs(:,:)*(-2d0)

  end subroutine projeflux_op_onee

    
end subroutine projeflux_double_time_int




module projbiomod
  use biorthotypemod
  implicit none
  type(biorthotype),target :: projbiovar
end module projbiomod

subroutine projeflux_single0(ifile,nt,alreadystate,outnumstate)
  use parameters    !! catavectorfiles and others
  use projbiomod
  use biorthomod
  use configmod
  use projefluxmod
  use mpimod
  use mpisubmod
  implicit none
!! necessary working variables
  integer,intent(in) :: nt,alreadystate,ifile
  integer,intent(out) :: outnumstate
  integer :: tau, i,ir,tnum2part,tnspf,tnumr,istate,ierr, cgflag, &
       spfcomplex, acomplex, tdims(3),imc,ioffset,myiostat, &
       targetms, targetrestrictflag, targetspinproject, targetspinval
  real*8 :: cgfac,doubleclebschsq,aa,bb,cc
  DATATYPE, allocatable :: &
       tmo(:,:),tavec(:,:,:),tmotemp(:,:),readta(:,:,:),&
       mobio(:,:,:),abio(:,:,:),mymo(:,:),myavec(:,:,:), &
       readmo(:,:),readavec(:,:,:)
  DATATYPE :: nullvector(numr)
!! mcscf specific read variables
  DATATYPE,target :: smo(nspf,nspf)

  nullvector(:)=0

!! read in the data from mcscf for our target cation state

  if (myrank.eq.1) then
     open(909,file=catspffiles(ifile),status="unknown",form="unformatted",iostat=myiostat)
     call checkiostat(myiostat,"opening "//catspffiles(ifile))
     open(910,file=catavectorfiles(ifile),status="unknown",form="unformatted",iostat=myiostat)
     call checkiostat(myiostat,"opening "//catavectorfiles(ifile))
     call avector_header_read(910,outnumstate,tnum2part,tnumr,tnumconfig,targetrestrictflag,targetms,&
          targetspinproject,targetspinval,acomplex,ierr)
  endif
  call mympiibcastone(outnumstate,1); call mympiibcastone(tnum2part,1); 
  call mympiibcastone(tnumr,1); call mympiibcastone(tnumconfig,1); 
  call mympiibcastone(targetrestrictflag,1);  call mympiibcastone(targetms,1);  
  call mympiibcastone(targetspinproject,1);  call mympiibcastone(targetspinval,1); 
  call mympiibcastone(acomplex,1);  call mympiibcastone(ierr,1)

  if (ierr.ne.0) then
     OFLWR "Avector header read error; redone 09/29/2014; recompute vector on disk."; CFLST
  endif

  if (myrank.eq.1) then
     call spf_header_read(909,tdims,tnspf,spfcomplex)
  endif
  call mympiibcast(tdims,1,3);
  call mympiibcastone(tnspf,1);
  call mympiibcastone(spfcomplex,1)

!! have to project on BO wfns.  Otherwise doesn't make sense in prolate.  
!! Not supported anymore, will support
!!   again, mcscf mode w/many r's on file

  if (tnspf.gt.nspf+numfrozen) then
     OFLWR "ERROR, for now can't do more orbs in projection than in calculation",tnspf,nspf+numfrozen; CFLST
  endif
  if ( (holeflag.eq.0.and.tnum2part.ne.num2part-2) .or. (holeflag.ne.0.and.tnum2part.ne.num2part+2) ) then
        OFLWR "Vectors read are not n-1 electron functions"; CFLST
  endif

  cgflag=0
  if (restrictflag.eq.0) then
     OFLWR " WARNING: propagated state is not spin projection restricted.  No CG algebra."
     CFL
  else
     if (targetrestrictflag.eq.0) then
        OFLWR "WARNING: N-1 electron state is not spin projection restricted (restrictflag=0).  No CG algebra."
        CFL
     else
        if (abs(targetms-restrict_ms).ne.1) then
           OFLWR "Targetms should differ from restrictms by 1.  Cation and neutral have no dyson orbital"
           WRFL  "    Targetms=", targetms, " restrictms=",restrict_ms; CFLST
        endif
        if (all_spinproject.eq.0) then
           OFLWR "WARNING: Wave function is not spin restricted (allspinproject=0).  No CG algebra.";CFL
        else
           if (targetspinproject.eq.0) then
              OFLWR "WARNING: N-1 electron state is not spin restricted (allspinproject=0).  No CG algebra."
              CFL
           else
              if (abs(targetspinval-spin_restrictval).gt.1) then
                 OFLWR "ERROR: Spin value of wave function and N-1 electron state differ by more than 1/2.",&
                      targetspinval,spin_restrictval; CFLST
              else
                 OFLWR "Spin adapted wave functions for propagation and projection. doing CG algebra.";CFL
                 cgflag=1
              endif
           endif
        endif
     endif
  endif


!!!!  100912 NEED CG COEFFICIENT RATIO 
              
!! programmed this with restrict_ms, targetms.  
!! Do have option to change spinrestrictval; should use this.
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
!!    if restrict_ms > targetms  then we couple high spin to high spin to get high spin, one CG coef = 1.
!!    
!!        otherwise we couple high spin target down to high spin neutral:
!!
!!         projection is spin down  x high spin target   e.g.   < s' s' |  s s  1/2 -1/2 > = A
!!         also have < s' s' | s s-1 1/2 1/2 > = B
!!  
!!         so mult cross section by ( A^2 + B^2 / A^2 )
!!    
!!   BUT if restrictflag is off, can't do this.  Must realize that you are only projecting one 
!!     component of C.G. sum.  that would be a problem.  But if you were to calculate say the 
!!     3 lowest states for a triplet, you'd get all components.
!!     otherwise cross section is meaningless.
!! all arguments are x 2, integers for half spin.

!!$                 if (restrict_ms.lt.targetms.and.restrictflag.eq.1) then
!!$  !! what we're calculating  targetms-restrict_ms=+/-1
!!$                    aa= doubleclebschsq(targetms,1,targetms,restrict_ms-targetms,restrict_ms) 
!!$  !! only one could be nonzero obviously
!!$                   bb= doubleclebschsq(targetms,1,targetms-2,restrict_ms-targetms+2,restrict_ms)  
!!$                    cc= doubleclebschsq(targetms,1,targetms+2,restrict_ms-targetms-2,restrict_ms)      
!!$                    cgfac = (aa+bb+cc)/aa
!!$                 else
!!$                    cgfac=1
!!$                 endif

  if (cgflag.eq.0) then
     cgfac=1
  else
!! what we're calculating  targetms-restrict_ms=+/-1
     aa= doubleclebschsq(targetspinval,1,targetms,restrict_ms-targetms,spin_restrictval)
!! only one could be nonzero obviously
     bb= doubleclebschsq(targetspinval,1,targetms-2,restrict_ms-targetms+2,spin_restrictval)
     cc= doubleclebschsq(targetspinval,1,targetms+2,restrict_ms-targetms-2,spin_restrictval)      
     cgfac = (aa+bb+cc)/aa
  endif

  allocate(tmo(spfsize,nspf),tavec(tnumconfig,outnumstate,numr),&
       tmotemp(spfsize,nspf+numfrozen),readta(tnumr,tnumconfig,outnumstate))
  tmo=0d0;  tavec=0d0; tmotemp=0d0; readta=0d0

  OFLWR "Reading", outnumstate," Born-Oppenheimer states."; CFL

  if (myrank.eq.1) then
     call simple_load_avectors(910,acomplex,readta(:,:,:),tnum2part,tnumr,tnumconfig,outnumstate)
     do ir=1,min(tnumr,numr)
        tavec(:,:,ir)=readta(ir,:,:)
     enddo
     do ir=min(tnumr,numr)+1,numr
        call staticvector(tavec(:,:,ir),tnumconfig*outnumstate)
     enddo
     do ir=1,numr

!! projecting on normalized electronic wfn at each R
!! we go to war with the army we've got

        do istate=1,outnumstate
           tavec(:,istate,ir)=tavec(:,istate,ir)/&
                sqrt(dot(tavec(:,istate,ir),tavec(:,istate,ir),tnumconfig)) !! no * bondweights(ir)
        enddo
     enddo
     call spf_read0(909,nspf+numfrozen,spfdims,tnspf,tdims,spfcomplex,spfdimtype,&
          tmotemp(:,:),(/0,0,0/))
     close(909);  close(910)
  else
     call spf_read0(-42,nspf+numfrozen,spfdims,tnspf,tdims,spfcomplex,spfdimtype,&
          tmotemp(:,:),(/0,0,0/))
  endif
  call mympibcast(tavec(:,:,:),1,tnumconfig*outnumstate*numr)

  tmo(:,1:tnspf-numfrozen) = tmotemp(:,numfrozen+1:tnspf)
  tnspf=tnspf-numfrozen

  do i=tnspf+1,nspf
     tmo(:,i)=0d0
     call staticvector(tmo(:,i),spfsize)
     if (parorbsplit.eq.3) then
        call gramschmidt(spfsize,i-1,spfsize,tmo(:,:),tmo(:,i),.true.)
     else
        call gramschmidt(spfsize,i-1,spfsize,tmo(:,:),tmo(:,i),.false.)
     endif
  enddo

  allocate(tconfiglist(tnum2part,tnumconfig));    tconfiglist=0

  if (myrank.eq.1) then
     open(910,file=catavectorfiles(ifile),status="unknown",form="unformatted",iostat=myiostat)
     call checkiostat(myiostat,"opening "//catavectorfiles(ifile))
     call avector_header_read_simple(910,outnumstate,tnum2part,tnumr,tnumconfig,acomplex)
     call get_avectorfile_configlist(910,acomplex,tconfiglist,tnum2part,tnumr,tnumconfig)
     close(910)
  endif

  call mympiibcast(tconfiglist,1,tnum2part*tnumconfig)

!! do the walks from the target state into our final state

  call projeflux_singlewalks()

!! allocate all necessary extra memory and io params to do this looping business

  OFLWR "Allocating arrays for projected flux"; CFL

  allocate(mobio(spfsize,nspf,numr),abio(first_config:last_config,mcscfnum,numr),&
       mymo(spfsize,nspf),  myavec(numr,first_config:last_config,mcscfnum))
  mobio=0; mymo=0
  if (last_config.ge.first_config) then
     abio=0; myavec=0
  endif
  if (myrank.eq.1) then
     if (parorbsplit.eq.3) then
        allocate(readmo(spfsize*nprocs,nspf))
     else
        allocate(readmo(spfsize,nspf))
     endif
     allocate(readavec(numr,num_config,mcscfnum))
  else
     allocate(readmo(1,nspf),readavec(1,numr,mcscfnum))
  endif
  readmo=0; readavec=0

  call mpibarrier()
  OFLWR "    ...allocated arrays for projected flux.  go time loop."; CFL

  if (myrank.eq.1) then
     inquire (iolength=i) readmo
     open(1001,file=fluxmofile,status="unknown",form="unformatted",&
          access="direct",recl=i,iostat=myiostat)
     call checkiostat(myiostat,"opening "//fluxmofile)
     inquire (iolength=i) readavec
     open(1002,file=fluxafile,status="unknown",form="unformatted",&
          access="direct",recl=i,iostat=myiostat)
     call checkiostat(myiostat,"opening "//fluxafile)
  endif

!! do the loop over all time, ALL TIME and now in parallel-o-vision

  do tau=0,nt    
!! read in this time's wavefucntion

     if (myrank.eq.1) then
        read(1001,rec=FluxSkipMult*tau+1,iostat=myiostat) readmo(:,:)
        call checkiostat(myiostat,"reading fluxmofile")
        read(1002,rec=FluxSkipMult*tau+1,iostat=myiostat) readavec(:,:,:)
        call checkiostat(myiostat,"reading fluxafile")
     endif
     if (parorbsplit.ne.3) then
        if (myrank.eq.1) then
           mymo(:,:)=readmo(:,:)
        endif
        call mympibcast(mymo(:,:),1,totspfdim)
     else
        do i=1,nspf
           call splitscatterv(readmo(:,i),mymo(:,i))
        enddo
     endif
     if (par_consplit.eq.0) then
        if (myrank.eq.1) then
           myavec(:,:,:)=readavec(:,:,:)
        endif
        call mympibcast(myavec(:,:,:),1,numr*num_config*mcscfnum)
     else
        do i=1,mcscfnum
           if (tot_adim.gt.0) then
              call myscatterv(readavec(:,:,i),myavec(:,:,i),configs_perproc(:)*numr)
           else
              call myscatterv(readavec(:,:,i),nullvector(:),configs_perproc(:)*numr)
           endif
        enddo
     endif

!! do biortho and construct the single particle function


!!$       if (1==0) then
!!$          do ir=1,numr
!!$             abio(:,:,ir)=myavec(ir,:,:)
!!$          enddo
!!$          do ir=1,numr
!!$             call bioset(projbiovar,smo,1,bioww); 
!!$             call biortho(mymo,tmo(:,:),mobio(:,:,ir),abio(:,1,ir),projbiovar)
!!$             do imc=2,mcscfnum
!!$                call biotransform(mymo,tmo(:,:),abio(:,imc,ir),projbiovar)
!!$             enddo
!!$          enddo
!!$       else

!! this should do the same and it looks like it does
     
     call bioset(projbiovar,smo,numr,bioww)

     if (tot_adim.gt.0) then
        call biortho(mymo,tmo(:,:),mobio(:,:,1),myavec(:,:,1),projbiovar)
        do imc=2,mcscfnum
           call biotransform(mymo,tmo(:,:),myavec(:,:,imc),projbiovar)
        enddo
        do ir=1,numr
           mobio(:,:,ir)=mobio(:,:,1)
           abio(:,:,ir)=myavec(ir,:,:)
        enddo
     else
        call biortho(mymo,tmo(:,:),mobio(:,:,1),nullvector(:),projbiovar)
        do imc=2,mcscfnum
           call biotransform(mymo,tmo(:,:),nullvector(:),projbiovar)
        enddo
     endif

     do ir=1,numr
        do imc=1,mcscfnum
           do istate=1,outnumstate

   ioffset=(istate-1+alreadystate)*mcscfnum*(nt+1)*numr + (imc-1)*(nt+1)*numr + tau*numr + ir

              if (tot_adim.gt.0) then
                 call projeflux_doproj(tavec(:,istate,ir),abio(:,imc,ir),mobio(:,:,ir),ioffset,cgfac*catfacs(ifile))
              else
                 call projeflux_doproj(tavec(:,istate,ir),nullvector(:),mobio(:,:,ir),ioffset,cgfac*catfacs(ifile))
              endif
           enddo
        enddo
     enddo

!! write out times cause we're bored
!    if(mod(tau,100).eq.0.or.tau.eq.nt) then
!      call openfile
!      write(mpifileptr,'(A52,F10.4)') " Timing Statistics: producing 1e- proj wfn as of T= ",tau*dt
!      write(mpifileptr,'(100A10)') "Times: ", "All", "IO", "Walks", "Biorth", "1e Bld", "1e Op", &
!          "Fluxeval", "FT g!tau"
!      write(mpifileptr,'(A10,100I10)') " ", times(1:8)/100
!    endif
  enddo

!! clean up
  close(1001);  close(1002)
  deallocate(readmo,readavec)
  deallocate(mobio,abio,mymo,myavec)
  deallocate(tavec,tmotemp,tmo,readta)
  deallocate(tconfiglist,numpwalk1,pwalk1,pspf1,pphase1)


end subroutine projeflux_single0



subroutine projeflux_single(mem)
  use parameters    !! numcatfiles, fluxinterval, par_timestep, etc.
  use projefluxmod  !! nstate, eachstate
  use mpimod   !! myrank
  implicit none
  integer,intent(in) :: mem
  integer :: nt,ifile,myiostat

  OFLWR ;  WRFL   "   *** DOING PROJECTED FLUX. ***    ";  WRFL; CFL

  nt=floor(real(numpropsteps,8)/fluxinterval/fluxskipmult)

  nstate=0
  allocate(eachstate(numcatfiles))
  eachstate=0

  do ifile=1,numcatfiles
     call projeflux_single0(ifile,nt,nstate,eachstate(ifile))
     nstate=nstate+eachstate(ifile)
  enddo

  call mpibarrier()
  OFLWR "Go double time integral";CFL

!! do the double time integral
  call projeflux_double_time_int(mem,nstate,nt)

  call mpibarrier()
  OFLWR "Cross Section acquired, cleaning and closing";CFL

  if (myrank.eq.1) then
     open(1003,file=projfluxfile,status="unknown",iostat=myiostat)
     write(1003,*) "Projfluxfile deleted"
     close(1003)
  endif
  call mpibarrier()
  call mpistop()

end subroutine projeflux_single

