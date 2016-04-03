

!! ELECTRON FLUX ROUTINES - ACTION 16 analysis total photoionization, 
!! ACTION 15 save file, saved file may be used for actions 16 17 23)


#include "Definitions.INC"

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


module fluxgtaubiomod
  use biorthotypemod
  implicit none
  type(biorthotype),target :: fluxgtaubiovar
end module fluxgtaubiomod


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


module fluxgtau0mod
contains

  subroutine fluxgtau0(alg,www,bioww)
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
    implicit none
!! FluxOpType:
!! 0       = use one-e potential and two-e contribution routines  (exact treatment)
!! 1       = replace one-e potential + two-e potential with halfnium one-e potential  (recommended)
!! 2       = use full one-e potential, no two-e 
!! other   = only KE

    type(walktype),target :: www,bioww
    integer,intent(in) :: alg
    integer :: curtime,oldtime,k,nt,i,molength,alength,  BatchSize,NBat,brabat,brareadsize, &
         bratime,ketbat,ketreadsize,kettime,bratop, atime,btime,itime,jtime,times(1:7)=0, &
         imc, tau, ispf, myiostat
    real*8 :: MemTot,MemVal,dt, myfac,wfi,estep,windowfunct
    complex*16, allocatable :: FTgtau(:,:), pulseft(:,:)
    real*8, allocatable :: pulseftsq(:)
    DATATYPE :: fluxevalval(mcscfnum), pots1(3)=0d0  !!$, pots2(3)=0d0
    DATATYPE, allocatable :: gtau(:,:), mobio(:,:),abio(:,:,:)
    DATATYPE, allocatable, target :: bramo(:,:,:),braavec(:,:,:,:), ketmo(:,:,:),ketavec(:,:,:,:)
    DATATYPE,allocatable :: read_bramo(:,:,:),read_braavec(:,:,:,:),read_ketmo(:,:,:),read_ketavec(:,:,:,:)
    DATATYPE, pointer :: moket(:,:),mobra(:,:),aket(:,:,:)
    DATATYPE, allocatable :: imke(:,:),impe(:,:),imV2(:,:,:,:), imyderiv(:,:)
    DATATYPE, allocatable :: imkeop(:,:),impeop(:,:),  imyop(:,:)
    DATATYPE, allocatable :: reke(:,:),repe(:,:),reV2(:,:,:,:), reyderiv(:,:)
    DATATYPE, allocatable :: rekeop(:,:),repeop(:,:),  reyop(:,:)
    DATATYPE,target :: smo(nspf,nspf)
    DATATYPE :: nullvector(numr)

    if (ceground.eq.(0d0,0d0)) then
       OFLWR "Eground is ZERO.  Are you sure?  If want zero just make it small. \n     Otherwise need eground: initial state energy."; CFLST
    endif

    nullvector(:)=0

!! initial setup

!!  dt=real(FluxInterval*FluxSkipMult,8)*par_timestep;  nt=floor(final time/dt)

    dt=real(FluxInterval*FluxSkipMult,8)*par_timestep;  nt=floor(real(numpropsteps,8)/fluxinterval/fluxskipmult)

    allocate(gtau(0:nt,mcscfnum))
    gtau(:,:)=0d0

    allocate(imke(nspf,nspf),impe(nspf,nspf),imV2(nspf,nspf,nspf,nspf),imyderiv(nspf,nspf),&
         reke(nspf,nspf),repe(nspf,nspf),reV2(nspf,nspf,nspf,nspf),reyderiv(nspf,nspf),&
         rekeop(spfsize,nspf),repeop(spfsize,nspf),reyop(spfsize,nspf),&
       mobio(spfsize,nspf),abio(numr,www%firstconfig:www%lastconfig,mcscfnum), &
       imkeop(spfsize,nspf),impeop(spfsize,nspf),imyop(spfsize,nspf))

    imke=0; impe=0; imV2=0; imyderiv=0; imkeop=0; impeop=0; imyop=0;
    reke=0; repe=0; reV2=0; reyderiv=0; rekeop=0; repeop=0; reyop=0
    mobio=0;
    if (www%lastconfig.ge.www%firstconfig) then
       abio=0
    endif

!! determine if we should do batching or not
!! 250,000 words/MB, real*8 2words/#, complex*16 4words/#

#ifdef REALGO
    MemVal = 1.25d5
#else
    MemVal = 6.25d4
#endif
    call openfile()
    write(mpifileptr,'(A30,F9.3,A3)') " Guess at necessary memory is ",&
         2d0*real((nt+1)*(www%localnconfig*numr*mcscfnum+spfsize*nspf),8)/MemVal," MB"
    if(alg.eq.0) then
       write(mpifileptr,*) "g(tau) will be computed with all of psi in core"
       BatchSize=nt+1
    else
       MemTot=real(alg,8)    
       write(mpifileptr,*) "g(tau) will be computed with all psi being read in batches"
       write(mpifileptr,'(A33,F9.3,A3)') " Desired amount of memory to use ",MemTot," MB"
       BatchSize=floor(MemTot * MemVal / (2d0*real(www%localnconfig*numr*mcscfnum+spfsize*nspf,8)))
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

    allocate(ketmo(spfsize,nspf,BatchSize),ketavec(numr,www%firstconfig:www%lastconfig,mcscfnum,BatchSize),&
         bramo(spfsize,nspf,BatchSize),braavec(numr,www%firstconfig:www%lastconfig,mcscfnum,BatchSize))
    ketmo=0; bramo=0
    if (www%lastconfig.ge.www%firstconfig) then
       ketavec=0; braavec=0
    endif
    if (myrank.eq.1) then
       if (parorbsplit.eq.3) then
          allocate(read_bramo(spfsize*nprocs,nspf,BatchSize),read_ketmo(spfsize*nprocs,nspf,BatchSize))
       else
          allocate(read_bramo(spfsize,nspf,BatchSize),read_ketmo(spfsize,nspf,BatchSize))
       endif
    else
       allocate(read_bramo(1,nspf,BatchSize),read_ketmo(1,nspf,BatchSize))
    endif
    read_bramo=0; read_ketmo=0

    if (myrank.eq.1) then
       allocate(read_braavec(numr,www%numconfig,mcscfnum,BatchSize),&
            read_ketavec(numr,www%numconfig,mcscfnum,BatchSize))
    else
       allocate(read_braavec(1,1,mcscfnum,BatchSize),read_ketavec(1,1,mcscfnum,BatchSize))
    endif
    read_braavec=0; read_ketavec=0

    NBat=ceiling(real(nt+1)/real(BatchSize))
    ketreadsize=0;  brareadsize=0

    if (myrank.eq.1) then
       inquire (iolength=molength) read_ketmo(:,:,1);  inquire (iolength=alength) read_ketavec(:,:,:,1)
    endif
    call mympiibcastone(molength,1); call mympiibcastone(alength,1)

    call openfile()
    write(mpifileptr,*) "MO record length is ",molength
    write(mpifileptr,*) "AVEC record length is ",alength
    call closefile()

    if (myrank.eq.1) then
       open(454, file="Dat/KVLsum.dat", status="unknown",iostat=myiostat)
       call checkiostat(myiostat,"opening kvlsum file")
       write(454,*,iostat=myiostat) "#KVL flux sum: itime, time, flux sum"
       call checkiostat(myiostat,"writing kvl sum file")
       write(454,*);  close(454)
    endif

!! begin the ket batch read loop
    do ketbat=1,NBat
       call system_clock(atime)
       OFLWR "Reading ket batch ", ketbat, " of ", NBat; CFL
       ketreadsize=min(BatchSize,nt+1-(ketbat-1)*BatchSize)
       if(myrank.eq.1) then
          open(1001,file=fluxmofile,status="old",form="unformatted",&
               access="direct",recl=molength,iostat=myiostat)
          call checkiostat(myiostat,"opening "//fluxmofile)
          open(1002,file=fluxafile,status="old",form="unformatted",&
               access="direct",recl=alength,iostat=myiostat)
          call checkiostat(myiostat,"opening "//fluxafile)
          do i=1,ketreadsize
             k=FluxSkipMult*((ketbat-1)*BatchSize+i-1)+1
             read(1001,rec=k,iostat=myiostat) read_ketmo(:,:,i) 
          enddo
          call checkiostat(myiostat,"reading "//fluxmofile)
          do i=1,ketreadsize
             k=FluxSkipMult*((ketbat-1)*BatchSize+i-1)+1
             read(1002,rec=k,iostat=myiostat) read_ketavec(:,:,:,i) 
          enddo
          call checkiostat(myiostat,"reading "//fluxafile)
          close(1001);      close(1002)
       endif

       if (parorbsplit.ne.3) then
          if (myrank.eq.1) then
             ketmo(:,:,1:ketreadsize)=read_ketmo(:,:,1:ketreadsize)
          endif
          call mympibcast(ketmo(:,:,:),1,totspfdim*ketreadsize)
       else
          do i=1,ketreadsize
             do ispf=1,nspf
                call splitscatterv(read_ketmo(:,ispf,i),ketmo(:,ispf,i))
             enddo
          enddo
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

       call system_clock(btime)
       times(1)=times(1)+btime-atime;    times(2)=times(2)+btime-atime
     
!! begin the bra batch read loop
       do brabat=1,ketbat
          call system_clock(atime)
          OFLWR "Reading bra batch ", brabat, " of ", ketbat; CFL
          brareadsize=min(BatchSize,nt+1-(brabat-1)*BatchSize)
          if(brabat.eq.ketbat) then
             bramo=ketmo;        
             if (tot_adim.gt.0) then
                braavec=ketavec
             endif
          else 
             if(myrank.eq.1) then
                open(1001,file=fluxmofile,status="old",form="unformatted",&
                     access="direct",recl=molength,iostat=myiostat)
                call checkiostat(myiostat,"opening "//fluxmofile)
                open(1002,file=fluxafile,status="old",form="unformatted",&
                     access="direct",recl=alength,iostat=myiostat)
                call checkiostat(myiostat,"opening "//fluxafile)
                do i=1,brareadsize
                   k=FluxSkipMult*((brabat-1)*BatchSize+i-1)+1
                   read(1001,rec=k,iostat=myiostat) read_bramo(:,:,i)
                enddo
                call checkiostat(myiostat,"reading "//fluxmofile)
                do i=1,brareadsize
                   k=FluxSkipMult*((brabat-1)*BatchSize+i-1)+1
                   read(1002,rec=k,iostat=myiostat) read_braavec(:,:,:,i) 
                enddo
                call checkiostat(myiostat,"reading "//fluxafile)
                close(1001);          close(1002)
             endif

             if (parorbsplit.ne.3) then
                if (myrank.eq.1) then
                   bramo(:,:,1:brareadsize)=read_bramo(:,:,1:brareadsize)
                endif
                call mympibcast(bramo(:,:,:),1,totspfdim*brareadsize)
             else
                do i=1,brareadsize
                   do ispf=1,nspf
                      call splitscatterv(read_bramo(:,ispf,i),bramo(:,ispf,i))
                   enddo
                enddo
             endif
             if (par_consplit.eq.0) then
                if (myrank.eq.1) then
                   braavec(:,:,:,1:brareadsize)=read_braavec(:,:,:,1:brareadsize)
                endif
                call mympibcast(braavec(:,:,:,1:brareadsize),1,&
                     numr*www%numconfig*mcscfnum*brareadsize)
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
          endif

          call system_clock(btime)
          times(1)=times(1)+btime-atime;      times(2)=times(2)+btime-atime
        
!! loop over all time for the ket of the flux integral
          do kettime=1,ketreadsize

!! get the one-e half transformed matrix elements for this ket time
             call system_clock(atime)
             curtime=(ketbat-1)*BatchSize+kettime-1 
             moket=>ketmo(:,:,kettime);        aket=>ketavec(:,:,:,kettime)

             call flux_op_onee(moket,imkeop,impeop,1)
             if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0)then
                call flux_op_onee(moket,rekeop,repeop,2)
             endif
             if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0) then
                call flux_op_nuc(moket,imyop,1)
                call flux_op_nuc(moket,reyop,2)
             endif

!! determine bounds of loop over bras and setup doing in parallel with mpi!
             if(brabat.lt.ketbat) then
                bratop=brareadsize
             else
                bratop=kettime
             endif
             call system_clock(btime);        times(4)=times(4)+btime-atime
           
!! loop over all previous time for the bra of the flux integral
             do bratime=1,bratop

                oldtime=(brabat-1)*BatchSize+bratime-1
              
!! biortho this pair of times!        
                call system_clock(itime)
                mobra=>bramo(:,:,bratime)
                if (tot_adim.gt.0) then
                   abio(:,:,:)=braavec(:,:,:,bratime)
                endif
                call bioset(fluxgtaubiovar,smo,numr,bioww)
#ifdef CNORMFLAG
                fluxgtaubiovar%hermonly=.true.
#endif
                if (tot_adim.gt.0) then
                   call biortho(mobra,moket,mobio,abio(:,:,1),fluxgtaubiovar)
                   do imc=2,mcscfnum
                      call biotransform(mobra,moket,abio(:,:,imc),fluxgtaubiovar)
                   enddo
#ifdef CNORMFLAG
                   abio(:,:,:)=conjg(abio(:,:,:))
#endif
                else
                   call biortho(mobra,moket,mobio,nullvector(:),fluxgtaubiovar)
                   do imc=2,mcscfnum
                      call biotransform(mobra,moket,nullvector(:),fluxgtaubiovar)
                   enddo
                endif
                fluxgtaubiovar%hermonly=.false.   !! in case fluxgtaubiovar is reused

                call system_clock(jtime);          times(3)=times(3)+jtime-itime

!! complete the one-e potential and kinetic energy matrix elements           
!! DO BLAS !!

                do i=1,nspf
                   do k=1,nspf
                      imke(k,i) = hermdot(mobio(:,k),imkeop(:,i),spfsize)
                      impe(k,i) = hermdot(mobio(:,k),impeop(:,i),spfsize)
                      if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0)then
                         reke(k,i) = hermdot(mobio(:,k),rekeop(:,i),spfsize)
                         repe(k,i) = hermdot(mobio(:,k),repeop(:,i),spfsize)
                         imyderiv(k,i) = hermdot(mobio(:,k),imyop(:,i),spfsize)
                         reyderiv(k,i) = hermdot(mobio(:,k),reyop(:,i),spfsize)
                      endif
                   enddo
                enddo

                if (parorbsplit.eq.3) then
                   call mympireduce(imke,nspf**2)
                   call mympireduce(impe,nspf**2)
                   if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0)then
                      call mympireduce(reke,nspf**2)
                      call mympireduce(repe,nspf**2)
                      call mympireduce(imyderiv,nspf**2)
                      call mympireduce(reyderiv,nspf**2)
                   endif
                endif

                call system_clock(itime);          times(4)=times(4)+itime-jtime

!! get the two-e contribution for exact formula (fluxoptype=0)

                imV2=0d0; reV2=0
                if(FluxOpType.eq.0) then
                   call flux_op_twoe(mobio,moket,imV2,1)
                   if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0)then
                      call flux_op_twoe(mobio,moket,reV2,2)
                   endif
                endif
                call system_clock(jtime)
                times(5)=times(5)+jtime-itime

!! evaluate the actual g(tau) expression

                do imc=1,mcscfnum
                   if (tot_adim.gt.0) then
                      fluxevalval(imc) = fluxeval(abio(:,:,imc),aket(:,:,imc),&
                           imke,impe,imV2,imyderiv,reke,repe,reV2,reyderiv)
                   else
                      fluxevalval(imc) = fluxeval(nullvector(:),nullvector(:),&
                           imke,impe,imV2,imyderiv,reke,repe,reV2,reyderiv)
                   endif
                enddo

                oldtime=(brabat-1)*BatchSize+bratime-1
                tau=curtime-oldtime; 
                gtau(tau,:) = gtau(tau,:) + fluxevalval(:)

                nullify(mobra)
                call system_clock(itime);          times(6)=times(6)+itime-jtime

             enddo !! end loop over specific bras

             nullify(moket,aket)
             call system_clock(jtime);        times(6)=times(6)+jtime-itime

!! only do this after we are sure we've gone through every bra
             if (brabat.eq.ketbat) then
                if (myrank.eq.1) then
                   open(454, file="Dat/KVLsum.dat", status="old", &
                        position="append",iostat=myiostat)
                   call checkiostat(myiostat,"opening kvlsum file")
                   write(454,'(I5,100F18.12)',iostat=myiostat) curtime, curtime*dt, gtau(0,:);    
                   call checkiostat(myiostat,"writing kvlsum file")
                   close(454)
                endif
             endif

          enddo !! do kettime

          if (brabat.eq.ketbat) then
             if (notiming.ne.2) then
                OFL
                write(mpifileptr,'(A28,F10.4)') " Timing statistics as of T= ",real(curtime,8)*dt
                write(mpifileptr,'(100A10)') "Times: ", "All", "Read",&
                     "Biorth", "One-e", "Two-e", "Fluxeval", "FT gtau"
                write(mpifileptr,'(A10,100I10)') " ", times(1:7)/100; CFL
             endif
          endif
       enddo !! do brabat
    enddo !! do ketbat

    if (curtime.ne.nt) then
       OFLWR "DOOG CURTIME NT ERR",curtime,nt; CFLST
    endif

!! create the xsec as it should be now

    OFLWR " Taking the FT of g(tau) to get xsection at T= ",curtime*dt; CFL

    allocate(ftgtau(-curtime:curtime,mcscfnum), pulseft(-curtime:curtime,3), &
         pulseftsq(-curtime:curtime))
    ftgtau(:,:)=0d0; pulseft(:,:)=0d0; pulseftsq(:)=0d0

    do i=0,curtime

       ftgtau(i,:) = ALLCON(gtau(i,:))   * windowfunct(i,curtime) * &
            exp((0.d0,-1.d0)*ALLCON(ceground)*par_timestep*FluxInterval*FluxSkipMult*i)

       call vectdpot(i*par_timestep*fluxinterval*fluxskipmult,0,pots1,-1)   !! LENGTH GAUGE.
       if (pulsewindowtoo == 0) then
          pulseft(i,:)=pots1(:)
       else
          pulseft(i,:)=pots1(:) * windowfunct(i,curtime)
       endif
    enddo

    do i=1,curtime
       ftgtau(-i,:) = ALLCON(ftgtau(i,:))
    enddo
  
    if (myrank.eq.1) then
       open(171,file=gtaufile,status="unknown",iostat=myiostat)
       call checkiostat(myiostat,"opening gtaufile")
       write(171,*,iostat=myiostat) "#   ", curtime
       call checkiostat(myiostat,"writing gtaufile")
       do i=0,curtime
          write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat)  &
               i*par_timestep*FluxInterval*FluxSkipMult, ftgtau(i,:)
       enddo
       call checkiostat(myiostat,"writing gtaufile")
       close(171)
    endif

    OFLWR "   ....Go ft...."; CFL
    do imc=1,mcscfnum
       call zfftf_wrap_diff(2*curtime+1,ftgtau(-curtime:curtime,imc),ftdiff)
    enddo
    OFLWR "   ....Go ft pulse...."; CFL
    do i=1,3
       call zfftf_wrap(2*curtime+1,pulseft(-curtime:curtime,i))
    enddo
    OFLWR "   ....Done with ft...."; CFL

    ftgtau(:,:)=ftgtau(:,:)*par_timestep*FluxInterval*FluxSkipMult
    pulseft(:,:)=pulseft(:,:)*par_timestep*FluxInterval*FluxSkipMult

    do i=-curtime,curtime
       ftgtau(i,:)=ftgtau(i,:)*exp((0.d0,1.d0)*(curtime+i)*curtime*2*pi/real(2*curtime+1))
       pulseft(i,:)=pulseft(i,:)*exp((0.d0,1.d0)*(curtime+i)*curtime*2*pi/real(2*curtime+1))
    enddo

    pulseftsq(:) = abs(pulseft(:,1)**2) + abs(pulseft(:,2)**2) + abs(pulseft(:,3)**2)

    estep=2*pi/par_timestep/fluxinterval/fluxskipmult/(2*curtime+1)

    if(myrank.eq.1) then
       open(1004,file=spifile,status="replace",action="readwrite",&
            position="rewind",iostat=myiostat)
       call checkiostat(myiostat,"opening "//spifile)
       write(1004,*,iostat=myiostat)
       call checkiostat(myiostat,"writing "//spifile)
       write(1004,*) "# six columns."
       write(1004,*) "# Omega (column 1); |pulse ft|^2 (2); cross section (Mb) (3); flux (column 5)"
       write(1004,*)
       do i=-curtime,curtime
          wfi=(i+curtime)*Estep

!! LENGTH GAUGE WAS FT'ed multiply by wfi dont divide
!! NEVERMIND FACTOR OF 1/3
!!        myfac = 5.291772108d0**2 / 3d0 * 2d0 * PI / 1.37036d2 * wfi

!! WITH THIS FACTOR, NOW THE QUANTUM MECHANICAL CROSS SECTION IN 
!! MEGABARNS (10^-18 cm^2) IS IN COLUMN 3 REAL PART
          myfac = 5.291772108d0**2 * 2d0 * PI / 1.37036d2 * wfi

          write(1004,'(F8.4,100E18.6)',iostat=myiostat) wfi, pulseftsq(i), &
               FTgtau(i,:)/pulseftsq(i) * myfac, ftgtau(i,:)
       enddo
       call checkiostat(myiostat,"writing "//spifile)
       close(1004)
    endif


    call mpibarrier()
    OFLWR "  finished electronflux, stopping."; CFL
    call mpibarrier()
    call mpistop()

    deallocate(ftgtau,pulseft,pulseftsq)
    deallocate(bramo,ketmo,braavec,ketavec)
    deallocate(read_bramo,read_ketmo,read_braavec,read_ketavec)
    deallocate(gtau)
    deallocate(imke,impe,imV2,imyderiv,imkeop,impeop,imyop)
    deallocate(reke,repe,reV2,reyderiv,rekeop,repeop,reyop)
    deallocate(mobio,abio)

  contains

!! begin the flux matrix element and contraction routine section
!! flag=1, imag part; flag=2, real part; flag=0, all

    subroutine flux_op_onee(inspfs,keop,peop,flag)
      use parameters
      use orbgathersubmod
      use orbmultsubmod
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

! no, later
! scale correctly and clear memory
!         if(flag.ne.0) then
!            peop(:,lowspf:highspf)=peop(:,lowspf:highspf)*(-2d0)
!            keop(:,lowspf:highspf)=keop(:,lowspf:highspf)*(-2d0)
!         endif

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

! no, later
! scale correctly and clear memory
!         if(flag.ne.0) then
!            yop(:,lowspf:highspf)=yop(:,lowspf:highspf)*(-2d0)
!         endif

      endif

      if (parorbsplit.eq.1) then
         call mpiorbgather(yop,spfsize)
      endif

    end subroutine flux_op_nuc


    subroutine flux_op_twoe(mobra,moket,V2,flag)
      use parameters
      use orbgathersubmod
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


    function fluxeval(abra,aket,imke,impe,imV2,imyderiv,reke,repe,reV2,reyderiv)
      use parameters 
      use walkmod
      use configptrmod
      use sparseptrmod
      use sparsemultmod
      use asssubmod
      implicit none
      DATATYPE,intent(in) :: abra(numr,first_config:last_config),aket(numr,first_config:last_config),&
           imke(nspf,nspf),impe(nspf,nspf),imV2(nspf,nspf,nspf,nspf),imyderiv(nspf,nspf),&
           reke(nspf,nspf),repe(nspf,nspf),reV2(nspf,nspf,nspf,nspf),reyderiv(nspf,nspf)
      DATATYPE, allocatable :: multket(:,:), ketwork(:,:), conjgket(:,:)
      DATATYPE :: fluxeval, outsum
      type(CONFIGPTR) :: matrix_ptr
      type(SPARSEPTR) :: sparse_ptr
      integer :: ispf,jspf

!! imaginary part of electronic operators, real part of nuclear operators

      call configptralloc(matrix_ptr,www)
      matrix_ptr%kefac = 0d0
      matrix_ptr%constfac = 0d0

      allocate(multket(numr,first_config:last_config), ketwork(numr,first_config:last_config),&
           conjgket(numr,first_config:last_config))

      matrix_ptr%xpotmatel(:,:) = impe(:,:)
      matrix_ptr%xopmatel(:,:)  = imke(:,:)
      matrix_ptr%xymatel(:,:)   = imyderiv(:,:)
      matrix_ptr%xtwoematel(:,:,:,:) = imV2(:,:,:,:)

      if (sparseopt.ne.0) then
         call sparseptralloc(sparse_ptr,www)
         call assemble_sparsemats(www,matrix_ptr,sparse_ptr,1,min(nucfluxopt,1),0,0)
      endif

      call sparseconfigmult(www,aket,ketwork,matrix_ptr,sparse_ptr,1,min(nucfluxopt,1),0,0,0d0,-1)

      matrix_ptr%xpotmatel(:,:)      = ALLCON(matrix_ptr%xpotmatel(:,:))
      matrix_ptr%xopmatel(:,:)       = ALLCON(matrix_ptr%xopmatel(:,:))
      matrix_ptr%xymatel(:,:)        = ALLCON(matrix_ptr%xymatel(:,:))
      matrix_ptr%xtwoematel(:,:,:,:) = ALLCON(matrix_ptr%xtwoematel(:,:,:,:))

      if (sparseopt.ne.0) then
         call assemble_sparsemats(www,matrix_ptr,sparse_ptr,1,min(nucfluxopt,1),0,0)
      endif

      if (tot_adim.gt.0) then
         conjgket(:,:) = ALLCON(aket(:,:))
      endif
      call sparseconfigmult(www,conjgket,multket,matrix_ptr,sparse_ptr,1,min(nucfluxopt,1),0,0,0d0,-1)
      if (tot_adim.gt.0) then
         multket(:,:)=ALLCON(multket(:,:))
      endif
   
      outsum=0d0

      if (tot_adim.gt.0.and.nucfluxopt.ne.2) then
         ketwork(:,:)=(ketwork(:,:)+multket(:,:)) / (-1d0)    !! (-2)xREAL PART OF SPARSECONFIGMULT
         outsum = hermdot(abra,ketwork,tot_adim)
      endif

      if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0) then

!! real part of electronic operators, imaginary part of nuclear operators

         matrix_ptr%kefac = 1d0

         matrix_ptr%xpotmatel(:,:) = repe(:,:)
         matrix_ptr%xopmatel(:,:)  = reke(:,:)
         matrix_ptr%xymatel(:,:)   = reyderiv(:,:)
         matrix_ptr%xtwoematel(:,:,:,:) = reV2(:,:,:,:)

         if (sparseopt.ne.0) then
            call assemble_sparsemats(www,matrix_ptr,sparse_ptr,1,min(nucfluxopt,1),0,0)
         endif
         call sparseconfigmult(www,aket,ketwork,matrix_ptr,sparse_ptr,1,min(nucfluxopt,1),0,0,0d0,-1)

         matrix_ptr%xpotmatel(:,:)      = ALLCON(matrix_ptr%xpotmatel(:,:))
         matrix_ptr%xopmatel(:,:)       = ALLCON(matrix_ptr%xopmatel(:,:))
         matrix_ptr%xymatel(:,:)        = ALLCON(matrix_ptr%xymatel(:,:))
         matrix_ptr%xtwoematel(:,:,:,:) = ALLCON(matrix_ptr%xtwoematel(:,:,:,:))

         if (sparseopt.ne.0) then
            call assemble_sparsemats(www,matrix_ptr,sparse_ptr,1,min(nucfluxopt,1),0,0)
         endif

         if (tot_adim.gt.0) then
            conjgket(:,:) = ALLCON(aket(:,:))
         endif
         call sparseconfigmult(www,conjgket,multket,matrix_ptr,sparse_ptr,1,min(nucfluxopt,1),0,0,0d0,-1)
         if (tot_adim.gt.0) then
            multket(:,:)=ALLCON(multket(:,:))
         endif

         if (tot_adim.gt.0) then
            ketwork(:,:)=(ketwork(:,:)-multket(:,:)) / (0d0,-1d0)   !! (-2)xIMAG PART OF SPARSECONFIGMULT
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

  end subroutine fluxgtau0

end module fluxgtau0mod



subroutine fluxgtau(alg)
  use fluxgtau0mod
  use configmod
  implicit none
  integer,intent(in) :: alg

  call fluxgtau0(alg,www,bwwptr)

end subroutine fluxgtau


