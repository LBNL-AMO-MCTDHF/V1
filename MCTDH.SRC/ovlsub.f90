
!! ALL ONE MODULE

!! ACTIONS 20 and 26 - overlaps 
!! (with time-dependent wfn, action 20, and between eigenfunctions, MCSCF-matel action 26)

#include "Definitions.INC"

module ovlsubmod
  use biorthotypemod
  implicit none
  DATATYPE, allocatable, private :: orig_spfs(:,:,:), orig_avectors(:,:), overlaps(:,:,:)
  integer, private :: numovl, calledflag=0 ,  xcalledflag=0 
  integer,private :: ovlbioalloc=0
  type(biorthotype),target,allocatable,private :: ovlbiovar(:,:)

contains

!! first draft; reads only one r point.   I.E. vib ovls not implemented.
!! numprop=1 only implemented

subroutine ovl_initial()
  use parameters
  use mpimod
  use readavectormod
  use mpisubmod
  use utilmod
  use configloadmod
  use loadstuffmod
  implicit none
  integer :: jnumovl, ifile,acomplex,spfcomplex,nstate,i,kk,tdims(3),&
       tnum2part,tnumconfig,tnumr,tnspf,myiostat
  DATATYPE, allocatable :: read_avectors(:,:), read_spfs(:,:)
  DATATYPE :: nullvector(numr)

  nullvector(:)=0

!! read in the data from mcscf for our target cation state

  numovl=0
  do ifile=1,numovlfiles
     if (myrank.eq.1) then
        open(909,file=ovlspffiles(ifile),status="unknown",&
             form="unformatted",iostat=myiostat)
        call checkiostat(myiostat,"opening "//ovlspffiles(ifile))
        open(910,file=ovlavectorfiles(ifile),status="unknown",&
             form="unformatted",iostat=myiostat)
        call checkiostat(myiostat,"opening "//ovlavectorfiles(ifile))
        call avector_header_read_simple(910,nstate,tnum2part,tnumr,tnumconfig,acomplex)
        call spf_header_read(909,tdims,tnspf,spfcomplex)
        close(909);     close(910)
     endif

     call mympiibcastone(nstate,1)
     call mympiibcastone(tnspf,1); call mympiibcastone(tnumconfig,1); 
     call mympiibcastone(tnumr,1)

     numovl=numovl+nstate

     if (tnspf.gt.nspf+numfrozen) then
        OFLWR " *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING *** "
        WRFL "  number of spfs for overlap states is LARGER than in calculation ", &
             tnspf,nspf+numfrozen," REMOVING THOSE ORBITALS!!"; CFL
     endif
     if (tnumconfig.gt.num_config)  then
        OFLWR " *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING *** "
        WRFL "  number of configs on overlap file greater than for calculation!"; CFL
     endif
     if(tnumr.gt.numr) then
        OFLWR "numr on disk too big. on file, current:",tnumr,numr; CFLST
     endif
  enddo

  allocate(overlaps(numovl,0:autosize,mcscfnum),&
       orig_spfs(spfsize,nspf,numovl),orig_avectors(tot_adim,numovl))
  overlaps=0; orig_spfs=0
  if (tot_adim.gt.0) then
     orig_avectors=0
  endif
  allocate(read_spfs(spfsize,nspf+numfrozen));    read_spfs=0

  if (myrank.eq.1) then
     allocate(read_avectors(numr*num_config,numovl))
  else
     allocate(read_avectors(1,numovl))
  endif
  read_avectors=0d0

  jnumovl=0

  do ifile=1,numovlfiles

     read_spfs(:,:)=0d0
     call load_spfs0(read_spfs(:,:), spfdims, nspf+numfrozen, &
          spfdimtype, ovlspffiles(ifile), tnspf, (/0,0,0/))
     orig_spfs(:,:,jnumovl+1)=read_spfs(:,1+numfrozen:nspf+numfrozen)

     do i=tnspf-numfrozen+1,nspf
        call staticvector(orig_spfs(:,i,jnumovl+1),spfsize)
        if (parorbsplit.eq.3) then
           call gramschmidt(spfsize,i-1,spfsize,&
                orig_spfs(:,:,jnumovl+1),orig_spfs(:,i,jnumovl+1),.true.)
        else
           call gramschmidt(spfsize,i-1,spfsize,&
                orig_spfs(:,:,jnumovl+1),orig_spfs(:,i,jnumovl+1),.false.)
        endif
     enddo

     if (myrank.eq.1) then
        open(910,file=ovlavectorfiles(ifile),status="unknown",&
             form="unformatted",iostat=myiostat)
        call checkiostat(myiostat,"opening "//ovlavectorfiles(ifile))
        call avector_header_read_simple(910,nstate,tnum2part,tnumr,tnumconfig,acomplex)
     endif
     call mympiibcastone(nstate,1); call mympiibcastone(tnum2part,1); 
     call mympiibcastone(tnumr,1);
     call mympiibcastone(tnumconfig,1); call mympiibcastone(acomplex,1)
     if (myrank.eq.1) then
        call easy_load_avectors(910,acomplex,read_avectors(:,jnumovl+1),&
             tnumr,tnumconfig,nstate)
        close(910)
     endif

     if (par_consplit.eq.0) then
        if (myrank.eq.1) then
           orig_avectors(:,jnumovl+1:jnumovl+nstate) = &
                read_avectors(:,jnumovl+1:jnumovl+nstate)
        endif
        call mympibcast(orig_avectors(:,jnumovl+1),1,num_config*numr*nstate)
     else
        do kk=1,nstate
           if (tot_adim.gt.0) then
              call myscatterv(read_avectors(:,jnumovl+kk),&
                   orig_avectors(:,jnumovl+kk),configs_perproc(:)*numr)
           else
              call myscatterv(read_avectors(:,jnumovl+kk),&
                   nullvector(:),configs_perproc(:)*numr)
           endif
        enddo
     endif

     do kk=2,nstate
        orig_spfs(:,:,jnumovl+kk)=orig_spfs(:,:,jnumovl+1)
     enddo

     jnumovl=jnumovl+nstate
     
  enddo
  
  deallocate(read_avectors, read_spfs)

end subroutine ovl_initial



subroutine getoverlaps(forceflag)
  use parameters
  use configmod
  use xxxmod
  use mpimod   !! myrank
  use autocorrelate_one_mod
  implicit none
  integer ::  i,imc,forceflag,myiostat
  DATATYPE :: nullvector1(numr),nullvector2(numr)

  nullvector1(:)=0;   nullvector2(:)=0

  calledflag = calledflag+1

  if (ovlbioalloc.eq.0) then
     allocate(ovlbiovar(numovl,mcscfnum))
  endif
  ovlbioalloc=1
     
  if (mod(calledflag-1,autosteps).eq.0) then
     do imc=1,mcscfnum
        do i=1,numovl
           if (tot_adim.gt.0) then
              call autocorrelate_one(www,bioww,yyy%cmfavec(:,imc,0),&
                   yyy%cmfspfs(:,0),orig_spfs(:,:,i), &
                   orig_avectors(:,i), overlaps(i,xcalledflag,imc),numr,ovlbiovar(i,imc))
           else
              call autocorrelate_one(www,bioww,nullvector1(:),&
                   yyy%cmfspfs(:,0),orig_spfs(:,:,i), &
                   nullvector2(:), overlaps(i,xcalledflag,imc),numr,ovlbiovar(i,imc))
           endif
        enddo
        xcalledflag=xcalledflag+1
     enddo
  endif
  if (mod(calledflag-1,autosteps).eq.0.or.forceflag.ne.0.and.myrank.eq.1) then
     open(881,file=outovl, status="unknown",iostat=myiostat)
     call checkiostat(myiostat,"opening "//outovl)
     do i=0,xcalledflag-1
        write(881,'(F12.3, 2000E20.10)',iostat=myiostat) &
             i*autotimestep,abs(overlaps(:,i,:))**2,overlaps(:,i,:)
     enddo
     call checkiostat(myiostat,"writing "//outovl)
     close(881)
  endif

end subroutine getoverlaps


subroutine mcscf_matel()
  use parameters
  use configmod
  use xxxmod
  use mpimod   !! myrank
  use biorthotypemod
  use autocorrelate_one_mod
  implicit none
  integer ::  i,j,myiostat
  DATATYPE :: myovl(numovl,numovl),nullvector1(numr),nullvector2(numr)
  type(biorthotype),target :: mcbiovar(numovl,numovl)

  nullvector1(:)=0;nullvector2(:)=0

!! REPLACE THIS - ADAPT FINALSTATS0.
!!   Make dual-purpose subroutine (finalstats00) contatining the guts of finalstats0 
!!   and call finalstats00 here and in finalstats0.  For now just overlaps

  do j=1,numovl
     do i=1,numovl
        if (tot_adim.gt.0) then
           call autocorrelate_one(www,bioww, orig_avectors(:,i), orig_spfs(:,:,i), &
                orig_spfs(:,:,j), orig_avectors(:,j), myovl(i,j), numr, mcbiovar(i,j))
        else
           call autocorrelate_one(www,bioww, nullvector1(:), orig_spfs(:,:,i), &
                orig_spfs(:,:,j), nullvector2(:), myovl(i,j), numr, mcbiovar(i,j))
        endif
     enddo
  enddo

  if (myrank.eq.1) then
     open(881,file=outmatel, status="unknown",iostat=myiostat)
     call checkiostat(myiostat,"opening "//outmatel)
     do i=1,numovl
        write(881,'(1000F17.10)',iostat=myiostat) abs(myovl(i,:))**2, myovl(i,:)
     enddo
     call checkiostat(myiostat,"writing "//outmatel)
     close(881)
  endif

  call mpibarrier()
  OFLWR "OK done mcscf_matel, stopping"; CFLST

end subroutine mcscf_matel

 
subroutine wfnovl()
  use parameters
  use configmod
  use mpimod
  use biorthotypemod
  use autocorrelate_one_mod
  use mpisubmod
  implicit none
  type(biorthotype),target :: wfnbiovar(mcscfnum)  
  integer :: k,molength,alength,nt,ketbat,imc,ispf,myiostat
  real*8 :: piover2,dt,angle(mcscfnum)
  DATATYPE :: myovl(mcscfnum) , bradot,phase,ketdot,blah,nullvector1(numr),&
       nullvector2(numr)
  DATATYPE, allocatable :: read_bramo(:,:), read_braavec(:,:), &
       read_ketmo(:,:), read_ketavec(:,:),&
       bramo(:,:),braavec(:,:),ketmo(:,:),ketavec(:,:)

  nullvector1(:)=0;nullvector2(:)=0
  allocate(bramo(spfsize,nspf),braavec(tot_adim,mcscfnum),&
       ketmo(spfsize,nspf),ketavec(tot_adim,mcscfnum))
  bramo=0; ketmo=0

  if (tot_adim.gt.0) then
     braavec=0; ketavec=0
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
  if (myrank.eq.1) then
     allocate(read_braavec(numr*num_config,mcscfnum),read_ketavec(numr*num_config,mcscfnum))
  else
     allocate(read_braavec(1,mcscfnum),read_ketavec(1,mcscfnum))
  endif
  read_bramo=0; read_ketmo=0; read_braavec=0; read_ketavec=0

!!  dt=real(FluxInterval*FluxSkipMult,8)*par_timestep;  nt=floor(final time/dt)

  dt=real(FluxInterval*FluxSkipMult,8)*par_timestep;  
  nt=floor(real(numpropsteps,8)/fluxinterval/fluxskipmult)

  piover2=atan2(1d0,1d0)*2

!! initial setup

  if (myrank.eq.1) then
     inquire (iolength=molength) read_ketmo(:,:)
     inquire (iolength=alength) read_ketavec(:,:)
  endif
  call mympiibcastone(molength,1); call mympiibcastone(alength,1)

  OFL
  write(mpifileptr,*) "MO record length is ",molength
  write(mpifileptr,*) "AVEC record length is ",alength
  CFL

  phase=(1d0,0d0)

  do ketbat=1,nt+1
     
     OFLWR "Reading ket batch ", ketbat, " of ", nt+1; CFL
     if (myrank.eq.1) then
        open(11001,file=fluxmofile2,status="old",form="unformatted",&
             access="direct",recl=molength,iostat=myiostat)
        call checkiostat(myiostat,"opening "//fluxmofile2)
        open(11002,file=fluxafile2,status="old",form="unformatted",&
             access="direct",recl=alength,iostat=myiostat)
        call checkiostat(myiostat,"opening "//fluxafile2)
        open(1001,file=fluxmofile,status="old",form="unformatted",&
             access="direct",recl=molength,iostat=myiostat)
        call checkiostat(myiostat,"opening "//fluxmofile)
        open(1002,file=fluxafile,status="old",form="unformatted",&
             access="direct",recl=alength,iostat=myiostat)
        call checkiostat(myiostat,"opening "//fluxafile)
        
        k=FluxSkipMult*(ketbat-1)+1
        read(1001,rec=k,iostat=myiostat) read_ketmo(:,:) 
        call checkiostat(myiostat,"reading "//fluxmofile2)
        read(1002,rec=k,iostat=myiostat) read_ketavec(:,:) 
        call checkiostat(myiostat,"reading "//fluxafile2)
        read(11001,rec=k,iostat=myiostat) read_bramo(:,:) 
        call checkiostat(myiostat,"reading "//fluxmofile)
        read(11002,rec=k,iostat=myiostat) read_braavec(:,:) 
        call checkiostat(myiostat,"reading "//fluxafile)
        
        close(1001);    close(1002);    close(11001);    close(11002)
     endif
     
     if (parorbsplit.ne.3) then
        if (myrank.eq.1) then
           bramo(:,:)=read_bramo(:,:)
           ketmo(:,:)=read_ketmo(:,:)
        endif
        call mympibcast(bramo(:,:),1,totspfdim)
        call mympibcast(ketmo(:,:),1,totspfdim)
     else
        do ispf=1,nspf
           call splitscatterv(read_bramo(:,ispf),bramo(:,ispf))
           call splitscatterv(read_ketmo(:,ispf),ketmo(:,ispf))
        enddo
     endif
     if (par_consplit.eq.0) then
        if (myrank.eq.1) then
           braavec(:,:)=read_braavec(:,:)
           ketavec(:,:)=read_ketavec(:,:)
        endif
        call mympibcast(braavec(:,:),1,numr*num_config*mcscfnum)
        call mympibcast(ketavec(:,:),1,numr*num_config*mcscfnum)
     else
        if (tot_adim.gt.0) then
           do imc=1,mcscfnum
              call myscatterv(read_braavec(:,imc),braavec(:,imc),configs_perproc(:)*numr)
              call myscatterv(read_ketavec(:,imc),ketavec(:,imc),configs_perproc(:)*numr)
           enddo
        else
           do imc=1,mcscfnum
              call myscatterv(read_braavec(:,imc),nullvector1(:),configs_perproc(:)*numr)
              call myscatterv(read_ketavec(:,imc),nullvector2(:),configs_perproc(:)*numr)
           enddo
        endif
     endif


     do imc=1,mcscfnum
        bradot=0; ketdot=0
        if (tot_adim.gt.0) then
           bradot=dot(braavec(:,imc),braavec(:,imc),tot_adim)
           ketdot=dot(ketavec(:,imc),ketavec(:,imc),tot_adim)
        endif
        if (par_consplit.ne.0) then
           call mympireduceone(bradot); call mympireduceone(ketdot)
        endif
        
        if (tot_adim.gt.0) then
           call autocorrelate_one(www,bioww,braavec(:,imc),bramo,ketmo,&
                ketavec(:,imc),myovl(imc),numr,wfnbiovar(imc))
        else
           call autocorrelate_one(www,bioww,nullvector1(:),bramo,ketmo,&
                nullvector2(:),myovl(imc),numr,wfnbiovar(imc))
        endif

        blah=myovl(imc)/sqrt(bradot*ketdot)
        angle(imc)=acos(abs(blah))
        myovl(imc)=bradot+ketdot-myovl(imc)-CONJUGATE(myovl(imc))
     enddo
     
     OFL; write(mpifileptr,'(A30,1000F18.10)') "ERRDOT,ABSERRDOT,ANGLE T= ",&
          dt*ketbat,myovl,(abs(myovl(imc)),angle(imc),imc=1,mcscfnum); CFL
  enddo

  deallocate(bramo,braavec,ketmo,ketavec)

  call mpibarrier()
  OFLWR "OK done wfnovl, stopping"; CFLST

end subroutine wfnovl

end module ovlsubmod


