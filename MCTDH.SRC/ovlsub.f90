
#include "Definitions.INC"

module ovlmod
  implicit none
  DATATYPE, allocatable :: orig_spfs(:,:,:), orig_avectors(:,:), overlaps(:,:,:)
  integer :: numovl, calledflag=0 ,  xcalledflag=0 
end module ovlmod

subroutine staticvector(vector,size)
  implicit none
  integer :: i,size
  DATATYPE :: vector(size)
  real*8 :: nextran
!! don't init... keep as is
  do i=1,size
     vector(i)=nextran() + (0d0,1d0) * nextran()
  enddo
end subroutine staticvector

!! first draft; reads only one r point.   I.E. vib ovls not implemented.
!! numprop=1 only implemented

subroutine ovl_initial()
  use ovlmod
  use parameters
  use mpimod
  implicit none

  integer :: jnumovl, ifile,acomplex,spfcomplex,nstate,i,kk,tdims(3),tndof,tnumconfig,tnumr,tnspf
  external :: readavectorsubsimple
  DATATYPE, allocatable :: read_avectors(:,:), read_spfs(:,:)

  if (numr.gt.1) then
     OFLWR "Need numr=1 for projone at this time TEMP CONTINUE.";CFL
  endif

!! read in the data from mcscf for our target cation state

  numovl=0
  do ifile=1,numovlfiles
     if (myrank.eq.1) then
        open(909,file=ovlspffiles(ifile),status="unknown",form="unformatted")
        open(910,file=ovlavectorfiles(ifile),status="unknown",form="unformatted")

        call avector_header_read_simple(910,nstate,tndof,tnumr,tnumconfig,acomplex)
        call spf_header_read(909,tdims,tnspf,spfcomplex)
        close(909);     close(910)
     endif

     call mympiibcastone(nstate,1); call mympiibcastone(tnspf,1); call mympiibcastone(tnumconfig,1); 
     call mympiibcastone(tnumr,1)

     numovl=numovl+nstate

     if (tnspf.gt.nspf+numfrozen) then
        OFLWR " *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING *** "
        WRFL "  number of spfs for overlap states is LARGER than in calculation ", tnspf,nspf+numfrozen," REMOVING THOSE ORBITALS!!"; CFL
     endif
     if (tnumconfig.gt.num_config)  then
        OFLWR " *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING *** "
        WRFL "  number of configs on overlap file greater than for calculation!"; CFL
     endif
     if(tnumr.gt.numr) then
        OFLWR "numr on disk too big. on file, current:",tnumr,numr; CFLST
     endif
  enddo

  allocate(overlaps(numovl,0:autosize,mcscfnum),orig_spfs(spfsize,nspf,numovl),orig_avectors(tot_adim,numovl))

  allocate(read_spfs(spfsize,nspf+numfrozen))

  if (myrank.eq.1) then
     allocate(read_avectors(numr*num_config,numovl))
  else
     allocate(read_avectors(1,numovl))
  endif

  orig_spfs=0d0;  orig_avectors=0d0; overlaps=0d0; read_avectors=0d0
  jnumovl=0

  do ifile=1,numovlfiles

!!     call load_spfs0(orig_spfs(:,:,jnumovl+1), spfdims, nspf, spfdimtype, ovlspffiles(ifile), tnspf, (/0,0,0/))

     read_spfs(:,:)=0d0
     call load_spfs0(read_spfs(:,:), spfdims, nspf+numfrozen, spfdimtype, ovlspffiles(ifile), tnspf, (/0,0,0/))
     orig_spfs(:,:,jnumovl+1)=read_spfs(:,1+numfrozen:nspf+numfrozen)

     do i=tnspf-numfrozen+1,nspf
        call staticvector(orig_spfs(:,i,jnumovl+1),spfsize)
        if (parorbsplit.eq.3) then
           call gramschmidt(spfsize,i-1,spfsize,orig_spfs(:,:,jnumovl+1),orig_spfs(:,i,jnumovl+1),.true.)
        else
           call gramschmidt(spfsize,i-1,spfsize,orig_spfs(:,:,jnumovl+1),orig_spfs(:,i,jnumovl+1),.false.)
        endif
     enddo

     if (myrank.eq.1) then
        open(910,file=ovlavectorfiles(ifile),status="unknown",form="unformatted")
        call avector_header_read_simple(910,nstate,tndof,tnumr,tnumconfig,acomplex)
     endif
     call mympiibcastone(nstate,1); call mympiibcastone(tndof,1); call mympiibcastone(tnumr,1);
     call mympiibcastone(tnumconfig,1); call mympiibcastone(acomplex,1)
     if (myrank.eq.1) then
        call easy_load_avectors(910,acomplex,read_avectors(:,jnumovl+1),tnumr,tnumconfig,nstate)
        close(910)
     endif

     if (par_consplit.eq.0) then
        if (myrank.eq.1) then
           orig_avectors(:,jnumovl+1:jnumovl+nstate)=read_avectors(:,jnumovl+1:jnumovl+nstate)
        endif
        call mympibcast(orig_avectors(:,jnumovl+1),1,num_config*numr*nstate)
     else
        do kk=1,nstate
           call myscatterv(read_avectors(:,jnumovl+kk),orig_avectors(:,jnumovl+kk),configs_perproc(:)*numr)
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
  use ovlmod
  use parameters
  use configmod
  use xxxmod
  implicit none
  integer ::  i,imc,forceflag

  calledflag = calledflag+1

  if (mod(calledflag-1,autosteps).eq.0) then
     do imc=1,mcscfnum
        do i=1,numovl

              call autocorrelate_one(www,bioww,yyy%cmfpsivec(astart(imc),0),yyy%cmfpsivec(spfstart,0),orig_spfs(:,:,i), &
                   orig_avectors(:,i), overlaps(i,xcalledflag,imc),numr)

        enddo
        xcalledflag=xcalledflag+1
     enddo
  endif
  if (mod(calledflag-1,autosteps).eq.0.or.forceflag.ne.0) then
     open(881,file=outovl, status="unknown")
     do i=0,xcalledflag-1
        write(881,'(F12.3, 1000E20.10)') i*autotimestep,abs(overlaps(:,i,:))**2
     enddo
     close(881)

  endif

end subroutine getoverlaps
 
subroutine wfnovl()
  use parameters
  use configmod
  use mpimod
  implicit none
  
  integer :: k,molength,alength,nt,ketbat,imc,ispf
  real*8 :: piover2,dt,angle(mcscfnum)
  DATATYPE :: dot,myovl(mcscfnum) , bradot,phase,ketdot,blah
  DATATYPE :: bramo(spfsize,nspf),braavec(tot_adim,mcscfnum),ketmo(spfsize,nspf),ketavec(tot_adim,mcscfnum)
  DATATYPE, allocatable :: read_bramo(:,:), read_braavec(:,:), read_ketmo(:,:), read_ketavec(:,:)

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

!!  dt=real(FluxInterval*FluxSkipMult,8)*par_timestep;  nt=floor(final time/dt)

  dt=real(FluxInterval*FluxSkipMult,8)*par_timestep;  nt=floor(real(numpropsteps,8)/fluxinterval/fluxskipmult)

  piover2=atan2(1d0,1d0)*2

!! initial setup

  if (myrank.eq.1) then
     inquire (iolength=molength) read_ketmo(:,:);  inquire (iolength=alength) read_ketavec(:,:)
  endif
  call mympiibcastone(molength,1); call mympiibcastone(alength,1)

  OFL
  write(mpifileptr,*) "MO record length is ",molength;  write(mpifileptr,*) "AVEC record length is ",alength
  CFL

  phase=(1d0,0d0)

  do ketbat=1,nt+1
     
     OFLWR "Reading ket batch ", ketbat, " of ", nt+1; CFL
     if (myrank.eq.1) then
        open(11001,file=fluxmofile2,status="old",form="unformatted",access="direct",recl=molength)
        open(11002,file=fluxafile2,status="old",form="unformatted",access="direct",recl=alength)
        open(1001,file=fluxmofile,status="old",form="unformatted",access="direct",recl=molength)
        open(1002,file=fluxafile,status="old",form="unformatted",access="direct",recl=alength)
        
        k=FluxSkipMult*(ketbat-1)+1
        read(1001,rec=k) read_ketmo(:,:) ;    read(1002,rec=k) read_ketavec(:,:) 
        read(11001,rec=k) read_bramo(:,:) ;    read(11002,rec=k) read_braavec(:,:) 
        
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
        do imc=1,mcscfnum
           call myscatterv(read_braavec(:,imc),braavec(:,imc),configs_perproc(:)*numr)
           call myscatterv(read_ketavec(:,imc),ketavec(:,imc),configs_perproc(:)*numr)
        enddo
     endif


     do imc=1,mcscfnum
        bradot=dot(braavec(:,imc),braavec(:,imc),tot_adim)
        ketdot=dot(ketavec(:,imc),ketavec(:,imc),tot_adim)
        if (par_consplit.ne.0) then
           call mympireduceone(bradot); call mympireduceone(ketdot)
        endif
        
           call autocorrelate_one(www,bioww,braavec(:,imc),bramo,ketmo,ketavec(:,imc),myovl(imc),numr)
       
        blah=myovl(imc)/sqrt(bradot*ketdot)
        angle(imc)=acos(abs(blah))
        myovl(imc)=bradot+ketdot-myovl(imc)-CONJUGATE(myovl(imc))
     enddo
     
     OFL; write(mpifileptr,'(A30,1000F18.10)') "ERRDOT,ABSERRDOT,ANGLE T= ",dt*ketbat,myovl,(abs(myovl(imc)),angle(imc),imc=1,mcscfnum); CFL
  enddo

end subroutine wfnovl

