
#include "Definitions.INC"

module ovlmod
  implicit none
  DATATYPE, allocatable :: orig_spfs(:,:,:), orig_avectors(:,:,:), overlaps(:,:,:)
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
  implicit none

  integer :: jnumovl, ifile,acomplex,spfcomplex,nstate,i,kk,tdims(3),tndof,tnumconfig,tnumr,tnspf
  external :: readavectorsubsimple

  if (numr.gt.1) then
     OFLWR "Need numr=1 for projone at this time TEMP CONTINUE.";CFL
  endif

!! read in the data from mcscf for our target cation state

  numovl=0
  do ifile=1,numovlfiles
     open(909,file=ovlspffiles(ifile),status="unknown",form="unformatted")
     open(910,file=ovlavectorfiles(ifile),status="unknown",form="unformatted")

     call avector_header_read_simple(910,nstate,tndof,tnumr,tnumconfig,acomplex)

     numovl=numovl+nstate

     call spf_header_read(909,tdims,tnspf,spfcomplex)

     if (tnspf.gt.nspf) then
        OFLWR " *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING *** "
        WRFL "  number of spfs for overlap states is LARGER than in calculation ", tnspf,nspf," REMOVING THOSE ORBITALS!!"; CFL
     endif
     if (tnumconfig.gt.numconfig)  then
        OFLWR " *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING *** "
        WRFL "  number of configs on overlap file greater than for calculation!"; CFL
     endif
     if(tnumr.gt.numr) then
        OFLWR "numr on disk too big. on file, current:",tnumr,numr; CFLST
     endif
     close(909);     close(910)
  enddo

  allocate(overlaps(numovl,0:autosize,mcscfnum),orig_spfs(spfsize,nspf,numovl),orig_avectors(numconfig,numr,numovl))
  orig_spfs=0d0;  orig_avectors=0d0
  jnumovl=0

!! no should be ok  call noparorbsupport("in ovl_initial")

  do ifile=1,numovlfiles
     open(909,file=ovlspffiles(ifile),status="unknown",form="unformatted")
     open(910,file=ovlavectorfiles(ifile),status="unknown",form="unformatted")
     call avector_header_read_simple(910,nstate,tndof,tnumr,tnumconfig,acomplex)
     call spf_header_read(909,tdims,tnspf,spfcomplex)
     call spf_read0(909,nspf,spfdims,tnspf,tdims,spfcomplex,spfdimtype,orig_spfs(:,:,jnumovl+1),(/0,0,0/))
     do i=tnspf+1,nspf
        call staticvector(orig_spfs(:,i,jnumovl+1),spfsize)
        if (parorbsplit.eq.3) then
           call gramschmidt(spfsize,i-1,spfsize,orig_spfs(:,:,jnumovl+1),orig_spfs(:,i,jnumovl+1),.true.)
        else
           call gramschmidt(spfsize,i-1,spfsize,orig_spfs(:,:,jnumovl+1),orig_spfs(:,i,jnumovl+1),.false.)
        endif
     enddo
     do kk=2,nstate
        orig_spfs(:,:,jnumovl+kk)=orig_spfs(:,:,jnumovl+1)
     enddo
     call load_avectors0(910,acomplex,orig_avectors(:,:,jnumovl+1),numr,numconfig,ndof,tnumr,tnumconfig,readavectorsubsimple, nstate)
     jnumovl=jnumovl+nstate
     
  enddo
  
end subroutine ovl_initial


subroutine getoverlaps(forceflag)
  use ovlmod
  use parameters
  use xxxmod
  implicit none
  integer ::  i,imc,forceflag

  calledflag = calledflag+1

  if (mod(calledflag-1,autosteps).eq.0) then
     do imc=1,mcscfnum
        do i=1,numovl
           call autocorrelate_one(yyy%cmfpsivec(astart(imc),0),yyy%cmfpsivec(spfstart,0),orig_spfs(:,:,i), &
                orig_avectors(:,:,i), overlaps(i,xcalledflag,imc),numr)
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
  implicit none

  integer :: k,molength,alength,nt,ketbat,imc
  real*8 :: piover2,dt,angle(mcscfnum)
  DATATYPE :: dot,myovl(mcscfnum) , bradot,phase,ketdot,blah
  DATATYPE :: bramo(spfsize,nspf),braavec(numconfig,numr,mcscfnum),ketmo(spfsize,nspf),ketavec(numconfig,numr,mcscfnum)

!!  dt=real(FluxInterval*FluxSkipMult,8)*par_timestep;  nt=floor(final time/dt)

  dt=real(FluxInterval*FluxSkipMult,8)*par_timestep;  nt=floor(real(numpropsteps,8)/fluxinterval/fluxskipmult)

  piover2=atan2(1d0,1d0)*2

!! initial setup
     
  inquire (iolength=molength) ketmo(:,:);  inquire (iolength=alength) ketavec(:,:,:)
  OFL
  write(mpifileptr,*) "MO record length is ",molength;  write(mpifileptr,*) "AVEC record length is ",alength
  CFL

  phase=(1d0,0d0)

  do ketbat=1,nt+1

    OFLWR "Reading ket batch ", ketbat, " of ", nt+1; CFL
    open(11001,file=fluxmofile2,status="old",form="unformatted",access="direct",recl=molength)
    open(11002,file=fluxafile2,status="old",form="unformatted",access="direct",recl=alength)
    open(1001,file=fluxmofile,status="old",form="unformatted",access="direct",recl=molength)
    open(1002,file=fluxafile,status="old",form="unformatted",access="direct",recl=alength)

    k=FluxSkipMult*(ketbat-1)+1
    read(1001,rec=k) ketmo(:,:) ;    read(1002,rec=k) ketavec(:,:,:) 
    read(11001,rec=k) bramo(:,:) ;    read(11002,rec=k) braavec(:,:,:) 

    close(1001);    close(1002);    close(11001);    close(11002)

    do imc=1,mcscfnum
       bradot=dot(braavec(:,:,imc),braavec(:,:,imc),numconfig*numr)
       ketdot=dot(ketavec(:,:,imc),ketavec(:,:,imc),numconfig*numr)
       
       call autocorrelate_one(braavec(:,:,imc),bramo,ketmo,ketavec(:,:,imc),myovl(imc),numr)
       
       blah=myovl(imc)/sqrt(bradot*ketdot)
       angle(imc)=acos(abs(blah))
       myovl(imc)=bradot+ketdot-myovl(imc)-CONJUGATE(myovl(imc))
    enddo

    OFL; write(mpifileptr,'(A30,1000F18.10)') "ERRDOT,ABSERRDOT,ANGLE T= ",dt*ketbat,myovl,(abs(myovl(imc)),angle(imc),imc=1,mcscfnum); CFL
 enddo

end subroutine wfnovl

