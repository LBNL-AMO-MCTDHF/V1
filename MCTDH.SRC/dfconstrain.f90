
#include "Definitions.INC"

!! ALL MODULES

!! SUBROUTINES FOR RESTRICTED CONFIG LIST (dfrestrictflag>0)
!! and constraint (constraintflag>0, for proper propagation according to variational
!! principle when using restricted configuration spaces
!!
!!   if constraintflag=2, dirac-frenkel.
!!   if constraintflag=1, density matrix.
!!   
!!  g_ij = < phi_i| d/dt | phi_j >    , that's variable conmatel, should be hermitian
!!
!! where i is virtual ("excluded") and j is active ("included") (other g_ij undefined, 0=0)
!! hermitian g matrix can then be immediattely gotten by fiat
!!

module dfconsubmod
contains

subroutine get_constraint(time)
  use fileptrmod
  use ham_parameters
  implicit none
  real*8,intent(in) :: time

  if (constraintflag.eq.1) then
     call get_denconstraint(time)
  else if (constraintflag.eq.2) then
  if (drivingflag.ne.0) then
     OFLWR "Driving with dfconstraint not implemented yet"; CFLST
  endif
     call get_dfconstraint(time)
  else 
     OFLWR "CONSTRAINTFLAG ERROR  ", constraintflag; CFLST
  endif

end subroutine get_constraint


subroutine dferror(www,cptr,sptr,avector,numvects,outerror,time)
  use configptrmod
  use sparseptrmod
  use walkmod
  use r_parameters
  use fileptrmod
  use dotmod
  use sparsemultmod
  use mpisubmod
  implicit none
  integer,intent(in) :: numvects
  type(walktype),intent(in) :: www
  type(configptr),intent(in) :: cptr
  type(sparseptr),intent(in) :: sptr
  real*8,intent(in) :: time
  DATATYPE,intent(in) :: avector(numr,www%firstconfig:www%lastconfig,numvects)
  CNORMTYPE,intent(out) :: outerror
  DATATYPE,allocatable :: temp1(:,:),temp2(:,:)
  integer :: i,imc

  allocate(temp1(numr,www%firstconfig:www%lastconfig), &
       temp2(numr,www%firstconfig:www%lastconfig))
  if (www%lastconfig.ge.www%firstconfig) then
     temp1=0; temp2=0
  endif

  outerror=0d0
  do imc=1,numvects
     if (www%numdfconfigs.gt.0) then

        temp1(:,:)=avector(:,:,imc)
        call sparseconfigmult(www,temp1,temp2,cptr,sptr,1,1,1,1,time,imc)
        do i=www%firstconfig,www%lastconfig
           temp1(:,i)=temp2(:,i)*(1-www%ddd%dfincludedmask(i))
        enddo
        
!! reporting norm of error appropriate to the variational principle.

        outerror = outerror + dot(temp1,temp1,www%totadim) !! ok conversion
     else if (www%numdfconfigs.lt.0) then
        OFLWR "error, dfconstrain not allocated"; CFLST
     endif
  enddo
  if (www%parconsplit.ne.0) then
#ifndef REALGO
#ifndef CNORMFLAG
     call mympirealreduceone(outerror)
#else
     call mympireduceone(outerror)
#endif
#else
     call mympireduceone(outerror)
#endif
  endif

  deallocate(temp1,temp2)

end subroutine dferror


subroutine get_smallwalkvects(www,avector, smallwalkvects,nblock,howmany)
  use fileptrmod
  use walkmod
  use mpisubmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: nblock,howmany
  DATATYPE,intent(in) ::  avector(nblock,www%firstconfig:www%lastconfig,howmany)
  DATATYPE,intent(out) :: smallwalkvects(nblock,www%configstart:www%configend,&
       howmany,www%nspf,www%nspf)
  DATATYPE, allocatable ::  bigavector(:,:,:)
  integer ::    i,ii,ix

  if (www%dfrestrictflag.le.www%dflevel) then
     OFLWR "WTF DFRESTRICT IS  ", www%dfrestrictflag,www%dflevel; CFLST
  endif
  
!! DO SUMMA (parconsplit.ne.0 and sparsesummaflag.eq.2, "circ")

  allocate(bigavector(nblock,www%numconfig,howmany))
  bigavector(:,:,:)=0d0
  if (www%lastconfig.ge.www%firstconfig) then
     bigavector(:,www%firstconfig:www%lastconfig,:)=avector(:,:,:)
  endif
  if (www%parconsplit.ne.0) then
     do ii=1,howmany
        call mpiallgather(bigavector(:,:,ii),nblock*www%numconfig,&
             nblock*www%configsperproc(:),nblock*www%maxconfigsperproc)
     enddo
  endif
  if (www%configend.ge.www%configstart) then
     smallwalkvects(:,:,:,:,:)=0d0
  endif
  do i=1,www%ddd%numdfwalks

     ii=www%ddd%includedorb(i);     ix=www%ddd%excludedorb(i)

     smallwalkvects(:,www%ddd%dfwalkto(i),:,ii,ix) = &
          smallwalkvects(:,www%ddd%dfwalkto(i),:,ii,ix) + &
          bigavector(:,www%ddd%dfwalkfrom(i),:) * www%ddd%dfwalkphase(i)
  enddo

  deallocate(bigavector)

end subroutine get_smallwalkvects



subroutine get_rhomat(www,avector, rhomat,nblock,howmany)
  use fileptrmod
  use sparse_parameters
  use spfsize_parameters !! parorbsplit
  use orbgathersubmod
  use walkmod
  use mpisubmod
  use utilmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: nblock,howmany
  DATATYPE,intent(in) ::  avector(nblock,www%firstconfig:www%lastconfig,howmany)
  DATATYPE,intent(out) :: rhomat(www%nspf,www%nspf,www%nspf,www%nspf)
!! nblock,configstart:configend,howmany,alpha,beta  alpha=included beta=excluded  :
  DATATYPE, allocatable ::  smallwalkvects(:,:,:,:,:)   
  integer :: ii,lowspf,highspf,numspf,flag

  if (www%dfrestrictflag.le.www%dflevel) then
     OFLWR "WTF DFRESTRICT IS  ", www%dfrestrictflag,www%dflevel; CFLST
  endif
  
  allocate(smallwalkvects(nblock,www%configstart:www%configend,&
       howmany,www%nspf,www%nspf))
  if (www%configend.ge.www%configstart) then
     smallwalkvects=0
  endif

!! call always
  call get_smallwalkvects(www,avector,smallwalkvects,nblock,howmany)

!! THIS TAKES A LONG TIME.

  lowspf=1; highspf=www%nspf
  if (parorbsplit.eq.1.and.sparseconfigflag.eq.0) then
     call checkorbsetrange(www%nspf,flag)
     if (flag.ne.0) then
        OFLWR "error exit, can't do get_rhomat with ",www%nspf,&
             " orbitals sparse mode"; CFLST
     endif
     call getOrbSetRange(lowspf,highspf)
  endif
  numspf=highspf-lowspf+1
  
  ii=(www%configend-www%configstart+1)*nblock*howmany
  if (ii.gt.0.and.numspf.gt.0) then
     call MYGEMM(CNORMCHAR,'N',www%nspf**2,www%nspf*numspf,ii,DATAONE,&
          smallwalkvects,ii,smallwalkvects(:,:,:,:,lowspf),ii,DATAZERO,&
          rhomat(:,:,:,lowspf),www%nspf**2)
  else if (numspf.gt.0) then
     rhomat(:,:,:,lowspf:highspf)=0d0
  endif
  if (parorbsplit.eq.1.and.sparseconfigflag.eq.0) then
     call mpiorbgather(rhomat,www%nspf**3)
  endif
!! BAD REDUCE.  REDUCE BEFORE ORBGATHER (need new communicators)
  if (sparseconfigflag.ne.0) then
     call mympireduce(rhomat,www%nspf**4) 
  endif

  deallocate(smallwalkvects)

end subroutine get_rhomat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   DIRAC FRENKEL (constraintflag=2) !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine get_dfconstraint0(inavectors,numvects,cptr,sptr,www,time)
  use fileptrmod
  use timing_parameters
  use sparse_parameters
  use spfsize_parameters   !! parorbsplit
  use basis_parameters
  use r_parameters
  use constraint_parameters
  use ham_parameters
  use configptrmod
  use sparseptrmod
  use walkmod
  use mpimod
  use dotmod
  use sparsemultmod
  use basissubmod
  use invsubmod
  use orbgathersubmod
  use mpisubmod
  use utilmod
  use clockmod
  implicit none
  integer,intent(in) :: numvects
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: inavectors(numr,www%firstconfig:www%lastconfig,numvects)
  type(CONFIGPTR),intent(inout) :: cptr
  type(SPARSEPTR),intent(in) :: sptr
  DATATYPE :: tempmatel(www%nspf,www%nspf)
  integer, save :: times(20)=0, icalled=0
  integer ::    i,     j,  ii, lwork,isize,ishell,iiyy,maxii,imc,itime,jtime,getlen,&
       lowspf,highspf,numspf,flag,myiostat
  integer :: ipairs(2,www%nspf*(www%nspf-1))
  DATATYPE, allocatable :: avectorp(:,:), rhs(:,:), avector(:,:), rhomat(:,:,:,:), &
!! numconfig,numr,alpha,beta  alpha=included beta=excluded
       smallwalkvects(:,:,:,:)  
  real*8, allocatable :: rhomatpairs(:,:,:,:),rhspairs(:,:), rhspairstemp(:,:), &
       rhspairsbig(:,:), projector(:,:,:,:,:), realrhomat(:,:,:,:,:,:), &
       bigprojector(:,:,:,:,:), rhomatpairsbig(:,:,:,:), &
       realrhs(:,:,:), sing(:), rwork(:), rhomatpairscopy(:,:,:,:), &
       rhomatpairsbigcopy(:,:,:,:), rhomatpairscopy2(:,:,:,:), rhspairstemp2(:,:),&
       mattemp(:,:),pseudoinv(:,:), rhspairsbigtemp(:,:), temppairs(:,:), &
       temppairsbig(:,:)
  complex*16, allocatable :: work(:)
#ifdef REALGO
  integer, parameter :: zzz=1
#else
  integer, parameter :: zzz=2
#endif
  real*8 :: time, myconjg(zzz,zzz), oneone(zzz,zzz), iimag(zzz,zzz), myperm(zzz,zzz)
  real*8 :: rrcond

  if (drivingflag.ne.0) then
     OFLWR "Driving with dfconstraint not implemented yet"; CFLST
  endif

#ifdef REALGO
  myconjg(:,:) = 1;  oneone(:,:) = 1;  iimag(:,:)=0; myperm(:,:)=1

#else
  myconjg(:,:)=RESHAPE((/1,0,0,-1/),(/2,2/));  oneone(:,:)=RESHAPE((/1,0,0,1/),(/2,2/))
  iimag(:,:)=RESHAPE((/0,1,-1,0/),(/2,2/));  myperm(:,:)=RESHAPE((/0,1, 1,0/),(/2,2/))

#endif

  call myclock(itime)

  if (www%dfrestrictflag.le.www%dflevel) then
     OFLWR "WTF DFRESTRICT IS  ", www%dfrestrictflag,www%dflevel; CFLST
  endif

  icalled=icalled+1

  if (notiming.eq.0.and.myrank.eq.1) then
     if (icalled.eq.1) then
        open(8712,file=timingdir(1:getlen(timingdir))//"/dfconstrain.dat",&
             status="unknown",iostat=myiostat)
        call checkiostat(myiostat,"opening dfconstrain timing file")
        write(8712,*,iostat=myiostat) "DF timings."
        call checkiostat(myiostat,"writing dfconstrain timing file")
        close(8712)
     else if (mod(icalled,30).eq.0) then
        open(8712,file=timingdir(1:getlen(timingdir))//"/dfconstrain.dat",&
             status="old", position="append",iostat=myiostat)
        call checkiostat(myiostat,"opening dfconstrain timing file")
        write(8712,'(100I12)',iostat=myiostat) times(1:14)/1000
        call checkiostat(myiostat,"writing dfconstrain timing file")
        close(8712)
     endif
  endif

  cptr%xconmatel(:,:)=0.d0;   cptr%xconmatelxx(:,:)=0.d0;   
  cptr%xconmatelyy(:,:)=0.d0;   cptr%xconmatelzz(:,:)=0.d0

  allocate(avectorp(numr,www%firstconfig:www%lastconfig), &
       avector(numr,www%firstconfig:www%lastconfig), &
       smallwalkvects(numr,www%configstart:www%configend,www%nspf,www%nspf), &
       rhs(www%nspf,www%nspf),  rhomat(www%nspf,www%nspf,www%nspf,www%nspf))
  if (www%lastconfig.ge.www%firstconfig) then
     avectorp=0; avector=0; 
  endif
  if (www%configend.ge.www%configstart) then
     smallwalkvects=0;
  endif
  rhs=0; rhomat=0

  isize=0
  do ishell=1,numshells-1
     do j=allshelltop(ishell)+1,www%nspf
        do i=allshelltop(ishell-1)+1,allshelltop(ishell)
           isize=isize+1
           ipairs(:,isize)=(/ i,j /)
        enddo
     enddo
  enddo
     
  allocate(projector(zzz,isize,zzz,www%nspf,www%nspf), &
       bigprojector(zzz,2*isize,zzz,www%nspf,www%nspf), &
       mattemp(zzz*isize,zzz*isize),pseudoinv(zzz*isize,2*zzz*isize))
  projector=0d0;   bigprojector=0d0;; mattemp=0; pseudoinv=0


  do ii=1,isize
     projector(:,ii,:,ipairs(1,ii),ipairs(2,ii))=oneone(:,:)
     projector(:,ii,:,ipairs(2,ii),ipairs(1,ii))=(-1)*myconjg(:,:)

     if (conway.ne.3) then
        bigprojector(:,ii,:,ipairs(1,ii),ipairs(2,ii))=oneone(:,:)
        bigprojector(:,ii+isize,:,ipairs(2,ii),ipairs(1,ii))=oneone(:,:)
     else
        bigprojector(:,ii,:,ipairs(1,ii),ipairs(2,ii))=oneone(:,:)        !! lagrangian
        bigprojector(:,ii,:,ipairs(2,ii),ipairs(1,ii))=(-1)*myconjg(:,:)  
        bigprojector(:,ii+isize,:,ipairs(1,ii),ipairs(2,ii))=oneone(:,:)*conprop  !!mclachlan
        bigprojector(:,ii+isize,:,ipairs(2,ii),ipairs(1,ii))=myconjg(:,:)*conprop
     endif

  enddo

  allocate(realrhs(zzz,www%nspf,www%nspf),&
       realrhomat(zzz,www%nspf,www%nspf,zzz,www%nspf,www%nspf), &
       rhomatpairs(zzz,isize,zzz,isize),rhspairs(zzz,isize),&
       rhspairstemp(zzz,isize), rhomatpairscopy(zzz,isize,zzz,isize), &
       rhomatpairscopy2(zzz,isize,zzz,isize),rhspairstemp2(zzz,isize), &
       rhomatpairsbig(zzz,2*isize,zzz,isize),rhspairsbig(zzz,2*isize),&
       rhspairsbigtemp(zzz,2*isize),rhomatpairsbigcopy(zzz,2*isize,zzz,isize),&
       temppairs(zzz*isize,zzz*www%nspf**2),temppairsbig(2*zzz*isize,zzz*www%nspf**2))
  realrhs=0; realrhomat=0; rhomatpairs=0; rhspairstemp=0;
  rhomatpairscopy=0; rhomatpairscopy2=0; rhspairstemp2=0; rhomatpairsbig=0;
  rhomatpairsbig=0; rhspairsbig=0; rhspairsbigtemp=0; rhomatpairsbigcopy=0;
  temppairs=0; temppairsbig=0

  lwork=100*isize*zzz
  allocate(sing(2*isize*zzz),work(lwork), rwork(lwork))
  sing=0; work=0; rwork=0

  rrcond=lioreg

  call myclock(jtime);   times(1)=times(1)+jtime-itime

  lowspf=1; highspf=www%nspf
  if (parorbsplit.eq.1.and.sparseconfigflag.eq.0) then
     call checkorbsetrange(www%nspf,flag)
     if (flag.ne.0) then
        OFLWR "error exit, can't do get_dfconstraint0 with ",www%nspf,&
             " orbitals sparse mode"; CFLST
     endif
     call getOrbSetRange(lowspf,highspf)
  endif
  numspf=highspf-lowspf+1

  maxii=1
  if (tdflag.ne.0) then
     maxii=4
  endif

  do iiyy=1,maxii

     call myclock(itime)

     rhomatpairscopy(:,:,:,:)=0d0
     rhspairstemp(:,:)=0d0
     rhomatpairsbigcopy(:,:,:,:)=0d0
     rhspairsbigtemp(:,:)=0d0

     call myclock(jtime);     times(1)=times(1)+jtime-itime

     do imc=1,numvects
        
        call myclock(itime)
        
        avector(:,:)=inavectors(:,:,imc)

        call basis_project(www,numr,avector)

        call get_smallwalkvects(www,avector,smallwalkvects,numr,1)
        
        call myclock(jtime);        times(3)=times(3)+jtime-itime;     itime=jtime

!! THIS TAKES A LONG TIME.
!! should only need to do off diagonal blocks but doing all to check
!!   ...no for general case do all I think (Except for diag, whatever)

        ii=(www%configend-www%configstart+1)*numr
        if (ii.gt.0.and.numspf.gt.0) then
           call MYGEMM(CNORMCHAR,'N',www%nspf**2,www%nspf*numspf,ii,DATAONE,&
                smallwalkvects,ii,smallwalkvects(:,:,:,lowspf),ii,DATAZERO,&
                rhomat(:,:,:,lowspf),www%nspf**2)
        else if (numspf.gt.0) then
           rhomat(:,:,:,lowspf:highspf)=0d0
        endif
        call myclock(jtime);    times(4)=times(4)+jtime-itime;     itime=jtime
        if (parorbsplit.eq.1.and.sparseconfigflag.eq.0) then
           call mpiorbgather(rhomat,www%nspf**3)
        endif
!! BAD REDUCE.  REDUCE BEFORE ORBGATHER (need new communicators)
        if (sparseconfigflag.ne.0) then
           call mympireduce(rhomat,www%nspf**4)  
        endif
        call myclock(jtime);     times(5)=times(5)+jtime-itime;     itime=jtime
        
        if (conway.eq.2.or.conway.eq.3) then
           rhomat(:,:,:,:)=rhomat(:,:,:,:)/timefac
        endif
        call myclock(jtime);     times(6)=times(6)+jtime-itime;     itime=jtime
        
        call assigncomplexmat(realrhomat,rhomat, www%nspf**2,www%nspf**2)
     
        call myclock(jtime);     times(7)=times(7)+jtime-itime;     itime=jtime

        call DGEMM('N','N',zzz*isize,zzz*www%nspf**2,zzz*www%nspf**2,1d0,projector,&
             zzz*isize,realrhomat,zzz*www%nspf**2,0d0,temppairs,zzz*isize)
        call DGEMM('N','T',zzz*isize,zzz*isize,zzz*www%nspf**2,1d0,temppairs,zzz*isize,&
             projector,zzz*isize,0d0,rhomatpairs,zzz*isize)

        call DGEMM('N','N',2*zzz*isize,zzz*www%nspf**2,zzz*www%nspf**2,1d0,bigprojector,&
             2*zzz*isize,realrhomat,zzz*www%nspf**2,0d0,temppairsbig,2*zzz*isize)
        call DGEMM('N','T',2*zzz*isize,zzz*isize,zzz*www%nspf**2,1d0,temppairsbig,&
             2*zzz*isize,projector,zzz*isize,0d0,rhomatpairsbig,2*zzz*isize)

        call myclock(jtime);     times(8)=times(8)+jtime-itime;     itime=jtime

        if (iiyy.eq.1) then
           call sparseconfigmult(www,avector,avectorp,cptr,sptr,1,1,0,0,time,imc)
        else
           call sparseconfigpulsemult(www,avector,avectorp,cptr,sptr,iiyy-1,imc)
        endif

        call myclock(jtime);    times(9)=times(9)+jtime-itime;     itime=jtime
        
        if (conway.ne.2.and.conway.ne.3) then
           avectorp(:,:)=avectorp(:,:)*timefac   
        endif

        rhs(:,:)=0d0
        if (www%configend.ge.www%configstart) then
           do j=lowspf,highspf
              do i=1,www%nspf
                 rhs(i,j)=dot(smallwalkvects(:,:,i,j),&
                      avectorp(:,www%configstart:www%configend),&
                      (www%configend-www%configstart+1)*numr)
              enddo
           enddo
        endif
        call myclock(jtime);    times(10)=times(10)+jtime-itime;   itime=jtime
        if (parorbsplit.eq.1.and.sparseconfigflag.eq.0) then
           call mpiorbgather(rhs,www%nspf)
        endif
        if (sparseconfigflag.ne.0) then
           call mympireduce(rhs,www%nspf**2)
        endif

        if (www%holeflag.ne.0) then
           rhs(:,:) = (-1) * rhs(:,:)
        endif

        call myclock(jtime);    times(11)=times(11)+jtime-itime;   itime=jtime

        call assigncomplexvec(realrhs(:,:,:),rhs(:,:), www%nspf**2)
        
        rhspairs(:,:)=RESHAPE(MATMUL(RESHAPE(projector,(/zzz*isize,zzz*www%nspf**2/)), &
             RESHAPE(realrhs(:,:,:),(/zzz*www%nspf**2,1/))),(/zzz,isize/))
        rhspairsbig(:,:)=RESHAPE(MATMUL(RESHAPE(bigprojector,&
             (/2*zzz*isize,zzz*www%nspf**2/)), &
             RESHAPE(realrhs(:,:,:),(/zzz*www%nspf**2,1/))),(/zzz,2*isize/))
        
        if (conway.ne.1.and.conway.ne.3) then
           rhomatpairscopy(:,:,:,:)=rhomatpairscopy(:,:,:,:)+rhomatpairs(:,:,:,:)
           rhspairstemp(:,:)=rhspairstemp(:,:)+rhspairs(:,:)
        else
           rhomatpairsbigcopy(:,:,:,:)=rhomatpairsbigcopy(:,:,:,:) + &
                rhomatpairsbig(:,:,:,:)
           rhspairsbigtemp(:,:)=rhspairsbigtemp(:,:)+rhspairsbig(:,:)
        endif

        call myclock(jtime);        times(12)=times(12)+jtime-itime;     
        
     enddo    !!! IMC

     call myclock(itime)

     if (conway.ne.1.and.conway.ne.3) then

!! rhomatpairscopy(zzz,isize,zzz,isize)

        if (conway.eq.2) then   !! make a symmetric matrix equation

#ifdef REALGO
           OFLWR "CONWAY=2 DOESN't MAKE SENSE FOR RELAXATION"; CFLST
#else
           rhomatpairscopy2(1,:,:,:)=rhomatpairscopy(2,:,:,:)*(-1)
           rhomatpairscopy2(2,:,:,:)=rhomatpairscopy(1,:,:,:)
           rhomatpairscopy(:,:,:,:)=rhomatpairscopy2(:,:,:,:)
           rhspairstemp2(1,:)=rhspairstemp(2,:)*(-1)
           rhspairstemp2(2,:)=rhspairstemp(1,:)
           rhspairstemp(:,:)=rhspairstemp2(:,:)
#endif

        endif

        call checksym(rhomatpairscopy,isize*zzz)

        call realinvmatsmooth(rhomatpairscopy,isize*zzz,lioreg,.false.)
        call dgemv('N',isize*zzz,isize*zzz,1d0,rhomatpairscopy,&
             isize*zzz,rhspairstemp,1,0d0,rhspairs,1)

     else
        
        mattemp(:,:)= &
             MATMUL(TRANSPOSE(RESHAPE(rhomatpairsbigcopy(:,:,:,:),&
             (/2*isize*zzz,isize*zzz/))),&
             RESHAPE(rhomatpairsbigcopy(:,:,:,:),(/2*isize*zzz,isize*zzz/)))
        call realinvmatsmooth(mattemp,isize*zzz,lioreg,.false.)
        pseudoinv(:,:)=MATMUL(mattemp,TRANSPOSE(RESHAPE(rhomatpairsbigcopy(:,:,:,:),&
             (/2*isize*zzz,isize*zzz/))))
        call DGEMV('N',zzz*isize,2*zzz*isize,1d0,pseudoinv(:,:),isize*zzz,&
             rhspairsbigtemp,1,0d0,rhspairs(:,:),1)

     endif

     call myclock(jtime);   times(13)=times(13)+jtime-itime;    itime=jtime

     tempmatel(:,:)=0d0

     do ii=1,isize
        i=ipairs(1,ii);     j=ipairs(2,ii)
#ifdef REALGO
        tempmatel(j,i) = rhspairs(1,ii)    
#else
        tempmatel(j,i) = rhspairs(1,ii) + (0d0,1d0) * rhspairs(2,ii)
#endif
     enddo

    tempmatel(:,:)=  ( tempmatel(:,:) -TRANSPOSE(CONJUGATE(tempmatel(:,:)))  ) &
         / timefac  * condamp

     select case(iiyy)
     case(1)
        cptr%xconmatel(:,:)=  tempmatel(:,:)
     case(2)
        cptr%xconmatelxx(:,:)=  tempmatel(:,:)
     case(3)
        cptr%xconmatelyy(:,:)=  tempmatel(:,:)
     case(4)
        cptr%xconmatelzz(:,:)=  tempmatel(:,:)
     end select

     call myclock(jtime);        times(14)=times(14)+jtime-itime;  

  enddo

  deallocate( avectorp,   avector,    &
       smallwalkvects,         rhs,          rhomat,  &
        projector,  bigprojector, mattemp,pseudoinv, &
          realrhs,  realrhomat,  &
       rhomatpairs,  rhspairs,      rhspairstemp,  rhomatpairscopy, &
       rhomatpairscopy2,rhspairstemp2, &
       rhomatpairsbig,  rhspairsbig, rhspairsbigtemp, &
       rhomatpairsbigcopy, &
       sing,  work,  rwork,temppairs,temppairsbig)

end subroutine get_dfconstraint0


subroutine get_dfconstraint(time)
  use parameters
  use xxxmod
  use configmod
  implicit none
  real*8, intent(in) :: time
  DATATYPE :: nullvectors(numr,mcscfnum)
  nullvectors=0
  if (tot_adim.gt.0) then
     call get_dfconstraint0(yyy%cmfavec(:,:,0),mcscfnum,yyy%cptr(0),yyysptr(0),www,time)
  else
     call get_dfconstraint0(nullvectors(:,:),mcscfnum,yyy%cptr(0),yyysptr(0),www,time)
  endif
end subroutine get_dfconstraint


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!  DENSITY MATRIX (constraintflag=1) !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine get_denconstraint(time)
  use ham_parameters
  implicit none
  real*8,intent(in) :: time

!! assume nothing, keep constant off block diag (lio solve)

  if (denmatfciflag.ne.0) then
!! "denmat FCI", wrong equation, as per restricted configuration paper
     call get_denconstraint1(2,time)    
  else
!! correct formula
     call new_get_denconstraint1(time)    
  endif

end subroutine get_denconstraint


!! lioden.    111510 WAS C-order.  1)  change refs to denmat because 
!!                                     denmat is now denmat not transpose denmat
!!                                 2)  change C-order because think that also has 
!!        to do with denmat transpose and not just internal to these subroutines
!!                                   
!! 111510 WAS function lind(ispf,jspf)

function lind(jspf,ispf)
  use parameters
  implicit none
  integer,intent(in) :: ispf,jspf
  integer :: lind
  print *, "CHECK LIODEN CODE. (SEE COMMENTS)"
  stop
  if (jspf.gt.ispf) then
     lind=(ispf-1)*(nspf-1)+jspf-1
  else
     lind=(ispf-1)*(nspf-1)+jspf
  endif
end function lind


function llind(ispf,jspf)
  use parameters
  implicit none
  integer,intent(in) :: ispf,jspf
  integer :: llind, lll, i

  lll=0
  do i=1,shells(ispf)-1
     lll=lll+ (nspf - allshelltop(i)+allshelltop(i-1)) &
          * (allshelltop(i)-allshelltop(i-1))
  enddo

  lll=lll + (nspf-allshelltop(shells(ispf))+allshelltop(shells(ispf)-1)) &
       * (ispf-allshelltop(shells(ispf)-1)-1) + jspf

  if (shells(jspf).gt.shells(ispf)) then
     lll=lll-(allshelltop(shells(ispf))-allshelltop(shells(ispf)-1))
  endif
  llind=lll

end function llind



subroutine get_denconstraint1_0(www,cptr,sptr,numvects,avector,drivingavectorsxx, &
     drivingavectorsyy,drivingavectorszz,denmat,iwhich,time)
  use fileptrmod
  use r_parameters
  use basis_parameters
  use ham_parameters
  use constraint_parameters
  use walkmod
  use spfsize_parameters   !! parorbsplit
  use configptrmod
  use sparseptrmod
  use dotmod
  use sparsemultmod
  use invsubmod
  use orbgathersubmod
  use mpisubmod
  use utilmod
  implicit none
  integer,intent(in) :: numvects
  type(walktype),intent(in) :: www
  type(CONFIGPTR) :: cptr
  type(SPARSEPTR) :: sptr
  DATATYPE,intent(in) :: denmat(www%nspf,www%nspf),&
       avector(numr,www%firstconfig:www%lastconfig,numvects),&
       drivingavectorsxx(numr,www%firstconfig:www%lastconfig,numvects),&
       drivingavectorsyy(numr,www%firstconfig:www%lastconfig,numvects),&
       drivingavectorszz(numr,www%firstconfig:www%lastconfig,numvects)
  DATATYPE ::  a1(numr,numvects), a2(numr,numvects), a1p(numr,numvects), &
       a2p(numr,numvects), tempconmatels(www%nspf,www%nspf), csum
  DATATYPE,allocatable :: bigavector(:,:,:),bigavectorp(:,:,:), avectorp(:,:,:)
  integer ::  config1,config2,   ispf,jspf,  dirphase,  i,   &
       iwalk, info, kspf, lspf, ind, jind,  lowspf,highspf, &
       flag, isize, iwhich,  iiyy,maxii,imc,ihop, liosize
  integer :: ipiv(www%nspf*(www%nspf-1))
  real*8 :: denom,time,rsum,rsum2,maxval,maxanti
  DATATYPE :: liosolve(www%nspf*(www%nspf-1)),lioden(www%nspf*(www%nspf-1), www%nspf*(www%nspf-1)),&
       liodencopy(www%nspf*(www%nspf-1),www%nspf*(www%nspf-1)),liosolvetemp(www%nspf*(www%nspf-1)),&
       myliosolve(www%nspf*(www%nspf-1)),nullvector1(numr),nullvector2(numr)

  nullvector1(:)=0; nullvector2(:)=0
  cptr%xconmatel(:,:)=0.d0;   cptr%xconmatelxx(:,:)=0.d0;   
  cptr%xconmatelyy(:,:)=0.d0;   cptr%xconmatelzz(:,:)=0.d0

  liosize=www%nspf*(www%nspf-1)

  if ((iwhich.eq.2).and.(numshells.eq.1)) then
     return
  endif
  if (real(timefac) /= 0.0) then
     OFLWR "Err, get_denconstraint1 (with lioden) can only be used for forward time propagation. "
     CFLST
  endif

  allocate(bigavector(numr,www%numconfig,numvects), &
       bigavectorp(numr,www%numconfig,numvects), &
       avectorp(numr,www%firstconfig:www%lastconfig,numvects))
  bigavector(:,:,:)=0d0; bigavectorp(:,:,:)=0d0; 
  if (www%lastconfig.ge.www%firstconfig) then
     avectorp(:,:,:)=0d0
  endif

  lioden=0.d0

  lowspf=1; highspf=www%nspf
  if (parorbsplit.eq.1) then
     call checkorbsetrange(www%nspf,flag)
     if (flag.ne.0) then
        OFLWR "error exit, can't do get_denconstraint1_0 with ",www%nspf,&
             " orbitals sparse mode"; CFLST
     endif
     call getOrbSetRange(lowspf,highspf)
  endif

  do ispf=lowspf,highspf
     do kspf=1,www%nspf
        do lspf=1,www%nspf
           flag=0
           select case (iwhich)
           case (1)
              if ((ispf/=lspf).and.(kspf/=lspf)) then
                 ind=lind(ispf,lspf);      jind=lind(kspf,lspf);    flag=1
              endif
           case (2)
              if ( (shells(ispf).ne.shells(lspf) ) .and. &
                   (shells(kspf).ne.shells(lspf) ) ) then
                 ind=llind(ispf,lspf);    jind=llind(kspf,lspf);       flag=1
              endif
           case default
              ind=0;     jind=0
              call openfile(); write(mpifileptr,*) 
              OFLWR "get_denconstraint1 error"; CFLST
           end select

           if (flag==1) then
              lioden(ind, jind) = lioden(ind, jind) + &
                   denmat(ispf,kspf) * CONJUGATE(timefac)
           endif
        enddo
     enddo
  enddo

  do jspf=lowspf,highspf
     do kspf=1,www%nspf
        do lspf=1,www%nspf
           flag=0
           select case (iwhich)
           case (1)
              if ((kspf/=jspf).and.(kspf/=lspf)) then
                 ind=lind(kspf,jspf);    jind=lind(kspf,lspf);    flag=1
              endif
           case (2)
              if ( (shells(kspf).ne.shells(jspf) ) .and. &
                   (shells(kspf).ne.shells(lspf) ) ) then
                 ind=llind(kspf,jspf);   jind=llind(kspf,lspf);     flag=1
              endif
           case default
              ind=0
              jind=0
              OFLWR "get_denconstraint1 error"; CFLST
           end select

           if (flag==1) then
              lioden(ind,jind) = lioden(ind,jind) + &
                   denmat(lspf,jspf) * timefac
           endif
        enddo
     enddo
  enddo

  if (parorbsplit.eq.1) then
     call mpiorbreduce(lioden,liosize**2)
  endif

  if (www%holeflag.ne.0) then
     lioden(:,:)=(-1)*lioden(:,:)
  endif

  maxii=1
  if (tdflag.ne.0) then
     maxii=4
  endif

  do iiyy=1,maxii
     do imc=1,numvects
        if (www%totadim.gt.0) then
           select case(iiyy)
           case(1)
              call sparseconfigmult(www,avector(:,:,imc),avectorp(:,:,imc),cptr,sptr,&
                   1,1,0,0,time,imc)
           case default
              call sparseconfigpulsemult(www,avector(:,:,imc),avectorp(:,:,imc),cptr,&
                   sptr,iiyy-1,imc)
              if (drivingflag.ne.0) then
                 if (iiyy.eq.2) then
                    avectorp(:,:,imc)=avectorp(:,:,imc)+drivingavectorsxx(:,:,imc)
                 else if (iiyy.eq.3) then
                    avectorp(:,:,imc)=avectorp(:,:,imc)+drivingavectorsyy(:,:,imc)
                 else if (iiyy.eq.4) then
                    avectorp(:,:,imc)=avectorp(:,:,imc)+drivingavectorszz(:,:,imc)
                 endif
              endif
           end select
        else
           select case(iiyy)
           case(1)
              call sparseconfigmult(www,nullvector1(:),nullvector2(:),cptr,sptr,&
                   1,1,0,0,time,imc)
           case default
              call sparseconfigpulsemult(www,nullvector1(:),nullvector2(:),cptr,&
                   sptr,iiyy-1,imc)
           end select
        endif
     enddo

!! DO SUMMA (parconsplit.ne.0 and sparsesummaflag.eq.2, "circ")

     if (www%lastconfig.ge.www%firstconfig) then
        bigavector(:,www%firstconfig:www%lastconfig,:)=avector(:,:,:)
        bigavectorp(:,www%firstconfig:www%lastconfig,:)=avectorp(:,:,:)
     endif

     if (www%parconsplit.ne.0) then
        do i=1,numvects
           call mpiallgather(bigavector(:,:,i),www%numconfig*numr,&
                www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
           call mpiallgather(bigavectorp(:,:,i),www%numconfig*numr,&
                www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
        enddo
     endif

     liosolve(:)=0.d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(config1,a1,a1p,ihop,iwalk,csum,config2,dirphase,a2,a2p,ispf,jspf,flag,ind, myliosolve)

     myliosolve(:)=0d0

!$OMP DO SCHEDULE(DYNAMIC)
     do config1=www%botconfig,www%topconfig

        a1(:,:)=avector(:,config1,:)
        a1p(:,:)=avectorp(:,config1,:)

        do ihop=1,www%numsinglehops(config1)
           config2=www%singlehop(ihop,config1)

           a2(:,:)=bigavector(:,config2,:)
           a2p(:,:)=bigavectorp(:,config2,:)

           csum = dot(a2p(:,:)*timefac,a1(:,:),numr*numvects) + &
                dot(a2(:,:),a1p(:,:)*timefac,numr*numvects)

           do iwalk=www%singlehopwalkstart(ihop,config1),www%singlehopwalkend(ihop,config1)

              dirphase=www%singlewalkdirphase(iwalk+www%scol(config1))

              ispf=www%singlewalkopspf(1,iwalk+www%scol(config1))
              jspf=www%singlewalkopspf(2,iwalk+www%scol(config1))
              
              flag=0
              select case (iwhich)
              case (1)
                 ind=lind(ispf,jspf) ; flag=1  !!problem?
              case (2)
                 if (shells(ispf).ne.shells(jspf)) then
!! FIRSTWAY
                    ind=llind(ispf,jspf);                 flag=1
                 endif
              case default
                 ind=0
                 OFLWR "get_denconstraint1 error"; CFLST
              end select

              if (flag==1) then
                 myliosolve(ind) = myliosolve(ind) +     dirphase*csum
              endif
           enddo  !! iwalk
        enddo     !! ihop
     enddo
!$OMP END DO
!$OMP CRITICAL
     liosolve(:)=liosolve(:)+myliosolve(:)
!$OMP END CRITICAL
!$OMP END PARALLEL

     call mympireduce(liosolve,liosize)

     select case (iwhich)
     case (1)
        isize=liosize
     case (2)
        isize=www%nspf**2
        do i=1,numshells
           isize=isize-(allshelltop(i)-allshelltop(i-1))**2
        enddo
     case default
        OFLWR "get_denconstraint1 error"; CFLST
     end select

     liodencopy(:,:)=lioden(:,:)
     if (lioreg.le.0.d0) then
        call MYGESV(isize, 1, liodencopy, liosize, ipiv, liosolve, liosize, info)
        if (info/=0) then
           OFLWR "Errx  mygesv lioden", info, " size ", liosize, isize ; CFLST
        endif
     else
        call invmatsmooth(liodencopy,liosize,liosize,lioreg,.false.)
        call MYGEMV('N',liosize,liosize,DATAONE,liodencopy,liosize,liosolve,1,&
             DATAZERO,liosolvetemp,1)
        liosolve(:)=liosolvetemp(:)
     endif

     tempconmatels(:,:)=0d0
     do ispf=1,www%nspf
        do jspf=1,www%nspf
           select case (iwhich)
           case (1)
              if (ispf/=jspf) then
                 ind=lind(ispf,jspf)   
                 tempconmatels(ispf,jspf)=liosolve(ind)
              endif
           case (2)
              if (shells(ispf)/=shells(jspf)) then
                 ind=llind(ispf,jspf)
                 tempconmatels(ispf,jspf)=liosolve(ind)
              endif
           case default
              OFLWR "get_denconstraint1 error"; CFLST
           end select
        enddo
     enddo

!! 111510   REGARDLESS!  #ifndef ECSFLAG
!! may require you to define constraint via hermitian part of hamiltonian for chmctdh.  
!! Don't want conmatels to be non-antiherm (chmctdh) or non real (cmctdh)

     maxval=0d0
     maxanti=0d0

     do ispf=1,www%nspf
        do jspf=ispf+1,www%nspf
           
           rsum=abs(timefac*tempconmatels(ispf,jspf)+CONJUGATE(timefac*tempconmatels(jspf,ispf)))

           rsum2=max(abs(tempconmatels(ispf,jspf)),abs(tempconmatels(jspf,ispf)))

           denom=  max(1d-5,rsum2)

           if (rsum / denom .gt.1.d-6) then
              OFLWR "Err herm incmatel temp continue"
              WRFL ispf, jspf, rsum, denom,iiyy; WRFL; CFL !!ST
           endif
           if (rsum.gt.maxanti) then
              maxanti=rsum
           endif
           if (rsum2.gt.maxval) then
              maxval=rsum2
           endif
        enddo
     enddo

!! 070414     
     tempconmatels(:,:)=0.5d0*(tempconmatels(:,:)+TRANSPOSE(CONJUGATE(tempconmatels(:,:))))

     select case(iiyy)
     case(1)
        cptr%xconmatel(:,:)=tempconmatels(:,:)
     case(2)
        cptr%xconmatelxx(:,:)=tempconmatels(:,:)
     case(3)
        cptr%xconmatelyy(:,:)=tempconmatels(:,:)
     case(4)
        cptr%xconmatelzz(:,:)=tempconmatels(:,:)
     end select
     
  end do

  deallocate(bigavector,bigavectorp,avectorp)

end subroutine get_denconstraint1_0


subroutine get_denconstraint1(iwhich,time)
  use parameters
  use xxxmod
  use configmod
  implicit none
  integer,intent(in) :: iwhich
  real*8,intent(in) :: time
  DATATYPE :: nullvectors(numr,mcscfnum)

  nullvectors(:,:)=0

  if (www%totadim.gt.0) then
     call get_denconstraint1_0(www,yyy%cptr(0),yyysptr(0),mcscfnum,&
          yyy%cmfavec(:,:,0),yyy%drivingavectorsxx(:,:,:,0),&
          yyy%drivingavectorsyy(:,:,:,0),yyy%drivingavectorszz(:,:,:,0),&
          yyy%denmat(:,:,0),iwhich,time)
  else
     call get_denconstraint1_0(www,yyy%cptr(0),yyysptr(0),mcscfnum,&
          nullvectors(:,:),nullvectors(:,:),&
          nullvectors(:,:),nullvectors(:,:),&
          yyy%denmat(:,:,0),iwhich,time)
  endif

end subroutine get_denconstraint1




subroutine new_get_denconstraint1_0(www,cptr,sptr,numvects,avector,drivingavectorsxx, &
     drivingavectorsyy,drivingavectorszz,denmat,time)
  use fileptrmod
  use constraint_parameters
  use ham_parameters
  use basis_parameters
  use sparse_parameters
  use spfsize_parameters !! parorbsplit
  use r_parameters
  use walkmod
  use configptrmod
  use sparseptrmod
  use dotmod
  use sparsemultmod
  use basissubmod
  use invsubmod
  use orbgathersubmod
  use mpisubmod
  use utilmod
  implicit none
  integer,intent(in) :: numvects
  type(walktype),intent(in) :: www
  type(CONFIGPTR) :: cptr
  type(SPARSEPTR) :: sptr
  DATATYPE,intent(in) :: denmat(www%nspf,www%nspf),&
       avector(numr,www%firstconfig:www%lastconfig,numvects),&
       drivingavectorsxx(numr,www%firstconfig:www%lastconfig,numvects),&
       drivingavectorsyy(numr,www%firstconfig:www%lastconfig,numvects),&
       drivingavectorszz(numr,www%firstconfig:www%lastconfig,numvects)
  integer :: ipairs(2,www%nspf*(www%nspf-1))
  DATATYPE ::  a1(numr,numvects), a2(numr,numvects), a1p(numr,numvects), &
       a2p(numr,numvects), csum, nullvector1(numr),nullvector2(numr)
  DATATYPE,allocatable :: bigavector(:,:,:), bigavectorp(:,:,:),avectorp(:,:,:),&
       rhomat(:,:,:,:),tempconmatels(:,:) 
  integer ::  config1,config2,   ispf,jspf,  dirphase,  i,  iwalk,  info, kspf, &
       lspf, ind, jind, flag, isize,   &
       iiyy,maxii,imc,j,ihop,lowspf,highspf,lowsize,highsize
  integer,allocatable :: ipiv(:)
  real*8 :: denom,time,rsum,rsum2,maxval,maxanti
  DATATYPE, allocatable :: liosolve(:),lioden(:, :),liodencopy(:,:),&
       liosolvetemp(:),myliosolve(:)

  nullvector1(:)=0;nullvector2(:)=0
  cptr%xconmatel(:,:)=0.d0;   cptr%xconmatelxx(:,:)=0.d0;   
  cptr%xconmatelyy(:,:)=0.d0;   cptr%xconmatelzz(:,:)=0.d0

  if ((numshells.eq.1)) then
     return
  endif

  allocate(bigavector(numr,www%numconfig,numvects), &
       bigavectorp(numr,www%numconfig,numvects),&
       avectorp(numr,www%firstconfig:www%lastconfig,numvects),&
       rhomat(www%nspf,www%nspf,www%nspf,www%nspf),&
       tempconmatels(www%nspf,www%nspf))
  bigavector(:,:,:)=0d0; bigavectorp(:,:,:)=0d0;   rhomat(:,:,:,:)=0d0;
  tempconmatels=0
  if (www%lastconfig.ge.www%firstconfig) then
     avectorp(:,:,:)=0d0
  endif

  if (www%dfrestrictflag.gt.www%dflevel) then
     call get_rhomat(www,avector,rhomat,numr,numvects)
  endif

  lowspf=1; highspf=www%nspf
  if (parorbsplit.eq.1) then
     call checkorbsetrange(www%nspf,flag)
     if (flag.ne.0) then
        OFLWR "error exit, can't do get_denconstraint1_0 with ",www%nspf,&
             " orbitals sparse mode"; CFLST
     endif
     call getOrbSetRange(lowspf,highspf)
  endif

  isize=0
  do i=1,www%nspf
     do j=1,www%nspf
        if (shells(i).ne.shells(j)) then
           isize=isize+1
           ipairs(:,isize)=(/ i,j /)
        endif
     enddo
  enddo

  highsize=0
  lowsize=isize+1

  isize=0
  do i=1,www%nspf
     if (i.eq.lowspf) then
        lowsize=isize+1
     endif
     do j=1,www%nspf
        if (shells(i).ne.shells(j)) then
           isize=isize+1
           ipairs(:,isize)=(/ i,j /)
        endif
     enddo
     if (i.eq.highspf) then
        highsize=isize
     endif
  enddo

  allocate(liosolve(isize),lioden(isize, isize),liodencopy(isize,isize),&
       liosolvetemp(isize),ipiv(isize))  
  liosolve=0; liodencopy=0; liosolvetemp=0;  lioden=0.d0

  do ind=lowsize,highsize
     ispf=ipairs(1,ind)
     jspf=ipairs(2,ind)
     do jind=1,isize
        kspf=ipairs(1,jind)
        lspf=ipairs(2,jind)

        rsum= abs(rhomat(ispf,jspf,kspf,lspf)-CONJUGATE(rhomat(kspf,lspf,ispf,jspf)))
        if (rsum.gt.1d-10) then
           OFLWR "DOOG",rsum,rhomat(ispf,jspf,kspf,lspf),&
                CONJUGATE(rhomat(kspf,lspf,ispf,jspf)); CFLST
        endif

!! rhomat included,excluded

        lioden(ind,jind)=  rhomat(kspf,lspf,ispf,jspf) - rhomat(jspf,ispf,lspf,kspf)   

        if (jspf.eq.lspf) then
           lioden(ind, jind) = lioden(ind, jind) - &
                denmat(ispf,kspf)  
        endif
        if (ispf.eq.kspf) then
           lioden(ind,jind) = lioden(ind,jind) + &
                denmat(lspf,jspf)  
        endif
     enddo
  enddo

  if (parorbsplit.eq.1) then
     call mpiorbreduce(lioden,isize**2)
  endif

  if (www%holeflag.ne.0) then
     lioden(:,:)=(-1)*lioden(:,:)
  endif

  maxii=1
  if (tdflag.ne.0) then
     maxii=4
  endif

  do iiyy=1,maxii
     do imc=1,numvects
        if (www%totadim.gt.0) then
           select case(iiyy)
           case(1)
              call sparseconfigmult(www,avector(:,:,imc),avectorp(:,:,imc),cptr,sptr,&
                   1,1,0,0,time,imc)
           case default
              call sparseconfigpulsemult(www,avector(:,:,imc),avectorp(:,:,imc),cptr,&
                   sptr,iiyy-1,imc)
              if (drivingflag.ne.0) then
                 if (iiyy.eq.2) then
                    avectorp(:,:,imc)=avectorp(:,:,imc)+drivingavectorsxx(:,:,imc)
                 else if (iiyy.eq.3) then
                    avectorp(:,:,imc)=avectorp(:,:,imc)+drivingavectorsyy(:,:,imc)
                 else if (iiyy.eq.4) then
                    avectorp(:,:,imc)=avectorp(:,:,imc)+drivingavectorszz(:,:,imc)
                 endif
              endif
           end select
        else
           select case(iiyy)
           case(1)
              call sparseconfigmult(www,nullvector1(:),nullvector2(:),cptr,sptr,&
                   1,1,0,0,time,imc)
           case default
              call sparseconfigpulsemult(www,nullvector1(:),nullvector2(:),cptr,&
                   sptr,iiyy-1,imc)
           end select
        endif
     enddo

     if (www%dfrestrictflag.gt.www%dflevel) then
        do imc=1,numvects
           if (www%totadim.gt.0) then
              call df_project(www,numr,avectorp(:,:,imc))
           else
              call df_project(www,numr,nullvector1(:))
           endif
        enddo
     endif

!! DO SUMMA (parconsplit.ne.0 and sparsesummaflag.eq.2, "circ")

     if (www%lastconfig.ge.www%firstconfig) then
        bigavector(:,www%firstconfig:www%lastconfig,:)=avector(:,:,:)
        bigavectorp(:,www%firstconfig:www%lastconfig,:)=avectorp(:,:,:)
     endif

     if (www%parconsplit.ne.0) then
        do i=1,numvects
           call mpiallgather(bigavector(:,:,i),www%numconfig*numr,&
                www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
           call mpiallgather(bigavectorp(:,:,i),www%numconfig*numr,&
                www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
        enddo
     endif

     liosolve(:)=0d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(config1,a1,a1p,ihop,iwalk,csum,config2,dirphase,a2,a2p,ispf,jspf,flag,ind,myliosolve)

     allocate(myliosolve(isize));     myliosolve(:)=0d0

!$OMP DO SCHEDULE(DYNAMIC)
     do config1=www%botconfig,www%topconfig

        a1(:,:)=avector(:,config1,:)
        a1p(:,:)=avectorp(:,config1,:)
           
        do ihop=1,www%numsinglehops(config1)
           config2=www%singlehop(ihop,config1)

           a2(:,:)=bigavector(:,config2,:)
           a2p(:,:)=bigavectorp(:,config2,:)

           csum = dot(a2p(:,:)*timefac,a1(:,:),numvects*numr) + &
                dot(a2(:,:),a1p(:,:)*timefac,numvects*numr)

           do iwalk=www%singlehopwalkstart(ihop,config1),www%singlehopwalkend(ihop,config1)
              
              dirphase=www%singlewalkdirphase(iwalk+www%scol(config1))
              ispf=www%singlewalkopspf(1,iwalk+www%scol(config1))
              jspf=www%singlewalkopspf(2,iwalk+www%scol(config1))

              flag=0
              if (shells(ispf).ne.shells(jspf)) then
!! FIRSTWAY
                 ind=llind(ispf,jspf);                 flag=1
              endif
              
              if (flag==1) then
                 myliosolve(ind)=myliosolve(ind) + dirphase*csum
              endif
           enddo  !! iwalk
        enddo     !! ihop
     enddo
!$OMP END DO
!$OMP CRITICAL
     liosolve(:)=liosolve(:) + myliosolve(:)
!$OMP END CRITICAL
     deallocate(myliosolve)
!$OMP END PARALLEL

     call mympireduce(liosolve,isize)

     liodencopy(:,:)=lioden(:,:)       

     if (lioreg.le.0.d0) then
        call MYGESV(isize, 1, liodencopy, isize, ipiv, liosolve, isize, info)
        if (info/=0) then
           OFLWR "Errx  mygesv lioden", info, " size ", isize ; CFLST
        endif
     else
        call invmatsmooth(liodencopy,isize,isize,lioreg,.false.)
        call MYGEMV('N',isize,isize,DATAONE,liodencopy,isize,&
             liosolve,1,DATAZERO,liosolvetemp,1)
        liosolve(:)=liosolvetemp(:)
     endif

     tempconmatels(:,:)=0d0
     do ispf=1,www%nspf
        do jspf=1,www%nspf
           if (shells(ispf)/=shells(jspf)) then
              ind=llind(ispf,jspf)
              tempconmatels(ispf,jspf)=liosolve(ind)/timefac 
           endif
        enddo
     enddo
     
!! 111510   REGARDLESS!  #ifndef ECSFLAG
!! may require you to define constraint via hermitian part of hamiltonian for chmctdh.  
!! Don't want conmatels to be non-antiherm (chmctdh) or non real (cmctdh)

     maxval=0d0
     maxanti=0d0

     do ispf=1,www%nspf
        do jspf=ispf+1,www%nspf
           
           rsum=abs(tempconmatels(ispf,jspf)-CONJUGATE(tempconmatels(jspf,ispf)))

           rsum2=max(abs(tempconmatels(ispf,jspf)),abs(tempconmatels(jspf,ispf)))

           denom=  max(1d-5,rsum2)

           if (rsum / denom .gt.1.d-6) then
              OFLWR "Err herm incmatel temp continue"
              WRFL ispf, jspf, rsum, denom,iiyy; WRFL; CFL !!ST
           endif

           if (rsum.gt.maxanti) then
              maxanti=rsum
           endif
           if (rsum2.gt.maxval) then
              maxval=rsum2
           endif
        enddo
     enddo

!! 070414     
     tempconmatels(:,:)=0.5d0*(tempconmatels(:,:)+TRANSPOSE(CONJUGATE(tempconmatels(:,:))))

     select case(iiyy)
     case(1)
        cptr%xconmatel(:,:)=tempconmatels(:,:)
     case(2)
        cptr%xconmatelxx(:,:)=tempconmatels(:,:)
     case(3)
        cptr%xconmatelyy(:,:)=tempconmatels(:,:)
     case(4)
        cptr%xconmatelzz(:,:)=tempconmatels(:,:)
     end select
     
  end do

  deallocate(bigavector,bigavectorp,avectorp,rhomat,tempconmatels)

end subroutine new_get_denconstraint1_0


subroutine new_get_denconstraint1(time)
  use parameters
  use xxxmod
  use configmod
  implicit none
  real*8,intent(in) :: time
  DATATYPE :: nullvectors(numr,mcscfnum)

  nullvectors(:,:)=0

  if (www%totadim.gt.0) then
  call new_get_denconstraint1_0(www,yyy%cptr(0),yyysptr(0),mcscfnum,&
       yyy%cmfavec(:,:,0),yyy%drivingavectorsxx(:,:,:,0),&
       yyy%drivingavectorsyy(:,:,:,0),yyy%drivingavectorszz(:,:,:,0),&
       yyy%denmat(:,:,0),time)
  else
     call new_get_denconstraint1_0(www,yyy%cptr(0),yyysptr(0),mcscfnum,&
          nullvectors(:,:),nullvectors(:,:),&
          nullvectors(:,:),nullvectors(:,:),&
          yyy%denmat(:,:,0),time)
  endif

end subroutine new_get_denconstraint1

end module dfconsubmod
