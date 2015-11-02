
#include "Definitions.INC"

!! INDICATES WHETHER A GIVEN CONFIG is included in the restricted configuration list.

!! it is only included if all double excitations from it are included in the total configuration list.



function dfallowed(thisconfig)
  use parameters
  implicit none
  integer,intent(in) :: thisconfig(ndof)
  logical :: allowedconfig0, dfallowed

  dfallowed = allowedconfig0(thisconfig,dfrestrictflag)

end function dfallowed

function getdfindex(config1)
  use parameters
  use dfconmod
  implicit none
  integer :: getdfindex, config1,j,flag, dir,newdir, step

  getdfindex=-1

  if (dfincludedmask(config1).eq.0) then
     return
  endif
  
  dir=1;  j=1;  step=max(1,numdfconfigs/4);  flag=0

  do while (flag.eq.0)
     flag=1

     if (dfincludedconfigs(j) .ne. config1) then
        flag=0
        if (dfincludedconfigs(j).lt.config1) then
           newdir=1
        else
           newdir=-1
        endif
        if (newdir.ne.dir) then
           step=max(1,step/2)
        endif
        dir=newdir;           j=j+dir*step
        if (j.le.1) then
           j=1;              dir=1
        endif
        if (j.ge.numdfconfigs) then
           j=numdfconfigs;              dir=-1
        endif
     endif
  enddo

  getdfindex=j
  
end function getdfindex



! "included" configurations are those kept in the calculation.
! there are numdfconfigs INCLUDED configurations.

subroutine init_dfcon()
  use parameters
  use dfconmod
  use configmod
  use walkmod
  use mpimod
  implicit none

  integer :: i,j,iconfig,nondfconfigs,dfrank,ii
  logical :: dfallowed
  integer,allocatable :: dfnotconfigs(:)

  allocate(dfincludedmask(numconfig), dfincludedconfigs(numconfig), &
       dfnotconfigs(numconfig),  dfconfsperproc(nprocs), &
       allbotdfconfigs(nprocs),alltopdfconfigs(nprocs))

  if (dfrestrictflag.eq.0) then
     dfincludedmask(:)=1; dfnotconfigs(:)=(-1)
     do i=1,numconfig
        dfincludedconfigs(i)=i
     enddo
     numdfconfigs=numconfig; nondfconfigs=0
     dfconfsperproc(:)=configsperproc(:)
     allbotdfconfigs(:)=allbotconfigs(:)
     alltopdfconfigs(:)=alltopconfigs(:)
     maxdfconfsperproc=maxconfigsperproc
     botdfconfig=botconfig;  topdfconfig=topconfig

  else

     dfincludedmask(:)=0;  dfincludedconfigs(:)=(-1);  dfnotconfigs(:)=(-1)
     numdfconfigs=0;  nondfconfigs=0

     do i=1,numconfig
        if (dfallowed(configlist(:,i))) then 
           numdfconfigs=numdfconfigs+1
           dfincludedmask(i)=1;        dfincludedconfigs(numdfconfigs)=i
        else
           nondfconfigs=nondfconfigs+1;        dfnotconfigs(nondfconfigs)=i
        endif
     enddo

     if (numdfconfigs.eq.0) then
        OFLWR "Error, no included DF configs!"; CFLST
     endif

     OFLWR "Number of included configurations ", numdfconfigs
     WRFL  "Number of not included ", nondfconfigs;  WRFL; CFL

     dfrank=0
     do i=botconfig,topconfig
        if (dfallowed(configlist(:,i))) then 
           dfrank=dfrank+1
        endif
     enddo

     dfconfsperproc(myrank)=dfrank
     maxdfconfsperproc=0
     ii=0
     do i=1,nprocs
        allbotdfconfigs(i)=ii+1
        call mympiibcastone(dfconfsperproc(i),i)
        ii=ii+dfconfsperproc(i)
        alltopdfconfigs(i)=ii
        if (dfconfsperproc(i).gt.maxdfconfsperproc) then
           maxdfconfsperproc=dfconfsperproc(i)
        endif
     enddo
     botdfconfig=allbotdfconfigs(myrank)
     topdfconfig=alltopdfconfigs(myrank)

  endif

  numdfwalks=0

  if (dfrestrictflag.eq.0) then

     allocate(dfwalkfrom(1), dfwalkto(1), includedorb(1), &
          excludedorb(1), dfwalkphase(1))
     dfwalkfrom(:)=(-1); dfwalkto(:)=(-1); includedorb(:)=(-1)
     excludedorb(:)=(-1); dfwalkphase(:)=(-798)

  else

     do i=1,NONdfconfigs
        iconfig=dfNOTconfigs(i)
        if (iconfig.ge.configstart.and.iconfig.le.configend) then
           do j=1,numsinglewalks(iconfig)
              if (dfallowed(configlist(:,singlewalk(j,iconfig)))) then
                 numdfwalks=numdfwalks+1
              endif
           enddo
        endif
     enddo

     OFLWR "numdfwalks:  ", numdfwalks," is total number DF single excitations on this processor"; WRFL; CFL

     allocate(dfwalkfrom(numdfwalks), dfwalkto(numdfwalks), includedorb(numdfwalks), excludedorb(numdfwalks), dfwalkphase(numdfwalks))

     numdfwalks=0
     do i=1,NONdfconfigs
        iconfig=dfNOTconfigs(i)
        if (iconfig.ge.configstart.and.iconfig.le.configend) then
           do j=1,numsinglewalks(iconfig)
              if (dfallowed(configlist(:,singlewalk(j,iconfig)))) then
                 numdfwalks=numdfwalks+1
              
                 dfwalkTO(numdfwalks)=iconfig
                 dfwalkFROM(numdfwalks)=singlewalk(j,iconfig)
                 EXcludedorb(numdfwalks)=singlewalkopspf(1,j,iconfig)
                 INcludedorb(numdfwalks)=singlewalkopspf(2,j,iconfig)
              
                 dfwalkphase(numdfwalks)=singlewalkdirphase(j,iconfig)
              endif
           enddo
        endif
     enddo

     OFLWR "numdfwalks:  ", numdfwalks," is really number DF single excitations on this processor"; WRFL; CFL
     
  endif

  deallocate(dfnotconfigs)

end subroutine init_dfcon

subroutine dfcondealloc()
  use parameters
  use dfconmod
  implicit none
  deallocate(dfincludedmask,dfincludedconfigs)
  deallocate(dfwalkfrom, dfwalkto, includedorb, excludedorb,dfwalkphase)
  numdfwalks=-1;  numdfconfigs=-1
end subroutine dfcondealloc


subroutine df_project(howmany,avector)
  use parameters
  implicit none
  integer :: howmany
  DATATYPE :: avector(howmany,firstconfig:lastconfig)

  if (parconsplit.eq.0) then
     call df_project_all(howmany,avector)
  else
     call df_project_local(howmany,avector)
  endif

end subroutine df_project


subroutine df_project_all(howmany,avector)
  use parameters
  use dfconmod
  implicit none
  integer :: howmany,i
  DATATYPE :: avector(howmany,numconfig)
  do i=1,numconfig
     avector(:,i)=avector(:,i)*dfincludedmask(i)
  enddo
end subroutine df_project_all


subroutine df_project_local(howmany,avector)
  use parameters
  use dfconmod
  implicit none
  integer :: howmany,i
  DATATYPE :: avector(howmany,botconfig:topconfig)
  do i=botconfig,topconfig
     avector(:,i)=avector(:,i)*dfincludedmask(i)
  enddo
end subroutine df_project_local


subroutine df_transformto_all(howmany,avectorin,avectorout)
  use parameters
  use dfconmod
  implicit none
  integer :: howmany,i
  DATATYPE,intent(in) :: avectorin(howmany,numconfig)
  DATATYPE,intent(out) :: avectorout(howmany,numdfconfigs)

  do i=1,numdfconfigs
     avectorout(:,i)=avectorin(:,dfincludedconfigs(i))
  enddo

end subroutine df_transformto_all


subroutine df_transformto_local(howmany,avectorin,avectorout)
  use parameters
  use dfconmod
  implicit none
  integer :: howmany,i
  DATATYPE,intent(in) :: avectorin(howmany,botconfig:topconfig)
  DATATYPE,intent(out) :: avectorout(howmany,botdfconfig:topdfconfig)

  do i=botdfconfig,topdfconfig
     avectorout(:,i)=avectorin(:,dfincludedconfigs(i))
  enddo

end subroutine df_transformto_local



subroutine df_transformfrom_all(howmany,avectorin,avectorout)
  use parameters
  use dfconmod
  implicit none
  integer :: howmany,i
  DATATYPE,intent(in) :: avectorin(howmany,numdfconfigs)
  DATATYPE,intent(out) :: avectorout(howmany,numconfig)

  avectorout(:,:)=0d0

  do i=1,numdfconfigs
     avectorout(:,dfincludedconfigs(i))=avectorin(:,i)
  enddo

end subroutine df_transformfrom_all



subroutine df_transformfrom_local(howmany,avectorin,avectorout)
  use parameters
  use dfconmod
  implicit none
  integer :: howmany,i
  DATATYPE,intent(in) :: avectorin(howmany,botdfconfig:topdfconfig)
  DATATYPE,intent(out) :: avectorout(howmany,botconfig:topconfig)

  avectorout(:,:)=0d0

  do i=botdfconfig,topdfconfig
     avectorout(:,dfincludedconfigs(i))=avectorin(:,i)
  enddo

end subroutine df_transformfrom_local


subroutine basis_transformto_all(howmany,avectorin,avectorout)
  use parameters
  implicit none
  integer :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,numconfig)
  DATATYPE,intent(out) :: avectorout(howmany,numdfbasis)
  DATATYPE :: workvec(howmany,numdfconfigs)

  if (dfrestrictflag.ne.0) then
     if (allspinproject.ne.0) then
        call df_transformto_all(howmany,avectorin(:,:),workvec(:,:))
        call dfspin_transformto_all(howmany,workvec(:,:),avectorout(:,:))
     else
        call df_transformto_all(howmany,avectorin(:,:),avectorout(:,:))
     endif
  else
     if (allspinproject.ne.0) then
        call configspin_transformto_all(howmany,avectorin(:,:),avectorout(:,:))
     else
        avectorout(:,:)=avectorin(:,:)
     endif
  endif

end subroutine basis_transformto_all



subroutine basis_transformto_local(howmany,avectorin,avectorout)
  use parameters
  implicit none
  integer :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,botconfig:topconfig)
  DATATYPE,intent(out) :: avectorout(howmany,botdfbasis:topdfbasis)
  DATATYPE :: workvec(howmany,botdfconfig:topdfconfig)

  if (topdfbasis-botdfbasis+1.ne.0) then
     if (dfrestrictflag.ne.0) then
        if (allspinproject.ne.0) then
           call df_transformto_local(howmany,avectorin(:,:),workvec(:,:))
           call dfspin_transformto_local(howmany,workvec(:,:),avectorout(:,:))
        else
           call df_transformto_local(howmany,avectorin(:,:),avectorout(:,:))
        endif
     else
        if (allspinproject.ne.0) then
           call configspin_transformto_local(howmany,avectorin(:,:),avectorout(:,:))
        else
           avectorout(:,:)=avectorin(:,:)
        endif
     endif
  endif

end subroutine basis_transformto_local



subroutine basis_transformfrom_all(howmany,avectorin,avectorout)
  use parameters
  implicit none
  integer :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,numdfbasis)
  DATATYPE,intent(out) :: avectorout(howmany,numconfig)
  DATATYPE :: workvec(howmany,numdfconfigs)

  if (dfrestrictflag.ne.0) then
     if (allspinproject.ne.0) then
        call dfspin_transformfrom_all(howmany,avectorin(:,:),workvec(:,:))
        call df_transformfrom_all(howmany,workvec(:,:),avectorout(:,:))
     else
        call df_transformfrom_all(howmany,avectorin(:,:),avectorout(:,:))
     endif
  else
     if (allspinproject.ne.0) then
        call configspin_transformfrom_all(howmany,avectorin(:,:),avectorout(:,:))
     else
        avectorout(:,:)=avectorin(:,:)
     endif
  endif

end subroutine basis_transformfrom_all



subroutine basis_transformfrom_local(howmany,avectorin,avectorout)
  use parameters
  implicit none
  integer :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,botdfbasis:topdfbasis)
  DATATYPE,intent(out) :: avectorout(howmany,botconfig:topconfig)
  DATATYPE :: workvec(howmany,botdfconfig:topdfconfig)

  if (topconfig-botconfig+1.ne.0) then
     avectorout(:,:)=0d0
  endif

  if (topdfbasis-botdfbasis+1.ne.0) then
     if (dfrestrictflag.ne.0) then
        if (allspinproject.ne.0) then
           call dfspin_transformfrom_local(howmany,avectorin(:,:),workvec(:,:))
           call df_transformfrom_local(howmany,workvec(:,:),avectorout(:,:))
        else
           call df_transformfrom_local(howmany,avectorin(:,:),avectorout(:,:))
        endif
     else
        if (allspinproject.ne.0) then
           call configspin_transformfrom_local(howmany,avectorin(:,:),avectorout(:,:))
        else
           avectorout(:,:)=avectorin(:,:)
        endif
     endif
  endif

end subroutine basis_transformfrom_local


subroutine fullbasis_transformto_local(howmany,avectorin,avectorout)
  use parameters
  implicit none
  integer :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,botconfig:topconfig)
  DATATYPE,intent(out) :: avectorout(howmany,botbasis:topbasis)

  if (topbasis-botbasis+1.ne.0) then
     if (allspinproject.ne.0) then
        call configspin_transformto_local(howmany,avectorin(:,:),avectorout(:,:))
     else
        avectorout(:,:)=avectorin(:,:)
     endif
  endif

end subroutine fullbasis_transformto_local
  


subroutine fullbasis_transformfrom_local(howmany,avectorin,avectorout)
  use parameters
  implicit none
  integer :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,botbasis:topbasis)
  DATATYPE,intent(out) :: avectorout(howmany,botconfig:topconfig)

  if (topconfig-botconfig+1.ne.0) then
     avectorout(:,:)=0d0
  endif

  if (topbasis-botbasis+1.ne.0) then
     if (allspinproject.ne.0) then
        call configspin_transformfrom_local(howmany,avectorin(:,:),avectorout(:,:))
     else
        avectorout(:,:)=avectorin(:,:)
     endif
  endif

end subroutine fullbasis_transformfrom_local



subroutine checksym(mat,dim)
  use fileptrmod
  implicit none
  integer :: dim,i,j
  real*8 :: mat(dim,dim),sym,asym,tot
  integer, save :: icalled=0

  sym=0; asym=0;tot=0

  do i=1,dim
     do j=1,dim
        tot=tot+ (2*mat(i,j))**2
        sym=sym+ (mat(i,j)+mat(j,i))**2
        asym=asym+ (mat(i,j)-mat(j,i))**2
     enddo
  enddo
  if (tot.gt.1d-10) then
     if (asym.gt.sym*1d-10) then
        icalled=icalled+1
!        if (icalled.lt.10) then
           OFLWR "SYM,ASYM,MAG ", sym/tot,asym/tot,tot; CFLST
!        endif
     endif
  endif

end subroutine checksym


subroutine dferror(avector,outerror,time)
  use xxxmod
  use parameters
  use dfconmod
  implicit none

  integer :: i,imc
  DATATYPE :: avector(numr,firstconfig:lastconfig,mcscfnum), dot
  DATATYPE :: temp1(numr,firstconfig:lastconfig), temp2(numr,firstconfig:lastconfig)
  CNORMTYPE :: outerror
  real *8 :: time

  outerror=0d0
  do imc=1,mcscfnum
     if (numdfconfigs.gt.0) then

        temp1(:,:)=avector(:,:,imc)
        call sparseconfigmult(temp1,temp2,yyy%cptr(0),yyy%sptr(0),1,1,1,1,time)
        do i=firstconfig,lastconfig
           temp1(:,i)=temp2(:,i)*(1-dfincludedmask(i))
        enddo
        
!! reporting norm of error appropriate to the variational principle.

        outerror=   outerror +  dot(temp1,temp1,localnconfig*numr)
     else if (numdfconfigs.lt.0) then
        OFLWR "error, dfconstrain not allocated"; CFLST
     endif
  enddo
  if (parconsplit.ne.0) then
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

end subroutine dferror


#ifdef REALGO
subroutine assigncomplex(realmat,complexf)
  implicit none
  real*8 :: realmat,complexf
  realmat=complexf
end subroutine assigncomplex
subroutine assigncomplexmat(realmat,complexf,m,n)
  implicit none
  integer :: n,m
  real*8 :: realmat(m,n),complexf(m,n)
  realmat(:,:)=complexf(:,:)
end subroutine assigncomplexmat
subroutine assigncomplexvec(realmat,complexf,m)
  implicit none
  integer :: m
  real*8 :: realmat(m),complexf(m)
  realmat(:)=complexf(:)
end subroutine assigncomplexvec
subroutine assignrealvec(complexf,realmat,m)
  implicit none
  integer :: m
  real*8 :: realmat(m),complexf(m)
  complexf(:)=realmat(:)
end subroutine assignrealvec

#else

subroutine assigncomplex(realmat,complexf)
  implicit none
  complex*16 :: complexf
  real*8 :: realmat(2,2)
  realmat(1,1)=real(complexf,8);  realmat(2,2)=real(complexf,8)
  realmat(2,1)=imag(complexf);  realmat(1,2)=(-1)*imag(complexf)
end subroutine assigncomplex

subroutine assigncomplexmat(realmat,complexf,m,n)
  implicit none
  integer :: n,m
  complex*16 :: complexf(m,n)
  real*8 :: realmat(2,m,2,n)
  realmat(1,:,1,:)=real(complexf(:,:),8);  realmat(2,:,2,:)=real(complexf(:,:),8)
  realmat(2,:,1,:)=imag(complexf(:,:));  realmat(1,:,2,:)=(-1)*imag(complexf(:,:))
end subroutine assigncomplexmat

subroutine assigncomplexvec(realmat,complexf,m)
  implicit none
  integer :: m
  complex*16 :: complexf(m)
  real*8 :: realmat(2,m)
  realmat(1,:)=real(complexf(:),8);  realmat(2,:)=imag(complexf(:))
end subroutine assigncomplexvec
subroutine assignrealvec(complexf,realmat,m)
  implicit none
  integer :: m
  complex*16 :: complexf(m)
  real*8 :: realmat(2,m)
  complexf(:)=realmat(1,:)+realmat(2,:)*(0d0,1d0)
end subroutine assignrealvec

#endif

!!  g_ij = < phi_i| d/dt | phi_j >    , that's conmatel, should be hermitian
!!
!!   where i is virtual ("excluded") and j is active ("included")   (other g_ij undefined, 0=0)
!! hermitian g matrix can then be immediattely gotten by fiat
!!
!!

subroutine get_dfconstraint(time)
  use parameters
  use xxxmod
  use mpimod
  use dfconmod
  implicit none

  DATATYPE ::  dot,tempmatel(nspf,nspf)
  integer, save :: times(20)=0, icalled=0
  integer ::    i,     j,  ii, lwork,isize,ishell,iiyy,maxii,imc,itime,jtime,getlen
  integer :: ipairs(2,nspf*(nspf-1))
  DATATYPE, allocatable :: avectorp(:,:), rhs(:,:), avector(:,:), rhomat(:,:,:,:), &
       smallwalkvects(:,:,:,:)   !! numconfig,numr,alpha,beta  alpha=included beta=excluded
  real*8, allocatable :: rhomatpairs(:,:,:,:),rhspairs(:,:), rhspairstemp(:,:), rhspairsbig(:,:), &
       projector(:,:,:,:,:), realrhomat(:,:,:,:,:,:), bigprojector(:,:,:,:,:), rhomatpairsbig(:,:,:,:), &
       realrhs(:,:,:), sing(:), rwork(:), rhomatpairscopy(:,:,:,:), rhomatpairsbigcopy(:,:,:,:), &
       rhomatpairscopy2(:,:,:,:), rhspairstemp2(:,:),&
       mattemp(:,:),pseudoinv(:,:), rhspairsbigtemp(:,:), walknorm(:), temppairs(:,:),temppairsbig(:,:)
  complex*16, allocatable :: work(:)
#ifdef REALGO
  integer, parameter :: zzz=1
#else
  integer, parameter :: zzz=2
#endif
  real*8 :: time,   myconjg(zzz,zzz)  , oneone(zzz,zzz), iimag(zzz,zzz), myperm(zzz,zzz)
  real*8 :: rrcond,rsum

  if (drivingflag.ne.0) then
     OFLWR "Driving with dfconstraint not implemented yet"; CFLST
  endif

#ifdef REALGO
  myconjg(:,:) = 1;  oneone(:,:) = 1;  iimag(:,:)=0; myperm(:,:)=1

#else
  myconjg(:,:)=RESHAPE((/1,0,0,-1/),(/2,2/));  oneone(:,:)=RESHAPE((/1,0,0,1/),(/2,2/))
  iimag(:,:)=RESHAPE((/0,1,-1,0/),(/2,2/));  myperm(:,:)=RESHAPE((/0,1, 1,0/),(/2,2/))

#endif

  call system_clock(itime)

  if (dfrestrictflag.lt.1) then
     OFLWR "WTF DFRESTRICT IS  ", dfrestrictflag; CFLST
  endif

  icalled=icalled+1

  if (notiming.eq.0.and.myrank.eq.1) then
     if (icalled.eq.1) then
        open(8712,file=timingdir(1:getlen(timingdir)-1)//"/dfconstrain.dat",status="unknown")
        write(8712,*) "DF timings."; close(8712)
     else if (mod(icalled,30).eq.0) then
        open(8712,file=timingdir(1:getlen(timingdir)-1)//"/dfconstrain.dat",status="old", position="append")
        write(8712,'(100I12)') times(1:14)/1000; close(8712)
     endif
  endif

  yyy%cptr(0)%xconmatel(:,:)=0.d0;   yyy%cptr(0)%xconmatelxx(:,:)=0.d0;   yyy%cptr(0)%xconmatelyy(:,:)=0.d0;   yyy%cptr(0)%xconmatelzz(:,:)=0.d0

  allocate(avectorp(numr,firstconfig:lastconfig), avector(numr,firstconfig:lastconfig), &
       smallwalkvects(numr,configstart:configend,nspf,nspf), rhs(nspf,nspf),  rhomat(nspf,nspf,nspf,nspf))

  isize=0
  do ishell=1,numshells-1
     do j=allshelltop(ishell)+1,nspf
        do i=allshelltop(ishell-1)+1,allshelltop(ishell)
           isize=isize+1
           ipairs(:,isize)=(/ i,j /)
        enddo
     enddo
  enddo
     
  allocate(projector(zzz,isize,zzz,nspf,nspf), &
       bigprojector(zzz,2*isize,zzz,nspf,nspf), &
       mattemp(zzz*isize,zzz*isize),pseudoinv(zzz*isize,2*zzz*isize))
  projector=0d0;   bigprojector=0d0;


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

  allocate(realrhs(zzz,nspf,nspf),realrhomat(zzz,nspf,nspf,zzz,nspf,nspf), &
       rhomatpairs(zzz,isize,zzz,isize),rhspairs(zzz,isize),rhspairstemp(zzz,isize), rhomatpairscopy(zzz,isize,zzz,isize), &
       rhomatpairscopy2(zzz,isize,zzz,isize),rhspairstemp2(zzz,isize), &
       rhomatpairsbig(zzz,2*isize,zzz,isize),rhspairsbig(zzz,2*isize),rhspairsbigtemp(zzz,2*isize), &
       rhomatpairsbigcopy(zzz,2*isize,zzz,isize))
  allocate(walknorm(2*isize))
  allocate(temppairs(zzz*isize,zzz*nspf**2),temppairsbig(2*zzz*isize,zzz*nspf**2))


  lwork=100*isize*zzz
  allocate(sing(2*isize*zzz),work(lwork), rwork(lwork))
  rrcond=lioreg

  call system_clock(jtime);   times(1)=times(1)+jtime-itime

  maxii=1
  if (tdflag.ne.0) then
     maxii=4
  endif

  do iiyy=1,maxii

     call system_clock(itime)

     rhomatpairscopy(:,:,:,:)=0d0
     rhspairstemp(:,:)=0d0
     rhomatpairsbigcopy(:,:,:,:)=0d0
     rhspairsbigtemp(:,:)=0d0
     walknorm(:)=0d0

     call system_clock(jtime);     times(1)=times(1)+jtime-itime

     do imc=1,mcscfnum
        
        call system_clock(itime)
        
        avector(:,:)=RESHAPE(yyy%cmfpsivec(astart(imc):aend(imc),0),(/numr,localnconfig/))
        
        call df_project(numr,avector)
        if (allspinproject.ne.0) then
           call configspin_project(numr,avector)  !! should commute
        endif

        call get_smallwalkvects(avector,smallwalkvects,numr,1)
        
        call system_clock(jtime);        times(3)=times(3)+jtime-itime;     itime=jtime

!! THIS TAKES A LONG TIME.        !! should only need to do off diagonal blocks but doing all to check
!!   ...no for general case do all I think (Except for diag, whatever)

        ii=(configend-configstart+1)*numr
        call MYGEMM(CNORMCHAR,'N',nspf**2,nspf**2,ii,DATAONE,smallwalkvects,ii,smallwalkvects,ii,DATAZERO,rhomat,nspf**2)

        call system_clock(jtime);        times(4)=times(4)+jtime-itime;     itime=jtime

        if (sparseconfigflag.ne.0) then
           call mympireduce(rhomat,nspf**4)  !! BAD REDUCE
        endif

        call system_clock(jtime);        times(5)=times(5)+jtime-itime;     itime=jtime
        
        if (conway.eq.2.or.conway.eq.3) then
           rhomat(:,:,:,:)=rhomat(:,:,:,:)/timefac
        endif

        if (conway.eq.1) then
           do ii=1,isize
              walknorm(ii)=walknorm(ii)+sqrt(abs(DOT_PRODUCT(&
                   RESHAPE(smallwalkvects(:,:,ipairs(1,ii),ipairs(2,ii)),(/(configend-configstart+1)*numr/)), &
                   RESHAPE(smallwalkvects(:,:,ipairs(1,ii),ipairs(2,ii)),(/(configend-configstart+1)*numr/))))) 
              walknorm(ii+isize)=walknorm(ii+isize)+sqrt(abs(DOT_PRODUCT(&
                   RESHAPE(smallwalkvects(:,:,ipairs(2,ii),ipairs(1,ii)),(/(configend-configstart+1)*numr/)), &
                   RESHAPE(smallwalkvects(:,:,ipairs(2,ii),ipairs(1,ii)),(/(configend-configstart+1)*numr/))))) 
           enddo
        else
           do ii=1,isize
              rsum=sqrt(abs(DOT_PRODUCT(&
                   RESHAPE(smallwalkvects(:,:,ipairs(1,ii),ipairs(2,ii)),(/(configend-configstart+1)*numr/)), &
                   RESHAPE(smallwalkvects(:,:,ipairs(1,ii),ipairs(2,ii)),(/(configend-configstart+1)*numr/)))) &
                   +abs(DOT_PRODUCT(&
                   RESHAPE(smallwalkvects(:,:,ipairs(2,ii),ipairs(1,ii)),(/(configend-configstart+1)*numr/)), &
                   RESHAPE(smallwalkvects(:,:,ipairs(2,ii),ipairs(1,ii)),(/(configend-configstart+1)*numr/)))))
              walknorm(ii)=walknorm(ii)+rsum
              walknorm(ii+isize)=walknorm(ii+isize)+rsum
           enddo
        endif


        call system_clock(jtime);        times(6)=times(6)+jtime-itime;     itime=jtime
        
        call assigncomplexmat(realrhomat,rhomat, nspf**2,nspf**2)
     
        call system_clock(jtime);        times(7)=times(7)+jtime-itime;     itime=jtime

        call DGEMM('N','N',zzz*isize,zzz*nspf**2,zzz*nspf**2,1d0,projector,zzz*isize,&
             realrhomat,zzz*nspf**2,0d0,temppairs,zzz*isize)
        call DGEMM('N','T',zzz*isize,zzz*isize,zzz*nspf**2,1d0,temppairs,zzz*isize,&
             projector,zzz*isize,0d0,rhomatpairs,zzz*isize)

        call DGEMM('N','N',2*zzz*isize,zzz*nspf**2,zzz*nspf**2,1d0,bigprojector,2*zzz*isize,&
             realrhomat,zzz*nspf**2,0d0,temppairsbig,2*zzz*isize)
        call DGEMM('N','T',2*zzz*isize,zzz*isize,zzz*nspf**2,1d0,temppairsbig,2*zzz*isize,&
             projector,zzz*isize,0d0,rhomatpairsbig,2*zzz*isize)

        call system_clock(jtime);        times(8)=times(8)+jtime-itime;     itime=jtime

!! DO SUMMA IN SPARSECONFIGMULT
        
        if (iiyy.eq.1) then
           call sparseconfigmult(avector,avectorp,yyy%cptr(0),yyy%sptr(0),1,1,0,0,time)
        else
           call sparseconfigpulsemult(avector,avectorp,yyy%cptr(0),yyy%sptr(0),iiyy-1)
        endif

        call system_clock(jtime);        times(9)=times(9)+jtime-itime;     itime=jtime
        
        if (conway.ne.2.and.conway.ne.3) then
           avectorp(:,:)=avectorp(:,:)*timefac   
        endif

        rhs(:,:)=0d0
        do j=1,nspf
           do i=1,nspf
              rhs(i,j)=dot(smallwalkvects(:,:,i,j),avectorp(:,configstart:configend),(configend-configstart+1)*numr)
           enddo
        enddo

        call system_clock(jtime);        times(10)=times(10)+jtime-itime;     itime=jtime

        if (sparseconfigflag.ne.0) then
           call mympireduce(rhs,nspf**2)
        endif

        call system_clock(jtime);        times(11)=times(11)+jtime-itime;     itime=jtime

        call assigncomplexvec(realrhs(:,:,:),rhs(:,:), nspf**2)
        
        rhspairs(:,:)=RESHAPE(MATMUL(RESHAPE(projector,(/zzz*isize,zzz*nspf**2/)), &
             RESHAPE(realrhs(:,:,:),(/zzz*nspf**2,1/))),(/zzz,isize/))
        rhspairsbig(:,:)=RESHAPE(MATMUL(RESHAPE(bigprojector,(/2*zzz*isize,zzz*nspf**2/)), &
             RESHAPE(realrhs(:,:,:),(/zzz*nspf**2,1/))),(/zzz,2*isize/))
        
        if (conway.ne.1.and.conway.ne.3) then
           
           rhomatpairscopy(:,:,:,:)=rhomatpairscopy(:,:,:,:)+rhomatpairs(:,:,:,:)
           rhspairstemp(:,:)=rhspairstemp(:,:)+rhspairs(:,:)
        else
           rhomatpairsbigcopy(:,:,:,:)=rhomatpairsbigcopy(:,:,:,:)+rhomatpairsbig(:,:,:,:)
           rhspairsbigtemp(:,:)=rhspairsbigtemp(:,:)+rhspairsbig(:,:)
           
        endif

        call system_clock(jtime);        times(12)=times(12)+jtime-itime;     
        
     enddo    !!! IMC

     call system_clock(itime)

     walknorm(:)=sqrt(walknorm(:))

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

        call realinvmatsmooth(rhomatpairscopy,isize*zzz,lioreg)
        call dgemv('N',isize*zzz,isize*zzz,1d0,rhomatpairscopy,isize*zzz,rhspairstemp,1,0d0,rhspairs,1)

     else
        
        if (connormflag.ne.0) then
           do ii=1,2*isize
              rhomatpairsbigcopy(:,ii,:,:)=rhomatpairsbigcopy(:,ii,:,:)/max(1d-12,walknorm(ii))
              rhspairsbigtemp(:,ii)=rhspairsbigtemp(:,ii)/max(1d-12,walknorm(ii))
           enddo
        endif

        mattemp(:,:)= &
             MATMUL(TRANSPOSE(RESHAPE(rhomatpairsbigcopy(:,:,:,:),(/2*isize*zzz,isize*zzz/))),RESHAPE(rhomatpairsbigcopy(:,:,:,:),(/2*isize*zzz,isize*zzz/)))

        call realinvmatsmooth(mattemp,isize*zzz,lioreg)
        pseudoinv(:,:)=MATMUL(mattemp,TRANSPOSE(RESHAPE(rhomatpairsbigcopy(:,:,:,:),(/2*isize*zzz,isize*zzz/))))
        call DGEMV('N',zzz*isize,2*zzz*isize,1d0,pseudoinv(:,:),isize*zzz,rhspairsbigtemp,1,0d0,rhspairs(:,:),1)

     endif

     call system_clock(jtime);        times(13)=times(13)+jtime-itime;    itime=jtime

     tempmatel(:,:)=0d0

     do ii=1,isize
        i=ipairs(1,ii);     j=ipairs(2,ii)
#ifdef REALGO
        tempmatel(j,i) = rhspairs(1,ii)    
#else
        tempmatel(j,i) = rhspairs(1,ii) + (0d0,1d0) * rhspairs(2,ii)
#endif
     enddo

    tempmatel(:,:)=  ( tempmatel(:,:) -TRANSPOSE(CONJUGATE(tempmatel(:,:)))  ) / timefac  * condamp

     select case(iiyy)
     case(1)
        yyy%cptr(0)%xconmatel(:,:)=  tempmatel(:,:)
     case(2)
        yyy%cptr(0)%xconmatelxx(:,:)=  tempmatel(:,:)
     case(3)
        yyy%cptr(0)%xconmatelyy(:,:)=  tempmatel(:,:)
     case(4)
        yyy%cptr(0)%xconmatelzz(:,:)=  tempmatel(:,:)
     end select

     call system_clock(jtime);        times(14)=times(14)+jtime-itime;  

  enddo

  deallocate(&
       avectorp,   avector,    &
       smallwalkvects,         rhs,          rhomat,  &
        projector,  bigprojector, mattemp,pseudoinv, &
          realrhs,  realrhomat,  &
       rhomatpairs,  rhspairs,      rhspairstemp,  rhomatpairscopy, &
       rhomatpairscopy2,rhspairstemp2, &
       rhomatpairsbig,  rhspairsbig, rhspairsbigtemp, &
       rhomatpairsbigcopy, &
       walknorm, &
       sing,  work,  rwork,temppairs,temppairsbig)

end subroutine get_dfconstraint





subroutine get_smallwalkvects(avector, smallwalkvects,nblock,howmany)
  use parameters
  use dfconmod
  implicit none

  integer ::    i,ii,ix,nblock,howmany
  DATATYPE,intent(in) ::  avector(nblock,firstconfig:lastconfig,howmany)
  DATATYPE,intent(out) :: smallwalkvects(nblock,configstart:configend,howmany,nspf,nspf)
  DATATYPE, allocatable ::  bigavector(:,:,:)

  if (dfrestrictflag.lt.1) then
     OFLWR "WTF DFRESTRICT IS  ", dfrestrictflag; CFLST
  endif
  
!! DO SUMMA
  allocate(bigavector(nblock,numconfig,howmany))
  bigavector(:,:,:)=0d0
  bigavector(:,firstconfig:lastconfig,:)=avector(:,:,:)

  if (parconsplit.ne.0) then
     do ii=1,howmany
        call mpiallgather(bigavector(:,:,ii),nblock*numconfig,nblock*configsperproc(:),nblock*maxconfigsperproc)
     enddo
  endif

  smallwalkvects(:,:,:,:,:)=0d0
  
  do i=1,numdfwalks
     ii=includedorb(i);     ix=excludedorb(i)
     smallwalkvects(:,dfwalkto(i),:,ii,ix)=smallwalkvects(:,dfwalkto(i),:,ii,ix) + bigavector(:,dfwalkfrom(i),:) * dfwalkphase(i)
  enddo

  deallocate(bigavector)

end subroutine get_smallwalkvects



subroutine get_rhomat(avector, rhomat,nblock,howmany)
  use parameters
  implicit none
  integer ::    ii,nblock,howmany
  DATATYPE,intent(in) ::  avector(nblock,firstconfig:lastconfig,howmany)
  DATATYPE,intent(out) :: rhomat(nspf,nspf,nspf,nspf)
  DATATYPE, allocatable ::  &
       smallwalkvects(:,:,:,:,:)   !! nblock,configstart:configend,howmany,alpha,beta  alpha=included beta=excluded

  if (dfrestrictflag.lt.1) then
     OFLWR "WTF DFRESTRICT IS  ", dfrestrictflag; CFLST
  endif
  
  allocate(smallwalkvects(nblock,configstart:configend,howmany,nspf,nspf))

  call get_smallwalkvects(avector,smallwalkvects,nblock,howmany)

!! THIS TAKES A LONG TIME.
  
  ii=(configend-configstart+1)*nblock*howmany
  call MYGEMM(CNORMCHAR,'N',nspf**2,nspf**2,ii,DATAONE,smallwalkvects,ii,smallwalkvects,ii,DATAZERO,rhomat,nspf**2)

  if (sparseconfigflag.ne.0) then
     call mympireduce(rhomat,nspf**4)  !! BAD REDUCE
  endif

  deallocate(smallwalkvects)

end subroutine get_rhomat

