
#include "Definitions.INC"

!! INDICATES WHETHER A GIVEN CONFIG is included in the restricted configuration list.

!! it is only included if all double excitations from it are included in the total configuration list.



function getdfindex(www,config1)
  use walkmod
  implicit none
  type(walktype) :: www
  integer :: getdfindex, config1,j,flag, dir,newdir, step

  getdfindex=-1

  if (www%ddd%dfincludedmask(config1).eq.0) then
     return
  endif
  
  dir=1;  j=1;  step=max(1,www%numdfconfigs/4);  flag=0

  do while (flag.eq.0)
     flag=1

     if (www%ddd%dfincludedconfigs(j) .ne. config1) then
        flag=0
        if (www%ddd%dfincludedconfigs(j).lt.config1) then
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
        if (j.ge.www%numdfconfigs) then
           j=www%numdfconfigs;              dir=-1
        endif
     endif
  enddo

  getdfindex=j
  
end function getdfindex


subroutine basis_set(www)
  use r_parameters
  use walkmod
  use mpimod !! nprocs
  implicit none
  type(walktype),intent(inout) :: www

  allocate(www%basisperproc(nprocs),www%dfbasisperproc(nprocs))

  if (www%allspinproject.eq.0) then
     www%numbasis=www%numconfig
     www%basisperproc(:)=www%configsperproc(:)
     www%maxbasisperproc=www%maxconfigsperproc
     www%botbasis=www%botconfig;  www%topbasis=www%topconfig

     www%numdfbasis=www%numdfconfigs
     www%dfbasisperproc(:)=www%dfconfsperproc(:)
     www%maxdfbasisperproc=www%maxdfconfsperproc
     www%botdfbasis=www%botdfconfig;  www%topdfbasis=www%topdfconfig
        
  else
     www%numbasis=www%sss%numspinconfig
     www%basisperproc(:)=www%sss%spinsperproc(:)
     www%maxbasisperproc=www%sss%maxspinsperproc
     www%botbasis=www%sss%botspin;  www%topbasis=www%sss%topspin

     www%numdfbasis=www%sss%numspindfconfig
     www%dfbasisperproc(:)=www%sss%spindfsperproc(:)
     www%maxdfbasisperproc=www%sss%maxspindfsperproc
     www%botdfbasis=www%sss%botdfspin;  www%topdfbasis=www%sss%topdfspin
  endif

  www%totadim=www%localnconfig*numr
end subroutine basis_set


! "included" configurations are those kept in the calculation.
! there are numdfconfigs INCLUDED configurations.

subroutine init_dfcon(www)
  use fileptrmod
  use walkmod
  use mpimod
  implicit none
  type(walktype) :: www
  integer :: i,j,iconfig,nondfconfigs,dfrank,ii
  logical :: allowedconfig0
  integer,allocatable :: dfnotconfigs(:)

  allocate(www%ddd%dfincludedmask(www%numconfig), www%ddd%dfincludedconfigs(www%numconfig), &
       www%ddd%dfincludedindex(www%numconfig),&
       dfnotconfigs(www%numconfig),  www%dfconfsperproc(nprocs), &
       www%allbotdfconfigs(nprocs),www%alltopdfconfigs(nprocs))

  if (www%dfrestrictflag.lt.www%dflevel) then
     OFLWR "error, set dfrestrictflag .ge. dflevel",www%dfrestrictflag,www%dflevel; CFLST
  endif

  if (www%dfrestrictflag.eq.www%dflevel) then
     www%ddd%dfincludedmask(:)=1; dfnotconfigs(:)=(-1)
     do i=1,www%numconfig
        www%ddd%dfincludedconfigs(i)=i
        www%ddd%dfincludedindex(i)=i
     enddo
     www%numdfconfigs=www%numconfig; nondfconfigs=0
     www%dfconfsperproc(:)=www%configsperproc(:)
     www%allbotdfconfigs(:)=www%allbotconfigs(:)
     www%alltopdfconfigs(:)=www%alltopconfigs(:)
     www%maxdfconfsperproc=www%maxconfigsperproc
     www%botdfconfig=www%botconfig;  www%topdfconfig=www%topconfig

  else

     www%ddd%dfincludedmask(:)=0;  www%ddd%dfincludedconfigs(:)=(-1);  dfnotconfigs(:)=(-1)
     www%ddd%dfincludedindex(:)=(-999)

     www%numdfconfigs=0;  nondfconfigs=0

     do i=1,www%numconfig
        if (allowedconfig0(www,www%configlist(:,i),www%dfrestrictflag)) then 
           www%numdfconfigs=www%numdfconfigs+1
           www%ddd%dfincludedmask(i)=1;        www%ddd%dfincludedconfigs(www%numdfconfigs)=i
           www%ddd%dfincludedindex(i)=www%numdfconfigs
        else
           nondfconfigs=nondfconfigs+1;        dfnotconfigs(nondfconfigs)=i
        endif
     enddo

     if (www%numdfconfigs.eq.0) then
        OFLWR "Error, no included DF configs!"; CFLST
     endif

     OFLWR "Number of included configurations ", www%numdfconfigs
     WRFL  "Number of not included ", nondfconfigs;  WRFL; CFL

     dfrank=0
     do i=www%botconfig,www%topconfig
        if (allowedconfig0(www,www%configlist(:,i),www%dfrestrictflag)) then 
           dfrank=dfrank+1
        endif
     enddo

     www%dfconfsperproc(myrank)=dfrank
     www%maxdfconfsperproc=0
     ii=0
     do i=1,nprocs
        www%allbotdfconfigs(i)=ii+1
        call mympiibcastone(www%dfconfsperproc(i),i)
        ii=ii+www%dfconfsperproc(i)
        www%alltopdfconfigs(i)=ii
        if (www%dfconfsperproc(i).gt.www%maxdfconfsperproc) then
           www%maxdfconfsperproc=www%dfconfsperproc(i)
        endif
     enddo
     www%botdfconfig=www%allbotdfconfigs(myrank)
     www%topdfconfig=www%alltopdfconfigs(myrank)

  endif

  www%ddd%numdfwalks=0

  if (www%dfrestrictflag.eq.www%dflevel) then

     allocate(www%ddd%dfwalkfrom(1), www%ddd%dfwalkto(1), www%ddd%includedorb(1), &
          www%ddd%excludedorb(1), www%ddd%dfwalkphase(1))
     www%ddd%dfwalkfrom(:)=(-1); www%ddd%dfwalkto(:)=(-1); www%ddd%includedorb(:)=(-1)
     www%ddd%excludedorb(:)=(-1); www%ddd%dfwalkphase(:)=(-798)

  else

     do i=1,NONdfconfigs
        iconfig=dfNOTconfigs(i)
        if (iconfig.ge.www%configstart.and.iconfig.le.www%configend) then
           do j=1,www%numsinglewalks(iconfig)
              if (allowedconfig0(www,www%configlist(:,www%singlewalk(j,iconfig)),www%dfrestrictflag)) then
                 www%ddd%numdfwalks=www%ddd%numdfwalks+1
              endif
           enddo
        endif
     enddo

     OFLWR "numdfwalks:  ", www%ddd%numdfwalks," is total number DF single excitations on this processor"; WRFL; CFL

     allocate(www%ddd%dfwalkfrom(www%ddd%numdfwalks), www%ddd%dfwalkto(www%ddd%numdfwalks),&
          www%ddd%includedorb(www%ddd%numdfwalks), www%ddd%excludedorb(www%ddd%numdfwalks), &
          www%ddd%dfwalkphase(www%ddd%numdfwalks))

     www%ddd%numdfwalks=0
     do i=1,NONdfconfigs
        iconfig=dfNOTconfigs(i)
        if (iconfig.ge.www%configstart.and.iconfig.le.www%configend) then
           do j=1,www%numsinglewalks(iconfig)
              if (allowedconfig0(www,www%configlist(:,www%singlewalk(j,iconfig)),www%dfrestrictflag)) then
                 www%ddd%numdfwalks=www%ddd%numdfwalks+1
              
                 www%ddd%dfwalkTO(www%ddd%numdfwalks)=iconfig
                 www%ddd%dfwalkFROM(www%ddd%numdfwalks)=www%singlewalk(j,iconfig)
                 www%ddd%EXcludedorb(www%ddd%numdfwalks)=www%singlewalkopspf(1,j,iconfig)
                 www%ddd%INcludedorb(www%ddd%numdfwalks)=www%singlewalkopspf(2,j,iconfig)
              
                 www%ddd%dfwalkphase(www%ddd%numdfwalks)=www%singlewalkdirphase(j,iconfig)
              endif
           enddo
        endif
     enddo

     OFLWR "numdfwalks:  ", www%ddd%numdfwalks," is really number DF single excitations on this processor"; WRFL; CFL
     
  endif

  deallocate(dfnotconfigs)

end subroutine init_dfcon

subroutine dfcondealloc(www)
  use walkmod
  implicit none
  type(walktype) :: www
  deallocate(www%ddd%dfincludedmask,www%ddd%dfincludedconfigs)
  deallocate(www%ddd%dfwalkfrom, www%ddd%dfwalkto, www%ddd%includedorb, &
       www%ddd%excludedorb,www%ddd%dfwalkphase)
  www%ddd%numdfwalks=-1;  www%numdfconfigs=-1
end subroutine dfcondealloc


subroutine df_project(www,howmany,avector)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: howmany
  DATATYPE :: avector(howmany,www%firstconfig:www%lastconfig)

  if (www%parconsplit.eq.0) then
     call df_project_all(www,howmany,avector)
  else
     call df_project_local(www,howmany,avector)
  endif

end subroutine df_project


subroutine df_project_all(www,howmany,avector)
  use walkmod
  implicit none
  integer :: howmany,i
  type(walktype),intent(in) :: www
  DATATYPE :: avector(howmany,www%numconfig)
  do i=1,www%numconfig
     avector(:,i)=avector(:,i)*www%ddd%dfincludedmask(i)
  enddo
end subroutine df_project_all


subroutine df_project_local(www,howmany,avector)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: howmany,i
  DATATYPE :: avector(howmany,www%botconfig:www%topconfig)
  do i=www%botconfig,www%topconfig
     avector(:,i)=avector(:,i)*www%ddd%dfincludedmask(i)
  enddo
end subroutine df_project_local


subroutine df_transformto_local(www,howmany,avectorin,avectorout)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: howmany,i
  DATATYPE,intent(in) :: avectorin(howmany,www%botconfig:www%topconfig)
  DATATYPE,intent(out) :: avectorout(howmany,www%botdfconfig:www%topdfconfig)

  do i=www%botdfconfig,www%topdfconfig
     avectorout(:,i)=avectorin(:,www%ddd%dfincludedconfigs(i))
  enddo

end subroutine df_transformto_local



subroutine df_transformfrom_local(www,howmany,avectorin,avectorout)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: howmany,i
  DATATYPE,intent(in) :: avectorin(howmany,www%botdfconfig:www%topdfconfig)
  DATATYPE,intent(out) :: avectorout(howmany,www%botconfig:www%topconfig)

  avectorout(:,:)=0d0

  do i=www%botdfconfig,www%topdfconfig
     avectorout(:,www%ddd%dfincludedconfigs(i))=avectorin(:,i)
  enddo

end subroutine df_transformfrom_local


subroutine basis_project(www,howmany,avector)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: howmany
  DATATYPE,intent(inout) :: avector(howmany,www%firstconfig:www%lastconfig)

  if (www%allspinproject==1) then
     call configspin_project(www,howmany,avector)
  endif
  if (www%dfrestrictflag.ne.www%dflevel) then
     call df_project(www,howmany,avector)
  endif

end subroutine basis_project



subroutine basis_transformto_all(www,howmany,avectorin,avectorout)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,www%numconfig)
  DATATYPE,intent(out) :: avectorout(howmany,www%numdfbasis)

  call basis_transformto_local(www,howmany,avectorin(:,www%botconfig),avectorout(:,www%botdfbasis))

  call mpiallgather(avectorout(:,:),www%numdfbasis*howmany,www%dfbasisperproc(:)*howmany,&
       www%maxdfbasisperproc*howmany)

end subroutine basis_transformto_all



subroutine basis_transformto_local(www,howmany,avectorin,avectorout)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,www%botconfig:www%topconfig)
  DATATYPE,intent(out) :: avectorout(howmany,www%botdfbasis:www%topdfbasis)
  DATATYPE :: workvec(howmany,www%botdfconfig:www%topdfconfig)

  if (www%topdfbasis-www%botdfbasis+1.ne.0) then
     if (www%dfrestrictflag.ne.www%dflevel) then
        if (www%allspinproject.ne.0) then
           call df_transformto_local(www,howmany,avectorin(:,:),workvec(:,:))
           call dfspin_transformto_local(www,howmany,workvec(:,:),avectorout(:,:))
        else
           call df_transformto_local(www,howmany,avectorin(:,:),avectorout(:,:))
        endif
     else
        if (www%allspinproject.ne.0) then
           call configspin_transformto_local(www,howmany,avectorin(:,:),avectorout(:,:))
        else
           avectorout(:,:)=avectorin(:,:)
        endif
     endif
  endif

end subroutine basis_transformto_local



subroutine basis_transformfrom_all(www,howmany,avectorin,avectorout)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,www%numdfbasis)
  DATATYPE,intent(out) :: avectorout(howmany,www%numconfig)

  call basis_transformfrom_local(www,howmany,avectorin(:,www%botdfbasis),avectorout(:,www%botconfig))

  call mpiallgather(avectorout(:,:),www%numconfig*howmany,www%configsperproc(:)*howmany,&
       www%maxconfigsperproc*howmany)

end subroutine basis_transformfrom_all



subroutine basis_transformfrom_local(www,howmany,avectorin,avectorout)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,www%botdfbasis:www%topdfbasis)
  DATATYPE,intent(out) :: avectorout(howmany,www%botconfig:www%topconfig)
  DATATYPE :: workvec(howmany,www%botdfconfig:www%topdfconfig)

  if (www%topconfig-www%botconfig+1.ne.0) then
     avectorout(:,:)=0d0
  endif

  if (www%topdfbasis-www%botdfbasis+1.ne.0) then
     if (www%dfrestrictflag.ne.www%dflevel) then
        if (www%allspinproject.ne.0) then
           call dfspin_transformfrom_local(www,howmany,avectorin(:,:),workvec(:,:))
           call df_transformfrom_local(www,howmany,workvec(:,:),avectorout(:,:))
        else
           call df_transformfrom_local(www,howmany,avectorin(:,:),avectorout(:,:))
        endif
     else
        if (www%allspinproject.ne.0) then
           call configspin_transformfrom_local(www,howmany,avectorin(:,:),avectorout(:,:))
        else
           avectorout(:,:)=avectorin(:,:)
        endif
     endif
  endif

end subroutine basis_transformfrom_local




subroutine fullbasis_transformto_local(www,howmany,avectorin,avectorout)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,www%botconfig:www%topconfig)
  DATATYPE,intent(out) :: avectorout(howmany,www%botbasis:www%topbasis)

  if (www%topbasis-www%botbasis+1.ne.0) then
     if (www%allspinproject.ne.0) then
        call configspin_transformto_local(www,howmany,avectorin(:,:),avectorout(:,:))
     else
        avectorout(:,:)=avectorin(:,:)
        endif
  endif

end subroutine fullbasis_transformto_local


subroutine fullbasis_transformfrom_local(www,howmany,avectorin,avectorout)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer :: howmany
  DATATYPE,intent(in) :: avectorin(howmany,www%botbasis:www%topbasis)
  DATATYPE,intent(out) :: avectorout(howmany,www%botconfig:www%topconfig)

  if (www%topconfig-www%botconfig+1.ne.0) then
     avectorout(:,:)=0d0
  endif

  if (www%topbasis-www%botbasis+1.ne.0) then
     if (www%allspinproject.ne.0) then
        call configspin_transformfrom_local(www,howmany,avectorin(:,:),avectorout(:,:))
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


subroutine dferror(www,cptr,sptr,avector,numvects,outerror,time)
  use configptrmod
  use sparseptrmod
  use walkmod
  use r_parameters
  use fileptrmod
  implicit none
  integer,intent(in) :: numvects
  type(walktype),intent(in) :: www
  type(configptr),intent(in) :: cptr
  type(sparseptr),intent(in) :: sptr
  integer :: i,imc
  DATATYPE :: avector(numr,www%firstconfig:www%lastconfig,numvects), dot
  DATATYPE :: temp1(numr,www%firstconfig:www%lastconfig), temp2(numr,www%firstconfig:www%lastconfig)
  CNORMTYPE :: outerror
  real *8 :: time

  outerror=0d0
  do imc=1,numvects
     if (www%numdfconfigs.gt.0) then

        temp1(:,:)=avector(:,:,imc)
        call sparseconfigmult(www,temp1,temp2,cptr,sptr,1,1,1,1,time,imc)
        do i=www%firstconfig,www%lastconfig
           temp1(:,i)=temp2(:,i)*(1-www%ddd%dfincludedmask(i))
        enddo
        
!! reporting norm of error appropriate to the variational principle.

        outerror=   outerror +  dot(temp1,temp1,www%totadim)
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
  use configmod
  implicit none
  real*8, intent(in) :: time
  call get_dfconstraint0(yyy%cmfpsivec(astart(1):aend(mcscfnum),0),mcscfnum,yyy%cptr(0),yyy%sptr(0),www,time)
end subroutine get_dfconstraint


subroutine get_dfconstraint0(inavectors,numvects,cptr,sptr,www,time)
  use fileptrmod
  use timing_parameters
  use sparse_parameters
  use basis_parameters
  use r_parameters
  use constraint_parameters
  use ham_parameters
  use configptrmod
  use sparseptrmod
  use walkmod
  use mpimod
  implicit none
  integer,intent(in) :: numvects
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: inavectors(numr,www%firstconfig:www%lastconfig,numvects)
  type(CONFIGPTR),intent(inout) :: cptr
  type(SPARSEPTR),intent(in) :: sptr
  DATATYPE ::  dot,tempmatel(www%nspf,www%nspf)
  integer, save :: times(20)=0, icalled=0
  integer ::    i,     j,  ii, lwork,isize,ishell,iiyy,maxii,imc,itime,jtime,getlen
  integer :: ipairs(2,www%nspf*(www%nspf-1))
  DATATYPE, allocatable :: avectorp(:,:), rhs(:,:), avector(:,:), rhomat(:,:,:,:), &
       smallwalkvects(:,:,:,:)   !! numconfig,numr,alpha,beta  alpha=included beta=excluded
  real*8, allocatable :: rhomatpairs(:,:,:,:),rhspairs(:,:), rhspairstemp(:,:), rhspairsbig(:,:), &
       projector(:,:,:,:,:), realrhomat(:,:,:,:,:,:), bigprojector(:,:,:,:,:), rhomatpairsbig(:,:,:,:), &
       realrhs(:,:,:), sing(:), rwork(:), rhomatpairscopy(:,:,:,:), rhomatpairsbigcopy(:,:,:,:), &
       rhomatpairscopy2(:,:,:,:), rhspairstemp2(:,:),&
       mattemp(:,:),pseudoinv(:,:), rhspairsbigtemp(:,:), temppairs(:,:),temppairsbig(:,:)
  complex*16, allocatable :: work(:)
#ifdef REALGO
  integer, parameter :: zzz=1
#else
  integer, parameter :: zzz=2
#endif
  real*8 :: time,   myconjg(zzz,zzz)  , oneone(zzz,zzz), iimag(zzz,zzz), myperm(zzz,zzz)
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

  call system_clock(itime)

!  OFLWR "*** GO CON ****"
!  WRFL real(inavectors)
!  WRFL;  CFL


  if (www%dfrestrictflag.le.www%dflevel) then
     OFLWR "WTF DFRESTRICT IS  ", www%dfrestrictflag,www%dflevel; CFLST
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

  cptr%xconmatel(:,:)=0.d0;   cptr%xconmatelxx(:,:)=0.d0;   cptr%xconmatelyy(:,:)=0.d0;   cptr%xconmatelzz(:,:)=0.d0

  allocate(avectorp(numr,www%firstconfig:www%lastconfig), avector(numr,www%firstconfig:www%lastconfig), &
       smallwalkvects(numr,www%configstart:www%configend,www%nspf,www%nspf), rhs(www%nspf,www%nspf),  rhomat(www%nspf,www%nspf,www%nspf,www%nspf))

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

  allocate(realrhs(zzz,www%nspf,www%nspf),realrhomat(zzz,www%nspf,www%nspf,zzz,www%nspf,www%nspf), &
       rhomatpairs(zzz,isize,zzz,isize),rhspairs(zzz,isize),rhspairstemp(zzz,isize), rhomatpairscopy(zzz,isize,zzz,isize), &
       rhomatpairscopy2(zzz,isize,zzz,isize),rhspairstemp2(zzz,isize), &
       rhomatpairsbig(zzz,2*isize,zzz,isize),rhspairsbig(zzz,2*isize),rhspairsbigtemp(zzz,2*isize), &
       rhomatpairsbigcopy(zzz,2*isize,zzz,isize))
  allocate(temppairs(zzz*isize,zzz*www%nspf**2),temppairsbig(2*zzz*isize,zzz*www%nspf**2))


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

     call system_clock(jtime);     times(1)=times(1)+jtime-itime

     do imc=1,numvects
        
        call system_clock(itime)
        
        avector(:,:)=inavectors(:,:,imc)

        call basis_project(www,numr,avector)

        call get_smallwalkvects(www,avector,smallwalkvects,numr,1)
        
        call system_clock(jtime);        times(3)=times(3)+jtime-itime;     itime=jtime

!! THIS TAKES A LONG TIME.        !! should only need to do off diagonal blocks but doing all to check
!!   ...no for general case do all I think (Except for diag, whatever)

        ii=(www%configend-www%configstart+1)*numr
        call MYGEMM(CNORMCHAR,'N',www%nspf**2,www%nspf**2,ii,DATAONE,smallwalkvects,ii,smallwalkvects,ii,DATAZERO,rhomat,www%nspf**2)

        call system_clock(jtime);        times(4)=times(4)+jtime-itime;     itime=jtime

        if (sparseconfigflag.ne.0) then
           call mympireduce(rhomat,www%nspf**4)  !! BAD REDUCE
        endif

        call system_clock(jtime);        times(5)=times(5)+jtime-itime;     itime=jtime
        
        if (conway.eq.2.or.conway.eq.3) then
           rhomat(:,:,:,:)=rhomat(:,:,:,:)/timefac
        endif

        call system_clock(jtime);        times(6)=times(6)+jtime-itime;     itime=jtime
        
        call assigncomplexmat(realrhomat,rhomat, www%nspf**2,www%nspf**2)
     
        call system_clock(jtime);        times(7)=times(7)+jtime-itime;     itime=jtime

        call DGEMM('N','N',zzz*isize,zzz*www%nspf**2,zzz*www%nspf**2,1d0,projector,zzz*isize,&
             realrhomat,zzz*www%nspf**2,0d0,temppairs,zzz*isize)
        call DGEMM('N','T',zzz*isize,zzz*isize,zzz*www%nspf**2,1d0,temppairs,zzz*isize,&
             projector,zzz*isize,0d0,rhomatpairs,zzz*isize)

        call DGEMM('N','N',2*zzz*isize,zzz*www%nspf**2,zzz*www%nspf**2,1d0,bigprojector,2*zzz*isize,&
             realrhomat,zzz*www%nspf**2,0d0,temppairsbig,2*zzz*isize)
        call DGEMM('N','T',2*zzz*isize,zzz*isize,zzz*www%nspf**2,1d0,temppairsbig,2*zzz*isize,&
             projector,zzz*isize,0d0,rhomatpairsbig,2*zzz*isize)

        call system_clock(jtime);        times(8)=times(8)+jtime-itime;     itime=jtime

        if (iiyy.eq.1) then
           call sparseconfigmult(www,avector,avectorp,cptr,sptr,1,1,0,0,time,imc)
        else
           call sparseconfigpulsemult(www,avector,avectorp,cptr,sptr,iiyy-1,imc)
        endif

        call system_clock(jtime);        times(9)=times(9)+jtime-itime;     itime=jtime
        
        if (conway.ne.2.and.conway.ne.3) then
           avectorp(:,:)=avectorp(:,:)*timefac   
        endif

        rhs(:,:)=0d0
        do j=1,www%nspf
           do i=1,www%nspf
              rhs(i,j)=dot(smallwalkvects(:,:,i,j),avectorp(:,www%configstart:www%configend),(www%configend-www%configstart+1)*numr)
           enddo
        enddo

        call system_clock(jtime);        times(10)=times(10)+jtime-itime;     itime=jtime

        if (sparseconfigflag.ne.0) then
           call mympireduce(rhs,www%nspf**2)
        endif

!        OFLWR "RHS"
!        WRFL real(rhs)
!        WRFL; CFL

        call system_clock(jtime);        times(11)=times(11)+jtime-itime;     itime=jtime

        call assigncomplexvec(realrhs(:,:,:),rhs(:,:), www%nspf**2)
        
        rhspairs(:,:)=RESHAPE(MATMUL(RESHAPE(projector,(/zzz*isize,zzz*www%nspf**2/)), &
             RESHAPE(realrhs(:,:,:),(/zzz*www%nspf**2,1/))),(/zzz,isize/))
        rhspairsbig(:,:)=RESHAPE(MATMUL(RESHAPE(bigprojector,(/2*zzz*isize,zzz*www%nspf**2/)), &
             RESHAPE(realrhs(:,:,:),(/zzz*www%nspf**2,1/))),(/zzz,2*isize/))
        
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
        
        mattemp(:,:)= &
             MATMUL(TRANSPOSE(RESHAPE(rhomatpairsbigcopy(:,:,:,:),(/2*isize*zzz,isize*zzz/))),RESHAPE(rhomatpairsbigcopy(:,:,:,:),(/2*isize*zzz,isize*zzz/)))

!        OFLWR "MATTEMP",lioreg
!        WRFL real(mattemp)
!        WRFL; CFL

        call realinvmatsmooth(mattemp,isize*zzz,lioreg)

!        OFLWR "MATTEMPnow",lioreg
!        WRFL real(mattemp)
!        WRFL; CFL

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
        cptr%xconmatel(:,:)=  tempmatel(:,:)

!        OFLWR "****CON",cptr%xconmatel(1,3); CFL

     case(2)
        cptr%xconmatelxx(:,:)=  tempmatel(:,:)
     case(3)
        cptr%xconmatelyy(:,:)=  tempmatel(:,:)
     case(4)
        cptr%xconmatelzz(:,:)=  tempmatel(:,:)
     end select

     call system_clock(jtime);        times(14)=times(14)+jtime-itime;  

  enddo

!!  OFLWR "TEMPSTOP" ; CFLST

  deallocate(&
       avectorp,   avector,    &
       smallwalkvects,         rhs,          rhomat,  &
        projector,  bigprojector, mattemp,pseudoinv, &
          realrhs,  realrhomat,  &
       rhomatpairs,  rhspairs,      rhspairstemp,  rhomatpairscopy, &
       rhomatpairscopy2,rhspairstemp2, &
       rhomatpairsbig,  rhspairsbig, rhspairsbigtemp, &
       rhomatpairsbigcopy, &
       sing,  work,  rwork,temppairs,temppairsbig)

end subroutine get_dfconstraint0





subroutine get_smallwalkvects(www,avector, smallwalkvects,nblock,howmany)
  use fileptrmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer ::    i,ii,ix,nblock,howmany
  DATATYPE,intent(in) ::  avector(nblock,www%firstconfig:www%lastconfig,howmany)
  DATATYPE,intent(out) :: smallwalkvects(nblock,www%configstart:www%configend,howmany,www%nspf,www%nspf)
  DATATYPE, allocatable ::  bigavector(:,:,:)

  if (www%dfrestrictflag.le.www%dflevel) then
     OFLWR "WTF DFRESTRICT IS  ", www%dfrestrictflag,www%dflevel; CFLST
  endif
  
!! DO SUMMA (parconsplit.ne.0 and sparsesummaflag.eq.2, "circ")
!! AND DO HOPS

  allocate(bigavector(nblock,www%numconfig,howmany))
  bigavector(:,:,:)=0d0
  bigavector(:,www%firstconfig:www%lastconfig,:)=avector(:,:,:)

  if (www%parconsplit.ne.0) then
     do ii=1,howmany
        call mpiallgather(bigavector(:,:,ii),nblock*www%numconfig,nblock*www%configsperproc(:),nblock*www%maxconfigsperproc)
     enddo
  endif

  smallwalkvects(:,:,:,:,:)=0d0
  
  do i=1,www%ddd%numdfwalks
     ii=www%ddd%includedorb(i);     ix=www%ddd%excludedorb(i)
     smallwalkvects(:,www%ddd%dfwalkto(i),:,ii,ix)=smallwalkvects(:,www%ddd%dfwalkto(i),:,ii,ix) + bigavector(:,www%ddd%dfwalkfrom(i),:) * www%ddd%dfwalkphase(i)
  enddo

  deallocate(bigavector)

end subroutine get_smallwalkvects



subroutine get_rhomat(www,avector, rhomat,nblock,howmany)
  use fileptrmod
  use sparse_parameters
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer ::    ii,nblock,howmany
  DATATYPE,intent(in) ::  avector(nblock,www%firstconfig:www%lastconfig,howmany)
  DATATYPE,intent(out) :: rhomat(www%nspf,www%nspf,www%nspf,www%nspf)
  DATATYPE, allocatable ::  &
       smallwalkvects(:,:,:,:,:)   !! nblock,configstart:configend,howmany,alpha,beta  alpha=included beta=excluded

  if (www%dfrestrictflag.le.www%dflevel) then
     OFLWR "WTF DFRESTRICT IS  ", www%dfrestrictflag,www%dflevel; CFLST
  endif
  
  allocate(smallwalkvects(nblock,www%configstart:www%configend,howmany,www%nspf,www%nspf))

  call get_smallwalkvects(www,avector,smallwalkvects,nblock,howmany)

!! THIS TAKES A LONG TIME.
  
  ii=(www%configend-www%configstart+1)*nblock*howmany
  call MYGEMM(CNORMCHAR,'N',www%nspf**2,www%nspf**2,ii,DATAONE,smallwalkvects,ii,smallwalkvects,ii,DATAZERO,rhomat,www%nspf**2)

  if (sparseconfigflag.ne.0) then
     call mympireduce(rhomat,www%nspf**4)  !! BAD REDUCE
  endif

  deallocate(smallwalkvects)

end subroutine get_rhomat

