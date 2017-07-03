
!! ALL MODULES except for all_matel()

!! ORBITAL AND CONFIGURATION MATRIX ELEMENTS.  TWO ELECTRON REDUCED MATRIX ELEMENTS
!! ARE CALCULATED HERE ALONG THE WAY.  GET_REDUCEDHAM MUST BE CALLED AFTER TWOE_MATEL
!! BECAUSE IT NEEDS THOSE REDUCED MATRIX ELEMENTS

#include "Definitions.INC"

module allmat0mod
contains

subroutine all_matel0(inholeflag,matrix_ptr,inspfs1,inspfs2,twoereduced,times,firstspf,lastspf)
  use parameters
  use configptrmod
  implicit none
  Type(CONFIGPTR),intent(inout) :: matrix_ptr
  integer,intent(in) :: firstspf,lastspf,inholeflag
  DATATYPE,intent(in) :: inspfs1(spfsize,nspf), inspfs2(spfsize,nspf)
  DATATYPE,intent(out) :: twoereduced(reducedpotsize,nspf,firstspf:lastspf)
  integer :: times(*), i,j,ispf

  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "     ...sparseops_matel"; CFL; call mpibarrier()
  endif
  call myclock(i)
  call sparseops_matel(matrix_ptr,inspfs1,inspfs2)
  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "     ...pot_matel"; CFL; call mpibarrier()
  endif
  call myclock(j);  times(1)=times(1)+j-i; i=j
  call pot_matel(matrix_ptr,inspfs1,inspfs2)
  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "     ...pulse_matel"; CFL; call mpibarrier()
  endif
  call myclock(j);  times(2)=times(2)+j-i; i=j
  call pulse_matel(matrix_ptr,inspfs1,inspfs2)
  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "     ...twoe_matel"; CFL; call mpibarrier()
  endif
  call myclock(j);  times(3)=times(3)+j-i; i=j
  call twoe_matel(matrix_ptr,inspfs1,inspfs2,twoereduced,firstspf,lastspf)
  call myclock(j);  times(4)=times(4)+j-i
  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "     ...done all_matel0"; CFL; call mpibarrier()
  endif

!!  OFLWR "CHECKREDUCEDxxx"
!!  write(mpifileptr,'(2F10.5)') twoereduced(1,:,1)
!!  OFLWR "REDSTOPxxx"; CFLST

  matrix_ptr%xtwoe_htrace(:,:) = 0d0     

  if (inholeflag.ne.0) then
     do ispf=1,nspf
        matrix_ptr%xtwoe_htrace(:,:) =  matrix_ptr%xtwoe_htrace(:,:) + &
             ( -2*matrix_ptr%xtwoematel(:,:,ispf,ispf) &
             + matrix_ptr%xtwoematel(:,ispf,ispf,:) )
     enddo
     matrix_ptr%xtwoe_htrace(:,:) =  TRANSPOSE(matrix_ptr%xtwoe_htrace(:,:))
  endif

  if (debugflag.eq.4242) then
     OFLWR "zeroing two-electron"; CFL
     twoereduced(:,:,:)=0d0; matrix_ptr%xtwoematel(:,:,:,:)=0d0
  endif

!!  OFLWR "HTRACEMAT", matrix_ptr%xtwoe_htrace(1,1); CFL

end subroutine all_matel0

subroutine twoe_matel(matrix_ptr,inspfs1,inspfs2,twoereduced,firstspf,lastspf)
  use parameters
  use mpimod
  use configptrmod
  use orbgathersubmod
  use mpisubmod
  implicit none
  integer,intent(in) :: firstspf,lastspf
  DATATYPE,intent(in) :: inspfs1(spfsize,nspf),inspfs2(spfsize,nspf)
  Type(CONFIGPTR),intent(inout) :: matrix_ptr
  DATATYPE,intent(out) :: twoereduced(reducedpotsize,nspf,firstspf:lastspf)
  integer :: lowspf,highspf,numspf

  lowspf=1; highspf=nspf
  if (parorbsplit.eq.1) then
     call getorbsetrange(lowspf,highspf)
  endif
  numspf=highspf-lowspf+1

  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "         In twoe_matel.  Calling call_twoe_matel"; CFL; call mpibarrier()
  endif
  if (numspf.gt.0) then
     call call_twoe_matel00(lowspf,highspf,inspfs1,inspfs2,matrix_ptr%xtwoematel(:,:,:,lowspf:highspf),&
          twoereduced(:,:,lowspf:highspf),timingdir,notiming)
  endif
  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "         Done call_twoe_matel.  Reduce"; CFL; call mpibarrier()
  endif
  if (parorbsplit.eq.1) then
     call mpiorbgather(matrix_ptr%xtwoematel(:,:,:,:),nspf**3)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(matrix_ptr%xtwoematel(:,:,:,:),nspf**4)
  endif
  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "         Done reduce, done twoe_matel."; CFL; call mpibarrier()
  endif

end subroutine twoe_matel

subroutine pot_matel(matrix_ptr,inspfs1,inspfs2)
  use parameters
  use configptrmod
  use orbgathersubmod
  use mpisubmod
  use orbmultsubmod
  use opmod   !! frozenreduced, hatomreduced
  implicit none
  DATATYPE,intent(in) :: inspfs1(spfsize,nspf), inspfs2(spfsize,nspf)
  Type(CONFIGPTR),intent(inout) :: matrix_ptr
  DATATYPE,allocatable :: ttempspfs(:,:),workspfs(:,:)
  integer :: lowspf,highspf,numspf

  matrix_ptr%xpotmatel(:,:)=0d0

  lowspf=1; highspf=nspf
  if (parorbsplit.eq.1) then
     call getorbsetrange(lowspf,highspf)
  endif
  numspf=highspf-lowspf+1

  if (numspf.gt.0) then

     allocate(ttempspfs(spfsize,lowspf:highspf),workspfs(spfsize,lowspf:highspf))
     ttempspfs=0; workspfs=0

     call mult_pot(numspf,inspfs2(:,lowspf:highspf),workspfs(:,lowspf:highspf))

!! usually just adding zeroes here
     call hatom_op(numspf,inspfs2(:,lowspf:highspf),ttempspfs(:,lowspf:highspf),hatomreduced(:))
     workspfs(:,lowspf:highspf)=workspfs(:,lowspf:highspf)+ttempspfs(:,lowspf:highspf)

!! potmatel is proper ordering.  fast index is conjg.

     if (numfrozen.gt.0) then
        call op_frozen_exchange(lowspf,highspf,inspfs2(:,lowspf:highspf),ttempspfs(:,lowspf:highspf))
        workspfs(:,lowspf:highspf)=workspfs(:,lowspf:highspf)+ttempspfs(:,lowspf:highspf)

        call op_frozenreduced(numspf,inspfs2(:,lowspf:highspf),ttempspfs(:,lowspf:highspf),frozenreduced(:))
        workspfs(:,lowspf:highspf)=workspfs(:,lowspf:highspf)+ttempspfs(:,lowspf:highspf)
     endif

     call MYGEMM(CNORMCHAR,'N',nspf,numspf, spfsize, DATAONE, inspfs1, spfsize, workspfs(:,lowspf:highspf), &
          spfsize, DATAZERO, matrix_ptr%xpotmatel(:,lowspf:highspf), nspf)

     deallocate(ttempspfs,workspfs)

  endif

  if (parorbsplit.eq.1) then
     call mpiorbgather(matrix_ptr%xpotmatel(:,:),nspf)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(matrix_ptr%xpotmatel(:,:),nspf**2)
  endif

end subroutine pot_matel



subroutine pulse_matel(matrix_ptr,inspfs1,inspfs2)
  use parameters
  use configptrmod
  use orbgathersubmod
  use mpisubmod
  use orbmultsubmod
  implicit none
  DATATYPE,intent(in) :: inspfs1(spfsize,nspf), inspfs2(spfsize,nspf)
  Type(CONFIGPTR),intent(inout) :: matrix_ptr
  DATATYPE :: nullcomplex(1), dipoles(3)
  DATATYPE,allocatable :: ttempspfsxx(:,:), ttempspfsyy(:,:), ttempspfszz(:,:)
  integer :: lowspf,highspf,numspf

  matrix_ptr%xpulsematelxx=0.d0;     matrix_ptr%xpulsematelyy=0.d0;     matrix_ptr%xpulsematelzz=0.d0

  lowspf=1; highspf=nspf
  if (parorbsplit.eq.1) then
     call getorbsetrange(lowspf,highspf)
  endif
  numspf=highspf-lowspf+1

  dipoles(:)=0d0
  if (velflag.eq.0) then
     call nucdipvalue(nullcomplex,dipoles)
  endif
  matrix_ptr%xpulsenuc(:)=dipoles(:)

  if (numspf.gt.0) then

     allocate(ttempspfsxx(spfsize,lowspf:highspf), ttempspfsyy(spfsize,lowspf:highspf), ttempspfszz(spfsize,lowspf:highspf))
     ttempspfsxx=0; ttempspfsyy=0; ttempspfszz=0

     if (velflag.eq.0) then
        call lenmultiply(numspf,inspfs2(:,lowspf:highspf),ttempspfsxx(:,lowspf:highspf), DATAONE,DATAZERO,DATAZERO)
        call lenmultiply(numspf,inspfs2(:,lowspf:highspf),ttempspfsyy(:,lowspf:highspf), DATAZERO,DATAONE,DATAZERO)
        call lenmultiply(numspf,inspfs2(:,lowspf:highspf),ttempspfszz(:,lowspf:highspf), DATAZERO,DATAZERO,DATAONE)
     else
        call velmultiply(numspf,inspfs2(:,lowspf:highspf),ttempspfsxx(:,lowspf:highspf), DATAONE,DATAZERO,DATAZERO)
        call velmultiply(numspf,inspfs2(:,lowspf:highspf),ttempspfsyy(:,lowspf:highspf), DATAZERO,DATAONE,DATAZERO)
        call velmultiply(numspf,inspfs2(:,lowspf:highspf),ttempspfszz(:,lowspf:highspf), DATAZERO,DATAZERO,DATAONE)
     endif
!! pulsematel is proper order.

     call MYGEMM(CNORMCHAR,'N',nspf,numspf, spfsize, DATAONE, inspfs1, spfsize, ttempspfsxx(:,lowspf:highspf), &
          spfsize, DATAONE, matrix_ptr%xpulsematelxx(:,lowspf:highspf) ,nspf)
     call MYGEMM(CNORMCHAR,'N',nspf,numspf, spfsize, DATAONE, inspfs1, spfsize, ttempspfsyy(:,lowspf:highspf), &
          spfsize, DATAONE, matrix_ptr%xpulsematelyy(:,lowspf:highspf) ,nspf)
     call MYGEMM(CNORMCHAR,'N',nspf,numspf, spfsize, DATAONE, inspfs1, spfsize, ttempspfszz(:,lowspf:highspf), &
          spfsize, DATAONE, matrix_ptr%xpulsematelzz(:,lowspf:highspf) ,nspf)

     deallocate(ttempspfsxx, ttempspfsyy, ttempspfszz)

  endif

  if (parorbsplit.eq.1) then
     call mpiorbgather(matrix_ptr%xpulsematelxx(:,:),nspf)
     call mpiorbgather(matrix_ptr%xpulsematelyy(:,:),nspf)
     call mpiorbgather(matrix_ptr%xpulsematelzz(:,:),nspf)
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(matrix_ptr%xpulsematelxx(:,:),nspf**2)
     call mympireduce(matrix_ptr%xpulsematelyy(:,:),nspf**2)
     call mympireduce(matrix_ptr%xpulsematelzz(:,:),nspf**2)
  endif

end subroutine pulse_matel


subroutine sparseops_matel(matrix_ptr,inspfs1,inspfs2)
  use parameters
  use configptrmod
  use orbgathersubmod
  use mpisubmod
  implicit none
  DATATYPE,intent(in) :: inspfs1(spfsize,nspf), inspfs2(spfsize,nspf)
  Type(CONFIGPTR),intent(inout) :: matrix_ptr
  DATATYPE,allocatable :: ttempspfs(:,:)
  integer :: lowspf,highspf,numspf

  matrix_ptr%xopmatel=0d0;  matrix_ptr%xymatel=0d0;  

  lowspf=1; highspf=nspf
  if (parorbsplit.eq.1) then
     call getorbsetrange(lowspf,highspf)
  endif
  numspf=highspf-lowspf+1

  if (numspf.gt.0) then

     allocate(ttempspfs(spfsize,lowspf:highspf)); ttempspfs=0

     call mult_ke(inspfs2(:,lowspf:highspf),ttempspfs(:,lowspf:highspf),numspf,timingdir,notiming)

     call MYGEMM(CNORMCHAR,'N',nspf,numspf, spfsize, DATAONE, inspfs1, spfsize, ttempspfs(:,lowspf:highspf), &
          spfsize, DATAZERO, matrix_ptr%xopmatel(:,lowspf:highspf) ,nspf)

     if (nonuc_checkflag/=1) then
        call op_yderiv(numspf,inspfs2(:,lowspf:highspf),ttempspfs(:,lowspf:highspf))
     
        call MYGEMM(CNORMCHAR,'N',nspf,numspf, spfsize, DATAONE, inspfs1, spfsize, ttempspfs(:,lowspf:highspf), &
             spfsize, DATAZERO, matrix_ptr%xymatel(:,lowspf:highspf) ,nspf)
     endif

     deallocate(ttempspfs)

  endif

  if (parorbsplit.eq.1) then
     call mpiorbgather(matrix_ptr%xopmatel(:,:),nspf)
     if (nonuc_checkflag/=1) then
        call mpiorbgather(matrix_ptr%xymatel(:,:),nspf)
     endif
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(matrix_ptr%xopmatel(:,:),nspf**2)
     if (nonuc_checkflag/=1) then
        call mympireduce(matrix_ptr%xymatel(:,:),nspf**2)
     endif
  endif

end subroutine sparseops_matel

end module allmat0mod



module cfgsubmod
contains

subroutine arbitraryconfig_matel_singles00transpose_hhh(www,in_onebodymat, smallmatrixtr,inholeflag)
  use fileptrmod
  use sparse_parameters
  use walkmod
  use mpisubmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: inholeflag
  DATATYPE,intent(in) :: in_onebodymat(www%nspf,www%nspf)
  DATATYPE,intent(out) :: smallmatrixtr(www%singlematsize,www%configstart:www%configend)
  DATATYPE :: onebodymat(www%nspf,www%nspf)        !! AUTOMATIC
  DATATYPE :: csum, myvec(www%singlematsize),holetrace
  integer ::    config1,  iwalk,myind,ihop,ispf

  if (www%configend.ge.www%configstart) then
     smallmatrixtr(:,:)=0d0
  endif

  onebodymat(:,:)=in_onebodymat(:,:)

!! HOLEFLAG

  if (www%holeflag.ne.0.and.inholeflag.eq.1) then
     holetrace=0
     do ispf=1,www%nspf
        holetrace=holetrace + onebodymat(ispf,ispf) * 2
     enddo

     onebodymat(:,:) = (-1) * in_onebodymat(:,:)

     do config1=www%botconfig,www%topconfig
        if (www%singlehopdiagflag(config1).ne.0) then
           if (sparseconfigflag.eq.0) then
              myind=config1
           else
              myind=www%singlediaghop(config1)
           endif
           smallmatrixtr(myind,config1) = holetrace
        endif
     enddo

  endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(config1,ihop,iwalk,myind,csum,myvec)
  
!$OMP DO SCHEDULE(DYNAMIC)
  do config1=www%botconfig,www%topconfig
     myvec(:)=smallmatrixtr(:,config1)
     do ihop=1,www%numsinglehops(config1)
        if (sparseconfigflag.eq.0) then
           myind=www%singlehop(ihop,config1)
        else
           myind=ihop
        endif

        csum=0d0
        do iwalk=www%singlehopwalkstart(ihop,config1),www%singlehopwalkend(ihop,config1)
           csum=csum+&
                onebodymat(www%singlewalkopspf(1,iwalk+www%scol(config1)), &
                www%singlewalkopspf(2,iwalk+www%scol(config1))) *  &
                www%singlewalkdirphase(iwalk+www%scol(config1))
        enddo
        myvec(myind)=myvec(myind)+csum
     enddo
     smallmatrixtr(:,config1)=myvec(:)
  enddo
!$OMP END DO
!$OMP END PARALLEL
  
  if (sparseconfigflag.eq.0) then
     call mpiallgather(smallmatrixtr(:,:),www%numconfig**2,www%configsperproc(:)*www%numconfig,&
          www%maxconfigsperproc*www%numconfig)
  endif
    
end subroutine arbitraryconfig_matel_singles00transpose_hhh
  

subroutine arbitraryconfig_matel_singles00transpose(www,onebodymat, smallmatrixtr)
  use fileptrmod
  use sparse_parameters
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: onebodymat(www%nspf,www%nspf)
  DATATYPE,intent(out) :: smallmatrixtr(www%singlematsize,www%configstart:www%configend)

  call arbitraryconfig_matel_singles00transpose_hhh(www,onebodymat, smallmatrixtr,1)

end subroutine arbitraryconfig_matel_singles00transpose


subroutine arbitraryconfig_matel_doubles00transpose(www,in_twobodymat, smallmatrixtr)
  use sparse_parameters
  use walkmod
  use mpisubmod
  implicit none
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: in_twobodymat(www%nspf,www%nspf,www%nspf,www%nspf)
  DATATYPE,intent(out) :: smallmatrixtr(www%doublematsize,www%configstart:www%configend)
  DATATYPE,allocatable :: twobodymat(:,:,:,:)
  DATATYPE :: csum, myvec(www%doublematsize),holetrace
  integer ::    config1, iwalk,myind,ihop,ispf,jspf

  if (www%configend.ge.www%configstart) then
     smallmatrixtr=0d0;
  endif

  allocate(twobodymat(www%nspf,www%nspf,www%nspf,www%nspf))
  twobodymat(:,:,:,:) = in_twobodymat(:,:,:,:)

!! HOLEFLAG
  if (www%holeflag.ne.0) then

     holetrace=0
     do ispf=1,www%nspf
        do jspf=1,www%nspf
           holetrace=holetrace + twobodymat(ispf,ispf,jspf,jspf) * 2 &
                - twobodymat(ispf,jspf,jspf,ispf)
        enddo
     enddo

!!     OFLWR "HTRACE001", holetrace,www%nspf; CFL

!! No!  Two-electron operator      twobodymat(:,:,:,:) = (-1) * in_twobodymat(:,:,:,:)

     do config1=www%botconfig,www%topconfig
        if (www%doublehopdiagflag(config1).ne.0) then
           if (sparseconfigflag.eq.0) then
              myind=config1
           else
              myind=www%doublediaghop(config1)
           endif
           smallmatrixtr(myind,config1) = holetrace
        endif
     enddo
  endif


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(config1,ihop,iwalk,myind,csum,myvec) 

!$OMP DO SCHEDULE(DYNAMIC)
  do config1=www%botconfig,www%topconfig

     myvec(:)=smallmatrixtr(:,config1)

     do ihop=1,www%numdoublehops(config1)
        if (sparseconfigflag.eq.0) then
           myind=www%doublehop(ihop,config1)
        else
           myind=ihop
        endif

        csum=0d0
        do iwalk=www%doublehopwalkstart(ihop,config1),www%doublehopwalkend(ihop,config1)
           csum=csum+&
                twobodymat(www%doublewalkdirspf(1,iwalk+www%dcol(config1)), &
                www%doublewalkdirspf(2,iwalk+www%dcol(config1)),   &
                www%doublewalkdirspf(3,iwalk+www%dcol(config1)),   &
                www%doublewalkdirspf(4,iwalk+www%dcol(config1)))* &
                www%doublewalkdirphase(iwalk+www%dcol(config1))
        enddo
        myvec(myind)=myvec(myind)+csum
     enddo
     smallmatrixtr(:,config1)=myvec(:)
  enddo
!$OMP END DO
!$OMP END PARALLEL

  if (sparseconfigflag.eq.0) then
     call mpiallgather(smallmatrixtr(:,:),www%numconfig**2,www%configsperproc(:)*www%numconfig,&
          www%maxconfigsperproc*www%numconfig)
  endif

end subroutine arbitraryconfig_matel_doubles00transpose

end module cfgsubmod


module asssubmod
contains

subroutine assemble_dfbasismat(www,outmatrix, matrix_ptr, boflag, nucflag, pulseflag, conflag, time,imc)
  use sparse_parameters
  use r_parameters
  use fileptrmod
  use configptrmod
  use walkmod
  use basissubmod
  implicit none
  type(walktype),intent(in) :: www
  Type(configptr),intent(in) :: matrix_ptr
  integer,intent(in) :: nucflag, pulseflag, conflag, boflag,imc
  DATATYPE,intent(out) ::        outmatrix(www%numdfbasis*numr,www%numdfbasis*numr)
  real*8,intent(in) :: time
  DATATYPE, allocatable ::        bigmatrix(:,:),halfmatrix(:,:),halfmatrix2(:,:),donematrix(:,:)
!  integer :: iconfig

  if (sparseconfigflag/=0) then
     OFLWR "Error, assemble_full_config_matel called when sparseconfigflag/=0";CFLST
  endif
  
  allocate(  bigmatrix(www%numconfig*numr,www%numconfig*numr), halfmatrix(www%numdfbasis*numr,www%numconfig*numr),&
       halfmatrix2(www%numconfig*numr,www%numdfbasis*numr), donematrix(www%numdfbasis*numr,www%numdfbasis*numr))
  bigmatrix(:,:)=0d0; halfmatrix(:,:)=0d0; halfmatrix2(:,:)=0d0; donematrix(:,:)=0d0; outmatrix(:,:)=0d0

  call assemble_configmat(www,bigmatrix,matrix_ptr,boflag,nucflag,pulseflag,conflag,time,imc)

!! obviously not good to transform for each r; should transform earlier. 

  call basis_transformto_all(www,www%numconfig*numr*numr,bigmatrix,halfmatrix2)
  halfmatrix=TRANSPOSE(halfmatrix2)
  call basis_transformto_all(www,www%numdfbasis*numr*numr,halfmatrix,donematrix)
  outmatrix=transpose(donematrix)

  deallocate(bigmatrix,halfmatrix,halfmatrix2,donematrix)

!!  OFLWR "OUTMATRIX"
!!  write(mpifileptr,'(4F20.10)')  (outmatrix(iconfig,1),outmatrix(iconfig,www%numdfbasis*numr),&
!!       iconfig=1,www%numdfbasis*numr)
!!  CFLST

end subroutine assemble_dfbasismat



!! FOR NONSPARSE (HERE, ASSEMBLE_CONFIGMAT) INCLUDE A-SQUARED TERM.

subroutine assemble_configmat(www,bigconfigmat,matrix_ptr, boflag, nucflag, pulseflag, conflag,time,imc)
  use fileptrmod
  use sparse_parameters
  use ham_parameters
  use r_parameters
  use configptrmod
  use walkmod
  use opmod   !! rkemod, proderivmod, frozenkediag, frozenpotdiag, bondpoints, bondweights
  use cfgsubmod
  use pulsesubmod
  implicit none
  type(walktype),intent(in) :: www
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  integer,intent(in) :: conflag,boflag,nucflag,pulseflag,imc
  real*8,intent(in) :: time
  DATATYPE,intent(out) :: bigconfigmat(numr,www%numconfig,numr,www%numconfig)
  DATATYPE :: tempmatel(www%nspf,www%nspf)   !! AUTOMATIC
  DATATYPE, allocatable :: tempconfigmat(:,:),tempconfigmat2(:,:),diagmat(:,:,:)
  real*8 :: gg
  DATATYPE :: facs(3), csum0,csum
  integer :: ir,jr,i

  if (sparseconfigflag.ne.0) then
     OFLWR "BADDDCALLL"; CFLST
  endif

  allocate(tempconfigmat(www%numconfig,www%numconfig),tempconfigmat2(www%numconfig,www%numconfig), &
       diagmat(www%numconfig,www%numconfig,numr))

  diagmat(:,:,:)=0d0

  call vectdpot(time,velflag,facs,imc)

  if (boflag==1) then
     do ir=1,numr
        do i=1,www%numconfig
         diagmat(i,i,ir)=( frozenkediag/bondpoints(ir)**2 + &
              (nucrepulsion+frozenpotdiag)/bondpoints(ir) + energyshift ) * matrix_ptr%constfac 
        enddo
     enddo

     call arbitraryconfig_matel_doubles00transpose(www,matrix_ptr%xtwoematel(:,:,:,:),tempconfigmat)

     if (www%holeflag.ne.0) then
        call arbitraryconfig_matel_singles00transpose_hhh(www,matrix_ptr%xtwoe_htrace(:,:), tempconfigmat2, 0)
        tempconfigmat(:,:)=tempconfigmat(:,:) + tempconfigmat2(:,:)
     endif

     call arbitraryconfig_matel_singles00transpose(www,matrix_ptr%xpotmatel(:,:), tempconfigmat2)

     tempconfigmat(:,:) = TRANSPOSE( tempconfigmat(:,:) + tempconfigmat2(:,:) )

     do ir=1,numr
        diagmat(:,:,ir)=diagmat(:,:,ir)+tempconfigmat(:,:)/bondpoints(ir)
     enddo
     call arbitraryconfig_matel_singles00transpose(www,matrix_ptr%xopmatel(:,:),  tempconfigmat)
     tempconfigmat(:,:)=TRANSPOSE(tempconfigmat(:,:))
     do ir=1,numr
        diagmat(:,:,ir)=diagmat(:,:,ir)+tempconfigmat(:,:)/bondpoints(ir)**2
     enddo
  endif

  if (conflag==1.and.constraintflag.ne.0) then
     tempmatel(:,:)=matrix_ptr%xconmatel(:,:)
     if (tdflag.ne.0) then
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelxx(:,:)*facs(1)
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelyy(:,:)*facs(2)
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelzz(:,:)*facs(3)
     endif
     call arbitraryconfig_matel_singles00transpose(www,tempmatel,tempconfigmat)
     tempconfigmat(:,:)=TRANSPOSE(tempconfigmat(:,:))
     do ir=1,numr
        diagmat(:,:,ir)=diagmat(:,:,ir)-tempconfigmat(:,:)
     enddo
  endif

  if ((pulseflag.eq.1.and.tdflag.ne.0)) then
     tempmatel(:,:)= matrix_ptr%xpulsematelxx*facs(1)+ &
          matrix_ptr%xpulsematelyy*facs(2)+ &
          matrix_ptr%xpulsematelzz*facs(3)
     call arbitraryconfig_matel_singles00transpose(www,tempmatel,tempconfigmat)
     tempconfigmat(:,:)=TRANSPOSE(tempconfigmat(:,:))
     do ir=1,numr
        if (velflag.eq.0) then
           diagmat(:,:,ir)=diagmat(:,:,ir)+tempconfigmat(:,:)*bondpoints(ir)
        else
           diagmat(:,:,ir)=diagmat(:,:,ir)+tempconfigmat(:,:)/bondpoints(ir)
        endif
     enddo

!! A-SQUARED TERM...  looked ok for z polarization.  number of electrons times a-squared.  fac of 2 needed apparently.
!!   for length... need nuclear dipole too.

!!   for velocity ...  also need derivative operator in bond length for hetero !!!!  TO DO !!!!


        gg=0.25d0


     if (velflag.ne.0) then 
        csum0 = matrix_ptr%constfac * www%numelec * (facs(1)**2 + facs(2)**2 + facs(3)**2) * 2  * gg   !! a-squared term.
     else
!! xpulsenuc is nuclear dipole, divided by R.    from nucdipvalue.

        csum0 = (matrix_ptr%xpulsenuc(1) * facs(1) + matrix_ptr%xpulsenuc(2) * facs(2) + matrix_ptr%xpulsenuc(3) * facs(3) )
     endif
     
     do ir=1,numr
        if (velflag.eq.0) then
           csum=bondpoints(ir) *csum0
        else
           csum=   csum0    !! A-squared !!    no R factor !!
        endif
        do i=1,www%numconfig
           diagmat(i,i,ir)=diagmat(i,i,ir) + csum 
        enddo
     enddo
  endif   !! pulse

  bigconfigmat(:,:,:,:)=0d0

  do ir=1,numr
     bigconfigmat(ir,:,ir,:)=diagmat(:,:,ir)
  enddo

  if (nonuc_checkflag.eq.0.and.nucflag.eq.1) then
       call arbitraryconfig_matel_singles00transpose(www,matrix_ptr%xymatel, tempconfigmat)
       tempconfigmat(:,:)=TRANSPOSE(tempconfigmat(:,:))

     do ir=1,numr
        do jr=1,numr
           bigconfigmat(ir,:,jr,:)=bigconfigmat(ir,:,jr,:)+tempconfigmat(:,:)*proderivmod(ir,jr)
        enddo
     enddo
     do i=1,www%numconfig
        bigconfigmat(:,i,:,i)=bigconfigmat(:,i,:,i) + rkemod(:,:)*matrix_ptr%kefac
     enddo
  endif

  deallocate(tempconfigmat,tempconfigmat2,diagmat)

end subroutine assemble_configmat

!! NO A-SQUARED TERM HERE; TIME IS NOT AN ARGUMENT.  A-squared in sparsesparsemult_byproc.

subroutine assemble_sparsemats(www,matrix_ptr, sparse_ptr,boflag, nucflag, pulseflag, conflag)
  use fileptrmod
  use r_parameters
  use sparse_parameters
  use ham_parameters
  use configptrmod
  use walkmod
  use sparseptrmod
  use opmod   !! rkemod, proderivmod
  use cfgsubmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: conflag,boflag,nucflag,pulseflag
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(inout) :: sparse_ptr

  if (sparseconfigflag.eq.0) then
     OFLWR "BADDDCAL555LL"; CFLST
  endif

  sparse_ptr%kefac = matrix_ptr%kefac
  sparse_ptr%constfac = matrix_ptr%constfac
  sparse_ptr%xpulsenuc(:) = matrix_ptr%xpulsenuc(:)

  if (boflag==1) then

     call arbitraryconfig_matel_doubles00transpose(www,matrix_ptr%xtwoematel(:,:,:,:),sparse_ptr%xpotsparsemattr(:,:))

     sparse_ptr%xpottracemattr(:,:)=0d0
     if (www%holeflag.ne.0) then
        call arbitraryconfig_matel_singles00transpose_hhh(www,matrix_ptr%xtwoe_htrace(:,:), sparse_ptr%xpottracemattr(:,:), 0)
     endif

     call arbitraryconfig_matel_singles00transpose(www,matrix_ptr%xpotmatel(:,:),sparse_ptr%xonepotsparsemattr(:,:))
     call arbitraryconfig_matel_singles00transpose(www,matrix_ptr%xopmatel(:,:),sparse_ptr%xopsparsemattr(:,:))
  endif

  if (conflag==1.and.constraintflag.ne.0) then
     call arbitraryconfig_matel_singles00transpose(www,matrix_ptr%xconmatel,  sparse_ptr%xconsparsemattr)
     if (tdflag.ne.0) then
        call arbitraryconfig_matel_singles00transpose(www,matrix_ptr%xconmatelxx,sparse_ptr%xconsparsemattrxx)
        call arbitraryconfig_matel_singles00transpose(www,matrix_ptr%xconmatelyy,sparse_ptr%xconsparsemattryy)
        call arbitraryconfig_matel_singles00transpose(www,matrix_ptr%xconmatelzz,sparse_ptr%xconsparsemattrzz)
     endif
  endif

  if ((pulseflag.eq.1.and.tdflag.ne.0)) then
     call arbitraryconfig_matel_singles00transpose(www,matrix_ptr%xpulsematelxx, sparse_ptr%xpulsesparsemattrxx)
     call arbitraryconfig_matel_singles00transpose(www,matrix_ptr%xpulsematelyy, sparse_ptr%xpulsesparsemattryy)
     call arbitraryconfig_matel_singles00transpose(www,matrix_ptr%xpulsematelzz, sparse_ptr%xpulsesparsemattrzz)
  endif
  if (nonuc_checkflag.eq.0.and.nucflag.eq.1) then
       call arbitraryconfig_matel_singles00transpose(www,matrix_ptr%xymatel, sparse_ptr%xysparsemattr)
  endif

end subroutine assemble_sparsemats

end module asssubmod


subroutine all_matel()
  use parameters
  use xxxmod
  use mpimod
  use configmod
  use opmod
  use mpi_orbsetmod
  use asssubmod
  use cfgsubmod
  use allmat0mod
  implicit none
  integer,save :: times(10)=0,xcalled=0
  integer :: itime,jtime,getlen,firstspf,lastspf
  xcalled=xcalled+1

  firstspf=1;lastspf=nspf
  if (parorbsplit.eq.1) then
     firstspf=firstmpiorb
     lastspf=firstmpiorb+orbsperproc-1
  endif

  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "     Call all_matel0 in all_matel"; CFL; call mpibarrier()
  endif
  call all_matel0(www%holeflag,yyy%cptr(0), yyy%cmfspfs(:,0), yyy%cmfspfs(:,0), &
       twoereduced,times,firstspf,lastspf)
  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "     ...called all_matel0 in all_matel"; CFL; call mpibarrier()
  endif

  call myclock(itime)
  if (sparseopt.ne.0) then
     call assemble_sparsemats(www,yyy%cptr(0),yyysptr(0),1,1,1,1)
     if (use_dfwalktype) then
        if (shuffle_dfwalktype) then
           call assemble_sparsemats(fdww,yyy%cptr(0),yyysfdptr(0),1,1,1,1)
        else
           call assemble_sparsemats(dfww,yyy%cptr(0),yyysdfptr(0),1,1,1,1)
        endif
     endif
  endif
  call myclock(jtime); times(5)=times(5)+jtime-itime

  if ((myrank.eq.1).and.(notiming.eq.0)) then
     if (xcalled==1) then
        open(853, file=timingdir(1:getlen(timingdir))//"/matel.time.dat", status="unknown")
        write(853,'(100A11)')   "op", "pot", "pulse", "two", "assemble"
        close(853)
     endif

     open(853, file=timingdir(1:getlen(timingdir))//"/matel.time.dat", &
          status="unknown", position="append")
     write(853,'(100I11)')  times(1:5)/1000;        close(853)
     close(853)
  endif

  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "     ...finished all_matel"; CFL; call mpibarrier()
  endif

end subroutine all_matel


