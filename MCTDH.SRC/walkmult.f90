
#include "Definitions.INC"

!! All purpose subroutine.

recursive subroutine sparseconfigmult(www,invector,outvector,matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time)
  use r_parameters
  use configptrmod
  use sparseptrmod
  use walkmod
  use sparse_parameters  !! sparsesummaflag
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: conflag,boflag,nucflag,pulseflag
  DATATYPE,intent(in) :: invector(www%localnconfig*numr)
  DATATYPE,intent(out) :: outvector(www%localnconfig*numr)
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr
  real*8,intent(in) :: time

  if (sparsesummaflag.eq.0) then
     call sparseconfigmultxxx_gather(www,invector,outvector,matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,0)
  else
     call sparseconfigmultxxx_summa(www,invector,outvector,matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,0)
  endif

end subroutine sparseconfigmult


!! For one bond length, born oppenheimer Hamiltonian

recursive subroutine sparseconfigmultone(www,invector,outvector,matrix_ptr,sparse_ptr, boflag,pulseflag,conflag,isplit,time)
  use configptrmod
  use sparseptrmod
  use sparse_parameters  !! sparsesummaflag
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: conflag, boflag, pulseflag, isplit
  DATATYPE,intent(in) :: invector(www%firstconfig:www%lastconfig)
  DATATYPE,intent(out) :: outvector(www%firstconfig:www%lastconfig)
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr
  real*8,intent(in) :: time

  if (sparsesummaflag.eq.0) then
     call sparseconfigmultone_gather(www,invector,outvector,matrix_ptr,sparse_ptr, boflag,pulseflag,conflag,isplit,time)
  else
     call sparseconfigmultone_summa(www,invector,outvector,matrix_ptr,sparse_ptr, boflag,pulseflag,conflag,isplit,time)
  endif
end subroutine sparseconfigmultone


recursive subroutine sparseconfigmultone_gather(www,invector,outvector,matrix_ptr,sparse_ptr, boflag,pulseflag,conflag,isplit,time)
  use configptrmod
  use sparseptrmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: conflag, boflag, pulseflag, isplit
  DATATYPE,intent(in) :: invector(www%firstconfig:www%lastconfig)
  DATATYPE,intent(out) :: outvector(www%firstconfig:www%lastconfig)
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr
  real*8,intent(in) :: time
  DATATYPE :: workvector(www%numconfig)      !! AUTOMATIC

  workvector(www%firstconfig:www%lastconfig)=invector(:)
  if (www%parconsplit.ne.0) then
     call mpiallgather(workvector,www%numconfig,www%configsperproc(:),www%maxconfigsperproc)
  endif
  call sparseconfigmult_nompi(www,workvector,outvector(www%botconfig),matrix_ptr,sparse_ptr, boflag, 0, pulseflag, conflag,time,0,isplit,isplit,0)
  if (www%parconsplit.eq.0) then
     call mpiallgather(outvector,www%numconfig,www%configsperproc(:),www%maxconfigsperproc)
  endif

end subroutine sparseconfigmultone_gather


recursive subroutine sparseconfigmultone_summa(www,invector,outvector,matrix_ptr,sparse_ptr, boflag,pulseflag,conflag,isplit,time)
  use configptrmod
  use sparseptrmod
  use walkmod
  use mpimod    !! nprocs,myrank
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: conflag, boflag, pulseflag, isplit
  DATATYPE,intent(in) :: invector(www%firstconfig:www%lastconfig)
  DATATYPE,intent(out) :: outvector(www%firstconfig:www%lastconfig)
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr
  real*8,intent(in) :: time
  DATATYPE :: workvector(www%maxconfigsperproc),workvector2(www%botconfig:www%topconfig)      !! AUTOMATIC
  integer :: iproc

  outvector(www%botconfig:www%topconfig)=0d0

  do iproc=1,nprocs

     if (www%parconsplit.ne.0) then     

        if (myrank.eq.iproc) then
           workvector(1:www%topconfig-www%botconfig+1)=invector(:)
        endif
        call mympibcast(workvector,iproc,(www%alltopconfigs(iproc)-www%allbotconfigs(iproc)+1))
        call sparseconfigmult_byproc(iproc,www,workvector,workvector2,matrix_ptr,sparse_ptr, boflag, 0, pulseflag, conflag,time,0,isplit,isplit,0)

     else

        call sparseconfigmult_byproc(iproc,www,invector(www%allbotconfigs(iproc):www%alltopconfigs(iproc)),&
             workvector2,matrix_ptr,sparse_ptr, boflag, 0, pulseflag, conflag,time,0,isplit,isplit,0)

     endif

     outvector(www%botconfig:www%topconfig)=outvector(www%botconfig:www%topconfig)+workvector2(:)

  enddo

  if (www%parconsplit.eq.0) then
     call mpiallgather(outvector,www%numconfig,www%configsperproc(:),www%maxconfigsperproc)
  endif

end subroutine sparseconfigmultone_summa



!! All purpose subroutine, just the pulse.

recursive subroutine sparseconfigpulsemult(www,invector,outvector,matrix_ptr, sparse_ptr,which)
  use r_parameters
  use configptrmod
  use sparseptrmod
  use walkmod
  use sparse_parameters  !! sparsesummaflag
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: which
  DATATYPE,intent(in) :: invector(www%totadim)
  DATATYPE,intent(out) :: outvector(www%totadim)
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr

  if (sparsesummaflag.eq.0) then
     call sparseconfigmultxxx_gather(www,invector,outvector,matrix_ptr,sparse_ptr, 0,0,1,0,-1d0,which)
  else
     call sparseconfigmultxxx_summa(www,invector,outvector,matrix_ptr,sparse_ptr, 0,0,1,0,-1d0,which)
  endif

end subroutine sparseconfigpulsemult


!! MPI subroutine

recursive subroutine sparseconfigmultxxx_gather(www,invector,outvector,matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag)
  use r_parameters
  use configptrmod
  use sparseptrmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: conflag,boflag,nucflag,pulseflag,onlytdflag
  DATATYPE,intent(in) :: invector(numr,www%firstconfig:www%lastconfig)
  DATATYPE,intent(out) :: outvector(numr,www%firstconfig:www%lastconfig)
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr
  real*8,intent(in) :: time
  DATATYPE :: workvector(numr,www%numconfig)    !!AUTOMATIC

  workvector(:,www%firstconfig:www%lastconfig) = invector(:,:)
  if (www%parconsplit.ne.0) then
     call mpiallgather(workvector,www%numconfig*numr,www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
  endif

  call sparseconfigmult_nompi(www,workvector,outvector(:,www%botconfig),matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag,1,numr,0)

  if (www%parconsplit.eq.0) then
     call mpiallgather(outvector,www%numconfig*numr,www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
  endif
     
end subroutine sparseconfigmultxxx_gather



recursive subroutine sparseconfigmultxxx_summa(www,invector,outvector,matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag)
  use r_parameters
  use configptrmod
  use sparseptrmod
  use mpimod    !! myrank,nprocs
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: conflag,boflag,nucflag,pulseflag,onlytdflag
  DATATYPE,intent(in) :: invector(numr,www%firstconfig:www%lastconfig)
  DATATYPE,intent(out) :: outvector(numr,www%firstconfig:www%lastconfig)
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr
  real*8,intent(in) :: time
  DATATYPE :: workvector(numr,www%maxconfigsperproc), workvector2(numr,www%botconfig:www%topconfig)    !!AUTOMATIC
  integer :: iproc

  outvector(:,www%botconfig:www%topconfig)=0d0

  do iproc=1,nprocs

     if (www%parconsplit.ne.0) then

        if (myrank.eq.iproc) then
           workvector(:,1:www%topconfig-www%botconfig+1) = invector(:,:)
        endif
        call mympibcast(workvector,iproc,(www%alltopconfigs(iproc)-www%allbotconfigs(iproc)+1)*numr)
        call sparseconfigmult_byproc(iproc,www,workvector,workvector2,matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag,1,numr,0)

     else

        call sparseconfigmult_byproc(iproc,www,invector(:,www%allbotconfigs(iproc):www%alltopconfigs(iproc)),&
             workvector2,matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag,1,numr,0)

     endif

     outvector(:,www%botconfig:www%topconfig)=outvector(:,www%botconfig:www%topconfig) + workvector2(:,:)

  enddo

  if (www%parconsplit.eq.0) then
     call mpiallgather(outvector,www%numconfig*numr,www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
  endif
     
end subroutine sparseconfigmultxxx_summa





recursive subroutine parsparseconfigpreconmult(www,invector,outvector,matrix_ptr, sparse_ptr, boflag, nucflag, pulseflag, conflag, inenergy,time)
  use r_parameters
  use configptrmod
  use sparseptrmod
  use sparse_parameters  !! sparsesummaflag
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: invector(numr,www%botconfig:www%topconfig)
  DATATYPE,intent(out) :: outvector(numr,www%botconfig:www%topconfig)
  DATATYPE,intent(in) :: inenergy
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr
  integer,intent(in) ::  conflag,boflag,nucflag,pulseflag
  real*8,intent(in) :: time

  if (sparsesummaflag.eq.0) then
     call parsparseconfigpreconmult_gather(www,invector,outvector,matrix_ptr, sparse_ptr, boflag, nucflag, pulseflag, conflag, inenergy,time)
  else
     call parsparseconfigpreconmult_summa(www,invector,outvector,matrix_ptr, sparse_ptr, boflag, nucflag, pulseflag, conflag, inenergy,time)
  endif

end subroutine parsparseconfigpreconmult


recursive subroutine parsparseconfigpreconmult_gather(www,invector,outvector,matrix_ptr, sparse_ptr, boflag, nucflag, pulseflag, conflag, inenergy,time)
  use r_parameters
  use configptrmod
  use sparseptrmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: invector(numr,www%botconfig:www%topconfig)
  DATATYPE,intent(out) :: outvector(numr,www%botconfig:www%topconfig)
  DATATYPE,intent(in) :: inenergy
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr
  integer,intent(in) ::  conflag,boflag,nucflag,pulseflag
  real*8,intent(in) :: time
  DATATYPE :: workvector(numr,www%numconfig),tempvector(numr,www%botconfig:www%topconfig)     !!AUTOMATIC

  if (www%topconfig-www%botconfig+1.eq.0) then
     return
  endif

  workvector(:,:)=1d0; tempvector(:,:)=0d0
  call sparseconfigmult_nompi(www,workvector,tempvector,matrix_ptr, sparse_ptr, boflag, nucflag, pulseflag, conflag,time,0,1,numr,1)
  
  outvector(:,:)=invector(:,:)/(tempvector(:,:)-inenergy)

end subroutine parsparseconfigpreconmult_gather


recursive subroutine parsparseconfigpreconmult_summa(www,invector,outvector,matrix_ptr, sparse_ptr, boflag, nucflag, pulseflag, conflag, inenergy,time)
  use r_parameters
  use configptrmod
  use sparseptrmod
  use walkmod
  use mpimod   !! myrank
  implicit none
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: invector(numr,www%botconfig:www%topconfig)
  DATATYPE,intent(out) :: outvector(numr,www%botconfig:www%topconfig)
  DATATYPE,intent(in) :: inenergy
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr
  integer,intent(in) ::  conflag,boflag,nucflag,pulseflag
  real*8,intent(in) :: time
  DATATYPE :: workvector(numr,www%botconfig:www%topconfig),tempvector(numr,www%botconfig:www%topconfig)     !!AUTOMATIC

  if (www%topconfig-www%botconfig+1.eq.0) then
     return
  endif

  workvector(:,:)=1d0; tempvector(:,:)=0d0
  call sparseconfigmult_byproc(myrank,www,workvector,tempvector,matrix_ptr, sparse_ptr, boflag, nucflag, pulseflag, conflag,time,0,1,numr,1)
  
  outvector(:,:)=invector(:,:)/(tempvector(:,:)-inenergy)

end subroutine parsparseconfigpreconmult_summa



recursive subroutine arbitraryconfig_mult_singles(www,onebodymat, rvector, avectorin, avectorout,inrnum)   
  use walkmod
  use sparse_parameters  !! sparsesummaflag
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: inrnum
  DATATYPE,intent(in) :: onebodymat(www%nspf,www%nspf), avectorin(inrnum,www%firstconfig:www%lastconfig)
  DATATYPE,intent(out) :: avectorout(inrnum,www%firstconfig:www%lastconfig)
  DATAECS,intent(in) :: rvector(inrnum)

  if (sparsesummaflag.eq.0) then
     call arbitraryconfig_mult_singles_gather(www,onebodymat, rvector, avectorin, avectorout,inrnum)   
  else
     call arbitraryconfig_mult_singles_summa(www,onebodymat, rvector, avectorin, avectorout,inrnum)   
  endif

end subroutine arbitraryconfig_mult_singles


recursive subroutine arbitraryconfig_mult_singles_gather(www,onebodymat, rvector, avectorin, avectorout,inrnum)   
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: inrnum
  DATATYPE,intent(in) :: onebodymat(www%nspf,www%nspf), avectorin(inrnum,www%firstconfig:www%lastconfig)
  DATATYPE,intent(out) :: avectorout(inrnum,www%firstconfig:www%lastconfig)
  DATAECS,intent(in) :: rvector(inrnum)
  DATATYPE :: workvector(inrnum,www%numconfig)         !! AUTOMATIC

  workvector(:,:)=0d0
  workvector(:,www%firstconfig:www%lastconfig)=avectorin(:,:)

  if (www%parconsplit.ne.0) then
     call mpiallgather(workvector,www%numconfig*inrnum,www%configsperproc(:)*inrnum,www%maxconfigsperproc*inrnum)
  endif

  call arbitraryconfig_mult_singles_nompi(www,onebodymat, rvector, workvector, avectorout(:,www%botconfig),inrnum,0)

  if (www%parconsplit.eq.0) then
     call mpiallgather(avectorout,www%numconfig*inrnum,www%configsperproc(:)*inrnum,www%maxconfigsperproc*inrnum)
  endif

end subroutine arbitraryconfig_mult_singles_gather


recursive subroutine arbitraryconfig_mult_singles_summa(www,onebodymat, rvector, avectorin, avectorout,inrnum)   
  use walkmod
  use mpimod   !! myrank,nprocs
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: inrnum
  DATATYPE,intent(in) :: onebodymat(www%nspf,www%nspf), avectorin(inrnum,www%firstconfig:www%lastconfig)
  DATATYPE,intent(out) :: avectorout(inrnum,www%firstconfig:www%lastconfig)
  DATAECS,intent(in) :: rvector(inrnum)
  DATATYPE :: workvector(inrnum,www%maxconfigsperproc),workvector2(inrnum,www%botconfig:www%topconfig)         !! AUTOMATIC
  integer :: iproc

  avectorout(:,www%botconfig:www%topconfig)=0d0

  do iproc=1,nprocs

     if (www%parconsplit.ne.0) then

        if (myrank.eq.iproc) then
           workvector(:,1:www%topconfig-www%botconfig+1) = avectorin(:,:)
        endif
        call mympibcast(workvector,iproc,(www%alltopconfigs(iproc)-www%allbotconfigs(iproc)+1)*inrnum)
        call arbitraryconfig_mult_singles_byproc(iproc,www,onebodymat, rvector, workvector, workvector2,inrnum,0)

     else
     
        call arbitraryconfig_mult_singles_byproc(iproc,www,onebodymat, rvector, avectorin(:,www%allbotconfigs(iproc):www%alltopconfigs(iproc)),&
             workvector2,inrnum,0)

     endif

     avectorout(:,www%botconfig:www%topconfig) = avectorout(:,www%botconfig:www%topconfig) + workvector2(:,:)

  enddo

  if (www%parconsplit.eq.0) then
     call mpiallgather(avectorout,www%numconfig*inrnum,www%configsperproc(:)*inrnum,www%maxconfigsperproc*inrnum)
  endif

end subroutine arbitraryconfig_mult_singles_summa




!!!!!!!    NOMPI SUBROUTINES    !!!!!!!

recursive subroutine sparseconfigmult_nompi(www,invector,outvector,matrix_ptr,sparse_ptr,&
     boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag)
  use sparse_parameters
  use configptrmod
  use sparseptrmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: conflag,boflag,nucflag,pulseflag,onlytdflag,botr,topr,diagflag
  DATATYPE,intent(in) :: invector(botr:topr,www%numconfig)
  DATATYPE,intent(out) :: outvector(botr:topr,www%botconfig:www%topconfig)
  real*8,intent(in) :: time
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr

  if (www%topconfig-www%botconfig+1.eq.0) then
     return
  endif
  if (sparseopt.eq.0) then
     call direct_sparseconfigmult_nompi(www,invector,outvector,matrix_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag)
  else
     call sparsesparsemult_nompi(www,invector,outvector,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag)
  endif

end subroutine sparseconfigmult_nompi


!! DIRECT MATVEC (with matrix_ptr not sparse_ptr)

recursive subroutine direct_sparseconfigmult_nompi(www,invector,outvector,matrix_ptr, &
     boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag)
  use ham_parameters
  use r_parameters
  use configptrmod
  use opmod   !! rkemod, proderivmod  
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: conflag,boflag,nucflag,pulseflag,onlytdflag,botr,topr,diagflag
  DATATYPE,intent(in) :: invector(botr:topr,www%numconfig)
  DATATYPE,intent(out) :: outvector(botr:topr,www%botconfig:www%topconfig)
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  real*8, intent(in) :: time
  integer :: mynumr,ir,flag,ii
  DATATYPE :: tempvector(botr:topr,www%botconfig:www%topconfig) !! AUTOMATIC
  DATATYPE :: tempmatel(www%nspf,www%nspf),facs(3),csum
  DATAECS :: rvector(botr:topr)
  real*8 :: gg

  if (www%topconfig-www%botconfig+1.eq.0) then
     return
  endif

  outvector(:,:)=0d0

  if (onlytdflag.ne.0) then
     facs(:)=0;     facs(onlytdflag)=1d0
  else
     call vectdpot(time,velflag,facs)
  endif
  
  mynumr=topr-botr+1
  
  if (boflag==1) then
     
!! easy nuclear repulsion, hardwired like this for now

     do ir=botr,topr

        outvector(ir,:)=outvector(ir,:) + invector(ir,www%botconfig:www%topconfig) * ( &
             energyshift + (nucrepulsion+frozenpotdiag)/bondpoints(ir) + &
             frozenkediag/bondpoints(ir)**2 ) * matrix_ptr%kefac

     enddo
     
     rvector(:)=1/bondpoints(botr:topr)
     call arbitraryconfig_mult_doubles_nompi(www,matrix_ptr%xtwoematel(:,:,:,:),rvector,invector,tempvector,mynumr,diagflag)
     outvector(:,:)=outvector(:,:)+tempvector(:,:)
     
     call arbitraryconfig_mult_singles_nompi(www,matrix_ptr%xpotmatel(:,:),rvector,invector,tempvector,mynumr,diagflag)
     outvector(:,:)=outvector(:,:)+tempvector(:,:)
     
     rvector(:)=1/bondpoints(botr:topr)**2
     call arbitraryconfig_mult_singles_nompi(www,matrix_ptr%xopmatel(:,:),rvector,invector,tempvector,mynumr,diagflag)
     outvector(:,:)=outvector(:,:)+tempvector(:,:)
     
  endif
  
  if (conflag==1.and.constraintflag.ne.0) then
     rvector(:)=1
     
     tempmatel(:,:)=matrix_ptr%xconmatel(:,:)
     if (tdflag.ne.0) then
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelxx(:,:)*facs(1)
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelyy(:,:)*facs(2)
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelzz(:,:)*facs(3)
     endif
     call arbitraryconfig_mult_singles_nompi(www,tempmatel,rvector,invector,tempvector,mynumr,diagflag)
     outvector(:,:)=outvector(:,:)-tempvector(:,:)

  endif
  if (pulseflag.eq.1.and.tdflag.eq.1.or.onlytdflag.ne.0) then
     if (velflag.eq.0) then
        rvector(:)=bondpoints(botr:topr)
     else
        rvector(:)=1/bondpoints(botr:topr)
     endif

     tempmatel(:,:)=0d0; flag=0
     if (onlytdflag.eq.0.or.onlytdflag.eq.1) then
        tempmatel(:,:)= tempmatel(:,:) + matrix_ptr%xpulsematelxx*facs(1); flag=1
     endif
     if (onlytdflag.eq.0.or.onlytdflag.eq.2) then
        tempmatel(:,:)= tempmatel(:,:) + matrix_ptr%xpulsematelyy*facs(2); flag=1
     endif
     if (onlytdflag.eq.0.or.onlytdflag.eq.3) then
        tempmatel(:,:)= tempmatel(:,:) + matrix_ptr%xpulsematelzz*facs(3); flag=1
     endif
     if (flag.eq.1) then
        call arbitraryconfig_mult_singles_nompi(www,tempmatel,rvector,invector,tempvector,mynumr,diagflag)
        outvector(:,:)=outvector(:,:)+tempvector(:,:)
     endif

!! STILL NEED A-VECTOR TERM IN BOND LENGTH FOR HETERO!! TO DO!!

        gg=0.25d0

!!  A-squared term and nuclear dipole:  nuclear dipole always
     if (velflag.eq.0) then 
        csum=0
        do ii=1,3
           if (onlytdflag.eq.0.or.onlytdflag.eq.ii) then
              csum=csum+facs(ii)*matrix_ptr%xpulsenuc(ii)
           endif
        enddo
        do ir=botr,topr
           outvector(ir,:)=outvector(ir,:)+invector(ir,www%botconfig:www%topconfig)*csum*bondpoints(ir)
        enddo
     else if (onlytdflag.eq.0) then
        csum=matrix_ptr%kefac * www%numelec * (facs(1)**2 + facs(2)**2 + facs(3)**2 ) * 2  * gg   !! a-squared
        do ir=botr,topr
           outvector(ir,:)=outvector(ir,:)+invector(ir,www%botconfig:www%topconfig)*csum     !! NO R FACTOR !!
        enddo
     endif

  endif   !! PULSE

  if (nonuc_checkflag.eq.0.and.nucflag.eq.1.and.mynumr.eq.numr) then
     if (diagflag.eq.0) then
        rvector(:)=1d0
        call arbitraryconfig_mult_singles_nompi(www,matrix_ptr%xymatel,rvector,invector,tempvector,numr,0)
        call MYGEMM('N','N',numr,www%topconfig-www%botconfig+1,numr,DATAONE,            proderivmod,numr,     tempvector,           numr,DATAONE,outvector,numr)
        call MYGEMM('N','N',numr,www%topconfig-www%botconfig+1,numr,DATAONE*matrix_ptr%kefac,rkemod,numr,     invector(:,www%botconfig),numr,DATAONE,outvector,numr)
!!$     else
!!$        do ir=1,numr
!!$           outvector(ir,:)=outvector(ir,:)+invector(ir,www%botconfig:www%topconfig)*matrix_ptr%kefac*rkemod(ir,ir)
!!$        enddo
     endif
  endif

end subroutine direct_sparseconfigmult_nompi


!! SPARSE MATVEC (with sparse_ptr not matrix_ptr)

recursive subroutine sparsesparsemult_nompi(www,invector,outvector,sparse_ptr,&
     boflag,nucflag,pulseflag,conflag,time,onlytdflag,botr,topr,diagflag)
  use sparse_parameters
  use ham_parameters
  use r_parameters
  use sparseptrmod
  use walkmod
  use opmod   !! rkemod, proderivmod
  use fileptrmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: boflag,nucflag,pulseflag,conflag,onlytdflag,botr,topr,diagflag
  DATATYPE,intent(in)  :: invector(botr:topr,www%numconfig)
  DATATYPE,intent(out) :: outvector(botr:topr,www%botconfig:www%topconfig)
  Type(SPARSEPTR),intent(in) :: sparse_ptr
  real*8,intent(in) :: time
  DATATYPE :: facs(3),csum
  DATATYPE :: tempvector(botr:topr,www%botconfig:www%topconfig), &
       tempsparsemattr(www%singlematsize,www%botconfig:www%topconfig) !! AUTOMATIC
  real*8 :: gg
  DATAECS :: rvector(botr:topr)
  integer :: ir,mynumr,ii,flag

  if (sparseconfigflag.eq.0) then
     OFLWR "error no sparsesparsemult if sparseconfigflag.eq.0"; CFLST
  endif

  if (www%topconfig-www%botconfig+1.eq.0) then
     return
  endif

  if (botr.ne.1.or.topr.ne.numr) then
     OFLWR "REPROGRAMMMM."; CFLST
  endif

  outvector(:,:)=0d0

  if (onlytdflag.ne.0) then
     facs(:)=0;     facs(onlytdflag)=1d0
  else
     call vectdpot(time,velflag,facs)
  endif

  mynumr=topr-botr+1

  if (boflag==1) then

    do ir=botr,topr

       outvector(ir,:)=outvector(ir,:) + invector(ir,www%botconfig:www%topconfig) * ( &
            energyshift + (nucrepulsion+frozenpotdiag)/bondpoints(ir) + &
            frozenkediag/bondpoints(ir)**2 ) * sparse_ptr%kefac
       
    enddo

     rvector(:)=1/bondpoints(botr:topr)
     call arbitrary_sparsemult_doubles_nompi(www,sparse_ptr%xpotsparsemattr(:,:),rvector,invector,tempvector,mynumr,diagflag)
     outvector(:,:)=outvector(:,:)+tempvector(:,:)

     call arbitrary_sparsemult_singles_nompi(www,sparse_ptr%xonepotsparsemattr(:,:),rvector,invector,tempvector,mynumr,diagflag)
     outvector(:,:)=outvector(:,:)+tempvector(:,:)

     rvector(:)=1/bondpoints(botr:topr)**2
     call arbitrary_sparsemult_singles_nompi(www,sparse_ptr%xopsparsemattr(:,:),rvector,invector,tempvector,mynumr,diagflag)
     outvector(:,:)=outvector(:,:)+tempvector(:,:)

  endif

  if (conflag==1.and.constraintflag.ne.0) then
     rvector(:)=1

     tempsparsemattr(:,:)=sparse_ptr%xconsparsemattr(:,:)
     if (tdflag.ne.0) then
        tempsparsemattr(:,:)=tempsparsemattr(:,:)+sparse_ptr%xconsparsemattrxx(:,:)*facs(1)
        tempsparsemattr(:,:)=tempsparsemattr(:,:)+sparse_ptr%xconsparsemattryy(:,:)*facs(2)
        tempsparsemattr(:,:)=tempsparsemattr(:,:)+sparse_ptr%xconsparsemattrzz(:,:)*facs(3)
     endif
     call arbitrary_sparsemult_singles_nompi(www,tempsparsemattr,rvector,invector,tempvector,mynumr,diagflag)
     outvector(:,:)=outvector(:,:)-tempvector(:,:)

  endif
  if (pulseflag.eq.1.and.tdflag.eq.1.or.onlytdflag.ne.0) then
     if (velflag.eq.0) then
        rvector(:)=bondpoints(botr:topr)
     else
        rvector(:)=1/bondpoints(botr:topr)
     endif

     tempsparsemattr(:,:)=0d0; flag=0
     if (onlytdflag.eq.0.or.onlytdflag.eq.1) then
        tempsparsemattr(:,:)= tempsparsemattr(:,:) + sparse_ptr%xpulsesparsemattrxx*facs(1); flag=1
     endif
     if (onlytdflag.eq.0.or.onlytdflag.eq.2) then
        tempsparsemattr(:,:)= tempsparsemattr(:,:) + sparse_ptr%xpulsesparsemattryy*facs(2); flag=1
     endif
     if (onlytdflag.eq.0.or.onlytdflag.eq.3) then
        tempsparsemattr(:,:)= tempsparsemattr(:,:) + sparse_ptr%xpulsesparsemattrzz*facs(3); flag=1
     endif
     if (flag.eq.1) then
        call arbitrary_sparsemult_singles_nompi(www,tempsparsemattr,rvector,invector,tempvector,mynumr,diagflag)
        outvector(:,:)=outvector(:,:)+tempvector(:,:)
     endif

!! STILL NEED A-VECTOR TERM IN BOND LENGTH FOR HETERO !! TO DO !!

!! A-squared term and nuclear dipole : nuclear dipole always

     gg=0.25d0
     
     if (velflag.eq.0) then 
        csum=0
        do ii=1,3
           if (onlytdflag.eq.0.or.onlytdflag.eq.ii) then
              csum=csum+facs(ii)*sparse_ptr%xpulsenuc(ii)
           endif
        enddo
        do ir=botr,topr
           outvector(ir,:)=outvector(ir,:)+invector(ir,www%botconfig:www%topconfig)*csum*bondpoints(ir)
        enddo
     else if (onlytdflag.eq.0) then
        csum= sparse_ptr%kefac * www%numelec * (facs(1)**2 + facs(2)**2 + facs(3)**2 ) * 2 * gg  !! a-squared
        do ir=botr,topr
           outvector(ir,:)=outvector(ir,:)+invector(ir,www%botconfig:www%topconfig)*csum     !! NO R FACTOR !!
        enddo
     endif

  endif !! PULSE

  if (nonuc_checkflag.eq.0.and.nucflag.eq.1.and.mynumr.eq.numr) then
     if (diagflag.eq.0) then
        rvector(:)=1d0
        call arbitrary_sparsemult_singles_nompi(www,sparse_ptr%xysparsemattr,rvector,invector,tempvector,numr,0)
        call MYGEMM('N','N',numr,www%topconfig-www%botconfig+1,numr,DATAONE,            proderivmod,numr,     tempvector,           numr,DATAONE,outvector,numr)
        call MYGEMM('N','N',numr,www%topconfig-www%botconfig+1,numr,DATAONE*sparse_ptr%kefac,rkemod,numr,     invector(:,www%botconfig),numr,DATAONE,outvector,numr)
!!$     else
!!$        do ir=1,numr
!!$           outvector(ir,:)=outvector(ir,:)+invector(ir,www%botconfig:www%topconfig)*sparse_ptr%kefac*rkemod(ir,ir)
!!$        enddo
     endif
  endif

end subroutine sparsesparsemult_nompi



recursive subroutine arbitrary_sparsemult_singles_nompi(www,mattrans, rvector,inbigvector,outsmallvector,mynumr,diagflag)
  use walkmod
  use sparse_parameters
  use fileptrmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: mynumr,diagflag
  DATATYPE,intent(in) :: inbigvector(mynumr,www%numconfig), mattrans(www%singlematsize,www%botconfig:www%topconfig)
  DATATYPE,intent(out) :: outsmallvector(mynumr,www%botconfig:www%topconfig)
  DATAECS,intent(in) :: rvector(mynumr)
  integer :: iwalk,config1,idiag,ihop

  if (sparseconfigflag.eq.0) then
     OFLWR "error no arbitrary_sparsemult if sparseconfigflag.eq.0"; CFLST
  endif

  if (www%topconfig-www%botconfig+1.eq.0) then
     return
  endif

  outsmallvector(:,:)=0d0

  if (sparsehopflag.eq.0) then

     if (diagflag.eq.0) then
        do config1=www%botconfig,www%topconfig
           do iwalk=1,www%numsinglewalks(config1)
              outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(iwalk,config1) * inbigvector(:,www%singlewalk(iwalk,config1)) * rvector(:)
           enddo
        enddo
     else
        do config1=www%botconfig,www%topconfig
           do idiag=1,www%numsinglediagwalks(config1)
              outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(www%singlediag(idiag,config1),config1) * inbigvector(:,config1) * rvector(:)
           enddo
        enddo
     endif

  else

     if (diagflag.eq.0) then
        do config1=www%botconfig,www%topconfig
           do ihop=1,www%numsinglehops(config1)
              outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(ihop,config1) * inbigvector(:,www%singlehop(ihop,config1)) * rvector(:)
           enddo
        enddo
     else
        do config1=www%botconfig,www%topconfig
           outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(www%singlediaghop(config1),config1) * inbigvector(:,config1) * rvector(:)
        enddo
     endif

  endif


end subroutine arbitrary_sparsemult_singles_nompi


recursive subroutine arbitrary_sparsemult_doubles_nompi(www,mattrans,rvector,inbigvector,outsmallvector,mynumr,diagflag)
  use fileptrmod
  use sparse_parameters
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: mynumr,diagflag
  DATATYPE,intent(in) :: inbigvector(mynumr,www%numconfig),  mattrans(www%doublematsize,www%botconfig:www%topconfig)
  DATATYPE,intent(out) :: outsmallvector(mynumr,www%botconfig:www%topconfig)
  DATAECS,intent(in) :: rvector(mynumr)
  integer :: iwalk,config1,idiag,ihop

  if (sparseconfigflag.eq.0) then
     OFLWR "error no arbitrary_sparsemult_doubles if sparseconfigflag.eq.0"; CFLST
  endif

  if (www%topconfig-www%botconfig+1.eq.0) then
     return
  endif

  outsmallvector(:,:)=0d0

  if (sparsehopflag.eq.0) then

     if (diagflag.eq.0) then
        do config1=www%botconfig,www%topconfig
           do iwalk=1,www%numdoublewalks(config1)
              outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(iwalk,config1) * inbigvector(:,www%doublewalk(iwalk,config1)) * rvector(:)
           enddo
        enddo
     else
        do config1=www%botconfig,www%topconfig
           do idiag=1,www%numdoublediagwalks(config1)
              outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(www%doublediag(idiag,config1),config1) * inbigvector(:,config1) * rvector(:)
           enddo
        enddo
     endif

  else

     if (diagflag.eq.0) then
        do config1=www%botconfig,www%topconfig
           do ihop=1,www%numdoublehops(config1)
              outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(ihop,config1) * inbigvector(:,www%doublehop(ihop,config1)) * rvector(:)
           enddo
        enddo
     else
        do config1=www%botconfig,www%topconfig
           outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(www%doublediaghop(config1),config1) * inbigvector(:,config1) * rvector(:)
        enddo
     endif

  endif

end subroutine arbitrary_sparsemult_doubles_nompi



recursive subroutine arbitraryconfig_mult_singles_nompi(www,onebodymat, rvector, avectorin, avectorout,inrnum,diagflag)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: inrnum,diagflag
  DATATYPE,intent(in) :: onebodymat(www%nspf,www%nspf), avectorin(inrnum,www%numconfig)
  DATATYPE,intent(out) :: avectorout(inrnum,www%botconfig:www%topconfig)
  DATAECS,intent(in) :: rvector(inrnum)
  DATATYPE :: csum
  integer ::    config2, config1, iwalk, idiag,ihop

  if (www%topconfig-www%botconfig+1.eq.0) then
     return
  endif

  avectorout(:,:)=0.d0

  if (diagflag.eq.0) then

     do config1=www%botconfig,www%topconfig

!!$        do iwalk=1,www%numsinglewalks(config1)
!!$           config2=www%singlewalk(iwalk,config1)
!!$           avectorout(:,config1)=avectorout(:,config1)+avectorin(:,config2)*&
!!$                onebodymat(www%singlewalkopspf(1,iwalk,config1), &
!!$                www%singlewalkopspf(2,iwalk,config1)) *  &
!!$                www%singlewalkdirphase(iwalk,config1) * rvector(:)
!!$        enddo

        do ihop=1,www%numsinglehops(config1)
           config2=www%singlehop(ihop,config1)
           csum=0d0
           do iwalk=www%singlehopwalkstart(ihop,config1),www%singlehopwalkend(ihop,config1)
              csum=csum +  onebodymat(www%singlewalkopspf(1,iwalk,config1), &
                   www%singlewalkopspf(2,iwalk,config1)) *  &
                   www%singlewalkdirphase(iwalk,config1)
           enddo
           avectorout(:,config1)=avectorout(:,config1)+avectorin(:,config2) * csum * rvector(:)
        enddo
     enddo

  else

     do config1=www%botconfig,www%topconfig
        csum=0d0
        do idiag=1,www%numsinglediagwalks(config1)
           iwalk=www%singlediag(idiag,config1)
           csum=csum+&
                onebodymat(www%singlewalkopspf(1,iwalk,config1), &
                www%singlewalkopspf(2,iwalk,config1)) *  &
                www%singlewalkdirphase(iwalk,config1)
        enddo
        avectorout(:,config1)=avectorout(:,config1)+avectorin(:,config1) * csum * rvector(:)
     enddo

  endif

end subroutine arbitraryconfig_mult_singles_nompi



recursive subroutine arbitraryconfig_mult_doubles_nompi(www,twobodymat, rvector, avectorin, avectorout,inrnum,diagflag)   
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: inrnum,diagflag
  DATATYPE,intent(in) :: twobodymat(www%nspf,www%nspf,www%nspf,www%nspf), avectorin(inrnum,www%numconfig)
  DATATYPE,intent(out) :: avectorout(inrnum,www%botconfig:www%topconfig)
  DATAECS,intent(in) :: rvector(inrnum)
  DATATYPE :: csum
  integer ::   config2, config1, iwalk, idiag, ihop

  avectorout(:,:)=0.d0

  if (diagflag.eq.0) then

     do config1=www%botconfig,www%topconfig

!!$        do iwalk=1,www%numdoublewalks(config1)
!!$           config2=www%doublewalk(iwalk,config1)
!!$           avectorout(:,config1)=avectorout(:,config1)+avectorin(:,config2)*&
!!$                twobodymat(www%doublewalkdirspf(1,iwalk,config1), &
!!$                www%doublewalkdirspf(2,iwalk,config1),   &
!!$                www%doublewalkdirspf(3,iwalk,config1),   &
!!$                www%doublewalkdirspf(4,iwalk,config1))* &
!!$                www%doublewalkdirphase(iwalk,config1) * rvector(:)
!!$        enddo

        do ihop=1,www%numdoublehops(config1)
           config2=www%doublehop(ihop,config1)
           csum=0d0
           do iwalk=www%doublehopwalkstart(ihop,config1),www%doublehopwalkend(ihop,config1)
              csum=csum+&
                twobodymat(www%doublewalkdirspf(1,iwalk,config1), &
                www%doublewalkdirspf(2,iwalk,config1),   &
                www%doublewalkdirspf(3,iwalk,config1),   &
                www%doublewalkdirspf(4,iwalk,config1))* &
                www%doublewalkdirphase(iwalk,config1)
           enddo
           avectorout(:,config1)=avectorout(:,config1)+avectorin(:,config2) *  csum * rvector(:)
        enddo
     enddo

  else

     do config1=www%botconfig,www%topconfig
        csum=0d0
        do idiag=1,www%numdoublediagwalks(config1)
           iwalk=www%doublediag(idiag,config1)
           csum=csum+&
                twobodymat(www%doublewalkdirspf(1,iwalk,config1), &
                www%doublewalkdirspf(2,iwalk,config1),   &
                www%doublewalkdirspf(3,iwalk,config1),   &
                www%doublewalkdirspf(4,iwalk,config1))* &
                www%doublewalkdirphase(iwalk,config1)
        enddo
        avectorout(:,config1)=avectorout(:,config1)+avectorin(:,config1) * csum * rvector(:)
     enddo

  endif

end subroutine arbitraryconfig_mult_doubles_nompi



!!!!!!      BYPROC SUBROUTINES   !!!!!!!


recursive subroutine sparseconfigmult_byproc(whichproc,www,invector,outvector,matrix_ptr,sparse_ptr,&
     boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag)
  use fileptrmod   !! TEMP
  use sparse_parameters
  use configptrmod
  use sparseptrmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: whichproc
  integer,intent(in) :: conflag,boflag,nucflag,pulseflag,onlytdflag,botr,topr,diagflag
  DATATYPE,intent(in) :: invector(botr:topr,www%allbotconfigs(whichproc):www%alltopconfigs(whichproc))
  DATATYPE,intent(out) :: outvector(botr:topr,www%botconfig:www%topconfig)
  real*8,intent(in) :: time
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr

  if (www%topconfig-www%botconfig+1.eq.0) then
     return
  endif

  if (sparseopt.eq.0) then
     call direct_sparseconfigmult_byproc(whichproc,www,invector,outvector,matrix_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag)
  else
     call sparsesparsemult_byproc(whichproc,www,invector,outvector,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag)
  endif

end subroutine sparseconfigmult_byproc


!! DIRECT MATVEC (with matrix_ptr not sparse_ptr)

recursive subroutine direct_sparseconfigmult_byproc(whichproc,www,invector,outvector,matrix_ptr, &
     boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag)
  use ham_parameters
  use r_parameters
  use configptrmod
  use opmod   !! rkemod, proderivmod  
  use mpimod  !! myrank
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: whichproc
  integer,intent(in) :: conflag,boflag,nucflag,pulseflag,onlytdflag,botr,topr,diagflag
  DATATYPE,intent(in) :: invector(botr:topr,www%allbotconfigs(whichproc):www%alltopconfigs(whichproc))
  DATATYPE,intent(out) :: outvector(botr:topr,www%botconfig:www%topconfig)
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  real*8, intent(in) :: time
  integer :: mynumr,ir,flag,ii
  DATATYPE :: tempvector(botr:topr,www%botconfig:www%topconfig) !! AUTOMATIC
  DATATYPE :: tempmatel(www%nspf,www%nspf),facs(3),csum
  DATAECS :: rvector(botr:topr)
  real*8 :: gg

  if (www%topconfig-www%botconfig+1.eq.0) then
     return
  endif

  outvector(:,:)=0d0

  if (onlytdflag.ne.0) then
     facs(:)=0;     facs(onlytdflag)=1d0
  else
     call vectdpot(time,velflag,facs)
  endif
  
  mynumr=topr-botr+1
  
  if (boflag==1) then
     
!! easy nuclear repulsion, hardwired like this for now

     if (myrank.eq.whichproc) then
        do ir=botr,topr
           outvector(ir,:)=outvector(ir,:) + invector(ir,www%botconfig:www%topconfig) * ( &
                energyshift + (nucrepulsion+frozenpotdiag)/bondpoints(ir) + &
                frozenkediag/bondpoints(ir)**2 ) * matrix_ptr%kefac
        enddo
     endif

     rvector(:)=1/bondpoints(botr:topr)
     call arbitraryconfig_mult_doubles_byproc(whichproc,www,matrix_ptr%xtwoematel(:,:,:,:),rvector,invector,tempvector,mynumr,diagflag)
     outvector(:,:)=outvector(:,:)+tempvector(:,:)
     
     call arbitraryconfig_mult_singles_byproc(whichproc,www,matrix_ptr%xpotmatel(:,:),rvector,invector,tempvector,mynumr,diagflag)
     outvector(:,:)=outvector(:,:)+tempvector(:,:)
     
     rvector(:)=1/bondpoints(botr:topr)**2
     call arbitraryconfig_mult_singles_byproc(whichproc,www,matrix_ptr%xopmatel(:,:),rvector,invector,tempvector,mynumr,diagflag)
     outvector(:,:)=outvector(:,:)+tempvector(:,:)
     
  endif
  
  if (conflag==1.and.constraintflag.ne.0) then
     rvector(:)=1
     
     tempmatel(:,:)=matrix_ptr%xconmatel(:,:)
     if (tdflag.ne.0) then
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelxx(:,:)*facs(1)
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelyy(:,:)*facs(2)
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelzz(:,:)*facs(3)
     endif
     call arbitraryconfig_mult_singles_byproc(whichproc,www,tempmatel,rvector,invector,tempvector,mynumr,diagflag)
     outvector(:,:)=outvector(:,:)-tempvector(:,:)

  endif
  if (pulseflag.eq.1.and.tdflag.eq.1.or.onlytdflag.ne.0) then
     if (velflag.eq.0) then
        rvector(:)=bondpoints(botr:topr)
     else
        rvector(:)=1/bondpoints(botr:topr)
     endif

     tempmatel(:,:)=0d0; flag=0
     if (onlytdflag.eq.0.or.onlytdflag.eq.1) then
        tempmatel(:,:)= tempmatel(:,:) + matrix_ptr%xpulsematelxx*facs(1); flag=1
     endif
     if (onlytdflag.eq.0.or.onlytdflag.eq.2) then
        tempmatel(:,:)= tempmatel(:,:) + matrix_ptr%xpulsematelyy*facs(2); flag=1
     endif
     if (onlytdflag.eq.0.or.onlytdflag.eq.3) then
        tempmatel(:,:)= tempmatel(:,:) + matrix_ptr%xpulsematelzz*facs(3); flag=1
     endif
     if (flag.eq.1) then
        call arbitraryconfig_mult_singles_byproc(whichproc,www,tempmatel,rvector,invector,tempvector,mynumr,diagflag)
        outvector(:,:)=outvector(:,:)+tempvector(:,:)
     endif

!! STILL NEED A-VECTOR TERM IN BOND LENGTH FOR HETERO!! TO DO!!

     gg=0.25d0

     if (myrank.eq.whichproc) then

!!  A-squared term and nuclear dipole:  nuclear dipole always
        if (velflag.eq.0) then 
           csum=0
           do ii=1,3
              if (onlytdflag.eq.0.or.onlytdflag.eq.ii) then
                 csum=csum+facs(ii)*matrix_ptr%xpulsenuc(ii)
              endif
           enddo
           do ir=botr,topr
              outvector(ir,:)=outvector(ir,:)+invector(ir,www%botconfig:www%topconfig)*csum*bondpoints(ir)
           enddo
        else if (onlytdflag.eq.0) then
           csum=matrix_ptr%kefac * www%numelec * (facs(1)**2 + facs(2)**2 + facs(3)**2 ) * 2  * gg   !! a-squared
           do ir=botr,topr
              outvector(ir,:)=outvector(ir,:)+invector(ir,www%botconfig:www%topconfig)*csum     !! NO R FACTOR !!
           enddo
        endif

     endif

  endif   !! PULSE

  if (nonuc_checkflag.eq.0.and.nucflag.eq.1.and.mynumr.eq.numr) then
     if (diagflag.eq.0) then
        rvector(:)=1d0
        call arbitraryconfig_mult_singles_byproc(whichproc,www,matrix_ptr%xymatel,rvector,invector,tempvector,numr,0)
        call MYGEMM('N','N',numr,www%topconfig-www%botconfig+1,numr,DATAONE,            proderivmod,numr,     tempvector,           numr,DATAONE,outvector,numr)
        if (myrank.eq.whichproc) then
           call MYGEMM('N','N',numr,www%topconfig-www%botconfig+1,numr,DATAONE*matrix_ptr%kefac,rkemod,numr,     invector(:,www%botconfig),numr,DATAONE,outvector,numr)
        endif
!!$     else
!!$        do ir=1,numr
!!$           outvector(ir,:)=outvector(ir,:)+invector(ir,www%botconfig:www%topconfig)*matrix_ptr%kefac*rkemod(ir,ir)
!!$        enddo
     endif
  endif

end subroutine direct_sparseconfigmult_byproc


!! SPARSE MATVEC (with sparse_ptr not matrix_ptr)

recursive subroutine sparsesparsemult_byproc(whichproc,www,invector,outvector,sparse_ptr,&
     boflag,nucflag,pulseflag,conflag,time,onlytdflag,botr,topr,diagflag)
  use sparse_parameters
  use ham_parameters
  use r_parameters
  use sparseptrmod
  use walkmod
  use opmod   !! rkemod, proderivmod
  use mpimod  !! myrank 
  use fileptrmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: whichproc
  integer,intent(in) :: boflag,nucflag,pulseflag,conflag,onlytdflag,botr,topr,diagflag
  DATATYPE,intent(in)  :: invector(botr:topr,www%allbotconfigs(whichproc):www%alltopconfigs(whichproc))
  DATATYPE,intent(out) :: outvector(botr:topr,www%botconfig:www%topconfig)
  Type(SPARSEPTR),intent(in) :: sparse_ptr
  real*8,intent(in) :: time
  DATATYPE :: facs(3),csum
  DATATYPE :: tempvector(botr:topr,www%botconfig:www%topconfig), &
       tempsparsemattr(www%singlematsize,www%botconfig:www%topconfig) !! AUTOMATIC
  real*8 :: gg
  DATAECS :: rvector(botr:topr)
  integer :: ir,mynumr,ii,flag

  if (sparseconfigflag.eq.0) then
     OFLWR "error no sparsesparsemult if sparseconfigflag.eq.0"; CFLST
  endif

  if (www%topconfig-www%botconfig+1.eq.0) then
     return
  endif

  if (botr.ne.1.or.topr.ne.numr) then
     OFLWR "REPROGRAMMMM."; CFLST
  endif

  outvector(:,:)=0d0

  if (onlytdflag.ne.0) then
     facs(:)=0;     facs(onlytdflag)=1d0
  else
     call vectdpot(time,velflag,facs)
  endif

  mynumr=topr-botr+1

  if (boflag==1) then

     if (myrank.eq.whichproc) then
        do ir=botr,topr
           outvector(ir,:)=outvector(ir,:) + invector(ir,www%botconfig:www%topconfig) * ( &
                energyshift + (nucrepulsion+frozenpotdiag)/bondpoints(ir) + &
                frozenkediag/bondpoints(ir)**2 ) * sparse_ptr%kefac
        enddo
     endif

     rvector(:)=1/bondpoints(botr:topr)
     call arbitrary_sparsemult_doubles_byproc(whichproc,www,sparse_ptr%xpotsparsemattr(:,:),rvector,invector,tempvector,mynumr,diagflag)
     outvector(:,:)=outvector(:,:)+tempvector(:,:)

     call arbitrary_sparsemult_singles_byproc(whichproc,www,sparse_ptr%xonepotsparsemattr(:,:),rvector,invector,tempvector,mynumr,diagflag)
     outvector(:,:)=outvector(:,:)+tempvector(:,:)

     rvector(:)=1/bondpoints(botr:topr)**2
     call arbitrary_sparsemult_singles_byproc(whichproc,www,sparse_ptr%xopsparsemattr(:,:),rvector,invector,tempvector,mynumr,diagflag)
     outvector(:,:)=outvector(:,:)+tempvector(:,:)

  endif

  if (conflag==1.and.constraintflag.ne.0) then
     rvector(:)=1

     tempsparsemattr(:,:)=sparse_ptr%xconsparsemattr(:,:)
     if (tdflag.ne.0) then
        tempsparsemattr(:,:)=tempsparsemattr(:,:)+sparse_ptr%xconsparsemattrxx(:,:)*facs(1)
        tempsparsemattr(:,:)=tempsparsemattr(:,:)+sparse_ptr%xconsparsemattryy(:,:)*facs(2)
        tempsparsemattr(:,:)=tempsparsemattr(:,:)+sparse_ptr%xconsparsemattrzz(:,:)*facs(3)
     endif
     call arbitrary_sparsemult_singles_byproc(whichproc,www,tempsparsemattr,rvector,invector,tempvector,mynumr,diagflag)
     outvector(:,:)=outvector(:,:)-tempvector(:,:)

  endif
  if (pulseflag.eq.1.and.tdflag.eq.1.or.onlytdflag.ne.0) then
     if (velflag.eq.0) then
        rvector(:)=bondpoints(botr:topr)
     else
        rvector(:)=1/bondpoints(botr:topr)
     endif

     tempsparsemattr(:,:)=0d0; flag=0
     if (onlytdflag.eq.0.or.onlytdflag.eq.1) then
        tempsparsemattr(:,:)= tempsparsemattr(:,:) + sparse_ptr%xpulsesparsemattrxx*facs(1); flag=1
     endif
     if (onlytdflag.eq.0.or.onlytdflag.eq.2) then
        tempsparsemattr(:,:)= tempsparsemattr(:,:) + sparse_ptr%xpulsesparsemattryy*facs(2); flag=1
     endif
     if (onlytdflag.eq.0.or.onlytdflag.eq.3) then
        tempsparsemattr(:,:)= tempsparsemattr(:,:) + sparse_ptr%xpulsesparsemattrzz*facs(3); flag=1
     endif
     if (flag.eq.1) then
        call arbitrary_sparsemult_singles_byproc(whichproc,www,tempsparsemattr,rvector,invector,tempvector,mynumr,diagflag)
        outvector(:,:)=outvector(:,:)+tempvector(:,:)
     endif

!! STILL NEED A-VECTOR TERM IN BOND LENGTH FOR HETERO !! TO DO !!

!! A-squared term and nuclear dipole : nuclear dipole always

     gg=0.25d0
     
     if (myrank.eq.whichproc) then

        if (velflag.eq.0) then 
           csum=0
           do ii=1,3
              if (onlytdflag.eq.0.or.onlytdflag.eq.ii) then
                 csum=csum+facs(ii)*sparse_ptr%xpulsenuc(ii)
              endif
           enddo
           do ir=botr,topr
              outvector(ir,:)=outvector(ir,:)+invector(ir,www%botconfig:www%topconfig)*csum*bondpoints(ir)
           enddo
        else if (onlytdflag.eq.0) then
           csum= sparse_ptr%kefac * www%numelec * (facs(1)**2 + facs(2)**2 + facs(3)**2 ) * 2 * gg  !! a-squared
           do ir=botr,topr
              outvector(ir,:)=outvector(ir,:)+invector(ir,www%botconfig:www%topconfig)*csum     !! NO R FACTOR !!
           enddo
        endif

     endif

  endif !! PULSE

  if (nonuc_checkflag.eq.0.and.nucflag.eq.1.and.mynumr.eq.numr) then
     if (diagflag.eq.0) then
        rvector(:)=1d0
        call arbitrary_sparsemult_singles_byproc(whichproc,www,sparse_ptr%xysparsemattr,rvector,invector,tempvector,numr,0)
        call MYGEMM('N','N',numr,www%topconfig-www%botconfig+1,numr,DATAONE,            proderivmod,numr,     tempvector,           numr,DATAONE,outvector,numr)
        if (myrank.eq.whichproc) then
           call MYGEMM('N','N',numr,www%topconfig-www%botconfig+1,numr,DATAONE*sparse_ptr%kefac,rkemod,numr,     invector(:,www%botconfig),numr,DATAONE,outvector,numr)
        endif
!!$     else
!!$        do ir=1,numr
!!$           outvector(ir,:)=outvector(ir,:)+invector(ir,www%botconfig:www%topconfig)*sparse_ptr%kefac*rkemod(ir,ir)
!!$        enddo
     endif
  endif

end subroutine sparsesparsemult_byproc



recursive subroutine arbitrary_sparsemult_singles_byproc(whichproc,www,mattrans, rvector,insmallvector,outsmallvector,mynumr,diagflag)
  use walkmod
  use sparse_parameters
  use fileptrmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: whichproc
  integer,intent(in) :: mynumr,diagflag
  DATATYPE,intent(in) :: insmallvector(mynumr,www%allbotconfigs(whichproc):www%alltopconfigs(whichproc)), &
       mattrans(www%singlematsize,www%botconfig:www%topconfig)
  DATATYPE,intent(out) :: outsmallvector(mynumr,www%botconfig:www%topconfig)
  DATAECS,intent(in) :: rvector(mynumr)
  integer :: iwalk,config1,idiag,ihop

  if (sparseconfigflag.eq.0) then
     OFLWR "error no arbitrary_sparsemult if sparseconfigflag.eq.0"; CFLST
  endif

  if (www%topconfig-www%botconfig+1.eq.0) then
     return
  endif

  outsmallvector(:,:)=0d0

  if (sparsehopflag.eq.0) then

     if (diagflag.eq.0) then
        do config1=www%botconfig,www%topconfig
           do iwalk=www%firstsinglewalkbyproc(whichproc,config1),www%lastsinglewalkbyproc(whichproc,config1)
              outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(iwalk,config1) * insmallvector(:,www%singlewalk(iwalk,config1)) * rvector(:)
           enddo
        enddo
     else
        do config1=www%botconfig,www%topconfig
           do idiag=1,www%numsinglediagwalks(config1)
              outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(www%singlediag(idiag,config1),config1) * insmallvector(:,config1) * rvector(:)
           enddo
        enddo
     endif

  else

     if (diagflag.eq.0) then
        do config1=www%botconfig,www%topconfig
           do ihop=www%firstsinglehopbyproc(whichproc,config1),www%lastsinglehopbyproc(whichproc,config1)
              outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(ihop,config1) * insmallvector(:,www%singlehop(ihop,config1)) * rvector(:)
           enddo
        enddo
     else
        do config1=www%botconfig,www%topconfig
           outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(www%singlediaghop(config1),config1) * insmallvector(:,config1) * rvector(:)
        enddo
     endif

  endif


end subroutine arbitrary_sparsemult_singles_byproc


recursive subroutine arbitrary_sparsemult_doubles_byproc(whichproc,www,mattrans,rvector,insmallvector,outsmallvector,mynumr,diagflag)
  use fileptrmod
  use sparse_parameters
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: whichproc
  integer,intent(in) :: mynumr,diagflag
  DATATYPE,intent(in) :: insmallvector(mynumr,www%allbotconfigs(whichproc):www%alltopconfigs(whichproc)), &
       mattrans(www%doublematsize,www%botconfig:www%topconfig)
  DATATYPE,intent(out) :: outsmallvector(mynumr,www%botconfig:www%topconfig)
  DATAECS,intent(in) :: rvector(mynumr)
  integer :: iwalk,config1,idiag,ihop

  if (sparseconfigflag.eq.0) then
     OFLWR "error no arbitrary_sparsemult_doubles if sparseconfigflag.eq.0"; CFLST
  endif

  if (www%topconfig-www%botconfig+1.eq.0) then
     return
  endif

  outsmallvector(:,:)=0d0

  if (sparsehopflag.eq.0) then

     if (diagflag.eq.0) then
        do config1=www%botconfig,www%topconfig
           do iwalk=www%firstdoublewalkbyproc(whichproc,config1),www%lastdoublewalkbyproc(whichproc,config1)
              outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(iwalk,config1) * insmallvector(:,www%doublewalk(iwalk,config1)) * rvector(:)
           enddo
        enddo
     else
        do config1=www%botconfig,www%topconfig
           do idiag=1,www%numdoublediagwalks(config1)
              outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(www%doublediag(idiag,config1),config1) * insmallvector(:,config1) * rvector(:)
           enddo
        enddo
     endif

  else

     if (diagflag.eq.0) then
        do config1=www%botconfig,www%topconfig
           do ihop=www%firstdoublehopbyproc(whichproc,config1),www%lastdoublehopbyproc(whichproc,config1)
              outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(ihop,config1) * insmallvector(:,www%doublehop(ihop,config1)) * rvector(:)
           enddo
        enddo
     else
        do config1=www%botconfig,www%topconfig
           outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(www%doublediaghop(config1),config1) * insmallvector(:,config1) * rvector(:)
        enddo
     endif

  endif

end subroutine arbitrary_sparsemult_doubles_byproc


recursive subroutine arbitraryconfig_mult_singles_byproc(whichproc,www,onebodymat, rvector, avectorin, avectorout,inrnum,diagflag)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: whichproc
  integer,intent(in) :: inrnum,diagflag
  DATATYPE,intent(in) :: onebodymat(www%nspf,www%nspf), avectorin(inrnum,www%allbotconfigs(whichproc):www%alltopconfigs(whichproc))
  DATATYPE,intent(out) :: avectorout(inrnum,www%botconfig:www%topconfig)
  DATAECS,intent(in) :: rvector(inrnum)
  DATATYPE :: csum
  integer ::    config2, config1, iwalk, idiag,ihop

  if (www%topconfig-www%botconfig+1.eq.0) then
     return
  endif

  avectorout(:,:)=0.d0

  if (diagflag.eq.0) then

     do config1=www%botconfig,www%topconfig

        do ihop=www%firstsinglehopbyproc(whichproc,config1),www%lastsinglehopbyproc(whichproc,config1)
           config2=www%singlehop(ihop,config1)
           csum=0d0
           do iwalk=www%singlehopwalkstart(ihop,config1),www%singlehopwalkend(ihop,config1)
              csum=csum +  onebodymat(www%singlewalkopspf(1,iwalk,config1), &
                   www%singlewalkopspf(2,iwalk,config1)) *  &
                   www%singlewalkdirphase(iwalk,config1)
           enddo
           avectorout(:,config1)=avectorout(:,config1)+avectorin(:,config2) * csum * rvector(:)
        enddo
     enddo

  else

     do config1=www%botconfig,www%topconfig
        csum=0d0
        do idiag=1,www%numsinglediagwalks(config1)
           iwalk=www%singlediag(idiag,config1)
           csum=csum+&
                onebodymat(www%singlewalkopspf(1,iwalk,config1), &
                www%singlewalkopspf(2,iwalk,config1)) *  &
                www%singlewalkdirphase(iwalk,config1)
        enddo
        avectorout(:,config1)=avectorout(:,config1)+avectorin(:,config1) * csum * rvector(:)
     enddo

  endif

end subroutine arbitraryconfig_mult_singles_byproc



recursive subroutine arbitraryconfig_mult_doubles_byproc(whichproc,www,twobodymat, rvector, avectorin, avectorout,inrnum,diagflag)   
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: whichproc
  integer,intent(in) :: inrnum,diagflag
  DATATYPE,intent(in) :: twobodymat(www%nspf,www%nspf,www%nspf,www%nspf), avectorin(inrnum,www%allbotconfigs(whichproc):www%alltopconfigs(whichproc))
  DATATYPE,intent(out) :: avectorout(inrnum,www%botconfig:www%topconfig)
  DATAECS,intent(in) :: rvector(inrnum)
  DATATYPE :: csum
  integer ::   config2, config1, iwalk, idiag, ihop

  avectorout(:,:)=0.d0

  if (diagflag.eq.0) then

     do config1=www%botconfig,www%topconfig

        do ihop=www%firstdoublehopbyproc(whichproc,config1),www%lastdoublehopbyproc(whichproc,config1)
           config2=www%doublehop(ihop,config1)
           csum=0d0
           do iwalk=www%doublehopwalkstart(ihop,config1),www%doublehopwalkend(ihop,config1)
              csum=csum+&
                twobodymat(www%doublewalkdirspf(1,iwalk,config1), &
                www%doublewalkdirspf(2,iwalk,config1),   &
                www%doublewalkdirspf(3,iwalk,config1),   &
                www%doublewalkdirspf(4,iwalk,config1))* &
                www%doublewalkdirphase(iwalk,config1)
           enddo
           avectorout(:,config1)=avectorout(:,config1)+avectorin(:,config2) *  csum * rvector(:)
        enddo
     enddo

  else

     do config1=www%botconfig,www%topconfig
        csum=0d0
        do idiag=1,www%numdoublediagwalks(config1)
           iwalk=www%doublediag(idiag,config1)
           csum=csum+&
                twobodymat(www%doublewalkdirspf(1,iwalk,config1), &
                www%doublewalkdirspf(2,iwalk,config1),   &
                www%doublewalkdirspf(3,iwalk,config1),   &
                www%doublewalkdirspf(4,iwalk,config1))* &
                www%doublewalkdirphase(iwalk,config1)
        enddo
        avectorout(:,config1)=avectorout(:,config1)+avectorin(:,config1) * csum * rvector(:)
     enddo

  endif

end subroutine arbitraryconfig_mult_doubles_byproc

