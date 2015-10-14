
#include "Definitions.INC"

!! All purpose subroutine.

recursive subroutine sparseconfigmult(myinvector,myoutvector,matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time)
  use parameters
  use configptrmod
  use sparseptrmod
  implicit none
  integer,intent(in) :: conflag,boflag,nucflag,pulseflag
  DATATYPE,intent(in) :: myinvector(totadim)
  DATATYPE,intent(out) :: myoutvector(totadim)
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr
  real*8,intent(in) :: time

  call sparseconfigmultxxx(myinvector,myoutvector,matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,0)

end subroutine sparseconfigmult


!! For one bond length, born oppenheimer Hamiltonian

recursive subroutine sparseconfigmultone(invector,outvector,matrix_ptr,sparse_ptr, boflag,pulseflag,conflag,isplit,time)
  use parameters
  use configptrmod
  use sparseptrmod
  implicit none
  integer,intent(in) :: conflag, boflag, pulseflag, isplit
  DATATYPE,intent(in) :: invector(firstconfig:lastconfig)
  DATATYPE,intent(out) :: outvector(firstconfig:lastconfig)
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr
  real*8,intent(in) :: time
  DATATYPE :: workvector(numconfig)      !! AUTOMATIC

  workvector(firstconfig:lastconfig)=invector(:)
  if (parconsplit.ne.0) then
     call mpiallgather(workvector,numconfig,configsperproc(:),maxconfigsperproc)
  endif
  call sparseconfigmult_nompi(workvector,outvector(botconfig),matrix_ptr,sparse_ptr, boflag, 0, pulseflag, conflag,time,0,isplit,isplit,0)
  if (parconsplit.eq.0) then
     call mpiallgather(outvector,numconfig,configsperproc,maxconfigsperproc)
  endif

end subroutine sparseconfigmultone


!! All purpose subroutine, just the pulse.

recursive subroutine sparseconfigpulsemult(myinvector,myoutvector,matrix_ptr, sparse_ptr,which)
  use parameters
  use configptrmod
  use sparseptrmod
  implicit none
  integer,intent(in) :: which
  DATATYPE,intent(in) :: myinvector(totadim)
  DATATYPE,intent(out) :: myoutvector(totadim)
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr

  call sparseconfigmultxxx(myinvector,myoutvector,matrix_ptr,sparse_ptr, 0,0,1,0,-1d0,which)

end subroutine sparseconfigpulsemult


!! MPI subroutine

recursive subroutine sparseconfigmultxxx(myinvector,myoutvector,matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag)
  use parameters
  use configptrmod
  use sparseptrmod
  implicit none
  integer,intent(in) :: conflag,boflag,nucflag,pulseflag,onlytdflag
  DATATYPE,intent(in) :: myinvector(numr,firstconfig:lastconfig)
  DATATYPE,intent(out) :: myoutvector(numr,firstconfig:lastconfig)
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr
  real*8,intent(in) :: time
  DATATYPE :: invector(numr,numconfig)    !!AUTOMATIC

  invector(:,firstconfig:lastconfig) = myinvector(:,:)
  if (parconsplit.ne.0) then
     call mpiallgather(invector,numconfig*numr,configsperproc(:)*numr,maxconfigsperproc*numr)
  endif

  call sparseconfigmult_nompi(invector,myoutvector(:,botconfig),matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag,1,numr,0)

  if (parconsplit.eq.0) then
     call mpiallgather(myoutvector,numconfig*numr,configsperproc*numr,maxconfigsperproc*numr)
  endif
     
end subroutine sparseconfigmultxxx



recursive subroutine parsparseconfigpreconmult(invector,outvector,matrix_ptr, sparse_ptr, boflag, nucflag, pulseflag, conflag, inenergy,time)
  use parameters
  use configptrmod
  use sparseptrmod
  implicit none
  DATATYPE,intent(in) :: invector(numr,botconfig:topconfig)
  DATATYPE,intent(out) :: outvector(numr,botconfig:topconfig)
  DATATYPE,intent(in) :: inenergy
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr
  integer,intent(in) ::  conflag,boflag,nucflag,pulseflag
  real*8,intent(in) :: time
  DATATYPE :: workvector(numr,numconfig),tempvector(numr,botconfig:topconfig)     !!AUTOMATIC

  if (topconfig-botconfig+1.eq.0) then
     return
  endif

  workvector(:,:)=1d0; tempvector(:,:)=0d0
  call sparseconfigmult_nompi(workvector,tempvector,matrix_ptr, sparse_ptr, boflag, nucflag, pulseflag, conflag,time,0,1,numr,1)
  
  outvector(:,:)=invector(:,:)/(tempvector(:,:)-inenergy)

end subroutine parsparseconfigpreconmult



!! NOW WRAPPER - SPARSEOPT HERE

recursive subroutine sparseconfigmult_nompi(myinvector,myoutvector,matrix_ptr,sparse_ptr,&
     boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag)
  use parameters
  use configptrmod
  use sparseptrmod
  implicit none
  integer,intent(in) :: conflag,boflag,nucflag,pulseflag,onlytdflag,botr,topr,diagflag
  DATATYPE,intent(in) :: myinvector(botr:topr,numconfig)
  DATATYPE,intent(out) :: myoutvector(botr:topr,botconfig:topconfig)
  real*8,intent(in) :: time
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  Type(SPARSEPTR),intent(in) :: sparse_ptr

  if (topconfig-botconfig+1.eq.0) then
     return
  endif
  if (sparseopt.eq.0) then
     call direct_sparseconfigmult_nompi(myinvector,myoutvector,matrix_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag)
  else
     call sparsesparsemult_nompi(myinvector,myoutvector,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag)
  endif

end subroutine sparseconfigmult_nompi


!! SPARSE MATVEC (with matrix_ptr not sparse_ptr)

recursive subroutine direct_sparseconfigmult_nompi(myinvector,myoutvector,matrix_ptr, &
     boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag)
  use parameters
  use configptrmod
  use opmod   !! rkemod, proderivmod  
  implicit none
  integer,intent(in) :: conflag,boflag,nucflag,pulseflag,onlytdflag,botr,topr,diagflag
  DATATYPE,intent(in) :: myinvector(botr:topr,numconfig)
  DATATYPE,intent(out) :: myoutvector(botr:topr,botconfig:topconfig)
  Type(CONFIGPTR),intent(in) :: matrix_ptr
  real*8, intent(in) :: time
  integer :: mynumr,ir,flag,ii
  DATATYPE :: tempvector(botr:topr,botconfig:topconfig) !! AUTOMATIC
  DATATYPE :: tempmatel(nspf,nspf),facs(3),csum
  DATAECS :: rvector(botr:topr)
  real*8 :: gg

  if (topconfig-botconfig+1.eq.0) then
     return
  endif

  myoutvector(:,:)=0d0

  if (onlytdflag.ne.0) then
     facs(:)=0;     facs(onlytdflag)=1d0
  else
     call vectdpot(time,velflag,facs)
  endif
  
  mynumr=topr-botr+1
  
  if (boflag==1) then
     
!! easy nuclear repulsion, hardwired like this for now

     do ir=botr,topr

        myoutvector(ir,:)=myoutvector(ir,:) + myinvector(ir,botconfig:topconfig) * ( &
             energyshift + (nucrepulsion+frozenpotdiag)/bondpoints(ir) + &
             frozenkediag/bondpoints(ir)**2 ) * matrix_ptr%kefac

     enddo
     
     rvector(:)=1/bondpoints(botr:topr)
     call arbitraryconfig_mult_doubles_nompi(matrix_ptr%xtwoematel(:,:,:,:),rvector,myinvector,tempvector,mynumr,diagflag)
     myoutvector(:,:)=myoutvector(:,:)+tempvector(:,:)
     
     call arbitraryconfig_mult_singles_nompi(matrix_ptr%xpotmatel(:,:),rvector,myinvector,tempvector,mynumr,diagflag)
     myoutvector(:,:)=myoutvector(:,:)+tempvector(:,:)
     
     rvector(:)=1/bondpoints(botr:topr)**2
     call arbitraryconfig_mult_singles_nompi(matrix_ptr%xopmatel(:,:),rvector,myinvector,tempvector,mynumr,diagflag)
     myoutvector(:,:)=myoutvector(:,:)+tempvector(:,:)
     
  endif
  
  if (conflag==1.and.constraintflag.ne.0) then
     rvector(:)=1
     
     tempmatel(:,:)=matrix_ptr%xconmatel(:,:)
     if (tdflag.ne.0) then
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelxx(:,:)*facs(1)
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelyy(:,:)*facs(2)
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelzz(:,:)*facs(3)
     endif
     call arbitraryconfig_mult_singles_nompi(tempmatel,rvector,myinvector,tempvector,mynumr,diagflag)
     myoutvector(:,:)=myoutvector(:,:)-tempvector(:,:)

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
        call arbitraryconfig_mult_singles_nompi(tempmatel,rvector,myinvector,tempvector,mynumr,diagflag)
        myoutvector(:,:)=myoutvector(:,:)+tempvector(:,:)
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
           myoutvector(ir,:)=myoutvector(ir,:)+myinvector(ir,botconfig:topconfig)*csum*bondpoints(ir)
        enddo
     else if (onlytdflag.eq.0) then
        csum=matrix_ptr%kefac * numelec * (facs(1)**2 + facs(2)**2 + facs(3)**2 ) * 2  * gg   !! a-squared
        do ir=botr,topr
           myoutvector(ir,:)=myoutvector(ir,:)+myinvector(ir,botconfig:topconfig)*csum     !! NO R FACTOR !!
        enddo
     endif

  endif   !! PULSE

  if (nonuc_checkflag.eq.0.and.nucflag.eq.1.and.mynumr.eq.numr) then
     if (diagflag.eq.0) then
        rvector(:)=1d0
        call arbitraryconfig_mult_singles_nompi(matrix_ptr%xymatel,rvector,myinvector,tempvector,numr,0)
        call MYGEMM('N','N',numr,topconfig-botconfig+1,numr,DATAONE,            proderivmod,numr,     tempvector,           numr,DATAONE,myoutvector,numr)
        call MYGEMM('N','N',numr,topconfig-botconfig+1,numr,DATAONE*matrix_ptr%kefac,rkemod,numr,     myinvector(:,botconfig),numr,DATAONE,myoutvector,numr)
!!$     else
!!$        do ir=1,numr
!!$           myoutvector(ir,:)=myoutvector(ir,:)+myinvector(ir,botconfig:topconfig)*matrix_ptr%kefac*rkemod(ir,ir)
!!$        enddo
     endif
  endif

end subroutine direct_sparseconfigmult_nompi


!! DIRECT MATVEC (with sparse_ptr not matrix_ptr)

recursive subroutine sparsesparsemult_nompi(myinvector,myoutvector,sparse_ptr,&
     boflag,nucflag,pulseflag,conflag,time,onlytdflag,botr,topr,diagflag)
  use parameters
  use sparseptrmod
  use walkmod
  use opmod   !! rkemod, proderivmod
  implicit none
  integer,intent(in) :: boflag,nucflag,pulseflag,conflag,onlytdflag,botr,topr,diagflag
  DATATYPE,intent(in)  :: myinvector(botr:topr,numconfig)
  DATATYPE,intent(out) :: myoutvector(botr:topr,botconfig:topconfig)
  Type(SPARSEPTR),intent(in) :: sparse_ptr
  real*8,intent(in) :: time
  DATATYPE :: facs(3),csum
  DATATYPE :: tempvector(botr:topr,botconfig:topconfig), tempsparsemattr(maxsinglewalks,botconfig:topconfig) !! AUTOMATIC
  real*8 :: gg
  DATAECS :: rvector(botr:topr)
  integer :: ir,mynumr,ii,flag

  if (topconfig-botconfig+1.eq.0) then
     return
  endif

  if (botr.ne.1.or.topr.ne.numr) then
     OFLWR "REPROGRAMMMM."; CFLST
  endif

  myoutvector(:,:)=0d0

  if (onlytdflag.ne.0) then
     facs(:)=0;     facs(onlytdflag)=1d0
  else
     call vectdpot(time,velflag,facs)
  endif

  mynumr=topr-botr+1

  if (boflag==1) then

    do ir=botr,topr

       myoutvector(ir,:)=myoutvector(ir,:) + myinvector(ir,botconfig:topconfig) * ( &
            energyshift + (nucrepulsion+frozenpotdiag)/bondpoints(ir) + &
            frozenkediag/bondpoints(ir)**2 ) * sparse_ptr%kefac
       
    enddo

     rvector(:)=1/bondpoints(botr:topr)
     call arbitrary_sparsemult_nompi_doubles(sparse_ptr%xpotsparsemattr(:,:),rvector,myinvector,tempvector,mynumr,diagflag)
     myoutvector(:,:)=myoutvector(:,:)+tempvector(:,:)

     call arbitrary_sparsemult_nompi_singles(sparse_ptr%xonepotsparsemattr(:,:),rvector,myinvector,tempvector,mynumr,diagflag)
     myoutvector(:,:)=myoutvector(:,:)+tempvector(:,:)

     rvector(:)=1/bondpoints(botr:topr)**2
     call arbitrary_sparsemult_nompi_singles(sparse_ptr%xopsparsemattr(:,:),rvector,myinvector,tempvector,mynumr,diagflag)
     myoutvector(:,:)=myoutvector(:,:)+tempvector(:,:)

  endif

  if (conflag==1.and.constraintflag.ne.0) then
     rvector(:)=1

     tempsparsemattr(:,:)=sparse_ptr%xconsparsemattr(:,:)
     if (tdflag.ne.0) then
        tempsparsemattr(:,:)=tempsparsemattr(:,:)+sparse_ptr%xconsparsemattrxx(:,:)*facs(1)
        tempsparsemattr(:,:)=tempsparsemattr(:,:)+sparse_ptr%xconsparsemattryy(:,:)*facs(2)
        tempsparsemattr(:,:)=tempsparsemattr(:,:)+sparse_ptr%xconsparsemattrzz(:,:)*facs(3)
     endif
     call arbitrary_sparsemult_nompi_singles(tempsparsemattr,rvector,myinvector,tempvector,mynumr,diagflag)
     myoutvector(:,:)=myoutvector(:,:)-tempvector(:,:)

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
        call arbitrary_sparsemult_nompi_singles(tempsparsemattr,rvector,myinvector,tempvector,mynumr,diagflag)
        myoutvector(:,:)=myoutvector(:,:)+tempvector(:,:)
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
           myoutvector(ir,:)=myoutvector(ir,:)+myinvector(ir,botconfig:topconfig)*csum*bondpoints(ir)
        enddo
     else if (onlytdflag.eq.0) then
        csum= sparse_ptr%kefac * numelec * (facs(1)**2 + facs(2)**2 + facs(3)**2 ) * 2 * gg  !! a-squared
        do ir=botr,topr
           myoutvector(ir,:)=myoutvector(ir,:)+myinvector(ir,botconfig:topconfig)*csum     !! NO R FACTOR !!
        enddo
     endif

  endif !! PULSE

  if (nonuc_checkflag.eq.0.and.nucflag.eq.1.and.mynumr.eq.numr) then
     if (diagflag.eq.0) then
        rvector(:)=1d0
        call arbitrary_sparsemult_nompi_singles(sparse_ptr%xysparsemattr,rvector,myinvector,tempvector,numr,0)
        call MYGEMM('N','N',numr,topconfig-botconfig+1,numr,DATAONE,            proderivmod,numr,     tempvector,           numr,DATAONE,myoutvector,numr)
        call MYGEMM('N','N',numr,topconfig-botconfig+1,numr,DATAONE*sparse_ptr%kefac,rkemod,numr,     myinvector(:,botconfig),numr,DATAONE,myoutvector,numr)
!!$     else
!!$        do ir=1,numr
!!$           myoutvector(ir,:)=myoutvector(ir,:)+myinvector(ir,botconfig:topconfig)*sparse_ptr%kefac*rkemod(ir,ir)
!!$        enddo
     endif
  endif

end subroutine sparsesparsemult_nompi



recursive subroutine arbitrary_sparsemult_nompi_singles(mattrans, rvector,inbigvector,outsmallvector,mynumr,diagflag)
  use parameters
  use walkmod
  implicit none
  integer,intent(in) :: mynumr,diagflag
  DATATYPE,intent(in) :: inbigvector(mynumr,numconfig), mattrans(maxsinglewalks,botconfig:topconfig)
  DATATYPE,intent(out) :: outsmallvector(mynumr,botconfig:topconfig)
  DATAECS :: rvector(mynumr)
  integer :: iwalk,config1,idiag

  if (topconfig-botconfig+1.eq.0) then
     return
  endif

  outsmallvector(:,:)=0d0

  if (diagflag.eq.0) then
     do config1=botconfig,topconfig
        do iwalk=1,numsinglewalks(config1)
           outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(iwalk,config1) * inbigvector(:,singlewalk(iwalk,config1))
        enddo
        outsmallvector(:,config1)=outsmallvector(:,config1)*rvector(:)
     enddo
  else
     do config1=botconfig,topconfig
        do idiag=1,numsinglediagwalks(config1)
           outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(singlediag(idiag,config1),config1) * inbigvector(:,config1)
        enddo
        outsmallvector(:,config1)=outsmallvector(:,config1)*rvector(:)
     enddo
  endif

end subroutine arbitrary_sparsemult_nompi_singles


recursive subroutine arbitrary_sparsemult_nompi_doubles(mattrans,rvector,inbigvector,outsmallvector,mynumr,diagflag)
  use parameters
  use walkmod
  implicit none
  integer,intent(in) :: mynumr,diagflag
  DATATYPE,intent(in) :: inbigvector(mynumr,numconfig),  mattrans(maxdoublewalks,botconfig:topconfig)
  DATATYPE,intent(out) :: outsmallvector(mynumr,botconfig:topconfig)
  DATAECS :: rvector(mynumr)
  integer :: iwalk,config1,idiag

  if (topconfig-botconfig+1.eq.0) then
     return
  endif

  outsmallvector(:,:)=0d0

  do config1=botconfig,topconfig
     if (diagflag.eq.0) then
        do iwalk=1,numdoublewalks(config1)
           outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(iwalk,config1) * inbigvector(:,doublewalk(iwalk,config1))
        enddo
     else
        do idiag=1,numdoublediagwalks(config1)
           outsmallvector(:,config1)=outsmallvector(:,config1)+mattrans(doublediag(idiag,config1),config1) * inbigvector(:,config1)
        enddo
     endif
     outsmallvector(:,config1)=outsmallvector(:,config1)*rvector(:)
  enddo

end subroutine arbitrary_sparsemult_nompi_doubles


recursive subroutine arbitraryconfig_mult_singles(onebodymat, rvector, avectorin, avectorout,inrnum)   
  use parameters
  implicit none
  integer,intent(in) :: inrnum
  DATATYPE,intent(in) :: onebodymat(nspf,nspf), avectorin(inrnum,firstconfig:lastconfig)
  DATATYPE,intent(out) :: avectorout(inrnum,firstconfig:lastconfig)
  DATATYPE :: workvector(inrnum,numconfig)         !! AUTOMATIC
  DATAECS :: rvector(inrnum)

  workvector(:,:)=0d0
  workvector(:,firstconfig:lastconfig)=avectorin(:,:)

!! DO SUMMA
  if (parconsplit.ne.0) then
     call mpiallgather(workvector,numconfig*inrnum,configsperproc(:)*inrnum,maxconfigsperproc*inrnum)
  endif

  call arbitraryconfig_mult_singles_nompi(onebodymat, rvector, workvector, avectorout(:,botconfig),inrnum,0)
  if (parconsplit.eq.0) then
     call mpiallgather(avectorout,numconfig*inrnum,configsperproc(:)*inrnum,maxconfigsperproc*inrnum)
  endif

end subroutine arbitraryconfig_mult_singles


recursive subroutine arbitraryconfig_mult_singles_nompi(onebodymat, rvector, avectorin, avectorout,inrnum,diagflag)
   use walkmod
  use parameters
  implicit none
  integer,intent(in) :: inrnum,diagflag
  DATATYPE,intent(in) :: onebodymat(nspf,nspf), avectorin(inrnum,numconfig)
  DATATYPE,intent(out) :: avectorout(inrnum,botconfig:topconfig)
  DATAECS :: rvector(inrnum)
  integer ::    config2, config1, iwalk, idiag

  if (topconfig-botconfig+1.eq.0) then
     return
  endif

  avectorout(:,:)=0.d0

  if (diagflag.eq.0) then
     do config1=botconfig,topconfig
        do iwalk=1,numsinglewalks(config1)
           config2=singlewalk(iwalk,config1)
           avectorout(:,config1)=avectorout(:,config1)+avectorin(:,config2)*&
                onebodymat(singlewalkopspf(1,iwalk,config1), &
                singlewalkopspf(2,iwalk,config1)) *  &
                singlewalkdirphase(iwalk,config1) * rvector(:)
        enddo
     enddo
  else
     do config1=botconfig,topconfig
        do idiag=1,numsinglediagwalks(config1)
           iwalk=singlediag(idiag,config1)
           avectorout(:,config1)=avectorout(:,config1)+avectorin(:,config1)*&
                onebodymat(singlewalkopspf(1,iwalk,config1), &
                singlewalkopspf(2,iwalk,config1)) *  &
                singlewalkdirphase(iwalk,config1) * rvector(:)
        enddo
     enddo
  endif
end subroutine arbitraryconfig_mult_singles_nompi



recursive subroutine arbitraryconfig_mult_doubles_nompi(twobodymat, rvector, avectorin, avectorout,inrnum,diagflag)   
  use walkmod
  use parameters
  implicit none
  integer,intent(in) :: inrnum,diagflag
  DATATYPE,intent(in) :: twobodymat(nspf,nspf,nspf,nspf), avectorin(inrnum,numconfig)
  DATATYPE,intent(out) :: avectorout(inrnum,botconfig:topconfig)
  DATAECS :: rvector(inrnum)
  integer ::   config2, config1, iwalk, idiag

  avectorout(:,:)=0.d0

  if (diagflag.eq.0) then
     do config1=botconfig,topconfig

        do iwalk=1,numdoublewalks(config1)
           config2=doublewalk(iwalk,config1)
           avectorout(:,config1)=avectorout(:,config1)+avectorin(:,config2)*&
                twobodymat(doublewalkdirspf(1,iwalk,config1), &
                doublewalkdirspf(2,iwalk,config1),   &
                doublewalkdirspf(3,iwalk,config1),   &
                doublewalkdirspf(4,iwalk,config1))* &
                doublewalkdirphase(iwalk,config1) * rvector(:)
        enddo
     enddo
  else
     do config1=botconfig,topconfig
        do idiag=1,numdoublediagwalks(config1)
           iwalk=doublediag(idiag,config1)
           config2=doublewalk(iwalk,config1)
           avectorout(:,config1)=avectorout(:,config1)+avectorin(:,config2)*&
                twobodymat(doublewalkdirspf(1,iwalk,config1), &
                doublewalkdirspf(2,iwalk,config1),   &
                doublewalkdirspf(3,iwalk,config1),   &
                doublewalkdirspf(4,iwalk,config1))* &
                doublewalkdirphase(iwalk,config1) * rvector(:)
        enddo
     enddo
  endif

end subroutine arbitraryconfig_mult_doubles_nompi
