
#include "Definitions.INC"

!! All purpose subroutine.

recursive subroutine sparseconfigmult(myinvector,myoutvector,matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time)
  use parameters
  use configptrmod
  use sparseptrmod
  implicit none
  DATATYPE :: myinvector(numconfig,numr), myoutvector(numconfig,numr)
  Type(CONFIGPTR) :: matrix_ptr
  Type(SPARSEPTR) :: sparse_ptr
  integer :: conflag,boflag,nucflag,pulseflag
  real*8 :: time

  call sparseconfigmultxxx(myinvector,myoutvector,matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,0)

end subroutine sparseconfigmult


!! For one bond length, born oppenheimer Hamiltonian

recursive subroutine sparseconfigmultone(invector,outvector,matrix_ptr,sparse_ptr, boflag,pulseflag,conflag,isplit,time)
  use parameters
  use configptrmod
  use sparseptrmod
  implicit none
  DATATYPE :: invector(numconfig), outvector(numconfig)
  Type(CONFIGPTR) :: matrix_ptr
  Type(SPARSEPTR) :: sparse_ptr
  integer :: conflag, boflag, pulseflag, isplit
  real*8 :: time

  call sparseconfigmultnew_transpose_nompi(invector,outvector(botwalk),matrix_ptr,sparse_ptr, boflag, 0, pulseflag, conflag,time,0,isplit,isplit,0)
  if (sparseconfigflag.ne.0) then
     call mpiallgather(outvector,numconfig,configsperproc,maxconfigsperproc)
  endif

end subroutine sparseconfigmultone


!! All purpose subroutine, just the pulse.

recursive subroutine sparseconfigpulsemult(myinvector,myoutvector,matrix_ptr, sparse_ptr,which)
  use parameters
  use configptrmod
  use sparseptrmod
  implicit none
  DATATYPE :: myinvector(numconfig,numr), myoutvector(numconfig,numr)
  Type(CONFIGPTR) :: matrix_ptr
  Type(SPARSEPTR) :: sparse_ptr
  integer :: which

  call sparseconfigmultxxx(myinvector,myoutvector,matrix_ptr,sparse_ptr, 0,0,1,0,-1d0,which)

end subroutine sparseconfigpulsemult


!! MPI subroutine

recursive subroutine sparseconfigmultxxx(myinvector,myoutvector,matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag)
  use parameters
  use configptrmod
  use sparseptrmod
  implicit none
  DATATYPE :: myinvector(numconfig,numr), myoutvector(numconfig,numr),&
       myinvectortr(numr,numconfig),myoutvectortr(numr,numconfig)  !!AUTOMATIC
  Type(CONFIGPTR) :: matrix_ptr
  Type(SPARSEPTR) :: sparse_ptr
  integer :: conflag,boflag,nucflag,pulseflag,onlytdflag
  real*8 :: time

  myinvectortr(:,:)=TRANSPOSE(myinvector(:,:)); myoutvectortr(:,:)=0d0
  call sparseconfigmultnew_transpose_nompi(myinvectortr,myoutvectortr(:,botwalk),matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag,1,numr,0)
  if (sparseconfigflag.ne.0) then
     call mpiallgather(myoutvectortr,numconfig*numr,configsperproc*numr,maxconfigsperproc*numr)
  endif
  myoutvector(:,:)=TRANSPOSE(myoutvectortr)
     
end subroutine sparseconfigmultxxx


!! So, counterintuitively, this routine is used within MPI propagation with
!!   parconfigsplit=1. 
!!   
!! Is a subroutine without mpi communication, that can be called, if vector is distributed
!!   (parconfigsplit=1)
!!   after communication to gather invector has been done.  Input whole vector; output only
!!   that part for which I have rows in the Hamiltonian.


recursive subroutine sparseconfigmult_nompi(myinvector,myoutvector,matrix_ptr, sparse_ptr,boflag, nucflag, pulseflag, conflag,time)
  use parameters
  use configptrmod
  use sparseptrmod
  implicit none
  DATATYPE :: myinvector(numconfig,numr), myoutvector(botwalk:topwalk,numr), &
       myinvectortr(numr,numconfig), myoutvectortr(numr,botwalk:topwalk)  !!AUTOMATIC
  real*8 :: time
  Type(CONFIGPTR) :: matrix_ptr
  Type(SPARSEPTR) :: sparse_ptr
  integer ::  conflag,boflag,nucflag,pulseflag

  if (topwalk-botwalk+1.eq.0) then
     return
  endif
  myinvectortr(:,:)=TRANSPOSE(myinvector(:,:)); myoutvectortr(:,:)=0d0
  call sparseconfigmultnew_transpose_nompi(myinvectortr,myoutvectortr,matrix_ptr, sparse_ptr, boflag, nucflag, pulseflag, conflag,time,0,1,numr,0)
  myoutvector(:,:)=TRANSPOSE(myoutvectortr(:,:))

end subroutine sparseconfigmult_nompi


recursive subroutine sparseconfigdiagmult(invector,outvector,matrix_ptr, sparse_ptr, boflag, nucflag, pulseflag, conflag, inenergy,time)
  use parameters
  use configptrmod
  use sparseptrmod
  implicit none
  DATATYPE :: invector(numconfig,numr), outvector(numconfig,numr), inenergy, &
       myinvectortr(numr,numconfig),myoutvectortr(numr,numconfig),facs(3)
  Type(CONFIGPTR) :: matrix_ptr
  Type(SPARSEPTR) :: sparse_ptr
  integer ::  conflag,boflag,nucflag,pulseflag
  real*8 :: time

  call vectdpot(time,velflag,facs)

  myinvectortr(:,:)=1d0; myoutvectortr(:,:)=0d0
  call sparseconfigmultnew_transpose_nompi(myinvectortr,myoutvectortr(:,botwalk),matrix_ptr, sparse_ptr, boflag, nucflag, pulseflag, conflag,time,0,1,numr,1)
  if (sparseconfigflag.ne.0) then
     call mpiallgather(myoutvectortr,numconfig*numr,configsperproc*numr,maxconfigsperproc*numr)
  endif
  outvector(:,:)=TRANSPOSE(myoutvectortr)

  if (quadpreconshift.eq.0d0) then
     outvector(:,:)=invector(:,:)/(outvector(:,:)-inenergy)
  else
     outvector(:,:)=invector(:,:)*0.5d0* (1d0/(outvector(:,:)-inenergy+(0d0,1d0)*quadpreconshift) + 1d0/(outvector(:,:)-inenergy-(0d0,1d0)*quadpreconshift))
  endif

end subroutine sparseconfigdiagmult


recursive subroutine parsparseconfigdiagmult_transpose(myinvectortr,myoutvectortr,matrix_ptr, sparse_ptr, boflag, nucflag, pulseflag, conflag, inenergy,time)
  use parameters
  use configptrmod
  use sparseptrmod
  implicit none
  DATATYPE :: inenergy,tempvectortr(numr,numconfig), &  !!AUTOMATIC
       myinvectortr(numr,numconfig),myoutvectortr(numr,botwalk:topwalk)
  Type(CONFIGPTR) :: matrix_ptr
  Type(SPARSEPTR) :: sparse_ptr
  integer ::  conflag,boflag,nucflag,pulseflag
  real*8 :: time

  if (topwalk-botwalk+1.eq.0) then
     return
  endif

  tempvectortr(:,:)=1d0; myoutvectortr(:,:)=0d0
  call sparseconfigmultnew_transpose_nompi(tempvectortr,myoutvectortr,matrix_ptr, sparse_ptr, boflag, nucflag, pulseflag, conflag,time,0,1,numr,1)
  
  if (quadpreconshift.eq.0d0) then
     myoutvectortr(:,:)=myinvectortr(:,:)/(myoutvectortr(:,:)-inenergy)
  else
     myoutvectortr(:,:)=myinvectortr(:,:)*0.5d0* (1d0/(myoutvectortr(:,:)-inenergy+(0d0,1d0)*quadpreconshift) + 1d0/(myoutvectortr(:,:)-inenergy-(0d0,1d0)*quadpreconshift))
  endif
end subroutine parsparseconfigdiagmult_transpose



!! NOW WRAPPER - SPARSEOPT HERE

recursive subroutine sparseconfigmultnew_transpose_nompi(myinvectortr,myoutvectortr,matrix_ptr,sparse_ptr,&
     boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag)
  use parameters
  use configptrmod
  use sparseptrmod
  implicit none
  integer :: conflag,boflag,nucflag,pulseflag,onlytdflag,botr,topr,diagflag
  real*8 :: time
  DATATYPE :: myinvectortr(botr:topr,numconfig), myoutvectortr(botr:topr,botwalk:topwalk)
  Type(CONFIGPTR) :: matrix_ptr
  Type(SPARSEPTR) :: sparse_ptr

  if (topwalk-botwalk+1.eq.0) then
     return
  endif
  if (sparseopt.eq.0) then
     call direct_sparseconfigmultnew_transpose_nompi(myinvectortr,myoutvectortr,matrix_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag)
  else
     call sparsesparsemult_transpose_nompi(myinvectortr,myoutvectortr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag)
  endif

end subroutine sparseconfigmultnew_transpose_nompi


!! SPARSE MATVEC (with matrix_ptr not sparse_ptr)

recursive subroutine direct_sparseconfigmultnew_transpose_nompi(myinvectortr,myoutvectortr,matrix_ptr, &
     boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag)
  use parameters
  use configptrmod
  use opmod   !! rkemod, proderivmod  
  implicit none

  integer :: conflag,boflag,nucflag,pulseflag,onlytdflag,botr,topr,mynumr,diagflag,ir,flag,ii
  DATATYPE :: myinvectortr(botr:topr,numconfig), myoutvectortr(botr:topr,botwalk:topwalk),tempvectortr(botr:topr,botwalk:topwalk), &
       tempmatel(nspf,nspf),facs(3),csum
  DATAECS :: rvector(botr:topr)
  Type(CONFIGPTR) :: matrix_ptr
  real*8 :: time,gg

  if (topwalk-botwalk+1.eq.0) then
     return
  endif

  myoutvectortr(:,:)=0d0

  if (onlytdflag.ne.0) then
     facs(:)=0;     facs(onlytdflag)=1d0
  else
     call vectdpot(time,velflag,facs)
  endif
  
  mynumr=topr-botr+1
  
  if (boflag==1) then
     
!! easy nuclear repulsion, hardwired like this for now

     do ir=botr,topr

        myoutvectortr(ir,:)=myoutvectortr(ir,:) + myinvectortr(ir,botwalk:topwalk) * ( &
             energyshift + (nucrepulsion+frozenpotdiag)/bondpoints(ir) + &
             frozenkediag/bondpoints(ir)**2 ) * matrix_ptr%kefac

     enddo
     
     rvector(:)=1/bondpoints(botr:topr)
     call arbitraryconfig_mult_doubles_transpose_nompi(matrix_ptr%xtwoematel(:,:,:,:),rvector,myinvectortr,tempvectortr,mynumr,diagflag)
     myoutvectortr(:,:)=myoutvectortr(:,:)+tempvectortr(:,:)
     
     call arbitraryconfig_mult_transpose_nompi(matrix_ptr%xpotmatel(:,:),rvector,myinvectortr,tempvectortr,mynumr,diagflag)
     myoutvectortr(:,:)=myoutvectortr(:,:)+tempvectortr(:,:)
     
     rvector(:)=1/bondpoints(botr:topr)**2
     call arbitraryconfig_mult_transpose_nompi(matrix_ptr%xopmatel(:,:),rvector,myinvectortr,tempvectortr,mynumr,diagflag)
     myoutvectortr(:,:)=myoutvectortr(:,:)+tempvectortr(:,:)
     
  endif
  
  if (conflag==1.and.constraintflag.ne.0) then
     rvector(:)=1
     
     tempmatel(:,:)=matrix_ptr%xconmatel(:,:)
     if (tdflag.ne.0) then
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelxx(:,:)*facs(1)
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelyy(:,:)*facs(2)
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelzz(:,:)*facs(3)
     endif
     call arbitraryconfig_mult_transpose_nompi(tempmatel,rvector,myinvectortr,tempvectortr,mynumr,diagflag)
     myoutvectortr(:,:)=myoutvectortr(:,:)-tempvectortr(:,:)

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
        call arbitraryconfig_mult_transpose_nompi(tempmatel,rvector,myinvectortr,tempvectortr,mynumr,diagflag)
        myoutvectortr(:,:)=myoutvectortr(:,:)+tempvectortr(:,:)
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
           myoutvectortr(ir,:)=myoutvectortr(ir,:)+myinvectortr(ir,botwalk:topwalk)*csum*bondpoints(ir)
        enddo
     else if (onlytdflag.eq.0) then
        csum=matrix_ptr%kefac * numelec * (facs(1)**2 + facs(2)**2 + facs(3)**2 ) * 2  * gg   !! a-squared
        do ir=botr,topr
           myoutvectortr(ir,:)=myoutvectortr(ir,:)+myinvectortr(ir,botwalk:topwalk)*csum     !! NO R FACTOR !!
        enddo
     endif

  endif   !! PULSE

  if (nonuc_checkflag.eq.0.and.nucflag.eq.1.and.mynumr.eq.numr) then
     if (diagflag.eq.0) then
        rvector(:)=1d0
        call arbitraryconfig_mult_transpose_nompi(matrix_ptr%xymatel,rvector,myinvectortr,tempvectortr,numr,0)
        call MYGEMM('N','N',numr,topwalk-botwalk+1,numr,DATAONE,            proderivmod,numr,     tempvectortr,           numr,DATAONE,myoutvectortr,numr)
        call MYGEMM('N','N',numr,topwalk-botwalk+1,numr,DATAONE*matrix_ptr%kefac,rkemod,numr,     myinvectortr(:,botwalk),numr,DATAONE,myoutvectortr,numr)
!!$     else
!!$        do ir=1,numr
!!$           myoutvectortr(ir,:)=myoutvectortr(ir,:)+myinvectortr(ir,botwalk:topwalk)*matrix_ptr%kefac*rkemod(ir,ir)
!!$        enddo
     endif
  endif

end subroutine direct_sparseconfigmultnew_transpose_nompi


!! DIRECT MATVEC (with sparse_ptr not matrix_ptr)

recursive subroutine sparsesparsemult_transpose_nompi(myinvectortr,myoutvectortr,sparse_ptr,&
     boflag,nucflag,pulseflag,conflag,time,onlytdflag,botr,topr,diagflag)
  use parameters
  use sparseptrmod
  use walkmod
  use opmod   !! rkemod, proderivmod
  implicit none
  integer :: ir,botr,topr,mynumr,ii
  integer :: boflag,nucflag,pulseflag,conflag,onlytdflag,diagflag,flag
  real*8 :: time,gg
  DATAECS :: rvector(botr:topr)
  Type(SPARSEPTR) :: sparse_ptr
  DATATYPE ::        myinvectortr(botr:topr,numconfig),myoutvectortr(botr:topr,botwalk:topwalk), facs(3),csum
  DATATYPE :: tempvectortr(botr:topr,botwalk:topwalk), tempsparsemattr(maxsinglewalks,botwalk:topwalk) !! AUTOMATIC

!!$  integer, save:: allochere=0
!!$  DATATYPE,save,allocatable :: tempvectortr(:,:), tempsparsemattr(:,:)
!!$  if (allochere.eq.0) then
!!$     allocate(tempvectortr(botr:topr,botwalk:topwalk), tempsparsemattr(maxsinglewalks,botwalk:topwalk))
!!$  endif
!!$  allochere=1

  if (topwalk-botwalk+1.eq.0) then
     return
  endif

  if (botr.ne.1.or.topr.ne.numr) then
     OFLWR "REPROGRAMMMM."; CFLST
  endif
  if (sparseconfigflag.eq.0) then
     OFLWR "BADDDCALLL6666666"; CFLST
  endif

  myoutvectortr(:,:)=0d0

  if (onlytdflag.ne.0) then
     facs(:)=0;     facs(onlytdflag)=1d0
  else
     call vectdpot(time,velflag,facs)
  endif

  mynumr=topr-botr+1

  if (boflag==1) then

    do ir=botr,topr

       myoutvectortr(ir,:)=myoutvectortr(ir,:) + myinvectortr(ir,botwalk:topwalk) * ( &
            energyshift + (nucrepulsion+frozenpotdiag)/bondpoints(ir) + &
            frozenkediag/bondpoints(ir)**2 ) * sparse_ptr%kefac
       
    enddo

     rvector(:)=1/bondpoints(botr:topr)
     call arbitrary_sparsemult_transpose_nompi_doubles(sparse_ptr%xpotsparsemattr(:,:),rvector,myinvectortr,tempvectortr,mynumr,diagflag)
     myoutvectortr(:,:)=myoutvectortr(:,:)+tempvectortr(:,:)

     call arbitrary_sparsemult_transpose_nompi_singles(sparse_ptr%xonepotsparsemattr(:,:),rvector,myinvectortr,tempvectortr,mynumr,diagflag)
     myoutvectortr(:,:)=myoutvectortr(:,:)+tempvectortr(:,:)

     rvector(:)=1/bondpoints(botr:topr)**2
     call arbitrary_sparsemult_transpose_nompi_singles(sparse_ptr%xopsparsemattr(:,:),rvector,myinvectortr,tempvectortr,mynumr,diagflag)
     myoutvectortr(:,:)=myoutvectortr(:,:)+tempvectortr(:,:)

  endif

  if (conflag==1.and.constraintflag.ne.0) then
     rvector(:)=1

     tempsparsemattr(:,:)=sparse_ptr%xconsparsemattr(:,:)
     if (tdflag.ne.0) then
        tempsparsemattr(:,:)=tempsparsemattr(:,:)+sparse_ptr%xconsparsemattrxx(:,:)*facs(1)
        tempsparsemattr(:,:)=tempsparsemattr(:,:)+sparse_ptr%xconsparsemattryy(:,:)*facs(2)
        tempsparsemattr(:,:)=tempsparsemattr(:,:)+sparse_ptr%xconsparsemattrzz(:,:)*facs(3)
     endif
     call arbitrary_sparsemult_transpose_nompi_singles(tempsparsemattr,rvector,myinvectortr,tempvectortr,mynumr,diagflag)
     myoutvectortr(:,:)=myoutvectortr(:,:)-tempvectortr(:,:)

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
        call arbitrary_sparsemult_transpose_nompi_singles(tempsparsemattr,rvector,myinvectortr,tempvectortr,mynumr,diagflag)
        myoutvectortr(:,:)=myoutvectortr(:,:)+tempvectortr(:,:)
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
           myoutvectortr(ir,:)=myoutvectortr(ir,:)+myinvectortr(ir,botwalk:topwalk)*csum*bondpoints(ir)
        enddo
     else if (onlytdflag.eq.0) then
        csum= sparse_ptr%kefac * numelec * (facs(1)**2 + facs(2)**2 + facs(3)**2 ) * 2 * gg  !! a-squared
        do ir=botr,topr
           myoutvectortr(ir,:)=myoutvectortr(ir,:)+myinvectortr(ir,botwalk:topwalk)*csum     !! NO R FACTOR !!
        enddo
     endif

  endif !! PULSE

  if (nonuc_checkflag.eq.0.and.nucflag.eq.1.and.mynumr.eq.numr) then
     if (diagflag.eq.0) then
        rvector(:)=1d0
        call arbitrary_sparsemult_transpose_nompi_singles(sparse_ptr%xysparsemattr,rvector,myinvectortr,tempvectortr,numr,0)
        call MYGEMM('N','N',numr,topwalk-botwalk+1,numr,DATAONE,            proderivmod,numr,     tempvectortr,           numr,DATAONE,myoutvectortr,numr)
        call MYGEMM('N','N',numr,topwalk-botwalk+1,numr,DATAONE*sparse_ptr%kefac,rkemod,numr,     myinvectortr(:,botwalk),numr,DATAONE,myoutvectortr,numr)
!!$     else
!!$        do ir=1,numr
!!$           myoutvectortr(ir,:)=myoutvectortr(ir,:)+myinvectortr(ir,botwalk:topwalk)*sparse_ptr%kefac*rkemod(ir,ir)
!!$        enddo
     endif
  endif

end subroutine sparsesparsemult_transpose_nompi



recursive subroutine arbitrary_sparsemult_transpose_nompi_singles(mattrans, rvector,inbigvectortr,outsmallvectortr,mynumr,diagflag)
  use parameters
  use walkmod
  implicit none
  integer :: iwalk,config1,mynumr,diagflag,idiag
  DATATYPE ::  inbigvectortr(mynumr,numconfig),outsmallvectortr(mynumr,botwalk:topwalk), mattrans(maxsinglewalks,botwalk:topwalk)
  DATAECS :: rvector(mynumr)

  if (sparseconfigflag.eq.0) then
     OFLWR "BADDDCALLL6666666"; CFLST
  endif

  if (topwalk-botwalk+1.eq.0) then
     return
  endif

  outsmallvectortr(:,:)=0d0

  if (diagflag.eq.0) then
     do config1=botwalk,topwalk
        do iwalk=1,numsinglewalks(config1)
           outsmallvectortr(:,config1)=outsmallvectortr(:,config1)+mattrans(iwalk,config1) * inbigvectortr(:,singlewalk(iwalk,config1))
        enddo
        outsmallvectortr(:,config1)=outsmallvectortr(:,config1)*rvector(:)
     enddo
  else
     if (sortwalks.ne.0) then
        OFLWR " ??? NOT DEBUGGED ??? sortwalks singlediag.  programmer remove me to try"; CFLST
     endif
     do config1=botwalk,topwalk
        do idiag=1,numsinglediagwalks(config1)
           outsmallvectortr(:,config1)=outsmallvectortr(:,config1)+mattrans(singlediag(idiag,config1),config1) * inbigvectortr(:,config1)
        enddo
        outsmallvectortr(:,config1)=outsmallvectortr(:,config1)*rvector(:)
     enddo
  endif

end subroutine arbitrary_sparsemult_transpose_nompi_singles


recursive subroutine arbitrary_sparsemult_transpose_nompi_doubles(mattrans,rvector,inbigvectortr,outsmallvectortr,mynumr,diagflag)
  use parameters
  use walkmod
  implicit none
  integer :: iwalk,config1,mynumr,diagflag,idiag
  DATATYPE :: inbigvectortr(mynumr,numconfig),outsmallvectortr(mynumr,botwalk:topwalk),  mattrans(maxdoublewalks,botwalk:topwalk)
  DATAECS :: rvector(mynumr)

  if (sparseconfigflag.eq.0) then
     OFLWR "BADDDCALLL6666666"; CFLST
  endif
  if (topwalk-botwalk+1.eq.0) then
     return
  endif

  outsmallvectortr(:,:)=0d0
  do config1=botwalk,topwalk
     if (diagflag.eq.0) then
        do iwalk=1,numdoublewalks(config1)
           outsmallvectortr(:,config1)=outsmallvectortr(:,config1)+mattrans(iwalk,config1) * inbigvectortr(:,doublewalk(iwalk,config1))
        enddo
     else
        if (sortwalks.ne.0) then
           OFLWR " ??? NOT DEBUGGED ??? sortwalks doublediag.  programmer remove me to try"; CFLST
        endif
        do idiag=1,numdoublediagwalks(config1)
           outsmallvectortr(:,config1)=outsmallvectortr(:,config1)+mattrans(doublediag(idiag,config1),config1) * inbigvectortr(:,config1)
        enddo
     endif
     outsmallvectortr(:,config1)=outsmallvectortr(:,config1)*rvector(:)
  enddo

end subroutine arbitrary_sparsemult_transpose_nompi_doubles


!! TAKES TRANSPOSED VECTOR (numr,numconfig).  temptransposeflag irrelevant because not constructing configmatels

recursive subroutine arbitraryconfig_mult(onebodymat, rvector, avectorin, avectorout,inrnum)   
  use parameters
  implicit none
  integer :: inrnum
  DATATYPE :: onebodymat(nspf,nspf), avectorin(numconfig,inrnum), avectorout(numconfig,inrnum), &
       avectorintr(inrnum,numconfig),avectorouttr(inrnum,numconfig)
  DATAECS :: rvector(inrnum)

  avectorintr(:,:)=TRANSPOSE(avectorin(:,:)); avectorouttr(:,:)=0d0
  call arbitraryconfig_mult_transpose_nompi(onebodymat, rvector, avectorintr, avectorouttr(:,botwalk),inrnum,0)
  if (sparseconfigflag.ne.0) then
     call mpiallgather(avectorouttr,numconfig*inrnum,configsperproc*inrnum,maxconfigsperproc*inrnum)
  endif
  avectorout(:,:)=TRANSPOSE(avectorouttr)

end subroutine arbitraryconfig_mult


recursive subroutine arbitraryconfig_mult_transpose_nompi(onebodymat, rvector, avectorintr, avectorouttr,inrnum,diagflag)
   use walkmod
  use parameters
  implicit none
  integer ::    config2, config1, iwalk, inrnum,diagflag,idiag
  DATATYPE :: onebodymat(nspf,nspf), avectorintr(inrnum,numconfig), avectorouttr(inrnum,botwalk:topwalk)
  DATAECS :: rvector(inrnum)

  if (topwalk-botwalk+1.eq.0) then
     return
  endif

  avectorouttr=0.d0

!! temptransposeflag doesn't matter; not assembling then operating with matrix.
!!       switching previous order for allg ather

  if (diagflag.eq.0) then
     do config1=botwalk,topwalk
        do iwalk=1,numsinglewalks(config1)
           config2=singlewalk(iwalk,config1)
           avectorouttr(:,config1)=avectorouttr(:,config1)+avectorintr(:,config2)*&
                onebodymat(singlewalkopspf(1,iwalk,config1), &
                singlewalkopspf(2,iwalk,config1)) *  &
                singlewalkdirphase(iwalk,config1) * rvector(:)
        enddo
     enddo
  else
     if (sortwalks.ne.0) then
        OFLWR " ??? NOT DEBUGGED ??? sortwalks singlediag.  programmer remove me to try"; CFLST
     endif
     do config1=botwalk,topwalk
        do idiag=1,numsinglediagwalks(config1)
           iwalk=singlediag(idiag,config1)
           avectorouttr(:,config1)=avectorouttr(:,config1)+avectorintr(:,config1)*&
                onebodymat(singlewalkopspf(1,iwalk,config1), &
                singlewalkopspf(2,iwalk,config1)) *  &
                singlewalkdirphase(iwalk,config1) * rvector(:)
        enddo
     enddo
  endif
end subroutine arbitraryconfig_mult_transpose_nompi



recursive subroutine arbitraryconfig_mult_doubles_transpose_nompi(twobodymat, rvector, avectorintr, avectorouttr,inrnum,diagflag)   
  use walkmod
  use parameters
  implicit none
  integer ::   config2, config1, iwalk, inrnum,diagflag, idiag
  DATATYPE :: twobodymat(nspf,nspf,nspf,nspf), avectorintr(inrnum,numconfig), avectorouttr(inrnum,botwalk:topwalk)
  DATAECS :: rvector(inrnum)

  avectorouttr=0.d0

  if (diagflag.eq.0) then
     do config1=botwalk,topwalk

        do iwalk=1,numdoublewalks(config1)
           config2=doublewalk(iwalk,config1)
           avectorouttr(:,config1)=avectorouttr(:,config1)+avectorintr(:,config2)*&
                twobodymat(doublewalkdirspf(1,iwalk,config1), &
                doublewalkdirspf(2,iwalk,config1),   &
                doublewalkdirspf(3,iwalk,config1),   &
                doublewalkdirspf(4,iwalk,config1))* &
                doublewalkdirphase(iwalk,config1) * rvector(:)
        enddo
     enddo
  else
     if (sortwalks.ne.0) then
        OFLWR " ??? NOT DEBUGGED ??? sortwalks doublediag.  programmer remove me to try"; CFLST
     endif
     do config1=botwalk,topwalk
        do idiag=1,numdoublediagwalks(config1)
           iwalk=doublediag(idiag,config1)
           config2=doublewalk(iwalk,config1)
           avectorouttr(:,config1)=avectorouttr(:,config1)+avectorintr(:,config2)*&
                twobodymat(doublewalkdirspf(1,iwalk,config1), &
                doublewalkdirspf(2,iwalk,config1),   &
                doublewalkdirspf(3,iwalk,config1),   &
                doublewalkdirspf(4,iwalk,config1))* &
                doublewalkdirphase(iwalk,config1) * rvector(:)
        enddo
     enddo
  endif

end subroutine arbitraryconfig_mult_doubles_transpose_nompi
