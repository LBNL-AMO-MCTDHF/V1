
#include "Definitions.INC"

subroutine printconfig(thisconfig)
  use parameters
  implicit none
  
  integer :: thisconfig(ndof), i
  character (len=4) :: mslabels(2) =["a ","b "]
  
  write(mpifileptr,'(100(I3,A2))') (thisconfig((i-1)*2+1), mslabels(thisconfig(i*2)), i=1,numelec)

end subroutine printconfig


!! NOT - this uses stopthresh as threshold to pass to laneigen.  is called for mcscf, and improved.


subroutine myconfigeig(thisconfigvects,thisconfigvals,order,printflag, guessflag,time,numshift)
  use parameters
  use mpimod
  use xxxmod
  implicit none

  integer :: order, printflag,i, guessflag,  l,k,numshift,flag
  DATATYPE :: thisconfigvects(totadim,order),lastval,dot
  DATAECS :: thisconfigvals(order), tempconfigvals(order+numshift)
  real*8 :: realconfigvect(totadim)
  real*8 :: sum,time
  DATATYPE, allocatable :: fullconfigmatel(:,:), fullconfigvects(:,:), &
       tempconfigvects(:,:)
  DATAECS, allocatable :: fullconfigvals(:)

  if (order+numshift.gt.totadim) then
     OFLWR "Error, want ",order," plus ",numshift," vectors but totadim= ",totadim;CFLST
  endif
  if (numshift.lt.0) then
     OFLWR "GG ERROR.", numshift,guessflag; CFLST
  endif

  if (sparseconfigflag/=0) then
     
     allocate(tempconfigvects(totadim,order+numshift))
     tempconfigvects(:,:)=0d0
     if (guessflag.ne.0) then
        do i=1,numshift
           call RANDOM_NUMBER(realconfigvect(:))
           tempconfigvects(:,i)=realconfigvect(:)
        enddo
        tempconfigvects(:,numshift+1:order+numshift)=thisconfigvects(:,1:order)
     endif
     call blocklanczos(order+numshift,tempconfigvects,tempconfigvals,printflag,guessflag)
     
        !! so if guessflag=1, numshift is 0.

     thisconfigvals(1:order)=tempconfigvals(numshift+1:numshift+order)
     thisconfigvects(:,1:order)=tempconfigvects(:,numshift+1:numshift+order)
     deallocate(tempconfigvects)

  else

     if (printflag.ne.0) then
        OFLWR "Construct big matrix:  ", numdfbasis*numr; CFL
     endif
     allocate(fullconfigvals(numdfbasis*numr))
     allocate(fullconfigmatel(numdfbasis*numr,numdfbasis*numr), fullconfigvects(numconfig*numr,numdfbasis*numr))
     allocate(tempconfigvects(numdfbasis*numr,numdfbasis*numr))
     fullconfigmatel=0.d0; fullconfigvects=0d0; fullconfigvals=0d0; tempconfigvects=0d0

     call assemble_dfbasismat(fullconfigmatel,yyy%cptr(0),1,1,0,0, time)

     if (printflag.ne.0) then
        OFLWR " Call eig ",numdfbasis*numr; CFL
     endif
     if (myrank.eq.1) then
        call CONFIGEIG(fullconfigmatel(:,:),numdfbasis*numr,numdfbasis*numr, tempconfigvects, fullconfigvals)
        do i=1,numdfbasis*numr
           call basis_transformfrom_all(numr,tempconfigvects(:,i),fullconfigvects(:,i))
        enddo
     endif

     call mympibcast(fullconfigvects,1,numdfbasis*numr*numconfig*numr)
     call mympibcast(fullconfigvals,1,numdfbasis*numr)

     if (printflag.ne.0) then
        OFLWR "  -- Nonsparse eigenvals --"
        do i=1,min(numdfbasis*numr,numshift+order+10)
           WRFL "   ", fullconfigvals(i)
        enddo
        CFL
     endif

     thisconfigvects(:,1:order) = fullconfigvects(:,1+numshift:order+numshift)
     thisconfigvals(1:order) = fullconfigvals(1+numshift:order+numshift)

     lastval = -99999999d0
     flag=0
     do i=1,order
        thisconfigvects(:,i)=thisconfigvects(:,i) / sqrt(dot(thisconfigvects(:,i),thisconfigvects(:,i),totadim))
        if (abs(thisconfigvals(i)-lastval).lt.1d-7) then
           flag=flag+1
           call gramschmidt(totadim,flag,totadim,thisconfigvects(:,i-flag),thisconfigvects(:,i),.false.)
        else
           flag=0
        endif
        lastval=thisconfigvals(i)
     enddo

     deallocate(fullconfigmatel, fullconfigvects, fullconfigvals,tempconfigvects)

  endif

!! this is called from either eigs (for which is irrelevant) or from 
!! config_eigen, which is used to start avectors as improved relax or at
!! start of prop with eigenflag.  Leave as c-normed; in config_eigen, 
!! herm-norm.
!! fix phase for chmctdh/pmctdh debug check

  l=-99
  do i=1,order
     sum=(-999d0)
     do k=1,totadim
        if (abs(thisconfigvects(k,i)).gt.sum) then
           sum=abs(thisconfigvects(k,i));           l=k
        endif
     enddo
     thisconfigvects(:,i)=thisconfigvects(:,i)*abs(thisconfigvects(l,i))/thisconfigvects(l,i)
  enddo

end subroutine myconfigeig



!! PROPAGATE A-VECTOR 

subroutine myconfigprop(avectorin,avectorout,time)
  use parameters
  implicit none
  real*8 :: time
  DATATYPE :: avectorin(totadim), avectorout(totadim)

  if (sparseconfigflag/=0) then
     call exposparseprop(avectorin,avectorout,time)
  else
     call nonsparseprop(avectorin,avectorout,time)
  endif

  if (allspinproject.ne.0) then
     call configspin_project(numr,avectorout)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(numr,avectorout)
  endif

end subroutine myconfigprop


subroutine nonsparseprop(avectorin,avectorout,time)
  use parameters
  use mpimod
  use configpropmod
  implicit none

  DATATYPE :: avectorin(totadim), avectorout(totadim),  avectortemp(numdfbasis*numr), avectortemp2(numdfbasis*numr)
  DATATYPE, allocatable :: bigconfigmatel(:,:), bigconfigvects(:,:) !!,bigconfigvals(:)
#ifndef REALGO
  real*8, allocatable :: realbigconfigmatel(:,:,:,:)
#endif
  integer, allocatable :: iiwork(:)
  integer :: iflag
  real*8 :: time

  if (sparseconfigflag/=0) then
     OFLWR "Error, nonsparseprop called but sparseconfigflag= ", sparseconfigflag; CFLST
  endif
  if (drivingflag.ne.0) then
     OFLWR "Driving flag not implemented for nonsparse"; CFLST
  endif
  allocate(iiwork(numdfbasis*numr*2))

  allocate(bigconfigmatel(numdfbasis*numr,numdfbasis*numr), bigconfigvects(numdfbasis*numr,2*(numdfbasis*numr+2)))

  call assemble_dfbasismat(bigconfigmatel, workconfigpointer,1,1,1,1, time)

  bigconfigmatel=bigconfigmatel*timefac

  call basis_transformto_all(numr,avectorin,avectortemp)

  if (myrank.eq.1) then
     if (nonsparsepropmode.eq.1) then
        avectortemp2(:)=avectortemp(:)
        call expmat(bigconfigmatel,numdfbasis*numr) 
        call MYGEMV('N',numdfbasis*numr,numdfbasis*numr,DATAONE,bigconfigmatel,numdfbasis*numr,avectortemp2,1,DATAZERO,avectortemp,1)
#ifndef REALGO
     else if (nonsparsepropmode.eq.2) then
        allocate(realbigconfigmatel(2,numdfbasis*numr,2,numdfbasis*numr))
        call assigncomplexmat(realbigconfigmatel,bigconfigmatel,numdfbasis*numr,numdfbasis*numr)
        call DGCHBV(numdfbasis*numr*2, 1.d0, realbigconfigmatel, numdfbasis*numr*2, avectortemp, bigconfigvects, iiwork, iflag)           
        deallocate(realbigconfigmatel)
#endif
     else
        call EXPFULL(numdfbasis*numr, 1.d0, bigconfigmatel, numdfbasis*numr, avectortemp, bigconfigvects, iiwork, iflag)
        
        if (iflag.ne.0) then
           OFLWR "Expo error A ", iflag; CFLST
        endif
     endif
  endif

  call mympibcast(avectortemp,1,numdfbasis*numr)

  call basis_transformfrom_all(numr,avectortemp,avectorout)

  deallocate(iiwork,bigconfigmatel,bigconfigvects)

end subroutine nonsparseprop



!! MULTIPLY BY A-VECTOR HAMILTONIAN MATRIX

!! NOTE BOUNDS !!  PADDED

subroutine parconfigexpomult_padded(inavectorspin,outavectorspin)
  use parameters
  use mpimod
  use configexpotimemod
  use configpropmod
  implicit none

  DATATYPE,intent(in) :: inavectorspin(numr,botdfbasis:botdfbasis+maxdfbasisperproc-1)
  DATATYPE,intent(out) :: outavectorspin(numr,botdfbasis:botdfbasis+maxdfbasisperproc-1)
  DATATYPE :: intemp(numr,numconfig), outavector(numr,botconfig:topconfig)

  call avectortime(3)

  if (sparseconfigflag.eq.0) then
     OFLWR "error, must use sparse for parconfigexpomult"; CFLST
  endif

  intemp(:,:)=0d0;  

  call basis_transformfrom_local(numr,inavectorspin,intemp(:,botconfig))

!! DO SUMMA  
  call mpiallgather(intemp,numconfig*numr,configsperproc*numr,maxconfigsperproc*numr)

  call sparseconfigmult_nompi(intemp(:,:),outavector(:,:), workconfigpointer, worksparsepointer, 1,1,1,1,configexpotime,0,1,numr,0)

  outavectorspin(:,:)=0d0

  call basis_transformto_local(numr,outavector,outavectorspin)

  outavectorspin=outavectorspin*timefac
  
  call avectortime(2)

end subroutine parconfigexpomult_padded





