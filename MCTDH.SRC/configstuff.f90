
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
        OFLWR "Construct big matrix:  ", numconfig*numr; CFL
     endif
     allocate(fullconfigvals(numbasis*numr))
     allocate(fullconfigmatel(numbasis*numr,numbasis*numr), fullconfigvects(numconfig*numr,numbasis*numr))
     allocate(tempconfigvects(numbasis*numr,numbasis*numr))
     fullconfigmatel=0.d0; fullconfigvects=0d0; fullconfigvals=0d0; tempconfigvects=0d0

     if (allspinproject.ne.0) then
        call assemble_spinconfigmat(fullconfigmatel,yyy%cptr(0),1,1,0,0,  time)
     else
        call assemble_configmat(fullconfigmatel,yyy%cptr(0),1,1,0,0, time)
     endif
     if (printflag.ne.0) then
        OFLWR " Call eig ",numbasis*numr; CFL
     endif
     if (myrank.eq.1) then
        call CONFIGEIG(fullconfigmatel(:,:),numbasis*numr,numbasis*numr, tempconfigvects, fullconfigvals)
        if (allspinproject.eq.0) then
           fullconfigvects(:,:)=tempconfigvects(:,:)
        else
           do i=1,numspinconfig*numr
              call configspin_transformfrom(numr,tempconfigvects(:,i),fullconfigvects(:,i))
           enddo
        endif
     endif

     call mympibcast(fullconfigvects,1,numbasis*numr*numconfig*numr)
     call mympibcast(fullconfigvals,1,numbasis*numr)

     if (printflag.ne.0) then
        OFLWR "  -- Nonsparse eigenvals --"
        do i=1,min(numbasis*numr,numshift+order+10)
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
  if (allspinproject.ne.0) then
     call configspin_project(avectorin,0)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(avectorin,numr)
  endif
  if (sparseconfigflag/=0) then
     call exposparseprop(avectorin,avectorout,time)
  else
     call nonsparseprop(avectorin,avectorout,time)
  endif
  if (allspinproject.ne.0) then
     call configspin_project(avectorout,0)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(avectorout,numr)
  endif
end subroutine myconfigprop


subroutine nonsparseprop(avectorin,avectorout,time)
  use parameters
  use mpimod
  use configpropmod
  implicit none

  DATATYPE :: avectorin(totadim), avectorout(totadim),  avectortemp(numbasis*numr), avectortemp2(numbasis*numr)
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
  allocate(iiwork(numbasis*numr*2))

  allocate(bigconfigmatel(numbasis*numr,numbasis*numr), bigconfigvects(numbasis*numr,2*(numbasis*numr+2)))

  if (allspinproject.eq.0) then
     call assemble_configmat(bigconfigmatel, workconfigpointer,1,1,1,1, time)
  else
     call assemble_spinconfigmat(bigconfigmatel, workconfigpointer,1,1,1,1, time)
  endif
  
  bigconfigmatel=bigconfigmatel*timefac

  if (allspinproject.eq.0) then
     avectortemp(:)=avectorin(:)
  else
     call configspin_transformto(numr,avectorin,avectortemp)
  endif

  if (myrank.eq.1) then
     if (nonsparsepropmode.eq.1) then
        avectortemp2(:)=avectortemp(:)
        call expmat(bigconfigmatel,numbasis*numr) 
        call MYGEMV('N',numbasis*numr,numbasis*numr,DATAONE,bigconfigmatel,numbasis*numr,avectortemp2,1,DATAZERO,avectortemp,1)
#ifndef REALGO
     else if (nonsparsepropmode.eq.2) then
        allocate(realbigconfigmatel(2,numbasis*numr,2,numbasis*numr))
        call assigncomplexmat(realbigconfigmatel,bigconfigmatel,numbasis*numr,numbasis*numr)
        call DGCHBV(numbasis*numr*2, 1.d0, realbigconfigmatel, numbasis*numr*2, avectortemp, bigconfigvects, iiwork, iflag)           
        deallocate(realbigconfigmatel)
#endif
     else
        call EXPFULL(numbasis*numr, 1.d0, bigconfigmatel, numbasis*numr, avectortemp, bigconfigvects, iiwork, iflag)
        
        if (iflag.ne.0) then
           OFLWR "Expo error A ", iflag; CFLST
        endif
     endif
  endif

  call mympibcast(avectortemp,1,numbasis*numr)

  if (allspinproject.eq.0) then
     avectorout(:)=avectortemp(:)
  else
     call configspin_transformfrom(numr,avectortemp,avectorout)
  endif

  deallocate(iiwork,bigconfigmatel,bigconfigvects)

end subroutine nonsparseprop



!! SUBROUTINE FOR LANCZOS PROPAGATION (CALLED BY EXPOKIT)

subroutine configexpotimeinit(intime)
  use parameters
  use configexpotimemod
  implicit none
  real*8 :: intime
  configexpotime=intime;  initialized=.true.
end subroutine configexpotimeinit




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





