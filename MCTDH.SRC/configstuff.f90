
#include "Definitions.INC"

subroutine printconfig(thisconfig)
  use parameters
  implicit none
  
  integer :: thisconfig(ndof), i
  character (len=4) :: mslabels(2) =["a ","b "]
  
!!  call openfile()

!  write(mpifileptr,'(100(I5,A3,A3,I3))') (thisconfig((i-1)*2+1), mslabels(thisconfig(i*2)),"M=",spfmvals(thisconfig((i-1)*2+1)), i=1,numelec)

  write(mpifileptr,'(100(I3,A2))') (thisconfig((i-1)*2+1), mslabels(thisconfig(i*2)), i=1,numelec)

!!  call closefile()
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
       spinconfigmatel(:,:), spinconfigvects(:,:), tempconfigvects(:,:)
  DATAECS, allocatable :: fullconfigvals(:)

!  if (numshift.lt.0.or.numshift.ne.0.and.guessflag.ne.0) then
!     OFLWR "GG ERROR.", numshift,guessflag; CFLST
!  endif

  if (order.gt.totadim) then
     OFLWR "Error, want ",order," vectors but totadim= ",totadim;CFLST
  endif
  if (numshift.lt.0) then
     OFLWR "GG ERROR.", numshift,guessflag; CFLST
  endif

  if (sparseconfigflag/=0) then
     
!!     OFLWR "    ... go myconfigeig, sparse"; CFL

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
     allocate(fullconfigvals(totadim))
     allocate(fullconfigmatel(numconfig*numr,numconfig*numr), fullconfigvects(totadim,totadim))
     fullconfigmatel=0.d0; fullconfigvects=0d0; fullconfigvals=0d0

     if (spineigflag.ne.0.and.allspinproject.ne.0) then
        allocate(spinconfigmatel(spintotrank*numr,spintotrank*numr), spinconfigvects(spintotrank*numr,spintotrank*numr))
        spinconfigmatel=0.d0

        call assemble_spinconfigmat(spinconfigmatel,yyy%cptr(0),1,1,0,0,  time)
        if (printflag.ne.0) then
           OFLWR " Call spin eig ",spintotrank*numr; CFL
        endif
        if (myrank.eq.1) then
           call CONFIGEIG(spinconfigmatel(:,:),spintotrank*numr,spintotrank*numr, spinconfigvects, fullconfigvals)
           call configspin_transformfrom(spintotrank*numr*numr,spinconfigvects(:,:),fullconfigvects(:,:))
        endif
        deallocate(spinconfigmatel,spinconfigvects)
     else
        call assemble_configmat(fullconfigmatel,yyy%cptr(0),1,1,0,0, time)
        if (printflag.ne.0) then
           OFLWR " Call big eig ",numconfig*numr; CFL
        endif
        if (myrank.eq.1) then
           call CONFIGEIG(fullconfigmatel(:,:),totadim,totadim, fullconfigvects, fullconfigvals)
        endif
     endif

     call mympibcast(fullconfigvects,1,totadim**2)
     call mympibcast(fullconfigvals,1,totadim)

     if (printflag.ne.0) then
        OFLWR "  -- Nonsparse eigenvals --"
        do i=1,min(numconfig*numr,numshift+order+10)
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


     deallocate(fullconfigmatel, fullconfigvects, fullconfigvals)

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
     call configspin_projectall(avectorin,0)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict(avectorin,numr)
  endif
  if (sparseconfigflag/=0) then
     call expoconfigprop(avectorin,avectorout,time)
  else
     call nonsparseprop(avectorin,avectorout,time)
  endif
  if (allspinproject.ne.0) then
     call configspin_projectall(avectorout,0)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict(avectorout,numr)
  endif
end subroutine myconfigprop


subroutine nonsparseprop(avectorin,avectorout,time)
  use parameters
  use mpimod
  use configpropmod
  implicit none

  DATATYPE :: avectorin(totadim), avectorout(totadim),  avectortemp(totadim)
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
  allocate(iiwork(totadim*2))


  if (allspinproject.eq.0) then
     allocate(bigconfigmatel(totadim,totadim), bigconfigvects(totadim,2*(totadim+2)))
     call assemble_configmat(bigconfigmatel, workconfigpointer,1,1,1,1, time)

     bigconfigmatel=bigconfigmatel*timefac

     avectorout = avectorin

     if (myrank.eq.1) then

        if (nonsparsepropmode.eq.1) then
           call expmat(bigconfigmatel,numconfig*numr); call MYGEMV('N',numconfig*numr,numconfig*numr,DATAONE,bigconfigmatel,numconfig*numr,avectorin,1,DATAZERO,avectorout,1)
#ifndef REALGO
        else if (nonsparsepropmode.eq.2) then
           allocate(realbigconfigmatel(2,numconfig*numr,2,numconfig*numr))
           call assigncomplexmat(realbigconfigmatel,bigconfigmatel,numconfig*numr,numconfig*numr)
           call DGCHBV(numconfig*numr*2, 1.d0, realbigconfigmatel, numconfig*numr*2, avectorout, bigconfigvects, iiwork, iflag)           
           deallocate(realbigconfigmatel)
#endif
        else
           call EXPFULL(numconfig*numr, 1.d0, bigconfigmatel, numconfig*numr, avectorout, bigconfigvects, iiwork, iflag)

           if (iflag.ne.0) then
              OFLWR "Expo error A ", iflag; CFLST
           endif
        endif

     endif
     call mympibcast(avectorout,1,numconfig*numr)
  else

     allocate(bigconfigmatel(spintotrank*numr,spintotrank*numr), bigconfigvects(spintotrank*numr,2*(spintotrank*numr+2)))


     call assemble_spinconfigmat(bigconfigmatel, workconfigpointer,1,1,1,1, time)



     bigconfigmatel=bigconfigmatel*timefac
     call configspin_transformto(1,avectorin,avectortemp)
     if (myrank.eq.1) then

        if (nonsparsepropmode.eq.1) then
           avectorout=avectortemp  !! in spin basis
           call expmat(bigconfigmatel,spintotrank*numr); call MYGEMV('N',spintotrank*numr,spintotrank*numr,DATAONE,bigconfigmatel,spintotrank*numr,avectorout,1,DATAZERO,avectortemp,1)
#ifndef REALGO
        else if (nonsparsepropmode.eq.2) then
           allocate(realbigconfigmatel(2,spintotrank*numr,2,spintotrank*numr))
           call assigncomplexmat(realbigconfigmatel,bigconfigmatel,spintotrank*numr,spintotrank*numr)
           call DGCHBV(spintotrank*numr*2, 1.d0, realbigconfigmatel, spintotrank*numr*2, avectortemp, bigconfigvects, iiwork, iflag)           
           deallocate(realbigconfigmatel)
#endif
        else
           call EXPFULL(spintotrank*numr, 1.d0, bigconfigmatel, spintotrank*numr, avectortemp, bigconfigvects, iiwork, iflag)
           
           if (iflag.ne.0) then
              OFLWR "Expo error A ", iflag; CFLST
           endif
        endif
     endif
     call mympibcast(avectortemp,1,spintotrank*numr)
     call configspin_transformfrom(1,avectortemp,avectorout)
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

!! TAKES TRANSPOSES AS INPUT AND OUTPUT

!! NOTE BOUNDS !!


subroutine parconfigexpomult_transpose(inavectortr,outavectortr)
  use parameters
  use mpimod
  use configexpotimemod
  use configpropmod
  implicit none

  DATATYPE :: inavectortr(numr,botwalk:botwalk+maxconfigsperproc-1), outavectortr(numr,botwalk:botwalk+maxconfigsperproc-1)
  integer, save :: allochere=0
  DATATYPE,save,allocatable :: intemptr(:,:), ttvector(:,:), ttvector2(:,:)

  if (allochere.eq.0) then
     allocate(intemptr(numr,numconfig), ttvector(numconfig,numr), ttvector2(botwalk:topwalk,numr))
  endif
  allochere=1

  call avectortime(3)

  if (sparseconfigflag.eq.0) then
     OFLWR "error, must use sparse for parconfigexpomult_transpose_direct"; CFLST
  endif
  intemptr(:,:)=0d0;  intemptr(:,botwalk:topwalk)=inavectortr(:,botwalk:topwalk)
  
  call mpiallgather(intemptr,numconfig*numr,configsperproc*numr,maxconfigsperproc*numr)

  ttvector(:,:)=TRANSPOSE(intemptr(:,:))

  if (dfrestrictflag.ne.0) then
     call dfrestrict(ttvector,numr)
  endif

  call sparseconfigmult_nompi(ttvector,ttvector2, workconfigpointer, worksparsepointer, 1,1,1,1,configexpotime)


  if (dfrestrictflag.ne.0) then
     call dfrestrict_par(ttvector2,numr)
  endif
  outavectortr(:,:)=0d0
  outavectortr(:,botwalk:topwalk)=TRANSPOSE(ttvector2(botwalk:topwalk,:))
  outavectortr=outavectortr*timefac
  
  call avectortime(2)

end subroutine parconfigexpomult_transpose



subroutine parconfigexpomult_transpose_spin(inavectortrspin,outavectortrspin)
  use parameters
  use mpimod
  use configexpotimemod
  implicit none

  DATATYPE :: inavectortrspin(numr,spinstart:spinstart+maxspinrank-1), outavectortrspin(numr,spinstart:spinstart+maxspinrank-1)
  integer,save :: allochere=0
  DATATYPE,save,allocatable :: inavectortr(:,:),outavectortr(:,:)

  if (allochere.eq.0) then
     allocate(inavectortr(numr,botwalk:botwalk+maxconfigsperproc-1),outavectortr(numr,botwalk:botwalk+maxconfigsperproc-1))
  endif
  allochere=1

  call avectortime(3)

  if (spinwalkflag.eq.0) then
     OFLWR "WTF SPIN"; CFLST
  endif

  call configspin_transformfrom_mine_transpose(inavectortrspin,inavectortr)

  call avectortime(2)

  call parconfigexpomult_transpose(inavectortr,outavectortr)

  outavectortrspin(:,:)=0d0

  call configspin_transformto_mine_transpose(outavectortr,outavectortrspin)

  call avectortime(2)
  
end subroutine parconfigexpomult_transpose_spin



