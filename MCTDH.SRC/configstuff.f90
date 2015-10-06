
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
       spinconfigmatel(:,:), spinconfigvects(:,:), tempconfigvects(:,:)
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
     allocate(fullconfigvals(totadim))
     allocate(fullconfigmatel(numconfig*numr,numconfig*numr), fullconfigvects(totadim,totadim))
     fullconfigmatel=0.d0; fullconfigvects=0d0; fullconfigvals=0d0

     if (spineigflag.ne.0.and.allspinproject.ne.0) then
        allocate(spinconfigmatel(numspinconfig*numr,numspinconfig*numr), spinconfigvects(numspinconfig*numr,numspinconfig*numr))
        spinconfigmatel=0.d0

        call assemble_spinconfigmat(spinconfigmatel,yyy%cptr(0),1,1,0,0,  time)
        if (printflag.ne.0) then
           OFLWR " Call spin eig ",numspinconfig*numr; CFL
        endif
        if (myrank.eq.1) then
           call CONFIGEIG(spinconfigmatel(:,:),numspinconfig*numr,numspinconfig*numr, spinconfigvects, fullconfigvals)
           do i=1,numspinconfig*numr
              call configspin_transformfrom(numr,spinconfigvects(:,i),fullconfigvects(:,i))
           enddo
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

     allocate(bigconfigmatel(numspinconfig*numr,numspinconfig*numr), bigconfigvects(numspinconfig*numr,2*(numspinconfig*numr+2)))


     call assemble_spinconfigmat(bigconfigmatel, workconfigpointer,1,1,1,1, time)



     bigconfigmatel=bigconfigmatel*timefac
     call configspin_transformto(numr,avectorin,avectortemp)

     if (myrank.eq.1) then
        if (nonsparsepropmode.eq.1) then
           avectorout=avectortemp  !! in spin basis
           call expmat(bigconfigmatel,numspinconfig*numr); call MYGEMV('N',numspinconfig*numr,numspinconfig*numr,DATAONE,bigconfigmatel,numspinconfig*numr,avectorout,1,DATAZERO,avectortemp,1)
#ifndef REALGO
        else if (nonsparsepropmode.eq.2) then
           allocate(realbigconfigmatel(2,numspinconfig*numr,2,numspinconfig*numr))
           call assigncomplexmat(realbigconfigmatel,bigconfigmatel,numspinconfig*numr,numspinconfig*numr)
           call DGCHBV(numspinconfig*numr*2, 1.d0, realbigconfigmatel, numspinconfig*numr*2, avectortemp, bigconfigvects, iiwork, iflag)           
           deallocate(realbigconfigmatel)
#endif
        else
           call EXPFULL(numspinconfig*numr, 1.d0, bigconfigmatel, numspinconfig*numr, avectortemp, bigconfigvects, iiwork, iflag)
           
           if (iflag.ne.0) then
              OFLWR "Expo error A ", iflag; CFLST
           endif
        endif
     endif
     call mympibcast(avectortemp,1,numspinconfig*numr)
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

!! TAKES TRANSPOSES AS INPUT AND OUTPUT

!! NOTE BOUNDS !!


subroutine parconfigexpomult(inavector,outavector)
  use parameters
  use mpimod
  use configexpotimemod
  use configpropmod
  implicit none

  DATATYPE,intent(in) :: inavector(numr,botwalk:botwalk+maxconfigsperproc-1)
  DATATYPE,intent(out) :: outavector(numr,botwalk:botwalk+maxconfigsperproc-1)
  DATATYPE :: intemp(numr,numconfig)

  call avectortime(3)

  if (sparseconfigflag.eq.0) then
     OFLWR "error, must use sparse for parconfigexpomult"; CFLST
  endif
  intemp(:,:)=0d0;  intemp(:,botwalk:topwalk)=inavector(:,botwalk:topwalk)

  if (dfrestrictflag.ne.0) then
     call df_project_local(intemp(:,botwalk),numr)
  endif

!! DO SUMMA  
  call mpiallgather(intemp,numconfig*numr,configsperproc*numr,maxconfigsperproc*numr)

  outavector(:,:)=0d0
  call sparseconfigmult_nompi(intemp(:,:),outavector(:,:), workconfigpointer, worksparsepointer, 1,1,1,1,configexpotime,0,1,numr,0)


  if (dfrestrictflag.ne.0) then
     call df_project_local(outavector,numr)
  endif

  outavector=outavector*timefac
  
  call avectortime(2)

end subroutine parconfigexpomult



subroutine parconfigexpomult_spin(inavectorspin,outavectorspin)
  use parameters
  use mpimod
  use configexpotimemod
  implicit none

  DATATYPE,intent(in) :: inavectorspin(numr,spinstart:spinstart+maxspinsperproc-1)
  DATATYPE,intent(out) :: outavectorspin(numr,spinstart:spinstart+maxspinsperproc-1)
  DATATYPE :: inavector(numr,botwalk:botwalk+maxconfigsperproc-1),&
       outavector(numr,botwalk:botwalk+maxconfigsperproc-1)

  call avectortime(3)

  if (spinwalkflag.eq.0) then
     OFLWR "WTF SPIN"; CFLST
  endif

  call configspin_transformfrom_local(numr,inavectorspin,inavector)

  call avectortime(2)

  call parconfigexpomult(inavector,outavector)

  outavectorspin(:,:)=0d0

  call configspin_transformto_local(numr,outavector,outavectorspin)

  call avectortime(2)
  
end subroutine parconfigexpomult_spin



