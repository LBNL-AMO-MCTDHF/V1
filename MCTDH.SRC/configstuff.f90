
#include "Definitions.INC"

subroutine printconfig(thisconfig,www)
  use walkmod
  use fileptrmod
  implicit none

  type(walktype),intent(in) :: www
  integer :: thisconfig(www%ndof),i
  character (len=4) :: mslabels(2) =["a ","b "]
  
  write(mpifileptr,'(100(I3,A2))') (thisconfig((i-1)*2+1), mslabels(thisconfig(i*2)), i=1,www%numelec)

end subroutine printconfig


!! NOT - this uses stopthresh as threshold to pass to laneigen.  is called for mcscf, and improved.


subroutine myconfigeig(www,dfww,cptr,thisconfigvects,thisconfigvals,order,printflag, guessflag,time,numshift)
  use fileptrmod
  use r_parameters
  use sparse_parameters
  use mpimod
  use configptrmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www,dfww
  type(CONFIGPTR),intent(in) :: cptr
  integer :: order, printflag,i, guessflag,  l,k,numshift,flag
  DATATYPE :: thisconfigvects(www%totadim,order),lastval,dot
  DATAECS :: thisconfigvals(order), tempconfigvals(order+numshift)
  real*8 :: realconfigvect(www%totadim)
  real*8 :: sum,time
  DATATYPE, allocatable :: fullconfigmatel(:,:), fullconfigvects(:,:), &
       tempconfigvects(:,:)
  DATAECS, allocatable :: fullconfigvals(:)

  if (order+numshift.gt.www%totadim) then
     OFLWR "Error, want ",order," plus ",numshift," vectors but totadim= ",www%totadim;CFLST
  endif
  if (numshift.lt.0) then
     OFLWR "GG ERROR.", numshift,guessflag; CFLST
  endif

  if (www%numdfbasis.ne.dfww%numdfbasis) then
     OFLWR "ERROR DF SETS MYCONFIGEIG",www%numdfbasis,dfww%numdfbasis; CFLST
  endif

  if (sparseconfigflag/=0) then
     
     allocate(tempconfigvects(www%totadim,order+numshift))
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
        OFLWR "Construct big matrix:  ", www%numdfbasis*numr; CFL
     endif
     allocate(fullconfigvals(www%numdfbasis*numr))
     allocate(fullconfigmatel(www%numdfbasis*numr,www%numdfbasis*numr), fullconfigvects(www%numconfig*numr,www%numdfbasis*numr))
     allocate(tempconfigvects(www%numdfbasis*numr,www%numdfbasis*numr))
     fullconfigmatel=0.d0; fullconfigvects=0d0; fullconfigvals=0d0; tempconfigvects=0d0

     call assemble_dfbasismat(dfww,fullconfigmatel,cptr,1,1,0,0, time)

     if (printflag.ne.0) then
        OFLWR " Call eig ",www%numdfbasis*numr; CFL
     endif
     if (myrank.eq.1) then
        call CONFIGEIG(fullconfigmatel(:,:),www%numdfbasis*numr,www%numdfbasis*numr, tempconfigvects, fullconfigvals)
        do i=1,www%numdfbasis*numr
           call basis_transformfrom_all(www,numr,tempconfigvects(:,i),fullconfigvects(:,i))
        enddo
     endif

     call mympibcast(fullconfigvects,1,www%numdfbasis*numr*www%numconfig*numr)
     call mympibcast(fullconfigvals,1,www%numdfbasis*numr)

     if (printflag.ne.0) then
        OFLWR "  -- Nonsparse eigenvals --"
        do i=1,min(www%numdfbasis*numr,numshift+order+10)
           WRFL "   ", fullconfigvals(i)
        enddo
        CFL
     endif

     thisconfigvects(:,1:order) = fullconfigvects(:,1+numshift:order+numshift)
     thisconfigvals(1:order) = fullconfigvals(1+numshift:order+numshift)

     lastval = -99999999d0
     flag=0
     do i=1,order
        thisconfigvects(:,i)=thisconfigvects(:,i) / sqrt(dot(thisconfigvects(:,i),thisconfigvects(:,i),www%totadim))
        if (abs(thisconfigvals(i)-lastval).lt.1d-7) then
           flag=flag+1
           call gramschmidt(www%totadim,flag,www%totadim,thisconfigvects(:,i-flag),thisconfigvects(:,i),.false.)
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
     do k=1,www%totadim
        if (abs(thisconfigvects(k,i)).gt.sum) then
           sum=abs(thisconfigvects(k,i));           l=k
        endif
     enddo
     thisconfigvects(:,i)=thisconfigvects(:,i)*abs(thisconfigvects(l,i))/thisconfigvects(l,i)
  enddo

end subroutine myconfigeig



!! PROPAGATE A-VECTOR .  CALLED WITHIN LITTLESTEPS LOOP

subroutine myconfigprop(www,dfww,avectorin,avectorout,time)
  use r_parameters
  use sparse_parameters
  use walkmod 
  implicit none
  type(walktype),intent(in) :: www,dfww
  real*8 :: time
  DATATYPE :: avectorin(www%totadim), avectorout(www%totadim)

  if (sparseconfigflag/=0) then
     call exposparseprop(www,avectorin,avectorout,time)
  else
     call nonsparseprop(www,dfww,avectorin,avectorout,time)
  endif

  call basis_project(www,numr,avectorout)

end subroutine myconfigprop


subroutine nonsparseprop(www,dfww,avectorin,avectorout,time)
  use fileptrmod
  use sparse_parameters
  use ham_parameters
  use r_parameters
  use mpimod
  use configpropmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www,dfww
  DATATYPE :: avectorin(www%totadim), avectorout(www%totadim),  &
       avectortemp(www%numdfbasis*numr), avectortemp2(www%numdfbasis*numr)
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

  if (www%numdfbasis.ne.dfww%numdfbasis) then
     OFLWR "ERROR DF SETS NONSPARSEOPROP",www%numdfbasis,dfww%numdfbasis; CFLST
  endif

  allocate(iiwork(www%numdfbasis*numr*2))

  allocate(bigconfigmatel(www%numdfbasis*numr,www%numdfbasis*numr), bigconfigvects(www%numdfbasis*numr,2*(www%numdfbasis*numr+2)))

  call assemble_dfbasismat(dfww,bigconfigmatel, workconfigpointer,1,1,1,1, time)

  bigconfigmatel=bigconfigmatel*timefac

  call basis_transformto_all(www,numr,avectorin,avectortemp)

  if (myrank.eq.1) then
     if (nonsparsepropmode.eq.1) then
        avectortemp2(:)=avectortemp(:)
        call expmat(bigconfigmatel,www%numdfbasis*numr) 
        call MYGEMV('N',www%numdfbasis*numr,www%numdfbasis*numr,DATAONE,bigconfigmatel,www%numdfbasis*numr,avectortemp2,1,DATAZERO,avectortemp,1)
#ifndef REALGO
     else if (nonsparsepropmode.eq.2) then
        allocate(realbigconfigmatel(2,www%numdfbasis*numr,2,www%numdfbasis*numr))
        call assigncomplexmat(realbigconfigmatel,bigconfigmatel,www%numdfbasis*numr,www%numdfbasis*numr)
        call DGCHBV(www%numdfbasis*numr*2, 1.d0, realbigconfigmatel, www%numdfbasis*numr*2, avectortemp, bigconfigvects, iiwork, iflag)           
        deallocate(realbigconfigmatel)
#endif
     else
        call EXPFULL(www%numdfbasis*numr, DATAONE, bigconfigmatel, www%numdfbasis*numr, avectortemp, bigconfigvects, iiwork, iflag)
        
        if (iflag.ne.0) then
           OFLWR "Expo error A ", iflag; CFLST
        endif
     endif
  endif

  call mympibcast(avectortemp,1,www%numdfbasis*numr)

  call basis_transformfrom_all(www,numr,avectortemp,avectorout)

  deallocate(iiwork,bigconfigmatel,bigconfigvects)

end subroutine nonsparseprop




!! MULTIPLY BY A-VECTOR HAMILTONIAN MATRIX

!! NOTE BOUNDS !!  PADDED

subroutine parconfigexpomult_padded(inavector,outavector)
  use fileptrmod
  use r_parameters
  use sparse_parameters
  use ham_parameters   !! timefac
  use mpimod
  use configexpotimemod
  use configpropmod
  use configmod
  implicit none
  DATATYPE,intent(in) :: inavector(numr,www%botdfbasis:www%botdfbasis+www%maxdfbasisperproc-1)
  DATATYPE,intent(out) :: outavector(numr,www%botdfbasis:www%botdfbasis+www%maxdfbasisperproc-1)

  if (sparsesummaflag.eq.0) then
     if (www%dfrestrictflag.eq.0.or.sparsedfflag.eq.0) then
        call parconfigexpomult_padded0_gather(www,workconfigpointer,worksparsepointer,inavector,outavector)
     else
        call parconfigexpomult_padded0_gather(dfww,workconfigpointer,workdfsparsepointer,inavector,outavector)
     endif
  else
     if (www%dfrestrictflag.eq.0.or.sparsedfflag.eq.0) then
        call parconfigexpomult_padded0_summa(www,workconfigpointer,worksparsepointer,inavector,outavector)
     else
        call parconfigexpomult_padded0_summa(dfww,workconfigpointer,workdfsparsepointer,inavector,outavector)
     endif
  endif

end subroutine parconfigexpomult_padded



subroutine parconfigexpomult_padded0_gather(www,workconfigpointer,worksparsepointer,inavector,outavector)
  use fileptrmod
  use r_parameters
  use sparse_parameters
  use ham_parameters   !! timefac
  use mpimod
  use configexpotimemod
  use walkmod
  use configptrmod
  use sparseptrmod
  implicit none
  type(walktype),intent(in) :: www
  type(CONFIGPTR),intent(in) :: workconfigpointer
  type(SPARSEPTR),intent(in) :: worksparsepointer
  DATATYPE,intent(in) :: inavector(numr,www%botdfbasis:www%botdfbasis+www%maxdfbasisperproc-1)
  DATATYPE,intent(out) :: outavector(numr,www%botdfbasis:www%botdfbasis+www%maxdfbasisperproc-1)
  DATATYPE,allocatable :: intemp(:,:)
  DATATYPE :: outtemp(numr,www%botconfig:www%topconfig)

  call avectortime(3)

  if (sparseconfigflag.eq.0) then
     OFLWR "error, must use sparse for parconfigexpomult"; CFLST
  endif
  
  allocate(intemp(numr,www%numconfig))

!! TRANSFORM SECOND TO REDUCE COMMUNICATION?

  if (www%topconfig-www%botconfig+1 .ne. 0) then
     call basis_transformfrom_local(www,numr,inavector,intemp(:,www%botconfig:www%topconfig))
  endif

  call mpiallgather(intemp,www%numconfig*numr,www%configsperproc(:)*numr,www%maxconfigsperproc*numr)

  call sparseconfigmult_byproc(1,nprocs,www,intemp,outtemp, workconfigpointer, worksparsepointer, 1,1,1,1,configexpotime,0,1,numr,0)
     
  outavector(:,:)=0d0   !! PADDED

  call basis_transformto_local(www,numr,outtemp,outavector)

  outavector=outavector*timefac

  deallocate(intemp)

  call avectortime(2)

end subroutine parconfigexpomult_padded0_gather



subroutine parconfigexpomult_padded0_summa(www,workconfigpointer,worksparsepointer,inavector,outavector)
  use fileptrmod
  use r_parameters
  use sparse_parameters
  use ham_parameters   !! timefac
  use mpimod
  use configexpotimemod
  use walkmod
  use configptrmod
  use sparseptrmod
  implicit none
  type(walktype),intent(in) :: www
  type(CONFIGPTR),intent(in) :: workconfigpointer
  type(SPARSEPTR),intent(in) :: worksparsepointer
  DATATYPE,intent(in) :: inavector(numr,www%botdfbasis:www%botdfbasis+www%maxdfbasisperproc-1)
  DATATYPE,intent(out) :: outavector(numr,www%botdfbasis:www%botdfbasis+www%maxdfbasisperproc-1)
  DATATYPE :: intemp(numr,www%maxconfigsperproc), outwork(numr,www%botconfig:www%topconfig),&
       outtemp(numr,www%botconfig:www%topconfig)
  integer :: iproc

  call avectortime(3)

  if (sparseconfigflag.eq.0) then
     OFLWR "error, must use sparse for parconfigexpomult"; CFLST
  endif
  
  outwork(:,:)=0d0

  do iproc=1,nprocs

!! TRANSFORM SECOND TO REDUCE COMMUNICATION?

     if (myrank.eq.iproc) then
        call basis_transformfrom_local(www,numr,inavector,intemp)
     endif

     call mympibcast(intemp,iproc,(www%alltopconfigs(iproc)-www%allbotconfigs(iproc)+1)*numr)

     call sparseconfigmult_byproc(iproc,iproc,www,intemp,outtemp, workconfigpointer, worksparsepointer, 1,1,1,1,configexpotime,0,1,numr,0)
     
     outwork(:,:)=outwork(:,:)+outtemp(:,:)

  enddo

  outavector(:,:)=0d0   !! PADDED

  call basis_transformto_local(www,numr,outwork,outavector)

  outavector=outavector*timefac
  
  call avectortime(2)

end subroutine parconfigexpomult_padded0_summa





