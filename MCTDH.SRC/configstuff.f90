
#include "Definitions.INC"

!! ALL MODULES

module configstuffmod
contains

!! guessflag=1  use thisconfigsvects for initial krylov vectors
!!          =2  do that and also follow, return eigvects with max overlap

subroutine myconfigeig(cptr,thisconfigvects,thisconfigvals,order,printflag, &
     guessflag,time,numshift)
  use fileptrmod
  use r_parameters
  use sparse_parameters
  use mpimod
  use configptrmod
  use configmod
  use basissubmod
  use mpisubmod
  use asssubmod
  use lanblockmod
  use eigenmod
  implicit none
  type(CONFIGPTR),intent(in) :: cptr
  integer,intent(in) :: order,printflag,guessflag,numshift
  DATATYPE,intent(out) :: thisconfigvects(www%totadim,order)
  DATAECS,intent(out) :: thisconfigvals(order)
  real*8,intent(in) :: time
  real*8,allocatable :: realconfigvect(:)
  DATATYPE, allocatable :: fullconfigmatel(:,:), fullconfigvects(:,:), &
       tempconfigvects(:,:)
  DATAECS, allocatable :: fullconfigvals(:)
  DATATYPE, allocatable :: followdots(:,:), tempeigvect(:)
  DATAECS ::  tempconfigvals(order+numshift), tempeigval
!!$  DATATYPE :: lastval,dot,csum
  real*8 :: rmaxdot
  integer :: i,j,ifollow

  if (order+numshift.gt.www%numdfbasis*numr) then
     OFLWR "Error, want ",order," plus ",numshift,&
          " vectors but total dimension= ",www%numdfbasis*numr;CFLST
  endif
  if (numshift.lt.0) then
     OFLWR "GG ERROR.", numshift,guessflag; CFLST
  endif

  if (www%numdfbasis.ne.dwwptr%numdfbasis) then
     OFLWR "ERROR DF SETS MYCONFIGEIG",www%numdfbasis,dwwptr%numdfbasis; CFLST
  endif

  if (numshift.gt.0.and.guessflag.eq.2) then
     OFLWR "Programmer fail, numshift.gt.0 not available with followflag",numshift,guessflag; CFLST
  endif
  
  if (sparseconfigflag/=0) then
     
     allocate(tempconfigvects(www%totadim,order+numshift))
     if (www%totadim.gt.0) then
        tempconfigvects(:,:)=0d0
     endif
     if (guessflag.ne.0) then
!! never mind consistency regardless of number of processors.  could scatter anyway.  
!!        allocate(realconfigvect(www%numconfig*numr));    realconfigvect=0
!!        do i=1,numshift
!!           call RANDOM_NUMBER(realconfigvect(:))
!!           if (www%totadim.gt.0) then
!!              tempconfigvects(:,i)=realconfigvect((www%firstconfig-1)*numr+1:www%lastconfig*numr)
!!           endif
!!        enddo

        if (www%totadim.gt.0) then
           if (numshift.gt.0) then
              allocate(realconfigvect(www%totadim*numr))
              do i=1,numshift
                 call RANDOM_NUMBER(realconfigvect(:))
                 tempconfigvects(:,i)=realconfigvect(:)
              enddo
              deallocate(realconfigvect)
           endif
           tempconfigvects(:,numshift+1:order+numshift)=thisconfigvects(:,1:order)
        endif
     endif
     call blocklanczos(order+numshift,tempconfigvects,tempconfigvals,printflag,guessflag)
     
     thisconfigvals(1:order)=tempconfigvals(numshift+1:numshift+order)
     if (www%totadim.gt.0) then
        thisconfigvects(:,1:order)=tempconfigvects(:,numshift+1:numshift+order)
     endif
     deallocate(tempconfigvects)

  else

     if (printflag.ne.0) then
        OFLWR "Construct big matrix:  ", www%numdfbasis*numr; CFL
     endif
     allocate(fullconfigvals(www%numdfbasis*numr),&
          fullconfigmatel(www%numdfbasis*numr,www%numdfbasis*numr), &
          fullconfigvects(www%numconfig*numr,www%numdfbasis*numr),&
          tempconfigvects(www%numdfbasis*numr,www%numdfbasis*numr))
     fullconfigvects=0d0; tempconfigvects=0d0; fullconfigmatel=0.d0; fullconfigvals=0d0

     call assemble_dfbasismat(dwwptr,fullconfigmatel,cptr,1,1,0,0, time,-1)

     if (printflag.ne.0) then
        OFLWR " Call eig ",www%numdfbasis*numr; CFL
     endif
     if (myrank.eq.1) then
        call CONFIGEIG(fullconfigmatel(:,:),www%numdfbasis*numr,www%numdfbasis*numr, &
             tempconfigvects, fullconfigvals)
     endif

     call mympibcast(tempconfigvects,1,www%numdfbasis*numr*www%numdfbasis*numr)
#ifndef ECSFLAG
#ifndef REALGO
     call mympirealbcast(fullconfigvals,1,www%numdfbasis*numr)
#else
     call mympibcast(fullconfigvals,1,www%numdfbasis*numr)
#endif
#else
     call mympibcast(fullconfigvals,1,www%numdfbasis*numr)
#endif
     do i=1,www%numdfbasis*numr
        call basis_transformfrom_all(www,numr,tempconfigvects(:,i),fullconfigvects(:,i))
     enddo

     call mpibarrier()

     if (guessflag.eq.2) then !! followflag

        allocate(followdots(order,www%numdfbasis*numr), tempeigvect(www%numconfig*numr))

        followdots=0; tempeigvect=0

        call MYGEMM(CNORMCHAR,'N',order,www%numdfbasis*numr,www%numconfig*numr, DATAONE, &
             thisconfigvects(:,:), www%numconfig*numr, fullconfigvects(:,:), www%numconfig*numr, &
             DATAZERO, followdots(:,:), order)


        do i=1,order
           ifollow=(-1)
           rmaxdot=0d0
           do j=i,www%numdfbasis*numr                   !! do j=i, followtaken not needed
              if (abs(followdots(i,j)).gt.rmaxdot) then !! .and.followtaken(j).eq.0) then
                 rmaxdot=abs(followdots(i,j))
                 ifollow=j
              endif
           enddo

           if (ifollow.ne.i) then
              tempeigvect(:) = fullconfigvects(:,i)
              fullconfigvects(:,i) = fullconfigvects(:,ifollow)
              fullconfigvects(:,ifollow) = tempeigvect(:)

              tempeigval = fullconfigvals(i)
              fullconfigvals(i) = fullconfigvals(ifollow)
              fullconfigvals(ifollow) = tempeigval
           endif
        enddo

        deallocate(followdots, tempeigvect)

     endif

     if (printflag.ne.0) then
        OFLWR "  -- Nonsparse eigenvals --"
        do i=1,min(www%numdfbasis*numr,numshift+order+10)
           WRFL "   ", fullconfigvals(i)
        enddo
        CFL
     endif

     thisconfigvects(:,1:order) = fullconfigvects(:,1+numshift:order+numshift)
     thisconfigvals(1:order) = fullconfigvals(1+numshift:order+numshift)


!!$  REMOVED THIS 2-2016 was buggy par_conplit.ne.0 logic was not implemented
!!$   implemented now but commented out.  only valid for pmchtdhf and cmctdhf not chmctdhf
!!$   if degeneracy is accidental (unlikely) as opposed to due to symmetry
!!$   in which case this orthogonality enforcement is fine.  But eigenvectors should
!!$   come out of diagonalization routines properly orthogonalized.
!!$  
!!$       lastval = -99999999d0
!!$       flag=0
!!$       do i=1,order
!!$  
!!$          csum=dot(thisconfigvects(:,i),thisconfigvects(:,i),www%totadim)
!!$          if (par_consplit.ne.0) then
!!$             call mympireduceone(csum)
!!$          endif
!!$          thisconfigvects(:,i)=thisconfigvects(:,i) / sqrt(csum)
!!$  
!!$          if (abs(thisconfigvals(i)-lastval).lt.1d-7) then
!!$             flag=flag+1
!!$             if (par_consplit.ne.0) then
!!$                call gramschmidt(www%totadim,flag,www%totadim,&
!!$                   thisconfigvects(:,i-flag),thisconfigvects(:,i),.true.)
!!$             else
!!$                call gramschmidt(www%totadim,flag,www%totadim,&
!!$                   thisconfigvects(:,i-flag),thisconfigvects(:,i),.false.)
!!$             endif
!!$  
!!$          else
!!$             flag=0
!!$          endif
!!$          lastval=thisconfigvals(i)
!!$       enddo

     deallocate(fullconfigmatel, fullconfigvects, fullconfigvals,tempconfigvects)

  endif   !! sparseconfigflag

!!$ REMOVED THIS FOR PARCONSPLIT 12-2015 v1.16 REINSTATE?
!!$  
!!$  !! this is called from either eigs (for which is irrelevant) or from 
!!$  !! config_eigen, which is used to start avectors as improved relax or at
!!$  !! start of prop with eigenflag.  Leave as c-normed; in config_eigen, 
!!$  !! herm-norm.
!!$  !! fix phase for chmctdh/pmctdh debug check
!!$  
!!$    l=-99
!!$    do i=1,order
!!$       sum=(-999d0)
!!$       do k=1,www%totadim
!!$          if (abs(thisconfigvects(k,i)).gt.sum) then
!!$             sum=abs(thisconfigvects(k,i));           l=k
!!$          endif
!!$       enddo
!!$       thisconfigvects(:,i)=thisconfigvects(:,i)*abs(thisconfigvects(l,i))/thisconfigvects(l,i)
!!$    enddo

end subroutine myconfigeig



!! PROPAGATE A-VECTOR .  CALLED WITHIN LITTLESTEPS LOOP

subroutine myconfigprop(avectorin,avectorout,time,imc,numiters)
  use r_parameters
  use sparse_parameters
  use configmod
  use basissubmod
  use expoavecpropmod
  implicit none
  integer,intent(in) :: imc
  integer,intent(out) :: numiters
  real*8,intent(in) :: time
  DATATYPE,intent(in) :: avectorin(www%totadim)
  DATATYPE,intent(out) :: avectorout(www%totadim)
  DATATYPE :: nullvector1(numr),nullvector2(numr)

  numiters=0
  if (sparseconfigflag/=0) then
     if (www%totadim.gt.0) then
        call expoavecprop(avectorin,avectorout,time,imc,numiters)
     else
        call expoavecprop(nullvector1,nullvector2,time,imc,numiters)
     endif
  else
     call nonsparseprop(www,dwwptr,avectorin,avectorout,time,imc)
  endif

  call basis_project(www,numr,avectorout)

end subroutine myconfigprop


subroutine nonsparseprop(wwin,dwin,avectorin,avectorout,time,imc)
  use fileptrmod
  use sparse_parameters
  use ham_parameters
  use r_parameters
  use mpimod
  use configpropmod
  use walkmod
  use basissubmod
  use expsubmod
  use mpisubmod
  use asssubmod
  use expokitmod, only: EXPFULL, DGCHBVQ
  use utilmod
  implicit none
  type(walktype),intent(in) :: wwin,dwin
  integer, intent(in) :: imc
  DATATYPE,intent(in) :: avectorin(wwin%totadim)
  DATATYPE,intent(out) :: avectorout(wwin%totadim)
  DATATYPE :: avectortemp(wwin%numdfbasis*numr), &
       avectortemp2(wwin%numdfbasis*numr) !! AUTOMATIC
  DATATYPE, allocatable :: bigconfigmatel(:,:)
#ifndef REALGO
  real*8, allocatable :: realbigconfigmatel(:,:,:,:)
#endif
  integer :: iflag
  real*8 :: time

  if (sparseconfigflag/=0) then
     OFLWR "Error, nonsparseprop called but sparseconfigflag= ", sparseconfigflag; CFLST
  endif
  if (drivingflag.ne.0) then
     OFLWR "Driving flag not implemented for nonsparse"; CFLST
  endif

  if (wwin%numdfbasis.ne.dwin%numdfbasis) then
     OFLWR "ERROR DF SETS NONSPARSEOPROP",wwin%numdfbasis,dwin%numdfbasis; CFLST
  endif

  allocate(bigconfigmatel(wwin%numdfbasis*numr,wwin%numdfbasis*numr))
  bigconfigmatel=0

  call assemble_dfbasismat(dwin,bigconfigmatel, workconfigpointer,1,1,1,1, time,imc)

  bigconfigmatel=bigconfigmatel*timefac

  call basis_transformto_all(wwin,numr,avectorin,avectortemp)

  if (myrank.eq.1) then
     if (nonsparsepropmode.eq.1) then
        avectortemp2(:)=avectortemp(:)
        call expmat(bigconfigmatel,wwin%numdfbasis*numr) 
        call MYGEMV('N',wwin%numdfbasis*numr,wwin%numdfbasis*numr,DATAONE,&
             bigconfigmatel,wwin%numdfbasis*numr,avectortemp2,1,DATAZERO,avectortemp,1)
#ifndef REALGO
     else if (nonsparsepropmode.eq.2) then
        allocate(realbigconfigmatel(2,wwin%numdfbasis*numr,2,wwin%numdfbasis*numr))
        realbigconfigmatel=0
        call assigncomplexmat(realbigconfigmatel,bigconfigmatel,&
             wwin%numdfbasis*numr,wwin%numdfbasis*numr)
        call DGCHBVQ(wwin%numdfbasis*numr*2, 1.d0, realbigconfigmatel, &
             wwin%numdfbasis*numr*2, avectortemp, iflag)           
        deallocate(realbigconfigmatel)
#endif
     else
        call EXPFULL(wwin%numdfbasis*numr, 1d0, bigconfigmatel, &
             wwin%numdfbasis*numr, avectortemp, iflag)
        
        if (iflag.ne.0) then
           OFLWR "Expo error A ", iflag; CFLST
        endif
     endif
  endif

  call mympibcast(avectortemp,1,wwin%numdfbasis*numr)

  call basis_transformfrom_all(wwin,numr,avectortemp,avectorout)

  deallocate(bigconfigmatel)

end subroutine nonsparseprop

end module configstuffmod
