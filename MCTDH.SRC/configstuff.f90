
#include "Definitions.INC"

subroutine printconfig(thisconfig,www)
  use walkmod
  use fileptrmod
  implicit none

  type(walktype),intent(in) :: www
  integer,intent(in) :: thisconfig(www%ndof)
  character (len=4) :: mslabels(2) =["a ","b "]
  integer :: i

  write(mpifileptr,'(100(I3,A2))') (thisconfig((i-1)*2+1), &
       mslabels(thisconfig(i*2)), i=1,www%numelec)

end subroutine printconfig


subroutine myconfigeig(cptr,thisconfigvects,thisconfigvals,order,printflag, &
     guessflag,time,numshift)
  use fileptrmod
  use r_parameters
  use sparse_parameters
  use mpimod
  use configptrmod
  use configmod
  use basissubmod
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
  DATAECS ::  tempconfigvals(order+numshift)
!!$  DATATYPE :: lastval,dot,csum
  integer :: i

  if (order+numshift.gt.www%numconfig*numr) then
     OFLWR "Error, want ",order," plus ",numshift,&
          " vectors but totadim= ",www%numconfig*numr;CFLST
  endif
  if (numshift.lt.0) then
     OFLWR "GG ERROR.", numshift,guessflag; CFLST
  endif

  if (www%numdfbasis.ne.dwwptr%numdfbasis) then
     OFLWR "ERROR DF SETS MYCONFIGEIG",www%numdfbasis,dwwptr%numdfbasis; CFLST
  endif

  if (sparseconfigflag/=0) then
     
     allocate(tempconfigvects(www%totadim,order+numshift))
     if (www%totadim.gt.0) then
        tempconfigvects(:,:)=0d0
     endif
     if (guessflag.ne.0) then
        allocate(realconfigvect(www%numconfig));    realconfigvect=0
        do i=1,numshift
           call RANDOM_NUMBER(realconfigvect(:))
           if (www%totadim.gt.0) then
              tempconfigvects(:,i)=realconfigvect(www%firstconfig*numr:www%lastconfig*numr)
           endif
        enddo
        deallocate(realconfigvect)
        if (www%totadim.gt.0) then
           tempconfigvects(:,numshift+1:order+numshift)=thisconfigvects(:,1:order)
        endif
     endif
     call blocklanczos(order+numshift,tempconfigvects,tempconfigvals,printflag,guessflag)
     
        !! so if guessflag=1, numshift is 0.

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
        call exposparseprop(avectorin,avectorout,time,imc,numiters)
     else
        call exposparseprop(nullvector1,nullvector2,time,imc,numiters)
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
  implicit none
  type(walktype),intent(in) :: wwin,dwin
  integer, intent(in) :: imc
  DATATYPE,intent(in) :: avectorin(wwin%totadim)
  DATATYPE,intent(out) :: avectorout(wwin%totadim)
  DATATYPE :: avectortemp(wwin%numdfbasis*numr), &
       avectortemp2(wwin%numdfbasis*numr) !! AUTOMATIC
  DATATYPE, allocatable :: bigconfigmatel(:,:), bigconfigvects(:,:)
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

  if (wwin%numdfbasis.ne.dwin%numdfbasis) then
     OFLWR "ERROR DF SETS NONSPARSEOPROP",wwin%numdfbasis,dwin%numdfbasis; CFLST
  endif

  allocate(iiwork(wwin%numdfbasis*numr*2));     iiwork=0

  allocate(bigconfigmatel(wwin%numdfbasis*numr,wwin%numdfbasis*numr), &
       bigconfigvects(wwin%numdfbasis*numr,2*(wwin%numdfbasis*numr+2)))
  bigconfigmatel=0; bigconfigvects=0

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
        call DGCHBV(wwin%numdfbasis*numr*2, 1.d0, realbigconfigmatel, &
             wwin%numdfbasis*numr*2, avectortemp, bigconfigvects, iiwork, iflag)           
        deallocate(realbigconfigmatel)
#endif
     else
        call EXPFULL(wwin%numdfbasis*numr, DATAONE, bigconfigmatel, &
             wwin%numdfbasis*numr, avectortemp, bigconfigvects, iiwork, iflag)
        
        if (iflag.ne.0) then
           OFLWR "Expo error A ", iflag; CFLST
        endif
     endif
  endif

  call mympibcast(avectortemp,1,wwin%numdfbasis*numr)

  call basis_transformfrom_all(wwin,numr,avectortemp,avectorout)

  deallocate(iiwork,bigconfigmatel,bigconfigvects)

end subroutine nonsparseprop


