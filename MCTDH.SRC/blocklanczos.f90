
#include "Definitions.INC"

!! ALL MODULES

module parblocklanmod
contains

  subroutine parblockconfigmult0_gather(wwin,cptr,sptr,inavector,outavector)
    use fileptrmod
    use r_parameters
    use sparse_parameters
    use ham_parameters
    use mpimod
    use walkmod
    use sparseptrmod
    use configptrmod
    use sparsemultmod
    use basissubmod
    use mpisubmod
    implicit none
    type(walktype),intent(in) :: wwin
    type(CONFIGPTR),intent(in) :: cptr
    type(SPARSEPTR),intent(in) :: sptr
    integer :: ii
    DATATYPE,intent(in) :: inavector(numr,wwin%botdfbasis:wwin%topdfbasis)
    DATATYPE,intent(out) :: outavector(numr,wwin%botdfbasis:wwin%topdfbasis)
    DATATYPE,allocatable :: intemp(:,:)
    DATATYPE :: outtemp(numr,wwin%botconfig:wwin%topconfig+1)   !! AUTOMATIC

    if (sparseconfigflag.eq.0) then
       OFLWR "error, must use sparse for parblockconfigmult"; CFLST
    endif

    allocate(intemp(numr,wwin%numconfig));    intemp(:,:)=0d0; outtemp=0

!! transform second to reduce communication?
!!   no, spin transformations done locally now.

    if (wwin%topdfbasis.ge.wwin%botdfbasis) then
       call basis_transformfrom_local(wwin,numr,inavector,&
            intemp(:,wwin%botconfig:wwin%topconfig))
    endif

    call mpiallgather_local(intemp,wwin%numconfig*numr,&
         wwin%nzconfsperproc(:)*numr,wwin%maxconfigsperproc*numr,&
         wwin%NZ_COMM,wwin%nzprocs,wwin%nzrank)

    if (wwin%topconfig.ge.wwin%botconfig) then
       call sparseconfigmult_byproc(1,nprocs,wwin,intemp,outtemp, cptr, sptr, &
            1,1,1,0,0d0,0,1,numr,0,-1)
       if (mshift.ne.0d0) then 
          do ii=wwin%botconfig,wwin%topconfig
             outtemp(:,ii)=outtemp(:,ii)+ intemp(:,ii)*wwin%configmvals(ii)*mshift
          enddo
       endif
    endif

    if (wwin%topdfbasis.ge.wwin%botdfbasis) then
       call basis_transformto_local(wwin,numr,outtemp,outavector)
    endif

    deallocate(intemp)

  end subroutine parblockconfigmult0_gather

#ifdef MPIFLAG

  subroutine parblockconfigmult0_summa(wwin,cptr,sptr,inavector,outavector)
    use fileptrmod
    use r_parameters
    use sparse_parameters
    use ham_parameters
    use mpimod
    use walkmod
    use sparseptrmod
    use configptrmod
    use sparsemultmod
    use basissubmod
    use mpisubmod
    implicit none
    type(walktype),intent(in) :: wwin
    type(CONFIGPTR),intent(in) :: cptr
    type(SPARSEPTR),intent(in) :: sptr
    integer :: ii,iproc,iiproc
    DATATYPE,intent(in) :: inavector(numr,wwin%botdfbasis:wwin%topdfbasis)
    DATATYPE,intent(out) :: outavector(numr,wwin%botdfbasis:wwin%topdfbasis)
    DATATYPE :: intemp(numr,wwin%maxconfigsperproc),&
         outwork(numr,wwin%botconfig:wwin%topconfig+1),&
         outtemp(numr,wwin%botconfig:wwin%topconfig+1)  !! AUTOMATIC

    if (sparseconfigflag.eq.0) then
       OFLWR "error, must use sparse for parblockconfigmult"; CFLST
    endif

    outwork(:,:)=0d0; outtemp=0; intemp=0

    do iiproc=1,wwin%nzprocs

       iproc=wwin%nzproclist(iiproc)

!! transform second to reduce communication?
!!   no, spin transformations done locally now.

       if (wwin%alltopconfigs(iproc).ge.wwin%allbotconfigs(iproc)) then
          if (myrank.eq.iproc) then
             intemp=0d0
             if (wwin%topdfbasis.ge.wwin%botdfbasis) then
                call basis_transformfrom_local(wwin,numr,inavector,intemp)
             endif
          endif
          call mympibcast_local(intemp,iiproc,wwin%nzconfsperproc(iiproc)*numr,wwin%NZ_COMM)

          call sparseconfigmult_byproc(iproc,iproc,wwin,intemp,outtemp, cptr, sptr, &
               1,1,1,0,0d0,0,1,numr,0,-1)

          if (myrank.eq.iproc) then
             if (mshift.ne.0d0) then 
                do ii=wwin%botconfig,wwin%topconfig
                   outtemp(:,ii)=outtemp(:,ii) + &
                        intemp(:,ii-wwin%botconfig+1)*wwin%configmvals(ii)*mshift
                enddo
             endif
          endif
          outwork(:,:)=outwork(:,:) + outtemp(:,:)
       endif
    enddo

    if (wwin%topdfbasis.ge.wwin%botdfbasis) then
       call basis_transformto_local(wwin,numr,outwork,outavector)
    endif

  end subroutine parblockconfigmult0_summa


  subroutine parblockconfigmult0_circ(wwin,cptr,sptr,inavector,outavector)
    use fileptrmod
    use r_parameters
    use sparse_parameters
    use ham_parameters
    use mpimod
    use walkmod
    use sparseptrmod
    use configptrmod
    use sparsemultmod
    use basissubmod
    use mpisubmod
    implicit none
    type(walktype),intent(in) :: wwin
    type(CONFIGPTR),intent(in) :: cptr
    type(SPARSEPTR),intent(in) :: sptr
    integer :: ii,iproc,prevproc,nextproc,deltaproc,iiproc
    DATATYPE,intent(in) :: inavector(numr,wwin%botdfbasis:wwin%topdfbasis)
    DATATYPE,intent(out) :: outavector(numr,wwin%botdfbasis:wwin%topdfbasis)
    DATATYPE :: workvector(numr,wwin%maxconfigsperproc), &
         workvector2(numr,wwin%maxconfigsperproc), &
         outwork(numr,wwin%botconfig:wwin%topconfig+1), &
         outtemp(numr,wwin%botconfig:wwin%topconfig+1)  !! AUTOMATIC

    if (sparseconfigflag.eq.0) then
       OFLWR "error, must use sparse for parblockconfigmult"; CFLST
    endif

!! doing circ mult slightly different than e.g. SINCDVR/coreproject.f90 
!!     and ftcore.f90, holding hands in a circle, prevproc and nextproc, 
!!     each chunk gets passed around the circle

    prevproc=mod(wwin%nzprocs+wwin%nzrank-2,wwin%nzprocs)+1
    nextproc=mod(wwin%nzrank,wwin%nzprocs)+1

    outwork(:,:)=0d0; outtemp=0; workvector=0; workvector2=0

    if (wwin%topdfbasis.ge.wwin%botdfbasis) then
       call basis_transformfrom_local(wwin,numr,inavector,workvector)
       if (mshift.ne.0d0) then 
          do ii=wwin%botconfig,wwin%topconfig
             outwork(:,ii)=outwork(:,ii) + &
                  workvector(:,ii-wwin%botconfig+1)*wwin%configmvals(ii)*mshift
          enddo
       endif
    endif
    
    do deltaproc=0,wwin%nzprocs-1

!! PASSING BACKWARD (plus deltaproc)

       iiproc=mod(wwin%nzrank-1+deltaproc,wwin%nzprocs)+1
       iproc=wwin%nzproclist(iiproc)

       if (wwin%alltopconfigs(iproc).ge.wwin%allbotconfigs(iproc)) then
          call sparseconfigmult_byproc(iproc,iproc,wwin,workvector,outtemp, cptr, sptr, &
               1,1,1,0,0d0,0,1,numr,0,-1)

          outwork(:,:)=outwork(:,:) + outtemp(:,:)
       endif

!! PASSING BACKWARD
!! mympisendrecv(sendbuf,recvbuf,dest,source,...)

       call mympisendrecv_local(workvector,workvector2,prevproc,nextproc,deltaproc,&
            numr * wwin%maxconfigsperproc,wwin%NZ_COMM)
       workvector(:,:)=workvector2(:,:)

    enddo

    if (wwin%topdfbasis.ge.wwin%botdfbasis) then
       call basis_transformto_local(wwin,numr,outwork,outavector)
    endif

  end subroutine parblockconfigmult0_circ

#endif

  subroutine parblockconfigmult(inavector,outavector)
    use r_parameters
    use sparse_parameters
    use configmod
    use xxxmod
    use fileptrmod
    implicit none
    DATATYPE,intent(in) :: inavector(numr,dwwptr%botdfbasis:dwwptr%topdfbasis)
    DATATYPE,intent(out) :: outavector(numr,dwwptr%botdfbasis:dwwptr%topdfbasis)

#ifdef MPIFLAG
    select case (sparsesummaflag)
    case(0)
#endif
       call parblockconfigmult0_gather(dwwptr,yyy%cptr(0),yyy%sptrptr(0),inavector,outavector)
#ifdef MPIFLAG
    case(1)
       call parblockconfigmult0_summa(dwwptr,yyy%cptr(0),yyy%sptrptr(0),inavector,outavector)
    case(2)
       call parblockconfigmult0_circ(dwwptr,yyy%cptr(0),yyy%sptrptr(0),inavector,outavector)
    case default
       OFLWR "Error sparsesummaflag ", sparsesummaflag; CFLST
    end select
#endif

  end subroutine parblockconfigmult

end module parblocklanmod


module lanblockmod
contains

subroutine blocklanczos0( lanblocknum, numout, lansize,maxlansize,order,maxiter, &
     outvectors,outvectorlda, outvalues,inprintflag,guessflag,lancheckmod,lanthresh,&
     multsub,logpar,targetflag,etarget)
  use mpimod
  implicit none 
  logical,intent(in) :: logpar
  integer,intent(in) :: lansize,maxlansize,maxiter,lanblocknum,inprintflag,order,&
       lancheckmod,outvectorlda,numout,targetflag,guessflag
  DATATYPE, intent(in) :: etarget
  DATAECS, intent(out) :: outvalues(numout)  
  DATATYPE, intent(inout) :: outvectors(outvectorlda,numout)
  real*8, intent(in) :: lanthresh
#ifndef MPIFLAG
  integer,parameter :: MPI_COMM_WORLD = (-798)
#endif
  external :: multsub

  call blocklanczos0_local( lanblocknum, numout, lansize,maxlansize,order,maxiter, &
       outvectors,outvectorlda, outvalues,inprintflag,guessflag,lancheckmod,lanthresh,&
       multsub,logpar,targetflag,etarget,MPI_COMM_WORLD)

end subroutine blocklanczos0



!! guessflag=1  use outvectors for initial krylov vectors
!!          =2  do that and also follow, return eigvects with max overlap

subroutine blocklanczos0_local( lanblocknum, numout, lansize,maxlansize,order,maxiter, &
     outvectors,outvectorlda, outvalues,inprintflag,guessflag,lancheckmod,lanthresh,&
     multsub,logpar,targetflag,etarget,IN_COMM)
  use fileptrmod
  use mpisubmod
  use eigenmod
  implicit none 
  logical,intent(in) :: logpar
  integer,intent(in) :: lansize,maxlansize,maxiter,lanblocknum,inprintflag,order,&
       lancheckmod,outvectorlda,numout,targetflag,guessflag,IN_COMM
  DATATYPE, intent(in) :: etarget
  DATAECS, intent(out) :: outvalues(numout)  
  DATATYPE, intent(inout) :: outvectors(outvectorlda,numout)
  real*8, intent(in) :: lanthresh
  external :: multsub
  integer :: printflag, iorder,k,flag,j,id,nn,i,thislanblocknum, thisdim,thisout
  real*8 :: error(numout), stopsum,rsum,nextran
  DATATYPE :: alpha(lanblocknum,lanblocknum),beta(lanblocknum,lanblocknum), &
       nullvector1(100)=0, nullvector2(100)=0, lastvalue(numout), thisvalue(numout), &
       valdot(numout),normsq(numout),sqdot(numout)
!! made these allocatable to fix lawrencium segfault 04-15
  DATATYPE, allocatable :: lanham(:,:,:,:), laneigvects(:,:,:), templanham(:,:),&
       betas(:,:,:),betastr(:,:,:),       invector(:), multvectors(:,:), &
       lanvects(:,:,:), tempvectors(:,:),  lanmultvects(:,:,:), tempvectors2(:,:), &
       followvects(:,:), followlandots(:,:,:)
  DATAECS, allocatable :: laneigvals(:)
  
  if (numout.lt.lanblocknum) then
     OFLWR "numout >= lanblocknum please",numout,lanblocknum; CFLST
  endif
  if (numout.gt.order*lanblocknum) then
     OFLWR "numout gt order*lanblocknum is impossible.",numout,order,lanblocknum; CFLST
  endif

  allocate(&
       invector(maxlansize), multvectors(maxlansize,lanblocknum), &
       lanvects(maxlansize,lanblocknum,order), tempvectors(maxlansize,numout), &
       lanmultvects(maxlansize,lanblocknum,order), tempvectors2(maxlansize,numout))

  invector=0d0; multvectors=0d0; lanvects=0d0; tempvectors=0d0; 
  lanmultvects=0d0; tempvectors2=0d0

  if (guessflag.eq.2) then
     allocate(followvects(maxlansize,numout),followlandots(numout,lanblocknum,order))
  else
     allocate(followvects(1,numout),followlandots(1,1,order))
  endif
  followvects=0; followlandots=0;

  allocate( &
       lanham(lanblocknum,order,lanblocknum,order),&
       laneigvects(lanblocknum,order,order*lanblocknum), &
       laneigvals(order*lanblocknum), &
       templanham(lanblocknum*order,lanblocknum*order))
  lanham=0d0; laneigvects=0d0; laneigvals=0d0; templanham=0d0

  printflag=inprintflag

  lastvalue(:)=1d10;  thisvalue(:)=1d9;   alpha=0; beta=0;

  call rand_init(1731.d0) 

  if (guessflag==0) then
     outvectors=0.0d0
     do k=1,numout
        do i=1,lansize
           outvectors(i,k)=nextran() + (0d0,1d0) * nextran()
        enddo
     enddo
  endif

!! follow vectors
  if (guessflag.eq.2) then
     followvects(1:lansize,:) = outvectors(1:lansize,:)
  endif

  outvectors(lansize+1:,:)=0d0  ! why not

  flag=0
  do while (flag.ne.1)
     flag=0

     do j=1,numout
        if (lansize.eq.0) then
           call multsub(nullvector1(:),nullvector2(:))
           call multsub(nullvector1(:),nullvector2(:))
        else
           call multsub(outvectors(:,j),tempvectors(:,j))
           call multsub(tempvectors(:,j),tempvectors2(:,j))
        endif
     enddo

     error(:)=1000d0

     if (lansize.eq.0) then
        call vecnulldots(numout,normsq,logpar)
        call vecnulldots(numout,valdot,logpar)
        call vecnulldots(numout,sqdot,logpar)
     else
        call vechdots(outvectors(:,:),outvectors(:,:),lansize,outvectorlda,&
             outvectorlda,numout,normsq,logpar)
        call vechdots(outvectors(:,:),tempvectors(:,:),lansize,outvectorlda,&
             maxlansize,numout,valdot,logpar)
        call vechdots(outvectors(:,:),tempvectors2(:,:),lansize,outvectorlda,&
             maxlansize,numout,sqdot,logpar)
     endif

     do j=1,numout
        error(j)=abs(&
             valdot(j)**2 / normsq(j)**2 - &
             sqdot(j)/normsq(j))
     enddo
     if (printflag.ne.0) then
        OFL; write(mpifileptr,'(A10,1000E8.1)') " FIRST ERRORS ", error(1:numout); CFL
     endif

     stopsum=0d0
     do nn=1,numout
        if (error(nn).gt.stopsum) then
           stopsum=error(nn)
        endif
     enddo
     if (stopsum.lt.lanthresh) then
        OFLWR "Vectors converged on read",stopsum,lanthresh; CFL
        flag=1
     else
        if (printflag.ne.0) then
           OFLWR "MAX ERROR : ", stopsum, lanthresh; CFL
        endif
     endif

     if (lansize.gt.0) then
        lanvects(1:lansize,:,1)=outvectors(1:lansize,1:lanblocknum)
     endif

     if (guessflag.eq.2) then
        if (lansize.eq.0) then
           call allnulldots(lanblocknum,lanblocknum,followlandots(:,:,1),logpar)
        else
           call allhdots(followvects(:,:),lanvects(:,:,1),lansize,&
                maxlansize,numout,lanblocknum,followlandots(:,:,1),logpar)
        endif
     endif

     lanham(:,:,:,:)=0.0d0
     do i=1,lanblocknum
        if (lansize.eq.0) then
           call nullgramschmidt_fast(i-1, logpar)
        else
           call myhgramschmidt_fast(lansize, i-1, maxlansize, &
                lanvects(:,:,1),lanvects(:,i,1),logpar)
        endif
     enddo
     do j=1,lanblocknum
        if (lansize.eq.0) then
           call multsub(nullvector1(:),nullvector2(:))
        else
           call multsub(lanvects(:,j,1),multvectors(:,j))
        endif
     enddo

     lanmultvects(:,:,1)=multvectors(:,:)

     if (lansize.eq.0) then
        call allnulldots(lanblocknum,lanblocknum,alpha,logpar)
     else
        call allhdots(lanvects(:,:,1),multvectors(:,:),lansize,&
             maxlansize,lanblocknum,lanblocknum,alpha,logpar)
     endif

     lanham(:,1,:,1)= alpha(:,:)

     if (printflag.ne.0) then
        OFLWR "FIRST ALPHA ", alpha(1,1); CFL
     endif


     iorder=1
     do while ( flag==0 )

        iorder=iorder+1
        if (iorder.eq.order) then
           flag=2
        endif
! each loop: multvector is h onto previous normalized vector

        if (lansize.gt.0) then
           lanvects(1:lansize,:,iorder)=multvectors(1:lansize,:)
        endif

        thislanblocknum=lanblocknum
        if (iorder*lanblocknum.ge.maxiter) then
           OFLWR "Max iters ",maxiter," reached, setting converged flag to 2"; CFL
           flag=2 
           thislanblocknum=maxiter-(iorder-1)*lanblocknum
        endif

        do i=1,thislanblocknum
           if (lansize.eq.0) then
              call nullgramschmidt_fast((iorder-1)*lanblocknum+i-1,logpar)
           else
              call myhgramschmidt_fast(lansize, (iorder-1)*lanblocknum+i-1, &
                   maxlansize, lanvects,lanvects(:,i,iorder),logpar)
           endif
        enddo

        if (guessflag.eq.2) then
           if (lansize.eq.0) then
              call allnulldots(lanblocknum,lanblocknum,followlandots(:,:,iorder),logpar)
           else
              call allhdots(followvects(:,:),lanvects(:,:,iorder),lansize,&
                   maxlansize,numout,lanblocknum,followlandots(:,:,iorder),logpar)
           endif
        endif

        multvectors(:,:)=0d0
        if (lansize.eq.0) then
           do j=1,thislanblocknum
              call multsub(nullvector2(:),nullvector1(:))
           enddo
        else
           do j=1,thislanblocknum
              call multsub(lanvects(:,j,iorder),multvectors(:,j))
           enddo
        endif
        lanmultvects(:,:,iorder)=multvectors(:,:)
        
        if (lansize.eq.0) then
           call allnulldots(lanblocknum,lanblocknum,alpha,logpar)
        else
           call allhdots(lanvects(:,:,iorder),multvectors(:,:),lansize,&
                maxlansize,lanblocknum,lanblocknum,alpha,logpar)
        endif
        lanham(:,iorder,:,iorder)=alpha(:,:)

        allocate(betas(lanblocknum,iorder-1,lanblocknum), &
             betastr(lanblocknum,lanblocknum,iorder-1))

        if (lansize.eq.0) then
           call allnulldots(lanblocknum*(iorder-1),lanblocknum,betas(:,:,:),logpar)
           call allnulldots(lanblocknum,lanblocknum*(iorder-1),betastr(:,:,:),logpar)
        else
           call allhdots(lanvects(:,:,:),lanmultvects(:,:,iorder),lansize,maxlansize,&
                lanblocknum*(iorder-1),lanblocknum,betas(:,:,:),logpar)
           call allhdots(lanvects(:,:,iorder),lanmultvects(:,:,:),lansize,maxlansize,&
                lanblocknum,lanblocknum*(iorder-1),betastr(:,:,:),logpar)
        endif

!!$        call allhdots(lanvects(:,:,iorder),lanmultvects(:,:,iorder-1),lansize,&
!!$            maxlansize,lanblocknum,lanblocknum,betastr(:,:,iorder-1),logpar)

        do i=1,iorder-1
           lanham(:,i,:,iorder)=betas(:,i,:)
        enddo
        do i=1,iorder-1                !!$  not needed on paper (arnoldi), just iorder-1
!!$        do i=iorder-1,iorder-1
           lanham(:,iorder,:,i)=betastr(:,:,i)
        enddo
        deallocate(betas,betastr)

!!! flag.ne.0 for maxiter, iorder
        if (mod(iorder,lancheckmod).eq.0.or.flag.ne.0) then

           templanham(:,:)=0d0
           templanham(1:lanblocknum*iorder,1:lanblocknum*iorder)=&
                RESHAPE(lanham(:,1:iorder,:,1:iorder),&
                (/lanblocknum*iorder,lanblocknum*iorder/))
           thisdim=min(maxiter,iorder*lanblocknum)
           call CONFIGEIG(templanham,thisdim,order*lanblocknum,laneigvects,laneigvals)
           if (targetflag.ne.0) then
              call mysort2(laneigvects,laneigvals,thisdim,order*lanblocknum)
           endif

           if (guessflag.eq.2) then
              call followsort(thisdim,followlandots)
           endif

           thisout=min(numout,thisdim)

           outvalues(1:thisout)=laneigvals(1:thisout)

           if (printflag.ne.0) then
              OFL; write(mpifileptr,'(A10,1000F19.12)') " ENERGIES ", outvalues(1:thisout); CFL
           endif

           thisvalue(1:thisout)=outvalues(1:thisout)
           stopsum=0d0
           do nn=1,thisout
              rsum=abs(thisvalue(nn)-lastvalue(nn))
              if (rsum.gt.stopsum) then
                 stopsum=rsum
              endif
           enddo
           lastvalue(1:thisout)=thisvalue(1:thisout)

! so lanthresh is for error of HPsi...  empirically make guess (divided by 4) as to 
!  when it might be converged based on change in energy

!!! flag=1 for maxiter,maxorder
           if ((thisout.eq.numout.and.stopsum.lt.lanthresh/4).or.flag.ne.0) then

              if (printflag.ne.0) then
                 OFL; write(mpifileptr,'(A25,2E12.5,I10)')  "checking convergence.",&
                      stopsum,lanthresh/10,thisdim; CFL
              endif

              outvectors = 0.0d0
              if (lansize.gt.0) then

                 do  j=1, numout
                    do k=1, iorder
                       do id=1,lanblocknum
                          if ((k-1)*lanblocknum+id.le.maxiter) then
                             outvectors(1:lansize,j) = outvectors(1:lansize,j) + &
                                  laneigvects(id,k,j) * lanvects(1:lansize,id,k)
                          endif
                       enddo
                    enddo
                 enddo

              endif  !! lansize 

              if (lansize.eq.0) then
                 call vecnulldots(numout,normsq,logpar)
              else
                 call vecthisdots(outvectors(:,:),outvectors(:,:),lansize,&
                      outvectorlda,outvectorlda,numout,normsq,logpar)
              endif
              do  j=1, numout
                 outvectors(:,j)=outvectors(:,j)/sqrt(normsq(j))
              enddo

              tempvectors=0d0
              do j=1,numout
                 if (lansize.eq.0) then
                    call multsub(nullvector2(:),nullvector1(:))
                    call multsub(nullvector1(:),nullvector2(:))
                 else
                    call multsub(outvectors(:,j),tempvectors(:,j))
                    call multsub(tempvectors(:,j),tempvectors2(:,j))
                 endif
              enddo

              if (lansize.eq.0) then
                 call vecnulldots(numout,normsq,logpar)
                 call vecnulldots(numout,valdot,logpar)
                 call vecnulldots(numout,sqdot,logpar)
              else
                 call vechdots(outvectors(:,:),outvectors(:,:),lansize,outvectorlda,&
                      outvectorlda,numout,normsq,logpar)
                 call vechdots(outvectors(:,:),tempvectors(:,:),lansize,outvectorlda,&
                      maxlansize,numout,valdot,logpar)
                 call vechdots(outvectors(:,:),tempvectors2(:,:),lansize,outvectorlda,&
                      maxlansize,numout,sqdot,logpar)
              endif

              do j=1,numout
                 error(j)=abs(&
                      valdot(j)**2 / normsq(j)**2 - &
                      sqdot(j)/normsq(j))
              enddo
              if (printflag.ne.0) then
                 OFL; write(mpifileptr,'(A10,1000E9.2)') " ERRORS ", error(1:numout); CFL
              endif
              stopsum=0d0
              do nn=1,numout
                 if (error(nn).gt.stopsum) then
                    stopsum=error(nn)
                 endif
              enddo
              if (stopsum.gt.lanthresh) then
                 if (printflag.ne.0) then
                    OFLWR "Not converged", stopsum, lanthresh; CFL
                 endif
                 if (thisdim.ge.maxiter) then
                    OFLWR "MAX DIM REACHED, NO CONVERGENCE -- THRESHOLDS TOO HIGH? BUG?",&
                         stopsum,lanthresh; CFLST
                 endif
              else
                 flag=1
              endif
           endif
        endif
     enddo

     if (flag==1) then
        if (printflag.ne.0) then
           OFLWR "Converged. ",stopsum,lanthresh  !!, "HDOTS:"
!!$           do i=1,numout
!!$              write(mpifileptr,'(10000F8.3)') &
!!$                  (hdot(outvectors(:,i),outvectors(:,j),lansize,logpar),j=1,numout)
!!$           enddo
           CFL
        endif
     else
        flag=0
        OFLWR "  Not converged. restarting.", stopsum,lanthresh; CFL
     endif
  enddo

  deallocate(  invector,multvectors, lanvects, tempvectors, lanmultvects,tempvectors2)
  deallocate(followvects, followlandots )
  deallocate( lanham,       laneigvects,       laneigvals,       templanham)

contains

  function nulldot(logpar)
    implicit none
    logical,intent(in) :: logpar
    DATATYPE :: nulldot
    if (logpar) then
       nulldot=0d0
       call mympireduceone_local(nulldot,IN_COMM)
    else
       OFLWR "WHAT? nulldot called but not doing parallel calc... dimension zero???"; CFLST
    endif
  end function nulldot

  function hdot(bra,ket,n,logpar)
    implicit none
    integer,intent(in) :: n
    logical,intent(in) :: logpar
    DATATYPE,intent(in) :: bra(n),ket(n)
    DATATYPE :: hdot
    hdot=DOT_PRODUCT(bra,ket)
    if (logpar) then
       call mympireduceone_local(hdot,IN_COMM)
    endif
  end function hdot

  function thisdot(bra,ket,n,logpar)
    implicit none
    integer,intent(in) :: n
    logical,intent(in) :: logpar
    DATATYPE,intent(in) :: bra(n), ket(n)
    DATATYPE :: thisdot,csum,csum2
    integer :: ii
    csum2=0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,csum)
    csum=0
!$OMP DO SCHEDULE(DYNAMIC)
    do ii=1,n
       csum=csum+CONJUGATE(bra(ii))*ket(ii)
    enddo
!$OMP END DO
!$OMP CRITICAL
    csum2=csum2+csum
!$OMP END CRITICAL
!$OMP END PARALLEL
    if (logpar) then
       call mympireduceone_local(csum2,IN_COMM)
    endif
    thisdot=csum2
  end function thisdot


  subroutine allhdots(bravectors,ketvectors,n,lda,num1,num2,outdots,logpar)
    implicit none
    integer,intent(in) :: num1,num2,lda,n
    logical,intent(in) :: logpar
    integer :: id,jd

    DATATYPE,intent(in) :: bravectors(lda,num1), ketvectors(lda,num2)
    DATATYPE,intent(out) :: outdots(num1,num2)

    if (n.le.0) then
       print *, "progerrrr"; stop
    endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(id,jd)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do id=1,num1
       do jd=1,num2
          outdots(id,jd)= DOT_PRODUCT(bravectors(1:n,id),ketvectors(1:n,jd))
       enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL
    if (logpar) then
       call mympireduce_local(outdots,num1*num2,IN_COMM)
    endif
  end subroutine allhdots

  subroutine allnulldots(num1,num2,outdots,logpar)
    implicit none
    integer,intent(in) :: num1,num2
    logical,intent(in) :: logpar
    DATATYPE,intent(out) :: outdots(num1,num2)
    
    if (logpar) then
       outdots(:,:)=0d0
       call mympireduce_local(outdots,num1*num2,IN_COMM)
    else
       OFLWR "WHAT? allnulldots called but not doing parallel calc... dimension zero???"; CFLST
    endif
  end subroutine allnulldots


  subroutine vechdots(bravectors,ketvectors,n,ldabra,ldaket,num,outdots,logpar)
    implicit none
    integer,intent(in) :: num,ldabra,ldaket,n
    logical,intent(in) :: logpar
    DATATYPE,intent(in) :: bravectors(ldabra,num), ketvectors(ldaket,num)
    DATATYPE,intent(out) :: outdots(num)
    integer :: id

    if (n.le.0) then
       print *, "progerrrrOR"; stop
    endif
    
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(id)
!$OMP DO SCHEDULE(STATIC)
    do id=1,num
       outdots(id)= DOT_PRODUCT(bravectors(1:n,id),ketvectors(1:n,id))
    enddo
!$OMP END DO
!$OMP END PARALLEL
    if (logpar) then
       call mympireduce_local(outdots,num,IN_COMM)
    endif
  end subroutine vechdots


  subroutine vecthisdots(bravectors,ketvectors,n,ldabra,ldaket,num,outdots,logpar)
    implicit none
    integer,intent(in) :: num,ldabra,ldaket,n
    logical,intent(in) :: logpar
    DATATYPE,intent(in) :: bravectors(ldabra,num), ketvectors(ldaket,num)
    DATATYPE,intent(out) :: outdots(num)
    DATATYPE :: tempdots(num)
    integer :: ii,id

    if (n.le.0) then
       print *, "progerrrrxxx"; stop
    endif

    outdots(:)=0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,id,tempdots)
    tempdots(:)=0
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2)
    do id=1,num
       do ii=1,n
          tempdots(id)=tempdots(id)+CONJUGATE(bravectors(ii,id))*ketvectors(ii,id)
       enddo
    enddo
!$OMP END DO
!$OMP CRITICAL
    outdots(:)=outdots(:)+tempdots(:)
!$OMP END CRITICAL
!$OMP END PARALLEL
    if (logpar) then
       call mympireduce_local(outdots,num,IN_COMM)
    endif
  end subroutine vecthisdots


  subroutine vecnulldots(num,outdots,logpar)
    implicit none
    integer,intent(in) :: num
    logical,intent(in) :: logpar
    DATATYPE,intent(out) :: outdots(num)

    if (logpar) then
       outdots(:)=0d0
       call mympireduce_local(outdots,num,IN_COMM)
    else
       OFLWR "WHAT? vecnulldots called but not doing parallel calc... dimension zero???"; CFLST
    endif
  end subroutine vecnulldots


! n is the length of the vectors; m is how many to orthogonalize to
  subroutine myhgramschmidt_fast(n, m, lda, previous, vector,logpar)
    implicit none
    integer,intent(in) :: n,m,lda
    DATATYPE,intent(in) :: previous(lda,m)
    DATATYPE,intent(inout) :: vector(n)
    DATATYPE :: norm,myhdots(m)
    logical,intent(in) :: logpar
    integer :: i,j

    if (n.le.0) then
       print *, "progerrrreeee"; stop
    endif

    do j=1,2

!!$ separate dot products for less numerical error; allhdots for speed

       if (m.ne.0) then
          call allhdots(previous(:,:),vector,n,lda,m,1,myhdots,logpar)
       endif

       do i=1,m
!!$        vector=vector-previous(1:n,i)* hdot(previous(1:n,i),vector(1:n),n,logpar) 
!!$  only the same with previous perfectly orth
          vector=vector-previous(1:n,i)*myhdots(i)                            
       enddo
       norm=sqrt(hdot(vector,vector,n,logpar))
       vector=vector/norm
       if (abs(norm).lt.1e-7) then
          OFLWR "Gram schmidt norm",norm,m; CFL
       endif
    enddo
  end subroutine myhgramschmidt_fast

  subroutine nullgramschmidt_fast(m,logpar)
    implicit none
    integer,intent(in) :: m
    logical,intent(in) :: logpar
    DATATYPE :: norm,mynulldots(m)
    integer :: j

    do j=1,2

       if (m.ne.0) then
          call allnulldots(m,1,mynulldots,logpar)
       endif

!!$     do i=1,m
!!$        norm=nulldot(logpar)
!!$     enddo

       norm=sqrt(nulldot(logpar))
    enddo
  end subroutine nullgramschmidt_fast


!! DATAECS data type for values
  subroutine mysort2(inout, values,n,lda)
    implicit none
    integer,intent(in) :: lda, n
    DATATYPE,intent(inout) :: inout(lda,n)
    DATAECS,intent(inout) :: values(lda)
    DATAECS,allocatable ::  newvals(:)
    DATATYPE,allocatable :: out(:,:)
    integer,allocatable :: taken(:), order(:)
    real*8 :: lowval 
    integer :: i,j,whichlowest, flag

    allocate(out(lda,n), newvals(lda), taken(lda), order(lda))

    taken=0;  order=-1
    do j=1,n
       whichlowest=-1; flag=0;     lowval=1.d+30  !! is not used (see flag)
       do i=1,n
          if ( taken(i) .eq. 0 ) then
             if ((flag.eq.0) .or.(abs(values(i)-etarget) .le. lowval)) then
                flag=1;              lowval=abs(values(i)-etarget); whichlowest=i
             endif
          endif
       enddo
       if ((whichlowest.gt.n).or.(whichlowest.lt.1)) then
          write(*,*) taken;        write(*,*) 
          write(*,*) "WHICHLOWEST ERROR, J=",j," WHICHLOWEST=", whichlowest
          call mpistop()
       endif
       if (taken(whichlowest).ne.0) then
          write(*,*) "TAKEN ERROR.";        call mpistop()
       endif
       taken(whichlowest)=1;     order(j)=whichlowest
       out(:,j)=inout(:,order(j));     newvals(j)=values(order(j))
    enddo

    values(1:n)=newvals(1:n)
    inout(1:n,1:n)=out(1:n,1:n)

    deallocate(out, newvals, taken, order)

  end subroutine mysort2

  subroutine followsort(indim,inlandots)
    implicit none
    integer,intent(in) :: indim
    DATATYPE,intent(in) :: inlandots(numout,indim)
    DATATYPE,allocatable :: followdots(:,:), tempeigvect(:,:)
    DATAECS :: tempeigval
    real*8 :: rmaxdot
    integer :: i,j,ifollow

!!    OFLWR "Blocklan ifollow "; CFL

    allocate(followdots(numout,indim), tempeigvect(lanblocknum,order))

    followdots=0; tempeigvect=0; 

    call MYGEMM('N','N',numout,indim,indim,DATAONE,inlandots(:,:),numout, &
         laneigvects(:,:,:),order*lanblocknum,DATAZERO,followdots(:,:),numout)

    do i=1,numout
       ifollow=(-1)
       rmaxdot=0d0
       do j=i,indim
          if (abs(followdots(i,j)).gt.rmaxdot) &
!!            do j=i, resorting, not needed:      .and.followtaken(j).eq.0) then
               then
             rmaxdot=abs(followdots(i,j))
             ifollow=j
          endif
       enddo

       if (ifollow.ne.i) then
          tempeigvect(:,:) = laneigvects(:,:,i)
          laneigvects(:,:,i) = laneigvects(:,:,ifollow)
          laneigvects(:,:,ifollow) = tempeigvect(:,:)
          
          tempeigval = laneigvals(i)
          laneigvals(i) = laneigvals(ifollow)
          laneigvals(ifollow) = tempeigval
       endif

!!       OFLWR "FF ",i,ifollow; CFL

    enddo

    deallocate(followdots, tempeigvect)

  end subroutine followsort

end subroutine blocklanczos0_local



subroutine blocklanczos(order,outvectors, outvalues,inprintflag,guessflag)
  use r_parameters
  use sparse_parameters
  use lanparameters
  use fileptrmod
  use configmod
  use parblocklanmod
  use basissubmod
  use mpisubmod
  implicit none 
  integer, intent(in) :: order,inprintflag,guessflag
  integer :: printflag,maxdim,vdim,ii
  DATAECS, intent(out) :: outvalues(order)
  DATATYPE, intent(inout) :: outvectors(numr,www%firstconfig:www%lastconfig,order)
  DATATYPE,allocatable :: workvectorsspin(:,:,:), shufflevectors(:,:,:)
  DATATYPE :: nullvectors(numr,order)

  printflag=max(lanprintflag,inprintflag);    nullvectors=0

  if (sparseconfigflag.eq.0) then
     OFLWR "error, can't use blocklanczos for a-vector with sparseconfigflag=0"; CFLST
  endif

  allocate(workvectorsspin(numr,www%botdfbasis:www%topdfbasis,order))

  if (www%topdfbasis-www%botdfbasis+1.ne.0) then
     workvectorsspin(:,:,:)=0d0
     if (guessflag.ne.0) then
        do ii=1,order
           call basis_transformto_local(www,numr,outvectors(:,www%botconfig,ii),&
                workvectorsspin(:,:,ii))
        enddo
     endif
  endif

  maxdim=dwwptr%numdfbasis*numr
  vdim=(dwwptr%topdfbasis-dwwptr%botdfbasis+1)*numr

  call mpibarrier()
  if (use_dfwalktype.and.shuffle_dfwalktype) then
     allocate(shufflevectors(numr,dwwptr%botdfbasis:dwwptr%topdfbasis,order))
     if (dwwptr%topdfbasis.ge.dwwptr%botdfbasis) then
        shufflevectors=0
     endif
     call basis_shuffle_several(numr,www,workvectorsspin,dwwptr,shufflevectors,order)
  else
     allocate(shufflevectors(numr,1,order))
  endif

  if (use_dfwalktype.and.shuffle_dfwalktype) then
     if (nzflag.eq.0.or.dwwptr%nzrank.gt.0) then
        call blocklanczos0_local(order,order,vdim,vdim,lanczosorder,maxdim,shufflevectors,vdim,&
             outvalues,printflag,guessflag,lancheckstep,lanthresh,parblockconfigmult,&
             .true.,0,DATAZERO,dwwptr%NZ_COMM)
     else
        if (dwwptr%topdfbasis.ge.dwwptr%botdfbasis) then
           shufflevectors=0
        endif
     endif
  else
     if (nzflag.eq.0.or.dwwptr%nzrank.gt.0) then
        call blocklanczos0_local(order,order,vdim,vdim,lanczosorder,maxdim,workvectorsspin,vdim,&
             outvalues,printflag,guessflag,lancheckstep,lanthresh,parblockconfigmult,&
             .true.,0,DATAZERO,dwwptr%NZ_COMM)
     else
        if (dwwptr%topdfbasis.ge.dwwptr%botdfbasis) then
           workvectorsspin=0d0
        endif
     endif
  endif

  call mpibarrier()

  if (use_dfwalktype.and.shuffle_dfwalktype) then
     call basis_shuffle_several(numr,dwwptr,shufflevectors,www,workvectorsspin,order)
  endif
  deallocate(shufflevectors)

!  if (printflag.ne.0) then
!     OFLWR "Done with blocklanczos0_local."; CFL
!  endif

  if (www%lastconfig.ge.www%firstconfig) then
     outvectors(:,:,:)=0d0
  endif

  if (www%topdfbasis-www%botdfbasis+1.ne.0) then
     do ii=1,order
        call basis_transformfrom_local(www,numr,workvectorsspin(:,:,ii),&
             outvectors(:,www%botconfig,ii))
     enddo
  endif

  deallocate(workvectorsspin)

  if (www%parconsplit.eq.0) then
     do ii=1,order
        call mpiallgather(outvectors(:,:,ii),www%numconfig*numr,&
             www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
     enddo
  endif

end subroutine blocklanczos

end module lanblockmod

