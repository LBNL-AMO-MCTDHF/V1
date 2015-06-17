
#include "Definitions.INC"

subroutine blocklanczos( order,outvectors, outvalues,inprintflag,guessflag)
  use parameters
  implicit none 
  integer, intent(in) :: order,inprintflag,guessflag
  integer :: printflag,maxdim,vdim,i,mytop,mybot, calcsize
  DATATYPE, intent(out) :: outvalues(order)
  DATATYPE :: outvectors(numconfig,numr,order)
  DATATYPE, allocatable :: tempoutvectorstr(:,:,:),outvectorstr(:,:,:), &
       initvectors(:,:,:), outvectorsspin(:,:,:)
  external :: parblockconfigmult_transpose, parblockconfigmult_transpose_spin

  printflag=max(lanprintflag,inprintflag)

  if (allspinproject.ne.0) then
     calcsize=spintotrank;   mytop=spinend; mybot=spinstart
  else
     calcsize=numconfig;    mytop=topwalk; mybot=botwalk
  endif
  if (dfrestrictflag.eq.0) then
     maxdim=calcsize*numr
  else
     if (allspinproject.ne.0) then
        maxdim=spintotdfrank*numr
     else
        maxdim=numdfconfigs*numr
     endif
  endif

  vdim=(mytop-mybot+1)*numr

  allocate(tempoutvectorstr(numr,mybot:mytop,order),outvectorstr(numr,order,calcsize), &
       initvectors(calcsize,numr,order),        outvectorsspin(calcsize,numr,order))

  tempoutvectorstr(:,:,:)=0d0

  if (guessflag.ne.0) then
     if (allspinproject.eq.0) then
        initvectors(:,:,:)=outvectors(:,:,:)  
     else
        call configspin_transformto(numr*order,outvectors,initvectors)
     endif
     do i=1,order
        tempoutvectorstr(:,:,i)=TRANSPOSE(initvectors(mybot:mytop,:,i))
     enddo
  endif
  if (allspinproject.eq.0) then
     call blocklanczos0(order,order,vdim,vdim,lanczosorder,maxdim,tempoutvectorstr,vdim,outvalues,printflag,guessflag,lancheckstep,lanthresh,parblockconfigmult_transpose,.true.,0,DATAZERO)
  else
     call blocklanczos0(order,order,vdim,vdim,lanczosorder,maxdim,tempoutvectorstr,vdim,outvalues,printflag,guessflag,lancheckstep,lanthresh,parblockconfigmult_transpose_spin,.true.,0,DATAZERO)
  endif


  outvectorstr(:,:,:)=0d0
  do i=1,order
     outvectorstr(:,i,mybot:mytop)=tempoutvectorstr(:,:,i)
  enddo
  if (allspinproject.eq.0) then
     call mpiallgather(outvectorstr,calcsize*numr*order,configsperproc*numr*order,&
          maxconfigsperproc*numr*order)
     do i=1,order
        outvectors(:,:,i)=TRANSPOSE(outvectorstr(:,i,:))
     enddo
  else
     call mpiallgather(outvectorstr,calcsize*numr*order,allspinranks*numr*order,maxspinrank*numr*order)
     do i=1,order
        outvectorsspin(:,:,i)=TRANSPOSE(outvectorstr(:,i,:))
     enddo
     outvectors(:,:,:)=0d0
     call configspin_transformfrom(numr*order,outvectorsspin,outvectors)
  endif

  deallocate(tempoutvectorstr,outvectorstr, initvectors,       outvectorsspin)

end subroutine blocklanczos

function nulldot(logpar)
  use parameters
  implicit none
  logical :: logpar
  DATATYPE :: nulldot
  if (logpar) then
     nulldot=0d0
     call mympireduceone(nulldot)
  else
     OFLWR "WHAT? nulldot called but not doing parallel calc... dimension zero???"; CFLST
  endif
end function nulldot

function hdot(in,out,n,logpar)
  use parameters
  implicit none
  integer :: n
  logical :: logpar
  DATATYPE :: in(n),out(n),hdot
  hdot=DOT_PRODUCT(in,out)
  if (logpar) then
     call mympireduceone(hdot)
  endif
end function hdot

function thisdot(in,out,n,logpar)
  use parameters
  implicit none
  integer :: n
  logical :: logpar
  DATATYPE :: in(n),out(n),thisdot,dot
  thisdot=dot(in,out,n)
  if (logpar) then
     call mympireduceone(thisdot)
  endif
end function thisdot

subroutine allhdots(bravectors,ketvectors,n,lda,num1,num2,outdots,logpar)
  use parameters
  implicit none
  integer :: id,jd,num1,num2,lda,n
  logical :: logpar
  DATATYPE :: bravectors(lda,num1), ketvectors(lda,num2), outdots(num1,num2)
  do id=1,num1
     do jd=1,num2
        outdots(id,jd)= DOT_PRODUCT(bravectors(1:n,id),ketvectors(1:n,jd))
     enddo
  enddo
  if (logpar) then
     call mympireduce(outdots,num1*num2)
  endif
end subroutine allhdots

subroutine nulldots(num1,num2,outdots,logpar)
  use parameters
  implicit none
  integer :: num1,num2
  logical :: logpar
  DATATYPE :: outdots(num1,num2)
  if (logpar) then
     outdots(:,:)=0d0
     call mympireduce(outdots,num1*num2)
  else
     OFLWR "WHAT? nulldots called but not doing parallel calc... dimension zero???"; CFLST
  endif
end subroutine nulldots

! n is the length of the vectors; m is how many to orthogonalize to
subroutine myhgramschmidt_fast(n, m, lda, previous, vector,logpar)
  use fileptrmod
  implicit none
  integer :: n,m,lda, i,j
  DATATYPE :: previous(lda,m), vector(n),hdot, myhdots(m), norm
  logical :: logpar
  do j=1,2
     if (m.ne.0) then
        call allhdots(previous(:,:),vector,n,lda,m,1,myhdots,logpar)
     endif
     do i=1,m
!!$?       vector=vector-previous(1:n,i)* hdot(previous(1:n,i),vector,n) 
!!$? only the same with previous perfectly orth
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
  use fileptrmod
  implicit none
  integer :: m,j
  DATATYPE :: myhdots(m),norm,nulldot
  logical :: logpar
  do j=1,2
     if (m.ne.0) then
        call nulldots(m,1,myhdots,logpar)
     endif
     norm=sqrt(nulldot(logpar))
  enddo
end subroutine nullgramschmidt_fast


!! TAKES TRANSPOSES AS INPUT AND OUTPUT

recursive subroutine parblockconfigmult_transpose(inavectortr,outavectortr)
  use parameters
  use mpimod
  use xxxmod
  implicit none
  integer :: ir
  DATATYPE :: inavectortr(numr,botwalk:topwalk), outavectortr(numr,botwalk:topwalk)
  DATATYPE :: intemptr(numr,numconfig), ttvector(numconfig,numr), ttvector2(botwalk:topwalk,numr)

!!$  integer,save :: allochere=0
!!$  DATATYPE,save,allocatable :: intemptr(:,:), ttvector(:,:), ttvector2(:,:)
!!$  if (allochere.eq.0) then
!!$     allocate(intemptr(numr,numconfig), ttvector(numconfig,numr), ttvector2(botwalk:topwalk,numr))
!!$  endif
!!$  allochere=1

  if (sparseconfigflag.eq.0) then
     OFLWR "error, must use sparse for parblockconfigmult_transpose"; CFLST
  endif

  intemptr(:,:)=0d0;   intemptr(:,botwalk:topwalk)=inavectortr(:,:)
  call mpiallgather(intemptr,numconfig*numr,configsperproc*numr,maxconfigsperproc*numr)
  ttvector(:,:)=TRANSPOSE(intemptr(:,:))
  if (dfrestrictflag.ne.0) then
     call dfrestrict(ttvector,numr)
  endif

  call sparseconfigmult_nompi(ttvector,ttvector2, yyy%cptr(0), yyy%sptr(0), 1,1,1,0,0d0)
  if (mshift.ne.0d0) then 
     do ir=1,numr
        ttvector2(:,ir)=ttvector2(:,ir)+ ttvector(botwalk:topwalk,ir)*configmvals(botwalk:topwalk)*mshift
     enddo
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict_par(ttvector2,numr)
  endif

  outavectortr(:,:)=TRANSPOSE(ttvector2(:,:))

end subroutine parblockconfigmult_transpose


recursive subroutine parblockconfigmult_transpose_spin(inavectortrspin,outavectortrspin)
  use parameters
  use mpimod
  use xxxmod
  implicit none
  integer :: ir
  DATATYPE :: inavectortrspin(numr,spinrank), outavectortrspin(numr,spinrank)
  DATATYPE :: intemptr(numr,spintotrank), ttvector(numconfig,numr), ttvector2(botwalk:topwalk,numr), &
       ttvectorspin(spintotrank,numr), ttvector2spin(spinrank,numr)

!!$  integer,save :: allochere=0
!!$  DATATYPE,save,allocatable :: intemptr(:,:), ttvector(:,:), ttvector2(:,:), &
!!$       ttvectorspin(:,:), ttvector2spin(:,:)
!!$  if (allochere.eq.0) then
!!$     allocate(intemptr(numr,spintotrank), ttvector(numconfig,numr), ttvector2(botwalk:topwalk,numr), &
!!$       ttvectorspin(spintotrank,numr), ttvector2spin(spinrank,numr))
!!$  endif
!!$  allochere=1

  if (sparseconfigflag.eq.0) then
     OFLWR "error, must use sparse for parblockconfigmult_transpose"; CFLST
  endif

  intemptr(:,:)=0d0;   intemptr(:,spinstart:spinend)=inavectortrspin(:,:)
  call mpiallgather(intemptr,spintotrank*numr,allspinranks*numr,maxspinrank*numr)
  ttvectorspin(:,:)=TRANSPOSE(intemptr(:,:))
  call configspin_transformfrom(numr,ttvectorspin,ttvector)
  if (dfrestrictflag.ne.0) then
     call dfrestrict(ttvector,numr)
  endif
  call sparseconfigmult_nompi(ttvector,ttvector2, yyy%cptr(0), yyy%sptr(0), 1,1,1,0,0d0)
  if (mshift.ne.0d0) then 
     do ir=1,numr
        ttvector2(:,ir)=ttvector2(:,ir)+ ttvector(botwalk:topwalk,ir)*configmvals(botwalk:topwalk)*mshift
     enddo
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict_par(ttvector2,numr)
  endif
  call configspin_transformto_mine(numr,ttvector2,ttvector2spin)
  
  outavectortrspin(:,:)=TRANSPOSE(ttvector2spin(:,:))

end subroutine parblockconfigmult_transpose_spin


subroutine blocklanczos0( lanblocknum, numout, lansize,maxlansize,order,maxiter,  outvectors,outvectorlda, outvalues,inprintflag,guessflag,lancheckmod,lanthresh,multsub,logpar,targetflag,etarget)
  use fileptrmod
  implicit none 
  logical,intent(in) :: logpar
  integer,intent(in) :: lansize,maxlansize,maxiter,lanblocknum,inprintflag,order,&
       lancheckmod,outvectorlda,numout,targetflag,guessflag
  DATATYPE, intent(in) :: etarget
  DATATYPE, intent(out) :: outvalues(numout)  
  DATATYPE, intent(inout) :: outvectors(outvectorlda,numout)
  real*8, intent(in) :: lanthresh
  integer :: printflag, iorder,k,flag,j,id,nn,i,&
       thislanblocknum, thisdim,ii,nfirst,nlast,myrank,nprocs
  real*8 :: error(numout), stopsum,rsum,nextran
  DATATYPE :: alpha(lanblocknum,lanblocknum),beta(lanblocknum,lanblocknum), csum, &
       thisdot ,nulldot , hdot,        nullvector1(1), nullvector2(1) , &
       lastvalue(numout), thisvalue(numout), valdot(numout),normsq(numout)
!! made these allocatable to fix lawrencium segfault 04-15
  DATATYPE, allocatable  ::       lanham(:,:,:,:),      laneigvects(:,:,:),&
       values(:),       templanham(:,:)
  DATATYPE, allocatable :: betas(:,:,:),betastr(:,:,:)
!! made these allocatable to fix carver segfault djh 03-30-15
  DATATYPE, allocatable :: &   
       initvectors(:,:),  invector(:), multvectors(:,:), &
       lanvects(:,:,:), tempvectors(:,:), &
       lanmultvects(:,:,:), tempvectors2(:,:)
  external :: multsub

  if (numout.lt.lanblocknum) then
     OFLWR "numout >= lanblocknum please",numout,lanblocknum; CFLST
  endif
  if (numout.gt.order*lanblocknum) then
     OFLWR "numout gt order*lanblocknum is impossible.",numout,order,lanblocknum; CFLST
  endif

  allocate(&
       initvectors(maxlansize,lanblocknum),  invector(maxlansize), multvectors(maxlansize,lanblocknum), &
       lanvects(maxlansize,lanblocknum,order), tempvectors(maxlansize,numout), &
       lanmultvects(maxlansize,lanblocknum,order), tempvectors2(maxlansize,numout))
  initvectors=0d0; invector=0d0; multvectors=0d0; lanvects=0d0; tempvectors=0d0; 
  lanmultvects=0d0; tempvectors2=0d0
  allocate( &
       lanham(lanblocknum,order,lanblocknum,order),&
       laneigvects(lanblocknum,order,order*lanblocknum), &
       values(order*lanblocknum), &
       templanham(lanblocknum*order,lanblocknum*order))
  lanham=0d0; laneigvects=0d0; values=0d0; templanham=0d0

  printflag=inprintflag

  lastvalue(:)=1d10;  thisvalue(:)=1d9;   alpha=0; beta=0;

  call rand_init(1731.d0) 

  if (guessflag==0) then

     initvectors=0.0d0

     do k=1,lanblocknum

!! want to be sure to have exactly the same calc regardless of nprocs.
!!    but doesn't seem to work
!!

!! attempt to make runs exactly the same, for parorbsplit=3 (orbparflag in sincproject)
!!   regardless of nprocs

        if (logpar) then   
           call getmyranknprocs(myrank,nprocs)
           nfirst=myrank-1;           nlast=nprocs-myrank
        else
           nfirst=0; nlast=0
        endif

        do ii=1,nfirst
           do i=1,lansize
              csum=nextran() + (0d0,1d0) * nextran()
           enddo
        enddo
        do i=1,lansize
           initvectors(i,k)=nextran() + (0d0,1d0) * nextran()
        enddo
        do ii=1,nlast
           do i=1,lansize
              csum=nextran() + (0d0,1d0) * nextran()
           enddo
        enddo
     enddo
  else
     initvectors(1:lansize,:)=outvectors(1:lansize,1:lanblocknum)
  endif

  initvectors(lansize+1:,:)=0d0  ! why not

  flag=0
!  do while (flag==0)
  do while (flag.ne.1)
     flag=0
     if (lansize.gt.0) then
        lanvects(1:lansize,:,1)=initvectors(1:lansize,:)
     endif
     lanham(:,:,:,:)=0.0d0
     do i=1,lanblocknum
        if (lansize.eq.0) then
           call nullgramschmidt_fast(i-1, logpar)
        else
           call myhgramschmidt_fast(lansize, i-1, maxlansize, lanvects(:,:,1),lanvects(:,i,1),logpar)
        endif
     enddo
     do j=1,lanblocknum
        if (lansize.eq.0) then
           call multsub(nullvector1(:),nullvector2(:))
           call multsub(nullvector1(:),nullvector2(:))
        else
           call multsub(lanvects(:,j,1),multvectors(:,j))
           call multsub(multvectors(:,j),tempvectors2(:,j))
        endif
     enddo
     lanmultvects(:,:,1)=multvectors(:,:)

     if (lansize.eq.0) then
        call nulldots(lanblocknum,lanblocknum,alpha,logpar)
     else
        call allhdots(lanvects(:,:,1),multvectors(:,:),lansize,maxlansize,lanblocknum,lanblocknum,alpha,logpar)
     endif

     lanham(:,1,:,1)= alpha(:,:)

     if (printflag.ne.0) then
        OFLWR "FIRST ALPHA ", alpha(1,1); CFL
     endif

     error(:)=1000d0

     do j=1,lanblocknum  !! not numout, don't have them yet.
        if (lansize.eq.0) then
           normsq(j)=nulldot(logpar)
           valdot(j)=nulldot(logpar)
        else
           normsq(j)=hdot(lanvects(:,j,1),lanvects(:,j,1),lansize,logpar)  !! yes should be normed, whatever
           valdot(j)=hdot(lanvects(:,j,1),multvectors(:,j),lansize,logpar)
        endif
        values(j)=valdot(j)/normsq(j)

        if (lansize.eq.0) then
           error(j)=abs(&
                valdot(j)**2 / normsq(j)**2 - &
                nulldot(logpar)/normsq(j))
        else
           error(j)=abs(&
                valdot(j)**2 / normsq(j)**2 - &
                hdot(lanvects(:,j,1),tempvectors2(:,j),lansize,logpar)/normsq(j))
        endif
     enddo
     if (printflag.ne.0) then
        OFL; write(mpifileptr,'(A10,100E8.1)') " FIRST ERRORS ", error(1:numout); CFL
     endif
     if (numout.eq.lanblocknum) then
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
              call myhgramschmidt_fast(lansize, (iorder-1)*lanblocknum+i-1, maxlansize, lanvects,lanvects(:,i,iorder),logpar)
           endif
        enddo

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
           call nulldots(lanblocknum,lanblocknum,alpha,logpar)
        else
           call allhdots(lanvects(:,:,iorder),multvectors(:,:),lansize,maxlansize,lanblocknum,lanblocknum,alpha,logpar)
        endif
        lanham(:,iorder,:,iorder)=alpha(:,:)

        allocate(betas(lanblocknum,iorder-1,lanblocknum), betastr(lanblocknum,lanblocknum,iorder-1))

        if (lansize.eq.0) then
           call nulldots(lanblocknum*(iorder-1),lanblocknum,betas(:,:,:),logpar)
           call nulldots(lanblocknum,lanblocknum*(iorder-1),betastr(:,:,:),logpar)
        else
           call allhdots(lanvects(:,:,:),lanmultvects(:,:,iorder),lansize,maxlansize,&
                lanblocknum*(iorder-1),lanblocknum,betas(:,:,:),logpar)
           call allhdots(lanvects(:,:,iorder),lanmultvects(:,:,:),lansize,maxlansize,&
                lanblocknum,lanblocknum*(iorder-1),betastr(:,:,:),logpar)
        endif

!!$        call allhdots(lanvects(:,:,iorder),lanmultvects(:,:,iorder-1),lansize,maxlansize,lanblocknum,lanblocknum,betastr(:,:,iorder-1),logpar)

        do i=1,iorder-1
           lanham(:,i,:,iorder)=betas(:,i,:)
        enddo
        do i=1,iorder-1                !!$  not needed on paper (arnoldi), just iorder-1
!!$        do i=iorder-1,iorder-1
           lanham(:,iorder,:,i)=betastr(:,:,i)
        enddo
        deallocate(betas,betastr)

        if (mod(iorder,lancheckmod).eq.0.or.flag.ne.0) then  ! flag.ne.0 for maxiter, iorder

           templanham(:,:)=0d0
           templanham(1:lanblocknum*iorder,1:lanblocknum*iorder)=&
                RESHAPE(lanham(:,1:iorder,:,1:iorder),(/lanblocknum*iorder,lanblocknum*iorder/))
           thisdim=min(maxiter,iorder*lanblocknum)
           call CONFIGEIG(templanham,thisdim,order*lanblocknum,laneigvects,values)
           if (targetflag.ne.0) then
              call mysort2(laneigvects,values,thisdim,order*lanblocknum,etarget)
           endif
           outvalues(:)=values(1:numout)
           if (printflag.ne.0) then
              OFL; write(mpifileptr,'(A10,1000F19.12)') " ENERIGES ", values(1:numout); CFL
           endif
           thisvalue(:)=values(1:numout)
           stopsum=0d0
           do nn=1,numout
              rsum=abs(thisvalue(nn)-lastvalue(nn))
              if (rsum.gt.stopsum) then
                 stopsum=rsum
              endif
           enddo
           lastvalue(:)=thisvalue(:)

! so lanthresh is for error of HPsi...  empirically make guess as to 
!  when it might be converged based on change in energy

!!$           if (stopsum.lt.lanthresh/10.or.flag.ne.0) then   ! flag=1 for maxiter,maxorder
           if (stopsum.lt.lanthresh/4.or.flag.ne.0) then   ! flag=1 for maxiter,maxorder
              if (printflag.ne.0) then
                 OFL; write(mpifileptr,'(A25,2E12.5,I10)')  "checking convergence.",stopsum,lanthresh/10,thisdim; CFL
              endif

              outvectors = 0.0d0
              do  j=1, numout
                 if (lansize.gt.0) then
                    do k=1, iorder
                       do id=1,lanblocknum
                          if ((k-1)*lanblocknum+id.le.maxiter) then
                             outvectors(1:lansize,j) = outvectors(1:lansize,j) + laneigvects(id,k,j) * lanvects(1:lansize,id,k)
                          endif
                       enddo
                    enddo
                    outvectors(:,j)=outvectors(:,j)/sqrt(thisdot(outvectors(:,j),outvectors(:,j),lansize,logpar))
                 else
                    outvectors(:,j)=outvectors(:,j)/sqrt(nulldot(logpar))
                 endif
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

              do j=1,numout
                 if (lansize.eq.0) then
                    normsq(j)=nulldot(logpar)
                    valdot(j)=nulldot(logpar)
                 else
                    normsq(j)=hdot(outvectors(:,j),outvectors(:,j),lansize,logpar)
                    valdot(j)=hdot(outvectors(:,j),tempvectors(:,j),lansize,logpar)
                 endif
                 values(j)=valdot(j)/normsq(j)

                 if (lansize.eq.0) then
                    error(j)=abs(&
                         valdot(j)**2 / normsq(j)**2 - &
                         nulldot(logpar)/normsq(j))
                 else
                    error(j)=abs(&
                         valdot(j)**2 / normsq(j)**2 - &
                         hdot(outvectors(:,j),tempvectors2(:,j),lansize,logpar)/normsq(j))
                 endif

!!$                 tempvectors(1:lansize,j) = tempvectors(1:lansize,j) - values(j) * outvectors(1:lansize,j)
!!$                 error(j)=sqrt(abs(hdot(tempvectors(:,j),tempvectors(:,j),lansize,logpar)))/ &
!!$                      sqrt(abs(hdot(outvectors(:,j),outvectors(:,j),lansize,logpar)))

              enddo
              if (printflag.ne.0) then
                 OFL; write(mpifileptr,'(A10,100E8.1)') " ERRORS ", error(1:numout); CFL
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
                    OFLWR "MAX DIM REACHED, NO CONVERGENCE -- THRESHOLDS TOO HIGH? BUG?",stopsum,lanthresh; CFLST
                 endif
              else
!!$                 if (printflag.ne.0) then
!!$                    OFLWR "Converged", stopsum, lanthresh; CFL
!!$                 endif
                 flag=1
              endif
           endif
        endif
     enddo

     if (flag==1) then
        if (printflag.ne.0) then
           OFLWR "Converged. ",stopsum,lanthresh  !!, "HDOTS:"
!!$           do i=1,numout
!!$              write(mpifileptr,'(10000F8.3)') (hdot(outvectors(:,i),outvectors(:,j),lansize,logpar),j=1,numout)
!!$           enddo
           CFL
        endif
     else
        flag=0
        OFLWR "  Not converged. restarting.", stopsum,lanthresh; CFL
     endif
     initvectors(:,:)=0d0
     if (lansize.gt.0) then
        initvectors(1:lansize,1:lanblocknum)=outvectors(1:lansize,1:lanblocknum)
     endif
  enddo

  deallocate(  initvectors,invector,multvectors, lanvects, tempvectors, lanmultvects,tempvectors2)

  deallocate( lanham,       laneigvects,       values,       templanham)

end subroutine blocklanczos0





subroutine mysort2(inout, values,n,lda,etarget)
  implicit none
  integer,intent(in) :: lda, n
  DATATYPE, intent(in) :: etarget
  DATATYPE :: inout(lda,n), out(lda,n), values(lda), newvals(lda)
  real*8 :: lowval 
  integer :: taken(lda), order(lda),i,j,whichlowest, flag

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
         write(*,*) "WHICHLOWEST ERROR, J=",j," WHICHLOWEST=", whichlowest;        call mpistop()
     endif
     if (taken(whichlowest).ne.0) then
        write(*,*) "TAKEN ERROR.";        call mpistop()
     endif
     taken(whichlowest)=1;     order(j)=whichlowest
     out(:,j)=inout(:,order(j));     newvals(j)=values(order(j))
  enddo

  values(1:n)=newvals(1:n)
  inout(1:n,1:n)=out(1:n,1:n)

end subroutine mysort2
