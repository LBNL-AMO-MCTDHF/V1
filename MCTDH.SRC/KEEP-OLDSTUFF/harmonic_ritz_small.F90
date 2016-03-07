
!! include this file in another to define what's below
!! (files blocklanczos_real.F90 and blocklanczos_complex.F90 provided for convenience)
!! or uncomment the following and compile this file as complex or real version
! 
!  uncomment
!#define DATATYPE complex*16
!#define REALFLAG 0
!     or
!#define DATATYPE real
!#define REALFLAG 1
!   ,  and
!#define OFLWR write(*,*)
!#define CFL
!#define mpifileptr *
!#define CFLST stop
!#define OFL


module hr_grammod
contains

!! allows real valued indefinite (possibly negative) inner product.
!!  norm of vector (1 or -1) output in mysign.
!!  vectors 1:m need to be orthonormal
!!  their norms are in prevsigns(1:m)

subroutine myhgramschmidt_many(n, howmany, m, previous, vector,inalldotsub,prevsigns,mysign,myfactors)
  implicit none
  
! n is the length of the vectors; m is how many to orthogonalize to

  integer :: n,m,mysign,prevsigns(m),howmany
  DATATYPE :: previous(n,howmany,m), vector(n,howmany),tempdots(m,howmany),myhdots(m,1),&
       tempprev(n,m,howmany)
  real*8 :: myfactors(howmany)
  integer :: i,j,k
  DATATYPE :: norm,temp(1,1)
  external :: inalldotsub

  do i=1,howmany
     tempprev(:,:,i)=previous(:,i,:)
  enddo
  do j=1,2 
     if (m.ne.0) then
        myhdots(:,:)=0d0
        do k=1,howmany
           call inalldotsub(tempprev(:,:,k),vector(:,k),n,m,1,tempdots(:,k))
           myhdots(:,1)=myhdots(:,1)+tempdots(:,k)*myfactors(k)
        enddo
     endif
     do i=1,m
       vector(:,:)=vector-previous(:,:,i)* myhdots(i,1) * prevsigns(i)
     enddo
     norm=0d0
     do k=1,howmany
        call inalldotsub(vector(:,k),vector(:,k),n,1,1,temp)
        norm=norm+temp(1,1)*myfactors(k)
     enddo
     if (abs(imag(norm)).gt.abs(real(norm,8))*1d-8) then
        OFLWR "AAUG IMAGNORM ", norm
     endif
     if (real(norm,8).ge.0d0) then
        mysign= 1
     else
        mysign= (-1)
     endif
     norm=sqrt(abs(norm))
     vector=vector/norm
     if (abs(norm).lt.1e-7) then
        OFLWR "Gram schmidt norm",norm,m; CFL
     endif
  enddo

end subroutine myhgramschmidt_many

end module hr_grammod




module hr_setmod

contains

subroutine fast_harmritz_setup( &
     lanblocknum, lansize,order,maxiter,initvectors,etarget,&
     inprintflag,&
     ingenmultsub,inh0multsub,indotsub,&
     lanham,lansimpleham,lansimpleovl,lansimplemultovl,lanvects)
  use hr_grammod
  implicit none 

  integer,intent(in) :: lansize,maxiter,lanblocknum,inprintflag,order
  external :: ingenmultsub,inh0multsub,indotsub
  integer, parameter :: manynum=2
  DATATYPE, intent(in) ::       initvectors(lansize,lanblocknum),etarget
  DATATYPE ::       lanham(lanblocknum,order,lanblocknum,order),lansimpleham(lanblocknum,order,lanblocknum,order),&
       lansimpleovl(lanblocknum,order,lanblocknum,order),lansimplemultovl(lanblocknum,order,lanblocknum,order),&
       lanvects(lansize,lanblocknum,order), &
       lansimplemultvects(lansize,lanblocknum,order), lanmanyvects(lansize,manynum,lanblocknum,order)
  real*8 :: factors(manynum)
  integer ::  iorder,i,thislanblocknum, printflag,lansigns(lanblocknum,order)

!! should be okay actually
  if (order*lanblocknum.gt.maxiter) then
     OFLWR "What in the world?  Order*lanblocknum is greater than lansize.", order,lanblocknum,lansize; CFLST
  endif

  printflag=inprintflag

  do iorder=1,order

!! DOES PARTIAL BLOCK IF AT THE END (leaving this)
     thislanblocknum=lanblocknum
     if (iorder*lanblocknum.ge.maxiter) then
        thislanblocknum=maxiter-(iorder-1)*lanblocknum
     endif
     if (iorder.eq.1) then
        lanvects(:,:,1)=initvectors(:,:)
     else
        call ingenmultsub(lanvects(:,:,iorder-1),lanvects(:,:,iorder),thislanblocknum,lansize,indotsub,inh0multsub)
     endif

     call inh0multsub(lanvects(:,:,iorder),lansimplemultvects(:,:,iorder),thislanblocknum,lansize)

     lansimplemultvects(:,:,iorder)=lansimplemultvects(:,:,iorder)-etarget*lanvects(:,:,iorder)
     lanmanyvects(:,1,:,:)=lansimplemultvects(:,:,:)
     lanmanyvects(:,2,:,:)=lanvects(:,:,:)
     factors(1)=1d0
     factors(2)=0d0
!!$     factors(2)=(-1)*abs(etarget)**2

     do i=1,thislanblocknum
        call myhgramschmidt_many(lansize, manynum, (iorder-1)*lanblocknum+i-1,  &
             lanmanyvects,lanmanyvects(:,:,i,iorder),&
             indotsub,lansigns,lansigns(i,iorder),factors)
     enddo

     lansimplemultvects(:,:,:)=lanmanyvects(:,1,:,:)
     lanvects(:,:,:)=lanmanyvects(:,2,:,:)

  enddo

  call indotsub(lansimplemultvects(:,:,:),lanvects(:,:,:),lansize,order*lanblocknum,order*lanblocknum,lanham)

  do iorder=1,order
     do i=1,lanblocknum
        lanham(:,:,i,iorder)=lanham(:,:,i,iorder)*lansigns(:,:)
     enddo
  enddo

  call indotsub(lanvects(:,:,:),lanvects(:,:,:),lansize,order*lanblocknum,order*lanblocknum,lansimpleovl)
  lansimplemultvects(:,:,:)=lansimplemultvects(:,:,:)+etarget*lanvects(:,:,:)
  call indotsub(lanvects(:,:,:),lansimplemultvects(:,:,:),lansize,order*lanblocknum,order*lanblocknum,lansimpleham)
  call indotsub(lansimplemultvects(:,:,:),lansimplemultvects(:,:,:),lansize,order*lanblocknum,order*lanblocknum,lansimplemultovl)

end subroutine fast_harmritz_setup


subroutine harm_orthog_oneblock(invects, &
     lanblocknum, lansize,etarget,&
     inprintflag,&
     inh0multsub,indotsub,&
     lantrans)
  use hr_grammod
  use hr_eigenmod
  implicit none 

  integer,intent(in) :: lansize,lanblocknum,inprintflag
  external :: inh0multsub,indotsub
  integer, parameter :: manynum=2
  DATATYPE, intent(in) :: etarget
  DATATYPE ::    invects(lansize,lanblocknum), tempvects(lansize,lanblocknum), &
       lantrans(lanblocknum,lanblocknum),&
       lantransinv(lanblocknum,lanblocknum),&
       tempdots(lanblocknum,lanblocknum),&
       indots(lanblocknum,lanblocknum),&
       lansimplemultvects(lansize,lanblocknum)
  integer ::  i,printflag,j,k

  printflag=inprintflag

  call inh0multsub(invects(:,:),lansimplemultvects(:,:),lanblocknum,lansize)

  lansimplemultvects(:,:)=lansimplemultvects(:,:)-etarget*invects(:,:)
  call indotsub(lansimplemultvects(:,:),lansimplemultvects(:,:),lansize,lanblocknum,lanblocknum,indots)
  lantransinv(:,:)=indots(:,:)

!! BUGGY ON EDISON WITH GNU COMPILERS (OR DJH MAKEFILE OR MODULES ERROR)
!!   zgesdd, zgesvd, zgeev

#if REALFLAG == 0

!!$ print *, "eeeeg"
!!$  call invsqrtmat(lantransinv,lanblocknum)

    call hermsqrtmatsmooth(lantransinv,lantrans,lanblocknum,lanblocknum,0d-10)

!!!  print *, "ooof";stop

  call ZGEMM('N','C',lansize,lanblocknum,lanblocknum,(1d0,0d0),invects,lansize,lantransinv,lanblocknum,&
       (0d0,0d0),tempvects,lansize)

!!! print *, "laaapf"

!!$ print *, "ooguuu"
!!$    call symorthogmat(lantransinv,lanblocknum)
!!$ print *, "fsfds88";stop
!!$    call ZGEMM('N','N',lansize,lanblocknum,lanblocknum,(1d0,0d0),invects,lansize,lantransinv,lanblocknum,&
!!$        (0d0,0d0),tempvects,lansize)

#else
  call realsqrtmatsmooth(lantransinv,lantrans,lanblocknum,lanblocknum,0d-10)
  call DGEMM('N','T',lansize,lanblocknum,lanblocknum,1d0,invects,lansize,lantransinv,lanblocknum,&
       0d0,tempvects,lansize)
#endif

  invects(:,:)=tempvects(:,:)

!! check

  call inh0multsub(invects(:,:),lansimplemultvects(:,:),lanblocknum,lansize)

  lansimplemultvects(:,:)=lansimplemultvects(:,:)-etarget*invects(:,:)

  call indotsub(lansimplemultvects(:,:),lansimplemultvects(:,:),lansize,lanblocknum,lanblocknum,tempdots)
  do i=1,lanblocknum
      tempdots(i,i)=tempdots(i,i)-1d0
  enddo

  do i=1,lanblocknum
     do j=1,lanblocknum
        if (abs(tempdots(i,j)).gt.1d-10) then
           print *, "AACK DOTS:"
           do k=1,lanblocknum
              write(*,'(100E20.3)') tempdots(:,k)
           enddo
           stop
        endif
     enddo
  enddo

end subroutine harm_orthog_oneblock

end module hr_setmod


!! GMRES LINEAR SOLVE.

!! fixed order, no restart, not for solving things per se, just a finite-order approximation
!!  to the inverse, error unchecked.

module hr_basicsolvemod

contains

subroutine basicblocklansolve(&
     lanblocknum, lansize,order, maxiter, invectors, outvectors, &
     inprintflag,&
     ingenmultsub,inh0multsub,indotsub,     etarget)
  use hr_setmod
  use hr_eigenmod
  implicit none 

  integer,intent(in) :: lansize,maxiter,lanblocknum,inprintflag,order
  DATATYPE, intent(in) :: etarget,invectors(lansize,lanblocknum)
  DATATYPE, intent(out) :: outvectors(lansize,lanblocknum)
  DATATYPE ::   initvectors(lansize,lanblocknum), &
       lanham(lanblocknum*order,lanblocknum*order),lansimpleovl(lanblocknum*order,lanblocknum*order),&
       lansimplemultovl(lanblocknum,order,lanblocknum,order),&
       lansimpleham(lanblocknum,order,lanblocknum,order),&
       lanvects(lansize,lanblocknum,order),        lantrans(lanblocknum,lanblocknum), &
       solvevecs(lanblocknum,order,lanblocknum),&
       rhs(lanblocknum,order,lanblocknum)
  integer ::  thisdim,printflag
  external :: ingenmultsub,inh0multsub,indotsub

  if (order.le.0) then
     print *, "ORDERERROR", order; stop
  endif

!! should be okay actually
  if (order*lanblocknum.gt.maxiter) then
     OFLWR "What in the world?  Order*lanblocknum is greater than lansize.", order,lanblocknum,lansize; CFLST
  endif

  printflag=inprintflag
  thisdim=min(maxiter,order*lanblocknum)

  initvectors(:,:)=invectors(:,:)

!!$  print *, "HHARMorth"

  call harm_orthog_oneblock(initvectors, &
       lanblocknum, lansize,etarget,&
       inprintflag,&
       inh0multsub,indotsub,&
       lantrans)

!!$  print *, "okaay..."; stop

  call fast_harmritz_setup( &
       lanblocknum, lansize,order,maxiter,initvectors,etarget,&
       inprintflag,&
       ingenmultsub,inh0multsub,indotsub,&
       lanham,lansimpleham,lansimpleovl,lansimplemultovl,lanvects)

  rhs(:,:,:)=0d0
  rhs(:,1,:)=lantrans(:,:)

#if REALFLAG == 0
  call ZGEMM('N','N',thisdim,lanblocknum,thisdim,(1d0,0d0),lanham(:,:),order*lanblocknum,rhs,order*lanblocknum,&
       (0d0,0d0),solvevecs,order*lanblocknum)
  call ZGEMM('N','N',lansize,lanblocknum,thisdim,(1d0,0d0),lanvects(:,:,:),lansize,solvevecs,order*lanblocknum,&
       (0d0,0d0),outvectors,lansize)
#else
  call DGEMM('N','N',thisdim,lanblocknum,thisdim,1d0,lanham(:,:),order*lanblocknum,rhs,order*lanblocknum,&
       0d0,solvevecs,order*lanblocknum)
  call DGEMM('N','N',lansize,lanblocknum,thisdim,1d0,lanvects(:,:,:),lansize,solvevecs,order*lanblocknum,&
       0d0,outvectors,lansize)
#endif

end subroutine basicblocklansolve

end module hr_basicsolvemod
