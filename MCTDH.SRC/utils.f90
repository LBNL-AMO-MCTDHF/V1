
!!
!!  UTILITIES.
!!  

#include "Definitions.INC"

subroutine waitawhile()
  implicit none
  integer :: i,j
  character (len=10) :: mytext
  j=1
  do i=1,10000000
     j=mod(j*i,1777)
  enddo
  write(mytext,'(I10)') j
  call system("echo "//mytext//" >> /dev/null")
end subroutine waitawhile


subroutine checkiostat(iniostat,intext)
  use fileptrmod
  implicit none
  integer,intent(in) :: iniostat
  character*(*),intent(in) :: intext
  if (iniostat /=0 ) then
     print *, "MCTDHF ABORT: I/O error ", iniostat,intext
     OFLWR "MCTDHF ABORT: I/O error ", iniostat,intext; CFL
     call waitawhile()
     stop
     stop   !!   STOP.   !!
     stop
  endif
end subroutine checkiostat

!! v1.27 getlen now reports length of string not length of string plus one

function getlen(buffer)
  implicit none
  character buffer*(*)
  integer :: j, getlen, mylen
  mylen=LEN(buffer)
  j=1
  do while (j.le.mylen)
     if (buffer(j:j) .eq. " ") then
        getlen=j-1
        return
     else
        j=j+1
     endif
  enddo
  getlen=mylen
end function getlen

function getlen2(buffer)
  implicit none
  character buffer*(*)
  integer :: j, getlen2, nn
  nn=LEN(buffer)-4
  j=1
  do while ((j.lt.nn).and..not.(buffer(j:j+3) .eq. "    "))
     j=j+1
  enddo
  getlen2=j-1
end function getlen2


module utilmod
implicit none
contains

subroutine staticvector(vector,size)
  implicit none
  integer :: i,size
  DATATYPE :: vector(size)
  real*8 :: nextran
!! don't init... keep as is
  do i=1,size
     vector(i)=nextran() + (0d0,1d0) * nextran()
  enddo
end subroutine staticvector

subroutine checkorbsetrange(checknspf,flag)
  use parameters
  implicit none
  integer,intent(in) :: checknspf
  integer,intent(out) :: flag
  flag=0
  if (nspf.ne.checknspf) then
     flag=1
  endif
end subroutine checkorbsetrange


!! difflevel:
!! 0=circular boundary conditions  1=interval bc  2=interval bc, subtract average after

subroutine complexdiff(size,in,out,difflevel)  
  use fileptrmod
  implicit none
  integer, intent(in) :: size, difflevel
  complex*16, intent(in) :: in(size)
  complex*16, intent(out) :: out(size)
  integer :: i,jj

  if (difflevel==0) then
     OFLWR "error complexdiff level=0"; CFLST
  elseif (difflevel==1) then

!! guarantees F.T. at zero is zero, right?

     do i=1,size
        out(i)= &
             1d0/280d0 * in(cindex(i-4)) &
             - 4d0/105d0 * in(cindex(i-3)) &
             + 1d0/5d0 * in(cindex(i-2)) &
             - 4d0/5d0 * in(cindex(i-1)) &
             + 4d0/5d0 * in(cindex(i+1)) &
             - 1d0/5d0 * in(cindex(i+2)) &
             + 4d0/105d0 * in(cindex(i+3)) &
             - 1d0/280d0 * in(cindex(i+4))
     enddo

  else

     out(:)=0d0
     do i=2,size-1
        jj=min(min(i-1,size-i),4)
        select case(jj)
        case(1)
           out(i)= &
                - 1d0/2d0 * in(i-1) &
                + 1d0/2d0 * in(i+1)
        case(2)
           out(i)= &
                1d0/12d0 * in(i-2) &
                - 2d0/3d0 * in(i-1) &
                + 2d0/3d0 * in(i+1) &
                - 1d0/12d0 * in(i+2)
        case(3)
           out(i)= &
                - 1d0/60d0 * in(i-3) &
                + 3d0/20d0 * in(i-2) &
                - 3d0/4d0 * in(i-1) &
                + 3d0/4d0 * in(i+1) &
                - 3d0/20d0 * in(i+2) &
                + 1d0/60d0 * in(i+3)
        case(4)
           out(i)= &
                1d0/280d0 * in(i-4) &
                - 4d0/105d0 * in(i-3) &
                + 1d0/5d0 * in(i-2) &
                - 4d0/5d0 * in(i-1) &
                + 4d0/5d0 * in(i+1) &
                - 1d0/5d0 * in(i+2) &
                + 4d0/105d0 * in(i+3) &
                - 1d0/280d0 * in(i+4)
        end select
     end do

     if (difflevel.eq.2) then
        out(:) = out(:) - sum(out(:))/size
     endif

  endif  !! docirc

contains
  function cindex(inindex)
    integer :: cindex,inindex
    cindex=mod(2*size+inindex-1,size)+1
  end function cindex

end subroutine complexdiff


subroutine zfftf_wrap_diff(size,inout,diffdflag)
  implicit none
  integer, intent(in) :: size,diffdflag
  complex*16, intent(inout) :: inout(size)
  complex*16,allocatable :: work(:)
  integer :: i

  if (diffdflag.eq.0) then
     call zfftf_wrap(size,inout)
  else
     allocate(work(size)); work=0

     call complexdiff(size,inout,work,diffdflag)

     call zfftf_wrap(size,work)

     do i=1,size
        inout(i)=work(i)*facfunct(i-1,size-1,1)
     enddo

     deallocate(work)
  endif

contains
  function facfunct(myindex,numdata,diffdflag)
    use fileptrmod
    implicit none
    integer, intent(in) :: myindex,diffdflag,numdata
    complex*16 :: facfunct,ccsum
    real*8, parameter :: twopi = 6.28318530717958647688d0
    if (myindex.lt.0.or.myindex.gt.numdata) then
       OFLWR "FACFUNCT ERR", myindex,0,numdata; CFLST
    endif
    ccsum=1d0
    if (diffdflag.ne.0) then
       if (myindex.ne.0) then
          ccsum= 1d0 / ((0d0,1d0)*myindex) / twopi * (numdata+1)
       else
          ccsum=0d0
       endif
    endif
    facfunct=ccsum
  end function facfunct

end subroutine zfftf_wrap_diff


subroutine zfftf_wrap(size,inout)
  implicit none
  integer, intent(in) :: size
  complex*16, intent(inout) :: inout(size)
  complex*16,allocatable :: wsave(:)
  allocate(wsave(4*size+15))
  wsave(:)=0d0
  call zffti(size,wsave)
  call zfftf(size,inout,wsave)
  deallocate(wsave)
end subroutine zfftf_wrap


subroutine zfftb_wrap(size,inout)
  implicit none
  integer, intent(in) :: size
  complex*16, intent(inout) :: inout(size)
  complex*16,allocatable :: wsave(:)
  allocate(wsave(4*size+15))
  wsave(:)=0d0
  call zffti(size,wsave)
  call zfftb(size,inout,wsave)
  deallocate(wsave)
end subroutine zfftb_wrap


function floatfac(in)
  implicit none
  integer,intent(in) :: in
  integer ::  i
  real*8 :: floatfac, sum
  sum=1.d0
  do i=1,in
     sum=sum*i
  enddo
  floatfac=sum
end function floatfac


function myisnan(input)
  implicit none
  real*8,intent(in) :: input
  logical :: myisnan
  if ((input+1.0d0.eq.input)) then
     myisnan=.true.
  else
     myisnan=.false.
  endif
end function myisnan



subroutine gramschmidt(n, m, lda, previous, vector,parflag)
  implicit none
  ! n is the length of the vectors; m is how many to orthogonalize to
  integer,intent(in) :: n,m,lda
  logical,intent(in) :: parflag
  DATATYPE,intent(in) :: previous(lda,m)
  DATATYPE,intent(inout) :: vector(n)
  CNORMTYPE :: norm
  integer ::  i,j

  do j=1,2
     do i=1,m
        vector=vector-previous(1:n,i)* heredot(previous(1:n,i),vector,n) 
     enddo
     norm=sqrt(heredot(vector,vector,n))  !! ok conversion
     vector=vector/norm
  enddo

contains

  function heredot(bra,ket,size)
    use mpisubmod
    implicit none
    integer,intent(in) :: size
    DATATYPE,intent(in) :: ket(size),bra(size)
    DATATYPE :: heredot,csum,csum2
    integer :: ii
    csum2=0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,csum)
    csum=0
!$OMP DO SCHEDULE(DYNAMIC)
    do ii=1,size
       csum=csum+CONJUGATE(bra(ii))*ket(ii)
    enddo
!$OMP END DO
!$OMP CRITICAL
    csum2=csum2+csum
!$OMP END CRITICAL
!$OMP END PARALLEL
    if (parflag) then
       call mympireduceone(csum2)
    endif
    heredot=csum2
  end function heredot

end subroutine gramschmidt


subroutine nullgramschmidt(m,parflag)
  implicit none
  ! n is the length of the vectors; m is how many to orthogonalize to
  integer,intent(in) :: m
  logical,intent(in) :: parflag
  DATATYPE :: norm
  integer ::  i,j

  if (.not.parflag) then
     print *, "parflag error dude"; stop
  endif
  do j=1,2
     do i=1,m
        norm=mynulldot()
     enddo
     norm=sqrt(mynulldot())
  enddo

contains

  function mynulldot()
    use mpisubmod
    implicit none
    DATATYPE :: mynulldot,csum
    if (parflag) then
       csum=0d0
       call mympireduceone(csum)
       mynulldot=csum
    else
       print *, "heredot error"; stop
    endif
  end function mynulldot

end subroutine nullgramschmidt

! enter as l1 l2 m1 m2 l3    TIMES TWO!!   INTEGER ARGUMENTS FOR HALF SPIN

function doubleclebschsq (l2,l1,m2,m1,l3)
  implicit none
  integer,intent(in) :: l1,l2,m1,m2,l3
  integer :: ierr,l3min,l3max
  real*8 :: dl1, dl2, dm1,dm2, dl3min, dl3max,dl3,dm3
  real*8 ::  doubleclebschsq, thrcof(100), sum
  integer, save :: ndim=100

!! real valued variables are m_s, not 2x m_s as in the arguments to this subroutine

  dl1=l1/2.d0
  dl2=l2/2.d0
  dm1=m1/2.d0
  dm2=m2/2.d0

  dm3=dm1+dm2

  call DRC3JJ(dl1,dl2,dm1,dm2,dl3min, dl3max, thrcof, ndim, ierr)

  l3min = int(2*dl3min + 0.0000001d0)
  l3max = int(2*dl3max + 0.0000001d0)

!! l3 and l3min are x 2    

  if (.not.(ierr .eq. 0)) then
     sum=0.0d0
  elseif ( (l3-l3min) < 0 ) then
     sum = 0.0d0
  elseif ( (l3-l3max) > 0 ) then
     sum = 0.0d0
  else
     !! cg coef squared.  avoids sqrt(-1) what do you do in this sitch  
     !!    is sqrt branch related the condon shortley phase convention

!!     sum = (-1)**(l1+l2+l3) * thrcof(l3-l3min+1)   myclebsch l's not doubled

!XXXXXARGH modulus function stupidly defined    if (mod(l3-l3min,2).eq.1) then

     if (mod(l3+l3min,2).eq.1) then
        print *, "L3ERROR",l3,l3min; stop
     endif

!! we are doing  clebsch^h.c.  clebsch  real valued (we are doing dot product squared in projeflux)

     sum = thrcof((l3-l3min)/2+1)**2

  endif

     dl3 = l3

!! again hermitian squared
!!  sum = sum * (-1)**(l1+l2-m3)
!!  myclebsch = sqrt(2*dl3+1) * sum

  doubleclebschsq = (2*dl3+1) * sum

end function doubleclebschsq

  

#ifdef REALGO

subroutine assigncomplex(realmat,complexf)
  implicit none
  real*8,intent(out) :: realmat
  real*8,intent(in) :: complexf
  realmat=complexf
end subroutine assigncomplex

subroutine assigncomplexmat(realmat,complexf,m,n)
  implicit none
  integer,intent(in) :: n,m
  real*8,intent(out) :: realmat(m,n)
  real*8,intent(in) :: complexf(m,n)
  realmat(:,:)=complexf(:,:)
end subroutine assigncomplexmat

subroutine assigncomplexvec(realmat,complexf,m)
  implicit none
  integer,intent(in) :: m
  real*8,intent(out) :: realmat(m)
  real*8,intent(in) :: complexf(m)
  realmat(:)=complexf(:)
end subroutine assigncomplexvec

subroutine assignrealvec(complexf,realmat,m)
  implicit none
  integer,intent(in) :: m
  real*8,intent(in) :: realmat(m)
  real*8,intent(out) :: complexf(m)
  complexf(:)=realmat(:)
end subroutine assignrealvec

#else

subroutine assigncomplex(realmat,complexf)
  implicit none
  complex*16,intent(in) :: complexf
  real*8,intent(out) :: realmat(2,2)
  realmat(1,1)=real(complexf,8);  realmat(2,2)=real(complexf,8)
  realmat(2,1)=imag(complexf);  realmat(1,2)=(-1)*imag(complexf)
end subroutine assigncomplex

subroutine assigncomplexmat(realmat,complexf,m,n)
  implicit none
  integer,intent(in) :: n,m
  complex*16,intent(in) :: complexf(m,n)
  real*8,intent(out) :: realmat(2,m,2,n)
  realmat(1,:,1,:)=real(complexf(:,:),8);  realmat(2,:,2,:)=real(complexf(:,:),8)
  realmat(2,:,1,:)=imag(complexf(:,:));  realmat(1,:,2,:)=(-1)*imag(complexf(:,:))
end subroutine assigncomplexmat

subroutine assigncomplexvec(realmat,complexf,m)
  implicit none
  integer,intent(in) :: m
  complex*16,intent(in) :: complexf(m)
  real*8,intent(out) :: realmat(2,m)
  realmat(1,:)=real(complexf(:),8);  realmat(2,:)=imag(complexf(:))
end subroutine assigncomplexvec

subroutine assignrealvec(complexf,realmat,m)
  implicit none
  integer,intent(in) :: m
  complex*16,intent(out) :: complexf(m)
  real*8,intent(in) :: realmat(2,m)
  complexf(:)=realmat(1,:)+realmat(2,:)*(0d0,1d0)
end subroutine assignrealvec

#endif



subroutine checksym(mat,dim)
  use fileptrmod
  implicit none
  integer,intent(in) :: dim
  real*8,intent(in) :: mat(dim,dim)
  integer :: i,j
  real*8 :: sym,asym,tot
  integer, save :: icalled=0

  sym=0; asym=0;tot=0

  do i=1,dim
     do j=1,dim
        tot=tot+ (2*mat(i,j))**2
        sym=sym+ (mat(i,j)+mat(j,i))**2
        asym=asym+ (mat(i,j)-mat(j,i))**2
     enddo
  enddo
  if (tot.gt.1d-10) then
     if (asym.gt.sym*1d-10) then
        icalled=icalled+1
        OFLWR "SYM,ASYM,MAG ", sym/tot,asym/tot,tot; CFLST
     endif
  endif

end subroutine checksym

end module utilmod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      MODULE INVSUBMOD    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module invsubmod
contains

!!$  if relflag then
!!$  tol is the maximum ratio of smallest to biggest eigenvalue
!!$    i.e. tol is relative not absolute
!!$  it .not.relflag then tol is absolute

subroutine invmatsmooth(A,N,LDA,tol,relflag)  !! inverse of ANY matrix.
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=(A_in)**-1
  implicit none
  logical,intent(in) :: relflag
  integer,intent(in) :: N,LDA
  real*8,intent(in) :: tol
  real*8 :: mytol
  DATATYPE,intent(inout) :: A(LDA,N)
  integer :: lwork,i,j,k
  real*8,allocatable :: SV(:)
  DATATYPE,allocatable :: SAVEA(:,:),U(:,:),VT(:,:),work(:)
#ifndef REALGO
  DATATYPE,allocatable :: zwork(:)

  allocate(zwork(5*N))
#endif
  allocate(SV(N),SAVEA(LDA,N),U(N,N),VT(N,N),work(5*N))

  lwork=5*N

  SAVEA=A

#ifdef REALGO 
  call dgesvd('A','A',N,N,A,LDA,SV,U,N,VT,N,work,lwork,i)
#else
  call zgesvd('A','A',N,N,A,LDA,SV,U,N,VT,N,zwork,lwork,work,i)
#endif

  if (i.ne.0) then
     print *, "ERR SVD",i;
     do i=1,N
        print *, SAVEA(:,i)
     enddo
     stop
  endif

  if (relflag) then 
     mytol=tol*SV(1)
  else
     mytol=tol
  endif

!! apply the function
  do k=1,N
     if (SV(k).ne.0d0) then
        if(SV(k).lt.mytol) then    !! it is positive ,whatever

!!$ NAAAH    SV(k)= 1d0 / mytol * (3 * (SV(k) / mytol / ) - 2 *( SV(k) / mytol )**2)

!!$    SVD so SV is real, keeping old regularization for now v1.19

           SV(k)= 1d0 / mytol

        else
           SV(k) = 1d0 / SV(k)
        endif
     endif
  enddo

!! rebuild the matrix
  do j=1,N
    do i=1,N
      A(i,j) = 0d0
      do k=1,N
         A(i,j) = A(i,j) + ALLCON(VT(k,i)) * SV(k) * ALLCON(U(j,k))
       enddo
    enddo
  enddo

#ifndef REALGO
  deallocate(zwork)
#endif
  deallocate(SV,SAVEA,U,VT,work)

end subroutine invmatsmooth



!!$  if relflag then
!!$  tol is the maximum ratio of smallest to biggest eigenvalue
!!$    i.e. tol is relative not absolute
!!$  it .not.relflag then tol is absolute


subroutine realinvmatsmooth(A,N,tol,relflag)  !! inverse of ANY matrix.
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=(A_in)**-1
  implicit none
  logical,intent(in) :: relflag
  integer,intent(in) :: N
  real*8,intent(in) :: tol 
  real*8 :: mytol
  real*8,intent(inout) :: A(N,N)
  real*8,allocatable :: U(:,:),VT(:,:),work(:),SV(:)
  integer :: lwork,i,j,k

  allocate( U(N,N),VT(N,N),work(5*N),SV(N) )

  lwork=5*N

!! do the svd

  call dgesvd('A','A',N,N,A,N,SV,U,N,VT,N,work,lwork,i)
  if (i.ne.0) then
     print *, "ERR SVD";     stop
  endif

  if (relflag) then 
     mytol=tol*SV(1)
  else
     mytol=tol
  endif

!! apply the function
  do k=1,N
     if (SV(k).ne.0d0) then
        if(abs(SV(k)).lt.mytol) then

!!$ NAAAH         SV(k)= 1d0 / mytol * (3 * (SV(k) / mytol) - 2 *( SV(k) / mytol )**2)

!!$    SVD so SV is real, keeping old regularization for now v1.19

           SV(k) = 1d0 / mytol 
        else
           SV(k) = 1d0 / SV(k)
        endif
     endif
  enddo
!! rebuild the matrix
  do j=1,N
    do i=1,N
      A(i,j) = 0d0
      do k=1,N
         A(i,j) = A(i,j) + VT(k,i) * SV(k) * U(j,k)
       enddo
    enddo
  enddo

  deallocate( U,VT,work,SV)

end subroutine realinvmatsmooth

end module invsubmod


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      MODULE LNSUBMOD     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module lnsubmod
contains

subroutine neglnmat(A,N)
  implicit none
  integer,intent(in) :: N
  DATATYPE,intent(inout) :: A(N,N)
  call bothlnmat(A,N,-1)
end subroutine neglnmat


subroutine lnmat(A,N)
  implicit none
  integer,intent(in) :: N
  DATATYPE,intent(inout) :: A(N,N)
  call bothlnmat(A,N,+1)
end subroutine lnmat


function djhlog(incomplex)
  use bio_parameters
  implicit none
  complex*16 :: djhlog,incomplex
  select case(logbranch)
  case(0)
     djhlog = log(incomplex)
  case(1)
     djhlog = log((0d0,1d0)*incomplex)-log((0d0,1d0))
  case (2)
     djhlog = log((-0.8d0,-0.6d0)*incomplex)-log((-0.8d0,-0.6d0))
  case default
     djhlog = log((-0.8d0,0.6d0)*incomplex)-log((-0.8d0,0.6d0))
  end select
end function djhlog


subroutine bothlnmat(A,N, which) 
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=-ln(A_in)
  use tol_parameters
  use invsubmod
  implicit none
  integer,intent(in) :: N,which
  DATATYPE,intent(inout) :: A(N,N)
  integer :: lwork,i,j,k
  complex*16 :: eig(N)
  DATATYPE,allocatable :: VL(:,:),VR(:,:),work(:),rwork(:)

  allocate(VL(N,N),VR(N,N),work(8*N),rwork(8*N))

  lwork=8*N

#ifdef REALGO 
  call dgeev('V','V',N,A,N,rwork(1),rwork(1+N),VL,N,VR,N,work,lwork,i)
  do j=1,N
    eig(j)= (1d0,0d0)*rwork(j) + (0d0,1d0)*rwork(j+N)
  enddo
#else
  call zgeev('V','V',N,A,N,eig,VL,N,VR,N,work,lwork,rwork,i)
  VL(1:N,1:N)=ALLCON(VL(1:N,1:N))
#endif


!! MAY 2014
  VL(:,:)=TRANSPOSE(VR(:,:))
  call invmatsmooth(VL,N,N,invtol,.true.)

!! apply the function
  do k=1,N
     if (eig(k).ne.0d0) then
        if (real(djhlog(eig(k)**which),8).gt.log(1d0/lntol)) then
           eig(k)= djhlog((eig(k)/abs(eig(k)))**which/lntol)
        elseif (real(djhlog(eig(k)**which),8).lt.log(invtol)) then
           eig(k)= djhlog(invtol*(eig(k)/abs(eig(k)))**which)
        else
           eig(k) = djhlog( eig(k)**which )
        endif
     else
        if (which.eq.1) then
           eig(k) = log(invtol)
        elseif (which.eq.-1) then
           eig(k) = log(1d0/lntol)
        else
           print *, "oogablah"; stop
        endif
     endif
  enddo

!! rebuild the matrix
  do j=1,N
    do i=1,N
      A(i,j) = 0d0
      do k=1,N
!!FOR real valued (mctdh) this is still ok.  imaginary cancels.
!!
        A(i,j) = A(i,j) + VR(i,k) * &   !! ok imp conv mctdh
             eig(k) * VL(j,k)           !! ok imp conv mctdh
      enddo
    enddo
  enddo

  deallocate(VL,VR,work,rwork)

end subroutine bothlnmat

end module lnsubmod

!!$ removing check because we are regularizing and the error message is confusing
!!$  
!!$  !! djh: check to make sure (degenerate case).  training wheels?
!!$  !!  now doing proper spectral expansion so can probably remove but... leave for now
!!$  
!!$    VL=(A)*which; time=1.d0;lwsp=6*N*N; ideg=6
!!$    wsp=0.d0
!!$    iexph=1
!!$    call MYGPADM( ideg, N, time, VL, N,wsp,lwsp,ipiv,iexph,nscale,iflag)
!!$    
!!$    tempmat=RESHAPE(wsp(iexph:iexph+N*N-1),(/N,N/))-saveA(:,:)
!!$  
!!$    sum=0;sum2=0;sum3=0
!!$    do i=1,N
!!$       do j=1,N
!!$          sum=sum+abs(tempmat(i,j)**2)
!!$          sum2=sum2+abs(A(i,j)**2)
!!$          sum3=sum3+abs(tempmat(i,j)**2)
!!$       enddo
!!$    enddo
!!$  
!!$  !!$ it is being regularized, this is to be expected
!!$    if (ierr.lt.10.and.((abs(sum/sum2).gt.1.d-7.and.abs(sum2).gt.1d-22).or.(iflag.ne.0))) then
!!$       print *, "NEGLN ERR", iflag,abs(sum/sum2),iexph,lwsp,sum,sum2,sum3
!!$       if (ierr.le.1) then
!!$          do i=1,min(10,N)
!!$             write(*,'(100F12.6)') abs(saveA(1:min(10,N),i))
!!$          enddo
!!$          print *, "XX"
!!$          tempmat=RESHAPE(wsp(iexph:iexph+N*N-1),(/N,N/))
!!$          do i=1,min(10,N)
!!$             write(*,'(100F12.6)') abs(tempmat(1:min(10,N),i))
!!$          enddo
!!$          print *
!!$       endif
!!$             
!!$       ierr=ierr+1
!!$       if (ierr.eq.10) then
!!$          print *, "FURTHER NEGLN ERRS SUPPRESSED"
!!$       endif
!!$    endif
!!$  
!!$    deallocate(VL,VR,work,rwork,tempmat,wsp,saveA,ipiv)
!!$
!!$  end subroutine bothlnmat
!!$
!!$  end module lnsubmod


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      MODULE EXPSUBMOD    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module expsubmod
contains

subroutine expmat(A,N)
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=exp(A_in)
  use tol_parameters
  use invsubmod
  implicit none
  integer,intent(in) :: N
  DATATYPE,intent(inout) :: A(N,N)
  integer :: lwork,i,k,j
  complex*16,allocatable :: eig(:), CVL(:,:),CVR(:,:),CA(:,:)
  DATATYPE,allocatable :: VL(:,:),VR(:,:),work(:),rwork(:)

  allocate( eig(N), CVL(N,N),CVR(N,N),CA(N,N), VL(N,N),VR(N,N),work(8*N),rwork(4*N) )
  lwork=8*N

#ifdef REALGO 
  call dgeev('V','V',N,A,N,rwork(1),rwork(1+N),VL,N,VR,N,work,lwork,i)
  do j=1,N
    eig(j)= (1d0,0d0)*rwork(j) + (0d0,1d0)*rwork(j+N)
  enddo
#else
  j=0
  call zgeev('V','V',N,A,N,eig,VL,N,VR,N,work,lwork,rwork,i)
#endif

  VL(:,:)=TRANSPOSE(VR(:,:))
  call invmatsmooth(VL,N,N,invtol,.true.)

!! apply the function
  do k=1,N
     eig(k) = exp( eig(k) )
     CVR(:,k)=VR(:,k)*eig(k)
  enddo
  CVL(:,:)=VL(:,:)

!! rebuild the matrix

  call ZGEMM('N','T',N,N,N,(1d0,0d0),CVR,N,CVL,N,(0d0,0d0),CA,N)
  A(:,:)=CA(:,:)  !! OK IMP CONV

  deallocate( eig,cvl,cvr,ca,vl,vr,work,rwork )

end subroutine expmat

end module expsubmod


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      MODULE MATSUBMOD    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! orthog routines:
!!  FLAG 1:  if on input A(i,j) = <phi_i|phi_j>
!!              then on output varphi_j = sum_i A(i,j) |phi_i>  are orthonormal.
!!  FLAG 2:  if on input A(i,j) = <phi_i|phi_j>
!!              on output varphi_j = sum_i A(i,j) |phi_i> are biorthonormal to phi:
!!            : <varphi_i|phi_j>=delta_ij


module matsubmod
contains

!!$subroutine biorthogmat(A,N)
!!$  implicit none
!!$  integer,intent(in) :: N
!!$  DATATYPE,intent(inout) :: A(N,N)
!!$#ifdef CNORMFLAG
!!$  call symorthogmat(A,N,2)
!!$#else
!!$  call symorthogmat(A,N,2)
!!$#endif
!!$end subroutine biorthogmat


subroutine orthogmat(A,N)
  implicit none
  integer,intent(in) :: N
  DATATYPE,intent(inout) :: A(N,N)
#ifdef CNORMFLAG
  call cnormorthogmat(A,N)
#else
  call symorthogmat(A,N,1)
#endif
end subroutine orthogmat


!! A must be hermitian (positive definite too right?) !!

subroutine symorthogmat(A,N,flag)  
  implicit none
  integer,intent(in) :: N,flag
  DATATYPE,intent(inout) :: A(N,N)
  real*8, allocatable :: SV(:),rwork(:)
  DATATYPE, allocatable :: U(:,:),VT(:,:),work(:),zwork(:),Asave(:,:)
  integer :: lwork,i,j,k

  if ((flag.ne.1)) then   !!  flag 2 not used. not sure about it.
     !!                        and.(flag.ne.2)) then
     print *, "bad flag",flag
     stop
  endif

#ifdef CNORMFLAG
  if (flag/=2) then
     print *, "FAIL."; stop
  endif
#endif

  if (flag.ne.2) then !! for sqrt, using eigen; must be herm.
     do i=1,N;  do j=1,i
        if (abs((CONJUGATE(A(i,j)))-A(j,i)).gt.1.d-7) then
           print *, "SYM ERR ORTHOG", abs((A(i,j))-A(j,i));           stop
        endif
     enddo;  enddo
  endif

  lwork=5*N
  allocate(U(N,N),VT(N,N),SV(N),work(lwork),zwork(lwork),rwork(lwork),Asave(N,N))
  Asave(:,:)=A(:,:)

#ifdef REALGO 
  call dgesvd('A','A',N,N,A,N,SV,U,N,VT,N,work,lwork,i)
#else
  call zgesvd('A','A',N,N,A,N,SV,U,N,VT,N,zwork,lwork,work,i)
#endif

  if (i.ne.0) then
     print *, "ERR herm SVD",i,N
     do k=1,N
        writE(*,'(10000E9.2)') Asave(:,k)
     enddo
     stop
  endif

!! apply the function
  do k=1,N
    if(abs(SV(k)).lt.1d-10) then
      print *, "SING symorthogmat!", sv(k)
      SV(k) = 0d0
    else if (flag.eq.2) then
      SV(k) = 1d0 / SV(k)
    else
      SV(k) = 1d0 / sqrt(SV(k))     
    endif
  enddo

!! rebuild the matrix
  do j=1,N
     do i=1,N
        A(i,j) = 0d0
        do k=1,N

!! TRANSPOSE OF INVERSE FOR CMCTDH    um what?
#ifdef CNORMFLAG
           A(i,j) = A(i,j) + ALLCON(U(i,k)) * SV(k) * ALLCON(VT(k,j))  
#else
           A(i,j) = A(i,j) + U(i,k) * SV(k) * VT(k,j)                  !HC OF INVERSE
#endif
        enddo
     enddo
  enddo
  deallocate(U,VT,SV,work,zwork,rwork,Asave)

end subroutine symorthogmat


subroutine cnormorthogmat(A,N)
  use fileptrmod
  implicit none
  integer,intent(in) :: N
  DATATYPE,intent(inout) :: A(N,N)

!!  OFLWR "cnormorthogmat CHECKME"; CFLST

  call allpurposemat(A,N,1)
end subroutine cnormorthogmat



subroutine allpurposemat(A,N,flag)
  use tol_parameters
  use fileptrmod
  use invsubmod
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=sqrt^-1(A_in) if flag=1 or sqrt(A_in) if flag=2
  implicit none
  integer,intent(in) :: N,flag
  DATATYPE,intent(inout) :: A(N,N)
  complex*16,allocatable :: eig(:), CVL(:,:),CVR(:,:),CA(:,:)
  DATATYPE,allocatable :: VL(:,:),VR(:,:),work(:),rwork(:)
  integer :: lwork,i,k,j

  allocate(eig(N), CVL(N,N),CVR(N,N),CA(N,N), VL(N,N),VR(N,N),work(8*N),rwork(4*N))
  lwork=8*N

#ifdef REALGO 
  call dgeev('V','V',N,A,N,rwork(1),rwork(1+N),VL,N,VR,N,work,lwork,i)
  do j=1,N
    eig(j)= (1d0,0d0)*rwork(j) + (0d0,1d0)*rwork(j+N)
  enddo
#else
  j=0
  call zgeev('V','V',N,A,N,eig,VL,N,VR,N,work,lwork,rwork,i)
#endif

  VL(:,:)=TRANSPOSE(VR(:,:))
  call invmatsmooth(VL,N,N,invtol,.true.)

!! apply the function
  do k=1,N
     select case(flag)
     case(1)
        eig(k) = 1d0 / sqrt( eig(k) )
     case(2)
        eig(k) = sqrt( eig(k) )
     case default
        OFLWR "OOGA ALLPURPOSEMAT ",flag; CFLST
     end select
     CVR(:,k)=VR(:,k)*eig(k)
  enddo
  CVL(:,:)=VL(:,:)

!! rebuild the matrix

  call ZGEMM('N','T',N,N,N,(1d0,0d0),CVR,N,CVL,N,(0d0,0d0),CA,N)
  A(:,:)=CA(:,:)  !! OK IMP CONV

  deallocate(eig, CVL,CVR,CA, VL,VR,work,rwork)

end subroutine allpurposemat

end module matsubmod


