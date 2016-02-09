
!!
!!  UTILITIES.
!!  

#include "Definitions.INC"

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


subroutine zfftf_wrap_diff(size,inout,diffdflag)
  implicit none
  integer, intent(in) :: size,diffdflag
  complex*16, intent(inout) :: inout(size)
  complex*16,allocatable :: work(:)
  complex*16 :: facfunct
  integer :: i,jj


  if (diffdflag.eq.0) then
     call zfftf_wrap(size,inout)
  else

     allocate(work(size)); work=0

#define DOCIRC
#ifdef DOCIRC
!! guarantees F.T. at zero is zero, right?

     do i=1,size
        work(i)= &
             1d0/280d0 * inout(cindex(i-4)) &
             - 4d0/105d0 * inout(cindex(i-3)) &
             + 1d0/5d0 * inout(cindex(i-2)) &
             - 4d0/5d0 * inout(cindex(i-1)) &
             + 4d0/5d0 * inout(cindex(i+1)) &
             - 1d0/5d0 * inout(cindex(i+2)) &
             + 4d0/105d0 * inout(cindex(i+3)) &
             - 1d0/280d0 * inout(cindex(i+4))
     enddo

#else

     work(:)=0d0
     do i=2,size-1
        jj=min(min(i-1,size-i),4)
        select case(jj)
        case(1)
           work(i)= &
                - 1d0/2d0 * inout(i-1) &
                + 1d0/2d0 * inout(i+1)
        case(2)
           work(i)= &
                1d0/12d0 * inout(i-2) &
                - 2d0/3d0 * inout(i-1) &
                + 2d0/3d0 * inout(i+1) &
                - 1d0/12d0 * inout(i+2)
        case(3)
           work(i)= &
                - 1d0/60d0 * inout(i-3) &
                + 3d0/20d0 * inout(i-2) &
                - 3d0/4d0 * inout(i-1) &
                + 3d0/4d0 * inout(i+1) &
                - 3d0/20d0 * inout(i+2) &
                + 1d0/60d0 * inout(i+3)
        case(4)
           work(i)= &
                1d0/280d0 * inout(i-4) &
                - 4d0/105d0 * inout(i-3) &
                + 1d0/5d0 * inout(i-2) &
                - 4d0/5d0 * inout(i-1) &
                + 4d0/5d0 * inout(i+1) &
                - 1d0/5d0 * inout(i+2) &
                + 4d0/105d0 * inout(i+3) &
                - 1d0/280d0 * inout(i+4)
        end select
     end do

#endif
     
     call zfftf_wrap(size,work)

     do i=1,size
        inout(i)=work(i)*facfunct(i-1,size-1,diffdflag)
     enddo

     deallocate(work)

  endif

contains
  function cindex(inindex)
    integer :: cindex,inindex
    cindex=mod(2*size+inindex-1,size)+1
  end function cindex

end subroutine zfftf_wrap_diff


subroutine zfftf_wrap(size,inout)
  implicit none
  integer, intent(in) :: size
  complex*16, intent(inout) :: inout(size)
  complex*16,allocatable :: wsave(:)
  allocate(wsave(4*size+15))
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
  call zffti(size,wsave)
  call zfftb(size,inout,wsave)
  deallocate(wsave)
end subroutine zfftb_wrap


subroutine checkiostat(iniostat)
  use fileptrmod
  implicit none
  integer :: iniostat
  if (iniostat /=0 ) then
     OFLWR "IOSTAT = ", iniostat; CFLST
  endif
end subroutine checkiostat

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


function hermdot(one,two,n)
  implicit none
  integer,intent(in) :: n
  DATATYPE,intent(in) :: one(n), two(n)
  DATATYPE :: hermdot, sum
  integer :: i
  sum=0.d0
  do i=1,n
     sum =   sum + ALLCON(one(i)) *  two(i) 
  enddo
  hermdot=sum
end function

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

!! USE THIS FOR C-NORM DOT
function cdot(one,two,n)
  implicit none
  integer,intent(in) :: n
  DATATYPE,intent(in) :: one(n), two(n)
  DATATYPE :: cdot, sum
  integer :: i
  sum=0.d0
  do i=1,n
     sum = sum + one(i) * two(i) 
  enddo
  cdot=sum
end function cdot

!! USE THIS FOR CNORM DOT OF DATAECS TYPE
function ecsdot(one,two,n)
  implicit none
  integer,intent(in) :: n
  DATAECS,intent(in) :: one(n), two(n)
  DATAECS :: ecsdot, sum
  integer :: i
  sum=0.d0
  do i=1,n
     sum = sum + one(i) * two(i) 
  enddo
  ecsdot=sum
end function ecsdot

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

subroutine dummysub()
end subroutine dummysub


!! CHOSEN DOT PRODUCT!! USE FOR <A|B> if A and B are orbs or A-vectors

function pardot(one,two,n)
  implicit none
  integer,intent(in) :: n
  DATATYPE,intent(in) :: one(n),two(n)
  DATATYPE :: csum,dot,pardot
  csum=dot(one,two,n)
  call mympireduceone(csum)
  pardot=csum
end function


function realpardot(one,two,n)
  implicit none
  integer,intent(in) :: n
  real*8,intent(in) :: one(n),two(n)
  real*8 :: sum,realdot,realpardot
  sum=realdot(one,two,n)
  call mympirealreduceone(sum)
  realpardot=sum
end function

subroutine realpardotsub(one,two,n,out)
  implicit none
  integer,intent(in) :: n
  real*8,intent(in) :: one(n),two(n)
  real*8,intent(out) :: out
  real*8 :: sum,realdot
  sum=realdot(one,two,n)
  call mympirealreduceone(sum)
  out=sum
end subroutine realpardotsub

subroutine pardotsub(one,two,n,out)
  implicit none
  integer,intent(in) :: n
  DATATYPE,intent(in) :: one(n),two(n)
  DATATYPE,intent(out) :: out
  DATATYPE :: sum,dot
  sum=dot(one,two,n)
  call mympireduceone(sum)
  out=sum
end subroutine pardotsub

!! CALLED INSIDE OMP LOOPS
recursive function dot(one,two,n)
  implicit none
  integer,intent(in) :: n
  DATATYPE,intent(in) :: one(n), two(n)
  DATATYPE :: dot, sum
  integer :: i
  sum=0.d0
  do i=1,n
     sum = sum + CONJUGATE(one(i)) * two(i) 
  enddo
  dot=sum
end function dot

function realdot(one,two,n)
  implicit none
  integer,intent(in) :: n
  real*8,intent(in) :: one(n), two(n)
  real*8 :: realdot, sum
  integer :: i
  sum=0.d0
  do i=1,n
     sum = sum + one(i) * two(i) 
  enddo
  realdot=sum
end function realdot


subroutine realdotsub(one,two,n,out)
  implicit none
  integer,intent(in) :: n
  real*8,intent(in) :: one(n), two(n)
  real*8 :: out, sum
  integer :: i
  sum=0.d0
  do i=1,n
     sum = sum + one(i) * two(i) 
  enddo
  out=sum
end subroutine realdotsub


subroutine realgramschmidt(n, m, lda, previous, vector)
  use fileptrmod
  implicit none
  ! n is the length of the vectors; m is how many to orthogonalize to
  integer,intent(in) :: n,m,lda
  real*8,intent(in) :: previous(lda,m)
  real*8,intent(inout) :: vector(n)
  real*8 :: norm, realdot
  integer :: i

  do i=1,m
     vector=vector-previous(1:n,i)* realdot(previous(1:n,i),vector,n) 
  enddo
  norm=real(sqrt(realdot(vector,vector,n)),8)  !! ok for imp conv (c/p/ch)
  if (abs(norm).lt.1.d-8) then
     call openfile();    write(mpifileptr,*) "Warning, small norm in realgramschmidt: ", norm;     call closefile()
  endif
  vector=vector/norm
end subroutine realgramschmidt

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
     norm=sqrt(heredot(vector,vector,n))  !! ok for impconv (chmctdh,pmctdh)
     vector=vector/norm
!!     if (nancheck.and.((.not.(abs(norm).gt.1e-6)).or.(.not.(abs(norm).lt.1e+10)))) then
!!        call openfile();    write(mpifileptr, *) "TE MP NOSTOP  Gram schmidt norm",norm, heredot(vector,vector,n), m; !!    call closefile()
!!     endif
  enddo
contains

  function heredot(bra,ket,size)
    implicit none
    integer,intent(in) :: size
    DATATYPE,intent(in) :: ket(size),bra(size)
    DATATYPE :: heredot,dot,pardot
    if (parflag) then
       heredot=pardot(bra,ket,size)
    else
       heredot=dot(bra,ket,size)
    endif
  end function heredot
end subroutine gramschmidt


subroutine ecsgramschmidt(n, m, lda, previous, vector, kdoneflag)
  implicit none
  ! n is the length of the vectors; m is how many to orthogonalize to
  integer,intent(in) :: n,m,lda
  integer,intent(out) :: kdoneflag
  DATAECS,intent(in) :: previous(lda,m)
  DATAECS,intent(inout) :: vector(n)
  DATAECS :: norm, ecsdot
  integer :: i

  do i=1,m
     vector=vector-previous(1:n,i)* ecsdot(previous(1:n,i),vector,n) 
  enddo
  norm=sqrt(ecsdot(vector,vector,n));  vector=vector/norm
  if ((abs(norm).lt.1e-7)) then
     kdoneflag=1
  endif
end subroutine ecsgramschmidt



!!$function hermnormsq(one,n)
!!$  implicit none
!!$  integer :: n,i
!!$  DATATYPE :: one(n)
!!$  real*8 :: hermnormsq, sum
!!$  sum=0.d0
!!$  do i=1,n
!!$     sum = sum + abs(one(i)**2)
!!$  enddo
!!$  hermnormsq=sum
!!$end function hermnormsq


!! Begin KVL routines

subroutine neglnmat(A,N,lntol)
  implicit none
  integer,intent(in) :: N
  real*8,intent(in) :: lntol
  DATATYPE,intent(inout) :: A(N,N)
  call bothlnmat(A,N,-1,lntol)
end subroutine neglnmat


subroutine lnmat(A,N,lntol)
  implicit none
  integer,intent(in) :: N
  real*8,intent(in) :: lntol
  DATATYPE,intent(inout) :: A(N,N)
  call bothlnmat(A,N,+1,lntol)
end subroutine lnmat

function djhlog(incomplex)
  implicit none
  complex*16 :: djhlog,incomplex
  djhlog = log((0d0,1d0)*incomplex)-log((0d0,1d0))
end function djhlog


subroutine bothlnmat(A,N, which ,lntol) 
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=-ln(A_in)
  implicit none
  integer,intent(in) :: N,which
  real*8,intent(in) :: lntol
  DATATYPE,intent(inout) :: A(N,N)
  integer, save :: ierr=0
  real*8 :: time
  integer :: lwork,i,j,k,nscale,iflag
  complex*16 :: eig(N),djhlog 
  DATATYPE :: sum,sum2,sum3
  integer :: lwsp,ideg,iexph
  DATATYPE,allocatable :: VL(:,:),VR(:,:),work(:),rwork(:), tempmat(:,:),wsp(:),saveA(:,:)
  integer,allocatable :: ipiv(:)

  allocate(VL(N,N),VR(N,N),work(8*N),rwork(8*N), tempmat(N,N),wsp(6*N*N),saveA(N,N),ipiv(2*N*N))

  lwork=8*N

  saveA(:,:)=A(:,:)

!! identity check  (Why?  forget if ever was needed for any reason...
!  rsum=0d0
!  do i=1,N
!    do j=1,N
!      if (i==j) then
!        rsum=rsum+abs(A(i,i)-1d0)
!      else
!        rsum=rsum+abs(A(i,j))
!      endif
!    enddo
!  enddo
!  if(rsum.lt.1d-8) then
!    A=0d0
!!!PREV 043012    A=0d0
!!    A=0d0
!!    do i=1,N
!!       A(i,i)=-log(1d-8)
!!    enddo
!    return
!  endif


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
  call invmatsmooth(VL,N,N,0d0)


!! apply the function
  do k=1,N
     if (eig(k).ne.0d0) then
        if (abs(eig(k)).lt.lntol) then
           eig(k)= which * djhlog(lntol*eig(k)/abs(eig(k)))
        else
           eig(k) = which * djhlog( eig(k) )
        endif
     else
        print *,  "BAD! ZERO EIG LN.  FIXME. TEMP CONTINUE"
        eig(k) = which * djhlog((0d0,1d0)*lntol)
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

!! djh: check to make sure (degenerate case).  training wheels?
!!  now doing proper spectral expansion so can probably remove but... leave for now

  VL=(A)*which; time=1.d0;lwsp=6*N*N; ideg=6
  wsp=0.d0
  iexph=1
  call MYGPADM( ideg, N, time, VL, N,wsp,lwsp,ipiv,iexph,nscale,iflag)
  
  tempmat=RESHAPE(wsp(iexph:iexph+N*N-1),(/N,N/))-saveA(:,:)

  sum=0;sum2=0;sum3=0
  do i=1,N
     do j=1,N
        sum=sum+abs(tempmat(i,j)**2);        sum2=sum2+abs(A(i,j)**2); sum3=sum3+abs(tempmat(i,j)**2)
     enddo
  enddo
  if (ierr.lt.10.and.((abs(sum/sum2).gt.1.d-7.and.abs(sum2).gt.1d-22).or.(iflag.ne.0))) then
     print *, "NEGLN ERR", iflag,abs(sum/sum2),iexph,lwsp,sum,sum2,sum3
     if (ierr.le.1) then
        do i=1,min(10,N)
           write(*,'(100F12.6)') abs(saveA(1:min(10,N),i))
        enddo
        print *, "XX"
        tempmat=RESHAPE(wsp(iexph:iexph+N*N-1),(/N,N/))
        do i=1,min(10,N)
           write(*,'(100F12.6)') abs(tempmat(1:min(10,N),i))
        enddo
        print *
     endif
           
     ierr=ierr+1
     if (ierr.eq.10) then
        print *, "FURTHER NEGLN ERRS SUPPRESSED"
     endif
  endif

  deallocate(VL,VR,work,rwork,tempmat,wsp,saveA,ipiv)

end subroutine bothlnmat





subroutine expmat(A,N)
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=exp(A_in)
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

  call invmatsmooth(VL,N,N,0d0)


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


!!$  For invmatsmooth, tol is the maximum ratio of smallest to biggest eigenvalue
!!$   i.e. tol is relative not absolute

subroutine invmatsmooth(A,N,LDA,tol)  !! inverse of ANY matrix.
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=(A_in)**-1
  implicit none
  integer,intent(in) :: N,LDA
  real*8,intent(in) :: tol
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

!! do the svd


#ifdef REALGO 
  call dgesvd('A','A',N,N,A,LDA,SV,U,N,VT,N,work,lwork,i)
#else
  call zgesvd('A','A',N,N,A,LDA,SV,U,N,VT,N,zwork,lwork,work,i)
#endif

!  print *, SV, "SV";

  if (i.ne.0) then
     print *, "ERR SVD",i;
     do i=1,N
        print *, SAVEA(:,i)
     enddo
     stop
  endif
!! apply the function
  do k=1,N
     if (SV(k).ne.0d0) then
        if(SV(k).lt.tol*SV(1)) then    !! it is positive ,whatever

!!$ NAAAH              SV(k)= 1d0 / tol / SV(1) * (3 * (SV(k) / tol / SV(1)) - 2 *( SV(k) / tol / SV(1))**2)

!!$    SVD so SV is real, keeping old regularization for now v1.19

           SV(k)= 1d0 / tol / SV(1)

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



!!$  For realinvmatsmooth, tol is absolute not relative (refers to absolute mag of eigenvalue)

subroutine realinvmatsmooth(A,N,tol)  !! inverse of ANY matrix.
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=(A_in)**-1
  implicit none
  integer,intent(in) :: N
  real*8,intent(in) :: tol 
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
!! apply the function
  do k=1,N
     if (SV(k).ne.0d0) then
        if(abs(SV(k)).lt.tol) then

!!$ NAAAH           SV(k)= 1d0 / tol * (3 * (SV(k) / tol) - 2 *( SV(k) / tol )**2)

!!$    SVD so SV is real, keeping old regularization for now v1.19

           SV(k) = 1d0 / tol 

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


!! orthog routines:
!!  FLAG 1:  if on input A(i,j) = <phi_i|phi_j>
!!              then on output varphi_j = sum_i A(i,j) |phi_i>  are orthonormal.
!!  FLAG 2:  if on input A(i,j) = <phi_i|phi_j>
!!              on output varphi_j = sum_i A(i,j) |phi_i> are biorthonormal to phi:
!!            : <varphi_i|phi_j>=delta_ij

subroutine biorthogmat(A,N)
  implicit none
  integer,intent(in) :: N
  DATATYPE,intent(inout) :: A(N,N)
#ifdef CNORMFLAG
  call symorthogmat(A,N,2)
#else
  call symorthogmat(A,N,2)
#endif
end subroutine biorthogmat


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


subroutine symorthogmat(A,N,flag)  
  implicit none
  integer,intent(in) :: N,flag
  DATATYPE,intent(inout) :: A(N,N)
  real*8, allocatable :: SV(:),rwork(:)
  DATATYPE, allocatable :: U(:,:),VT(:,:),work(:),zwork(:),Asave(:,:)
  integer :: lwork,i,j,k

  if ((flag.ne.1).and.(flag.ne.2)) then
     print *, "bad flag",flag
     stop
  endif
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
#ifdef CNORMFLAG
  if (flag/=2) then
     print *, "FAIL."; stop
  endif
#endif
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
  use fileptrmod
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
  call invmatsmooth(VL,N,N,0d0)

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
!        if (icalled.lt.10) then
           OFLWR "SYM,ASYM,MAG ", sym/tot,asym/tot,tot; CFLST
!        endif
     endif
  endif

end subroutine checksym
