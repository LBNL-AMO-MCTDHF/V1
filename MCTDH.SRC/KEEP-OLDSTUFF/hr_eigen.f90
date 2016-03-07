
!! LAPACK WRAPPERS

!! REMOVED ALL LEFT RIGHT EIGENVECTOR ROUTINES, NOT NEEDED, AND USE POOR SPECTRAL EXPANSION

module hr_eigenmod

contains


subroutine mycsort_getorder(values,n,etarget,order)
  implicit none
  integer :: n,i,j,whichlowest, flag,taken(n), order(n)
  complex*16,intent(in) :: values(n),etarget
  real*8 :: lowval 

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
         write(*,*) "WHICHLOWEST ERROR, J=",j," WHICHLOWEST=", whichlowest;        stop
     endif
     if (taken(whichlowest).ne.0) then
        write(*,*) "TAKEN ERROR.";        stop
     endif
     taken(whichlowest)=1;     order(j)=whichlowest
  enddo

end subroutine mycsort_getorder




subroutine mycsort_vecs(inout,n,lda,order)
  implicit none
  complex*16 :: inout(lda,n), out(lda,n)
  integer :: lda, n,j
  integer :: order(n)

  do j=1,n
     out(:,j)=inout(:,order(j))
  enddo
  inout(:,:)=out(:,:)

end subroutine mycsort_vecs

subroutine mycsort_vals(values,n,order)
  implicit none
  complex*16 :: values(n),newvals(n)
  integer ::  n,j
  integer :: order(n)

  do j=1,n
     newvals(j)=values(order(j))
  enddo
  values(1:n)=newvals(1:n)

end subroutine mycsort_vals




subroutine myrsort_getorder(values,n,etarget,order)
  implicit none
  integer :: n,i,j,whichlowest, flag,taken(n), order(n)
  real*8,intent(in) :: values(n),etarget
  real*8 :: lowval 

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
         write(*,*) "WHICHLOWEST ERROR, J=",j," WHICHLOWEST=", whichlowest;        stop
     endif
     if (taken(whichlowest).ne.0) then
        write(*,*) "TAKEN ERROR.";        stop
     endif
     taken(whichlowest)=1;     order(j)=whichlowest
  enddo

end subroutine myrsort_getorder




subroutine myrsort_vecs(inout,n,lda,order)
  implicit none
  real*8 :: inout(lda,n), out(lda,n)
  integer :: lda, n,j
  integer :: order(n)

  do j=1,n
     out(:,j)=inout(:,order(j))
  enddo
  inout(:,:)=out(:,:)

end subroutine myrsort_vecs

subroutine myrsort_vals(values,n,order)
  implicit none
  real*8 :: values(n),newvals(n)
  integer ::  n,j
  integer :: order(n)

  do j=1,n
     newvals(j)=values(order(j))
  enddo
  values(1:n)=newvals(1:n)

end subroutine myrsort_vals



!! leave herm normed as given by zgeev:

subroutine get_eigen_c_plain(inmat, n, lda, outmat,values)
  implicit none
  
  complex*16 :: inmat(lda,n), outmat(lda,n), values(lda), vr(lda,n),vl(lda,n),tempmat(lda,n)
  integer :: n,lda, info, lwork
  character*1 :: jobvl="V", jobvr="V"
  complex*16 :: work(5*lda)
  real*8 :: rwork(5*lda)

  lwork=5*lda;  tempmat(1:n,1:n)=inmat(1:n,1:n)

  call zgeev( jobvl, jobvr, n, tempmat, lda, values, vl,lda,vr,lda,work, lwork, rwork, info)
  if (info /= 0) then
     write(*, *) "AAAUGH, info in get_eigen_c:",info;     stop
  endif
  outmat(:,:)=vr(:,:)

end subroutine get_eigen_c_plain




subroutine get_eigen_two(inmat, n, lda, outvects, outvals)
  implicit none
  real*8 :: inmat(lda,n), outvects(lda,n), outvals(n)
  integer :: n,lda, info, lwork
  character*1 :: job="V", uplo="U"
  real*8 :: work(5*lda)

  lwork=5*lda;  outvects = inmat
  call dsyev( job, uplo, n, outvects, lda, outvals, work, lwork, info)
  if (info /= 0) then
       write(*, *)"AAAUGH, in get_eigen_two dsyev info = ",info;     stop
  endif

end subroutine 


subroutine get_eigen_tri(indiag,insubdiag, n, lda, outvects, outvals)
  implicit none
  real*8 :: indiag(n), insubdiag(n-1), outvals(n), outvects(lda,n)
  integer :: n,lda, info
  character*1 :: job="V"
  real*8 :: work(3*n), subdiag(n-1) 
  outvals = indiag;  subdiag = insubdiag
  call dstev( job, n, outvals, subdiag, outvects, lda, work, info )
  if (info /= 0) then
    write(*, *)"AAAUGH, in get_eigen_tri dstev info = ",info;     stop
  endif
end subroutine 


subroutine get_eigen_herm(inmat, n, lda, outmat,values)
  implicit none
  complex*16 :: inmat(lda,n), outmat(lda,n)
  integer :: info, n,lda, lwork
  real*8 :: values(lda)
  character*1 :: jobz="V", uplo="U"
  complex*16 :: work(5*lda)
  real*8 :: rwork(5*lda)
  lwork=5*lda;  outmat=inmat
  call zheev( jobz, uplo, n, outmat, lda, values, work, lwork, rwork, info)
  if (info /= 0) then
    write(*, *) "AAAUGH, info in get_eigen_c:",info;     stop
  endif
end subroutine 

!! inverse of ANY matrix.
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=(A_in)**-1
subroutine cinvmatsmooth(A,N,lda,tol)
  implicit none
  integer :: N,lda
  real*8 :: tol
  complex*16 :: A(lda,N),B(lda,N)
  call gen_cmatsub(A,B,N,lda,tol,1)
end subroutine cinvmatsmooth


subroutine hermsqrtmatsmooth(A,B,N,lda,tol)
  implicit none
  integer :: N,lda
  real*8 :: tol
  complex*16 :: A(lda,N),B(lda,N)
  call gen_cmatsub(A,B,N,lda,tol,2)
end subroutine hermsqrtmatsmooth

!! inverse of any matrix for flag=1; sqrt of Hermitian matrix if flag=2.
subroutine gen_cmatsub(A,B,N,lda,tol,flag) 
  implicit none
  integer :: N,lwork,i,j,k,lda,flag
  real*8 :: SV(N),tol,SVsave(N)
  complex*16 :: A(lda,N),B(lda,N)
  real*8 :: work(20*N**2+100) 
  complex*16 :: U(N,N),VT(N,N),zwork(20*N**2+100)
!  integer :: iwork(8*N)
  integer :: iwork(100*N)
  lwork=10*N**2+100
!  lwork=20*N**2+100

!! do the svd... buggy on edison with gnu compilers I think

!  call zgesvd('A','A',N,N,A,lda,SV,U,N,VT,N,zwork,lwork,work,i)

!print *, "GOSVD",N,lda
!print *, A(:,:)
!print *, "AA";stop

  call zgesdd('A',N,N,A,lda,SV,U,N,VT,N,zwork,lwork,work,iwork,i)

!print *, "OKSVD";stop

  SVsave(:)=SV(:)

  if (i.ne.0) then
     print *, "ERR SVD cmatsub TEMP CONTINUE",i,n;    !! stop
  endif
!! apply the function
  do k=1,N
     if (SV(k).ne.0d0) then
        if(SV(k).lt.tol*SV(1)) then    !! it is positive ,whatever
           SV(k)= 1d0 / tol / SV(1)
        else
           SV(k) = 1d0 / SV(k)
        endif
     endif
  enddo
  select case(flag)
  case(1)
  case(2)
     SV(:)=SQRT(SV(:))
     SVsave(:)=SQRT(SVsave(:))
  case default
     print *, "AAAUGH!!! flagcheck",flag; stop
  end select

!! rebuild the matrix

  A(:,:) = 0d0
  do j=1,N
    do i=1,N
      do k=1,N
         A(i,j) = A(i,j) + CONJG(VT(k,i)) * SV(k) * CONJG(U(j,k))
       enddo
    enddo
  enddo
  if (flag.eq.2) then
     B(:,:) = 0d0
     do j=1,N
        do i=1,N
           do k=1,N
              B(i,j) = B(i,j) + CONJG(VT(k,i)) * SVsave(k) * CONJG(U(j,k))
           enddo
        enddo
     enddo
  endif
end subroutine gen_cmatsub


!! inverse of ANY matrix.
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=(A_in)**-1
subroutine realinvmatsmooth(A,N,lda,tol)
  implicit none
  integer :: N,lda
  real*8 :: A(lda,N),tol,B(lda,N)
  call gen_rmatsub(A,B,N,lda,tol,1)
end subroutine realinvmatsmooth

!! on input : matrix A
!! on output : A-> A^(-1/2); B->A^(1/2)

subroutine realsqrtmatsmooth(A,B,N,lda,tol)
  implicit none
  integer :: N,lda
  real*8 :: A(lda,N),tol,B(lda,N)
  call gen_rmatsub(A,B,N,lda,tol,2)
end subroutine realsqrtmatsmooth


!! inverse of any matrix for flag=1; sqrt of Hermitian matrix if flag=2.
subroutine gen_rmatsub(A,B,N,lda,tol,flag)  
  implicit none
  integer :: N,lwork,i,j,k,lda,flag
  real*8 :: SV(N),tol,SVsave(N)
  real*8 :: A(lda,N),B(lda,N)
  real*8 :: U(N,N),VT(N,N),work(5*N)
  lwork=5*N

!! do the svd

  call dgesvd('A','A',N,N,A,lda,SV,U,N,VT,N,work,lwork,i)
  SVsave(:)=SV(:)

  if (i.ne.0) then
     print *, "ERR SVD rmatsub",i;     stop
  endif
!! apply the function
  do k=1,N
     if (SV(k).ne.0d0) then
        if(SV(k).lt.tol*SV(1)) then    !! it is positive ,whatever
           SV(k)= 1d0 / tol / SV(1)
        else
           SV(k) = 1d0 / SV(k)
        endif
     endif
  enddo

  select case(flag)
  case(1)
  case(2)
     SV(:)=SQRT(SV(:))
     SVsave(:)=SQRT(SVsave(:))
  case default
     print *, "AAAUGH!!! flagcheck",flag; stop
  end select

!! rebuild the matrix
  A(:,:) = 0d0
  do j=1,N
    do i=1,N
      do k=1,N
         A(i,j) = A(i,j) + VT(k,i) * SV(k) * U(j,k)
       enddo
    enddo
  enddo
  if (flag.eq.2) then
     B(:,:) = 0d0
     do j=1,N
        do i=1,N
           do k=1,N
              B(i,j) = B(i,j) + VT(k,i) * SVsave(k) * U(j,k)
           enddo
        enddo
     enddo
  endif

end subroutine gen_rmatsub


!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=exp(A_in)
!!
subroutine cexpmat(A,N,lda)
  implicit none
  integer :: N,lwork,i,k,lda
  complex*16 :: eig(N)  ,CVL(N,N),CVR(N,N),CA(lda,N)
  complex*16 :: A(lda,N)
  complex*16 :: VL(N,N),VR(N,N),work(8*N),rwork(4*N)

  lwork=8*N

  call zgeev('V','V',N,A,lda,eig,VL,N,VR,N,work,lwork,rwork,i)

  VL(:,:)=TRANSPOSE(VR(:,:))

  call cinvmatsmooth(VL,N,N,0d0)

!! apply the function
  do k=1,N
     eig(k) = exp( eig(k) )
     CVR(:,k)=VR(:,k)*eig(k)
  enddo
  CVL(:,:)=VL(:,:)

!! rebuild the matrix

  call ZGEMM('N','T',N,N,N,(1d0,0d0),CVR,N,CVL,N,(0d0,0d0),CA,lda)
  A(:,:)=CA(:,:)

end subroutine cexpmat




subroutine rexpmat(A,N,lda)
  implicit none
  integer :: N,lwork,i,k,j,lda
  complex*16 :: eig(N)  ,CVL(N,N),CVR(N,N),CA(lda,N)
  real*8 :: A(lda,N)
  real*8 :: VL(N,N),VR(N,N),work(8*N),rwork(4*N)

  lwork=8*N

  call dgeev('V','V',N,A,lda,rwork(1),rwork(1+N),VL,N,VR,N,work,lwork,i)
  do j=1,N
    eig(j)= (1d0,0d0)*rwork(j) + (0d0,1d0)*rwork(j+N)
  enddo

  VL(:,:)=TRANSPOSE(VR(:,:))

  call realinvmatsmooth(VL,N,N,0d0)

!! apply the function
  do k=1,N
     eig(k) = exp( eig(k) )
     CVR(:,k)=VR(:,k)*eig(k)
  enddo
  CVL(:,:)=VL(:,:)

!! rebuild the matrix.  Is real.

  call ZGEMM('N','T',N,N,N,(1d0,0d0),CVR,N,CVL,N,(0d0,0d0),CA,lda)
  A(:,:)=CA(:,:)  !! Ok conversion

end subroutine rexpmat

!! not used, not debugged
subroutine complexsolve(mat,size,lda,rhs,num,solution)
  implicit none
  integer, intent(in) :: size,lda,num
  complex*16,intent(in) :: mat(lda,size),rhs(lda,num)
  complex*16,intent(out) :: solution(lda,num)
  complex*16 :: tempmat(lda,size)
  integer :: ipiv(size),info

  tempmat(:,:)=mat(:,:)
  solution(:,:)=rhs(:,:)
  call zgesv(size,num,tempmat,lda,ipiv,solution,lda,info)
  if (info.ne.0) then
     print *, "AAAUGH ZGESV INFO ",info
     stop
  endif
end subroutine complexsolve

end module hr_eigenmod


