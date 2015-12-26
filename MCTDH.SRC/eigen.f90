
!! LAPACK WRAPPERS

!! REMOVED ALL LEFT RIGHT EIGENVECTOR ROUTINES, NOT NEEDED, AND USE POOR SPECTRAL EXPANSION

#include "Definitions.INC"

subroutine mysort(in, values,n,lda,out)
  implicit none
  complex*16 :: in(lda,n), out(lda,n), values(lda)
  integer :: lda, n,i,j,whichlowest, flag
  real*8 :: lowval 
  integer :: taken(lda), order(lda)
  complex*16 :: newvals(lda)

  taken=0;  order=-1
  do j=1,n
     whichlowest=-1; flag=0;     lowval=1.d+30  !! is not used (see flag)
     do i=1,n
        if ( taken(i) .eq. 0 ) then
           if ((flag.eq.0) .or.(real(values(i)) .le. lowval)) then
              flag=1;              lowval=real(values(i)); whichlowest=i
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
     out(:,j)=in(:,order(j));     newvals(j)=values(order(j))
  enddo

  values(1:n)=newvals(1:n)

end subroutine mysort

!! this wrapper orthogonalizes WRT dot(): it is for complex symmetric cmctdh, or complex unsymmetric 
!!  chmctdh matrices like configuration matrices.  
!!  for cmctdh, may be linearly dep if degenerate.  could gs-orthog somehow.

subroutine get_eigen_c(inmat, n, lda, outmat,values)
  implicit none

  integer :: n,lda
  complex*16 :: inmat(lda,n), outmat(lda,n), values(n)
!  complex*16, allocatable :: onemat(:,:)

!! 080711 naah  082911 back again
!! 112110 now sending to cnorm for better orthog

#ifdef CNORMFLAG
  call get_eigen_cnorm(inmat,n,lda,outmat,values)
#else

!!! should be left hermitian normed in c0
  call get_eigen_c0(inmat, n, lda, outmat,values)

!  allocatE(onemat(lda,n)); onemat(:,:)=0d0
!  do i=1,n
!     onemat(i,i)=1d0
!  enddo
!  call get_eigen_gen(inmat, onemat, n, lda, outmat, values)
!  deallocate(onemat)
#endif
end subroutine get_eigen_c


!! C-normed.

subroutine get_eigen_cnorm(inmat, n, lda, outmat,values)
  implicit none

  integer :: n,lda,i, low
  complex*16 :: inmat(lda,n), outmat(lda,n), values(n) !!, mydot
  real*8, parameter :: lowthresh=1.d-8
  real*8 :: valnorm, rrealdot

  call get_eigen_c0(inmat, n, lda, outmat,values)

!! for degeneracy, partial gram schmidt. should be sorted by real value. 

  valnorm=max(abs(values(1)),abs(values(n)))

  low=1
  do i=1,n

! this is screwing stuff up  DANGER DANGER .... IS NEEDED sometimes
!     do j=low,i-1
!        if (real(values(j)-values(i),8).gt.lowthresh*valnorm) then
!           low=low+1
!           if (low.gt.j) then
!!              print *, "LOWERROR!!"; stop  !! temp training wheels
!           endif
!        else
!           mydot=DOT_PRODUCT(conjg(outmat(:,j)),outmat(:,i))
!           outmat(:,i)=outmat(:,i)-mydot*outmat(:,j)
!        endif
!     enddo

     outmat(:,i)=outmat(:,i)/sqrt(DOT_PRODUCT(CONJG(outmat(:,i)),outmat(:,i)))
     rrealdot=sqrt(real(DOT_PRODUCT(outmat(:,i),outmat(:,i))))
!     if (rrealdot.gt.1.d5) then
!        print *, "large!", rrealdot
!        outmat(:,i)=outmat(:,i)*1.d3/rrealdot
!     endif
  enddo


end subroutine get_eigen_cnorm




!! leave herm normed as given by zgeev:

subroutine get_eigen_c0(inmat, n, lda, outmat,values)
  implicit none
  
  complex*16 :: inmat(lda,n), outmat(lda,n), values(lda), vr(lda,n),vl(lda,n),tempmat(lda,n)
  integer :: n,lda, info, lwork
  character*1 :: jobvl="V", jobvr="V"
  complex*16 :: work(5*lda)
  real*8 :: rwork(5*lda)

  lwork=5*lda;  tempmat(1:n,1:n)=inmat(1:n,1:n)

  call zgeev( jobvl, jobvr, n, tempmat, lda, values, vl,lda,vr,lda,work, lwork, rwork, info)
  if (info /= 0) then
     write(*, *) "AAAUGH, info in get_eigen_c:",info;     call mpistop()
  endif

  call mysort(vr, values ,n,lda,outmat)

end subroutine get_eigen_c0




subroutine get_eigen_two(inmat, n, lda, outvects, outvals)
  implicit none
  real*8 :: inmat(lda,n), outvects(lda,n), outvals(n)
  integer :: n,lda, info, lwork
  character*1 :: job="V", uplo="U"
  real*8 :: work(5*lda)

  lwork=5*lda;  outvects = inmat
  call dsyev( job, uplo, n, outvects, lda, outvals, work, lwork, info)
  if (info /= 0) then
       write(*, *)"AAAUGH, in get_eigen_two dsyev info = ",info;     call mpistop()
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
    write(*, *)"AAAUGH, in get_eigen_tri dstev info = ",info;     call mpistop()
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
    write(*, *) "AAAUGH, info in get_eigen_c:",info;     call mpistop()
  endif
end subroutine 



!! for case with no complex eigvals

subroutine get_eigen_asym(inmat, n, lda, leftvects, rightvects,values)
  implicit none
  real*8 :: inmat(lda,n), leftvects(lda,n), rightvects(lda,n), values(lda)
  real*8 :: tempmat(lda,n), ivalues(lda), work(10*n)
  integer :: n,lda,i, info, lwork

  lwork=10*n
  
  tempmat(1:n,1:n)=inmat(1:n,1:n)

!      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
!     $                  LDVR, WORK, LWORK, INFO )

  call dgeev( 'V', 'V', n, tempmat, lda, values, ivalues, leftvects,lda,rightvects,lda,work, lwork, info)
  if (info /= 0) then
     print *, "AAAUGH, info in get_eigen_asym:",info;     stop
  endif
  do i=1,n
     if (ivalues(i)/=0.d0) then
        print *, "IMAG EIGEN!"; stop
     endif
  enddo

end subroutine get_eigen_asym


subroutine get_eigen_gen(inmat, ovlmat, n, lda, outvects, outvals)
  implicit none
  complex*16 :: inmat(lda,n), outvects(lda,n), outvals(n), ovlmat(lda,n)
  integer :: n,lda, info, lwork
  character*1 :: jobvr="N", jobvl="V" !! uplo="U"
  complex*16 :: work(5*lda),tempovl(lda,n), alpha(n), beta(n), tempvects(lda,n)
  real*8 :: rwork(8*n)

  lwork=5*lda;  outvects = inmat;  tempovl=ovlmat

  call zggev( jobvl,jobvr, n, outvects, lda, tempovl,lda, alpha,beta, tempvects, lda, outvects,lda, work, lwork, rwork, info)

  outvals=alpha/beta

  if (info /= 0) then
     print *, "AAAUGH, in get_eigen_gen zggev info = ",info;stop
  endif

end subroutine 






