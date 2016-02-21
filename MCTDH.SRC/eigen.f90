
!! LAPACK WRAPPERS

!! REMOVED ALL LEFT RIGHT EIGENVECTOR ROUTINES, NOT NEEDED, AND USE POOR SPECTRAL EXPANSION

#include "Definitions.INC"

subroutine mysort(in, values,n,lda,out)
  implicit none
  complex*16,intent(in) :: in(lda,n)
  complex*16,intent(out) :: out(lda,n)
  complex*16,intent(inout) :: values(lda)
  integer :: lda, n,i,j,whichlowest, flag
  real*8 :: lowval 
  integer,allocatable :: taken(:), order(:)
  complex*16,allocatable :: newvals(:)

  allocate(taken(lda), order(lda),newvals(lda))

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
         write(*,*) "WHICHLOWEST ERROR, J=",j," WHICHLOWEST=", whichlowest;        
         call mpistop()
     endif
     if (taken(whichlowest).ne.0) then
        write(*,*) "TAKEN ERROR.";        call mpistop()
     endif
     taken(whichlowest)=1;     order(j)=whichlowest
     out(:,j)=in(:,order(j));     newvals(j)=values(order(j))
  enddo

  values(1:n)=newvals(1:n)

  deallocate(taken,order,newvals)

end subroutine mysort

!! this wrapper orthogonalizes WRT dot(): it is for 
!!  complex symmetric cmctdh, or complex unsymmetric 
!!  chmctdh matrices like configuration matrices.  
!!  for cmctdh, may be linearly dep if degenerate.  could gs-orthog somehow.

subroutine get_eigen_c(inmat, n, lda, outmat,values)
  implicit none
  integer,intent(in) :: n,lda
  complex*16,intent(in) :: inmat(lda,n)
  complex*16,intent(out) :: outmat(lda,n), values(n)

!! 080711 naah  082911 back again
!! 112110 now sending to cnorm for better orthog

#ifdef CNORMFLAG
  call get_eigen_cnorm(inmat,n,lda,outmat,values)
#else

!!! should be left hermitian normed in c0
  call get_eigen_c0(inmat, n, lda, outmat,values)

#endif
end subroutine get_eigen_c


!! C-normed.

subroutine get_eigen_cnorm(inmat, n, lda, outmat,values)
  implicit none
  integer,intent(in) :: n,lda
  complex*16,intent(in) :: inmat(lda,n)
  complex*16,intent(out) :: outmat(lda,n), values(n)
  real*8, parameter :: lowthresh=1.d-8
  real*8 :: valnorm, rrealdot
  integer :: i,low

  call get_eigen_c0(inmat, n, lda, outmat,values)

!! for degeneracy, partial gram schmidt. should be sorted by real value. 

  valnorm=max(abs(values(1)),abs(values(n)))

  low=1
  do i=1,n

! this is screwing stuff up  DANGER DANGER .... IS NEEDED sometimes
!  (attempt to deal with linear dependence in degenerate eigs)
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
  integer,intent(in) :: n,lda
  complex*16,intent(in) :: inmat(lda,n)
  complex*16,intent(out) :: outmat(lda,n), values(lda)
  complex*16,allocatable :: vr(:,:),vl(:,:),tempmat(:,:),work(:)
  real*8,allocatable :: rwork(:)
  integer :: info, lwork
  character*1 :: jobvl="V", jobvr="V"

  allocate(vr(lda,n),vl(lda,n),tempmat(lda,n),work(5*lda),rwork(5*lda))
  vr=0; vl=0; tempmat=0; work=0; rwork=0

  lwork=5*lda;  tempmat(1:n,1:n)=inmat(1:n,1:n)

  call zgeev( jobvl, jobvr, n, tempmat, lda, values, vl,lda,vr,lda,&
       work, lwork, rwork, info)
  if (info /= 0) then
     write(*, *) "AAAUGH, info in get_eigen_c:",info;     
     call mpistop()
  endif

  call mysort(vr, values ,n,lda,outmat)

  deallocate(vr,vl,tempmat,work,rwork)

end subroutine get_eigen_c0




subroutine get_eigen_two(inmat, n, lda, outvects, outvals)
  implicit none
  integer,intent(in) :: n,lda
  real*8,intent(in) :: inmat(lda,n)
  real*8,intent(out) :: outvects(lda,n), outvals(n)
  integer :: info, lwork
  character*1 :: job="V", uplo="U"
  real*8,allocatable :: work(:)

  allocate(work(5*lda)); work=0
  lwork=5*lda;  outvects(:,:) = inmat(:,:)
  call dsyev( job, uplo, n, outvects, lda, outvals, work, lwork, info)
  if (info /= 0) then
       write(*, *)"AAAUGH, in get_eigen_two dsyev info = ",info;    
       call mpistop()
  endif
  deallocate(work)

end subroutine 


subroutine get_eigen_tri(indiag,insubdiag, n, lda, outvects, outvals)
  implicit none
  integer,intent(in) :: n,lda
  real*8,intent(in) :: indiag(n), insubdiag(n-1)
  real*8,intent(out) :: outvals(n), outvects(lda,n)
  integer :: info
  character*1 :: job="V"
  real*8 :: work(3*n), subdiag(n-1)    !! AUTOMATIC
  outvals = indiag;  subdiag = insubdiag
  call dstev( job, n, outvals, subdiag, outvects, lda, work, info )
  if (info /= 0) then
    write(*, *)"AAAUGH, in get_eigen_tri dstev info = ",info;    
    call mpistop()
  endif
end subroutine 


subroutine get_eigen_herm(inmat, n, lda, outmat,values)
  implicit none
  integer,intent(in) :: n,lda
  complex*16,intent(in) :: inmat(lda,n)
  complex*16,intent(out) :: outmat(lda,n)
  real*8,intent(out) :: values(lda)
  integer :: info, lwork
  character*1 :: jobz="V", uplo="U"
  complex*16,allocatable :: work(:)
  real*8,allocatable :: rwork(:)

  allocate(work(5*lda),rwork(5*lda))
  work=0; rwork=0
  lwork=5*lda;  outmat(:,:)=inmat(:,:)
  call zheev( jobz, uplo, n, outmat, lda, values, work, lwork, rwork, info)
  if (info /= 0) then
    write(*, *) "AAAUGH, info in get_eigen_herm:",info;    
    call mpistop()
  endif
  deallocate(work,rwork)

end subroutine 



!! for case with no complex eigvals

subroutine get_eigen_asym(inmat, n, lda, leftvects, rightvects,values)
  implicit none
  integer,intent(in) :: n,lda
  real*8,intent(in) :: inmat(lda,n)
  real*8,intent(out) :: leftvects(lda,n), rightvects(lda,n), values(lda)
  real*8,allocatable :: tempmat(:,:), ivalues(:), work(:)
  integer :: i, info, lwork

  lwork=10*n
  allocate( tempmat(lda,n), ivalues(lda), work(10*n))
  ivalues=0; work=0
  tempmat(1:n,1:n)=inmat(1:n,1:n)

!      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
!     $                  LDVR, WORK, LWORK, INFO )

  call dgeev( 'V', 'V', n, tempmat, lda, values, ivalues, leftvects,&
       lda,rightvects,lda,work, lwork, info)
  if (info /= 0) then
     print *, "AAAUGH, info in get_eigen_asym:",info;     stop
  endif
  do i=1,n
     if (ivalues(i)/=0.d0) then
        print *, "IMAG EIGEN!"; stop
     endif
  enddo

  deallocate( tempmat, ivalues, work)

end subroutine get_eigen_asym


subroutine get_eigen_gen(inmat, ovlmat, n, lda, outvects, outvals)
  implicit none
  integer,intent(in) :: n,lda
  complex*16,intent(in) :: inmat(lda,n), ovlmat(lda,n)
  complex*16,intent(out) :: outvects(lda,n), outvals(n)
  integer :: info, lwork
  character*1 :: jobvr="N", jobvl="V" !! uplo="U"
  complex*16,allocatable :: work(:),tempovl(:,:), &
       alpha(:), beta(:), tempvects(:,:)
  real*8,allocatable :: rwork(:)

  allocate(work(5*lda),tempovl(lda,n), alpha(n), beta(n), &
       tempvects(lda,n),rwork(8*n))
  lwork=5*lda;  outvects(:,:) = inmat(:,:);  tempovl(:,:)=ovlmat(:,:)
  work=0; alpha=0; beta=0; tempvects=0; rwork=0
 
  call zggev( jobvl,jobvr, n, outvects, lda, tempovl,lda, alpha,beta, &
       tempvects, lda, outvects,lda, work, lwork, rwork, info)

  outvals=alpha/beta

  if (info /= 0) then
     print *, "AAAUGH, in get_eigen_gen zggev info = ",info;stop
  endif

  deallocate(work,tempovl, alpha, beta, tempvects,rwork)

end subroutine 






