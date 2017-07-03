
#include "Definitions.INC"

!! these are meant to be #included not compiled directly!  need definitions header.  
!! don't include definitions here.  
!! uses only DATAECS, MYGETRF, MYGETRI, possibly ECSFLAG.


module myprojectmod
  implicit none

  real*8, allocatable :: jacobideriv(:, :, :) , jacobirhoderiv(:,:),jacobiKE(:, :, :)  , &
       jacobiweights(:), jacobipoints(:),glpoints2d(:,:), glweights2d(:,:),ylmvals(:,:,:)

!! for YLMS - eigenvalues of each T_r^(l)
!! or for DVR - first entry only used; eigenvalues of T_r

!! for DVR - eigvals of jacobi for each m

  DATAECS, allocatable, target :: ddrhopot(:,:)

  DATAECS, allocatable :: glke(:,:,:), glfirstdertot(:,:,:),glrhoderivs(:,:)
  DATAECS, allocatable :: glpoints(:), glweights(:),zdipole(:,:), xydipole(:,:)
  DATATYPE, allocatable, target :: &
       sparseops_xi_banded(  :,   :,    :,   :), &
       sparseops_eta(  :,   :, :,   :),       sparseops_diag(   :,  :,  : )

  DATATYPE, allocatable, target :: &
       sparseddz_xi_banded(  :,   :,    :,     :), &
       sparseddz_eta(  :,   :,    :,     :),        sparseddz_diag(   :,    :,     :), &
       sparseddrho_xi_banded(  :,   :,    :,     :),        sparseddrho_eta(  :,   :,    :,     :), &
       sparseddrho_diag(   :,    :,     :) !!, &

  DATAECS, allocatable :: rmatrix(:,:,:,:)

!!! Eq.68, McCurdy Baertschy Rescigno
!!! for atom
!  DATAECS :: rmatrix(numerad,numerad,mseriesmax+1,lseriesmax+1)   

  real*8, parameter :: rpoints(1)=1d0

end module myprojectmod


subroutine myprojectalloc()
  use myparams
  use myprojectmod
  implicit none
  
  allocate( jacobiKE(lbig+1, lbig+1, 0:mbig) , jacobideriv(lbig+1, lbig+1, 0:mbig) ,jacobirhoderiv(lbig+1,lbig+1), &
       jacobiweights(lbig+1), jacobipoints(lbig+1) ,glke(hegridpoints,hegridpoints,0:1) , &
       glfirstdertot(hegridpoints,hegridpoints,0:1), glrhoderivs(hegridpoints,hegridpoints), &
       glpoints(hegridpoints), glweights(hegridpoints), glpoints2d(henumpoints,2), glweights2d(henumpoints,2) , &
       xydipole(numerad,lbig+1),zdipole(numerad,lbig+1), ddrhopot(numerad,lbig+1));
  jacobike=0; jacobideriv=0; jacobirhoderiv=0; jacobiweights=0; jacobipoints=0; glke=0; glfirstdertot=0; 
  glrhoderivs=0; glpoints=0; glweights=0; glpoints2d=0; glweights2d=0; xydipole=0; zdipole=0; ddrhopot=0

  allocate(sparseddz_xi_banded(2*bandwidth+1,numerad,lbig+1,mbig+1), &
       sparseddz_eta(  lbig+1,lbig+1,numerad,mbig+1 ) ,  &
       sparseddz_diag(  numerad,lbig+1,mbig+1 )  , &
       sparseddrho_xi_banded(2*bandwidth+1,numerad,lbig+1,mbig+1), &
       sparseddrho_eta(  lbig+1,lbig+1,numerad,mbig+1 ) ,  &
       sparseddrho_diag(  numerad,lbig+1,mbig+1 ) , &
       sparseops_xi_banded(2*bandwidth+1,numerad,lbig+1,mbig+1), &
       sparseops_eta(  lbig+1,lbig+1,numerad,mbig+1), &
       sparseops_diag(  numerad,lbig+1,mbig+1),&
       rmatrix(numerad,numerad,mseriesmax+1,lseriesmax+1)  , &
       ylmvals(0:2*mbig, 1:lbig+1, lseriesmax+1))

  sparseddz_xi_banded = 0;       sparseddz_eta = 0;       sparseddz_diag = 0;  
  sparseddrho_xi_banded = 0;         sparseddrho_eta = 0;         sparseddrho_diag = 0;  
  sparseops_xi_banded = 0;         sparseops_eta = 0;         sparseops_diag = 0;  
  rmatrix=0    ;   ylmvals=0

end subroutine myprojectalloc


subroutine myprojectdealloc()
  implicit none
end subroutine myprojectdealloc

module gettwoemod
contains

subroutine get_twoe_new()
  use myparams
  use myprojectmod  

  implicit none
  DATAECS, allocatable :: kearray(:,:,:), invkearray(:,:,:)
  integer :: lwork, info,i,k,ii,j, j1a, lsum, deltam
  DATAECS, allocatable :: work(:)
  integer, allocatable :: ipiv(:)
  real*8 :: sum

  lwork=hegridpoints*(lbig+1)*(2*mbig+1) * 10
  allocate(work(lwork),ipiv(lwork),kearray(numerad,numerad,lseriesmax+1), &
       invkearray(numerad,numerad,lseriesmax+1))
  work=0; ipiv=0; kearray=0; invkearray=0

  do i=0,lseriesmax
    ii=numerad
    kearray(:,:,i+1) = glke(2:hegridpoints-1,2:hegridpoints-1,0)
    invkearray(:,:,i+1) = kearray(:,:,i+1)
    do j=1,numerad
       invkearray(j,j,i+1) = invkearray(j,j,i+1) - i*(i+1) / glpoints(j+1)**2
    enddo

    call MYGETRF(ii,ii,invkearray(:,:,i+1),ii,ipiv,info)
    if (info/=0) then
       print *, "dgetrf info ", info;       stop
    endif
    call MYGETRI(ii,invkearray(:,:,i+1),ii,ipiv,work,lwork,info)
    if (info/=0) then
       print *, "dgetri info ", info;       stop
    endif

    do j=1,numerad
    do k=1,numerad
       rmatrix(j,k,1,i+1) = (2*i+1) / (glpoints(j+1)*glpoints(k+1)*&
            sqrt(glweights(j+1)*glweights(k+1))) * invkearray(j,k,i+1) &
            - (glpoints(j+1)*glpoints(k+1))**i / glpoints(hegridpoints)**(2*i+1)
    enddo
    enddo
  enddo


  ylmvals=0.d0
  do deltam   = 0,2*mbig
     do j1a=0,lbig   
        lsum = abs(deltam)-1
        do while (lsum .lt. jacobisummax)
           lsum=lsum+1
           
           call jacobinormed(lsum-abs(deltam),abs(deltam),abs(deltam),jacobipoints(j1a+1),sum)
           ylmvals(deltam,j1a+1,lsum+1-abs(deltam)) = sum * &
                sqrt( ( 2.d0 ) / ( 2.0d0 * lsum + 1.d0 ) ) * &
                sqrt( ( 1.d0 - jacobipoints(j1a+1)**2 ) )**abs(deltam)
           
        enddo
     enddo
  enddo

  deallocate(work,ipiv, kearray, invkearray)

end subroutine get_twoe_new

end module


subroutine op_yderiv(howmanynotused,notusedin,notusedout)
  use myparams
  implicit none
  integer :: howmanynotused
  DATATYPE :: notusedin(*), notusedout(*)
  OFLWR "WHAT! no op_yderiv atom."; CFLST
  notusedout(1)=notusedin(1)*howmanynotused
end subroutine op_yderiv

