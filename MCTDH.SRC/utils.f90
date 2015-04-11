
!!
!!  UTILITIES.
!!  

#include "Definitions.INC"

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
  integer :: n,i
  DATATYPE :: one(n), two(n), hermdot, sum
  sum=0.d0
  do i=1,n
     sum =   sum + ALLCON(one(i)) *  two(i) 
  enddo
  hermdot=sum
end function

function floatfac(in)
  implicit none
  integer ::  in, i
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
  integer :: n,i
  DATATYPE :: one(n), two(n), cdot, sum
  sum=0.d0
  do i=1,n
     sum = sum + one(i) * two(i) 
  enddo
  cdot=sum
end function cdot

!! USE THIS FOR CNORM DOT OF DATAECS TYPE
function ecsdot(one,two,n)
  implicit none
  integer :: n,i
  DATAECS :: one(n), two(n), ecsdot, sum
  sum=0.d0
  do i=1,n
     sum = sum + one(i) * two(i) 
  enddo
  ecsdot=sum
end function ecsdot

function myisnan(input)
  implicit none
  real*8 :: input
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
  DATATYPE :: csum,dot,pardot,one(n),two(n)
  integer :: n
  csum=dot(one,two,n)
  call mympireduceone(csum)
  pardot=csum
end function


function realpardot(one,two,n)
  implicit none
  real*8 :: sum,realdot,realpardot,one(n),two(n)
  integer :: n
  sum=realdot(one,two,n)
  call mympirealreduceone(sum)
  realpardot=sum
end function

subroutine realpardotsub(one,two,n,out)
  implicit none
  real*8 :: sum,realdot,out,one(n),two(n)
  integer :: n
  sum=realdot(one,two,n)
  call mympirealreduceone(sum)
  out=sum
end subroutine realpardotsub

subroutine pardotsub(one,two,n,out)
  implicit none
  DATATYPE :: sum,dot,out,one(n),two(n)
  integer :: n
  sum=dot(one,two,n)
  call mympireduceone(sum)
  out=sum
end subroutine pardotsub

function dot(one,two,n)
  implicit none
  integer :: n,i
  DATATYPE :: one(n), two(n), dot, sum
  sum=0.d0
  do i=1,n
     sum = sum + CONJUGATE(one(i)) * two(i) 
  enddo
  dot=sum
end function dot

function realdot(one,two,n)
  implicit none
  integer :: n,i
  real*8 :: one(n), two(n), realdot, sum
  sum=0.d0
  do i=1,n
     sum = sum + one(i) * two(i) 
  enddo
  realdot=sum
end function realdot


subroutine realdotsub(one,two,n,out)
  implicit none
  integer :: n,i
  real*8 :: one(n), two(n), out, sum
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

  integer :: n,m,lda, i
  real*8 :: norm,previous(lda,m), vector(n), realdot
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

  integer :: n,m,lda, i,j
  logical :: parflag
  DATATYPE :: previous(lda,m), vector(n)
  CNORMTYPE :: norm
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
    integer :: size
    DATATYPE :: ket(size),bra(size),heredot,dot,pardot
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

  integer :: n,m,lda, i, kdoneflag
  DATAECS :: previous(lda,m), vector(n), norm, ecsdot
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

function matdet(N,A)
!subroutine getdet(N,A,matdet)
!! this gets a determinant based off of LU decomposition of square matrix A
!! A = PLU : det(A)=det(P)*det(L)*det(U)
!! ?getrf uses Doolittle factorization to get L and U
!! input :
!! N - dimension of A
!! A - an N by N matrix
!! output :
!! matdet - the determinat of A, if anything goes bad, it just returns a value of 0d0
  implicit none
  integer :: N,info,i,pow
  integer :: ipiv(N)
  DATATYPE :: A(N,N),matdet
!! rank 0 or less is not sufficient for a matrix, boo to your input!
  if(N.le.0) then
    matdet=0d0
!! if A is 1x1, ?getrf will just return A, so save the effort
  else if(N.eq.1) then
    matdet=A(1,1)
  else
!! go the actual ?getrf because we need to 

    call MYGETRF(N,N,A,N,ipiv,info)
!! this is the regular case where the factorization worked
!! in Doolittle factorization the diagonal of L=1 => det(L)=1
!! off-diagonal of L and whole of U are stored in A => det(U)=Prod_i{A(i,i)} 
!! det(P)=(-1)**num_row_swaps
!! if ipiv(i).ne.i then row was swapped, last row can't be swapped
    if(info.eq.0) then
      pow=0
      matdet=A(N,N)
      do i=1,N-1
        matdet=matdet*A(i,i)
        if(ipiv(i).ne.i) pow=pow+1
      enddo
      matdet=matdet*((-1d0)**pow)
!! this case is just the matrix is singular
    else if(info.gt.0) then
      matdet=0d0
!! this case is that the ?getrf had an error, so out comes zero
    else 
      matdet=0d0
    endif
  endif
end function matdet


subroutine get_petite_mat(M,N,A,B,left,right)
!! input :
!! M - the dimension of A
!! N - the dimension of B
!! A - an M by M matrix
!! left - an integer list containing all the LHS elements of A to keep in B (N long)
!! right - an integer list containing all the RHS elements of A to keep in B (N long)
!! output : 
!! B - an N by N matrix that is a subset of A
  implicit none
  integer :: i,j,M,N,left(N),right(N)
  DATATYPE :: A(M,M),B(N,N)
  do i=1,N
    do j=1,N
      B(i,j)=A(left(i),right(j))
    enddo
  enddo 
end subroutine get_petite_mat


subroutine neglnmat(A,N,lntol)
  implicit none
  integer :: N
  real*8 :: lntol
  DATATYPE :: A(N,N)
  call bothlnmat(A,N,-1,lntol)
end subroutine neglnmat


subroutine lnmat(A,N,lntol)
  implicit none
  integer :: N
  real*8 :: lntol
  DATATYPE :: A(N,N)
  call bothlnmat(A,N,+1,lntol)
end subroutine lnmat


subroutine bothlnmat(A,N, which ,lntol) 
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=-ln(A_in)
  implicit none
  real*8 :: time,lntol
  integer :: N,lwork,i,j,k,nscale,iflag,ipiv(2*N*N),which
  complex*16 :: eig(N)   !!,saveeig(N)
  DATATYPE :: A(N,N),sum, tempmat(N,N), sum2,sum3,wsp(6*N*N),saveA(N,N)
  integer :: lwsp,ideg,iexph
  DATATYPE :: VL(N,N),VR(N,N),work(8*N),rwork(8*N)
  integer, save :: ierr=0

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

!! get the eigenvalues / eigenvectors
#ifdef REALGO 
  call dgeev('V','V',N,A,N,rwork(1),rwork(1+N),VL,N,VR,N,work,lwork,i)
  do j=1,N
    eig(j)= (1d0,0d0)*rwork(j) + (0d0,1d0)*rwork(j+N)
  enddo
#else
  call zgeev('V','V',N,A,N,eig,VL,N,VR,N,work,lwork,rwork,i)
  VL(1:N,1:N)=ALLCON(VL(1:N,1:N))
#endif

!!  saveeig(:)=eig(:)

!!! normalize the eigenvectors to each other    NO NEED BETTER PROPER SPECTRAL EXPANSION, THE RIGHT WAY, BUT CHECKING.
if (1==0) then
  do j=1,N
    sum=0d0
    do i=1,N
      sum=sum+VL(i,j)*VR(i,j)
    enddo
    VL(:,j)=VL(:,j)/sqrt(sum)
    VR(:,j)=VR(:,j)/sqrt(sum)
  enddo

else

!! MAY 2014
  VL(:,:)=TRANSPOSE(VR(:,:))
  call invmatsmooth(VL,N,N,0d0)

endif



!! apply the function
  do k=1,N
     if (eig(k).ne.0d0) then
        if (abs(eig(k)).lt.lntol) then
           eig(k)= which * log(lntol*eig(k)/abs(eig(k)))
        else
           eig(k) = which * log( eig(k) )
        endif
     else
!        print *,  "BAD! ZERO EIG LN.  FIXME."; call mpistop()

        print *,  "BAD! ZERO EIG LN.  FIXME. TEMP CONTINUE"

        eig(k) = which * log(lntol)

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
  


end subroutine bothlnmat





subroutine expmat(A,N)
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=exp(A_in)
  implicit none
  integer :: N,lwork,i,k,j
  complex*16 :: eig(N)  ,CVL(N,N),CVR(N,N),CA(N,N)
  DATATYPE :: A(N,N)
  DATATYPE :: VL(N,N),VR(N,N),work(8*N),rwork(4*N)

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

end subroutine expmat



subroutine invmatsmooth(A,N,LDA,tol)  !! inverse of ANY matrix.
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=(A_in)**-1
  implicit none
  integer :: N,lwork,i,j,k,LDA
  real*8 :: SV(N),tol
  DATATYPE :: A(LDA,N), SAVEA(LDA,N)
  DATATYPE :: U(N,N),VT(N,N),work(5*N),zwork(5*N)
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

end subroutine invmatsmooth





subroutine realinvmatsmooth(A,N,tol)  !! inverse of ANY matrix.
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=(A_in)**-1
  implicit none
  integer :: N,lwork,i,j,k
  real*8 :: SV(N),tol
  real*8 :: A(N,N)
  real*8 :: U(N,N),VT(N,N),work(5*N)
  lwork=5*N

!! do the svd

  call dgesvd('A','A',N,N,A,N,SV,U,N,VT,N,work,lwork,i)

  if (i.ne.0) then
     print *, "ERR SVD";     stop
  endif
!! apply the function
  do k=1,N
     if (SV(k).eq.0d0) then
        SV(k)=0d0
     else
        if(abs(SV(k)).lt.tol) then
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

end subroutine realinvmatsmooth


!! orthog routines:
!!  FLAG 1:  if on input A(i,j) = <phi_i|phi_j>
!!              then on output varphi_j = sum_i A(i,j) |phi_i>  are orthonormal.
!!  FLAG 2:  if on input A(i,j) = <phi_i|phi_j>
!!              on output varphi_j = sum_i A(i,j) |phi_i> are biorthonormal to phi:
!!            : <varphi_i|phi_j>=delta_ij

subroutine biorthogmat(A,N)
  implicit none
  integer :: N
  DATATYPE :: A(N,N)
#ifdef CNORMFLAG
!!  call cnormorthogmat(A,N,2,0)
  call symorthogmat(A,N,2)
#else
  call symorthogmat(A,N,2)
#endif
end subroutine biorthogmat


subroutine orthogmat(A,N)
  implicit none
  integer :: N
  DATATYPE :: A(N,N)
#ifdef CNORMFLAG
  call cnormorthogmat(A,N,1,0)
!!  call symorthogmat(A,N,1)   no such alg?
#else
  call symorthogmat(A,N,1)
#endif
end subroutine orthogmat


subroutine symorthogmat(A,N,flag)  
  implicit none
  integer :: N,lwork,i,j,k,flag
  real*8, allocatable :: SV(:)
  DATATYPE :: A(N,N),Asave(N,N)
  DATATYPE, allocatable :: U(:,:),VT(:,:),work(:),zwork(:)
  real*8, allocatable :: rwork(:)
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

Asave(:,:)=A(:,:)

!! do the svd
  lwork=5*N
  allocate(U(N,N),VT(N,N),SV(N),work(lwork),zwork(lwork),rwork(lwork))
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
  deallocate(U,VT,SV,work,zwork,rwork)

end subroutine symorthogmat


subroutine cnormorthogmat(ovl,num,flag,printflag)
  implicit none
  DATATYPE :: ovl(num,num), dot
  integer :: num,i,j,printflag,k, flag
  DATATYPE :: ovlvals(num),ovlvects(num,num)
  if ((flag.ne.1).and.(flag.ne.2)) then
     print *, "bad flag",flag;     stop
  endif
#ifndef CNORMFLAG
  print *, "ERR why did you call cnormorthogmat?  eigen_cnorm won't do what you want it to"
  stop
#endif
  do i=1,num
     do j=1,i
        if (abs((ovl(i,j))-ovl(j,i)).gt.1.d-7) then
           print *, "CNORM ERR ORTHOG", abs((ovl(i,j))-ovl(j,i))
           stop
        endif
     enddo
  enddo

print *,  "REDO SPECTRAL EXPANSION."; call mpistop()

#ifdef REALGO
  call get_eigen_two(ovl,num,num,ovlvects,ovlvals)
#else
  call get_eigen_cnorm(ovl,num,num,ovlvects,ovlvals)
#endif
  do i=1,num
     do j=1,i-1
        if (abs(dot(ovlvects(:,i),ovlvects(:,j),num)).gt.1.d-7) then
           print *, " ERR CNORM ", abs(dot(ovlvects(:,i),ovlvects(:,j),num)),i,j,ovlvals(i),ovlvals(j)
        endif
     enddo
  enddo
  if (printflag.eq.1) then
     do i=1,num
        if (abs(sqrt(ovlvals(i))-1.d0).gt.0.01d0) then
            write(*,*) "Large ovl!! ", ovlvals(i)
        endif
     enddo
  endif
  ovl=0.d0
  do k=1,num
     if (abs(ovlvals(k)).lt.1d-10) then
        ovlvals(k)=0.d0
        print *, "CSING!"
     else if (flag.eq.1) then
        ovlvals(k)=1.d0/sqrt(ovlvals(k))
     else
        ovlvals(k)=1.d0/ovlvals(k)
     endif
     do i=1,num
        do j=1,num
!! inverse square root or inverse
           ovl(i,j) = ovl(i,j) + ovlvects(i,k)*ovlvects(j,k) * ovlvals(k)
        enddo
     enddo
  enddo


end subroutine cnormorthogmat



! enter as l1 l2 m1 m2 l3    TIMES TWO!!   INTEGER ARGUMENTS FOR HALF SPIN

function doubleclebschsq (l2,l1,m2,m1,l3)
  implicit none

  integer :: l1,l2,m1,m2,l3,ierr,l3min,l3max

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

  




