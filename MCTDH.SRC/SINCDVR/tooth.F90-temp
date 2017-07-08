
#include "Definitions.INC"


MODULE TOOTH
IMPLICIT NONE

CONTAINS

!--------------------------------------------------------------------------------------

RECURSIVE SUBROUTINE GETINVERSE(TINVOUT,N,NBIG,NSMALL,DELTA)
  USE pfileptrmod
  IMPLICIT NONE
! COMPUTES TINV, THE INVERSE OF T, ON THE BIG GRID
! EACH DIMENSION OF TINV IS INDEXED FROM -NBIG TO NBIG
! LIKEWISE FOR THE SMALL GRID

INTEGER, INTENT(IN) :: N, NBIG, NSMALL
REAL *8, INTENT(IN) :: DELTA
REAL *8, INTENT(OUT) :: TINVOUT(-N:N,-N:N,-N:N)
REAL *8 :: PI, MYCON
REAL *8, ALLOCATABLE :: WORK(:), W(:,:), WT(:,:), E(:), EIGS(:,:,:), QQ(:,:,:),&
     TINV(:,:,:),TR(:,:,:),Q(:,:,:),T(:,:)
INTEGER :: I, J, K, INFO, LWORK, NS

OFLWR "        ..allocating for getinverse"; CFL

allocate(TINV(-NBIG:NBIG,-NBIG:NBIG,-NBIG:NBIG), TR(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL),&
     Q(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL), T(-NSMALL:NSMALL,-NSMALL:NSMALL))
TINV=0; TR=0; Q=0; T=0

OFLWR "      ..ok allocated"; CFL

! COMPUTE RHS = 2*PI - T*R;    STORE R IN TINV

PI = 4D0*ATAN(1D0)
MYCON = (PI/DELTA)**2/6D0

OFLWR "        ..calling rvector"; CFL
CALL RVECTOR(TINV,TR,NBIG,NSMALL,DELTA)
OFLWR "        ..called rvector"; CFL

TR(0,0,0) = TR(0,0,0) + 2*PI

! FILL IN ENTRIES OF T_SMALL

DO J = -NSMALL,NSMALL
   DO I = -NSMALL,NSMALL
      T(I,J) = (-1D0)**(I-J)/(DELTA**2*(I-J)**2)
   END DO
   T(J,J) = MYCON
END DO

! SOLVE T*Q = 2*PI - T*R FOR Q

NS = 2*NSMALL+1
LWORK = 6*NS
ALLOCATE(WORK(LWORK),W(NS,NS),WT(NS,NS),E(NS))
ALLOCATE(EIGS(NS,NS,NS),QQ(NS,NS,NS))
WORK=0; W=0; WT=0; E=0; EIGS=0; QQ=0

W = T
CALL DSYEV('V','U',NS,W,NS,E,WORK,LWORK,INFO)
IF (INFO.NE.0) THEN
   OFLWR "ERROR DSYEV GETINVERSE"; CFLST
ENDIF

DO K = 1,NS
DO J = 1,NS
   EIGS(:,J,K) = E(:) + E(J) + E(K)
END DO
END DO

OFLWR "        ..calling eigmult"; CFL
WT = TRANSPOSE(W)
CALL TOOTHEIGMULT(TR,QQ,WT,NS)
QQ = QQ/EIGS
CALL TOOTHEIGMULT(QQ,Q,W,NS)

OFLWR "        ..called eigmult"; CFL

TINV(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL) = &
     TINV(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL) + &
     Q(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL)

DEALLOCATE(WORK,W,WT,E,EIGS,QQ)

TINVOUT(:,:,:)=TINV(-N:N,-N:N,-N:N)

deallocate(TINV,TR,Q,T)

OFLWR "        ..done getinverse"; CFL

RETURN

!--------------------------------------------------------------------------------------

CONTAINS
RECURSIVE SUBROUTINE RVECTOR(R,TR,NBIG,NSMALL,DELTA)
  IMPLICIT NONE
! COMPUTES THE VECTORS R(I,J,K) = 1/|R_{IJK}|
! AND TR = -T_{BIG}*R

  INTEGER, INTENT(IN) :: NBIG, NSMALL
  REAL *8, INTENT(IN) :: DELTA
  REAL *8, INTENT(OUT) :: R(-NBIG:NBIG,-NBIG:NBIG,-NBIG:NBIG), &
       TR(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL)
  INTEGER :: I,J,K
  REAL *8 :: PI,T1,T2,T3,DELTACUBED,DELTASQUAREDINV
  REAL *8, ALLOCATABLE :: PPOS(:),       ALLTT(:), MYR(:,:,:)
  INTEGER, ALLOCATABLE :: IPOS(:)

  if (NBIG.LE.NSMALL) THEN
     print *, "NBIG is less than nsmall.  error", nbig,nsmall
     call mpistop()
  endif

  R=0; TR=0
  
  PI = 4D0*ATAN(1D0)

! COMPUTE R VECTOR

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,PPOS,IPOS,ALLTT,DELTACUBED,DELTASQUAREDINV,MYR)

  ALLOCATE( PPOS(-NSMALL-NBIG:NSMALL+NBIG), ALLTT(-NSMALL-NBIG:NSMALL+NBIG), &
       IPOS(-NSMALL-NBIG:NSMALL+NBIG), MYR(-NBIG:NBIG,-NBIG:NBIG,-NBIG:NBIG) )

  PPOS=0; ALLTT=0; IPOS=0 

  DELTACUBED=DELTA**3
  DO I = -NSMALL-NBIG,NSMALL+NBIG
     IPOS(I) = I
     PPOS(I) = I*DELTA
  ENDDO

!$OMP DO SCHEDULE(STATIC)
  DO K = -NBIG,NBIG
  DO J = -NBIG,NBIG
     R(:,J,K) = DELTACUBED/SQRT(PPOS(-NBIG:NBIG)**2 + PPOS(J)**2 + PPOS(K)**2)
  END DO
  END DO
!$OMP END DO

  R(0,0,0) = 0d0

  MYR(:,:,:) = R(:,:,:);    !! LOCAL COPY...HMM

! COMPUTE TR VECTOR

  DELTASQUAREDINV=1d0/DELTA**2

  ALLTT(:) = (-1)**(IPOS(:)) * DELTASQUAREDINV/(IPOS(:))**2

  ALLTT(0) = PI**2/(6D0*DELTA**2)

!$OMP DO SCHEDULE(STATIC)
  DO K = -NSMALL,NSMALL
  DO J = -NSMALL,NSMALL
  DO I = -NSMALL,NSMALL

     TR(I,J,K) = (-1) * SUM( &
          ALLTT(-NBIG-I:NBIG-I)*MYR(:,J,K) + &
          ALLTT(-NBIG-J:NBIG-K)*MYR(I,:,K) + &
          ALLTT(-NBIG-K:NBIG-K)*MYR(I,J,:))

!     DO P = -NBIG,NBIG
!        T1 = AMOUNT
!        T2 = AMOUNT
!        T3 = AMOUNT
!        IF(I /= P) T1 = (-1D0)**(I-P)*DELTASQUAREDINV/(I-P)**2
!        IF(J /= P) T2 = (-1D0)**(J-P)*DELTASQUAREDINV/(J-P)**2
!        IF(K /= P) T3 = (-1D0)**(K-P)*DELTASQUAREDINV/(K-P)**2
!        TR(I,J,K) = TR(I,J,K) - (T1*R(P,J,K) + T2*R(I,P,K)+ T3*R(I,J,P))
!     END DO

  END DO
  END DO
  END DO
!$OMP END DO

  DEALLOCATE(ALLTT,PPOS,IPOS,MYR)

!$OMP END PARALLEL

END SUBROUTINE RVECTOR

!--------------------------------------------------------------------------------------

RECURSIVE SUBROUTINE TOOTHEIGMULT(X,Y,W,NN)
IMPLICIT NONE
INTEGER, INTENT(IN) :: NN
REAL *8, INTENT(IN) :: X(1:NN,1:NN,1:NN), W(1:NN,1:NN)
REAL *8, INTENT(OUT) :: Y(1:NN,1:NN,1:NN)
REAL *8,allocatable :: X1(:,:,:), X2(:,:,:)
REAL *8 :: MYW(1:NN,1:NN)                   !! AUTOMATIC
INTEGER :: I,J,K

allocate(X1(1:NN,1:NN,1:NN), X2(1:NN,1:NN,1:NN))
X1=0; X2=0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,MYW)

MYW(:,:) = W(:,:)

DO J = 1,NN
DO I = 1,NN
X1(I,J,:) = MATMUL(MYW,X(I,J,:))
END DO
END DO
!$OMP DO SCHEDULE(STATIC)
DO K = 1,NN
DO I = 1,NN
X2(I,:,K) = MATMUL(MYW,X1(I,:,K))
END DO
END DO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
DO K = 1,NN
DO J = 1,NN
Y(:,J,K) = MATMUL(MYW,X2(:,J,K))
END DO
END DO
!$OMP END DO

!$OMP END PARALLEL

deallocate(X1,X2)

END SUBROUTINE TOOTHEIGMULT

END SUBROUTINE GETINVERSE

!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------

SUBROUTINE xFOURIER_INVERSE(TINV,NB)

INTEGER, INTENT(IN) :: NB
REAL *8, INTENT(OUT) :: TINV(-NB:NB,-NB:NB,-NB:NB)

INTEGER :: I,J,K
REAL *8 :: PI
COMPLEX *16,allocatable :: X(:,:,:), Y(:,:,:), T(:)
COMPLEX *16 :: T1,T2,T3

allocate( X(-NB:NB,-NB:NB,-NB:NB), Y(-NB:NB,-NB:NB,-NB:NB),T(-NB:NB) )
X=0; Y=0; T=0

PI = 4D0*ATAN(1D0)

DO J = -NB,NB
   T(J) = ((-1D0)**J)/J**2
   IF(J == 0) T(J) = PI**2/6D0
END DO

DO K = -NB,NB
DO J = -NB,NB
DO I = -NB,NB
   T1 = 0; T2 = 0; T3 = 0
   IF(J == 0 .AND. K == 0) T1 = T(I)
   IF(I == 0 .AND. K == 0) T2 = T(J)
   IF(I == 0 .AND. J == 0) T3 = T(K)
   X(I,J,K) = T1+T2+T3
END DO
END DO
END DO


CALL DFT3(X,Y,NB)
CALL IDFT3(1D0/Y,X,NB)
TINV = 2D0*PI*DREAL(X)

deallocate( X,Y,T )

CONTAINS

! ---------------------------------------------------------------

SUBROUTINE DFT1(X,Y,N)

INTEGER, INTENT(IN) :: N
COMPLEX *16, INTENT(IN) :: X(-N:N)
COMPLEX *16, INTENT(OUT) :: Y(-N:N)
INTEGER :: J,K,NN
COMPLEX *16 :: A
REAL *8 :: PI

NN = 2*N+1
PI = 4D0*ATAN(1D0)
A = -2D0*PI*DCMPLX(0D0,1D0)/NN

DO K = -N,N
   Y(K) = 0
   DO J = -N,N
      Y(K) = Y(K) + X(J)*EXP(A*J*K)
   END DO
END DO

END SUBROUTINE DFT1

! ---------------------------------------------------------------

SUBROUTINE IDFT1(X,Y,N)

INTEGER, INTENT(IN) :: N
COMPLEX *16, INTENT(IN) :: X(-N:N)
COMPLEX *16, INTENT(OUT) :: Y(-N:N)
INTEGER :: J,K,NN
COMPLEX *16 :: A
REAL *8 :: PI

NN = 2*N+1
PI = 4D0*ATAN(1D0)
A = 2D0*PI*DCMPLX(0D0,1D0)/NN

DO K = -N,N
   Y(K) = 0
   DO J = -N,N
      Y(K) = Y(K) + X(J)*EXP(A*J*K)/NN
   END DO
END DO

END SUBROUTINE IDFT1

! ---------------------------------------------------------------

SUBROUTINE DFT3(X,Y,N)

INTEGER, INTENT(IN) :: N
COMPLEX *16, INTENT(IN) :: X(-N:N,-N:N,-N:N)
COMPLEX *16, INTENT(OUT) :: Y(-N:N,-N:N,-N:N)
COMPLEX *16 :: Y1(-N:N,-N:N,-N:N)    !! AUTOMATIC
INTEGER :: I,J

DO J = -N,N
DO I = -N,N
   CALL DFT1(X(:,I,J),Y(:,I,J),N)
END DO
END DO

DO J = -N,N
DO I = -N,N
   CALL DFT1(Y(I,:,J),Y1(I,:,J),N)
END DO
END DO

DO J = -N,N
DO I = -N,N
   CALL DFT1(Y1(I,J,:),Y(I,J,:),N)
END DO
END DO

END SUBROUTINE DFT3

! ---------------------------------------------------------------

SUBROUTINE IDFT3(X,Y,N)

INTEGER, INTENT(IN) :: N
COMPLEX *16, INTENT(IN) :: X(-N:N,-N:N,-N:N)
COMPLEX *16, INTENT(OUT) :: Y(-N:N,-N:N,-N:N)
COMPLEX *16 :: Y1(-N:N,-N:N,-N:N)   !! AUTOMATIC
INTEGER :: I,J

DO J = -N,N
DO I = -N,N
   CALL IDFT1(X(:,I,J),Y(:,I,J),N)
END DO
END DO

DO J = -N,N
DO I = -N,N
   CALL IDFT1(Y(I,:,J),Y1(I,:,J),N)
END DO
END DO

DO J = -N,N
DO I = -N,N
   CALL IDFT1(Y1(I,J,:),Y(I,J,:),N)
END DO
END DO

END SUBROUTINE IDFT3

END SUBROUTINE xFOURIER_INVERSE

! ---------------------------------------------------------------

END MODULE TOOTH




!!$
!!$SUBROUTINE GETFIRSTINVERSE(FINV,N,NBIG,NSMALL,DELTA)
!!$
!!$INTEGER, INTENT(IN) :: N, NBIG, NSMALL
!!$REAL *8, INTENT(IN) :: DELTA
!!$REAL *8, INTENT(OUT) :: FINV(-NBIG:NBIG,-NBIG:NBIG,-NBIG:NBIG)
!!$
!!$INTEGER :: I, J, K, INFO, LWORK, NS
!!$REAL *8 :: Q(-NBIG:NBIG,-NBIG:NBIG,-NBIG:NBIG)
!!$REAL *8 :: T(-NBIG:NBIG,-NBIG:NBIG)
!!$REAL *8 :: PI
!!$REAL *8, ALLOCATABLE :: WORK(:), W(:,:), E(:)
!!$
!!$
!!$PI = 4D0*ATAN(1D0)
!!$
!!$DO J = -NBIG,NBIG
!!$DO I = -NBIG,NBIG
!!$
!!$   IF(I == J) THEN
!!$   T(I,J) = (PI/DELTA)**2/6D0
!!$   ELSE
!!$   T(I,J) = (-1D0)**(I-J)/(DELTA**2*(I-J)**2)
!!$   END IF
!!$
!!$END DO
!!$END DO
!!$
!!$
!!$NS = 2*NBIG+1
!!$LWORK = 6*NS
!!$ALLOCATE(WORK(LWORK),W(-NBIG:NBIG,-NBIG:NBIG),E(-NBIG:NBIG))
!!$
!!$W = T
!!$CALL DSYEV('V','U',NS,T,NS,E,WORK,LWORK,INFO)
!!$W = TRANSPOSE(T)
!!$
!!$DO K = -NBIG,NBIG
!!$DO J = -NBIG,NBIG
!!$DO I = -NBIG,NBIG
!!$
!!$   Q(I,J,K) = W(0,I)*W(0,J)*W(0,K)/SQRT(E(I)+E(J)+E(K))
!!$   
!!$END DO
!!$END DO
!!$END DO
!!$
!!$
!!$CALL TOOTHEIGMULT(Q,FINV,W,NS)
!!$
!!$
!!$
!!$END SUBROUTINE GETFIRSTINVERSE
!!$
!!$!--------------------------------------------------------------------------------------
!!$
!!$SUBROUTINE GETINVERSE3(TINV,N,NBIG,NSMALL,DELTA)
!!$
!!$! COMPUTES TINV, THE INVERSE OF T, ON THE BIG GRID
!!$! EACH DIMENSION OF TINV IS INDEXED FROM -NBIG TO NBIG
!!$! LIKEWISE FOR THE SMALL GRID
!!$
!!$INTEGER, INTENT(IN) :: N, NBIG, NSMALL
!!$REAL *8, INTENT(IN) :: DELTA
!!$REAL *8, INTENT(OUT) :: TINV(-NBIG:NBIG,-NBIG:NBIG,-NBIG:NBIG)
!!$
!!$INTEGER :: I, J, K, NS, NSUM
!!$REAL *8 :: Q(-NBIG:NBIG,-NBIG:NBIG,-NBIG:NBIG)
!!$REAL *8 :: X(-NBIG:NBIG,-NBIG:NBIG), E(-NBIG:NBIG)
!!$REAL *8 :: PI,A,S
!!$
!!$NS = 2*NBIG+1
!!$PI = 4D0*ATAN(1D0)
!!$A = PI/NS
!!$NSUM = 1000
!!$
!!$DO J = -NBIG,NBIG
!!$DO I = -NBIG,NBIG
!!$   X(I,J) = SIN(A*(I+NBIG+1)*(J+NBIG+1))
!!$END DO
!!$END DO
!!$
!!$X = X*SQRT(2D0/NS)
!!$
!!$DO J = -NBIG,NBIG
!!$   S = 0
!!$   DO I = 1,NBIG
!!$      S = S + (-1)**I*COS(A*I*(J+NBIG+1))/I**2
!!$   END DO
!!$   E(J) = PI**2/6D0 + 2D0*S
!!$END DO
!!$
!!$DO K = -NBIG,NBIG
!!$DO J = -NBIG,NBIG
!!$DO I = -NBIG,NBIG
!!$ 
!!$   Q(I,J,K) = X(I,0)*X(J,0)*X(K,0)/(E(I)+E(J)+E(K))
!!$
!!$END DO
!!$END DO
!!$END DO
!!$
!!$CALL TOOTHEIGMULT(Q,TINV,X,NS)
!!$TINV = 2D0*PI*DELTA**2*TINV
!!$
!!$END SUBROUTINE GETINVERSE3
!!$




!!$
!!$!--------------------------------------------------------------------------------------
!!$
!!$SUBROUTINE CONJGRAD(T,U,V,NN)
!!$
!!$! SOLVES THE SYSTEM T*U = V WITH CONJUGATE GRADIENT
!!$
!!$INTEGER, INTENT(IN) :: NN
!!$REAL *8, INTENT(IN) :: T(1:NN,1:NN)
!!$REAL *8, INTENT(INOUT) :: V(1:NN,1:NN,1:NN)
!!$REAL *8, INTENT(OUT) :: U(1:NN,1:NN,1:NN)
!!$
!!$INTEGER :: K, MAXITER
!!$REAL *8 :: R(NN,NN,NN), ROLD(NN,NN,NN), P(NN,NN,NN), Q(NN,NN,NN), A, B, TOL, RNORM
!!$REAL *8 :: U1(NN,NN,NN)
!!$
!!$! INITIALIZE
!!$
!!$MAXITER = 10000
!!$TOL = 1D-12
!!$U = 0
!!$
!!$CALL KEMULTREAL(T,U,P,NN)
!!$R = V - P
!!$RNORM = SQRT(SUM(R*R))
!!$P = R
!!$K = 0
!!$
!!$! CG ITERATION
!!$
!!$DO WHILE(RNORM > TOL)
!!$
!!$   CALL KEMULTREAL(T,P,Q,NN)
!!$   A = SUM(R*R)/SUM(P*Q)
!!$   U = U + A*P
!!$   ROLD = R
!!$   R = R - A*Q
!!$   RNORM = SUM(R*R)
!!$   B = RNORM/SUM(ROLD*ROLD)
!!$   P = R + B*P
!!$   RNORM = SQRT(RNORM)
!!$   K = K+1
!!$
!!$END DO
!!$
!!$END SUBROUTINE CONJGRAD
!!$
!!$!--------------------------------------------------------------------------------------
!!$
!!$SUBROUTINE KEMULTREAL(T,U,V,NN)
!!$
!!$! COMPUTES THE STRICTLY REAL MATRIX-VECTOR PRODUCT V = T*U
!!$
!!$INTEGER, INTENT(IN) :: NN
!!$REAL *8, INTENT(IN) :: T(1:NN,1:NN), U(1:NN,1:NN,1:NN)
!!$REAL *8, INTENT(OUT) :: V(1:NN,1:NN,1:NN)
!!$INTEGER :: I,J,K
!!$
!!$DO K = 1,NN
!!$DO J = 1,NN
!!$DO I = 1,NN
!!$
!!$   V(I,J,K) = SUM(T(:,I)*U(:,J,K)) &
!!$            + SUM(T(:,J)*U(I,:,K)) &
!!$            + SUM(T(:,K)*U(I,J,:))
!!$
!!$END DO
!!$END DO
!!$END DO
!!$
!!$END SUBROUTINE KEMULTREAL
!!$
!!$!--------------------------------------------------------------------------------------
!!$
!!$SUBROUTINE KETOOTHCLUSTER(T,N,NS,W,DELTA)
!!$
!!$INTEGER, INTENT(IN) :: N, NS
!!$REAL *8, INTENT(OUT) :: T(-NS:NS,-NS:NS)
!!$REAL *8, INTENT(IN) :: W, DELTA
!!$INTEGER :: I,J
!!$REAL *8 :: X(-NS:NS), D(-NS:NS,-NS:NS), R1(-NS:NS),A,C
!!$
!!$A = N*DELTA
!!$C = A/SINH(W*A)
!!$
!!$DO I = -NS,NS
!!$   X(I) = C*SINH(W*I*DELTA)
!!$END DO
!!$
!!$R1 = 1D0/(C*W*SQRT(1+X**2/C**2))
!!$
!!$D = 0
!!$DO J = -NS,NS
!!$DO I = -NS,NS
!!$
!!$   IF(I /= J) D(I,J) = (-1)**(I-J)/(DELTA*(I-J)) 
!!$
!!$END DO
!!$END DO
!!$
!!$DO J = -NS,NS
!!$DO I = -NS,NS
!!$
!!$   T(I,J) = 0.5*SUM(R1**2*D(:,I)*D(:,J))
!!$
!!$END DO
!!$END DO
!!$
!!$END SUBROUTINE KETOOTHCLUSTER




!!$ This never made sense
!!$SUBROUTINE FVECTOR(TINV,FINV,NBIG,NSMALL,DELTA)
!!$
!!$  INTEGER, INTENT(IN) :: NBIG, NSMALL
!!$  REAL *8, INTENT(IN) :: DELTA
!!$  REAL *8, INTENT(IN) :: TINV(-NBIG:NBIG,-NBIG:NBIG,-NBIG:NBIG)
!!$  REAL *8, INTENT(OUT) :: FINV(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL,3)
!!$  INTEGER :: I,J,K,P
!!$  REAL *8 :: T1,T2,T3,DELTAINV
!!$
!!$  if (NBIG.LE.NSMALL) THEN
!!$     print *, "NBIG is less than nsmall.  error", nbig,nsmall
!!$     call mpistop()
!!$  endif
!!$  
!!$  FINV(:,:,:,:) = 0
!!$
!!$!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,P,T1,T2,T3,DELTAINV)
!!$  
!!$ 
!!$  DELTAINV=1d0/DELTA
!!$  
!!$  DO P = -NBIG,NBIG
!!$
!!$!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
!!$     DO K = -NSMALL,NSMALL
!!$     DO J = -NSMALL,NSMALL
!!$     DO I = -NSMALL,NSMALL
!!$
!!$        T1 = 0d0
!!$        T2 = 0d0
!!$        T3 = 0d0
!!$        
!!$        IF(I /= P) T1 = (-1D0)**(I-P)*DELTAINV/(I-P)
!!$        IF(J /= P) T2 = (-1D0)**(J-P)*DELTAINV/(J-P)
!!$        IF(K /= P) T3 = (-1D0)**(K-P)*DELTAINV/(K-P)
!!$        
!!$        FINV(I,J,K,1) = FINV(I,J,K,1) + T1*TINV(P,J,K)
!!$        FINV(I,J,K,2) = FINV(I,J,K,2) + T2*TINV(I,P,K)
!!$        FINV(I,J,K,3) = FINV(I,J,K,3) + T3*TINV(I,J,P)
!!$        
!!$     END DO
!!$     END DO
!!$     END DO
!!$!$OMP END DO
!!$  END DO
!!$
!!$!$OMP BARRIER
!!$!$OMP END PARALLEL
!!$
!!$
!!$END SUBROUTINE FVECTOR
!!$
!!$
!!$
!!$!--------------------------------------------------------------------------------------
!!$
!!$SUBROUTINE FVECTOR2(TINV,FINVOUT,FINVSQOUT,NBIG,NSMALL,DELTA)
!!$
!!$  INTEGER, INTENT(IN) :: NBIG, NSMALL
!!$  REAL *8, INTENT(IN) :: DELTA
!!$  REAL *8, INTENT(IN) :: TINV(-NBIG:NBIG,-NBIG:NBIG,-NBIG:NBIG)
!!$  REAL *8, INTENT(OUT) :: FINVOUT(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL,3)
!!$  REAL *8, INTENT(OUT) :: FINVSQOUT(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL,3)
!!$  REAL *8 :: WORK(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL)
!!$  REAL *8 :: FINV(-NBIG:NBIG,-NBIG:NBIG,-NBIG:NBIG,3)  !! AUTOMATIC
!!$  INTEGER :: I,J,K,P,myrank,nprocs
!!$  REAL *8 :: T1,T2,T3,DELTAINV
!!$  REAL*8 :: aaa,bbb,ccc
!!$
!!$  if (NBIG.LE.NSMALL) THEN
!!$     print *, "NBIG is less than nsmall.  error", nbig,nsmall
!!$     call mpistop()
!!$  endif
!!$  
!!$  FINV(:,:,:,:) = 0
!!$
!!$!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,P,T1,T2,T3,DELTAINV)
!!$  
!!$  DELTAINV=1d0/DELTA
!!$  
!!$  DO P = -NBIG,NBIG
!!$
!!$!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
!!$     DO K = -NBIG,NBIG
!!$     DO J = -NBIG,NBIG
!!$     DO I = -NBIG,NBIG
!!$
!!$        T1 = 0d0
!!$        T2 = 0d0
!!$        T3 = 0d0
!!$        
!!$        IF(I /= P) T1 = (-1D0)**(I-P)*DELTAINV/(I-P)
!!$        IF(J /= P) T2 = (-1D0)**(J-P)*DELTAINV/(J-P)
!!$        IF(K /= P) T3 = (-1D0)**(K-P)*DELTAINV/(K-P)
!!$        
!!$        FINV(I,J,K,1) = FINV(I,J,K,1) + T1*TINV(P,J,K)
!!$        FINV(I,J,K,2) = FINV(I,J,K,2) + T2*TINV(I,P,K)
!!$        FINV(I,J,K,3) = FINV(I,J,K,3) + T3*TINV(I,J,P)
!!$        
!!$     END DO
!!$     END DO
!!$     END DO
!!$!$OMP END DO
!!$  END DO
!!$
!!$!$OMP BARRIER
!!$!$OMP END PARALLEL
!!$
!!$  CALL CONVOLVE_3D(FINV,FINV,FINVSQOUT,NBIG,NSMALL,3)
!!$
!!$  FINVOUT(:,:,:,:)=FINV(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL,:)
!!$
!!$  call getmyranknprocs(myrank,nprocs)
!!$  if (myrank.eq.1) then
!!$     print *, "CHECKING FVECTOR2."
!!$     WORK(:,:,:) = FINVSQOUT(:,:,:,1)+FINVSQOUT(:,:,:,1)+FINVSQOUT(:,:,:,3)
!!$
!!$     aaa= DOT_PRODUCT(RESHAPE(WORK(:,:,:),(/(2*nsmall+1)**3/)),&
!!$          RESHAPE(WORK(:,:,:),(/(2*nsmall+1)**3/)))
!!$     bbb= DOT_PRODUCT(RESHAPE(TINV(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL),(/(2*nsmall+1)**3/)),&
!!$          RESHAPE(TINV(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL),(/(2*nsmall+1)**3/)))
!!$     ccc= DOT_PRODUCT(RESHAPE(TINV(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL),(/(2*nsmall+1)**3/)),&
!!$          RESHAPE(WORK(:,:,:),(/(2*nsmall+1)**3/)))
!!$     print *, aaa,bbb,ccc,aaa*bbb/ccc/ccc
!!$     print *, "TEMPSTOP."
!!$  endif
!!$  call mpibarrier(); call mpistop()
!!$
!!$END SUBROUTINE FVECTOR2
!!$

