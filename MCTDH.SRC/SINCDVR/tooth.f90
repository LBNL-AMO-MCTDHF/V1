MODULE TOOTH
IMPLICIT NONE

CONTAINS

!--------------------------------------------------------------------------------------

SUBROUTINE GETINVERSE(TINVOUT,N,NBIG,NSMALL,DELTA)

! COMPUTES TINV, THE INVERSE OF T, ON THE BIG GRID
! EACH DIMENSION OF TINV IS INDEXED FROM -NBIG TO NBIG
! LIKEWISE FOR THE SMALL GRID

INTEGER, INTENT(IN) :: N, NBIG, NSMALL
REAL *8, INTENT(IN) :: DELTA
REAL *8, INTENT(OUT) :: TINVOUT(-N:N,-N:N,-N:N)
REAL *8 :: TINV(-NBIG:NBIG,-NBIG:NBIG,-NBIG:NBIG)
INTEGER :: I, J, K, INFO, LWORK, NS
REAL *8 :: TR(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL)
REAL *8 :: Q(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL)
REAL *8 :: T(-NSMALL:NSMALL,-NSMALL:NSMALL)
!REAL *8 :: TBIG(-NBIG:NBIG,-NBIG:NBIG)
!REAL *8 :: TRBIG(-NBIG:NBIG,-NBIG:NBIG,-NBIG:NBIG)
REAL *8 :: PI
REAL *8, ALLOCATABLE :: WORK(:), W(:,:), WT(:,:)
REAL *8, ALLOCATABLE :: E(:), EIGS(:,:,:), QQ(:,:,:)

! COMPUTE RHS = 2*PI - T*R
! STORE R IN TINV

PI = 4D0*ATAN(1D0)

!!$  CALL RVECTOR(TINV,TR,N,NBIG,NSMALL,DELTA)

CALL RVECTOR(TINV,TR,NBIG,NSMALL,DELTA)

TR(0,0,0) = TR(0,0,0) + 2*PI

! FILL IN ENTRIES OF T_SMALL

DO J = -NSMALL,NSMALL
DO I = -NSMALL,NSMALL

   IF(I == J) THEN
   T(I,J) = (PI/DELTA)**2/6D0
   ELSE
   T(I,J) = (-1D0)**(I-J)/(DELTA**2*(I-J)**2)
   END IF

END DO
END DO

! SOLVE T*Q = 2*PI - T*R FOR Q

!GOTO 40

NS = 2*NSMALL+1
LWORK = 6*NS
ALLOCATE(WORK(LWORK),W(NS,NS),WT(NS,NS),E(NS))
ALLOCATE(EIGS(NS,NS,NS),QQ(NS,NS,NS))

W = T
CALL DSYEV('V','U',NS,W,NS,E,WORK,LWORK,INFO)
WT = TRANSPOSE(W)

DO K = 1,NS
DO J = 1,NS
DO I = 1,NS

   EIGS(I,J,K) = E(I) + E(J) + E(K)

END DO
END DO
END DO

CALL TOOTHEIGMULT(TR,QQ,WT,NS)
QQ = QQ/EIGS
CALL TOOTHEIGMULT(QQ,Q,W,NS)

DO K = -NSMALL,NSMALL
DO J = -NSMALL,NSMALL
DO I = -NSMALL,NSMALL

   TINV(I,J,K) = TINV(I,J,K) + Q(I,J,K)

END DO
END DO
END DO

DEALLOCATE(WORK,W,WT,E,EIGS,QQ)

TINVOUT=TINV(-N:N,-N:N,-N:N)

END SUBROUTINE GETINVERSE

!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------

!!SUBROUTINE RVECTOR(R,TR,N,NBIG,NSMALL,DELTA)

SUBROUTINE RVECTOR(R,TR,NBIG,NSMALL,DELTA)

! COMPUTES THE VECTORS R(I,J,K) = 1/|R_{IJK}|
! AND TR = -T_{BIG}*R

INTEGER, INTENT(IN) :: NBIG, NSMALL
REAL *8, INTENT(IN) :: DELTA
REAL *8, INTENT(OUT) :: R(-NBIG:NBIG,-NBIG:NBIG,-NBIG:NBIG)
REAL *8, INTENT(OUT) :: TR(-NSMALL:NSMALL,-NSMALL:NSMALL,-NSMALL:NSMALL)

INTEGER :: I,J,K,P
REAL *8 :: T1,T2,T3,RAD
REAL *8 :: PI,X,Y,Z !!,A,C,XX(-NBIG:NBIG)

PI = 4D0*ATAN(1D0)

! COMPUTE R VECTOR

!$OMP PARALLEL DO PRIVATE(I,J,K,RAD) SHARED(NBIG,DELTA,R)

DO K = -NBIG,NBIG
DO J = -NBIG,NBIG
DO I = -NBIG,NBIG

   X = I*DELTA; Y = J*DELTA; Z = K*DELTA 
   RAD = SQRT(X**2 + Y**2 + Z**2)
   
   IF(RAD == 0) THEN
      R(I,J,K) = 0
   ELSE
      R(I,J,K) = DELTA**3/RAD
   END IF

END DO
END DO
END DO
!$OMP END PARALLEL DO

! COMPUTE TR VECTOR

TR = 0
!$OMP PARALLEL DO PRIVATE(I,J,K,P,T1,T2,T3) SHARED(DELTA,PI,TR)
DO K = -NSMALL,NSMALL
DO J = -NSMALL,NSMALL
DO I = -NSMALL,NSMALL

   DO P = -NBIG,NBIG
   T1 = PI**2/(6D0*DELTA**2)
   T2 = T1
   T3 = T1
   
   IF(I /= P) T1 = (-1D0)**(I-P)/(DELTA**2*(I-P)**2)
   IF(J /= P) T2 = (-1D0)**(J-P)/(DELTA**2*(J-P)**2)
   IF(K /= P) T3 = (-1D0)**(K-P)/(DELTA**2*(K-P)**2)
   
   TR(I,J,K) = TR(I,J,K) - (T1*R(P,J,K) + T2*R(I,P,K)+ T3*R(I,J,P))

   END DO

END DO
END DO
END DO
!$OMP END PARALLEL DO

END SUBROUTINE RVECTOR


!--------------------------------------------------------------------------------------

SUBROUTINE TOOTHEIGMULT(X,Y,W,NN)

INTEGER, INTENT(IN) :: NN
REAL *8, INTENT(IN) :: X(1:NN,1:NN,1:NN), W(1:NN,1:NN)
REAL *8, INTENT(OUT) :: Y(1:NN,1:NN,1:NN)
REAL *8 :: X1(1:NN,1:NN,1:NN), X2(1:NN,1:NN,1:NN)
INTEGER :: I,J,K

DO J = 1,NN
DO I = 1,NN
X1(I,J,:) = MATMUL(W,X(I,J,:))
END DO
END DO

DO K = 1,NN
DO I = 1,NN
X2(I,:,K) = MATMUL(W,X1(I,:,K))
END DO
END DO

DO K = 1,NN
DO J = 1,NN
Y(:,J,K) = MATMUL(W,X2(:,J,K))
END DO
END DO

END SUBROUTINE TOOTHEIGMULT

!--------------------------------------------------------------------------------------

SUBROUTINE FOURIER_INVERSE(TINV,NB)

INTEGER, INTENT(IN) :: NB
REAL *8, INTENT(OUT) :: TINV(-NB:NB,-NB:NB,-NB:NB)

INTEGER :: I,J,K
REAL *8 :: PI
COMPLEX *16 :: X(-NB:NB,-NB:NB,-NB:NB), Y(-NB:NB,-NB:NB,-NB:NB)
COMPLEX *16 :: T1,T2,T3,T(-NB:NB)

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

END SUBROUTINE FOURIER_INVERSE

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

COMPLEX *16 :: Y1(-N:N,-N:N,-N:N)
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

COMPLEX *16 :: Y1(-N:N,-N:N,-N:N)
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



