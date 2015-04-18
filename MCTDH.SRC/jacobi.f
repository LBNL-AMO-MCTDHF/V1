

      function fac ( in )
      integer in, i
      real*8 fac

      fac = 1.0d0

      do i=1,in
         fac = fac*i
      enddo


      end function

      subroutine jacobinormed( n, myalpha, mybeta, x, cx )
      implicit none

      integer*4 n,alpha,beta,myalpha,mybeta
      real*8 fac

      double precision  alphadouble, betadouble
      double precision  cx,x
!!$      logical checknan2real

      alpha=abs(myalpha)
      beta=abs(mybeta)
      alphadouble=real(alpha)
      betadouble=real(beta)

c      write(*,'(10F17.13)'), 69 ,n,alpha,fac(4)
      call jacobi_poly( n, alphadouble, betadouble, x, cx )

!!$      if (checknan2real(cx,1)) then
!!$         print *, "NAN jacobi!"
!!$         stop
!!$      endif


!      cx = cx/  sqrt(real(2**(alpha+beta+1) * fac( n + alpha ) * 
!     & fac ( n + beta ),8) / real(( 2 * n + alpha + beta+1)  * 
!     & fac( n ) * fac ( n + alpha + beta  ),8))

      cx = cx/  sqrt(real(2.d0**(alpha+beta+1) * fac( n + alpha ) * 
     & fac ( n + beta ),8) / real(( 2 * n + alpha + beta+1)  * 
     & fac( n ) * fac ( n + alpha + beta  ),8))

!!$      if (checknan2real(cx,1)) then
!!$         print *, "NAN jacobi! xx"
!!$         print *, alpha+beta+1
!!$         print *, 2**(alpha+beta+1)
!!$         print *,  fac( n + alpha ) 
!!$         print *, fac ( n + beta )
!!$         print *, ( 2 * n + alpha + beta+1) 
!!$         print *, fac( n ) 
!!$         print *,  fac ( n + alpha + beta  )
!!$         stop
!!$      endif

      end subroutine

      
      subroutine jacobi_poly ( n, alpha, beta, x, cx )
      
!******************************************************************************
!     
!     ! JACOBI_POLY evaluates the Jacobi polynomials at X.
!     
!     Differential equation:
!     
!     (1-X*X) Y'' + (BETA-ALPHA-(ALPHA+BETA+2) X) Y' + N (N+ALPHA+BETA+1) Y = 0
!     
!     Recursion:
!     
!     P(0,ALPHA,BETA,X) = 1,
!     
!     P(1,ALPHA,BETA,X) = ( (2+ALPHA+BETA)*X + (ALPHA-BETA) ) / 2
!     
!     P(N,ALPHA,BETA,X)  = 
!     ( 
!     (2*N+ALPHA+BETA-1) 
!     * ((ALPHA**2-BETA**2)+(2*N+ALPHA+BETA)*(2*N+ALPHA+BETA-2)*X) 
!     * P(N-1,ALPHA,BETA,X)
!     -2*(N-1+ALPHA)*(N-1+BETA)*(2*N+ALPHA+BETA) * P(N-2,ALPHA,BETA,X)
!     ) / 2*N*(N+ALPHA+BETA)*(2*N-2+ALPHA+BETA)
!     
!     Restrictions:
!     
!     -1 < ALPHA
!     -1 < BETA
!     
!     Norm:
!     
!     Integral ( -1 <= X <= 1 ) ( 1 - X )**ALPHA * ( 1 + X )**BETA 
!     * P(N,ALPHA,BETA,X)**2 dX 
!     = 2**(ALPHA+BETA+1) * Gamma ( N + ALPHA + 1 ) * Gamma ( N + BETA + 1 ) /
!     ( 2 * N + ALPHA + BETA ) * N! * Gamma ( N + ALPHA + BETA + 1 )
!     
!     Special values:
!     
!     P(N,ALPHA,BETA)(1) = (N+ALPHA)!/(N!*ALPHA!) for integer ALPHA.
!     
!     Modified:
!     
!     01 October 2002
!     
!     Author:
!     
!     John Burkardt
!     
!     Reference:
!     
!     Milton Abramowitz, Irene Stegun,
!     Handbook of Mathematical Functions,
!     US Department of Commerce, 1964.
!     
!     Parameters:
!     
!     Input, integer N, the highest order polynomial to compute.  Note
!     that polynomials 0 through N will be computed.
!     
!     Input, double precision ( kind = 8 ) ALPHA, one of the 
!     parameters defining the Jacobi
!     polynomials, ALPHA must be greater than -1.
!     
!     Input, double precision ( kind = 8 ) BETA, the second 
!     parameter defining the Jacobi
!     polynomials, BETA must be greater than -1.
!     
!     Input, double precision ( kind = 8 ) X, the point at which 
!     the polynomials are 
!     to be evaluated.
!     
!     Output, double precision ( kind = 8 ) CX, the value
!     of the Jacobi polynomial of order N at the point X.
!     
      implicit none

      integer n

      double precision  alpha
      double precision  beta
      double precision  cx,cxtemp1,cxtemp2
      double precision  c1
      double precision  c2
      double precision  c3
      double precision  c4
      integer i
      double precision  r_i
      double precision  x

      if ( alpha <= -1.0D+00 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'JACOBI_POLY - Fatal error!'
         write ( *, '(a,g14.6)' ) '  Illegal input value of ALPHA = '
     1        , alpha
         write ( *, '(a)' ) '  But ALPHA must be greater than -1.'
         stop
      end if
      
      if ( beta <= -1.0D+00 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'JACOBI_POLY - Fatal error!'
         write ( *, '(a,g14.6)' ) '  Illegal input value of BETA = '
     1        , beta
         write ( *, '(a)' ) '  But BETA must be greater than -1.'
         stop
      end if
      
      if ( n < 0 ) then
         return
      end if

      cxtemp1 = 1.0D+00

      if ( n == 0 ) then
         cx = cxtemp1      
         return
      end if
 
      cxtemp2 = ( 1.0D+00 + 0.5D+00 * ( alpha + beta ) ) * x 
     1     + 0.5D+00 * ( alpha - beta )

      if(n == 1) then
         cx = cxtemp2
         return
      end if
      
      do i = 2, n

         r_i =dble(i)

         c1 = 2.0D0 * r_i * ( r_i + alpha + beta ) 
     1        * ( 2.0D+00 * r_i - 2.0D+00 + alpha + beta )

         c2 = ( 2.0D+00 * r_i - 1.0D+00 + alpha + beta ) 
     1        * ( 2.0D+00 * r_i  + alpha + beta ) 
     2        * ( 2.0D+00 * r_i - 2.0D+00 + alpha + beta )

         c3 = ( 2.0D+00 * r_i - 1.0D+00 + alpha + beta ) 
     1        * ( alpha + beta ) * ( alpha - beta )

         c4 = 2.0D+00 * ( r_i - 1.0D+00 + alpha ) 
     1        * ( r_i - 1.0D+00 + beta )  * ( 2.0D+00 * r_i + alpha 
     2        + beta )

         cx = ( ( c3 + c2 * x ) * cxtemp2 - c4 * cxtemp1 ) / c1

         cxtemp1 = cxtemp2
         cxtemp2 = cx
      end do

      return
      end



