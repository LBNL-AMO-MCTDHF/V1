      DOUBLE PRECISION FUNCTION D1MACH (I)
      IMPLICIT NONE
      INTEGER :: I
      DOUBLE PRECISION :: B, X
      X = 1.0D0
      
      B = RADIX(X)
      SELECT CASE (I)
      CASE (1)
         D1MACH = B**(MINEXPONENT(X)-1) ! the smallest positive magnitude.
      CASE (2)
         D1MACH = HUGE(X)       ! the largest magnitude.
      CASE (3)
         D1MACH = B**(-DIGITS(X)) ! the smallest relative spacing.
      CASE (4)
         D1MACH = B**(1-DIGITS(X)) ! the largest relative spacing.
      CASE (5)
         D1MACH = log10(B)
      CASE DEFAULT
         WRITE (*, FMT = 9000)
 9000    FORMAT ('1ERROR    1 IN D1MACH - I OUT OF BOUNDS')
         STOP
      END SELECT

!!      write(*,*)I,  D1MACH

      RETURN
      
      end function d1mach
