
#define DATATYPE complex*16
#define REALFLAG 0

#define OFLWR write(*,*)
#define CFL
#define mpifileptr *
#define CFLST stop
#define OFL

#include "harmonic_ritz_small.F90"

!
!  CALL THIS SUBROUTINE (AT END OF FILE)
!
!subroutine harmonic_ritz( &
!     lanblocknum, numvects, lansize,order, maxiter, outvectors, outvalues,&
!     inprintflag,guessflag,lanthresh,&
!     multsub,indotsub,   etarget,range, &
!     invmult_order, invexpmult_order, expmult_order)




!!! ON INPUT
!!
!!   lanblocknum        size of krylov block; number of eigenvectors calculated
!!   numvects           number of those to converge within lanthresh tolerance
!!   lansize            size of problem per processor
!!   order              Krylov space dimension
!!   maxiter            size of problem
!!   outvectors         if guessflag=1, guess vectors    outvectors(lansize,lanblocknum)
!!   inprintflag        print flag
!!   guessflag          use outvectors for initial krylov block?
!!   lanthresh          convergence criterion
!!   multsub            matrix vector multiplication subroutine.  arbitrary complex nonsymmetric.
!!                      multiples NN vectors.  must be programmed as 
!!                           multsub(in(lansize,NN),out(lansize,NN),NN,lansize)
!!   indotsub           dot product subroutine.
!!                      computes HERMITIAN dot products between a set of NN vectors and one of MM vectors.
!!                           indotsub(bravectors(lansize,NN),ketvectors(lansize,MM),lansize,NN,MM,result(NN,MM)
!!   etarget            TARGET ENERGY. 
!!   range              RANGE (in units of energy).  small range will cause breakdown of gram schmidt proceedure.
!!   invmult_order      NOT USED  (please set to -1)
!!   invexpmult_order   Krylov space dimension for subroutine that generates each member of main Krylov space
!!   expmult_order      Krylov space dimension for subroutine that generates initial Krylov vector
!!                        if guessflag=0
!!! ON OUTPUT
!!
!!   outvectors    eigenvectors
!!   outvalues     eigenvalues
