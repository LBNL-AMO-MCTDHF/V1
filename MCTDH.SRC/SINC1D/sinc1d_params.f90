
!! PARAMETERS FILE FOR 1-D SINC DVR BASIS AND HAMILTONIAN:
!!            sinc1dparinp namelist input.

#include "Definitions.INC"

module myparams
implicit none

!! HAMILTONIAN PARAMS

real*8  :: twostrength=1d0   !!  2-electron interaction coef.
real*8  :: nucstrength=1d0   !!  2-nuclei interaction coef.
integer :: twotype = 1       !!  Two-particle interaction:
                             !!    1 = potential interaction  0 = constant interaction
integer :: twomode = 0       !!  One and (if twotype=1) two particle interaction:
!                            !!    0 = sech^2 potential  1 = soft coulomb  2 = linear
integer :: coulmode = -1     !!  For soft coulomb: -1: plain ad-hoc softened coulomb with 
!                            !!    unphysical even parity, softness determines groundstate
!                            !!    0,1: 1d 3d atomic analog with nonlocal centrifugal
!                            !!    potential even parity but molecular not figured out
!                            !!    0: integer quantum numbers even parity (3d analog)
!                            !!    1: half-integer
real*8 :: softness=1d0       !!  scale parameter for sech and softening param for coulomb
                             !!    softcoul = 1/sqrt(x^2+softness^2)
real*8 :: softnesstwoe=1d0,& !!  for coulmode = (-1), additional parameters for ee and nn
     softnessnuc=1d0         !!
! For twomode = 0 : sech-squared 1d model for atoms
integer :: sechmode = 1      !!  0:  n=0 "s" energy is correct
!                            !!      softness determines range of potential.
!                            !!      electron-nucleus potentials simply add
!                            !!  1:  orbital energies n = 0 and 1 ("s" and "p") are
!                            !!      correct; softness variable not used
integer :: combinesech = 1   !!  Add potentials to get united atom limit correct
integer :: nucmode = 0       !!  For sech-squared: 0=physically motivated repulsion
!                            !!  1=ad hoc repulsion to set R=0 energy for hypothetical
!                            !!      noninteracting boson problem, all electrons in
!                            !!      lowest orbital equal to asymptote energy
integer :: nucrangefac = 1d0 !! For nucmode=1, ad hoc, 
! Nuclei
integer :: numcenters=1
integer :: centershift(100)=0       !! grid point index for each center
real*8 :: nuccharges(100)=2d0       !! nuclear charges

!! BASIS PARAMS

integer :: numpoints=15
real*8 :: spacing=0.25d0

DATATYPE :: harmstrength=0d0   !! harmonic potential

!!    CAP & Exterior Complex Scaling

integer :: capflag=0               !! Number of complex absorbing potentials
integer :: capmode=0               !! Capmode=1 is 
integer :: cappower(100)=2         !!   v_i(r)= capstrength_i*(r/capstart_i)^cappower_i
real*8 :: capstart(100)=0.001d0    !! Capmode=0 is    
real*8 :: capstrength(100)=0.01d0  !!   v_i(r)= capstrength_i*max(0,r-capstart_i)^cappower_i
real*8 :: mincap=0d0 , maxcap=1d30 !! V_CAP = -i* max(mincap,min(maxcap,sum_i v_i))

integer :: scalingflag=0           !! 1 = SMOOTH EXTERIOR COMPLEX SCALING
real*8 :: scalingdistance=10000d0  !! atomic units (bohr)
real*8 :: smoothness=5             !! atomic units (bohr)
real*8 :: scalingtheta=0d0         !! scaling angle
real*8 :: scalingstretch=1d0         !! stretching factor

!! NUMERICAL PARAMS

integer :: eigmode = 0       !!  0 = exact lapack diag 1 = block lanczos

!!    FOR KE MULT IN 1D, toepflag.ne.0 does single toeplitz matvec
integer :: toepflag=0

!!    FOR TOTAL ORBITAL PARALLELIZATION, SET orbparflag=.true., AND parorbsplit=3 in &parinp
logical :: orbparflag=.false.

!!    THE FOLLOWING FLAG IS THEN RELEVANT.  Option for parallel KE matvec, rate limiting step.
integer :: zke_paropt=1   !! 0=sendrecv 1=SUMMA (bcast before) 2=reduce after

!!    fft_batchdim: determines batch size for matrix elements and 
!!    fft_circbatchdim: determines sub batch size for FFT 
!!       defaults set small (less memory, more MPI messages) to avoid MPI problems when doing large 
!!       calculations.  Otherwise bigger values will be faster.  There is a message size sweet spot
!!       on many machines.
integer :: fft_batchdim=1     !! 1 = do nspf matrix elements in nspf batches (less memory)
                              !! 2 = do nspf^2 in one batch (faster unless MPI problems)
integer :: fft_circbatchdim=1 !! 0,1,2, circbatchdim < batchdim; larger faster unless MPI problems
integer :: fft_ct_paropt=1           !! fft_ct_paropt,
                                     !!    like zke_paropt: 0 = sendrecv 1 = summa
integer :: num_skip_orbs=0
integer :: orb_skip(200)=-1

!!  for block lanzos:

integer :: orblanorder=500         !! krylov order for block lanczos orbital calculation
integer :: orblancheckmod=10       !! check every
real*8 :: orblanthresh=1d-4

!! construct virtual orbtials using density based on loaded orbitals
integer :: ivoflag=0
real*8 :: loadedocc(200)=2d0
!! XXSNIPXX 
!! INTERNAL

integer :: orbtargetflag=0
DATATYPE :: orbtarget=DATAZERO

integer :: numspf=-1
integer :: nbox=1  !! BOXES FOR PAR.
integer :: qbox=1

integer :: gridpoints=10

integer, parameter :: numr=1,bornopflag=1

!! internal

integer :: totpoints=(-1)
integer :: nonucrepflag=0
real*8 :: sumcharge, nucrepulsion
integer :: debugflag=0


end module myparams


