!! PARAMETERS FILE FOR SINC DVR POLYATOMIC BASIS AND HAMILTONIAN:
!!            sincparinp namelist input.

module myparams
implicit none

!! FOR TOTAL ORBITAL PARALLELIZATION, SET orbparflag=.true., AND parorbsplit=3 in &parinp
logical :: orbparflag=.false.
!! integer :: orbparlevel=3  !!  in namelist, but not in myparams

!! THE FOLLOWING FLAG IS THEN RELEVANT.  Option for parallel KE matvec, rate limiting step.
integer :: zke_paropt=1   !! 0=sendrecv 1=SUMMA (bcast before) 2=reduce after

!! fft_batchdim: determines batch size for matrix elements and 
!! fft_circbatchdim: determines sub batch size for FFT 
!!    defaults set small (less memory, more MPI messages) to avoid MPI problems when doing large 
!!    calculations.  Otherwise bigger values will be faster.  There is a message size sweet spot
!!    on many machines.
integer :: fft_batchdim=1     !! 1 = do nspf matrix elements in nspf batches (less memory)
                              !! 2 = do nspf^2 in one batch (faster unless MPI problems)
integer :: fft_circbatchdim=1 !! 0,1,2, circbatchdim < batchdim; larger faster unless MPI problems
integer :: fft_mpi_inplaceflag=1     !! fft_mpi_inplaceflag:
                                     !!  0 = out-of-place fft, out-of-place fft inverse
                                     !!      3d FFT + (summa/circ C.T. depending on fft_ctflag)
                                     !!  1 = 3 x (1d FFT , all-to-all index transposition)
integer :: fft_ct_paropt=1           !! fft_ct_paropt, relevant if fft_mpi_inplaceflag=0
                                     !!    like zke_paropt: 0 = sendrecv 1 = summa
integer :: num_skip_orbs=0
integer :: orb_skip(200)=-1

integer :: toothnsmall=40
integer :: toothnbig=240

integer :: numcenters=1
integer :: centershift(3,100)=0    !! grid point index for each center
real*8 :: nuccharges(100)=2d0

integer :: numpoints(100)=15
real*8 :: spacing=0.25d0

integer :: orblanorder=500         !! krylov order for block lanczos orbital calculation
integer :: orblancheckmod=10       !! check every
real*8 :: orblanthresh=1d-4

integer :: capflag=0               !! Number of complex absorbing potentials
integer :: capmode=0               !! Capmode=1 is 
integer :: cappower(100)=2         !!   v_i(r)= capstrength_i*(r/capstart_i)^cappower_i
real*8 :: capstart(100)=0.001d0    !! Capmode=0 is    
real*8 :: capstrength(100)=0.01d0  !!   v_i(r)= capstrength_i*max(0,r-capstart_i)^cappower_i
real*8 :: mincap=0d0 , maxcap=1d30 !! V_CAP = -i* max(mincap,min(maxcap,sum_i v_i))

integer :: maskflag=0
integer :: masknumpoints=0

integer :: scalingflag=0           !! 1 = SMOOTH EXTERIOR COMPLEX SCALING
real*8 :: scalingdistance=10000d0  !! atomic units (bohr)
real*8 :: smoothness=5             !! atomic units (bohr)
real*8 :: scalingtheta=0d0         !! scaling angle
integer :: scalingorder=2          !! should be 2 or greater!
 real*8 :: tinv_tol=1d-3
!! XXSNIPXX 
!! INTERNAL

!!$ integer :: scalingorders(3)=0      !!   put nodes where nuclei are, that's the idea
!!$ real*8 :: scalingterms(100,3)=0d0  !!   doesn't work well yet
!!$ integer :: scalingdflag=0
!!$ real*8 :: scalingdconst(3)=1d0
!!$ real*8 :: tinv_tol=1d-2

!!real*8 :: ecstheta=0d0             !! X(x)= x + e^(i ecstheta)*sum_j scalingterms(j,1) x^(j-1)



integer :: notwoflag=0
  integer :: numspf=-1
integer :: nbox(100)=1  !! BOXES FOR PAR.

integer :: griddim=3, ccc=0


integer :: gridpoints(100)=10,totpoints=-1,maxgridpoints=-99

integer, parameter :: numr=1,bornopflag=1



!! internal

integer :: nonucrepflag=0
real*8 :: sumcharge, nucrepulsion
integer :: debugflag=0


end module myparams


