
!! PARAMETERS FILE FOR SINC DVR POLYATOMIC BASIS AND HAMILTONIAN:
!!            sincparinp namelist input.

module myparams
implicit none

!! FOR TOTAL ORBITAL PARALLELIZATION, SET THIS TO TRUE, ALSO
!! PARORBSPLIT=3 in &parinp
logical :: orbparflag=.false.

!! THE FOLLOWING FLAG IS THEN RELEVANT.  Option for parallel KE matvec,
!!   rate limiting step in algorithm.
integer :: zke_paropt=1   !! 0=sendrecv 1=SUMMA (bcast before) 2=reduce after

!! toepflag=1, toeplitz T^-1.  toepflag=2, T^-1 and T.   toepflag=2 slower; just do 1.
integer :: toepflag=1

!! fft_mpi_inplaceflag:
!!  0 = out-of-place fft, out-of-place fft inverse
!!      3d FFT + (summa/circ C.T. depending on fft_ctflag)
!!  1 = 3 x (1d FFT , all-to-all index transposition)
integer :: fft_mpi_inplaceflag=1    

!! fft_ct_paropt: if fft_mpi_inplaceflag=0:
!!    like zke_paropt: 0 = sendrecv 1 = summa
integer :: fft_ct_paropt=1

integer :: fft_batchopt=1    !! 1 = do nspf matrix elements in nspf batches (less memory)
                             !! 2 = do nspf^2 matrix elements in one batch (faster)

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
integer :: orblancheckmod=10         !! check every
real*8 :: orblanthresh=1d-4

integer :: capflag=0            !! add complex absorbing potential to one-elec Hamiltonian
real*8  :: capstart=30d0        !! radius in bohr
real*8  :: capstrength=0.001d0  !! atomic units (hartree)
integer :: cappower=2           !! V_CAP = (-i) * capstrength * r^cappower  (r>capstart)

!! XXSNIPXX 
!! INTERNAL

  integer :: myrank=(-1),nprocs=(-1)

  integer :: gridoffset,gridlow,gridhigh,gridsize(100)=(-1)
  logical :: localflag



integer :: notwoflag=0
  integer :: numspf=-1
integer :: nbox(100)=1  !! BOXES FOR PAR.

integer :: griddim=3, ccc=0

integer :: mpifileptr=-1

integer :: gridpoints(100)=10,totpoints=-1,maxgridpoints=-99

integer, parameter :: numr=1,bornopflag=1



!! internal

integer :: coulflag=1
integer :: nonucrepflag=0
real*8 :: sumcharge, nucrepulsion
integer :: debugflag=0


end module myparams


