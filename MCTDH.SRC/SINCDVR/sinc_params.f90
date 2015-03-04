
!! PARAMETERS FILE FOR SINC DVR POLYATOMIC BASIS AND HAMILTONIAN:
!!            sincparinp namelist input.

module myparams
implicit none

!! FOR TOTAL ORBITAL PARALLELIZATION, SET THIS TO TRUE, ALSO
!! PARORBSPLIT=3 in &parinp
logical :: orbparflag=.false.

!! toepflag=1, toeplitz T^-1.  toepflag=2, T^-1 and T.   toepflag=2 slower; just do 1.
integer :: toepflag=1

!! T^-1 must be either read from file or computed with JRJ routine.  Still ironing it out.
!!  runtoothflag=1  computes;  runtoothflag=0 reads in which case you need the invke2file

integer :: runtoothflag=0

integer :: num_skip_orbs=0
integer :: orb_skip(200)=-1

character(len=80) :: invke2file="Invke2.Bin"
integer :: toothnbig=57
integer :: toothnsmall=19
integer :: toothngrid=19

integer :: numcenters=1
integer :: centershift(3,100)=0    !! grid point index for each center
real*8 :: nuccharges(100)=2d0

integer :: numpoints(100)=(-1)
real*8 :: spacing=0.25d0

integer :: orblanorder=500         !! krylov order for block lanczos orbital calculation
integer :: orblancheckmod=10         !! check every
real*8 :: orblanthresh=1d-4

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


