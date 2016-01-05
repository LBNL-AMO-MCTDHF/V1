
!! PARAMETERS FILE FOR ATOM DVR INPUT:  heparinp NAMELIST input.

module myparams
implicit none

integer :: debugflag=0

!! HAMILTONIAN AND BASIS

integer :: hecelement=100,&      !! Celement=    !! Frist element that is scaled
     henumelements=2,&           !! Numel=       !! Number of elements in r
     henumpoints=14              !! Numpts=      !! Number of points per elements 
real*8 :: heelementsizes(100)=5,&!!              !! Sizes of elements in r
     heecstheta=0.157d0          !!              !! ECS scaling angle
integer :: lbig=3,mbig=0
real*8 :: nuccharge1=1d0

!! ORBITAL INITIALIZATION

integer :: ivoflag=0
real*8 :: loadedocc(200)=2d0
integer :: num_skip_orbs=0                       !! For skipping core orbitals for initial diagonalization.
integer :: orb_skip_mvalue(100)= 99              !! Input how many you want to skip, their index (orb_skip)
integer :: orb_skip(20)= -1                      !!   (order in energy, each m value) and m-value

!! XXSNIPXX
!! MISC AND INTERNAL

 !! dummy of main program variable
  integer, allocatable :: spfmvals(:),spfugvals(:)
  integer :: numspf=-1, spfrestrictflag=0,spfugrestrict=0             


integer :: temp_glflag=1


integer :: numhatoms=0        !!  number of h atoms (only for prolate) ( variables not used here )
                             !!! old way: hlocs puts hatoms on gridpoints
integer :: hlocs(3,100)=1    !! dimension 3, numhatoms: first index
                             !!  xi, eta gridpoints; -1 or 1 for left or right
                             !!  i.e. all h's in plane for now
        !!! new way: turn on hlocrealflag; hlocreal is position r,theta (h2 or he)  phi not yet
integer :: hlocrealflag=0
real*8 :: hlocreal(2,100)=0d0

integer :: bandwidth = -1,numerad = -1, edim, mpifileptr = -1, mseriesmax = -1,hegridpoints
integer :: lseriesmax = -1,jacobisummax = -1
integer, parameter :: atomflag=1,numr=1, bornopflag=1

end module myparams

