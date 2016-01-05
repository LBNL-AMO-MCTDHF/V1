!! PARAMETERS FILE FOR PROLATE SPHEROIDAL DVR BASIS and HAMILTONIAN :
!!            h2parinp NAMELIST input.
!!
module myparams
implicit none

integer :: debugflag=0

!! HAMILTONIAN
real*8 :: nuccharge1=1, &        !!              !! Nuclear charges
     nuccharge2=1                !!              !!
integer :: twoeattractflag=0     !!              !! make 1/r12 attractive
integer :: reducedflag=0         !!              !! Use reduced e^- masses?  Default 1 (yes) with bornopflag=0
integer :: bornopflag=1          !!              !! Born-Op calculation, or with nuclear KE?
                                                 !!   (Sets &parinp nonuc_checkflag=1.)
real*8 :: pro_Hmass=1836.152701d0!!              !! Masses of nuclei: pro_Hmass is the mass of nucleus one and 
!! real*8 :: pro_Dmass=3670.483014d0 
real*8 :: pro_Dmass=1836.152701d0!!              !!   pro_Dmass is the mass of nuc 2.
integer :: JVALUE=0              !!              !! J value for improved adiabatic

!! Additional atoms - hardwired Z=1 (hydrogen) for now
!! old way: hlocs puts hatoms on gridpoints, hlocrealflag=0
!! new way: turn on hlocrealflag; hlocreal is position r,theta (h2 or he)  phi not yet

integer :: numhatoms=0           !!              !! number of h atoms 
integer :: hlocs(3,100)=1        !!              !! dimension (3, numhatoms): first 2 indices xi,eta gridpoint
                                 !!              !! then -1 or 1 for left or right -- all h's in plane for now
integer :: hlocrealflag=0
real*8 :: hlocreal(2,100)=0d0

!! BASIS

integer :: lbig=3,&              !!              !! Number of points in eta minus 1
     mbig=0                      !!              !! Number of m-values exp(imphi)

!! DVR / FEM-DVR for R
!! R
integer :: rnumelements=1        !!              !!  Number of elements in R
integer :: rcelement=100         !!              !!  First element at which ecs begins
real*8 :: rthetaecs=0.d0         !!              !!  ECS scaling angle
integer :: rnumpoints=3          !!              !!  Number of points per element including endpoints in R
real*8 :: relementsize=2d-16     !!              !!  Size of R elements
real*8 :: rstart=1.4d0           !!              !!  First (excluded) gridpoint in R
integer :: capflag=0             !!              !!
real*8 :: capstrength=0.d0       !!              !!  
integer :: cappower=2            !!              !!  

!! XI / r  radial electronic DOF

integer :: xinumpoints=14        !! Numpts=      !!  Number of points per element including endpoints in xi
integer :: xinumelements=2       !! Numel=       !!  Number of elements in xi
real*8 :: xielementsizes(100)=5  !!              !!  Sizes of elements
integer :: xicelement=100        !! Celement=    !!  First element which is scaled
real*8 :: xiecstheta=0.157d0     !!              !!  Scaling angle, radians

!! ORBITAL INITIALIZATION

integer :: ivoflag=0
real*8 :: loadedocc(200)=2d0
integer :: num_skip_orbs=0                       !! For skipping core orbitals for initial diagonalization.
integer :: orb_skip_mvalue(100)= 99              !! Input how many you want to skip, their index (orb_skip)
integer :: orb_skip(20)= -1                      !!   (order in energy, each m value) and m-value
!! XXSNIPXX
!! INTERNAL

 !! dummy of main program variable
  integer, allocatable :: spfmvals(:),spfugvals(:)
  integer :: numspf=-1, spfrestrictflag=0,spfugrestrict=0             


integer :: edim, numerad =  -1,    numr= -1

integer, parameter :: atomflag=0
integer :: bandwidth = -1
integer :: mpifileptr = -1
integer :: mseriesmax = -1
integer :: lseriesmax = -1
integer :: jacobisummax = -1


logical :: temppot=.false.
integer :: lobattoopt=0          !!              !! turns on G-L for first element, not gauss-radau


real*8 :: capstart   !! set to be first ecs gridpoint (uses xicelement)
integer :: xigridpoints
real*8 :: Rmass 
real*8 :: totalmass
real*8 :: mass_asym
integer :: rgridpoints
integer :: numeta

end module myparams

