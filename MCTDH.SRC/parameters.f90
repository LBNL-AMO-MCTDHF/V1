
#include "Definitions.INC"

module fileptrmod
integer :: mpifileptr
integer,parameter :: nullfileptr=798  !! MUST MATCH COOLEYTUKEY_SHARED.F90
end module

module littleparmod
integer, parameter :: MXF=5
end module littleparmod

module r_parameters
integer :: numr
real*8 :: nucrepulsion
DATATYPE, allocatable :: bondpoints(:),bondweights(:)
end module r_parameters

module class_parameters
integer :: numclasses=1
integer, allocatable :: classorb(:,:),nperclass(:), orbclass(:)
end module class_parameters

module dotmod
contains

!! CALLED INSIDE OMP LOOPS

recursive function dot(one,two,n)
  implicit none
  integer,intent(in) :: n
  DATATYPE,intent(in) :: one(n), two(n)
  DATATYPE :: dot, sum
  integer :: i
  sum=0.d0
  do i=1,n
     sum = sum + CONJUGATE(one(i)) * two(i) 
  enddo
  dot=sum
end function dot


!! CALLED INSIDE OMP LOOPS

recursive function hermdot(one,two,n)
  implicit none
  integer,intent(in) :: n
  DATATYPE,intent(in) :: one(n), two(n)
  DATATYPE :: hermdot, sum
  integer :: i
  sum=0.d0
  do i=1,n
     sum =   sum + ALLCON(one(i)) *  two(i) 
  enddo
  hermdot=sum
end function hermdot
end module dotmod

!!YYSNIPYY
!!BB
!! *********************************************************************************************************** !!
!!   Parameters for MCTDHF calculation; parinp NAMELIST input from Input.Inp (default)
!! *********************************************************************************************************** !!
!!
!! Type, variable, default value !! Command-line !! Description 
!!                               !!  option      !! 
!!EE
!!\textbf{\qquad SPARSE - if sparseconfigflag .ne. 0}
!!BB
module lan_parameters
integer :: lanprintflag=0
integer :: lanczosorder=200      !!              !!   lanczos order used in A-vector eigen.
integer :: lancheckstep=20       !!              !! lanczos eigen routine checks for convergence every this # steps
real*8 :: lanthresh=1.d-9        !!              !! convergence criterion.
end module lan_parameters
module sparse_parameters
integer :: sparseprime=1                         !! For reordering config list (experimental)
integer :: sparsesummaflag=0                     !! 0=gather 1=summa(bad) 2=sendrecv(ok)
integer :: sparsedfflag=1                        !! If zero disable separate restricted config walk list
integer :: sparseconfigflag=0    !! Sparse       !! Sparse configuration routines on or off (for large # configs)
integer :: sparseopt =1                          !! 0= direct CI  1= sparse matrix algebra (faster, more memory)
integer :: nonsparsepropmode=1                   !! 0 = ZGCHBV expokit; 1 = mine expmat
integer :: nzflag=1                              !! use only processors with nonzero number of Slaters for a-vec
logical :: shuffle_dfwalktype=.true.             !! redistribute restricted configs evenly over processors
logical :: use_dfwalktype=.false.                !! internal, IGNORE me
end module sparse_parameters
module ham_parameters
!!EE
!!\textbf{\qquad HAMILTONIAN PARAMETERS}
!!BB
integer :: nonuc_checkflag=1     !!              !! Turn off deriv operators in nuclear dofs.
integer :: tdflag=0              !! Pulse        !! Use pulse?
integer :: velflag=0             !!              !!  Length (V(t)) or velocity (A(t))      
!!  Constraintflag=1: Density matrix constraint: assume nothing, keep constant off block diag
!!  2: Dirac-Frenkel (McLachlan/Lagrangian) variational principle.
integer :: constraintflag=0      !! Constraint=  !! As described below (see CONSTRAINT)
integer :: denmatfciflag=0                       !! If .ne. 0 then does denmat constrant as programmed
                                                 !!  before Miyagi's help
DATATYPE :: energyshift=0d0      !!              !! complex shift for making energy real for imperfect CAP/ECS
integer :: drivingflag=0                         !!  Solve for the change in the wave function not wave function 
real*8 :: drivingproportion=0.999999999999d0     !!   -- "psi-prime" treatment.
DATATYPE :: timefac=&            !! Prop/        !! d/dt psi = timefac * H * psi
        DATANEGONE               !!  Relax       !!
real*8 :: mshift=0d0                             !! shift configurations based on m-value.. to break 
                                                 !!  degeneracy for state averaged sym restricted
                                                 !!  (mrestrictmin, mrestrictmax) mcscf; good idea.
integer :: offaxispulseflag=0                    !! internal (not namelist), IGNORE me
end module ham_parameters
module basis_parameters
!!EE
!! {\large \quad For restricted configuration lists (not full CI): SEE MANUAL about dfrestrictflag}
!!BB
integer :: numshells=1           !!              !! number of shells.  greater than one: possibly not full CI. 
!!$ integer :: shelltop(100)=-1  !! Numfrozen=   !! shelltop is namelist input in parinp; the internal variable 
                                                 !!    is allshelltop.   shelltop(1) only may be assigned via 
                                                 !!    command line, with Numfrozen.
integer :: numexcite(100)=999    !! Numexcite=   !! excitations from core shells (i.e. defined for shells 1 
                                                 !!    through numshells-1).  Only numexcite(1) may be 
                                                 !!    assigned via command line input.
integer :: minocc(100)=-999      !!              !! minimum occupation, each shell
integer :: maxocc(100)=999       !!              !!    maximum
integer :: vexcite=999           !!              !! excitations INTO last shell. Use to restrict to doubles, etc.
integer :: shells(100)=1
integer :: allshelltop(0:100)=0           
!!EE
!!{\large \quad CONFIGURATIONS}
!!BB
integer :: mrestrictflag=0       !!              !! If spfrestrictflag=1, restrict wfn to given total M.
integer :: mrestrictval=0        !!              !!    This is the value.
integer :: mrestrictmax= 99999   !!              !! If doing state averaged MCSCF, can include a range of m vals;
integer :: mrestrictmin=-99999   !!              !!    set these variables, with mrestrictflag=0, spfrestrictflag=1
integer :: ugrestrictflag=0      !!              !! like mrestrictflag but for parity
integer :: ugrestrictval=1       !!              !!    like mrestrictval but for parity (1=even,-1=odd)
integer :: restrictflag=1        !!              !! Restrict spin projection of determinants?
!!EE
!!{\large \quad ORBITALS (SINGLE PARTICLE FUNCTIONS, SPFS)}
!!BB
integer :: spfrestrictflag=0     !!              !! Restrict m values of orbitals?  
integer :: spfmvals(1000)=0      !!              !!   M-values of orbitals 
integer :: spfugrestrict=0       !!              !! Restrict parity of orbitals? 
integer :: spfugvals(1000)=0     !!              !!   Parity (+/-1; 0=either) of orbitals (ungerade/gerade)
end module basis_parameters
module output_parameters
integer :: iprintconfiglist=0
end module output_parameters
module timing_parameters
!!EE
!!{\large \quad Timing}
!!BB
integer :: notiming=2            !!NoTiming=0,1,2!! 0=write all 1=write some 2= write none
                                 !!  Timing=2,1,0!!     controls writing of all timing and some info files
integer :: timingout=499         !!              !! various routines output to file (timing info) every this 
                                 !!              !!   # of calls
character(len=200):: timingdir="timing"                       !!
end module timing_parameters
module tol_parameters
real*8 :: lntol=1d-4
real*8 :: invtol=1d-8
end module tol_parameters
module bio_parameters
!!EE
!!{\large \quad Biorthogonalization }
!!BB
integer ::      maxbiodim=100, &
     biodim=10                   !! Krylov dim for biorthogonalization
real*8 ::     biotol=1.d-6
integer :: auto_biortho=1        !! do we want to use biorthonormalization or permutation overlaps? 0 perm overlaps, 1 biortho
integer :: logbranch=1           !! branch of logarithm 3 options 0,1,2
end module bio_parameters
module spfsize_parameters
integer :: parorbsplit=1                         !!  Parallelize orbital calculation.  Might speed up, might
                                                 !!   slow down; check timing.
integer :: spfsize,spfsmallsize                  !! internal, IGNORE
integer :: reducedpotsize = -1                   !! internal, iGNORE
end module spfsize_parameters
module constraint_parameters
!!EE
!!{\large \quad CONSTRAINT: With constraintflag. NEED FOR RESTRICTED CONFIG LIST.}
!!BB
real*8 :: lioreg= 1d-9           !!              !! Regularization for linear solve constraint
integer :: conway=0                              !! for constraintflag=2, dirac frenkel constraint
                                                 !!   0=McLachlan 1=50/50 mix 2=Lagrangian
                                                 !!   3=Lagrangian with epsilon times McLachlan
real*8 :: conprop=1d-1                           !! epsilon for conway=3
real*8 :: condamp=1
integer :: liosize
end module constraint_parameters
module denreg_parameters
real*8 :: denreg=1d-10           !! Denreg=      !! density matrix regularization parameter.
end module denreg_parameters
module df_parameters
!!EE
!!{\large \quad CONFIGURATIONS}
!!BB
integer :: df_restrictflag=0      !!              !! apply constraint to configuration list?  Must use this
                                                 !!  option if constraintflag /= 0.  1 is sufficient;
                                                 !!  dfrestrictflag=2 necessary for action 22. 
                                                 !!  SEE MANUAL FOR PROPER USE OF dfrestrictflag/shell options.
end module df_parameters
module parameters
  use littleparmod;  use fileptrmod;  use r_parameters; use sparse_parameters; use tol_parameters
  use ham_parameters;  use basis_parameters;  use timing_parameters; use spfsize_parameters;
  use df_parameters; use dotmod
  implicit none
!!EE
!!{\large \quad MAIN PARAMETERS }
!!BB
integer :: mcscfnum=1            !! MCSCF=       !! Number of A-vectors (state avgd mcscf or prop)
integer :: numelec=2             !!              !! NUMBER OF ELECTRONS
integer :: orbcompact=1          !!              !! Compact orbitals for expo prop with spfrestrictflag?  Probably ok.
integer :: saveflag=1                            !! if zero does not save wave function at the end
integer :: save_every=0                          !! if nonzero saves wave function every save_every mean field steps
integer :: walkwriteflag=0                       !! Turning OFF writing of walks by default
integer :: spf_flag=1            !!              !! IF ZERO, FREEZE SPFS. (for debugging, or TDCI)
integer :: avector_flag=1        !!              !! IF ZERO, FREEZE AVECTOR. (for debugging)
!! FOR TOTAL ORBITAL PARALLELIZATION with SINC DVR, SET PARORBSPLIT=3
!!   and orbparflag=.true. in &sinc_params.  parorbsplit=3 not supported for atom or diatom.
integer :: par_consplit=0
character (len=200) :: &         !!              !! MAY BE SET BY COMMAND LINE OPTION ONLY: not namelist
  inpfile="Input.Inp        "    !! Inp=filename !!  input.  (=name of input file where namelist input is)
!!EE
!!{\large \quad PROPAGATION/RELAXATION}
!!BB
real*8 :: par_timestep=0.1d0     !! Step=        !! MEAN FIELD TIMESTEP
integer :: improvedrelaxflag=0   !! Relax        !! For improved versus regular relaxtion.   
integer :: threshflag=0          !!              !! Set to 1 for regular relaxation
real*8 :: expotol=1d-8           !!              !! Orbital krylov convergence parameter
integer :: maxexpodim=100        !!              !! Orbital maximum kry dimension OR DGMRES DIM improvedquad=2,3
real*8  :: expostepfac=1.2d0     !!              !! Miscellaneous algorithm parameter
!!EE
!!\textbf{\qquad SPARSE - if sparseconfigflag .ne. 0}
!!BB
integer :: maxaorder=100         !!              !!   lanczos order for sparse a-vector prop and improvedquad=1,3
!!EE
!!{\large \quad PROPAGATION}
!!BB
integer :: littlesteps=1         !!              !! Sub intervals of mean field time step for avector prop
real*8 :: finaltime=4d4          !! T=           !! length of prop.  Overridden for pulse and relax.  
integer :: expodim=10            !!              !! Starting krylov size for orbital propagation (expokit)
!!EE
!!\textbf{\qquad SPARSE - if sparseconfigflag .ne. 0}
!!BB
integer :: aorder=30             !!              !!   
real*8 :: aerror=1d-9            !!              !! lanczos error criterion for sparse a-vector CMF propagation
                                                 !!  within aerror to stop.
!!EE
!!{\large \quad RELAXATION}
!!BB
integer :: improvednatflag=0     !!              !! If improved relax, replace with natorbs every iteration
real*8 :: stopthresh=1d-5        !!              !! Spf error tolerance for relaxation convergence (PRIMARY)
real*8 :: astoptol=1d-7                          !! Avector error tolerance for relax  (BACKUP - WAS STOPTHRESH)
real*8 :: timestepfac=1d0                        !! accelerate relax. multiply par_timestep by this each time
real*8 :: max_timestep=1d10                      !!    maximum time step (limit on exponential growth)
integer :: improvedquadflag=0    !!              !! Use newton iteration not diagonalization for improvedrelax.
                                                 !!     (1 = A-vector, 2 = orbitals, 3 = both)
real*8 :: quadstarttime=-1d0                     !! Waits to turn on orbital quad (2 or 3) until this time
real*8 :: aquadstarttime=-1d0                    !! Waits to turn on avector quad (1 or 3) until this time
real*8 :: maxquadnorm=1d10                       !! brakes to use if improvedquadflag=2 or 3 is diverging
real*8 :: quadtol=1d-1           !!              !! Threshold for solution of Newton solve iterations orbitals.
integer :: quadprecon=1          !!              !! Precondition newton iterations for A-vector?
integer :: quadorthflag=0        !!              !! If eigenfunctions are becoming linearly dependent try this
integer :: normboflag=0          !!              !! Enforce norm at each r value
!!EE
!!{\large \quad ORBITALS (SINGLE PARTICLE FUNCTIONS, SPFS)}
!!BB
integer :: nspf=1                !! Nspf=        !! number of orbitals
integer :: numfrozen=0           !!              !! number of doubly occ orbs (removed from calculation)
!!EE
!!{\large \quad CONFIGURATIONS}
!!BB
integer :: restrict_ms=0         !!              !! For restrictflag=1: 2*m_s: 2x total m_s (multiplicity of 
                                                 !!    lowest included spin states minus one)
integer :: spin_restrictval=0    !!              !! For allspinproject=1: determines spin. Default high spin S=M_s.
                                                 !!  To override use this variable. Equals 2S if S^2 eigval is S(S+1)
integer :: all_spinproject=1     !!              !! Constrain S(S+1) for propagation?
!!EE
!!{\large \quad INITIALIZATION}
!!BB
integer :: loadavectorflag=0     !! A=file       !! load avector to start calculation?
integer :: numavectorfiles=1     !!              !! number of avector files containing a-vectors to load into the
                                                 !! mcscfnum available slots, or if load_avector_product.ne.0,
                                                 !! number of one-wfn files for product, w/ total numelec=numelec
integer :: load_avector_product=0                !! make product wave function for multiple-molecule load
integer :: loadspfflag=0         !! Spf=file     !! load spfs to start calculation?  
                                                 !!    (Otherwise, core eigenfunctions.)
integer :: reinterp_orbflag=0                    !! sinc dvr only, half spacing interpolation for orb load
integer :: spf_gridshift(3,100)=0                !! sinc dvr only, shift orbitals on read, slow index spffile
integer :: numspffiles=1         !!              !! for multiple-molecule (e.g. chemistry) calcs, load many
integer :: numskiporbs=0         !!              !! Reading orbs on file(s), skips members of combined set.
integer :: orbskip(1000)=0       !!              !! Which to skip
character (len=200) :: &         !! A=file       !! A-vector binary file to read.  Can have different configs
     avectorfile(MXF)="Bin/avector.bin"          !!   but should have same number of electrons.   
character (len=200) :: &         !! Spf=file     !! Spf file to read.  Can have fewer m vals, smaller radial 
     spffile(MXF)="Bin/spfs.bin"                 !!   grid, or fewer than nspf total orbitals. 
integer :: avecloadskip(100)=0
integer :: numholes=0                            !! Load a-vector with this many more electrons and annihilate
integer :: numholecombo=1                        !! Number of (products of) annihilation operators to combine 
                                                 !!   (for spin adapt)
integer :: numloadfrozen=0                       !! For loading a vector with orbitals to be frozen (dangerous)
integer, allocatable :: myavectorhole(:,:,:)     !! Namelist input is avectorhole.  Fast index numholes (#
                                                 !!   annihilation operators to multiply together); then    
                                                 !!   numholecombo, number of such products to combine; then
                                                 !!   mcscfnum, wfn of current propagation.
integer :: excitations=0                         !! Similar to holes: number of products of excitation ops
integer :: excitecombos=1                        !!   number of products to linearly combine
integer,allocatable:: myavectorexcitefrom(:,:,:) !! Similar to avectorhole.  Namelist input avectorexcitefrom, etc.
integer,allocatable:: myavectorexciteto(:,:,:)   !!
                                                 !! For both excite and hole: value is spin orbital index
                                                 !! 1=1alpha, 2=1beta, 3=1alpha, etc.
                                                 !! negative input -> negative coefficent
!!EE
!!{\large \quad INPUT / OUTPUT }
!!BB
character(len=200):: finalstatsfile="Dat/finalstats.dat"      !! output for relaxation calculation
real*8 :: pulseft_estep=0.01d0                   !! energy step in hartree for Pulseft.Dat output files
character (len=200) ::      avectoroutfile="Bin/avector.bin"  !! A-vector output file.        
character (len=200) ::      spfoutfile="Bin/spfs.bin"         !! Spf output file.
character(len=200):: psistatsfile="Dat/psistats.dat"          !! for action 25
integer :: psistatfreq=1                                      !!  "
character(len=200):: dendatfile="Dat/denmat.eigs.dat"         !! if notiming=0
character(len=200):: denrotfile="Dat/denmat.rotate.dat"       !!  "
character(len=200):: rdendatfile="Dat/rdenmat.eigs.dat"       !! deprecated
character (len=200) :: ovlspffiles(50)="Bin/ovl.spfs.bin"     !! for actions 20 and 26 
character (len=200) :: ovlavectorfiles(50)="Bin/ovl.avector.bin"!!      (see numovlfiles in ACTIONS)
character(len=200):: outovl="Dat/Overlaps.dat"                !! for action 20
character(len=200):: outmatel="Dat/Matel.dat"                 !! for action 20
character(len=200):: zdipfile="Dat/ZDipoleexpect.Dat"         !! for action 21
character(len=200):: zdftfile="Dat/ZDipoleft.Dat"             !!  "
character(len=200):: ydipfile="Dat/YDipoleexpect.Dat"         !!  "
character(len=200):: ydftfile="Dat/YDipoleft.Dat"             !!  "
character(len=200):: xdipfile="Dat/XDipoleexpect.Dat"         !!  "
character(len=200):: xdftfile="Dat/XDipoleft.Dat"             !!  "
character(len=200):: corrdatfile="Dat/Correlation.Dat"        !! for action 1
character(len=200):: corrftfile="Dat/Corrft.Dat"              !!  "
character(len=200):: fluxmofile="Flux/flux.mo.bin"            !! for actions 15,16,17,23
character(len=200):: fluxafile="Flux/flux.avec.bin"           !!  "
character(len=200):: spifile="Dat/xsec.spi.dat"               !! for action 16 (cross section)
character(len=200):: gtaufile="Dat/gtau.dat"                  !!  " (total flux(t) without e_ke resolution)
character(len=200):: projspifile="Dat/xsec.proj.spi"          !! for action 17 (partial cross section)
character(len=200):: projgtaufile="Dat/gtau.dat"              !!  " (projected flux(t))
character(len=200):: projfluxfile="Flux/proj.flux.wfn.bin"    !!  "
character (len=200):: catspffiles(50)="Bin/cation.spfs.bin"       !!  " (see numcatfiles in ACTIONS)
character (len=200):: catavectorfiles(50)="Bin/cation.avector.bin"!!  "
character(len=200):: fluxafile2="Flux/flux.avec.bin"          !! for action 23
character(len=200):: fluxmofile2="Flux/flux.mo.bin"           !!  "
character(len=200):: natplotbin="Bin/Natlorb.bin"             !! for actions 2,8
character(len=200):: spfplotbin="Bin/Spfplot.bin"             !! for actions 3,9
character(len=200):: denplotbin="Bin/Density.bin"             !! for actions 4,10
character(len=200):: rnatplotbin="Bin/RNatorb.bin"            !! for actions 5,11
character(len=200):: denprojplotbin="Bin/Denproj.bin"         !! for actions 6,12
character(len=200):: natprojplotbin="Bin/Natproj.bin"         !!  "
!!EE
!!{\large \quad PULSE.  (If tdflag=1) }
!!BB
integer :: numpulses=1
integer ::  pulsetype(100)=1      !!              !!  Pulsetype=1:  A(t) = pulsestrength * sin(w t)^2,
real*8  :: omega(100)=1.d0        !!              !!  2:  A(t) = strength * sin(w t)^2 
real*8 :: omega2(100)=1.d0        !!              !!             * sin(w2 t + phaseshift),
real*8 :: pulsestart(100)=0.1d0   !!              !!   
real*8 :: phaseshift(100)=0.d0    !!              !!    pulsestart < t < pulsestart + pi/w; 0 otherwise
real*8 :: chirp(100)=0d0          !!              !!
real*8 :: ramp(100)=0d0
real*8 :: longstep(100)=1d0       !!              !!  Pulsetype 3 available: monochromatic, sinesq start+end
!! NOW COMPLEX
DATATYPE :: pulsestrength(100)=.5d0 !!            !!  A_0 = E_0/omega (strength of field)  
real*8 :: intensity(100)= -1.d0   !!              !! overrides pulse strength.  Intensity, 10^16 W cm^-2 
real*8 :: pulsetheta(100)=0.d0    !!              !!  angle between polarization and bond axis (radians)
real*8 :: pulsephi(100)=0.d0      !!              !!  polarization in xy plane
real*8 :: maxpulsetime=1.d20     !!              !!  
real*8 :: minpulsetime=0.d0      !!              !!  By default calc stops after pulse (overrides finaltime,
                                                 !!   numpropsteps); this will enforce minimum duration
!!EE 
!!{\large \quad ACTIONS} \verb# may also be specified by Act=X where X is an integer on the command line #
!!BB
integer :: numactions=0          !! 
integer :: actions(100)=0        !!              !! ACTIONS
!!   Act=1    Autocorrelation; set corrflag=1 for fourier transform
!!   Act=2    Save natorbs
!!   Act=3    Save spfs
!!   Act=4    Save density
!!   Act=5    Save R-natorbs
!!   Act=6    Save projections of natural configurations  (with Mathematica data in NatCurves/)
!!   Act=7    Save curve data files in LanCurves/ for Mathematica plotting.
!!   Act=8    Enter plotting mode (do not run calculation) and plot natorbs 
!!   Act=9    Enter plotting mode and plot spfs
!!   Act=10   Enter plotting mode and plot density
!!   Act=11   Enter plotting mode and plot natorbs in R
!!   Act=12   Enter plotting mode and plot projections from act=6
!!   Act=13   Nuclear FLUX 
!!   Act=14   Enter plotting mode and analyze nuclear flux
!!   Act=15   Save ELECTRONFLUX
!!   Act=16   Enter plotting mode and analyze ELECTRONFLUX 
!!   Act=17   Enter plotting mode and analyze ELECTRONFLUX (projected)
!!   Act=18   Plot denproj from act=6
!!   Act=19   Enforce natorbs between steps (experimental)
!!   Act=20   Overlaps with supplied eigenfunctions
!!   Act=21   Fourier transform dipole moment with pulse for emission/absorption
!!   Act=22   With Dirac Frenkel restrition, constraintflag=2, check norm of error - NEEDS
!!                 dfrestrictflag=2 not 1
!!   Act=23   Enter plotting/analysis mode and read flux.bin files from Act=15 for overlaps
!!                 between two time dependent wave functions
!!   Act=24   keprojector
!!   Act=25   make psistats.dat
!!   Act=26   mcscf_matel with supplied eigenfunctions 
!!EE
!!{\large \quad ACTION VARIABLES (also see filenames in INPUT/OUTPUT above)}
!!BB
integer :: numovlfiles=1         !!  see ovlspffiles and ovlavectorfiles in INPUT/OUTPUT
integer :: nkeproj=200           !!  For keprojector ACTION 24
real*8 :: keprojminenergy=0.04d0 !!   "
real*8 :: keprojenergystep=0.04d0!!   "
real*8 :: keprojminrad=30        !!   "
real*8 :: keprojmaxrad=40        !!   "
real*8 :: eground=0.d0           !! Eground=     !! energy to shift fourier transform for ACTIONS 1 16 17
                                                 !!   (AUTOCORRELATION AND PHOTOIONIZATION)
complex*16 :: ceground=(0.d0,0d0)!!              !! input as complex-valued instead if you like
real*8 :: autotimestep=1.d0      !! ACTIONS 1 and 21 (autocorrelation and emission/absorption):
                                 !!   time step for fourier transform
!!EE
!!{\large \quad AUTOCORRELATION, PHOTOIONIZATION and EMISSION/ABSORPTION (actions 1,16,17,21)}
!!BB
integer :: ftwindowlength=-99 !! FOR ACTIONS 1,16,17,21 options for damping function of time to be transformed
                              !!    if .ge.0, only damp at end (high frequency cutoff): last ftwindowlength 
                              !!    points in fourier transform are damped by cosine function
integer :: ftwindowpower=1    !! if ftwindowlength not set, use previous (v1.16) windowing function: multiply
                              !!    all points by cos(pi t / 2 / tmax)**ftwindowpower 
integer :: fttriwindow=1      !! If nonzero override other two options, do straight linear damping function
                              !!    NOW DEFAULT linear damping looks best v1.16
integer :: ftdiff=0           !! fourier transform derivative of dipole moment not dipole moment
!!EE
!!{\large \quad EMISSION/ABSORPTION (action 21)}
!!BB
integer :: hanningflag=0         !! for hanning window set nonzero action 1 autocorr
integer :: diptime=100           !! For act=20, outputs copies every diptime atomic units
integer :: dipmodtime=100        !! do ft every autotimestep*dipmodtime
real*8 :: dipolesumstart=1d10,&  !! range for integration of oscillator strength (e.g. for sum rule)
     dipolesumend=1d9            !!    photon energy, atomic units (Hartree), start and end
!!EE
!!{\large \quad PHOTOIONIZATION (actions 15,16,17)}
!!BB
integer :: computeFlux=500, &    ! 0=All in memory other: MBs to allocate
     FluxInterval=50,&           !! Multiple of par_timestep at which to save flux
     FluxSkipMult=1              !! Read every this number of time points.  Step=FluxInterval*FluxSkipMult
integer :: nucfluxopt=0          !! Include imaginary part of hamiltonian from nuc ke 
integer :: FluxOpType=1          !! 0=Full ham 1=halfnium 
integer :: numcatfiles=1         !! see catspffiles and catavectorfiles in INPUT/OUTPUT for action 17
!!$ IMPLEMENT ME (DEPRECATE fluxinterval as namelist input) 
!!$ real*8 :: fluxtimestep=0.1d0
!!EE
!!{\large \quad PLOTTING OPTIONS }
!!BB
integer :: plotmodulus=10        !! PlotModulus= !! For saving nat/spf (Act=2, 6), par_timestep interval
real*8 :: plotpause=0.25d0       !! PlotPause=   !!    for saving natorbs.  plotskip is for stepping over
real*8 :: plotrange=0.2d0        !! PlotZ=       !!    the saved natorbs on read.  others are dimensions
real*8 :: plotcbrange=0.001d0    !!              !!    
real*8 :: plotxyrange=2.d0       !! PlotXY=      !!    
real*8 :: plotview1=70.d0        !!              !! viewing angle, degrees
real*8 :: plotview2=70.d0        !!              !! viewing angle, degrees
integer :: plotnum=10            !! PlotNum=     !! Max number of plots
integer :: plotterm=0            !!              !! 0=x11, 1=aqua
integer :: pm3d=1                !! PM3D         !! Turn pm3d on when plotting
integer :: plotres=50            !!              !! Resolution of plot
integer :: plotskip=1            !! PlotSkip=    !! For plotting (Act=3,5,7), number to skip over
real*8 :: povmult=1d0            !!  Mult df3 data by factor.  For small part of orbs.
integer :: povres=10             !!              !! Povray resolution
integer :: numpovranges=1        !!              !! number of magnifications to plot
real*8 :: povrange(10)=(/ 5,15,& !!              !! Povray plotting ranges (unitless - each magnification)
  80,80,80,80,80,80,80,80 /)
real*8  :: povsparse=1.d-3       !!              !! Sparsity threshold for transformation matrix in povray
!!EE
!!{\large \quad MISC AND EXPERIMENTAL}
!!BB
integer :: conjgpropflag=0       !! for complex Domcke
integer :: pulsewindowtoo=1      !! use window function for E(omega) as well
integer :: nucfluxflag=0         !! 0 = both 1 =electronic 2= nuclear  NOT nuclear flux action 13,14
logical :: readfullvector=.true.
logical :: walksinturn=.false.   !!              !! if you have problems with MPI i/o, maybe try this
integer :: turnbatchsize=5
integer :: nosparseforce=0       !!              !! to override exit with large number of configs, no sparse
integer :: noftflag=0            !!              !! turns off f.t. for flux. use for e.g. core hole propag'n.
integer :: timefacforce=0        !!              !!  override defaults
integer :: timedepexpect=0  !! expectation value of H_0(t) or H(t) reported
integer :: cmf_flag=1            !! CMF/VMF      !! CMF/LMF/QMF or VMF?
integer :: intopt=3              !! RK, GBS      !! SPF/VMF Integrator: 0, RK; 1, GBS, 2, DLSODPK  
                                                 !!  for CMF: 3=expo 4=verlet
integer :: verletnum=80          !!              !! Number of verlet steps per CMF step
integer :: jacprojorth=0         !! 1: projector = sum_i |phi_i> <phi_i|phi_i>^-1 <phi_i|
                                 !! 0:             sum_i |phi_i> <phi_i|    default
integer :: jacsymflag=0          !! 1:  use WP - PW  not (1-P)W   0: default (1-P)W
integer :: jacgmatthird=0        !! 0: default g (constraintflag.ne.0) is linear operator
                                 !! 1: g |phi_c> -> sum_ab |phi_a> g_ab <phi_b|phi_c>
real*8 :: autopermthresh=0.001d0 !! Autoperm=
real*8 :: autonormthresh=0.d0    !! 
integer :: debugflag=0
real*8 :: debugfac=1d0
!!EE
!! XXSNIPXX

integer, parameter :: nodgexpthirdflag=1  !! =1 HARDWIRE 10-2015 not sure about dgexpthird

!integer :: noorthogflag=1        !! 082010 NOW TURNING THIS ON !!
!     hardwire.  eliminated realproject for quad which didn't have the 
!     call to orthog.  don't remember if that was purposeful.
integer :: eigprintflag=0
!integer :: intopt=3              !! RK, GBS      !! SPF/VMF Integrator: 0, RK; 1, GBS, 2, DLSODPK  
!                                                 !!  for CMF: 3=expo 4=verlet
real*8 :: relerr=1.d-10          !!              !! relative error for integrator (RK/GBS/DLSODPK)
                                                 !!   norm error for expo (controls exposteps)
real*8 :: myrelerr=1.d-10        !!              !! absolute error set to norm*myrelerr
                                 !!  or for implicit: uses quadtol (dgmres parameter), myrelerr=orb change,
                                 !!     aerror=a change iteration thresholds; relerr=dgmres call threshold
!! number of MCSCF files
integer :: numfluxfiles=1
integer :: numfluxcurves(20)
!character (len=200) :: fluxfilenames(20)
!!integer :: whichsideproj=0       !! 0 if want to project on N-1 e- state, 1 if on N e- state

integer :: drivingmethod=0   !! 0 = V(t)  1 = H(t)-E
integer :: nonatrotate=0

!!!!!  integer :: whichquad=1           !!     0=my eqn 1=rayleigh quotient  WHICHQUAD 1 HARDWIRE

!!$ integer,parameter :: orderflag=0           !! Order=       !! ordering of configs. 1= first all alphas,  then betas;  
!!$                                                 !!    0= 1a1b2a2b etc.

!! INTERNAL

!!$ IMPLEMENT ME (DEPRECATE fluxinterval as namelist input)   
!!$ integer :: fluxsteps=1

integer :: eachloaded(MXF)=(-99)

integer :: numpropsteps=1000000  !!              !! length of prop.  Overridden for pulse and relax.
integer, parameter :: maxmc=100
DATATYPE :: drivingenergies(maxmc)=0d0
integer :: spfdims(3)=0   !! general - not numerad,lbig+1,2*mbig+1
integer :: spfdimtype(3)=(/0,2,1/)    !! 0 = [0,infty];  1=[-infty,infty];     2=[-A,A] must match on read
integer :: numcurves=1
DATATYPE, allocatable :: elecweights(:,:,:), elecradii(:)
real*8 :: langramthresh=1d-9
integer :: numreduced=1
integer :: headersize=200

integer :: tot_adim

integer :: num_config=-1
integer, allocatable :: configs_perproc(:)
integer :: first_config,last_config,local_nconfig

integer ::  ndof
integer :: totspfdim

integer :: autosize, autosteps
integer :: cdenflag=0            !!              !! Calulate natconfig?  Not necessary if Act=6.
integer :: rdenflag=0            !!              !! Calulate denmat in R? 
real*8 :: globaltime=0.d0   !! for ease of output
!!! real*8 :: rcond=1.d-4      !! singular value for split operator cranck-nicholson
real*8 :: pi
integer :: povplotflag=0
integer, parameter :: natorbfile=765
integer, parameter :: rnatorbfile=769
integer, parameter :: spfplotfile=303
integer, parameter :: denfile=904
integer, parameter :: denprojfile=597
integer, parameter :: natprojfile=488

real*8 :: abserr
integer :: messflag=0
real*8 :: messamount=1.d-2
integer :: skipflag=0

character (len=200) :: nullbuff="                                                                                                                                                                                                        ";;;;;;
end module parameters


module mpi_orbsetmod
  implicit none
  integer :: mpi_orbset_init=0
  integer :: orbsperproc=(-1),  norbsets=(-1),  maxprocsperset=(-1)
  integer :: myorbset=-1, firstmpiorb=-1  !! specific to each processor
  integer, allocatable :: MPI_GROUP_ORB(:),MPI_COMM_ORB(:)
end module

subroutine getOrbSetRange(out_lowspf,out_highspf)
  use parameters
  use mpi_orbsetmod
  implicit none
  integer, intent(out) :: out_lowspf, out_highspf
  if (parorbsplit.ne.1) then
     out_lowspf=1
     out_highspf=nspf
  else
     if (mpi_orbset_init.ne.1) then
        OFLWR "error, mpiorbset init.ne.1",mpi_orbset_init; CFLST
     endif
     out_lowspf=min(firstmpiorb,nspf+1)
     out_highspf=min(nspf,firstmpiorb+orbsperproc-1)
  endif
end subroutine getOrbSetRange



module mpimod
  implicit none

#ifdef MPIFLAG
  include "mpif.h"
#endif
  integer :: MPI_GROUP_WORLD
  integer :: MY_COMM_ORB=(-798)
  integer :: nprocs=1, myrank
  character(len=200), parameter :: mpioutfilebase="MPIOUTS/MPI.Out."
  character(len=200) :: mpioutfile
  integer, parameter :: mpioutfilelen=16
  integer :: stdoutflag=0
  integer :: mpitime=0
  integer :: nonmpitime=0
  integer :: mpiatime=0
  integer :: mpibtime=0

end module mpimod



