
#include "Definitions.INC"

module constant_parameters
real*8 :: pi
end module constant_parameters

module fileptrmod
integer :: mpifileptr
integer,parameter :: nullfileptr=798  !! MUST MATCH COOLEYTUKEY_SHARED.F90
end module

module r_parameters
integer :: numr
real*8 :: nucrepulsion
DATATYPE, allocatable :: bondpoints(:),bondweights(:)
end module r_parameters

module class_parameters
integer :: numclasses=1
integer, allocatable :: classorb(:,:),nperclass(:), orbclass(:)
end module class_parameters

module output_parameters
integer :: iprintconfiglist=0
end module output_parameters

module clockmod
contains
  subroutine myclock(mytime)
    integer :: values(10),mytime
    integer, parameter :: fac(5:8)=(/60*60*1000,60*1000,1000,1/)  !! hour,minute,second,millisecond
    call date_and_time(values=values)
    mytime=values(8)+values(7)*fac(7)+values(6)*fac(6)+values(5)*fac(5)
  end subroutine myclock
end module clockmod

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
!!   Input variables for MCTDHF: can be specified in namelist &parinp or namelist &pulse in input file.
!! *********************************************************************************************************** !!
!!
!! Type, variable, default value !! Command-line !! Description 
!!                               !!  option      !! 
!!EE
!!{\large \quad HAMILTONIAN PARAMETERS}
!!BB
module ham_parameters
integer :: nonuc_checkflag=1     !!              !! Turn off deriv operators in nuclear dofs.
integer :: tdflag=0              !! Pulse        !! Use pulse?
integer :: velflag=0             !!              !!  Length (V(t)) or velocity (A(t))NAMELIST PULSE
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
!!EE
!!{\large \quad PULSE.  (If tdflag=1) }
!!BB
module pulse_parameters           !!      NAMELIST PULSE
integer :: numpulses=1            !!  number of pulses, enter pulsetype, omega, etc. for each
integer :: reference_pulses=0     !!  for domcke method action 21,29. Enter numpulses+reference_pulses pulses, then
                                  !!    output (eg absorption) is calc'd with e-field for pulses numpules+1 and up
integer ::  pulsetype(100)=1      !!              !!  Pulsetype=1:  A(t) = pulsestrength * sin(w t)^2,
real*8  :: omega(100)=1.d0        !!              !!  2:  A(t) = strength * sin(w t)^2 
real*8 :: omega2(100)=1.d0        !!              !!             * sin(w2 t + phaseshift),
real*8 :: pulsestart(100)=0d0     !!              !!   
real*8 :: phaseshift(100)=0d0     !!              !!    pulsestart < t < pulsestart + pi/w; 0 otherwise
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
end module pulse_parameters
!!EE
!!{\large \quad CONFIGURATIONS / SLATER DETERMINANTS}
!!BB
module df_parameters
integer :: df_restrictflag=0      !!              !! apply constraint to configuration list?  Must use this
                                                 !!  option if constraintflag /= 0.  1 is sufficient;
                                                 !!  dfrestrictflag=2 necessary for action 22. 
                                                 !!  SEE MANUAL FOR PROPER USE OF dfrestrictflag/shell options.
end module df_parameters
module basis_parameters
integer :: holeflag=0            !! Holes        !! index holes not electrons
integer :: mrestrictflag=0       !!              !! If spfrestrictflag=1, restrict wfn to given total M.
integer :: mrestrictval=0        !!              !!    This is the value.
integer :: mrestrictmax= 99999   !!              !! If doing state averaged MCSCF, can include a range of m vals;
integer :: mrestrictmin=-99999   !!              !!    set these variables, with mrestrictflag=0, spfrestrictflag=1
integer :: ugrestrictflag=0      !!              !! like mrestrictflag but for parity
integer :: ugrestrictval=1       !!              !!    like mrestrictval but for parity (1=even,-1=odd)
integer :: restrictflag=1        !!              !! Restrict spin projection of determinants?
!!EE
!! \textbf{\qquad For restricted configuration lists (not full CI) (see manual)}
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
integer :: shells(100)=1
integer :: allshelltop(0:100)=0           
!!EE
!!{\large \quad ORBITALS (SINGLE PARTICLE FUNCTIONS, SPFS)}
!!BB
integer :: spfrestrictflag=0     !!              !! Restrict m values of orbitals?  
integer :: spfmvals(1000)=0      !!              !!   M-values of orbitals 
integer :: spfugrestrict=0       !!              !! Restrict parity of orbitals? 
integer :: spfugvals(1000)=0     !!              !!   Parity (+/-1; 0=either) of orbitals (ungerade/gerade)
integer :: tot_adim              !! = numr*localnconfig     INTERNAL, IGNORE ME
integer :: num_config=-1                         !!         INTERNAL, IGNORE ME
integer, allocatable :: configs_perproc(:)       !!         INTERNAL, IGNORE ME
integer :: first_config,last_config,local_nconfig!!         INTERNAL, IGNORE ME
integer :: numpart,num2part    !! particles or holes        INTERNAL, IGNORE ME
end module basis_parameters
!!EE
!!{\large \quad for SPARSE CONFIGURATION ROUTINES - if sparseconfigflag .ne. 0}
!!BB
module lanparameters
integer :: lanprintflag=0
integer :: lanczosorder=300      !!              !! lanczos order used in A-vector eigen routine.
integer :: lancheckstep=20       !!              !! lanczos eigen routine checks for convergence every this # steps
real*8 :: lanthresh=1.d-9        !!              !! convergence criterion.
end module lanparameters
module sparse_parameters
integer :: sparseconfigflag=0    !! Sparse       !! Sparse configuration routines on or off (for large # configs)
integer :: sparseopt =1                          !! 0= direct CI  1= sparse matrix algebra (faster, more memory)
integer :: par_consplit=1                        !! Fully distribute a-vector
integer :: sparsesummaflag=0                     !! a-vector matvec mode 0=gather 1=summa(bad) 2=sendrecv(ok)
integer :: sparsedfflag=1                        !! If zero disable separate restricted config walk list
integer :: nonsparsepropmode=1                   !! 0 = ZGCHBV expokit; 1 = mine expmat
integer :: nzflag=1                              !! use only processors with nonzero number of Slaters for a-vec
logical :: shuffle_dfwalktype=.true.             !! redistribute restricted configs evenly over processors
integer :: sparseprime=1                         !! For reordering config list (experimental)
logical :: use_dfwalktype=.false.                !! internal, IGNORE me
logical :: use_fockmatrix=.false.                !! internal, IGNORE me
end module sparse_parameters
!!EE
!!{\large \quad Orbital parallelization}
!!BB
module spfsize_parameters
integer :: parorbsplit=1                !! 1=keep all orbs, prop in parallel   3=polyatomic, full MPI
integer :: spfsize,spfsmallsize         !! size of orbitals each processor                internal, IGNORE
integer :: reducedpotsize = -1          !! size of array used to store reduced potential  internal, IGNORE
integer :: totspfdim                    !! nspf * spfsize (local size only)               internal, IGNORE
integer :: spfdims(3)=0                 !! general - not numerad,lbig+1,2*mbig+1          internal, IGNORE
integer :: spfdimtype(3)=(/0,2,1/)      !! 0 = [0,infty];  1=[-infty,infty]; 2=[-A,A]     internal, IGNORE
end module spfsize_parameters
!!EE
!!{\large \quad CONSTRAINT: With constraintflag. NEED FOR RESTRICTED CONFIG LIST.}
!!BB
module constraint_parameters
real*8 :: lioreg= 1d-9           !!              !! Regularization for linear solve constraint
integer :: conway=0                              !! for constraintflag=2, dirac frenkel constraint
                                                 !!   0=McLachlan 1=50/50 mix 2=Lagrangian
                                                 !!   3=Lagrangian with epsilon times McLachlan
real*8 :: conprop=1d-1                           !! epsilon for conway=3
real*8 :: condamp=1
end module constraint_parameters
!!EE
!!{\large \quad Timing}
!!BB
module timing_parameters
integer :: notiming=2            !!NoTiming=0,1,2!! 0=write all 1=write some 2= write none
                                 !!  Timing=2,1,0!!     controls writing of all timing and some info files
integer :: timingout=499         !!              !! various routines output to file (timing info) every this 
                                 !!              !!   # of calls
character(len=SLN):: timingdir="timing"                       !!
end module timing_parameters
!!EE
!!{\large \quad Biorthogonalization }
!!BB
module bio_parameters
integer ::      maxbiodim=100, & !! Max krylov dim for biorthogonalization
     biodim=10                   !! Starting krylov dim for biorthogonalization
real*8 ::     biotol=1.d-7       !! Expokit tolerance parameter for biorthogonalization
integer :: logbranch=1           !! branch of logarithm 4 options 0,1,2,3
integer :: auto_biortho=1        !! internal, IGNORE me
end module bio_parameters
!!EE
!!{\large \quad Regularization}
!!BB
module tol_parameters
real*8 :: lntol=1d-4             !! Regularization of natural logarithm
real*8 :: invtol=1d-8            !! Regularization of inverses especially for spectral expansions
end module tol_parameters
module denreg_parameters
real*8 :: denreg=1d-10           !! Denreg=      !! density matrix regularization parameter.
end module denreg_parameters
module parameters
  use fileptrmod;  use r_parameters; use sparse_parameters; use tol_parameters
  use ham_parameters;  use basis_parameters;  use timing_parameters; use spfsize_parameters;
  use df_parameters; use dotmod;    use constant_parameters; use clockmod
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
character (len=SLN) :: &         !!              !! MAY BE SET BY COMMAND LINE OPTION ONLY: not namelist
  inpfile="Input.Inp        "    !! Inp=filename !!  input.  (=name of input file where namelist input is)
!!EE
!!{\large \quad PROPAGATION/RELAXATION}
!!BB
real*8 :: par_timestep=0.1d0     !! Step=        !! MEAN FIELD TIMESTEP
integer :: improvedrelaxflag=0   !! Relax        !! For improved versus regular relaxtion.   
integer :: threshflag=0          !!              !! Set to 1 for regular relaxation
real*8 :: expotol=1d-8           !!              !! Orbital krylov convergence parameter
integer :: expodim=10            !!              !! Starting krylov dimension for expokit orbital propagation
integer :: maxexpodim=100        !!              !! Maximum krylov dimension for expokit orbital propagation
integer :: scalarflag=0          !!              !! Flag for inclusion of scalar terms in reduced hamiltonian
real*8  :: expostepfac=1.2d0     !!              !! Miscellaneous algorithm parameter
integer :: exact_exchange=0      !!              !! Exact exchange available for debug
integer :: exchange_mode=0       !!              !! For default inexact: 0= previous 1,2=new
!!EE
!!\textbf{\qquad SPARSE - if sparseconfigflag .ne. 0}
!!BB
real*8 :: aerror=1d-9            !!              !! error criterion for sparse a-vector prop
                                                 !! and eigenfunctions with improvedquad=1,3
integer :: maxaorder=100         !!              !! max lanczos order for sparse a-vector expokit prop or
                                                 !!    dgmres solve for avector improvedquadflag=1,3
!!EE
!!{\large \quad PROPAGATION}
!!BB
integer :: step_flag=2           !!              !! Main propagation steps, CMF or LMF depending on pre/propflag
integer :: prepropflag=1         !!              !! Option for pre-propagation step
integer :: postpropflag=2        !!              !! Option for LMF post-propagation 1=orbitals 2=A-vector
integer :: littlesteps=1         !!              !! Sub intervals of mean field time step for avector prop
real*8 :: finaltime=4d4          !! T=           !! length of prop.  Overridden for pulse and relax.  
!!EE
!!\textbf{\qquad SPARSE - if sparseconfigflag .ne. 0}
!!BB
integer :: aorder=30             !!              !! Starting lanczos order for sparse a-vector expokit prop
!!EE
!!{\large \quad RELAXATION - if improvedrelaxflag.ne.0}
!!BB
integer :: improvednatflag=0     !!              !! Replace with natorbs every iteration
integer :: improvedfockflag=0    !!              !! Alternatively, replace with generalized fock eigfuncts
integer :: improvedquadflag=0    !!              !! Use newton iteration not diagonalization for improvedrelax.
                                                 !!     (1 = A-vector, 2 = orbitals, 3 = both)
integer :: followflag=0                          !! use overlaps to follow vectors
real*8 :: stopthresh=1d-5        !!              !! Spf error tolerance for relaxation convergence
real*8 :: astoptol=1d-7                          !! Avector error tolerance for relax
real*8 :: timestepfac=1d0                        !! accelerate relax. multiply par_timestep by this each time
real*8 :: max_timestep=1d10                      !!    maximum time step (limit on exponential growth)
integer :: jacsymquad=1                          !! Turn on jacsymflag for orbital quad 2,3
real*8 :: quadstarttime=-1d0                     !! Waits to turn on orbital quad (2 or 3) until this time
real*8 :: aquadstarttime=-1d0                    !! Waits to turn on avector quad (1 or 3) until this time
real*8 :: maxquadnorm=1d10                       !! brakes to use if improvedquadflag=2 or 3 is diverging
real*8 :: quadtol=0.316d0        !!              !! Threshold for solution of Newton solve iterations orbitals.
integer :: quadprecon=1          !!              !! Precondition newton iterations for A-vector?
integer :: quadorthflag=0        !!              !! If eigenfunctions are becoming linearly dependent try this
integer :: normboflag=0          !!              !! Enforce norm at each r value
!!EE
!!\textbf{\qquad SPARSE - for relaxation if sparseconfigflag .ne. 0 see module lanparameters above, and}
!!BB
integer :: maxdgdim = 800                        !! Max dimension for improvedquadflag options
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
character (len=SLN) :: &         !! A=file       !! A-vector binary file to read.  Can have different configs
     avectorfile(MXF)="Bin/avector.bin"          !!   but should have same number of electrons.   
character (len=SLN) :: &         !! Spf=file     !! Spf file to read.  Can have fewer m vals, smaller radial 
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
integer :: messflag=0                            !! add static to initial orbitals
real*8 :: messamount=1.d-4                       !!    how much
integer :: messavec=0                            !! add static to initial a-vector
real*8 :: messaamount=1.d-4                      !!    how much
!!EE
!!{\large \quad INPUT / OUTPUT }
!!BB
character(len=SLN):: finalstatsfile="Dat/finalstats.dat"         !! output for relaxation calculation
real*8 :: pulseft_estep=0.01d0                                   !! step in hartree for Pulseft.Dat output files
character (len=SLN) ::      avectoroutfile="Bin/avector.bin"     !! A-vector output file.        
character (len=SLN) ::      spfoutfile="Bin/spfs.bin"            !! Spf output file.
character(len=SLN):: psistatsfile="Dat/psistats.dat"             !! for action 25
integer :: psistatfreq=1                                         !!  "
character(len=SLN):: dendatfile="Dat/denmat.eigs.dat"            !! if notiming=0
character(len=SLN):: denrotfile="Dat/denmat.rotate.dat"          !!  "
character(len=SLN):: rdendatfile="Dat/rdenmat.eigs.dat"          !! deprecated
character (len=SLN) :: ovlspffiles(50)="Bin/ovl.spfs.bin"        !! for actions 20 and 26 
character (len=SLN) :: ovlavectorfiles(50)="Bin/ovl.avector.bin" !!      (see numovlfiles in ACTIONS)
character(len=SLN):: outovl="Dat/Overlaps.dat"                   !! for action 20
character(len=SLN):: outmatel="Dat/Matel.dat"                    !! for action 20
integer :: act21circ=0                                           !! set nonzero to enable circular polarization 
                                                                 !!     output for action 21,29
character(len=SLN):: xdipfile="Dat/XDipoleexpect.Dat",&          !! for action 21,29 dipole moments D(t)
 ydipfile="Dat/YDipoleexpect.Dat",&
 zdipfile="Dat/ZDipoleexpect.Dat"                 
character(len=SLN):: xtworkfile="Dat/XTWork.Dat",&               !! for action 21,29 work per pulse integral dt
 ytworkfile="Dat/YTWork.Dat",ztworkfile="Dat/ZTWork.Dat",&       !!  "
 xytworkfile="Dat/XYTWork.Dat",yztworkfile="Dat/YZTWork.Dat",&   !!  "  circ polarization: three torque components
 zxtworkfile="Dat/ZXTWork.Dat"                                   !!  "  XY, YZ, ZX.
character(len=SLN):: xdftfile="Dat/XDipoleft.Dat",&              !! for action 21,29 dipole moments D(omega)
 ydftfile="Dat/YDipoleft.Dat", zdftfile="Dat/ZDipoleft.Dat",&    !! and emission/absorption
 xydftfile="Dat/XYDipoleft.Dat",xzdftfile="Dat/XZDipoleft.Dat",& !!  "  circ polarization: six components
 yxdftfile="Dat/YXDipoleft.Dat",yzdftfile="Dat/YZDipoleft.Dat",& !!  "  XY, XZ, YX, YZ, ZX, ZY
 zxdftfile="Dat/ZXDipoleft.Dat",zydftfile="Dat/ZYDipoleft.Dat"   !!  "  e.g. XY: X->Y, Y->(-X)
character(len=SLN):: xoworkfile="Dat/XOWork.Dat",&               !! for action 21,29 work per pulse integral domega
 yoworkfile="Dat/YOWork.Dat",zoworkfile="Dat/ZOWork.Dat",&       !!  "
 xyoworkfile="Dat/XYOWork.Dat",xzoworkfile="Dat/XZOWork.Dat",&   !!  "  circ polarization, six components
 yxoworkfile="Dat/YXOWork.Dat",yzoworkfile="Dat/YZOWork.Dat",&   !!  "
 zxoworkfile="Dat/ZXOWork.Dat",zyoworkfile="Dat/ZYOWork.Dat"     !!  "
character(len=SLN):: xophotonfile="Dat/XOPhoton.Dat",&           !! for action 21,29 work per pulse integral domega
 yophotonfile="Dat/YOPhoton.Dat",zophotonfile="Dat/ZOPhoton.Dat",&     !!  "
 xyophotonfile="Dat/XYOPhoton.Dat",xzophotonfile="Dat/XZOPhoton.Dat",& !!  "  circ polarization, six components
 yxophotonfile="Dat/YXOPhoton.Dat",yzophotonfile="Dat/YZOPhoton.Dat",& !!  "
 zxophotonfile="Dat/ZXOPhoton.Dat",zyophotonfile="Dat/ZYOPhoton.Dat"   !!  "
character(len=SLN):: corrdatfile="Dat/Correlation.Dat"        !! for action 1
character(len=SLN):: corrftfile="Dat/Corrft.Dat"              !!  "
character(len=SLN):: fluxmofile="Flux/flux.mo.bin"            !! for actions 15,16,17,23
character(len=SLN):: fluxafile="Flux/flux.avec.bin"           !!  "
character(len=SLN):: spifile="Dat/xsec.spi.dat"               !! for action 16 (cross section)
character(len=SLN):: fluxtsumfile="Dat/fluxtsum.dat"          !!  " total flux integral dt action 28 too
character(len=SLN):: projspifile="Dat/xsec.proj.spi"          !! for action 17 (partial cross section)
character(len=SLN):: projfluxtsumfile="Dat/projfluxtsum"      !!  " total projected flux integral dt action 28 too
character(len=SLN):: projfluxfile="Flux/proj.flux.wfn.bin"    !!  " main output action 17
character (len=SLN):: catspffiles(50)="Bin/cation.spfs.bin"   !!  " (see numcatfiles in ACTIONS)
character (len=SLN):: catavectorfiles(50)="Bin/cation.avector.bin"!!  "
character (len=SLN):: strongcatspffiles(50)=&                 !! for action 17 if strongcatflag.ne.0, use back-
     "Flux/cation.mo.bin"                                     !!  propagated cation states stored with action 15
character (len=SLN):: strongcatavectorfiles(50)=&             !!  "
     "Flux/cation.avec.bin"                                   !!  "
character(len=SLN):: angprojspifile="Dat/xsec.angproj.spi"    !!  " if angularflag.ne.0, output file
character(len=SLN):: angprojfluxtsumfile="Dat/angprojfluxtsum"!!  " and angular integral dt action 28 too
character(len=SLN):: fluxafile2="Flux/flux.avec.bin"          !! for action 23
character(len=SLN):: fluxmofile2="Flux/flux.mo.bin"           !!  "
character(len=SLN):: natplotbin="Bin/Natlorb.bin"             !! for actions 2,8
character(len=SLN):: spfplotbin="Bin/Spfplot.bin"             !! for actions 3,9
character(len=SLN):: denplotbin="Bin/Density.bin"             !! for actions 4,10
character(len=SLN):: rnatplotbin="Bin/RNatorb.bin"            !! for actions 5,11
character(len=SLN):: denprojplotbin="Bin/Denproj.bin"         !! for actions 6,12
character(len=SLN):: natprojplotbin="Bin/Natproj.bin"         !!  "
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
!!   Act=27   like action 16 total photoionization during calculation just integral dt
!!   Act=28   like action 17 partial photoionization during calculation just integral dt
!!   Act=29   redo action 21 with file saved from flux.bin files from Act=15
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
real*8 :: autostart=0d0          !! time zero for action 1. negative setting -> end of pulse.
real*8 :: autotimestep=1.d0      !! ACTIONS 1,21,29 (autocorrelation and emission/absorption):
                                 !!   time step for fourier transform
integer :: auto_subtract=0       !! if nonzero subtract tent function for autocorrelation (action 1)
!!EE
!!{\large \quad FTs for AUTOCORRELATION, PHOTOIONIZATION and EMISSION/ABSORPTION (actions 1,16,17,21,29)}
!!BB
!! now defined for each action.  Variables in actions.f90.  Defaults in getparams.f90.
!!
!! int fttriwindow(MAXACTIONS)=1    !! If nonzero, window function is ((tmax-t)/tmax)**ftwindowpower
!! int ftwindowpower(MAXACTIONS)=1  !! If fttriwindow=0, 
!!                                  !!    window function is cos(pi t / 2 / tmax)**ftwindowpower 
integer :: ftdiff=1                 !! fourier transform derivative of dipole moment not dipole moment
!!EE
!!{\large \quad EMISSION/ABSORPTION (action 21,29)}
!!BB
integer :: hanningflag=0         !! for hanning window set nonzero action 1 autocorr
integer :: diptime=100           !! For act=20, outputs copies every diptime atomic units
integer :: dipmodtime=100        !! do ft every autotimestep*dipmodtime
real*8 :: dipolesumstart=0d0,&   !! range for integration of oscillator strength (e.g. for sum rule), photon
     dipolesumend=0d0            !!    energy atomic units (Hartree) start and end, dipolesumend default set
                                 !!    based on &pulse namelist input (depending on omegas)
logical :: redobra=.false.       !! for complex Domcke action 29
!!EE
!!{\large \quad PHOTOIONIZATION (actions 15,16,17,27,28)}
!!BB
integer :: computeFlux=500, &    ! 0=All in memory other: MBs to allocate
     FluxInterval=50,&           !! Multiple of par_timestep at which to save flux
     FluxSkipMult=1,&            !! Read every this number of time points.  Step=FluxInterval*FluxSkipMult
     flux_subtract=1             !! if nonzero subtract t=0 flux matrix element for decaying initial state
integer :: gaugefluxflag=0       !! Transform length to velocity (1) or vice versa (2) for 16,17,27,28
integer :: nucfluxopt=1          !! NucFlux=     !! Include imag part of hamiltonian from nuc ke 2=only that
integer :: FluxOpType=1          !! 0=Full ham 1=halfnium 
integer :: numcatfiles=1         !! see catspffiles and catavectorfiles in INPUT/OUTPUT for action 17,28
integer :: strongcatflag=0       !! for strong field ionization during the pulse, use back-propagated cation
DATATYPE :: catenergies(50)=0d0  !! need to input cation energies for strongcatflag
real*8 :: catfacs(50)=1d0       !! weighting of partial cross sections for proj.xsec.spi_all_001.dat, etc.
integer :: angularflag=0         !! for action 17,28 calculate fully differential partial ionization,
                                 !!   angular output in angprojspifile (atom/diatom only)
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
integer :: pulsewindowtoo=1      !! use window function for E(omega) as well
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

integer :: tentmode=0            !! tentmode=1 attempt not good

integer :: jacquaddir=0

integer, parameter :: nodgexpthirdflag=1  !! =1 HARDWIRE 10-2015 not sure about dgexpthird

!integer :: noorthogflag=1        !! 082010 NOW TURNING THIS ON !!
!     hardwire.  eliminated realproject for quad which didn't have the 
!     call to orthog.  don't remember if that was purposeful.
integer :: eigprintflag=0
!integer :: intopt=3              !! RK, GBS      !! SPF/VMF Integrator: 0, RK; 1, GBS, 2, DLSODPK  
!                                                 !!  for CMF: 3=expo 4=verlet
real*8 :: abserr
real*8 :: relerr=1.d-10          !!              !! relative error for integrator (RK/GBS/DLSODPK)
                                                 !!   norm error for expo (controls exposteps)
real*8 :: myrelerr=1.d-10        !!              !! absolute error set to norm*myrelerr
                                 !!  or for implicit: uses quadtol (dgmres parameter), myrelerr=orb change,
                                 !!     aerror=a change iteration thresholds; relerr=dgmres call threshold
!! number of MCSCF files
integer :: numfluxfiles=1
integer :: numfluxcurves(20)
!character (len=SLN) :: fluxfilenames(20)
!!integer :: whichsideproj=0       !! 0 if want to project on N-1 e- state, 1 if on N e- state

integer :: drivingmethod=0   !! 0 = V(t)  1 = H(t)-E
integer :: nonatrotate=0

!!!!!  integer :: whichquad=1           !!     0=my eqn 1=rayleigh quotient  WHICHQUAD 1 HARDWIRE

!!$ integer,parameter :: orderflag=0           !! Order=       !! ordering of configs. 1= first all alphas,  then betas;  
!!$                                                 !!    0= 1a1b2a2b etc.

!! INTERNAL

!!$ IMPLEMENT ME (DEPRECATE fluxinterval as namelist input)   
!!$ integer :: fluxsteps=1

real*8 :: lastfinish=0d0  !! final time of pulse

integer :: eachloaded(MXF)=(-99)

integer :: numpropsteps=1000000  !!              !! length of prop.  Overridden for pulse and relax.
integer, parameter :: maxmc=100
DATATYPE :: drivingenergies(maxmc)=0d0
integer :: numcurves=1
DATATYPE, allocatable :: elecweights(:,:,:,:), elecradii(:)
real*8 :: langramthresh=1d-9
integer :: numreduced=1
integer :: headersize=200

integer :: autosize, autosteps
integer :: cdenflag=0            !!              !! Calulate natconfig?  Not necessary if Act=6.
integer :: rdenflag=0            !!              !! Calulate denmat in R? 
real*8 :: globaltime=0.d0   !! for ease of output
!!! real*8 :: rcond=1.d-4      !! singular value for split operator cranck-nicholson
integer :: povplotflag=0
integer, parameter :: natorbfile=765
integer, parameter :: rnatorbfile=769
integer, parameter :: spfplotfile=303
integer, parameter :: denfile=904
integer, parameter :: denprojfile=597
integer, parameter :: natprojfile=488

integer :: skipflag=0

character (len=SLN) :: nullbuff

end module parameters


module mpi_orbsetmod
  implicit none
  integer :: mpi_orbset_init=0
  integer :: orbsperproc=(-1),  norbsets=(-1),  maxprocsperset=(-1),&
       nzprocsperset=(-1)
!! specific to each processor: 
  integer :: myorbset=-1, firstmpiorb=-1, orbrank=(-1)
  integer, allocatable :: MPI_COMM_ORB(:), NZ_COMM_ORB(:)
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
  integer :: nprocs=1, myrank
  character(len=SLN), parameter :: mpioutfilebase="MPIOUTS/MPI.Out."
  character(len=SLN) :: mpioutfile
  integer, parameter :: mpioutfilelen=16
  integer :: stdoutflag=0
  integer :: mpitime=0
  integer :: nonmpitime=0
  integer :: mpiatime=0
  integer :: mpibtime=0

end module mpimod



