
!! FILE CONTAINING **ALL** MODULES THAT ARE SHARED AMONG .F90 FILES !! 
!!
!!  THIS MAKES IT EASY TO COMPILE -- ALL THINGS DEPEND ON THIS AND MAYBE JUST A FEW OTHER FILES
!!

#include "Definitions.INC"


!! DATA TYPE FOR CONSTANT MEAN FIELD PROPAGATION.  TYPE XARR IS WHERE A WAVEFUNCTION AND ASSOCIATED STUFF
!!   IS STORED -- FOR INSTANCE MEAN FIELDS.  FOR CMF/LMF WE NEED VARIOUS THINGS FOR PREVIOUS TIMESTEPS AND
!!   THIS KEEPS THAT ALL SORTED OUT.  THERE ARE POSSIBLY SEVERAL CONCURRENTLY PROPAGATED WAVE FUNCTIONS
!!   (NUMPROP .GT. 1, E.G. FOR MCSCF)
!! I USE A VARIABLE XXX, AN ARRAY OF TYPE XARR DATA STRUCTURES, WHICH IS IN XARR.F90
!!   XXX IS HIDDEN TO MOST ROUTINES.  POINTER YYY POINTS TO CURRENT INDEX IN XXX VIA XXXPTR_INIT
!!   THE CURRENT INDEX IS IGLOBALPROP, INTEGER IN IGLOBALMOD WHICH IS A MODULE USED ONLY SPARINGLY 
!!   IN HIGHER LEVELS OF PROGRAM.  YYY CAN BE USED WHENEVER NEEDED AND IS IN YYYMOD BELOW
!! SEE XARR.F90 FOR FURTHER EXPLANATION

module configptrmod
  type CONFIGPTR

     DATATYPE, allocatable :: &
          xpotmatel(:,:),xopmatel(:,:),xconmatel(:,:),xymatel(:,:), xpulsematelxx(:,:), &
          xpulsematelyy(:,:), xpulsematelzz(:,:), xtwoematel(:,:,:,:), &
          xconmatelxx(:,:),xconmatelyy(:,:),xconmatelzz(:,:)

     DATATYPE :: xpulsenuc(3)=0
     DATATYPE :: kefac=1          !! scalar by which R ke is mutiplied
     DATATYPE :: constfac=1       !! scalar by constant terms are multiplied

!! for holes
     DATATYPE,allocatable :: xtwoe_htrace(:,:)

  end type CONFIGPTR
end module configptrmod



module sparseptrmod
  type SPARSEPTR

     DATATYPE, allocatable :: &
          xpotsparsemattr(:,:),          xpottracemattr(:,:),&
          xopsparsemattr(:,:),xonepotsparsemattr(:,:), &
          xconsparsemattr(:,:),xysparsemattr(:,:), xpulsesparsemattrxx(:,:), &
          xpulsesparsemattryy(:,:), xpulsesparsemattrzz(:,:), &
          xconsparsemattrxx(:,:),xconsparsemattryy(:,:),xconsparsemattrzz(:,:)

     DATATYPE :: xpulsenuc(3)=0
     DATATYPE :: kefac=1          !! scalar by which R ke is mutiplied
     DATATYPE :: constfac=1       !!                 constant terms are multiplied

  end type SPARSEPTR
end module sparseptrmod



module xmod
  use configptrmod
  use sparseptrmod
  type xarr
     DATATYPE, allocatable :: &  

          denmat(:,:,:)  , &           !! DENSITY MATRIX!
          invdenmat(:,:,:), &  !!          !! INVERSE !!
          reducedinvr(:,:,:),  &
          reducedr(:,:,:),   &
          reducedinvrsq(:,:,:), &
          reducedproderiv(:,:,:),  &
          reducedpot(:,:,:,:), &
          reducedpottally(:,:,:,:,:)

     Type(CONFIGPTR), allocatable :: cptr(:)
     Type(SPARSEPTR), pointer :: sptrptr(:)

     CNORMTYPE, allocatable :: &
          denvals(:)
     DATATYPE, allocatable ::    &
          cmfspfs(:,:), &
          cmfavec(:,:,:),&
          denvects(:,:) 
     DATATYPE, allocatable :: frozenexchinvr(:,:,:), frozenexchmat(:,:,:)

     DATATYPE, allocatable :: fockmatrix(:,:,:)

  end type xarr

end module xmod

module xxxmod
  use xmod
  implicit none
  Type(xarr) :: yyy
  Type(SPARSEPTR), target, allocatable :: yyysptr(:), yyysdfptr(:), yyysfdptr(:)
end module xxxmod


!! USE THIS MODULE FOR ACCESS TO WAVE FUNCTION AND ASSOCIATED DATA IN DATA TYPE XARR, VIA YYY%CMFSPFS() FOR INSTANCE

!! CONFIGPTRs are used to pass pieces of configuration matrix elements, for sparseconfigflag=1, 
!!  or whole R x elec configuration matel for sparseconfigflag=0


!! AARR returns the 2D space, spin index from condensed spinorbital index (which determines slater determinant order)
!! It is an array valued function and so must be included in a module that is included in any routines that use it

module aarrmod
contains
  function aarr(xind)
    implicit none
    integer, dimension(2) :: aarr
    integer, save :: temp(2)
    integer :: ind,q, xind
    ind=xind-1
!!$    if (orderflag==1) then
!!$       temp(1)           = mod(ind,numspf)+1
!!$       q=(ind-mod(ind,numspf))/numspf
!!$       temp(2)           = q+1
!!$       aarr=temp
!!$    else

       temp(2)           = mod(ind,2)+1
       q=(ind-mod(ind,2))/2
       temp(1)           = q+1
       aarr=temp

!!$    endif
  end function aarr
end module aarrmod



module jactimingmod
  implicit none
  integer :: times(20), numcalledhere=0
end module jactimingmod



module linearmod
  implicit none
  real*8 :: lasttime, firsttime
  integer :: effective_cmf_linearflag

end module linearmod

!! CONTAINS ORBITAL MATRIX ELEMENTS AND ASSOCIATED STUFF
!! POINTERS ARE NECESSARY BECAUSE HAVE SAME NAMED ARRAYS IN BOTH H2PROJECTMOD AND HEPROJECTMOD
!! THUS VIA HESPARSEASSIGN OR H2SPARSEASSIGN I POINT THE POINTERS TO THE APPROPRIATE ARRAYS

module opmod
use configptrmod
  implicit none

  DATATYPE,  allocatable :: &
       rkemod(:,:), &                         !! KE (diagonal) in bond length R
       proderivmod(:,:) !, &                    !! Prolate nuclear off diagonal

  DATATYPE, allocatable :: &
       pot(:), &                     !!  One-e potential
       halfniumpot(:)                     !!  One-e potential

  DATATYPE, allocatable :: frozenspfs(:,:)
  DATATYPE :: frozenpotdiag=0d0  !! potential matrix element
  DATATYPE :: frozenkediag=0d0  !! ke matrix element

  DATATYPE, allocatable ::   twoereduced(:,:,:),&        !! numerad,lbig+1,-2*mbig:2*mbig,nspf,nspf 
                                                         !! that's reducedpotsize,nspf,nspf
       hatomreduced(:), frozenreduced(:)

end module opmod



!!PREVIOUS:
!! the sparse spin projector is done as follows.  For the configurations config1 which do 
!! not have S^2 matrix elements with any other configuration, if their diagonal S^2 
!! eigenvalue is equal to the desired one, then iispinset(config1) is set to 1.  
!! iispinset is otherwise set to zero.  Then sum over all contributing spinsets using
!! spinsetprojector.
!!NOW KEEPING 1x1 , no iispinset

module spinwalkmod
  implicit none

  type settype
     DATATYPE, allocatable :: mat(:,:), vects(:,:)
  end type settype

  type spintype
     integer, allocatable ::  &
          spinsetsize(:,:),spinsets(:,:,:), spinsetrank(:,:), spindfsets(:,:,:), spindfsetindex(:,:)
     integer :: maxspinsetsize=0,maxnumspinsets=0,maxnumspindfsets=0
     integer,allocatable :: numspinsets(:),numspindfsets(:)
     type(settype), allocatable :: spinsetprojector(:,:)

     integer :: spinrestrictval=0

     integer :: numspinconfig=0
     integer, allocatable :: spinsperproc(:), alltopspins(:), allbotspins(:)
     integer ::       maxspinsperproc
     integer :: botspin=-1,topspin=-1
     integer :: firstspinconfig,lastspinconfig,localnspin

     integer :: numspindfconfig=0
     integer, allocatable :: spindfsperproc(:), alltopspindfs(:), allbotspindfs(:)
     integer ::       maxspindfsperproc
     integer :: botdfspin=-1,topdfspin=-1

  end type spintype

end module spinwalkmod



module dfconmod
  implicit none

type dfcontype
  integer, allocatable :: dfincludedmask(:), dfincludedconfigs(:)
!!$SP, dfincludedindex(:)

!! careful!  parconfigsplit, configstart, configend, etc.

  integer :: numdfwalks=-1    !! all walks, for all configurations
  integer, allocatable :: dfwalkto(:)    !! excluded  
  integer, allocatable :: dfwalkfrom(:)  !! included, out of numdfconfigs, but 
                                         !!   only those on this processor if parconfigsplit=1.
  integer, allocatable :: dfwalkphase(:) 
  integer, allocatable :: includedorb(:), excludedorb(:)   !! this is the excitation
end type dfcontype

end module


module walkmod
  use dfconmod
  use spinwalkmod
  implicit none

  type walktype

     integer :: m_restrictflag=0, m_restrictval=0, m_restrictmin=-99999, m_restrictmax=99999,&
          ug_restrictflag=0, ug_restrictval=0

     integer :: holeflag=0   !! index holes not electrons

     integer :: dflevel=0, dfwalklevel=0, dfrestrictflag=0
     integer :: singlewalkflag=0, doublewalkflag=0

     integer :: parconsplit=0

     integer, allocatable :: &   !! LIST OF CONFIGURATIONS num2part, numconfig
       configlist(  :  ,   : )   !!  num2part = 2x numpart = 2x number of elec or holes
     integer, allocatable :: &
       configmvals(  : ), &      !! if spfrestrict, their mvalues : numconfig
       configugvals(  : ), &     !! ugvalues based on info in spfugvals
       configtypes( : )          !! Hamiltonian matrix elements only among same type
!!$SP       configorder(:)

     integer :: nspf=0, numelec=0, numpart=0, num2part=0
     integer :: allspinproject=0
     integer :: restrictms=0

     integer :: sparseconfigflag=0

     integer :: numconfig=-1
     integer,allocatable :: configsperproc(:), alltopconfigs(:), allbotconfigs(:),&
          nzconfsperproc(:), nzproclist(:)
     integer :: maxconfigsperproc
     integer :: botconfig=-1,topconfig=-1
     integer :: configstart=-1,configend=-1, startrank,endrank
     integer :: firstconfig,lastconfig,localnconfig
     integer :: totadim=-1
     integer :: NZ_COMM=(-1), nzprocs=(-1), nzrank=(-798)

     integer :: numdfconfigs=-1
     integer, allocatable :: dfconfsperproc(:), alltopdfconfigs(:), allbotdfconfigs(:)
     integer ::       maxdfconfsperproc
     integer :: botdfconfig=-1,topdfconfig=-1

     integer :: numbasis=0
     integer, allocatable :: basisperproc(:), allbotbasis(:), alltopbasis(:)
     integer ::       maxbasisperproc
     integer :: botbasis=-1,topbasis=-1
     
     integer :: numdfbasis=0
     integer, allocatable :: dfbasisperproc(:), allbotdfbasis(:), alltopdfbasis(:)
     integer ::       maxdfbasisperproc
     integer :: botdfbasis=-1,topdfbasis=-1

!!       WALKS         !!
!!       WALKS         !!

     integer :: maxtotsinglewalks=0,maxtotdoublewalks=0 !! total for all configs on each proc
     integer :: singlemaxwalks=0, doublemaxwalks=0      !! old maxsinglewalks and maxdoublewalks
                                                        !! only used for listorder
     ! gives config2 for (iwalk,config1)       for sparse indexing
     integer, allocatable :: singlewalk(:),  scol(:)  
     ! gives walk number for (1:numsinglediagwalks,config1), max numsinglediagwalks=numelec
     integer, allocatable :: singlediag(:,:)  
     
     ! numconfig, numconfig - (iwalk, config1) - 
     integer, allocatable :: singlewalkdirphase(:)
     
     ! indices on matrix element - 1:2, numconfig,numconfig
     integer, allocatable :: singlewalkopspf(:,:)

     integer, allocatable :: numsinglewalks(:) ! numconfig
     integer, allocatable :: numsinglediagwalks(:) ! numconfig

     ! 1:4, numconfig, numconfig - (1:4, iwalk, config1)
     integer, allocatable :: doublewalkdirspf(:,:)

     ! iwalk, config1
     integer, allocatable :: doublewalkdirphase(:)

     integer, allocatable :: doublewalk(:),   dcol(:)
     ! gives walk number for (1:numdoublediagwalks,config1), max numdoublediagwalks=numelec*(numelec-1)
     integer, allocatable :: doublediag(:,:)
     integer, allocatable :: numdoublewalks(:)
     integer, allocatable :: numdoublediagwalks(:)

     type(spintype) :: sss
     type(dfcontype) :: ddd

     integer :: maxnumsinglehops=0
     integer,allocatable :: numsinglehops(:)
     integer,allocatable :: singlehop(:,:), singlehopwalkstart(:,:), singlehopwalkend(:,:)
     integer,allocatable :: singlediaghop(:)
     integer,allocatable :: singlehopdiagflag(:)
     integer,allocatable :: firstsinglehopbyproc(:,:), lastsinglehopbyproc(:,:)

     integer :: maxnumdoublehops=0
     integer,allocatable :: numdoublehops(:)
     integer,allocatable :: doublehop(:,:), doublehopwalkstart(:,:), doublehopwalkend(:,:)
     integer,allocatable :: doublediaghop(:)
     integer,allocatable :: doublehopdiagflag(:)
     integer,allocatable :: firstdoublehopbyproc(:,:), lastdoublehopbyproc(:,:)

     integer :: singlematsize=0,doublematsize=0

  end type walktype

end module walkmod


module configmod
  use walkmod
  implicit none
  type(walktype),target :: www
!!$  type(walktype),pointer:: bwwptr
  type(walktype),target :: bioww
  type(walktype),pointer:: dwwptr
  type(walktype),target :: dfww
  type(walktype),target :: fdww
  type(walktype),target :: catww
!! cation arrays dimensioned nspf,www%botdfconfig:www%topdfconfig.  
!! If value is (-1) then
!! there is no cation hop... (this should not happen given if df_restrictflag
!! is greater than zero given logic max(df_restrictflag-1,0) in mctdhf.f90...)
!! else value is cation configuration index
  integer :: numcathops = -1
  integer, allocatable :: cathopfrom(:), cathopto(:)
!! catphase(nspf,idfconfig)=0 no walk else +/-1
  integer,allocatable :: catphase(:)
  integer,allocatable :: catspf(:)
end module configmod


module configpropmod
  use configptrmod
  use sparseptrmod
  implicit none
  Type(CONFIGPTR) :: workconfigpointer
  Type(SPARSEPTR),target :: worksparsepointer, workdfsparsepointer, workfdsparsepointer
  Type(SPARSEPTR),pointer :: worksparsepointerptr

end module configpropmod

subroutine configpropalloc()
  use parameters
  use configpropmod
  use configmod
  implicit none

  call configptralloc(workconfigpointer,www)  !! nspf not used for regular walks (not walks2)

  if (sparseopt.ne.0) then  !! (sparseconfigflag is also 0, see getparams)
     call sparseptralloc(worksparsepointer,www)  !! nspf not used for regular walks (not walks2)

     if (use_dfwalktype) then
        if (shuffle_dfwalktype) then
           call sparseptralloc(workfdsparsepointer,fdww)
        else
           call sparseptralloc(workdfsparsepointer,dfww)
        endif
     endif
  endif


end subroutine configpropalloc


subroutine configpropdealloc()
  use parameters
  use configpropmod
  implicit none
  call configptrdealloc(workconfigpointer)
  if (sparseopt.ne.0) then
     call sparseptrdealloc(worksparsepointer)
     if (use_dfwalktype) then
        if (shuffle_dfwalktype) then
           call sparseptrdealloc(workfdsparsepointer)
        else
           call sparseptrdealloc(workdfsparsepointer)
        endif
     endif
  endif

end subroutine configpropdealloc





module orblabelmod
contains
  function orblabel(index)
    use aarrmod
    use parameters
    implicit none
    integer :: aarray(2), index
    character (len=2) :: mslabels(2) =["a ","b "]
    character (len=6) :: orblabel, label
    label="      "
    aarray=aarr(index)

!! does not work on FRANKLIN.
#ifndef PGFFLAG
    write(label,'(I4,A2)') aarray(1), mslabels(aarray(2))
#endif
    orblabel=label
  end function orblabel
end module orblabelmod


module natprojmod
  implicit none
  DATATYPE, allocatable :: natproj(:,:,:), natconfigs(:,:), natderiv(:)
  CNORMTYPE, allocatable :: natvals(:), natdot(:,:), curves(:,:)
  integer :: nnalloc=0
end module natprojmod



subroutine add_cptr(aptr,bptr,sumptr,afac,bfac)
  use configptrmod
  use parameters
  implicit none
  DATATYPE,intent(in) :: afac,bfac
  Type(CONFIGPTR),intent(in) :: aptr,bptr
  Type(CONFIGPTR),intent(inout) :: sumptr
  call add_cptr0(aptr,bptr,sumptr, afac,bfac, afac,bfac, afac,bfac, afac,bfac )
end subroutine



subroutine assign_cptr(outptr,inptr,fac)
  use configptrmod
  use parameters
  implicit none
  DATATYPE,intent(in) :: fac
  Type(CONFIGPTR),intent(in) :: inptr
  Type(CONFIGPTR),intent(inout) :: outptr

  outptr%kefac = inptr%kefac   * fac
  outptr%constfac = inptr%constfac   * fac
  outptr%xpulsenuc(:) = inptr%xpulsenuc(:) * fac

  outptr%xtwoe_htrace(:,:) = inptr%xtwoe_htrace(:,:) * fac

  outptr%xpotmatel(:,:)=inptr%xpotmatel(:,:)   * fac
  outptr%xopmatel(:,:)=inptr%xopmatel(:,:)   * fac
  outptr%xconmatel(:,:)=inptr%xconmatel(:,:)   * fac
  outptr%xconmatelxx(:,:)=inptr%xconmatelxx(:,:)   * fac
  outptr%xconmatelyy(:,:)=inptr%xconmatelyy(:,:)   * fac
  outptr%xconmatelzz(:,:)=inptr%xconmatelzz(:,:)   * fac
  outptr%xymatel(:,:)=inptr%xymatel(:,:)   * fac
  outptr%xpulsematelxx(:,:)=inptr%xpulsematelxx(:,:)   * fac
  outptr%xpulsematelyy(:,:)=inptr%xpulsematelyy(:,:)   * fac
  outptr%xpulsematelzz(:,:)=inptr%xpulsematelzz(:,:)   * fac
  outptr%xpulsenuc=inptr%xpulsenuc   * fac
  outptr%xtwoematel(:,:,:,:)=inptr%xtwoematel(:,:,:,:)   * fac

end subroutine assign_cptr


subroutine add_cptr0(aptr,bptr,sumptr,afacbo,bfacbo,  afacnuc,bfacnuc,   afacpulse,bfacpulse, afaccon,bfaccon)
  use configptrmod
  use parameters
  implicit none
  Type(CONFIGPTR),intent(in) :: aptr,bptr
  Type(CONFIGPTR),intent(inout) :: sumptr
  DATATYPE,intent(in) :: afacbo,bfacbo,  afacnuc,bfacnuc,   afacpulse,bfacpulse, afaccon,bfaccon

  sumptr%kefac = afacnuc*aptr%kefac + bfacnuc*bptr%kefac
  sumptr%constfac = afacnuc*aptr%constfac + bfacnuc*bptr%constfac
  sumptr%xpulsenuc(:) = afacnuc*aptr%xpulsenuc(:) + bfacnuc*bptr%xpulsenuc(:)

  sumptr%xtwoe_htrace(:,:) = afacbo*aptr%xtwoe_htrace(:,:) + bfacbo*bptr%xtwoe_htrace(:,:)

  sumptr%xpotmatel(:,:)=aptr%xpotmatel(:,:)*afacbo                +bptr%xpotmatel(:,:)*bfacbo
  sumptr%xopmatel(:,:)=aptr%xopmatel(:,:)*afacbo                  +bptr%xopmatel(:,:)*bfacbo
  sumptr%xconmatel(:,:)=aptr%xconmatel(:,:)*afaccon               +bptr%xconmatel(:,:)*bfaccon
  sumptr%xconmatelxx(:,:)=aptr%xconmatelxx(:,:)*afaccon               +bptr%xconmatelxx(:,:)*bfaccon
  sumptr%xconmatelyy(:,:)=aptr%xconmatelyy(:,:)*afaccon               +bptr%xconmatelyy(:,:)*bfaccon
  sumptr%xconmatelzz(:,:)=aptr%xconmatelzz(:,:)*afaccon               +bptr%xconmatelzz(:,:)*bfaccon
  sumptr%xymatel(:,:)=aptr%xymatel(:,:)*afacnuc                   +bptr%xymatel(:,:)*bfacnuc
  sumptr%xpulsematelxx(:,:)=aptr%xpulsematelxx(:,:)*afacpulse         +bptr%xpulsematelxx(:,:)*bfacpulse
  sumptr%xpulsematelyy(:,:)=aptr%xpulsematelyy(:,:)*afacpulse         +bptr%xpulsematelyy(:,:)*bfacpulse
  sumptr%xpulsematelzz(:,:)=aptr%xpulsematelzz(:,:)*afacpulse         +bptr%xpulsematelzz(:,:)*bfacpulse
  sumptr%xpulsenuc=aptr%xpulsenuc*afacpulse         +bptr%xpulsenuc*bfacpulse
  sumptr%xtwoematel(:,:,:,:)=aptr%xtwoematel(:,:,:,:)*afacbo      +bptr%xtwoematel(:,:,:,:)*bfacbo


end subroutine add_cptr0


subroutine zero_cptr(outptr)
  use configptrmod
  use parameters
  implicit none
  Type(CONFIGPTR),intent(inout) :: outptr

  outptr%kefac=0d0
  outptr%constfac=0d0
  outptr%xpulsenuc(:)=0

  outptr%xtwoe_htrace(:,:) = 0

  outptr%xpotmatel(:,:)=0d0
  outptr%xopmatel(:,:)=0d0
  outptr%xconmatel(:,:)=0d0
  outptr%xconmatelxx(:,:)=0d0
  outptr%xconmatelyy(:,:)=0d0
  outptr%xconmatelzz(:,:)=0d0
  outptr%xymatel(:,:)=0d0
  outptr%xpulsematelxx(:,:)=0d0
  outptr%xpulsematelyy(:,:)=0d0
  outptr%xpulsematelzz(:,:)=0d0
  outptr%xpulsenuc=0d0
  outptr%xtwoematel(:,:,:,:)=0d0

end subroutine zero_cptr


subroutine configptralloc(inptr,www)
  use configptrmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  Type(CONFIGPTR) :: inptr

  allocate(inptr%xpotmatel(www%nspf,www%nspf))
  allocate(inptr%xopmatel(www%nspf,www%nspf))
  allocate(inptr%xpulsematelxx(www%nspf,www%nspf))
  allocate(inptr%xpulsematelyy(www%nspf,www%nspf))
  allocate(inptr%xpulsematelzz(www%nspf,www%nspf))
  allocate(inptr%xconmatel(www%nspf,www%nspf))
  allocate(inptr%xconmatelxx(www%nspf,www%nspf))
  allocate(inptr%xconmatelyy(www%nspf,www%nspf))
  allocate(inptr%xconmatelzz(www%nspf,www%nspf))
  allocate(inptr%xymatel(www%nspf,www%nspf))
  allocate(inptr%xtwoematel(www%nspf,www%nspf,www%nspf,www%nspf))
  inptr%xpotmatel=0;  inptr%xopmatel=0;  inptr%xpulsematelxx=0;  inptr%xpulsematelyy=0;
  inptr%xpulsematelzz=0;   inptr%xconmatel=0;  inptr%xconmatelxx=0;  inptr%xconmatelyy=0;
  inptr%xconmatelzz=0;   inptr%xymatel=0;  inptr%xtwoematel=0; 

  allocate(inptr%xtwoe_htrace(www%nspf,www%nspf))
  inptr%xtwoe_htrace=0

  inptr%kefac=1
  inptr%constfac=1
  inptr%xpulsenuc(:)=0

end subroutine configptralloc

subroutine configptrdealloc(inptr)
  use configptrmod
  Type(CONFIGPTR) :: inptr

  deallocate(inptr%xpotmatel)
  deallocate(inptr%xopmatel)
  deallocate(inptr%xpulsematelxx)
  deallocate(inptr%xpulsematelyy)
  deallocate(inptr%xpulsematelzz)
  deallocate(inptr%xconmatel)
  deallocate(inptr%xconmatelxx)
  deallocate(inptr%xconmatelyy)
  deallocate(inptr%xconmatelzz)
  deallocate(inptr%xymatel)
  deallocate(inptr%xtwoematel)
  deallocate(inptr%xtwoe_htrace)

  inptr%kefac=0
  inptr%constfac=0
  inptr%xpulsenuc(:)=0

end subroutine configptrdealloc





subroutine add_sptr(aptr,bptr,sumptr,afac,bfac,www)
  use sparseptrmod
  use parameters
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: afac,bfac
  Type(SPARSEPTR),intent(in) :: aptr,bptr
  Type(SPARSEPTR),intent(inout) :: sumptr
  call add_sptr0(aptr,bptr,sumptr, afac,bfac, afac,bfac, afac,bfac, afac,bfac,www )
end subroutine



subroutine assign_sptr(outptr,inptr,fac,www)
  use sparseptrmod
  use parameters
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: fac
  Type(SPARSEPTR),intent(in) :: inptr
  Type(SPARSEPTR),intent(inout) :: outptr

  outptr%kefac = inptr%kefac   * fac
  outptr%constfac = inptr%constfac   * fac
  outptr%xpulsenuc=inptr%xpulsenuc   * fac

  if (www%configend.ge.www%configstart) then
     outptr%xpotsparsemattr(:,:)=inptr%xpotsparsemattr(:,:)   * fac
     outptr%xpottracemattr(:,:)=inptr%xpottracemattr(:,:)   * fac

     outptr%xopsparsemattr(:,:)=inptr%xopsparsemattr(:,:)   * fac
     outptr%xonepotsparsemattr(:,:)=inptr%xonepotsparsemattr(:,:)   * fac
     outptr%xconsparsemattr(:,:)=inptr%xconsparsemattr(:,:)   * fac
     outptr%xconsparsemattrxx(:,:)=inptr%xconsparsemattrxx(:,:)   * fac
     outptr%xconsparsemattryy(:,:)=inptr%xconsparsemattryy(:,:)   * fac
     outptr%xconsparsemattrzz(:,:)=inptr%xconsparsemattrzz(:,:)   * fac
     outptr%xysparsemattr(:,:)=inptr%xysparsemattr(:,:)   * fac
     outptr%xpulsesparsemattrxx(:,:)=inptr%xpulsesparsemattrxx(:,:)   * fac
     outptr%xpulsesparsemattryy(:,:)=inptr%xpulsesparsemattryy(:,:)   * fac
     outptr%xpulsesparsemattrzz(:,:)=inptr%xpulsesparsemattrzz(:,:)   * fac
  endif

end subroutine assign_sptr


subroutine add_sptr0(aptr,bptr,sumptr,afacbo,bfacbo,  afacnuc,bfacnuc,   afacpulse,bfacpulse, afaccon,bfaccon,www)
  use sparseptrmod
  use parameters
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  Type(SPARSEPTR),intent(in) :: aptr,bptr
  Type(SPARSEPTR),intent(inout) :: sumptr
  DATATYPE,intent(in) :: afacbo,bfacbo,  afacnuc,bfacnuc,   afacpulse,bfacpulse, afaccon,bfaccon

  sumptr%kefac = afacnuc*aptr%kefac + bfacnuc*bptr%kefac
  sumptr%constfac = afacnuc*aptr%constfac + bfacnuc*bptr%constfac
  sumptr%xpulsenuc=aptr%xpulsenuc*afacpulse         +bptr%xpulsenuc*bfacpulse

  if (www%configend.ge.www%configstart) then
     sumptr%xpotsparsemattr(:,:)=aptr%xpotsparsemattr(:,:)*afacbo    +bptr%xpotsparsemattr(:,:)*bfacbo
     sumptr%xpottracemattr(:,:)=aptr%xpottracemattr(:,:)*afacbo    +bptr%xpottracemattr(:,:)*bfacbo

     sumptr%xopsparsemattr(:,:)=aptr%xopsparsemattr(:,:)*afacbo      +bptr%xopsparsemattr(:,:)*bfacbo
     sumptr%xonepotsparsemattr(:,:)=aptr%xonepotsparsemattr(:,:)*afacbo      +bptr%xonepotsparsemattr(:,:)*bfacbo
     sumptr%xconsparsemattr(:,:)=aptr%xconsparsemattr(:,:)*afaccon   +bptr%xconsparsemattr(:,:)*bfaccon
     sumptr%xconsparsemattrxx(:,:)=aptr%xconsparsemattrxx(:,:)*afaccon   +bptr%xconsparsemattrxx(:,:)*bfaccon
     sumptr%xconsparsemattryy(:,:)=aptr%xconsparsemattryy(:,:)*afaccon   +bptr%xconsparsemattryy(:,:)*bfaccon
     sumptr%xconsparsemattrzz(:,:)=aptr%xconsparsemattrzz(:,:)*afaccon   +bptr%xconsparsemattrzz(:,:)*bfaccon
     sumptr%xysparsemattr(:,:)=aptr%xysparsemattr(:,:)*afacnuc       +bptr%xysparsemattr(:,:)*bfacnuc
     sumptr%xpulsesparsemattrxx(:,:)=aptr%xpulsesparsemattrxx(:,:)*afacpulse  +bptr%xpulsesparsemattrxx(:,:)*bfacpulse
     sumptr%xpulsesparsemattryy(:,:)=aptr%xpulsesparsemattryy(:,:)*afacpulse  +bptr%xpulsesparsemattryy(:,:)*bfacpulse
     sumptr%xpulsesparsemattrzz(:,:)=aptr%xpulsesparsemattrzz(:,:)*afacpulse  +bptr%xpulsesparsemattrzz(:,:)*bfacpulse
  endif

end subroutine add_sptr0


subroutine zero_sptr(outptr,www)
  use sparseptrmod
  use parameters
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  Type(SPARSEPTR) :: outptr

  outptr%kefac=0d0
  outptr%constfac=0d0
  outptr%xpulsenuc=0d0

  if (www%configend.ge.www%configstart) then
     outptr%xpotsparsemattr(:,:)=0d0
     outptr%xpottracemattr(:,:)=0d0

     outptr%xopsparsemattr(:,:)=0d0
     outptr%xonepotsparsemattr(:,:)=0d0
     outptr%xconsparsemattr(:,:)=0d0
     outptr%xconsparsemattrxx(:,:)=0d0
     outptr%xconsparsemattryy(:,:)=0d0
     outptr%xconsparsemattrzz(:,:)=0d0
     outptr%xysparsemattr(:,:)=0d0
     outptr%xpulsesparsemattrxx(:,:)=0d0
     outptr%xpulsesparsemattryy(:,:)=0d0
     outptr%xpulsesparsemattrzz(:,:)=0d0
  endif


end subroutine zero_sptr


subroutine sparseptralloc(inptr,www)
  use fileptrmod
  use sparse_parameters
  use sparseptrmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  Type(SPARSEPTR) :: inptr

  if (sparseconfigflag.eq.0) then
     OFLWR "error, sparseptralloc called when sparseconfigflag.eq.0"; CFLST
  endif

!!  OFLWR "  ... go alloc sparsopt=1 .. "; CFL
!!  call mpibarrier()

  allocate(inptr%xpotsparsemattr    (www%doublematsize,www%configstart:www%configend))
  allocate(inptr%xpottracemattr     (www%singlematsize,www%configstart:www%configend))

  allocate(inptr%xopsparsemattr     (www%singlematsize,www%configstart:www%configend))
  allocate(inptr%xonepotsparsemattr (www%singlematsize,www%configstart:www%configend))
  allocate(inptr%xpulsesparsemattrxx(www%singlematsize,www%configstart:www%configend))
  allocate(inptr%xpulsesparsemattryy(www%singlematsize,www%configstart:www%configend))
  allocate(inptr%xpulsesparsemattrzz(www%singlematsize,www%configstart:www%configend))
  allocate(inptr%xconsparsemattr    (www%singlematsize,www%configstart:www%configend))
  allocate(inptr%xconsparsemattrxx  (www%singlematsize,www%configstart:www%configend))
  allocate(inptr%xconsparsemattryy  (www%singlematsize,www%configstart:www%configend))
  allocate(inptr%xconsparsemattrzz  (www%singlematsize,www%configstart:www%configend))
  allocate(inptr%xysparsemattr      (www%singlematsize,www%configstart:www%configend))

  if (www%configend.ge.www%configstart) then
     inptr%xpotsparsemattr=0; inptr%xpottracemattr=0;
     inptr%xopsparsemattr=0; inptr%xonepotsparsemattr=0; 
     inptr%xpulsesparsemattrxx=0; inptr%xpulsesparsemattryy=0; inptr%xpulsesparsemattrzz=0; 
     inptr%xconsparsemattr=0; inptr%xconsparsemattrxx=0; inptr%xconsparsemattryy=0; 
     inptr%xconsparsemattrzz=0; inptr%xysparsemattr=0
  endif

!!  call mpibarrier()
!!  OFLWR "   ... ok alloc sparseopt=1"; CFL

end subroutine sparseptralloc

subroutine sparseptrdealloc(inptr)
  use sparseptrmod
  Type(SPARSEPTR) :: inptr

  deallocate(inptr%xpotsparsemattr)
  deallocate(inptr%xpottracemattr)
  deallocate(inptr%xopsparsemattr)
  deallocate(inptr%xonepotsparsemattr)
  deallocate(inptr%xpulsesparsemattrxx)
  deallocate(inptr%xpulsesparsemattryy)
  deallocate(inptr%xpulsesparsemattrzz)
  deallocate(inptr%xconsparsemattr)
  deallocate(inptr%xconsparsemattrxx)
  deallocate(inptr%xconsparsemattryy)
  deallocate(inptr%xconsparsemattrzz)
  deallocate(inptr%xysparsemattr)

end subroutine sparseptrdealloc





subroutine opalloc()
  use parameters
  use opmod
  use configmod
  use mpi_orbsetmod
  implicit none

  allocate(rkemod(numr,numr), proderivmod(numr,numr))   
  rkemod=0.d0; proderivmod=0.d0; 

  if (parorbsplit.eq.1) then
     allocate(twoereduced(reducedpotsize,nspf,firstmpiorb:firstmpiorb+orbsperproc-1))
  else
     allocate(twoereduced(reducedpotsize,nspf,nspf))
  endif
  twoereduced(:,:,:)=0d0

  allocate(   pot(spfsize),  halfniumpot(spfsize) )
  pot(:)=0; halfniumpot(:)=0

  allocate(frozenspfs(spfsize,max(1,numfrozen)))
  frozenspfs(:,:)=0
  allocate(hatomreduced(reducedpotsize),frozenreduced(reducedpotsize))
  hatomreduced=0; frozenreduced=0

end subroutine opalloc

subroutine opdealloc()
  use opmod
  use parameters
  implicit none
  deallocate(rkemod,proderivmod,   pot, halfniumpot)
  deallocate(frozenspfs,hatomreduced,frozenreduced)

end subroutine opdealloc

subroutine natprojalloc
  use parameters
  use natprojmod
  use configmod
  implicit none
  if (par_consplit.ne.0) then
     OFLWR "Error, no natproj parconsplit.ne.0"; CFLST
  endif
  if (nnalloc==0) then
     allocate( natproj(numr,numr,mcscfnum), natconfigs(num_config,numr), natderiv(num_config),&
          natvals(num_config), natdot(numr,mcscfnum), curves(numr,numr))
     natproj=0; natconfigs=0; natderiv=0; natvals=0; natdot=0; curves=0
  endif
  nnalloc=1
end subroutine

subroutine natprojdealloc
  use natprojmod
  use parameters
  implicit none
  if (nnalloc.ne.0) then
     deallocate( natproj, natconfigs, natderiv, natvals, natdot, curves)
  endif
  nnalloc=0
end subroutine

module configexpotimemod
  implicit none
  real*8 :: configexpotime = 0d0
  logical :: initialized = .false.
  integer :: imc = -99
end module configexpotimemod

subroutine configexpoinit(intime,inmc)
  use configexpotimemod
  implicit none
  real*8,intent(in) :: intime
  integer,intent(in) :: inmc
  configexpotime=intime;  initialized=.true.
  imc=inmc
end subroutine configexpoinit




!!  NEVER CHANGE IGLOBALPROP WITHOUT CALLING XXXPTR_INIT !!!!!

subroutine xalloc()
  use parameters
  use configmod
  use xxxmod
  use mpi_orbsetmod
  implicit none
  integer :: ii

  allocate(yyy%denvals(nspf), yyy%denvects(nspf,nspf),yyy%cptr(0:numreduced))
  yyy%denvals=0; yyy%denvects=0

  do ii=0,numreduced
     call configptralloc(yyy%cptr(ii),www)
  enddo
  allocate(yyysptr(0:numreduced))
  if (sparseopt.ne.0) then
     do ii=0,numreduced
        call sparseptralloc(yyysptr(ii),www)
     enddo
  endif
  if (use_dfwalktype) then
     if (shuffle_dfwalktype) then
        allocate(yyysfdptr(0:numreduced))
     else
        allocate(yyysdfptr(0:numreduced))
     endif
     if (sparseopt.ne.0) then
        do ii=0,numreduced
           if (shuffle_dfwalktype) then
              call sparseptralloc(yyysfdptr(ii),fdww)
           else
              call sparseptralloc(yyysdfptr(ii),dfww)
           endif
        enddo
     endif
  endif

  allocate( yyy%denmat( nspf, nspf,0:numreduced), &
       yyy%invdenmat( nspf, nspf,0:numreduced), &
       yyy%reducedinvr( nspf, nspf, 0:numreduced), &
       yyy%reducedr( nspf, nspf, 0:numreduced),&
       yyy%reducedinvrsq( nspf, nspf, 0:numreduced),&
       yyy%reducedproderiv(nspf,nspf, 0:numreduced),  &
       yyy%cmfspfs(totspfdim,0:numreduced),&
       yyy%cmfavec(tot_adim,mcscfnum,0:numreduced),&
       yyy%reducedpottally(nspf,nspf, nspf, nspf, 0:numreduced))

  if (parorbsplit.eq.1) then
     allocate(yyy%reducedpot(reducedpotsize,  nspf, &
          firstmpiorb:firstmpiorb+orbsperproc-1, 0:numreduced))
  else
     allocate(yyy%reducedpot(reducedpotsize,  nspf, nspf, 0:numreduced))
  endif

  yyy%denmat=0; yyy%invdenmat=0; yyy%reducedr=0; yyy%reducedinvrsq=0; yyy%reducedpot=0;
  yyy%reducedpottally=0; yyy%reducedproderiv=0; yyy%cmfspfs=0; 
  if (tot_adim.gt.0) then
     yyy%cmfavec=0
  endif
  if (numfrozen.gt.0) then
     allocate(yyy%frozenexchinvr(spfsize,nspf,0:numreduced),yyy%frozenexchmat(nspf,nspf,0:numreduced))
     yyy%frozenexchinvr=0; yyy%frozenexchmat=0
  endif

  if (use_fockmatrix) then
     allocate(yyy%fockmatrix(nspf,nspf,0:numreduced))
     yyy%fockmatrix=0
  endif

end subroutine xalloc


subroutine xdealloc()
  use parameters
  use mpimod
  use xxxmod
  implicit none
  integer :: ii

  deallocate(yyy%denvects,  yyy%denvals)
  do ii=0,numreduced
     call configptrdealloc(yyy%cptr(ii))
  enddo
  if (sparseopt.ne.0) then
     do ii=0,numreduced
        call sparseptrdealloc(yyysptr(ii))
     enddo
  endif
  if (use_dfwalktype) then
     if (sparseopt.ne.0) then
        do ii=0,numreduced
           if (shuffle_dfwalktype) then
              call sparseptrdealloc(yyysfdptr(ii))
           else
              call sparseptrdealloc(yyysdfptr(ii))
           endif
        enddo
     endif
  endif

  deallocate(yyy%denmat, yyy%invdenmat, yyy%reducedinvr, yyy%reducedr,&
       yyy%reducedinvrsq, yyy%reducedpot, yyy%reducedpottally, yyy%reducedproderiv,&
       yyy%cmfspfs,yyy%cmfavec)

  if (numfrozen.gt.0) then
     deallocate(yyy%frozenexchinvr, yyy%frozenexchmat)
  endif

  if (use_fockmatrix) then
     deallocate(yyy%fockmatrix)
  endif

end subroutine xdealloc


module savenormmod
  implicit none
  DATATYPE,allocatable :: savenorms(:,:)

contains

  subroutine get_bonorms(howmany,vectors,outnorms)
    use sparse_parameters
    use basis_parameters
    use r_parameters
    use dotmod
    use mpisubmod
    implicit none
    integer,intent(in) :: howmany
    DATATYPE,intent(in) :: vectors(numr,first_config:last_config,howmany)
    DATATYPE,intent(out) :: outnorms(numr,howmany)
    DATATYPE :: tempnorms(numr,howmany),tempvec(first_config:last_config)   !! AUTOMATIC
    integer :: kk,mm
    do mm=1,howmany
       do kk=1,numr
          tempvec(:)=vectors(kk,:,mm)
          tempnorms(kk,mm)=dot(tempvec,tempvec,local_nconfig)
       enddo
    enddo
    if (par_consplit.ne.0) then
       call mympireduce(tempnorms,howmany*numr)
    endif
    outnorms(:,:)=sqrt(tempnorms)
  end subroutine get_bonorms

  subroutine enforce_bonorms(howmany,vectors,innorms)
    use sparse_parameters
    use basis_parameters
    use r_parameters
    use dotmod
    implicit none
    integer,intent(in) :: howmany
    DATATYPE,intent(inout) :: vectors(numr,first_config:last_config,howmany)
    DATATYPE,intent(in) :: innorms(numr,howmany)
    DATATYPE :: tempnorms(numr,howmany),tempvec(first_config:last_config)   !! AUTOMATIC
    integer :: kk,mm
    call get_bonorms(howmany,vectors,tempnorms)
    do mm=1,howmany
       do kk=1,numr
          tempvec(:)=vectors(kk,:,mm)
          vectors(kk,:,mm)=tempvec(:)/tempnorms(kk,mm)*innorms(kk,mm)
       enddo
    enddo
  end subroutine enforce_bonorms

end module savenormmod

