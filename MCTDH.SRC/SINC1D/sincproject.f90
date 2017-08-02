
#include "Definitions.INC"

module myprojectmod
  implicit none

  DATATYPE, allocatable ::  threed_two(:)

  type fourmat
     DATATYPE, allocatable :: mat(:,:,:,:), tam(:,:,:,:)
  end type fourmat

  type twomat
     real*8, allocatable :: mat(:,:)
  end type twomat

  type onemat
     real*8, allocatable :: rmat(:)
     DATATYPE, allocatable :: cmat(:)
  end type onemat

  type(fourmat) :: ketot,fdtot
  type(twomat) ::  sinepoints
  type(onemat) :: kevect,fdvect

  DATATYPE, allocatable :: dipoles(:),&

!! WAS:  e.g. X(x) = x + i scalefunction(x,1)


       jacobian(:),&         !! jacobian(:,1) should only be a function of x, etc.
       invjacobian(:),&
       invsqrtjacobian(:),&
       scalediag(:),&
       invsqrtscaleweights(:)

end module myprojectmod


subroutine myprojectalloc()
  use myparams
  use pmpimod
  use pfileptrmod
  use myprojectmod
  implicit none

  allocate(dipoles(totpoints))
  dipoles=0
  
  if (scalingflag.ne.0) then
     allocate(          jacobian(totpoints),invjacobian(totpoints), &
          invsqrtjacobian(totpoints), &
          scalediag(totpoints),&
          invsqrtscaleweights(totpoints))
     jacobian=0; invjacobian=0; invsqrtjacobian=0; scalediag=0;
     invsqrtscaleweights=0;
  endif

  if (toepflag.eq.0) then
  
     if (numpoints*nbox.gt.10000) then
        OFLWR "WOW THAT'S BIG!  Are you sure you don't want to try toepflag?", numpoints*nbox; CFL
        OFLWR "WOW THAT'S BIG!  Are you sure you don't want to try toepflag?", numpoints*nbox; CFL
        OFLWR "WOW THAT'S BIG!  Are you sure you don't want to try toepflag?", numpoints*nbox; CFL
     endif

!! not needed now with tam.  old comment:
!! !! Allocating extra here for fdtot%mat and ketot%mat (+1's) --
!! !!   see Z/GEMM calls in coreproject.f90... leading dimension not
!! !!   allocated as passed to Z/GEMM without extra

     if (orbparflag) then
        allocate( &
!!             fdtot%mat(numpoints,nbox,numpoints,myrank:myrank +1), &
!!             ketot%mat(numpoints,nbox,numpoints,myrank:myrank +1))
             fdtot%mat(numpoints,nbox,numpoints,myrank:myrank), &
             ketot%mat(numpoints,nbox,numpoints,myrank:myrank),&
             fdtot%tam(numpoints,numpoints,nbox,myrank:myrank), &
             ketot%tam(numpoints,numpoints,nbox,myrank:myrank))
     else
        allocate( &
             fdtot%mat(numpoints,nbox,numpoints,nbox), &
             ketot%mat(numpoints,nbox,numpoints,nbox), &
             fdtot%tam(numpoints,numpoints,nbox,nbox), &
             ketot%tam(numpoints,numpoints,nbox,nbox))

     endif

     fdtot%mat=0; ketot%mat=0; fdtot%tam=0; ketot%tam=0
  endif

  allocate( &
       kevect%rmat(0-gridpoints:gridpoints-1),&
       kevect%cmat(0-gridpoints:gridpoints-1),&
       fdvect%rmat(0-gridpoints:gridpoints-1),&
       fdvect%cmat(0-gridpoints:gridpoints-1),&
       sinepoints%mat(numpoints,nbox))

  kevect%rmat=0; kevect%cmat=0;
  fdvect%rmat=0; kevect%cmat=0

  allocate(threed_two(0-numpoints:numpoints-1))
  threed_two=0

end subroutine myprojectalloc


module onedfunmod
contains
  
  function onedfun(inarray,num,incharge1,incharge2)
    use myparams
    implicit none
    integer,intent(in) :: num
    DATATYPE,intent(in) :: inarray(num)
    real*8, intent(in) :: incharge1, incharge2    
    DATATYPE :: onedfun(num)

    if (twomode==0) then
       onedfun = sechsq(inarray,num,incharge1,incharge2)
    elseif (twomode==1) then
       onedfun = softcoul(inarray,num,incharge1,incharge2)
    else
       onedfun = linearfun(inarray,num,incharge1,incharge2)
    endif
    
  end function onedfun

  function sechsq(inarray,num,incharge1,incharge2)
    use myparams
    use pfileptrmod
    implicit none
    integer,intent(in) :: num
    DATATYPE,intent(in) :: inarray(num)
    real*8, intent(in) :: incharge1, incharge2
    DATATYPE :: sechsq(num)
    real*8 :: xchargeprod, esoft
      
    xchargeprod = incharge1*incharge2

    if (sechmode.eq.0) then    !! all sech potentials have same range.
       esoft = softness
    else                       !! sech potentials are summed then squared to get united atom limit
       esoft = 2d0/xchargeprod
    endif

    !! with secmode.ne.0, potential is 3/4 Z^2 sech^2(x)
    
    sechsq(:) = ( 2d0/(exp(inarray(:)/esoft) + exp((-1)*inarray(:)/esoft)) )**2 * &
         0.5d0 * xchargeprod*(xchargeprod+1d0/esoft)

  end function sechsq


  function onlysech(inarray,num,icenter)
    use myparams
    use pfileptrmod
    implicit none
    integer,intent(in) :: num,icenter
    DATATYPE,intent(in) :: inarray(num)
    DATATYPE :: onlysech(num)
    real*8 :: esoft
    
    if (sechmode.eq.0) then
       esoft = softness
    else
       esoft = 2d0/nuccharges(icenter)
    endif

    onlysech(:) = &
         ( 2d0/(exp(inarray(:)/esoft) + exp((-1)*inarray(:)/esoft)) )
      
  end function onlysech

  function softcoul(inarray,num,incharge1,incharge2)
    use myparams
    use pfileptrmod
    implicit none
    integer,intent(in) :: num
    DATATYPE,intent(in) :: inarray(num)
    real*8, intent(in) :: incharge1, incharge2
    real*8 :: xp
    DATATYPE :: softcoul(num)
      
    xp = incharge1*incharge2;

    softcoul(:) = xp / sqrt(softness**2 + inarray**2)
      
  end function softcoul

    function linearfun(inarray,num,incharge1,incharge2)
    use myparams
    use pfileptrmod
    implicit none
    integer,intent(in) :: num
    DATATYPE,intent(in) :: inarray(num)
    real*8, intent(in) :: incharge1, incharge2
    real*8 :: xp
    DATATYPE :: linearfun(num)
      
    xp = incharge1*incharge2;

    linearfun(:) = (-0.5) * xp * (sqrt(inarray(:)**2 + softness**2) - softness)
      
  end function linearfun

  !! returns positive number
  function elecpot(inarray, insize, incharge, iskipcenter)
    use myparams
    use pfileptrmod
    implicit none
    integer,intent(in) :: insize, iskipcenter
    real*8, intent(in) :: incharge
    DATATYPE, intent(in) :: inarray(insize)
    DATATYPE :: elecpot(insize), pot(insize), myarray(insize), pot2(insize), potA(insize)
    real*8 :: esoft
    integer :: icenter
    
    pot(:)=0;    pot2(:)=0;   potA(:)=0
    
    if (twomode.ne.0.or.combinesech.eq.0) then
       !! coulomb or linear, twomode.ne.0 or sechmode.eq.0, old version

       do icenter=1,numcenters
          if (icenter.ne.iskipcenter) then
             myarray(:)=( inarray(:) - centershift(icenter)*spacing/2d0 )
             pot(:)=pot(:) +  onedfun(myarray(:), insize, nuccharges(icenter), incharge)
          endif
       enddo

    else  !! combinesech.ne.0 ::

       do icenter=1,numcenters
          esoft = 2d0/nuccharges(icenter)
          myarray(:)=( inarray(:) - centershift(icenter)*spacing/2d0 )
          potA(:)=potA(:) +  onlysech(myarray(:), insize, icenter) * nuccharges(icenter)
          if (icenter.ne.iskipcenter) then
             pot(:)=pot(:) +  onlysech(myarray(:), insize, icenter) * nuccharges(icenter)
             pot2(:)=pot2(:) +  onlysech(myarray(:), insize, icenter)**2 * nuccharges(icenter)**2
          endif
       enddo
          
       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !! THIS PRESCRIPTION GETS THE UNITED ATOM LIMIT CORRECT AS BEST POSSIBLE
       !! (analogous to Coulomb; charges add) WITHOUT ANY OTHER AD HOC CHOICE
       
       pot = incharge * 0.25 * ( 2*pot(:)*potA(:) + pot2(:) )

       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    endif
       
    elecpot(:) = pot(:)

  end function elecpot
  
end module onedfunmod


module gettwoemod
contains

subroutine get_twoe_new(pot)
  use myparams
  use pfileptrmod
  use myprojectmod  
  use pmpimod
  use onedfunmod
  implicit none
  DATATYPE,intent(out) :: pot(totpoints)
  DATATYPE,allocatable :: myarray(:)
  integer :: ii,jj,pp,gridoffset,istart
  real*8 :: xfac
  
  istart=1
  if (orbparflag.and.myrank.ne.1) then
     istart=0
  endif
  gridoffset=0; pp=1-gridpoints
  if (orbparflag) then
     gridoffset=(myrank-1)*numpoints
     pp=istart+2*gridoffset-gridpoints
  endif

  allocate(myarray(istart-numpoints:numpoints-1)); myarray=0

  jj=pp
  do ii=istart-numpoints,numpoints-1
     myarray(ii)=jj*spacing
     jj=jj+1
  enddo

  threed_two(:)=0d0

  if (twotype.eq.0) then
     !$$ infinite range     threed_two(istart-numpoints:numpoints-1) = twostrength
     if (istart-numpoints<0 .and. numpoints-1>0) then
        threed_two(0) = twostrength/spacing
     endif
  else
     xfac=1;
     if (twomode.eq.(-1)) then
       xfac = softness/softnesstwoe
     endif
     threed_two(istart-numpoints:numpoints-1) = twostrength * &
          onedfun(xfac*myarray(:),2*numpoints-istart, 1d0, 1d0) * xfac
  endif

  deallocate(myarray); allocate(myarray(numpoints))

  if (numcenters.eq.0) then
     pot(:) = 0.5d0 * harmstrength * dipoles(:)**2
  else
!! Add in after get orbs  
     pot(:)=0d0
  endif

  pot(:)=pot(:) - elecpot(dipoles(:), numpoints, 1d0, 0)
  
  deallocate(myarray)

end subroutine get_twoe_new

end module gettwoemod


subroutine op_yderiv(notint,notused1,notused2)
  use pfileptrmod
  implicit none
  integer :: notint
  DATATYPE :: notused1(*), notused2(*)
  OFLWR "WHAT! no op_yderiv sincdvr, not yet."; CFLST
  notused1(1)=0*notused2(1)
end subroutine op_yderiv

