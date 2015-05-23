
#include "Definitions.INC"

module myprojectmod
  implicit none

  integer :: threedtwosize=0
  DATATYPE, allocatable ::  threed_two(:,:,:,:)

  type fourmat
     DATATYPE, allocatable :: mat(:,:,:,:)
  end type fourmat

  type twomat
     real*8, allocatable :: mat(:,:)
  end type twomat

  type onemat
     real*8, allocatable :: rmat(:)
     DATATYPE, allocatable :: cmat(:)
  end type onemat

  type(fourmat), allocatable :: ketot(:),fdtot(:)
  type(twomat), allocatable ::  sinepoints(:)
  type(onemat),allocatable :: kevect(:),fdvect(:)

  DATATYPE, allocatable :: dipoles(:,:), &

!! smooth complex scaling:  e.g. X(x) = x + i scalefunction(x,1)
       scalefunction(:,:), &  !! scalefunction(:,1) and 
       jacobian(:,:),&        !! jacobian(:,1) should only be a function of x, etc.
       sqrtjacobian(:,:),&
       invsqrtjacobian(:,:),&
       invjacobian(:,:),&
       scaleweights(:),&
       invsqrtscaleweights(:)
  

end module myprojectmod


subroutine myprojectalloc()
  use myparams
  use myprojectmod
  implicit none
  integer :: idim

  allocate(dipoles(totpoints,griddim))
  if (scalingflag.ne.0) then
     allocate(scalefunction(totpoints,3),jacobian(totpoints,3), &
          invsqrtjacobian(totpoints,3), sqrtjacobian(totpoints,3), invjacobian(totpoints,3),&
          scaleweights(totpoints), invsqrtscaleweights(totpoints))
  endif
  allocate(ketot(griddim),sinepoints(griddim),kevect(griddim),fdtot(griddim),fdvect(griddim))
  do idim=1,griddim
     allocate( &

!! Allocating extra here for fdtot%mat and ketot%mat (+1's) --
!!   see Z/GEMM calls in coreproject.f90... leading dimension not
!!   allocated as passed to Z/GEMM without extra

          fdtot(idim)%mat(numpoints(idim),nbox(idim),numpoints(idim),nbox(idim)   +1), &
          ketot(idim)%mat(numpoints(idim),nbox(idim),numpoints(idim),nbox(idim)   +1), &
          kevect(idim)%rmat(1-gridpoints(idim):gridpoints(idim)-1),&
          kevect(idim)%cmat(1-gridpoints(idim):gridpoints(idim)-1),&
          fdvect(idim)%rmat(1-gridpoints(idim):gridpoints(idim)-1),&
          fdvect(idim)%cmat(1-gridpoints(idim):gridpoints(idim)-1),&
          sinepoints(idim)%mat(numpoints(idim),nbox(idim)))
  enddo

  if (griddim.ne.3) then
     OFLWR "griddim.ne.3 not supported no mo"; CFLST
  endif

  if (scalingflag.ne.0) then
     threedtwosize=4
  else
     threedtwosize=1
  endif

  allocate(threed_two(0-gridsize(1):gridsize(1)-1,0-gridsize(2):gridsize(2)-1,0-gridsize(3):gridsize(3)-1,threedtwosize))

end subroutine myprojectalloc



subroutine get_twoe_new(pot)
  use myparams
  use myprojectmod  
  implicit none
  real*8 :: realpot(totpoints)
  DATATYPE :: pot(totpoints)

  if (griddim.ne.3) then
     OFLWR "griddim.ne.3 not supported get_twoe_new"; CFLST
  endif

  call get_3dpoisson(realpot)        
  pot(:)=realpot(:)

  if (notwoflag.eq.1) then
     threed_two(:,:,:,:)=0d0
  endif

  if (debugflag.eq.33.or.scalingflag.ne.0) then
     call get_3dpoisson_scaledoption(pot)
  endif

end subroutine get_twoe_new



subroutine op_yderiv(notused1,notused2)
  use myparams
  implicit none
  DATATYPE :: notused1(*), notused2(*)
  OFLWR "WHAT! no op_yderiv sincdvr, not yet."; CFLST
  notused1(1)=0*notused2(1)
end subroutine op_yderiv
subroutine op_reyderiv(notused1,notused2)
  use myparams
  implicit none
  DATATYPE :: notused1(*), notused2(*)
  OFLWR "WHAT! no op_yderiv sincdvr, not yet."; CFLST
  notused1(1)=0*notused2(1)
end subroutine op_reyderiv
subroutine op_imyderiv(notused1,notused2)
  use myparams
  implicit none
  DATATYPE :: notused1(*), notused2(*)
  OFLWR "WHAT! no op_yderiv sincdvr, not yet."; CFLST
  notused1(1)=0*notused2(1)
end subroutine op_imyderiv



