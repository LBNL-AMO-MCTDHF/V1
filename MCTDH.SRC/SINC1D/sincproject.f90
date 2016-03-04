
#include "Definitions.INC"

module myprojectmod
  implicit none

  DATATYPE, allocatable ::  threed_two(:)

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

  type(fourmat) :: ketot,fdtot
  type(twomat) ::  sinepoints
  type(onemat) :: kevect,fdvect

  DATATYPE, allocatable :: dipoles(:),&

!! WAS:  e.g. X(x) = x + i scalefunction(x,1)


       jacobian(:),&         !! jacobian(:,1) should only be a function of x, etc.
       invjacobian(:),&
       invsqrtjacobian(:),&
       scalediag(:),&
       invsqrtscaleweights(:),&
       scaleweights13(:), &
       invscaleweights13(:),&
       scaleweights16(:), &
       invscaleweights16(:)

end module myprojectmod


subroutine myprojectalloc()
  use myparams
  use pfileptrmod
  use myprojectmod
  implicit none

  allocate(dipoles(totpoints))
  dipoles=0
  
  if (scalingflag.ne.0) then
     allocate(          jacobian(totpoints),invjacobian(totpoints), &
          invsqrtjacobian(totpoints), &
          scalediag(totpoints),&
          invsqrtscaleweights(totpoints),scaleweights13(totpoints),&
          invscaleweights13(totpoints),scaleweights16(totpoints),&
          invscaleweights16(totpoints))
     jacobian=0; invjacobian=0; invsqrtjacobian=0; scalediag=0;
     invsqrtscaleweights=0;scaleweights13=0;invscaleweights13=0;
     scaleweights16=0;invscaleweights16=0;
  endif
  

  allocate( &

!! Allocating extra here for fdtot%mat and ketot%mat (+1's) --
!!   see Z/GEMM calls in coreproject.f90... leading dimension not
!!   allocated as passed to Z/GEMM without extra

       fdtot%mat(numpoints,nbox,numpoints,nbox   +1), &
       ketot%mat(numpoints,nbox,numpoints,nbox   +1), &
       kevect%rmat(1-gridpoints:gridpoints-1),&
       kevect%cmat(1-gridpoints:gridpoints-1),&
       fdvect%rmat(1-gridpoints:gridpoints-1),&
       fdvect%cmat(1-gridpoints:gridpoints-1),&
       sinepoints%mat(numpoints,nbox))
  fdtot%mat=0; ketot%mat=0; kevect%rmat=0; kevect%cmat=0;
  fdvect%rmat=0; kevect%cmat=0

  allocate(threed_two(0-numpoints:numpoints-1))
  threed_two=0

end subroutine myprojectalloc


module twoemod
  implicit none
  DATATYPE, allocatable :: frozenreduced(:)
end module twoemod


subroutine get_twoe_new(pot)
  use myparams
  use pfileptrmod
  use myprojectmod  
  implicit none
  DATATYPE,intent(out) :: pot(totpoints)

  pot(:)=dipoles(:)**2 * 0.5d0

  if (notwoflag.eq.1) then
     threed_two(:)=0d0
  endif

end subroutine get_twoe_new



subroutine op_yderiv(notint,notused1,notused2)
  use pfileptrmod
  implicit none
  integer :: notint
  DATATYPE :: notused1(*), notused2(*)
  OFLWR "WHAT! no op_yderiv sincdvr, not yet."; CFLST
  notused1(1)=0*notused2(1)
end subroutine op_yderiv

