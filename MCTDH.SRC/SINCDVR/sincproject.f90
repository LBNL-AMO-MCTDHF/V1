
#include "Definitions.INC"

module myprojectmod
  implicit none

  DATATYPE, allocatable ::  threed_two(:,:,:)

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

  type(onemat), allocatable :: maskfunction(:)

  DATATYPE, allocatable :: dipoles(:,:),&

!! WAS:  e.g. X(x) = x + i scalefunction(x,1)


       jacobian(:,:),&         !! jacobian(:,1) should only be a function of x, etc.
       invjacobian(:,:),&
       invsqrtjacobian(:,:),&
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
  integer :: idim

  allocate(dipoles(totpoints,griddim))
  if (maskflag.ne.0) then
     if (masknumpoints.lt.0.or.&
          masknumpoints.gt.gridpoints(1).or.&
          masknumpoints.gt.gridpoints(2).or.&
          masknumpoints.gt.gridpoints(3)) then
        OFLWR "masknumpoints not allowed",masknumpoints; CFLST
     endif
     allocate(maskfunction(3))
     do idim=1,griddim
        allocate(maskfunction(idim)%rmat(numpoints(idim)))
     enddo
  endif
  
  if (scalingflag.ne.0) then
     allocate(          jacobian(totpoints,3),invjacobian(totpoints,3), &
          invsqrtjacobian(totpoints,3), &
          scalediag(totpoints),&
          invsqrtscaleweights(totpoints),scaleweights13(totpoints),&
          invscaleweights13(totpoints),scaleweights16(totpoints),&
          invscaleweights16(totpoints))

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

  allocate(threed_two(0-numpoints(1):numpoints(1)-1,0-numpoints(2):numpoints(2)-1,0-numpoints(3):numpoints(3)-1))

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
  real*8,allocatable :: realpot(:)

  if (griddim.ne.3) then
     OFLWR "griddim.ne.3 not supported get_twoe_new"; CFLST
  endif

  allocate(realpot(totpoints))
  call get_3dpoisson(realpot)        
  pot(:)=realpot(:)
  deallocate(realpot)

  if (notwoflag.eq.1) then
     threed_two(:,:,:)=0d0
  endif

!!$TINV  if (debugflag.eq.33.or.scalingflag.ne.0) then
!!$TINV     call get_3dpoisson_scaledoption(pot)
!!$TINV   endif

end subroutine get_twoe_new



subroutine op_yderiv(notint,notused1,notused2)
  use pfileptrmod
  implicit none
  integer :: notint
  DATATYPE :: notused1(*), notused2(*)
  OFLWR "WHAT! no op_yderiv sincdvr, not yet."; CFLST
  notused1(1)=0*notused2(1)
end subroutine op_yderiv

