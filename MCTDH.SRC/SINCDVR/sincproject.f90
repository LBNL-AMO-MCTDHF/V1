
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

  type(onemat), allocatable :: maskfunction(:)

  DATATYPE, allocatable :: dipoles(:,:),&

!! WAS:  e.g. X(x) = x + i scalefunction(x,1)


       jacobian(:,:),&         !! jacobian(:,1) should only be a function of x, etc.
       invjacobian(:,:),&
       scalediag(:),&      !!sum  (-1/4) J^-4 (d/dx J(x))^2 = (-1)* sum_i=1..3 scaleder(:,i)**2
       scaleder(:,:),&       !!     (1/2) J^-2 (d/dx J(x))   
       invsqrtscaleweights(:)

end module myprojectmod


subroutine myprojectalloc()
  use myparams
  use fileptrmod
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
          scalediag(totpoints),scaleder(totpoints,3),&
          invsqrtscaleweights(totpoints))

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

!!$  if (scalingflag.ne.0) then
!!$!!$     threedtwosize=7
!!$     threedtwosize=4
!!$  else

     threedtwosize=1

!!$  endif

  allocate(threed_two(0-numpoints(1):numpoints(1)-1,0-numpoints(2):numpoints(2)-1,0-numpoints(3):numpoints(3)-1,threedtwosize))

end subroutine myprojectalloc



subroutine get_twoe_new(pot)
  use myparams
  use fileptrmod
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
  use fileptrmod
  implicit none
  DATATYPE :: notused1(*), notused2(*)
  OFLWR "WHAT! no op_yderiv sincdvr, not yet."; CFLST
  notused1(1)=0*notused2(1)
end subroutine op_yderiv
subroutine op_reyderiv(notused1,notused2)
  use fileptrmod
  implicit none
  DATATYPE :: notused1(*), notused2(*)
  OFLWR "WHAT! no op_yderiv sincdvr, not yet."; CFLST
  notused1(1)=0*notused2(1)
end subroutine op_reyderiv
subroutine op_imyderiv(notused1,notused2)
  use fileptrmod
  implicit none
  DATATYPE :: notused1(*), notused2(*)
  OFLWR "WHAT! no op_yderiv sincdvr, not yet."; CFLST
  notused1(1)=0*notused2(1)
end subroutine op_imyderiv



