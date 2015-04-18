
#include "Definitions.INC"

module myprojectmod
  implicit none

  DATATYPE, allocatable ::  oned_two(:),twod_two(:,:),threed_two(:,:,:)

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
  type(twomat), allocatable :: littlepot(:), sinepoints(:)
  type(onemat),allocatable :: kevect(:),fdvect(:)

  DATATYPE, allocatable :: dipoles(:,:)

end module myprojectmod


subroutine myprojectalloc()
  use myparams
  use myprojectmod
  implicit none
  integer :: idim

  allocate(dipoles(totpoints,griddim))

  allocate(ketot(griddim),littlepot(griddim),sinepoints(griddim),kevect(griddim),fdtot(griddim),fdvect(griddim))
  do idim=1,griddim
     allocate( &
          fdtot(idim)%mat(numpoints(idim),nbox(idim),numpoints(idim),nbox(idim)), &
          ketot(idim)%mat(numpoints(idim),nbox(idim),numpoints(idim),nbox(idim)), &
          littlepot(idim)%mat(numpoints(idim),nbox(idim)), &
          kevect(idim)%rmat(1-gridpoints(idim):gridpoints(idim)-1),&
          kevect(idim)%cmat(1-gridpoints(idim):gridpoints(idim)-1),&
          fdvect(idim)%rmat(1-gridpoints(idim):gridpoints(idim)-1),&
          fdvect(idim)%cmat(1-gridpoints(idim):gridpoints(idim)-1),&
          sinepoints(idim)%mat(numpoints(idim),nbox(idim)))
  enddo

  if (griddim.ne.3) then
     OFLWR "griddim.ne.3 not supported no mo"; CFLST
  endif

  allocate(threed_two(0-gridsize(1):gridsize(1)-1,0-gridsize(2):gridsize(2)-1,0-gridsize(3):gridsize(3)-1))

end subroutine myprojectalloc



subroutine get_twoe_new(pot)
  use myparams
  use myprojectmod  
  implicit none
  real*8 :: realpot(totpoints)
  DATATYPE :: pot(totpoints)

  if (griddim.ne.3.or.coulflag.eq.0) then
     OFLWR "griddim.ne.3 or coulflag.eq.0 not supported get_twoe_new"; CFLST
  endif

  call get_3dpoisson(realpot)        

  if (toepflag.ne.0) then
     call setblock(nprocs,myrank,gridpoints(:)*2) !! cube only
  endif

  pot(:)=realpot(:)
        
  if (notwoflag.eq.1) then
     threed_two(:,:,:)=0d0
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



