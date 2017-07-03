
#include "Definitions.INC"

module myprojectmod
  implicit none

  DATATYPE, allocatable ::  threed_two(:,:,:)

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
  dipoles=0
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
        maskfunction(idim)%rmat=0
     enddo
  endif
  
  if (scalingflag.ne.0) then
     allocate(          jacobian(totpoints,3),invjacobian(totpoints,3), &
          invsqrtjacobian(totpoints,3), &
          scalediag(totpoints),&
          invsqrtscaleweights(totpoints),scaleweights13(totpoints),&
          invscaleweights13(totpoints),scaleweights16(totpoints),&
          invscaleweights16(totpoints))
     jacobian=0; invjacobian=0; invsqrtjacobian=0; scalediag=0;
     invsqrtscaleweights=0;scaleweights13=0;invscaleweights13=0;
     scaleweights16=0;invscaleweights16=0;
  endif
  
  allocate(ketot(griddim),sinepoints(griddim),kevect(griddim),fdtot(griddim),fdvect(griddim))

  do idim=1,griddim
     allocate( &

!! not needed with tam.  old comment:
!! !! Allocating extra here for fdtot%mat and ketot%mat (+1's) --
!! !!   see Z/GEMM calls in coreproject.f90... leading dimension not
!! !!   allocated as passed to Z/GEMM without extra

!!          fdtot(idim)%mat(numpoints(idim),nbox(idim),numpoints(idim),nbox(idim)   +1), &
!!          ketot(idim)%mat(numpoints(idim),nbox(idim),numpoints(idim),nbox(idim)   +1), &
          fdtot(idim)%mat(numpoints(idim),nbox(idim),numpoints(idim),nbox(idim)), &
          ketot(idim)%mat(numpoints(idim),nbox(idim),numpoints(idim),nbox(idim)), &
          fdtot(idim)%tam(numpoints(idim),numpoints(idim),nbox(idim),nbox(idim)), &
          ketot(idim)%tam(numpoints(idim),numpoints(idim),nbox(idim),nbox(idim)), &
          kevect(idim)%rmat(1-gridpoints(idim):gridpoints(idim)-1),&
          kevect(idim)%cmat(1-gridpoints(idim):gridpoints(idim)-1),&
          fdvect(idim)%rmat(1-gridpoints(idim):gridpoints(idim)-1),&
          fdvect(idim)%cmat(1-gridpoints(idim):gridpoints(idim)-1),&
          sinepoints(idim)%mat(numpoints(idim),nbox(idim)))
     fdtot(idim)%mat=0; ketot(idim)%mat=0; kevect(idim)%rmat=0; kevect(idim)%cmat=0;
     fdvect(idim)%rmat=0; kevect(idim)%cmat=0; fdtot(idim)%tam=0; ketot(idim)%tam=0
  enddo

  if (griddim.ne.3) then
     OFLWR "griddim.ne.3 not supported no mo"; CFLST
  endif

  allocate(threed_two(0-numpoints(1):numpoints(1)-1,0-numpoints(2):numpoints(2)-1,0-numpoints(3):numpoints(3)-1))
  threed_two=0

end subroutine myprojectalloc


module poissonmod
contains

subroutine get_3dpoisson(pot)
  use tooth
  use myparams
  use pmpimod
  use pfileptrmod
  use myprojectmod
  use mpisubmod   !! IN PARENT DIRECTORY
  implicit none
  real*8,intent(out) :: pot(totpoints)
  real*8 :: rval
  real*8, allocatable :: tempcoulomb(:,:,:)
  real*8, allocatable :: xcoulomb(:,:,:,:,:,:)
  real*8, allocatable :: threed_two_big(:,:,:)
  integer :: ldims(3),udims(3),istart(3),gridoffset(3)
  integer :: i,k,j,ii,jj,kk,qq(3),pp(3)

  if (griddim.ne.3) then
     OFLWR "WTF! should not happen."; CFLST
  endif

  OFLWR "GO GET_3DPOISSON. calling getinverse.";CFL

  allocate(tempcoulomb(1-gridpoints(1):gridpoints(1)-1,&
       1-gridpoints(2):gridpoints(2)-1,&
       1-gridpoints(3):gridpoints(3)-1))

  tempcoulomb=0

  if (gridpoints(1).ne.gridpoints(2).or.gridpoints(1).ne.gridpoints(3)) then
     OFLWR" NEED CUBE"; CFLST
  endif

  if (toothnbig.lt.gridpoints(1)-1.or.toothnbig.lt.gridpoints(2)-1.or.toothnbig.lt.gridpoints(3)-1) then
     OFLWR "Error, need toothnbig > gridpoints-1",toothnbig, gridpoints(1:3); CFLST
  endif
  if (myrank.eq.1) then
     call getinverse(tempcoulomb,gridpoints(1)-1,toothnbig,toothnsmall,spacing)
     tempcoulomb(:,:,:)=tempcoulomb(:,:,:)/spacing**3
  endif
  call mympirealbcast(tempcoulomb,1,(2*gridpoints(1)-1)*(2*gridpoints(2)-1)*(2*gridpoints(3)-1))

  istart(:)=1
  if (orbparflag) then
     do ii=orbparlevel,3
        if (boxrank(ii).ne.1) then
           istart(ii)=0
        endif
     enddo
  endif

  if (orbparflag) then

     gridoffset(1:3)=(boxrank(1:3)-1)*numpoints(1:3)

     pp(1:3)=istart(1:3)+2*gridoffset(1:3)-gridpoints(1:3)
     qq(1:3)=2*(gridoffset(1:3)+numpoints(1:3))-gridpoints(1:3)-1

     threed_two(:,:,:)=0
     threed_two(istart(1)-numpoints(1):numpoints(1)-1,&
          istart(2)-numpoints(2):numpoints(2)-1,&
          istart(3)-numpoints(3):numpoints(3)-1)=&
          tempcoulomb(pp(1):qq(1),pp(2):qq(2),pp(3):qq(3))
  else
     threed_two(:,:,:)=0
     threed_two(1-gridpoints(1):gridpoints(1)-1,1-gridpoints(2):gridpoints(2)-1,1-gridpoints(3):gridpoints(3)-1)=&
          tempcoulomb(&
          1-gridpoints(1):gridpoints(1)-1,&
          1-gridpoints(2):gridpoints(2)-1,&
          1-gridpoints(3):gridpoints(3)-1)
  endif

  OFLWR "     ...getting potential..."; CFL

  ldims(1:3)=1-gridpoints(1:3)
  udims(1:3)=gridpoints(1:3)-1

!! centershift entered times two!!  will be divided by two...
  do i=1,numcenters
     pp(:) = -gridpoints(1:3)+centershift(:,i)+1 
     qq(:) =  gridpoints(1:3)+centershift(:,i)-1   
     do ii=1,3
        if (mod(qq(ii),2).eq.1.or.mod(pp(ii),2).eq.1) then
           OFLWR "WTF even!!!",i,ii,gridpoints(ii),centershift(ii,i); 
           WRFL pp(:),qq(:); CFLST
        endif
     enddo
     pp(:)=pp(:)/2; qq(:)=qq(:)/2

     ldims(:)=min(ldims(:),pp(:))
     udims(:)=max(udims(:),qq(:))
  enddo

  allocate(threed_two_big(ldims(1):udims(1),ldims(2):udims(2),ldims(3):udims(3)))

  threed_two_big(:,:,:)=0d0

  OFLWR "filling in...."; CFL

  do k=ldims(3),udims(3)
     do j=ldims(2),udims(2)
        do i=ldims(1),udims(1)
           rval=spacing*sqrt(real(i**2+j**2+k**2))
           threed_two_big(i,j,k)=1d0/rval
        enddo
     enddo
  enddo
  
  OFLWR "    ...filled."; CFL
     
  ii=gridpoints(1)-1
  jj=gridpoints(2)-1
  kk=gridpoints(3)-1

  threed_two_big(-ii:ii,-jj:jj,-kk:kk)=tempcoulomb(-ii:ii,-jj:jj,-kk:kk)

  deallocate(tempcoulomb)

  allocate(xcoulomb(numpoints(1),nbox(1),numpoints(2),nbox(2),numpoints(3),nbox(3)))
  xcoulomb=0d0

  do i=1,numcenters

     pp(:) = -gridpoints(1:3)+centershift(:,i)+1 
     qq(:) =  gridpoints(1:3)+centershift(:,i)-1   
     pp(:)=pp(:)/2; qq(:)=qq(:)/2

     xcoulomb(:,:,:,:,:,:)=xcoulomb(:,:,:,:,:,:)   -1d0 * nuccharges(i) *RESHAPE(&
          threed_two_big(pp(1):qq(1), pp(2):qq(2), pp(3):qq(3)), &
          (/numpoints(1),nbox(1),numpoints(2),nbox(2),numpoints(3),nbox(3)/))
  enddo

  pot(:)=RESHAPE(xcoulomb(:,qbox(1),:,qbox(2),:,qbox(3)),(/totpoints/))

  deallocate(xcoulomb,threed_two_big)

  OFLWR "    ... done get_3dpoisson."; CFL

end subroutine get_3dpoisson


!!$TINV subroutine get_3dpoisson_scaledoption(cpot)
!!$TINV   use tooth
!!$TINV   use myparams
!!$TINV   use pmpimod
!!$TINV   use pfileptrmod
!!$TINV   use myprojectmod
!!$TINV   use tinvsubmod
!!$TINV   implicit none
!!$TINV   DATATYPE,intent(out) :: cpot(numpoints(1),numpoints(2),numpoints(3))
!!$TINV   integer :: i,jj,qq(3)
!!$TINV   DATATYPE,allocatable :: sourceterm(:,:,:)
!!$TINV   integer :: null1,null2,null3,null4,nullft(10)
!!$TINV 
!!$TINV   OFLWR "    ... go scaled option ..."; CFL
!!$TINV 
!!$TINV   allocate(sourceterm(numpoints(1),numpoints(2),numpoints(3)))
!!$TINV 
!!$TINV   sourceterm(:,:,:)=0d0; cpot(:,:,:)=0d0
!!$TINV 
!!$TINV   do i=1,numcenters
!!$TINV 
!!$TINV      qq(:) = ( gridpoints(1:3)+centershift(:,i)+1 ) / 2
!!$TINV 
!!$TINV      do jj=1,3
!!$TINV         if (qq(jj).lt.1.or.qq(jj).gt.gridpoints(jj)) then
!!$TINV            OFLWR "error, can't used scaledoption if nuclei are off grid",&
!!$TINV                 qq(1:3),gridpoints(1:3); CFLST
!!$TINV         endif
!!$TINV      enddo
!!$TINV 
!!$TINV      qq(1:3) = qq(1:3) - numpoints(1:3)*(qbox(1:3)-1)
!!$TINV      
!!$TINV      if (qq(1).ge.1.and.qq(1).le.numpoints(1) .and. &
!!$TINV           qq(2).ge.1.and.qq(2).le.numpoints(2) .and. &
!!$TINV           qq(3).ge.1.and.qq(3).le.numpoints(3) ) then
!!$TINV         sourceterm(qq(1),qq(2),qq(3)) = sourceterm(qq(1),qq(2),qq(3)) + nuccharges(i) * (-1)
!!$TINV      endif
!!$TINV   enddo
!!$TINV 
!!$TINV   call op_tinv(sourceterm,cpot,1,1,null1,null2,null3,null4,nullft)
!!$TINV 
!!$TINV   deallocate(sourceterm)
!!$TINV 
!!$TINV   OFLWR "    ... done scaledoption."; CFL
!!$TINV 
!!$TINV end subroutine get_3dpoisson_scaledoption

end module poissonmod

module gettwoemod
contains

subroutine get_twoe_new(pot)
  use myparams
  use pfileptrmod
  use myprojectmod  
  use poissonmod
  implicit none
  DATATYPE,intent(out) :: pot(totpoints)
  real*8,allocatable :: realpot(:)

  if (griddim.ne.3) then
     OFLWR "griddim.ne.3 not supported get_twoe_new"; CFLST
  endif

  allocate(realpot(totpoints))
  realpot=0
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

end module


subroutine op_yderiv(notint,notused1,notused2)
  use pfileptrmod
  implicit none
  integer :: notint
  DATATYPE :: notused1(*), notused2(*)
  OFLWR "WHAT! no op_yderiv sincdvr, not yet."; CFLST
  notused1(1)=0*notused2(1)
end subroutine op_yderiv

