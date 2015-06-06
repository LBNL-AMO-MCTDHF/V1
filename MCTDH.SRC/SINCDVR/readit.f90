
#include "Definitions.INC"


subroutine get_3dpoisson(pot)
  use tooth
  use myparams
  use myprojectmod
  implicit none
  integer :: i,k,j,ii,jj,kk,qq(3),pp(3),qbox
  real*8 :: rval,pot(totpoints)
  real*8, allocatable :: tempcoulomb(:,:,:,:)
  real*8, allocatable :: xcoulomb(:,:,:,:,:,:)
  real*8, allocatable :: threed_two_big(:,:,:)
  integer :: ldims(3),udims(3),istart,gridoffset

  if (griddim.ne.3) then
     OFLWR "WTF! should not happen."; CFLST
  endif

  OFLWR "GO GET_3DPOISSON. calling getinverse.",threedtwosize;CFL

  allocate(tempcoulomb(1-gridpoints(1):gridpoints(1)-1,&
       1-gridpoints(2):gridpoints(2)-1,&
       1-gridpoints(3):gridpoints(3)-1,      threedtwosize))

  if (gridpoints(1).ne.gridpoints(2).or.gridpoints(1).ne.gridpoints(3)) then
     OFLWR" NEED CUBE"; CFLST
  endif

!!$  if (scalingflag.ne.0) then
!!$     call getinverse(tempcoulomb,1,gridpoints(1)-1,toothnbig,toothnsmall,spacing)
!!$!!$     call getinverse(tempcoulomb,2,gridpoints(1)-1,toothnbig,toothnsmall,spacing)
!!$  else
     call getinverse(tempcoulomb,0,gridpoints(1)-1,toothnbig,toothnsmall,spacing)
!!$  endif

  tempcoulomb(:,:,:,:)=tempcoulomb(:,:,:,:)/spacing**3

  istart=0
  if (myrank.eq.1.or.(.not.orbparflag)) then
     istart=1
  endif

  if (orbparflag) then
     gridoffset=(myrank-1)*numpoints(3)

     pp(1:2)=1-gridpoints(1:2)
     qq(1:2)=gridpoints(1:2)-1
     pp(3)=istart+2*gridoffset-gridpoints(3)
     qq(3)=2*(gridoffset+numpoints(3))-gridpoints(3)-1

     threed_two(:,:,:,:)=0
     threed_two(1-numpoints(1):numpoints(1)-1,1-numpoints(2):numpoints(2)-1,istart-numpoints(3):numpoints(3)-1,:)=&
          tempcoulomb(pp(1):qq(1),pp(2):qq(2),pp(3):qq(3),:)
  else
     threed_two(:,:,:,:)=0
     threed_two(1-gridpoints(1):gridpoints(1)-1,1-gridpoints(2):gridpoints(2)-1,1-gridpoints(3):gridpoints(3)-1,:)=&
          tempcoulomb(&
          1-gridpoints(1):gridpoints(1)-1,&
          1-gridpoints(2):gridpoints(2)-1,&
          1-gridpoints(3):gridpoints(3)-1,:)
  endif

  if (orbparflag.and.myrank.gt.nbox(3)) then
     OFLWR "DOODODFODxxx",myrank,nbox(3); CFLST
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

  threed_two_big(-ii:ii,-jj:jj,-kk:kk)=tempcoulomb(-ii:ii,-jj:jj,-kk:kk,1)

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

  pot(:)=RESHAPE(xcoulomb(:,1,:,1,:,qbox(3)),(/totpoints/))

  deallocate(xcoulomb,threed_two_big)

  OFLWR "    ... done get_3dpoisson."; CFL

end subroutine get_3dpoisson





subroutine get_3dpoisson_scaledoption(cpot)
  use tooth
  use myparams
  use myprojectmod
  implicit none
  integer :: i,jj,qq(3)
  DATATYPE :: cpot(numpoints(1),numpoints(2),numpoints(3))
  DATATYPE :: sourceterm(numpoints(1),numpoints(2),numpoints(3))
  integer :: null1,null2,null3,null4,nullft(10)

  OFLWR "    ... go scaled option ..."; CFL

  sourceterm(:,:,:)=0d0; cpot(:,:,:)=0d0

  do i=1,numcenters

     qq(:) = ( gridpoints(1:3)+centershift(:,i)+1 ) / 2

     do jj=1,3
        if (qq(jj).lt.1.or.qq(jj).gt.gridpoints(jj)) then
           OFLWR "error, can't used scaledoption if nuclei are off grid",&
                qq(1:3),gridpoints(1:3); CFLST
        endif
     enddo
     qq(3) = qq(3) - numpoints(3)*(myrank-1)
     
     if (qq(3).ge.1.and.qq(3).le.numpoints(3)) then
        sourceterm(qq(1),qq(2),qq(3)) = sourceterm(qq(1),qq(2),qq(3)) + nuccharges(i) * (-1)
     endif
  enddo

  call op_tinv(sourceterm,cpot,1,1,null1,null2,null3,null4,nullft)

  OFLWR "    ... done scaledoption."; CFL

end subroutine get_3dpoisson_scaledoption


