
#include "Definitions.INC"


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





subroutine get_3dpoisson_scaledoption(cpot)
  use tooth
  use myparams
  use pmpimod
  use pfileptrmod
  use myprojectmod
  implicit none
  DATATYPE,intent(out) :: cpot(numpoints(1),numpoints(2),numpoints(3))
  integer :: i,jj,qq(3)
  DATATYPE,allocatable :: sourceterm(:,:,:)
  integer :: null1,null2,null3,null4,nullft(10)

  OFLWR "    ... go scaled option ..."; CFL

  allocate(sourceterm(numpoints(1),numpoints(2),numpoints(3)))

  sourceterm(:,:,:)=0d0; cpot(:,:,:)=0d0

  do i=1,numcenters

     qq(:) = ( gridpoints(1:3)+centershift(:,i)+1 ) / 2

     do jj=1,3
        if (qq(jj).lt.1.or.qq(jj).gt.gridpoints(jj)) then
           OFLWR "error, can't used scaledoption if nuclei are off grid",&
                qq(1:3),gridpoints(1:3); CFLST
        endif
     enddo

     qq(1:3) = qq(1:3) - numpoints(1:3)*(qbox(1:3)-1)
     
     if (qq(1).ge.1.and.qq(1).le.numpoints(1) .and. &
          qq(2).ge.1.and.qq(2).le.numpoints(2) .and. &
          qq(3).ge.1.and.qq(3).le.numpoints(3) ) then
        sourceterm(qq(1),qq(2),qq(3)) = sourceterm(qq(1),qq(2),qq(3)) + nuccharges(i) * (-1)
     endif
  enddo

  call op_tinv(sourceterm,cpot,1,1,null1,null2,null3,null4,nullft)

  deallocate(sourceterm)

  OFLWR "    ... done scaledoption."; CFL

end subroutine get_3dpoisson_scaledoption


