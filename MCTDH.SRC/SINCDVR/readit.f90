
#include "Definitions.INC"


subroutine get_3dpoisson(pot)
  use tooth
  use myparams
  use myprojectmod
  implicit none
  integer :: ngrid,i,k,myiostat,j,ii,jj,kk,qq(3),pp(3),qbox
  real*8 :: inspacing,rval,pot(totpoints)
  real*8, allocatable :: tempcoulomb(:,:,:)
  real*8, allocatable :: xcoulomb(:,:,:,:,:,:)
  real*8, allocatable :: threed_two_big(:,:,:)
  integer :: ldims(3),udims(3),istart


  if (griddim.ne.3) then
     OFLWR "WTF! should not happen."; CFLST
  endif


!!$    if (runtoothflag.eq.0) then
!!$       if (myrank.eq.1) then
!!$          open(5679,file=invke2file, form="unformatted", status="old",iostat=myiostat)
!!$          if (myiostat.ne.0) then  
!!$             OFLWR "Invke2.Bin file couldn't be read.  ";
!!$             WRFL "   File=", invke2file; CFLST
!!$          endif
!!$          read(5679) ngrid,inspacing  
!!$       endif
!!$       call mympiibcastone(ngrid,1)
!!$       call mympirealbcastone(inspacing,1)
!!$       OFLWR "GO get_3dpoisson. Gridpoints here, on file=",gridpoints(1:3),ngrid*2+1; CFL
!!$       if (abs(inspacing-spacing).gt.1d-7) then
!!$          OFLWR "Spacing does not agree for Invke2.Bin file.  Here, on file=", spacing,inspacing; 
!!$          WRFL "   File = ", invke2file; CFLST
!!$       endif
!!$    else

     ngrid=toothnsmall
     inspacing=spacing

     OFLWR "GO GET TWOCOULOMB. calling tooth.";CFL

!!$    endif
  
  allocate(tempcoulomb(-ngrid:ngrid,-ngrid:ngrid,-ngrid:ngrid));   tempcoulomb(:,:,:)=0d0

!!$    if (runtoothflag.ne.0) then

     call getinverse(tempcoulomb,ngrid,toothnbig,toothnsmall,spacing)
     tempcoulomb(:,:,:)=tempcoulomb(:,:,:)/spacing**3

!!$    else
!!$       if (myrank.eq.1) then
!!$          do i=-ngrid,ngrid
!!$             read(5679) tempcoulomb(:,:,i)
!!$          enddo
!!$          close(5679)
!!$       endif
!!$       call mympirealbcast(tempcoulomb(:,:,:),1,(2*ngrid+1)**3)
!!$    endif


  istart=0
  if (myrank.eq.1.or.(.not.localflag)) then
     istart=1
  endif

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
     
  ii=min(ngrid,gridpoints(1)-1)
  jj=min(ngrid,gridpoints(2)-1)
  kk=min(ngrid,gridpoints(3)-1)

  threed_two_big(-ii:ii,-jj:jj,-kk:kk)=tempcoulomb(-ii:ii,-jj:jj,-kk:kk)

  deallocate(tempcoulomb)

  if (localflag) then
     pp(1:2)=1-gridpoints(1:2)
     qq(1:2)=gridpoints(1:2)-1
     pp(3)=istart+2*gridoffset-gridpoints(3)
     qq(3)=2*(gridoffset+numpoints(3))-gridpoints(3)-1

     threed_two(:,:,:)=0
     threed_two(1-numpoints(1):numpoints(1)-1,1-numpoints(2):numpoints(2)-1,istart-numpoints(3):numpoints(3)-1)=&
          threed_two_big(pp(1):qq(1),pp(2):qq(2),pp(3):qq(3))
  else
     threed_two(:,:,:)=0
     threed_two(1-gridpoints(1):gridpoints(1)-1,1-gridpoints(2):gridpoints(2)-1,1-gridpoints(3):gridpoints(3)-1)=&
          threed_two_big(&
          1-gridpoints(1):gridpoints(1)-1,&
          1-gridpoints(2):gridpoints(2)-1,&
          1-gridpoints(3):gridpoints(3)-1)
  endif

  if (orbparflag.and.myrank.gt.nbox(3)) then
     OFLWR "DOODODFODxxx",myrank,nbox(3); CFLST
  endif

  OFLWR "     ...getting potential..."; CFL

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


