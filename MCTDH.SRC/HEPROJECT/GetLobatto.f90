
!! RADIAL QUADRATURE

#include "Definitions.INC"


subroutine getlobatto(points,weights,points2d,weights2d, ketot, numpoints,&
     numelements,elementsizes, gridpoints, celement, ecstheta,  xi_derivs, xi_rhoderivs, xi_cent, evenodd)
  implicit none
  integer,intent(in) :: numelements,gridpoints,celement, numpoints, evenodd
  real*8, intent(in) :: elementsizes(numelements), ecstheta
  DATAECS,intent(out) :: points(gridpoints),  weights(gridpoints),  xi_derivs(gridpoints,gridpoints), &
       xi_rhoderivs(gridpoints,gridpoints), ketot(gridpoints,gridpoints), xi_cent(gridpoints,gridpoints)
  integer ::  one=1, izero=0, two=2, i,j,k,l,jj, qq, extraorder
  real*8 :: endpoints(2), zero = 0.0,rsum
  integer, parameter :: numextra=11
  real*8 :: extrapoints0(numpoints+numextra), extraweights0(numpoints+numextra),&
       firstder(numpoints+numextra,numpoints),  points2d(numpoints), weights2d(numpoints), &
       scratch(2*numpoints+numextra), xivals0(numpoints+numextra,numpoints)
  DATAECS :: extrapoints(numpoints+numextra,numelements), extraweights(numextra+numpoints,numelements), &
       firstdertot(numpoints+numextra,numelements,gridpoints), &
       xivals(numpoints+numextra,numelements,gridpoints),  cweight, sum, sum1, sum2
  DATAECS :: xi_ovl(gridpoints,gridpoints), xi_coul(gridpoints,gridpoints), ppoints(gridpoints)   ! temp debug? AUTOMATIC
  
  extrapoints0=0; extraweights0=0; firstder=0; points2d=0; weights2d=0;
  scratch=0; xivals0=0; extrapoints=0; extraweights=0; firstdertot=0; xivals=0

  i=celement; sum=ecstheta !! avoid warn unused

  if (evenodd.ne.0.and.evenodd.ne.1) then
     print *, "evenodd error", evenodd
  endif

  extraorder=numextra+numpoints
  if (gridpoints.ne.(numelements*(numpoints-1)+1)) then
     print *, "gridpoints error", gridpoints, numelements, numpoints;     stop
  endif

  !! for exact quad.  points used on [-1,1] first then changed.
  call gaussq(one,extraorder,zero,zero,izero,endpoints,scratch,extrapoints0(:),extraweights0(:))
     
  !! for G-L elements
  endpoints= [ -1.0d0 , 1.d0 ]
  call gaussq(one,numpoints,zero,zero,two,endpoints,scratch,points2d(:),weights2d(:))

  xivals0=0.d0
  do k=1,extraorder  ! point k
     do j=1,numpoints  ! which basis function
        rsum=1
        do i=1,numpoints
           if (j/=i) then
              rsum=rsum * (extrapoints0(k)-points2d(i))/(points2d(j)-points2d(i))
           endif
        enddo
        xivals0(k,j) = rsum
     enddo
  enddo
  if (numextra.ne.0) then
     do k=1,extraorder  ! point k
        do j=1,numpoints
           rsum=0.d0
           do i=1,numpoints
              if (j/=i) then
                 rsum=rsum+ xivals0(k,j)/(extrapoints0(k)-points2d(i))
              endif
           enddo
           firstder(k,j)=rsum
        enddo
     enddo
  else
     print *, "Numextra not supported.", numextra
     stop
  endif

  weights=0.0d0;  j=0;  jj=0;  sum=0.d0;  firstdertot=0.d0;  xivals=0d0

  do l=1,numelements
     cweight=1.0d0
#ifdef ECSFLAG
     if (l.ge.celement) then
        cweight=exp((0.d0, 1.d0)*ecstheta)
     endif
#endif
     qq = (numpoints-1)*(l-1)

     !! for now value of derivative of interp function at point

     firstdertot(:,l,qq+1:qq+numpoints)=firstdertot(:,l,qq+1:qq+numpoints) + &
          firstder(:,:)/cweight*2.d0/elementsizes(l)
     xivals(:,l,qq+1:qq+numpoints)=xivals(:,l,qq+1:qq+numpoints)+xivals0(:,:)
     extraweights(:,l)=extraweights0(:)* elementsizes(l)/2.0d0 * cweight
     extrapoints(:,l)= cweight*elementsizes(l)/2.0d0 * (extrapoints0(:)+1d0) + sum

     do k=1,numpoints
        jj=jj+1
        if ((l==1) .or. (k/=1)) then    
           j = j+1
        endif
        points(j) = cweight*elementsizes(l)/2.0d0 * (points2d(k)+1.0d0) + sum
        weights(j) = weights(j) + weights2d(k) * elementsizes(l)/2.0d0 * cweight
     enddo
     sum=sum+elementsizes(l)*cweight
  enddo

  do jj=1,gridpoints
     firstdertot(:,:,jj) = firstdertot(:,:,jj) / sqrt(weights(jj)) 
     xivals(:,:,jj) = xivals(:,:,jj) / sqrt(weights(jj)) 
  enddo

  do i=1,gridpoints
     do j=1,gridpoints
        sum=0.d0
        if (evenodd.eq.0) then
           do l=1,numelements
              do k=1,extraorder
                 sum=sum + extraweights(k,l) * xivals(k,l,i) * firstdertot(k,l,j)
              enddo
           enddo
        else
           do l=1,numelements
              do k=1,extraorder
                 sum=sum + extraweights(k,l) * xivals(k,l,i) * &
                      (extrapoints(k,l)**2*firstdertot(k,l,j) + extrapoints(k,l)*xivals(k,l,j))
              enddo
           enddo
           sum=sum/(points(i)*points(j))
        endif
        xi_derivs(i,j)=sum
        sum=0
        do l=1,numelements
           do k=1,extraorder
              sum=sum + extraweights(k,l) * xivals(k,l,i) * &
                   ( extrapoints(k,l) * firstdertot(k,l,j) + 0.5d0 * xivals(k,l,j) )
           enddo
        enddo
        xi_rhoderivs(i,j)=sum
     enddo
  enddo

!! negative definite operators for xi.

  do j=1,gridpoints
     do i=1,gridpoints
        if (evenodd.eq.0) then !! even, d/dx  d/dx
           sum=0.0d0
           do l=1,numelements
              do k=1,extraorder
                 sum  = sum  - firstdertot(k,l,j) * firstdertot(k,l,i) * extraweights(k,l) 
              enddo
           enddo
        else
           sum=0.0d0
           do l=1,numelements
              do k=1,extraorder   
                 sum= sum - (firstdertot(k,l,j)*extrapoints(k,l) + xivals(k,l,j) )  &
                      * (firstdertot(k,l,i)*extrapoints(k,l) + xivals(k,l,i)  )  &
                      * extraweights(k,l) 
              enddo
           enddo
           sum=sum/(points(i)*points(j))
        endif
        ketot(i,j)=sum   !! symmetric kinetic energy
     enddo
  enddo

  ! DVR expression for centrifugal potential saved here from init_HE_new.f90
  !
  !     do i=1,hegridpoints-2
  !        do j=1,hegridpoints-2
  !           bigham(i,:,j,:)=bigham(i,:,j,:) + jacobike(:,:,abs(imvalue)) * &
  !                glfirstdertot(1,i+1,0) *  glfirstdertot(1,j+1,0)    
  !        enddo
  !     enddo
  !
  !! matrix elements for centrifugal operator...
  !! for acceleration operator (dipole acceleration) and centrifugal potential in hamiltonian
  !! prior attempt above (first derivative term addition to 1/r^2 according to DVR principles)
  !!    was not successful for hamiltonian
  !! current dvr expression for 1/r^2 works well for hamiltonian,
  !!    very accurate but not variational (is like collocation)
  !! using xi_cent as representation of 1/r^2 in hamiltonian leads to more variational
  !!    hydrogenic eigenvalues (upper bounds) 
  !!    it is not exactly variational because overlap matrix is approximated
  !!    but keep in mind, with one finite element the DVR approx to coulomb operator 1/r is exact

  do i=1,gridpoints
     do j=1,gridpoints
        sum=0.d0
        sum1=0.d0
        sum2=0.d0

        do l=1,numelements
           do k=1,extraorder
              sum2=sum2 + extraweights(k,l) * xivals(k,l,i) * xivals(k,l,j) / (extrapoints(k,l))**2 
              sum1=sum1 + extraweights(k,l) * xivals(k,l,i) * xivals(k,l,j) / (extrapoints(k,l))
              sum=sum + extraweights(k,l) * xivals(k,l,i) * xivals(k,l,j)
           enddo
        enddo
        
        xi_cent(i,j)=sum2 
        xi_coul(i,j)=sum1 
        xi_ovl(i,j)=sum
     enddo

     ! print '(13F10.3)',  real(xi_cent(i,1:13) * points(i) * points(1:13),8)
     ! print '(1F11.5)',  real( xi_cent(i,i) * points(i)**2, 8)
     ! print '(1F11.5)',  real( xi_cent(i,i) * points(i), 8)
     
  enddo

  if (1==0) then
     ! temp output for debug
     ppoints(:) = points(:)
     ppoints(1) = 1;  
     print *, ' '
     print *, ' ovl '
     do i=1,13
        print '(13F10.3)',  real(xi_ovl(i,1:13),8)
     enddo
     print *, ' coul '
     do i=1,13
        print '(13F10.3)',  real(xi_coul(i,1:13)*sqrt(ppoints(i)*ppoints(1:13)),8)
     enddo
     print *, ' cent '
     do i=1,13
        print '(13F10.3)',  real(xi_cent(i,1:13)*ppoints(i)*ppoints(1:13),8)
     enddo
     print *, ' '
  endif
  
end subroutine getlobatto


