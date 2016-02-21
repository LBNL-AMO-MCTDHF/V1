

!! eta derivs is (1-eta^2) d/deta - eta (E1) or O1
!! eta_rhoderivs is d/deta
!! adapted from eta_init
subroutine GetJacobiKE(points,weights,ketot, lbig, eta_derivs,eta_rhoderivs,evenodd)
  implicit none
  integer,intent(in) :: lbig, evenodd
  real*8,intent(out) :: points((lbig+1)), weights((lbig+1)),ketot((lbig+1),(lbig+1)),&
       eta_derivs((lbig+1),(lbig+1)),eta_rhoderivs(lbig+1,lbig+1)
  integer, parameter :: numextra=6
  integer :: numpoints,  one=1, two=2, izero=0, i,j,k,jj , extraorder
  real*8 :: extrapoints((lbig+1)+numextra), extraweights(numextra+(lbig+1)), &
       firstdertot((lbig+1)+numextra,(lbig+1)), scratch(numextra+(lbig+1)), &
       etavals((lbig+1)+numextra,(lbig+1)), endpoints(2)=[-1.0d0,1.0d0] , zero = 0.0, sum

  extrapoints=0; extraweights=0; firstdertot=0; scratch=0; etavals=0

  numpoints=lbig+1
  one=1;two=2
  extraorder=numextra+numpoints

  call gaussq(one,numpoints,zero,zero,izero,endpoints,scratch,points,weights)
  call gaussq(one,extraorder,zero,zero,izero,endpoints,scratch,extrapoints,extraweights)
 
 ! firstderiv(i,j)  First derivative matrix.  Deriv of ith dvr at x_j

  do k=1,extraorder  ! point k
     do j=1,numpoints  ! which basis function
        sum=1.d0
        do i=1,numpoints
           if (j/=i) then
              sum=sum * (extrapoints(k)-points(i))/(points(j)-points(i))
           endif
        enddo
        etavals(k,j)=sum
     enddo
  enddo
  
  if (numextra.ne.0) then
     do k=1,extraorder  ! point k
        do j=1,numpoints
           sum=0.d0
           do i=1,numpoints
              if (j/=i) then
                 sum=sum+ etavals(k,j)/(extrapoints(k)-points(i))
              endif
           enddo
           firstdertot(k,j)=sum
        enddo
     enddo
  else
     print *, "Numextra not supported.", numextra;     stop
  endif

  do jj=1,numpoints 
     firstdertot(:,jj) = firstdertot(:,jj) / sqrt(weights(jj)) 
     etavals(:,jj) = etavals(:,jj) / sqrt(weights(jj)) 
  enddo
  
  do i=1,numpoints
     do j=1,numpoints
        sum=0
        do k=1,extraorder
           sum=sum + extraweights(k) * etavals(k,i) * ( extrapoints(k)*(1-extrapoints(k)**2) * &
                firstdertot(k,j) + (0.5d0 - 3d0/2d0 * extrapoints(k)**2) * etavals(k,j) )
        enddo
        eta_rhoderivs(i,j)=sum
        
        sum=0.0d0
        if (evenodd.eq.0) then 
           do k=1,extraorder
              sum= sum + etavals(k,i) &
              * (firstdertot(k,j)* (1-extrapoints(k)**2) - etavals(k,j) * extrapoints(k) )  &
                   * extraweights(k) 
           enddo
        else
           do k=1,extraorder
              sum= sum + etavals(k,i) &
                   * (firstdertot(k,j)* (1-extrapoints(k)**2)**2 - 2*etavals(k,j) * &
                   (1-extrapoints(k)**2)*extrapoints(k) )  &
                   * extraweights(k) 
           enddo
           sum=sum/sqrt((1-points(i)**2)*(1-points(j)**2))
        endif
        eta_derivs(i,j)=sum
     enddo
  enddo

!! negative definite operators for eta.
  
  do j=1,numpoints
     do i=1,numpoints
        if (evenodd.eq.0) then !! even, d/dx (1-x^2)**2 d/dx
           sum=0.0d0
           do k=1,extraorder
              sum = sum - firstdertot(k,j) * firstdertot(k,i) * extraweights(k) * (1-extrapoints(k)**2)
           enddo
        else
           sum=0.0d0
           do k=1,extraorder
              sum= sum - (firstdertot(k,j)* (1-extrapoints(k)**2) - etavals(k,j) * extrapoints(k) )  &
                   * (firstdertot(k,i)* (1-extrapoints(k)**2) - etavals(k,i) * extrapoints(k) )  &
                   * extraweights(k) 
           enddo
           sum=sum/sqrt((1-points(i)**2)*(1-points(j)**2))
        endif
        ketot(i,j)= sum   !! symmetric kinetic energy
     enddo
  enddo

end subroutine GetJacobiKE






