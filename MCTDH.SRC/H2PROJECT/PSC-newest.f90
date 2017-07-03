
#include "Definitions.INC"

module pscmod
contains

subroutine PSC()
  use myparams
  use myprojectmod  
  use myprojecttempmod  
  use lpmsubmod
  implicit none
  integer :: i,k,m,iii,  j, thismval
  real*8 :: eta1,eta2,xi1,xi2, dummypoints(rgridpoints), dummyweights(rgridpoints) 

  dummypoints=0; dummyweights=0

  call PSCtempalloc()
  call r_init(dummypoints,dummyweights, usualke, rnumpoints,rnumelements,relementsize,rgridpoints,rstart)
  usualke=usualke*(-0.5d0)/Rmass

  call prolate_init(rpoints,rweights, rnumpoints,rnumelements,relementsize,rgridpoints,rstart, &
       prolate_derivs, rketot, rcelement, rthetaecs)
  rketot=rketot*(-0.5d0/Rmass)

  if (capflag.eq.1) then
     OFLWR "Doing CAP. "; CFL
     do i=1,rgridpoints
        
        if (abs(imag(rpoints(i)+(0.d0,0.d0))).gt.1.d-10) then
           OFLWR "Error, looks like you have scaling on with cap. "
           write(mpifileptr,*) i,rpoints(i);  CFLST
        endif
        
        if (real(rpoints(i)).gt.capstart) then
           rketot(i,i) = rketot(i,i) + (0.d0,-1.d0) * capstrength*(rpoints(i)-capstart)**cappower
        endif
     enddo
  endif

!!!      ELECTRONIC      !!!
  
  call xi_init(xipoints,xiweights,xipoints2d,xike(:,:,0),xinumpoints,xinumelements,xielementsizes,&
       xigridpoints,xicelement,xiecstheta, lobattoopt, xideg1(:,:,0),xideg2(:,:,0),xi_derivs(:,:,0),&
       xi_rhoderivs(:,:),xi_fourth(:,:,0), xi_third(:,:,0), xi_second(:,:,0), 0)
  call eta_init(etapoints,etaweights,etake(:,:,0),numeta, etadeg1(:,:,0),etadeg2(:,:,0),&
       eta_derivs(:,:,0),eta_rhoderivs(:,:),eta_fourth(:,:,0), eta_third(:,:,0),eta_second(:,:,0),0)
  
  call xi_init(xipoints,xiweights,xipoints2d,xike(:,:,1),xinumpoints,xinumelements,xielementsizes,&
       xigridpoints,xicelement,xiecstheta, lobattoopt, xideg1(:,:,1),xideg2(:,:,1),xi_derivs(:,:,1),&
       xi_rhoderivs(:,:),xi_fourth(:,:,1), xi_third(:,:,1), xi_second(:,:,1), 1)
  call eta_init(etapoints,etaweights,etake(:,:,1),numeta, etadeg1(:,:,1),etadeg2(:,:,1),&
       eta_derivs(:,:,1),eta_rhoderivs(:,:),eta_fourth(:,:,1), eta_third(:,:,1),eta_second(:,:,1),1)

! valid for odd because of cancellation of inverse rho squared; valid for M=0 due to zero;
!    valid for higher, no l'hospital term.

  do thismval=mbig+1,0,-1    !!!mbig+1
     etake(:,:,thismval)=etake(:,:,mod(thismval,2))
     if (thismval.ne.0) then
        do i=1,numeta
           etake(i,i,thismval) = etake(i,i,thismval) - thismval**2/(1-etapoints(i)**2)   
        enddo
     endif
  enddo

  do thismval=mseriesmax,0,-1
     xike(:,:,thismval)=xike(:,:,mod(thismval,2))
     if (thismval.ne.0) then
        do i=1,xigridpoints 
         xike(i,i,thismval) = xike(i,i,thismval) - thismval**2/(xipoints(i)**2-1)   
        enddo
     endif
  enddo

  etake=etake*(-0.5d0);  xike=xike*(-0.5d0)

  call pro_ham()

  if (bornopflag.ne.1) then
     call pro_lpm();     call get_both()

     yderivs=0.d0
     do iii=0,1
           if (iii.eq.0) then
              eta2=1d0;   eta1=1d0;   xi2=-1d0;   xi1=-1d0
           else
              eta2=-1d0;   eta1=-1d0;   xi2=-1d0;   xi1=-1d0  !! likely; starred in notebook
!!              eta2=1d0;   eta1=-1d0;   xi2=1d0;   xi1=-1d0
           endif
           do m=1,xigridpoints
              Yderivs(m,:,m,:,iii) = eta2 *  etadeg2(:,:,iii) + eta1 * etadeg1(:,:,iii)*mass_asym*xipoints(m) 
!! YEAH!
!              Yderivs(m,:,m,:,iii) =  etadeg2(:,:,iii) + etadeg1(:,:,iii)*mass_asym*xipoints(m) 
           enddo
             
           do i=1,numeta
              Yderivs(:,i,:,i,iii) = Yderivs(:,i,:,i,iii) + xi2* xideg2(:,:,iii) + xi1* xideg1(:,:,iii)*mass_asym*etapoints(i) 
!! YEAH!
!              Yderivs(:,i,:,i,iii) = Yderivs(:,i,:,i,iii) - xideg2(:,:,iii) - xideg1(:,:,iii)*mass_asym*etapoints(i) 
           enddo
     end do !! iii=0,1

     do m=1,xigridpoints
        do i=1,numeta
           Yderivs(m,i,:,:,:) = Yderivs(m,i,:,:,:) / sqrt( (xipoints(m)**2 - etapoints(i)**2)  )
           Yderivs(:,:,m,i,:) = Yderivs(:,:,m,i,:) / sqrt( (xipoints(m)**2 - etapoints(i)**2)  )
        enddo
     enddo
     do thismval=0,mbig
        do i=1,xigridpoints-1
           do j=1,xigridpoints-1
              myYderivs(:,:,:,:,thismval+1) = Yderivs(1:xigridpoints-1,:,1:xigridpoints-1,:,mod(thismval,2))  * (-0.5d0) / Rmass
           enddo
        enddo
     enddo

     do thismval=0,mbig
        proham(:,:,:,:,thismval+1)=proham(:,:,:,:,thismval+1) + pro_both(:,:,:,:,mod(thismval,2))*(-0.5d0/Rmass)
        do i=1,numeta
           do j=1,numerad
              proham(j,i,j,i,thismval+1) = proham(j,i,j,i,thismval+1) + 0.5d0/Rmass * thismval**2 * &
                   (1+ (etapoints(i)*xipoints(j)+mass_asym)**2/(xipoints(j)**2-1)/(1-etapoints(i)**2))
           enddo
        enddo
     enddo
  endif !! bornopflag

  sparseops_xi_banded(:,:,:,:)=0.d0;  sparseops_eta(:,:,:,:)=0.d0;  sparseops_diag(:,:,:)=0.d0
  sparseddz_xi_banded(:,:,:,:)=0.d0;  sparseddz_eta(:,:,:,:)=0.d0;  sparseddz_diag(:,:,:)=0.d0
  sparseddrho_xi_banded(:,:,:,:)=0.d0;  sparseddrho_eta(:,:,:,:)=0.d0;  sparseddrho_diag(:,:,:)=0.d0
  
  if (bornopflag.eq.0) then
     sparseyops_xi_banded(:,:,:,:)=0.d0;  sparseyops_eta(:,:,:,:)=0.d0;  sparseyops_diag(:,:,:)=0.d0
  endif
  do thismval=0,mbig
     do j=1,lbig+1
        do i=1,numerad
           do k=max(1,i-bandwidth),min(numerad,i+bandwidth)
              sparseddz_xi_banded(k-i+bandwidth+1,i,j,thismval+1) = proddz(k,j,i,j,thismval)
              sparseddrho_xi_banded(k-i+bandwidth+1,i,j,thismval+1) = proddrho(k,j,i,j,thismval)
           enddo
           do k=1,lbig+1
              sparseddz_eta(k,j,i,thismval+1) = proddz(i,k,i,j,thismval)
              sparseddrho_eta(k,j,i,thismval+1) = proddrho(i,k,i,j,thismval)
           enddo
           sparseddz_diag(i,j,thismval+1) = proddz(i,j,i,j,thismval)
           sparseddrho_diag(i,j,thismval+1) = proddrho(i,j,i,j,thismval)
        enddo
     enddo
     do j=1,lbig+1
        do i=1,numerad
           do k=max(1,i-bandwidth),min(numerad,i+bandwidth)
              sparseops_xi_banded(k-i+bandwidth+1,i,j,thismval+1) = proham(k,j,i,j,thismval+1)
           enddo
           do k=1,lbig+1
              sparseops_eta(k,j,i,thismval+1) = proham(i,k,i,j,thismval+1)
           enddo
           sparseops_diag(i,j,thismval+1) = proham(i,j,i,j,thismval+1)
        enddo
     enddo
     if (bornopflag.eq.0) then
        do j=1,lbig+1
           do i=1,numerad
              do k=max(1,i-bandwidth),min(numerad,i+bandwidth)
                 sparseyops_xi_banded(k-i+bandwidth+1,i,j,thismval+1) = myYderivs(k,j,i,j,thismval+1)
              enddo
              do k=1,lbig+1
                 sparseyops_eta(k,j,i,thismval+1) = myYderivs(i,k,i,j,thismval+1)
              enddo
              sparseyops_diag(i,j,thismval+1) = myYderivs(i,j,i,j,thismval+1)
           enddo
        enddo
     endif
  enddo   !! thismval

  call psctempdealloc()

contains

subroutine PSCtempalloc()
  implicit none

  allocate( xi_derivs(xigridpoints,xigridpoints,0:1), xi_rhoderivs(xigridpoints,xigridpoints), &
       xi_elderivs(xinumpoints*xinumelements,xigridpoints),       eta_derivs(numeta,numeta,0:1), &
       eta_rhoderivs(numeta,numeta))
  xi_derivs=0; xi_rhoderivs=0; xi_elderivs=0; eta_derivs=0; eta_rhoderivs=0

  if (bornopflag.eq.0) then
     allocate( &
       Yderivs(xigridpoints,numeta,xigridpoints,numeta,0:1) , &
       Ydiag(xigridpoints,numeta,xigridpoints,numeta,0:1), &
       Ytemp2(numeta,xigridpoints,numeta,xigridpoints),&
       pro_lop(xigridpoints,numeta,xigridpoints,numeta,0:1),&
       pro_ldiag(xigridpoints,numeta),pro_both(numerad,numeta,numerad,numeta,0:1))
     yderivs=0; ydiag=0; ytemp2=0; pro_lop=0; pro_ldiag=0
  endif
  allocate( etadeg1(numeta,numeta, 0:1),  etadeg2(numeta,numeta, 0:1), eta_second(numeta,numeta, 0:1), &
       eta_third(numeta,numeta, 0:1),  eta_fourth(numeta,numeta, 0:1),&
       xideg1(xigridpoints, xigridpoints,0:1),  xideg2(xigridpoints, xigridpoints,0:1),&
       xi_fourth(xigridpoints, xigridpoints,0:1),  xi_third(xigridpoints, xigridpoints,0:1),  &
       xi_second(xigridpoints, xigridpoints,0:1))
  etadeg1=0; etadeg2=0; eta_second=0; eta_third=0; eta_fourth=0; 
  xideg1=0; xideg2=0; xi_fourth=0; xi_third=0; xi_second=0

end subroutine PSCtempalloc

subroutine PSCtempdealloc()
  implicit none
  if (bornopflag.eq.0) then
     deallocate(Yderivs, Ydiag, ytemp2,pro_lop, pro_ldiag, pro_both)
  endif
  deallocate(  xi_derivs, eta_derivs,eta_rhoderivs,xi_rhoderivs, xi_elderivs, xi_fourth, xi_third, xi_second, &
       etadeg1, etadeg2, eta_fourth, eta_third, eta_second, xideg1, xideg2)
end subroutine PSCtempdealloc

! returns matrices which should be added
!
!   1/R * propot  +  1/R^2  * proham

subroutine pro_ham()
  implicit none
  DATAECS :: csum, csum2
  integer :: i,j,k,m,n, mvalue
  real*8 :: fac !, nucfac

  proham(:,:,:,:,:)=0.0d0
  do j=1,numerad
     do k=1,numeta
        
!! this multiplies 1/r to be the potential

        propot(j,k) = (nuccharge1+nuccharge2)/2.d0 *  &
             (-4.d0)* xipoints(j) / (xipoints(j)**2 - etapoints(k)**2) + &
             (nuccharge1-nuccharge2)/2.d0 * &
             (-4.d0)* etapoints(k) / (xipoints(j)**2 - etapoints(k)**2) 

! similarly  - but halfpot depends on numelec; multiply this by
!   (nuccharge1+nuccharge2-numelec+1)

        halfpot(j,k) = (-2.d0)* xipoints(j) / (xipoints(j)**2 - etapoints(k)**2) 
     enddo
  enddo

!!!    1/R^2 terms. 

!! for multiple electrons, we do not include mass polarization nor coriolis terms nor 
!! orbit-orbit two-electron interactions (l^+_1 l^-_2).
!! Hamiltonian is exact for one-electron diatomic molecules except for coriolis 
!! coupling (lambda doubling, coupling between electronic m-values due to nuclear KE) 
!! terms. Exact nonrelativistic one-electron diatomic hamiltonian for J=0 only.

  fac=1d0
  if (reducedflag.eq.1) then
     fac = ( totalmass + numelec ) / ( totalmass + numelec - 1 )
  endif

  do mvalue=0,mbig
     do j=1,numeta
        do m=1,numerad
           do n=1,numerad
              proham(n,j,m,j,mvalue+1) = proham(n,j,m,j,mvalue+1) + xike(n,m,mvalue) * 4 * fac
           enddo
        enddo
     enddo
     do j=1,numerad
        do m=1,numeta
           do n=1,numeta
              proham(j,n,j,m,mvalue+1) = proham(j,n,j,m,mvalue+1) + etake(n,m,mvalue) * 4 * fac
           enddo
        enddo
     enddo

! coriolis couplings between m-values not included, only this diagonal part

     if (bornopflag /= 1) then
        do j=1,numerad
           do k=1,numeta
              proham(j,k,j,k,mvalue+1) =  proham(j,k,j,k,mvalue+1) &
                   + mvalue**2 * ( xipoints(j)**2 - etapoints(k)**2) *  (-2.d0) *  0.5d0/Rmass
              proham(j,k,j,k,mvalue+1) =  proham(j,k,j,k,mvalue+1) &
                   + Jvalue*(Jvalue+1.d0) * ( xipoints(j)**2 - etapoints(k)**2) *  0.5d0/Rmass
           enddo
        enddo
     endif
  enddo

     proddz(:,:,:,:,:) = 0.d0;      proddrho(:,:,:,:,:) = 0.d0
     do mvalue=0,mbig   
        do m=1,numerad
           do i=1,numeta
              n=m
              do j=1,numeta
                 csum = eta_derivs(i,j,mod(mvalue,2)) * ( xipoints(m) + mass_asym*etapoints(i) ) 
                 proddz(m,i,n,j,mvalue) = proddz(m,i,n,j,mvalue) + 0.5d0 * csum
                 proddz(n,j,m,i,mvalue) = proddz(n,j,m,i,mvalue) - 0.5d0 * csum

                 csum2 = eta_rhoderivs(i,j) * sqrt((xipoints(m)**2-1)/(1-etapoints(i)**2))
                 proddrho(m,i,n,j,mvalue) = proddrho(m,i,n,j,mvalue) -  csum2  
              enddo
           enddo
        enddo
        do m=1,numerad
           do i=1,numeta
              do n=1,numerad
                 j=i
                 
                 csum = xi_derivs(m,n,mod(mvalue,2)) * ( etapoints(i) + mass_asym*xipoints(m) ) 
                 proddz(m,i,n,j,mvalue) = proddz(m,i,n,j,mvalue) + 0.5d0 * csum
                 proddz(n,j,m,i,mvalue) = proddz(n,j,m,i,mvalue) - 0.5d0 * csum

                 csum2 = xi_rhoderivs(m,n) * sqrt((1-etapoints(i)**2)/(xipoints(m)**2-1))
                 proddrho(m,i,n,j,mvalue) = proddrho(m,i,n,j,mvalue) +  csum2 
              enddo
           enddo
        enddo
     enddo  !! mvalue
     
     do m=1,numerad
        do i=1,numeta
           ddrhopot(m,i) =  1.d0/sqrt((xipoints(m)**2-1.d0)*(1.d0-etapoints(i)**2))  
        enddo
     enddo
     
     proddrho(:,:,:,:,:) = proddrho(:,:,:,:,:)*2.d0 
     ddrhopot(:,:) = ddrhopot(:,:)*2.d0 
     
!! 060811:
     proddz=proddz*2d0
!  endif


  do m=1,numerad
     do i=1,numeta
        proham(m,i,:,:,:) = proham(m,i,:,:,:) / sqrt( (xipoints(m)**2 - etapoints(i)**2)  )
        proham(:,:,m,i,:) = proham(:,:,m,i,:) / sqrt( (xipoints(m)**2 - etapoints(i)**2)  )

        proddz(m,i,:,:,:) = proddz(m,i,:,:,:) / sqrt( (xipoints(m)**2 - etapoints(i)**2)  )
        proddrho(m,i,:,:,:) = proddrho(m,i,:,:,:) / sqrt( (xipoints(m)**2 - etapoints(i)**2)  )
        
        proddz(:,:,m,i,:) = proddz(:,:,m,i,:) / sqrt( (xipoints(m)**2 - etapoints(i)**2)  )
        proddrho(:,:,m,i,:) = proddrho(:,:,m,i,:) / sqrt( (xipoints(m)**2 - etapoints(i)**2)  )
     enddo
  enddo

  do j=1,numeta
     do i=1,numerad

!! so if particle 1 (Hmass) is much heavier, this is zero at xi=1, eta=1
        zdipole(i,j) =             0.5d0 * (etapoints(j) * xipoints(i) + mass_asym)  
        xydipole(i,j) =             0.5d0 * sqrt( (1.d0-etapoints(j)**2) * (xipoints(i)**2 - 1.d0))
     enddo
  enddo
  
end subroutine pro_ham


subroutine eta_init(points,weights,ketot, numpoints, etadeg1, etadeg2, eta_derivs, &
     eta_rhoderivs, eta_fourth, eta_third, eta_second, evenodd)   !! 0=even 1=odd
  implicit none

  integer, parameter :: numextra=6
  integer,intent(in) :: numpoints,evenodd
  integer :: one=1, two=2, izero=0, i,j,k,jj , extraorder
  real*8,intent(out) :: ketot(numpoints,numpoints)  , points(numpoints),&
       etadeg1(numpoints,numpoints),        eta_second(numpoints,numpoints), &
       etadeg2(numpoints,numpoints), eta_fourth(numpoints,numpoints), &
       eta_third(numpoints,numpoints),       eta_derivs(numpoints,numpoints),&
       eta_rhoderivs(numpoints,numpoints), weights(numpoints)

  real*8 :: extrapoints(numpoints+numextra), extraweights(numextra+numpoints), &
       firstdertot(numpoints+numextra,numpoints),&
        scratch(numextra+numpoints), &
       etavals(numpoints+numextra, numpoints), &
       endpoints(2)=[-1.0d0,1.0d0] , zero = 0.0, sum, sum2, sum3, sum4

  extrapoints=0; extraweights=0; firstdertot=0; weights=0; scratch=0; eta_rhoderivs=0; etavals=0;

  if (evenodd.ne.0.and.evenodd.ne.1) then
     print *, "evenodd error", evenodd
  endif

  one=1;two=2;  extraorder=numextra+numpoints

  call gaussq(one,numpoints,zero,zero,izero,endpoints,scratch,points,weights)

!print *, "TEMPQUAD ETA!!"  !! (gl)
!  call gaussq(one,numpoints,zero,zero,two,endpoints,scratch,points,weights)

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

if (1==0) then
         sum=0.d0
         if (evenodd.eq.0) then 
            do k=1,extraorder
               sum=sum + extraweights(k) * etavals(k,i) * firstdertot(k,j)
            enddo
         else
            do k=1,extraorder
               sum= sum + etavals(k,i) &
                    * (firstdertot(k,j)* (1-extrapoints(k)**2) - etavals(k,j) * extrapoints(k) )  &
                    * extraweights(k) 
            enddo
            sum=sum/sqrt((1-points(i)**2)*(1-points(j)**2))
         endif
         eta_derivs(i,j)=sum
else
         sum=0.d0
         if (evenodd.eq.0) then 
            do k=1,extraorder
               sum= sum + etavals(k,i) &
                    * (firstdertot(k,j)* (1-extrapoints(k)**2) - etavals(k,j) * extrapoints(k) )  &
                    * extraweights(k) 
            enddo
         else
            do k=1,extraorder
               sum= sum + etavals(k,i) &
                    * (firstdertot(k,j)* (1-extrapoints(k)**2)**2 - &
                    2*etavals(k,j)*(1-extrapoints(k)**2) * extrapoints(k) ) * extraweights(k) 
            enddo
            sum=sum/sqrt((1-points(i)**2)*(1-points(j)**2))
         endif
         eta_derivs(i,j)=sum
endif

         sum=0
         do k=1,extraorder
            sum=sum + extraweights(k) * etavals(k,i) * ( extrapoints(k)*(1-extrapoints(k)**2) * &
                 firstdertot(k,j) + (0.5d0 - 3d0/2d0 * extrapoints(k)**2) * etavals(k,j) )
         enddo
         eta_rhoderivs(i,j)=sum

      enddo
   enddo

!! negative definite operators for eta.

   do j=1,numpoints
      do i=1,numpoints
         if (evenodd.eq.0) then !! even, d/dx (1-x^2)**2 d/dx
            sum=0.0d0;            sum2=0.0d0;            sum3=0.0d0;            sum4=0.0d0
            do k=1,extraorder
               sum = sum - firstdertot(k,j) * firstdertot(k,i) * extraweights(k) * (1-extrapoints(k)**2)
               sum2 = sum2 - firstdertot(k,j) * firstdertot(k,i) * extraweights(k) * (1-extrapoints(k)**2)**2
               sum3 = sum3 - firstdertot(k,j) * firstdertot(k,i) * extraweights(k) * (extrapoints(k)**3-extrapoints(k))
               sum4 = sum4 - firstdertot(k,j) * firstdertot(k,i) * extraweights(k) * (1-extrapoints(k)**2)
            enddo
         else
            sum=0.0d0;            sum2=0.0d0;            sum3=0.0d0;            sum4=0.0d0
            do k=1,extraorder
               sum= sum - (firstdertot(k,j)* (1-extrapoints(k)**2) - etavals(k,j) * extrapoints(k) )  &
                    * (firstdertot(k,i)* (1-extrapoints(k)**2) - etavals(k,i) * extrapoints(k) )  &
                    * extraweights(k) 

!! (1-x^2)**2 
               sum2= sum2 - (firstdertot(k,j)* (1-extrapoints(k)**2) - etavals(k,j) * extrapoints(k) ) &
                    * (1-extrapoints(k)**2) &
                    * (firstdertot(k,i)* (1-extrapoints(k)**2) - etavals(k,i) * extrapoints(k) )  &
                    * extraweights(k) 

!! (x-x^3)
               sum3= sum3 + (firstdertot(k,j)* (1-extrapoints(k)**2) - etavals(k,j) * extrapoints(k) ) &
                    * extrapoints(k) &
                    * (firstdertot(k,i)* (1-extrapoints(k)**2) - etavals(k,i) * extrapoints(k) )  &
                    * extraweights(k)  

!! (1-x^2)
               sum4= sum4 - (firstdertot(k,j)* (1-extrapoints(k)**2) - etavals(k,j) * extrapoints(k) )  &
                    * (firstdertot(k,i)* (1-extrapoints(k)**2) - etavals(k,i) * extrapoints(k) )  &
                    * extraweights(k)  
               
               
!               print *, "ETA SUM ", sum
            enddo
            sum=sum/sqrt((1-points(i)**2)*(1-points(j)**2))
            sum2=sum2/sqrt((1-points(i)**2)*(1-points(j)**2))
            sum3=sum3/sqrt((1-points(i)**2)*(1-points(j)**2))
            sum4=sum4/sqrt((1-points(i)**2)*(1-points(j)**2))
!            stop
         endif
         ketot(i,j)=sum   !! symmetric kinetic energy
         eta_fourth(i,j)=sum2 
         eta_third(i,j)=sum3
         eta_second(i,j)=sum4
      enddo
   enddo


   do i=1,numpoints
      do j=1,numpoints
         sum=0.d0;         sum2=0.d0
         if (evenodd.eq.0) then
            do k=1,extraorder
               sum=sum + extraweights(k) * etavals(k,i) * &
                    ( (extrapoints(k)**2-1d0) * firstdertot(k,j) + extrapoints(k) * etavals(k,j) )
               sum2=sum2 + extraweights(k) * etavals(k,i) * &
                    ( extrapoints(k)*(extrapoints(k)**2-1d0) * firstdertot(k,j) + &
                    ((-0.5d0) + 3d0/2d0 * extrapoints(k)**2) * etavals(k,j) )
            enddo
         else
            do k=1,extraorder
               sum=sum + extraweights(k) * etavals(k,i) * &
                    ( (extrapoints(k)**2-1d0)**2 * firstdertot(k,j) &
                    - 2*extrapoints(k)*(1-extrapoints(k)**2) * etavals(k,j) )
               sum2=sum2 + extraweights(k) * etavals(k,i) * &
                    ( extrapoints(k)*(extrapoints(k)**2-1d0)**2 * firstdertot(k,j) + &
                    (1d0/2d0 - 3 * extrapoints(k)**2 + 5./2.*extrapoints(k)**4) * etavals(k,j) )
            enddo
            sum=sum/sqrt((1-points(i)**2)*(1-points(j)**2))
            sum2=sum2/sqrt((1-points(i)**2)*(1-points(j)**2))
         endif
         etadeg1(i,j)=sum;         etadeg2(i,j)=sum2
      enddo
   enddo

end subroutine eta_init


subroutine xi_init(points,weights,points2d,ketot, numpoints,numelements,elementsizes, gridpoints, celement,&
     ecstheta, lobattoopt, xideg1,xideg2, xi_derivs, xi_rhoderivs, xi_fourth, xi_third, xi_second, evenodd)
  implicit none
  integer,intent(in) :: numpoints,numelements,gridpoints,celement,lobattoopt,evenodd
  real*8,intent(in) ::  elementsizes(numelements), ecstheta
  DATAECS, intent(out) ::        ketot(gridpoints,gridpoints),points(gridpoints), weights(gridpoints),&
       xi_rhoderivs(gridpoints,gridpoints),&
       xi_derivs(gridpoints,gridpoints), xideg1(gridpoints,gridpoints),xideg2(gridpoints,gridpoints), &
       xi_fourth(gridpoints,gridpoints), xi_third(gridpoints,gridpoints), xi_second(gridpoints,gridpoints)
  real*8, intent(out) :: points2d(numpoints,2)
  integer ::   one=1, izero=0, two=2, i,j,k,l,jj, qq, extraorder
  integer, parameter :: numextra=11
  real*8 :: extrapoints0(numpoints+numextra), extraweights0(numpoints+numextra), &
       firstder(numpoints+numextra,numpoints,2),  weights2d(numpoints,2), &
       scratch(2*numpoints+numextra), xivals0(numpoints+numextra,numpoints,2),  endpoints(2), &
       zero = 0.0,rsum,rsum2
  DATAECS :: extrapoints(numpoints+numextra,numelements), extraweights(numextra+numpoints,numelements), &
       firstdertot(numpoints+numextra,numelements,gridpoints), &
       xivals(numpoints+numextra,numelements,gridpoints), xi_test(gridpoints,gridpoints), &
       xi_vects(gridpoints,gridpoints), xi_vals(gridpoints) 
  DATAECS :: cweight, sum, sum2, sum3, sum4, sum5
  i=celement; sum=ecstheta !! avoid warn unused

  extrapoints0=0; extraweights0=0; firstder=0; points2d=0; weights2d=0; scratch=0;
  xivals0=0; extrapoints=0; extraweights=0; firstdertot=0; xivals=0; xi_test=0; 
  xi_vals=0; xi_vects=0

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
  call gaussq(one,numpoints,zero,zero,two,endpoints,scratch,points2d(:,2),weights2d(:,2))

  select case(lobattoopt)
  case (0)
     endpoints= [ 1.0d0 , 0.d0 ]
     call gaussq(one,numpoints,zero,zero,one,endpoints,scratch,points2d(:,1),weights2d(:,1))
  case (1)
     endpoints= [ -1.0d0 , 1.d0 ]
     call gaussq(one,numpoints,zero,zero,two,endpoints,scratch,points2d(:,1),weights2d(:,1))
  case (2)
     print *, "GAUSS LOBATTO TRUE DVR FIRST ELEMENT   BETA !!!!"
!     endpoints= [ -1.0d0 , 1.d0 ]
     print *, "not supported";     stop
!     call gaussq(one,numpoints*2,zero,zero,two,endpoints,scratch,points2d(:,1),weights2d(:,1))
!     points2d=points2d*2d0-1d0   !! points2d -3 to 1
!     weights2d=weights2d*2
  case default
     print *, "AAAAAAUGH what is this lobatto opt!!!", lobattoopt;     stop
  end select

!  print *, "TEMPQUAD XI!!!"
!  endpoints= [ -1.0d0 , 1.d0 ]
!  call gaussq(one,numpoints,zero,zero,two,endpoints,scratch,points2d(:,1),weights2d(:,1))
!  call gaussq(one,numpoints,zero,zero,izero,endpoints,scratch,points2d(:,1),weights2d(:,1))

  xivals0=0.d0
  do k=1,extraorder  ! point k
     do j=1,numpoints  ! which basis function
        rsum=1;        rsum2=1
        do i=1,numpoints
           if (j/=i) then
              rsum2=rsum2 * (extrapoints0(k)-points2d(i,2))/(points2d(j,2)-points2d(i,2)) 
              rsum=rsum * (extrapoints0(k)-points2d(i,1))/(points2d(j,1)-points2d(i,1))
           endif
        enddo
        xivals0(k,j,1) = rsum;        xivals0(k,j,2) = rsum2
     enddo
  enddo

  if (numextra.ne.0) then
     do k=1,extraorder  ! point k
        do j=1,numpoints
           rsum=0.d0;           rsum2=0.d0
           do i=1,numpoints
              if (j/=i) then
                 rsum=rsum+ xivals0(k,j,1)/(extrapoints0(k)-points2d(i,1))
                 rsum2=rsum2+ xivals0(k,j,2)/(extrapoints0(k)-points2d(i,2))
              endif
           enddo
           firstder(k,j,1)=rsum;           firstder(k,j,2)=rsum2
        enddo
     enddo
  else
     print *, "Numextra not supported.", numextra;     stop
  endif

  weights=0.0d0;  j=0;  jj=0;  sum=1.d0;  firstdertot=0.d0;  xivals=0d0

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
!!BUG? 
!!         firstder(:,:,min(l,2))*cweight*2.d0/elementsizes(l)
!! FIX?
          firstder(:,:,min(l,2))/cweight*2.d0/elementsizes(l)

     xivals(:,l,qq+1:qq+numpoints)=xivals(:,l,qq+1:qq+numpoints)+xivals0(:,:,min(l,2))
     extraweights(:,l)=extraweights0(:)* elementsizes(l)/2.0d0 * cweight
     extrapoints(:,l)= cweight*elementsizes(l)/2.0d0 * (extrapoints0(:)+1d0) + sum

     do k=1,numpoints
        jj=jj+1
        if ((l==1) .or. (k/=1)) then    
           j = j+1
        endif
!!        if (lobattoopt.ne.2.or.l.ne.1) then  ! DELETED LOBATTOOPT 2
        points(j) = cweight*elementsizes(l)/2.0d0 * (points2d(k,min(2,l))+1.0d0) + sum
        weights(j) = weights(j) + weights2d(k,min(2,l)) * elementsizes(l)/2.0d0 * cweight
     enddo
     sum=sum+elementsizes(l)*cweight
  enddo

  !! now deriv and value are deriv and value of normalized FEM DVR function.

  do jj=1,gridpoints
     firstdertot(:,:,jj) = firstdertot(:,:,jj) / sqrt(weights(jj)) 
     xivals(:,:,jj) = xivals(:,:,jj) / sqrt(weights(jj)) 
  enddo
  do i=1,gridpoints
     do j=1,gridpoints
        
        if (1==0) then
           
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
                    sum= sum + xivals(k,l,i) &
                         * (firstdertot(k,l,j)* (extrapoints(k,l)**2-1) + xivals(k,l,j) * extrapoints(k,l) )  &
                         * extraweights(k,l) 
                 enddo
              enddo
              sum=sum/sqrt((1-points(i)**2)*(1-points(j)**2))
           endif
           xi_derivs(i,j)=sum
           
        else
           
           sum=0.d0
           if (evenodd.eq.0) then
              do l=1,numelements
                 do k=1,extraorder
                    sum= sum + xivals(k,l,i) &
                         * (firstdertot(k,l,j)* (extrapoints(k,l)**2-1) + xivals(k,l,j) * extrapoints(k,l) )  &
                         * extraweights(k,l) 
                 enddo
              enddo
           else
              do l=1,numelements
                 do k=1,extraorder
                    sum= sum + xivals(k,l,i) &
                         * (firstdertot(k,l,j)* (extrapoints(k,l)**2-1)**2 + &
                         2*xivals(k,l,j)*(extrapoints(k,l)**2-1) * extrapoints(k,l) )  &
                         * extraweights(k,l) 
                 enddo
              enddo
              sum=sum/sqrt((points(i)**2-1)*(points(j)**2-1))
           endif
           xi_derivs(i,j)=sum
        endif
        
        sum=0
        do l=1,numelements
           do k=1,extraorder
              sum=sum + extraweights(k,l) * xivals(k,l,i) * &
                   ( extrapoints(k,l)*(extrapoints(k,l)**2-1d0) * firstdertot(k,l,j) + &
                   ((-0.5d0) + 3d0/2d0 * extrapoints(k,l)**2) * xivals(k,l,j) )
           enddo
        enddo
        xi_rhoderivs(i,j)=sum
     enddo
  enddo
  
!! negative definite operators for xi.

  do j=1,gridpoints
     do i=1,gridpoints
        if (evenodd.eq.0) then !! even, d/dx (1-x^2)**2 d/dx
           sum=0.0d0;           sum2=0.d0;     sum3=0.d0;    sum4=0.d0;      sum5=0.d0
           do l=1,numelements
              do k=1,extraorder

                 sum  = sum  - firstdertot(k,l,j) * firstdertot(k,l,i) * &
                      extraweights(k,l) * (extrapoints(k,l)**2-1)
!! (x^2-1)^2                 
                 sum2 = sum2 - firstdertot(k,l,j) * firstdertot(k,l,i) * &
                      extraweights(k,l) * (extrapoints(k,l)**2-1)**2
!! (x^3-x)
                 sum3 = sum3 - firstdertot(k,l,j) * firstdertot(k,l,i) * &
                      extraweights(k,l) * (extrapoints(k,l)**3-extrapoints(k,l)) 
!! (x^2-1)
                 sum4 = sum4 - firstdertot(k,l,j) * firstdertot(k,l,i) * &
                      extraweights(k,l) * (extrapoints(k,l)**2-1)

                 sum5  = sum5  - firstdertot(k,l,j) * firstdertot(k,l,i) * extraweights(k,l) 
              enddo
           enddo
           xi_test(i,j)=sum5
        else

           sum=0.0d0;           sum2=0.d0;           sum3=0.d0;           sum4=0.d0
           do l=1,numelements
              do k=1,extraorder   
                 sum= sum - (firstdertot(k,l,j)* (extrapoints(k,l)**2-1) + xivals(k,l,j) * extrapoints(k,l) )  &
                      * (firstdertot(k,l,i)* (extrapoints(k,l)**2-1) + xivals(k,l,i) * extrapoints(k,l) )  &
                      * extraweights(k,l) 

!! (x^2-1)**2 
                 sum2= sum2 - (firstdertot(k,l,j)* (extrapoints(k,l)**2-1) + xivals(k,l,j) * extrapoints(k,l) ) &
                      * (extrapoints(k,l)**2 -1) &
                      * (firstdertot(k,l,i)* (extrapoints(k,l)**2-1) + xivals(k,l,i) * extrapoints(k,l) )  &
                      * extraweights(k,l) 
!!  (x^3-x)
                 sum3= sum3 - (firstdertot(k,l,j)* (extrapoints(k,l)**2-1) + xivals(k,l,j) * extrapoints(k,l) ) &
                      * extrapoints(k,l) &
                      * (firstdertot(k,l,i)* (extrapoints(k,l)**2-1) + xivals(k,l,i) * extrapoints(k,l) )  &
                      * extraweights(k,l) 
!!  (x^2-1)
                 sum4= sum4 - (firstdertot(k,l,j)* (extrapoints(k,l)**2-1) + xivals(k,l,j) * extrapoints(k,l) )  &
                      * (firstdertot(k,l,i)* (extrapoints(k,l)**2-1) + xivals(k,l,i) * extrapoints(k,l) )  &
                      * extraweights(k,l) 


!                 print *, "XI SUM ", sum
              enddo
           enddo
           sum=sum/sqrt((points(i)**2-1)*(points(j)**2-1))
           sum2=sum2/sqrt((points(i)**2-1)*(points(j)**2-1))
           sum3=sum3/sqrt((points(i)**2-1)*(points(j)**2-1))
           sum4=sum4/sqrt((points(i)**2-1)*(points(j)**2-1))
!           stop
        endif
        ketot(i,j)=sum   !! symmetric kinetic energy
        xi_fourth(i,j)=sum2  ;        xi_third(i,j)=sum3   ;        xi_second(i,j)=sum4   
     enddo
  enddo

  do i=1,gridpoints
     do j=1,gridpoints
        sum=0.d0;        sum2=0.d0
        if (evenodd.eq.0) then
           do l=1,numelements
              do k=1,extraorder
                 sum=sum + extraweights(k,l) * xivals(k,l,i) * ( (extrapoints(k,l)**2-1d0) * &
                      firstdertot(k,l,j) + extrapoints(k,l) * xivals(k,l,j) )
                 sum2=sum2 + extraweights(k,l) * xivals(k,l,i) * ( extrapoints(k,l)*(extrapoints(k,l)**2-1d0) * &
                      firstdertot(k,l,j) - (0.5d0 - 3d0/2d0 * extrapoints(k,l)**2) * xivals(k,l,j) )
              enddo
           enddo
        else
           do l=1,numelements
              do k=1,extraorder
                 sum=sum + extraweights(k,l) * xivals(k,l,i) * ( (extrapoints(k,l)**2-1d0)**2 * &
                      firstdertot(k,l,j) - 2*extrapoints(k,l)*(1-extrapoints(k,l)**2) * xivals(k,l,j) )
                 sum2=sum2 + extraweights(k,l) * xivals(k,l,i) * ( extrapoints(k,l)*(extrapoints(k,l)**2-1d0)**2 * &
                      firstdertot(k,l,j) + &
                      (1d0/2d0 - 3 * extrapoints(k,l)**2 + 5./2.*extrapoints(k,l)**4) * xivals(k,l,j) )
              enddo
           enddo
           sum=sum/sqrt((1-points(i)**2)*(1-points(j)**2))
           sum2=sum2/sqrt((1-points(i)**2)*(1-points(j)**2))
        endif
        xideg1(i,j)=sum;        xideg2(i,j)=sum2
     enddo
  enddo

end subroutine xi_init


subroutine prolate_init(points,weights, numpoints,numelements,elementsize,gridpoints,start, &
     prolate_derivs, ketot, celement, ecstheta)
  implicit none
  integer,intent(in) :: numpoints,numelements,gridpoints, celement
  DATAECS,intent(out) ::  prolate_derivs(gridpoints,gridpoints), points(gridpoints), weights(gridpoints), &
       ketot(gridpoints,gridpoints)
  real*8,intent(in) ::  elementsize,start, ecstheta
  call prolate_init_new(points,weights, numpoints,numelements,elementsize,gridpoints,start, &
       prolate_derivs, ketot, celement, ecstheta)
end subroutine prolate_init


subroutine prolate_init_new(points,weights, numpoints,numelements,elementsize,gridpoints,start, &
     prolate_derivs, ketot, celement, ecstheta)
  use myparams
  implicit none
  integer, parameter :: numextra=20
  integer,intent(in) :: numpoints,numelements,gridpoints,celement
  real*8,intent(in) :: elementsize,start,ecstheta
  integer ::one=1, two=2, izero=0, i,j,k,l,jj,kk, qq, extraorder
  DATAECS,intent(out) ::  prolate_derivs(gridpoints,gridpoints),&
       points(gridpoints), weights(gridpoints), ketot(gridpoints,gridpoints)
  DATAECS :: sum, cweight
  real*8 ::  endpoints(2)=[-1.0d0,1.0d0],  zero = 0.0, rsum, rsum2
  real*8, allocatable :: firstderiv(:,:), ke(:,:), points2d(:), weights2d(:), scratch(:), &  
       extrapoints0(:), extraweights0(:), rvals0(:,:), firstder(:,:)
  DATAECS, allocatable :: extrapoints(:,:), extraweights(:,:), rvals(:,:,:), firstdertot(:,:,:), &
       cextrapoints(:,:), cextraweights(:,:)
  i=celement; sum=ecstheta !! avoid warn unused

  extraorder=numpoints+numextra

  allocate(points2d(numpoints), weights2d(numpoints), scratch(extraorder),  &
       firstderiv(numpoints,numpoints), ke(numpoints,numpoints),&
       extrapoints(extraorder,numelements), extraweights(extraorder,numelements), rvals0(extraorder,numpoints), &
       rvals(extraorder,numelements,gridpoints), firstder(extraorder,numpoints), &
       firstdertot(extraorder,numelements,gridpoints), extrapoints0(extraorder), extraweights0(extraorder),&
       cextrapoints(extraorder,numelements), cextraweights(extraorder,numelements))

  points2d=0; weights2d=0; firstderiv=0; ke=0; extrapoints=0; extraweights=0; rvals0=0;
  rvals=0; firstder=0; firstdertot=0; extrapoints0=0; extraweights0=0;
  cextrapoints=0; cextraweights=0

  call gaussq(one,numpoints,zero,zero,two,endpoints,scratch,points2d,weights2d)
  call gaussq(one,extraorder,zero,zero,izero,endpoints,scratch,extrapoints0(:),extraweights0(:))

  do j=1,numpoints
     do i=1,numpoints
        if (i/=j) then
           firstderiv(i,j) = 1.0d0 / (points2d(i)-points2d(j))
           do k=1,numpoints  
              if ( (k/=i) .and. (k/=j) ) then
                 firstderiv(i,j) = firstderiv(i,j) * (points2d(j)-points2d(k))/(points2d(i)-points2d(k))
              endif
           enddo
        else
           firstderiv(i,i) = 0.0d0
           do k=1,numpoints
              if (k/=i) then
                 firstderiv(i,i) = firstderiv(i,i) + 1.d0 / (points2d(i) - points2d(k))
              endif
           enddo
        endif
     enddo
  enddo
  do j=1,numpoints
     do i=1,numpoints
        ke(i,j)=0.0d0
        do k=1,numpoints
           ke(i,j) = ke(i,j) - firstderiv(j,k) * firstderiv(i,k) *weights2d(k)
        enddo
     enddo
  enddo

  weights=0.0d0;  j=0;  jj=0;  sum=start

  do l=1,numelements
     cweight=1.d0
#ifdef ECSFLAG
     if (l.ge.celement) then
        cweight=exp((0.d0, 1.d0)*ecstheta)
     endif
#endif
     
     do k=1,numpoints
        jj=jj+1
        if ((l==1) .or. (k/=1)) then    
           j = j+1
        endif
        points(j) = sum + cweight*elementsize/2.0d0 * (points2d(k)+1.0d0)
        weights(j) = weights(j) + weights2d(k) * elementsize/2.0d0 * cweight
     enddo
     sum=sum+elementsize*cweight
  enddo

  ketot = 0.0d0
  do l=1,numelements
     cweight=1.d0
#ifdef ECSFLAG
     if (l.ge.celement) then
        cweight=exp((0.d0, 1.d0)*ecstheta)
     endif
#endif
     do j=1,numpoints
        jj = (l-1)*(numpoints-1)+j
        do k=1,numpoints
           kk = (l-1)*(numpoints-1)+k
           ketot(jj,kk) = ketot(jj,kk) + ke(j,k) / cweight
        enddo
     enddo
  enddo
  do jj=1,gridpoints
     do kk=1,gridpoints
        ketot(jj,kk) = ketot(jj,kk) / sqrt(weights(jj)*weights(kk))
     enddo
  enddo

  ketot = ketot * (2.0/elementsize)    !! DONE WITH D^2/DR^2   (BILL's WAY!)

!! NEW- prolate deriv operator 2/R d/dR + 1/R

  rvals0=0.d0
  do k=1,extraorder  ! point k
     do j=1,numpoints  ! which basis function
        rsum=1
        do i=1,numpoints
           if (j/=i) then
              rsum=rsum * (extrapoints0(k)-points2d(i))/(points2d(j)-points2d(i))
           endif
        enddo
        rvals0(k,j) = rsum
     enddo
  enddo
  do k=1,extraorder  ! point k
     do j=1,numpoints
        rsum=0.d0
        do i=1,numpoints
           if (j/=i) then
              rsum=rsum+ rvals0(k,j)/(extrapoints0(k)-points2d(i))
           endif
        enddo
        firstder(k,j)=rsum
     enddo
  enddo

  firstdertot=0.d0;  sum=start

  do l=1,numelements
     cweight=1.d0
#ifdef ECSFLAG
     if (l.ge.celement) then
        cweight=exp((0.d0, 1.d0)*ecstheta)
     endif
#endif
     qq = (numpoints-1)*(l-1)

     !! for now value of derivative of interp function at point
     firstdertot(:,l,qq+1:qq+numpoints)=firstdertot(:,l,qq+1:qq+numpoints) + &
          firstder(:,:)*2.d0/elementsize       /cweight

     rvals(:,l,qq+1:qq+numpoints)=rvals(:,l,qq+1:qq+numpoints)+rvals0(:,:)
     extraweights(:,l)=extraweights0(:)* elementsize/2.0d0 
     extrapoints(:,l)= elementsize/2.0d0 * (extrapoints0(:)+1d0) + sum

     cextraweights(:,l)=extraweights0(:)* elementsize/2.0d0 *cweight
     cextrapoints(:,l)= elementsize/2.0d0 * (extrapoints0(:)+1d0)*cweight + sum

     sum=sum+elementsize * cweight
  enddo

  !! now deriv and value are deriv and value of normalized FEM DVR function.

  do jj=1,gridpoints
     firstdertot(:,:,jj) = firstdertot(:,:,jj) / sqrt(weights(jj)) 
     rvals(:,:,jj) = rvals(:,:,jj) / sqrt(weights(jj)) 
  enddo

!! hmm, apparently firstdertot is (-1) what it should be...

  do j=1,gridpoints
     do i=1,gridpoints
        sum=0.0d0
        do l=1,numelements
           do k=1,extraorder
              
!!!              sum = sum -  rvals(k,l,j) * ( 2.d0/extrapoints(k,l) * firstdertot(k,l,i) - &
!!!!               1.d0/extrapoints(k,l)**2 * rvals(k,l,i) ) * extraweights(k,l) 

              sum = sum -  rvals(k,l,j) * ( 2.d0/cextrapoints(k,l) * firstdertot(k,l,i) - &
                   1.d0/cextrapoints(k,l)**2 * rvals(k,l,i) ) * cextraweights(k,l) 

!!              sum = sum -  rvals(k,l,j) * firstdertot(k,l,i)  * extraweights(k,l) 

           enddo
        enddo
        prolate_derivs(i,j)=sum
     enddo
  enddo

  rsum=0.d0
  do i=2,gridpoints-1
     do j=i,gridpoints-1
        rsum2=abs(prolate_derivs(i,j)+prolate_derivs(j,i))
        if (rsum2.gt.rsum) then
           rsum=rsum2
        endif
     enddo
  enddo
  if (rsum.gt.1d-5) then
     OFLWR "MAX prolate_derivs asym :  ", rsum; write(mpifileptr,*);CFL
  endif

  prolate_derivs= (0.5d0) * (prolate_derivs - transpose(prolate_derivs))
  deallocate(points2d,weights2d,scratch,  ke)

end subroutine prolate_init_new


subroutine r_init(points,weights,ketot,  numpoints,numelements,elementsize,gridpoints,start)
  implicit none
  integer,intent(in) :: numpoints,numelements, gridpoints
  integer :: one=1, two=2, i,j,k,l,jj,kk
  real*8,intent(in) :: start, elementsize
  real*8, intent(out) :: points(gridpoints), weights(gridpoints),ketot(gridpoints,gridpoints)
  real*8 :: points2d(numpoints), weights2d(numpoints), &
       scratch(numpoints),endpoints(2)=[-1.0d0,1.0d0],zero = 0.0
 ! ke(i,j,k) is the i,j-th matrix element of the ke 
 !              matrix for the kth finite element
  real*8 :: ke(numpoints,numpoints), firstderiv(numpoints,numpoints)

  points2d=0; weights2d=0; scratch=0; ke=0; firstderiv=0

  call gaussq(one,numpoints,zero,zero,two,endpoints,scratch,points2d,weights2d)

 ! firstderiv(i,j)  First derivative matrix.  Deriv of ith dvr at x_j

  do j=1,numpoints
     do i=1,numpoints
        if (i/=j) then
           firstderiv(i,j) = 1.0d0 / (points2d(i)-points2d(j))
           do k=1,numpoints  
              if ( (k/=i) .and. (k/=j) ) then
                 firstderiv(i,j) = firstderiv(i,j) * (points2d(j)-points2d(k))/(points2d(i)-points2d(k))
              endif
           enddo
        else

!  I FORGET HOW THIS GOES - THE TWO SECTIONS BELOW ARE EQUIV FOR G-L but not for GAUSS.

!           if (i.eq.1) then
!              firstderiv(i,i) = -0.5d0/weights2d(i)
!           elseif (i.eq.numpoints) then
!              firstderiv(i,i) = 0.5d0/weights2d(i)
!           else
!              firstderiv(i,i) = 0
!           endif

         firstderiv(i,i) = 0.0d0
         do k=1,numpoints
            if (k/=i) then
               firstderiv(i,i) = firstderiv(i,i) + 1.d0 / (points2d(i) - points2d(k))
            endif
           enddo
        endif
     enddo
  enddo

  do j=1,numpoints
     do i=1,numpoints
        ke(i,j)=0.0d0
        do k=1,numpoints
           ke(i,j) = ke(i,j) - firstderiv(j,k) * firstderiv(i,k) *weights2d(k)
        enddo
     enddo
  enddo

  

! Now I have kinetic energy matrices for each set of Gauss-Lobatto points including
! endpoints, with respect to basis functions NOT normalized on (-1,1) (0.01 - amplitude 1).  
! Put them together
! into FEM KE matrices.

  weights=0.0d0
  
  j=0
  do l=1,numelements
     do k=1,numpoints
        if ((l==1) .or. (k/=1)) then    
           j = j+1
        endif
        points(j) = (l-1)*elementsize + elementsize/2.0d0 * (points2d(k)+1.0d0) + start
        weights(j) = weights(j) + weights2d(k) * elementsize/2.0d0
     enddo
  enddo
  ketot = 0.0d0

  do l=1,numelements
     do j=1,numpoints
        jj = (l-1)*(numpoints-1)+j
        do k=1,numpoints
           kk = (l-1)*(numpoints-1)+k
           ketot(jj,kk) = ketot(jj,kk) + ke(j,k)
        enddo
     enddo
  enddo
  do jj=1,gridpoints
     do kk=1,gridpoints
        ketot(jj,kk) = ketot(jj,kk) / sqrt(weights(jj)*weights(kk))
     enddo
  enddo
  ketot = ketot * 2.0/elementsize

end subroutine r_init

end subroutine PSC

end module


!! NOT A FUNCTION OF rvalue YET - still using easy prolate coords only - need to work on it

subroutine nucdipvalue(rvaluenotused,dipoles)
  use myparams
  implicit none
  real*8 :: nucfac

  DATATYPE,intent(in) :: rvaluenotused(1)
  DATATYPE,intent(out) :: dipoles(3)

!! nuccharge1 goes with pro_Hmass, nuccharge2 with pro_Dmass.
!! nucleus 1 (H) is on top (positive z).
!! if nucleus1 (H) is heavier, mass_asym is negative; if it is
!! heavier, it is closer to the origin.  (1+mass_asym) multiplies
!! nuccharge1.
!! This reports (-1) times the dipole because that is how the 
!! rest of the program is done (the program is correct for antimatter,
!! the dipole operator is programmed as +z)

  nucfac=( (1.d0-mass_asym)*NucCharge2 - (1.d0+mass_asym)*NucCharge1 )/2.d0
  dipoles(:)=0d0
!!  dipoles(3)=nucfac*rvalue(1)
  dipoles(3)=nucfac              !! will be multiplied by rvalue (length gauge only)

end subroutine nucdipvalue

