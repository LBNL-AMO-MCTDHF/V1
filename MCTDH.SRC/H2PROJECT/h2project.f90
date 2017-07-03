
#include "Definitions.INC"


module myprojecttempmod
  implicit none

  real*8, allocatable :: eta_derivs(:,:,:), eta_rhoderivs(:,:)
  DATAECS, allocatable ::  Ytemp2(:,:,:,:), Yderivs(:,:,:,:,:), Ydiag(:,:,:,:,:), xi_derivs(:,:,:), &
       xi_elderivs(:,:), Yallderivs(:,:,:,:,:),xi_rhoderivs(:,:)
  DATAECS, allocatable :: xi_fourth(:,:,:), xi_third(:,:,:), xi_second(:,:,:), xideg1(:,:,:), xideg2(:,:,:)
  real*8, allocatable :: etadeg1(:,:,:), etadeg2(:,:,:), eta_second(:,:,:), eta_third(:,:,:), eta_fourth(:,:,:)
  DATAECS, allocatable  :: pro_lop(:,:,:,:,  :), pro_both(:,:,:,:,  :)  !! last index 0:1 even or odd
  DATAECS, allocatable  :: pro_ldiag(:,:)

end module myprojecttempmod


module myprojectmod
  implicit none
  real*8, allocatable :: usualke(:,:)
  DATAECS, allocatable :: usualke2(:,:)
  DATAECS, allocatable :: usualvects(:,:,:)
  DATAECS, allocatable :: usualvals(:,:)
  DATAECS, allocatable :: myYderivs(:,:,:,:, :) ! for operating with mvalue (streamline this later)
  DATAECS, allocatable  :: proham(:,:,:,:,   :),    propot(:,:), halfpot(:,:)

!! dimensioning these xigridpoints now!!  check nonsparse mult!! 

  DATAECS, allocatable  :: pro_lplus(:,:,:,:,  :)  !! (-mbig-1):(mbig+1)
  DATAECS, allocatable  :: pro_lminus(:,:,:,:,  :)  
  DATAECS, allocatable  :: pro_lsquared(:,:,:,:,   :) !! (-mbig-1):(mbig+1)

!! NOW ECS:
  DATAECS, allocatable :: rketot(:,:)
  DATAECS, allocatable ::    rpoints(:),   prolate_derivs(:,:), rweights(:)
  DATAECS, allocatable :: xipoints(:),  xike(:,:,:)
  real*8, allocatable :: xipoints2d(:,:)
  DATAECS, allocatable :: xiweights(:)
  DATAECS, allocatable :: proddz(:,:,:,:,:) ! for velocity.  Mult by 2/R
  DATAECS, allocatable :: proddrho(:,:,:,:,:) ! for velocity.  Mult by 2/R
  DATAECS, allocatable,target :: ddrhopot(:,:) ! for velocity.  Mult by 2/R
  DATAECS, allocatable :: xydipole(:,:),zdipole(:,:)
  real*8, allocatable ::  etapoints(:),  etake(:,:,:)
  real*8, allocatable :: etaweights(:)

  DATATYPE, allocatable, target :: &
       sparseops_xi_banded(:,:,:,:),sparseyops_xi_banded(:,:,:,:),&
       sparseops_eta(:,:,:,:),sparseyops_eta(:,:,:,:),&
       sparseops_diag(:,:,:),sparseyops_diag(:,:,:),&
       sparselplus_xi_banded( :,  :,  :,  :), sparselminus_xi_banded( :, :,  :,   :), &
       sparselplus_eta(  :,   :, :,   :    ), sparselminus_eta(:,:,:,:), &
       sparselplus_diag(   :,  :,  :), sparselminus_diag(:,:,:), &
       sparseddz_xi_banded(  :,   :,    :,     :), &
       sparseddz_eta(  :,   :,    :,     :), &
       sparseddz_diag(   :,    :,     :), &
       sparseddrho_xi_banded(  :,   :,    :,     :), &
       sparseddrho_eta(  :,   :,    :,     :), &
       sparseddrho_diag(   :,    :,     :)

  DATAECS, allocatable :: rmatrix(:,:,:,:)

!!! Eq.68, McCurdy Baertschy Rescigno
!  DATAECS :: rmatrix(numerad,numerad,mseriesmax+1,lseriesmax+1)   

  real*8, allocatable :: ylmvals(:,:,:)

end module myprojectmod

subroutine myprojectalloc()
  use myparams
  use myprojectmod  
  implicit none

  call mpibarrier()
  OFLWR "allocating for h2 project.... "; CFL
  call mpibarrier()

  allocate( usualke(rgridpoints,rgridpoints) , usualke2(rgridpoints,rgridpoints) )
  allocate( usualvects(rgridpoints,rgridpoints,0:mbig) ,usualvals(rgridpoints,0:mbig) )
  allocate( xiKE(xigridpoints, xigridpoints,0:max(mbig+2,mseriesmax+1)) )
  allocate(proham(numerad,numeta,numerad,numeta,mbig+1),     &
       propot(numerad,numeta), halfpot(numerad,numeta) )
  usualke=0; usualke2=0;  usualvects=0; usualvals=0; xike=0; proham=0; propot=0; halfpot=0

  if (bornopflag.eq.0) then
     allocate(pro_lplus(numerad,numeta,numerad,numeta,-mbig-1:mbig+1))
     allocate(pro_lminus(numerad,numeta,numerad,numeta,-mbig-1:mbig+1))
     allocate(pro_lsquared(numerad,numeta,numerad,numeta,-mbig-1:mbig+1))
     allocate( myYderivs(numerad,numeta,numerad,numeta, mbig+1) )
     pro_lplus=0; pro_lminus=0; myyderivs=0
  endif
  allocate( rketot(rgridpoints,rgridpoints),    rpoints(rgridpoints),  &
       prolate_derivs(rgridpoints,rgridpoints) , rweights(rgridpoints))
  allocate( etaweights(numeta), etapoints(numeta))
  allocate( etaKE(numeta, numeta,0:mbig+1) )
  allocate( xiweights(xigridpoints), xipoints(xigridpoints), xipoints2d(xinumpoints,2) )
  allocate(proddz(numerad,numeta,numerad,numeta,0:mbig))
  allocate(proddrho(numerad,numeta,numerad,numeta,0:mbig))
  allocate(ddrhopot(numerad,numeta))
  allocate(xydipole(numerad,numeta),zdipole(numerad,numeta))
  allocate(sparseddz_xi_banded(2*bandwidth+1,numerad,lbig+1,mbig+1), &
       sparseddz_eta(  lbig+1,lbig+1,numerad,mbig+1 ) ,  &
       sparseddz_diag(  numerad,lbig+1,mbig+1 )  , &
       sparseddrho_xi_banded(2*bandwidth+1,numerad,lbig+1,mbig+1), &
       sparseddrho_eta(  lbig+1,lbig+1,numerad,mbig+1 ) ,  &
       sparseddrho_diag(  numerad,lbig+1,mbig+1 ) )
  allocate(sparseops_xi_banded(2*bandwidth+1,numerad,lbig+1,mbig+1), &
       sparseops_eta(  lbig+1,lbig+1,numerad,mbig+1), &
       sparseops_diag(  numerad,lbig+1,mbig+1))
  rketot=0; rpoints=0; prolate_derivs=0; rweights=0; etaweights=0; etapoints=0;
  etake=0; xiweights=0; xipoints=0; xipoints2d=0; proddz=0; proddrho=0; ddrhopot=0;
  xydipole=0; sparseddz_xi_banded=0; sparseddz_eta=0; sparseddrho_xi_banded=0;
  sparseddrho_eta=0; sparseddrho_diag=0; sparseops_xi_banded=0; 
  sparseops_eta=0; sparseops_diag=0

  if (bornopflag==0) then
     allocate(&   
          sparseyops_xi_banded(2*bandwidth+1,numerad,lbig+1,mbig+1),   &
          sparseyops_eta(lbig+1,lbig+1,numerad,mbig+1),  &
          sparseyops_diag(numerad,lbig+1,mbig+1),&
          sparselplus_xi_banded(2*bandwidth+1,numerad,lbig+1,-mbig-1:mbig+1),   &
          sparselplus_eta(lbig+1,lbig+1,numerad,-mbig-1:mbig+1),  &
          sparselplus_diag(numerad,lbig+1,-mbig-1:mbig+1),&
          sparselminus_xi_banded(2*bandwidth+1,numerad,lbig+1,-mbig-1:mbig+1),   &
          sparselminus_eta(lbig+1,lbig+1,numerad,-mbig-1:mbig+1),  &
          sparselminus_diag(numerad,lbig+1,-mbig-1:mbig+1))
     sparseyops_xi_banded=0; sparseyops_eta=0; sparseyops_diag=0; 
     sparselplus_xi_banded=0; sparselplus_eta=0; sparselplus_diag=0;
     sparselminus_xi_banded=0; sparselminus_eta=0; sparselminus_diag=0
  endif
  allocate( rmatrix(numerad,numerad,mseriesmax+1,lseriesmax+1)   )
  allocate(ylmvals(0:2*mbig, 1:lbig+1, lseriesmax+1))
  rmatrix=0; ylmvals=0

  call mpibarrier()
  OFLWR "    ...allocated for h2 project. "; CFL
  call mpibarrier()

end subroutine myprojectalloc


subroutine myprojectdealloc()
implicit none
end subroutine myprojectdealloc


module gettwoemod
contains

subroutine get_twoe_new()
  use myparams
  use myprojectmod  
  use utilmod   !! IN PARENT DIRECTORY
  implicit none

  DATAECS, allocatable :: work(:),kearray(:,:,:,:), invkearray(:,:,:,:)
  integer, allocatable :: ipiv(:)
  integer ::  j1a,lsum,deltam,  i,k,ii,m,j, lwork, info
  complex*16,allocatable :: pval(:,:,:), qval(:,:,:), pder(:,:),qder(:,:)
  DATAECS :: surface, wronsk
  real*8 :: dgamma, rsum

  OFLWR "   Calc two electron.";CFL

  allocate(pval(mseriesmax+1,lseriesmax+mseriesmax+1, xigridpoints), &
       qval(mseriesmax+1,lseriesmax+mseriesmax+1, xigridpoints), &
       pder(mseriesmax+1,lseriesmax+mseriesmax+1), &
       qder(mseriesmax+1,lseriesmax+mseriesmax+1))
  pval=0; qval=0; pder=0; qder=0

  lwork=xigridpoints*(lbig+1)*(2*mbig+1) * 10
  allocate(work(lwork), ipiv(lwork))
  work=0; ipiv=0

  allocate( kearray(numerad,numerad,mseriesmax+1,lseriesmax+1), &
       invkearray(numerad,numerad,mseriesmax+1, lseriesmax+1))
  kearray=0; invkearray=0

  do j=1,xigridpoints

#ifdef ECSFLAG
     call clpmn(mseriesmax,mseriesmax,lseriesmax+mseriesmax,real(xipoints(j),8), &
          imag(xipoints(j)-(0.d0,0.000000000000001d0)), pval(:,:,j), pder(:,:))
     call clqmn(mseriesmax,mseriesmax,lseriesmax+mseriesmax,real(xipoints(j),8), &
          imag(xipoints(j)-(0.d0,0.000000000000001d0)), qval(:,:,j),qder(:,:))
#else
     call clpmn(mseriesmax,mseriesmax,lseriesmax+mseriesmax,xipoints(j), 0.d0, pval(:,:,j), pder(:,:))
     call clqmn(mseriesmax,mseriesmax,lseriesmax+mseriesmax,xipoints(j), 0.d0, qval(:,:,j),qder(:,:))
#endif
  enddo

  do i=0,mseriesmax
     do j=0,lseriesmax

! this would be normed
!        pval(i+1,j+i+1,:) = pval(i+1,j+i+1,:) * sqrt(real(2*(j+i)+1)/2.d0*real(floatfac(j+i-i)) / real(floatfac(j+i+i)))
!        qval(i+1,j+i+1,:) = qval(i+1,j+i+1,:) * sqrt(real(2*(j+i)+1)/2.d0*real(floatfac(j+i-i)) / real(floatfac(j+i+i)))

        pval(i+1,j+i+1,:) = pval(i+1,j+i+1,:) * sqrt(floatfac(j+i-i) / floatfac(j+i+i))
        qval(i+1,j+i+1,:) = qval(i+1,j+i+1,:) * sqrt(floatfac(j+i-i) / floatfac(j+i+i))
     enddo
  enddo
  
  do i=0,lseriesmax
     do m=0,mseriesmax
        
        ii=numerad
        kearray(:,:,m+1,i+1) = xike(1:xigridpoints-1,1:xigridpoints-1,m) * (-2.0d0)

!! factor removed from xike.
!    if (reducedflag.eq.1) then
!       kearray(:,:,m+1,i+1) = kearray(:,:,m+1,i+1) / ( totalmass + 2 ) * ( totalmass + 1 )
!    endif

        invkearray(:,:,m+1,i+1) = kearray(:,:,m+1,i+1)
        
        do j=1,numerad
           invkearray(j,j,m+1,i+1) = invkearray(j,j,m+1,i+1) - (i+m)*(i+m+1) 
        enddo
!!$        if (checknan2(invkearray(:,:,m+1,i+1),ii*ii)) then
!!$           print *, "NAN KE!!  Lindex=",i+1," Mindex=", m;call mpistop()
!!$        endif
        
        call MYGETRF(ii,ii,invkearray(:,:,m+1,i+1),ii,ipiv,info)
        if (info/=0) then
           OFLWR "dgetrf info ", info;CFLST
        endif
!!$        if (checknan2(invkearray(:,:,m+1,i+1),ii*ii)) then
!!$           print *, "NAN INVKE 1!!  Lindex=",i+1," Mindex=", m+1;       stop
!!$        endif

        call MYGETRI(ii,invkearray(:,:,m+1,i+1),ii,ipiv,work,lwork,info)
        if (info/=0) then
           OFLWR "dgetri info ", info;CFLST
        endif
!!$        if (checknan2(invkearray(:,:,m+1,i+1),ii*ii)) then
!!$           print *, "NAN INVKE 2!!  Lindex=",i+1," Mindex=", m+1
!!$           stop
!!$        endif

        do j=1,numerad
           do k=1,numerad
              surface =  (-1.d0)**(m+1) * pval(m+1,i+m+1,j) * pval(m+1,i+m+1,k) * qval(m+1,i+m+1,xigridpoints) &
                   / pval(m+1,i+m+1,xigridpoints)
!!$              if (checknan2(surface,1)) then
!!$                 print *, "NAN SURFACE!! ";       stop
!!$              endif
              wronsk = 2**(2*m) * dgamma( real( i+m + m + 2 ,8)/2.d0 ) * dgamma( real( i+m + m + 1 ,8)/2.d0 ) / &
                   dgamma( real( i+m - m + 2 ,8)/2.d0 ) / dgamma( real( i+m - m + 1 ,8)/2.d0 )
              wronsk=wronsk*floatfac(i+m - m)/floatfac(i+m+m)
              rmatrix(j,k,m+1,i+1) = wronsk / (sqrt(xiweights(j)*xiweights(k))) * invkearray(j,k,m+1,i+1) + surface
           enddo
        enddo
!!$        if (checknan2(rmatrix(:,:,m+1,i+1),numerad*numerad)) then
!!$           print *, "NAN RMATRIX!!  Lindex=",i+1," Mindex=", m+1
!!$           stop
!!$        endif
     enddo
  enddo
  
  ylmvals=0.d0
  do deltam   = 0,2*mbig
     do j1a=0,lbig   
        lsum = abs(deltam)-1
        do while (lsum .lt. jacobisummax)
           lsum=lsum+1
           call jacobinormed(lsum-abs(deltam),abs(deltam),abs(deltam),etapoints(j1a+1),rsum)
           ylmvals(deltam,j1a+1,lsum+1-abs(deltam)) = rsum * &
                sqrt( ( 1.d0 - etapoints(j1a+1)**2 ) )**abs(deltam)
        enddo
     enddo
  enddo
!!$  if (checknan2real(ylmvals,(2*mbig+1)*(lbig+1)*lseriesmax+1)) then
!!$     print *, "NAN ylmvals!";     stop
!!$  endif
  deallocate(ipiv,work,kearray,invkearray)
  OFLWR "   ...done.";CFL


  if (twoeattractflag.ne.0) then
     rmatrix(:,:,:,:)=(-1)*rmatrix(:,:,:,:)
  endif

  deallocate(pval,qval,pder,qder)

end subroutine get_twoe_new

end module gettwoemod


subroutine op_lsquaredone(in,out,m2val)   !! right now used just to check operators.
  use myparams
  implicit none
  DATATYPE,intent(out) :: out(numerad,lbig+1)
  DATATYPE,intent(in) :: in(numerad,lbig+1)
  DATATYPE :: temp2(numerad,lbig+1),temp3(numerad,lbig+1)  !!AUTOMATIC
  integer :: m2val

  out=0.d0; temp2=0; temp3=0

  call op_lplusminus_one(in,temp2,1,m2val)      !!lplus
  call op_lplusminus_one(temp2,temp3,2,m2val+1)   !!lminus
  out=out+temp3(:,:)*(0.5d0)
  call op_lplusminus_one(in,temp2,2,m2val)  
  call op_lplusminus_one(temp2,temp3,1,m2val-1)
  out=out+temp3(:,:)*(0.5d0)
  out(:,:)=out(:,:)+m2val**2*in(:,:)
end subroutine


!! kind=1: lplus; 2, lminus 

#ifndef REALGO
#define XXMVXX zgemv 
#define XXBBXX zgbmv 
#else
#define XXMVXX dgemv 
#define XXBBXX dgbmv 
#endif

subroutine op_lplusminus_one(in,out,inkind,m2val)
  use myparams
  use myprojectmod
  implicit none
  integer, intent(in) :: m2val,inkind
  DATATYPE,intent(in) :: in(numerad,lbig+1)
  DATATYPE,intent(out) :: out(numerad,lbig+1)
  DATATYPE :: work(lbig+1) 
  integer :: ieta,i,ixi

  out=0.d0; work=0

  if ( (inkind/=1) .and. (inkind/=2) ) then
     OFLWR "kind error", inkind; CFLST
  endif

  do ixi=1,numerad
     select case(inkind)  
     case (1)           
        call XXMVXX('N',lbig+1,lbig+1,(1.d0,0.d0),sparselplus_eta(:,:,ixi,m2val),lbig+1,&
             in(ixi,:),1,(0.d0,0.d0), work, 1)
     case (2)
        call XXMVXX('N',lbig+1,lbig+1,(1.d0,0.d0),sparselminus_eta(:,:,ixi,m2val),lbig+1,&
             in(ixi,:),1,(0.d0,0.d0), work, 1)
     case default
        print *, "AAUGH!!!"
        call mpistop()
     end select
     out(ixi,:)= out(ixi,:) + work(1:lbig+1) 
  enddo

  i=2*bandwidth+1
  do ieta=1,lbig+1
     select case(inkind)
     case (1)           
        call XXBBXX('N',numerad,numerad,bandwidth,bandwidth,(1.d0,0.d0),&
             sparselplus_xi_banded(:,:,ieta,m2val),i, in(:,ieta),1,(1.d0,0.d0), out(:,ieta), 1)
     case (2)
        call XXBBXX('N',numerad,numerad,bandwidth,bandwidth,(1.d0,0.d0),&
             sparselminus_xi_banded(:,:,ieta,m2val),i, in(:,ieta),1,(1.d0,0.d0), out(:,ieta), 1)
     case default
        print *, "AAUGH!!!";        call mpistop()
     end select
  enddo
  
  select case(inkind)
  case (1)
     out(:,:) = out(:,:) - sparselplus_diag(:,:,m2val) * in(:,:) 
  case (2)
     out(:,:) = out(:,:) - sparselminus_diag(:,:,m2val) * in(:,:)
  case default
     print *, "AAUGH!!!";     call mpistop()
  end select
  
end subroutine op_lplusminus_one
  

subroutine op_yderiv(howmany,in,out)
  use myprojectmod
  use myparams
  implicit none
  integer,intent(in) :: howmany
  DATATYPE, intent(in) :: in(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE, intent(out) :: out(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE :: work(lbig+1),work2(lbig+1)  !!AUTOMATIC
  integer ::  ixi, ieta, i,m2val,ii

  out=0.d0;

  if (bornopflag==1) then
     return
  endif

 work=0; work2=0

  do ii=1,howmany
  do m2val=-mbig,mbig
     do ixi=1,numerad
        work2=in(ixi,:,m2val,ii)
        call XXMVXX('N',lbig+1,lbig+1,(1.d0,0.d0),sparseyops_eta(:,:,ixi,abs(m2val)+1),lbig+1,&
             work2,1,(0.d0,0.d0), work, 1)
        out(ixi,:,m2val,ii)= out(ixi,:,m2val,ii) + work(1:lbig+1)  
     enddo

     i=2*bandwidth+1
     do ieta=1,lbig+1
        call XXBBXX('N',numerad,numerad,bandwidth,bandwidth,(1.d0,0.d0),&
             sparseyops_xi_banded(:,:,ieta,abs(m2val)+1),i, &
             in(:,ieta,m2val,ii),1,(1.d0,0.d0), out(:,ieta,m2val,ii), 1)
     enddo
     out(:,:,m2val,ii) = out(:,:,m2val,ii) - &
          sparseyops_diag(:,:,abs(m2val)+1) * in(:,:,m2val,ii) 
  enddo
  enddo

end subroutine op_yderiv


