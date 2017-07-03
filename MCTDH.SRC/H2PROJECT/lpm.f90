
#include "Definitions.INC"

module lpmsubmod
contains

subroutine pro_lpm()
  use myparams
  use myprojectmod
  use myprojecttempmod
  implicit none

  integer :: i,j,k,m, mvalue,ii,l
  DATAECS :: sum
#define CHECxxKL
#ifdef CHECKL
  complex*16, allocatable :: orb1(:,:), orb2(:,:)
  DATAECS, allocatable :: pro_ltemp(:,:,:), pro_tempv(:)
  DATAECS :: csum, simplecdot
#endif

  pro_lop=0.d0;  pro_lplus=0.d0;  pro_lminus=0.d0

  if (bornopflag == 1) then
     return
  endif
  do ii=0,1
     do i=1,numerad
        do j=1,numerad
           if ((abs(xideg1(i,j,ii)+xideg1(j,i,ii))/( abs(xideg1(i,j,ii)+xideg1(j,i,ii)) +1 )).gt.1d-4) then
              print *, "xideg1 warn ", i,j, xideg1(i,j,ii),xideg1(j,i,ii)
!           stop
           endif
        enddo
     end do
     
     do i=1,numerad
        do j=1,numerad
           if (( ( abs(xideg2(i,j,ii)+xideg2(j,i,ii))/(abs(xideg2(i,j,ii)-xideg2(j,i,ii))+1)) .gt.1d-4)) then
              print *, "xideg2 warn ", i,j, xideg2(i,j,ii),xideg2(j,i,ii)
!           stop
           endif
        enddo
     end do
     
     do i=1,numeta
        do j=1,numeta
           if (( ( abs(etadeg1(i,j,ii)+etadeg1(j,i,ii))/(abs(etadeg1(i,j,ii)-etadeg1(j,i,ii))+1).gt.1d-4))) then
              print *, "etadeg1 warn ", i,j, etadeg1(i,j,ii),etadeg1(j,i,ii)
!           stop
           endif
        enddo
     end do
     
     do i=1,numeta
        do j=1,numeta
           if ((abs(etadeg2(i,j,ii)+etadeg2(j,i,ii))/(abs(etadeg2(i,j,ii)-etadeg2(j,i,ii))+1).gt.1d-4)) then
              print *, "etadeg2 warn ", i,j, etadeg2(i,j,ii),etadeg2(j,i,ii)
!           stop
           endif
        enddo
     end do
  enddo
  
  do i=1,xigridpoints
     do k=1,xigridpoints
        do j=1,numeta
           l=j
           sum= sqrt((1-etapoints(j)**2)/(xipoints(i)**2-1))* &
                ( etapoints(j) * xideg1(i,k,0) + mass_asym*xideg2(i,k,0) )
           pro_lop(i,j,k,l,0) = sum
           sum= (1-etapoints(j)**2)*sqrt(1d0/(1-etapoints(l)**2)/(xipoints(k)**2-1))* &
                ( etapoints(j) * xideg1(i,k,0) + mass_asym*xideg2(i,k,0) )
           pro_lop(i,j,k,l,1) = sum
        enddo
     enddo
  enddo

  do i=1,xigridpoints
     k=i
     do j=1,numeta
        do l=1,numeta
           sum=  sqrt((xipoints(i)**2-1)/(1-etapoints(j)**2))* &
                ( xipoints(i) * etadeg1(j,l,0) + mass_asym*etadeg2(j,l,0) )
           pro_lop(i,j,k,l,0) = pro_lop(i,j,k,l,0) + sum
           sum=  (xipoints(i)**2-1)*sqrt(1d0/(xipoints(k)**2-1)/(1-etapoints(l)**2))* &
                ( xipoints(i) * etadeg1(j,l,0) + mass_asym*etadeg2(j,l,0) )
           pro_lop(i,j,k,l,1) = pro_lop(i,j,k,l,1) + sum
        enddo
     enddo
  enddo

  do ii=0,1
     do i=1,numerad !!xigridpoints
        do j=1,numeta
           do k=1,numerad !!xigridpoints
              do l=1,numeta
                 if (abs(pro_lop(i,j,k,l,ii)+pro_lop(k,l,i,j,mod(ii+1,2))).gt.1.d-7) then
                    print *, "LOP ERR", i,j,k,l,ii, pro_lop(i,j,k,l,ii),pro_lop(k,l,i,j,mod(ii+1,2))
                    stop
                 endif
              enddo
           enddo
        enddo
     enddo
  enddo
  
  do m=1,xigridpoints
     do i=1,numeta
        sum=  sqrt( (xipoints(m)**2-etapoints(i)**2) ) 
        pro_lop(:,:,m,i,:)=pro_lop(:,:,m,i,:)/sum
        pro_lop(m,i,:,:,:)=pro_lop(m,i,:,:,:)/sum
     enddo
  enddo

  do j=1,numerad
     do k=1,numeta
        pro_ldiag(j,k)= (xipoints(j) * etapoints(k)+ mass_asym) / sqrt( ( xipoints(j)**2 - 1.d0 ) * (1.d0 - etapoints(k)**2) ) 
     enddo
  enddo
  
!  do ii=0,1
!     do i=1,xigridpoints
!        do j=1,numeta
!           do k=1,xigridpoints
!              do l=1,numeta
!                 if (abs(imag(pro_lop(i,j,k,l,ii))).gt.1d-10) then
!                    print *, "err lop", pro_lop(i,j,k,l,ii), i,j,k,l
!                    stop
!                 endif
!              enddo
!           enddo
!        enddo
!     enddo
!  enddo
  

!! pro_lplus and minus

   do mvalue = -(mbig+1),(mbig+1)
      pro_lplus(:,:,:,:,mvalue)=pro_lop(1:numerad,:,1:numerad,:,mod(mvalue+100,2))
      pro_lminus(:,:,:,:,mvalue)=(-1)*pro_lop(1:numerad,:,1:numerad,:,mod(mvalue+100,2))
      do j=1,numeta
         do i=1,numerad
            pro_lplus(i,j,i,j,mvalue)  = pro_lplus(i,j,i,j,mvalue)  - (mvalue+mod(mvalue+99,2)) * pro_ldiag(i,j)    !! lplus
            pro_lminus(i,j,i,j,mvalue) = pro_lminus(i,j,i,j,mvalue) - (mvalue-mod(mvalue+99,2)) * pro_ldiag(i,j)    !! lminus
         enddo
      enddo
   enddo

   pro_lsquared=0.d0
   do mvalue=-mbig,mbig
      pro_lsquared(1:numerad,:,1:numerad,:,mvalue) = RESHAPE(matmul( &
           RESHAPE(pro_lminus(1:numerad,:,1:numerad,:,mvalue+1),(/numeta*numerad,numeta*numerad/)), &
           RESHAPE(pro_lplus(1:numerad,:,1:numerad,:,mvalue),(/numeta*numerad,numeta*numerad/)) &
           ),(/ numerad,numeta,numerad,numeta /))*(0.5d0)
      pro_lsquared(1:numerad,:,1:numerad,:,mvalue) = pro_lsquared(1:numerad,:,1:numerad,:,mvalue) + RESHAPE(matmul( &
           RESHAPE(pro_lplus(1:numerad,:,1:numerad,:,mvalue-1),(/numeta*numerad,numeta*numerad/)), &
           RESHAPE(pro_lminus(1:numerad,:,1:numerad,:,mvalue),(/numeta*numerad,numeta*numerad/)) &
           ),(/ numerad,numeta,numerad,numeta /))*(0.5d0)

!      do i=1,numerad
!         do j=1,numeta
!            do k=1,numerad
!               do l=1,numeta
!                  if (abs(imag(pro_lsquared(i,j,k,l,mvalue))).gt.1d-10) then
!                     print *, "err lsq", pro_lsquared(i,j,k,l,mvalue), i,j,k,l
!                     stop
!                  endif
!               enddo
!            enddo
!         enddo
!      enddo
      
      do j=1,numeta
         do i=1,numerad
            pro_lsquared(i,j,i,j,mvalue) = pro_lsquared(i,j,i,j,mvalue) + mvalue**2 
         enddo
      enddo
   enddo  !! MVALUE

   if (bornopflag.eq.0) then
      sparselplus_xi_banded(:,:,:,:)=0.d0;  sparselplus_eta(:,:,:,:)=0.d0;  sparselplus_diag(:,:,:)=0.d0
      sparselminus_xi_banded(:,:,:,:)=0.d0;  sparselminus_eta(:,:,:,:)=0.d0;  sparselminus_diag(:,:,:)=0.d0
   endif
   if (bornopflag/=1) then
      do mvalue=-mbig-1,mbig+1
         do j=1,lbig+1
            do i=1,numerad
               do k=max(1,i-bandwidth),min(numerad,i+bandwidth)
                  sparselplus_xi_banded(k-i+bandwidth+1,i,j,mvalue) = pro_lplus(k,j,i,j,mvalue)
                  sparselminus_xi_banded(k-i+bandwidth+1,i,j,mvalue) = pro_lminus(k,j,i,j,mvalue) 
               enddo
               do k=1,lbig+1
                  sparselplus_eta(k,j,i,mvalue) = pro_lplus(i,k,i,j,mvalue)
                  sparselminus_eta(k,j,i,mvalue) = pro_lminus(i,k,i,j,mvalue) 
               enddo
               sparselplus_diag(i,j,mvalue) = pro_lplus(i,j,i,j,mvalue) 
               sparselminus_diag(i,j,mvalue) = pro_lminus(i,j,i,j,mvalue) 
            enddo
         enddo
      enddo
   endif

!! CHECK L OPERATORS!!    new lplus, lminus make better l-squared. indeed, replace, not just here temporarily to check.

#ifdef CHECKL
   
   print *;         print *, "Doing LSQ CHECK...", numerad*numeta;         print *
   allocate(orb1(numerad,lbig+1), orb2(numerad,lbig+1))
   allocate(pro_ltemp(numerad,numeta,numeta*numerad),  pro_tempv(numeta*numerad))
   do ii=-mbig,mbig
      print *
      print *, " *****   MVALUE ", ii , " *****"
      print *
#ifdef ECSFLAG
      call get_eigen_cnorm(pro_lsquared(:,:,:,:,ii),numerad*numeta,numerad*numeta,pro_ltemp, pro_tempv)
#else
      call get_eigen_two(pro_lsquared(:,:,:,:,ii),numerad*numeta,numerad*numeta,pro_ltemp, pro_tempv)
#endif
      do i=1,min(90,numerad*numeta)
         
         orb1(:,:)=pro_ltemp(1:numerad,:,i)
         call op_lsquaredone(orb1,orb2,ii)
         csum=simplecdot(orb1,orb2,edim)
         
!!            csum=DOT_PRODUCT(conjg(RESHAPE(orb1(:,:,:),(/edim*(2*mbig-mval+1)/))), RESHAPE(orb2(:,:,:),(/edim*(2*mbig-mval+1)/)))   !! c-product
!            print *, "lap"
         print *, pro_tempv(i),csum 
      enddo
      print *;         print *, "DONE LSQ CHECK...";         print *
   enddo
   deallocate(orb1,orb2);     deallocate(pro_ltemp,pro_tempv)
   call mpistop()
#endif
   
end subroutine pro_lpm


subroutine get_both()
  use myparams
  use myprojectmod
  use myprojecttempmod
  implicit none
  integer :: i,j,k,ii,l

#define CHECKxxY
#ifdef CHECKY
  DATAECS, allocatable :: pro_ltemp(:,:,:,:), pro_lvect(:)
#endif

  pro_both=0.d0
  if ((bornopflag == 1)) then
     return
  endif

!! operators are negative definite...  in the end pro_both is subtracted.

  do j=1,numeta
     pro_both(:,j,:,j,:) = pro_both(:,j,:,j,:) &
          + xi_fourth(1:numerad,1:numerad,:) &
          + (etapoints(j)**2 + mass_asym**2) * xi_second(1:numerad,1:numerad,:) &
          + 2 * mass_asym * etapoints(j) * xi_third(1:numerad,1:numerad,:)
  enddo
  do j=1,numerad
     pro_both(j,:,j,:,:) = pro_both(j,:,j,:,:) &
          - eta_fourth(:,:,:) &
          + (xipoints(j)**2 + mass_asym**2) * eta_second(:,:,:) &
          - 2 * mass_asym * xipoints(j) * eta_third(:,:,:)
  enddo
  do i=1,numeta
     do j=1,numerad
        pro_both(j,i,:,:,:) = pro_both(j,i,:,:,:)/sqrt(xipoints(j)**2 - etapoints(i)**2)
        pro_both(:,:,j,i,:) = pro_both(:,:,j,i,:)/sqrt(xipoints(j)**2 - etapoints(i)**2)
     enddo
  enddo
  do i=1,numeta
     do j=1,numerad
        pro_both(j,i,j,i,:) = pro_both(j,i,j,i,:) + 9d0/4d0
     enddo
  enddo
  do ii=0,1
     do i=1,numerad
        do j=1,numeta
           do k=1,numerad
              do l=1,numeta
                 if (abs(pro_both(i,j,k,l,ii)-pro_both(k,l,i,j,ii)).gt.1.d-5) then
                    print *, "BOTH ERR", i,j,k,l,ii, pro_both(i,j,k,l,ii),pro_both(k,l,i,j,ii)
                    stop
                 endif
              enddo
           enddo
        enddo
     enddo
  enddo
  
#ifdef CHECKY

  allocate(pro_ltemp(numerad,numeta,numerad,numeta),pro_lvect(numerad*numeta))
  do ii=0,1
     print *;         print *, "Doing both CHECK...", numerad*numeta, ii ;         print *

#ifdef ECSFLAG
     call get_eigen_cnorm(pro_both(:,:,:,:,ii),numerad*numeta,numerad*numeta,pro_ltemp, pro_lvect)
#else
     call get_eigen_two(pro_both(:,:,:,:,ii),numerad*numeta,numerad*numeta,pro_ltemp, pro_lvect)
#endif
     
     do i=1,min(90,numerad*numeta)
        print *, pro_lvect(i)
     enddo
     print *;         print *, "DONE both CHECK...";         print *
  enddo
  deallocate(pro_ltemp, pro_lvect)
  call mpistop()

#endif

end subroutine get_both

end module lpmsubmod


