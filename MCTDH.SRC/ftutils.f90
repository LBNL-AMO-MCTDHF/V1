
#include "Definitions.INC"


module half_ft_mod
  ! report halfft_fac times as many pos-frequency results as input data
  integer, parameter :: halfft_fac   = 3
  ! oversample by this factor
  integer, parameter :: NN           = 7
  ! halfft_fac should be less than NN/2
end module half_ft_mod


module ftutilmod
contains
  
  subroutine complexdiff_prev(size,in,out,diffoption)  
    use fileptrmod
    implicit none
    integer, intent(in) :: size, diffoption ! diffoption = ftdiff hardwired to 3
    complex*16, intent(in) :: in(size)
    complex*16, intent(out) :: out(size)
    integer :: i,jj

    out(:)=0d0
    do i=2,size-1
       jj=min(min(i-1,size-i),4)
       select case(jj)
       case(1)
          out(i)= &
               - 1d0/2d0 * in(i-1) &
               + 1d0/2d0 * in(i+1)
       case(2)
          out(i)= &
               1d0/12d0 * in(i-2) &
               - 2d0/3d0 * in(i-1) &
               + 2d0/3d0 * in(i+1) &
               - 1d0/12d0 * in(i+2)
       case(3)
          out(i)= &
               - 1d0/60d0 * in(i-3) &
               + 3d0/20d0 * in(i-2) &
               - 3d0/4d0 * in(i-1) &
               + 3d0/4d0 * in(i+1) &
               - 3d0/20d0 * in(i+2) &
               + 1d0/60d0 * in(i+3)
       case(4)
          out(i)= &
               1d0/280d0 * in(i-4) &
               - 4d0/105d0 * in(i-3) &
               + 1d0/5d0 * in(i-2) &
               - 4d0/5d0 * in(i-1) &
               + 4d0/5d0 * in(i+1) &
               - 1d0/5d0 * in(i+2) &
               + 4d0/105d0 * in(i+3) &
               - 1d0/280d0 * in(i+4)
       end select
    end do

  end subroutine complexdiff_prev


  subroutine banded_diff_mat(size,diffpower,bw,diffmat,lda)
    use fileptrmod
    implicit none
    integer, intent(in) :: size, diffpower, bw, lda
    real*8, intent(out) :: diffmat(-bw:bw,size)
    integer :: ipiv(2*bw+1)  !! AUTOMATIC
    real*8, allocatable :: xmat(:,:),vec(:),ymat(:,:),xmat0(:,:),ymat2(:,:),&
         smat(:,:), pts(:)
    real*8 :: fac
    integer :: tbw,ii,jj,neqn,factorial,imin,imax,info,nn

    ! note:  bw = difforder
    if (diffpower==0) then
       OFLWR "Programmer failure: diffpower=0 banded_diff_mat not necessary"; CFLST
    endif
    if (diffpower<0 .or. bw.lt.1) then
       OFLWR "ERROR diffpower or difforder in banded_diff_mat ",diffpower,bw; CFLST
    endif

    tbw = 2*bw + 1
    neqn = tbw

    if (lda.ne.tbw) then
       OFLWR "ERROR banded_diff_mat programmer failure",lda,tbw; CFLST
    endif

    fac = bw**((bw-1d0)/bw)
    factorial = 1
    do ii=2,diffpower
       factorial = factorial * ii
    enddo

    fac = 1  ! TEMP

    allocate(vec(neqn), xmat(neqn,-bw:bw), ymat(-bw:bw,-bw:bw), &
         xmat0(neqn,-bw:bw), ymat2(-bw:bw,-bw:bw), &
         smat(-bw:bw,-bw:bw), pts(-bw:bw) )

    diffmat(:,:) = 0  ; vec(:) = 0;
    xmat(:,:) = 0; ymat(:,:) = 0; xmat0(:,:) = 0
    ymat2(:,:) = 0; smat(:,:) = 0; pts(:) = 0

    do ii = -bw,bw
       pts(ii) = ii / fac
    enddo
    do ii = 1, neqn
       xmat0(ii,:) = pts(:)**(ii-1)
    enddo
    do ii = -bw , bw
       xmat(:,:) = xmat0(:,:)
       imin = max(-bw,-bw-ii)
       imax = min(bw,bw-ii)
       nn = imax-imin+1
       vec(:) = 0
       vec(diffpower+1) = 1
       call dgesv(nn,1,xmat(1,imin),tbw,ipiv,vec,neqn,info)
       if (info.ne.0)then
          OFLWR "INFO DGESV IN banded_diff_mat ",info; CFLST
       endif
       vec = vec * factorial / fac**diffpower
       ymat(ii,imin:imax) = vec(1:nn)
       imin = max(-bw,-bw+ii)
       imax = min(bw,bw+ii)
       ymat2(ii,imin:imax) = vec(1:nn)
    enddo

    if (1==0) then
       print *, "YMAT2:  ", fac, diffpower
       do ii = -bw,bw
          print '(100F10.5)', ymat2(ii,:)
       enddo
       print *, "-----------end diffmat"
       print '(9F10.5)', matmul(ymat,transpose(xmat0))
       !print *, "TEMPSTOP"
       !CFLST
    endif

    do jj = -bw, -1
       ii = jj + bw + 1
       diffmat(-bw-jj:bw,ii) = ymat2(jj,-bw:bw+jj)
    enddo
    do ii=bw+1,size-bw
       diffmat(:,ii) = ymat2(0,:)
    enddo
    do jj= 1,bw
       ii = jj - bw + size
       diffmat(-bw:bw-jj,ii) = ymat2(jj,-bw+jj:bw)
    enddo

    if (1==0) then
       print *, "DMAT - T:  ", fac, diffpower
       do ii = 1,size
          print '(100F10.5)', diffmat(:,ii)
       enddo

       !print *, "DMAT:  ", fac, diffpower
       !do ii = -bw, bw
       !   print '(100F10.5)', diffmat(ii,1:10)
       !enddo

       print *, "-----------end diffmat"
       !print *, "TEMPSTOP"
       !CFLST
    endif

    deallocate(vec,xmat,ymat,xmat0,ymat2,smat,pts)

  end subroutine banded_diff_mat


  !! NOW DOES OTHER DERIVATIVES.. BEFORE WAS ONLY 1ST DERIVATIVE
  !!   and ftdiff selected different options for it..  now
  !!   we select the derivative order with diffpower
  !!   difforder is the order of the finite difference approximation
  !!
  subroutine complexdiff(size,in,out,howmany,diffpower,difforder)
    use fileptrmod
    implicit none
    integer, intent(in) :: size, howmany,diffpower, difforder
    complex*16, intent(in) :: in(size,howmany)
    complex*16, intent(out) :: out(size,howmany)
    real*8, allocatable :: bandmat(:,:)
    real*8 :: rein(size),reout(size),imin(size),imout(size) ! AUTOMATIC
    integer :: ii,tbw,bw

    if (diffpower<0) then
       OFLWR "error diffpower >= 0 please ", diffpower; CFLST
    endif
    if (diffpower==0) then
       out(:,:) = in(:,:)
       return
    endif
    if (difforder<1) then
       OFLWR "error, difforder > 0 please ",difforder; CFLST
    endif
    bw      = difforder
    tbw     = 2*difforder+1
    allocate(bandmat(tbw,size))

    call banded_diff_mat(size,diffpower,bw,bandmat,tbw)

    out(:,:) = 0d0

    do ii = 1,howmany
       rein(:) = real(in(:,ii),8)
       imin(:) = imag(in(:,ii))

       !   do jj=1,size
       !      rein(jj) = (jj)**3
       !   enddo

       call DGBMV('T',size,size,bw,bw,1d0,bandmat,tbw,rein,1,0d0,reout,1)

       !   print *, "DIFFVEC"
       !   print '(3F14.0)', reout(:)
       !   print *, "TEMPSTOP"
       !   CFLST

       call DGBMV('T',size,size,bw,bw,0d0,bandmat,tbw,imin,1,0d0,imout,1)
       out(:,ii)  = reout(:) + (0d0,1d0) * imout(:)
    enddo
    deallocate(bandmat)

  end subroutine complexdiff


!!!!!!!!!!  Half FT  !!!!!!!!!!!

  function half_ft_size(size)
    use half_ft_mod
    implicit none
    integer,intent(in) :: size
    integer :: half_ft_size
    half_ft_size = halfft_fac * size;
  end function half_ft_size


  subroutine half_ft_estep(size,estart,estep)
    use half_ft_mod
    implicit none
    integer,intent(in) :: size
    real*8, intent(out) :: estart,estep
    real*8, parameter :: twopi = 6.28318530717958647688d0
    estep = twopi / size / 2 / NN
    estart = estep / 2
  end subroutine half_ft_estep


  subroutine half_ft_evals(size,evals,outsize)
    use fileptrmod
    use half_ft_mod
    implicit none
    integer,intent(in) :: size,outsize
    real*8, intent(out) :: evals(outsize)
    real*8 :: estart,estep
    integer :: i
    if (outsize.ne.half_ft_size(size)) then
       OFLWR "ack half_ft_evals",outsize,half_ft_size(size); CFLST
    endif
    call half_ft_estep(size,estart,estep)
    evals(:) = 0
    do i=1,outsize
       evals(i) = estart + estep*(i-1)
    enddo
  end subroutine half_ft_evals


  subroutine half_ft_wrap_diff(size,in,out,outsize,howmany,diffpower,difforder)
    use fileptrmod
    use half_ft_mod
    implicit none  
    integer, intent(in) :: size,outsize,diffpower,difforder,howmany
    complex*16, intent(in) :: in(size,howmany)
    complex*16, intent(out) :: out(outsize,howmany)
    complex*16,allocatable :: work(:,:)
    real*8, allocatable :: evals(:)
    integer :: i

    if (outsize.ne.half_ft_size(size)) then
       OFLWR "progammer error 99801"; CFLST
    endif

    if (diffpower.eq.0) then
       call half_ft_wrap(size,in,out,outsize,howmany)
    else
       allocate(work(size,howmany),evals(outsize));
       work=0; evals=0

       call half_ft_evals(size,evals,outsize)

       call complexdiff(size,in,work,howmany,diffpower,difforder)

       call half_ft_wrap(size,work,out,outsize,howmany)

       do i=1,howmany
          out(:,i) = out(:,i) * ( (0d0,-1d0) / evals(:) )**diffpower
       enddo

       deallocate(work,evals)
    endif

  end subroutine half_ft_wrap_diff



  subroutine half_ft_wrap(size,in,out,outsize,howmany)
    use fileptrmod
    use half_ft_mod
    implicit none
    integer, intent(in) :: size,outsize,howmany
    complex*16, intent(in) :: in(size,howmany)
    complex*16, intent(out) :: out(outsize,howmany)
    complex*16,allocatable :: vv(:,:), csave(:)
    real*8, allocatable :: rsave(:)
    integer :: ifac(15)  !! AUTOMATIC
    integer :: bigsize,ii

    if (size.gt.6000) then
       OFLWR "ARE YOU SURE?  big size, programmer checkme ft ",size; CFLST
    endif
    if (outsize > 2*NN*size) then
       OFLWR "ack bad outsize,NN,size: ",outsize,NN,size; CFLST
    endif
    if (half_ft_size(size) .ne. outsize) then
       OFLWR "ack programmer failure 543"; CFLST
    endif

    bigsize           = 8*NN*size

    allocate(vv(2,bigsize/2))
    vv(:,:) = 0; 
    allocate(rsave(2*bigsize),csave(2*bigsize))
    rsave(:) = 0; csave(:) = 0
    ifac(:) = 0

    call cffti1(bigsize,rsave,ifac)

    out(:,:) = 0;
    do ii             = 1,howmany
       vv(:,:)        = 0;
       vv(2,1:size)   = in(:,ii);
       call cfftf1(bigsize,vv,csave,rsave,ifac)
       out(:,ii)    = vv(2,1:outsize);
    enddo

    deallocate(rsave,csave)
    deallocate(vv)

  end subroutine half_ft_wrap


!!!!!!!!!  ZFFT (DFFTPACK) FT subroutines   !!!!!!!!!!



  function zfftf_size(size)
    implicit none
    integer,intent(in) :: size
    integer :: zfftf_size
    zfftf_size = size;
  end function zfftf_size


  subroutine zfftf_evals(size,evals,outsize)
    use fileptrmod
    implicit none
    integer,intent(in) :: size,outsize
    real*8, intent(out) :: evals(size)
    real*8 :: estart,estep
    integer :: i
    if (outsize.ne.zfftf_size(size)) then
       OFLWR "ack half_ft_evals",outsize,half_ft_size(size); CFLST
    endif
    call zfftf_estep(size,estart,estep)
    evals(:) = 0
    do i=1,size
       evals(i) = estart + estep*(i-1)
    enddo
  end subroutine zfftf_evals


  subroutine zfftf_estep(size,estart,estep)
    implicit none
    integer,intent(in) :: size
    real*8, intent(out) :: estart,estep
    real*8, parameter :: twopi = 6.28318530717958647688d0
    estep = twopi / size
    estart = 0
  end subroutine zfftf_estep


  subroutine zfftf_wrap_diff(size,inout,diffpower,difforder)
    implicit none
    integer, intent(in) :: size,diffpower,difforder
    complex*16, intent(inout) :: inout(size)
    complex*16,allocatable :: work(:)
    real*8, allocatable :: evals(:)

    if (diffpower.eq.0) then
       call zfftf_wrap(size,inout)
    else
       allocate(work(size),evals(size));
       work=0; evals=0
       call zfftf_evals(size,evals,size)

       !! NOW DOES OTHER DERIVATIVES.. BEFORE WAS ONLY 1ST DERIVATIVE
       !!   and ftdiff selected different options for it..  now
       !!   we select the derivative order with diffpower
       !!   difforder is the order of the finite difference approximation
       !!
       call complexdiff(size,inout,work,1,diffpower,difforder)

       call zfftf_wrap(size,work)

       inout(:) = work(:) * ( (0d0,-1d0) / evals(:) )**diffpower

       deallocate(work,evals)
    endif

  end subroutine zfftf_wrap_diff


  !  *********       edits 08-2021                  *****
  !  zfft functions use bad casting of wsave
  !  they are just little wrappers for cfft1 functions
  !  *********   now using those instead 08-2021    *****

  subroutine zfftf_wrap(size,inout)
    implicit none
    integer, intent(in) :: size
    complex*16, intent(inout) :: inout(size)
    complex*16,allocatable :: csave(:)
    real*8, allocatable :: rsave(:)
    integer :: ifac(15)  !! AUTOMATIC
    allocate(rsave(2*size),csave(2*size))
    rsave(:) = 0; csave(:) = 0; ifac(:) = 0
    call cffti1(size,rsave,ifac)
    call cfftf1(size,inout,csave,rsave,ifac)
    deallocate(rsave,csave)
  end subroutine zfftf_wrap


  subroutine zfftb_wrap(size,inout)
    implicit none
    integer, intent(in) :: size
    complex*16, intent(inout) :: inout(size)
    complex*16,allocatable :: csave(:)
    real*8, allocatable :: rsave(:)
    integer :: ifac(15)  !! AUTOMATIC
    allocate(rsave(2*size),csave(2*size))
    rsave(:) = 0; csave(:) = 0; ifac(:) = 0
    call cffti1(size,rsave,ifac)
    call cfftb1(size,inout,csave,rsave,ifac)
    deallocate(rsave,csave)
  end subroutine zfftb_wrap


  !subroutine zfftb_wrap(size,inout)
  !  implicit none
  !  integer, intent(in) :: size
  !  complex*16, intent(inout) :: inout(size)
  !  complex*16,allocatable :: wsave(:)
  !  allocate(wsave(4*size+15))
  !  wsave(:)=0d0
  !  call zffti(size,wsave)          ! no bad casting of wsave 08-2021
  !  call zfftb(size,inout,wsave)
  !  deallocate(wsave)
  !end subroutine zfftb_wrap

end module ftutilmod
