
module ct_mpimod
  implicit none
  integer :: myrank = -1
  integer :: nprocs = -1
  integer :: mpifileptr=6
  integer :: ctparopt=999    !! 0 = sendrecv  1 = SUMMA
end module ct_mpimod


subroutine ctset(in_ctparopt)
  use ct_mpimod
  implicit none
  integer :: in_ctparopt
  if (in_ctparopt.ne.0.and.in_ctparopt.ne.1) then
     print *, "COOLEY TUKEY OPTION ERROR (fft_ct_paropt)", in_ctparopt; call mpistop()
  endif
  ctparopt=in_ctparopt
  call getmyranknprocs(myrank,nprocs)
  if (myrank.eq.1) then
     mpifileptr=6
  else
     mpifileptr=987
     open(mpifileptr,file="/dev/null", status="unknown")
  endif
end subroutine ctset


subroutine ctcheck(int1)
  use ct_mpimod
  implicit none
  integer :: int1,int2
  int2=(int1/nprocs)*nprocs
  if (int2.ne.int1) then
     write(mpifileptr,*) "ACK DIVISOR", int1,nprocs,int2; call mpistop()
  endif
end subroutine ctcheck


!! INVERSE OF cooleytukey_outofplace_mpi  EXCEPT FOR DIVISON

recursive subroutine cooleytukey3d_outofplace_backward_mpi(intranspose,out,dim,times,howmany)
  use ct_mpimod
  implicit none
  integer, intent(in) :: dim,howmany
  integer, intent(inout) :: times(8)
  integer :: atime,btime
  complex*16, intent(in) :: intranspose(dim,dim,1,dim/nprocs,howmany)
  complex*16, intent(out) :: out(dim,dim,dim/nprocs,1,howmany)
  complex*16 ::  intransconjg(dim,dim,1,dim/nprocs,howmany), & !! AUTOMATIC
       outconjg(dim,dim,dim/nprocs,1,howmany)

  call ctcheck(dim)

  call myclock(atime)
  intransconjg(:,:,:,:,:)=CONJG(intranspose(:,:,:,:,:))
  call myclock(btime); times(1)=times(1)+btime-atime

  call cooleytukey3d_outofplaceinput_mpi(intransconjg,outconjg,dim,times,howmany)

  call myclock(atime)
  out(:,:,:,:,:)=CONJG(outconjg(:,:,:,:,:))     !!! IN CIRC3D_SUB  /dim**3
  call myclock(btime); times(1)=times(1)+btime-atime

end subroutine cooleytukey3d_outofplace_backward_mpi


!! fourier transform with OUT-OF-PLACE OUTPUT. 

!! times(1) saved for conjugate, backward, not used here
!! times(2) gettwiddle
!! times(3) 1d ft-slowindex
!! times(4) raise twiddle factor to power
!! times(5) multiply
!! times(6) 3d f.t.

recursive subroutine cooleytukey3d_outofplace_mpi(in,outtrans,dim,times,howmany)
  use ct_mpimod
  implicit none
  integer, intent(in) :: dim,howmany
  integer, intent(inout) :: times(8)
    complex*16, intent(in) :: in(dim,dim,dim/nprocs,1,howmany)
  complex*16, intent(out) :: outtrans(dim,dim,1,dim/nprocs,howmany)
  complex*16 :: twiddle1(dim/nprocs),tt1(dim/nprocs), &  !! AUTOMATIC
       tempout(dim,dim,dim/nprocs,1,howmany),      outtemp(dim,dim,dim/nprocs,1,howmany)
  integer :: ii,jj,atime,btime

  call ctcheck(dim)

  call myclock(atime)
  call gettwiddlesmall(twiddle1(:),dim/nprocs,nprocs)
  call myclock(btime); times(2)=times(2)+btime-atime; atime=btime

  call myzfft1d_slowindex_mpi(dim**2,in,tempout,dim,howmany)

  call myclock(btime); times(3)=times(3)+btime-atime; atime=btime

  tt1(:)=twiddle1(:)**(myrank-1)

  call myclock(btime); times(4)=times(4)+btime-atime; atime=btime

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,jj)
!$OMP DO SCHEDULE(STATIC)
  do ii=1,howmany
     do jj=1,dim/nprocs
        outtemp(:,:,jj,1,ii) = tempout(:,:,jj,1,ii) * tt1(jj)
     enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

  call myclock(btime); times(5)=times(5)+btime-atime; atime=btime

  call myzfft3d(outtemp,outtrans,dim,dim,dim/nprocs,howmany)

  call myclock(btime); times(6)=times(6)+btime-atime

end subroutine cooleytukey3d_outofplace_mpi



recursive subroutine cooleytukey3d_outofplaceinput_mpi(intranspose,out,dim,times,howmany)
  use ct_mpimod
  implicit none
  integer, intent(in) :: dim,howmany
  integer, intent(inout) :: times(8)
  complex*16, intent(in) :: intranspose(dim,dim,dim/nprocs,1,howmany)
  complex*16, intent(out) :: out(dim,dim,1,dim/nprocs,howmany)
  complex*16 :: twiddlefacs(dim/nprocs), tt(dim/nprocs),&  !! AUTOMATIC
       temptrans(dim,dim,dim/nprocs,1,howmany),outtrans(dim,dim,dim/nprocs,1,howmany)
  integer :: ii,jj,atime,btime

  call ctcheck(dim)

  call myclock(atime)

  call gettwiddlesmall(twiddlefacs(:),dim/nprocs,nprocs)

  call myclock(btime); times(2)=times(2)+btime-atime; atime=btime

  call myzfft3d(intranspose,temptrans,dim,dim,dim/nprocs,howmany)

  call myclock(btime); times(6)=times(6)+btime-atime; atime=btime

  tt(:)=twiddlefacs(:)**(myrank-1)

  call myclock(btime); times(4)=times(4)+btime-atime; atime=btime

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,jj)
!$OMP DO SCHEDULE(STATIC)
  do ii=1,howmany
     do jj=1,dim/nprocs
        outtrans(:,:,jj,1,ii) = temptrans(:,:,jj,1,ii) * tt(jj)
     enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

  call myclock(btime); times(5)=times(5)+btime-atime; atime=btime

  call myzfft1d_slowindex_mpi(dim**2,outtrans,out,dim,howmany)

  call myclock(btime); times(3)=times(3)+btime-atime; atime=btime

end subroutine cooleytukey3d_outofplaceinput_mpi


recursive subroutine myzfft1d_slowindex_mpi(bunchsize,in,out,dim,howmany)
  use ct_mpimod
  implicit none
  integer, intent(in) :: dim,howmany,bunchsize
  complex*16, intent(in) :: in(bunchsize,dim/nprocs,howmany)
  complex*16, intent(out) :: out(bunchsize,dim/nprocs,howmany)
  complex*16 :: fouriermatrix(nprocs,nprocs),twiddle(nprocs)
  integer :: ii
  call ctcheck(dim)
  call gettwiddlesmall(twiddle,nprocs,1)
  do ii=1,nprocs
     fouriermatrix(:,ii)=twiddle(:)**(ii-1)
  enddo
  select case (ctparopt)
  case(0)
     call simple_circ(in,out,fouriermatrix,bunchsize*dim/nprocs,howmany,myrank,nprocs)
  case(1)
     call simple_summa(in,out,fouriermatrix,bunchsize*dim/nprocs,howmany,myrank,nprocs)
  case default
     print *, "How did this happen?  was ctset() not called? internal ctparopt=",ctparopt; call mpistop()
  end select
end subroutine myzfft1d_slowindex_mpi


recursive subroutine simple_circ(in, out,mat,sizeper,howmany,myrank,nprocs)
  implicit none
  integer :: nnn,howmany,totsize,sizeper,myrank,nprocs
  complex*16 :: in(sizeper,howmany),     out(sizeper,howmany),&
      work2(sizeper,howmany),mat(nprocs,nprocs),work(sizeper,howmany)
  integer :: ibox,jbox,deltabox
  nnn=1
  totsize=sizeper*howmany
  out(:,:)=0
  do deltabox=0,nprocs-1
     ibox=mod(nprocs+myrank-1+deltabox,nprocs)+1
     jbox=mod(nprocs+myrank-1-deltabox,nprocs)+1
     work(:,:)=in(:,:)*mat(ibox,myrank)
     if (deltabox.ne.0) then
        call mympisendrecv(work(:,:),work2(:,:),ibox,jbox,deltabox,totsize)
        out(:,:)=out(:,:)+work2(:,:)
     else
        out(:,:)=out(:,:)+work(:,:)
     endif
  enddo
end subroutine simple_circ


recursive subroutine simple_summa(in, out,mat,sizeper,howmany,myrank,nprocs)
  implicit none
  integer :: nnn,sizeper,howmany,totsize,ibox,myrank,nprocs
  complex*16 :: in(sizeper,howmany),     out(sizeper,howmany),&
       work(sizeper,howmany),mat(nprocs,nprocs)
  nnn=1
  totsize=sizeper*howmany
  out(:,:)=0d0
  do ibox=1,nprocs
     if (myrank.eq.ibox) then
        work(:,:)=in(:,:)
     endif
     call mympibcast(work(:,:),ibox,totsize)
     out(:,:)=out(:,:)+work(:,:)*mat(myrank,ibox)
  enddo
end subroutine simple_summa


subroutine gettwiddlesmall(twiddlefacs,dim2,dim1)
  implicit none
  integer, intent(in) :: dim2,dim1
  complex*16, intent(out) :: twiddlefacs(dim2)
  complex*16 :: phi
  integer :: k2, itwiddle(dim2)
  real*8, parameter :: pi=3.14159265358979323846264338327950d0
  phi=exp((0d0,-2d0) * pi / (dim2*dim1))
  do k2=1,dim2
     itwiddle(k2)=(k2-1)
  enddo
  twiddlefacs(:)=phi**itwiddle(:)
end subroutine gettwiddlesmall
