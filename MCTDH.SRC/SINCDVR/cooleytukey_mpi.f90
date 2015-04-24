
module ct_mpimod
  implicit none
  integer :: myrank = -1
  integer :: nprocs = -1
  integer :: mpifileptr=6
end module ct_mpimod




subroutine ctset()
  use ct_mpimod
  implicit none
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


!! INVERSE OF cooleytukey_outofplace_mpi

subroutine cooleytukey_outofplace_inverse_mpi(intranspose,out,dim,howmany)
  use ct_mpimod
  implicit none
  integer, intent(in) :: dim,howmany
  complex*16, intent(in) :: intranspose(1,dim/nprocs,howmany)
  complex*16, intent(out) :: out(dim/nprocs,1,howmany)
  complex*16 ::  intransconjg(1,dim/nprocs,howmany), & !! AUTOMATIC
       outconjg(dim/nprocs,1,howmany)

  call ctcheck(dim)

  intransconjg(:,:,:)=CONJG(intranspose(:,:,:))
  call cooleytukey_outofplaceinput_mpi(intransconjg,outconjg,dim,howmany)
  out(:,:,:)=CONJG(outconjg(:,:,:))/dim

end subroutine cooleytukey_outofplace_inverse_mpi






!! fourier transform with OUT-OF-PLACE OUTPUT. 

subroutine cooleytukey_outofplace_mpi(in,outtrans,dim,howmany)
  use ct_mpimod
  implicit none
  integer :: ii
  integer, intent(in) :: dim,howmany
  complex*16, intent(in) :: in(dim/nprocs,1,howmany)
  complex*16, intent(out) :: outtrans(1,dim/nprocs,howmany)
  complex*16 :: twiddle1(dim/nprocs),tt1(dim/nprocs), &  !! AUTOMATIC
       tempout(dim/nprocs,1,howmany),      outtemp(dim/nprocs,1,howmany)

  call ctcheck(dim)

  call gettwiddlesmall(twiddle1(:),dim/nprocs,nprocs)

  call myzfft1d_slowindex_mpi(in,tempout,dim,howmany)

  tt1(:)=twiddle1(:)**(myrank-1)
  do ii=1,howmany
     outtemp(:,1,ii) = tempout(:,1,ii) * tt1(:)
  enddo

  call myzfft1d(outtemp,outtrans,dim/nprocs,howmany)

end subroutine cooleytukey_outofplace_mpi



subroutine cooleytukey_outofplaceinput_mpi(intranspose,out,dim,howmany)
  use ct_mpimod
  implicit none
  integer :: ii
  integer, intent(in) :: dim,howmany
  complex*16, intent(in) :: intranspose(dim/nprocs,1,howmany)
  complex*16, intent(out) :: out(1,dim/nprocs,howmany)
  complex*16 :: twiddlefacs(dim/nprocs), tt(dim/nprocs),&  !! AUTOMATIC
       temptrans(dim/nprocs,1,howmany),outtrans(dim/nprocs,1,howmany)

  call ctcheck(dim)

  call gettwiddlesmall(twiddlefacs(:),dim/nprocs,nprocs)

  call myzfft1d(intranspose,temptrans,dim/nprocs,howmany)

  tt(:)=twiddlefacs(:)**(myrank-1)
  do ii=1,howmany
     outtrans(:,1,ii) = temptrans(:,1,ii) * tt(:)
  enddo

  call myzfft1d_slowindex_mpi(outtrans,out,dim,howmany)

end subroutine cooleytukey_outofplaceinput_mpi








subroutine myzfft1d_slowindex_mpi(in,out,dim,howmany)
  use ct_mpimod
  implicit none
  integer, intent(in) :: dim,howmany
  complex*16, intent(in) :: in(dim/nprocs,howmany)
  complex*16, intent(out) :: out(dim/nprocs,howmany)
  complex*16 :: fouriermatrix(nprocs,nprocs),twiddle(nprocs)
  integer :: ii

  call ctcheck(dim)

  call gettwiddlesmall(twiddle,nprocs,1)

  do ii=1,nprocs
     fouriermatrix(:,ii)=twiddle(:)**(ii-1)
  enddo

  call simple_circ(in,out,fouriermatrix,dim/nprocs,howmany,myrank,nprocs)

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



!! NON-MPI SUBROUTINES

subroutine cooleytukey_putinplace(intranspose,out,dim1,dim2,howmany)
  implicit none
  integer, intent(in) :: dim1,dim2,howmany
  integer :: ii
  complex*16, intent(in) :: intranspose(dim2,dim1,howmany)
  complex*16, intent(out) :: out(dim1,dim2,howmany)
  do ii=1,howmany
     out(:,:,ii)=RESHAPE(TRANSPOSE(RESHAPE(intranspose(:,:,ii),(/dim1,dim2/))),(/dim1,dim2/))
  enddo
end subroutine cooleytukey_putinplace



  

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


