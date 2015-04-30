

module ct_fileptrmod
  implicit none
  integer :: mpifileptr = 6
end module ct_fileptrmod

module ct_mpimod
  implicit none
  integer :: myrank = -1
  integer :: nprocs = -1
  integer :: CT_COMM_WORLD = -1
  integer :: CT_GROUP_WORLD = -1
end module ct_mpimod


subroutine ctset()
  use ct_fileptrmod
  use ct_mpimod
  implicit none
  call getmyranknprocs(myrank,nprocs)
  call getworldcommgroup(CT_COMM_WORLD,CT_GROUP_WORLD)
  if (myrank.eq.1) then
     mpifileptr=6
  else
     mpifileptr=989
     open(mpifileptr,file="/dev/null", status="unknown")
  endif
end subroutine ctset


subroutine twiddlemult_mpi(blocksize,in,out,dim1,numfactored,myfactor,localnumprocs,ctrank,howmany)
  use ct_fileptrmod
  use ct_mpimod   !! nprocs check
  implicit none
  integer, intent(in) :: blocksize,dim1,howmany,localnumprocs,numfactored,myfactor,ctrank
  complex*16, intent(in) :: in(blocksize,dim1,howmany)
  complex*16, intent(out) :: out(blocksize,dim1,howmany)
  complex*16 :: twiddle1(dim1,numfactored),tt1(dim1)
  integer :: ii,n1

  if (myfactor.lt.0.or.myfactor.gt.numfactored) then
     write (mpifileptr,*) "bad factor",myfactor,numfactored; call mpistop()
  endif

  call gettwiddlesmall(twiddle1(:,:),dim1*numfactored,localnumprocs)

  tt1(:)=twiddle1(:,myfactor)**(ctrank-1)   !!! ???? CRUX 
  do ii=1,howmany
     do n1=1,dim1
        out(:,n1,ii) = in(:,n1,ii) * tt1(n1)
     enddo
  enddo

end subroutine twiddlemult_mpi


subroutine myzfft1d_slowindex_mpi(in,out,localnumprocs,ctrank,proclist,totsize)
  implicit none
  integer, intent(in) :: totsize,localnumprocs,ctrank,proclist(localnumprocs)
  complex*16, intent(in) :: in(totsize)
  complex*16, intent(out) :: out(totsize)
  complex*16 :: fouriermatrix(localnumprocs,localnumprocs),twiddle(localnumprocs)
  integer :: ii

  call gettwiddlesmall(twiddle,localnumprocs,1)
  do ii=1,localnumprocs
     fouriermatrix(:,ii)=twiddle(:)**(ii-1)
  enddo
  call simple_summa(in,out,fouriermatrix,totsize,ctrank,localnumprocs,proclist)

end subroutine myzfft1d_slowindex_mpi




subroutine simple_circ(in, out,mat,howmany,ctrank,localnumprocs,proclist)
  use ct_fileptrmod
  use ct_mpimod
  implicit none
  integer, intent(in) :: howmany,ctrank,localnumprocs,proclist(localnumprocs)
  complex*16, intent(in) :: in(howmany), mat(localnumprocs,localnumprocs)
  complex*16, intent(out) :: out(howmany)
  complex*16 :: work2(howmany),work(howmany)
  integer :: ibox,jbox,deltabox,nnn,CT_GROUP_LOCAL,CT_COMM_LOCAL,ierr,procshift(localnumprocs)
  
  ierr=798
  if (localnumprocs.gt.nprocs.or.localnumprocs.le.1.or.ctrank.lt.1.or.ctrank.gt.localnumprocs) then
     write(mpifileptr,*) "local error", ctrank,localnumprocs,nprocs; call mpistop()
  endif

  procshift(:)=proclist(:)-1
#ifndef MPIFLAG
  if (ctrank.ne.1) then
     write(mpifileptr,*) "Error non-mpi rank ne 1", ctrank; call mpistop()
  endif
  CT_GROUP_LOCAL=(-42)
  CT_COMM_LOCAL=798
#else
  call mpi_group_incl(CT_GROUP_WORLD,localnumprocs,procshift,CT_GROUP_LOCAL,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "Error group incl simple_circ",ierr; call mpistop()
  endif
  call mpi_comm_create(CT_COMM_WORLD, CT_GROUP_LOCAL, CT_COMM_LOCAL,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "Error comm create simple_circ",ierr; call mpistop()
  endif
#endif

  nnn=1

  out(:)=0

  do deltabox=0,localnumprocs-1

     ibox=mod(localnumprocs+ctrank-1+deltabox,localnumprocs)+1
     jbox=mod(localnumprocs+ctrank-1-deltabox,localnumprocs)+1

     work(:)=in(:)*mat(ibox,ctrank)

     if (deltabox.ne.0) then
        call mympisendrecv(work(:),work2(:),ibox,jbox,deltabox,howmany,CT_COMM_LOCAL)
        out(:)=out(:)+work2(:)
     else
        out(:)=out(:)+work(:)
     endif
  enddo

#ifdef MPIFLAG
  call mpi_comm_free(CT_COMM_LOCAL,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "Error comm destroy simple_circ",ierr; call mpistop()
  endif
  call mpi_group_free(CT_GROUP_LOCAL,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "Error group destroy simple_circ",ierr; call mpistop()
  endif
#endif

end subroutine simple_circ



subroutine simple_summa(in, out,mat,howmany,ctrank,localnumprocs,proclist)
  use ct_fileptrmod
  use ct_mpimod
  implicit none
  integer, intent(in) :: howmany,ctrank,localnumprocs,proclist(localnumprocs)
  complex*16, intent(in) :: in(howmany), mat(localnumprocs,localnumprocs)
  complex*16, intent(out) :: out(howmany)
  complex*16 :: work(howmany)
  integer :: ibox,nnn,CT_GROUP_LOCAL,CT_COMM_LOCAL,ierr,procshift(localnumprocs)

  ierr=(-798)
  if (localnumprocs.gt.nprocs.or.localnumprocs.le.1.or.ctrank.lt.1.or.ctrank.gt.localnumprocs) then
     write(mpifileptr,*) "local error", ctrank,localnumprocs,nprocs; call mpistop()
  endif

  procshift(:)=proclist(:)-1

#ifndef MPIFLAG
  if (ctrank.ne.1) then
     write(mpifileptr,*) "Error non-mpi rank ne 1", ctrank; call mpistop()
  endif
  CT_GROUP_LOCAL=(-42)
  CT_COMM_LOCAL=798
#else
  call mpi_group_incl(CT_GROUP_WORLD,localnumprocs,procshift,CT_GROUP_LOCAL,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "Error group incl simple_summa",ierr; call mpistop()
  endif
  call mpi_comm_create(CT_COMM_WORLD, CT_GROUP_LOCAL, CT_COMM_LOCAL,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "Error comm create simple_summa",ierr; call mpistop()
  endif
#endif

  nnn=1

  out(:)=0d0
  do ibox=1,localnumprocs
     if (ctrank.eq.ibox) then
        work(:)=in(:)
     endif
     call mympibcast(work(:),ibox,howmany,CT_COMM_LOCAL)
     out(:)=out(:)+work(:)*mat(ctrank,ibox)
  enddo

#ifdef MPIFLAG
  call mpi_comm_free(CT_COMM_LOCAL,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "Error comm destroy simple_summa",ierr; call mpistop()
  endif
  call mpi_group_free(CT_GROUP_LOCAL,ierr)
  if (ierr.ne.0) then
     write(mpifileptr,*) "Error group destroy simple_summa",ierr; call mpistop()
  endif
#endif

end subroutine simple_summa




subroutine myzfft1d_slowindex_local(in,out,dim1,dim2,howmany)
  implicit none
  integer, intent(in) :: dim1,dim2,howmany
  complex*16, intent(in) :: in(dim1,dim2,howmany)
  complex*16, intent(out) :: out(dim1,dim2,howmany)
  complex*16 :: intrans(dim2,dim1,howmany),outtrans(dim2,dim1,howmany)
  integer :: ii

  do ii=1,howmany
     intrans(:,:,ii)=TRANSPOSE(in(:,:,ii))
  enddo
  call myzfft1d(intrans,outtrans,dim2,dim1*howmany)

  do ii=1,howmany
     out(:,:,ii)=TRANSPOSE(outtrans(:,:,ii))
  enddo

end subroutine myzfft1d_slowindex_local
  


subroutine getprimefactor(dim,myfactor)
  implicit none
  integer, intent(in) :: dim
  integer, intent(out) :: myfactor
  integer :: iprime
  integer, parameter :: numprimes=31
  integer, parameter :: primelist(numprimes)=&
       (/  2,  3,  5,  7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,&
          73, 79, 83, 89, 97,101,103,107,109,113,127 /)  ! no need to go remotely this high

  myfactor=dim

  do iprime=1,numprimes
     if (mod(dim,primelist(iprime)).eq.0) then
        myfactor=primelist(iprime)
        return
     endif
  enddo

end subroutine getprimefactor

subroutine getallprimefactors(dim,numfactors,allfactors)
  implicit none
  integer, intent(in) :: dim
  integer, intent(out) :: allfactors(7),numfactors
  integer :: thisdim,flag
  allfactors(:)=1
  numfactors=1
  thisdim=dim
  flag=0

798 if (numfactors.eq.7) then
     allfactors(7)=thisdim
     return
  else
     call getprimefactor(thisdim,allfactors(numfactors))
     if (allfactors(numfactors).eq.thisdim) then
        return
     endif
     thisdim=thisdim/allfactors(numfactors)
     numfactors=numfactors+1
  endif
go to 798
end subroutine getallprimefactors
  


!! haven't done proper version with recursion
!
!subroutine cooleytukey_replace(blocksize,intranspose,out,dim1,dim2,howmany)
!  implicit none
!  integer, intent(in) :: dim1,dim2,howmany,blocksize
!  complex*16, intent(in) :: intranspose(blocksize,dim1,dim2,howmany)
!  complex*16, intent(out) :: out(blocksize,dim1,dim2,howmany)
!  integer :: ii,jj
!
!  do ii=1,howmany
!     do jj=1,blocksize
!        out(jj,:,:,ii)=RESHAPE(TRANSPOSE(intranspose(jj,:,:,ii)),(/dim1,dim2/))
!     enddo
!  enddo
!
!end subroutine cooleytukey_replace



  

subroutine gettwiddlesmall(twiddlefacs,dim1,dim2)
  implicit none
  integer, intent(in) :: dim1,dim2
  complex*16, intent(out) :: twiddlefacs(dim1)
  complex*16 :: phi
  integer :: k1, itwiddle(dim1)
  real*8, parameter :: pi=3.14159265358979323846264338327950d0
  phi=exp((0d0,-2d0) * pi / (dim1*dim2))
  do k1=1,dim1
     itwiddle(k1)=(k1-1)
  enddo
  twiddlefacs(:)=phi**itwiddle(:)
end subroutine gettwiddlesmall



