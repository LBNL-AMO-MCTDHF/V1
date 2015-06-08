
module pfileptrmod
  implicit none
  integer :: mpifileptr=-1
  integer,parameter :: nullfileptr=798   !! HARDWIRE TO MATCH PARAMETERS.F90
end module pfileptrmod

module pmpimod
  implicit none
  integer :: myrank=(-1),nprocs=(-1),sqnprocs=(-1),cbnprocs=(-1), &
       PROJ_COMM_WORLD=(-1), PROJ_GROUP_WORLD=(-1)
  integer :: boxrank(3)=1
  integer, allocatable :: BOX_COMM(:,:,:), BOX_GROUP(:,:,:), rankbybox(:,:,:)
!!  allocated BOX_COMM(procsplit(1),procsplit(2),orbparlevel:3)

  integer :: orbparlevel=3  !!  in namelist, but not in myparams

  integer :: procsplit(3)

end module pmpimod


subroutine bigdimsub(localdims,alldims)
  use pmpimod
  implicit none
  integer, intent(in) :: localdims(3)
  integer, intent(out) :: alldims(3)
  alldims(:)=localdims(:)*procsplit(:)
end subroutine bigdimsub


!! not crucially needed yet

subroutine splitgatherv_real(inlocal,outbig,bcastflag)
  use pmpimod
  use myparams
  implicit none
  logical,intent(in) :: bcastflag
  real*8, intent(in) :: inlocal(numpoints(1),numpoints(2),numpoints(3))
  real*8, allocatable :: ingather(:,:,:,:,:,:)
  real*8, intent(out) :: outbig(numpoints(1),procsplit(1),numpoints(2),procsplit(2),numpoints(3),procsplit(3))
  integer :: qqblocks(nprocs)
  integer :: ii,jj,kk

  if (myrank.eq.1) then
     allocate(ingather(numpoints(1),numpoints(2),numpoints(3),procsplit(1),procsplit(2),procsplit(3)))
  else
     allocate(ingather(1,1,1,1,1,1))
  endif

  qqblocks(:)=totpoints
  call mygatherv_real(inlocal,ingather,qqblocks,.false.)

  if (myrank.eq.1) then
     do ii=1,procsplit(3)
     do jj=1,procsplit(2)
     do kk=1,procsplit(1)
        outbig(:,kk,:,jj,:,ii) = ingather(:,:,:,kk,jj,ii)
     enddo
     enddo
     enddo
  endif

  deallocate(ingather)

  if (bcastflag) then
     call mympirealbcast(outbig,1,totpoints*nprocs)
  endif


end subroutine splitgatherv_real


subroutine splitgatherv_complex(inlocal,outbig,bcastflag)
  use pmpimod
  use myparams
  implicit none
  logical,intent(in) :: bcastflag
  complex*16, intent(in) :: inlocal(numpoints(1),numpoints(2),numpoints(3))
  complex*16, allocatable :: ingather(:,:,:,:,:,:)
  complex*16, intent(out) :: outbig(numpoints(1),procsplit(1),numpoints(2),procsplit(2),numpoints(3),procsplit(3))
  integer :: qqblocks(nprocs)
  integer :: ii,jj,kk

  if (myrank.eq.1) then
     allocate(ingather(numpoints(1),numpoints(2),numpoints(3),procsplit(1),procsplit(2),procsplit(3)))
  else
     allocate(ingather(1,1,1,1,1,1))
  endif

  qqblocks(:)=totpoints
  call mygatherv_complex(inlocal,ingather,qqblocks,.false.)

  if (myrank.eq.1) then
     do ii=1,procsplit(3)
     do jj=1,procsplit(2)
     do kk=1,procsplit(1)
        outbig(:,kk,:,jj,:,ii) = ingather(:,:,:,kk,jj,ii)
     enddo
     enddo
     enddo
  endif

  deallocate(ingather)

  if (bcastflag) then
     call mympicomplexbcast(outbig,1,totpoints*nprocs)
  endif


end subroutine splitgatherv_complex



subroutine splitscatterv_real(inbig,outlocal)
  use pmpimod
  use myparams
  implicit none
  real*8, intent(out) :: outlocal(numpoints(1),numpoints(2),numpoints(3))
  real*8, allocatable :: inscatter(:,:,:,:,:,:)
  real*8, intent(in) :: inbig(numpoints(1),procsplit(1),numpoints(2),procsplit(2),numpoints(3),procsplit(3))
  integer :: qqblocks(nprocs)
  integer :: ii,jj,kk

  if (myrank.eq.1) then
     allocate(inscatter(numpoints(1),numpoints(2),numpoints(3),procsplit(1),procsplit(2),procsplit(3)))

     do ii=1,procsplit(3)
     do jj=1,procsplit(2)
     do kk=1,procsplit(1)
        inscatter(:,:,:,kk,jj,ii)=inbig(:,kk,:,jj,:,ii)
     enddo
     enddo
     enddo

  else
     allocate(inscatter(1,1,1,1,1,1))
  endif

  qqblocks(:)=totpoints
  call myscatterv_real(inscatter,outlocal,qqblocks)


  deallocate(inscatter)

end subroutine splitscatterv_real




subroutine splitscatterv_complex(inbig,outlocal)
  use pmpimod
  use myparams
  implicit none
  complex*16, intent(out) :: outlocal(numpoints(1),numpoints(2),numpoints(3))
  complex*16, allocatable :: inscatter(:,:,:,:,:,:)
  complex*16, intent(in) :: inbig(numpoints(1),procsplit(1),numpoints(2),procsplit(2),numpoints(3),procsplit(3))
  integer :: qqblocks(nprocs)
  integer :: ii,jj,kk

  if (myrank.eq.1) then
     allocate(inscatter(numpoints(1),numpoints(2),numpoints(3),procsplit(1),procsplit(2),procsplit(3)))

     do ii=1,procsplit(3)
     do jj=1,procsplit(2)
     do kk=1,procsplit(1)
        inscatter(:,:,:,kk,jj,ii)=inbig(:,kk,:,jj,:,ii)
     enddo
     enddo
     enddo

  else
     allocate(inscatter(1,1,1,1,1,1))
  endif

  qqblocks(:)=totpoints
  call myscatterv_complex(inscatter,outlocal,qqblocks)


  deallocate(inscatter)

end subroutine splitscatterv_complex




