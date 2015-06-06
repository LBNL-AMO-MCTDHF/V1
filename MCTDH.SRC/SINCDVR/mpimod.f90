
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
!!  allocated BOX_COMM(nbox(1),nbox(2),orbparlevel:3)

  integer :: orbparlevel=3  !!  in namelist, but not in myparams

end module pmpimod


subroutine getProjCommGroup(outcommunicator,outgroup)
  use pmpimod
  implicit none
  integer, intent(out) :: outcommunicator,outgroup
  outcommunicator=PROJ_COMM_WORLD
  outgroup=PROJ_GROUP_WORLD
end subroutine getProjCommGroup

