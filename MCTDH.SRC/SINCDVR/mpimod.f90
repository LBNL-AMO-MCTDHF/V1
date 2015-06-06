
module fileptrmod
  implicit none
  integer :: mpifileptr=-1
  integer,parameter :: nullfileptr=798   !! HARDWIRE TO MATCH PARAMETERS.F90
end module fileptrmod

module mpimod
  implicit none
  integer :: myrank=(-1),nprocs=(-1)
end module mpimod


