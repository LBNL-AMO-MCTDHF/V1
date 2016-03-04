
module pfileptrmod
  implicit none
  integer :: mpifileptr=-1
  integer,parameter :: nullfileptr=798   !! HARDWIRE TO MATCH PARAMETERS.F90
end module pfileptrmod

module pmpimod
  implicit none
  integer :: myrank=(-1),nprocs=(-1),&
       PROJ_COMM_WORLD=(-1), PROJ_GROUP_WORLD=(-1)

end module pmpimod


subroutine bigdimsub(localdims,alldims)
  use myparams
  use pmpimod
  implicit none
  integer, intent(in) :: localdims(3)
  integer, intent(out) :: alldims(3)
  if (localdims(1).ne.1.or.localdims(2).ne.1) then
     call waitawhile()
     print *, "bigdimsub error!",localdims,myrank
     call waitawhile()
     stop
  endif
  alldims(1:2)=1
  if (orbparflag) then
     alldims(3)=localdims(3)*nprocs
  else
     alldims(3)=localdims(3)
  endif
end subroutine bigdimsub






