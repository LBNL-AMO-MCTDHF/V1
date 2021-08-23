
module pfileptrmod
  implicit none
  integer :: mpifileptr=-1
  integer,parameter :: nullfileptr=798   !! HARDWIRE TO MATCH PARAMETERS.F90
end module pfileptrmod

module pmpimod
  implicit none
  integer :: myrank=(-1),nprocs=(-1),&
       PROJ_COMM_WORLD=(-1), PROJ_GROUP_WORLD=(-1)
  integer,parameter :: boxrank(3) = (/ 1,1,1 /), &  !! not used 1d (dummies)
       procsplit(3) = (/ -888,-888,-888 /), &       !! not used 1d (dummies)
       box_comm(1,1,3) = -798,   orbparlevel=1      !! not used 1d (dummies)

end module pmpimod


subroutine bigdimsub(localdims,alldims)
  use myparams
  use pmpimod
  use miscmod    !! IN PARENT DIRECTORY
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






