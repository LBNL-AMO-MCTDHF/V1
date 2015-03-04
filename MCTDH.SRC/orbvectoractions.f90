
#include "Definitions.INC"


module orbvecmod
  implicit none
  integer, parameter :: maxsavefile=2000
  logical, save :: myopened(maxsavefile)=.false.
end module

subroutine close_orbvector(ifile)
  use orbvecmod
  implicit none
  integer :: ifile
  if (myopened(ifile)) then
     close(ifile)
     myopened(ifile)=.false.
  endif
end subroutine

subroutine read_orbvector(returnval,outvec, isize, ifile, filename, header)
  use orbvecmod
  use parameters
  use mpimod
  implicit none
  integer :: returnval, isize, ifile
  character :: filename*(*)
  character (len=headersize) :: header
  DATATYPE :: outvec(isize)

  if (ifile.gt.maxsavefile) then
     OFLWR "Error: programmer, please use file numbers less than ", maxsavefile;CFLST
  endif
  if (.not.myopened(ifile)) then
     open(ifile,file=filename, status="unknown", form="unformatted")
     myopened(ifile)=.true.
  endif
  read(ifile,iostat=returnval) header
  if (returnval/=0) then
     OFLWR "Done with vectors on file for header, iostat= ", returnval, " Filename= ", filename;CFL
     myopened(ifile)=.false.;     close(ifile);     return
  endif
  read(ifile,iostat=returnval) outvec
  if (returnval/=0) then
     OFLWR "Done with vectors on file, iostat= ", returnval, " Filename= ", filename, "  Size=", isize
     call closefile();     myopened(ifile)=.false.;     close(ifile);     return
  endif

end subroutine read_orbvector

subroutine save_orbvector(outvec, isize, ifile, filename, header)
  use orbvecmod
  use parameters
  use mpimod
  implicit none

  integer :: isize, ifile  !!, returnval
  character :: filename*(*)
  character (len=headersize) :: header
  DATATYPE :: outvec(isize)

  call noparorbsupport("save_orbvector")

  if (ifile.gt.maxsavefile) then
     OFLWR "Error: programmer, please use file numbers less than ", maxsavefile;CFLST
  endif
  if (.not.myopened(ifile)) then
     open(ifile,file=filename, form="unformatted");     myopened(ifile)=.true.
  endif
  write(ifile) header

!  if (returnval/=0) then
!     call openfile()
!     write(mpifileptr,*) "Error on write for header, iostat= ", returnval, "  Filename= ", filename
!     call closefile()
!     myopened(ifile)=.false.
!     close(ifile)
!     call mpistop()
!  endif

  write(ifile) outvec

end subroutine save_orbvector


subroutine get_spf_title(titleline, inspf, imvalue, thistime)
  implicit none
  character :: titleline*(*)
  integer :: inspf, imvalue
  real*8 :: thistime
  write(titleline,'(A15,I2,A6,I3,A6,F12.4,A3)') "set title 'Spf=",inspf,"  Mval=",imvalue,"  T= ", thistime, " ' "
end subroutine get_spf_title

subroutine get_natorb_title(titleline, inspf, imvalue, thistime, mydenval)
  implicit none
  character :: titleline*(*)
  integer :: inspf, imvalue
  real*8 :: thistime
  DATATYPE :: mydenval
#ifdef REALGO
  write(titleline,'(A15,I2,A6,I3,A6,F12.4,A6,E12.3,A3)') "set title 'Nat=",inspf,"  Mval=",imvalue,"  T= ", thistime, " Occ= ",mydenval," ' "
#else
  write(titleline,'(A15,I2,A6,I3,A6,F12.4,A6,2E12.3,A3)') "set title 'Nat=",inspf,"  Mval=",imvalue,"  T= ", thistime, " Occ= ",mydenval," ' "
#endif
end subroutine get_natorb_title

subroutine get_density_title(titleline, thistime)
  implicit none
  character :: titleline*(*)
  real*8 :: thistime
  write(titleline,'(A23,F12.4,A3)') "set title 'Density T= ", thistime, " ' "
end subroutine get_density_title

subroutine get_rnat_title(titleline, inspf, thistime, mydenval,iwhich)
  implicit none
  character :: titleline*(*)
  real*8 :: thistime
  integer :: inspf,iwhich
  DATATYPE :: mydenval
  character (len=5) :: labels(2) = (/ "Nat= ", "Proj=" /)
#ifdef REALGO
  write(titleline,'(A12,A5,I2,A6,F12.4,A10,1F13.8,A3)') "set title '",labels(iwhich),inspf,"  T= ", thistime, "  Eigval= ", mydenval," ' "
#else
  write(titleline,'(A12,A5,I2,A6,F12.4,A10,2F13.8,A3)') "set title '",labels(iwhich),inspf,"  T= ", thistime, "  Eigval= ", mydenval," ' "
#endif
end subroutine get_rnat_title

