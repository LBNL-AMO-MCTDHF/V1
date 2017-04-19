
!! FOR PLOTTING ACTIONS (2-12)

#include "Definitions.INC"


module orbvecmod
  implicit none
  integer, parameter :: maxsavefile=2000
  logical, save :: myopened(maxsavefile)=.false.
end module

subroutine close_orbvector(ifile)
  use orbvecmod
  use mpimod
  implicit none
  integer :: ifile
  if (myrank.eq.1) then
     if (myopened(ifile)) then
        close(ifile)
        myopened(ifile)=.false.
     endif
  endif
end subroutine

subroutine read_orbvector(returnval,outvec, isize, ifile, filename, header)
  use orbvecmod
  use parameters
  use mpimod
  use mpisubmod
  implicit none
  integer :: returnval, isize, ifile
  character :: filename*(*)
  character (len=headersize) :: header
  DATATYPE,intent(out) :: outvec(isize)

  if (ifile.gt.maxsavefile) then
     OFLWR "Error: programmer, please use file numbers less than ", maxsavefile;CFLST
  endif
  if (myrank.eq.1) then
     if (.not.myopened(ifile)) then
        open(ifile,file=filename, status="unknown", form="unformatted")
        myopened(ifile)=.true.
     endif
     read(ifile,iostat=returnval) header
  endif
  call mympiibcastone(returnval,1)

  if (returnval/=0) then
     OFLWR "Done with vectors on file for header, iostat= ", returnval, " Filename= ", filename;CFL
     if (myrank.eq.1) then
        myopened(ifile)=.false.;     
        close(ifile);
     endif
     return
  endif
  if (myrank.eq.1) then
     read(ifile,iostat=returnval) outvec
  endif
  call mympibcast(outvec,1,isize)
  call mympiibcastone(returnval,1)

  if (returnval/=0) then
     OFLWR "Done with vectors on file, iostat= ", returnval, " Filename= ", filename, "  Size=", isize
     call closefile(); 
     if (myrank.eq.1) then
        myopened(ifile)=.false.;
        close(ifile);
     endif
  endif

end subroutine read_orbvector

subroutine save_orbvector(outvec, isize, ifile, filename, header)
  use orbvecmod
  use parameters
  use mpimod
  implicit none

  integer :: isize, ifile
  character :: filename*(*)
  character (len=headersize) :: header
  DATATYPE :: outvec(isize)

  if (ifile.gt.maxsavefile) then
     OFLWR "Error: programmer, please use file numbers less than ", maxsavefile;CFLST
  endif
  if (myrank.eq.1) then
     if (.not.myopened(ifile)) then
        open(ifile,file=filename, form="unformatted")
        myopened(ifile)=.true.
     endif
     write(ifile) header
     write(ifile) outvec
  endif
!! BUGFIX 041817  call mpibarrier()

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

