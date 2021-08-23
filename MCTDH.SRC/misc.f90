
#include "Definitions.INC"


module miscmod
contains
  
subroutine waitawhile()
  implicit none
  integer :: i,j
  character (len=10) :: mytext
  j=1
  do i=1,10000000
     j=mod(j*i,1777)
  enddo
  write(mytext,'(I10)') j
  call system("echo "//mytext//" >> /dev/null")
end subroutine waitawhile


subroutine checkiostat(iniostat,intext)
  use fileptrmod
  implicit none
  integer,intent(in) :: iniostat
  character*(*),intent(in) :: intext
  if (iniostat /=0 ) then
     print *, "MCTDHF ABORT: I/O error ", iniostat,intext
     OFLWR "MCTDHF ABORT: I/O error ", iniostat,intext; CFL
     call waitawhile()
     stop
     stop   !!   STOP.   !!
     stop
  endif
end subroutine checkiostat

!! v1.27 getlen now reports length of string not length of string plus one

function getlen(buffer)
  implicit none
  character buffer*(*)
  integer :: j, getlen, mylen
  mylen=LEN(buffer)
  j=1
  do while (j.le.mylen)
     if (buffer(j:j) .eq. " ") then
        getlen=j-1
        return
     else
        j=j+1
     endif
  enddo
  getlen=mylen
end function getlen

function getlen2(buffer)
  implicit none
  character buffer*(*)
  integer :: j, getlen2, nn
  nn=LEN(buffer)-4
  j=1
  do while ((j.lt.nn).and..not.(buffer(j:j+3) .eq. "    "))
     j=j+1
  enddo
  getlen2=j-1
end function getlen2

end module miscmod

