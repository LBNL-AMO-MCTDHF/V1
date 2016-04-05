
!! LOTS OF ACTIONS FOR READING ORBITALS FOR PLOTTING ROUTINES

#include "Definitions.INC"

!! iwhich=4, denproj, not debugged

subroutine read_orb_initial(iwhich)
  use parameters
  use mpimod
  implicit none
  integer ::  ispf,imvalue=0, iflag, jflag, iwhich, xxflag=0, jspf, iallflag
  character (len=10) :: labels(4) = (/ "Natorb    ", "Spfs      ", "Density   ", "Denproj   " /)
  character (len=10) :: inchar

  if (stdoutflag/=1) then
     OFLWR "Need to run the program with standard output to see ", labels(iwhich),"."; CFLST
  endif

  select case (iwhich)
  case (1)
     print *, " ****************************** " 
     print *, "    VIEWING NATURAL ORBITALS    "
     print *, " ****************************** "
  case (2)
     print *, " ****************************** " 
     print *, "           VIEWING SPFS"
     print *, " ****************************** "
  case (3)
     print *, " ****************************** " 
     print *, "         VIEWING DENSITY"
     print *, " ****************************** "
  case (4)
     print *, " ****************************** " 
     print *, "         VIEWING DENPROJ"
     print *, " ****************************** "
  case default
     print *, " Error, iwhich not supported in read_orb_initial, ", iwhich
     call mpistop()
  end select
  print *

  xxflag=xxflag+1
  if (xxflag==1) then
     print *, "USE POVRAY? y for yes."
     read(*,*) inchar
     if (inchar(1:1).eq.'y') then
        povplotflag=1
        print *, "OK, using povray."
     else
        povplotflag=0
        print *, "OK, just using gnuplot."
     endif
  endif
  print *

  iflag=1
  do while (iflag==1) 
     jflag=0
     do while ( jflag == 0 ) 
        call whatsub(jflag)
     enddo
     if (iwhich.ne.3) then
        print *, "Enter ", labels(iwhich), " number.  Zero to change options.  Negative to stop."
        read(*,*) ispf
     else
        print *, "Enter positive number to continue to plot.  Zero to change options.  Negative to stop."
        read(*,*) jspf
        ispf=1
     endif
!     iprop=1
!     if (numprop.gt.2) then
!        print *, "Which calculation?  There are ", numprop," concurrent."
!        read(*,*) iprop
!     endif
     if (ispf<0) then
        print *, "OK! Done."
        iflag=0
     else if (ispf>0) then
        iallflag=0
        if (povplotflag/=1) then
           if (iwhich.ne.3) then
              imvalue=spfmvals(ispf);              print *, "this one is ", imvalue
              print *, "Enter m projection";              read(*,*) imvalue
           else
              print *, "Enter m projection";              read(*,*) imvalue
           endif
        endif
        if (ispf>nspf) then
           print *, "Will do all!";           iallflag=1
        else
           if (ispf>nspf) then
              print *, "Will do all!";              iallflag=1
           endif
        endif
        if (iallflag==1) then
           do jspf=1,nspf
              if (spfrestrictflag==1) then
                 imvalue=spfmvals(jspf)
              endif
              print *, "imvalue is ", imvalue
              call read_orb(jspf,imvalue, 1, iwhich)
!! somewhat kloogey... here hardwired for m-value third index
!!              call read_orb(jspf,imvalue+(spfdims(3)+1)/2, 1, iwhich)
           enddo
        else
           call read_orb(ispf,imvalue, 1, iwhich)           
!!           call read_orb(ispf,imvalue+(spfdims(3)+1)/2, 1, iwhich)           
        endif
     endif
  enddo
  call mpistop()

end subroutine read_orb_initial


subroutine read_Rorb_initial(iwhich)
  use parameters
  use mpimod
  implicit none
  integer ::  ispf, iflag,jflag, iwhich
  character (len=10) :: labels(2) = (/ "R-Natorbs ", "Natprojs  " /)

  if (stdoutflag/=1) then
     OFLWR "Need to run the program with standard output to see natorbs."; CFLST
  endif
  print *, " *********************************************** " 
  select case (iwhich)
  case (1)
     print *, "    VIEWING NATURAL ORBITALS in R    "
  case (2)
     print *, "    VIEWING PROJECTIONS OF NATURAL CONFIGURATIONS."
  case default
     print *, "Error, iwhich not supported in read_Rorb, ", iwhich
     call mpistop()
  end select
  print *, " *********************************************** ";  print *;  print *

  iflag=1
  do while (iflag==1) 
     jflag=0
     do while ( jflag == 0 )
        call whatsub(jflag)
     enddo
     print *, "Enter ", labels(iwhich)," number.  Zero to change options.  Negative to stop."
     read(*,*) ispf
!     iprop=1
!     if (numprop.gt.2) then
!        print *, "Which calculation?  There are ", numprop," concurrent."
!        read(*,*) iprop
!     endif
     if (ispf.gt.numr) then
        print *, "Um, that is greater than numr. goodbye.";        stop
     endif
     if (ispf<0) then
        print *, "OK! Done.";        iflag=0
     else if (ispf>0) then
!!        call read_Rorb(ispf, iprop,iwhich)
        call read_Rorb(ispf, 1,iwhich)
     endif
  enddo
  call mpistop()
end subroutine read_Rorb_initial


module zerorecmod
contains
function zerorecord(inint)
  implicit none
  character (len=4) :: zerorecord, work
  integer :: inint
  work="0000"
  if (inint.lt.10) then
     write(work(4:),'(I1)') inint
  else if (inint.lt.100) then
     write(work(3:),'(I2)') inint
  else if (inint.lt.1000) then
     write(work(2:),'(I3)') inint
  else if (inint.lt.10000) then
     write(work(1:),'(I4)') inint
  else
     print *, "CHAR ERROR."
     stop
  endif
  zerorecord=work
end function zerorecord
end module zerorecmod


!! iwhich: 1=spf 2=natorb 3=density  4=denproj   denproj not debugged

subroutine read_orb(inspf,imvalue, iprop, iwhich)   
  use parameters
  use mpimod
  use zerorecmod
  implicit none

  integer :: readprop,  returnval, ixi, ieta, inspf, imvalue, flag, xflag, irecord, first,  iprop, &
       getlen, iwhich, ispf, iii, itable
  integer, parameter :: ifilenums(4)=(/  natorbfile, spfplotfile, denfile, denprojfile /)
  character (len=SLN) :: ifilenames(4)
  character (len=8), parameter :: ignufile(4)=(/ "Natplot_", "Spfplot_", "Denplot_", "Denproj_" /)
  character (len=10) :: povdirs(4) = (/ "Natorb    ", "Spfs      ", "Density   ", "Denproj   " /)
  character(len=10) :: mypovdir
  real*8 :: thistime, zval, xval,  costheta, rho
  DATATYPE :: sum, cylindricalvalue, mydenval
  DATATYPE :: ttempspf(spfdims(1),spfdims(2),spfdims(3))
  character (len=8) :: tableexts(4)=(/'.table0"','.table1"','.table2"','.table3"' /)
  character (len=100) :: filename,sysline, titleline
  character(len=20) :: giffilename
  character(len=40) :: epsfilename,epstable(4)
  character(len=100) :: filex="                                                                                                    "
  character (len=headersize) :: header

  ifilenames(1)= natplotbin;  ifilenames(2)=  spfplotbin 
  ifilenames(3)= denplotbin;  ifilenames(4)= denprojplotbin 
  irecord=1

  OFLWR "READ IMVALUE ", imvalue;CFL

  if ((iwhich.lt.1).or.(iwhich.gt.4)) then
     OFLWR "Improper iwhich in read_orb : ", iwhich; CFLST
  endif

  first=1;  flag=0
  do while (flag.lt.2)
     flag=0
     do while ( flag.eq.0 )
        call read_orbvector(returnval,ttempspf(:,:,:),spfsize,ifilenums(iwhich),ifilenames(iwhich),header)
        select case (iwhich)
        case (1)
           call read_nat_header(header,thistime,readprop,ispf,mydenval)
        case (2)
           call read_spf_header(header,thistime,readprop,ispf)
        case (3)
           call read_den_header(header,thistime,readprop)
           ispf=1
        case (4)
           call read_nat_header(header,thistime,readprop,ispf,mydenval)
        end select
        if (iwhich.eq.1) then  !! natorb
           call fixphase0(ttempspf(:,:,:),spfsize)
        endif
        if (returnval/=0) then
           call openfile();   write(mpifileptr,*) "Done with vectors on file.";  call closefile()
           if (povplotflag.eq.0) then
              close(871)   
           endif
           flag=3
        else
           xflag=0
           select case (iwhich)
           case (1,2,4)
              if ((inspf==ispf).and.(readprop.eq.iprop)) then
                 xflag=1
              endif
           case (3)
              if (readprop.eq.iprop) then
                 xflag=1
              endif
           end select
           if (xflag==1) then
              if (mod(irecord-1,plotskip).eq.0) then 
                 flag=1
                 print *,  "Good read at record ",irecord
              endif
              if (irecord.ge.plotnum*plotskip) then
                 print *,  "Stopping at ", plotnum, " plots."
                 flag=2
              else
                 irecord=irecord+1
              endif
           endif
        endif
     enddo
     if ((flag==1).or.(flag==2)) then
        if (povplotflag==1) then
           mypovdir=povdirs(iwhich)
           call read_povother( thistime, ttempspf(:,:,:), iprop, mypovdir(1:getlen(mypovdir)), iwhich, ispf, irecord-1)
        else
           if (first==1) then
              first=0;              filename=filex
              giffilename="                    "
              if (iwhich.eq.3) then
                 write(filename,'(A8,I1,A4)') ignufile(iwhich),readprop,".gnu"
                 write(giffilename,'(A1,A8,I1,A5)') '"',ignufile(iwhich),readprop,'.gif"'
              else
                 if (inspf.lt.10) then
                    write(filename,'(A8,I1,A1,I1,A4)') ignufile(iwhich),readprop,"_",inspf,".gnu"
if (plotterm.eq.3) then
                    write(giffilename,'(A1,A8,I1,A1,I1,A5)') '"',ignufile(iwhich),readprop,"_",inspf,'.eps"'
else
                    write(giffilename,'(A1,A8,I1,A1,I1,A5)') '"',ignufile(iwhich),readprop,"_",inspf,'.gif"'
endif
                 else
                    write(filename,'(A8,I1,A1,I2,A4)') ignufile(iwhich),readprop,"_",inspf,".gnu"
                    write(giffilename,'(A1,A8,I1,A1,I2,A5)') '"',ignufile(iwhich),readprop,"_",inspf,'.gif"'
                 endif
              endif
              open(871,file=filename(1:getlen(filename)),status="unknown")
              call writegnuoptions(871, giffilename)
           endif

           if (plotterm.eq.3.and.pm3d.eq.0) then
#ifdef REALGO
              do itable=1,2
#else
              do itable=1,4
#endif
                 epstable(itable)="                                        "
                 if (iwhich.eq.3) then
                    write(epstable(itable),'(A1,A8,I1,A8)') '"',ignufile(iwhich),readprop,tableexts(itable)
                 else
                    if (inspf.lt.10) then
                       write(epstable(itable),'(A1,A8,I1,A1,I1,A1,A4,A8)') '"',ignufile(iwhich),readprop,"_",inspf,'_', zerorecord(irecord),tableexts(itable)
                    else
                       write(epstable(itable),'(A1,A8,I1,A1,I2,A1,A4,A8)') '"',ignufile(iwhich),readprop,"_",inspf,'_', zerorecord(irecord),tableexts(itable)
                    endif
                 endif
                 write(871,*) "unset surf"
                 write(871,*) "set out"
                 write(871,'(A11,A40)') " set table  ", epstable(itable)(1:getlen(epstable(itable)))
                 if (mod(itable,2).eq.0) then
                    write(871,*) " set zrange [0:*]"
                 else
                    write(871,*) " set zrange [*:0]"
                 endif
#ifdef REALGO
                 write(871,*) " splot '-' using 2:1:3"
#else
                 write(871,*) " splot '-' using 1:2:3"
#endif
                 do ixi=(-1)*plotres,plotres-1
                    do ieta=(-1)*plotres,plotres-1
                       xval=(ixi+0.5d0)*plotxyrange/plotres
                       zval=(ieta+0.5d0)*plotxyrange/plotres
                       rho=sqrt(xval**2+zval**2)
                       costheta=cos(atan2(zval,xval))
                       
                       sum=cylindricalvalue(rho,costheta,1.d0,imvalue, ttempspf(:,:,imvalue+(spfdims(3)+1)/2))
                       if (zval.lt.0.d0) then
                          sum=sum*(-1)**abs(imvalue)
                       endif
#ifndef REALGO
                       if (itable.lt.3) then
#endif
                          write(871,'(2F8.3, F14.8)') zval,xval,real(sum)
#ifndef REALGO
                       else
                          write(871,'(2F8.3, F14.8)') zval,xval,imag(sum)
                       endif
#endif                       
                    enddo
                    write(871,*)
                 enddo
                 write(871,*) "e";                 write(871,*) ;                 write(871,*) "unset table"
                 write(871,*) 
              enddo !! itable
           endif

           select case(iwhich)
           case (1)
              call get_natorb_title(titleline, inspf, imvalue, thistime, mydenval)
           case (2)
              call get_spf_title(titleline, inspf, imvalue, thistime)
           case (3)
              call get_density_title(titleline, thistime)
           case (4)
              call get_natorb_title(titleline, inspf, imvalue, thistime, mydenval)
           end select

           write(871,*) titleline

           if (plotterm.eq.3) then
              epsfilename="                                        "
              if (iwhich.eq.3) then
                 write(epsfilename,'(A1,A8,I1,A5)') '"',ignufile(iwhich),readprop,'.eps"'
              else
                 if (inspf.lt.10) then
                    write(epsfilename,'(A1,A8,I1,A1,I1,A1,A4,A5)') '"',ignufile(iwhich),readprop,"_",inspf,'_', zerorecord(irecord),'.eps"'
                 else
                    write(epsfilename,'(A1,A8,I1,A1,I2,A1,A4,A5)') '"',ignufile(iwhich),readprop,"_",inspf,'_', zerorecord(irecord),'.eps"'
                 endif
              endif
              
              write(871,'(A11,A40)') " set out  ", epsfilename(1:getlen(epsfilename))

              if (pm3d.eq.1) then
                 write(871,*) " unset surf "
                 write(871,*) " splot '-' using 2:1:(sqrt($3**2+$4**2)):(atan($3/$4)) lt -1 lw 4, 'points.dat' using 1:2:(0) with points pt 7 ps 2 lt -1"
              else
                 write(871,*) " set surf"
                 write(871,*) " set zrange [*:*]"
                 write(871,*) " set style line 1 lt -1 lw 3"
                 write(871,*) " set term post color enhanced solid  24 eps"

#ifdef REALGO
                 write(871,*) " plot ", epstable(1)(1:getlen(epstable(1)))," lw 4, \\"
                 write(871,*)  epstable(2)(1:getlen(epstable(2))), " lw 4 , 'points.dat' with points pt 7 ps 2 lt -1"
#else
                 write(871,*) " plot ", epstable(1)(1:getlen(epstable(1)))," lw 4, ", epstable(2)(1:getlen(epstable(2))), " lw 4, ", &
                      epstable(3)(1:getlen(epstable(3)))," using ($1+2*",plotxyrange,"):2:3 lw 4, ", epstable(4)(1:getlen(epstable(4))), &
                      " using ($1+2*",plotxyrange,"):2 lw 4 , 'points.dat' with points pt 7 ps 2 lt -1, 'points.dat' using ($1+2*",plotxyrange,"):2 with points pt 7 ps 2 lt -1"
#endif
              endif  !! pm3d
           else !! plotterm eq 3
              if (pm3d==1) then
                 write(871,*) " splot '-' using 1:($2+2*",plotxyrange,"):3:3,'-' using 1:2:3:3"
              else
                 write(871,*) " splot '-' using 1:($2+2*",plotxyrange,"):3,'-' using 1:2:3"
              endif
           endif
           write(871,*)
           
           if (plotterm.ne.3.or.pm3d.ne.0) then  !! we want to plot out some data otherwise term 3, 3d 0 and we have tables and we are done
              do ixi=(-1)*plotres,plotres-1
                 do ieta=(-1)*plotres,plotres-1
                    xval=(ixi+0.5d0)*plotxyrange/plotres
                    zval=(ieta+0.5d0)*plotxyrange/plotres
                    rho=sqrt(xval**2+zval**2)
                    costheta=cos(atan2(zval,xval))
                    
                    sum=cylindricalvalue(rho,costheta,1.d0,imvalue, ttempspf(:,:,imvalue+(spfdims(3)+1)/2))
                    
                    if (zval.lt.0.d0) then
                       sum=sum*(-1)**abs(imvalue)
                    endif
                    if (plotterm.ne.3) then
                       write(871,'(2F8.3, F14.8)') zval,xval,imag(sum+(0.d0,0.d0))
                    else
                       write(871,'(2F8.3, 2F14.8)') zval,xval,sum
                    endif
                 enddo
                 write(871,*)
              enddo
              write(871,*) "e"
           
              if ((plotterm.ne.3)) then  ! plot out the imaginary set as well unless term 3
                 
                 do ixi=(-1)*plotres,plotres-1
                    do ieta=(-1)*plotres,plotres-1
                       xval=(ixi+0.5d0)*plotxyrange/plotres
                       zval=(ieta+0.5d0)*plotxyrange/plotres
                       rho=sqrt(xval**2+zval**2)
                       costheta=cos(atan2(zval,xval))
                       
                       sum=cylindricalvalue(rho,costheta,1.d0,imvalue, ttempspf(:,:,imvalue+(spfdims(3)+1)/2))
                       if (zval.lt.0.d0) then
                          sum=sum*(-1)**abs(imvalue)
                       endif
                       write(871,'(2F8.3, F18.12)') zval,xval,real(sum+(0.d0,0.d0))
                    enddo
                    write(871,*)
                 enddo
                 write(871,*) "e"
              endif  !! plot out second set
           endif  !! plotting out the data
           write(871,*) "pause ", plotpause
        endif  !! povplotflag
     endif !! flag==1 or 2
     open(222,file="stop",status="old", iostat=iii)
     if (iii==0) then
        close(222)
        print *, "Stopping due to stopfile!"
        stop
     endif
  enddo

  call close_orbvector(ifilenums(iwhich))
  if (povplotflag.eq.0) then
     close(871)

     if (myrank.eq.1) then
        print *, "Plotted to file.  CAlling GNU!!!"
     endif
     sysline(1:100)=filex
     write(sysline,'(A18,A20,A1)') "gnuplot -persist ",filename,"&";     call system(sysline)
  endif

end subroutine read_orb


subroutine read_Rorb(inspf, iprop,iwhich)
  use parameters
  use mpimod
  implicit none

  integer :: readprop,   returnval, ir, inspf, flag, irecord, first, iprop, getlen, iwhich, ispf
  real*8 :: thistime
  DATATYPE :: mydenval
  character (len=100) :: filename,sysline , titleline
  character(len=20) :: giffilename
  character(len=20) :: epsfilename
  character(len=100) :: filex="                                                                                                    "
  character (len=headersize) :: header
  character (len=SLN) :: ifilenames(2)
  integer, parameter :: ifilenums(2) = (/ rnatorbfile, natprojfile /)
  character (len=9), parameter :: ignufile(2)=(/ "RNatplot_", "Projplot_" /)
  DATATYPE :: inrdenvect(numr)  

  ifilenames(1)= rnatplotbin
  ifilenames(2)= natprojplotbin 

  irecord=1;  first=1;  flag=0
  do while (flag.lt.2)
     flag=0
     do while ( flag.eq.0 )
        call read_orbvector(returnval,inrdenvect,numr,ifilenums(iwhich),ifilenames(iwhich),header)
        call read_nat_header(header,thistime,readprop,ispf, mydenval)
        if (iwhich.eq.1) then  !! natorb
           call fixphase0(inrdenvect,numr)
        endif
        if (returnval/=0) then
           print *,  "Done with vectors on file for rden= ", inspf
           close(871);           flag=3
        else
           if ((inspf==ispf).and.(readprop.eq.iprop)) then
              if (mod(irecord-1,plotskip).eq.0) then;                 flag=1
                 print *, "Good read at record ",irecord
              endif
              if (irecord.ge.plotnum*plotskip) then
                 print *,  "Stopping at ",plotnum," plots:", irecord, plotnum, plotskip;              flag=2
              else
                 irecord=irecord+1
              endif
           endif
        endif
     enddo
     if ((flag==1).or.(flag==2)) then
        if (first==1) then
           first=0;           filename=filex
           giffilename="                    "
           if (inspf.lt.10) then
              write(filename,'(A9,I1,A1,I1,A4)') ignufile(iwhich),readprop,"_",inspf,".gnu"
              write(giffilename,'(A1,A9,I1,A1,I1,A5)') '"', ignufile(iwhich),readprop,"_",inspf,'.gif"'
           else
              write(filename,'(A9,I1,A1,I2,A4)') ignufile(iwhich),readprop,"_",inspf,".gnu"
              write(giffilename,'(A1,A9,I1,A1,I2,A5)') '"', ignufile(iwhich),readprop,"_",inspf,'.gif"'
           endif
           open(871,file=filename(1:getlen(filename)),status="unknown")
           call writergnuoptions(871, giffilename)
        endif
        if (plotterm.eq.3) then
           epsfilename="                    "
           if (inspf.lt.10) then
              write(epsfilename,'(A1,A9,I1,A1,I1,A5)') '"', ignufile(iwhich),readprop,"_",inspf,'.gif"'
           else
              write(epsfilename,'(A1,A9,I1,A1,I2,A5)') '"', ignufile(iwhich),readprop,"_",inspf,'.gif"'
           endif
           write(871,'(A11,A20)') " set out  ", epsfilename(1:getlen(epsfilename))
        endif
        call get_rnat_title(titleline, inspf, thistime, mydenval,iwhich)
        write(871,*) titleline
        write(871,*) " plot '-' using 1:2, '-' using 1:2"
        write(871,*)
        do ir=1,numr
           write(871,'(F8.3, F14.8)') real(bondpoints(ir),8),imag((0.d0,0.d0)+inrdenvect(ir)/sqrt(bondweights(ir)))
        enddo
        write(871,*) "e"
        do ir=1,numr
           write(871,'(F8.3, F14.8)') real(bondpoints(ir),8),real(inrdenvect(ir)/sqrt(bondweights(ir)),8)
        enddo
        write(871,*) "e";        write(871,*) "pause ", plotpause
     endif
  enddo
  close(871)

  call close_orbvector(ifilenums(iwhich))
  if (myrank.eq.1) then
     print *, "Plotted to file.  CAlling GNU!!!"
  endif

  sysline(1:100)=filex
  write(sysline,'(A18,A20,A1)') "gnuplot -persist ",filename(1:getlen(filename)),"&";  call system(sysline)
  close(rnatorbfile)

end subroutine read_Rorb


subroutine writegnuoptions(myunit, giffilename)
  use parameters
  implicit none
  integer :: myunit, getlen
  character (len=20) :: giffilename

  write(myunit,*) 
  write(myunit,*) " set contour"
  write(myunit,*) " set view ", plotview1," , ", plotview2 
  write(myunit,*) " set style data lines           "
  write(myunit,*) " set border 0"

  if (pm3d==1) then
     write(myunit,*) " set pm3d"
     if (plotterm.ne.3) then
        write(myunit,*) " unset surf                       "
     endif
  else
     write(myunit,*) " set surf                       "
  endif
  write(myunit,*) " set hidden3d"
  !  write(myunit,*) " set dgrid3d"
  write(myunit,*) " unset ztics; unset ytics; unset xtics"
  write(myunit,*) " set nokey                      "
  if (plotterm.ne.3) then
     write(myunit,*) " set label 'Real' at screen 0.2,0.8"
     write(myunit,*) " set label 'Imag' at screen 0.6,0.8"
     write(myunit,'(A14,F10.5,A1,F10.5,A1)') &
          " set zrange [",-plotrange,":",plotrange,"]"
  endif
  if (pm3d.ne.0) then
     if (plotterm.eq.3) then
        write(myunit,*) &
             " set cbrange [-3.14159 : 3.14159 ] "
     else
        write(myunit,'(A14,F10.5,A1,F10.5,A1)') &
             " set cbrange [",-plotcbrange,":",plotcbrange,"]"
     endif
  endif
  write(myunit,*) " set ticslevel 0"
  write(myunit,*) " unset clabel "
  write(myunit,'(A30,F17.10,A2,F17.10,A2,F17.10)') " set cntrparam levels incremental",-plotrange,",",plotrange/60,",",plotrange
  if (plotterm.ne.3) then
     write(myunit,*) " set xrange [",-plotxyrange,":",plotxyrange,"]"
     write(myunit,*) " set yrange [",-plotxyrange,":",3*plotxyrange,"]"
  else
#ifndef REALGO
     if (pm3d.ne.0) then
#endif
        write(myunit,*) " set xrange [",-plotxyrange,":",plotxyrange,"]"
        write(myunit,*) " set yrange [",-plotxyrange,":",plotxyrange,"]" 
#ifndef REALGO
     else
        write(myunit,*) " set xrange [",-plotxyrange,":",plotxyrange*3d0,"]"
        write(myunit,*) " set yrange [",-plotxyrange,":",plotxyrange,"]" 
     endif
#endif
  endif
  select case (plotterm)
  case (1)
     write(myunit,*) " set term aqua                   "
     write(myunit,*) " set out                        "
  case (2)
     write(myunit,*) " set term gif animate delay 10      "
     write(myunit,'(A11,A20)') " set out  ", giffilename(1:getlen(giffilename))
  case (3)
        write(myunit,*) "     set border 31 lw 2"
        write(myunit,*) "     unset ztics"
!! 032712 new version gnuplot        write(myunit,*) "      set term table"
!!        write(myunit,*) " set out 'points.dat'"
        write(myunit,*) "      set table 'points.dat'"
        write(myunit,*) " plot '-' "
#ifndef REALGO
if (pm3d.ne.0) then
#endif
        write(myunit,*) " -0.5  0.0"
        write(myunit,*) " 0.5 0.0"
        write(myunit,*) "e"
        write(myunit,*) "set size 1,1.25"
#ifndef REALGO
else
        write(myunit,*) "0.0 -0.5 "
        write(myunit,*) "0.0 0.5 "
        write(myunit,*) "e"
        write(myunit,*) "set size 2,1.4"
endif
#endif
     if (pm3d==1) then
        write(myunit,*) " set term post eps enhanced 20  color solid   "
        write(myunit,*) "     set view 0,0"
!!     else
!! 032712 dammit.  syntax change in gnuplot; 
!!      set out and set term together in set term table "filename".        
!!        write(myunit,*) " set term table"
!! ok!  not needed; redundant.
     endif
  case default
     write(myunit,*) " set term x11                   "
     write(myunit,*) " set out                        "
  end select
  write(myunit,*) 
  write(myunit,*) 

end subroutine writegnuoptions


subroutine writergnuoptions(myunit, giffilename)
  use parameters
  implicit none
  integer :: myunit, getlen
  character (len=20) :: giffilename
  write(myunit,*) 
  write(myunit,*) " set style data lines           "
  write(myunit,*) " set nokey                      "
  write(myunit,'(A20)') " set xrange [*:*]"
  write(myunit,'(A20)') " set yrange [*:*]"
  select case (plotterm)
  case (1)
     write(myunit,*) " set term aqua                   "
     write(myunit,*) " set out                        "
  case (2)
     write(myunit,*) " set term gif animate delay 10      "
     write(myunit,'(A11,A20)') " set out  ", giffilename(1:getlen(giffilename))
  case (3)
     write(myunit,*) " set term post eps enhanced 20  color solid   "
  case default
     write(myunit,*) " set term x11                   "
     write(myunit,*) " set out                        "
  end select
  write(myunit,*) 
  write(myunit,*) 
end subroutine writergnuoptions


subroutine whatsub(flag)
  use parameters
  implicit none
  character(len=6) :: termlabels(0:3) = (/ " x11  ", " aqua ", " .gif ", " .eps " /)
  integer :: flag
  character (len=2) :: inchar

  print *, "What should I do?"
  print *, "         s = change plotskip (now ",plotskip,")"
  print *, "         n = change plotnum (now ",plotnum,")"
     
  if (povplotflag /= 1) then
     print *, "         z = change zrange (now ",plotrange,")"
     print *, "         c = change cbrange (now ",plotcbrange,")"
     print *, "         x = change xyrange (now ",plotxyrange,")"
     print *, "         t = change terminal (now ",termlabels(plotterm),")"
     print *, "         d = change pm3d (now ",pm3d,")"
     print *, "         v1= change view rotation 1 (now ",plotview1,") degrees"
     print *, "         v2= change view rotation 2 (now ",plotview2,") degrees"
  endif
  print *, "   default = continue (plot with these options)"
  print *

  inchar="  ";  read(*,*) inchar
  if (povplotflag == 1) then
     select case (inchar(1:1))
     case ('s')
        print *, "   enter plotskip ( now ", plotskip, " ) "
        read(*,*) plotskip
     case ('n')
        print *, "   enter plotnum ( now ", plotnum, " ) "
        read(*,*) plotnum
     case default
        flag=1
     end select
  else
     select case (inchar(1:1))
     case ('s')
        print *, "   enter plotskip ( now ", plotskip, " ) "
        read(*,*) plotskip
     case ('n')
        print *, "   enter plotnum ( now ", plotnum, " ) "
        read(*,*) plotnum
     case ('z')
        print *, "   enter zrange ( now ", plotrange, " ) "
        read(*,*) plotrange
     case ('c')
        print *, "   enter cbrange ( now ", plotcbrange, " ) "
        read(*,*) plotcbrange
     case ('x')
        print *, "   enter xyrange ( now ", plotxyrange, " ) "
        read(*,*) plotxyrange
     case ('d')
        print *, "   enter pm3d ( now ", pm3d, " ) "
        read(*,*) pm3d
     case ('t')
        plotterm=mod(plotterm+1,4)
        print *, "Term now ", termlabels(plotterm)
     case ('v')
        select case (inchar(1:2))
        case ('v1')
           print *, "   enter view 1 ( now ", plotview1, " degrees ) "
           read(*,*) plotview1
        case ('v2')
           print *, "   enter view 2 ( now ", plotview2, " degrees ) "
           read(*,*) plotview2
        end select
     case default
        flag=1
     end select
  endif
end subroutine whatsub

