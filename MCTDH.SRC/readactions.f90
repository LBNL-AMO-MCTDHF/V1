
!! LOTS OF ACTIONS FOR READING ORBITALS FOR PLOTTING ROUTINES

#include "Definitions.INC"

!! iwhich=4, denproj, not debugged

subroutine read_orb_initial(iwhich)
  use parameters
  use mpimod
  implicit none
  integer ::  ispf,imvalue, iflag, jflag, iwhich, xxflag=0, jspf, iallflag
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
     print *, "                    PLOTTING MODE:"
     print *, "Enter -2 for 1D reduced; -1 for 1D slice; 0 for 2D reduced; 1 for 3D povray"
     read(*,*) plotmodeflag
     if (plotmodeflag==0) then
        print *, "OK, gnuplot, 2D slice"
     elseif (plotmodeflag==-1) then
        print *, "OK, gnuplot, 1D slice"
     elseif (plotmodeflag==-2) then
        print *, "Ok, gnuplot, 1D reduced"
     else
        print *, "OK, using povray."
        plotmodeflag=1
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
        print *, "Enter ", labels(iwhich), " number.  Zero to change options.  Negative number to stop."
        read(*,*) ispf
        
     else

        !!        print *, "Enter positive number to continue to plot.  Zero to change options.  Negative number to stop."
        !!        read(*,*) ispf
        !!        ispf=min(1,ispf)

        ispf = 1
     endif

!$$     iprop=1
!$$     if (numprop.gt.2) then
!$$        print *, "Which calculation?  There are ", numprop," concurrent."
!$$        read(*,*) iprop
!$$     endif

     imvalue=0
     if (ispf<0) then
        print *, "OK! Done."
        iflag=0
     else if (ispf>0) then
        iallflag=0
        if (plotmodeflag==0) then
           if (iwhich.ne.3) then
              if (spfrestrictflag.ne.0) then
                 imvalue=spfmvals(ispf);              
                 print *, "m-value of this orbital is ", imvalue
              else
                 print *, "Enter m projection of orbital to plot"
                 read(*,*) imvalue
              endif
           else
              print *, "Enter m projection of density to plot"
              read(*,*) imvalue
           endif
        endif

        if (ispf>nspf) then
           print *, "Will do all!";           iallflag=1
        endif

        if (iallflag==1) then
           do jspf=1,nspf
              if (spfrestrictflag==1) then
                 imvalue=spfmvals(jspf)
              endif
              print *, "imvalue is ", imvalue
              call read_orb(jspf,imvalue, 1, iwhich)
           enddo
        else
           call read_orb(ispf,imvalue, 1, iwhich)           
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

!$$     iprop=1
!$$     if (numprop.gt.2) then
!$$        print *, "Which calculation?  There are ", numprop," concurrent."
!$$        read(*,*) iprop
!$$     endif

     if (ispf.gt.numr) then
        print *, "Um, that is greater than numr. goodbye.";        stop
     endif
     if (ispf<0) then
        print *, "OK! Done.";        iflag=0
     else if (ispf>0) then
!$$        call read_Rorb(ispf, iprop,iwhich)
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

subroutine read_orb(whichspf,imvalue, iprop, iwhich)   
  use parameters
  use mpimod
  use zerorecmod
  implicit none

  integer, intent(in) :: iwhich,imvalue,iprop,whichspf
  integer :: readprop,  returnval, ixi, ieta, flag, xflag, irecord, first, &
       getlen, ispf, itable, tabmax, iz
  integer, parameter :: ifilenums(4)=(/  natorbfile, spfplotfile, denfile, denprojfile /)
  character (len=SLN) :: ifilenames(4)
  character (len=8), parameter :: ignufile(4)=(/ "Natplot_", "Spfplot_", "Denplot_", "Denproj_" /)
  character (len=10) :: povdirs(4) = (/ "Natorb    ", "Spfs      ", "Density   ", "Denproj   " /)
  character(len=10) :: mypovdir
  real*8 :: thistime, zval, xval,  costheta, rho, rsum
  DATATYPE :: sum, cylindricalvalue, mydenval
  DATATYPE :: ttempspf(spfdims(1),spfdims(2),spfdims(3))      !! AUTOMATIC
  DATATYPE :: savespf(spfdims(1),spfdims(2),spfdims(3))       !! AUTOMATIC
  character (len=8) :: tableexts(4)=(/'.table0"','.table1"','.table2"','.table3"' /)
  character (len=100) :: filename,sysline, titleline
  character(len=20) :: giffilename
  character(len=40) :: epsfilename,epstable(4)
  character(len=100) :: filex="                                                                                                    "
  character (len=headersize) :: header

  ttempspf=0; savespf=0;

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
        if (iwhich.eq.1) then  !! natorb
           call fixphase0(ttempspf(:,:,:),spfsize)
        endif
        if (plotsubtract.ne.0) then
           if (first==1) then
              savespf(:,:,:) = ttempspf(:,:,:)
              ttempspf=0
           else
              ttempspf(:,:,:) = ttempspf(:,:,:) - savespf(:,:,:)
           endif
        endif
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
        if (returnval/=0) then
           call openfile();   write(mpifileptr,*) "Done with vectors on file.";  call closefile()
           if (plotmodeflag.ne.1) then   !! if not povray
              close(871)   
           endif
           flag=3
        else
           xflag=0
           select case (iwhich)
           case (1,2,4)
              if ((whichspf==ispf).and.(readprop.eq.iprop)) then
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
        if (plotmodeflag==1) then
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
                 if (whichspf.lt.10) then
                    write(filename,'(A8,I1,A1,I1,A4)') ignufile(iwhich),readprop,"_",whichspf,".gnu"
                    if (plotterm.eq.3) then
                       write(giffilename,'(A1,A8,I1,A1,I1,A5)') '"',ignufile(iwhich),readprop,"_",whichspf,'.eps"'
                    else
                       write(giffilename,'(A1,A8,I1,A1,I1,A5)') '"',ignufile(iwhich),readprop,"_",whichspf,'.gif"'
                    endif
                 else
                    write(filename,'(A8,I1,A1,I2,A4)') ignufile(iwhich),readprop,"_",whichspf,".gnu"
                    write(giffilename,'(A1,A8,I1,A1,I2,A5)') '"',ignufile(iwhich),readprop,"_",whichspf,'.gif"'
                 endif
              endif
              open(871,file=filename(1:getlen(filename)),status="unknown")
              if (plotmodeflag==0) then
                 call writegnuoptions3d(871, giffilename)
              elseif (plotmodeflag==-1 .or. plotmodeflag==-2) then
                 call writegnuoptions1d(871, giffilename)
              else
                 OFLWR "OOPS error programmer fail gnuoptions"; CFLST
              endif

           endif  !! ifirst

           if (plotterm.eq.3.and.pm3d.eq.0.and.plotmodeflag.ge.0) then

              if (plotmodeflag.ne.0) then
                 OFLWR "PROGRAMMER FAAAAIL"; CFLST
              endif
#ifdef REALGO
              tabmax = 2
#else
              tabmax = 4
#endif              
              do itable=1,tabmax
                 epstable(itable)="                                        "
                 if (iwhich.eq.3) then
                    write(epstable(itable),'(A1,A8,I1,A8)') '"',ignufile(iwhich),readprop,tableexts(itable)
                 else
                    if (whichspf.lt.10) then
                       write(epstable(itable),'(A1,A8,I1,A1,I1,A1,A4,A8)') '"',ignufile(iwhich),readprop,"_",whichspf,'_', zerorecord(irecord),tableexts(itable)
                    else
                       write(epstable(itable),'(A1,A8,I1,A1,I2,A1,A4,A8)') '"',ignufile(iwhich),readprop,"_",whichspf,'_', zerorecord(irecord),tableexts(itable)
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
                       
                       !!                       sum=cylindricalvalue(rho,costheta,1.d0,imvalue, ttempspf(:,:,imvalue+(spfdims(3)+1)/2))
                       sum=cylindricalvalue(rho,costheta,1.d0,imvalue, ttempspf)
                       if (zval.lt.0.d0) then
                          sum=sum*(-1)**abs(imvalue)
                       endif
#ifndef REALGO
                       if (itable.lt.3) then
#endif
                          write(871,'(2F8.3, F14.8)') zval,xval,real(sum)
#ifndef REALGO
                       else
                          write(871,'(2F8.3, F14.8)') zval,xval,imag(sum+(0d0,0d0))
                       endif
#endif                       
                    enddo
                    write(871,*)
                 enddo  !! ixi
                 write(871,*) "e"
                 write(871,*)
                 write(871,*) "unset table"
                 write(871,*)
              enddo !! itable
           endif  !! ( plotterm.eq.3 .and. pm3d.eq.0 .and. plotmodeflag.ge.0 )

           select case(iwhich)
           case (1)
              call get_natorb_title(titleline, whichspf, imvalue, thistime, mydenval)
           case (2)
              call get_spf_title(titleline, whichspf, imvalue, thistime)
           case (3)
              call get_density_title(titleline, thistime)
           case (4)
              call get_natorb_title(titleline, whichspf, imvalue, thistime, mydenval)
           end select

           write(871,*) titleline

           if (plotterm.eq.3.and.plotmodeflag.ge.0) then
              if (plotmodeflag.ne.0) then
                 OFLWR "Rogrammer fail"; CFLST
              endif
              epsfilename="                                        "
              if (iwhich.eq.3) then
                 write(epsfilename,'(A1,A8,I1,A5)') '"',ignufile(iwhich),readprop,'.eps"'
              else
                 if (whichspf.lt.10) then
                    write(epsfilename,'(A1,A8,I1,A1,I1,A1,A4,A5)') '"',ignufile(iwhich),readprop,"_",whichspf,'_', zerorecord(irecord),'.eps"'
                 else
                    write(epsfilename,'(A1,A8,I1,A1,I2,A1,A4,A5)') '"',ignufile(iwhich),readprop,"_",whichspf,'_', zerorecord(irecord),'.eps"'
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
           else !! plotterm.eq.3 .and. plotmodeflag.ge.0
              if (plotmodeflag.eq.0) then
                 if (pm3d==1) then
                    write(871,*) " splot '-' using 1:($2+2*",plotxyrange,"):3:3,'-' using 1:2:3:3"
                 else
                    write(871,*) " splot '-' using 1:($2+2*",plotxyrange,"):3,'-' using 1:2:3"
                 endif
              elseif (plotmodeflag.eq.-1 .or. plotmodeflag.eq.-2) then
                 write(871,*) " plot '-' using 1:2 title 'real','-' using 1:2 title 'imag'"
              else
                 OFLWR "programmaerrrr fail"; CFLST
              endif
           endif !! plotterm eq 3 and plotmodeflag ne -1
           write(871,*)
           
           if (plotmodeflag.eq.0) then
              if (plotterm.ne.3.or.pm3d.ne.0) then  !! we want to plot out some data otherwise term 3, 3d 0 and we have tables and we are done
                 do ixi=(-1)*plotres,plotres-1
                    do ieta=(-1)*plotres,plotres-1
                       xval=(ixi+0.5d0)*plotxyrange/plotres
                       zval=(ieta+0.5d0)*plotxyrange/plotres
                       rho=sqrt(xval**2+zval**2)
                       costheta=cos(atan2(zval,xval))

                       sum=cylindricalvalue(rho,costheta,1.d0,imvalue, ttempspf)
                       
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

                          sum=cylindricalvalue(rho,costheta,1.d0,imvalue, ttempspf)
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

           elseif (plotmodeflag.eq.-2) then
              rsum=0
              do iz = (-1)*plotres,plotres-1
                 zval=(iz+0.5d0)*plotxyrange/plotres
                 sum=0d0
                 do ixi = 1,500             !! HARDWIRE 500 points for now
                    xval = (ixi-0.5)/10;    !! HARDWIRE 1/10 step for now

                    rho=sqrt(xval**2+zval**2)
                    costheta=zval / rho;  !!  cos(atan2(zval,xval))

                    sum = sum + 2 * pi * xval / 10 * cylindricalvalue(rho,costheta,1.d0,0, ttempspf)
                 enddo
                 rsum = rsum + plotxyrange/plotres * real(sum,8)
                 write(871,'(F8.3, 2F14.8)') zval,real(sum,8), rsum
              enddo
              write(871,*) "e"
              rsum=0
              do iz = (-1)*plotres,plotres-1
                 zval=(iz+0.5d0)*plotxyrange/plotres
                 rho=abs(zval)
                 costheta = sign(1d0,zval)

                 sum=0d0
                 do ixi = 1,500             !! HARDWIRE 500 points for now
                    xval = (ixi-0.5)/10;    !! HARDWIRE 1/10 step for now

                    rho=sqrt(xval**2+zval**2)
                    costheta=zval / rho;  !!  cos(atan2(zval,xval))

                    sum = sum + 2 * pi * xval / 10 * cylindricalvalue(rho,costheta,1.d0,0, ttempspf)
                 enddo
                 rsum = rsum + plotxyrange/plotres * imag(sum+(0d0,0d0))
                 write(871,'(F8.3, 2F14.8)') zval,imag(sum+(0d0,0d0)),rsum
              enddo
              write(871,*) "e"
              write(871,*) "pause ", plotpause


           elseif (plotmodeflag.eq.-1) then
              rsum=0d0
              do iz = (-1)*plotres,plotres-1
                 zval=(iz+0.5d0)*plotxyrange/plotres
                 rho=abs(zval)
                 costheta = sign(1d0,zval)

                 sum=cylindricalvalue(rho,costheta,1.d0,0, ttempspf)

                 rsum = rsum + plotxyrange/plotres * real(sum,8)
                 write(871,'(F8.3, 2F14.8)') zval,real(sum,8),rsum
              enddo
              write(871,*) "e"
              rsum=0d0
              do iz = (-1)*plotres,plotres-1
                 zval=(iz+0.5d0)*plotxyrange/plotres
                 rho=abs(zval)
                 costheta = sign(1d0,zval)

                 sum=cylindricalvalue(rho,costheta,1.d0,0, ttempspf)
                 rsum = rsum + plotxyrange/plotres * imag(sum+(0d0,0d0))
                 write(871,'(F8.3, 2F14.8)') zval,imag(sum+(0d0,0d0)),rsum
              enddo
              write(871,*) "e"
              write(871,*) "pause ", plotpause
           else
              OFLWR "PROGFAIL"; CFLST
           endif  !! plotmodeflag.eq.0 or -1
        endif  !! plotmodeflag.eq.1
     endif !! flag==1 or 2
     call checkstopfile()
  enddo

  call close_orbvector(ifilenums(iwhich))
  if (plotmodeflag.ne.1) then
     close(871)

     if (myrank.eq.1) then
        print *, "Plotted to file.  CAlling GNU!!!"
     endif
     sysline(1:100)=filex

     !!     write(sysline,'(A18,A20,A1)') "gnuplot -persist ",filename,"&"
     !!     call system(sysline)

     write(sysline,'(A18,A20)') "gnuplot -persist ",filename
     call system(sysline)

  endif

end subroutine read_orb


subroutine read_Rorb(whichspf, iprop,iwhich)
  use parameters
  use mpimod
  implicit none

  integer :: readprop,   returnval, ir, whichspf, flag, irecord, first, iprop, getlen, iwhich, ispf
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
           print *,  "Done with vectors on file for rden= ", whichspf
           close(871)
           flag=3
        else
           if ((whichspf==ispf).and.(readprop.eq.iprop)) then
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
           if (whichspf.lt.10) then
              write(filename,'(A9,I1,A1,I1,A4)') ignufile(iwhich),readprop,"_",whichspf,".gnu"
              write(giffilename,'(A1,A9,I1,A1,I1,A5)') '"', ignufile(iwhich),readprop,"_",whichspf,'.gif"'
           else
              write(filename,'(A9,I1,A1,I2,A4)') ignufile(iwhich),readprop,"_",whichspf,".gnu"
              write(giffilename,'(A1,A9,I1,A1,I2,A5)') '"', ignufile(iwhich),readprop,"_",whichspf,'.gif"'
           endif
           open(871,file=filename(1:getlen(filename)),status="unknown")
           call writergnuoptions(871, giffilename)
        endif
        if (plotterm.eq.3) then
           epsfilename="                    "
           if (whichspf.lt.10) then
              write(epsfilename,'(A1,A9,I1,A1,I1,A5)') '"', ignufile(iwhich),readprop,"_",whichspf,'.gif"'
           else
              write(epsfilename,'(A1,A9,I1,A1,I2,A5)') '"', ignufile(iwhich),readprop,"_",whichspf,'.gif"'
           endif
           write(871,'(A11,A20)') " set out  ", epsfilename(1:getlen(epsfilename))
        endif
        call get_rnat_title(titleline, whichspf, thistime, mydenval,iwhich)
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
        write(871,*) "e"
        write(871,*) "pause ", plotpause
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


subroutine writegnuoptions3d(myunit, giffilename)
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
     if (plotrange < 0) then
        write(myunit,*) " set zrange [*:*]"
     else
        write(myunit,'(A14,F10.5,A1,F10.5,A1)') " set zrange [",-plotrange,":",plotrange,"]"
     endif
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
  if (plotrange > 0) then
     write(myunit,'(A30,F17.10,A2,F17.10,A2,F17.10)') " set cntrparam levels incremental",-plotrange,",",plotrange/60,",",plotrange
  endif
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

end subroutine writegnuoptions3d


subroutine writegnuoptions1d(myunit, giffilename)
  use parameters
  implicit none
  integer :: myunit, getlen
  character (len=20) :: giffilename

  write(myunit,*) 
  write(myunit,*) " set style data lines           "
  write(myunit,*) "     set border 31 lw 2"

  write(myunit,*) " set xrange [",-plotxyrange,":",plotxyrange,"]"

  select case (plotterm)
  case (1)
     write(myunit,*) " set term aqua                   "
     write(myunit,*) " set out                        "
  case (2)
     write(myunit,*) " set term gif animate delay 10      "
     write(myunit,'(A11,A20)') " set out  ", giffilename(1:getlen(giffilename))
  case (3)
  case default
     write(myunit,*) " set term x11                   "
     write(myunit,*) " set out                        "
  end select
  write(myunit,*) 
  write(myunit,*) 

end subroutine writegnuoptions1d


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
  print *, "         q = quit"
  print *, "         s = change plotskip (now ",plotskip,")"
  print *, "         b = toggle plotsubtract (now ",plotsubtract,")"
  print *, "         n = change plotnum (now ",plotnum,")"
     
  if (plotmodeflag /= 1) then  !! not povray
     print *, "         r = change plotres (now ",plotres,")"
     print *, "         x = change xy (abscissa) range (now ",plotxyrange,")"
     if (plotrange<0) then
        print *, "         z = change zrange (now auto)"
     else
        print *, "         z = change zrange (now ",plotrange,")"
     endif
     print *, "         t = change terminal (now ",termlabels(plotterm),")"
  endif
  if (plotmodeflag == 0) then  !! 2D plot
     print *, "         d = change pm3d (now ",pm3d,")"
     if (pm3d.ne.0) then
        print *, "         c = change cbrange (now ",plotcbrange,")"
     endif
     print *, "         v1= change view rotation 1 (now ",plotview1,") degrees"
     print *, "         v2= change view rotation 2 (now ",plotview2,") degrees"
  endif
  print *,    "         e = continue (plot with these options)"
  print *

  inchar="  ";  read(*,*) inchar
  if (plotmodeflag == 1) then    !! povray
     select case (inchar(1:1))
     case('q')
        print *, "Ok, quitting"
        call mpistop()
     case ('s')
        print *, "   enter plotskip ( now ", plotskip, " ) "
        read(*,*) plotskip
     case ('b')
        if (plotsubtract.eq.0) then
           plotsubtract=1
        else
           plotsubtract=0
        endif
     case ('n')
        print *, "   enter plotnum ( now ", plotnum, " ) "
        read(*,*) plotnum
     case default
        flag=1
     end select
  else
     select case (inchar(1:1))
     case('q')
        print *, "Ok, quitting"
        call mpistop()
     case ('s')
        print *, "   enter plotskip ( now ", plotskip, " ) "
        read(*,*) plotskip
     case ('b')
        if (plotsubtract.eq.0) then
           plotsubtract=1
        else
           plotsubtract=0
        endif
     case ('n')
        print *, "   enter plotnum ( now ", plotnum, " ) "
        read(*,*) plotnum
     case ('r')
        print *, "   enter plotres ( now ", plotres, " ) "
        read(*,*) plotres
     case ('z')
        if (plotrange<0) then
           print *, "   enter zrange, negative for auto ( now auto ) "
        else
           print *, "   enter zrange, negative for auto ( now ", plotrange, " ) "
        endif
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

