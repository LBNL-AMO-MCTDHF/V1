
!! FOR POVRAY PLOTTING.

#include "Definitions.INC"

module povmod
  implicit none
  real*8, allocatable :: xvals(:,:), yvals(:,:), zvals(:,:)
  complex*16, allocatable :: sparsetransmat(:,:)
  integer, allocatable :: sparseend(:,:,:,:),sparsexyz(:,:,:),sparsestart(:,:,:,:)
  integer :: xxcall=0,maxsparse=0
end module

subroutine povinit()
  use povmod
  use parameters
  implicit none
  integer :: ix, isparse, irange
  real*8 :: sum
  if (xxcall==0) then
     allocate(sparsestart(spfdims(1),spfdims(2),-(spfdims(3)-1)/2:(spfdims(3)-1)/2,numpovranges))
     allocate(sparseend(spfdims(1),spfdims(2),-(spfdims(3)-1)/2:(spfdims(3)-1)/2,numpovranges))
     allocate( xvals(-povres:povres,numpovranges),yvals(-povres:povres,numpovranges))
     allocate( zvals(-povres:povres,numpovranges))
     do irange=1,numpovranges
        do ix=-povres,povres
           yvals(ix,irange)=ix+povres
        enddo
        sum=yvals(povres,irange)
        yvals(:,irange)=(yvals(:,irange)/sum*2.d0 -1.d0)*povrange(irange)
        yvals(:,irange)=yvals(:,irange) - (yvals(1,irange)-yvals(0,irange))/2.d0
        xvals(:,irange)=yvals(:,irange)
        zvals(:,irange)=yvals(:,irange)
     enddo

     OFLWR "  Get maxsparse for prolate spheroidal/spherical transformation matrix....";CFL
     maxsparse=0
     do irange=1,numpovranges
        call get_maxsparse(2*povres+1,2*povres+1,2*povres+1,zvals(:,irange),yvals(:,irange),xvals(:,irange), isparse,povsparse)
        OFL; write(mpifileptr,'(A18,I16,F30.0)') "Got maxsparse: ", isparse, real(2*povres+1,8)**3*spfsize;CFL
        if (isparse.gt.maxsparse) then
           maxsparse=isparse
        endif
     enddo
     OFLWR "  Get prolate spheroidal/spherical transformation matrix.";CFL
     allocate (sparsexyz(maxsparse,3, numpovranges),sparsetransmat(maxsparse, numpovranges))
     sparsexyz=-999
     do irange=1,numpovranges
        call get_sphericalsparse(2*povres+1,2*povres+1,2*povres+1,zvals(:,irange),yvals(:,irange),xvals(:,irange), maxsparse, sparsetransmat(:,irange), sparsestart(:,:,:,irange), sparseend(:,:,:,irange), sparsexyz(:,:,irange), povsparse)        
     enddo
     sparsexyz=sparsexyz-1-povres
     OFLWR "     ...done.";CFL
  endif
  xxcall=1

end subroutine povinit



!! iwhich : 1 for nat, 2 for spf, 3 for den, 4 for denproj,   denproj not debugged.

subroutine read_povother( thistime, inspf, iprop, dirlabel, iwhich, ispf, irecord)
  use parameters
  use povmod
  implicit none

  DATATYPE :: inspf(spfdims(1),spfdims(2),-(spfdims(3)-1)/2:(spfdims(3)-1)/2)
  character :: dirlabel*(*)
  integer ::  iprop, getlen2, iwhich, itotpov,   irange, irecord, imval, lvalue, ixi, kkk, ispf, yycall(4)=0, &
       ix,iy,iz, jjj,iii, iilen
  real*8 :: mmax(10),thistime
  integer, save :: maxflag(40,4,10)=0
  real*8, save :: plotmax(40,4,10) 
  character (len=100) :: filename, filename2, filename3
  character (len=100) :: filex="                                        "
  character(len=30) :: iinum, iinum2, iinum3, iinum1,tempchar, timechar
  character (len=2) ::  mychar(20)= (/ "1 ","2 ","3 ","4 ","5 ","6 ","7 ","8 ","9 ","10", "11","12","13","14","15","16","17","18","19","20" /)
  real*8, allocatable  :: threedensity3(:,:,:),threedensity1(:,:,:),threedensity2(:,:,:)
  complex*16, allocatable  :: cthreedensity(:,:,:)
  integer, parameter :: zerolen=4
  character (len=10) :: zeroes = "0000000000"

  allocate( &
       threedensity1(-povres:povres, -povres:povres, -povres:povres), &
       threedensity2(-povres:povres, -povres:povres, -povres:povres), &
       threedensity3(-povres:povres, -povres:povres, -povres:povres), &
       cthreedensity(-povres:povres, -povres:povres, -povres:povres))
  yycall(iwhich)=yycall(iwhich)+1

#ifdef NOCFLAG
  if (yycall(iwhich).eq.1) then
     OFLWR "NOCFLAG turned on in compilation; no povray output possible.";CFL
  endif
#else
  if (yycall(iwhich).eq.1) then
     call system("mkdir "//dirlabel);     call system("mkdir "//dirlabel//"/df3")
  endif
  itotpov=(2*povres+1)**3

  call povinit()

!  if ((numprop.gt.20).or.(nspf.gt.40).or.(iwhich.gt.4)) then
!     call openfile()
!     write(mpifileptr,*) "REDIM"
!     call closefile();     call mpistop()
!  endif

  do irange=1,numpovranges
     OFLWR " POV-plotting "//dirlabel//" !! ", irecord, ispf, povrange(irange);CFL
     filename=filex
     filename2=filex
     filename3=filex
     iinum="                              "
     iinum1="                              "
     iinum2="                              "
     iinum3="                              "

if (1==0) then
     
     write(filename,'(A'//mychar(LEN(dirlabel)+5)//',I'//mychar(iilen(ispf))//',A2,I'//mychar(iilen(irecord))//',A5,I1,A6,I1,A6)') dirlabel//"/df3/", ispf,"_T",irecord,"_prop",iprop,"_range",irange,"_r.df3"
     
     write(iinum,'(I'//mychar(iilen(ispf))//',A2,I'//mychar(iilen(irecord))//',A5,I1,A6,I1)')  ispf,"_T",irecord,"_prop",iprop,"_range",irange
     
     write(iinum1,'(I'//mychar(iilen(ispf))//',A2,I'//mychar(iilen(irecord))//',A5,I1,A6,I1,A2)')  ispf,"_T",irecord,"_prop",iprop,"_range",irange,"_r"

     
     write(filename2,'(A'//mychar(LEN(dirlabel)+5)//',I'//mychar(iilen(ispf))//',A2,I'//mychar(iilen(irecord))//',A5,I1,A6,I1,A7)') dirlabel//"/df3/", ispf,"_T",irecord,"_prop",iprop,"_range",irange,"_nr.df3"
     
     write(iinum2,'(I'//mychar(iilen(ispf))//',A2,I'//mychar(iilen(irecord))//',A5,I1,A6,I1,A3)')  ispf,"_T",irecord,"_prop",iprop,"_range",irange,"_nr"
     
     write(filename3,'(A'//mychar(LEN(dirlabel)+5)//',I'//mychar(iilen(ispf))//',A2,I'//mychar(iilen(irecord))//',A5,I1,A6,I1,A6)') dirlabel//"/df3/", ispf,"_T",irecord,"_prop",iprop,"_range",irange,"_i.df3"
     
     write(iinum3,'(I'//mychar(iilen(ispf))//',A2,I'//mychar(iilen(irecord))//',A5,I1,A6,I1,A2)')  ispf,"_T",irecord,"_prop",iprop,"_range",irange,"_i"
     

else

     write(filename,'(A'//mychar(LEN(dirlabel)+5)//',I'//mychar(iilen(ispf))//',A'//mychar(2+zerolen-iilen(irecord))//',I'//mychar(iilen(irecord))//',A5,I1,A6,I1,A6)') dirlabel//"/df3/", &
     ispf,"_T"//zeroes(1:zerolen-iilen(irecord)),irecord,"_prop",iprop,"_range",irange,"_r.df3"
     
     write(iinum,'(I'//mychar(iilen(ispf))//',A'//mychar(2+zerolen-iilen(irecord))//',I'//mychar(iilen(irecord))//',A5,I1,A6,I1)')  ispf,"_T"//zeroes(1:zerolen-iilen(irecord)),irecord,"_prop",iprop,"_range",irange
     
     write(iinum1,'(I'//mychar(iilen(ispf))//',A'//mychar(2+zerolen-iilen(irecord))//',I'//mychar(iilen(irecord))//',A5,I1,A6,I1,A2)')  ispf,"_T"//zeroes(1:zerolen-iilen(irecord)),irecord,"_prop",iprop,"_range",irange,"_r"
     
     write(filename2,'(A'//mychar(LEN(dirlabel)+5)//',I'//mychar(iilen(ispf))//',A'//mychar(2+zerolen-iilen(irecord))//',I'//mychar(iilen(irecord))//',A5,I1,A6,I1,A7)') dirlabel//"/df3/", &
     ispf,"_T"//zeroes(1:zerolen-iilen(irecord)),irecord,"_prop",iprop,"_range",irange,"_nr.df3"
     
     write(iinum2,'(I'//mychar(iilen(ispf))//',A'//mychar(2+zerolen-iilen(irecord))//',I'//mychar(iilen(irecord))//',A5,I1,A6,I1,A3)')  ispf,"_T"//zeroes(1:zerolen-iilen(irecord)),irecord,"_prop",iprop,"_range",irange,"_nr"
     
     write(filename3,'(A'//mychar(LEN(dirlabel)+5)//',I'//mychar(iilen(ispf))//',A'//mychar(2+zerolen-iilen(irecord))//',I'//mychar(iilen(irecord))//',A5,I1,A6,I1,A6)') dirlabel//"/df3/", &
     ispf,"_T"//zeroes(1:zerolen-iilen(irecord)),irecord,"_prop",iprop,"_range",irange,"_i.df3"
     
     write(iinum3,'(I'//mychar(iilen(ispf))//',A'//mychar(2+zerolen-iilen(irecord))//',I'//mychar(iilen(irecord))//',A5,I1,A6,I1,A2)')  ispf,"_T"//zeroes(1:zerolen-iilen(irecord)),irecord,"_prop",iprop,"_range",irange,"_i"
endif     

     cthreedensity=0.d0
     do imval=-(spfdims(3)-1)/2,(spfdims(3)-1)/2
        do lvalue=1,spfdims(2)
           do ixi=1,spfdims(1)
              do iii=sparsestart(ixi,lvalue,imval,irange), sparseend(ixi,lvalue,imval,irange)
                 cthreedensity(sparsexyz(iii,3,irange), sparsexyz(iii,2,irange), sparsexyz(iii,1,irange)) = &
                      cthreedensity(sparsexyz(iii,3,irange), sparsexyz(iii,2,irange), sparsexyz(iii,1,irange)) + &
                      sparsetransmat(iii,irange) * inspf(ixi,lvalue,imval)
              enddo
           enddo
        enddo
     enddo
     
     if (maxflag(ispf,iwhich,irange)==0) then
        mmax(1)=0.0;        kkk=0
        do ix=-povres,povres
           do iy=-povres,povres
              do iz=-povres,povres
                 if (abs(cthreedensity(iz,iy,ix)).gt.mmax(1)) then
                    mmax(2:10)=mmax(1:9)
                    mmax(1)=abs(cthreedensity(iz,iy,ix))
                    kkk=kkk+1
                 endif
              enddo
           enddo
        enddo
        plotmax(ispf,iwhich,irange)=mmax(min(kkk,1))
        maxflag(ispf,iwhich,irange)=1
     endif

     threedensity1=real(cthreedensity)/plotmax(ispf,iwhich,irange) * povmult
     threedensity2=real(cthreedensity*exp((0.d0,1.d0)*2.d0*pi/3.d0))/plotmax(ispf,iwhich,irange) * povmult
     threedensity3=real(cthreedensity*exp((0.d0,1.d0)*4.d0*pi/3.d0))/plotmax(ispf,iwhich,irange) * povmult

     jjj=2*povres+1
     call writepovray(jjj,threedensity1,filename(1:getlen2(filename)))
     call writepovray(jjj,threedensity2,filename2(1:getlen2(filename2)))
     call writepovray(jjj,threedensity3,filename3(1:getlen2(filename3)))

#ifndef NOPOVRAY
     write(tempchar,'(F10.4)') povrange(irange)
     write(timechar,'(F12.4)') thistime *0.02418884d0

     open(809, file="Density.Bat", status="old", iostat=iii)
     if ( iii /=0 ) then
        if (xxcall==1) then
           OFLWR; WRFL "  **  I did not find Density.Bat in the directory, thus povray will not be called to render."
           WRFL; CFL
        endif
     else
        close(809)
        if (iwhich.eq.3) then
           call system("./Density.Bat "//dirlabel//" "//iinum1//" "//dirlabel//"_"//iinum//" "//tempchar//" "//timechar)
           call system('echo "./Density.Bat '//dirlabel//' '//iinum1//' '//dirlabel//'_'//iinum//' '//tempchar//' '//timechar//'">lastcall.bat')

        else
           call system("./Density.Bat "//dirlabel//" "//iinum1//" "//iinum2//" "//iinum3//" "//dirlabel//"_"//iinum//" "//tempchar//" "//timechar)
           call system('echo "./Density.Bat '//dirlabel//' '//iinum1//' '//iinum2//' '//iinum3//' '//dirlabel//'_'//iinum//' '//tempchar//' '//timechar//'">lastcall2.bat')

        endif
     endif
#endif
  enddo
  OFLWR "     ...povother done";CFL
#endif
  deallocate(      threedensity1,      threedensity2,      threedensity3,      cthreedensity)

end subroutine read_povother
