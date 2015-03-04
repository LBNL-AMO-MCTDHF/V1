
#include "Definitions.INC"

subroutine save_natproj_final()
  use parameters
  implicit none
  call close_orbvector(natprojfile)
end subroutine save_natproj_final

subroutine save_natproj( thistime )
  use parameters
  use natprojmod
  use xxxmod !! denproj
  implicit none

  DATATYPE :: csum, dot
  integer :: xcalledhere=0,  iflag, korder, pulseflag, isplit, i, imc, getlen2, ilen, iii, nn, getlen, iorb
  real*8 :: thistime, rsum
  real*8, save :: mmax,mmin
  character (len=headersize) :: header
  character (len=40) :: filename,filename2   !!, sysline
  character(len=40) :: filex="                                        "

  character(len=256) :: cwdline  !!, cwd 
  character(len=256) :: nnbb="                                                                                                &
&                                                                                                                                                                "
  DATATYPE :: tempdenmat(nspf,nspf,4),tempinvden(nspf,nspf),tempdenvects(nspf,nspf)
  CNORMTYPE :: tempdenvals(nspf)

  xcalledhere=xcalledhere+1
  if (xcalledhere==1) then
     call system("mkdir NatCurves")
  endif

!  do iorb=1,numr
!     call fixphase0(natproj(:,iorb),numr)
!  enddo

  do imc=1,mcscfnum
     do iorb=1,numr
        csum=natvals(iorb)
        call write_nat_header(header,thistime,imc,iorb,csum)
        call save_orbvector(natproj(:,iorb,imc), numr, natprojfile,  natprojplotbin, header)
     enddo
  enddo

  do pulseflag=0,1  !pulse
     filename=filex;     filename2=filex
     if (pulseflag.eq.0) then
        iflag=0
        if (xcalledhere.lt.10) then
           write(filename,'(A10,I1,A5)') "NatCurves/",xcalledhere,".dat0"
           write(filename2,'(A10,I1,A6)') "NatCurves/",xcalledhere,"W.dat0"
           iflag=1
        else if (xcalledhere.lt.100) then
           write(filename,'(A10,I2,A5)') "NatCurves/",xcalledhere,".dat0"
           write(filename2,'(A10,I2,A6)') "NatCurves/",xcalledhere,"W.dat0"
           iflag=1
        else if (xcalledhere.lt.1000) then
           write(filename,'(A10,I3,A5)') "NatCurves/",xcalledhere,".dat0"
           write(filename2,'(A10,I3,A6)') "NatCurves/",xcalledhere,"W.dat0"
           iflag=1
        else if (xcalledhere.lt.10000) then
           write(filename,'(A10,I4,A5)') "NatCurves/",xcalledhere,".dat0"
           write(filename2,'(A10,I4,A6)') "NatCurves/",xcalledhere,"W.dat0"
           iflag=1
        endif
     else
        iflag=0
        if (xcalledhere.lt.10) then
           write(filename,'(A10,I1,A4)') "NatCurves/",xcalledhere,".dat"
           write(filename2,'(A10,I1,A5)') "NatCurves/",xcalledhere,"W.dat"
           iflag=1
        else if (xcalledhere.lt.100) then
           write(filename,'(A10,I2,A4)') "NatCurves/",xcalledhere,".dat"
           write(filename2,'(A10,I2,A5)') "NatCurves/",xcalledhere,"W.dat"
           iflag=1
        else if (xcalledhere.lt.1000) then
           write(filename,'(A10,I3,A4)') "NatCurves/",xcalledhere,".dat"
           write(filename2,'(A10,I3,A5)') "NatCurves/",xcalledhere,"W.dat"
           iflag=1
        else if (xcalledhere.lt.10000) then
           write(filename,'(A10,I4,A4)') "NatCurves/",xcalledhere,".dat"
           write(filename2,'(A10,I4,A5)') "NatCurves/",xcalledhere,"W.dat"
           iflag=1
        endif
     endif

     korder=min(numcurves,numr)

     if (iflag==1) then
        open(888,file=filename(1:getlen(filename)), status="unknown")
        open(889,file=filename2(1:getlen(filename2)), status="unknown")
        open(890,file="NatCurves/Points.dat", status="unknown")
        open(891,file="NatCurves/Params.dat", status="unknown")
        write(891,*) korder, numr, xcalledhere
        close(891) 
        do isplit=1,numr
           do i=1,korder
                 call sparseconfigmultone(natconfigs(:,i),natderiv,yyy%cptr(0), yyy%sptr(0), 1,pulseflag,0, isplit ,thistime)

              curves(isplit,i)=dot(natconfigs(:,i),natderiv,numconfig)/dot(natconfigs(:,i),natconfigs(:,i), &
                   numconfig) !! ok implicit.

           enddo

           write(888,'(100F12.6)') real(curves(isplit,1:korder))
!           OFLWR "program quick plotting routine in PROJECT"; CFLST
           write(889,'(100F12.6)') (abs(natproj(isplit,1:korder,imc)**2/bondweights(isplit)),imc=1,mcscfnum)
           write(890,'(100F12.6)') real(bondpoints(isplit),8)
        enddo
        write(888,*);        write(889,*);        write(890,*)
        close(888);        close(889);        close(890)
     endif
  enddo
  
  if (xcalledhere==1) then
     mmax=-1.d+10;     mmin=1.d+10
     do i=1,numr
        if (real(curves(i,1)).lt.mmin) then
           mmin=real(curves(i,1))
        endif
        if (real(curves(i,1)).gt.mmax) then
           mmax=real(curves(i,1))
        endif
     enddo
     mmax=mmin+1.d0;     mmin=mmin-0.2d0
  endif
  
#ifndef NOMATHEMATICA

  open(809, file="MC-LinesQQQ.txt", status="old", iostat=iii)
  if ( iii /=0 ) then
     if (xcalledhere.lt.2) then
        OFLWR
        write(mpifileptr,*) "  **  I did not find MC-LinesQQQ.txt in the directory, will not make MC-Nat.txt"
        WRFL; CFL
     endif
  else
     close(809)
     call system("cp MC-LinesQQQ.txt MC-Nat.txt")
     call system('perl -pi -e "s/KKK/Nat/" MC-Nat.txt')

     cwdline=nnbb
     write(cwdline,*) 'perl -pi -e "s/MMM/'
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A33)') "H(t) Natural Configuration Curves"
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A13)') '/" MC-Nat.txt'

     call system(cwdline)

     cwdline=nnbb
     write(cwdline,*) 'perl -pi -e "s/LLL/'

     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A31)') "H0 Natural Configuration Curves"

     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A13)') '/" MC-Nat.txt'

     call system(cwdline)

     cwdline=nnbb
     write(cwdline,*) 'perl -pi -e "s/RRR/'
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A17)') "MathNatCurves.gif"
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A13)') '/" MC-Nat.txt'

     call system(cwdline)

     cwdline=nnbb
     write(cwdline,*) 'perl -pi -e "s/RXR/'
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A18)') "MathNatCurves0.gif"
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A13)') '/" MC-Nat.txt'

     call system(cwdline)

     do nn=1,2
        cwdline=nnbb
        write(cwdline,*) 'perl -pi -e "s/XYZ/'
        ilen=getlen2(cwdline)+1
        write(cwdline(ilen:),'(F12.6)') mmax
        ilen=getlen2(cwdline)+1
        write(cwdline(ilen:),'(A13)') '/" MC-Nat.txt'
        
     call system(cwdline)

     cwdline=nnbb
     write(cwdline,*) 'perl -pi -e "s/ZYX/'
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(F12.6)') mmin
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A13)') '/" MC-Nat.txt'

     call system(cwdline)

     cwdline=nnbb
     write(cwdline,*) 'perl -pi -e "s/AAA/'
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(F9.6)') 999d0  !!! T E M P   !!!real(bondpoints(1),8)
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A13)') '/" MC-Nat.txt'

     call system(cwdline)

     cwdline=nnbb
     write(cwdline,*) 'perl -pi -e "s/BBB/'
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(F9.6)') 999d0   !!! T E M P   !!!real(bondpoints(numr),8)
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A13)') '/" MC-Nat.txt'

     call system(cwdline)

  enddo

!  cwdline=nnbb
!  write(cwdline,*) 'perl -pi -e "s/SSS/'
!  ilen=getlen2(cwdline)+2
!  write(cwdline(ilen:),'(I2)') iprop
!  ilen=getlen2(cwdline)+2
!  write(cwdline(ilen:),*) '/" MC-Nat.txt'
!  call system(cwdline)
  
  cwdline=nnbb
  write(cwdline,*) 'perl -pi -e "s/TTT/'
  ilen=getlen2(cwdline)+2
  if (numcurves.lt.10) then
     write(cwdline(ilen:),'(I1)') numcurves
  else 
     write(cwdline(ilen:),'(I2)') numcurves
  endif
  
  ilen=getlen2(cwdline)+2
  write(cwdline(ilen:),*) '/" MC-Nat.txt'
  call system(cwdline)
  
  cwdline=nnbb
  write(cwdline,*) 'perl -pi -e "s/UUU/'
  ilen=getlen2(cwdline)+2
  if (numr.lt.10) then
     write(cwdline(ilen:),'(I1)') numr
  else if (numr.lt.100) then
     write(cwdline(ilen:),'(I2)') numr
  else if (numr.lt.1000) then
     write(cwdline(ilen:),'(I3)') numr
  else
     write(cwdline(ilen:),'(I4)') numr
  endif
  
  ilen=getlen2(cwdline)+2
  write(cwdline(ilen:),*) '/" MC-Nat.txt'
  call system(cwdline)
  
  cwdline=nnbb
  write(cwdline,*) 'perl -pi -e "s/VVV/'
  ilen=getlen2(cwdline)+2
  if (xcalledhere.lt.10) then
     write(cwdline(ilen:),'(I1)') xcalledhere
  else if (xcalledhere.lt.100) then
     write(cwdline(ilen:),'(I2)') xcalledhere
  else if (xcalledhere.lt.1000) then
     write(cwdline(ilen:),'(I3)') xcalledhere
  else if (xcalledhere.lt.10000) then
     write(cwdline(ilen:),'(I4)') xcalledhere
  else
     write(cwdline(ilen:),'(I5)') xcalledhere
  endif
  
  ilen=getlen2(cwdline)+2
  write(cwdline(ilen:),*) '/" MC-Nat.txt'
  call system(cwdline)
  
  cwdline=nnbb
  write(cwdline,*) 'perl -pi -e "s/XXX/'
  ilen=getlen2(cwdline)+2
  
  rsum=plotmodulus*par_timestep*0.02418884d0
  
  if (rsum.lt.10.d0) then
     write(cwdline(ilen:),'(F7.5)') rsum
  else if (rsum.lt.100.d0) then
     write(cwdline(ilen:),'(F8.5)') rsum
  else if (rsum.lt.1000.d0) then
     write(cwdline(ilen:),'(F9.5)') rsum
  else if (rsum.lt.10000.d0) then
     write(cwdline(ilen:),'(F10.5)') rsum
  else if (rsum.lt.100000.d0) then
     write(cwdline(ilen:),'(F11.5)') rsum
  else if (rsum.lt.1000000.d0) then
     write(cwdline(ilen:),'(F12.5)') rsum
  else 
     write(cwdline(ilen:),'(F13.5)') rsum
  endif
  
  ilen=getlen2(cwdline)+2
  write(cwdline(ilen:),*) '/" MC-Nat.txt'
  
  call system(cwdline)

!     if (mod(xcalledhere,10).eq.1) then
!        call openfile()
!        write(mpifileptr,*) "not Calling Mathematica -- do it yourself: "
!        write(mpifileptr,*) '     math -noprompt -run "<<MC-Nat.txt"'
!        call closefile()
!!        call system('math -noprompt -run "<<MC-Nat.txt"')
!     endif
!     write(sysline,'(A24,I1)') "mv MC-Nat.txt MC-Nat.txt",iprop
!     call system(sysline(1:25))

  endif

#endif

  tempdenmat(:,:,:)=0d0
  do i=1,min(numr,min(numconfig,4))
     call getdenmat0(natconfigs(:,i),tempdenmat(:,:,i),tempinvden(:,:),tempdenvals(:),tempdenvects(:,:),1)
  enddo
  call save_denproj( 4, thistime, yyy%cmfpsivec(spfstart,0),  tempdenmat, denprojplotbin)
  
end subroutine save_natproj




subroutine save_denproj_final()
  use parameters
  implicit none
  call close_orbvector(denprojfile)
end subroutine save_denproj_final

subroutine save_denproj_initial(denfilename)
  use parameters
  implicit none
  character :: denfilename*(*)
  open(denprojfile,file=denfilename, status="replace", form="unformatted");  close(denprojfile)
end subroutine save_denproj_initial

subroutine save_denproj( nproj, thistime, inspfs, indenmats, denfilename)
  use parameters
  use natprojmod
  use mpimod
  implicit none

  DATATYPE :: inspfs(spfdims(1),spfdims(2),spfdims(3),nspf), indenmats(nspf,nspf,nproj), csum
  character :: denfilename*(*)
  integer :: i,j, imval, iproj, nproj,jj
  real*8 :: thistime
  character (len=headersize) :: header

  complex*16 ::        mtrans(spfdims(3),spfdims(3))
  DATATYPE :: mdensity(spfdims(1),spfdims(2),spfdims(3))
  complex*16 ::     cmdensity(spfdims(1),spfdims(2),spfdims(3)), density(spfdims(1),spfdims(2),spfdims(3))

  if (spfdimtype(3).eq.1) then  !! assume fourier basis (-mbig:mbig)
     if (mod(spfdims(3),2).ne.1) then
        OFLWR "FOURRR ERROR"; CFLST
     endif
     do i=1,spfdims(3)
        do j=1,spfdims(3)
           jj=j-(spfdims(3)+1)/2;           mtrans(i,j) = exp((0.d0,1.d0)*jj*2*pi*i/real(spfdims(3)))
        enddo
     enddo
  else
     mtrans(:,:)=0d0
     do i=1,spfdims(3)
        mtrans(i,i)=1d0
     enddo
  endif

  do iproj=1,nproj
     call getdensity(density, indenmats(:,:,iproj), inspfs,nspf)
     call zgemm('N', 'N', spfdims(1)*spfdims(2),spfdims(3),spfdims(3), (1.d0,0.d0), density, spfdims(1)*spfdims(2), mtrans,spfdims(3), (0.d0,0.d0), cmdensity, spfdims(1)*spfdims(2))

     do imval=1,spfdims(3)
        mdensity(:,:,imval)=  cmdensity(:,:,imval) / &
             sqrt(elecweights(:,:,imval))   !! ok imp conv mctdh
     enddo
     csum=natvals(iproj)
     call write_nat_header(header,thistime,0,iproj,csum)
     call save_orbvector(mdensity, spfsize, denprojfile, denfilename, header)
  end do

end subroutine save_denproj



