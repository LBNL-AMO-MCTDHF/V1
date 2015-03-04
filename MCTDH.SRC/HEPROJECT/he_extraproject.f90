

!!!  EXTRA COORINATE SPECIFIC ROUTINES INCLUDING DIATOMSTATS FOR EXPECTATION VALUES (DiatomStats.Dat)


#include "Definitions.INC"

subroutine getmyparams(inmpifileptr,inpfile,spfdims,spfdimtype,reducedpotsize,outnumr,nucrepulsion,nonuc_checkflag)   !! inpfile input, others output dummy variables
  use myparams
  implicit none

  integer :: spfdims(3),spfdimtype(3),inmpifileptr,nonuc_checkflag,myiostat, reducedpotsize, outnumr,&
       nargs, getlen, i, len
  real*8 :: nucrepulsion
  character (len=200) :: buffer
  character (len=200) :: nullbuff="                                                                                "
  character :: inpfile*(*)

  NAMELIST /heparinp/  &
       henumpoints,  henumelements,  hecelement,  heecstheta,   heelementsizes, numhatoms, hlocs,hlocrealflag,&
       hlocreal,lbig,nuccharge1,mbig, temp_glflag, num_skip_orbs, orb_skip_mvalue, orb_skip,debugflag
#ifdef PGFFLAG
  integer :: myiargc
  nargs=myiargc()
#else
  nargs=iargc()
#endif
  mpifileptr=inmpifileptr

  open(971,file=inpfile, status="old", iostat=myiostat)

  if (myiostat/=0) then

     OFLWR "No Input.Inp found, not reading heparams. iostat=",myiostat;CFL
  else
     read(971,nml=heparinp)
  endif
  close(971)

  call openfile()
  do i=1,nargs
     buffer=nullbuff
     call getarg(i,buffer)
     len=getlen(buffer)
     if (buffer(1:6) .eq. 'Numel=') then
        read(buffer(7:len),*) henumelements
        write(mpifileptr, *) "Radial electronic numelements set to ", henumelements
     endif
     if (buffer(1:7) .eq. 'Numpts=') then
        read(buffer(8:len),*) henumpoints
        write(mpifileptr, *) "Radial electronic numpoints set to ", henumpoints
     endif
     if (buffer(1:9) .eq. 'Celement=') then
        read(buffer(10:len),*) hecelement
        write(mpifileptr, *) "Radial electronic celement set to ", hecelement
     endif
  enddo
  call closefile()

  hegridpoints = henumelements*henumpoints-henumelements+1   
  bandwidth=2*henumpoints

  jacobisummax=lbig
  lseriesmax=jacobisummax  

  numerad = henumelements*henumpoints-henumelements-1
  nonuc_checkflag=1
  mseriesmax=0
  edim=numerad*(lbig+1)

  spfdims(1)=numerad;  spfdims(2)=lbig+1;  spfdims(3)=2*mbig+1

  spfdimtype(1)=0;  spfdimtype(2)=2;  spfdimtype(3)=1

!!!       xxx%reducedpot(numerad,lbig+1,-2*mbig:2*mbig,  nspf,nspf, 0:numreduced), &
  reducedpotsize=numerad*(lbig+1)*(4*mbig+1)

  outnumr=numr   !! which is 1
  nucrepulsion=0d0

end subroutine getmyparams

subroutine printmyopts()
  use myparams !! mpifileptr
  implicit none
  write(mpifileptr,'(A40,3I5)') "Henumpoints, henumelements:", henumpoints,  henumelements
  write(mpifileptr,'(A40,3I5)') "  Lbig,mbig                ",lbig,mbig
#ifdef ECSFLAG
  write(mpifileptr,'(A40)') "ECS is ON!"
  write(mpifileptr,'(A40,I5, F18.12)')   "hecelement,  heecstheta",    hecelement,  heecstheta
#else
  write(mpifileptr,'(A40)') "ECS is OFF!"
#endif
  write(mpifileptr,'(A20,100F10.5)') "He elementsizes:",  heelementsizes(1:henumelements)
end subroutine printmyopts

subroutine checkmyopts()
end subroutine checkmyopts


!!  PASS ARGUMENTS AS REAL (e.g. for xi) !!

function  radiallobatto(n,x, mvalue)
  use myparams
  use myprojectmod
  integer :: mvalue,   n
  real*8 :: x
  DATAECS :: xilobatto, radiallobatto
  radiallobatto=xilobatto(n,x,mvalue,henumpoints,glpoints2d(:,1),glpoints2d(:,2),glpoints,glweights,heelementsizes,henumelements,0)
end function radiallobatto

function  angularlobatto(n,x, mvalue)
  use myparams
  use myprojectmod
  implicit none
  integer :: mvalue,   n
  real*8 :: x,angularlobatto, etalobatto
     angularlobatto=etalobatto(n,x,mvalue,jacobipoints,jacobiweights)
end function angularlobatto

function cylindricalvalue(radpoint, thetapoint,nullrvalue,mvalue, invector)
  use myparams
  use myprojectmod
  implicit none
  DATATYPE ::  cylindricalvalue, sum,invector(numerad,lbig+1)
  real*8 :: radpoint,thetapoint,nullrvalue, angularlobatto
  DATAECS :: radiallobatto
  integer :: mvalue, ixi,lvalue
  sum=0.d0
  do ixi=1,numerad
     do lvalue=1,lbig+1
        sum=sum + &
             radiallobatto(ixi,radpoint, mvalue) * &
             angularlobatto(lvalue,thetapoint, mvalue) * invector(ixi,lvalue) * &
             1.d0/sqrt(2.d0*3.14159265d0)
     enddo
  enddo
  cylindricalvalue=sum
end function cylindricalvalue

subroutine get_maxsparse(nx,ny,nz,xvals,yvals,zvals, maxsparse,povsparse)
   use myparams
   use myprojectmod
   implicit none

   integer :: nx,ny,nz, maxsparse,  mvalue, ixi,lvalue, ix,iy,iz, iii, jjj
   real*8 :: xvals(nx),yvals(ny),zvals(nz),povsparse, xval,yval,zval, angularlobatto, costhetaval,phival, rhoval
   DATAECS :: radiallobatto
   complex*16 :: csum

   maxsparse=0;   jjj=0
   do mvalue=-mbig,mbig
   do lvalue=1,lbig+1
   do ixi=1,numerad

   iii=0
   do iz=1,nz
   do iy=1,ny
   do ix=1,nx

      xval=xvals(ix);      yval=yvals(iy);      zval=zvals(iz)
      phival=atan2(xval,yval);   rhoval=sqrt(xval**2+yval**2+zval**2)
      costhetaval=zval/rhoval

      if (.not.(costhetaval.gt.-.9999999d0)) then
         costhetaval=-0.9999999d0
      endif
      if (.not.(costhetaval.lt.0.9999999d0)) then
         costhetaval=0.9999999d0
      endif

      OFLWR "CHECK RADIALLOBATTO"; CFLST
      
      csum=          radiallobatto(ixi,rhoval, mvalue) * &
           angularlobatto(lvalue,costhetaval, mvalue) * 1.d0 / &
           sqrt(glpoints(ixi)**2 - jacobipoints(lvalue)**2) * &
           1.d0/sqrt(2.d0*3.14159265d0) * exp((0.d0,1.d0)*mvalue*phival) 
      if (abs(csum).gt.povsparse) then
         iii=iii+1
      endif
   enddo
   enddo
   enddo

   jjj=jjj+iii

  enddo
  enddo
  enddo

  maxsparse=jjj

end subroutine get_maxsparse



subroutine get_sphericalsparse(nx,ny,nz,xvals,yvals,zvals, maxsparse,sparsetransmat,sparsestart,sparseend,sparsexyz,povsparse)
   use myparams
   use myprojectmod
   implicit none
   integer :: nx,ny,nz, maxsparse, sparsexyz(maxsparse,3), iflag
   integer :: mvalue, ixi,lvalue, ix,iy,iz, iii, sparsestart(numerad,lbig+1,-mbig:mbig), &
        sparseend(numerad,lbig+1,-mbig:mbig)
   real*8 :: xvals(nx),yvals(ny),zvals(nz),povsparse, xval,yval,zval, angularlobatto, costhetaval,phival, rhoval
   DATAECS :: radiallobatto
   complex*16 :: sparsetransmat(maxsparse),  csum

   iii=0;   sparsestart=-1;   sparseend=-100
   do mvalue=-mbig,mbig
   do lvalue=1,lbig+1
   do ixi=1,numerad

      iflag=0

   do iz=1,nz
   do iy=1,ny
   do ix=1,nx

      xval=xvals(ix);      yval=yvals(iy);      zval=zvals(iz)
      phival=atan2(xval,yval);   rhoval=sqrt(xval**2+yval**2+zval**2)
      costhetaval=zval/rhoval
      if (.not.(costhetaval.gt.-.9999999d0)) then
         costhetaval=-0.9999999d0
      endif
      if (.not.(costhetaval.lt.0.9999999d0)) then
         costhetaval=0.9999999d0
      endif
      
OFLWR "CHECK RADIALLOBATTO"; CFLST

      csum=         radiallobatto(ixi,rhoval, mvalue) * &
           angularlobatto(lvalue,costhetaval, mvalue) * 1.d0 / &
           sqrt(glpoints(ixi)**2 - jacobipoints(lvalue)**2) * &
           1.d0/sqrt(2.d0*3.14159265d0) * exp((0.d0,1.d0)*mvalue*phival) 
   
      if (abs(csum).gt.povsparse) then
         iii=iii+1
         if (iflag==0) then
            sparsestart(ixi,lvalue,mvalue)=iii
            iflag=1
         endif
         sparseend(ixi,lvalue,mvalue)=iii
         if (iii.gt.maxsparse) then
            OFLWR"Sparse error!";         CFLST
            stop
      endif
      sparsexyz(iii,1) = ix;      sparsexyz(iii,2) = iy;      sparsexyz(iii,3) = iz;      sparsetransmat(iii) = csum
   endif
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo

end subroutine get_sphericalsparse


function interpolate(radpoint, thetapoint,nullrvalue,mvalue, ixi,lvalue)
  use myparams
  use myprojectmod
  implicit none
  DATATYPE ::  interpolate, sum
  real*8 :: radpoint,thetapoint, angularlobattoint,nullrvalue
  DATAECS :: radiallobattoint
  integer :: mvalue, ixi,lvalue
  sum=0.d0
  sum=sum +        radiallobattoint(ixi,radpoint, mvalue) *       angularlobattoint(lvalue,thetapoint, mvalue) 
  interpolate=sum
end function interpolate

!!  PASS ARGUMENTS AS REAL (e.g. for xi) !!

function  radiallobattoint(n,x, mvalue)
  use myparams
  use myprojectmod
  implicit none
  integer :: mvalue,   n
  real*8 :: x
  DATAECS :: xilobattoint, radiallobattoint
     radiallobattoint=xilobattoint(n,x,mvalue,henumpoints,glpoints2d(:,1),glpoints2d(:,2),glpoints,glweights,heelementsizes,henumelements,0)
end function radiallobattoint

function  angularlobattoint(n,x, mvalue)
  use myparams
  use myprojectmod
  implicit none
  integer :: mvalue,   n
  real*8 :: x,angularlobattoint, etalobattoint
     angularlobattoint=etalobattoint(n,x,mvalue,jacobipoints,jacobiweights)
end function angularlobattoint



!! rvalue not used for atom
function radiusvalue(spfindex,notused)
  use myparams
  use myprojectmod
  implicit none
  integer :: spfindex,ir
  DATATYPE :: radiusvalue, notused

  ir=mod(spfindex-1,numerad)+1

  radiusvalue=glpoints(ir+1)

end function radiusvalue
