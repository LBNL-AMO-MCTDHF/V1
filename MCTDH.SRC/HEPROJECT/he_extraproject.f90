

!!!  EXTRA COORINATE SPECIFIC ROUTINES INCLUDING DIATOMSTATS FOR EXPECTATION VALUES (DiatomStats.Dat)


#include "Definitions.INC"

!! inpfile input, others output dummy variables

subroutine getmyparams(inmpifileptr,inpfile,spfdims,spfdimtype,reducedpotsize,outnumr,nucrepulsion,&
     nonuc_checkflag)
  use myparams
  implicit none
  integer,intent(in) :: inmpifileptr
  character,intent(in) :: inpfile*(*)
  integer,intent(out) :: spfdims(3),spfdimtype(3),nonuc_checkflag,reducedpotsize, outnumr
  real*8,intent(out) :: nucrepulsion
  character (len=SLN) :: buffer
  character (len=SLN) :: nullbuff
  integer ::        nargs, getlen, i, len,myiostat

  NAMELIST /heparinp/  &
       henumpoints,  henumelements,  hecelement,  heecstheta,   heelementsizes, numhatoms, hlocs,hlocrealflag,&
       hlocreal,lbig,nuccharge1,mbig, temp_glflag, num_skip_orbs, orb_skip_mvalue, orb_skip,debugflag,&
       loadedocc,ivoflag

#ifdef PGFFLAG
  integer :: myiargc
  nargs=myiargc()
#else
  nargs=iargc()
#endif
  mpifileptr=inmpifileptr

  do i=1,SLN
     nullbuff(i:i)=" "
  enddo

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


module dvrvalmod
contains

!!  PASS ARGUMENTS AS REAL (e.g. for xi) !!

function  radiallobatto(n,x, mvalue)
  use myparams
  use myprojectmod
  use lobstuffmod
  implicit none
  integer,intent(in) :: mvalue,   n
  real*8,intent(in) :: x
  DATAECS :: radiallobatto
  radiallobatto=xilobatto(n,x,mvalue,henumpoints,glpoints2d(:,1),glpoints2d(:,2),&
       glpoints,glweights,heelementsizes,henumelements,0)
end function radiallobatto

function  angularlobatto(n,x, mvalue)
  use myparams
  use myprojectmod
  use lobstuffmod
  implicit none
  integer,intent(in) :: mvalue,   n
  real*8,intent(in) :: x
  real*8 :: angularlobatto
     angularlobatto=etalobatto(n,x,mvalue,jacobipoints,jacobiweights)
end function angularlobatto

!!  PASS ARGUMENTS AS REAL (e.g. for xi) !!

function  radiallobattoint(n,x, mvalue)
  use myparams
  use myprojectmod
  use lobstuffmod
  implicit none
  integer,intent(in) :: mvalue, n
  real*8,intent(in) :: x
  DATAECS :: radiallobattoint
  radiallobattoint=xilobattoint(n,x,mvalue,henumpoints,glpoints2d(:,1),glpoints2d(:,2),&
       glpoints,glweights,heelementsizes,henumelements,0)
end function radiallobattoint

function  angularlobattoint(n,x, mvalue)
  use myparams
  use myprojectmod
  use lobstuffmod
  implicit none
  integer,intent(in) :: mvalue,   n
  real*8,intent(in) :: x
  real*8 :: angularlobattoint
     angularlobattoint=etalobattoint(n,x,mvalue,jacobipoints,jacobiweights)
end function angularlobattoint

end module dvrvalmod


function cylindricalvalue(radpoint, thetapoint,nullrvalue,mvalue, invector)
  use myparams
  use myprojectmod
  use dvrvalmod
  implicit none
  integer,intent(in) :: mvalue
  DATATYPE,intent(in) :: invector(numerad,lbig+1)
  real*8,intent(in) :: radpoint, thetapoint, nullrvalue
  DATATYPE ::  cylindricalvalue, sum
  integer :: ixi,lvalue

  sum=0.d0*nullrvalue !! avoid warn unused
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
   use dvrvalmod
   implicit none
   integer,intent(in) :: nx,ny,nz
   integer,intent(out) :: maxsparse
   real*8, intent(in) :: xvals(nx),yvals(ny),zvals(nz), povsparse
   integer :: mvalue, ixi,lvalue, ix,iy,iz, iii, jjj
   real*8 :: xval,yval,zval, costhetaval,phival, rhoval
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


subroutine get_sphericalsparse(nx,ny,nz,xvals,yvals,zvals, maxsparse,sparsetransmat,&
     sparsestart,sparseend,sparsexyz,povsparse)
   use myparams
   use myprojectmod
   use dvrvalmod
   implicit none
   integer,intent(in) :: nx,ny,nz,maxsparse
   real*8,intent(in) :: xvals(nx),yvals(ny),zvals(nz),povsparse
   complex*16,intent(out) :: sparsetransmat(maxsparse)
   integer,intent(out) :: sparsexyz(maxsparse,3),  sparsestart(numerad,lbig+1,-mbig:mbig), &
        sparseend(numerad,lbig+1,-mbig:mbig)
   real*8 :: xval,yval,zval, costhetaval,phival, rhoval
   complex*16 :: csum
   integer :: mvalue, ixi,lvalue, ix,iy,iz, iii, iflag

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
      sparsexyz(iii,1) = ix;      sparsexyz(iii,2) = iy;     
      sparsexyz(iii,3) = iz;      sparsetransmat(iii) = csum
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
  use dvrvalmod
  implicit none
  integer,intent(in) :: mvalue, ixi,lvalue
  real*8,intent(in) :: radpoint,thetapoint
  DATATYPE ::  interpolate, sum
  real*8 :: nullrvalue

  sum=0.d0* nullrvalue !! avoid warn unused
  sum=sum +    radiallobattoint(ixi,radpoint, mvalue) * angularlobattoint(lvalue,thetapoint, mvalue)
  interpolate=sum
end function interpolate



!!$  !! rvalue not used for atom
!!$  function radiusvalue(spfindex,notused)
!!$    use myparams
!!$    use myprojectmod
!!$    implicit none
!!$    integer,intent(in) :: spfindex
!!$    integer :: ir
!!$    DATATYPE :: radiusvalue, notused,notused2
!!$    notused2=notused*0 !! avoid warn unused
!!$    ir=mod(spfindex-1,numerad)+1
!!$    radiusvalue=glpoints(ir+1)
!!$  end function radiusvalue
