

!!!  EXTRA COORINATE SPECIFIC ROUTINES INCLUDING DIATOMSTATS FOR EXPECTATION VALUES (DiatomStats.Dat)


#include "Definitions.INC"

subroutine getmyparams(inmpifileptr,inpfile,spfdims,spfdimtype,reducedpotsize,&
     outnumr,nucrepulsion,nonuc_checkflag)   !! inpfile input, others output dummy variables
  use myparams
  implicit none
  integer,intent(in) :: inmpifileptr
  character,intent(in) :: inpfile*(*)
  integer,intent(out) :: spfdims(3),spfdimtype(3), reducedpotsize, outnumr, nonuc_checkflag
  real*8,intent(out) :: nucrepulsion
  integer :: nargs, getlen, len, myiostat, i
  character (len=SLN) :: buffer
  character (len=SLN) :: nullbuff
#ifdef PGFFLAG
  integer ::    myiargc
#endif
  NAMELIST /h2parinp/  &
       JVALUE,    rcelement, rthetaecs,        xinumpoints,  xinumelements,  xicelement,  xiecstheta,  &
       pro_Hmass, pro_Dmass,  rnumelements,  rnumpoints,  numr, &
       relementsize ,  rstart ,xielementsizes, xiecstheta, bornopflag, capflag, capstrength, cappower, &
       numhatoms, hlocs,hlocrealflag,hlocreal,mbig,lbig,nuccharge1,nuccharge2, &
       num_skip_orbs, orb_skip_mvalue, orb_skip,twoeattractflag,debugflag,ivoflag,loadedocc

  do i=1,SLN
     nullbuff(i:i)=" "
  enddo

#ifdef PGFFLAG
  nargs=myiargc()
#else
  nargs=iargc()
#endif
  mpifileptr=inmpifileptr

  open(971,file=inpfile, status="old", iostat=myiostat)
  if (myiostat/=0) then
     OFLWR "No Input.Inp found, not reading h2params. iostat=",myiostat;CFL
  else
     read(971,nml=h2parinp)
  endif
  close(971)

  call openfile()
  do i=1,nargs
     buffer=nullbuff;     call getarg(i,buffer);     len=getlen(buffer)
     if (buffer(1:6) .eq. 'Numel=') then
        read(buffer(7:len),*) xinumelements
        write(mpifileptr, *) "Radial electronic numelements set to ", xinumelements
     endif
     if (buffer(1:7) .eq. 'RNumel=') then
        read(buffer(8:len),*) rnumelements
        write(mpifileptr, *) "Bond distance (R) numelements set to ", rnumelements
     endif

     if (buffer(1:8) .eq. 'RNumpts=') then
        read(buffer(9:len),*) rnumpoints
        write(mpifileptr, *) "Bond distance (R) numpoints set to ", rnumpoints
     endif
     if (buffer(1:7) .eq. 'Numpts=') then
        read(buffer(8:len),*) xinumpoints
        write(mpifileptr, *) "Radial electronic numpoints set to ", xinumpoints
     endif
     if (buffer(1:6) .eq. 'Theta=') then
        read(buffer(7:len),*) xiecstheta
        write(mpifileptr, *) "Scaling angle for xi set to ", xiecstheta," by command line option."
     endif
     if (buffer(1:9) .eq. 'Celement=') then
        read(buffer(10:len),*) xicelement
        write(mpifileptr, *) "Radial electronic celement set to ", xicelement
     endif
  enddo
  call closefile()

  Rmass = pro_Hmass*pro_Dmass/(pro_Hmass+pro_Dmass)
  totalmass=pro_Hmass+pro_Dmass

!! nuccharge1 goes with pro_Hmass, nuccharge2 with pro_Dmass.
!! nucleus 1 (H) is on top (positive z).
!! if nucleus 1 (H) is heavier, mass_asym is negative; if it is
!! heavier, it is closer to the origin.  The result
!! is that z is evaluated to be near zero when xi=1,eta=1.

  mass_asym = (pro_Dmass - pro_Hmass) / (pro_Dmass+pro_Hmass)

  xigridpoints = xinumelements*xinumpoints-xinumelements+1   
  bandwidth=2*xinumpoints
  rgridpoints = (rnumpoints-1)*rnumelements+1
  numr=rgridpoints-2

  jacobisummax=lbig
  lseriesmax=jacobisummax  !!lbig+mbig+4
  numerad = xinumelements*xinumpoints-xinumelements
  numeta =  lbig + 1
  mseriesmax=mbig*2+4
  edim=(lbig+1)*numerad
  spfdims(1)=numerad;  spfdims(2)=lbig+1;  spfdims(3)=2*mbig+1
  spfdimtype(1)=0;  spfdimtype(2)=2;  spfdimtype(3)=1

!!!       xxx%reducedpot(numerad,lbig+1,-2*mbig:2*mbig,  nspf,nspf, 0:numreduced), &
  reducedpotsize=numerad*(lbig+1)*(4*mbig+1)

  outnumr=numr
  nucrepulsion=nuccharge1*nuccharge2

  if (numhatoms.gt.1) then
     do i=1,numhatoms
        if ( (hlocs(1,i).lt.0).or.(hlocs(2,i).lt.0).or.(hlocs(1,i).gt.numerad).or.&
             (hlocs(2,i).gt.lbig+1) .or. (abs(hlocs(3,i)).ne.1) ) then
           OFLWR "Error, numhatoms set to ", numhatoms
           write(mpifileptr,*) "  but hlocs assigned wrong for atom ", i
           write(mpifileptr,*) "  ", hlocs(:,i);CFLST
        endif
     enddo
  endif
  if (capflag.eq.1) then
     rthetaecs=0.d0
     write(mpifileptr,*) "CAPflag is on, setting rthetaecs to zero"
  endif
  if (bornopflag.eq.1) then
     reducedflag=0
  else
     reducedflag=1   
  endif
  if (capflag.eq.1) then
        capstart=rstart + (rcelement-1)*relementsize
  endif
  nonuc_checkflag=bornopflag

end subroutine getmyparams


!subroutine checkmyopts()
!  use myparams
!  implicit none
!  if ((nonuc_checkflag.eq.0) .and. (tdflag.eq.1) .and. (velflag/=0) .and. &
!        (abs (NucCharge1*pro_dmass - NucCharge2*pro_hmass) .gt. 1.d-3 ) ) then
!     call openfile()
!     write(mpifileptr,*) " Velocity with nuclear motion not available for hetero, only homo   - TEMP CONTINUE"
!     call closefile()
!  endif
!end subroutine checkmyopts



subroutine printmyopts()
  use myparams
  implicit none

  write(mpifileptr,*) "  *****  prolate spheroidal coordinate options ****"
  write(mpifileptr,'(A40,3I5)') "bornopflag:     ",  bornopflag
!  if (bornopflag==1) then
!     if (nonuc_checkflag==0) then
!        write(mpifileptr, *) "flag err 555"
!        call closefile()
!        call  mpistop()
!     else
!        write(mpifileptr,*) "    --> Born Oppenheimer calculation "
!     endif
!  else
!     if (nonuc_checkflag==0) then     
!        write(mpifileptr,*) "    --> Improved adiabatic calculation "
!     else
!        write(mpifileptr,*) "    --> NONADIABATIC calc with full nuclear KE "
!     endif
!  endif

  write(mpifileptr,'(A40,3I5)') "Xinumpoints, xinumelements:", xinumpoints,  xinumelements
#ifdef ECSFLAG
  write(mpifileptr,'(A40)') "ECS is ON!"
  write(mpifileptr,'(A40,I5, F18.12)')   "xicelement,  xiecstheta",    xicelement,  xiecstheta
#else
  write(mpifileptr,'(A40)') "ECS is OFF!"
#endif
  write(mpifileptr,'(A20,100F10.5)') "Xi elementsizes:",  xielementsizes(1:xinumelements)
  write(mpifileptr,'(A20,100I10.5)') "  lbig,mbig     ",lbig,mbig

  write(mpifileptr, *)
  write(mpifileptr, *) "**************************   Parameters: nuclear    *************************   "
  write(mpifileptr, *)
  write(mpifileptr,'(A25,3F18.8)') "Hmass, Dmass, Reduced:", pro_Hmass, pro_Dmass, Rmass
  write(mpifileptr,'(A25,2I5)') "rnumelements, rnumpoints:",  rnumelements,  rnumpoints
  write(mpifileptr,'(A25,2F10.5)') "relementsize ,  rstart:", relementsize ,  rstart
#ifdef ECSFLAG
  write(mpifileptr,'(A40)') "ECS is ON!"
  write(mpifileptr,'(A40,I5, F18.12)')   "rcelement,  rthetaecs",  rcelement,  rthetaecs
  if (capflag==1) then
     write(mpifileptr,'(A40)') "CAP is ON!"     
     write(mpifileptr,'(A40, F10.5)') "Capstrength: ", capstrength
     write(mpifileptr,'(A40, I5)') "Cappower:    ", cappower
     write(mpifileptr,'(A40, F10.5)') "Capstart (from ecs):    ", capstart
  else
     write(mpifileptr,'(A40)') "CAP is Off!"     
  endif
#else
  write(mpifileptr,'(A40)') "ECS is OFF!"
#endif
  write(mpifileptr,'(A25,3I5)') "Jvalue:", JVALUE
  write(mpifileptr,'(A40,3I5)') " Reducedflag:",  reducedflag

end subroutine printmyopts

module dvrvalmod
contains

function  radiallobatto(n,x, mvalue)
  use myprojectmod
  use myparams
  use lobstuffmod
  integer,intent(in) :: mvalue,   n
  real*8,intent(in) :: x
  DATAECS :: radiallobatto
  radiallobatto=xilobatto(n,x,mvalue,xinumpoints,xipoints2d(:,1), xipoints2d(:,2), xipoints,&
       xiweights,xielementsizes,xinumelements,1)
end function radiallobatto

function  angularlobatto(n,x, mvalue)
  use myparams
  use myprojectmod
  use lobstuffmod
  implicit none
  integer,intent(in) :: mvalue,   n
  real*8,intent(in) :: x
  real*8 :: angularlobatto
  angularlobatto=etalobatto(n,x,mvalue,etapoints,etaweights)
end function angularlobatto

!!  PASS ARGUMENTS AS REAL (e.g. for xi) !!

function  radiallobattoint(n,x, mvalue)
  use myprojectmod
  use myparams
  use lobstuffmod
  implicit none
  integer,intent(in) :: mvalue,   n
  real*8,intent(in) :: x
  DATAECS :: radiallobattoint
     radiallobattoint=xilobattoint(n,x,mvalue,xinumpoints,xipoints2d(:,1), &
          xipoints2d(:,2), xipoints,xiweights,xielementsizes,xinumelements,1)
end function radiallobattoint

function  angularlobattoint(n,x, mvalue)
  use myparams
  use myprojectmod
  use lobstuffmod
  implicit none
  integer,intent(in) :: mvalue,   n
  real*8,intent(in) :: x
  real*8 :: angularlobattoint
  angularlobattoint=etalobattoint(n,x,mvalue,etapoints,etaweights)
end function angularlobattoint

function xifunct(radpoint,thetapoint,rvalue)
  use myparams
  use lobstuffmod
  implicit none
  real*8,intent(in) :: radpoint,thetapoint,rvalue
  real*8 ::  sum1, sum2,xifunct

  sum1 = 0.25d0 + ( radpoint / rvalue )**2
  sum2 = radpoint / rvalue * thetapoint     !! thetapoint is costheta
  xifunct=sqrt(sum1+sum2) + sqrt(sum1-sum2) + 0.0000000001d0
  if (xifunct .le. 1.d0) then
     call openfile();     write(mpifileptr,*) "Xifunct err!"
     write(mpifileptr,*)  "    ", radpoint,thetapoint, rvalue
     call closefile();     call mpistop()
  endif
end function xifunct

!! returns eta/xi as function of spherical coords

function etafunct(radpoint,thetapoint,rvalue)
  use myparams
  use lobstuffmod
  implicit none
  real*8,intent(in) :: radpoint,thetapoint,rvalue
  real*8 :: sum1, sum2, etafunct

  sum1 = 0.25d0 + ( radpoint / rvalue )**2
  sum2 = radpoint / rvalue * thetapoint     !! thetapoint is costheta
  etafunct=sqrt(sum1+sum2) - sqrt(sum1-sum2) 
  if (abs(etafunct) .gt. 1.d0) then
     OFLWR "Etafunct err!", etafunct, sum1, sum2, radpoint, thetapoint, rvalue;CFLST
  endif
end function etafunct

end module dvrvalmod


function cylindricalvalue(radpoint, thetapoint,rvalue,mvalue, invector)
  use myparams
  use myprojectmod
  use dvrvalmod
  implicit none
  integer,intent(in) :: mvalue
  DATATYPE,intent(in) :: invector(numerad,lbig+1)
  real*8,intent(in) :: radpoint,thetapoint,rvalue
  DATATYPE ::  cylindricalvalue, sum
  integer :: ixi,lvalue
  sum=0.d0

  do ixi=1,numerad
     do lvalue=1,lbig+1
        sum=sum + &
             radiallobatto(ixi,xifunct(radpoint,thetapoint,rvalue), mvalue) * &
             angularlobatto(lvalue,etafunct(radpoint,thetapoint,rvalue), mvalue) * 1.d0 / &
             sqrt(xipoints(ixi)**2 - etapoints(lvalue)**2) * invector(ixi,lvalue) * &
             1.d0/sqrt(2.d0*3.14159265d0)
     enddo
  enddo
  cylindricalvalue=sum
end function cylindricalvalue


function sphericalvalue(xval,yval,zval,  inspf)
   use myprojectmod
   use myparams
   use dvrvalmod
   implicit none
   real*8,intent(in) :: xval,yval,zval
   DATATYPE,intent(in) ::  inspf(numerad,lbig+1,-mbig:mbig)
   complex*16 :: sphericalvalue, csum,csum2
   real*8 :: xival, etaval,phival, rhoval
   integer :: mvalue, ixi,lvalue

   phival=atan2(xval,yval);   rhoval=sqrt(xval**2+yval**2)
   etaval= (sqrt(rhoval**2 + (zval+1.0)**2) - sqrt(rhoval**2 + (zval-1.0)**2))/2.d0
   xival= (sqrt(rhoval**2 + (zval+1.0)**2) + sqrt(rhoval**2 + (zval-1.0)**2))/2.d0

   if (.not.(etaval.gt.-.9999999d0)) then
      etaval=-0.9999999d0
   endif
   if (.not.(etaval.lt.0.9999999d0)) then
      etaval=0.9999999d0
   endif
   if (.not.(xival.gt.1.0000001d0)) then
      xival=1.0000001d0
   endif
   csum=0.d0
  do mvalue=-mbig,mbig
   do ixi=1,numerad
   do lvalue=1,lbig+1
         csum2= &
              radiallobatto(ixi,xival, mvalue) * &
              angularlobatto(lvalue,etaval, mvalue) * 1.d0 / &
              sqrt(xipoints(ixi)**2 - etapoints(lvalue)**2) * inspf(ixi,lvalue,mvalue) * &
              1.d0/sqrt(2.d0*3.14159265d0) * exp((0.d0,1.d0)*mvalue*phival)
         csum=csum+csum2
  enddo
  enddo
  enddo
  sphericalvalue=csum

end function sphericalvalue


subroutine get_maxsparse(nx,ny,nz,xvals,yvals,zvals, maxsparse,povsparse)
   use myprojectmod
   use myparams
   use dvrvalmod
   implicit none
   integer,intent(in) :: nx,ny,nz
   integer,intent(out) :: maxsparse
   real*8,intent(in) :: xvals(nx),yvals(ny),zvals(nz),povsparse
   integer :: ixi,lvalue, ix,iy,iz, iii, jjj, mvalue
   real*8 :: xval,yval,zval, xival, etaval,phival, rhoval
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
      phival=atan2(xval,yval);   rhoval=sqrt(xval**2+yval**2)
      etaval= (sqrt(rhoval**2 + (zval+1.0)**2) - sqrt(rhoval**2 + (zval-1.0)**2))/2.d0
      xival= (sqrt(rhoval**2 + (zval+1.0)**2) + sqrt(rhoval**2 + (zval-1.0)**2))/2.d0

      if (.not.(etaval.gt.-.9999999d0)) then
         etaval=-0.9999999d0
      endif
      if (.not.(etaval.lt.0.9999999d0)) then
         etaval=0.9999999d0
      endif
      
      if (.not.(xival.gt.1.0000001d0)) then
         xival=1.0000001d0
      endif
      
      csum=          radiallobatto(ixi,xival, mvalue) * &
           angularlobatto(lvalue,etaval, mvalue) * 1.d0 / &
           sqrt(xipoints(ixi)**2 - etapoints(lvalue)**2) * &
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
   use myprojectmod
   use myparams
   use dvrvalmod
   implicit none
   integer,intent(in) :: nx,ny,nz, maxsparse
   integer,intent(out) :: sparsexyz(maxsparse,3),  sparsestart(numerad,lbig+1,-mbig:mbig), &
        sparseend(numerad,lbig+1,-mbig:mbig)
   real*8,intent(in) :: xvals(nx),yvals(ny),zvals(nz), povsparse
   complex*16,intent(out) :: sparsetransmat(maxsparse)
   real*8 :: xval,yval,zval, xival, etaval,phival, rhoval
   integer :: mvalue, ixi,lvalue, ix,iy,iz, iii, iflag
   complex*16 :: csum

   iii=0;   sparsestart=-1;   sparseend=-100

   do mvalue=-mbig,mbig
   do lvalue=1,lbig+1
   do ixi=1,numerad

      iflag=0

   do iz=1,nz
   do iy=1,ny
   do ix=1,nx

      xval=xvals(ix);      yval=yvals(iy);      zval=zvals(iz)
      phival=atan2(xval,yval);   rhoval=sqrt(xval**2+yval**2)
      etaval= (sqrt(rhoval**2 + (zval+1.0)**2) - sqrt(rhoval**2 + (zval-1.0)**2))/2.d0
      xival= (sqrt(rhoval**2 + (zval+1.0)**2) + sqrt(rhoval**2 + (zval-1.0)**2))/2.d0

      if (.not.(etaval.gt.-.9999999d0)) then
         etaval=-0.9999999d0
      endif
      if (.not.(etaval.lt.0.9999999d0)) then
         etaval=0.9999999d0
      endif
      if (.not.(xival.gt.1.0000001d0)) then
         xival=1.0000001d0
      endif
   
      csum=         radiallobatto(ixi,xival, mvalue) * &
           angularlobatto(lvalue,etaval, mvalue) * 1.d0 / &
           sqrt(xipoints(ixi)**2 - etapoints(lvalue)**2) * &
           1.d0/sqrt(2.d0*3.14159265d0) * exp((0.d0,1.d0)*mvalue*phival) 
   
      if (abs(csum).gt.povsparse) then
         iii=iii+1
         if (iflag==0) then
            sparsestart(ixi,lvalue,mvalue)=iii
            iflag=1
         endif
         sparseend(ixi,lvalue,mvalue)=iii
         if (iii.gt.maxsparse) then
            call openfile();         write(mpifileptr,*) "Sparse error!";         call closefile()
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


function interpolate(radpoint, thetapoint,rvalue,mvalue,ixi,lvalue)
  use myprojectmod
  use myparams
  use dvrvalmod
  implicit none
  integer,intent(in) :: mvalue, ixi,lvalue
  real*8,intent(in) :: radpoint,thetapoint,rvalue
  DATATYPE ::  interpolate, sum

  sum=0.d0
  sum=sum + &
       radiallobattoint(ixi,xifunct(radpoint,thetapoint,rvalue), mvalue) * &
       angularlobattoint(lvalue,etafunct(radpoint,thetapoint,rvalue), mvalue) 
  interpolate=sum
end function interpolate



!!$  function radiusvalue(spfindex,rvalue)
!!$    use myparams
!!$    use myprojectmod
!!$    implicit none
!!$    integer,intent(in) :: spfindex
!!$    DATAECS,intent(in) :: rvalue
!!$    integer :: ixi,ieta,jindex
!!$    DATATYPE :: radiusvalue
!!$    DATAECS :: radiusfun
!!$  
!!$    ixi=mod(spfindex-1,numerad)+1
!!$    jindex=(spfindex-ixi)/numerad+1
!!$  
!!$    if (jindex*numerad.ne.spfindex-ixi) then
!!$       OFLWR "ERRRR"; CFLST
!!$    endif
!!$  
!!$    ieta=mod(jindex-1,numeta)+1
!!$  
!!$    radiusvalue=radiusfun(xipoints(ixi),etapoints(ieta),rvalue)
!!$  
!!$  end function radiusvalue
!!$    
!!$  function radiusfun(xivalue,etavalue,rvalue)
!!$    implicit none
!!$    DATAECS,intent(in) :: xivalue,rvalue
!!$    real*8,intent(in) :: etavalue
!!$    DATAECS :: radiusfun
!!$  
!!$    radiusfun=sqrt(etavalue**2+xivalue**2-1)*rvalue
!!$  end function radiusfun
