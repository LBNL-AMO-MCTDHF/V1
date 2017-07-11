
#include "Definitions.INC"

module jfunctmod
contains

function ffunct(xval)
  use myparams
  implicit none
  real*8, intent(in) :: xval
  DATATYPE :: ffunct, fac
  real*8 :: xx

  fac=exp((0d0,1d0)*scalingtheta)*scalingstretch

  if (xval.le.(-1)*(scalingdistance+smoothness)) then
     ffunct=xval + (fac-1)*(xval+scalingdistance+smoothness/2d0)
  else if (xval.ge.scalingdistance+smoothness) then
     ffunct=xval + (fac-1)*(xval-scalingdistance-smoothness/2d0)
  else if (abs(xval).le.scalingdistance) then
     ffunct=xval
  else if (xval.ge.scalingdistance) then
     xx = (xval-scalingdistance)/smoothness
     ffunct = xval + (fac-1) * smoothness * fscaled(xx)
  else if (xval.le.(-1)*scalingdistance) then
     xx = (-xval-scalingdistance)/smoothness
     ffunct = xval - (fac-1) * smoothness * fscaled(xx)
  else
     print *, "OOGABLAH0"; stop
  endif
end function ffunct


function jfunct(xval)
  use myparams
  implicit none
  real*8, intent(in) :: xval
  DATATYPE :: jfunct, fac
  real*8 :: xx

  fac=exp((0d0,1d0)*scalingtheta)*scalingstretch

  if (xval.le.(-1)*(scalingdistance+smoothness)) then
     jfunct=fac
  else if (xval.ge.scalingdistance+smoothness) then
     jfunct=fac
  else if (abs(xval).le.scalingdistance) then
     jfunct=1d0
  else if (xval.ge.scalingdistance) then
     xx = (xval-scalingdistance)/smoothness
     jfunct = 1d0 + (fac-1d0) * jscaled(xx) 
  else if (xval.le.(-1)*scalingdistance) then
     xx = (-xval-scalingdistance)/smoothness
     jfunct = 1d0 + (fac-1d0) * jscaled(xx) 
  else
     print *, "OOGABLAHd"; stop
  endif
end function jfunct


function djfunct(xval)
  use myparams
  implicit none
  real*8, intent(in) :: xval
  DATATYPE :: djfunct, fac
  real*8 :: xx

  fac=exp((0d0,1d0)*scalingtheta)*scalingstretch

  if (xval.le.(-1)*(scalingdistance+smoothness)) then
     djfunct=0
  else if (xval.ge.scalingdistance+smoothness) then
     djfunct=0
  else if (abs(xval).le.scalingdistance) then
     djfunct=0
  else if (xval.ge.scalingdistance) then
     xx = (xval-scalingdistance)/smoothness
     djfunct = (fac-1d0) * djscaled(xx) / smoothness
  else if (xval.le.(-1)*scalingdistance) then
     xx = (-xval-scalingdistance)/smoothness
     djfunct = (1d0-fac) * djscaled(xx) / smoothness
  else
     print *, "OOGABLAHd"; stop
  endif
end function djfunct


function ddjfunct(xval)
  use myparams
  implicit none
  real*8, intent(in) :: xval
  DATATYPE :: ddjfunct, fac
  real*8 :: xx

  fac=exp((0d0,1d0)*scalingtheta)*scalingstretch

  if (xval.le.(-1)*(scalingdistance+smoothness)) then
     ddjfunct=0
  else if (xval.ge.scalingdistance+smoothness) then
     ddjfunct=0
  else if (abs(xval).le.scalingdistance) then
     ddjfunct=0
  else if (xval.ge.scalingdistance) then
     xx = (xval-scalingdistance)/smoothness
     ddjfunct = (fac-1d0) * ddjscaled(xx) / smoothness**2
  else if (xval.le.(-1)*scalingdistance) then
     xx = (-xval-scalingdistance)/smoothness
     ddjfunct = (fac-1d0) * ddjscaled(xx) / smoothness**2
  else
     print *, "OOGABLAHd"; stop
  endif
end function ddjfunct

!! f happens to go from 0 to 1/2 on (0,1)

function fscaled(xval)
  implicit none
  real*8 :: fscaled,xval,pixval
  real*8, parameter :: pi = 3.14159265358979323844d0 !!TEMP

  if (xval.lt.0d0.or.xval.gt.1d0) then
     print *, "SCALEDERR F",xval; stop
  endif
  pixval=pi*xval
  fscaled = (-1d0)/12d0/pi * sin(pixval)**3 - 1d0/2d0/pi * sin(pixval) + xval/2d0
end function fscaled

function jscaled(xval)
  implicit none
  real*8 :: jscaled,xval,pixval
  real*8, parameter :: pi = 3.14159265358979323844d0 !!TEMP

  if (xval.lt.0d0.or.xval.gt.1d0) then
     print *, "SCALEDERR F",xval; stop
  endif
  pixval=pi*xval
  jscaled = 1d0/4d0 * cos(pixval)**3 - 3d0/4d0 * cos(pixval) + 0.5d0
end function jscaled

function djscaled(xval)
  implicit none
  real*8 :: djscaled,xval,pixval
  real*8, parameter :: pi = 3.14159265358979323844d0 !!TEMP

  if (xval.lt.0d0.or.xval.gt.1d0) then
     print *, "SCALEDERR F",xval; stop
  endif
  pixval=pi*xval
  djscaled = 0.75d0 * pi * sin(pixval)**3 

end function djscaled

function ddjscaled(xval)
  implicit none
  real*8 :: ddjscaled,xval,pixval
  real*8, parameter :: pi = 3.14159265358979323844d0 !!TEMP

  if (xval.lt.0d0.or.xval.gt.1d0) then
     print *, "SCALEDERR F",xval; stop
  endif
  pixval=pi*xval
  ddjscaled = 9d0/4d0*pi**2 * sin(pixval)**2 * cos(pixval)
end function ddjscaled

end module jfunctmod


module ivopotmod
  implicit none
  integer :: numocc=(-1)
  DATATYPE, allocatable :: ivopot(:),ivo_occupied(:,:)
end module ivopotmod


module initspfsmod
contains

subroutine init_spfs(inspfs,numloaded,numfrozen,frozenreduced)
  use myparams
  use pmpimod
  use pfileptrmod
  use ivopotmod
  use tinvsubmod
  use eigenmod      !! IN PARENT DIRECTORY
  use lanblockmod   !! IN PARENT DIRECTORY
  use utilmod !! IN PARENT DIRECTORY
  implicit none
  DATATYPE,intent(inout) :: inspfs(totpoints,numspf)
  integer, intent(in) :: numloaded,numfrozen
  DATATYPE,intent(in) :: frozenreduced(totpoints)
  DATATYPE,allocatable :: lanspfs(:,:),density(:)
  DATAECS,allocatable :: energies(:)
  !$$
  DATAECS,allocatable :: mybigspfham(:,:), bigvects(:,:), bigenergies(:)
  !! AUTOMATIC ::  DATAECS :: mybigspfham(totpoints,totpoints), bigvects(totpoints,totpoints), bigenergies(totpoints)
  integer :: ibig,iorder,ispf,ppfac,ii,jj,kk,olist(numspf),flag
  integer :: null1,null2,null3,null4,null10(10),numcompute

  if (ivoflag.ne.0) then
     if (numloaded.le.0) then
        OFLWR "error, for ivo must load orbitals",numloaded; CFLST
     endif
     numcompute = numspf + max(0,num_skip_orbs) - numloaded
  else
     numcompute = numspf + num_skip_orbs
  endif

  if (numcompute.lt.0) then
     OFLWR "error numcompute lt 0 ", numloaded,numspf; CFLST
  endif
  if (numcompute.eq.0) then
     OFLWR "numcompute eq 0 RETURN"; CFL
     return
  endif
  if (numspf-numloaded.le.0) then
     OFLWR "numget eq 0 return", numspf-numloaded; CFL
     return
  endif

  if (ivoflag.ne.0) then
     kk=0
  else
     if (num_skip_orbs.lt.0) then
        kk=numloaded+num_skip_orbs
     else
        kk=numloaded
     endif
  endif

  ii=1
  do while (ii.le.numspf-numloaded)
     kk=kk+1
     flag=0
     do jj=1,num_skip_orbs
        if (orb_skip(jj).eq.kk) then
           flag=1
           exit
        endif
     enddo
     if (flag.eq.0) then
        olist(ii)=kk
        ii=ii+1
     endif
  enddo

  if (kk.gt.numcompute) then
     OFLWR "FSFFD EEEEE 555",kk,numspf,num_skip_orbs; CFLST
  endif

  if (ivoflag.ne.0) then
     OFLWR "getting IVO pot.  occupations are "
     WRFL loadedocc(1:numloaded); CFL

     allocate(ivopot(totpoints), density(totpoints),ivo_occupied(totpoints,numloaded))
     ivopot(:)=0d0; density(:)=0d0; ivo_occupied=0d0

     numocc=numloaded
     ivo_occupied(:,:)=inspfs(:,1:numloaded)
     do ispf=1,numloaded

        call gramschmidt(totpoints,ispf-1,totpoints,ivo_occupied(:,:),ivo_occupied(:,ispf),orbparflag)

        density(:)=density(:)+ivo_occupied(:,ispf)*CONJUGATE(ivo_occupied(:,ispf))*loadedocc(ispf)
     enddo
     call op_tinv(density,ivopot,1,1,null1,null2,null3,null4,null10)
     deallocate(density)

     if (numfrozen.gt.0) then
        ivopot(:)=ivopot(:)+frozenreduced(:)
     endif
  endif

  allocate(lanspfs(totpoints,numcompute),energies(numcompute))
  lanspfs=0; energies=0

  ibig=totpoints
  ppfac=1
  if (orbparflag) then
     ppfac=nprocs
  endif
  iorder=min(ibig*ppfac,orblanorder)


  if (eigmode.ne.0) then
     OFLWR "CALL BLOCK LAN FOR ORBS, ",numcompute," VECTORS",orbparflag; CFL
     call blocklanczos0(1,numcompute,ibig,ibig,iorder,ibig*ppfac,lanspfs,ibig,&
          energies,1,0,orblancheckmod,orblanthresh,mult_bigspf,orbparflag,orbtargetflag,orbtarget)
     OFLWR "BLOCKLAN CALLED. ENERGIES: ";CFL
  else
     OFLWR "Exact Diagonalization with vector size ",totpoints; CFL
     !$$
     allocate(mybigspfham(totpoints,totpoints),bigvects(totpoints,totpoints), &
          bigenergies(totpoints))
     call getbigspfham(mybigspfham(:,:))
     bigvects=0; bigenergies=0
     call ECSEIG(mybigspfham,ibig,ibig,bigvects,bigenergies)
     lanspfs(:,:) = bigvects(:,1:numcompute)
     energies(:) = bigenergies(1:numcompute)
     !$$
     deallocate(mybigspfham,bigvects,bigenergies)
     OFLWR "DIRECT DIAG CALLED. ENERGIES: ";CFL
  endif
  if (ivoflag.ne.0) then
     deallocate(ivopot,ivo_occupied)
  endif
  
  do ispf=1,numcompute
     OFLWR ispf,energies(ispf); CFL
  enddo
  OFLWR;  CFL

  do ispf=1,numspf-numloaded
     inspfs(:,ispf+numloaded)=lanspfs(:,olist(ispf))
  enddo

  deallocate(lanspfs,energies)

contains

  subroutine ivo_project(inbigspf,outbigspf)
    use myparams
    implicit none
    DATATYPE,intent(in) :: inbigspf(totpoints)
    DATATYPE, intent(out) :: outbigspf(totpoints)
    integer :: ii
    outbigspf(:)=0d0
    do ii=1,numocc
       outbigspf(:)=outbigspf(:) + ivo_occupied(:,ii) * ivodot(ivo_occupied(:,ii),inbigspf(:),totpoints)
    enddo

  end subroutine ivo_project

  function ivodot(inbra,inket,size)
    use myparams
    use mpisubmod    !! IN PARENT DIRECTORY
    implicit none
    integer,intent(in) :: size 
    DATATYPE,intent(in) :: inbra(size),inket(size)
    DATATYPE :: ivodot,csum
    csum=DOT_PRODUCT(inbra,inket)
    if (orbparflag) then
       call mympireduceone(csum)
    endif
    ivodot=csum
  end function ivodot

  subroutine mult_bigspf_ivo(inbigspf,outbigspf)
    use myparams
    use orbmultsubmod   !! IN PARENT DIRECTORY
    use orbprojectmod
    implicit none
    DATATYPE,intent(in) :: inbigspf(totpoints)
    DATATYPE, intent(out) :: outbigspf(totpoints)
    DATATYPE :: tempspf(totpoints),inwork(totpoints),inwork2(totpoints),&
         workspf(totpoints)   !! AUTOMATIC

    tempspf=0; inwork=0; inwork2=0; workspf=0

    call ivo_project(inbigspf,outbigspf)
    call project_onfrozen(inbigspf,workspf)
    outbigspf=outbigspf+workspf

    inwork2(:)=inbigspf(:)-outbigspf(:)

    outbigspf(:)=outbigspf(:) * (1d2)

    call mult_ke(inwork2(:),inwork(:),1,"booga",2)

    call mult_pot(1,inwork2(:),tempspf(:))

    inwork(:) = inwork(:) + tempspf(:) + ivopot(:)*inwork2(:)

    call ivo_project(inwork,inwork2)
    call project_onfrozen(inwork,workspf)
    inwork2=inwork2+workspf

    outbigspf(:) = outbigspf(:) + inwork(:) - inwork2(:)

  end subroutine mult_bigspf_ivo

  subroutine mult_bigspf0(inbigspf,outbigspf)
    use myparams
    use orbmultsubmod   !! IN PARENT DIRECTORY
    implicit none
    DATATYPE,intent(in) :: inbigspf(totpoints)
    DATATYPE, intent(out) :: outbigspf(totpoints)
    DATATYPE :: tempspf(totpoints)   !! AUTOMATIC

    tempspf=0

    call mult_ke(inbigspf(:),outbigspf(:),1,"booga",2)

    call mult_pot(1,inbigspf(:),tempspf(:))
    outbigspf(:)=outbigspf(:)+tempspf(:)

  end subroutine mult_bigspf0

  subroutine mult_bigspf(inbigspf,outbigspf)
    use myparams
    implicit none
    DATATYPE,intent(in) :: inbigspf(totpoints)
    DATATYPE, intent(out) :: outbigspf(totpoints)
    if (ivoflag.eq.0) then
       call mult_bigspf0(inbigspf,outbigspf)
    else
       call mult_bigspf_ivo(inbigspf,outbigspf)
    endif
  end subroutine mult_bigspf

  subroutine getbigspfham(bigspfham)
    use myparams
    use myprojectmod
    use orbmultsubmod   !! IN PARENT DIRECTORY
    implicit none
    DATAECS,intent(out) :: bigspfham(totpoints,totpoints)
    DATATYPE :: temppot(totpoints) , pot(totpoints)  !! AUTOMATIC
    integer :: ii
    
    bigspfham = 0; temppot=0; pot=0

    temppot = 1;
    call mult_pot(1,temppot(:),pot(:))
    
    bigspfham(:,:) = RESHAPE(ketot%mat(:,:,:,:),(/totpoints,totpoints/));
    do ii=1,totpoints
       bigspfham(ii,ii) = bigspfham(ii,ii) + pot(ii)
    enddo
    
  end subroutine getbigspfham

end subroutine init_spfs

end module initspfsmod


subroutine init_project(inspfs,spfsloaded,pot,halfniumpot,rkemod,proderivmod,skipflag,&
     bondpoints,bondweights,elecweights,elecradii,notused,& 
     numfrozen, infrozens, frozenkediag, frozenpotdiag, frozenreduced, notusedreduced)
  use myparams
  use pmpimod
  use pfileptrmod
  use myprojectmod
  use jfunctmod
  use initspfsmod
  use gettwoemod
  implicit none
  integer, intent(in) :: skipflag,notused
  integer,intent(inout) :: spfsloaded
  DATATYPE,intent(inout) :: inspfs(totpoints, numspf)
  DATATYPE,intent(out) :: pot(totpoints),proderivmod(numr,numr),rkemod(numr,numr),&
       bondpoints(numr),bondweights(numr), elecweights(totpoints,3),elecradii(totpoints),&
       halfniumpot(totpoints)
  integer,intent(in) :: numfrozen
  DATATYPE,intent(in) :: infrozens(totpoints,numfrozen)
  DATATYPE,intent(out) :: frozenkediag, frozenpotdiag, &
       frozenreduced(totpoints)
  DATATYPE :: notusedreduced(totpoints)

#ifndef REALGO
  real*8,allocatable :: temppot(:)
#endif
  real*8 :: rsum,pi
  character (len=2) :: th(4)
  integer ::  i,   j,ii,jj,k,l,ilow,ihigh

!! smooth exterior scaling
  DATATYPE,allocatable :: scalefunction(:), djacobian(:), ddjacobian(:)

  pi=4d0*atan(1d0)

#ifndef REALGO
  allocate(temppot(totpoints))
  temppot=0
#endif
  allocate(scalefunction(totpoints), djacobian(totpoints), ddjacobian(totpoints))
  scalefunction=0; djacobian=0; ddjacobian=0

  rkemod(:,:)=0d0; proderivmod(:,:)=0d0; bondpoints(:)=1d0; bondweights(:)=1d0

  elecweights(:,3)=(1d0/spacing)
  elecweights(:,1:2)=1d0

  call sineDVR(kevect%rmat(1-gridpoints:gridpoints-1),&
       fdvect%rmat(1-gridpoints:gridpoints-1), sinepoints%mat(:,:),gridpoints,spacing)

  kevect%cmat(:)=kevect%rmat(:)
  fdvect%cmat(:)=fdvect%rmat(:)

  ilow=1; ihigh=nbox
  if (orbparflag) then
     ilow=myrank; ihigh=myrank
  endif
  ii=(ilow-1)*numpoints

  if (toepflag.eq.0) then
     do i=ilow,ihigh
        do j=1,numpoints
           ii=ii+1
           jj=0
           do k=1,nbox
              do l=1,numpoints
                 jj=jj+1
                 ketot%mat(l,k,j,i)=kevect%rmat(ii-jj)
                 fdtot%mat(l,k,j,i)=fdvect%rmat(ii-jj)
                 ketot%tam(l,j,k,i)=kevect%rmat(ii-jj)
                 fdtot%tam(l,j,k,i)=fdvect%rmat(ii-jj)
              enddo
           enddo
        enddo
     enddo

     !! (e.g. coulmode = -1 disables centrifugal)
     !!
     if (twomode .eq. 1 .and. coulmode > -1) then  !! soft coulomb fix
        if (numpoints.ne.totpoints) then
           OFLWR "Error, bad points",numpoints,totpoints; CFLST
        endif
        if (toepflag.ne.0) then
           OFLWR "ERROR, toep not supported with twmode.ne.0 (softcoul)"; CFLST
        endif
        OFLWR "ADDING CENTRIFUGAL"; CFL
        call addit(ketot%mat)    !!! not ketot%tam, can't use ketot%tam
        ketot%tam = 0
     endif  !! two coulmode >=0  
  endif

  th=(/ "st", "nd", "rd", "th" /)

!!$  if (skipflag.gt.1) then
!!$     return
!!$  endif

  call get_dipoles()


  if (scalingflag.ne.0) then

     if (debugflag.eq.54321) then
        OFLWR "#Testing scaling..."; CFL
        if (myrank.eq.1) then
           open(76333,file="Contour.Dat",status="unknown")
           do i=-100,100
              rsum=i/100d0*spacing*gridpoints/2d0
              write(76333,'(100F12.5)') rsum, ffunct(rsum), jfunct(rsum), djfunct(rsum), ddjfunct(rsum)
           enddo
           close(76333)
        endif
        call mpibarrier()
        OFLWR "#Done testing scaling"; CFLST
     endif

     OFLWR "Getting scaling..."; CFL

     do i=1,totpoints
        scalefunction(i)=ffunct(real(dipoles(i),8))
        jacobian(i)=jfunct(real(dipoles(i),8))
        djacobian(i)=djfunct(real(dipoles(i),8))
        ddjacobian(i)=ddjfunct(real(dipoles(i),8))
     enddo

     invjacobian(:)=1d0/jacobian(:)

     invsqrtjacobian(:)=sqrt(invjacobian(:))

     invsqrtscaleweights(:) = sqrt(invjacobian(:)) 

     scalediag(:) = 3d0/8d0 * invjacobian(:)**4 * djacobian(:)**2 &
          - 1d0/4d0 * invjacobian(:)**3 * ddjacobian(:)

     elecweights(:,3)=elecweights(:,3)*jacobian(:)

     dipoles(:)=scalefunction(:)
     OFLWR "    ....Ok got scaling."; CFL
  endif

  elecradii(:)=sqrt( dipoles(:)**2 )
   
  call get_twoe_new(pot)

  if (numfrozen.gt.0) then
     call call_frozen_matels_core(infrozens,numfrozen,frozenkediag,frozenpotdiag,frozenreduced)
  endif

#ifndef REALGO

  if (capflag.gt.0) then
     temppot(:)=0d0
     do i=1,capflag
        if (capmode.eq.1) then
           temppot(:)=temppot(:) + capstrength(i)*( real(elecradii(:),8)/capstart(i) )**cappower(i)
        else
           temppot(:)=temppot(:) + capstrength(i)*( max(0d0,real(elecradii(:),8)-capstart(i)) )**cappower(i)
        endif
     enddo
     temppot(:)= min(maxcap,max(mincap,temppot(:)))
     halfniumpot(:) = (0d0,-1d0) * temppot(:)
     pot(:)=pot(:) + (0d0,-1d0) * temppot(:)
  else
#endif
     halfniumpot(:)=0d0
#ifndef REALGO
  endif
#endif

  if (spfsloaded.lt.numspf) then
     call init_spfs(inspfs(:,:),spfsloaded,numfrozen,frozenreduced)
  endif

!! now add in harmonic

  if (numcenters.gt.0) then
     pot(:) = pot(:) + 0.5d0 * harmstrength * dipoles(:)**2
  endif

  spfsloaded=numspf   !! for mcscf... really just for BO curve to skip eigen

  if (debugflag.eq.90210) then
     OFLWR "TEMPSTOP init_project debugflag 90210"; CFLST
  endif

#ifndef REALGO
  deallocate(temppot)
#endif
  deallocate(scalefunction, djacobian, ddjacobian)

  !! END INIT_PROJECT !!
  
contains

  !!
  !! for soft coulomb to reproduce p-wave eigvals need to fix up even parity.
  !!    add P 1/x^2 P  where P projects onto even.. l(l+1)/2 = 1 for p-wave centrifugal potential
  !! even 1d functions = p wave   odd 1d functions = s wave
  !!
  subroutine addit(inoutmat)
    use myparams
    use pfileptrmod
    implicit none
    DATATYPE,intent(inout) :: inoutmat(totpoints,totpoints)
    real*8, allocatable :: pproj(:,:), pproj2(:,:), myarray(:), allarray(:)
    integer :: ii, ihalf, ibot, itop, icenter
    real*8 :: dcenter

    if (coulmode<0) then
       OFLWR "error, addit called when coulmode<0, programmer fail"; CFLST
    endif
    
    if (twomode.ne.1) then
       OFLWR "error, addit called when twomode.ne.1: ", twomode; CFLST
    endif
    
    allocate(pproj(totpoints,totpoints), pproj2(totpoints,totpoints), &
         myarray(totpoints), allarray(totpoints))
    pproj=0; pproj2=0; myarray=0; allarray=0;
    
    myarray(:) = 0d0
    dcenter = (1+totpoints)/2d0
    do ii=1,totpoints
       myarray(ii) = (ii-dcenter)*spacing
    enddo

    do icenter=1,numcenters
       
       if (coulmode==0) then  !! integer quantum numbers even functions
          allarray(:)= 1d0 / ( softness**2 + ( myarray(:) - centershift(icenter)*spacing/2d0 )**2 )
       else                   !! half-integer quantum numbers even functions
          allarray(:)= 0.375d0 / ( softness**2 + ( myarray(:) - centershift(icenter)*spacing/2d0 )**2 )
       endif
       
       ihalf = floor((totpoints + 1 +1d-8 + centershift(icenter))/2d0)

       itop = min(totpoints,totpoints+centershift(icenter))
       ibot = max(1,centershift(icenter))

       pproj = 0
       do ii=1,totpoints
          pproj(ii,ii) = 1d0
       enddo
       do ii=ibot,ihalf
          pproj(          ii,           ii) =  0.5d0
          pproj(itop+ibot-ii,           ii) =  0.5d0
          pproj(          ii, itop+ibot-ii) =  0.5d0
          pproj(itop+ibot-ii, itop+ibot-ii) =  0.5d0
       enddo

       do ii=1,totpoints
          pproj2(ii,:) = pproj(ii,:) * allarray
       enddo
       inoutmat = inoutmat + MATMUL( pproj2, pproj )
       
    enddo

    deallocate(myarray,pproj,pproj2,allarray)
    
  end subroutine addit
  
  subroutine get_dipoles()
    use myparams
    use myprojectmod  
    implicit none
    call get_one_dipole(dipoles(:),qbox,1,1)
  end subroutine get_dipoles

  subroutine get_one_dipole(out,whichbox,nnn,mmm)
    use myparams
    use myprojectmod
    implicit none
    integer,intent(in) :: mmm,nnn,whichbox
    DATATYPE,intent(out) :: out(nnn,numpoints,mmm)
    integer :: jj,ii
    do jj=1,mmm
       do ii=1,numpoints
          out(:,ii,jj)=sinepoints%mat(ii,whichbox)
       enddo
    enddo
  end subroutine get_one_dipole

end subroutine init_project


!! fixed nuclei only for now

subroutine nucdipvalue(notused,dipoles)
  use myparams
  implicit none
  DATATYPE,intent(in) :: notused(1)
  DATATYPE,intent(out) :: dipoles(3)
  integer :: i
  dipoles(:)=0d0 
  do i=1,numcenters
     dipoles(3)=dipoles(3) + nuccharges(i) * centershift(i)/2d0 * spacing
  enddo
end subroutine nucdipvalue
