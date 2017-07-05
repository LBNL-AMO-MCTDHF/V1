
#include "Definitions.INC"

module myprojectmod
  implicit none

  DATATYPE, allocatable ::  threed_two(:)

  type fourmat
     DATATYPE, allocatable :: mat(:,:,:,:), tam(:,:,:,:)
  end type fourmat

  type twomat
     real*8, allocatable :: mat(:,:)
  end type twomat

  type onemat
     real*8, allocatable :: rmat(:)
     DATATYPE, allocatable :: cmat(:)
  end type onemat

  type(fourmat) :: ketot,fdtot
  type(twomat) ::  sinepoints
  type(onemat) :: kevect,fdvect

  DATATYPE, allocatable :: dipoles(:),&

!! WAS:  e.g. X(x) = x + i scalefunction(x,1)


       jacobian(:),&         !! jacobian(:,1) should only be a function of x, etc.
       invjacobian(:),&
       invsqrtjacobian(:),&
       scalediag(:),&
       invsqrtscaleweights(:)

end module myprojectmod


subroutine myprojectalloc()
  use myparams
  use pmpimod
  use pfileptrmod
  use myprojectmod
  implicit none

  allocate(dipoles(totpoints))
  dipoles=0
  
  if (scalingflag.ne.0) then
     allocate(          jacobian(totpoints),invjacobian(totpoints), &
          invsqrtjacobian(totpoints), &
          scalediag(totpoints),&
          invsqrtscaleweights(totpoints))
     jacobian=0; invjacobian=0; invsqrtjacobian=0; scalediag=0;
     invsqrtscaleweights=0;
  endif

  if (toepflag.eq.0) then
  
     if (numpoints*nbox.gt.10000) then
        OFLWR "WOW THAT'S BIG!  Are you sure you don't want to try toepflag?", numpoints*nbox; CFL
        OFLWR "WOW THAT'S BIG!  Are you sure you don't want to try toepflag?", numpoints*nbox; CFL
        OFLWR "WOW THAT'S BIG!  Are you sure you don't want to try toepflag?", numpoints*nbox; CFL
     endif

!! not needed now with tam.  old comment:
!! !! Allocating extra here for fdtot%mat and ketot%mat (+1's) --
!! !!   see Z/GEMM calls in coreproject.f90... leading dimension not
!! !!   allocated as passed to Z/GEMM without extra

     if (orbparflag) then
        allocate( &
!!             fdtot%mat(numpoints,nbox,numpoints,myrank:myrank +1), &
!!             ketot%mat(numpoints,nbox,numpoints,myrank:myrank +1))
             fdtot%mat(numpoints,nbox,numpoints,myrank:myrank), &
             ketot%mat(numpoints,nbox,numpoints,myrank:myrank),&
             fdtot%tam(numpoints,numpoints,nbox,myrank:myrank), &
             ketot%tam(numpoints,numpoints,nbox,myrank:myrank))
     else
        allocate( &
             fdtot%mat(numpoints,nbox,numpoints,nbox), &
             ketot%mat(numpoints,nbox,numpoints,nbox), &
             fdtot%tam(numpoints,numpoints,nbox,nbox), &
             ketot%tam(numpoints,numpoints,nbox,nbox))

     endif

     fdtot%mat=0; ketot%mat=0; fdtot%tam=0; ketot%tam=0
  endif

  allocate( &
       kevect%rmat(0-gridpoints:gridpoints-1),&
       kevect%cmat(0-gridpoints:gridpoints-1),&
       fdvect%rmat(0-gridpoints:gridpoints-1),&
       fdvect%cmat(0-gridpoints:gridpoints-1),&
       sinepoints%mat(numpoints,nbox))

  kevect%rmat=0; kevect%cmat=0;
  fdvect%rmat=0; kevect%cmat=0

  allocate(threed_two(0-numpoints:numpoints-1))
  threed_two=0

end subroutine myprojectalloc


module onedfunmod
contains
  
  function onedfun(inarray,num,incharge1,incharge2)
    use myparams
    implicit none
    integer,intent(in) :: num
    DATATYPE,intent(in) :: inarray(num)
    real*8, intent(in) :: incharge1, incharge2    
    DATATYPE :: onedfun(num)

    if (twomode==0) then
       onedfun = sechsq(inarray,num,incharge1,incharge2)
    else
       onedfun = softcoul(inarray,num,incharge1,incharge2)
    endif
    
  contains  
    function sechsq(inarray,num,incharge1,incharge2)
      use myparams
      use pfileptrmod
      implicit none
      integer,intent(in) :: num
      DATATYPE,intent(in) :: inarray(num)
      real*8, intent(in) :: incharge1, incharge2
      DATATYPE :: sechsq(num)
      real*8 :: xchargeprod     , esoft 
      
      xchargeprod = incharge1*incharge2

      esoft = softness;
      
      !      if (nucflag==0) then
      
      sechsq(:) = ( 2d0/(exp(inarray(:)/esoft) + exp((-1)*inarray(:)/esoft)) )**2 * &
           0.5d0 * xchargeprod*(xchargeprod+1d0/esoft)
      
      !      else
      !         !!
      !         !! this gets it to zero.  correction term such that united-atom limit is correct
      !         !!
      !         sechsq(:) = ( 2d0/(exp(inarray(:)/esoft) + exp((-1)*inarray(:)/esoft)) )**2 * &
      !              ( xchargeprod - (incharge1+incharge2)/esoft/2d0 - 0.25d0/esoft**2 &
      !              + 0.5/esoft * &
      !              sqrt(0.25/esoft**2 + incharge1**2 + incharge2**2 + (incharge1+incharge2)/esoft) )
      !      end if
    end function sechsq
    
    function softcoul(inarray,num,incharge1,incharge2)
      use myparams
      use pfileptrmod
      implicit none
      integer,intent(in) :: num
      DATATYPE,intent(in) :: inarray(num)
      real*8, intent(in) :: incharge1, incharge2
      real*8 :: xp
      DATATYPE :: softcoul(num)
      
      xp = incharge1*incharge2;

      softcoul(:) = xp / sqrt(softness**2 + inarray**2)
      
!$$      softcoul(:) = inarray(:)
!$$      call getsoftcoul(softcoul)
!$$      softcoul(:) = xp * softcoul(:)
      
    end function softcoul

  end function onedfun

!$$  subroutine getsoftcoul(inoutarray)
!$$    use myparams
!$$    use myprojectmod
!$$    implicit none
!$$    DATATYPE, intent(out) :: inoutarray(totpoints)
!$$    DATATYPE, allocatable :: hmat(:,:), hvects(:,:)
!$$    integer :: minindex, ii !! hack
!$$    real*8 :: minval
!$$      
!$$    allocate(hmat(totpoints,totpoints),hvects(totpoints,totpoints))
!$$      
!$$    inoutarray(:) = softness**2 +inoutarray(:)**2
!$$         
!$$    hmat = RESHAPE(ketot%mat,(/totpoints,totpoints/))
!$$    minval = 10**10
!$$    minindex = -8      
!$$    do ii=1,totpoints
!$$       if (abs(inoutarray(ii)) < minval) then
!$$          minindex = ii
!$$          minval = abs(inoutarray(ii))
!$$       end if
!$$       hmat(ii,:) = hmat(ii,:) * inoutarray(:)
!$$       hmat(:,ii) = hmat(:,ii) * inoutarray(:)
!$$    enddo
!$$      
!$$    inoutarray=0
!$$    inoutarray( minval ) = 1d0;
!$$    call hermsolve(totpoints,hmat,inoutarray)
!$$    print *, inoutarray
!$$    print *, "tempstop getsoftcoul"
!$$    stop
!$$    
!$$  contains
!$$    subroutine hermsolve(npts, inmat, inoutvect)
!$$      use pfileptrmod
!$$      implicit none
!$$      integer,intent(in) :: npts
!$$      DATATYPE, intent(in) :: inmat(npts,npts)
!$$      DATATYPE, intent(inout) :: inoutvect(npts)
!$$      DATATYPE, allocatable :: work(:)
!$$      integer,allocatable :: ipiv(:)
!$$      integer :: lwork, info
!$$      
!$$      lwork = 10*totpoints
!$$      allocate(ipiv(totpoints),work(lwork))
!$$      
!$$      !      
!$$      !      call ZHESV('U',totpoints,1,inmat,totpoints,ipiv,inoutvect,totpoints,work,lwork,info)
!$$      !
!$$      
!$$      call MYGESV(totpoints,1,inmat,totpoints,ipiv,inoutvect,totpoints,info)
!$$      if (info.ne.0) then
!$$         OFLWR "error zhesv hermsolve ", info; CFLST
!$$      endif
!$$      deallocate(ipiv,work)
!$$    end subroutine hermsolve
!$$    
!$$  end subroutine getsoftcoul
  

  
  function getscalefac(icenter)
    use myparams
    use pfileptrmod
    implicit none
    integer,intent(in) :: icenter
    integer :: jcenter
    real*8 :: getscalefac, fac, z1, z2, eee, esoft, ns(1)
    DATATYPE, parameter :: myzero(1)=0

    if (twomode.ne.0) then  !! can't use this with coulomb
       OFLWR "oops twomode getscalefac"; CFLST
    endif
    
    esoft = softness;
       
    z1=nuccharges(icenter)
    fac = 1d0
    do jcenter=1,numcenters
       if (jcenter.ne.icenter) then
          !
          ns(:) = real(onedfun(spacing*DATAONE*(centershift(icenter:icenter)-centershift(jcenter:jcenter))/2,1,1d0,1d0) / &
               onedfun(myzero,1,1d0,1d0),8)

          z2=nuccharges(jcenter)
          eee = ( -(z1+z2)/esoft    +   sqrt((z1+z2)**2/esoft**2 + &
               4*(z1**2+z2**2)*(z1**2 + z2**2 + 2*z1*z2 + (z1+z2)/esoft)) ) / 2 / (z1**2 + z2**2)
          
          fac = fac + eee*ns(1) + (1-ns(1)) - 1d0
       endif
    enddo
    getscalefac = fac
    
  end function getscalefac

  function getscalefac0(icenter)
    use myparams
    implicit none
    integer,intent(in) :: icenter
    integer :: jcenter
    real*8 :: getscalefac0, fac, z1, z2, eee, esoft
    real*8, parameter :: myone(1)=0
    
    esoft = softness;
       
    z1=nuccharges(icenter)
    fac = 1d0
    do jcenter=1,numcenters
       if (jcenter.ne.icenter) then
          !
          z2=nuccharges(jcenter)
          eee = ( -(z1+z2)/esoft    +   sqrt((z1+z2)**2/esoft**2 + &
               4*(z1**2+z2**2)*(z1**2 + z2**2 + 2*z1*z2 + (z1+z2)/esoft)) ) / 2 / (z1**2 + z2**2)
          
          fac = fac + eee - 1d0
       endif
    enddo
    getscalefac0 = fac
    
  end function getscalefac0

end module onedfunmod


module gettwoemod
contains

subroutine get_twoe_new(pot)
  use myparams
  use pfileptrmod
  use myprojectmod  
  use pmpimod
  use onedfunmod
  implicit none
  DATATYPE,intent(out) :: pot(totpoints)
  DATATYPE,allocatable :: myarray(:)
  integer :: ii,jj,pp,gridoffset,istart,icenter
  real*8 :: fac
  
  istart=1
  if (orbparflag.and.myrank.ne.1) then
     istart=0
  endif
  gridoffset=0; pp=1-gridpoints
  if (orbparflag) then
     gridoffset=(myrank-1)*numpoints
     pp=istart+2*gridoffset-gridpoints
  endif

  allocate(myarray(istart-numpoints:numpoints-1)); myarray=0

  jj=pp
  do ii=istart-numpoints,numpoints-1
     myarray(ii)=jj*spacing
     jj=jj+1
  enddo

  threed_two(:)=0d0

  if (twotype.eq.0) then
     threed_two(istart-numpoints:numpoints-1) = twostrength
  else
     threed_two(istart-numpoints:numpoints-1) = twostrength * &
          onedfun(myarray(:),2*numpoints-istart, 1d0, 1d0)
  endif

  deallocate(myarray); allocate(myarray(numpoints))

  if (numcenters.eq.0) then
     pot(:) = 0.5d0 * harmstrength * dipoles(:)**2
  else
!! Add in after get orbs  
     pot(:)=0d0
  endif

  do icenter=1,numcenters
     myarray(:)=( dipoles(:) - centershift(icenter)*spacing/2d0 )

     fac = 1d0
     if (twomode.eq.0) then    !! sechsq
        fac = getscalefac(icenter)
     endif
     
     pot(:)=pot(:) -  onedfun(myarray(:), numpoints, nuccharges(icenter)*fac, 1d0)

  enddo

  deallocate(myarray)

end subroutine get_twoe_new

end module gettwoemod


subroutine op_yderiv(notint,notused1,notused2)
  use pfileptrmod
  implicit none
  integer :: notint
  DATATYPE :: notused1(*), notused2(*)
  OFLWR "WHAT! no op_yderiv sincdvr, not yet."; CFLST
  notused1(1)=0*notused2(1)
end subroutine op_yderiv

