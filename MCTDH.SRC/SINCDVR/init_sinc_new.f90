

#include "Definitions.INC"


subroutine init_project(inspfs,spfsloaded,pot,halfniumpot,rkemod,proderivmod,skipflag,&
     bondpoints,bondweights,elecweights,elecradii,notused )
  use myparams
  use myprojectmod
  implicit none
  integer, intent(in) :: skipflag,notused
  integer ::  i,   spfsloaded,j,ii,jj,k,l,idim,ilow(3),ihigh(3)
  DATATYPE :: inspfs(totpoints, numspf),pot(totpoints),proderivmod(numr,numr),rkemod(numr,numr),&
       bondpoints(numr),bondweights(numr), elecweights(totpoints),elecradii(totpoints),&
       halfniumpot(totpoints)
#ifndef REALGO
  real*8 :: temppot(totpoints) !! AUTOMATIC
#endif
  real*8 :: rsum,pi
  character (len=2) :: th(4)

!! smooth exterior scaling
  DATATYPE :: scalefunction(totpoints,3),scaleweights(totpoints), djacobian(totpoints,3)
  real*8 :: exteriorcoord(totpoints,3)
  DATATYPE :: fac

!!$ smooth, not exterior scaling
!!$#ifndef REALGO
!!$  DATATYPE :: scaledenom(totpoints,3),scaledenomderiv(totpoints,3),scaledenomsecderiv(totpoints,3),&
!!$       djacobian(totpoints,3), new_jacobian(totpoints,3),new_djacobian(totpoints,3),& !! AUTOMATIC
!!$       scalefunction(totpoints,3),scaleweights(totpoints)
!!$#endif


  if (fft_mpi_inplaceflag.eq.0) then
     call ct_init(fft_ct_paropt,mpifileptr)
  endif

  rkemod(:,:)=0d0; proderivmod(:,:)=0d0; bondpoints(:)=1d0; bondweights(:)=1d0

  elecweights(:)=(1d0/spacing)**griddim

  do idim=1,griddim

     call sineDVR(kevect(idim)%rmat(:),fdvect(idim)%rmat(:), sinepoints(idim)%mat(:,:),gridpoints(idim),spacing)

     kevect(idim)%cmat(:)=kevect(idim)%rmat(:)
     fdvect(idim)%cmat(:)=fdvect(idim)%rmat(:)

     ii=0
     do i=1,nbox(idim)
        do j=1,numpoints(idim)
           ii=ii+1
           jj=0
           do k=1,nbox(idim)
              do l=1,numpoints(idim)
                 jj=jj+1
                 ketot(idim)%mat(l,k,j,i)=kevect(idim)%rmat(ii-jj)
                 fdtot(idim)%mat(l,k,j,i)=fdvect(idim)%rmat(ii-jj)
              enddo
           enddo
        enddo
     enddo
  enddo

  th=(/ "st", "nd", "rd", "th" /)

  if (skipflag.gt.1) then
     return
  endif

  call get_rad(elecradii(:))
  call get_dipoles()


  fac=1d0

  if (scalingflag.ne.0) then
     fac=exp((0d0,1d0)*scalingtheta) !! ok imp conv mctdhf


     if (mod(scalingorder,2).eq.0) then
        OFLWR "Only odd scaling order supported right now", scalingorder; CFLST
     endif

     if (scalingorder.lt.3) then
        OFLWR "need scalingorder .ge. 3", scalingorder; CFLST
     endif

!! exterior smooth scaling
!! X(x,y,z) = x + scalefunction(x,1) etc.

     exteriorcoord(:,:)=&
          max(real(dipoles(:,:),8)-scalingdistance,0d0) + &
          min(real(dipoles(:,:),8)+scalingdistance,0d0)

     scalefunction(:,:) = fac * &
          exteriorcoord(:,:)**scalingorder / smoothness**scalingorder*scalingdistance

     jacobian(:,:) = 1 + fac * scalingorder * &
          exteriorcoord(:,:)**(scalingorder-1) / smoothness**scalingorder*scalingdistance

     djacobian(:,:) = fac * scalingorder * (scalingorder-1) * &
          exteriorcoord(:,:)**(scalingorder-2) / smoothness**scalingorder*scalingdistance

     invjacobian(:,:)=1d0/jacobian(:,:)
     scaleweights(:)=jacobian(:,1)*jacobian(:,2)*jacobian(:,3)

     invsqrtscaleweights(:)=sqrt(1d0/scaleweights(:))

     scaleder(:,:)=0.5d0 * invjacobian(:,:)**2 * djacobian(:,:)
     scalediag(:)=0d0
     do jj=1,3
        scalediag(:)=scalediag(:) - scaleder(:,jj)**2    !! MINUS  (negative 1 goes here, add scalediag to ke mult)
     enddo
  endif
  
   
!!$ smooth, not exterior scaling
!!$#ifndef REALGO
!!$  if (scalingflag.ne.0) then
!!$
!!$!! X(x,y,z) = x + scalefunction(x,1) etc.
!!$
!!$     djacobian(:,:)=0d0; jacobian(:,:)=0d0; scalefunction(:,:)=0d0
!!$
!!$     do jj=1,3
!!$        do i=0,scalingorders(jj)
!!$           scalefunction(:,jj)=scalefunction(:,jj) +   dipoles(:,jj)**i * scalingterms(i+1,jj) * (0d0,1d0)
!!$        enddo
!!$        do i=1,scalingorders(jj)
!!$           jacobian(:,jj)=jacobian(:,jj)       + i*dipoles(:,jj)**(i-1) * scalingterms(i+1,jj) * (0d0,1d0)
!!$        enddo
!!$        do i=2,scalingorders(jj)
!!$           djacobian(:,jj)=djacobian(:,jj)   + i*(i-1)*dipoles(:,jj)**(i-2) * scalingterms(i+1,jj) * (0d0,1d0)
!!$        enddo
!!$     enddo
!!$     
!!$     if (scalingdflag.ne.0) then
!!$        do jj=1,3
!!$           if (scalingorders(jj).lt.1) then
!!$              OFLWR "error, must have scalingorders>1 for scalingdflag", scalingorders(jj),jj,scalingdflag; CFLST
!!$           endif
!!$
!!$           scaledenomderiv(:,jj)=0d0; scaledenomsecderiv(:,jj) = 0d0
!!$
!!$           scaledenom(:,jj) = 1d0 +  scalingdconst(jj)* scalingterms(scalingorders(jj)+1,jj) * dipoles(:,jj)**(scalingorders(jj)-1)
!!$
!!$           if (scalingorders(jj).gt.1) then
!!$              scaledenomderiv(:,jj) = scalingdconst(jj)* scalingterms(scalingorders(jj)+1,jj) * dipoles(:,jj)**(scalingorders(jj)-2) * (scalingorders(jj)-1)
!!$           endif
!!$           if (scalingorders(jj).gt.2) then
!!$              scaledenomsecderiv(:,jj) = scalingdconst(jj)* scalingterms(scalingorders(jj)+1,jj) * dipoles(:,jj)**(scalingorders(jj)-3) * (scalingorders(jj)-1) * (scalingorders(jj)-2)
!!$           endif
!!$        enddo
!!$
!!$        new_jacobian(:,:) = ( jacobian(:,:) * scaledenom(:,:) - scalefunction(:,:) * scaledenomderiv(:,:) ) / scaledenom(:,:)**2
!!$
!!$        new_djacobian(:,:) = ( djacobian(:,:) * scaledenom(:,:) - 2* jacobian(:,:) * scaledenomderiv(:,:) - scalefunction(:,:) * scaledenomsecderiv(:,:) ) / scaledenom(:,:)**2 + &
!!$             2 * scalefunction(:,:) * scaledenomderiv(:,:)**2 / scaledenom(:,:)**3
!!$
!!$        scalefunction(:,:) = scalefunction(:,:) / scaledenom(:,:)
!!$
!!$        jacobian(:,:)=new_jacobian(:,:)
!!$        djacobian(:,:)=new_djacobian(:,:)
!!$
!!$     endif
!!$
!!$     scalefunction(:,:)=scalefunction(:,:) + dipoles(:,:)
!!$     jacobian(:,:) = jacobian + 1d0
!!$
!!$!! attempt at a kloodge
!!$     sumjacobian(:)=0d0
!!$     do jj=1,3
!!$        sumjacobian(:)=sumjacobian(:)+jacobian(:,jj) / 3d0
!!$     enddo
!!$
!!$
!!$     invjacobian(:,:)=1d0/jacobian(:,:)
!!$     scaleweights(:)=jacobian(:,1)*jacobian(:,2)*jacobian(:,3)
!!$     invsqrtscaleweights(:)=sqrt(1d0/scaleweights(:))
!!$
!!$     scaleder(:,:)=0.5d0 * invjacobian(:,:)**2 * djacobian(:,:)
!!$     scalediag(:)=0d0
!!$     do jj=1,3
!!$        scalediag(:)=scalediag(:) - scaleder(:,jj)**2    !! MINUS  (negative 1 goes here, add scalediag to ke mult)
!!$     enddo
!!$
!!$  endif
!!$
!!$#endif


  call get_twoe_new(pot)

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
     pot(:)=pot(:) + (0d0,-1d0) * temppot(:)
  endif
#endif

  ilow(1:3)=1
  ihigh(1:3)=gridpoints(1:3)

  if (orbparflag) then
     ilow(3)=(myrank-1)*numpoints(3)+1
     ihigh(3)=myrank*numpoints(3)
  endif

  pi=4d0*atan(1d0)

  if (maskflag.ne.0) then
     do jj=1,3
        maskfunction(jj)%rmat(:)=1d0
        do ii=1,masknumpoints

           rsum=cos(pi * ii / (masknumpoints+1) )

!!$           rsum=0.5d0 * ( cos( 2*pi * ii / (masknumpoints+1) ) + 1d0 )

           i=masknumpoints+1-ii

           if (i.ge.ilow(jj).and.i.le.ihigh(jj)) then
              maskfunction(jj)%rmat(i+1-ilow(jj)) = maskfunction(jj)%rmat(i+1-ilow(jj)) * rsum
           endif

           i=gridpoints(jj)-masknumpoints+ii

           if (i.ge.ilow(jj).and.i.le.ihigh(jj)) then
              maskfunction(jj)%rmat(i+1-ilow(jj)) = maskfunction(jj)%rmat(i+1-ilow(jj)) * rsum
           endif
        enddo
     enddo
  endif

  halfniumpot(:)=pot(:)/sumcharge

  if (spfsloaded.lt.numspf) then
     call init_spfs(inspfs(:,:),spfsloaded)
  endif

  spfsloaded=numspf   !! for mcscf... really just for BO curve to skip eigen

!!  OFLWR "TEMPSTOP init_project"; CFLST

end subroutine init_project


recursive subroutine mult_bigspf(inbigspf,outbigspf)
  use myparams
  implicit none
  DATATYPE :: inbigspf(totpoints),outbigspf(totpoints),tempspf(totpoints)

  call mult_ke(inbigspf(:),outbigspf(:),1,"booga",2)

  call mult_pot(inbigspf(:),tempspf(:))

  outbigspf(:)=outbigspf(:)+tempspf(:)
  
end subroutine mult_bigspf


subroutine init_spfs(inspfs,numloaded)
  use myparams
  implicit none
  DATATYPE :: inspfs(totpoints,numspf), energies(numspf+num_skip_orbs)
  DATATYPE,allocatable :: lanspfs(:,:)
  integer :: ibig,iorder,ispf,numloaded,ppfac,ii,jj,kk,olist(numspf),flag
  external :: mult_bigspf

  if (numloaded.ge.numspf) then
     OFLWR "already loaded ", numloaded,numspf; CFL
     return
  endif

!  if (num_skip_orbs.lt.0) then
!     OFLWR "WOOWOWWO"; CFLST
!  endif
  if (num_skip_orbs+numloaded.lt.0) then
     OFLWR "WOOWOWWOxxx",num_skip_orbs,numloaded; CFLST
  endif
  if (num_skip_orbs+numspf.le.0) then
     OFLWR "WOOWOWWOyyy",num_skip_orbs,numspf; CFLST
  endif

  do ii=1,num_skip_orbs
     if (orb_skip(ii).lt.1.or.orb_skip(ii).gt.numspf+num_skip_orbs) then
        OFLWR "TTTT DOG!"; CFLST
     endif
  enddo
  kk=0-min(num_skip_orbs,0)
  olist(:)=(-9999999)
  do ii=1,numspf+num_skip_orbs
     flag=0
     do jj=1,num_skip_orbs
        if (orb_skip(jj).eq.ii) then
           flag=1
           exit
        endif
     enddo
     if (flag.eq.0) then
        kk=kk+1
        olist(kk)=ii
     endif
  enddo
  if (kk.ne.numspf) then
     OFLWR "FSFFD EEEEE 555",kk,numspf,num_skip_orbs; CFLST
  endif

  allocate(lanspfs(totpoints,numspf+num_skip_orbs))

  ibig=totpoints
  iorder=min(ibig,orblanorder)
  ppfac=1
  if (orbparflag) then
     ppfac=nprocs
  endif

  OFLWR "CALL BLOCK LAN FOR ORBS, ",numspf+num_skip_orbs," VECTORS"; CFL

  call blocklanczos0(min(3,numspf),numspf+num_skip_orbs,ibig,ibig,iorder,ibig*ppfac,lanspfs,ibig,energies,1,0,orblancheckmod,orblanthresh,mult_bigspf,orbparflag)
  
  OFLWR "BL CALLED. ENERGIES: ";CFL
  do ispf=1,numspf+num_skip_orbs
     OFLWR ispf,energies(ispf); CFL
  enddo
  OFLWR;  CFL

  do ispf=numloaded+1,numspf
     inspfs(:,ispf)=lanspfs(:,olist(ispf))
  enddo

  deallocate(lanspfs)

end subroutine init_spfs

subroutine nucdipvalue(notused,dipoles)
  use myparams
  implicit none
  DATATYPE :: notused(1),dipoles(3)
  dipoles(:)=0d0
end subroutine nucdipvalue

