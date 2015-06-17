

#include "Definitions.INC"


subroutine init_project(inspfs,spfsloaded,pot,halfniumpot,rkemod,proderivmod,skipflag,&
     bondpoints,bondweights,elecweights,elecradii,notused )
  use myparams
  use pmpimod
  use pfileptrmod
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
     call ct_init(fft_ct_paropt)
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

!!$  later, just from dipoles
!!$  call get_rad(elecradii(:))

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
     invsqrtjacobian(:,:)=1d0/sqrt(jacobian(:,:))
     scaleder2(:,:)=1d0*invsqrtjacobian(:,:)**3 * djacobian(:,:)

     scaleweights(:)=jacobian(:,1)*jacobian(:,2)*jacobian(:,3)

!!in case angle gt 60 use next formula     invsqrtscaleweights(:)=sqrt(1d0/scaleweights(:))

     invsqrtscaleweights(:)=invsqrtjacobian(:,1)*invsqrtjacobian(:,2)*invsqrtjacobian(:,3)

     scaleder(:,:)=0.5d0 * invjacobian(:,:)**2 * djacobian(:,:)
     scalediag(:)=0d0
     do jj=1,3
        scalediag(:)=scalediag(:) - scaleder(:,jj)**2    !! MINUS  (negative 1 goes here, add scalediag to ke mult)
     enddo

     elecweights(:)=elecweights(:)*scaleweights(:)

     dipoles(:,:)=dipoles(:,:)+exteriorcoord(:,:)

  endif

  elecradii(:)=sqrt( dipoles(:,1)**2 + dipoles(:,2)**2 + dipoles(:,3)**2)
   
  call get_twoe_new(pot)

  if (debugflag .eq. 4040.or.debugflag.eq.3939) then
     pot(:) = 0.35d0 * elecradii(:)**2 * ( &
          exp((-0.13d0)*(dipoles(:,1)**2+dipoles(:,2)**2+(dipoles(:,3)-2d0)**2)) + &
          exp((-0.13d0)*(dipoles(:,1)**2+dipoles(:,2)**2+(dipoles(:,3)+2d0)**2)) )
     threed_two(:,:,:,:)=0d0
  else if (debugflag .eq. 4242.or.debugflag.eq.4141) then
     pot(:) = 7.5d0 * elecradii(:)**2 * exp((-1)*elecradii(:))
     threed_two(:,:,:,:)=0d0
  else if (debugflag.eq.4343) then
     pot(:) = 4.5d0 * elecradii(:)**2 
     threed_two(:,:,:,:)=0d0
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
     pot(:)=pot(:) + (0d0,-1d0) * temppot(:)
  endif
#endif

  ilow(1:3)=1
  ihigh(1:3)=gridpoints(1:3)

  if (orbparflag) then
     ilow(orbparlevel:3)=(boxrank(orbparlevel:3)-1)*numpoints(orbparlevel:3)+1
     ihigh(orbparlevel:3)=boxrank(orbparlevel:3)*numpoints(orbparlevel:3)
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
  DATATYPE,intent(in) :: inbigspf(totpoints)
  DATATYPE, intent(out) :: outbigspf(totpoints)
  if (ivoflag.eq.0) then
     call mult_bigspf0(inbigspf,outbigspf)
  else
     call mult_bigspf_ivo(inbigspf,outbigspf)
  endif
end subroutine mult_bigspf


recursive subroutine mult_bigspf0(inbigspf,outbigspf)
  use myparams
  implicit none
  DATATYPE,intent(in) :: inbigspf(totpoints)
  DATATYPE, intent(out) :: outbigspf(totpoints)
  DATATYPE :: tempspf(totpoints)

  call mult_ke(inbigspf(:),outbigspf(:),1,"booga",2)
  call mult_pot(inbigspf(:),tempspf(:))
  outbigspf(:)=outbigspf(:)+tempspf(:)

end subroutine mult_bigspf0


function ivodot(inbra,inket,size)
  use myparams
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


module ivopotmod
  implicit none
  integer :: numocc=(-1)
  DATATYPE, allocatable :: ivopot(:),ivo_occupied(:,:)
end module ivopotmod


subroutine ivo_project(inbigspf,outbigspf)
  use myparams
  use ivopotmod
  implicit none
  DATATYPE,intent(in) :: inbigspf(totpoints)
  DATATYPE, intent(out) :: outbigspf(totpoints)
  DATATYPE :: ivodot
  integer :: ii
  outbigspf(:)=0d0
  do ii=1,numocc
     outbigspf(:)=outbigspf(:) + ivo_occupied(:,ii) * ivodot(ivo_occupied(:,ii),inbigspf(:),totpoints)
  enddo
end subroutine ivo_project

  
recursive subroutine mult_bigspf_ivo(inbigspf,outbigspf)
  use myparams
  use ivopotmod
  implicit none
  DATATYPE,intent(in) :: inbigspf(totpoints)
  DATATYPE, intent(out) :: outbigspf(totpoints)
  DATATYPE :: tempspf(totpoints),inwork(totpoints),inwork2(totpoints)

  call ivo_project(inbigspf,outbigspf)

  inwork2(:)=inbigspf(:)-outbigspf(:)

  outbigspf(:)=outbigspf(:) * (1d2)

  call mult_ke(inwork2(:),inwork(:),1,"booga",2)

  call mult_pot(inwork2(:),tempspf(:))

  inwork(:) = inwork(:) + tempspf(:) + ivopot(:)*inwork2(:)

  call ivo_project(inwork,inwork2)

  outbigspf(:) = outbigspf(:) + inwork(:) - inwork2(:)


end subroutine mult_bigspf_ivo



subroutine init_spfs(inspfs,numloaded)
  use myparams
  use pmpimod
  use pfileptrmod
  use ivopotmod
  implicit none
  DATATYPE :: inspfs(totpoints,numspf)
  DATATYPE,allocatable :: lanspfs(:,:),density(:), energies(:)
  integer, intent(in) :: numloaded
  integer :: ibig,iorder,ispf,ppfac,ii,jj,kk,olist(numspf),flag
  integer :: null1,null2,null3,null4,null10(10),numcompute
  external :: mult_bigspf

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
     numocc=numloaded
     ivo_occupied(:,:)=inspfs(:,1:numloaded)
     do ispf=1,numloaded

        call myhgramschmidt_fast(totpoints,ispf-1,totpoints,ivo_occupied(:,:),ivo_occupied(:,ispf),orbparflag)

        density(:)=density(:)+abs(ivo_occupied(:,ispf)**2)*loadedocc(ispf)
     enddo
     call op_tinv(density,ivopot,1,1,null1,null2,null3,null4,null10)
     deallocate(density)
  endif

  allocate(lanspfs(totpoints,numcompute),energies(numcompute))

  ibig=totpoints
  iorder=min(ibig,orblanorder)
  ppfac=1
  if (orbparflag) then
     ppfac=nprocs
  endif

  OFLWR "CALL BLOCK LAN FOR ORBS, ",numcompute," VECTORS"; CFL

  call blocklanczos0(min(3,numspf),numcompute,ibig,ibig,iorder,ibig*ppfac,lanspfs,ibig,energies,1,0,orblancheckmod,orblanthresh,mult_bigspf,orbparflag,orbtargetflag,orbtarget)

  if (ivoflag.ne.0) then
     deallocate(ivopot,ivo_occupied)
  endif
  
  OFLWR "BL CALLED. ENERGIES: ";CFL
  do ispf=1,numcompute
     OFLWR ispf,energies(ispf); CFL
  enddo
  OFLWR;  CFL

  do ispf=1,numspf-numloaded
     inspfs(:,ispf+numloaded)=lanspfs(:,olist(ispf))
  enddo

  deallocate(lanspfs,energies)

end subroutine init_spfs

subroutine nucdipvalue(notused,dipoles)
  use myparams
  implicit none
  DATATYPE :: notused(1),dipoles(3)
  dipoles(:)=0d0
end subroutine nucdipvalue

