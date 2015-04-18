

!! INITIALIZATION ROUTINES FOR PROLATE AND ATOM: INIT_H2 AND INIT_HELIUM.  THESE CALL SUBROUTINES IN PROJECTS/ LIKE PSC.F90

#include "Definitions.INC"


!!subroutine init_project(inspfs,spfsloaded,  pot,myrank,nprocs)


subroutine init_project(inspfs,spfsloaded,pot,halfniumpot,rkemod,proderivmod,skipflag,&
     bondpoints,bondweights,elecweights,elecradii,notused )
  use myparams
  use myprojectmod
  implicit none

  integer, intent(in) :: skipflag,notused
  integer ::  i,   spfsloaded,j,ii,jj,k,l,idim
  DATATYPE :: inspfs(totpoints, numspf),pot(totpoints),proderivmod(numr,numr),rkemod(numr,numr),&
       bondpoints(numr),bondweights(numr), elecweights(totpoints),elecradii(totpoints),&
       halfniumpot(totpoints)
  character (len=2) :: th(4)

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


!  littlepot(:,:) = sinep oints(:,:)**2  *(41d0/4d0)

  do idim=1,griddim
     littlepot(idim)%mat(:,:) = (32d0/4d0) * sinepoints(idim)%mat(:,:) **2

!RESHAPE(sinep oints(1+jj/2:max gridpoints-jj/2),(/numpoints(idim),nbox(idim)/))  **2

  enddo

  pot(:)=0d0; elecradii(:)=0d0

  call get_pot(pot(:))
  call get_rad(elecradii(:))

#ifndef REALGO
  if (capflag.ne.0) then
     pot(:)=pot(:) + capstrength * (0d0,-1d0) * max(0d0,real(elecradii(:),8)-capstart)**cappower
  endif
#endif

  call get_dipoles()

  call get_twoe_new(pot)

  halfniumpot(:)=pot(:)/sumcharge

  if (spfsloaded.lt.numspf) then
     call init_spfs(inspfs(:,:),spfsloaded)
  endif

  spfsloaded=numspf   !! for mcscf... really just for BO curve to skip eigen

end subroutine init_project



subroutine mult_bigspf(inbigspf,outbigspf)
  use myparams
  implicit none
  DATATYPE :: inbigspf(totpoints),outbigspf(totpoints)
  DATATYPE, allocatable :: tempspf(:)

  allocate(tempspf(totpoints))
  call mult_ke(inbigspf(:),outbigspf(:),1,"booga",2)

  call mult_pot(inbigspf(:),tempspf(:))
  outbigspf(:)=outbigspf(:)+tempspf(:)
  
  deallocatE(tempspf)

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

!     inspfs(:,ispf)=lanspfs(:,ispf-numloaded)

!     inspfs(:,ispf)=lanspfs(:,olist(ispf-numloaded))

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

