

!! INITIALIZATION ROUTINES FOR PROLATE AND ATOM: INIT_H2 AND INIT_HELIUM.  THESE CALL SUBROUTINES IN PROJECTS/ LIKE PSC.F90

#include "Definitions.INC"

!! This is dumb; use reshape

subroutine init_project(inspfs,spfsloaded,pot,halfniumpot,rkemod,proderivmod,skipflag,&
     bondpoints,bondweights,elecweights,elecradii,numelec )
  use myparams
  use myprojectmod
  implicit none
  integer, intent(in) ::  skipflag, numelec
  integer,intent(inout) :: spfsloaded
  DATATYPE,intent(inout) ::       inspfs(numerad,lbig+1, -mbig:mbig, numspf)
  DATATYPE,intent(out) :: proderivmod(numr,numr),rkemod(numr,numr),bondpoints(numr),&
       bondweights(numr),  halfniumpot(numerad,lbig+1, -mbig:mbig),pot(numerad,lbig+1, -mbig:mbig), &
       elecweights(numerad,lbig+1, -mbig:mbig),elecradii(numerad,lbig+1, -mbig:mbig)
  character (len=2) :: th(4)
  DATAECS, allocatable :: bigham(:,:,:,:), bigvects(:,:,:,:), bigvals(:)
  integer ::  i,ii,k,j,   taken(200)=0, flag,xiug, iug, ugvalue(200,0:10), getsmallugvalue

  halfniumpot=0d0
  do i=1,numspf
     if (abs(spfmvals(i)).gt.mbig) then
        OFLWR "MBIG INSUFFICIENT FOR SPFMVAL ",i,spfmvals(i),mbig; CFLST
     endif
  enddo

  call getjacobiKE(jacobipoints,jacobiweights, jacobike(:,:,0), lbig, jacobideriv(:,:,0),jacobirhoderiv(:,:),0)
  if (mbig.ge.1) then
     call getjacobiKE(jacobipoints,jacobiweights, jacobike(:,:,1), lbig, jacobideriv(:,:,1),jacobirhoderiv(:,:),1)
  endif

  do i=1,mbig
     jacobike(:,:,i)=jacobike(:,:,mod(i,2))
     jacobideriv(:,:,i)=jacobideriv(:,:,mod(i,2))
     do j=1,lbig+1

!! OBSERVATION!  must do this.  Didn't work the math out why.
        if (mod(i,2).eq.1.and.i.ne.1) then
           jacobike(j,j,i)=jacobike(j,j,i) - (i**2-1)/(1-jacobipoints(j)**2)
        else
           jacobike(j,j,i)=jacobike(j,j,i) - i**2/(1-jacobipoints(j)**2)
        endif
     enddo
  enddo

  jacobike(:,:,:)=jacobike(:,:,:)*(-0.5d0)


  call getLobatto(glpoints,glweights,glpoints2d,glweights2d,glke(:,:,0), henumpoints, henumelements, heelementsizes, hegridpoints, hecelement, heecstheta, glfirstdertot(:,:,0),glrhoderivs(:,:),0)

  call getLobatto(glpoints,glweights,glpoints2d,glweights2d,glke(:,:,1), henumpoints, henumelements, heelementsizes, hegridpoints, hecelement, heecstheta, glfirstdertot(:,:,1),glrhoderivs(:,:),1)


  if (temp_glflag.ne.0) then
     glke(:,:,1)=glke(:,:,0)
  endif

  do i=1,hegridpoints-2
     elecradii(i,:,:)=glpoints(i+1)
  enddo

  glpoints2d(:,2)=glpoints2d(:,1);   glweights2d(:,2)=glweights2d(:,1)

  th=(/ "st", "nd", "rd", "th" /)

  bondpoints(:)=1d0
  bondweights(:)=1d0
  do i=1,lbig+1
     elecweights(:,i,:)=jacobiweights(i)
  enddo
  do i=2,hegridpoints-1
     elecweights(i-1,:,:)=elecweights(i-1,:,:)*glweights(i)
  enddo
  rkemod(:,:)=0d0
  proderivmod(:,:)=0d0


!!$  if (skipflag.gt.1) then
!!$     return
!!$  endif

  allocate(bigham(numerad, lbig+1, numerad, lbig+1), bigvects(numerad,lbig+1, edim,0:mbig), bigvals(edim))
  do j=1,lbig+1
     do i=1,hegridpoints-2
        zdipole(i,j) = glpoints(i+1) * jacobipoints(j)
        xydipole(i,j) = glpoints(i+1) * sqrt(1.0d0-jacobipoints(j)**2)   !!TWOFIX / 2d0
     enddo
  enddo


  do ii=0,mbig

!! jacobideriv is matrix elements of ( (1-q^2) ddq - q) Where q=cos theta.

!! r ddz f =  [ 1/r ( (1-q^2) ddq - q ) + q ddr ] rf

!! such that ( (1-q^2) ddq - q) is antihermitian.

     do j=1,lbig+1
        do i=1,numerad
           do k=max(1,i-bandwidth),min(numerad,i+bandwidth)
              sparseddz_xi_banded(k-i+bandwidth+1,i,j,ii+1) = (-1) * glfirstdertot(k+1,i+1,mod(ii,2)) * jacobipoints(j)  * (-1)

              sparseddrho_xi_banded(k-i+bandwidth+1,i,j,ii+1) = glrhoderivs(k+1,i+1 ) * sqrt(1-jacobipoints(j)**2) /glpoints(k+1)    !!TWOFIX  /2d0

           enddo
           do k=1,lbig+1
              sparseddz_eta(k,j,i,ii+1) = (-1) * 1.d0/glpoints(i+1) * jacobideriv(k,j,ii)  * (-1)
              sparseddrho_eta(k,j,i,ii+1) =  (-1)* 1.d0/glpoints(i+1) * jacobirhoderiv(k,j) / sqrt(1-jacobipoints(k)**2)     !!TWOFIX  /2d0
           enddo
           sparseddz_diag(i,j,ii+1) = 0d0
           sparseddrho_diag(i,j,ii+1) = 0d0
        enddo
     enddo
  enddo
  do j=1,lbig+1
     do i=1,numerad
       ddrhopot(i,j)=1d0 /sqrt(1d0-jacobipoints(j)**2)  /glpoints(i+1)      !!TWOFIX  /2d0
     enddo
  enddo

  do ii=0,mbig

!!$!! 062111 added this logic here for atom too
!!$     if ((skipflag.eq.0).and.(spfsloaded.lt.numspf)) then

     if (spfsloaded.lt.numspf) then

        bigham=0.d0
        do i=1,hegridpoints-2
           bigham(i,:,i,:) = jacobike(:,:,abs(ii)) / glpoints(i+1)**2
        enddo


!! for centrifugal term add value of derivative of bra and ket at zero multiplied together

!! with setup analogous to prolate I don't get degeneracies even/odd m for sparse radial grid.
!!   seeing if this (accurate 1/r^2 integration for even m) will fix.
!!   did okay but not exact.
!!   could just replace gkle(:,:,1) with glke(:,:,0).
!! whatever.  will leave this as improved for usual (temp_glflag=0) 
!!   with temp_glflag option still available
!NOT!  doesn't work with sparseops setup

!if (mod(ii,2).eq.0.or.temp_glflag.ne.0) then
!     do i=1,hegridpoints-2
!        do j=1,hegridpoints-2
!           bigham(i,:,j,:)=bigham(i,:,j,:) + jacobike(:,:,abs(ii)) * &
!                glfirstdertot(1,i+1,0) *  glfirstdertot(1,j+1,0)    
!        enddo
!     enddo
!endif



        do i=1,lbig+1
           bigham(:,i,:,i) = bigham(:,i,:,i) + (-0.5d0)*glke(2:hegridpoints-1, 2:hegridpoints-1,mod(ii,2))
        enddo
        do i=1,hegridpoints-2
           do k=1,lbig+1
              bigham(i,k,i,k) = bigham(i,k,i,k) - nuccharge1 /glpoints(i+1)
           enddo
        enddo
     
        OFLWR "Get spfs.  Electronic dim, mval  ",edim, numerad,lbig+1, ii;CFL
        call ECSEIG(bigham,edim,edim,bigvects(:,:,:,ii),bigvals)

        OFL;        write(mpifileptr,'(A12,100F18.12)') "  Eigvals:"
#ifdef ECSFLAG
        write(mpifileptr,'(2F18.12)') bigvals(1:10)
#else
        write(mpifileptr,'(1F18.12)') bigvals(1:10)
#endif
        call closefile()
     endif


     bigham=0.d0
     do i=1,hegridpoints-2
        bigham(i,:,i,:) = jacobike(:,:,abs(ii)) / glpoints(i+1)**2
     enddo


!if (mod(ii,2).eq.0.or.temp_glflag.ne.0) then
!     do i=1,hegridpoints-2
!        do j=1,hegridpoints-2
!           bigham(i,:,j,:)=bigham(i,:,j,:) + jacobike(:,:,abs(ii)) * &
!                glfirstdertot(1,i+1,0) *  glfirstdertot(1,j+1,0)    
!        enddo
!     enddo
!endif


     do i=1,lbig+1
        bigham(:,i,:,i) = bigham(:,i,:,i) + (-0.5d0)*glke(2:hegridpoints-1, 2:hegridpoints-1,mod(ii,2))
     enddo
     do j=1,lbig+1
        do i=1,numerad
           do k=max(1,i-bandwidth),min(numerad,i+bandwidth)
              sparseops_xi_banded(k-i+bandwidth+1,i,j,ii+1) = bigham(k,j,i,j)
           enddo
           do k=1,lbig+1
              sparseops_eta(k,j,i,ii+1) = bigham(i,k,i,j)
           enddo
           sparseops_diag(i,j,ii+1) = bigham(i,j,i,j)
        enddo
     enddo
  enddo


  do i=1,hegridpoints-2
     pot(i,:,:) = (-1.0d0) * nuccharge1 /glpoints(i+1)
     halfniumpot(i,:,:) = (-1d0) /glpoints(i+1) * (nuccharge1 - numelec + 1)
  enddo
!!$  if (skipflag.ne.0) then
!!$     return
!!$  endif


  call openfile()
  if (spfsloaded.lt.numspf) then
     do ii=0,mbig
        !! get ug value
        do j=1,min(4*numspf,edim)
           ugvalue(j,ii)=getsmallugvalue((bigvects(:,:,j,abs(ii))),ii)
        enddo
        !           print *, "ugvalues m= ", ii," : ", ugvalue(1:numspf,ii)
     enddo
     
     do ii=-mbig,mbig
        taken(:)=0
        do i=1,num_skip_orbs
           if (orb_skip_mvalue(i).eq.ii) then
              taken(orb_skip(i))=1
           endif
        enddo
        do xiug=-1,1    !  iug=-1,0,1:   u, unspecified, or g sym...  now need taken array
           !!              want to do -1,1,0 so
           iug=floor(0.5d0*xiug-1.5d0*xiug**2+1 + 1.d-9)
           j=0
           do i=1,numspf
              if ((spfmvals(i)==ii).and.(spfugvals(i).eq.iug)) then           !! going through orbitals i with specified mvalue and ugvalue
                 
                 flag=0
                 do while (flag==0)                                          !! move to next vector, haven't found a bigvect of proper ugvalue yet
                    flag=0
                    do while (flag==0)                                         !! move to next calculated vector of specified mvalue
                       j=j+1
                       if (taken(j).eq.0) then
                          flag=1
                       endif
                    enddo
                    flag=0
                    if ((iug.eq.0).or.(iug.eq.ugvalue(j,abs(ii)))) then                 !! this one will do
                       flag=1
                    endif
                 enddo
                 if (i.le.spfsloaded) then
                    write(mpifileptr, *) "WON'T assign spf ", i, " to ", j, th(min(j,4)), " eigval of m=",ii, " because it is already loaded "
                 else
                    taken(j)=1
                    if (iug.eq.0) then
                       write(mpifileptr, *) "Assigning spf ", i, " to ", j, th(min(j,4)), " eigval of m=",ii, " ; ugvalue not fixed, is ", ugvalue(j,abs(ii))
                    else
                       write(mpifileptr, *) "Assigning spf ", i, " to ", j, th(min(j,4)), " eigval of m=",ii, " ; has specified ugvalue= ", spfugvals(i)
                    endif
                    inspfs(:,:,:,i)=0d0
                    inspfs(:,:,ii,i) = bigvects(:,:,j,abs(ii))
                 endif
              endif !! spfmval
           enddo !! nspf
        enddo !! iug
     enddo !! ii=-mbig,mbig
  else
     write(mpifileptr,*) "Found all spfs I need on file."
  endif  !! spfsloaded


  spfsloaded=numspf   !! for mcscf... really just for BO curve to skip eigen
  WRFL; CFL
  deallocate(bigham, bigvects,  bigvals)   !! twoe

  OFLWR "Ok init project.  Call get_twoe new"; CFL
  call get_twoe_new()
  OFLWR "Ok called get twoe_new"; CFL

end subroutine init_project

subroutine nucdipvalue(nullrvalue,dipoles)
  use myparams
  implicit none
  DATATYPE,intent(in) :: nullrvalue(1)
  DATATYPE,intent(out) :: dipoles(3)
  dipoles(:)=0d0*nullrvalue(1)  !! avoid warn unused
end subroutine nucdipvalue
