

!! INITIALIZATION ROUTINES FOR PROLATE AND ATOM: INIT_H2 AND INIT_HELIUM.  THESE CALL SUBROUTINES IN PROJECTS/ LIKE PSC.F90

#include "Definitions.INC"


subroutine init_project(inspfs,spfsloaded,pot,halfniumpot,rkemod,proderivmod,skipflag,&
     bondpoints,bondweights,elecweights,elecradii,numelec)
  use myparams
  use myprojectmod
  implicit none

  integer,intent(in) :: skipflag, numelec
  
  DATATYPE :: inspfs(numerad,lbig+1, -mbig:mbig, numspf),proderivmod(numr,numr),rkemod(numr,numr), &
       bondpoints(numr),bondweights(numr), halfniumpot(numerad,lbig+1, -mbig:mbig),pot(numerad,lbig+1, -mbig:mbig), &
       elecweights(numerad,lbig+1, -mbig:mbig), elecradii(numerad,lbig+1, -mbig:mbig)
  integer :: i,ii,j,    taken(200)=0, flag, jj, jflag,spfsloaded, xiug, iug, ugvalue(200,0:30), getsmallugvalue
  DATAECS :: thisrvalue  
  character (len=2) :: th(4)
  DATAECS, allocatable :: bigham(:,:,:,:), bigvects(:,:,:,:), bigvals(:) !, regvects(:,:,:,:)
  th=(/ "st", "nd", "rd", "th" /)

  if ((mbig.gt.30).or.(numspf.gt.25)) then
     print *, "redim ugvalue";     stop
  endif
  do i=1,numspf
     if (abs(spfmvals(i)).gt.mbig) then
        OFLWR "MBIG INSUFFICIENT FOR SPFMVAL ",i,spfmvals(i),mbig; CFLST
     endif
  enddo

  call PSC()

  bondpoints(:)=rpoints(2:rgridpoints-1)
  bondweights(:)=rweights(2:rgridpoints-1)

  do i=1,lbig+1
     elecradii(:,i,:)=etapoints(i)**2
  enddo
  do i=2,xigridpoints-1
     elecradii(i-1,:,:)=elecradii(i-1,:,:) + xipoints(i)**2
  enddo


  do i=1,lbig+1
     elecweights(:,i,:)=etaweights(i)
  enddo
  do i=2,xigridpoints-1
     elecweights(i-1,:,:)=elecweights(i-1,:,:)*xiweights(i)
  enddo

  rkemod=rketot(2:numr+1,2:numr+1)
  proderivmod(:,:) = prolate_derivs(2:numr+1,2:numr+1)

  if (skipflag.gt.1) then
     return
  endif

!!MAY2014  thisrvalue=rpoints(2) 
  
  thisrvalue=(rpoints(2)+rpoints(rgridpoints-1))/2d0

  allocate(bigham(numerad, lbig+1, numerad, lbig+1),bigvects(numerad,lbig+1, edim,0:mbig), bigvals(edim))

  do ii=0,mbig
     jflag=0
     do jj=1,numspf
        if (abs(spfmvals(jj)).eq.ii) then
           jflag=1
        endif
     enddo

     if ((skipflag.eq.0).and.(spfsloaded.lt.numspf).and.(jflag==1)) then
        bigham=proham(:,:,:,:,ii+1)/thisrvalue**2
        do i=1,lbig+1
           do j=1,numerad
              bigham(j,i,j,i) = bigham(j,i,j,i) + propot(j,i) / thisrvalue
           enddo
        enddo

        OFLWR "Calculating orbitals.  Electronic dim, mval  ",edim, numerad,lbig+1, ii; CFL
        call ECSEIG(bigham,edim,edim,bigvects(:,:,:,ii),bigvals)

        !! fix phase, mostly for messflag to debug
        do j=1,min(edim,numspf)
           if (real(bigvects(1,1,j,ii)).lt.0.d0) then
              bigvects(:,:,j,ii)=(-1)*bigvects(:,:,j,ii)
           endif
        enddo

        OFLWR; WRFL "  Eigvals at R=", thisrvalue
#ifdef ECSFLAG
        write(mpifileptr,'(2F18.12)') bigvals(1:10) + NucCharge1*NucCharge2/thisrvalue
#else
        write(mpifileptr,'(F18.12)') bigvals(1:10) + NucCharge1*NucCharge2/thisrvalue
#endif
        call closefile()
     endif
  enddo

  do ii=-mbig,mbig
     pot(:,:,ii)= propot(:,:)
     halfniumpot(:,:,ii)= halfpot(:,:) * (nuccharge1+nuccharge2-numelec+1)
  enddo

  if (numhatoms.gt.0) then
     call hatomcalc()
  endif

  if (skipflag.ne.0) then
     return
  endif

  call openfile()
  if (spfsloaded.lt.numspf) then
     do ii=0,mbig
!! get ug value
        do j=1,min(4*numspf,edim)
           ugvalue(j,ii)=getsmallugvalue((bigvects(:,:,j,abs(ii))),ii)
        enddo
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
     enddo !! ii
  else
     write(mpifileptr,*) "Found all spfs I need on file."
  endif  !! spfsloaded

  spfsloaded=numspf   !! for mcscf... really just for BO curve to skip eigen
  WRFL; CFL

  deallocate(bigham );  deallocate( bigvects, bigvals)   !! twoe
  call openfile();  write(mpifileptr,*) "Done init_h2.";  write(mpifileptr, *);  call closefile()

  OFLWR "Ok init project.  Call get_twoe new"; CFL
  call get_twoe_new()
  OFLWR "Ok called get twoe_new"; CFL

end subroutine init_project


