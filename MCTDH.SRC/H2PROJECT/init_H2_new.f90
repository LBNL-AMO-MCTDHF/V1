

!! INITIALIZATION ROUTINES FOR PROLATE AND ATOM: INIT_H2 AND INIT_HELIUM.  THESE CALL SUBROUTINES IN PROJECTS/ LIKE PSC.F90

#include "Definitions.INC"


subroutine init_project(inspfs,spfsloaded,pot,halfniumpot,rkemod,proderivmod,skipflag,&
     bondpoints,bondweights,elecweights,elecradii,in_numelec,&
     numfrozen, infrozens, frozenkediag, frozenpotdiag, frozenreduced, hatomreduced)
  use myparams
  use myprojectmod
  use eigenmod    !! IN PARENT DIRECTORY
  use pscmod
  use tinvsubmod
  use gettwoemod
  implicit none
  integer,intent(in) :: skipflag, in_numelec
  integer,intent(inout) :: spfsloaded
  DATATYPE,intent(inout) :: inspfs(numerad,lbig+1, -mbig:mbig, numspf)
  DATATYPE,intent(out) :: proderivmod(numr,numr),rkemod(numr,numr), &
       bondpoints(numr),bondweights(numr), halfniumpot(numerad,lbig+1, -mbig:mbig),&
       pot(numerad,lbig+1, -mbig:mbig), &
       elecweights(numerad,lbig+1, -mbig:mbig,3), elecradii(numerad,lbig+1, -mbig:mbig)
  integer,intent(in) :: numfrozen
  DATATYPE,intent(in) :: infrozens(numerad,lbig+1, -mbig:mbig,numfrozen)
  DATATYPE,intent(out) :: frozenkediag, frozenpotdiag, &
       frozenreduced(numerad,lbig+1, -2*mbig:2*mbig),&
       hatomreduced(numerad,lbig+1, -2*mbig:2*mbig)

  integer :: i,ii,imvalue,j, taken(200)=0, flag, jj, jflag, xiug, iug, ugvalue(200,0:30), &
       getsmallugvalue, istart
  DATAECS :: thisrvalue  
  character (len=2) :: th(4)
  DATAECS, allocatable :: bigham(:,:,:,:), bigvects(:,:,:,:), bigvals(:)
  DATATYPE,allocatable  ::  mydensity(:,:),ivopot(:,:),ivoproj(:,:,:,:)

  numelec=in_numelec

  th=(/ "st", "nd", "rd", "th" /)

  if ((mbig.gt.30).or.(numspf.gt.25)) then
     print *, "redim ugvalue";     stop
  endif
  do i=1,numspf
     if (abs(spfmvals(i)).gt.mbig) then
        OFLWR "MBIG INSUFFICIENT FOR SPFMVAL ",i,spfmvals(i),mbig; CFLST
     endif
  enddo

  OFLWR "Go init_project diatomic.   PSC..."; CFL

  call PSC()

  bondpoints(:)=rpoints(2:rgridpoints-1)
  if (bornopflag.eq.0) then
     bondweights(:)=rweights(2:rgridpoints-1)
  else
     bondweights(:)=1d0
  endif

  do i=1,lbig+1
     elecradii(:,i,:)=etapoints(i)**2
  enddo
  do i=2,xigridpoints-1
     elecradii(i-1,:,:)=elecradii(i-1,:,:) + xipoints(i)**2
  enddo


  elecweights(:,:,:,3)=1d0
  do i=1,lbig+1
     elecweights(:,i,:,2)=etaweights(i)
  enddo
  do i=2,xigridpoints-1
     elecweights(i-1,:,:,1)=xiweights(i)
  enddo

!!$      elecweights(i-1,:,:)=elecweights(i-1,:,:)*xiweights(i)

  rkemod=rketot(2:numr+1,2:numr+1)
  proderivmod(:,:) = prolate_derivs(2:numr+1,2:numr+1)

!! 01-2016 now here for frozen
  OFLWR "   ... ok PSC. call get_twoe_new"; CFL
  call get_twoe_new()
  OFLWR "   ... ok called get twoe_new"; CFL

  if (numhatoms.gt.0) then
     call hatomcalc(hatomreduced)
  endif
  if (numfrozen.gt.0) then
     call call_frozen_matels_core(infrozens,numfrozen,frozenkediag,frozenpotdiag,&
          frozenreduced,hatomreduced)
  endif


!!$  if (skipflag.gt.1) then
!!$     return
!!$  endif


!!MAY2014  thisrvalue=rpoints(2) 
  
  thisrvalue=(rpoints(2)+rpoints(rgridpoints-1))/2d0

  allocate(bigham(numerad, lbig+1, numerad, lbig+1),bigvects(numerad,lbig+1, edim,0:mbig), bigvals(edim))
  bigham=0; bigvals=0

  do imvalue=0,mbig
     jflag=0
     do jj=1,numspf
        if (abs(spfmvals(jj)).eq.imvalue) then
           jflag=1
        endif
     enddo

!!$  if ((skipflag.eq.0).and.(spfsloaded.lt.numspf).and.(jflag==1)) then
     if ((spfsloaded.lt.numspf).and.(jflag==1)) then
        bigham=proham(:,:,:,:,imvalue+1)/thisrvalue**2
        do i=1,lbig+1
           do j=1,numerad
              bigham(j,i,j,i) = bigham(j,i,j,i) + propot(j,i) / thisrvalue
           enddo
        enddo

        if (numfrozen.gt.0) then
           do i=1,lbig+1
              do j=1,numerad
                 bigham(j,i,j,i) = bigham(j,i,j,i) + &   !! ok conversion
                      frozenreduced(j,i,0) / thisrvalue  !! ok conversion
              enddo
           enddo
        endif

        if (numhatoms.gt.0) then
           do i=1,lbig+1
              do j=1,numerad
                 bigham(j,i,j,i) = bigham(j,i,j,i) + &   !! ok conversion
                      hatomreduced(j,i,0) / thisrvalue  !! ok conversion
              enddo
           enddo
        endif

        if (ivoflag.ne.0) then
           OFLWR "getting IVO pot.  occupations are "
           WRFL loadedocc(1:spfsloaded); CFL
           allocate(mydensity(numerad,lbig+1),ivopot(numerad,lbig+1))
           mydensity(:,:)=0d0; ivopot(:,:)=0d0;
           do ii=1,spfsloaded
              do j=-mbig,mbig
                 mydensity(:,:)=mydensity(:,:)+inspfs(:,:,j,ii)*CONJUGATE(inspfs(:,:,j,ii))*loadedocc(ii)
              enddo
           enddo
           call op_tinv(0,0,1,mydensity,ivopot)
           allocate(ivoproj(numerad,lbig+1,numerad,lbig+1))
           ivoproj(:,:,:,:)=0d0           
           do i=1,lbig+1
              do j=1,numerad
                 bigham(j,i,j,i) = bigham(j,i,j,i) + ivopot(j,i) / thisrvalue  !! ok conversion
                 ivoproj(j,i,j,i) = ivoproj(j,i,j,i) + 1d0
                 do ii=1,spfsloaded
                    ivoproj(:,:,j,i) = ivoproj(:,:,j,i) - &
                         inspfs(:,:,imvalue,ii)*CONJUGATE(inspfs(j,i,imvalue,ii))
                 enddo
              enddo
           enddo

!! DO BLAS with temporary arrays and WATCH DATA TYPES !!

           bigham(:,:,:,:) = RESHAPE(MATMUL(MATMUL(&
                RESHAPE(ivoproj(:,:,:,:),(/edim,edim/)),&
                RESHAPE(bigham(:,:,:,:),(/edim,edim/))),&
                RESHAPE(ivoproj(:,:,:,:),(/edim,edim/))),&
                (/numerad,lbig+1,numerad,lbig+1/))

           deallocate(mydensity,ivopot,ivoproj)
        endif

        OFLWR "Calculating orbitals.  Electronic dim, mval  ",edim, numerad,lbig+1, imvalue; CFL
        call ECSEIG(bigham,edim,edim,bigvects(:,:,:,imvalue),bigvals)

        !! fix phase, mostly for messflag to debug
        do j=1,min(edim,numspf)
           if (real(bigvects(1,1,j,imvalue)).lt.0.d0) then
              bigvects(:,:,j,imvalue)=(-1)*bigvects(:,:,j,imvalue)
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

  do imvalue=-mbig,mbig
     pot(:,:,imvalue)= propot(:,:)
     halfniumpot(:,:,imvalue)= halfpot(:,:) * (nuccharge1+nuccharge2-numelec+1)
  enddo

!!$  if (skipflag.ne.0) then
!!$     return
!!$  endif

  call openfile()
  if (spfsloaded.lt.numspf) then
     do imvalue=0,mbig
!! get ug value
        do j=1,min(4*numspf,edim)
           ugvalue(j,imvalue)=getsmallugvalue((bigvects(:,:,j,abs(imvalue))),imvalue)
        enddo
     enddo
     do imvalue=-mbig,mbig
        taken(:)=0
        do i=1,num_skip_orbs
           if (orb_skip_mvalue(i).eq.imvalue) then
              taken(orb_skip(i))=1
           endif
        enddo

        do xiug=-1,1    !  iug=-1,0,1:   u, unspecified, or g sym...  now need taken array
!!              want to do -1,1,0 so
           iug=floor(0.5d0*xiug-1.5d0*xiug**2+1 + 1.d-9)
           j=0
           if (ivoflag.eq.0) then
              istart=1
           else
              istart=spfsloaded+1
           endif
           do i=istart,numspf
              if ((spfmvals(i)==imvalue).and.(spfugvals(i).eq.iug)) then           !! going through orbitals i with specified mvalue and ugvalue
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
                    if ((iug.eq.0).or.(iug.eq.ugvalue(j,abs(imvalue)))) then                 !! this one will do
                       flag=1
                    endif
                 enddo
                 if (i.le.spfsloaded) then
                    write(mpifileptr, *) "WON'T assign spf ", i, " to ", j, th(min(j,4)), " eigval of m=",imvalue, " because it is already loaded "
                 else
                    taken(j)=1
                    if (iug.eq.0) then
                       write(mpifileptr, *) "Assigning spf ", i, " to ", j, th(min(j,4)), " eigval of m=",imvalue, " ; ugvalue not fixed, is ", ugvalue(j,abs(imvalue))
                    else
                       write(mpifileptr, *) "Assigning spf ", i, " to ", j, th(min(j,4)), " eigval of m=",imvalue, " ; has specified ugvalue= ", spfugvals(i)
                    endif
                    inspfs(:,:,:,i)=0d0
                    inspfs(:,:,imvalue,i) = bigvects(:,:,j,abs(imvalue))
                 endif
              endif !! spfmval
           enddo !! nspf
        enddo !! iug
     enddo !! imvalue
  else
     write(mpifileptr,*) "Found all spfs I need on file."
  endif  !! spfsloaded

  spfsloaded=numspf   !! for mcscf... really just for BO curve to skip eigen
  WRFL; CFL

  deallocate(bigham );  deallocate( bigvects, bigvals)   !! twoe

  OFLWR "Done init h2 project."; WRFL; CFL

end subroutine init_project


