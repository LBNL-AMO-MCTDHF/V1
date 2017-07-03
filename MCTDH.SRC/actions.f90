
!! ACTION DRIVER SUBROUTINES AND SEVERAL ACTIONS LIKE PLOTTING

#include "Definitions.INC"

module actionlistmod
  
  integer, parameter :: MAXACTION=29

  character (len=18) :: action_list(MAXACTION) = (/ &
       "Autocorr          ", &      !! 1
       "Save nat          ", &      !! 2
       "Save spf          ", &      !! 3
       "Save den          ", &      !! 4
       "Save RNat         ", &      !! 5
       "Save NProj        ", &      !! 6
       "LanCurves         ", &      !! 7
       "Plot natorb       ", &      !! 8
       "Plot spf          ", &      !! 9
       "Plot density      ", &      !! 10
       "Plot RNat         ", &      !! 11
       "Plot NProj        ", &      !! 12
       "DissocFlux        ", &      !! 13
       "DissocFlux FT     ", &      !! 14
       "IonFlux Save      ", &      !! 15
       "IonFlux FT        ", &      !! 16
       "ProjIonFlux FT    ", &      !! 17
       "Plot denproj      ", &      !! 18
       "Enforce natorb    ", &      !! 19
       "Overlaps          ", &      !! 20
       "DipoleFT          ", &      !! 21
       "Check DF error    ", &      !! 22
       "Analyze overlaps  ", &      !! 23
       "KE projector      ", &      !! 24
       "psistats.dat      ", &      !! 25
       "Matrixelements    ", &      !! 26
       "IonFlux during    ", &      !! 27
       "ProjIon during    ", &      !! 28
       "Redo DipoleFT     " /)      !! 29

!! 0 = propagation action 1 = analysis action 2 = deprecated action
  integer :: action_type(MAXACTION) = (/ &
       0,&  !!       "Autocorr       ", &      !! 1
       0,&  !!       "Save nat       ", &      !! 2
       0,&  !!       "Save spf       ", &      !! 3
       0,&  !!       "Save den       ", &      !! 4
       0,&  !!       "Save RNat      ", &      !! 5
       0,&  !!       "Save NProj     ", &      !! 6
       0,&  !!       "LanCurves      ", &      !! 7
       1,&  !!       "Plot natorb    ", &      !! 8
       1,&  !!       "Plot spf       ", &      !! 9
       1,&  !!       "Plot density   ", &      !! 10
       2,&  !!       "Plot RNat      ", &      !! 11
       2,&  !!       "Plot NProj     ", &      !! 12
       2,&  !!       "DissocFlux     ", &      !! 13
       2,&  !!       "DissocFlux FT  ", &      !! 14
       0,&  !!       "IonFlux Save   ", &      !! 15
       1,&  !!       "IonFlux FT     ", &      !! 16
       1,&  !!       "ProjIonFlux FT ", &      !! 17
       1,&  !!       "Plot denproj   ", &      !! 18
       0,&  !!       "Enforce natorb ", &      !! 19
       0,&  !!       "Overlaps       ", &      !! 20
       0,&  !!       "DipoleFT       ", &      !! 21
       0,&  !!       "Check DF error ", &      !! 22
       1,&  !!       "Analyze overlap", &      !! 23
       0,&  !!       "KE projector   ", &      !! 24
       0,&  !!       "psistats.dat   ", &      !! 25
       1,&  !!       "Matrixelements ", &      !! 26
       0,&  !!       "IonFlux during ", &      !! 27
       0,&  !!       "ProjIon during " /)      !! 28
       1 /) !!    "Redo DipoleFT     " /)      !! 29

!! for windowing function.  &parinp namelist input.  Defaults in getparams.f90.

!! If nonzero, window function is ((tmax-t)/tmax)**ftwindowpower
 integer :: fttriwindow(MAXACTION)=1

!! If fttriwindow=0, window function is cos(pi t / 2 / tmax)**ftwindowpower 
 integer :: ftwindowpower(MAXACTION)=1

end module actionlistmod


subroutine get_skipflag_from_actions(outskip)
  use parameters
  use actionlistmod
  implicit none
  integer, intent(out) :: outskip
  integer :: i,j

  outskip=0
  j=0
  do while (j.lt.numactions)
     j=j+1

     if (action_type(actions(j)).lt.0.or.action_type(actions(j)).gt.2) then
        OFLWR "programmer error action_type ",actions(j),action_type(actions(j)); CFLST
     endif

     if (action_type(actions(j)).eq.2) then   !! analysis routine
        OFLWR "Error, action ", actions(j), action_list(actions(j))," is deprecated"; CFLST
     endif
        
     if (action_type(actions(j)).eq.1) then   !! analysis routine

!! only one analysis routine
        actions(1)=actions(j)
        numactions=1
        outskip=1
     endif

! never mind skipflag=2

!     if (((actions(j).gt.7).and.(actions(j).lt.13)).or.(actions(j).eq.14).or.(actions(j).eq.18)) then
!        outskip=2
!     endif
!     if ((actions(j).eq.16).or.(actions(j).eq.17).or.(actions(j).eq.23).or.(actions(j).eq.26)) then
!        outskip=1
!     endif

  enddo

  !! make sure no duplicates
  i=0
  do while (i.lt.numactions)
     i=i+1;     j=i
     do while ( j.lt. numactions)
        j=j+1
        if (actions(i)==actions(j)) then
           actions(j:numactions-1)=actions(j+1:numactions);           numactions=numactions-1
           j=j-1
        endif
     enddo
  enddo

end subroutine get_skipflag_from_actions


!! Windowfunct for actions 1 16 17 21.
!! now different windows for different actions 6-16 v1.31

function windowfunct(i,numdata,iaction)
  use parameters
  use actionlistmod
  implicit none
  real*8 :: windowfunct
  integer,intent(in) :: i,numdata,iaction

  if (iaction.gt.MAXACTION.or.iaction.lt.1) then
     OFLWR "windowfunct action error", iaction, MAXACTION; CFLST
  endif

  if (i.lt.0.or.i.gt.numdata) then
     OFLWR "ERROR, windowfunct ",i,numdata; CFLST
  endif

  if (ftwindowpower(iaction).eq.0) then
     windowfunct = 1d0
  else
     if (fttriwindow(iaction).ne.0) then
        windowfunct = ( real(numdata-i,8) / real(numdata,8) )**ftwindowpower(iaction)
     else
        windowfunct = cos( pi/2d0 * i / real(numdata,8) )**ftwindowpower(iaction)
     endif
  endif

end function windowfunct


subroutine actionsub(thistime)
  use parameters
  use mpimod
  use xxxmod
  use configmod
  use actionlistmod
  use autosubmod
  use dfconsubmod
  use dipsubmod
  use fluxgtaumod
  use keprojsubmod
  use ovlsubmod
  use projactionmod
  use psistatsubmod
  implicit none
  
  integer :: i, calledhere=0, atime, btime, times(MAXACTION)=0,getlen,myiostat
  real*8,intent(in) :: thistime
  CNORMTYPE :: error

  calledhere=calledhere+1

  do i=1,numactions
     select case (actions(i))
     case (1)    !! auto-correlation
        if (thistime.gt.autostart-1d-8) then
           call myclock(atime)
           call autocorrelate()
           call myclock(btime);        times(1)=times(1)+btime-atime
        endif
     case(2)
        call myclock(atime)
        if (myrank.eq.1) then
           if (mod(calledhere-1,plotmodulus).eq.0) then
                 call save_natorb( thistime, yyy%cmfspfs(:,0), yyy%denvects, yyy%denvals , 1)
           endif
        endif
        call myclock(btime);        times(2)=times(2)+btime-atime
     case(3)
        call myclock(atime)
        if (myrank.eq.1) then
           if (mod(calledhere-1,plotmodulus).eq.0) then
                 call save_spf( thistime, yyy%cmfspfs(:,0), 1)
           endif
        endif
        call myclock(btime);        times(3)=times(3)+btime-atime
     case(4)
        call myclock(atime)
        if (myrank.eq.1) then
           if (mod(calledhere-1,plotmodulus).eq.0) then
                 call save_density( thistime, yyy%cmfspfs(:,0),  yyy%denmat(:,:,0), 1, denplotbin)
           endif
        endif
        call myclock(btime);        times(4)=times(4)+btime-atime
     case(5)
        OFLWR "redo save rnatorb - need to calc"; CFLST
!        call myclock(atime)
!        if (myrank.eq.1) then
!           if (mod(calledhere-1,plotmodulus).eq.0) then
!                 call save_rnatorb(thistime, yyy%rdenvects(:,:), yyy%rdenvals(:), 1)
!           endif
!        endif
!        call myclock(btime);        times(5)=times(5)+btime-atime
     case(6)
        call myclock(atime)
        if (myrank.eq.1) then
           if (mod(calledhere-1,plotmodulus).eq.0) then
              call save_natproj( thistime )
           endif
        endif
        call myclock(btime);        times(6)=times(6)+btime-atime
     case(7)

        OFLWR "lancurves deprecated"; CFLST

     case(13)
        OFLWR "newflux not supported any more"; CFLST

     case(15)

        call myclock(atime)
        if(mod(calledhere-1,FluxInterval).eq.0) then  !! writes all mcscfnum
           call fluxwrite((calledhere-1)/FluxInterval,yyy%cmfspfs(:,0), yyy%cmfavec(:,:,0))
       endif
        call myclock(btime);        times(15)=times(15)+btime-atime

     case(19)
        call myclock(atime)
        call action_replacenat()
        call myclock(btime);        times(19)=times(19)+btime-atime
     case (20)    
        call myclock(atime)
        call getoverlaps(0)
        call myclock(btime);        times(20)=times(20)+btime-atime
     case (21)    
        call myclock(atime)
        call dipolesub()
        call myclock(btime);        times(21)=times(21)+btime-atime
     case (22)    
        call myclock(atime)
        call dferror(www,yyy%cptr(0),yyysptr(0),yyy%cmfavec(:,:,0),&
             mcscfnum,error,thistime)  !! does all mcscfnum
        OFL; write(mpifileptr,'(A15,2F25.10)') " DF error is ", error; CFL
        call myclock(btime);        times(22)=times(22)+btime-atime
     case (24)    
        if(mod(calledhere-1,FluxInterval).eq.0) then 
           call myclock(atime)
           call keprojector(yyy%cmfavec(:,:,0),yyy%cmfspfs(:,0),par_timestep*FluxInterval,www)
           call myclock(btime);        times(24)=times(24)+btime-atime
        endif
     case (25)
        if (mod(calledhere,psistatfreq).eq.0) then
           call myclock(atime)
           call psistats(thistime)
           call myclock(btime);        times(25)=times(25)+btime-atime
        endif
     case (27)    !! total ionization during calculation
        if (mod(calledhere-1,FluxInterval*FluxSkipMult).eq.0) then
           call fluxgtau_during(yyy%cmfspfs(:,0),yyy%cmfavec(:,:,0),par_timestep*FluxInterval*FluxSkipMult)
        endif
     case (28)
        if (mod(calledhere-1,FluxInterval*FluxSkipMult).eq.0) then
           call projeflux_during(yyy%cmfspfs(:,0),yyy%cmfavec(:,:,0),par_timestep*FluxInterval*FluxSkipMult)
        endif
     end select
  enddo

  if ((myrank.eq.1).and.(notiming.eq.0)) then
     if (calledhere.eq.1) then
        open(4132,file=timingdir(1:getlen(timingdir))//"/Actions.time.dat",&
             status="unknown",iostat=myiostat)
        call checkiostat(myiostat,"opening actions timing file")
        write(4132,'(500A18)',iostat=myiostat) (action_list(actions(i)),i=1,numactions)
        call checkiostat(myiostat,"writing actions timing file")
        close(4132)
     endif
     if (mod(calledhere,10).eq.1) then
        open(4132,file=timingdir(1:getlen(timingdir))//"/Actions.time.dat",&
             status="old",position="append",iostat=myiostat)
        call checkiostat(myiostat,"opening actions timing file")
        write(4132,'(500I15)',iostat=myiostat) (times(actions(i))/1000,i=1,numactions)
        call checkiostat(myiostat,"writing actions timing file")
        close(4132)
     endif
  endif
  
end subroutine actionsub


subroutine write_actions()
  use parameters
  use actionlistmod
  implicit none
  integer :: i,j,k

  do i=1,numactions
     if ((actions(i).gt.MAXACTION).or.(actions(i).lt.1)) then
        OFLWR "ACTION NOT SUPPORTED!!!! ", actions(i); CFLST
     endif
  enddo

  call openfile()
  if (numactions.gt.0) then
     write(mpifileptr,*) "ACTIONS:  ", (action_list(actions(i)),i=1,numactions)
  endif
  write(mpifileptr, *)

!  j=0
!  do i=1,numactions
!     if (actions(i).eq.13) then
!        j=1
!        exit
!     endif
!  enddo
!  if (j==1) then
!     write(mpifileptr,*) "FOR DISSOCIATIVE FLUX:"
!     write(mpifileptr,'(A20)') "Flux filenames:"
!     do i=1,numfluxfiles
!        write(mpifileptr,'(A20)') fluxfilenames(i)
!        write(mpifileptr,*) "   Curve # ", (whichfluxcurves(j,i), j=1,numfluxcurves(i))
!     enddo
!  endif

  j=0
  k=0
  do i=1,numactions
     if (actions(i).eq.1) then
        WRFL "For action 1 autostart= ",autostart," autotimestep= ",autotimestep
     endif
     if ((actions(i).eq.15).or.(actions(i).eq.16).or.(actions(i).eq.17)) then
        j=1
     endif
     if ((actions(i).eq.16).or.(actions(i).eq.17)) then
        k=1
     endif
     if (actions(i).eq.21) then
        WRFL "For emission/absorption action 21:"
        WRFL "   dipolesumstart/end ", dipolesumstart,dipolesumend
     endif
  enddo
  if (j==1) then
     write(mpifileptr,*) "FOR IONIZATION FLUX:"
     write(mpifileptr,*) "   Flux Interval ",FluxInterval,FluxInterval*par_timestep
     if (k==1) then
        write(mpifileptr,*) "   Fluxskipmult  ",FluxSkipMult
        if (nonuc_checkflag.eq.0) then
           write(mpifileptr,*) "   nucfluxopt  ", nucfluxopt
        endif
     endif
  endif

  call closefile()

end subroutine write_actions


subroutine actions_final()
  use parameters
  use mpimod
  use xxxmod
  use autosubmod
  use dipsubmod
  use fluxgtaumod
  use ovlsubmod
  implicit none
  integer :: i

  do i=1,numactions
     if (actions(i) == 1) then    !! auto-correlation
        call autocorrelate()
     endif
     if (actions(i) == 21) then   
        call dipolesub()
     endif
  enddo

  do i=1,numactions
     select case (actions(i))
     case(1)
        call autocorrelate_final()
     case(2)
        if (myrank.eq.1) then
           call save_natorb_final()
        endif
     case(3)
        if (myrank.eq.1) then
           call save_spf_final()
        endif
     case(4)
        if (myrank.eq.1) then
           call save_density_final()
        endif

     case(5)
        if (myrank.eq.1) then
           call save_rnatorb_final()
        endif

     case(6)
        if (myrank.eq.1) then
           call save_natproj_final()
        endif
     case(15)
        if(mod(numpropsteps,FluxInterval).eq.0) then
          call fluxwrite(numpropsteps/FluxInterval,yyy%cmfspfs(:,0), yyy%cmfavec(:,:,0))
       endif
    case(20)
        call getoverlaps(1)
     case(21)
        call dipolesub_final()
     end select
  enddo

end subroutine actions_final



!! FOR STARTUP OF ACTIONS ROUTINES (BEFORE PROPATAGATION)
!! ANALYSIS ROUTINES GO HERE.  NEVER GETS BEYOND ACTIONS_INITIAL.

subroutine actions_initial()
  use parameters
  use mpimod
  use xxxmod
  use autosubmod
  use dipsubmod
  use fluxgtaumod
  use ovlsubmod
  use projactionmod
  implicit none
  integer :: i, dflag

  dflag=0

  do i=1,numactions
     select case (actions(i))
     case (1)    !! auto-correlation
!!        call autocorrelate_initial() no, now in autocorrelate with autostart
     case (2)    
        if (myrank.eq.1) then
           call save_natorb_initial()
        endif
     case (3)    
        if (myrank.eq.1) then
           call save_spf_initial()
        endif
     case (4)
        dflag=1
        if (myrank.eq.1) then    
           call save_density_initial(denplotbin)
        endif
     case (5)    
        if (myrank.eq.1) then
           call save_rnatorb_initial()
        endif
     case (6)    
        call natprojalloc()
        dflag=1
        if (myrank.eq.1) then
           call save_denproj_initial(denprojplotbin)
           call save_natproj(-1.d0)
        endif
     case (7)
        OFLWR "lancurves not supported any more"; CFLST

     case (8)    
        call read_orb_initial(1)  !natorb
     case (9)    
        call read_orb_initial(2)  !spf
     case (10)    
        call read_orb_initial(3)  !den
     case (11)    
        call read_rorb_initial(1)  !! R-natorb
     case (12)    
        call read_rorb_initial(2)  !! Natproj
     case (13)
        OFLWR "newflux not supported any more"; CFLST
     case(14)
        OFLWR "newflux not supported any more"; CFLST
     case(15)
        OFLWR "Saving the wavefunction for reactive flux."; CFL
        call closefile()
     case(16)
        call fluxgtau(computeFlux)
     case(17)
          call projeflux_single(computeFlux)
     case (18)    
        call read_orb_initial(4)  !denproj
     case (19)    
     case (20)    !! overlaps onto final states
        call ovl_initial()
        call getoverlaps(0)
     case (21)    
        call dipolesub_initial()
     case (22)    
        if (df_restrictflag.eq.0) then
           OFLWR "Um, you are doing action 22 without dfrestrictflag..."; CFLST
        endif
        if (df_restrictflag.lt.2) then
           OFLWR "WARNING, usually need dfrestrictflag=2 or greater for action 22...."; CFL
        endif
     case (23)
        call wfnovl()
     case(24)
     case(25)
     case(26)
        call ovl_initial()
        call mcscf_matel()
     case (27)
     case (28)
     case (29)
        call dipolesub_initial()
        call redo_dipolesub(computeFlux,redobra)
     case default
        OFLWR "Action not supported: ", actions(i); CFLST
     end select
  enddo
end subroutine actions_initial


subroutine action_replacenat()
  use repnatmod
  use getstuffmod
  implicit none
  call replace_withnat(1)
  call get_stuff(0d0)
end subroutine action_replacenat



subroutine fixphase0(invector,isize)
  use parameters
  implicit none
  integer :: isize,ii
  DATATYPE :: invector(isize),csum
  real*8 :: rsum

  rsum=0d0; csum=0d0
  do ii=1,isize
     rsum=rsum+abs(invector(ii))**2
     csum=csum+invector(ii) * mod(floor(ii*45972.2203141d0+49.08d0/ii + ii**3),8)
  enddo
  if (rsum.eq.0) then
     return
  else
     if (abs(csum).eq.0) then
        return
     else
        invector(:)=invector(:)*abs(csum)/csum
     endif
  endif

end subroutine fixphase0

