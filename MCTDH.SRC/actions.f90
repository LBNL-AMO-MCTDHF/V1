
!! ACTION DRIVER SUBROUTINES AND SEVERAL ACTIONS LIKE PLOTTING

!! actions_initial may change yyy%cmfpsivec; actions may not!!!

#include "Definitions.INC"

module actionlistmod
  
  character (len=15) :: action_list(25) = (/ &
       "Autocorr       ", &      !! 1
       "Save nat       ", &      !! 2
       "Save spf       ", &      !! 3
       "Save den       ", &      !! 4
       "Save RNat      ", &      !! 5
       "Save NProj     ", &      !! 6
       "LanCurves      ", &      !! 7

       "Plot natorb    ", &      !! 8
       "Plot spf       ", &      !! 9
       "Plot density   ", &      !! 10
       "Plot RNat      ", &      !! 11
       "Plot NProj     ", &      !! 12
       "Dissoc.FLUX    ", &      !! 13
       "Dissoc.Flux  FT", &      !! 14
       "Ion.  Flux Save", &      !! 15
       "Ion.  Flux  FT ", &      !! 16
       "Proj Ion Flux  ", &      !! 17
       "Plot denproj   ", &      !! 18
       "Enforce natorb ", &      !! 19
       "Final state ovl", &      !! 20
       "DipoleFT       ", &      !! 21
       "Check D-F error", &      !! 22
       "Wfn ovls       ", &      !! 23
       "KE projector   ", &      !! 24
       "psistats.dat   " /)      !! 25

end module actionlistmod


subroutine actionsub(thistime)
  use parameters
  use mpimod
  use xxxmod
  use actionlistmod
  implicit none
  
  integer :: i, calledhere=0, atime, btime, times(25)=0,getlen
  real*8 :: thistime
  CNORMTYPE :: error

  calledhere=calledhere+1

  do i=1,numactions
     select case (actions(i))
     case (1)    !! auto-correlation
        call system_clock(atime)
        call autocorrelate()
        call system_clock(btime);        times(1)=times(1)+btime-atime
     case(2)
        call system_clock(atime)
        if (myrank.eq.1) then
           if (mod(calledhere-1,plotmodulus).eq.0) then
                 call save_natorb( thistime, yyy%cmfpsivec(spfstart,0), yyy%denvects, yyy%denvals , 1)
           endif
        endif
        call system_clock(btime);        times(2)=times(2)+btime-atime
     case(3)
        call system_clock(atime)
        if (myrank.eq.1) then
           if (mod(calledhere-1,plotmodulus).eq.0) then
                 call save_spf( thistime, yyy%cmfpsivec(spfstart,0), 1)
           endif
        endif
        call system_clock(btime);        times(3)=times(3)+btime-atime
     case(4)
        call system_clock(atime)
        if (myrank.eq.1) then
           if (mod(calledhere-1,plotmodulus).eq.0) then
                 call save_density( thistime, yyy%cmfpsivec(spfstart,0),  yyy%denmat(:,:,0), 1, denplotbin)
           endif
        endif
        call system_clock(btime);        times(4)=times(4)+btime-atime
     case(5)
        OFLWR "redo save rnatorb - need to calc"; CFLST
!        call system_clock(atime)
!        if (myrank.eq.1) then
!           if (mod(calledhere-1,plotmodulus).eq.0) then
!                 call save_rnatorb(thistime, yyy%rdenvects(:,:), yyy%rdenvals(:), 1)
!           endif
!        endif
!        call system_clock(btime);        times(5)=times(5)+btime-atime
     case(6)
        call system_clock(atime)
        if (myrank.eq.1) then
           if (mod(calledhere-1,plotmodulus).eq.0) then
              call save_natproj( thistime )  !!, yyy%cmfpsivec(astart(1),0))
           endif
        endif
        call system_clock(btime);        times(6)=times(6)+btime-atime
     case(7)

        OFLWR "lancurves deprecated"; CFLST

     case(13)
        OFLWR "newflux not supported any more"; CFLST

     case(15)

        call system_clock(atime)
        if(mod(calledhere-1,FluxInterval).eq.0) then  !! writes all mcscfnum
           call fluxwrite((calledhere-1)/FluxInterval,yyy%cmfpsivec(spfstart,0),yyy%cmfpsivec(astart(1),0))
       endif
        call system_clock(btime);        times(15)=times(15)+btime-atime

     case(19)
        call system_clock(atime)
        call action_replacenat()
        call system_clock(btime);        times(19)=times(19)+btime-atime
     case (20)    
        call system_clock(atime)
        call getoverlaps(0)
        call system_clock(btime);        times(20)=times(20)+btime-atime
     case (21)    
        call system_clock(atime)
        call dipolesub()
        call system_clock(btime);        times(21)=times(21)+btime-atime
     case (22)    
        call system_clock(atime)
        call dferror(yyy%cmfpsivec(astart(1),0),error,thistime)  !! does all mcscfnum
        OFL; write(mpifileptr,'(A15,2F25.10)') " DF error is ", error; CFL
        call system_clock(btime);        times(22)=times(22)+btime-atime
     case (24)    
        if(mod(calledhere-1,FluxInterval).eq.0) then 
           call system_clock(atime)
           call keprojector(yyy%cmfpsivec(astart(1),0),yyy%cmfpsivec(spfstart,0),par_timestep*FluxInterval)
           call system_clock(btime);        times(24)=times(24)+btime-atime
        endif
     case (25)
        if (mod(calledhere,psistatfreq).eq.0) then
           call system_clock(atime)
           call psistats(thistime)
           call system_clock(btime);        times(25)=times(25)+btime-atime
        endif
     end select

  enddo

  if ((myrank.eq.1).and.(notiming.eq.0)) then
     if (calledhere.eq.1) then
!!        open(4132,file="timing/Actions.time.dat",status="unknown")
        open(4132,file=timingdir(1:getlen(timingdir)-1)//"/Actions.time.dat",status="unknown")
        write(4132,'(500A15)') (action_list(actions(i)),i=1,numactions)
        close(4132)
     endif
     if (mod(calledhere,10).eq.1) then


!!        open(4132,file="timing/Actions.time.dat",status="old",position="append")
        open(4132,file=timingdir(1:getlen(timingdir)-1)//"/Actions.time.dat",status="old",position="append")
        write(4132,'(500I15)') (times(actions(i))/1000,i=1,numactions)
        close(4132)
     endif
  endif
  
end subroutine actionsub


subroutine write_actions()
  use parameters
  use actionlistmod
  implicit none
  integer :: i,j

  do i=1,numactions
     if ((actions(i).gt.25).or.(actions(i).lt.1)) then
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
  do i=1,numactions
     if ((actions(i).eq.15).or.(actions(i).eq.16).or.(actions(i).eq.17)) then
        j=1;        exit
     endif
  enddo
  if (j==1) then
     write(mpifileptr,*) "FOR IONIZATION FLUX:"
     write(mpifileptr,*) "   Flux Interval ",FluxInterval,FluxInterval*par_timestep
  endif

  call closefile()

end subroutine write_actions


subroutine actions_final()
  use parameters
  use mpimod
  use xxxmod
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
          call fluxwrite(numpropsteps/FluxInterval,yyy%cmfpsivec(spfstart,0),yyy%cmfpsivec(astart(1),0))
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
  implicit none
  integer :: i, dflag

  dflag=0

  do i=1,numactions
     select case (actions(i))
     case (1)    !! auto-correlation
        call autocorrelate_initial()
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
        if (dfrestrictflag.eq.0) then
           OFLWR "Um, you are doing action 22 without dfrestrictflag...wtf..."; CFLST
        endif
        if (dfrestrictflag.lt.2) then
           OFLWR "WARNING, usually need dfrestrictflag=2 or greater for action 22...."; CFL
        endif
        OFLWR "Getting D-F info for act=22"; CFL
        call getdfcon()
     case (23)
        call wfnovl()
     case(24)
     case(25)
     case default
        OFLWR "Action not supported: ", actions(i); CFLST
     end select
  enddo
end subroutine actions_initial


subroutine action_replacenat()
implicit none

  call replace_withnat(1)
  call get_stuff(0d0)

end subroutine action_replacenat



subroutine fixphase0(invector,isize)
  use parameters
  implicit none
  integer :: isize,ii
!  integer, save :: icount=0
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


!  icount=icount+1
!  if (icount.lt.10) then
!     OFLWR "FIXPHASE DEPRECATED"; CFL
!  endif

end subroutine fixphase0

