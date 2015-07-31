
!! FOR VERLET AND IMPLICIT PROPAGATION (EXPERIMENTAL)

#include "Definitions.INC"

!! derspfs0 is the first derivative; sdspfs is second derivative.

!! jacunitflag, jacsymflag not implementd

subroutine second_derivs(thistime,inspfs,sdspfs)
  use parameters
  use linearmod
  implicit none
  real*8 :: thistime, gridtime
  DATATYPE :: inspfs(spfsize,nspf), sdspfs(spfsize,nspf), nullspfs(1)
  integer,save :: allocated=0
  DATATYPE, save, allocatable :: workoutspfs0(:,:), workoutspfs2(:,:), workspfs2(:,:), &
       workspfs0(:,:), derspfs0(:,:)

  gridtime=(thistime-firsttime)/(lasttime-firsttime) 
  if (jacprojorth.ne.0) then
     OFLWR "program jacprojorth secondderivs???"; CFLST
  endif
  if (jacsymflag.ne.0) then
     OFLWR "program jacsym secondderivs"; CFLST
  endif
  if (jacgmatthird.ne.0) then
     OFLWR "program jacgmatthird secondderivs"; CFLST
  endif
  if (constraintflag.ne.0) then
     OFLWR "constraintflag with verlet not checked"; CFLST
  endif
  if (allocated==0) then
     allocate(workoutspfs0(spfsize,nspf),workoutspfs2(spfsize,nspf), &
          workspfs0(spfsize,nspf), workspfs2(spfsize,nspf),derspfs0(spfsize,nspf))
     allocated=1
  endif

!! First derivative.

  call actreduced0(thistime, inspfs, nullspfs, workspfs0, 1, 0,1)
  if (effective_cmf_linearflag.eq.0) then
     call oneminusproject(workspfs0, derspfs0, inspfs)
  else
     call actreduced0(thistime, inspfs, nullspfs, workspfs2, 0, 0,1)
     workoutspfs0=(1.d0-gridtime)*workspfs0 + gridtime*workspfs2
     call oneminusproject(workoutspfs0, derspfs0, inspfs)
  endif

!! d/dt phi = (1-P) Q phi

!! d^2/dt^2 phi = [ d/dt (1-P) ] Q phi + (1-P) Q [d/dt phi] for stype 1

 !! second term    

  call actreduced0(thistime, derspfs0, nullspfs, workoutspfs0, 1, 0,1) 

  if (effective_cmf_linearflag.eq.0) then
     call actreduced0(thistime, derspfs0, nullspfs, workoutspfs2, 0, 0,1)  !! second term
     workoutspfs0 = (1-gridtime)*workoutspfs0 + gridtime*workoutspfs2

!! d/dt Q term

     workoutspfs0 = workoutspfs0 + 1.d0/(lasttime-firsttime) * (workspfs2-workspfs0)
  endif
  call oneminusproject(workoutspfs0, sdspfs, inspfs)

!! first term.

  if (effective_cmf_linearflag.eq.0) then
     workspfs0 = (1-gridtime)*workspfs0 + gridtime*workspfs2
  endif

  !!   Projector is just sum_i |phi_i><phi_i| without regard to orthonorm (consistent with noorthogflag=1)

  call derproject(workspfs0, workoutspfs0, derspfs0, inspfs)
  sdspfs=sdspfs-workoutspfs0


end subroutine second_derivs


  
subroutine verlet(inspfs0, time1,time2,outiter)
  use parameters
  implicit none

  DATATYPE, allocatable, save :: prevspfs(:,:)
  integer, save :: jstep=0
  DATATYPE :: inspfs(spfsize,nspf), inspfs0(spfsize,nspf)
  DATATYPE :: tempspfs(spfsize,nspf),  sdspfs(spfsize,nspf)
  real*8 :: time1,time2, thistime, cstep, orthogerror
  integer :: istep,outiter
  real*8, external :: spf_linear_derivs

!  OFLWR "GO VERLET!"; CFL

  outiter=verletnum

  inspfs(:,:)=inspfs0(:,:)


  cstep = (time2-time1)/real(verletnum,8)

  do istep=1,verletnum

     jstep=jstep+1
     
     thistime=(istep-1)*cstep + time1
     call second_derivs(thistime,inspfs,sdspfs)

     if (jstep.eq.1) then
        allocate(prevspfs(spfsize,nspf))
     endif

     if (istep.eq.1) then

        tempspfs=inspfs
        call expoprop(time1,time2,tempspfs,outiter)
     else
        tempspfs = 2*inspfs - prevspfs + cstep**2 * sdspfs 
     endif
     prevspfs=inspfs
     inspfs=tempspfs


!     if (orthogerror.gt.1d-6) then
!        OFLWR "ORTHOGERROR VERLET : ", orthogerror; CFLST
!     endif
  end do

  call spf_orthogit(inspfs,orthogerror)

  inspfs0(:,:)=inspfs(:,:)

  OFLWR "final ORTHOGERROR VERLET : ", orthogerror; CFL



end subroutine verlet


