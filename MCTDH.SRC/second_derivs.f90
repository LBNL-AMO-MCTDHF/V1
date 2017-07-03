
!! ALL MODULES

!! FOR VERLET AND IMPLICIT PROPAGATION (EXPERIMENTAL)

#include "Definitions.INC"

!! derspfs0 is the first derivative; sdspfs is second derivative.

!! jacunitflag, jacsymflag not implementd

module secondsubmod
contains

subroutine second_derivs00(lowspf,highspf,thistime,inspfs,sdspfs)
  use parameters
  use linearmod
  use derivativemod
  use orbprojectmod
  use orbgathersubmod
  implicit none
  integer,intent(in) :: lowspf,highspf
  real*8,intent(in) :: thistime
  DATATYPE,intent(in) :: inspfs(spfsize,nspf)
  DATATYPE,intent(out) :: sdspfs(spfsize,lowspf:highspf)
  DATATYPE :: nullspfs(1)
  DATATYPE :: workoutspfs0(spfsize,lowspf:highspf+1), &       !! AUTOMATIC
          workspfs0(spfsize,lowspf:highspf+1), workspfs2(spfsize,lowspf:highspf+1),&
          workoutspfs2(spfsize,lowspf:highspf+1),  derspfs0(spfsize,nspf)
  real*8 :: gridtime
  integer :: numspf

  if (numfrozen.gt.0) then
     OFLWR "second_derivs not done for frozen orbitals",numfrozen; CFLST
  endif
  if (drivingflag.ne.0) then
     OFLWR "second_derivs not done for drivingflag"; CFLST
  endif

  numspf=highspf-lowspf+1
  if (numspf.le.0) then
     call waitawhile()
     print *, "ACKDOOG 5689"
     call waitawhile()
     stop
  endif

  gridtime=(thistime-firsttime)/(lasttime-firsttime) 
  if (jacprojorth.ne.0) then
     OFLWR "program jacprojorth secondderivs"; CFLST
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

!! First derivative.

  call actreduced00(lowspf,highspf,1,thistime, inspfs, nullspfs, workspfs0, 0, 0,1)

  if (effective_cmf_linearflag.eq.0) then
     call project00(lowspf,highspf,workspfs0(:,lowspf:highspf), &
          derspfs0(:,lowspf:highspf), inspfs)
     derspfs0(:,lowspf:highspf)=workspfs0(:,lowspf:highspf)-derspfs0(:,lowspf:highspf)
  else
     call actreduced00(lowspf,highspf,1,thistime, inspfs, nullspfs, workspfs2, 1, 0,1)
     workoutspfs0=(1.d0-gridtime)*workspfs2 + gridtime*workspfs0
     call project00(lowspf,highspf,workoutspfs0(:,lowspf:highspf), &
          derspfs0(:,lowspf:highspf), inspfs)
     derspfs0(:,lowspf:highspf)=workoutspfs0(:,lowspf:highspf)-derspfs0(:,lowspf:highspf)
  endif

  if (parorbsplit.eq.1) then
     call mpiorbgather_nz(derspfs0,spfsize)
  endif

!! d/dt phi = (1-P) Q phi

!! d^2/dt^2 phi = [ d/dt (1-P) ] Q phi + (1-P) Q [d/dt phi] for stype 1

 !! second term    

  call actreduced00(lowspf,highspf,1,thistime, derspfs0, nullspfs, workoutspfs0, 0, 0,1) 

  if (effective_cmf_linearflag.ne.0) then

     call actreduced00(lowspf,highspf,1,thistime, derspfs0, nullspfs, workoutspfs2, 1, 0,1)
     workoutspfs0 = (1-gridtime)*workoutspfs2 + gridtime*workoutspfs0

!! d/dt Q term
     workoutspfs0 = workoutspfs0 + 1.d0/(lasttime-firsttime) * (workspfs0-workspfs2)

  endif

  call project00(lowspf,highspf,workoutspfs0(:,lowspf:highspf), &
       sdspfs(:,lowspf:highspf), inspfs)
  sdspfs(:,lowspf:highspf)=workoutspfs0(:,lowspf:highspf)-sdspfs(:,lowspf:highspf)

!! first term.

  if (effective_cmf_linearflag.ne.0.and.numspf.gt.0) then
     workspfs0 = (1-gridtime)*workspfs2 + gridtime*workspfs0
  endif

!!   Projector is just sum_i |phi_i><phi_i| without regard to orthonorm 
!!      (consistent with noorthogflag=1)

  call derproject00(lowspf,highspf,workspfs0, workoutspfs0, derspfs0, inspfs)
  sdspfs(:,lowspf:highspf)=sdspfs(:,lowspf:highspf)-workoutspfs0(:,lowspf:highspf)

end subroutine second_derivs00


subroutine second_derivs(thistime,inspfs,sdspfs)
  use parameters
  use orbgathersubmod
  implicit none
  real*8,intent(in) :: thistime
  DATATYPE,intent(in) :: inspfs(spfsize,nspf)
  DATATYPE,intent(out) :: sdspfs(spfsize,nspf)
  integer :: lowspf,highspf

  lowspf=1; highspf=nspf
  if (parorbsplit.eq.1) then
     call getOrbSetRange(lowspf,highspf)
  endif

  if (highspf.ge.lowspf) then
     call second_derivs00(lowspf,highspf,thistime,inspfs,sdspfs(:,lowspf:highspf))
  endif
  if (parorbsplit.eq.1) then
     call mpiorbgather(sdspfs,spfsize)
  endif

end subroutine second_derivs

end module


module verletmod
  implicit none
  integer :: jstep=0
  DATATYPE, allocatable :: prevspfs(:,:)
  
contains
  
subroutine verlet(inspfs0, time1,time2,outiter)
  use parameters
  use orbgathersubmod
  use mpi_orbsetmod
  use expospfpropmod
  use secondsubmod
  use spfsubmod
  implicit none
  real*8,intent(in) :: time1,time2
  integer,intent(out) :: outiter
  DATATYPE,intent(inout) :: inspfs0(spfsize,nspf)
  DATATYPE,allocatable :: inspfs(:,:), tempspfs(:,:),  sdspfs(:,:)
  real*8 :: thistime, cstep, orthogerror
  integer :: istep,lowspf,highspf,numspf

  lowspf=1; highspf=nspf
  if (parorbsplit.eq.1) then
     call getOrbSetRange(lowspf,highspf)
  endif
  numspf=highspf-lowspf+1

  allocate(inspfs(spfsize,nspf), tempspfs(spfsize,nspf), &
       sdspfs(spfsize,firstmpiorb:firstmpiorb+orbsperproc-1))
  tempspfs=0; sdspfs=0
  inspfs(:,:)=inspfs0(:,:)

  outiter=verletnum

  cstep = (time2-time1)/real(verletnum,8)

  do istep=1,verletnum

     jstep=jstep+1
     
     thistime=(istep-1)*cstep + time1
     if (numspf.gt.0) then
        call second_derivs00(lowspf,highspf,thistime,inspfs,sdspfs)
     endif
     if (jstep.eq.1) then
        allocate(prevspfs(spfsize,nspf))
     endif

     if (istep.eq.1) then
        tempspfs(:,:)=inspfs0(:,:)
        call expospfprop(time1,time2,tempspfs,outiter)
     else
        if (numspf.gt.0) then
           tempspfs(:,lowspf:highspf) = 2*inspfs(:,lowspf:highspf) - prevspfs(:,lowspf:highspf) &
                + cstep**2 * sdspfs(:,lowspf:highspf)
        endif
        if (parorbsplit.eq.1) then
           call mpiorbgather(tempspfs,spfsize)
        endif
     endif

     prevspfs=inspfs
     inspfs=tempspfs

!     if (orthogerror.gt.1d-6) then
!        OFLWR "ORTHOGERROR VERLET : ", orthogerror; CFLST
!     endif
  end do

  call spf_orthogit(inspfs,orthogerror)

  inspfs0(:,:)=inspfs(:,:)

  deallocate(inspfs,tempspfs,sdspfs)
  
  OFLWR "final ORTHOGERROR VERLET : ", orthogerror; CFL

end subroutine verlet

end module verletmod
