

!! SUBROUTINES FOR IMPROVEDQUADFLAG, AVECTOR AND ORBITAL.  CALLS DGMRES ROUTINE.

#include "Definitions.INC"

!      SUBROUTINE DGMRES (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,
!     +   ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX, RGWK, LRGW,
!     +   IGWK, LIGW, RWORK, IWORK)

!! ON INPUT SOLUTION IS GUESS.

subroutine dgsolve0(rhs, solution, numiter, inmult, preconflag, inprecon, &
     intolerance, indimension, inkrydim,parflag,ierr)
  use fileptrmod
  implicit none
  integer,intent(in) :: indimension,inkrydim,parflag,preconflag
  integer,intent(out) :: ierr
  real*8,intent(in) :: intolerance
  DATATYPE,intent(in) :: rhs(indimension)
  integer,intent(out) :: numiter
  DATATYPE,intent(out) :: solution(indimension)
  integer :: itol, iunit, jjxx, rgwkdim
  external :: inmult, dummysub, inprecon
  integer :: igwk(20), ligwk=20, nullint1,nullint2,nullint3,nullint4, nullint5
  real*8 :: nulldouble1, errest, nulldouble2, nulldouble3, nulldouble4, tol
  real*8, allocatable :: rgwk(:)
#ifdef REALGO
  integer, parameter :: zzz=1
#else
  integer, parameter :: zzz=2
#endif
  integer, parameter :: jacwrite=0

!!$  looks ok now dgmres bugfix 10-2015
!!$  mindim=indimension; maxdim=indimension
!!$  call mympiimax(maxdim); call mympiimin(mindim)
!!$  if (mindim.ne.maxdim) then
!!$     OFLWR "Error, for dgsolve now must use same dimension all processors", &
!!$         mindim,maxdim; CFLST
!!$  endif

  jjxx = zzz * indimension

  rgwkdim= (42 + jjxx*(inkrydim*zzz + 6) + inkrydim*zzz*(inkrydim*zzz + 3)) / 5 * 6

  if (rgwkdim.lt.0) then
     OFLWR "OOPS, integer overflow, decrease maxexpodim or reduce size of orbitals (per processor)",&
          rgwkdim; CFLST
  endif

  allocate(rgwk(rgwkdim)); rgwk(:)=0d0

  itol=0;  iunit=0

  igwk=0
  igwk(1) = inkrydim*zzz    !! max krylov dim
  igwk(2) = inkrydim*zzz    !! max orthog
  igwk(3) = 0             !! scaling (nulldouble 2 and 3)
  igwk(4) = preconflag  !!  preconditioning
  igwk(5)=100  !! max restarts
  igwk(5)=-1 

  tol=intolerance

  if (parflag.eq.0) then
     call dgmres(jjxx, rhs, solution, nullint1, nullint2, nullint3, nulldouble1, nullint3, inmult, &
          inprecon,itol , tol, nullint4, numiter, errest, ierr, iunit, nulldouble2, nulldouble3, rgwk, &
          rgwkdim, igwk, ligwk, nulldouble4, nullint5)
  else
     call dgmrespar(jjxx, rhs, solution, nullint1, nullint2, nullint3, nulldouble1, nullint3, inmult, &
          inprecon,itol , tol, nullint4, numiter, errest, ierr, iunit, nulldouble2, nulldouble3, rgwk, &
          rgwkdim, igwk, ligwk, nulldouble4, nullint5)
  endif

!!$  if (debugflag.ne.0) then
!!$     OFLWR "     DGSOLVE - restarts,iters,errest ",igwk(5), numiter,errest; CFL
!!$  endif

  if (jacwrite==1) then
     OFLWR "   NEWTON: numiter ", numiter, " error ", errest;     call closefile()
  endif

!! AUTHOR COMMENTS ARE WRONG IN DGMRES.F ierr=1 denotes max restarts (zero) reached !! 
  if (ierr.eq.1) then 
     OFLWR "WARNING! dgmres did not converge. increase krylov order?",inkrydim; CFL
  else if (ierr.eq.2) then
     OFLWR "WARNING! dgmres stalled."; CFL
  else  if (ierr/=0) then
     OFLWR "Error dgmres ", ierr
     WRFL "igwk(6) = ",igwk(6)
     WRFL "rgwkdim = ",rgwkdim
     CFLST
  endif

  deallocate(rgwk)

end subroutine dgsolve0


!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NEWTON SOLVE (IMPROVEDQUADFLAG=1,2) FOR ORBITALS
!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine quadinit(inspfs, thistime)
  use parameters
  implicit none
  DATATYPE,intent(in) :: inspfs(spfsize,nspf) 
  real*8,intent(in) :: thistime
  call jacinit0(0,inspfs,thistime)
end subroutine quadinit


subroutine quadoperate(notusedint,inspfs,outspfs)
   use parameters
   implicit none
   integer :: notusedint
   DATATYPE,intent(in) ::  inspfs(totspfdim)
   DATATYPE,intent(out) :: outspfs(totspfdim)
   DATATYPE :: workspfs(totspfdim)              !! AUTOMATIC

   if (jacprojorth.ne.0) then   
      call jacorth(inspfs,outspfs)
      call jacoperate0(0,0,outspfs,workspfs)
      call jacorth(workspfs,outspfs)
   else
      call jacoperate0(0,0,inspfs,outspfs)
   endif

end subroutine quadoperate


subroutine quadopcompact(notusedint,com_inspfs,com_outspfs)
   use parameters
   implicit none
   integer :: notusedint
   DATATYPE,intent(in) ::  com_inspfs(spfsmallsize*nspf)
   DATATYPE,intent(out) :: com_outspfs(spfsmallsize*nspf)
   DATATYPE :: inspfs(spfsize*nspf),outspfs(spfsize*nspf)  !! AUTOMATIC

   call spfs_expand(com_inspfs,inspfs)
   if (jacprojorth.ne.0) then   
      call jacorth(inspfs,outspfs)
      call jacoperate0(0,0,outspfs,inspfs)
      call jacorth(inspfs,outspfs)
   else
      call jacoperate0(0,0,inspfs,outspfs)
   endif
   call spfs_compact(outspfs,com_outspfs)

end subroutine quadopcompact


subroutine quadspfs(inspfs,jjcalls)
  use parameters
  use mpimod
  use linearmod
  implicit none
  DATATYPE,intent(inout) :: inspfs(totspfdim)  
  integer,intent(out) :: jjcalls
  integer :: icount,maxdim,ierr
  real*8 :: orthogerror,dev,mynorm
  DATATYPE, allocatable ::  vector(:), vector2(:), vector3(:), com_vector2(:), com_vector3(:)
  EXTERNAL :: quadoperate,dummysub,quadopcompact

  allocate( vector(totspfdim), vector2(totspfdim), vector3(totspfdim) )
  vector=0; vector2=0; vector3=0
  if (orbcompact.ne.0) then
     allocate( com_vector2(spfsmallsize*nspf), com_vector3(spfsmallsize*nspf) )
     com_vector2=0; com_vector3=0
  endif

  if (jacsymflag.ne.1.and.jacprojorth.eq.0) then
     OFLWR "setting jacsymflag=1 for orbital quad"; CFL
     jacsymflag=1
  endif

  effective_cmf_linearflag=0

  vector=inspfs

  do icount=0,0

     call spf_orthogit(vector,orthogerror)

     call quadinit(vector,0d0)

     call spf_linear_derivs0(0,0,0d0,vector,vector2,1,1)

     call apply_spf_constraints(vector2)

     dev=abs(hermdot(vector2,vector2,totspfdim))
     if (parorbsplit.eq.3) then
        call mympirealreduceone(dev)
     endif
     dev=sqrt(dev)
  
     OFLWR "   Quad orbitals: Deviation is     ", dev; CFL

     vector3=vector   !! guess

     if (orbcompact.ne.0) then
        call spfs_compact(vector2,com_vector2)
        call spfs_compact(vector3,com_vector3)
        if (parorbsplit.eq.3) then
           maxdim=min(spfsmallsize*nspf*nprocs,maxexpodim)
           call dgsolve0( com_vector2, com_vector3, jjcalls, quadopcompact,0,dummysub, &
                quadtol,spfsmallsize*nspf,maxdim,1,ierr)
        else
           maxdim=min(spfsmallsize*nspf,maxexpodim)
           call dgsolve0( com_vector2, com_vector3, jjcalls, quadopcompact,0,dummysub, &
                quadtol,spfsmallsize*nspf,maxdim,0,ierr)
        endif
        call spfs_expand(com_vector3,vector3)
     else
        if (parorbsplit.eq.3) then
           maxdim=min(totspfdim*nprocs,maxexpodim)
           call dgsolve0( vector2, vector3, jjcalls, quadoperate,0,dummysub, &
                quadtol,totspfdim,maxdim,1,ierr)
        else
           maxdim=min(totspfdim,maxexpodim)
           call dgsolve0( vector2, vector3, jjcalls, quadoperate,0,dummysub, &
                quadtol,totspfdim,maxdim,0,ierr)
        endif
     endif

     if (ierr.eq.1) then
        vector3(:)=0d0
     endif

     mynorm=abs(hermdot(vector3,vector3,totspfdim))
     
     if (parorbsplit.eq.3) then
        call mympirealreduceone(mynorm)
     endif
     mynorm=sqrt(mynorm)

     if (mynorm.gt.maxquadnorm*nspf) then
        vector3=vector3*maxquadnorm*nspf/mynorm
     endif

     vector=vector+vector3

  enddo

  call spf_orthogit(vector,orthogerror)

  inspfs = vector

  call spf_linear_derivs0(0,0,0d0,vector,vector2,1,1)

  call apply_spf_constraints(vector2)

  dev=abs(hermdot(vector2,vector2,totspfdim))
  if (parorbsplit.eq.3) then
     call mympirealreduceone(dev)
  endif
  dev=sqrt(dev)
  
  OFLWR "   Quad orbitals: Deviation is now ", dev;
  WRFL  "                  Orthog error is  ", orthogerror; CFL

  deallocate( vector,vector2,vector3)
  
  if (orbcompact.ne.0) then
     deallocate(com_vector2,com_vector3)
  endif

end subroutine quadspfs


!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NEWTON SOLVE (IMPROVEDQUADFLAG=1,3) FOR AVECTOR
!!!!!!!!!!!!!!!!!!!!!!!!!!


module aaonedmod
  implicit none
  DATATYPE :: quadexpect=0d0
end module



subroutine aaonedinit(www,inavector) 
  use aaonedmod
  use xxxmod
  use r_parameters
  use walkmod
  use dotmod
  implicit none
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: inavector(www%totadim)
  DATATYPE,allocatable :: jacaa(:),jacaamult(:)
  DATATYPE :: csum2

  allocate(jacaa(www%totadim), jacaamult(www%totadim))
  if (www%totadim.gt.0) then
     jacaamult(:)=0
     jacaa(:)=inavector(:)
  endif
  call basis_project(www,numr,jacaa)

  csum2=0d0
  if (www%totadim.gt.0) then
     csum2=dot(jacaa,jacaa,www%totadim)
  endif
  if (www%parconsplit.ne.0) then
     call mympireduceone(csum2)
  endif
  if (www%totadim.gt.0) then
     jacaa(:)=jacaa(:)/sqrt(csum2)
  endif

!! 01-2016 setting conflag to zero here too.  Tau not used in avector quad equations.
!!   here (defining quadexpect) perhaps it could be set to 1
  call sparseconfigmult(www,jacaa, jacaamult, yyy%cptr(0), yyy%sptr(0),1,1,1,0,0d0,-1)

  call basis_project(www,numr,jacaamult)

  quadexpect=0
  if (www%totadim.gt.0) then
     quadexpect=dot(jacaa,jacaamult,www%totadim)
  endif
  if (www%parconsplit.ne.0) then
     call mympireduceone(quadexpect)
  endif

  deallocate(jacaa,jacaamult)

end subroutine aaonedinit



subroutine quadavector(inavector,jjcalls)
  use parameters
  use configmod
  implicit none
  DATATYPE,intent(inout) :: inavector(www%totadim,mcscfnum)
  DATATYPE :: nullvector(numr)
  integer,intent(out) :: jjcalls
  integer :: imc

  if (quadorthflag.ne.0) then
     do imc=1,mcscfnum
        if (www%parconsplit.eq.0) then
           call gramschmidt(www%totadim,imc-1,www%totadim,&
                inavector(:,:),inavector(:,imc),.false.)
        else
           if (www%totadim.gt.0) then
              call gramschmidt(www%totadim,imc-1,www%totadim,&
                   inavector(:,:),inavector(:,imc),.true.)
           else
              call nullgramschmidt(imc-1,.true.)
           endif
        endif
     enddo
  endif

  do imc=1,mcscfnum
     if (sparseconfigflag==0) then
        call nonsparsequadavector(www,inavector(:,imc))
        jjcalls=0
     else
        if (www%totadim.gt.0) then
           call sparsequadavector(www,inavector(:,imc),jjcalls)
        else
           call sparsequadavector(www,nullvector(:),jjcalls)
        endif
     endif
  enddo

end subroutine quadavector


subroutine sparsequadavector(www,inavector,jjcalls0)
  use parameters
  use mpimod
  use aaonedmod
  use xxxmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  EXTERNAL :: paraamult, parquadpreconsub
  DATATYPE, intent(inout) ::  inavector(numr,www%firstconfig:www%lastconfig)
  integer,intent(out) :: jjcalls0
  integer :: jjcalls, ss, maxdim, mysize,flag,ierr
  real*8 ::  dev,  thisaerror
  DATATYPE :: csum 
  DATATYPE, allocatable ::  vector(:,:), vector2(:,:), vector3(:,:), &
       smallvectorspin(:,:),smallvectorspin2(:,:)

  if (sparseconfigflag.eq.0) then
     OFLWR "Error, must use sparseconfigflag.ne.0 for sparsequadavector"; CFLST
  endif

  allocate(smallvectorspin(numr,www%botdfbasis:www%topdfbasis+1), &
       smallvectorspin2(numr,www%botdfbasis:www%topdfbasis+1))
  smallvectorspin=0; smallvectorspin2=0
  
  allocate( vector(numr,www%firstconfig:www%lastconfig+1), &
       vector2(numr,www%firstconfig:www%lastconfig+1), &
       vector3(numr,www%firstconfig:www%lastconfig+1))
  vector=0; vector2=0; vector3=0

  if (www%lastconfig.ge.www%firstconfig) then
     vector(:,www%firstconfig:www%lastconfig)=inavector(:,:)
  endif

  jjcalls0=0
  ss=0

  do while (42.eq.42)

     call basis_project(www,numr,vector)

     call aaonedinit(www,vector)

!! CONSISTENT WITH PARAAMULT/PARBLOCKCONFIGMULT yes pulse no tau
     call sparseconfigmult(www,vector, vector2, yyy%cptr(0), yyy%sptr(0), &
          1,1,1,0,0d0,-1)

     call basis_project(www,numr,vector2)

     vector3=vector2-quadexpect*vector              !! error term.

     dev=0
     if (www%totadim.gt.0) then
        dev=abs(hermdot(vector3,vector3,www%totadim))
     endif
     if (www%parconsplit.ne.0) then
        call mympirealreduceone(dev)
     endif
     dev=sqrt(dev)

!! 12-2015 sqrt
     thisaerror=min(0.1d0,sqrt(aerror/min(dev,1d0)))

     OFL;write(mpifileptr,'(A20,E12.5,A6,2E12.5,A7,100F14.7)') &
          "   SPARSEQUAD: DEV", dev, " TOL ",aerror,thisaerror,"ENERGY",quadexpect; CFL

     if (abs(dev).lt.aerror.or.ss.ge.10) then
        if (www%lastconfig.ge.www%firstconfig) then
           inavector(:,:) = vector(:,www%firstconfig:www%lastconfig)
        endif
        if (ss.ge.10) then
           OFLWR "10 iterations, quad aborting"; CFL
        endif
        deallocate(smallvectorspin,smallvectorspin2)
        deallocate(vector, vector2, vector3)
        call mpibarrier()
        return  !! RETURN
     endif

     smallvectorspin(:,:)=0d0
     if (www%topconfig.ge.www%botconfig) then
        call basis_transformto_local(www,numr,vector(:,www%botconfig:www%topconfig),smallvectorspin)
     endif
     smallvectorspin2(:,:)=smallvectorspin(:,:)    !! guess

     maxdim=min(maxaorder,numr*www%numdfbasis)
     mysize=numr*(www%topdfbasis-www%botdfbasis+1)

     flag=0
     if (mysize.eq.0) then
        print *, "ACK, CAN'T DO A-VECTOR QUAD WITH ZERO CONFIGS PER PROCESSOR RANK",myrank,mysize; 
        flag=99
     endif
     call mympiireduceone(flag)
     if (flag.ne.0) then
        call mpistop()
     endif

     call dgsolve0( smallvectorspin, smallvectorspin2, jjcalls, paraamult,&
          quadprecon,parquadpreconsub, thisaerror,mysize,maxdim,1,ierr)

!!$  call basicblocklansolve(1,mysize,maxdim,maxdim,smallvectorspin,smallvectorspin2,1,&
!!$     parhrmult,parhrdotsub,quadexpect)
!!$  jjcalls=0

     vector3(:,:)=0d0; 

     if (www%topconfig.ge.www%botconfig) then
        call basis_transformfrom_local(www,numr,smallvectorspin2,&
             vector3(:,www%botconfig:www%topconfig))
     endif
     if (www%parconsplit.eq.0) then
        call mpiallgather(vector3(:,:),www%numconfig*numr,&
             www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
     endif

     csum=0d0
     if (www%totadim.gt.0) then
        csum=dot(vector3,vector3,www%totadim)
     endif
     if (www%parconsplit.ne.0) then
        call mympireduceone(csum)
     endif
     vector=vector3/sqrt(csum)

!!  OFLWR "    ITERATIONS ***  ",jjcalls

     jjcalls0=jjcalls0+jjcalls

     ss=ss+1

  enddo  !! INFINITE LOOP

end subroutine sparsequadavector


subroutine nonsparsequadavector(www,avectorout)
  use mpimod
  use aaonedmod
  use parameters
  use xxxmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  DATATYPE,intent(inout) :: avectorout(www%totadim)
  DATATYPE :: csum
  CNORMTYPE :: norm
  real*8 :: dev
  integer :: k, info,ss
  DATATYPE,allocatable ::  spinmatel(:,:),spinavectorout(:),err(:)
  integer,allocatable :: ipiv(:)

  if (sparseconfigflag/=0) then
     OFLWR "Error, nonsparsequadavector called but sparseconfigflag= ", sparseconfigflag
     CFLST
  endif

  allocate(spinmatel(www%numdfbasis*numr,www%numdfbasis*numr),&
       spinavectorout(www%numdfbasis*numr), err(www%totadim),ipiv(www%numdfbasis*numr))
  spinmatel=0; spinavectorout=0; ipiv=0;     err=0; 

  ss=0
  do while (42.eq.42)

     call basis_project(www,numr,avectorout)

     call aaonedinit(www,avectorout)

!! CONSISTENT WITH PARAAMULT/PARBLOCKCONFIGMULT yes pulse no tau
     call sparseconfigmult(www,avectorout(:),err(:),yyy%cptr(0),yyy%sptr(0),1,1,1,0,0d0,-1)

     call basis_project(www,numr,err(:))

     err(:)=err(:)-quadexpect*avectorout(:)              !! error term.

     dev=sqrt(abs(hermdot(err(:),err(:),www%totadim)))

     OFL;write(mpifileptr,'(A22,E12.5,A6,E12.5,A7,100F14.7)') &
          "   NONSPARSEQUAD: DEV", dev, " TOL ",aerror,"ENERGY",quadexpect; CFL

     if (abs(dev).lt.aerror.or.(ss.gt.10)) then

        deallocate(spinmatel,spinavectorout,err,ipiv)

!#ifdef REALGO
!  OFL; write(mpifileptr,'(A40,F15.10, A5,1F18.10,A5,1E6.1)') &
!      "      Converged nonsparse quad: dev ", dev, "  E= ", quadexpect, " tol ", aerror; CFL
!#else
!  OFL; write(mpifileptr,'(A40,F15.10, A5,2F18.10,A5,1E6.1)') &
!      "      Converged nonsparse quad: dev ", dev, "  E= ", quadexpect, " tol ", aerror; CFL
!#endif

        call mpibarrier()
        return   !! RETURN
     endif

     spinmatel(:,:)=0d0

!! CONSISTENT WITH PARAAMULT/PARBLOCKCONFIGMULT yes pulse no tau
!! no constraint (tau) in eigenvalue equation
     call assemble_dfbasismat(www, spinmatel, yyy%cptr(0),1,1,1,0, 0d0,-1)

     do k=1,www%numdfbasis*numr
        spinmatel(k,k)=spinmatel(k,k)-quadexpect
     enddo

     spinavectorout(:)=0d0

     call basis_transformto_all(www,numr,avectorout(:),spinavectorout(:))


     if (myrank.eq.1) then
        call MYGESV(www%numdfbasis*numr,1,spinmatel,www%numdfbasis*numr,&
             ipiv,spinavectorout,www%numdfbasis*numr,info)
        spinavectorout=spinavectorout/&
             sqrt(dot(spinavectorout,spinavectorout,www%numdfbasis*numr))
        if (info/=0) then
           OFLWR "Info in dgesv nonsparsequadavector= ", info; CFLST
        endif
     endif

     call mpibarrier()

     call mympibcast(spinavectorout,1,www%numdfbasis*numr)

     call basis_transformfrom_all(www,numr,spinavectorout(:),avectorout(:))

     csum=0
     if (www%totadim.gt.0) then
        csum=dot(avectorout,avectorout,www%totadim)  !ok conv
     endif
     if (www%parconsplit.ne.0) then
        call mympireduceone(csum)
     endif
     norm=sqrt(csum)  !! ok conversion
     if (www%totadim.gt.0) then
        avectorout=avectorout/norm
     endif

     ss=ss+1

  enddo  !! INFINITE LOOOP

end subroutine nonsparsequadavector



subroutine paraamult(notusedint,inavectorspin,outavectorspin)
  use r_parameters
  use aaonedmod
  use configmod
  implicit none
  integer,intent(in) :: notusedint
  DATATYPE,intent(in) :: inavectorspin(numr,www%botdfbasis:www%topdfbasis)
  DATATYPE,intent(out) :: outavectorspin(numr,www%botdfbasis:www%topdfbasis)

  call parblockconfigmult(inavectorspin,outavectorspin)

  if (www%topdfbasis.ge.www%botdfbasis) then
     outavectorspin(:,:)= outavectorspin(:,:) - quadexpect*inavectorspin(:,:)
  endif

end subroutine paraamult



subroutine parquadpreconsub(notusedint, inavectorspin,outavectorspin)
  use fileptrmod
  use r_parameters
  use sparse_parameters
  use aaonedmod
  use xxxmod
  use configmod
  implicit none
  integer,intent(in) :: notusedint
  DATATYPE,intent(in) :: inavectorspin(numr,www%botdfbasis:www%topdfbasis)
  DATATYPE,intent(out) :: outavectorspin(numr,www%botdfbasis:www%topdfbasis)

  if (use_dfwalktype) then
     call parquadpreconsub0(dwwptr,yyy%cptr(0),yyy%sdfptr(0),&
          notusedint,inavectorspin,outavectorspin)
  else
     call parquadpreconsub0(dwwptr,yyy%cptr(0),yyy%sptr(0),&
          notusedint,inavectorspin,outavectorspin)
  endif

end subroutine parquadpreconsub


subroutine parquadpreconsub0(www,cptr,sptr,notusedint, inavectorspin,outavectorspin)
  use fileptrmod
  use r_parameters
  use sparse_parameters
  use aaonedmod
  use walkmod
  use configptrmod
  use sparseptrmod
  implicit none
  type(walktype),intent(in) :: www
  type(CONFIGPTR),intent(in) :: cptr
  type(SPARSEPTR),intent(in) :: sptr
  integer,intent(in) :: notusedint
  DATATYPE,intent(in) :: inavectorspin(numr,www%botdfbasis:www%topdfbasis)
  DATATYPE,intent(out) :: outavectorspin(numr,www%botdfbasis:www%topdfbasis)
  DATATYPE :: inavector(numr,www%botconfig:www%topconfig+1),&
       outavector(numr,www%botconfig:www%topconfig+1)  !!AUTOMATIC

  if (sparseconfigflag.eq.0) then
     OFLWR "error, no parquadpreconsub if sparseconfigflag=0"; CFLST
  endif

  call basis_transformfrom_local(www,numr,inavectorspin,inavector)

  call parsparseconfigpreconmult(www,inavector, outavector, cptr, sptr,&
       1,1,1,1, quadexpect,0d0)

  call basis_transformto_local(www,numr,outavector,outavectorspin)

end subroutine parquadpreconsub0



!!$
!!$subroutine parhrmult(lanblocknum,inavectorspin,outavectorspin)
!!$  use parameters
!!$  use aaonedmod
!!$  implicit none
!!$  integer,intent(in) :: lanblocknum
!!$  DATATYPE,intent(in) :: inavectorspin(numr,botbasis:topbasis)
!!$  DATATYPE,intent(out) :: outavectorspin(numr,botbasis:topbasis)
!!$
!!$  if (lanblocknum.ne.1) then
!!$     OFLWR "DOOGSNATCH",lanblocknum; CFLST
!!$  endif
!!$  call parblockconfigmult(inavectorspin,outavectorspin)
!!$
!!$end subroutine parhrmult
!!$
!!$subroutine parhrdotsub(bra,ket,size,numbra,numket,outdots)
!!$  implicit none
!!$  integer,intent(in) :: size,numbra,numket
!!$  DATATYPE,intent(in) :: bra(size,numbra),ket(size,numbra)
!!$  DATATYPE,intent(out) :: outdots(numbra,numket)
!!$  DATATYPE :: dot
!!$  integer :: ii,jj
!!$  do jj=1,numket
!!$     do ii=1,numbra
!!$        outdots(ii,jj)=dot(bra(:,ii),ket(:,jj),size)
!!$     enddo
!!$  enddo
!!$end subroutine parhrdotsub
!!$

