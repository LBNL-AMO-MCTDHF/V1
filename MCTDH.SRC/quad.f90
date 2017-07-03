
!! ALL MODULES

!! SUBROUTINES FOR IMPROVEDQUADFLAG, AVECTOR AND ORBITAL.  CALLS DGMRES ROUTINE.

#include "Definitions.INC"

!      SUBROUTINE DGMRES (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,
!     +   ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX, RGWK, LRGW,
!     +   IGWK, LIGW, RWORK, IWORK)

!! ON INPUT SOLUTION IS GUESS.


module dgsolvemod
contains

  subroutine dgsolve00(rhs, solution, numiter, inmult, preconflag, inprecon, &
       intolerance, indimension, inkrydim,parflag,ierr,incomm)
    use fileptrmod
    use dgm_comm_mod
    use dgmresmod
    implicit none
    integer,intent(in) :: indimension,inkrydim,parflag,preconflag,incomm
    integer,intent(out) :: ierr
    real*8,intent(in) :: intolerance
    DATATYPE,intent(in) :: rhs(indimension)
    integer,intent(out) :: numiter
    DATATYPE,intent(out) :: solution(indimension)
    external inmult, inprecon
    integer :: igwk(20), ligwk=20, nullint1,nullint2,nullint3,nullint4, nullint5, nullint6,&
         itol, iunit, jjxx, rgwkdim
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

    if (indimension.le.0) then
       print *, "ACK, DGSOLVE CAN'T DO ZERO PER PROCESSOR, STOPPING"
       print *, "ACK, DGSOLVE CAN'T DO ZERO PER PROCESSOR, STOPPING"
       print *, "ACK, DGSOLVE CAN'T DO ZERO PER PROCESSOR, STOPPING"
       stop
    endif

    dgmcomm=incomm

    jjxx = zzz * indimension

    rgwkdim= (42 + jjxx*(inkrydim*zzz + 6) + inkrydim*zzz*(inkrydim*zzz + 3)) / 5 * 6

    if (rgwkdim.lt.0) then
       OFLWR "OOPS, integer overflow, decrease kyrlov dimension in dgsolve";
       WRFL rgwkdim; CFLST
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
       call dgmres(jjxx, rhs, solution, inmult, &
            inprecon,itol , tol, numiter, errest, ierr, iunit, rgwk, &
            rgwkdim, igwk, ligwk)
    else
       call dgmrespar(jjxx, rhs, solution, inmult, &
            inprecon,itol , tol, numiter, errest, ierr, iunit, rgwk, &
            rgwkdim, igwk, ligwk)
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

  end subroutine dgsolve00

  subroutine dgsolve0(rhs, solution, numiter, inmult, preconflag, inprecon, &
       intolerance, indimension, inkrydim,parflag,ierr)
#ifdef MPIFLAG
    use mpimod     !! MPI_COMM_WORLD
#endif
    implicit none
    integer,intent(in) :: indimension,inkrydim,parflag,preconflag
    integer,intent(out) :: ierr
    real*8,intent(in) :: intolerance
    DATATYPE,intent(in) :: rhs(indimension)
    integer,intent(out) :: numiter
    DATATYPE,intent(out) :: solution(indimension)
    external inmult, inprecon
#ifndef MPIFLAG
  integer,parameter :: MPI_COMM_WORLD=(-798)
#endif

    call dgsolve00(rhs, solution, numiter, inmult, preconflag, inprecon, &
       intolerance, indimension, inkrydim,parflag,ierr,MPI_COMM_WORLD)
  end subroutine dgsolve0
    
end module dgsolvemod


!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NEWTON SOLVE (IMPROVEDQUADFLAG=1,2) FOR ORBITALS
!!!!!!!!!!!!!!!!!!!!!!!!!!

module quadspfsmod
contains

subroutine quadspfs(inspfs,jjcalls)
  use parameters
  use mpimod
  use linearmod
  use dgsolvemod
  use orbdermod
  use mpi_orbsetmod
  use orbgathersubmod
  use mpisubmod
  use jacinitmod
  use spfsubmod
  implicit none
  DATATYPE,intent(inout) :: inspfs(spfsize,nspf)
  integer,intent(out) :: jjcalls
  integer :: ierr,minerr,maxnorbs,lastmpiorb
  real*8 :: orthogerror,indev,dev,mynorm,tempquadtol
  DATATYPE, allocatable ::  invector(:,:), errvec(:,:), solution(:,:), &
       com_errvec(:,:), com_solution(:,:), workvec(:,:)

  lastmpiorb=firstmpiorb+orbsperproc-1

  if (parorbsplit.eq.1) then
!! TEMP? was     maxnorbs=nzprocsperset*orbsperproc
     maxnorbs=maxprocsperset*orbsperproc
  else
     maxnorbs=nspf
  endif

  allocate( invector(spfsize,maxnorbs), errvec(spfsize,maxnorbs), &
       solution(spfsize,maxnorbs), workvec(spfsize,maxnorbs) )
  invector=0; errvec=0; solution=0; workvec=0


!!$ 06-16 
!!$   Don't understand jacsymflag!  See    workvec=invector+solution    below
!!$
  if (jacsymflag.eq.0.and.jacsymquad.ne.0) then
     OFLWR "setting jacsymflag=1 for orbital quad"; CFL
     jacsymflag=1
  endif

  effective_cmf_linearflag=0

  invector(:,1:nspf)=inspfs(:,:)

  call spf_orthogit(invector,orthogerror)

  call jacinit0(0,invector,0d0)

  call spf_linear_derivs0(0,0,0d0,invector,errvec,1,1)

  call apply_spf_constraints(errvec)

  indev=abs(hermdot(errvec,errvec,totspfdim))
  if (parorbsplit.eq.3) then
     call mympirealreduceone(indev)
  endif
  indev=sqrt(indev)
  
  OFLWR "   Quad orbitals: Deviation is     ", indev; CFL

  ierr=1; minerr=1
  tempquadtol=quadtol

  do while (ierr.ne.0.or.minerr.ne.0)

     ierr=0; minerr=0

     solution=invector   !! guess

     if (orbcompact.ne.0) then
        allocate( com_errvec(spfsmallsize,maxnorbs), com_solution(spfsmallsize,maxnorbs) )
        com_errvec=0; com_solution=0
        call spfs_compact(errvec,com_errvec)
        call spfs_compact(solution,com_solution)
        if (parorbsplit.eq.3) then
           call dgsolve0( com_errvec, com_solution, jjcalls, quadopcompact,0,dummysub, &
                tempquadtol,spfsmallsize*nspf,min(maxdgdim,spfsmallsize*nspf*nprocs),1,ierr)
        elseif (parorbsplit.eq.1) then
           if (firstmpiorb.le.nspf) then
              call dgsolve00( com_errvec(:,firstmpiorb:lastmpiorb), &
                   com_solution(:,firstmpiorb:lastmpiorb), jjcalls, parquadopcompact,&
                   0,dummysub,tempquadtol,spfsmallsize*orbsperproc,&
                   min(maxdgdim,spfsmallsize*nspf),1,ierr,NZ_COMM_ORB(myorbset))
           endif
           call mpiorbgather(com_solution,spfsmallsize)
        else
           call dgsolve0( com_errvec, com_solution, jjcalls, quadopcompact,0,dummysub, &
                tempquadtol,spfsmallsize*nspf,min(maxdgdim,spfsmallsize*nspf),0,ierr)
        endif
        call spfs_expand(com_solution,solution)
        deallocate(com_errvec,com_solution)
     else
        if (parorbsplit.eq.3) then
           call dgsolve0( errvec, solution, jjcalls, quadoperate,0,dummysub, &
                tempquadtol,spfsize*nspf,min(maxdgdim,spfsize*nspf*nprocs),1,ierr)
        elseif (parorbsplit.eq.1) then
           if (firstmpiorb.le.nspf) then
              call dgsolve00( errvec(:,firstmpiorb:lastmpiorb), &
                   solution(:,firstmpiorb:lastmpiorb), jjcalls, parquadoperate,0,dummysub, &
                   tempquadtol,spfsize*orbsperproc,min(maxdgdim,spfsize*nspf),1,ierr,&
                   NZ_COMM_ORB(myorbset))
           endif
           call mpiorbgather(solution,spfsize)
        else
           call dgsolve0( errvec, solution, jjcalls, quadoperate,0,dummysub, &
                tempquadtol,spfsize*nspf,min(maxdgdim,spfsize*nspf),0,ierr)
        endif
     endif
     minerr=ierr
     call mympiimax(ierr)
     call mympiimin(minerr)
     if (ierr.ne.0.or.minerr.ne.0) then
        tempquadtol=sqrt(tempquadtol)
        OFLWR "  *** Error in dgsolve, trying with lower tolerance",ierr,minerr
        WRFL  "  *** tolerance now", tempquadtol; CFL
     else

        mynorm=abs(hermdot(solution,solution,totspfdim))
        if (parorbsplit.eq.3) then
           call mympirealreduceone(mynorm)
        endif
        mynorm=sqrt(mynorm)
  
        if (mynorm.gt.maxquadnorm*nspf) then
           solution=solution*maxquadnorm*nspf/mynorm
        endif
!! HOW IT SHOULD BE
!! logic now in getparams        if (jacsymflag.eq.0.or.(numfrozen.gt.0.and.exact_exchange.eq.0)) then

        select case(jacquaddir)
        case(1)

           workvec=invector-solution * 0.5d0   !! TEMP 0.5 ?  Works better

        case(-1)
!!
!! WHAT?  WHY?  This is what I've been doing, with jacsymflag.  It works,
!!              with plus not minus.  Where is the sign error?  Or something?
!!              What is going on with jacsymflag.ne.0 and numfrozen.eq.0 ?????
!!
           workvec=invector+solution * 0.5d0   !! TEMP 0.5 ?

        case default
           OFLWR "programmer fail jacquaddir=",jacquaddir; CFLST
        end select

        call spf_orthogit(workvec,orthogerror)

        call spf_linear_derivs0(0,0,0d0,workvec,errvec,1,1)

        call apply_spf_constraints(errvec)

        dev=abs(hermdot(errvec,errvec,totspfdim))
        if (parorbsplit.eq.3) then
           call mympirealreduceone(dev)
        endif
        dev=sqrt(dev)
  
        OFLWR "   Quad orbitals: Deviation is now ", dev; CFL
        if (dev.le.indev) then
           inspfs(:,:) = workvec(:,1:nspf)
        else
           tempquadtol=sqrt(tempquadtol)
           OFLWR "        ... rejecting step, trying with lower tolerance"
           WRFL  "  *** tolerance now", tempquadtol; CFL
           ierr=1
        endif
     endif
  enddo  !! do while ierr.ne.0.or.minerr.ne.0

  deallocate( invector,errvec,solution,workvec)

contains
  subroutine dummysub()
  end subroutine dummysub

!! SUBROUTINES PASSED TO DGSOLVE FOR ORBITAL NEWTON SOLVE IMPROVEDQUADFLAG 2 and 3

  subroutine quadoperate(notusedint,inspfs,outspfs)
    use parameters
    use jacopmod
    implicit none
    integer :: notusedint
    DATATYPE,intent(in) ::  inspfs(totspfdim)
    DATATYPE,intent(out) :: outspfs(totspfdim)
    DATATYPE :: workspfs(totspfdim)              !! AUTOMATIC
    workspfs=0
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
    use jacopmod
    implicit none
    integer :: notusedint
    DATATYPE,intent(in) ::  com_inspfs(spfsmallsize*nspf)
    DATATYPE,intent(out) :: com_outspfs(spfsmallsize*nspf)
    DATATYPE :: inspfs(spfsize*nspf),outspfs(spfsize*nspf)  !! AUTOMATIC
    inspfs=0; outspfs=0
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

  subroutine parquadopcompact(notusedint,com_inspfs,com_outspfs)
    use parameters
    use jacopmod
    use mpi_orbsetmod
    implicit none
    integer :: notusedint
    DATATYPE,intent(in) :: com_inspfs(spfsmallsize,firstmpiorb:firstmpiorb+orbsperproc-1)
    DATATYPE,intent(out) :: com_outspfs(spfsmallsize,firstmpiorb:firstmpiorb+orbsperproc-1)
    DATATYPE ::  inspfs(spfsize,firstmpiorb:firstmpiorb+orbsperproc-1),  &
         outspfs(spfsize,firstmpiorb:firstmpiorb+orbsperproc-1)  !! AUTOMATIC
    inspfs=0; outspfs=0
    call spfs_expand_local(com_inspfs,inspfs)
    call parjacoperate0(0,0,inspfs,outspfs)
    call spfs_compact_local(outspfs,com_outspfs)
  end subroutine parquadopcompact

  subroutine parquadoperate(notusedint,inspfs,outspfs)
    use parameters
    use jacopmod
    use mpi_orbsetmod
    implicit none
    integer :: notusedint
    DATATYPE,intent(in) :: inspfs(spfsize,firstmpiorb:firstmpiorb+orbsperproc-1)
    DATATYPE,intent(out) :: outspfs(spfsize,firstmpiorb:firstmpiorb+orbsperproc-1)
    call parjacoperate0(0,0,inspfs,outspfs)
  end subroutine parquadoperate

end subroutine quadspfs

end module quadspfsmod


!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NEWTON SOLVE (IMPROVEDQUADFLAG=1,3) FOR AVECTOR
!!!!!!!!!!!!!!!!!!!!!!!!!!


module quadavecmod
  implicit none
  DATATYPE :: quadexpect=0d0
  contains

subroutine aaonedinit(www,inavector) 
  use xxxmod
  use r_parameters
  use walkmod
  use dotmod
  use sparsemultmod
  use basissubmod
  use mpisubmod
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
  call sparseconfigmult(www,jacaa, jacaamult, yyy%cptr(0), yyysptr(0),1,1,1,0,0d0,-1)

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
  use utilmod
  implicit none
  DATATYPE,intent(inout) :: inavector(www%totadim,mcscfnum)
  DATATYPE :: nullvector(numr)
  integer,intent(out) :: jjcalls
  integer :: imc

  nullvector(:)=0

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
           call sparsequadavector(inavector(:,imc),jjcalls)
        else
           call sparsequadavector(nullvector(:),jjcalls)
        endif
     endif
  enddo

end subroutine quadavector


subroutine sparsequadavector(inavector,jjcalls0)
  use parameters
  use mpimod
  use xxxmod
  use configmod
  use dgsolvemod
  use sparsemultmod
  use basissubmod
  use mpisubmod
  implicit none
  DATATYPE, intent(inout) ::  inavector(numr,www%firstconfig:www%lastconfig)
  integer,intent(out) :: jjcalls0
  integer :: jjcalls, ss, maxdim, mysize,flag,ierr,minerr
  real*8 ::  dev,  thisaerror, nowaerror
  DATATYPE :: csum 
  DATATYPE, allocatable ::  vector(:,:), vector2(:,:), vector3(:,:), &
       smallvectorspin(:,:),smallvectorspin2(:,:),shufflevector(:,:),&
       shufflevector2(:,:)

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
     call sparseconfigmult(www,vector, vector2, yyy%cptr(0), yyysptr(0), &
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
!!     thisaerror=min(0.1d0,sqrt(aerror/min(dev,1d0)))
!! v1.28
!!     thisaerror=min(0.316d0,sqrt(aerror/min(dev,1d0)))

     thisaerror=min(0.316d0,sqrt(sqrt(aerror/min(dev,1d0))))

     OFL;write(mpifileptr,'(A20,E12.5,A6,2E12.5,A7,100F14.7)') &
          "   SPARSEQUAD: DEV", dev, " TOL ",aerror,thisaerror,"ENERGY",quadexpect; CFL

     if (abs(dev).lt.aerror.or.ss.ge.10) then
        if (ss.ge.10) then
           OFLWR "10 iterations, quad aborting"; CFL
        endif
        if (www%lastconfig.ge.www%firstconfig) then
           inavector(:,:) = vector(:,www%firstconfig:www%lastconfig)
        endif
        deallocate(smallvectorspin,smallvectorspin2)
        deallocate(vector, vector2, vector3)
        return
     endif

     smallvectorspin(:,:)=0d0
     if (www%topconfig.ge.www%botconfig) then
        call basis_transformto_local(www,numr,vector(:,www%botconfig:www%topconfig),smallvectorspin)
     endif
     smallvectorspin2(:,:)=smallvectorspin(:,:)    !! guess

!!$ v1.23     maxdim=min(maxdgdim,numr*www%numdfbasis)

!! back to maxaorder here so that can do different orbitals and avec improvedquadflag.eq.3
     maxdim=min(maxaorder,numr*www%numdfbasis)

     mysize=numr*(dwwptr%topdfbasis-dwwptr%botdfbasis+1)

     flag=0; jjcalls=0

     if (nzflag.eq.0) then
        if (mysize.eq.0) then
           print *, "ACK, CAN'T DO A-VECTOR QUAD WITH ZERO CONFIGS PER PROCESSOR RANK",myrank,mysize
           flag=99
        endif
        call mympiireduceone(flag)
        if (flag.ne.0) then
           call mpistop()
        endif
     endif

     call mpibarrier()
     if (use_dfwalktype.and.shuffle_dfwalktype) then
        allocate(shufflevector(numr,dwwptr%botdfbasis:dwwptr%topdfbasis),&
             shufflevector2(numr,dwwptr%botdfbasis:dwwptr%topdfbasis))
        if (dwwptr%topdfbasis.ge.dwwptr%botdfbasis) then
           shufflevector=0; shufflevector2=0
        endif
        call basis_shuffle(numr,www,smallvectorspin,dwwptr,shufflevector)
        call basis_shuffle(numr,www,smallvectorspin2,dwwptr,shufflevector2)
     else
        allocate(shufflevector(numr,1),shufflevector2(numr,1))
     endif

     nowaerror=thisaerror

     ierr=1; minerr=1

     do while (ierr.ne.0.or.minerr.ne.0)
        ierr=0; minerr=0

        if (nzflag.eq.0.or.dwwptr%nzrank.gt.0) then
           if (use_dfwalktype.and.shuffle_dfwalktype) then
              call dgsolve00( shufflevector, shufflevector2, jjcalls, paraamult,&
                   quadprecon,parquadpreconsub, nowaerror,mysize,maxdim,1,ierr,dwwptr%NZ_COMM)
           else
              call dgsolve00( smallvectorspin, smallvectorspin2, jjcalls, paraamult,&
                   quadprecon,parquadpreconsub, nowaerror,mysize,maxdim,1,ierr,dwwptr%NZ_COMM)
           endif
        else
           smallvectorspin2=0
           if (use_dfwalktype.and.shuffle_dfwalktype) then
              shufflevector2=0
           endif
        endif
        call mpibarrier()

        minerr=ierr
        call mympiimax(ierr)
        call mympiimin(minerr)
        if (ierr.ne.0.or.minerr.ne.0) then
           nowaerror=sqrt(nowaerror)
           OFLWR " ->Error in dgsolve, lowering tolerance",nowaerror; CFL
        endif

!           OFLWR "          Error in dgsolve, not changing a-vector",ierr,minerr; CFL
!           if (www%lastconfig.ge.www%firstconfig) then
!              inavector(:,:) = vector(:,www%firstconfig:www%lastconfig)
!           endif
!           deallocate(smallvectorspin,smallvectorspin2)
!           deallocate(vector, vector2, vector3)
!           deallocate(shufflevector,shufflevector2)
!           return
!        endif

     enddo  !! do while ierr 

     if (use_dfwalktype.and.shuffle_dfwalktype) then
        call basis_shuffle(numr,dwwptr,shufflevector2,www,smallvectorspin2)
     endif
     deallocate(shufflevector,shufflevector2)

!!$  call basicblocklansolve(1,mysize,maxdim,maxdim,smallvectorspin,smallvectorspin2,1,&
!!$     parhrmult,parhrdotsub,quadexpect)

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

contains

  subroutine paraamult(notusedint,inavectorspin,outavectorspin)
    use r_parameters
    use configmod
    use parblocklanmod
    implicit none
    integer,intent(in) :: notusedint
    DATATYPE,intent(in) :: inavectorspin(numr,dwwptr%botdfbasis:dwwptr%topdfbasis)
    DATATYPE,intent(out) :: outavectorspin(numr,dwwptr%botdfbasis:dwwptr%topdfbasis)
    call parblockconfigmult(inavectorspin,outavectorspin)
    if (dwwptr%topdfbasis.ge.dwwptr%botdfbasis) then
       outavectorspin(:,:)= outavectorspin(:,:) - quadexpect*inavectorspin(:,:)
    endif
  end subroutine paraamult

  subroutine parquadpreconsub(notusedint, inavectorspin,outavectorspin)
    use fileptrmod
    use r_parameters
    use sparse_parameters
    use xxxmod
    use configmod
    implicit none
    integer,intent(in) :: notusedint
    DATATYPE,intent(in) :: inavectorspin(numr,dwwptr%botdfbasis:dwwptr%topdfbasis)
    DATATYPE,intent(out) :: outavectorspin(numr,dwwptr%botdfbasis:dwwptr%topdfbasis)

    call parquadpreconsub0(dwwptr,yyy%cptr(0),yyy%sptrptr(0),&
         notusedint,inavectorspin,outavectorspin)

  end subroutine parquadpreconsub

  subroutine parquadpreconsub0(wwin,cptr,sptr,notusedint, inavectorspin,outavectorspin)
    use fileptrmod
    use r_parameters
    use sparse_parameters
    use walkmod
    use configptrmod
    use sparseptrmod
    use sparsemultmod
    use basissubmod
    implicit none
    type(walktype),intent(in) :: wwin
    type(CONFIGPTR),intent(in) :: cptr
    type(SPARSEPTR),intent(in) :: sptr
    integer,intent(in) :: notusedint
    DATATYPE,intent(in) :: inavectorspin(numr,wwin%botdfbasis:wwin%topdfbasis)
    DATATYPE,intent(out) :: outavectorspin(numr,wwin%botdfbasis:wwin%topdfbasis)
    DATATYPE :: inavector(numr,wwin%botconfig:wwin%topconfig+1),&
         outavector(numr,wwin%botconfig:wwin%topconfig+1)  !!AUTOMATIC
    
    if (sparseconfigflag.eq.0) then
       OFLWR "error, no parquadpreconsub if sparseconfigflag=0"; CFLST
    endif
    
    call basis_transformfrom_local(wwin,numr,inavectorspin,inavector)

    call parsparseconfigpreconmult(wwin,inavector, outavector, cptr, sptr,&
         1,1,1,1, quadexpect,0d0)

    call basis_transformto_local(wwin,numr,outavector,outavectorspin)

  end subroutine parquadpreconsub0

!!$
!!$subroutine parhrmult(lanblocknum,inavectorspin,outavectorspin)
!!$  use parameters
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

end subroutine sparsequadavector


subroutine nonsparsequadavector(www,avectorout)
  use mpimod
  use parameters
  use xxxmod
  use walkmod
  use sparsemultmod
  use basissubmod
  use mpisubmod
  use asssubmod
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
     call sparseconfigmult(www,avectorout(:),err(:),yyy%cptr(0),yyysptr(0),1,1,1,0,0d0,-1)

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

end module quadavecmod

