

!! SUBROUTINES FOR IMPROVEDQUADFLAG, AVECTOR AND ORBITAL.  CALLS DGMRES ROUTINE.

#include "Definitions.INC"

!      SUBROUTINE DGMRES (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,
!     +   ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX, RGWK, LRGW,
!     +   IGWK, LIGW, RWORK, IWORK)

!! ON INPUT SOLUTION IS GUESS.

subroutine dgsolve0(rhs, solution, numiter, inmult, preconflag, inprecon, intolerance, indimension, inkrydim,parflag)
  use parameters
  implicit none
  integer, parameter :: jacwrite=0
  integer :: indimension, inkrydim,parflag,itol, numiter, ierr, iunit, jjxx, preconflag,rgwkdim
  DATATYPE :: rhs(indimension), solution(indimension)
  external :: inmult, dummysub, inprecon
  integer :: igwk(20), ligwk=20, nullint1,nullint2,nullint3,nullint4, nullint5
  real*8 :: nulldouble1, errest, nulldouble2, nulldouble3, nulldouble4, tol, intolerance
  real*8, allocatable :: rgwk(:)
#ifdef REALGO
  integer, parameter :: zzz=1
#else
  integer, parameter :: zzz=2
#endif

!!$  looks ok now dgmres bugfix 10-2015
!!$  mindim=indimension; maxdim=indimension
!!$  call mympiimax(maxdim); call mympiimin(mindim)
!!$  if (mindim.ne.maxdim) then
!!$     OFLWR "Error, for dgsolve now must use same dimension all processors", mindim,maxdim; CFLST
!!$  endif

  jjxx = zzz * indimension

  rgwkdim= (42 + jjxx*(inkrydim*zzz + 6) + inkrydim*zzz*(inkrydim*zzz + 3)) / 5 * 6

  if (rgwkdim.lt.0) then
     OFLWR "OOPS, integer overflow, decrease maxexpodim or reduce size of orbitals (per processor)",rgwkdim; CFLST
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
     call dgmres(jjxx, rhs, solution, nullint1, nullint2, nullint3, nulldouble1, nullint3, inmult, inprecon,itol , tol, nullint4, numiter, errest, ierr, iunit, nulldouble2, nulldouble3, rgwk, rgwkdim, igwk, ligwk, nulldouble4, nullint5)
  else
     call dgmrespar(jjxx, rhs, solution, nullint1, nullint2, nullint3, nulldouble1, nullint3, inmult, inprecon,itol , tol, nullint4, numiter, errest, ierr, iunit, nulldouble2, nulldouble3, rgwk, rgwkdim, igwk, ligwk, nulldouble4, nullint5)
  endif

!!$  if (debugflag.ne.0) then
!!$     OFLWR "     DGSOLVE - restarts,iters,errest ",igwk(5), numiter,errest; CFL
!!$  endif

  if (jacwrite==1) then
     OFLWR "   NEWTON: numiter ", numiter, " error ", errest;     call closefile()
  endif
  if (ierr.eq.1) then     !! AUTHOR COMMENTS ARE WRONG IN DGMRES.F ierr=1 denotes max restarts (zero) reached !! 
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
!!! NEWTON SOLVE (IMPROVEDQUADFLAG) FOR ORBITALS
!!!!!!!!!!!!!!!!!!!!!!!!!!


recursive subroutine quadoperate(notusedint,inspfs,outspfs)
   use parameters
   implicit none
   integer :: notusedint
   DATATYPE ::  inspfs(totspfdim), outspfs(totspfdim)

   call jacoperate(inspfs,outspfs)

end subroutine quadoperate


subroutine quadspfs(inspfs,jjcalls)
  use parameters
  use mpimod
  use linearmod
  implicit none
  EXTERNAL :: quadoperate,dummysub
  integer :: jjcalls,icount,maxdim
  real*8 :: orthogerror,dev,mynorm
  DATATYPE ::  hermdot,inspfs(totspfdim)  
  DATATYPE ::       vector(totspfdim), vector2(totspfdim), vector3(totspfdim)

  if (jacsymflag.ne.1) then
     OFLWR "setting jacsymflag=1 for orbital quad"; CFL
     jacsymflag=1
  endif

  effective_cmf_linearflag=0
  effective_cmf_spfflag=1

  vector=inspfs

  icount=0

333 continue

  icount=icount+1

  call spf_orthogit(vector,orthogerror)

  call jacinit(vector,0d0)

  call spf_linear_derivs(0d0,vector,vector2)

  dev=abs(hermdot(vector2,vector2,totspfdim))
  if (parorbsplit.eq.3) then
     call mympirealreduceone(dev)
  endif
  dev=sqrt(dev)
  


  OFLWR "   Quad orbitals: Deviation is     ", dev
  WRFL  "                  Orthog error is  ", orthogerror; CFL

  if (dev.lt.stopthresh.or.icount.gt.1) then
     inspfs = vector
     OFLWR "    --> Converged newton"; CFL
     return
  else
     vector3=vector   !! guess
     if (parorbsplit.eq.3) then
        maxdim=min(totspfdim*nprocs,maxexpodim)
        call dgsolve0( vector2, vector3, jjcalls, quadoperate,0,dummysub, quadtol,totspfdim,maxdim,1)
     else
        maxdim=min(totspfdim,maxexpodim)
        call dgsolve0( vector2, vector3, jjcalls, quadoperate,0,dummysub, quadtol,totspfdim,maxdim,0)
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
  endif

  go to 333

end subroutine quadspfs


!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NEWTON SOLVE (IMPROVEDQUADFLAG) FOR AVECTOR
!!!!!!!!!!!!!!!!!!!!!!!!!!


module aaonedmod
  implicit none
  DATATYPE :: quadexpect=0d0
end module



subroutine aaonedinit(inavector) 
  use aaonedmod
  use xxxmod
  use parameters
  implicit none
  DATATYPE,intent(in) :: inavector(totadim)
  DATATYPE :: jacaa(totadim), jacaamult(totadim) !! AUTOMATIC
  DATATYPE :: dot, csum2

  jacaa(:)=inavector(:)

  if (allspinproject==1) then
     call configspin_project(jacaa,1)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(jacaa,numr)
  endif

  csum2=dot(jacaa,jacaa,totadim)
  if (parconsplit.ne.0) then
     call mympireduceone(csum2)
  endif
  jacaa(:)=jacaa(:)/sqrt(csum2)

  call sparseconfigmult(jacaa, jacaamult, yyy%cptr(0), yyy%sptr(0),1,1,0,1,0d0)
  if (allspinproject==1) then
     call configspin_project(jacaamult, 1)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(jacaamult, numr)
  endif

  quadexpect=dot(jacaa,jacaamult,totadim)
  if (parconsplit.ne.0) then
     call mympireduceone(quadexpect)
  endif

end subroutine aaonedinit



subroutine quadavector(inavector,jjcalls)
  use parameters
  implicit none
  DATATYPE :: inavector(totadim,mcscfnum)
  integer :: jjcalls,imc
  do imc=1,mcscfnum
     if (sparseconfigflag==0) then
        call nonsparsequadavector(inavector(:,imc))
        jjcalls=0
     else
        call sparsequadavector(inavector(:,imc),jjcalls)
     endif
  enddo
end subroutine quadavector


subroutine sparsequadavector(inavector,jjcalls0)
  use parameters
  use mpimod
  use aaonedmod
  use xxxmod
  implicit none
  EXTERNAL :: paraamult, parquadpreconsub
  DATATYPE, intent(inout) ::  inavector(numr,firstconfig:lastconfig)
  integer :: jjcalls, ss,jjcalls0,maxdim,mysize
  real*8 ::  dev,  thisaerror
  DATATYPE :: dot, hermdot,csum 
  DATATYPE, allocatable ::  vector(:,:), vector2(:,:), vector3(:,:), &
       smallvectorspin(:,:),smallvectorspin2(:,:)

  if (sparseconfigflag.eq.0) then
     OFLWR "Error, must use sparseconfigflag.ne.0 for sparsequadavector"; CFLST
  endif

  if (aerror.lt.1d-7) then
     OFLWR "Error, set tolerance parameter arror at least 1d-7 for reliable performance sparse quad a-vector."; CFLST
  endif

  allocate(smallvectorspin(numr,botdfbasis:topdfbasis), smallvectorspin2(numr,botdfbasis:topdfbasis))

  allocate( vector(numr,firstconfig:lastconfig), vector2(numr,firstconfig:lastconfig), vector3(numr,firstconfig:lastconfig))

  vector(:,:)=inavector(:,:)

  jjcalls0=0
  ss=0

333 continue

  if (allspinproject==1) then
     call configspin_project(vector, 1)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(vector, numr)
  endif

  call aaonedinit(vector)

  call sparseconfigmult(vector, vector2, yyy%cptr(0), yyy%sptr(0), 1,1,0,1,0d0)

  if (allspinproject.ne.0) then
     call configspin_project(vector2,0)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(vector2,numr)
  endif

  vector3=vector2-quadexpect*vector              !! error term.

  dev=abs(hermdot(vector3,vector3,totadim))
  if (parconsplit.ne.0) then
     call mympirealreduceone(dev)
  endif
  dev=sqrt(dev)

!! 10-2015
  thisaerror=min(0.1d0,aerror/min(dev,1d0))

  OFL;write(mpifileptr,'(A20,E12.5,A6,2E12.5,A7,100F14.7)') "   SPARSEQUAD: DEV", dev, " TOL ",aerror,thisaerror,"ENERGY",quadexpect; CFL

  if (abs(dev).lt.aerror.or.ss.ge.10) then
     inavector = vector
     if (ss.ge.10) then
        OFLWR "10 iterations, quad aborting"; CFL
     endif

     !OFLWR "     Converged sparse quad: dev ", dev, "  Expect ", quadexpect, " tol is ", thisaerror; CFL

     deallocate(smallvectorspin,smallvectorspin2)
     deallocate(vector, vector2, vector3)
     return
  endif

!! zero dimension will fail in dgsolve, but putting this logic in anyway
  if (topdfbasis-botdfbasis+1.ne.0) then
     smallvectorspin(:,:)=0d0
  endif
  call basis_transformto_local(numr,vector(:,botconfig),smallvectorspin(:,:))
  if (topdfbasis-botdfbasis+1.ne.0) then
     smallvectorspin2(:,:)=smallvectorspin(:,:)    !! guess
  endif

  maxdim=min(maxaorder,numr*numdfbasis)
  mysize=numr*(topdfbasis-botdfbasis+1)

  call dgsolve0( smallvectorspin(:,:), smallvectorspin2(:,:), jjcalls, paraamult,quadprecon,parquadpreconsub, thisaerror,mysize,maxdim,1)

!!$  call basicblocklansolve(1,mysize,maxdim,maxdim,smallvectorspin,smallvectorspin2,1,parhrmult,parhrdotsub,quadexpect)
!!$  jjcalls=0

  vector3(:,:)=0d0; 

  call basis_transformfrom_local(numr,smallvectorspin2(:,:),vector3(:,botconfig))

  if (parconsplit.eq.0) then
     call mpiallgather(vector3(:,:),numconfig*numr,configsperproc*numr,maxconfigsperproc*numr)
  endif

  csum=dot(vector3,vector3,totadim)
  if (parconsplit.ne.0) then
     call mympireduceone(csum)
  endif
  vector=vector3/sqrt(csum)
     
!!  OFLWR "    ITERATIONS ***  ",jjcalls

  jjcalls0=jjcalls0+jjcalls

  ss=ss+1
  go to 333

end subroutine sparsequadavector


subroutine nonsparsequadavector(avectorout)
  use mpimod
  use aaonedmod
  use parameters
  use xxxmod
  implicit none
  DATATYPE,intent(inout) :: avectorout(totadim)
  DATATYPE :: dot, hermdot
  CNORMTYPE :: norm
  real*8 :: dev
  integer :: k, info,ss
  DATATYPE, allocatable :: halfmatel(:,:),halfmatel2(:,:),&
       spinmatel(:,:),spinerr(:),spinavectorout(:)
  DATATYPE :: bigconfigmatel(totadim,totadim),err(totadim)
  integer :: ipiv(totadim)

  if (sparseconfigflag/=0) then
     OFLWR "Error, nonsparsequadavector called but sparseconfigflag= ", sparseconfigflag; CFLST
  endif

  ss=0
333 continue

  if (allspinproject==1) then
     call configspin_project(avectorout, 1)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(avectorout, numr)
  endif

  call aaonedinit(avectorout)
  call sparseconfigmult(avectorout(:),err(:),yyy%cptr(0),yyy%sptr(0),1,1,0,1,0d0)

  if (allspinproject==1) then
     call configspin_project(err(:), 0)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(err(:), numr)
  endif

  err(:)=err(:)-quadexpect*avectorout(:)              !! error term.

  dev=sqrt(abs(hermdot(err(:),err(:),totadim)))

  OFL;write(mpifileptr,'(A22,E12.5,A6,E12.5,A7,100F14.7)') "   NONSPARSEQUAD: DEV", dev, " TOL ",aerror,"ENERGY",quadexpect; CFL

  if (abs(dev).lt.aerror.or.(ss.gt.10)) then

!#ifdef REALGO
!  OFL; write(mpifileptr,'(A40,F15.10, A5,1F18.10,A5,1E6.1)') "      Converged nonsparse quad: dev ", dev, "  E= ", quadexpect, " tol ", aerror; CFL
!#else
!  OFL; write(mpifileptr,'(A40,F15.10, A5,2F18.10,A5,1E6.1)') "      Converged nonsparse quad: dev ", dev, "  E= ", quadexpect, " tol ", aerror; CFL
!#endif

     return
  endif

!! no spin project; doing transform below
  call assemble_configmat(bigconfigmatel, yyy%cptr(0),1,1,1,1, 0d0)

  do k=1,totadim
     bigconfigmatel(k,k)=bigconfigmatel(k,k)-quadexpect
  enddo

  if (allspinproject.ne.0) then

!! can make this complex or whatever; multiplying by real matrix
!!   but doing it carefully, only condensing leading index as complex

     allocate(halfmatel(numspinconfig*numr,numconfig*numr), &
          halfmatel2(numconfig*numr,numspinconfig*numr), &
          spinmatel(numspinconfig*numr,numspinconfig*numr), &
          spinerr(numspinconfig*numr),spinavectorout(numspinconfig*numr))

     call configspin_transformto(numr*numconfig*numr,bigconfigmatel,halfmatel2)
     
     halfmatel(:,:)=TRANSPOSE(halfmatel2(:,:))

     call configspin_transformto(numr*numspinconfig*numr,halfmatel,spinmatel)

     spinmatel(:,:)=TRANSPOSE(spinmatel(:,:))

     call configspin_transformto(numr,avectorout(:),spinavectorout(:))
     call configspin_transformto(numr,err(:),spinerr(:))

     if (myrank.eq.1) then
        call MYGESV(numspinconfig*numr,1,spinmatel,numspinconfig*numr,ipiv,spinavectorout,numspinconfig*numr,info)
        spinavectorout=spinavectorout/sqrt(dot(spinavectorout,spinavectorout,numspinconfig*numr))
        if (info/=0) then
           OFLWR "Info in dgesv nonsparsequadavector= ", info; CFLST
        endif
     endif

     call mympibcast(spinavectorout,1,numspinconfig*numr)
     call configspin_transformfrom(numr,spinavectorout(:),avectorout(:))

     deallocate(halfmatel, halfmatel2, spinmatel,     spinerr, spinavectorout)

  else

     if (myrank.eq.1) then

        call MYGESV(totadim,1,bigconfigmatel,totadim,ipiv,avectorout,totadim,info)
        avectorout=avectorout/sqrt(dot(avectorout,avectorout,totadim))

        if (info/=0) then
           OFLWR "Info in dgesv nonsparsequadavector= ", info; CFLST
        endif
     endif

     call mympibcast(avectorout,1,totadim)

  endif

  if (dfrestrictflag.ne.0) then
     call df_project(avectorout, numr)
  endif

  norm=sqrt(dot(avectorout,avectorout,totadim))  !ok conv
  avectorout=avectorout/norm

  ss=ss+1
  go to 333
  
end subroutine nonsparsequadavector



recursive subroutine paraamult(notusedint,inavectorspin,outavectorspin)
  use parameters
  use aaonedmod
  implicit none
  integer :: notusedint
  DATATYPE,intent(in) :: inavectorspin(numr,botdfbasis:topdfbasis)
  DATATYPE,intent(out) :: outavectorspin(numr,botdfbasis:topdfbasis)

  call parblockconfigmult(inavectorspin,outavectorspin)

  outavectorspin(:,:)= outavectorspin(:,:) &
       - quadexpect*inavectorspin(:,:)

end subroutine paraamult



recursive subroutine parquadpreconsub(notusedint, inavectorspin,outavectorspin)
  use parameters
  use aaonedmod
  use xxxmod
  implicit none
  integer :: notusedint
  DATATYPE,intent(in) :: inavectorspin(numr,botdfbasis:topdfbasis)
  DATATYPE,intent(out) :: outavectorspin(numr,botdfbasis:topdfbasis)
  DATATYPE :: inavector(numr,botconfig:topconfig),outavector(numr,botconfig:topconfig)

  if (sparseconfigflag.eq.0) then
     OFLWR "error, no parquadpreconsub if sparseconfigflag=0"; CFLST
  endif

  call basis_transformfrom_local(numr,inavectorspin,inavector)

  call parsparseconfigpreconmult(inavector, outavector, yyy%cptr(0), yyy%sptr(0),1,1,1,1, quadexpect,0d0)

  call basis_transformto_local(numr,outavector,outavectorspin)

end subroutine parquadpreconsub



!!$
!!$recursive subroutine parhrmult(lanblocknum,inavectorspin,outavectorspin)
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

