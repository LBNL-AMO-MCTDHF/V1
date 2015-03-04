

!! SUBROUTINES FOR IMPROVEDQUADFLAG, AVECTOR AND ORBITAL.  CALLS DGMRES ROUTINE.

#include "Definitions.INC"

!! THIS IS dgmres_driver from DEPRECATED, minimally changed, plus additional
!!  routines for 2nd order relax.

!      SUBROUTINE DGMRES (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,
!     +   ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX, RGWK, LRGW,
!     +   IGWK, LIGW, RWORK, IWORK)



!! For Marquardt-levenberg; using hessianoperate to solve AX=B where A is
!!   hessian.  So pass B= jacobian onto f(phi) and X will be guess
!!   at phi solution where f(phi) is zero.  (hessianoperate is inumult)


!! ON INPUT SOLUTION IS GUESS.


subroutine dgsolve0(rhs, solution, numiter, inmult, preconflag, inprecon, intolerance, indimension, inkrydim,parflag)
  use parameters
  implicit none
  integer, parameter :: jacwrite=0
  integer :: indimension, inkrydim,parflag,itol, numiter, ierr, iunit, jjxx, preconflag,rgwkdim
  DATATYPE :: rhs(indimension), solution(indimension)
  external :: inmult, dummysub, inprecon
  integer :: igwk(20), ligwk=20, nullint1,nullint2,nullint3,nullint4, nullint5 ,maxdim,mindim
  real*8 :: nulldouble1, errest, nulldouble2, nulldouble3, nulldouble4, tol, intolerance
  real*8, allocatable :: rgwk(:)

  mindim=indimension; maxdim=indimension
  call mympiimax(maxdim); call mympiimin(mindim)
  if (mindim.ne.maxdim) then
     OFLWR "Error, for dgsolve now must use same dimension all processors", mindim,maxdim; CFLST
  endif


#ifdef REALGO
  jjxx = indimension
#else
  jjxx = 2*indimension
#endif
  rgwkdim= 3*(1 + jjxx*(inkrydim+6) + inkrydim*(inkrydim+3)) * 6/5
  allocate(rgwk(rgwkdim))

  itol=0;  iunit=0

#ifdef REALGO
  jjxx = indimension
#else
  jjxx = 2*indimension
#endif
  igwk=0
  igwk(1) = inkrydim    !! max krylov dim
  igwk(2) = inkrydim    !! max orthog
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

  OFLWR "        DGSOLVE - restarts ",igwk(5), "iters", numiter; CFL

  if (jacwrite==1) then
     OFLWR "   NEWTON: numiter ", numiter, " error ", errest;     call closefile()
  endif
!  if (ierr.eq.1) then     !! AUTHOR COMMENTS ARE WRONG IN DGMRES.F ierr=1 denotes !! 
!     OFLWR "WARNING! dgmres did not converge. increase krylov order?"; CFL
!  else  
if (ierr/=0) then
     OFLWR "Error dgmres ", ierr," TEMP CONTINUE"; CFL  !!ST
  endif

  deallocate(rgwk)

end subroutine dgsolve0


!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! ONEDMOD SUBROUTINES REMOVED !! USE JACOPERATE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine quadoperate(nullint,inspfs,outspfs)
  use parameters
!!  use onedmod
  implicit none
  integer :: nullint
  DATATYPE ::  inspfs(totspfdim), outspfs(totspfdim)

!  quadcalled=quadcalled+1
!  if (mod(quadcalled,1000).eq.0) then
!     call openfile(); write(mpifileptr,*) "    Quad -iterations ", quadcalled; call closefile()
!  endif

  OFLWR "PROGRAM ME (CALL JACOPERATE)"; CFLST

!  if (whichquad.eq.1) then
!     outspfs=outspfs-jacexpect*inspfs
!  else
!     call realproject(outspfs,jacspftemp,jacspf)
!     outspfs=outspfs-jacspftemp
!     call realderproject(jacspfmult,jacspftemp,jacspf,inspfs)
!     outspfs=outspfs-jacspftemp
!
!  endif

end subroutine quadoperate



subroutine quadspfs(inspfs,jjcalls)
  use parameters
!  use onedmod
  implicit none
  EXTERNAL :: quadoperate,dummysub
  integer :: jjcalls, ss
  real*8 :: nulldouble
  DATATYPE :: dev, dot,inspfs(totspfdim)  
  DATATYPE, allocatable ::       vector(:), vector2(:), vector3(:)!!, work(:)

  allocate(       vector(totspfdim), vector2(totspfdim), vector3(totspfdim))

  vector=inspfs
  call spf_orthogit(vector,nulldouble)


  ss=0

333 continue

  call spf_orthogit(vector,nulldouble)

  OFLWR "PROGRAM ME (JAC INIT)"; CFLST

!  call onedinit(vector) 
!  call wmult(vector,vector2,ireduced)

!  if (whichquad.ne.1) then
!     call realproject(vector2,vector3,vector)
!     vector2=vector3-vector2                     !! error term.
!  else
!     vector2=(-1)*vector2                     
!  endif
!  dev=dot(vector2,vector2,totspfdim)
  dev=0

  call openfile(); write(mpifileptr,*)  "    NEWTON: Deviation is ", dev; call closefile()

  if (abs(dev).lt.stopthresh) then
     inspfs = vector
     OFLWR "    --> Converged newton, Deviation ", dev; CFL
     OFLWR "PROGME"; CFLST
!!!     call onedfinal()
     deallocate(vector, vector2, vector3);     return
  else
     vector3=inspfs
     call dgsolve0( vector2, vector3, jjcalls, quadoperate,0,dummysub, quadtol,totspfdim,maxexpodim,0)
     OFLWR "PROGME"; CFLST
!!     call onedfinal()
     if (whichquad.ne.1) then
        vector=vector+ vector3
     else
        vector=vector3
     endif
  endif


  ss=1
  go to 333

end subroutine quadspfs







module aaonedmod
  implicit none
  integer, parameter :: ireduced=1

  DATATYPE :: jacexpect=0.d0
  DATATYPE :: quadexpect=0d0 !! ONLY FOR REPORTING ENERGY, NOT USED IN EQNS    no used for rayleigh whichquad=1
  DATATYPE, allocatable :: jacaa(:,:), jacaamult(:,:)
  DATATYPE, allocatable :: bigconfigmatel(:,:), tempconfigmatel(:,:), jacaaproj(:,:), jacaaproj2(:,:), err(:)
  real*8, allocatable :: realconfigmatel(:,:,:,:)
  integer, allocatable :: ipiv(:)

#ifdef REALGO
  integer, parameter :: zzz=1
#else
  integer, parameter :: zzz=2
#endif

end module

subroutine quadpreconsq(nullint, inavector,outavector)
  use parameters
  use aaonedmod
  implicit none
  integer :: nullint
  DATATYPE :: inavector(totadim), outavector(totadim),jtempaa(totadim)
  call quadpreconsub(nullint, inavector,jtempaa);  jtempaa=conjg((0d0,0d0)+jtempaa)
  call quadpreconsub(nullint, jtempaa,outavector);  outavector=conjg((0d0,0d0)+outavector)
end subroutine quadpreconsq

subroutine quadpreconsub(nullint, inavector,outavector)
  use parameters
  use aaonedmod
  use xxxmod
  implicit none
  integer :: nullint
  DATATYPE :: inavector(totadim), outavector(totadim),tempvector(totadim)

  tempvector=inavector

  if (allspinproject==1) then
     call configspin_projectall(tempvector, 0)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict(tempvector, numr)
  endif

  
  if (whichquad.eq.0) then
     call sparseconfigdiagmult(tempvector, outavector, yyy%cptr(0),yyy%sptr(0), 1,1,1,1, jacexpect,0d0) !+quadpreconshift,0d0)
  else
     call sparseconfigdiagmult(tempvector, outavector, yyy%cptr(0),yyy%sptr(0), 1,1,1,1, quadexpect,0d0) !+quadpreconshift,0d0)
  endif
  if (allspinproject==1) then
     call configspin_projectall(outavector, 0)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict(outavector, numr)
  endif
end subroutine



subroutine aamultsq(nullint,inavector,outavector)
  use parameters
  use aaonedmod
  implicit none
  integer :: nullint
  DATATYPE :: inavector(totadim), outavector(totadim),jtempaa(totadim)
  call aamult(nullint,inavector,jtempaa)
  call aamultconjg(nullint,jtempaa,outavector)
end subroutine aamultsq

subroutine aamultconjg(nullint,inavector,outavector)
  use parameters
  use aaonedmod
  implicit none
  integer :: nullint
  DATATYPE :: inavector(totadim), outavector(totadim)
  inavector=conjg((0d0,0d0)+inavector)
  call aamult(nullint,inavector,outavector)
  inavector=conjg((0d0,0d0)+inavector)
  outavector=conjg((0d0,0d0)+outavector)  
end subroutine aamultconjg

subroutine aamult(nullint, inavector0,outavector)
  use parameters
  use xxxmod
  use aaonedmod
  implicit none
  integer :: nullint
  DATATYPE :: inavector(numconfig,numr), outavector(numconfig,numr), dot  ,jacaatemp(numconfig,numr), &
       inavector0(numconfig, numr)

  inavector(:,:)=inavector0(:,:)

!  call blockconfigmult(inavector0,jacaatemp,yyy%cptr(0))

  if (allspinproject==1) then
     call configspin_projectall(inavector, 0)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict(inavector, numr)
  endif

  call sparseconfigmult(inavector, jacaatemp, yyy%cptr(0), yyy%sptr(0), 1,1,0,1,0d0)

  if (allspinproject==1) then
     call configspin_projectall(jacaatemp, 0)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict(jacaatemp, numr)
  endif

  if (whichquad.eq.0) then
     outavector= jacaatemp - jacexpect*inavector - ( dot(jacaa,jacaatemp,totadim) + &
          dot(inavector,jacaamult,totadim) ) * jacaa 
  else
     outavector= jacaatemp - quadexpect*inavector 
  endif
!  quadcalled=quadcalled+1
!  if (mod(quadcalled,1000).eq.0) then
!     OFLWR "    Avector Quad - iterating!!  Is your maxaorder high enough??? "
!     WRFL "         iterations ", quadcalled, " maxaorder=",maxaorder; CFL
!  endif


end subroutine aamult




subroutine aaonedinit(inavector) 
  use aaonedmod
  use xxxmod
  use parameters
  implicit none
  real*8 :: norm
  DATATYPE :: inavector(numconfig,numr),  dot 

!  quadcalled=0
  if (sparseconfigflag.eq.0) then
     allocate(bigconfigmatel(totadim,totadim), tempconfigmatel(totadim,totadim), jacaaproj(totadim,totadim), jacaaproj2(totadim,totadim), realconfigmatel(zzz,totadim,zzz,totadim), ipiv(zzz*totadim), err(totadim))
  endif
  allocate(jacaa(numconfig,numr), jacaamult(numconfig,numr))

  jacaa=inavector
  if (allspinproject==1) then
     call configspin_projectall(jacaa,1)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict(jacaa,numr)
  endif
  norm=abs(sqrt(dot(jacaa,jacaa,totadim)))
  if (abs(norm)-1.d0.gt.1.d-5) then
     print *, "warning large norm error  avec ", norm, abs(norm)
  endif

  call sparseconfigmult(jacaa, jacaamult, yyy%cptr(0), yyy%sptr(0),1,1,0,1,0d0)
  if (allspinproject==1) then
     call configspin_projectall(jacaamult, 1)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict(jacaamult, numr)
  endif

  jacexpect=dot(jacaa,jacaamult,totadim);  quadexpect=dot(jacaa,jacaamult,totadim)/dot(jacaa,jacaa,totadim)

end subroutine aaonedinit

subroutine aaonedfinal
  use aaonedmod
  use parameters
  implicit none
  deallocate(jacaa, jacaamult) 
  if (sparseconfigflag.eq.0) then
     deallocate(bigconfigmatel, tempconfigmatel, jacaaproj, jacaaproj2, realconfigmatel, ipiv, err)
  endif
end subroutine

function  fquadenergy()
  use parameters
  use aaonedmod
  implicit none
  DATATYPE :: fquadenergy
  fquadenergy=quadexpect
end function fquadenergy

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

!!  OFLWR "ITERATIONS ***  ",jjcalls

     endif
  enddo


end subroutine quadavector

subroutine sparsequadavector(inavector,jjcalls0)
  use parameters
  use mpimod
  use aaonedmod
  use xxxmod
  implicit none

  EXTERNAL :: paraamult_transpose, parquadpreconsub_transpose
  EXTERNAL :: paraamult_transpose_spin, parquadpreconsub_transpose_spin
  EXTERNAL :: paraamult_transpose_padded, parquadpreconsub_transpose_padded
  EXTERNAL :: paraamult_transpose_spin_padded, parquadpreconsub_transpose_spin_padded

  integer :: jjcalls, ss,ii,jjcalls0
  real*8 ::  dev,  thisaerror
  DATATYPE :: dot, hermdot, inavector(numconfig,numr)
  DATATYPE, allocatable ::  vector(:,:), vector2(:,:), vector3(:,:), vectortr(:,:), vector2tr(:,:), vector3tr(:,:), &
       vector2spin(:,:), vector2trspin(:,:),vector3trspin(:,:), vectortrspin(:,:), vectorspin(:,:), vector3spin(:,:), &
       smallvectortr(:,:),smallvectortrspin(:,:),smallvectortr2(:,:),smallvectortrspin2(:,:)
  CNORMTYPE :: norm

  jjcalls0=0

  allocate(smallvectortr(numr,botwalk:botwalk+maxconfigsperproc-1),smallvectortrspin(numr,spinstart:spinstart+maxspinrank-1), &
       smallvectortr2(numr,botwalk:botwalk+maxconfigsperproc-1),smallvectortrspin2(numr,spinstart:spinstart+maxspinrank-1))

  allocate(       vector(numconfig,numr), vector2(numconfig,numr), vector3(numconfig,numr),vector2spin(spintotrank,numr), &
       vectorspin(spintotrank,numr), vector3spin(spintotrank,numr))

  allocate(       vectortr(numr,numconfig), vector2tr(numr,numconfig), vector3tr(numr,numconfig))
  allocate(       vectortrspin(numr,spintotrank), vector2trspin(numr,spintotrank), vector3trspin(numr,spintotrank))

  vector=inavector

#define NEWFLAG 1

if (NEWFLAG==1) then

  if (allspinproject.ne.0) then
     ii=maxspinrank*nprocs*numr
  else
     ii=maxconfigsperproc*nprocs*numr
  endif

else
  if (allspinproject.ne.0) then
     ii=spintotrank*numr
  else
     ii=totadim
  endif
endif

  if (maxaorder.gt.ii) then
     maxaorder=ii
  endif
  if (aorder.gt.maxaorder) then
     aorder=maxaorder
  endif
  if (aorder.lt.3) then
     OFLWR "Error, you must have like no configs, use nonsparse fool!!"; CFLST
  endif

  ss=0
333 continue

  if (whichquad.eq.0) then
     OFLWR "Something is screwy with sparsequadavector.  Dev not being calculated correctly. Run with whichquad=1 not zero. TEMP CONTINUE"; CFL    !!!  CFLST
  endif
  
  if (allspinproject==1) then
     call configspin_projectall(vector, 1)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict(vector, numr)
  endif

  call aaonedinit(vector)

  call sparseconfigmult(vector, vector2, yyy%cptr(0), yyy%sptr(0), 1,1,0,1,0d0)

  if (allspinproject.ne.0) then
     call configspin_projectall(vector2,0)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict(vector2,numr)
  endif

  if (allspinproject.ne.0) then
     call configspin_transformto(numr,vector2,vector2spin)
  endif


  if (whichquad.eq.0) then   !! mine
     vector3=vector2-jacexpect*vector              !! error term.
  else   !! rayleigh
     vector3=vector2-quadexpect*vector              !! error term.
  endif

  if (allspinproject.ne.0) then
     call configspin_transformto(numr,vector,vectorspin)
  endif

  dev=abs(sqrt(hermdot(vector3,vector3,totadim)))

!  OFLWR "     SPARSEQUAD: DEV,EXPECT", dev,quadexpect,jacexpect; CFL
  OFLWR "     SPARSEQUAD: DEV", dev; CFL

  thisaerror=min(0.1d0,max(sqrt(aerror),aerror/dev))

!061814
!  thisaerror=aerror/min(dev,1d0)

  if (abs(dev).lt.aerror.or.ss.ge.10) then
     inavector = vector
     if (ss.ge.10) then
        OFLWR "10 iterations, quad aborting"; CFL
     endif

     !OFLWR "     Converged sparse quad: dev ", dev, "  Expect ", quadexpect, " tol is ", thisaerror; CFL

     call aaonedfinal()
     deallocate(smallvectortr,smallvectortr2,smallvectortrspin,smallvectortrspin2)
     deallocate(vector, vector2, vector3,vector2spin,vectorspin,vector3spin,vectortr, vector2tr, vector3tr,vectortrspin,vector2trspin,vector3trspin)
     return
  endif

  vector3=0d0; 

  if (whichquad.eq.0) then

     if (allspinproject.ne.0) then

        vector2trspin(:,:)=TRANSPOSE(vector2spin(:,:))

if (NEWFLAG==0) then
        call dgsolve0( vector2trspin(:,spinstart), vector3trspin(:,spinstart), jjcalls, paraamult_transpose_spin,quadprecon,parquadpreconsub_transpose_spin, thisaerror,numr*(spinend-spinstart+1),maxaorder,1)
else
        smallvectortrspin(:,:)=0d0; smallvectortrspin(:,spinstart:spinend)=vector2trspin(:,spinstart:spinend)
        call dgsolve0( smallvectortrspin(:,:), smallvectortrspin2(:,:), jjcalls, paraamult_transpose_spin_padded,quadprecon,parquadpreconsub_transpose_spin_padded, thisaerror,numr*maxspinrank,maxaorder,1)
        vector3trspin(:,:)=0d0;  vector3trspin(:,spinstart:spinend)=smallvectortrspin2(:,spinstart:spinend)
endif

        call mpiallgather(vector3trspin(:,:),spintotrank*numr,allspinranks*numr,maxspinrank*numr)
        vector3spin=TRANSPOSE(vector3trspin)
        call configspin_transformfrom(numr,vector3spin,vector3)

     else

        vector2tr(:,:)=TRANSPOSE(vector2(:,:))

if (NEWFLAG==0) then
        call dgsolve0( vector2tr(:,botwalk), vector3tr(:,botwalk), jjcalls, paraamult_transpose,quadprecon,parquadpreconsub_transpose, thisaerror,numr*(topwalk-botwalk+1),maxaorder,1)
else
        smallvectortr(:,:)=0d0; smallvectortr(:,botwalk:topwalk)=vector2tr(:,botwalk:topwalk)
        call dgsolve0( smallvectortr(:,:), smallvectortr2(:,:), jjcalls, paraamult_transpose_padded,quadprecon,parquadpreconsub_transpose_padded, thisaerror,numr*maxconfigsperproc,maxaorder,1)
        vector3tr(:,:)=0d0; vector3tr(:,botwalk:topwalk)=smallvectortr2(:,botwalk:topwalk)
endif

        call mpiallgather(vector3tr(:,:),numconfig*numr,configsperproc*numr,maxconfigsperproc*numr)
        vector3=TRANSPOSE(vector3tr)

     endif

     vector=vector - vector3

  else

     if (allspinproject.ne.0) then

        vectortrspin(:,:)=TRANSPOSE(vectorspin(:,:))

if (NEWFLAG==0) then
        call dgsolve0( vectortrspin(:,spinstart), vector3trspin(:,spinstart), jjcalls, paraamult_transpose_spin,quadprecon,parquadpreconsub_transpose_spin, thisaerror,numr*(spinend-spinstart+1),maxaorder,1)
else
        smallvectortrspin(:,:)=0d0; smallvectortrspin(:,spinstart:spinend)=vectortrspin(:,spinstart:spinend)
        call dgsolve0( smallvectortrspin(:,:), smallvectortrspin2(:,:), jjcalls, paraamult_transpose_spin_padded,quadprecon,parquadpreconsub_transpose_spin_padded, thisaerror,numr*maxspinrank,maxaorder,1)
        vector3trspin(:,:)=0d0; vector3trspin(:,spinstart:spinend)=smallvectortrspin2(:,spinstart:spinend)
endif

        call mpiallgather(vector3trspin(:,:),spintotrank*numr,allspinranks*numr,maxspinrank*numr)
        vector3spin=TRANSPOSE(vector3trspin)
        call configspin_transformfrom(numr,vector3spin,vector3)
        
     else
        
        vectortr(:,:)=TRANSPOSE(vector(:,:))

if (NEWFLAG==0) then
        call dgsolve0( vectortr(:,botwalk), vector3tr(:,botwalk), jjcalls, paraamult_transpose,quadprecon,parquadpreconsub_transpose, thisaerror,numr*(topwalk-botwalk+1),maxaorder,1)
else
        smallvectortr(:,:)=0d0; smallvectortr(:,botwalk:topwalk)=vectortr(:,botwalk:topwalk)
        call dgsolve0( smallvectortr(:,:), smallvectortr2(:,:), jjcalls, paraamult_transpose_padded,quadprecon,parquadpreconsub_transpose_padded, thisaerror,numr*maxconfigsperproc,maxaorder,1)
        vector3tr(:,:)=0d0; vector3tr(:,botwalk:topwalk)=smallvectortr2(:,botwalk:topwalk)
endif

        call mpiallgather(vector3tr(:,:),numconfig*numr,configsperproc*numr,maxconfigsperproc*numr)
        vector3=TRANSPOSE(vector3tr)
        
     endif

     vector=vector3/sqrt(dot(vector3,vector3,totadim))

  endif

!!  OFLWR "    ITERATIONS ***  ",jjcalls

  jjcalls0=jjcalls0+jjcalls

  if (allspinproject==1) then
     call configspin_projectall(vector, 1)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict(vector, numr)
  endif

  norm=sqrt(dot(vector,vector,totadim)) ! ok conv
  vector=vector/norm

#ifndef CNORMFLAG
  vector=vector*abs(sum(vector))/sum(vector)
#endif

  call aaonedfinal()

  ss=ss+1
  go to 333

end subroutine sparsequadavector


subroutine nonsparsequadavector(avectorout)
  use mpimod
  use aaonedmod
  use parameters
  use xxxmod
  implicit none

  DATATYPE :: avectorout(totadim), dot, hermdot, flatjacaa(totadim),flatjacaamult(totadim)
  CNORMTYPE :: norm, dev  !! making dev hermdot now 030912
  integer :: k,i, info,ss
!! spin
  real*8, allocatable :: realhalfmatel(:,:,:,:),realhalfmatel2(:,:,:,:),realspinmatel(:,:,:,:)
  DATATYPE, allocatable :: crealconfigmatel(:,:,:),crealhalfmatel(:,:,:),crealhalfmatel2(:,:,:),&
       crealspinmatel(:,:,:),spinerr(:),spinavectorout(:)

  if (sparseconfigflag/=0) then
     OFLWR "Error, nonsparsequadavector called but sparseconfigflag= ", sparseconfigflag; CFLST
  endif

  ss=0
333 continue

  if (allspinproject==1) then
     call configspin_projectall(avectorout, 1)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict(avectorout, numr)
  endif

  call aaonedinit(avectorout)
  call sparseconfigmult(avectorout(:),err(:),yyy%cptr(0),yyy%sptr(0),1,1,0,1,0d0)

  if (allspinproject==1) then
     call configspin_projectall(err(:), 0)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict(err(:), numr)
  endif
  if (whichquad.eq.0) then
     err(:)=err(:)-jacexpect*avectorout(:)              !! error term.
  else
     err(:)=err(:)-quadexpect*avectorout(:)              !! error term.
  endif

  dev=abs(sqrt(hermdot(err(:),err(:),totadim)))

!!  OFLWR "  NONSPARSEQUAD: DEV,EXPECT", dev,quadexpect,jacexpect; CFL
  OFLWR "  NONSPARSEQUAD: DEV", dev; CFL

  if (abs(dev).lt.aerror.or.(ss.gt.10)) then

!#ifdef REALGO
!  OFL; write(mpifileptr,'(A40,F15.10, A5,1F18.10,A5,1E6.1)') "      Converged nonsparse quad: dev ", dev, "  E= ", quadexpect, " tol ", aerror; CFL
!#else
!  OFL; write(mpifileptr,'(A40,F15.10, A5,2F18.10,A5,1E6.1)') "      Converged nonsparse quad: dev ", dev, "  E= ", quadexpect, " tol ", aerror; CFL
!#endif

     call aaonedfinal()
     return
  endif

  flatjacaa(:)=RESHAPE(jacaa,(/totadim/))
  flatjacaamult(:)=RESHAPE(jacaamult,(/totadim/))

  do k=1,totadim
     jacaaproj(:,k)=flatjacaa(:)*CONJUGATE(flatjacaa(k))
  enddo

!! this operates on conjugate of vector
  do k=1,totadim
     jacaaproj2(:,k)=flatjacaamult(:)*(flatjacaa(k))
  enddo

!! no spin project; doing transform below
  call assemble_configmat(bigconfigmatel, yyy%cptr(0),1,1,1,1, 0d0)

  if (whichquad.eq.0) then
#ifdef REALGO
     OFLWR "CHECKME ZGEMM"; CFLST
#endif
     call zgemm('N', 'N', totadim, totadim, totadim, (1d0,0d0), jacaaproj, &
       totadim, bigconfigmatel, totadim, (0d0,0d0), tempconfigmatel, totadim)
     bigconfigmatel=bigconfigmatel - tempconfigmatel
  endif
  if (whichquad.eq.0) then
  do k=1,totadim
     bigconfigmatel(k,k)=bigconfigmatel(k,k)-jacexpect
  enddo
  else
  do k=1,totadim
     bigconfigmatel(k,k)=bigconfigmatel(k,k)-quadexpect
  enddo
  endif

  realconfigmatel=0d0
  do i=1,totadim
     realconfigmatel(1,:,1,i) = real(bigconfigmatel(:,i),8)
#ifndef REALGO     
     realconfigmatel(2,:,2,i) = real(bigconfigmatel(:,i),8)
     realconfigmatel(1,:,2,i) = (-1)*imag(bigconfigmatel(:,i))
     realconfigmatel(2,:,1,i) = imag(bigconfigmatel(:,i))
#endif
  enddo

  if (whichquad.eq.0) then
!! jacaaproj2 operates on conjugate of avectorin.
     do i=1,totadim
        realconfigmatel(1,:,1,i) =realconfigmatel(1,:,1,i) - real(jacaaproj2(:,i),8)
#ifndef REALGO
        realconfigmatel(2,:,2,i) =realconfigmatel(2,:,2,i)- (-1)*real(jacaaproj2(:,i),8)
        realconfigmatel(1,:,2,i) =realconfigmatel(1,:,2,i)- imag(jacaaproj2(:,i))
        realconfigmatel(2,:,1,i) =realconfigmatel(2,:,1,i)- imag(jacaaproj2(:,i))
#endif
     enddo
  endif



  if (allspinproject.ne.0) then

!! can make this complex or whatever; multiplying by real matrix
!!   but doing it carefully, only condensing leading index as complex

allocate(crealconfigmatel(numconfig*numr,zzz,numconfig*numr), &
     crealhalfmatel(spintotrank*numr,zzz,numconfig*numr), &
     realhalfmatel(zzz,spintotrank*numr,zzz,numconfig*numr), &
     realhalfmatel2(zzz,numconfig*numr,zzz,spintotrank*numr), &
     crealhalfmatel2(numconfig*numr,zzz,spintotrank*numr), &
     crealspinmatel(spintotrank*numr,zzz,spintotrank*numr), &
     realspinmatel(zzz,spintotrank*numr,zzz,spintotrank*numr), &
     spinerr(spintotrank*numr),spinavectorout(spintotrank*numr))

#ifdef REALGO
     crealconfigmatel=realconfigmatel(1,:,:,:)
#else
     crealconfigmatel=realconfigmatel(1,:,:,:) + realconfigmatel(2,:,:,:) * (0d0,1d0)
#endif

     call configspin_transformto(numconfig*numr*numr*zzz,crealconfigmatel,crealhalfmatel)

     realhalfmatel(1,:,:,:) = real( crealhalfmatel(:,:,:) ,8)

#ifndef REALGO
     realhalfmatel(2,:,:,:) = imag( crealhalfmatel(:,:,:) )
#endif

     realhalfmatel2(:,:,:,:)=RESHAPE(TRANSPOSE(RESHAPE(realhalfmatel(:,:,:,:),(/zzz*spintotrank*numr,zzz*numconfig*numr/))),(/zzz,numconfig*numr,zzz,spintotrank*numr/))

#ifdef REALGO
     crealhalfmatel2(:,:,:)=realhalfmatel2(1,:,:,:)
#else
     crealhalfmatel2(:,:,:)=realhalfmatel2(1,:,:,:) + realhalfmatel2(2,:,:,:) * (0d0,1d0) 
#endif

     call configspin_transformto(spintotrank*numr*numr*zzz,crealhalfmatel2,crealspinmatel)

     realspinmatel(1,:,:,:) = real( crealspinmatel(:,:,:) ,8)
#ifndef REALGO
     realspinmatel(2,:,:,:) = imag( crealspinmatel(:,:,:))
#endif

     realspinmatel(:,:,:,:)=RESHAPE(TRANSPOSE(RESHAPE(realspinmatel(:,:,:,:),(/zzz*spintotrank*numr,zzz*spintotrank*numr/))),(/zzz,spintotrank*numr,zzz,spintotrank*numr/))

     call configspin_transformto(numr,avectorout(:),spinavectorout(:))
     call configspin_transformto(numr,err(:),spinerr(:))

     if (myrank.eq.1) then
        if (whichquad.eq.0) then
           call dgesv(spintotrank*numr*zzz,1,realspinmatel,spintotrank*numr*zzz,ipiv,spinerr,spintotrank*numr*zzz,info)
           spinavectorout=spinavectorout-spinerr
        else
           call dgesv(spintotrank*numr*zzz,1,realspinmatel,spintotrank*numr*zzz,ipiv,spinavectorout,spintotrank*numr*zzz,info)
           spinavectorout=spinavectorout/sqrt(dot(spinavectorout,spinavectorout,spintotrank*numr))
        endif
        if (info/=0) then
           OFLWR "Info in dgesv nonsparsequadavector= ", info; CFLST
        endif
     endif

     call mympibcast(spinavectorout,1,spintotrank*numr)

     call configspin_transformfrom(numr,spinavectorout(:),avectorout(:))

     deallocate(crealconfigmatel,     crealhalfmatel,     realhalfmatel,     realhalfmatel2, &
          crealhalfmatel2,     crealspinmatel,      realspinmatel,     spinerr, spinavectorout)

  else

     if (myrank.eq.1) then
        if (whichquad.eq.0) then
           call dgesv(totadim*zzz,1,realconfigmatel,totadim*zzz,ipiv,err,totadim*zzz,info)
           avectorout=avectorout-err
        else
           call dgesv(totadim*zzz,1,realconfigmatel,totadim*zzz,ipiv,avectorout,totadim*zzz,info)
           avectorout=avectorout/sqrt(dot(avectorout,avectorout,totadim))
        endif
        if (info/=0) then
           OFLWR "Info in dgesv nonsparsequadavector= ", info; CFLST
        endif
     endif

     call mympibcast(avectorout,1,totadim)

  endif


  if (allspinproject==1) then
     call configspin_projectall(avectorout, 1)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict(avectorout, numr)
  endif

  norm=sqrt(dot(avectorout,avectorout,totadim))  !ok imp conv ch/pmctdh
  avectorout=avectorout/norm
#ifndef CNORMFLAG
  avectorout=avectorout*abs(sum(avectorout))/sum(avectorout)
#endif

  call aaonedfinal()
  ss=ss+1
  go to 333
  
end subroutine nonsparsequadavector


subroutine paraamult_transpose_spin_padded(nullint,inavectortrspin,outavectortrspin)
  use parameters
  implicit none

  integer :: nullint
  DATATYPE :: inavectortrspin(numr,maxspinrank), outavectortrspin(numr,maxspinrank)

  outavectortrspin(:,:)=0d0   !!inavectortrspin(:,:)
  call paraamult_transpose_spin(nullint,inavectortrspin,outavectortrspin)
end subroutine paraamult_transpose_spin_padded


subroutine paraamult_transpose_spin(nullint,inavectortrspin,outavectortrspin)
  use parameters
  implicit none

  integer :: nullint
  DATATYPE :: inavectortrspin(numr,spinstart:spinend), outavectortrspin(numr,spinstart:spinend)
  DATATYPE :: inavectortr(numr,botwalk:topwalk), outavectortr(numr,botwalk:topwalk)

!  DATATYPE :: inavectorspin(spinstart:spinend,numr),inavector(botwalk:topwalk,numr),outavector(botwalk:topwalk,numr)
!  DATATYPE :: outavectorspin(spinstart:spinend,numr)

  if (spinwalkflag.eq.0) then
     OFLWR "WTF SPIN"; CFLST
  endif

!  inavectorspin(:,:)=TRANSPOSE(inavectortrspin(:,spinstart:spinend))
!  call configspin_transformfrom_mine(numr,inavectorspin,inavector)
!  inavectortr(:,:)=TRANSPOSE(inavector(:,:))

  call configspin_transformfrom_mine_transpose(inavectortrspin,inavectortr)

  call paraamult_transpose(0,inavectortr,outavectortr)

  outavectortrspin(:,:)=0d0

!  outavector(:,:)=TRANSPOSE(outavectortr(:,:))
!  call configspin_transformto_mine(numr,outavector,outavectorspin)
!  outavectortrspin(:,spinstart:spinend)=TRANSPOSE(outavectorspin(:,:))

  call configspin_transformto_mine_transpose(outavectortr,outavectortrspin)
  
end subroutine paraamult_transpose_spin



subroutine paraamult_transpose_padded(nullint, inavectortr,outavectortr)
  use parameters
  use xxxmod
  use aaonedmod
  implicit none
  integer :: nullint
  DATATYPE :: inavectortr(numr,maxconfigsperproc), outavectortr(numr,maxconfigsperproc)

  outavectortr(:,:)=0d0    !!inavectortr(:,:)
  call paraamult_transpose(nullint, inavectortr,outavectortr)
end subroutine paraamult_transpose_padded


subroutine paraamult_transpose(nullint, inavectortr,outavectortr)
  use parameters
  use xxxmod
  use aaonedmod
  implicit none
  integer :: nullint
  DATATYPE :: inavectortr(numr,botwalk:topwalk), outavectortr(numr,botwalk:topwalk), pardot , &
       jacaatemptr(numr,botwalk:topwalk),jacaalittletr(numr,botwalk:topwalk), jacaamultlittletr(numr,botwalk:topwalk)

  jacaalittletr(:,:)=TRANSPOSE(jacaa(botwalk:topwalk,:)); jacaamultlittletr(:,:)=TRANSPOSE(jacaamult(botwalk:topwalk,:))

  call parblockconfigmult_transpose(inavectortr,jacaatemptr)

  if (whichquad.eq.0) then

     outavectortr(:,:)= jacaatemptr(:,:) - jacexpect*inavectortr(:,:) - ( pardot(jacaalittletr,jacaatemptr,numr*(topwalk-botwalk+1)) + &
          pardot(inavectortr,jacaamultlittletr,numr*(topwalk-botwalk+1)) ) * jacaalittletr(:,:)
  else
     outavectortr= jacaatemptr - quadexpect*inavectortr 
  endif
!  quadcalled=quadcalled+1
!  if (mod(quadcalled,1000).eq.0) then
!     OFLWR "    Avector Quad - iterating!!  Is your maxaorder high enough??? "
!     WRFL "         iterations ", quadcalled, " maxaorder=",maxaorder; CFL
!  endif


end subroutine paraamult_transpose



subroutine parquadpreconsub_transpose_spin_padded(nullint,inavectortrspin,outavectortrspin)
  use parameters
  implicit none

  integer :: nullint
  DATATYPE :: inavectortrspin(numr,maxspinrank), outavectortrspin(numr,maxspinrank)

  outavectortrspin(:,:)=inavectortrspin(:,:)
  call parquadpreconsub_transpose_spin(nullint,inavectortrspin,outavectortrspin)
end subroutine parquadpreconsub_transpose_spin_padded


subroutine parquadpreconsub_transpose_spin(nullint,inavectortrspin,outavectortrspin)
  use parameters
  implicit none

  integer :: nullint
  DATATYPE :: inavectortrspin(numr,spinstart:spinend), outavectortrspin(numr,spinstart:spinend)
  DATATYPE :: inavectortr(numr,botwalk:topwalk),outavectortr(numr,botwalk:topwalk)

  if (spinwalkflag.eq.0) then
     OFLWR "WTF SPIN"; CFLST
  endif

  call configspin_transformfrom_mine_transpose(inavectortrspin,inavectortr)

  call parquadpreconsub_transpose(0,inavectortr,outavectortr)

  call configspin_transformto_mine_transpose(outavectortr,outavectortrspin)

end subroutine parquadpreconsub_transpose_spin


subroutine parquadpreconsub_transpose_padded(nullint, inavectortr,outavectortr)
  use parameters
  use aaonedmod
  use xxxmod
  implicit none
  integer :: nullint
  DATATYPE :: inavectortr(numr,maxconfigsperproc), outavectortr(numr,maxconfigsperproc)

  outavectortr(:,:)=inavectortr(:,:)
  call parquadpreconsub_transpose(nullint, inavectortr,outavectortr)

end subroutine parquadpreconsub_transpose_padded


subroutine parquadpreconsub_transpose(nullint, inavectortr,outavectortr)
  use parameters
  use aaonedmod
  use xxxmod
  implicit none
  integer :: nullint
  DATATYPE :: inavectortr(numr,botwalk:topwalk), outavectortr(numr,botwalk:topwalk),tempvector(botwalk:topwalk,numr), &
       tempvectortr(numr,botwalk:topwalk),tempvectortr2(numr,botwalk:topwalk),outavector(botwalk:topwalk,numr)

  tempvector(:,:)=TRANSPOSE(inavectortr(:,:))

  if (allspinproject==1) then
     call configspin_projectmine(tempvector,0)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict_par(tempvector, numr)
  endif

  tempvectortr(:,:)=TRANSPOSE(tempvector(:,:))

  if (whichquad.eq.0) then
     call parsparseconfigdiagmult_transpose(tempvectortr, tempvectortr2, yyy%cptr(0), yyy%sptr(0),1,1,1,1, jacexpect,0d0) !+quadpreconshift,0d0)
  else
     call parsparseconfigdiagmult_transpose(tempvectortr, tempvectortr2, yyy%cptr(0), yyy%sptr(0),1,1,1,1, quadexpect,0d0) !+quadpreconshift,0d0)
  endif
  outavector=TRANSPOSE(tempvectortr2)
  if (allspinproject==1) then
     call configspin_projectmine(outavector, 0)
  endif
  if (dfrestrictflag.ne.0) then
     call dfrestrict_par(outavector, numr)
  endif
  outavectortr(:,:)=TRANSPOSE(outavector)

end subroutine


