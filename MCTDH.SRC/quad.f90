

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
  rgwkdim= (42 + jjxx*(inkrydim+6) + inkrydim*(inkrydim+3)) / 5 * 6

  if (rgwkdim.lt.0) then
     OFLWR "OOPS, integer overflow, decrease maxexpodim or reduce size of orbitals (per processor)",rgwkdim; CFLST
  endif

  allocate(rgwk(rgwkdim)); rgwk(:)=0d0

  itol=0;  iunit=0

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

!!$  if (debugflag.ne.0) then
!!$     OFLWR "        DGSOLVE - restarts ",igwk(5), "iters", numiter; CFL
!!$  endif

  if (jacwrite==1) then
     OFLWR "   NEWTON: numiter ", numiter, " error ", errest;     call closefile()
  endif
  if (ierr.eq.1) then     !! AUTHOR COMMENTS ARE WRONG IN DGMRES.F ierr=1 denotes max restarts (zero) reached !! 
     OFLWR "WARNING! dgmres did not converge. increase krylov order?"; CFL
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
  use linearmod
  implicit none
  EXTERNAL :: quadoperate,dummysub
  integer :: jjcalls,icount
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
        call dgsolve0( vector2, vector3, jjcalls, quadoperate,0,dummysub, quadtol,totspfdim,maxexpodim,1)
     else
        call dgsolve0( vector2, vector3, jjcalls, quadoperate,0,dummysub, quadtol,totspfdim,maxexpodim,0)
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
  integer, parameter :: ireduced=1
  DATATYPE :: quadexpect=0d0
  DATATYPE, allocatable :: jacaa(:,:), jacaamult(:,:)
  DATATYPE, allocatable :: bigconfigmatel(:,:), jacaaproj(:,:), jacaaproj2(:,:), err(:)
  real*8, allocatable :: realconfigmatel(:,:,:,:)
  integer, allocatable :: ipiv(:)
#ifdef REALGO
  integer, parameter :: zzz=1
#else
  integer, parameter :: zzz=2
#endif
end module


recursive subroutine quadpreconsq(nullint, inavector,outavector)
  use parameters
  use aaonedmod
  implicit none
  integer :: nullint
  DATATYPE :: inavector(totadim), outavector(totadim),jtempaa(totadim)
  call quadpreconsub(nullint, inavector,jtempaa);  jtempaa=conjg((0d0,0d0)+jtempaa)
  call quadpreconsub(nullint, jtempaa,outavector);  outavector=conjg((0d0,0d0)+outavector)
end subroutine quadpreconsq

recursive subroutine quadpreconsub(notusedint, inavector,outavector)
  use parameters
  use aaonedmod
  use xxxmod
  implicit none
  integer :: notusedint
  DATATYPE :: inavector(totadim), outavector(totadim),&
       tempvector(totadim)  !!AUTOMATIC

  tempvector=inavector

  if (allspinproject==1) then
     call configspin_project(tempvector, 0)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(tempvector, numr)
  endif

  call sparseconfigdiagmult(tempvector, outavector, yyy%cptr(0),yyy%sptr(0), 1,1,1,1, quadexpect,0d0)

  if (allspinproject==1) then
     call configspin_project(outavector, 0)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(outavector, numr)
  endif
end subroutine



recursive subroutine aamult(notusedint, inavector0,outavector)
  use parameters
  use xxxmod
  use aaonedmod
  implicit none
  integer :: notusedint
  DATATYPE, intent(in) :: inavector0(numr,firstconfig:lastconfig)
  DATATYPE, intent(out) :: outavector(numr,firstconfig:lastconfig)
  DATATYPE :: inavector(numr,firstconfig:lastconfig), jacaatemp(numr,firstconfig:lastconfig)

  inavector(:,:)=inavector0(:,:)

  if (allspinproject==1) then
     call configspin_project(inavector, 0)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(inavector, numr)
  endif
  call sparseconfigmult(inavector, jacaatemp, yyy%cptr(0), yyy%sptr(0), 1,1,0,1,0d0)
  if (allspinproject==1) then
     call configspin_project(jacaatemp, 0)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(jacaatemp, numr)
  endif

  outavector= jacaatemp - quadexpect*inavector 

end subroutine aamult


subroutine aaonedinit(inavector) 
  use aaonedmod
  use xxxmod
  use parameters
  implicit none
  real*8 :: norm
  DATATYPE :: inavector(numr,firstconfig:lastconfig),  dot, csum1, csum2

  if (sparseconfigflag.eq.0) then
     allocate(bigconfigmatel(totadim,totadim), jacaaproj(totadim,totadim), jacaaproj2(totadim,totadim), realconfigmatel(zzz,totadim,zzz,totadim), ipiv(zzz*totadim), err(totadim))
  endif
  allocate(jacaa(numr,firstconfig:lastconfig), jacaamult(numr,firstconfig:lastconfig))

  jacaa(:,:)=inavector(:,:)
  if (allspinproject==1) then
     call configspin_project(jacaa,1)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(jacaa,numr)
  endif
  norm=abs(dot(jacaa,jacaa,totadim))
  if (parconsplit.ne.0) then
     call mympirealreduceone(norm)
  endif
  norm=sqrt(norm)

  if (abs(norm)-1.d0.gt.1.d-5) then
     print *, "warning large norm error  avec ", norm, abs(norm)
  endif

  call sparseconfigmult(jacaa, jacaamult, yyy%cptr(0), yyy%sptr(0),1,1,0,1,0d0)
  if (allspinproject==1) then
     call configspin_project(jacaamult, 1)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(jacaamult, numr)
  endif

  csum1=dot(jacaa,jacaamult,totadim)
  csum2=dot(jacaa,jacaa,totadim)
  if (parconsplit.ne.0) then
     call mympireduceone(csum1); call mympireduceone(csum2)
  endif

  quadexpect=csum1/csum2

end subroutine aaonedinit


subroutine aaonedfinal
  use aaonedmod
  use parameters
  implicit none
  deallocate(jacaa, jacaamult) 
  if (sparseconfigflag.eq.0) then
     deallocate(bigconfigmatel, jacaaproj, jacaaproj2, realconfigmatel, ipiv, err)
  endif
end subroutine


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
  EXTERNAL :: paraamult_padded, parquadpreconsub_padded
  EXTERNAL :: paraamult_spin_padded, parquadpreconsub_spin_padded
  DATATYPE, intent(inout) ::  inavector(numr,firstconfig:lastconfig)
  integer :: jjcalls, ss,ii,jjcalls0
  real*8 ::  dev,  thisaerror
  DATATYPE :: dot, hermdot,csum 
  DATATYPE, allocatable ::  vector(:,:), vector2(:,:), vector3(:,:), &
       smallvector(:,:),smallvectorspin(:,:),smallvector2(:,:),smallvectorspin2(:,:)
  CNORMTYPE :: norm

  jjcalls0=0

  allocate(smallvector(numr,botwalk:botwalk+maxconfigsperproc-1),smallvectorspin(numr,spinstart:spinstart+maxspinrank-1), &
       smallvector2(numr,botwalk:botwalk+maxconfigsperproc-1),smallvectorspin2(numr,spinstart:spinstart+maxspinrank-1))

  allocate( vector(numr,firstconfig:lastconfig), vector2(numr,firstconfig:lastconfig), vector3(numr,firstconfig:lastconfig))

  vector(:,:)=inavector(:,:)

  if (allspinproject.ne.0) then
     ii=maxspinrank*nprocs*numr
  else
     ii=maxconfigsperproc*nprocs*numr
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
     deallocate(smallvector,smallvector2,smallvectorspin,smallvectorspin2)
     deallocate(vector, vector2, vector3)
     return
  endif

  vector3=0d0; 


  if (allspinproject.ne.0) then
     
     smallvectorspin(:,:)=0d0
     call configspin_transformto_local(numr,vector(:,botwalk),smallvectorspin(:,spinstart))
     smallvectorspin2(:,:)=smallvectorspin(:,:)    !! guess
     call dgsolve0( smallvectorspin(:,:), smallvectorspin2(:,:), jjcalls, paraamult_spin_padded,quadprecon,parquadpreconsub_spin_padded, thisaerror,numr*maxspinrank,maxaorder,1)

     call configspin_transformfrom_local(numr,smallvectorspin2(:,:),smallvector2(:,:))
     
  else
     
     smallvector(:,:)=0d0; 
     smallvector(:,botwalk:topwalk)=vector(:,botwalk:topwalk)
     smallvector2(:,:)=smallvector(:,:)       !! guess
     call dgsolve0( smallvector(:,:), smallvector2(:,:), jjcalls, paraamult_padded,quadprecon,parquadpreconsub_padded, thisaerror,numr*maxconfigsperproc,maxaorder,1)

  endif

  vector3(:,:)=0d0; 
  vector3(:,botwalk:topwalk)=smallvector2(:,botwalk:topwalk)
  
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

  if (allspinproject==1) then
     call configspin_project(vector, 1)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(vector, numr)
  endif

  csum=dot(vector,vector,totadim)
  call mympireduceone(csum)
  norm=sqrt(csum)         ! ok conv
  vector=vector/norm

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
  DATATYPE,intent(inout) :: avectorout(totadim)
  DATATYPE :: dot, hermdot
  DATATYPE :: flatjacaa(totadim),flatjacaamult(totadim)  !! AUTOMATIC
  CNORMTYPE :: norm
  real*8 :: dev
  integer :: k,i, info,ss
  real*8, allocatable :: realhalfmatel(:,:,:,:),realhalfmatel2(:,:,:,:),realspinmatel(:,:,:,:)
  DATATYPE, allocatable :: crealconfigmatel(:,:,:),crealhalfmatel(:,:,:),crealhalfmatel2(:,:,:),&
       crealspinmatel(:,:,:),spinerr(:),spinavectorout(:)

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

  do k=1,totadim
     bigconfigmatel(k,k)=bigconfigmatel(k,k)-quadexpect
  enddo

  realconfigmatel=0d0
  do i=1,totadim
     realconfigmatel(1,:,1,i) = real(bigconfigmatel(:,i),8)
#ifndef REALGO     
     realconfigmatel(2,:,2,i) = real(bigconfigmatel(:,i),8)
     realconfigmatel(1,:,2,i) = (-1)*imag(bigconfigmatel(:,i))
     realconfigmatel(2,:,1,i) = imag(bigconfigmatel(:,i))
#endif
  enddo

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

     call configspin_transformto(numr*numconfig*numr*zzz,crealconfigmatel,crealhalfmatel2)

     realhalfmatel2(1,:,:,:) = real( crealhalfmatel2(:,:,:) ,8)

#ifndef REALGO
     realhalfmatel2(2,:,:,:) = imag( crealhalfmatel2(:,:,:) )
#endif

     realhalfmatel(:,:,:,:)=RESHAPE(TRANSPOSE(RESHAPE(realhalfmatel2(:,:,:,:),(/zzz*numconfig*numr,zzz*spintotrank*numr/))),(/zzz,spintotrank*numr,zzz,numconfig*numr/))

#ifdef REALGO
     crealhalfmatel(:,:,:)=realhalfmatel(1,:,:,:)
#else
     crealhalfmatel(:,:,:)=realhalfmatel(1,:,:,:) + realhalfmatel(2,:,:,:) * (0d0,1d0) 
#endif

     call configspin_transformto(numr*spintotrank*numr*zzz,crealhalfmatel,crealspinmatel)

     realspinmatel(1,:,:,:) = real( crealspinmatel(:,:,:) ,8)
#ifndef REALGO
     realspinmatel(2,:,:,:) = imag( crealspinmatel(:,:,:))
#endif

     realspinmatel(:,:,:,:)=RESHAPE(TRANSPOSE(RESHAPE(realspinmatel(:,:,:,:),(/zzz*spintotrank*numr,zzz*spintotrank*numr/))),(/zzz,spintotrank*numr,zzz,spintotrank*numr/))

     call configspin_transformto(numr,avectorout(:),spinavectorout(:))
     call configspin_transformto(numr,err(:),spinerr(:))

     if (myrank.eq.1) then
        call dgesv(spintotrank*numr*zzz,1,realspinmatel,spintotrank*numr*zzz,ipiv,spinavectorout,spintotrank*numr*zzz,info)
        spinavectorout=spinavectorout/sqrt(dot(spinavectorout,spinavectorout,spintotrank*numr))
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

        call dgesv(totadim*zzz,1,realconfigmatel,totadim*zzz,ipiv,avectorout,totadim*zzz,info)
        avectorout=avectorout/sqrt(dot(avectorout,avectorout,totadim))

        if (info/=0) then
           OFLWR "Info in dgesv nonsparsequadavector= ", info; CFLST
        endif
     endif

     call mympibcast(avectorout,1,totadim)

  endif


  if (allspinproject==1) then
     call configspin_project(avectorout, 1)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project(avectorout, numr)
  endif

  norm=sqrt(dot(avectorout,avectorout,totadim))  !ok conv
  avectorout=avectorout/norm

  call aaonedfinal()
  ss=ss+1
  go to 333
  
end subroutine nonsparsequadavector


recursive subroutine paraamult_spin_padded(nullint,inavectorspin,outavectorspin)
  use parameters
  implicit none
  integer :: nullint
  DATATYPE,intent(in) :: inavectorspin(numr,maxspinrank)
  DATATYPE,intent(out) :: outavectorspin(numr,maxspinrank)
  outavectorspin(:,:)=0d0
  call paraamult_spin(nullint,inavectorspin,outavectorspin)
end subroutine paraamult_spin_padded


recursive subroutine paraamult_spin(notusedint,inavectorspin,outavectorspin)
  use parameters
  implicit none
  integer :: notusedint
  DATATYPE,intent(in) :: inavectorspin(numr,spinstart:spinend)
  DATATYPE,intent(out) :: outavectorspin(numr,spinstart:spinend)
  DATATYPE :: inavector(numr,botwalk:topwalk), outavector(numr,botwalk:topwalk)
  if (spinwalkflag.eq.0) then
     OFLWR "WTF SPIN"; CFLST
  endif
  call configspin_transformfrom_local(numr,inavectorspin,inavector)
  call paraamult(0,inavector,outavector)
  call configspin_transformto_local(numr,outavector,outavectorspin)
end subroutine paraamult_spin


recursive subroutine paraamult_padded(nullint, inavector,outavector)
  use parameters
  implicit none
  integer :: nullint
  DATATYPE,intent(in) :: inavector(numr,maxconfigsperproc)
  DATATYPE,intent(out) :: outavector(numr,maxconfigsperproc)
  outavector(:,:)=0d0   
  call paraamult(nullint, inavector,outavector)
end subroutine paraamult_padded


recursive subroutine paraamult(notusedint, inavector,outavector)
  use parameters
  use aaonedmod
  implicit none
  integer :: notusedint
  DATATYPE,intent(in) :: inavector(numr,botwalk:topwalk)
  DATATYPE,intent(out) :: outavector(numr,botwalk:topwalk)

  call parblockconfigmult(inavector,outavector)

  outavector= outavector - quadexpect*inavector 

end subroutine paraamult


recursive subroutine parquadpreconsub_spin_padded(nullint,inavectorspin,outavectorspin)
  use parameters
  implicit none
  integer :: nullint
  DATATYPE,intent(in) :: inavectorspin(numr,maxspinrank)
  DATATYPE,intent(out) :: outavectorspin(numr,maxspinrank)
  outavectorspin(:,:)=inavectorspin(:,:)  !! MUST BE FULL RANK
  call parquadpreconsub_spin(nullint,inavectorspin,outavectorspin)
end subroutine parquadpreconsub_spin_padded


recursive subroutine parquadpreconsub_spin(notusedint,inavectorspin,outavectorspin)
  use parameters
  implicit none
  integer :: notusedint
  DATATYPE,intent(in) :: inavectorspin(numr,spinstart:spinend)
  DATATYPE,intent(out) :: outavectorspin(numr,spinstart:spinend)
  DATATYPE :: inavector(numr,botwalk:topwalk),outavector(numr,botwalk:topwalk)
  if (spinwalkflag.eq.0) then
     OFLWR "WTF SPIN"; CFLST
  endif
  call configspin_transformfrom_local(numr,inavectorspin,inavector)
  call parquadpreconsub(0,inavector,outavector)
  call configspin_transformto_local(numr,outavector,outavectorspin)
end subroutine parquadpreconsub_spin


recursive subroutine parquadpreconsub_padded(nullint, inavector,outavector)
  use parameters
  implicit none
  integer :: nullint
  DATATYPE,intent(in) :: inavector(numr,maxconfigsperproc)
  DATATYPE,intent(out) :: outavector(numr,maxconfigsperproc)
  outavector(:,:)=inavector(:,:) !! MUST BE FULL RANK
  call parquadpreconsub(nullint, inavector,outavector)
end subroutine parquadpreconsub_padded


recursive subroutine parquadpreconsub(notusedint, inavector,outavector)
  use parameters
  use aaonedmod
  use xxxmod
  implicit none
  integer :: notusedint
  DATATYPE,intent(in) :: inavector(numr,botwalk:topwalk)
  DATATYPE,intent(out) :: outavector(numr,botwalk:topwalk)
  DATATYPE :: tempvector(numr,botwalk:topwalk)

  tempvector(:,:)=inavector(:,:)

  if (allspinproject==1) then
     call configspin_project_local(tempvector,0)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project_local(tempvector, numr)
  endif

  call parsparseconfigdiagmult(tempvector, outavector, yyy%cptr(0), yyy%sptr(0),1,1,1,1, quadexpect,0d0)

  if (allspinproject==1) then
     call configspin_project_local(outavector, 0)
  endif
  if (dfrestrictflag.ne.0) then
     call df_project_local(outavector, numr)
  endif

end subroutine parquadpreconsub
