
!! BIORTHOGONALZATION ROUTINES


#include "Definitions.INC"

!! this is the bad boy that does biorthonormalization
!! it gets fed in all the parameters it way need to eschew parameters and the mods
!! pointers are used to save on memory slightly for the external routine needed in expokit 
!! and the sparse algorithm     they are instanced only once for the call to abio_sparse and then dropped 
!! so we don't leave any unwanted dangling pointers
!! the main thing is verifying that biortho's self-re-initializer is functioning properly
!! it should be for all electron flux subroutine needs, but in case of added functionality that uses 
!! this routine double check this    kvl

module biorthotypemod
  type biorthotype
     DATATYPE, pointer :: smo(:,:)
     integer :: thisbiodim=100000000, biomaxdim=-1, bionr=-1
     logical :: started=.false.
  end type biorthotype
end module biorthotypemod

module biomatvecmod
  use biorthotypemod
  implicit none
  type(biorthotype), pointer :: biopointer
end module

module matvecsetmod
contains
subroutine biomatvecset(inbiotarget)
  use biomatvecmod
  implicit none
  type(biorthotype), target :: inbiotarget
  biopointer=>inbiotarget
end subroutine biomatvecset
end module matvecsetmod


module abiosparsemod
contains
subroutine abio_sparse(abio,aout,inbiovar)
  use parameters
  use mpimod
  use matvecsetmod
  use biorthotypemod
  implicit none
  type(biorthotype),target :: inbiovar
  integer :: liwsp,lwsp,itrace,iflag,ixx,getlen
  integer, allocatable :: iwsp(:)
  integer*8, save :: icalled=0
  integer :: biofileptr=6719
  real*8 :: t,anorm, tol
  real*8,save :: tempstepsize=-1d0
  DATATYPE :: abio(numconfig,inbiovar%bionr), aout(numconfig,inbiovar%bionr)
  DATATYPE, allocatable :: wsp(:), smallvector(:,:), smallvectorout(:,:),          aouttr(:,:)
  external biomatvec, parbiomatvec_transpose,realpardotsub


  t=1d0;  anorm=1d0 ;  itrace=0;   iflag=0; tol=biotol;   icalled=icalled+1
  if (icalled.eq.1) then
     tempstepsize=1d0
  endif

  if (inbiovar%thisbiodim.lt.inbiovar%biomaxdim) then
     tempstepsize=tempstepsize*4
  else
     tempstepsize=tempstepsize*1.1
  endif

  if ((myrank.eq.1).and.(notiming==0)) then
     if (icalled.eq.1) then
        open(biofileptr,file=timingdir(1:getlen(timingdir)-1)//"/biortho.dat",status="unknown")
        write(biofileptr,*) 
        write(biofileptr,*);        close(biofileptr)
     endif

!     open(biofileptr,file=timingdir(1:getlen(timingdir)-1)//"/biortho.dat",status="old", position="append")
!     write(biofileptr,*) " BIORTHO. numr",inbiovar%bionr," biodim ",inbiovar%thisbiodim, "step ",min(1d0,tempstepsize)
!     close(biofileptr)
  endif

  if (myrank.eq.1.and.notiming.eq.0) then
     open(biofileptr,file=timingdir(1:getlen(timingdir)-1)//"/biortho.dat",status="old", position="append")
  else
     open(biofileptr,file="/dev/null",status="unknown")
  endif

  call biomatvecset(inbiovar)

  ixx=maxconfigsperproc*inbiovar%bionr
  liwsp=max(12,inbiovar%thisbiodim+3);  lwsp=ixx*(liwsp+1) + 6*(liwsp)**2 + 6+1
  
  allocate(wsp(lwsp),iwsp(liwsp))
  allocate(smallvector(inbiovar%bionr,maxconfigsperproc),smallvectorout(inbiovar%bionr,maxconfigsperproc), &
          aouttr(inbiovar%bionr,numconfig))
  smallvector(:,:)=0; smallvectorout(:,:)=0
  smallvector(:,1:topwalk-botwalk+1)=TRANSPOSE(abio(botwalk:topwalk,:))

#ifdef REALGO
  call DGEXPVxxx2(ixx,inbiovar%thisbiodim,t,smallvector,smallvectorout,tol,anorm,wsp,lwsp,iwsp,liwsp,parbiomatvec_transpose,itrace,iflag,biofileptr,tempstepsize,realpardotsub,maxconfigsperproc*nprocs*inbiovar%bionr)     !!!! numconfig*inbiovar%bionr)
#else
  lwsp=2*lwsp
  call DGEXPVxxx2(2*ixx,inbiovar%thisbiodim,t,smallvector,smallvectorout,tol,anorm,wsp,lwsp,iwsp,liwsp,parbiomatvec_transpose,itrace,iflag,biofileptr,tempstepsize,realpardotsub,2*maxconfigsperproc*nprocs*inbiovar%bionr)    !! 2*numconfig*inbiovar%bionr)
  lwsp=lwsp/2
#endif
  if(iflag.ne.0) then
     print *, "Stopping due to bad iflag in sparsebiortho: ",iflag;     call mpistop()
  endif
  aouttr(:,:)=0d0
  aouttr(:,botwalk:topwalk)=smallvectorout(:,1:topwalk-botwalk+1)
  if (sparseconfigflag.ne.0) then
     call mpiallgather(aouttr,numconfig*inbiovar%bionr,configsperproc*inbiovar%bionr,maxconfigsperproc*inbiovar%bionr)
  endif
  aout(:,:)=TRANSPOSE(aouttr(:,:))
  deallocate(aouttr,smallvector,smallvectorout)

  if (myrank.eq.1.and.notiming.eq.0) then
     open(biofileptr,file=timingdir(1:getlen(timingdir)-1)//"/biortho.dat",status="old", position="append")
     write(biofileptr,*) " Bio: steps, iter, stepsize", iwsp(4),iwsp(1),tempstepsize; 
     if ((iwsp(4).gt.1)) then
        write(biofileptr,*) "Warning - biorthogonalization restarting. Ddim, tol: ", inbiovar%thisbiodim, biotol; 
     endif
     close(biofileptr)
  endif


!!  if (mod(icalled,4000)==0.and.iwsp(4).eq.1) then
  if (mod(icalled,20)==0.and.iwsp(4).eq.1) then
     inbiovar%thisbiodim=max(min(2,inbiovar%biomaxdim),inbiovar%thisbiodim-1)
  endif
  if (iwsp(4).gt.1) then
     inbiovar%thisbiodim=min(inbiovar%thisbiodim+2, inbiovar%biomaxdim)
  endif
  if (iwsp(4).gt.5) then
     inbiovar%thisbiodim=min(inbiovar%thisbiodim*2,inbiovar%biomaxdim)
  endif
  deallocate(wsp,iwsp)

end subroutine abio_sparse
end module abiosparsemod


!!!!!!!!!!!!!!!!!!!!!!!!!

!! these are the encapuslated routines necessary for bi-orthonoramlization as we've derived it
!! the program takes a ton of inputs so it needs not ever use a module other than it's own and mpi
!! the code is also set to be self-resetting if a sour input is passed to it
!! right now if the number of configurations or electrons changes it resets, 
!! this is because the only variables not instanced just once are the spin config lists (alpha, beta, etc) 
!! that are needed for nonsparse routines, and the depend on the number of congigs and electrons  
!! the sparse algorithm needs no lingering variables as everything is passed and pointed to 
!! for just the duration of the call to biortho
!! xoxo -KVL 

!! all the global variables needed to do biothonormalization without any mod save this one and mpi
!! this mod is also NOT in main_modules because we do not want these variables to be touchable 
!! by any other routines than the ones in this file if we can avoid it

module biorthomod
contains
  subroutine bioset(biotypevar,insmo,innumr)
    use biorthotypemod
    use parameters
    implicit none
    Type(biorthotype) :: biotypevar
    DATATYPE,target :: insmo(nspf,nspf)
    integer :: innumr

    biotypevar%smo=>insmo
    biotypevar%biomaxdim=min(maxbiodim,numconfig*innumr-1)

#ifndef REALGO
    if (biocomplex.eq.0) then
       biotypevar%biomaxdim=numconfig*innumr*2-1
    endif
#endif
    biotypevar%bionr=innumr
    if (.not.biotypevar%started) then
       biotypevar%thisbiodim=min(biodim,biotypevar%biomaxdim)
    endif
    biotypevar%started=.true.
  end subroutine bioset


!! on input:  origmo and abio wave function;  oppmo orbitals to biorthogonalize against
!!   output : origmos have been transformed to mobio, and abio transformed in place

  subroutine biortho(origmo,oppmo,mobio,abio,inbiovar)
    use parameters
    use biorthotypemod
    use abiosparsemod
    implicit none

    type(biorthotype),target :: inbiovar
    integer :: i,j
    integer, save :: icalled=0
    DATATYPE :: origmo(spfsize,nspf),oppmo(spfsize,nspf),mobio(spfsize,nspf),abio(numconfig,inbiovar%bionr),dot,data0,data1
    DATATYPE :: atmp(numconfig,inbiovar%bionr),smosave(nspf,nspf)
    
    icalled=icalled+1
    
    do i=1,nspf                 !! Start by finding S**-1 of the original orbitals
       do j=1,nspf
          inbiovar%smo(i,j)=dot(oppmo(:,i),origmo(:,j),spfsize)
       enddo
    enddo

    if (parorbsplit.eq.3) then
       call mympireduce(inbiovar%smo,nspf**2)
    endif

    
!! make the  bi-orthonormal orbitals 
    smosave(:,:)=inbiovar%smo(:,:)
    call invmatsmooth(inbiovar%smo,nspf,nspf,invtol);   data1=1d0;  data0=0d0   

    call MYGEMM('N','N',spfsize,nspf,nspf,data1,origmo,spfsize,inbiovar%smo,nspf,data0,mobio,spfsize)
    
!! check the biorthonormalization to be within hardwired tolerance 1d-10
!! training wheels turned off by KVL for speed on larger jobs
!do i=1,nspf
!  do j=1,nspf
!    data0=dot(oppmo(:,i),mobio(:,j),biospfsize);    val=0d0;   if(i.eq.j) val=1d0
!    if(abs(real(data0))-val      .gt.1d-10) print *, 'bad real ',i,j,data0 
!    if(abs(imag(data0+(0d0,0d0))).gt.1d-10) print *, 'bad imag ',i,j,data0 
!  enddo
!enddo
!! do the a-vector backtransform, sparse or nonsparse
    
    if(sparseconfigflag.eq.0) then
       call abio_nonsparse(inbiovar%smo,abio,atmp,inbiovar%bionr)
    else
!       call neglnmat(inbiovar%smo,nspf,lntol) !! transform s to -ln(s)
       inbiovar%smo(:,:)=smosave(:,:)
       call lnmat(inbiovar%smo(:,:),nspf,lntol) 
       call abio_sparse(abio,atmp,inbiovar)
    endif
    
! hangs sometimes (n2) uncertain why yet
!  if (mod(icalled,1000).eq.0) then
! call checkbio(origmo,mobio,abio,atmp,biospfsize,nspf,nr)
!  endif
    
    abio(:,:)=atmp
    
  end subroutine biortho

!! given origmo and abio wave function, transform abio such that it is now coefficients of oppmo.
!!   i.e. operate with S-inverse on a-vector from s computed from input orbs

  subroutine biotransform(origmo,oppmo,abio,inbiovar)
    use parameters
    use biorthotypemod
    use abiosparsemod
    implicit none

    type(biorthotype),target :: inbiovar
    integer :: i,j
    integer, save :: icalled=0
    DATATYPE :: origmo(spfsize,nspf),oppmo(spfsize,nspf),abio(numconfig,inbiovar%bionr),dot
    DATATYPE :: atmp(numconfig,inbiovar%bionr)
    
    icalled=icalled+1
    
    do i=1,nspf         
       do j=1,nspf
          inbiovar%smo(i,j)=dot(origmo(:,i),oppmo(:,j),spfsize)
       enddo
    enddo

    if (parorbsplit.eq.3) then
       call mympireduce(inbiovar%smo,nspf**2)
    endif


    if(sparseconfigflag.eq.0) then
       call abio_nonsparse(inbiovar%smo,abio,atmp,inbiovar%bionr)
    else
       call neglnmat(inbiovar%smo,nspf,lntol)    !! transform s to -ln(s)
       call abio_sparse(abio,atmp,inbiovar)
    endif
    
    abio(:,:)=atmp
    
  end subroutine biotransform


!! operate with S

  subroutine biooverlap(origmo,oppmo,abio,inbiovar)
    use parameters
    use biorthotypemod
    use abiosparsemod
    implicit none

    type(biorthotype),target :: inbiovar
    integer :: i,j
    integer, save :: icalled=0
    DATATYPE :: origmo(spfsize,nspf),oppmo(spfsize,nspf),abio(numconfig,inbiovar%bionr),dot
    DATATYPE :: atmp(numconfig,inbiovar%bionr)
    
    icalled=icalled+1
    
    do i=1,nspf         
       do j=1,nspf
          inbiovar%smo(i,j)=dot(origmo(:,i),oppmo(:,j),spfsize)
       enddo
    enddo

    if (parorbsplit.eq.3) then
       call mympireduce(inbiovar%smo,nspf**2)
    endif


    if(sparseconfigflag.eq.0) then
       call abio_nonsparse(inbiovar%smo,abio,atmp,inbiovar%bionr)
    else
       call lnmat(inbiovar%smo,nspf,lntol)    !! transform s to ln(s)
       call abio_sparse(abio,atmp,inbiovar)
    endif
    
    abio(:,:)=atmp
    
  end subroutine biooverlap


  
end module biorthomod



!! this is the non-sparse specific routine that does the back-transofrmation of the a-vector
subroutine abio_nonsparse(smo,abio,aout,nr)
  use parameters
  use configmod
  use mpimod
  use aarrmod

  implicit none
  integer :: nr,i,j,iflag,clow,chigh,jproc,cnum,nnn(2),iind,mmm(2)
  integer :: ipiv(numconfig),bioconfiglist(numelec,numconfig)
  DATATYPE :: abio(numconfig,nr),aout(numconfig,nr),matdet,smobig(nspf*2,nspf*2),Stmpbig(numelec,numelec)
  DATATYPE  :: Sconfig(numconfig,numconfig),smo(nspf,nspf)
  
!! for the nonsparse routine this builds the full nonsparse configuration overlap matrix
!! this relies on the unique properties of the Doolittle algorithm of LU factorization
!! to take the overlap matrices of the ith and jth configs alpha and beta orbitals and get their determinant that way 
!! this is very much so a brute force way to approach this problem

  aout(:,:)=abio(:,:);  smobig(:,:)=0d0

  do i=1,nspf*2
     mmm(:)=aarr(i,nspf)
     do j=1,nspf*2
        nnn(:)=aarr(j,nspf)
        if (mmm(2).eq.nnn(2)) then
           smobig(i,j)=smo(mmm(1),nnn(1))
        endif
     enddo
  enddo

  do i=1,numconfig
     do j=1,numelec
        bioconfiglist(j,i)=iind( configlist(j*2-1:j*2,i) )
     enddo
  enddo
  
  do j=1,numconfig
     do i=1,numconfig
        call get_petite_mat(nspf*2,numelec,Smobig,Stmpbig,bioconfiglist(:,i),bioconfiglist(:,j))
        sconfig(i,j) = matdet(numelec,Stmpbig)
     enddo
  enddo
  
!! this is where the linear equation solver is called to solve S*abio'=abio to get our  abio'

  clow = (myrank-1)*nr/nprocs+1;  chigh = myrank*nr/nprocs

  call MYGESV(numconfig,chigh-clow+1,Sconfig,numconfig,ipiv,aout(:,clow),numconfig,iflag)

  if(iflag.ne.0) then
     OFLWR "Stopping due to bad iflag in nonsparsebiortho: ",iflag; CFLST
  endif

!! if we aren't parallel this means several different calls to biortho maybe going on at once across all procs
!! so we REALLY do not want to broadcast the a-vector in that case
!! this is the case in electronflux where the parallelization is over the bras of the flux integral for a specific ket
!! this is also the case in projeflux where the parallelization is over time in forming the projected single particle function

  do jproc=1,nprocs
     clow = (jproc-1)*nr/nprocs+1;      chigh = jproc*nr/nprocs;      cnum = (chigh-clow+1)*numconfig
     call mympibcast(aout(:,clow:),jproc,cnum)
  enddo

end subroutine abio_nonsparse


!! NOTE BOUNDS !!

subroutine parbiomatvec_transpose(inavectortr,outavectortr)
  use parameters
  use mpimod
  use matvecsetmod
  use biomatvecmod
  implicit none

!!  DATATYPE :: inavectortr(biopointer%bionr,botwalk:topwalk), outavectortr(biopointer%bionr,botwalk:topwalk)
  DATATYPE :: inavectortr(biopointer%bionr,botwalk:botwalk+maxconfigsperproc-1), outavectortr(biopointer%bionr,botwalk:botwalk+maxconfigsperproc-1)

  DATATYPE :: intemptr(biopointer%bionr,numconfig), ttvector(numconfig,biopointer%bionr), ttvector2(botwalk:topwalk,biopointer%bionr)

  outavectortr(:,:)=0d0

  if (sparseconfigflag.eq.0) then
     OFLWR "error, must use sparse for parbiomatvec_transpose"; CFLST
  endif

  intemptr(:,:)=0d0
  intemptr(:,botwalk:topwalk)=inavectortr(:,botwalk:topwalk)

  if (sparseconfigflag.ne.0) then
     call mpiallgather(intemptr,numconfig*biopointer%bionr,configsperproc*biopointer%bionr,maxconfigsperproc*biopointer%bionr)
  endif

  ttvector(:,:)=TRANSPOSE(intemptr(:,:))

  call biomatvec_nompi(ttvector,ttvector2)

  outavectortr(:,botwalk:topwalk)=TRANSPOSE(ttvector2(botwalk:topwalk,:))
  
end subroutine parbiomatvec_transpose




!! this is the sparse specific routine that does the back-transofrmation of the a-vector

!! the external subroutine
!! does y=-ln(S)*x for the formation of the Krylov vectors
!! input:
!! x - previous krylov vector
!! output:
!! y - the  krylov vector y = -ln(S)x
!! data needed from biorthomod:
!! Smo - the overlap matrix containing -ln(s)
!! the config and walk mod objects are pointers to the configlist and singlewalk arrays 
!! that were passed to biortho when it was called this time

subroutine biomatvec(x,y)
  use matvecsetmod
  use biomatvecmod
  use walkmod
  use parameters
  implicit none
  integer :: i,j
  DATATYPE :: x(numconfig,biopointer%bionr),y(numconfig,biopointer%bionr),ytr(biopointer%bionr,numconfig)

  y(:,:)=0d0

  do i=botwalk,topwalk
!!$     tmpval=0d0
!!$removed 09-2016 simple walks
!!$     do j=1,numelec
!!$        tmpval = tmpval + biopointer%smo(configlist((j-1)*2+1,i),configlist((j-1)*2+1,i))
!!$     enddo
     y(i,:)=0    !!$tmpval*x(i,:)
     do j=1,numsinglewalks(i) !! summing over nonconjugated second index in s(:), good -- s(i,j) is ln(<i|j>) correct?
        y(i,:) = y(i,:) + biopointer%smo(singlewalkopspf(1,j,i),singlewalkopspf(2,j,i)) * singlewalkdirphase(j,i) * x(singlewalk(j,i),:) 
     enddo
  enddo

  if (sparseconfigflag.ne.0) then
     ytr(:,:)=TRANSPOSE(y)
     call mpiallgather(ytr,numconfig*biopointer%bionr,configsperproc*biopointer%bionr,maxconfigsperproc*biopointer%bionr)
     y=TRANSPOSE(ytr)
  endif

!!!!  call mympireduce(y,numconfig*biopointer%bionr)



end subroutine biomatvec





subroutine biomatvec_nompi(x,y)
  use matvecsetmod
  use biomatvecmod
  use walkmod
  use parameters
  implicit none
  integer :: i,j
  DATATYPE :: x(numconfig,biopointer%bionr),y(botwalk:topwalk,biopointer%bionr)

  y(:,:)=0d0

  do i=botwalk,topwalk
!!$removed 09-2016 simple walks     tmpval=0d0
!!$     do j=1,numelec
!!$        tmpval = tmpval + biopointer%smo(configlist((j-1)*2+1,i),configlist((j-1)*2+1,i))
!!$     enddo
     y(i,:)=0   !!$tmpval*x(i,:)
     do j=1,numsinglewalks(i) !! summing over nonconjugated second index in s(:), good -- s(i,j) is ln(<i|j>) correct?
        y(i,:) = y(i,:) + biopointer%smo(singlewalkopspf(1,j,i),singlewalkopspf(2,j,i)) * singlewalkdirphase(j,i) * x(singlewalk(j,i),:) 
     enddo
  enddo


end subroutine biomatvec_nompi




!! NOW for nat check, not biortho.  ONLY FOR DEBUG, HANGS.

subroutine checkbio(origmo,mobio,abio,atmp)
  use parameters
  use configmod
  implicit none
  DATATYPE :: mobio(spfsize,nspf), origmo(spfsize,nspf), atmp(numconfig,numr), abio(numconfig,numr), checkoverlap, check1,check2
  integer, save :: icalled=-1
  real*8 :: err
  icalled=icalled+1

  OFLWR "SOMETHING IS UP.  BIOTRANSFORM DEFINITELY WORKS BUT CHECKBIO DOES NOT AGREE. AUTOCORRELATE_ONE LIKELY BROKEN."; CFLST

  call openfile(); write(mpifileptr,*) "   .... check biortho perm " ; call closefile()

  call permoverlaps(numr,numelec,spfsize,origmo,mobio,abio,atmp,checkoverlap,0,0.d-8,0.d-8,nspf,nspf,numconfig,numconfig,configlist,numelec*2,configlist,numelec*2,0,parorbsplit)
  call permoverlaps(numr,numelec,spfsize,origmo,origmo,abio,abio,check1,0,0.d-8,0.d-8,nspf,nspf,numconfig,numconfig,configlist,numelec*2,configlist,numelec*2,0,parorbsplit)
  call permoverlaps(numr,numelec,spfsize,mobio,mobio,atmp,atmp,check2,0,0.d-8,0.d-8,nspf,nspf,numconfig,numconfig,configlist,numelec*2,configlist,numelec*2,0,parorbsplit)

  err=abs((checkoverlap-sqrt(check1*check2))/checkoverlap)
  if (err.gt.1.d-7) then
     OFLWR "Checkbio: FAIL ", checkoverlap,check1,check2,err,numconfig*numr; CFL
  else
     OFLWR "Checkbio: err",err; CFL
  endif

end subroutine checkbio


!! OLD OLD OLD
!! this is hanging sometimes.... disabling above
!!$ 
!!$ subroutine checkbio(origmo,mobio,abio,atmp,insize,norb,nr)
!!$   use biorthomod
!!$   use mpimod
!!$   implicit none
!!$   DATATYPE :: mobio(insize,norb), origmo(insize,norb), atmp(numconfig,nr), abio(numconfig,nr), checkoverlap, check1,check2
!!$   integer, save :: icalled=-1
!!$   integer :: nr,insize,norb
!!$   real*8 :: err
!!$   icalled=icalled+1
!!$ 
!!$ !!  call openfile(); write(mpifileptr,*) "   .... check biortho perm " ; call closefile()
!!$ 
!!$   call permoverlaps(nr,numelec,insize,origmo,mobio,abio,atmp,checkoverlap,0,0.d-8,0.d-8,norb,norb,numconfig,numconfig,configlist,numelec*2,configlist,numelec*2,0)
!!$   call permoverlaps(nr,numelec,insize,origmo,origmo,abio,abio,check1,0,0.d-8,0.d-8,norb,norb,numconfig,numconfig,configlist,numelec*2,configlist,numelec*2,0)
!!$   call permoverlaps(nr,numelec,insize,mobio,mobio,atmp,atmp,check2,0,0.d-8,0.d-8,norb,norb,numconfig,numconfig,configlist,numelec*2,configlist,numelec*2,0)
!!$ 
!!$   err=abs((checkoverlap-sqrt(check1*check2))/checkoverlap)
!!$   if (err.gt.1.d-6) then
!!$      OFLWR "PERM FAIL ", checkoverlap,check1,check2,err,numconfig*nr; CFL
!!$   else
!!$      OFLWR "Checkbio: err",err; CFL
!!$   endif
!!$ end subroutine checkbio






