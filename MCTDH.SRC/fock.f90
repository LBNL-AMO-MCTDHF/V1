
#include "Definitions.INC"


!! here we are doing walks from our BO target state to our wavefunction \Psi(t)
subroutine getcathops()
  use fileptrmod
  use aarrmod
  use configmod
  use mpimod
  use mpisubmod
  use configsubmod
  implicit none
  integer :: iconfig,jconfig,idof,iindex,flag,dirphase, idfconfig,&
       tempconfig(catww%num2part),temporb(2), ii, ihop, totcathops

!! get the number of single hops from the target N-1 e- state to our regular N e- state

  do ii=1,2

     if (ii.eq.1) then
        OFLWR "Counting cation single hops for fock"; CFL
     else
        OFLWR "Getting cation single hops for fock"; CFL

        allocate(cathopfrom(numcathops),cathopto(numcathops),catphase(numcathops),&
             catspf(numcathops))
!! catphase=0 no hop
        cathopfrom(:)=(-1);  cathopto(:)=(-1);  catphase(:)=0; catspf(:)=(-1)

     endif

     ihop=0
     do idfconfig=1,www%numdfconfigs

        iconfig=www%ddd%dfincludedconfigs(idfconfig)

!!     do iindex=1,2*www%nspf

        do iindex=1,2*www%nspf,2   !! spin down for now right?

           temporb=aarr(iindex)

           if (www%holeflag.ne.0) then
              tempconfig(1:2)=temporb(:)
              tempconfig(3:catww%num2part)=www%configlist(:,iconfig)
              flag=0
              do idof=2,catww%numpart 
                 if(iind(tempconfig(idof*2-1:idof*2)).eq.iindex) then
                    flag=1
                    exit
                 endif
              enddo
              dirphase=reorder(tempconfig,catww%numpart)
           else   !! holeflag.eq.0:
              flag=1
              do idof=1,www%numpart
                 dirphase=(-1)**idof
                 if(iind(www%configlist(idof*2-1:idof*2,iconfig)).eq.iindex) then
                    flag=0
                    tempconfig(1:(idof-1)*2) = www%configlist(1:(idof-1)*2,iconfig)
                    tempconfig(idof*2-1:catww%num2part) = www%configlist(idof*2+1:www%num2part,iconfig)
                    exit
                 endif
              enddo
           endif   !! holeflag
           if(flag.eq.0) then
              if(allowedconfig0(catww,tempconfig,catww%dfrestrictflag)) then
                 jconfig=getconfiguration(tempconfig,catww)
                 if (jconfig.ge.catww%botconfig.and.jconfig.le.catww%topconfig) then
                    ihop=ihop+1
                    if (ii.eq.2) then
                       cathopfrom(ihop)=iconfig
                       cathopto(ihop)=jconfig
                       catspf(ihop)=temporb(1)
                       catphase(ihop)=dirphase
                    endif
                 endif
              endif
           endif
        enddo   !! iindex (orbital)
     enddo  !! idfconfig

     if (ii.eq.1) then
        OFLWR "number of cation fock hops this processor",ihop; CFL
        numcathops=ihop
        totcathops=ihop
        call mympiireduceone(totcathops)
        OFLWR "         number of cation fock hops total",totcathops; CFL
     else
        if (numcathops.ne.ihop) then
           OFLWR "ERROR NUMCATHOPS", numcathops, ihop; CFLST
        endif
     endif
  enddo

  call mpibarrier()
  OFLWR "DONE getting fock cation single hops"; CFL


!!
!  print *, "CATHOPS"
!  do ihop=1,numcathops
!     print *, ihop, cathopfrom(ihop), cathopto(ihop), catspf(ihop), catphase(ihop)
!  enddo
!  print *, "TTSTOP"
!  call mpistop()

end subroutine getcathops


subroutine get_fockmatrix0(fockmatrix,fockden,avectors,cptr,neutsptr)
  use parameters
  use aarrmod
  use configmod
  use sparseptrmod
  use configptrmod
  use mpisubmod
  use sparsemultmod
  use asssubmod
  implicit none
  DATATYPE,intent(out) :: fockmatrix(www%nspf,www%nspf), fockden(www%nspf,www%nspf)
  DATATYPE,intent(in) :: avectors(numr,www%firstconfig:www%lastconfig,mcscfnum)
  Type(CONFIGPTR),intent(in) :: cptr
  Type(SPARSEPTR),intent(in) :: neutsptr
  Type(SPARSEPTR) :: catsptr
  DATATYPE,allocatable :: catvects(:,:,:), hugeavector(:,:), catmult(:,:),bigcatvector(:,:),&
       neutvects(:,:,:), neutmult(:,:), bigneutvector(:,:)
  DATATYPE :: tempfockmatrix(www%nspf),nullvector1(numr),nullvector2(numr),&          !! AUTOMATIC
       tempfockden(www%nspf,www%nspf)
  DATATYPE :: norm,csum,myexpect
  integer :: imc, ispf, jspf, ihop

#define FIRSxxxTFOCK

#ifdef FIRSTFOCK
  integer :: idof, iconfig, ispin
#endif

!! this isn't the best, lots of zeroes.  should switch the way it's done.

  if (.not.use_fockmatrix) then
     OFLWR "What? use_fockmatrix get_fockmatrix"; CFLST
  endif

  tempfockmatrix=0; tempfockden=0
  fockmatrix=0; fockden=0

  allocate( neutmult(numr,www%firstconfig:www%lastconfig))
  if (www%totadim.gt.0) then
     neutmult=0
  endif

  norm=0; csum=0
  do imc=1,mcscfnum
     if (www%totadim.gt.0) then
        norm=norm+dot(avectors(:,:,imc),avectors(:,:,imc),numr*www%totadim)
        call sparseconfigmult(www,avectors(:,:,imc),neutmult(:,:),cptr,neutsptr,&
             1,1,0,0,0d0,-1)
        csum=csum + dot(avectors(:,:,imc),neutmult(:,:),numr*www%totadim)
     else
        call sparseconfigmult(www,nullvector1,nullvector2,cptr,neutsptr,&
             1,1,0,0,0d0,-1)
     endif
  enddo

  if (www%parconsplit.ne.0) then
     call mympireduceone(norm)
     call mympireduceone(csum)
  endif
  myexpect=csum/norm

!!$  OFLWR " expectation val for fock", myexpect; CFL

  myexpect=real(myexpect,8)   !! making things hermitian

  allocate(catvects(numr,catww%botconfig:catww%topconfig,nspf))
  if (catww%topconfig.ge.catww%botconfig) then
     catvects=0
  endif

  allocate(neutvects(numr,www%botconfig:www%topconfig,nspf))
  if (www%topconfig.ge.www%botconfig) then
     neutvects=0
  endif

  allocate(bigcatvector(numr,catww%firstconfig:catww%lastconfig),&
       bigneutvector(numr,www%firstconfig:www%lastconfig),&
       catmult(numr,catww%firstconfig:catww%lastconfig))
  if (catww%totadim.gt.0) then
     bigcatvector=0; catmult=0
  endif
  if (www%totadim.gt.0) then
     bigneutvector=0
  endif

  allocate(hugeavector(numr,www%numconfig))

  if (sparseopt.ne.0) then
     call sparseptralloc(catsptr,catww)
     call assemble_sparsemats(catww,cptr,catsptr,1,1,0,0)
  endif

  do imc=1,mcscfnum

     hugeavector=0d0

     if (www%totadim.gt.0) then
        hugeavector(:,www%firstconfig:www%lastconfig) = avectors(:,www%firstconfig:www%lastconfig,imc)
     endif

     if (www%parconsplit.ne.0) then
        call mpiallgather(hugeavector,numr*www%numconfig,numr*www%configsperproc(:),&
             numr*www%maxconfigsperproc)
     endif

     if (catww%topconfig.ge.catww%botconfig) then
        catvects=0d0
     endif
     do ihop=1,numcathops
        catvects(:,cathopto(ihop),catspf(ihop))=hugeavector(:,cathopfrom(ihop))*catphase(ihop)
     enddo

#ifdef FIRSTFOCK

     if (www%topconfig.ge.www%botconfig) then
        neutvects=0d0
     endif
     do iconfig=www%botconfig,www%topconfig
        do idof=1,www%numpart
           ispf=www%configlist(idof*2-1,iconfig)
           ispin=www%configlist(idof*2,iconfig)
!!           if (ispin.eq.1)  then  !! OK?  remove ispin=1 spin down not ispin=2 spin up
              neutvects(:,iconfig,ispf) = hugeavector(:,iconfig) * (-1)**idof
!!           endif
        enddo
     enddo

#endif

     do ispf=1,www%nspf

        if (catww%totadim.gt.0) then
           bigcatvector=0
           if (catww%topconfig.ge.catww%botconfig) then
              bigcatvector(:,catww%botconfig:catww%topconfig)=catvects(:,:,ispf)
           endif
        endif
        if (catww%parconsplit.eq.0) then
           call mpiallgather(bigcatvector,numr*catww%numconfig,numr*catww%configsperproc(:),&
                numr*catww%maxconfigsperproc)
        endif

        if (catww%totadim.gt.0) then
           call sparseconfigmult(catww,bigcatvector,catmult,cptr,catsptr,&
                1,1,0,0,0d0,-1)
        else
           call sparseconfigmult(catww,nullvector1,nullvector2,cptr,catsptr,&
                1,1,0,0,0d0,-1)
        endif

        tempfockmatrix(:)=0d0

        if (catww%topconfig.ge.catww%botconfig) then

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ispf,jspf)
!$OMP DO SCHEDULE(DYNAMIC)
           do jspf=1,www%nspf
              tempfockmatrix(jspf) = &
                   dot(catvects(:,:,jspf),catmult(:,catww%botconfig:catww%topconfig),&
                   numr*(catww%topconfig-catww%botconfig+1))
           enddo
!$OMP END DO
!$OMP END PARALLEL

        endif

        fockmatrix(:,ispf)=fockmatrix(:,ispf)-tempfockmatrix(:)


#ifdef FIRSTFOCK

        if (www%totadim.gt.0) then
           bigneutvector=0
           if (www%topconfig.ge.www%botconfig) then
              bigneutvector(:,www%botconfig:www%topconfig)=neutvects(:,:,ispf)
           endif
        endif
        if (www%parconsplit.eq.0) then
           call mpiallgather(bigneutvector,numr*www%numconfig,numr*www%configsperproc(:),&
                numr*www%maxconfigsperproc)
        endif

        if (www%totadim.gt.0) then
           call sparseconfigmult(www,bigneutvector,neutmult,cptr,neutsptr,&
                1,1,0,0,0d0,-1)
        else
           call sparseconfigmult(www,nullvector1,nullvector2,cptr,neutsptr,&
                1,1,0,0,0d0,-1)
        endif


        tempfockmatrix(:)=0d0

        if (www%topconfig.ge.www%botconfig) then

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ispf,jspf)
!$OMP DO SCHEDULE(DYNAMIC)
           do jspf=1,www%nspf
              tempfockmatrix(jspf) = &
                   dot(neutvects(:,:,jspf),neutmult(:,www%botconfig:www%topconfig),&
                   numr*(www%topconfig-www%botconfig+1))
           enddo
!$OMP END DO
!$OMP END PARALLEL

        endif

        fockmatrix(:,ispf)=fockmatrix(:,ispf)+tempfockmatrix(:)

#endif

     enddo  !! ispf

     tempfockden=0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ispf,jspf)
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2)
     do ispf=1,www%nspf
        do jspf=1,www%nspf
           tempfockden(jspf,ispf)=dot(catvects(:,:,jspf),catvects(:,:,ispf),&
                numr*(catww%topconfig-catww%botconfig+1))
        enddo
     enddo
!$OMP END DO
!$OMP END PARALLEL

     fockden(:,:)=fockden(:,:)+tempfockden(:,:)

  end do  !! imc

  deallocate(catvects,neutvects,hugeavector,bigcatvector,bigneutvector,catmult,neutmult)
  if (sparseopt.ne.0) then
     call sparseptrdealloc(catsptr)
  endif

  call mympireduce(fockmatrix,www%nspf**2)
  call mympireduce(fockden,www%nspf**2)

#ifndef FIRSTFOCK
!! SECOND WAY  
  fockmatrix(:,:) = fockmatrix(:,:) + myexpect * fockden(:,:)
#endif

end subroutine get_fockmatrix0



subroutine get_fockmatrix()
  use xxxmod
  implicit none
  call get_fockmatrix0(yyy%fockmatrix(:,:,0),yyy%fockden(:,:,0),yyy%cmfavec(:,:,0),yyy%cptr(0),yyysptr(0))
end subroutine get_fockmatrix


module fockrepbiomod
  use biorthotypemod
  implicit none
  type(biorthotype),target :: fockrepbiovar
end module


subroutine replace_withfock(printflag)
  use parameters
  use class_parameters
  use configmod
  use biorthomod
  use fockrepbiomod
  use xxxmod
  use matsubmod
  implicit none
  integer,intent(in) :: printflag
  DATATYPE,allocatable :: fockeigvects(:,:), fockop(:,:), invsqrtden(:,:), myfock(:,:),&
       fockvects(:,:), outspfs(:,:)
  CNORMTYPE,allocatable :: fockeigvals(:),fockvals(:)
  DATATYPE,target :: smo(nspf,nspf)    !! AUTOMATIC
  real*8 :: errorval
  integer :: iclass, ispf,jspf, imc

  if (.not.use_fockmatrix) then
     OFLWR "What? use fockmatrix replace_withfock"; CFLST
  endif

  smo=0d0;  

!  OFLWR "FOCKMATRIX:"
!  do ispf=1,nspf
!     WRFL yyy%fockmatrix(:,ispf,0)
!  enddo
!  WRFL
!  WRFL "TEMPSSSTOP"; CFLST


  allocate(fockvals(nspf),fockvects(nspf,nspf))
  fockvals=0; fockvects=0

  do iclass=1,numclasses

     allocate(fockeigvals(nperclass(iclass)), fockeigvects(nperclass(iclass),nperclass(iclass)), &
          fockop(nperclass(iclass),nperclass(iclass)),invsqrtden(nperclass(iclass),nperclass(iclass)),&
          myfock(nperclass(iclass),nperclass(iclass)))

     fockeigvals=0; fockeigvects=0; fockop=0; invsqrtden=0; myfock=0

     do jspf=1,nperclass(iclass)
        do ispf=1,nperclass(iclass)
           invsqrtden(ispf,jspf)=yyy%fockden(classorb(ispf,iclass),classorb(jspf,iclass),0)
           myfock(ispf,jspf)=yyy%fockmatrix(classorb(ispf,iclass),classorb(jspf,iclass),0)
        enddo
        invsqrtden(jspf,jspf)=invsqrtden(jspf,jspf)+fockreg
     enddo

!! TAKE HERMITIAN PART IF HERM NORM (does nothing if c-norm, CONJUGATE is nothing)

     myfock(:,:) = 0.5d0 * ( myfock(:,:) + TRANSPOSE(CONJUGATE(myfock(:,:))) )

     call allpurposemat(invsqrtden,nperclass(iclass),1)

!! DO BLAS
     fockop(:,:)=MATMUL(MATMUL(invsqrtden(:,:),myfock(:,:)),invsqrtden(:,:)) !! * (-1)

!!     fockop(:,:)=myfock(:,:) !! * (-1)

!!TRY     fockop(:,:)=TRANSPOSE(myfock(:,:)) * (-1)

! FOR HERMITIAN if HERM NORM
     call EIGEN(fockop,nperclass(iclass),nperclass(iclass),fockeigvects,fockeigvals)

!!     fockeigvals(:) =  fockeigvals(:) * (-1)

     do ispf=1,nperclass(iclass)
        fockvals(classorb(ispf,iclass))=fockeigvals(ispf)
        do jspf=1,nperclass(iclass)
           fockvects(classorb(ispf,iclass),classorb(jspf,iclass)) = fockeigvects(ispf,jspf)
        enddo
     enddo

     deallocate(fockeigvals, fockeigvects, fockop,invsqrtden,myfock)

  enddo  !! do iclass


  allocate(outspfs(spfsize,nspf))
  outspfs=0d0

  do jspf=1,nspf  ! which natorb
     do ispf=1,nspf  ! which original
        outspfs(:,jspf)=outspfs(:,jspf)+ &
             yyy%cmfspfs((ispf-1)*spfsize+1:ispf*spfsize,0)*fockvects(ispf,jspf)
     enddo
  enddo

  call spf_orthogit(outspfs, errorval)
  if (errorval.gt.1d-7) then
     OFLWR "WTF!  ERROR IN REPLACE_WITHFOCK ", errorval; CFLST
  endif

  if (printflag==1) then
     OFLWR "REPLACING SPFS - GENERALIZED FOCK EIGS"
     do ispf=1,nspf
        write(mpifileptr,'(2E25.10)') fockvals(ispf)
     enddo
     WRFL; CFL
  endif

  call bioset(fockrepbiovar,smo,numr,bioww)

  do imc=1,mcscfnum
     call biotransform(yyy%cmfspfs(:,0),outspfs, yyy%cmfavec(:,imc,0),fockrepbiovar)
  enddo

  yyy%cmfspfs(:,0)=RESHAPE(outspfs,(/totspfdim/))

  deallocate(outspfs, fockvals,fockvects)

end subroutine replace_withfock
