
#include "Definitions.INC"

!! ALL MODULES except for get_fockmatrix

module cathopsubmod
contains

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

end subroutine getcathops

end module cathopsubmod


module focksubmod
contains

subroutine get_fockmatrix0(fockmatrix,avectors,cptr,sptr)
  use parameters
  use aarrmod
  use configmod
  use sparseptrmod
  use configptrmod
  use mpisubmod
  use sparsemultmod
  use asssubmod
  implicit none
  DATATYPE,intent(out) :: fockmatrix(www%nspf,www%nspf)
  DATATYPE,intent(in) :: avectors(numr,www%firstconfig:www%lastconfig,mcscfnum)
  Type(CONFIGPTR),intent(in) :: cptr
  Type(SPARSEPTR),intent(in) :: sptr
  Type(SPARSEPTR) :: catsptr
  DATATYPE,allocatable :: catvects(:,:,:), hugeavector(:,:), catmult(:,:),neutmult(:,:),bigcatvector(:,:)
  DATATYPE :: tempfockmatrix(www%nspf),nullvector1(numr),nullvector2(numr)          !! AUTOMATIC
  DATATYPE :: neutenergy
  integer :: imc, ispf, jspf, ihop

!! this isn't the best, lots of zeroes.  should switch the way it's done.

  if (.not.use_fockmatrix) then
     OFLWR "What? use_fockmatrix get_fockmatrix"; CFLST
  endif

  tempfockmatrix=0
  fockmatrix=0

  allocate(catvects(numr,catww%botconfig:catww%topconfig,nspf))
  if (catww%topconfig.ge.catww%botconfig) then
     catvects=0
  endif
  allocate(bigcatvector(numr,catww%firstconfig:catww%lastconfig),&
       catmult(numr,catww%firstconfig:catww%lastconfig),&
       neutmult(numr,www%firstconfig:www%lastconfig))
  if (catww%totadim.gt.0) then
     bigcatvector=0; catmult=0
  endif
  if (www%totadim.gt.0) then
     neutmult=0
  endif

  allocate(hugeavector(numr,www%numconfig))

  if (sparseopt.ne.0) then
     call sparseptralloc(catsptr,catww)
     call assemble_sparsemats(catww,cptr,catsptr,1,1,0,0)
  endif

  do imc=1,mcscfnum

     if (www%totadim.gt.0) then
        call sparseconfigmult(www,avectors(:,:,imc),neutmult,cptr,sptr,&
             1,1,0,0,0d0,-1)
     else
        call sparseconfigmult(www,nullvector1,nullvector2,cptr,sptr,&
             1,1,0,0,0d0,-1)
     endif
     neutenergy=0d0
     if (www%totadim.gt.0) then
        neutenergy = dot(avectors(:,:,imc),neutmult(:,:),www%totadim)
     endif
     if (www%parconsplit.ne.0) then
        call mympireduceone(neutenergy)
     endif

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

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(jspf)
!$OMP DO SCHEDULE(DYNAMIC)
           do jspf=1,www%nspf
              tempfockmatrix(jspf) = &
                   dot(catvects(:,:,jspf),catvects(:,:,ispf),&
                   numr*(catww%topconfig-catww%botconfig+1)) * neutenergy - &
                   dot(catvects(:,:,jspf),catmult(:,catww%botconfig:catww%topconfig),&
                   numr*(catww%topconfig-catww%botconfig+1))
           enddo
!$OMP END DO
!$OMP END PARALLEL

        endif

        fockmatrix(:,ispf)=fockmatrix(:,ispf)+tempfockmatrix(:)

     enddo  !! ispf

  end do  !! imc

  deallocate(catvects,hugeavector,bigcatvector,catmult,neutmult)
  if (sparseopt.ne.0) then
     call sparseptrdealloc(catsptr)
  endif

  call mympireduce(fockmatrix,www%nspf**2)

  fockmatrix(:,:)=fockmatrix(:,:)/mcscfnum

!! TRANSPOSE !!
  fockmatrix(:,:)=TRANSPOSE(fockmatrix(:,:))

end subroutine get_fockmatrix0

end module focksubmod


subroutine get_fockmatrix()
  use xxxmod
  use focksubmod
  implicit none
  call get_fockmatrix0(yyy%fockmatrix(:,:,0),yyy%cmfavec(:,:,0),yyy%cptr(0),yyysptr(0))
end subroutine get_fockmatrix



module fockrepsubmod
  use biorthotypemod
  implicit none
  type(biorthotype),target, private :: fockrepbiovar

contains

subroutine replace_withfock(printflag)
  use parameters
  use class_parameters
  use configmod
  use biorthomod
  use xxxmod
  use matsubmod
  use eigenmod
  use spfsubmod
  implicit none
  integer,intent(in) :: printflag
  DATATYPE,allocatable :: fockeigvects(:,:), fockop(:,:), myfock(:,:),&
       fockvects(:,:), outspfs(:,:)
  CNORMTYPE,allocatable :: fockeigvals(:),fockvals(:)
  DATATYPE,target :: smo(nspf,nspf)    !! AUTOMATIC
  DATATYPE :: csum
  real*8 :: errorval
  integer :: iclass, ispf,jspf, imc

  if (.not.use_fockmatrix) then
     OFLWR "What? use fockmatrix replace_withfock"; CFLST
  endif

  smo=0d0;  

  allocate(fockvals(nspf),fockvects(nspf,nspf))
  fockvals=0; fockvects=0

  do iclass=1,numclasses

     allocate(fockeigvals(nperclass(iclass)), fockeigvects(nperclass(iclass),nperclass(iclass)), &
          fockop(nperclass(iclass),nperclass(iclass)),myfock(nperclass(iclass),nperclass(iclass)))

     fockeigvals=0; fockeigvects=0; fockop=0; myfock=0

     do jspf=1,nperclass(iclass)
        do ispf=1,nperclass(iclass)
           myfock(ispf,jspf)=yyy%fockmatrix(classorb(ispf,iclass),classorb(jspf,iclass),0)
        enddo
     enddo

!! TAKE HERMITIAN PART IF HERM NORM (does nothing if c-norm, CONJUGATE is nothing)

     fockop(:,:) = (+0.5d0) * ( myfock(:,:) + TRANSPOSE(CONJUGATE(myfock(:,:))) )

! FOR HERMITIAN if HERM NORM
     call EIGEN(fockop,nperclass(iclass),nperclass(iclass),fockeigvects,fockeigvals)

     do ispf=1,nperclass(iclass)
        fockvals(classorb(ispf,iclass))=fockeigvals(ispf) * (+1)
        do jspf=1,nperclass(iclass)
           fockvects(classorb(ispf,iclass),classorb(jspf,iclass)) = fockeigvects(ispf,jspf)
        enddo
     enddo

     deallocate(fockeigvals, fockeigvects, fockop,myfock)

  enddo  !! do iclass

  do ispf=1,nspf
     csum=fockvects(ispf,ispf)
#ifdef CNORMFLAG
     if (real(csum,8).lt.0d0) then
        fockvects(:,ispf) = fockvects(:,ispf) * (-1)
     endif
#else
     if (csum.ne.0d0) then
        fockvects(:,ispf) = fockvects(:,ispf) * abs(csum) / csum
     endif
#endif
  enddo

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
     OFLWR "ERROR IN REPLACE_WITHFOCK ", errorval; CFLST
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

end module fockrepsubmod


