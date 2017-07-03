
!! DENSITY MATRIX SUBROUTINES !!  
  
#include "Definitions.INC"
  
!! ALL MODULES

module denutilmod
contains

!! ASSUMES denmat,denvects is ready. 

subroutine output_denmat( incmf, time)
  use parameters
  use mpimod
  use xxxmod
  implicit none
  integer :: i,j, incmf, numcalledhere=0, mytimingout
  real*8 ::  rsum, rsum2, time, rsum3

  if (myrank.ne.1) then
     return
  endif
  numcalledhere=numcalledhere+1

  if (incmf==1) then
     mytimingout=1
  else
     mytimingout=timingout
  endif
  if (numcalledhere==1) then
     open(853,file=dendatfile,status="unknown")
     write(853,*) ;           close(853)
     open(853,file=denrotfile,status="unknown")
     write(853,*) ;           close(853)
  endif
  
  if (mod(numcalledhere,mytimingout).eq.0)  then
     rsum=0.d0;  rsum2=0.d0;  rsum3=0.d0
     do i=1,nspf
        do j=1,nspf
           if (i /= j) then
              rsum=rsum+abs(yyy%denmat(i,j,0)**2)
           else
              rsum2=rsum2+abs(yyy%denmat(i,j,0)**2)
           endif
           if (shells(i) /= shells(j)) then
              rsum3=rsum3+abs(yyy%denmat(i,j,0)**2)
           endif
        enddo
     enddo

     open(853,file=dendatfile,status="unknown", position="append")
     write(853,'(F17.8, 100E17.7)') time, yyy%denvals(1:nspf); close(853)

     open(853,file=denrotfile,status="unknown", position="append")
     write(853,'(100F23.14)') time, rsum/rsum2, rsum3/rsum2;     close(853)
           
!!$     if (cdenflag.ne.0.and.rdenflag.ne.0) then
!!$        !!     call openfile(); write(mpifileptr,*) "Calling schmidt check"; call closefile()
!!$        call schmidtcheck(1)
!!$     endif
  endif
end subroutine output_denmat


subroutine denprint(iprop)
  use parameters
  use xxxmod
  implicit none
  integer :: i,iprop
  CNORMTYPE ::  sum2

  sum2=0d0
  do i=1,nspf
     sum2=sum2+yyy%denvals(i)
  enddo
  if (sum2.eq.0.d0) then
     print *, "SUM2 ERR!!"
  endif
  if (iprop.eq.1) then
     OFL; WRFL "|-------Natural orbital occupation numbers for each state-------------------|"
  endif
  write(mpifileptr,'(100F9.5)') (abs(yyy%denvals(i)/sum2*numelec),i=1,nspf); CFL

end subroutine denprint

end module denutilmod


module natrepbiomod
  use biorthotypemod
  implicit none
  type(biorthotype),target :: natrepbiovar
end module

module repnatmod
contains

subroutine replace_withnat(printflag)
  use natrepbiomod
  use biorthomod
  use parameters
  use configmod    !! bioww
  use xxxmod
  use spfsubmod
  implicit none

  DATATYPE,allocatable :: outspfs(:,:)
  integer :: i,j,printflag,imc
  real*8 :: errorval
  DATATYPE,target :: smo(nspf,nspf)

  allocate(outspfs(spfsize,nspf))
  outspfs=0d0

  do j=1,nspf  ! which natorb
     do i=1,nspf  ! which original
        outspfs(:,j)=outspfs(:,j)+ &
             yyy%cmfspfs((i-1)*spfsize+1:i*spfsize,0)*yyy%denvects(i,j)
     enddo
  enddo

  call spf_orthogit(outspfs, errorval)
  if (errorval.gt.1d-7) then
     OFLWR "ERROR IN REPLACENAT ", errorval; CFLST
  endif

  if (printflag==1) then
     OFLWR "REPLACING SPFS - DENMAT EIGS"
     do i=1,nspf
        write(mpifileptr,'(2E25.10)') yyy%denvals(i)/mcscfnum
     enddo
     WRFL; CFL
  endif

  call bioset(natrepbiovar,smo,numr,bioww)

  do imc=1,mcscfnum

!!$     tempavector(:)=yyy%cmfavec(:,imc,0)

     call biotransform(yyy%cmfspfs(:,0),outspfs, yyy%cmfavec(:,imc,0),natrepbiovar)

!!$! doesn't work.  permoverlaps broken presumably.
!!$#ifdef NOWAYDUDE
!!$     call autocorrelate_one(yyy%cmfavec(:,imc,0), outspfs, yyy%cmfspfs(:,0), tempavector, sum,1,needbiovar)
!!$!! CHECK
!!$     tempavector2(:)=yyy%cmfavec(:,imc,0)
!!$     call biotransform(outspfs,yyy%cmfspfs(:,0),tempavector2(:),natrepbiovar)
!!$     print *, DOT_PRODUCT(tempavector,tempavector),"AAAAAA"
!!$     tempavector2(:)=     tempavector2(:) - tempavector(:)
!!$     print *, DOT_PRODUCT(tempavector2,tempavector2),"BBBBB";
!!$     nullspfs(:,:)=0d0
!!$nullspfs(:,:)=nullspfs(:,:)-outspfs(:,:)
!!$print *, "BIOCHECK ", DOT_PRODUCT(RESHAPE(nullspfs, (/spfsize*nspf/)), RESHAPE(nullspfs, (/spfsize*nspf/)))
!!$!     call checkbio(yyy%cmfspfs(:,0),outspfs,tempavector,yyy%cmfavec(:,imc,0))
!!$!     call checkbio(yyy%cmfspfs(:,0),nullspfs,tempavector,yyy%cmfavec(:,imc,0))
!!$#endif

  enddo

  yyy%cmfspfs(:,0)=RESHAPE(outspfs,(/totspfdim/))

  deallocate(outspfs)

end subroutine replace_withnat

end module


!! denmat is the true denmat, not transposed.

module densubmod
contains

subroutine getdenmat00(inholeflag,wwin,avector1,in_avector2,rvector, denmat, numpoints,howmany)
  use walkmod
  use dotmod
  use mpisubmod
  implicit none
  type(walktype),intent(in) :: wwin
  integer,intent(in) ::  numpoints,howmany,inholeflag
  DATATYPE, intent(in) :: in_avector2(numpoints,wwin%firstconfig:wwin%lastconfig,howmany),&
       avector1(numpoints,wwin%firstconfig:wwin%lastconfig,howmany)
  DATAECS,intent(in) :: rvector(numpoints)
  DATATYPE,intent(out) :: denmat(wwin%nspf,wwin%nspf)
  DATATYPE :: a1(numpoints,howmany), a2(numpoints,howmany), &
       mydenmat(wwin%nspf,wwin%nspf), csum
  DATATYPE, allocatable :: avector2(:,:,:)
  integer :: config1,config2,  ispf,jspf,  dirphase, iwalk, ii, ihop

  allocate(avector2(numpoints,wwin%numconfig,howmany))
  avector2(:,:,:)=0d0
  if (wwin%lastconfig.ge.wwin%firstconfig) then
     avector2(:,wwin%firstconfig:wwin%lastconfig,:) = in_avector2(:,:,:)
  endif

!! DO SUMMA (parconsplit.ne.0 and sparsesummaflag.eq.2, "circ")

  if (wwin%parconsplit.ne.0) then
     do ii=1,howmany
        call mpiallgather(avector2(:,:,ii),wwin%numconfig*numpoints,&
             wwin%configsperproc(:)*numpoints,wwin%maxconfigsperproc*numpoints)
     enddo
  endif

  denmat(:,:)=0.d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(config1,ii,iwalk,config2,dirphase,ispf,jspf,mydenmat,a1,a2,csum,ihop)

!! removed  REDUCTION(+:denmat) after problems in get_twoereduced, OMP CRITICAL instead

  mydenmat(:,:)=0.d0

!$OMP DO SCHEDULE(DYNAMIC)
  do config1=wwin%botconfig,wwin%topconfig

     do ii=1,howmany
        a1(:,ii)=avector1(:,config1,ii) * rvector(:)
     enddo

     do ihop=1,wwin%numsinglehops(config1)
        config2=wwin%singlehop(ihop,config1)

        a2(:,:)=avector2(:,config2,:)
        csum=dot(a2(:,:),a1(:,:),numpoints*howmany)

        do iwalk=wwin%singlehopwalkstart(ihop,config1),wwin%singlehopwalkend(ihop,config1)

           dirphase=wwin%singlewalkdirphase(iwalk+wwin%scol(config1))
           ispf=wwin%singlewalkopspf(1,iwalk+wwin%scol(config1))  !! goes with config1
           jspf=wwin%singlewalkopspf(2,iwalk+wwin%scol(config1))  !! goes with config2

           mydenmat(ispf,jspf)=mydenmat(ispf,jspf)+ &
                csum*dirphase

!!$!!  ONCE AND FOR ALL 2014           :                     RIGHT MULTIPLYING!!!
!!$           denmat(ispf,jspf)=denmat(ispf,jspf)+ &
!!$                dirphase*a1*CONJUGATE(a2)
!!$!! NOT
!!$!!           denmat(jspf,ispf)=denmat(jspf,ispf)+ &
!!$!!                dirphase*a1*CONJUGATE(a2)

        enddo

     enddo
  enddo
!$OMP END DO
!$OMP CRITICAL
  denmat(:,:)=denmat(:,:)+mydenmat(:,:)
!$OMP END CRITICAL
!$OMP END PARALLEL

  call mympireduce(denmat,wwin%nspf**2)

  if (wwin%holeflag.ne.0.and.inholeflag.ne.0) then
     csum=0d0
     if (wwin%lastconfig.ge.wwin%firstconfig) then
        csum=dot(in_avector2,avector1,numpoints*wwin%localnconfig*howmany)
     endif
     if (wwin%parconsplit.ne.0) then
        call mympireduceone(csum)
     endif
     denmat(:,:) = (-1) * denmat(:,:)
     do ispf=1,wwin%nspf
        denmat(ispf,ispf)=denmat(ispf,ispf) + csum * 2d0
     enddo
  endif

  deallocate(avector2)


end subroutine getdenmat00

end module densubmod


module denmatxmod
contains

subroutine getdenmatx()
  use parameters
  use configmod
  use xxxmod
  implicit none

  call getdenmatstuff(www,yyy%cmfavec(:,:,0), yyy%denmat(:,:,0) , &
       yyy%invdenmat(:,:,0) , yyy%denvals(:) , yyy%denvects(:,:), numr, mcscfnum)

!!  OFLWR "DENCHECK"
!!#ifdef REALGO
!!  write(mpifileptr,'(F20.10)') yyy%denmat(:,1,0)
!!#else
!!  write(mpifileptr,'(2F20.10)') yyy%denmat(:,1,0)
!!#endif
!!  WRFL "DENSTOP!"; CFLST



!!$  if (rdenflag==1) then
!!$     call getrdenmat()
!!$  endif
!!$  if (cdenflag==1) then
!!$     call getnatconfig()
!!$  endif
!!$  if (cdenflag==1.and.rdenflag==1) then
!!$     call schmidtcheck(0)
!!$  endif

contains
subroutine getdenmatstuff(wwin,avector, denmat, invdenmat, denvals, &
     denvects, numpoints,howmany)
  use class_parameters
  use denreg_parameters
  use walkmod
  use invsubmod
  use densubmod
  use eigenmod
  implicit none
  type(walktype),intent(in) :: wwin
  integer,intent(in) ::  numpoints,howmany
  CNORMTYPE,intent(out) :: denvals(wwin%nspf)
  DATATYPE, intent(in) :: avector(numpoints,wwin%firstconfig:wwin%lastconfig,howmany)
  DATATYPE,intent(out) :: denmat(wwin%nspf,wwin%nspf),&
       invdenmat(wwin%nspf,wwin%nspf),denvects(wwin%nspf,wwin%nspf)
  CNORMTYPE :: tempdenvals(wwin%nspf,numclasses)
  DATATYPE :: tempinvden(wwin%nspf,wwin%nspf,numclasses),&
       tempdenvects(wwin%nspf,wwin%nspf,numclasses), csum
  integer :: ispf,jspf,iclass
  DATAECS :: rvector(numpoints)

  rvector(:)=1d0
  call getdenmat00(1,wwin,avector,avector,rvector,denmat,numpoints,howmany)

  denvects(:,:)=0d0; invdenmat(:,:)=0d0

  tempinvden(:,:,:)=0d0; tempdenvals(:,:)=0d0; tempdenvects(:,:,:)=0d0; 

  do iclass=1,numclasses
     do ispf=1,nperclass(iclass)
        do jspf=1,nperclass(iclass)
           tempinvden(ispf,jspf,iclass) = &
                denmat(classorb(ispf,iclass),classorb(jspf,iclass))
        enddo
     enddo
     tempinvden(:,:,iclass)= (-1) * tempinvden(:,:,iclass)
     call EIGEN(tempinvden(:,:,iclass),nperclass(iclass), wwin%nspf, &
          tempdenvects(:,:,iclass),tempdenvals(:,iclass))
     tempinvden(:,:,iclass)= (-1) * tempinvden(:,:,iclass)
     tempdenvals(:,iclass)= (-1) * tempdenvals(:,iclass)

     call invmatsmooth(tempinvden(:,:,iclass),nperclass(iclass),wwin%nspf,denreg,.false.)
     do ispf=1,nperclass(iclass)
        denvals(classorb(ispf,iclass))=tempdenvals(ispf,iclass)
        do jspf=1,nperclass(iclass)
           denvects(classorb(ispf,iclass),classorb(jspf,iclass)) = &
                tempdenvects(ispf,jspf,iclass)
           invdenmat(classorb(ispf,iclass),classorb(jspf,iclass)) = &
                tempinvden(ispf,jspf,iclass)
        enddo
     enddo
  enddo

  do ispf=1,wwin%nspf
     csum=denvects(ispf,ispf)
#ifdef CNORMFLAG
     if (real(csum,8).lt.0d0) then
        denvects(:,ispf) = denvects(:,ispf) * (-1)
     endif
#else
     if (csum.ne.0d0) then
        denvects(:,ispf) = denvects(:,ispf) * abs(csum) / csum
     endif
#endif
  enddo

end subroutine getdenmatstuff

end subroutine getdenmatx

end module denmatxmod



module getoccmod
contains

subroutine getoccupations(wwin,in_avector, numpoints, occupations)
  use walkmod
  use dotmod
  use mpisubmod
  implicit none
  type(walktype),intent(in) :: wwin
  integer,intent(in) ::  numpoints
  CNORMTYPE, intent(out) :: occupations(wwin%nspf)
  integer :: config1,  ispf,jspf,  iwalk,idiag
  DATATYPE, intent(in) :: in_avector(numpoints,wwin%firstconfig:wwin%lastconfig)
  DATATYPE, allocatable :: avector(:,:)

  allocate(avector(numpoints,wwin%numconfig))
  avector(:,:)=0d0

  if (wwin%lastconfig.ge.wwin%firstconfig) then
     avector(:,wwin%firstconfig:wwin%lastconfig) = in_avector(:,:)
  endif

!! DO SUMMA (parconsplit.ne.0 and sparsesummaflag.eq.2, "circ")

  if (wwin%parconsplit.ne.0) then
     call mpiallgather(avector,wwin%numconfig*numpoints,&
          wwin%configsperproc(:)*numpoints,wwin%maxconfigsperproc*numpoints)
  endif

  occupations(:)=0d0

  do config1=wwin%botconfig,wwin%topconfig
     do idiag=1,wwin%numsinglediagwalks(config1)
        iwalk=wwin%singlediag(idiag,config1)

        ispf=wwin%singlewalkopspf(1,iwalk+wwin%scol(config1))  !! goes with config1
        jspf=wwin%singlewalkopspf(2,iwalk+wwin%scol(config1))  !! goes with config2

        occupations(ispf)=occupations(ispf) + &  !! ok conversion
             dot(avector(:,config1),avector(:,config1),numpoints)

     enddo
  enddo

  deallocate(avector)

#ifndef REALGO
#ifndef CNORMFLAG
  call mympirealreduce(occupations,wwin%nspf)
#else
  call mympireduce(occupations,wwin%nspf)
#endif
#else
  call mympireduce(occupations,wwin%nspf)
#endif

  if (wwin%holeflag.ne.0) then
     occupations(:)=2d0-occupations(:)
  endif

end subroutine getoccupations

end module getoccmod

!!$subroutine getrdenmat()
!!$  use parameters
!!$  use xxxmod
!!$  implicit none
!!$
!!$  DATATYPE :: avector(numr,numconfig,mcscfnum), a1(mcscfnum), a2(mcscfnum), &
!!$       csum, csum2, sum,dot
!!$  integer ::     ii,jj, config1
!!$
!!$  avector(:,:,:)=RESHAPE(yyy%cmfavec(:,:,0),(/numr,numconfig,mcscfnum/))
!!$  do ii=1,numr
!!$     do jj=1,numr
!!$        sum=0.d0
!!$        do config1=1,numconfig  !! ok walks
!!$           a1(:)=avector(ii,config1,:);           
!!$           a2(:)=avector(jj,config1,:)
!!$           sum=sum+dot(a2(:),a1(:),mcscfnum)
!!$        enddo
!!$        yyy%rdenmat(ii,jj)=sum
!!$     enddo
!!$  enddo
!!$  
!!$  yyy%rdenmat=(-1)*yyy%rdenmat
!!$  call EIGEN(yyy%rdenmat,numr, numr, yyy%rdenvects(:,:), yyy%rdenvals(:))
!!$  yyy%rdenvals(:)=(-1)*yyy%rdenvals(:) ;  yyy%rdenmat=(-1)*yyy%rdenmat     
!!$  
!!$  do ii=1,numr
!!$     csum=0.d0;     csum2=0.d0
!!$     do jj=1,numr
!!$        
!!$!! 111510 this is just for output, no big deal, but rdenexpect should be 
!!$        CNORMTYPE, really DATATYPE if want complex expect for chmctdh
!!$
!!$        csum=csum+              ALLCON(yyy%rdenvects(jj,ii))*     yyy%rdenvects(jj,ii)  
!!$        csum2=csum2+ ALLCON(yyy%rdenvects(jj,ii))*yyy%rdenvects(jj,ii)*bondpoints(ii) 
!!$     enddo
!!$     yyy%rdenexpect(ii)=csum2/csum   !! ok for implicit conv (chmctdh, pmctdh, mctdh)
!!$  enddo
!!$
!!$end subroutine getrdenmat




!!$  KEEPME KEEPME KEEMPE
!!$
!!$!! checks whether projections of natural configurations are 
!!$!! proportional to natural orbitals in r
!!$!! only valid for mcscfnum=1....   projections are for each state,
!!$!! r natorbs are averaged
!!$
!!$subroutine schmidtcheck(schmidtprint)
!!$  use parameters
!!$  use xxxmod
!!$  use natprojmod
!!$  implicit none
!!$
!!$  DATATYPE :: dot
!!$  real*8 :: rsum
!!$  integer ::  i,j,ssize, schmidtprint
!!$  DATATYPE, allocatable :: schmidt(:,:)
!!$
!!$  if (numr.eq.1.or.mcscfnum.ne.1) then
!!$     return
!!$  endif
!!$  do i=1,numr
!!$     if (abs(yyy%rdenvals(i)).lt.1.d-8) then
!!$        exit
!!$     endif
!!$  enddo
!!$  ssize=i-1;  allocate(schmidt(ssize,ssize))
!!$  do i=1,ssize
!!$     do j=1,ssize
!!$        schmidt(i,j)=dot(natproj(:,i,1),yyy%rdenvects(:,j),numr)
!!$     enddo
!!$  enddo
!!$  if (schmidtprint.ne.0) then
!!$     call openfile()
!!$     write(mpifileptr,*) "------- Schmidt Overlaps ----------"
!!$     do i=1,ssize
!!$        write(mpifileptr,'(100E12.4)') schmidt(:,i)
!!$     enddo
!!$     write(mpifileptr,*) "-----------------------------------"
!!$  endif
!!$  do i=1,ssize
!!$     do j=1,ssize
!!$        rsum=abs(schmidt(i,j))**2
!!$        if (i==j) then 
!!$#ifndef CNORMFLAG
!!$           rsum=rsum-yyy%rdenvals(i)
!!$        endif
!!$#else
!!$        else
!!$#endif
!!$           if (abs(rsum).gt.1.d-12) then
!!$              OFLWR "WARNING!  Schmidt dot fail", rsum, i,j,schmidt(i,j),yyy%rdenvals(i);CFL
!!$           endif
!!$#ifdef CNORMFLAG
!!$        endif
!!$#endif
!!$     enddo
!!$  enddo
!!$
!!$  deallocate(schmidt)
!!$end subroutine schmidtcheck
!!$
!!$
!!$subroutine getnatconfig()
!!$  use parameters
!!$  use natprojmod
!!$  use xxxmod
!!$  implicit none
!!$
!!$  DATATYPE :: avector(numr,numconfig,mcscfnum), a1(mcscfnum), a2(mcscfnum),dot
!!$  integer ::  config1,config2,     ii,imc,i,ir,numvects
!!$  DATATYPE :: natmat(numconfig,numconfig), natvects(numconfig,numconfig),tempconfig(numconfig)
!!$
!!$  avector(:,:,:)=RESHAPE(yyy%cmfavec(:,:,0),(/numr,numconfig,numr,mcscfnum/))
!!$
!!$  if (numconfig.gt.2000) then
!!$     OFLWR "Error: getting natconfigs but numconfig gt 2000: ", numconfig; CFLST
!!$  endif
!!$  natmat=0.d0
!!$  do ii=1,numr
!!$     do config1=1,numconfig    !! ok walks
!!$        a1(:)=avector(ii,config1,:)
!!$        do config2=1,numconfig !! ok walks
!!$           a2(:)=avector(ii,config2,:)
!!$           natmat(config1,config2) = natmat(config1,config2) + dot(a2,a1,mcscfnum)
!!$        enddo
!!$     enddo
!!$  enddo
!!$
!!$  natmat=(-1)*natmat
!!$  call EIGEN(natmat,numconfig, numconfig, natvects, natvals)
!!$
!!$  natmat=(-1)*natmat;  natvals=(-1)*natvals
!!$  numvects=min(numr,numconfig);  natconfigs=0.d0
!!$
!!$  do imc=1,mcscfnum
!!$     do i=1,numvects
!!$        do ir=1,numr
!!$        tempvector(:)=avector(ir,:,imc)
!!$           natproj(ir,i,imc)=dot(natvects(:,i),  & !! ok imp conv (p,ch)
!!$                tempvector(:),numconfig)
!!$        enddo
!!$        natdot(i,imc)=dot(natproj(:,i,imc),       & !! ok imp conv (p,ch)
!!$             natproj(:,i,imc),numr)
!!$     enddo
!!$  enddo
!!$
!!$  natconfigs(:,1:numvects)=natvects(:,1:numvects)
!!$
!!$end subroutine getnatconfig



