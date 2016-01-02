
!! DENSITY MATRIX SUBROUTINES !!  
  
#include "Definitions.INC"
  

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
     open(853,file=dendatfile,status="unknown");   write(853,*) ;           close(853)
     open(853,file=denrotfile,status="unknown");     write(853,*) ;           close(853)
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

module natrepbiomod
  use biorthotypemod
  implicit none
  type(biorthotype),target :: natrepbiovar
end module

subroutine replace_withnat(printflag)
  use natrepbiomod
  use biorthomod
  use parameters
  use configmod
  use xxxmod
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
             yyy%cmfpsivec(spfstart+(i-1)*spfsize:spfstart+i*spfsize-1,0)*yyy%denvects(i,j)
     enddo
  enddo

  call spf_orthogit(outspfs, errorval)
  if (errorval.gt.1d-7) then
     OFLWR "WTF!  ERROR IN REPLACENAT ", errorval; CFLST
  endif

  if (printflag==1) then
     OFLWR "REPLACING SPFS - DENMAT EIGS"
     do i=1,nspf
        write(mpifileptr,'(2E25.10)') yyy%denvals(i)/mcscfnum  !!!*numelec
     enddo
     WRFL; CFL
  endif

     call bioset(natrepbiovar,smo,numr,bwwptr)

  do imc=1,mcscfnum

!!$     tempavector(:)=yyy%cmfpsivec(astart(imc):aend(imc),0)

     call biotransform(yyy%cmfpsivec(spfstart,0),outspfs,yyy%cmfpsivec(astart(imc),0),natrepbiovar)

!!$! doesn't work.  permoverlaps broken presumably.
!!$#ifdef NOWAYDUDE
!!$     call autocorrelate_one(yyy%cmfpsivec(astart(imc),0), outspfs, yyy%cmfpsivec(spfstart,0), tempavector, sum,1)
!!$!! CHECK
!!$     tempavector2(:)=yyy%cmfpsivec(astart(imc):aend(imc),0)
!!$     call biotransform(outspfs,yyy%cmfpsivec(spfstart,0),tempavector2(:),natrepbiovar)
!!$     print *, DOT_PRODUCT(tempavector,tempavector),"AAAAAA"
!!$     tempavector2(:)=     tempavector2(:) - tempavector(:)
!!$     print *, DOT_PRODUCT(tempavector2,tempavector2),"BBBBB";
!!$     nullspfs(:,:)=0d0
!!$nullspfs(:,:)=nullspfs(:,:)-outspfs(:,:)
!!$print *, "BIOCHECK ", DOT_PRODUCT(RESHAPE(nullspfs, (/spfsize*nspf/)), RESHAPE(nullspfs, (/spfsize*nspf/)))
!!$!     call checkbio(yyy%cmfpsivec(spfstart,0),outspfs,tempavector,yyy%cmfpsivec(astart(imc),0))
!!$!     call checkbio(yyy%cmfpsivec(spfstart,0),nullspfs,tempavector,yyy%cmfpsivec(astart(imc),0))
!!$#endif

  enddo

  yyy%cmfpsivec(spfstart:spfend,0)=RESHAPE(outspfs,(/totspfdim/))

  deallocate(outspfs)

end subroutine replace_withnat



subroutine getdenmatx()
  use parameters
  use configmod
  use xxxmod
  implicit none

  call getdenmatstuff(www,yyy%cmfpsivec(astart(1),0), yyy%denmat(:,:,0) , yyy%invdenmat(:,:,0) , yyy%denvals(:) , yyy%denvects(:,:), numr, mcscfnum)

!!$  if (rdenflag==1) then
!!$     call getrdenmat()
!!$  endif
!!$  if (cdenflag==1) then
!!$     call getnatconfig()
!!$  endif
!!$  if (cdenflag==1.and.rdenflag==1) then
!!$     call schmidtcheck(0)
!!$  endif

end subroutine getdenmatx

!! denmat is the true denmat, not transposed.

subroutine getdenmat00(www,avector1,in_avector2,rvector, denmat, numpoints,howmany)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) ::  numpoints,howmany
  DATATYPE, intent(in) :: in_avector2(numpoints,www%firstconfig:www%lastconfig,howmany),&
       avector1(numpoints,www%firstconfig:www%lastconfig,howmany)
  DATAECS,intent(in) :: rvector(numpoints)
  DATATYPE,intent(out) :: denmat(www%nspf,www%nspf)
  DATATYPE :: a1(numpoints,howmany), a2(numpoints,howmany), dot
  DATATYPE, allocatable :: avector2(:,:,:)
  integer :: config1,config2,  ispf,jspf,  dirphase, iwalk, ii

  allocate(avector2(numpoints,www%numconfig,howmany))
  avector2(:,:,:)=0d0

  avector2(:,www%firstconfig:www%lastconfig,:) = in_avector2(:,:,:)

!! DO SUMMA (parconsplit.ne.0 and sparsesummaflag.eq.2, "circ")
!! AND DO HOPS

  if (www%parconsplit.ne.0) then
     do ii=1,howmany
        call mpiallgather(avector2(:,:,ii),www%numconfig*numpoints,www%configsperproc(:)*numpoints,www%maxconfigsperproc*numpoints)
     enddo
  endif

  denmat(:,:)=0.d0

  !! single off diagonal walks

  do config1=www%botconfig,www%topconfig

     do ii=1,howmany
        a1(:,ii)=avector1(:,config1,ii) * rvector(:)
     enddo

     do iwalk=1,www%numsinglewalks(config1)
        config2=www%singlewalk(iwalk,config1)
        dirphase=www%singlewalkdirphase(iwalk,config1)
        a2(:,:)=avector2(:,config2,:)
        ispf=www%singlewalkopspf(1,iwalk,config1)  !! goes with config1
        jspf=www%singlewalkopspf(2,iwalk,config1)  !! goes with config2


        denmat(ispf,jspf)=denmat(ispf,jspf)+ &
             dot(a2(:,:),a1(:,:),numpoints*howmany)*dirphase

!!$!!  ONCE AND FOR ALL 2014           :                     RIGHT MULTIPLYING!!!
!!$           denmat(ispf,jspf)=denmat(ispf,jspf)+ &
!!$                dirphase*a1*CONJUGATE(a2)
!!$!! NOT
!!$!!           denmat(jspf,ispf)=denmat(jspf,ispf)+ &
!!$!!                dirphase*a1*CONJUGATE(a2)

     enddo
  enddo

  deallocate(avector2)

  call mympireduce(denmat,www%nspf**2)

end subroutine getdenmat00



subroutine getdenmatstuff(www,avector, denmat, invdenmat, denvals, denvects, numpoints,howmany)
  use class_parameters
  use denreg_parameters
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) ::  numpoints,howmany
  CNORMTYPE,intent(out) :: denvals(www%nspf)
  DATATYPE, intent(in) :: avector(numpoints,www%firstconfig:www%lastconfig,howmany)
  DATATYPE,intent(out) :: denmat(www%nspf,www%nspf),invdenmat(www%nspf,www%nspf),denvects(www%nspf,www%nspf)
  CNORMTYPE :: tempdenvals(www%nspf,numclasses)
  DATATYPE :: tempinvden(www%nspf,www%nspf,numclasses),tempdenvects(www%nspf,www%nspf,numclasses)
  integer :: ispf,jspf,iclass
  DATAECS :: rvector(numpoints)

  rvector(:)=1d0
  call getdenmat00(www,avector,avector,rvector,denmat,numpoints,howmany)

  denvects(:,:)=0d0; invdenmat(:,:)=0d0

  tempinvden(:,:,:)=0d0; tempdenvals(:,:)=0d0; tempdenvects(:,:,:)=0d0; 

  do iclass=1,numclasses
     do ispf=1,nperclass(iclass)
        do jspf=1,nperclass(iclass)
           tempinvden(ispf,jspf,iclass)=denmat(classorb(ispf,iclass),classorb(jspf,iclass))
        enddo
     enddo
     tempinvden(:,:,iclass)= (-1) * tempinvden(:,:,iclass)
     call EIGEN(tempinvden(:,:,iclass),nperclass(iclass), www%nspf, tempdenvects(:,:,iclass),tempdenvals(:,iclass))
     tempinvden(:,:,iclass)= (-1) * tempinvden(:,:,iclass)
     tempdenvals(:,iclass)= (-1) * tempdenvals(:,iclass)

     call invmatsmooth(tempinvden(:,:,iclass),nperclass(iclass),www%nspf,denreg)
     do ispf=1,nperclass(iclass)
        denvals(classorb(ispf,iclass))=tempdenvals(ispf,iclass)
        do jspf=1,nperclass(iclass)
           denvects(classorb(ispf,iclass),classorb(jspf,iclass))=tempdenvects(ispf,jspf,iclass)
           invdenmat(classorb(ispf,iclass),classorb(jspf,iclass))=tempinvden(ispf,jspf,iclass)
        enddo
     enddo
  enddo


end subroutine getdenmatstuff


subroutine getoccupations(www,in_avector, numpoints, occupations)
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) ::  numpoints
  CNORMTYPE, intent(out) :: occupations(www%nspf)
  integer :: config1,  ispf,jspf,  iwalk,idiag
  DATATYPE, intent(in) :: in_avector(numpoints,www%firstconfig:www%lastconfig)
  DATATYPE :: dot
  DATATYPE, allocatable :: avector(:,:)

  allocate(avector(numpoints,www%numconfig))
  avector(:,:)=0d0

  avector(:,www%firstconfig:www%lastconfig) = in_avector(:,:)

!! DO SUMMA (parconsplit.ne.0 and sparsesummaflag.eq.2, "circ")
!! AND DO HOPS

  if (www%parconsplit.ne.0) then
     call mpiallgather(avector,www%numconfig*numpoints,www%configsperproc(:)*numpoints,www%maxconfigsperproc*numpoints)
  endif

  occupations(:)=0d0

  do config1=www%botconfig,www%topconfig
     do idiag=1,www%numsinglediagwalks(config1)
        iwalk=www%singlediag(idiag,config1)

        ispf=www%singlewalkopspf(1,iwalk,config1)  !! goes with config1
        jspf=www%singlewalkopspf(2,iwalk,config1)  !! goes with config2

        occupations(ispf)=occupations(ispf) + dot(avector(:,config1),avector(:,config1),numpoints)

     enddo
  enddo

  deallocate(avector)

#ifndef REALGO
#ifndef CNORMFLAG
  call mympirealreduce(occupations,www%nspf)
#else
  call mympireduce(occupations,www%nspf)
#endif
#else
  call mympireduce(occupations,www%nspf)
#endif


end subroutine getoccupations


!!$subroutine getrdenmat()
!!$  use parameters
!!$  use xxxmod
!!$  implicit none
!!$
!!$  DATATYPE :: avector(numconfig,numr,mcscfnum), a1(mcscfnum), a2(mcscfnum), csum, csum2, sum,dot
!!$  integer ::     ii,jj, config1
!!$
!!$  avector(:,:,:)=RESHAPE(yyy%cmfpsivec(astart(1):aend(mcscfnum),0),(/numconfig,numr,mcscfnum/))
!!$  do ii=1,numr
!!$     do jj=1,numr
!!$        sum=0.d0
!!$        do config1=1,numconfig  !! ok walks
!!$           a1(:)=avector(config1,ii,:);           a2(:)=avector(config1,jj,:)
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
!!$!! 111510 this is just for output, no big deal, but rdenexpect should be CNORMTYPE, really DATATYPE if want complex expect for chmctdh
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
!!$!! checks whether projections of natural configurations are proportional to natural orbitals in r
!!$!! only valid for mcscfnum=1....   projections are for each state, r natorbs are averaged
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
!!$  DATATYPE :: avector(numconfig,numr,mcscfnum), a1(mcscfnum), a2(mcscfnum),dot
!!$  integer ::  config1,config2,     ii,imc,i,ir,numvects
!!$  DATATYPE :: natmat(numconfig,numconfig), natvects(numconfig,numconfig)
!!$
!!$  avector(:,:,:)=RESHAPE(yyy%cmfpsivec(astart(1):aend(mcscfnum),0),(/numconfig,numr,mcscfnum/))
!!$
!!$  if (numconfig.gt.2000) then
!!$     OFLWR "Error: getting natconfigs but numconfig gt 2000: ", numconfig; CFLST
!!$  endif
!!$  natmat=0.d0
!!$  do ii=1,numr
!!$     do config1=1,numconfig    !! ok walks
!!$        a1(:)=avector(config1,ii,:)
!!$        do config2=1,numconfig !! ok walks
!!$           a2(:)=avector(config2,ii,:)
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
!!$           natproj(ir,i,imc)=dot(natvects(:,i),  & !! ok imp conv (p,ch)
!!$                avector(:,ir,imc),numconfig)
!!$        enddo
!!$        natdot(i,imc)=dot(natproj(:,i,imc),       & !! ok imp conv (p,ch)
!!$             natproj(:,i,imc),numr)
!!$     enddo
!!$  enddo
!!$
!!$  natconfigs(:,1:numvects)=natvects(:,1:numvects)
!!$
!!$end subroutine getnatconfig



