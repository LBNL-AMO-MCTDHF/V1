
!! ORBITAL AND CONFIGURATION MATRIX ELEMENTS.  TWO ELECTRON REDUCED MATRIX ELEMENTS
!! ARE CALCULATED HERE ALONG THE WAY.  GET_REDUCEDHAM MUST BE CALLED AFTER TWOE_MATEL
!! BECAUSE IT NEEDS THOSE REDUCED MATRIX ELEMENTS

#include "Definitions.INC"


recursive subroutine all_matel()
  use parameters
  use xxxmod
  use mpimod
  use opmod
  implicit none
  integer,save :: times(10)=0,xcalled=0
  integer :: itime,jtime,getlen
  xcalled=xcalled+1

  call all_matel0(yyy%cptr(0), yyy%cmfpsivec(spfstart,0), yyy%cmfpsivec(spfstart,0), twoereduced,times)

  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "     ...called all_matel0 in all_matel"; CFL; call mpibarrier()
  endif

  call system_clock(itime)
  if (sparseopt.ne.0) then
     call assemble_sparsemats(yyy%cptr(0),yyy%sptr(0),1,1,1,1)
  endif
  call system_clock(jtime); times(5)=times(5)+jtime-itime

  if ((myrank.eq.1).and.(notiming.eq.0)) then
     if (xcalled==1) then
        open(853, file=timingdir(1:getlen(timingdir)-1)//"/matel.time.dat", status="unknown")
!!$        open(8853, file=timingdir(1:getlen(timingdir)-1)//"/matel.abs.time.dat", status="unknown")
        write(853,'(100A11)')   "op", "pot", "pulse", "two", "assemble"
!!$        write(8853,'(100A11)')   "matel absolute timing"
        close(853)
!!$        close(8853)
     endif

     open(853, file=timingdir(1:getlen(timingdir)-1)//"/matel.time.dat", status="unknown", position="append")
     write(853,'(100I11)')  times(1:5)/1000;        close(853)
     close(853)
!!$     call system("date >> "//timingdir(1:getlen(timingdir)-1)//"/matel.abs.time.dat")
  endif

  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "     ...finished all_matel"; CFL; call mpibarrier()
  endif

end subroutine all_matel

recursive subroutine all_matel0(matrix_ptr,inspfs1,inspfs2,twoereduced,times)
  use parameters
  use configptrmod
  implicit none
  Type(CONFIGPTR) :: matrix_ptr
  DATATYPE :: inspfs1(spfsize,nspf), inspfs2(spfsize,nspf), twoereduced(reducedpotsize,nspf,nspf)
  integer :: times(*), i,j

  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "     ...sparseops_matel"; CFL; call mpibarrier()
  endif
  call system_clock(i)
  call sparseops_matel(matrix_ptr,inspfs1,inspfs2)
  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "     ...pot_matel"; CFL; call mpibarrier()
  endif
  call system_clock(j);  times(1)=times(1)+j-i; i=j
  call pot_matel(matrix_ptr,inspfs1,inspfs2)
  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "     ...pulse_matel"; CFL; call mpibarrier()
  endif
  call system_clock(j);  times(2)=times(2)+j-i; i=j
  call pulse_matel(matrix_ptr,inspfs1,inspfs2)
  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "     ...twoe_matel"; CFL; call mpibarrier()
  endif
  call system_clock(j);  times(3)=times(3)+j-i; i=j
  call twoe_matel(matrix_ptr,inspfs1,inspfs2,twoereduced)
  call system_clock(j);  times(4)=times(4)+j-i
  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "     ...done all_matel0"; CFL; call mpibarrier()
  endif

end subroutine all_matel0

recursive subroutine twoe_matel(matrix_ptr,inspfs1,inspfs2,twoereduced)
  use parameters
  use mpimod
  use configptrmod
  implicit none
  Type(CONFIGPTR) :: matrix_ptr
  DATATYPE :: inspfs1(spfsize,nspf),inspfs2(spfsize,nspf),twoereduced(reducedpotsize,nspf,nspf)

  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "         In twoe_matel.  Calling call_twoe_matel"; CFL; call mpibarrier()
  endif
  call call_twoe_matel(inspfs1,inspfs2,matrix_ptr%xtwoematel(:,:,:,:),twoereduced,timingdir,notiming)
  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "         Done call_twoe_matel.  Reduce"; CFL; call mpibarrier()
  endif
  if (parorbsplit.eq.3) then
     call mympireduce(matrix_ptr%xtwoematel(:,:,:,:),nspf**4)
  endif
  if (debugflag.eq.42) then
     call mpibarrier();     OFLWR "         Done reduce, done twoe_matel."; CFL; call mpibarrier()
  endif

end subroutine twoe_matel

subroutine pot_matel(matrix_ptr,inspfs1,inspfs2)
  use parameters
  use configptrmod
  implicit none
  Type(CONFIGPTR) :: matrix_ptr
  integer :: j
  DATATYPE :: inspfs1(spfsize,nspf), inspfs2(spfsize,nspf),  ttempspfs(spfsize,nspf),tempmatel(nspf,nspf)

  matrix_ptr%xpotmatel(:,:)=0d0

  do j=1,nspf
     call mult_pot(inspfs2(:,j),ttempspfs(:,j))
  enddo

!! potmatel is proper ordering.  fast index is conjg.

  call MYGEMM(CNORMCHAR,'N',nspf,nspf, spfsize, DATAONE, inspfs1, spfsize, ttempspfs, spfsize, DATAZERO, matrix_ptr%xpotmatel(:,:), nspf)

  call hatom_matel(inspfs1,inspfs2,tempmatel, nspf)

  matrix_ptr%xpotmatel(:,:)=matrix_ptr%xpotmatel(:,:)+ tempmatel(:,:)

  if (numfrozen.gt.0) then
!! only call of op_frozen_exchange.

     ttempspfs=0d0
     call op_frozen_exchange(inspfs2,ttempspfs)
     do j=1,nspf
        call op_frozenreduced(inspfs2(:,j),ttempspfs(:,j))   !! ADDS to ttempspfs
     enddo

     call MYGEMM(CNORMCHAR,'N',nspf,nspf, spfsize, DATAONE, inspfs1, spfsize, ttempspfs, spfsize, DATAONE, matrix_ptr%xpotmatel(:,:) ,nspf)
  endif

  if (parorbsplit.eq.3) then
     call mympireduce(matrix_ptr%xpotmatel(:,:),nspf**2)
  endif

end subroutine pot_matel



subroutine pulse_matel(matrix_ptr,inspfs1,inspfs2)
  use parameters
  use configptrmod
  implicit none
  Type(CONFIGPTR) :: matrix_ptr
  integer :: ispf
  DATATYPE :: inspfs1(spfsize,nspf), inspfs2(spfsize,nspf), nullcomplex(1), dipoles(3)
  DATATYPE :: ttempspfsxx(spfsize,nspf), ttempspfsyy(spfsize,nspf), ttempspfszz(spfsize,nspf)

  matrix_ptr%xpulsematelxx=0.d0;     matrix_ptr%xpulsematelyy=0.d0;     matrix_ptr%xpulsematelzz=0.d0

  if (velflag.eq.0) then

     do ispf=1,nspf
        call lenmultiply(inspfs2(:,ispf),ttempspfsxx(:,ispf), DATAONE,DATAZERO,DATAZERO)
        call lenmultiply(inspfs2(:,ispf),ttempspfsyy(:,ispf), DATAZERO,DATAONE,DATAZERO)
        call lenmultiply(inspfs2(:,ispf),ttempspfszz(:,ispf), DATAZERO,DATAZERO,DATAONE)
     enddo
  else
     do ispf=1,nspf
        call velmultiply(inspfs2(:,ispf),ttempspfsxx(:,ispf), DATAONE,DATAZERO,DATAZERO)
        call velmultiply(inspfs2(:,ispf),ttempspfsyy(:,ispf), DATAZERO,DATAONE,DATAZERO)
        call velmultiply(inspfs2(:,ispf),ttempspfszz(:,ispf), DATAZERO,DATAZERO,DATAONE)
     enddo
  endif
!! pulsematel is proper order.

  dipoles(:)=0d0
  if (velflag.eq.0) then
     call nucdipvalue(nullcomplex,dipoles)
  endif

  call MYGEMM(CNORMCHAR,'N',nspf,nspf, spfsize, DATAONE, inspfs1, spfsize, ttempspfsxx, spfsize, DATAONE, matrix_ptr%xpulsematelxx(:,:) ,nspf)
  call MYGEMM(CNORMCHAR,'N',nspf,nspf, spfsize, DATAONE, inspfs1, spfsize, ttempspfsyy, spfsize, DATAONE, matrix_ptr%xpulsematelyy(:,:) ,nspf)
  call MYGEMM(CNORMCHAR,'N',nspf,nspf, spfsize, DATAONE, inspfs1, spfsize, ttempspfszz, spfsize, DATAONE, matrix_ptr%xpulsematelzz(:,:) ,nspf)
  matrix_ptr%xpulsenuc(:)=dipoles(:)

  if (parorbsplit.eq.3) then
     call mympireduce(matrix_ptr%xpulsematelxx(:,:),nspf**2)
     call mympireduce(matrix_ptr%xpulsematelyy(:,:),nspf**2)
     call mympireduce(matrix_ptr%xpulsematelzz(:,:),nspf**2)
  endif

end subroutine pulse_matel


subroutine sparseops_matel(matrix_ptr,inspfs1,inspfs2)
  use parameters
  use configptrmod
  implicit none
  Type(CONFIGPTR) :: matrix_ptr
  integer :: ispf  
  DATATYPE :: inspfs1(spfsize,nspf), inspfs2(spfsize,nspf)  , ttempspfs(spfsize,nspf)  

  matrix_ptr%xopmatel=0d0;  matrix_ptr%xymatel=0d0;  

  if (multmanyflag.ne.0) then
     call mult_ke(inspfs2(:,:),ttempspfs(:,:),nspf,timingdir,notiming)
  else
     do ispf=1,nspf
        call mult_ke(inspfs2(:,ispf),ttempspfs(:,ispf),1,timingdir,notiming)
     enddo
  endif

  call MYGEMM(CNORMCHAR,'N',nspf,nspf, spfsize, DATAONE, inspfs1, spfsize, ttempspfs, spfsize, DATAZERO, matrix_ptr%xopmatel(:,:) ,nspf)

  if (nonuc_checkflag/=1) then
     call noparorbsupport("op_yderiv")
     do ispf=1,nspf
        call op_yderiv(inspfs2(:,ispf),ttempspfs(:,ispf))
     enddo
     call MYGEMM(CNORMCHAR,'N',nspf,nspf, spfsize, DATAONE, inspfs1, spfsize, ttempspfs, &
          spfsize, DATAZERO, matrix_ptr%xymatel(:,:) ,nspf)
  endif

  if (parorbsplit.eq.3) then
     call mympireduce(matrix_ptr%xopmatel(:,:),nspf**2)
  endif

end subroutine sparseops_matel


!! direct (not exchange)

subroutine frozen_matels()
  use parameters
  implicit none

  call call_frozen_matels()

end subroutine frozen_matels


!! ADDS TO OUTSPFS

subroutine op_frozen_exchange(inspfs,outspfs)
  use parameters
  implicit none
  DATATYPE :: inspfs(spfsize,nspf)
  DATATYPE :: outspfs(spfsize,nspf)

  call call_frozen_exchange(inspfs,outspfs)

end subroutine op_frozen_exchange


subroutine arbitraryconfig_matel00transpose(onebodymat, smallmatrixtr,topdim)
  use walkmod
  use parameters
  implicit none

  integer ::    config1,  iwalk,topdim,myind
  DATATYPE :: onebodymat(nspf,nspf), smallmatrixtr(topdim,botwalk:topwalk)

  if ( (sparseconfigflag.ne.0.and.topdim.ne.maxsinglewalks) .or. &
       (sparseconfigflag.eq.0.and.topdim.ne.numconfig) ) then
     OFLWR "BAD CALL 0909 ",sparseconfigflag,topdim,maxsinglewalks,numconfig; CFLST
  endif

  smallmatrixtr=0d0; 

!! 06-2015   do config1=botwalk,topwalk

  do config1=botconfig,topconfig
     do iwalk=1,numsinglewalks(config1)
        if (sparseconfigflag.eq.0) then
           myind=singlewalk(iwalk,config1)
        else
           myind=iwalk
        endif
        smallmatrixtr(myind,config1)=smallmatrixtr(myind,config1)+   &
             onebodymat(singlewalkopspf(1,iwalk,config1), &
             singlewalkopspf(2,iwalk,config1)) *  &
             singlewalkdirphase(iwalk,config1)
     enddo
  enddo

!! 06-2015

  if (sparseconfigflag.eq.0) then
     call mpiallgather(smallmatrixtr(:,:),numconfig**2,configsperproc(:)*numconfig,maxconfigsperproc*numconfig)
  endif

end subroutine arbitraryconfig_matel00transpose



subroutine arbitraryconfig_matel_doubles00transpose(twobodymat, smallmatrixtr,topdim)   
  use walkmod
  use parameters
  implicit none
  integer ::    config1, iwalk,topdim,myind
  DATATYPE :: twobodymat(nspf,nspf,nspf,nspf), smallmatrixtr(topdim,botwalk:topwalk)

  if ( (sparseconfigflag.ne.0.and.topdim.ne.maxdoublewalks) .or. &
       (sparseconfigflag.eq.0.and.topdim.ne.numconfig) ) then
     OFLWR "BAD CALL 09119 ",sparseconfigflag,topdim,maxdoublewalks,numconfig; CFLST
  endif

  smallmatrixtr=0d0;

!! 06-2015   do config1=botwalk,topwalk

  do config1=botconfig,topconfig
     do iwalk=1,numdoublewalks(config1)
        if (sparseconfigflag.eq.0) then
           myind=doublewalk(iwalk,config1)
        else
           myind=iwalk  
        endif
        smallmatrixtr(myind,config1)=smallmatrixtr(myind,config1) + &           
             twobodymat(doublewalkdirspf(1,iwalk,config1), &
             doublewalkdirspf(2,iwalk,config1),   &
             doublewalkdirspf(3,iwalk,config1),   &
             doublewalkdirspf(4,iwalk,config1))* &
             doublewalkdirphase(iwalk,config1) 
     enddo
  enddo

!! 06-2015

  if (sparseconfigflag.eq.0) then
     call mpiallgather(smallmatrixtr(:,:),numconfig**2,configsperproc(:)*numconfig,maxconfigsperproc*numconfig)
  endif

end subroutine arbitraryconfig_matel_doubles00transpose



subroutine assemble_spinconfigmat(outmatrix, matrix_ptr, boflag, nucflag, pulseflag, conflag, time)
  use parameters
  use configptrmod
  implicit none
  Type(configptr) :: matrix_ptr
  DATATYPE ::        outmatrix(spintotrank*numr,spintotrank*numr)
  DATATYPE, allocatable ::        bigmatrix(:,:),halfmatrix(:,:),halfmatrix2(:,:)
  integer :: nucflag, pulseflag, conflag, boflag
  real*8 :: time

  if (sparseconfigflag/=0) then
     OFLWR "Error, assemble_full_config_matel called when sparseconfigflag/=0";CFLST
  endif
  if (spinwalkflag.eq.0) then
     OFLWR "Error, assemble spin matel called but spinwalkflag is zero"; CFLST
  endif
  
  allocate(  bigmatrix(numconfig*numr,numconfig*numr), halfmatrix(spintotrank*numr,numconfig*numr),halfmatrix2(numconfig*numr,spintotrank*numr))

  call assemble_configmat(bigmatrix,matrix_ptr,boflag,nucflag,pulseflag,conflag,time)

!! PROJECT ON SPIN !!

!! obviously not good to project on spin for each r; should project earlier. 

  call configspin_transformto(numconfig*numr*numr,bigmatrix,halfmatrix2)
  halfmatrix=TRANSPOSE(halfmatrix2)
  call configspin_transformto(spintotrank*numr*numr,halfmatrix,outmatrix)
  outmatrix=transpose(outmatrix)

  deallocate(bigmatrix,halfmatrix,halfmatrix2)

end subroutine assemble_spinconfigmat



!! FOR NONSPARSE (HERE, ASSEMBLE_CONFIGMAT) INCLUDE A-SQUARED TERM.

subroutine assemble_configmat(bigconfigmat,matrix_ptr, boflag, nucflag, pulseflag, conflag,time)
  use parameters
  use configptrmod
  use opmod   !! rkemod, proderivmod
  implicit none
  integer :: conflag,boflag,nucflag,pulseflag,ir,jr,i
  DATATYPE :: bigconfigmat(numr,numconfig,numr,numconfig), tempmatel(nspf,nspf)
  DATATYPE, allocatable :: tempconfigmat(:,:),tempconfigmat2(:,:),diagmat(:,:,:)
  Type(CONFIGPTR) :: matrix_ptr
  real*8 :: time,gg
  DATATYPE :: facs(3), csum0,csum

  if (sparseconfigflag.ne.0) then
     OFLWR "BADDDCALLL"; CFLST
  endif

  allocate(tempconfigmat(numconfig,numconfig),tempconfigmat2(numconfig,numconfig), &
       diagmat(numconfig,numconfig,numr))

  diagmat(:,:,:)=0d0

  call vectdpot(time,velflag,facs)

  if (boflag==1) then
     do ir=1,numr
        do i=1,numconfig

         diagmat(i,i,ir)=( frozenkediag/bondpoints(ir)**2 + (nucrepulsion+frozenpotdiag)/bondpoints(ir) + energyshift ) * matrix_ptr%kefac 

        enddo
     enddo

     call arbitraryconfig_matel_doubles00transpose(matrix_ptr%xtwoematel(:,:,:,:),tempconfigmat,numconfig)
     call arbitraryconfig_matel00transpose(matrix_ptr%xpotmatel(:,:),tempconfigmat2,numconfig)

     tempconfigmat(:,:) = TRANSPOSE( tempconfigmat(:,:) + tempconfigmat2(:,:) )
     do ir=1,numr
        diagmat(:,:,ir)=diagmat(:,:,ir)+tempconfigmat(:,:)/bondpoints(ir)
     enddo
     call arbitraryconfig_matel00transpose(matrix_ptr%xopmatel(:,:),tempconfigmat,numconfig)
     tempconfigmat(:,:)=TRANSPOSE(tempconfigmat(:,:))
     do ir=1,numr
        diagmat(:,:,ir)=diagmat(:,:,ir)+tempconfigmat(:,:)/bondpoints(ir)**2
     enddo
  endif

  if (conflag==1.and.constraintflag.ne.0) then
     tempmatel(:,:)=matrix_ptr%xconmatel(:,:)
     if (tdflag.ne.0) then
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelxx(:,:)*facs(1)
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelyy(:,:)*facs(2)
        tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelzz(:,:)*facs(3)
     endif
     call arbitraryconfig_matel00transpose(tempmatel,tempconfigmat,numconfig)
     tempconfigmat(:,:)=TRANSPOSE(tempconfigmat(:,:))
     do ir=1,numr
        diagmat(:,:,ir)=diagmat(:,:,ir)-tempconfigmat(:,:)
     enddo
  endif

  if ((pulseflag.eq.1.and.tdflag.eq.1)) then
     tempmatel(:,:)= matrix_ptr%xpulsematelxx*facs(1)+ &
          matrix_ptr%xpulsematelyy*facs(2)+ &
          matrix_ptr%xpulsematelzz*facs(3)
     call arbitraryconfig_matel00transpose(tempmatel,tempconfigmat,numconfig)
     tempconfigmat(:,:)=TRANSPOSE(tempconfigmat(:,:))
     do ir=1,numr
        if (velflag.eq.0) then
           diagmat(:,:,ir)=diagmat(:,:,ir)+tempconfigmat(:,:)*bondpoints(ir)
        else
           diagmat(:,:,ir)=diagmat(:,:,ir)+tempconfigmat(:,:)/bondpoints(ir)
        endif
     enddo

     !! A-SQUARED TERM...  looked ok for z polarization.  number of electrons times a-squared.  fac of 2 needed apparently.
     !!   for length... need nuclear dipole too.

     !!   for velocity ...  also need derivative operator in bond length for hetero !!!!  TO DO !!!!


        gg=0.25d0


     if (velflag.ne.0) then 
        csum0 = matrix_ptr%kefac * numelec * (facs(1)**2 + facs(2)**2 + facs(3)**2) * 2  * gg   !! a-squared term.
     else
        !! xpulsenuc is nuclear dipole, divided by R.    from nucdipvalue.

        csum0 = (matrix_ptr%xpulsenuc(1) * facs(1) + matrix_ptr%xpulsenuc(2) * facs(2) + matrix_ptr%xpulsenuc(3) * facs(3) )
     endif
     
     do ir=1,numr
        if (velflag.eq.0) then
           csum=bondpoints(ir) *csum0
        else
           csum=   csum0    !! A-squared !!    no R factor !!
        endif
        do i=1,numconfig
           diagmat(i,i,ir)=diagmat(i,i,ir) + csum 
        enddo
     enddo
  endif   !! pulse

  do ir=1,numr
     call dfrestrictmatrix(diagmat(:,:,ir))
  enddo

  bigconfigmat(:,:,:,:)=0d0

  do ir=1,numr
     bigconfigmat(ir,:,ir,:)=diagmat(:,:,ir)
  enddo

  if (nonuc_checkflag.eq.0.and.nucflag.eq.1) then
       call arbitraryconfig_matel00transpose(matrix_ptr%xymatel,tempconfigmat,numconfig)   !! fills in all, wasteful
       tempconfigmat(:,:)=TRANSPOSE(tempconfigmat(:,:))
       call dfrestrictmatrix(tempconfigmat)

     do ir=1,numr
        do jr=1,numr
           bigconfigmat(ir,:,jr,:)=bigconfigmat(ir,:,jr,:)+tempconfigmat(:,:)*proderivmod(ir,jr)
        enddo
     enddo
     do i=1,numconfig
        bigconfigmat(:,i,:,i)=bigconfigmat(:,i,:,i) + rkemod(:,:)*matrix_ptr%kefac
     enddo
  endif

  deallocate(tempconfigmat,tempconfigmat2,diagmat)

end subroutine assemble_configmat


!! NO A-SQUARED TERM HERE; TIME IS NOT AN ARGUMENT.  A-squared in sparsesparsemult_nompi.

subroutine assemble_sparsemats(matrix_ptr, sparse_ptr,boflag, nucflag, pulseflag, conflag)
  use parameters
  use configptrmod
  use sparseptrmod
  use opmod   !! rkemod, proderivmod
  implicit none
  integer :: conflag,boflag,nucflag,pulseflag
  Type(CONFIGPTR) :: matrix_ptr
  Type(SPARSEPTR) :: sparse_ptr

  if (sparseconfigflag.eq.0) then
     OFLWR "BADDDCAL555LL"; CFLST
  endif

  if (boflag==1) then
     call arbitraryconfig_matel_doubles00transpose(matrix_ptr%xtwoematel(:,:,:,:),sparse_ptr%xpotsparsemattr(:,:),maxdoublewalks)
     call arbitraryconfig_matel00transpose(matrix_ptr%xpotmatel(:,:),sparse_ptr%xonepotsparsemattr(:,:),maxsinglewalks)
     call arbitraryconfig_matel00transpose(matrix_ptr%xopmatel(:,:),sparse_ptr%xopsparsemattr(:,:),maxsinglewalks)
  endif

  if (conflag==1.and.constraintflag.ne.0) then
     call arbitraryconfig_matel00transpose(matrix_ptr%xconmatel,sparse_ptr%xconsparsemattr,maxsinglewalks)
     if (tdflag.ne.0) then
        call arbitraryconfig_matel00transpose(matrix_ptr%xconmatelxx,sparse_ptr%xconsparsemattrxx,maxsinglewalks)
        call arbitraryconfig_matel00transpose(matrix_ptr%xconmatelyy,sparse_ptr%xconsparsemattryy,maxsinglewalks)
        call arbitraryconfig_matel00transpose(matrix_ptr%xconmatelzz,sparse_ptr%xconsparsemattrzz,maxsinglewalks)
     endif
  endif

  if ((pulseflag.eq.1.and.tdflag.eq.1)) then
     call arbitraryconfig_matel00transpose(matrix_ptr%xpulsematelxx,sparse_ptr%xpulsesparsemattrxx,maxsinglewalks)
     call arbitraryconfig_matel00transpose(matrix_ptr%xpulsematelyy,sparse_ptr%xpulsesparsemattryy,maxsinglewalks)
     call arbitraryconfig_matel00transpose(matrix_ptr%xpulsematelzz,sparse_ptr%xpulsesparsemattrzz,maxsinglewalks)
  endif
  if (nonuc_checkflag.eq.0.and.nucflag.eq.1) then
       call arbitraryconfig_matel00transpose(matrix_ptr%xymatel,sparse_ptr%xysparsemattr,maxsinglewalks)
  endif

end subroutine assemble_sparsemats

