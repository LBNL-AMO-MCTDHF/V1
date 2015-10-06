
!! SUBROUTINES FOR SPIN EIGENFUNCTION PROJECTOR
!!   SEE SPINWALKS.F90 ALSO

#include "Definitions.INC"


subroutine configspin_project(vector, iprint)
  use parameters
  implicit none
  integer :: iprint
  DATATYPE :: vector(numr,firstconfig:lastconfig)

  if (parconsplit.eq.0) then
     call configspin_project_all(vector(:,:),iprint)
  else
     call configspin_project_local(vector(:,:),iprint)
  endif

end subroutine configspin_project


subroutine configspin_project_all(vector, iprint)
  use spinwalkmod
  use parameters
  implicit none
  integer :: iprint
  DATATYPE :: vector(numr,numconfig)

  call configspin_project_local(vector(:,botwalk),iprint)

  if (sparseconfigflag.ne.0) then
     call mpiallgather(vector,numconfig*numr, configsperproc*numr,maxconfigsperproc*numr)
  endif

end subroutine configspin_project_all



subroutine configspin_project_local(vector,iprint)
  use spinwalkmod
  use parameters
  implicit none
  integer :: iprint, iset, ii, isize
  real*8 :: normsq, normsq2
  DATATYPE :: hermdot, vector(numr,botwalk:topwalk)
  DATATYPE :: smallvect(numr,maxspinsetsize), smalltemp(numr,maxspinsetsize), &
       outvector(numr,botwalk:topwalk)

  if (spinwalkflag==0) then
     OFLWR "Error, configspin_projectmany called but spinwalkflag is 1"; CFLST
  endif

  isize=topwalk-botwalk+1

  normsq=real(hermdot(vector,vector,numr*isize))  !! ok hermdot
  if (sparseconfigflag.ne.0) then
     call mympirealreduceone(normsq)
  endif

  outvector(:,:) = 0.d0
  do iset=1,numspinsets
     do ii=1,spinsetsize(iset)
        smallvect(:,ii)=vector(:,spinsets(ii,iset))
     enddo
     call MYGEMM('N', 'T', numr, spinsetsize(iset), spinsetsize(iset), DATAONE,  smallvect, numr, spinsetprojector(iset)%mat, spinsetsize(iset),DATAZERO,smalltemp, numr)
     do ii=1,spinsetsize(iset)
        outvector(:,spinsets(ii,iset)) = outvector(:,spinsets(ii,iset)) + smalltemp(:,ii)
     enddo
  enddo

  vector(:,:)=outvector(:,:)

  normsq2=real(hermdot(vector,vector,numr*isize))  !! ok hermdot
  if (sparseconfigflag.ne.0) then
     call mympirealreduceone(normsq2)
  endif

  if (iprint/=0) then
     if (abs(normsq/normsq2-1.d0).gt.1.d-7) then
        OFLWR "Warning, in configspin_project_local I lost norm: ", normsq, normsq2; CFL
     endif
  endif

end subroutine configspin_project_local


subroutine configspin_transformto(nblock,invector,outvector)
  use spinwalkmod
  use parameters
  implicit none
  integer,intent(in) :: nblock
  DATATYPE,intent(in) :: invector(nblock,localnconfig)
  DATATYPE,intent(out) :: outvector(nblock,localnspin)
  if (parconsplit.eq.0) then
     call configspin_transformto_all(nblock,invector,outvector)
  else
     call configspin_transformto_local(nblock,invector,outvector)
  endif
end subroutine configspin_transformto


!! ALL NUMCONFIG.

subroutine configspin_transformto_all(nblock,invector,outvector)
  use spinwalkmod
  use parameters
  implicit none
  integer :: nblock
  DATATYPE,intent(in) :: invector(nblock,numconfig)
  DATATYPE,intent(out) :: outvector(nblock,numspinconfig)

  outvector(:,:)=0d0
  call configspin_transformto_local(nblock,invector(:,botwalk),outvector(:,spinstart))
  if (sparseconfigflag.ne.0) then
     call mpiallgather(outvector(:,:),numspinconfig*nblock, spinsperproc(:)*nblock,maxspinsperproc*nblock)
  endif

end subroutine configspin_transformto_all




subroutine configspin_transformto_local(nblock,invector,outvector)
  use spinwalkmod
  use parameters
  implicit none
  integer :: nblock, iset, iind,ii
  DATATYPE,intent(in) :: invector(nblock,botwalk:topwalk)
  DATATYPE,intent(out) :: outvector(nblock,spinrank)

  DATATYPE :: smallvect(nblock,maxspinsetsize), smalltemp(nblock,maxspinsetsize)

  outvector(:,:)=0d0

  iind=1
  do iset=1,numspinsets
     smallvect(:,:)=0d0
     do ii=1,spinsetsize(iset)
        smallvect(:,ii)=invector(:,spinsets(ii,iset))
     enddo
     call MYGEMM('N', 'N', nblock, spinsetrank(iset), spinsetsize(iset), DATAONE, smallvect,nblock, spinsetprojector(iset)%vects, spinsetsize(iset), DATAZERO,smalltemp, nblock)
     outvector(:,iind:iind+spinsetrank(iset)-1) = smalltemp(:,1:spinsetrank(iset))
     iind=iind+spinsetrank(iset)
  enddo

  if (iind.ne.spinrank+1) then
     OFLWR "IIND ERRO", iind,spinrank; CFLST
  endif

end subroutine configspin_transformto_local



subroutine configspin_transformfrom(nblock,invector,outvector)
  use parameters
  implicit none
  integer :: nblock
  DATATYPE,intent(out) :: outvector(nblock,localnconfig)
  DATATYPE,intent(in) :: invector(nblock,localnspin)
  if (parconsplit.eq.0) then
     call configspin_transformfrom_all(nblock,invector,outvector)
  else
     call configspin_transformfrom_local(nblock,invector,outvector)
  endif
end subroutine configspin_transformfrom


subroutine configspin_transformfrom_all(nblock,invector,outvector)
  use spinwalkmod
  use parameters
  implicit none
  integer :: nblock
  DATATYPE,intent(out) :: outvector(nblock,numconfig)
  DATATYPE,intent(in) :: invector(nblock,numspinconfig)

  outvector(:,:)=0d0
  call configspin_transformfrom_local(nblock,invector(:,spinstart),outvector(:,botwalk))
  if (sparseconfigflag.ne.0) then
        call mpiallgather(outvector(:,:),numconfig*nblock, configsperproc(:)*nblock,maxconfigsperproc*nblock)
  endif

end subroutine configspin_transformfrom_all



subroutine configspin_transformfrom_local(nblock,invector,outvector)
  use spinwalkmod
  use parameters
  implicit none
  integer :: nblock,iset, iind, ii
  DATATYPE,intent(out) :: outvector(nblock,botwalk:topwalk)
  DATATYPE,intent(in) :: invector(nblock,spinrank)
  DATATYPE :: smallvect(nblock,maxspinsetsize), smalltemp(nblock,maxspinsetsize)

  outvector(:,:)=0d0

  iind=1
  do iset=1,numspinsets
     smallvect(:,1:spinsetrank(iset))=invector(:,iind:iind+spinsetrank(iset)-1)
     call MYGEMM('N', 'T', nblock, spinsetsize(iset), spinsetrank(iset), DATAONE, smallvect, nblock, spinsetprojector(iset)%vects, spinsetsize(iset), DATAZERO,smalltemp, nblock)
     do ii=1,spinsetsize(iset)
        outvector(:,spinsets(ii,iset))=smalltemp(:,ii)
     enddo
     iind=iind+spinsetrank(iset)
  enddo

  if (iind.ne.spinrank+1) then
     OFLWR "IIND ERROxx", iind, spinrank; CFLST
  endif

end subroutine configspin_transformfrom_local





subroutine configspin_matel()   
  use spinwalkmod
  use parameters
  implicit none
  integer ::     config2, config1,   iwalk, myind,myiostat
  if (spinwalkflag == 0) then
     return
  endif

  if (walksonfile.ne.0) then

     if (walksinturn) then 
        call beforebarrier()
     endif
     
     OFLWR "   ...reading configspinmatel..."; CFL
     read(751,iostat=myiostat) configspinmatel
     OFLWR "   ...ok, done reading configspinmatel..."; CFL

     if (walksinturn) then 
        call afterbarrier()
     endif

     call mympiimax(myiostat)
     if (myiostat.ne.0) then
        OFLWR "Read error for savewalks.BIN!  Delete it to recompute walks. 662", myiostat;CFLST
     endif

  else   !WALKSONFILE

     configspinmatel(:,:)=0.d0
     do config1=botwalk,topwalk
        myind=1
        
!! msvalue is 2x ms quantum number

        configspinmatel(myind,config1) = msvalue(config1)**2/4.d0 + numunpaired(config1)/2.d0
        
        do iwalk=1,numspinwalks(config1)
           config2=spinwalk(iwalk,config1)
           myind=1+iwalk
           configspinmatel(myind,config1) = configspinmatel(myind,config1) + &
                spinwalkdirphase(iwalk,config1)
        enddo
     enddo
     if (walkwriteflag.ne.0) then
        if (walksinturn) then
           call beforebarrier()
        endif

        OFLWR "   ...writing configspinmatel..."; CFL
        write(751) configspinmatel
        OFLWR "   ...ok, done writing configspinmatel..."; CFL

        if (walksinturn) then
           call afterbarrier()
        endif
     endif
  endif

end subroutine configspin_matel


function spinallowed(spinval)
  use parameters
  implicit none
  real*8 :: spinval
  logical :: spinallowed
  if (abs(spinval-(spinrestrictval/2.d0*(spinrestrictval/2.d0+1))).lt.1.d-3) then
     spinallowed=.true.
  else
     spinallowed=.false.
  endif
end function

  
subroutine configspinset_projector()   
  use spinwalkmod
  use configmod   !! configlist for spintotdfrank
  use mpimod
  use parameters
  implicit none
  integer :: info, lwork,j,i,ii,iset,jj, elim, elimsets, flag, iwalk,myiostat,spindfrank
  real*8, allocatable :: spinvects(:,:), spinvals(:), work(:), realprojector(:,:)
  logical :: spinallowed,dfincluded
  integer :: allspinstart(nprocs)
!  DATATYPE :: doublevects(maxspinsetsize**2)
!  real*8 :: doublevals(maxspinsetsize)
  
  OFLWR "Getting SPARSE spinset projector.  Numspinsets is ", numspinsets, " maxspinsetsize is ", maxspinsetsize; CFL


  if (walksonfile.ne.0) then

     if (walksinturn) then
        call beforebarrier()
     endif

     OFLWR "Reading spin set projectors...";CFL
     read(751,iostat=myiostat) numspinsets,maxspinsetsize,spinrank,spindfrank

     if (walksinturn) then
        call afterbarrier()
     endif

     call mympiimax(myiostat)
     if (myiostat.ne.0) then
        OFLWR "Read error xx44 savewalks.BIN"; CFLST
     endif

     allocate(spinsets(maxspinsetsize,numspinsets),spinsetprojector(numspinsets))
     spinsetsize=0;spinsets=0; spinsetrank=0

     if (walksinturn) then
        call beforebarrier()
     endif
     
     OFLWR "   ...reading spinsets..."; CFL
     read(751,iostat=myiostat) spinsetsize(1:numspinsets),spinsetrank(1:numspinsets),spinsets(:,:)

     if (walksinturn) then
        call afterbarrier()
     endif

     call mympiimax(myiostat)
     if (myiostat.ne.0) then
        OFLWR "Read error xx45 savewalks.BIN"; CFLST
     endif

     do iset=1,numspinsets
        allocate(spinsetprojector(iset)%mat(spinsetsize(iset), spinsetsize(iset)))
        allocate(spinsetprojector(iset)%vects(spinsetsize(iset), spinsetrank(iset)))
     enddo

     if (walksinturn) then
        call beforebarrier()
     endif

     OFLWR "   ...reading spin set projector ... "; CFL
     read(751,iostat=myiostat) (spinsetprojector(i)%mat,spinsetprojector(i)%vects,i=1,numspinsets)     
     OFLWR "   ...ok, done reading spin set projector ..."; CFL

     if (walksinturn) then
        call afterbarrier()
     endif

     call mympiimax(myiostat)
     if (myiostat.ne.0) then
        OFLWR "Read error setprojector savewalks.BIN"; CFLST
     endif

  else

     allocate(spinsetprojector(numspinsets))

     allocate(spinvects(maxspinsetsize,maxspinsetsize), spinvals(maxspinsetsize), &
          realprojector(maxspinsetsize,maxspinsetsize))
     lwork=10*maxspinsetsize;  allocate(work(lwork))
     
     OFLWR "Getting spin set projectors...";CFL

     elim=0;  elimsets=0;  iset=1; spinrank=0; spindfrank=0
  
     do while (iset.le.numspinsets)
        spinvects=0.d0
        do ii=1,spinsetsize(iset)
           do jj=1,spinsetsize(iset)
              
              if (ii.eq.jj) then
                 spinvects(ii,jj)=configspinmatel(1, spinsets(jj,iset))
              else
                 flag=0
                 do iwalk=1,numspinwalks(spinsets(jj,iset))
                    if (spinwalk(iwalk,spinsets(jj,iset)).eq.spinsets(ii,iset)) then
                       spinvects(ii,jj)=configspinmatel(iwalk+1, spinsets(jj,iset))
                       flag=1
                       exit
                    endif
                 enddo
              endif
           enddo
        enddo
     
        call dsyev('V','U', spinsetsize(iset), spinvects, maxspinsetsize, spinvals, work, lwork, info)
        if (info/=0) then
           OFLWR  "INFO SSYEV", info; CFLST
        endif
        j=0; 
        do i=1,spinsetsize(iset)
           if (spinallowed(spinvals(i))) then
              j=j+1;           spinvects(:,j)=spinvects(:,i)
           endif
        enddo
        spinsetrank(iset)=j
        spinrank=spinrank+j

        if (dfrestrictflag.eq.0) then
           spindfrank=spindfrank+j
        else
           if (dfincluded(configlist(:,spinsets(1,iset)))) then
              spindfrank=spindfrank+j
           endif
        endif

        spinvects(:,j+1:maxspinsetsize)=0d0
        
        if (spinsetrank(iset)==0) then 
           elimsets=elimsets+1
           elim=elim+spinsetsize(iset)
           spinsetsize(iset:numspinsets-1)=spinsetsize(iset+1:numspinsets)
           spinsets(:,iset:numspinsets-1)=spinsets(:,iset+1:numspinsets)
           numspinsets=numspinsets-1
        else
           allocate(spinsetprojector(iset)%mat(spinsetsize(iset), spinsetsize(iset)))
           allocate(spinsetprojector(iset)%vects(spinsetsize(iset), spinsetrank(iset)))
           
           spinsetprojector(iset)%vects(:,:)=spinvects(1:spinsetsize(iset), 1:spinsetrank(iset))
           
           call dgemm('N', 'T', spinsetsize(iset), spinsetsize(iset),spinsetrank(iset),1.0d0, spinvects, maxspinsetsize, &
                spinvects, maxspinsetsize ,0.0d0,realprojector, maxspinsetsize)
           
           spinsetprojector(iset)%mat(:,:) = realprojector(1:spinsetsize(iset), 1:spinsetsize(iset))

!just checking right
!        call CONFIGHERM(spinsetprojector(iset)%mat,spinsetsize(iset),spinsetsize(iset), doublevects,doublevals)
!        do i=1,spinsetsize(iset)
!           if (abs(doublevals(i)*2-1.d0)-1.d0 .gt. 1.d-9) then
!              OFLWR "SPIN PROJECTOR ERROR", doublevals(i); CFLST
!           endif
!        enddo

           iset=iset+1

        endif
     enddo

     if (walkwriteflag.ne.0) then
        if (walksinturn) then
           call beforebarrier()
        endif

        OFLWR "   ...writing spinsets..."; CFL
        write(751) numspinsets,maxspinsetsize,spinrank,spindfrank
        write(751) spinsetsize(1:numspinsets),spinsetrank(1:numspinsets),spinsets(:,1:numspinsets)
        write(751) (spinsetprojector(i)%mat,spinsetprojector(i)%vects,i=1,numspinsets)
        OFLWR "   ...ok, done writing spinsets..."; CFL

        if (walksinturn) then
           call afterbarrier()
        endif
     endif

     deallocate(spinvects, spinvals, realprojector, work)
  
     OFLWR "...done.  Eliminated ", elimsets, " sets with total rank ", elim   ; CFL

  endif !! walksonfile


  i=0
  do ii=1,numspinsets
     i=i+spinsetrank(ii)
  enddo
  if (i.ne.spinrank) then
     OFLWR "CHECKMEERRR",i,spinrank; CFLST
  endif

  
  OFLWR  "Number of spin sets is now ", numspinsets
  WRFL "Total # of spinvects with S^2 = ", (spinrestrictval/2.d0*(spinrestrictval/2.d0+1)), " is ", spinrank, " out of ", topwalk-botwalk+1; CFL

  maxspinsperproc=spinrank

  allocate(spinsperproc(nprocs)); spinsperproc=(-1)


  if (sparseconfigflag.eq.0) then
     numspinconfig=spinrank;      spinstart=1; spinend=numspinconfig
     spinsperproc(1)=numspinconfig
     maxspinsperproc=numspinconfig
     firstspinconfig=1; lastspinconfig=numspinconfig; localnspin=numspinconfig
  else
     spinsperproc(myrank)=spinrank
     ii=0
     do i=1,nprocs
        allspinstart(i)=ii+1
        call mympiibcast(spinsperproc(i),i,1)
        ii=ii+spinsperproc(i)
        if (spinsperproc(i).gt.maxspinsperproc) then
           maxspinsperproc=spinsperproc(i)
        endif
     enddo
     numspinconfig=ii
     spinstart=allspinstart(myrank)
     spinend=spinstart+spinrank-1
     if (parconsplit.eq.0) then
        firstspinconfig=1
        lastspinconfig=numspinconfig
        localnspin=numspinconfig

     else
        firstspinconfig=spinstart
        lastspinconfig=spinend
        localnspin=spinrank
     endif
  endif

  spintotdfrank=spindfrank; 
  if (sparseconfigflag.ne.0) then !! 05-08-2015 Okay.. was only used if sparseconfigflag.ne.0.. 
     call mympiireduceone(spintotdfrank)
  endif
  

  OFLWR "TOTAL (all processors) rank",numspinconfig," out of ", numconfig
  WRFL "     Spin DF rank is ", spindfrank, "   all processors ", spintotdfrank
  WRFL "   This proc:  spinstart : ", spinstart
  WRFL; CFL
  

  
end subroutine configspinset_projector







