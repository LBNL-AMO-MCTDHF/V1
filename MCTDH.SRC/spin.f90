
!! SUBROUTINES FOR SPIN EIGENFUNCTION PROJECTOR
!!   SEE SPINWALKS.F90 ALSO

#include "Definitions.INC"


subroutine configspin_projectall(vector, iprint)
  use spinwalkmod
  use parameters
  implicit none
  integer :: iprint
  DATATYPE :: vector(numconfig,numr)
  integer,save :: allochere=0
  DATATYPE,save,allocatable :: smallvect(:,:), smalltemp(:,:),  configvector1(:,:), vectortr(:,:)

  if (allochere.eq.0) then
     allocate(smallvect(maxspinsetsize,numr), smalltemp(maxspinsetsize,numr), &
       configvector1(numconfig,numr), vectortr(numr,numconfig))
  endif
  allochere=1

  call configspin_project0(vector,iprint,1,numconfig,smallvect,smalltemp,configvector1,vectortr)

end subroutine configspin_projectall



subroutine configspin_projectmine(vector, iprint)
  use spinwalkmod
  use parameters
  implicit none
  integer :: iprint
  DATATYPE :: vector(botwalk:topwalk,numr)
  integer,save :: allochere=0
  DATATYPE,save,allocatable :: smallvect(:,:), smalltemp(:,:),  configvector1(:,:), vectortr(:,:)

  if (allochere.eq.0) then
     allocate(smallvect(maxspinsetsize,numr), smalltemp(maxspinsetsize,numr), &
       configvector1(botwalk:topwalk,numr), vectortr(numr,botwalk:topwalk))
  endif
  allochere=1

  call configspin_project0(vector,iprint,botwalk,topwalk,smallvect,smalltemp,configvector1,vectortr)

end subroutine configspin_projectmine


subroutine configspin_project0(vector,iprint,ibot,itop,smallvect,smalltemp,configvector1,vectortr)
  use spinwalkmod
  use parameters
  implicit none
  integer :: iprint, iset, ii, ibot,itop,isize
  real*8 :: normsq, normsq2
  DATATYPE :: hermdot        , vector(ibot:itop,numr)
  !! keep hermdot, want to see if wfn has gotten chopped; dot could increase with cmctdh

  DATATYPE :: smallvect(maxspinsetsize,numr), smalltemp(maxspinsetsize,numr), &
       configvector1(ibot:itop,numr), vectortr(numr,ibot:itop)

  if (.not.((ibot.eq.1.and.itop.eq.numconfig).or.(ibot.eq.botwalk.and.itop.eq.topwalk))) then
     OFLWR "ISIZE UNRECOGNIIIIZED"; CFLST
  endif

  isize=itop-ibot+1

  if (spinwalkflag==0) then
     OFLWR "Error, configspin_projectmany called but spinwalkflag is 1"; CFLST
  endif
  normsq=real(hermdot(vector,vector,numr*isize))  !! ok hermdot    !! just for check, so don't allreduce.

  configvector1(:,:) = 0.d0
  do iset=1,numspinsets
     do ii=1,spinsetsize(iset)
        smallvect(ii,:)=vector(spinsets(ii,iset),:)
     enddo
     call MYGEMM('N', 'N', spinsetsize(iset), numr, spinsetsize(iset), DATAONE, spinsetprojector(iset)%mat, spinsetsize(iset), smallvect, maxspinsetsize,DATAZERO,smalltemp, maxspinsetsize)
     do ii=1,spinsetsize(iset)
        configvector1(spinsets(ii,iset),1:numr) = configvector1(spinsets(ii,iset),1:numr) + smalltemp(ii,:)
     enddo
  enddo

  if (isize.eq.numconfig.and.sparseconfigflag.ne.0) then
     
     vectortr=TRANSPOSE(vector)
     call mpiallgather(vectortr,numconfig*numr, configsperproc*numr,maxconfigsperproc*numr)
     vector=TRANSPOSE(vectortr)
     
  endif
     
  normsq2=real(hermdot(vector,vector,numr*numconfig))  !! ok hermdot

  if (iprint/=0) then
     if (abs(normsq/normsq2-1.d0).gt.1.d-7) then
        OFLWR "Warning, in configspin_projectmany I lost norm: ", normsq, normsq2; CFL
     endif
  endif
end subroutine configspin_project0



subroutine configspin_transformto(howmany,invector,outvector)
  use spinwalkmod
  use parameters
  implicit none
  integer :: howmany, iset, iind,ii
  DATATYPE :: invector(numconfig,howmany),outvector(spintotrank,howmany), smallvect(maxspinsetsize,howmany), smalltemp(maxspinsetsize,howmany), &
       outvectortr(howmany,spintotrank)

  outvector(:,:)=0d0

  iind=spinstart
  do iset=1,numspinsets
     smallvect(:,:)=0d0
     do ii=1,spinsetsize(iset)
        smallvect(ii,:)=invector(spinsets(ii,iset),:)
     enddo
     call MYGEMM('T', 'N', spinsetrank(iset), howmany, spinsetsize(iset), DATAONE, spinsetprojector(iset)%vects, spinsetsize(iset), smallvect, maxspinsetsize,DATAZERO,smalltemp, maxspinsetsize)
     outvector(iind:iind+spinsetrank(iset)-1,:) = smalltemp(1:spinsetrank(iset),:)
     iind=iind+spinsetrank(iset)
  enddo

  if (iind.ne.spinend+1) then
     OFLWR "IIND ERRO", iind-spinstart,spinrank,iind,spinend; CFLST
  endif

  if (sparseconfigflag.ne.0) then
     outvectortr=TRANSPOSE(outvector)
     call mpiallgather(outvectortr,spintotrank*howmany, allspinranks*howmany,maxspinrank*howmany)
     outvector=TRANSPOSE(outvectortr)
  endif

end subroutine configspin_transformto



subroutine configspin_transformto_mine(howmany,invector,outvector)
  use spinwalkmod
  use parameters
  implicit none
  integer :: howmany, iset, iind,ii
  DATATYPE :: invector(botwalk:topwalk,howmany),outvector(spinrank,howmany), smallvect(maxspinsetsize,howmany), smalltemp(maxspinsetsize,howmany)


  outvector(:,:)=0d0

  iind=1
  do iset=1,numspinsets
     smallvect(:,:)=0d0
     do ii=1,spinsetsize(iset)
        smallvect(ii,:)=invector(spinsets(ii,iset),:)
     enddo
     call MYGEMM('T', 'N', spinsetrank(iset), howmany, spinsetsize(iset), DATAONE, spinsetprojector(iset)%vects, spinsetsize(iset), smallvect, maxspinsetsize,DATAZERO,smalltemp, maxspinsetsize)
     outvector(iind:iind+spinsetrank(iset)-1,:) = smalltemp(1:spinsetrank(iset),:)
     iind=iind+spinsetrank(iset)
  enddo

  if (iind.ne.spinrank+1) then
     OFLWR "IIND ERRO", iind,spinrank; CFLST
  endif

end subroutine configspin_transformto_mine




subroutine configspin_transformto_mine_transpose(invectortr,outvectortr)
  use spinwalkmod
  use parameters
  implicit none
  integer :: iset, iind,ii
  DATATYPE :: invectortr(numr,botwalk:topwalk),outvectortr(numr,spinrank)
  integer, save :: allochere=0
  DATATYPE,save,allocatable :: smallvecttr(:,:), smalltemptr(:,:)

  if (allochere.eq.0) then
     allocate(smallvecttr(numr,maxspinsetsize), smalltemptr(numr,maxspinsetsize))
  endif
  allochere=1

  outvectortr(:,:)=0d0

  iind=1
  do iset=1,numspinsets
     smallvecttr(:,:)=0d0
     do ii=1,spinsetsize(iset)
        smallvecttr(:,ii)=invectortr(:,spinsets(ii,iset))
     enddo
!     call MYGEMM('T', 'N', spinsetrank(iset), numr, spinsetsize(iset), DATAONE, spinsetprojector(iset)%vects, spinsetsize(iset), smallvect, maxspinsetsize,DATAZERO,smalltemp, maxspinsetsize)

     call MYGEMM('N', 'N', numr, spinsetrank(iset), spinsetsize(iset), DATAONE, smallvecttr, numr, spinsetprojector(iset)%vects, spinsetsize(iset), DATAZERO,smalltemptr, numr)

     outvectortr(:,iind:iind+spinsetrank(iset)-1) = smalltemptr(:,1:spinsetrank(iset))
     iind=iind+spinsetrank(iset)
  enddo

  if (iind.ne.spinrank+1) then
     OFLWR "IIND ERRO", iind,spinrank; CFLST
  endif

end subroutine configspin_transformto_mine_transpose



subroutine configspin_transformfrom(howmany,invector,outvector)
  use spinwalkmod
  use parameters
  implicit none
  integer :: howmany, iset, iind,ii
  DATATYPE :: outvector(numconfig,howmany),invector(spintotrank,howmany), smallvect(maxspinsetsize,howmany), smalltemp(maxspinsetsize,howmany), &
       outvectortr(howmany,numconfig)

  outvector(:,:)=0d0

  iind=spinstart
  do iset=1,numspinsets
     smallvect(:,:)=0d0
     smallvect(1:spinsetrank(iset),:)=invector(iind:iind+spinsetrank(iset)-1,:)
     call MYGEMM('N', 'N', spinsetsize(iset), howmany, spinsetrank(iset), DATAONE, spinsetprojector(iset)%vects, spinsetsize(iset), smallvect, maxspinsetsize,DATAZERO,smalltemp, maxspinsetsize)

     do ii=1,spinsetsize(iset)
        outvector(spinsets(ii,iset),:)=smalltemp(ii,:)
     enddo
     iind=iind+spinsetrank(iset)
  enddo

  if (iind.ne.spinend+1) then
     OFLWR "IIND ERRO", iind-spinstart,spinrank,iind,spinend; CFLST
  endif

  if (sparseconfigflag.ne.0) then
     outvectortr=TRANSPOSE(outvector)
     call mpiallgather(outvectortr,numconfig*howmany, configsperproc*howmany,maxconfigsperproc*howmany)
     outvector=TRANSPOSE(outvectortr)
  endif

end subroutine configspin_transformfrom



subroutine configspin_transformfrom_mine(howmany,invector,outvector)
  use spinwalkmod
  use parameters
  implicit none
  integer :: howmany, iset, iind,ii
  DATATYPE :: outvector(botwalk:topwalk,howmany),invector(spinrank,howmany), smallvect(maxspinsetsize,howmany), smalltemp(maxspinsetsize,howmany)

  outvector(:,:)=0d0

  iind=1
  do iset=1,numspinsets
     smallvect(:,:)=0d0
     smallvect(1:spinsetrank(iset),:)=invector(iind:iind+spinsetrank(iset)-1,:)
     call MYGEMM('N', 'N', spinsetsize(iset), howmany, spinsetrank(iset), DATAONE, spinsetprojector(iset)%vects, spinsetsize(iset), smallvect, maxspinsetsize,DATAZERO,smalltemp, maxspinsetsize)

     do ii=1,spinsetsize(iset)
        outvector(spinsets(ii,iset),:)=smalltemp(ii,:)
     enddo
     iind=iind+spinsetrank(iset)
  enddo

  if (iind.ne.spinrank+1) then
     OFLWR "IIND ERROxx", iind, spinrank; CFLST
  endif

end subroutine configspin_transformfrom_mine




subroutine configspin_transformfrom_mine_transpose(invectortr,outvectortr)
  use spinwalkmod
  use parameters
  implicit none
  integer :: iset, iind,ii
  DATATYPE :: outvectortr(numr,botwalk:topwalk),invectortr(numr,spinrank)
  integer, save :: allochere=0
  DATATYPE,save,allocatable :: smallvecttr(:,:), smalltemptr(:,:)

  if (allochere.eq.0) then
     allocate(smallvecttr(numr,maxspinsetsize), smalltemptr(numr,maxspinsetsize))
  endif
  allochere=1

  outvectortr(:,:)=0d0

  iind=1
  do iset=1,numspinsets
     smallvecttr(:,:)=0d0
     smallvecttr(:,1:spinsetrank(iset))=invectortr(:,iind:iind+spinsetrank(iset)-1)

!!     call MYGEMM('N', 'N', spinsetsize(iset), numr, spinsetrank(iset), DATAONE, spinsetprojector(iset)%vects, spinsetsize(iset), smallvect, maxspinsetsize,DATAZERO,smalltemp, maxspinsetsize)

     call MYGEMM('N', 'T', numr, spinsetsize(iset), spinsetrank(iset), DATAONE, smallvecttr, numr, spinsetprojector(iset)%vects, spinsetsize(iset), DATAZERO,smalltemptr, numr)

     do ii=1,spinsetsize(iset)
        outvectortr(:,spinsets(ii,iset))=smalltemptr(:,ii)
     enddo
     iind=iind+spinsetrank(iset)
  enddo

  if (iind.ne.spinrank+1) then
     OFLWR "IIND ERROxx", iind, spinrank; CFLST
  endif

end subroutine configspin_transformfrom_mine_transpose






subroutine configspin_matel()   
  use spinwalkmod
  use parameters
  implicit none
  integer ::     config2, config1,   iwalk, myind,myiostat
  if (spinwalkflag == 0) then
     return
  endif

  if (walksonfile.ne.0) then

     read(751,iostat=myiostat) configspinmatel

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
        write(751) configspinmatel
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
!  DATATYPE :: doublevects(maxspinsetsize**2)
!  real*8 :: doublevals(maxspinsetsize)
  
  OFLWR "Getting SPARSE spinset projector.  Numspinsets is ", numspinsets, " maxspinsetsize is ", maxspinsetsize; CFL


  if (walksonfile.ne.0) then

     OFLWR "Reading spin set projectors...";CFL

     read(751,iostat=myiostat) numspinsets,maxspinsetsize,spinrank,spindfrank
     if (myiostat.ne.0) then
        OFLWR "Read error xx44 savewalks.BIN"; CFLST
     endif

     allocate(spinsets(maxspinsetsize,numspinsets),spinsetprojector(numspinsets))
     spinsetsize=0;spinsets=0; spinsetrank=0

     read(751,iostat=myiostat) spinsetsize(1:numspinsets),spinsetrank(1:numspinsets),spinsets(:,:)
     if (myiostat.ne.0) then
        OFLWR "Read error xx45 savewalks.BIN"; CFLST
     endif

     do iset=1,numspinsets
        allocate(spinsetprojector(iset)%mat(spinsetsize(iset), spinsetsize(iset)))
        allocate(spinsetprojector(iset)%vects(spinsetsize(iset), spinsetrank(iset)))
     enddo

     read(751,iostat=myiostat) (spinsetprojector(i)%mat,spinsetprojector(i)%vects,i=1,numspinsets)     
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
           WRFL  "INFO SSYEV", info
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
        write(751) numspinsets,maxspinsetsize,spinrank,spindfrank
        write(751) spinsetsize(1:numspinsets),spinsetrank(1:numspinsets),spinsets(:,1:numspinsets)
        write(751) (spinsetprojector(i)%mat,spinsetprojector(i)%vects,i=1,numspinsets)
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

  maxspinrank=spinrank
  allocate(allspinstart(nprocs)); allspinstart=(-1); allocate(allspinranks(nprocs)); allspinranks=(-1)




  if (sparseconfigflag.eq.0) then
     spintotrank=spinrank;      spinstart=1; spinend=spintotrank
  else
     allspinranks(myrank)=spinrank
     ii=0
     do i=1,nprocs
        allspinstart(i)=ii+1
        call mympiibcast(allspinranks(i),i,1)
        ii=ii+allspinranks(i)
        if (allspinranks(i).gt.maxspinrank) then
           maxspinrank=allspinranks(i)
        endif
     enddo
     spintotrank=ii
     spinstart=allspinstart(myrank)
     spinend=spinstart+spinrank-1
  endif

  spintotdfrank=spindfrank; call mympiireduceone(spintotdfrank)
  

  OFLWR "TOTAL (all processors) rank",spintotrank," out of ", numconfig
  WRFL "     Spin DF rank is ", spindfrank, "   all processors ", spintotdfrank
  WRFL "   This proc:  spinstart : ", spinstart
  WRFL; CFL
  

  
end subroutine configspinset_projector







