
!! ALL MODULES

!! subroutines for sparse matrix multiplication of Slater determinant coefficient vectors
!! the matrices are represented by derived type configptr or sparseptr depending on sparseopt
!! these derived types are handled using derived type walktype - their indexing etc is stored there


#include "Definitions.INC"

!!!!!!      BYPROC SUBROUTINES   !!!!!!!

module arbitrarymultmod
  use mpisubmod

contains

  subroutine arbitrary_sparsemult_singles_byproc(firstproc,lastproc,www,mattrans, &
       rvector,insmallvector,outsmallvector,mynumr,diagflag)
    use walkmod
    use sparse_parameters
    use fileptrmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: firstproc,lastproc
    integer,intent(in) :: mynumr,diagflag
    DATATYPE,intent(in) :: insmallvector(mynumr,www%allbotconfigs(firstproc):www%alltopconfigs(lastproc)), &
         mattrans(www%singlematsize,www%botconfig:www%topconfig)
    DATATYPE,intent(out) :: outsmallvector(mynumr,www%botconfig:www%topconfig)
    DATAECS,intent(in) :: rvector(mynumr)
    DATATYPE :: outsum(mynumr)
    DATAECS :: myrvector(mynumr)
    integer :: config1,ihop

    if (sparseconfigflag.eq.0) then
       OFLWR "error no arbitrary_sparsemult if sparseconfigflag.eq.0"; CFLST
    endif

    if (www%topconfig-www%botconfig+1.eq.0) then
       return
    endif

    outsmallvector(:,:)=0d0

    if (www%alltopconfigs(lastproc).lt.www%allbotconfigs(firstproc)) then
       return
    endif

    if (diagflag.eq.0) then

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(config1,ihop,outsum,myrvector)
       myrvector(:)=rvector(:)
!$OMP DO SCHEDULE(DYNAMIC)
       do config1=www%botconfig,www%topconfig
          outsum(:)=0d0

          do ihop=www%firstsinglehopbyproc(firstproc,config1),www%lastsinglehopbyproc(lastproc,config1)
             outsum(:)=outsum(:)+mattrans(ihop,config1) * &
                  insmallvector(:,www%singlehop(ihop,config1)) * myrvector(:)
          enddo
          outsmallvector(:,config1)=outsum(:)
       enddo
!$OMP END DO
!$OMP END PARALLEL

    else

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(config1,myrvector)
       myrvector(:)=rvector(:)

!$OMP DO SCHEDULE(DYNAMIC)
       do config1=www%botconfig,www%topconfig
          if (www%singlehopdiagflag(config1).ne.0) then
             outsmallvector(:,config1)=mattrans(www%singlediaghop(config1),config1) * &
                  insmallvector(:,config1) * myrvector(:)
          endif
       enddo
!$OMP END DO
!$OMP END PARALLEL

    endif

  end subroutine arbitrary_sparsemult_singles_byproc


  subroutine arbitrary_sparsemult_doubles_byproc(firstproc,lastproc,www,mattrans,&
       rvector,insmallvector,outsmallvector,mynumr,diagflag)
    use fileptrmod
    use sparse_parameters
    use walkmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: firstproc,lastproc,mynumr,diagflag
    DATATYPE,intent(in) :: insmallvector(mynumr,www%allbotconfigs(firstproc):www%alltopconfigs(lastproc)), &
         mattrans(www%doublematsize,www%botconfig:www%topconfig)
    DATATYPE,intent(out) :: outsmallvector(mynumr,www%botconfig:www%topconfig)
    DATAECS,intent(in) :: rvector(mynumr)
    DATATYPE :: outsum(mynumr)
    DATAECS :: myrvector(mynumr)
    integer :: config1,ihop

    if (sparseconfigflag.eq.0) then
       OFLWR "error no arbitrary_sparsemult_doubles if sparseconfigflag.eq.0"; CFLST
    endif

    if (www%topconfig-www%botconfig+1.eq.0) then
       return
    endif

    outsmallvector(:,:)=0d0

    if (www%alltopconfigs(lastproc).lt.www%allbotconfigs(firstproc)) then
       return
    endif

    if (diagflag.eq.0) then

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(config1,ihop,outsum,myrvector)
       myrvector(:)=rvector(:)

!$OMP DO SCHEDULE(DYNAMIC)
       do config1=www%botconfig,www%topconfig
          outsum(:)=0d0
          do ihop=www%firstdoublehopbyproc(firstproc,config1),www%lastdoublehopbyproc(lastproc,config1)
             outsum(:)=outsum(:)+mattrans(ihop,config1) * &
                  insmallvector(:,www%doublehop(ihop,config1)) * myrvector(:)
          enddo
          outsmallvector(:,config1)=outsum(:)
       enddo
!$OMP END DO
!$OMP END PARALLEL

    else

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(config1,myrvector)
       myrvector(:)=rvector(:)

!$OMP DO SCHEDULE(DYNAMIC)
       do config1=www%botconfig,www%topconfig
          if (www%doublehopdiagflag(config1).ne.0) then
             outsmallvector(:,config1)=mattrans(www%doublediaghop(config1),config1) * &
                  insmallvector(:,config1) * myrvector(:)
          endif
       enddo
!$OMP END DO
!$OMP END PARALLEL

    endif

  end subroutine arbitrary_sparsemult_doubles_byproc



  subroutine arbitraryconfig_mult_singles_byproc_hhh(firstproc,lastproc,www,in_onebodymat, &
       rvector, avectorin, avectorout,inrnum,diagflag,inholeflag)
    use walkmod
    use mpimod   !! myrank
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: firstproc,lastproc,inholeflag,inrnum,diagflag
    DATATYPE,intent(in) :: in_onebodymat(www%nspf,www%nspf), &
         avectorin(inrnum,www%allbotconfigs(firstproc):www%alltopconfigs(lastproc))
    DATATYPE,intent(out) :: avectorout(inrnum,www%botconfig:www%topconfig)
    DATAECS,intent(in) :: rvector(inrnum)
    DATATYPE :: onebodymat(www%nspf,www%nspf)    !! AUTOMATIC
    DATATYPE :: csum, outsum(inrnum), holetrace
    DATAECS :: myrvector(inrnum)
    integer ::    config2, config1, iwalk, idiag,ihop,ispf

    if (www%topconfig-www%botconfig+1.eq.0) then
       return
    endif

    avectorout(:,:)=0d0

    if (www%alltopconfigs(lastproc).lt.www%allbotconfigs(firstproc)) then
       return
    endif

    onebodymat(:,:)=in_onebodymat(:,:)

!! HOLEFLAG

    if (www%holeflag.ne.0.and.inholeflag.eq.1) then
       holetrace=0
       do ispf=1,www%nspf
          holetrace=holetrace + onebodymat(ispf,ispf) * 2
       enddo

       onebodymat(:,:) = (-1) * in_onebodymat(:,:)

       if (myrank.ge.firstproc.and.myrank.le.lastproc) then
          do config1=www%botconfig,www%topconfig
             avectorout(:,config1) = holetrace*avectorin(:,config1)*rvector(:)
          enddo
       endif
    endif    !! holeflag

    if (diagflag.eq.0) then

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(config1,ihop,iwalk,config2,csum,outsum,myrvector)
       myrvector(:)=rvector(:)

!$OMP DO SCHEDULE(DYNAMIC)
       do config1=www%botconfig,www%topconfig
          outsum(:)=avectorout(:,config1)
          do ihop=www%firstsinglehopbyproc(firstproc,config1),www%lastsinglehopbyproc(lastproc,config1)
             config2=www%singlehop(ihop,config1)
             csum=0d0
             do iwalk=www%singlehopwalkstart(ihop,config1),www%singlehopwalkend(ihop,config1)
                csum=csum +  onebodymat(www%singlewalkopspf(1,iwalk+www%scol(config1)), &
                     www%singlewalkopspf(2,iwalk+www%scol(config1))) *  &
                     www%singlewalkdirphase(iwalk+www%scol(config1))
             enddo
             outsum(:)=outsum(:)+avectorin(:,config2) * csum * myrvector(:)
          enddo
          avectorout(:,config1)=outsum(:)

       enddo
!$OMP END DO
!$OMP END PARALLEL

    else

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(config1,idiag,iwalk,csum,myrvector)
       myrvector(:)=rvector(:)

!$OMP DO SCHEDULE(DYNAMIC)
       do config1=www%botconfig,www%topconfig
          csum=0d0
          do idiag=1,www%numsinglediagwalks(config1)
             iwalk=www%singlediag(idiag,config1)
             csum=csum+&
                  onebodymat(www%singlewalkopspf(1,iwalk+www%scol(config1)), &
                  www%singlewalkopspf(2,iwalk+www%scol(config1))) *  &
                  www%singlewalkdirphase(iwalk+www%scol(config1))
          enddo
          avectorout(:,config1)=avectorout(:,config1) + &
               avectorin(:,config1) * csum * myrvector(:)
       enddo
!$OMP END DO
!$OMP END PARALLEL

    endif

  end subroutine arbitraryconfig_mult_singles_byproc_hhh


  subroutine arbitraryconfig_mult_singles_byproc(firstproc,lastproc,www,onebodymat, &
       rvector, avectorin, avectorout,inrnum,diagflag)
    use walkmod
    use mpimod   !! myrank
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: firstproc,lastproc,inrnum,diagflag
    DATATYPE,intent(in) :: onebodymat(www%nspf,www%nspf), &
         avectorin(inrnum,www%allbotconfigs(firstproc):www%alltopconfigs(lastproc))
    DATATYPE,intent(out) :: avectorout(inrnum,www%botconfig:www%topconfig)
    DATAECS,intent(in) :: rvector(inrnum)

    call arbitraryconfig_mult_singles_byproc_hhh(firstproc,lastproc,www,onebodymat, &
         rvector, avectorin, avectorout,inrnum,diagflag,1)

  end subroutine arbitraryconfig_mult_singles_byproc


  subroutine arbitraryconfig_mult_doubles_byproc(firstproc,lastproc,www,in_twobodymat, &
       rvector, avectorin, avectorout,inrnum,diagflag)
    use walkmod
    use mpimod   !! myrank
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: firstproc,lastproc,inrnum,diagflag
    DATATYPE,intent(in) :: in_twobodymat(www%nspf,www%nspf,www%nspf,www%nspf), &
         avectorin(inrnum,www%allbotconfigs(firstproc):www%alltopconfigs(lastproc))
    DATATYPE,intent(out) :: avectorout(inrnum,www%botconfig:www%topconfig)
    DATAECS,intent(in) :: rvector(inrnum)
    DATATYPE, allocatable :: twobodymat(:,:,:,:)
    DATATYPE :: csum, outsum(inrnum),holetrace
    DATAECS :: myrvector(inrnum)
    integer ::   config2, config1, iwalk, idiag, ihop,ispf,jspf

    if (www%topconfig-www%botconfig+1.eq.0) then
       return
    endif

    avectorout(:,:)=0.d0

    if (www%alltopconfigs(lastproc).lt.www%allbotconfigs(firstproc)) then
       return
    endif

    allocate(twobodymat(www%nspf,www%nspf,www%nspf,www%nspf))
    twobodymat(:,:,:,:) = in_twobodymat(:,:,:,:)

    if (www%holeflag.ne.0) then

       holetrace=0
       do ispf=1,www%nspf
          do jspf=1,www%nspf
             holetrace=holetrace + twobodymat(ispf,ispf,jspf,jspf) * 2 &
                  - twobodymat(ispf,jspf,jspf,ispf)
          enddo
       enddo

!!       OFLWR "HTRACE002",holetrace,www%nspf; CFL

!! NO!  Two-electron operator       twobodymat(:,:,:,:) = (-1) * in_twobodymat(:,:,:,:)

       if (myrank.ge.firstproc.and.myrank.le.lastproc) then
          do config1=www%botconfig,www%topconfig
             avectorout(:,config1) = holetrace*avectorin(:,config1) * rvector(:)
          enddo
       endif
    endif    !! holeflag

    if (diagflag.eq.0) then

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(config1,ihop,iwalk,config2,csum,outsum,myrvector)
       myrvector(:)=rvector(:)
!$OMP DO SCHEDULE(DYNAMIC)
       do config1=www%botconfig,www%topconfig
          outsum(:)=0d0
          do ihop=www%firstdoublehopbyproc(firstproc,config1),www%lastdoublehopbyproc(lastproc,config1)
             config2=www%doublehop(ihop,config1)
             csum=0d0
             do iwalk=www%doublehopwalkstart(ihop,config1),www%doublehopwalkend(ihop,config1)
                csum=csum+&
                     twobodymat(www%doublewalkdirspf(1,iwalk+www%dcol(config1)), &
                     www%doublewalkdirspf(2,iwalk+www%dcol(config1)),   &
                     www%doublewalkdirspf(3,iwalk+www%dcol(config1)),   &
                     www%doublewalkdirspf(4,iwalk+www%dcol(config1)))* &
                     www%doublewalkdirphase(iwalk+www%dcol(config1))
             enddo
             outsum(:)=outsum(:)+avectorin(:,config2) *  csum * myrvector(:)
          enddo
          avectorout(:,config1)=avectorout(:,config1)+outsum(:)
       enddo
!$OMP END DO
!$OMP END PARALLEL

    else

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(config1,idiag,iwalk,csum,myrvector)
       myrvector(:)=rvector(:)
!$OMP DO SCHEDULE(DYNAMIC)
       do config1=www%botconfig,www%topconfig
          csum=0d0
          do idiag=1,www%numdoublediagwalks(config1)
             iwalk=www%doublediag(idiag,config1)
             csum=csum+&
                  twobodymat(www%doublewalkdirspf(1,iwalk+www%dcol(config1)), &
                  www%doublewalkdirspf(2,iwalk+www%dcol(config1)),   &
                  www%doublewalkdirspf(3,iwalk+www%dcol(config1)),   &
                  www%doublewalkdirspf(4,iwalk+www%dcol(config1)))* &
                  www%doublewalkdirphase(iwalk+www%dcol(config1))
          enddo
          avectorout(:,config1)=avectorout(:,config1) + &
               avectorin(:,config1) * csum * myrvector(:)
       enddo
!$OMP END DO
!$OMP END PARALLEL

    endif

    deallocate(twobodymat)

  end subroutine arbitraryconfig_mult_doubles_byproc



  subroutine arbitraryconfig_mult_singles_noparcon(www,onebodymat, rvector, avectorin, avectorout,inrnum)   
    use walkmod
    use mpimod   !! myrank,nprocs
    use fileptrmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: inrnum
    DATATYPE,intent(in) :: onebodymat(www%nspf,www%nspf), avectorin(inrnum,www%firstconfig:www%lastconfig)
    DATATYPE,intent(out) :: avectorout(inrnum,www%firstconfig:www%lastconfig)
    DATAECS,intent(in) :: rvector(inrnum)

    if (www%parconsplit.ne.0) then
       OFLWR "ERROR NOPARCON ARB"; CFLST
    endif   !! now botconfig:topconfig = firstconfig:lastconfig
    
    call arbitraryconfig_mult_singles_byproc(1,nprocs,www,onebodymat, &
         rvector, avectorin,avectorout(:,www%botconfig:www%topconfig),inrnum,0)
  
    call mpiallgather(avectorout,www%numconfig*inrnum,www%configsperproc(:)*inrnum,&
         www%maxconfigsperproc*inrnum)

  end subroutine arbitraryconfig_mult_singles_noparcon

#ifdef MPIFLAG

  subroutine arbitraryconfig_mult_singles_gather(www,onebodymat, rvector, avectorin, avectorout,inrnum)   
    use walkmod
    use mpimod   !! myrank,nprocs
    use fileptrmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: inrnum
    DATATYPE,intent(in) :: onebodymat(www%nspf,www%nspf), avectorin(inrnum,www%firstconfig:www%lastconfig)
    DATATYPE,intent(out) :: avectorout(inrnum,www%firstconfig:www%lastconfig)
    DATAECS,intent(in) :: rvector(inrnum)
    DATATYPE, allocatable :: workvector(:,:)

    if (www%parconsplit.eq.0) then
       OFLWR "ERROR PARCON ARB"; CFLST
    endif   !! now botconfig:topconfig = firstconfig:lastconfig

    allocate(workvector(inrnum,www%numconfig))
    workvector(:,:)=0d0

    if (www%topconfig.ge.www%botconfig) then
       workvector(:,www%botconfig:www%topconfig) = avectorin(:,:)
    endif
    call mpiallgather(workvector,www%numconfig*inrnum,www%configsperproc(:)*inrnum,&
         www%maxconfigsperproc*inrnum)

    if (www%lastconfig.ge.www%firstconfig) then
       call arbitraryconfig_mult_singles_byproc(1,nprocs,www,onebodymat,&
            rvector, workvector, avectorout,inrnum,0)
    endif
    deallocate(workvector)

  end subroutine arbitraryconfig_mult_singles_gather

  subroutine arbitraryconfig_mult_singles_summa(www,onebodymat, rvector, avectorin, avectorout,inrnum)   
    use walkmod
    use mpimod   !! myrank,nprocs
    use fileptrmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: inrnum
    DATATYPE,intent(in) :: onebodymat(www%nspf,www%nspf), &
         avectorin(inrnum,www%firstconfig:www%lastconfig)
    DATATYPE,intent(out) :: avectorout(inrnum,www%firstconfig:www%lastconfig)
    DATAECS,intent(in) :: rvector(inrnum)
    DATATYPE :: workvector(inrnum,www%maxconfigsperproc),&
         outtemp(inrnum,www%botconfig:www%topconfig+1)         !! AUTOMATIC
    integer :: iproc

    if (www%parconsplit.eq.0) then
       OFLWR "ERROR PARCON ARB"; CFLST
    endif   !! now botconfig:topconfig = firstconfig:lastconfig

    workvector=0; outtemp=0

    if (www%lastconfig.ge.www%firstconfig) then
       avectorout(:,:)=0d0
    endif

    do iproc=1,nprocs

       if (myrank.eq.iproc.and.www%lastconfig.ge.www%firstconfig) then
          workvector(:,1:www%topconfig-www%botconfig+1) = avectorin(:,:)
       endif
       if (www%alltopconfigs(iproc).ge.www%allbotconfigs(iproc)) then
          call mympibcast(workvector,iproc,&
               (www%alltopconfigs(iproc)-www%allbotconfigs(iproc)+1)*inrnum)
          if (www%lastconfig.ge.www%firstconfig) then
             call arbitraryconfig_mult_singles_byproc(iproc,iproc,www,onebodymat, &
                  rvector, workvector, outtemp,inrnum,0)           
             avectorout(:,:) = avectorout(:,:) + outtemp(:,www%botconfig:www%topconfig)
          endif
       endif
    enddo

  end subroutine arbitraryconfig_mult_singles_summa

  subroutine arbitraryconfig_mult_singles_circ(www,onebodymat, rvector, avectorin, avectorout,inrnum)   
    use walkmod
    use mpimod   !! myrank,nprocs
    use fileptrmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: inrnum
    DATATYPE,intent(in) :: onebodymat(www%nspf,www%nspf),&
         avectorin(inrnum,www%firstconfig:www%lastconfig)
    DATATYPE,intent(out) :: avectorout(inrnum,www%firstconfig:www%lastconfig)
    DATAECS,intent(in) :: rvector(inrnum)
    DATATYPE :: workvector(inrnum,www%maxconfigsperproc),workvector2(inrnum,www%maxconfigsperproc),&
         outtemp(inrnum,www%botconfig:www%topconfig+1)         !! AUTOMATIC
    integer :: iproc,prevproc,nextproc,deltaproc

    if (www%parconsplit.eq.0) then
       OFLWR "ERROR PARCON ARB"; CFLST
    endif   !! now botconfig:topconfig = firstconfig:lastconfig

!! doing circ mult slightly different than e.g. SINCDVR/coreproject.f90 
!!     and ftcore.f90, holding hands in a circle, prevproc and nextproc, 
!!     each chunk gets passed around the circle
    prevproc=mod(nprocs+myrank-2,nprocs)+1
    nextproc=mod(myrank,nprocs)+1

    workvector=0; workvector2=0; outtemp=0

    if (www%topconfig.ge.www%botconfig) then
       avectorout=0d0
       workvector(:,1:www%topconfig-www%botconfig+1) = avectorin(:,:)
    endif

    do deltaproc=0,nprocs-1

!! PASSING BACKWARD (plus deltaproc)
       iproc=mod(myrank-1+deltaproc,nprocs)+1

       if (www%topconfig.ge.www%botconfig.and.www%alltopconfigs(iproc).ge.www%allbotconfigs(iproc)) then
          call arbitraryconfig_mult_singles_byproc(iproc,iproc,www,onebodymat, &
               rvector, workvector, outtemp,inrnum,0)
          avectorout(:,:) = avectorout(:,:) + outtemp(:,www%botconfig:www%topconfig)
       endif

!! PASSING BACKWARD
!! mympisendrecv(sendbuf,recvbuf,dest,source,...)

       call mympisendrecv(workvector,workvector2,prevproc,nextproc,deltaproc,&
            inrnum * www%maxconfigsperproc)
       workvector(:,:)=workvector2(:,:)

    enddo

  end subroutine arbitraryconfig_mult_singles_circ

#endif

  subroutine arbitraryconfig_mult_singles(www,onebodymat, rvector, avectorin, avectorout,inrnum)   
    use walkmod
    use sparse_parameters
    use fileptrmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: inrnum
    DATATYPE,intent(in) :: onebodymat(www%nspf,www%nspf), &
         avectorin(inrnum,www%firstconfig:www%lastconfig)
    DATATYPE,intent(out) :: avectorout(inrnum,www%firstconfig:www%lastconfig)
    DATAECS,intent(in) :: rvector(inrnum)
    
#ifdef MPIFLAG
    if (www%parconsplit.eq.0) then
#endif
       call arbitraryconfig_mult_singles_noparcon(www,onebodymat, rvector, avectorin, avectorout,inrnum)
#ifdef MPIFLAG
    else
       select case (sparsesummaflag)
       case (0)
          call arbitraryconfig_mult_singles_gather(www,onebodymat, rvector, avectorin, avectorout,inrnum)
       case (1)
          call arbitraryconfig_mult_singles_summa(www,onebodymat, rvector, avectorin, avectorout,inrnum)
       case (2)
          call arbitraryconfig_mult_singles_circ(www,onebodymat, rvector, avectorin, avectorout,inrnum)
       case default
          OFLWR "Error sparsesummaflag ", sparsesummaflag; CFLST
       end select
    endif
#endif

  end subroutine arbitraryconfig_mult_singles

end module arbitrarymultmod


module sparsemultmod
use arbitrarymultmod
contains

!! DIRECT MATVEC (with matrix_ptr not sparse_ptr)

  subroutine direct_sparseconfigmult_byproc(firstproc,lastproc,www,invector,outvector,matrix_ptr, &
       boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag,imc)
    use ham_parameters
    use r_parameters
    use configptrmod
    use opmod   !! rkemod, proderivmod  
    use mpimod  !! myrank
    use walkmod
    use pulsesubmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: firstproc,lastproc,imc
    integer,intent(in) :: conflag,boflag,nucflag,pulseflag,onlytdflag,botr,topr,diagflag
    DATATYPE,intent(in) :: invector(botr:topr,www%allbotconfigs(firstproc):www%alltopconfigs(lastproc))
    DATATYPE,intent(out) :: outvector(botr:topr,www%botconfig:www%topconfig)
    Type(CONFIGPTR),intent(in) :: matrix_ptr
    real*8, intent(in) :: time
    integer :: mynumr,ir,flag,ii
    DATATYPE :: tempvector(botr:topr,www%botconfig:www%topconfig) !! AUTOMATIC
    DATATYPE :: tempmatel(www%nspf,www%nspf),facs(3),csum
    DATAECS :: rvector(botr:topr)
    real*8 :: gg

    if (www%topconfig.lt.www%botconfig) then
       return
    endif

    outvector(:,:)=0d0

    if (www%alltopconfigs(lastproc).lt.www%allbotconfigs(firstproc)) then
       return
    endif

    if (onlytdflag.ne.0) then
       facs(:)=0;     facs(onlytdflag)=1d0
    else
       call vectdpot(time,velflag,facs,imc)
    endif
  
    mynumr=topr-botr+1
  
    if (boflag==1) then
     
!! easy nuclear repulsion, hardwired like this for now

       if (myrank.ge.firstproc.and.myrank.le.lastproc) then
          do ir=botr,topr
             outvector(ir,:)=outvector(ir,:) + invector(ir,www%botconfig:www%topconfig) * ( &
                  energyshift + (nucrepulsion+frozenpotdiag)/bondpoints(ir) + &
                  frozenkediag/bondpoints(ir)**2 ) * matrix_ptr%constfac
          enddo
       endif

       rvector(:)=1/bondpoints(botr:topr)
       call arbitraryconfig_mult_doubles_byproc(firstproc,lastproc,www,matrix_ptr%xtwoematel(:,:,:,:),&
            rvector,invector,tempvector,mynumr,diagflag)
       outvector(:,:)=outvector(:,:)+tempvector(:,:)

       if (www%holeflag.ne.0) then
          call arbitraryconfig_mult_singles_byproc_hhh(firstproc,lastproc,www,matrix_ptr%xtwoe_htrace(:,:),&
               rvector,invector,tempvector,mynumr,diagflag,0)
          outvector(:,:)=outvector(:,:) + tempvector(:,:)
       endif

       call arbitraryconfig_mult_singles_byproc(firstproc,lastproc,www,matrix_ptr%xpotmatel(:,:),&
            rvector,invector,tempvector,mynumr,diagflag)
       outvector(:,:)=outvector(:,:)+tempvector(:,:)

       rvector(:)=1/bondpoints(botr:topr)**2
       call arbitraryconfig_mult_singles_byproc(firstproc,lastproc,www,matrix_ptr%xopmatel(:,:),&
            rvector,invector,tempvector,mynumr,diagflag)
       outvector(:,:)=outvector(:,:)+tempvector(:,:)
     
    endif
  
    if (conflag==1.and.constraintflag.ne.0) then
       rvector(:)=1
     
       tempmatel(:,:)=matrix_ptr%xconmatel(:,:)
      if (tdflag.ne.0) then
          tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelxx(:,:)*facs(1)
          tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelyy(:,:)*facs(2)
          tempmatel(:,:)=tempmatel(:,:)+matrix_ptr%xconmatelzz(:,:)*facs(3)
       endif
       call arbitraryconfig_mult_singles_byproc(firstproc,lastproc,www,tempmatel,&
            rvector,invector,tempvector,mynumr,diagflag)
       outvector(:,:)=outvector(:,:)-tempvector(:,:)

    endif
    if (pulseflag.eq.1.and.tdflag.ne.0.or.onlytdflag.ne.0) then
       if (velflag.eq.0) then
          rvector(:)=bondpoints(botr:topr)
       else
          rvector(:)=1/bondpoints(botr:topr)
       endif

       tempmatel(:,:)=0d0; flag=0
       if (onlytdflag.eq.0.or.onlytdflag.eq.1) then
          tempmatel(:,:)= tempmatel(:,:) + matrix_ptr%xpulsematelxx*facs(1); flag=1
       endif
       if (onlytdflag.eq.0.or.onlytdflag.eq.2) then
          tempmatel(:,:)= tempmatel(:,:) + matrix_ptr%xpulsematelyy*facs(2); flag=1
       endif
       if (onlytdflag.eq.0.or.onlytdflag.eq.3) then
          tempmatel(:,:)= tempmatel(:,:) + matrix_ptr%xpulsematelzz*facs(3); flag=1
       endif
       if (flag.eq.1) then
          call arbitraryconfig_mult_singles_byproc(firstproc,lastproc,www,tempmatel,&
               rvector,invector,tempvector,mynumr,diagflag)
          outvector(:,:)=outvector(:,:)+tempvector(:,:)
       endif

!! STILL NEED A-VECTOR TERM IN BOND LENGTH FOR HETERO!! TO DO!!

       gg=0.25d0

       if (myrank.ge.firstproc.and.myrank.le.lastproc) then

!!  A-squared term and nuclear dipole:  nuclear dipole always
          if (velflag.eq.0) then 
             csum=0
             do ii=1,3
                if (onlytdflag.eq.0.or.onlytdflag.eq.ii) then
                   csum=csum+facs(ii)*matrix_ptr%xpulsenuc(ii)
                endif
             enddo
             do ir=botr,topr
                outvector(ir,:)=outvector(ir,:)+invector(ir,www%botconfig:www%topconfig) * &
                     csum*bondpoints(ir)
             enddo
          else if (onlytdflag.eq.0) then
!! a-squared
             csum=matrix_ptr%constfac * www%numelec * (facs(1)**2 + facs(2)**2 + facs(3)**2 ) * 2  * gg
             do ir=botr,topr
!! NO R FACTOR !!
                outvector(ir,:)=outvector(ir,:)+invector(ir,www%botconfig:www%topconfig)*csum
             enddo
          endif
       endif
    endif   !! PULSE

    if (nonuc_checkflag.eq.0.and.nucflag.eq.1.and.mynumr.eq.numr) then
       if (diagflag.eq.0) then
          rvector(:)=1d0
          call arbitraryconfig_mult_singles_byproc(firstproc,lastproc,www,matrix_ptr%xymatel,&
               rvector,invector,tempvector,numr,0)
          call MYGEMM('N','N',numr,www%topconfig-www%botconfig+1,numr,DATAONE, &
               proderivmod,numr,     tempvector,           numr,DATAONE,outvector,numr)
          if (myrank.ge.firstproc.and.myrank.le.lastproc) then
             call MYGEMM('N','N',numr,www%topconfig-www%botconfig+1,numr,DATAONE*matrix_ptr%kefac,&
                  rkemod,numr,     invector(:,www%botconfig),numr,DATAONE,outvector,numr)
          endif
!!$     else
!!$        do ir=1,numr
!!$           outvector(ir,:)=outvector(ir,:)+&
!!$                invector(ir,www%botconfig:www%topconfig)*matrix_ptr%kefac*rkemod(ir,ir)
!!$        enddo
       endif
    endif

  end subroutine direct_sparseconfigmult_byproc

!! SPARSE MATVEC (with sparse_ptr not matrix_ptr)

  subroutine sparsesparsemult_byproc(firstproc,lastproc,www,invector,outvector,sparse_ptr,&
       boflag,nucflag,pulseflag,conflag,time,onlytdflag,botr,topr,diagflag,imc)
    use sparse_parameters
    use ham_parameters
    use r_parameters
    use sparseptrmod
    use walkmod
    use opmod   !! rkemod, proderivmod
    use mpimod  !! myrank 
    use fileptrmod
    use pulsesubmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: firstproc,lastproc,imc
    integer,intent(in) :: boflag,nucflag,pulseflag,conflag,onlytdflag,botr,topr,diagflag
    DATATYPE,intent(in)  :: invector(botr:topr,www%allbotconfigs(firstproc):www%alltopconfigs(lastproc))
    DATATYPE,intent(out) :: outvector(botr:topr,www%botconfig:www%topconfig)
    Type(SPARSEPTR),intent(in) :: sparse_ptr
    real*8,intent(in) :: time
    DATATYPE :: facs(3),csum
    DATATYPE :: tempvector(botr:topr,www%botconfig:www%topconfig), &
         tempsparsemattr(www%singlematsize,www%botconfig:www%topconfig) !! AUTOMATIC
    real*8 :: gg
    DATAECS :: rvector(botr:topr)
    integer :: ir,mynumr,ii,flag

    if (sparseconfigflag.eq.0) then
       OFLWR "error no sparsesparsemult if sparseconfigflag.eq.0"; CFLST
    endif

    if (www%topconfig.lt.www%botconfig) then
       return
    endif

    if (botr.ne.1.or.topr.ne.numr) then
       OFLWR "REPROGRAMMMM."; CFLST
    endif

    outvector(:,:)=0d0

    if (www%alltopconfigs(lastproc).lt.www%allbotconfigs(firstproc)) then
       return
    endif

    if (onlytdflag.ne.0) then
       facs(:)=0;     facs(onlytdflag)=1d0
    else
       call vectdpot(time,velflag,facs,imc)
    endif

    mynumr=topr-botr+1

    if (boflag==1) then

       if (myrank.ge.firstproc.and.myrank.le.lastproc) then
          do ir=botr,topr
             outvector(ir,:)=outvector(ir,:) + invector(ir,www%botconfig:www%topconfig) * ( &
                  energyshift + (nucrepulsion+frozenpotdiag)/bondpoints(ir) + &
                  frozenkediag/bondpoints(ir)**2 ) * sparse_ptr%constfac
          enddo
       endif

       rvector(:)=1/bondpoints(botr:topr)
       call arbitrary_sparsemult_doubles_byproc(firstproc,lastproc,www,sparse_ptr%xpotsparsemattr,&
            rvector,invector,tempvector,mynumr,diagflag)
       outvector(:,:)=outvector(:,:)+tempvector(:,:)

       if (www%holeflag.ne.0) then
          call arbitrary_sparsemult_singles_byproc(firstproc,lastproc,www,sparse_ptr%xpottracemattr,&
               rvector,invector,tempvector,mynumr,diagflag)
          outvector(:,:)=outvector(:,:) + tempvector(:,:)
       endif

       call arbitrary_sparsemult_singles_byproc(firstproc,lastproc,www,sparse_ptr%xonepotsparsemattr,&
            rvector,invector,tempvector,mynumr,diagflag)
       outvector(:,:)=outvector(:,:)+tempvector(:,:)

       rvector(:)=1/bondpoints(botr:topr)**2
       call arbitrary_sparsemult_singles_byproc(firstproc,lastproc,www,sparse_ptr%xopsparsemattr,&
            rvector,invector,tempvector,mynumr,diagflag)
       outvector(:,:)=outvector(:,:)+tempvector(:,:)

    endif

    if (conflag==1.and.constraintflag.ne.0) then

       rvector(:)=1

       tempsparsemattr(:,:)=sparse_ptr%xconsparsemattr(:,:)
       if (tdflag.ne.0) then
          tempsparsemattr(:,:)=tempsparsemattr(:,:)+sparse_ptr%xconsparsemattrxx(:,:)*facs(1)
          tempsparsemattr(:,:)=tempsparsemattr(:,:)+sparse_ptr%xconsparsemattryy(:,:)*facs(2)
          tempsparsemattr(:,:)=tempsparsemattr(:,:)+sparse_ptr%xconsparsemattrzz(:,:)*facs(3)
       endif
       call arbitrary_sparsemult_singles_byproc(firstproc,lastproc,www,tempsparsemattr,&
            rvector,invector,tempvector,mynumr,diagflag)
       outvector(:,:)=outvector(:,:)-tempvector(:,:)

    endif

    if (pulseflag.eq.1.and.tdflag.ne.0.or.onlytdflag.ne.0) then
       if (velflag.eq.0) then
          rvector(:)=bondpoints(botr:topr)
       else
          rvector(:)=1/bondpoints(botr:topr)
       endif

       tempsparsemattr(:,:)=0d0; flag=0
       if (onlytdflag.eq.0.or.onlytdflag.eq.1) then
          tempsparsemattr= tempsparsemattr + sparse_ptr%xpulsesparsemattrxx*facs(1); flag=1
       endif
       if (onlytdflag.eq.0.or.onlytdflag.eq.2) then
          tempsparsemattr= tempsparsemattr + sparse_ptr%xpulsesparsemattryy*facs(2); flag=1
       endif
       if (onlytdflag.eq.0.or.onlytdflag.eq.3) then
          tempsparsemattr= tempsparsemattr + sparse_ptr%xpulsesparsemattrzz*facs(3); flag=1
       endif
       if (flag.eq.1) then
          call arbitrary_sparsemult_singles_byproc(firstproc,lastproc,www,tempsparsemattr,&
               rvector,invector,tempvector,mynumr,diagflag)
          outvector(:,:)=outvector(:,:)+tempvector(:,:)
       endif

!! STILL NEED A-VECTOR TERM IN BOND LENGTH FOR HETERO !! TO DO !!

!! A-squared term and nuclear dipole : nuclear dipole always

       gg=0.25d0

       if (myrank.ge.firstproc.and.myrank.le.lastproc) then

          if (velflag.eq.0) then 
             csum=0
             do ii=1,3
                if (onlytdflag.eq.0.or.onlytdflag.eq.ii) then
                   csum=csum+facs(ii)*sparse_ptr%xpulsenuc(ii)
                endif
             enddo
             do ir=botr,topr
                outvector(ir,:)=outvector(ir,:)+invector(ir,www%botconfig:www%topconfig) * &
                     csum*bondpoints(ir)
             enddo
          else if (onlytdflag.eq.0) then
!! a-squared
             csum= sparse_ptr%constfac * www%numelec * (facs(1)**2 + facs(2)**2 + facs(3)**2 ) * 2 * gg
             do ir=botr,topr
!! NO R FACTOR !!
                outvector(ir,:)=outvector(ir,:)+invector(ir,www%botconfig:www%topconfig)*csum
             enddo
          endif
       endif
    endif !! PULSE

    if (nonuc_checkflag.eq.0.and.nucflag.eq.1.and.mynumr.eq.numr) then
       if (diagflag.eq.0) then
          rvector(:)=1d0
          call arbitrary_sparsemult_singles_byproc(firstproc,lastproc,www,sparse_ptr%xysparsemattr,&
               rvector,invector,tempvector,numr,0)
          call MYGEMM('N','N',numr,www%topconfig-www%botconfig+1,numr,DATAONE, &
               proderivmod,numr,     tempvector,           numr,DATAONE,outvector,numr)
          if (myrank.ge.firstproc.and.myrank.le.lastproc) then
             call MYGEMM('N','N',numr,www%topconfig-www%botconfig+1,numr,DATAONE*sparse_ptr%kefac,&
                  rkemod,numr,     invector(:,www%botconfig),numr,DATAONE,outvector,numr)
          endif
!!$     else
!!$        do ir=1,numr
!!$           outvector(ir,:)=outvector(ir,:)+&
!!$             invector(ir,www%botconfig:www%topconfig)*sparse_ptr%kefac*rkemod(ir,ir)
!!$        enddo
       endif
    endif

  end subroutine sparsesparsemult_byproc

  subroutine sparseconfigmult_byproc(firstproc,lastproc,www,invector,outvector,matrix_ptr,sparse_ptr,&
       boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag,imc)
    use sparse_parameters
    use configptrmod
    use sparseptrmod
    use walkmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: firstproc,lastproc
    integer,intent(in) :: conflag,boflag,nucflag,pulseflag,onlytdflag,botr,topr,diagflag,imc
    DATATYPE,intent(in) :: invector(botr:topr,www%allbotconfigs(firstproc):www%alltopconfigs(lastproc))
    DATATYPE,intent(out) :: outvector(botr:topr,www%botconfig:www%topconfig)
    real*8,intent(in) :: time
    Type(CONFIGPTR),intent(in) :: matrix_ptr
    Type(SPARSEPTR),intent(in) :: sparse_ptr
    
    if (www%topconfig-www%botconfig+1.eq.0) then
       return
    endif

    if (sparseopt.eq.0) then
       call direct_sparseconfigmult_byproc(firstproc,lastproc,www,invector,outvector,matrix_ptr, &
            boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag,imc)
    else
       call sparsesparsemult_byproc(firstproc,lastproc,www,invector,outvector,sparse_ptr, &
            boflag, nucflag, pulseflag, conflag,time,onlytdflag,botr,topr,diagflag,imc)
    endif

  end subroutine sparseconfigmult_byproc

!! sparseconfigmultone:
!! For one bond length, born oppenheimer Hamiltonian

  subroutine sparseconfigmultone_noparcon(www,invector,outvector,matrix_ptr,sparse_ptr, &
       boflag,pulseflag,conflag,isplit,time,imc)
    use configptrmod
    use sparseptrmod
    use walkmod
    use mpimod    !! nprocs,myrank
    use fileptrmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: conflag, boflag, pulseflag, isplit,imc
    DATATYPE,intent(in) :: invector(www%firstconfig:www%lastconfig)
    DATATYPE,intent(out) :: outvector(www%firstconfig:www%lastconfig)
    Type(CONFIGPTR),intent(in) :: matrix_ptr
    Type(SPARSEPTR),intent(in) :: sparse_ptr
    real*8,intent(in) :: time

    if (www%parconsplit.ne.0) then
       OFLWR "ERROR NOPARCON", www%parconsplit; CFLST
    endif

    if (www%topconfig.ge.www%botconfig) then
       call sparseconfigmult_byproc(1,nprocs,www,invector,outvector(www%botconfig:www%topconfig),&
            matrix_ptr,sparse_ptr, boflag, 0, pulseflag, conflag,time,0,isplit,isplit,0,imc)
    endif
    call mpiallgather(outvector,www%numconfig,www%configsperproc(:),www%maxconfigsperproc)

  end subroutine sparseconfigmultone_noparcon

#ifdef MPIFLAG

  subroutine sparseconfigmultone_gather(www,invector,outvector,matrix_ptr,sparse_ptr, &
       boflag,pulseflag,conflag,isplit,time,imc)
    use configptrmod
    use sparseptrmod
    use walkmod
    use mpimod    !! nprocs,myrank
    use fileptrmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: conflag, boflag, pulseflag, isplit,imc
    DATATYPE,intent(in) :: invector(www%firstconfig:www%lastconfig)
    DATATYPE,intent(out) :: outvector(www%firstconfig:www%lastconfig)
    Type(CONFIGPTR),intent(in) :: matrix_ptr
    Type(SPARSEPTR),intent(in) :: sparse_ptr
    real*8,intent(in) :: time
    DATATYPE,allocatable :: workvector(:)

    if (www%parconsplit.eq.0) then
       OFLWR "ERROR PARCON"; CFLST
    endif    !! now botconfig:topconfig = firstconfig:lastconfig

    allocate(workvector(www%numconfig))
    workvector(:)=0d0

    if (www%topconfig.ge.www%botconfig) then
       workvector(www%botconfig:www%topconfig)=invector(:)
    endif
    call mpiallgather(workvector,www%numconfig,www%configsperproc(:),www%maxconfigsperproc)
  
    if (www%topconfig.ge.www%botconfig) then
       call sparseconfigmult_byproc(1,nprocs,www,workvector,outvector,&
          matrix_ptr,sparse_ptr, boflag, 0, pulseflag, conflag,time,0,isplit,isplit,0,imc)
    endif
    deallocate(workvector)

  end subroutine sparseconfigmultone_gather

  subroutine sparseconfigmultone_summa(www,invector,outvector,matrix_ptr,sparse_ptr, &
       boflag,pulseflag,conflag,isplit,time,imc)
    use configptrmod
    use sparseptrmod
    use walkmod
    use mpimod    !! nprocs,myrank
    use fileptrmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: conflag, boflag, pulseflag, isplit,imc
    DATATYPE,intent(in) :: invector(www%firstconfig:www%lastconfig)
    DATATYPE,intent(out) :: outvector(www%firstconfig:www%lastconfig)
    Type(CONFIGPTR),intent(in) :: matrix_ptr
    Type(SPARSEPTR),intent(in) :: sparse_ptr
    real*8,intent(in) :: time
    DATATYPE :: workvector(www%maxconfigsperproc),outtemp(www%botconfig:www%topconfig+1)  !! AUTOMATIC
    integer :: iproc

    if (www%parconsplit.eq.0) then
       OFLWR "ERROR PARCON"; CFLST
    endif    !! now botconfig:topconfig = firstconfig:lastconfig

    workvector=0; outtemp=0

    if (www%lastconfig.ge.www%firstconfig) then
       outvector(:)=0d0
    endif

    do iproc=1,nprocs
       if (myrank.eq.iproc.and.www%lastconfig.ge.www%firstconfig) then
          workvector(1:www%topconfig-www%botconfig+1)=invector(:)
       endif
       if (www%alltopconfigs(iproc).ge.www%allbotconfigs(iproc)) then
          call mympibcast(workvector,iproc,(www%alltopconfigs(iproc)-www%allbotconfigs(iproc)+1))
          if (www%lastconfig.ge.www%firstconfig) then
             call sparseconfigmult_byproc(iproc,iproc,www,workvector,outtemp,matrix_ptr,sparse_ptr, &
                  boflag, 0, pulseflag, conflag,time,0,isplit,isplit,0,imc)
             outvector(:)=outvector(:)+outtemp(www%botconfig:www%topconfig)
          endif
       endif
    enddo

  end subroutine sparseconfigmultone_summa


  subroutine sparseconfigmultone_circ(www,invector,outvector,matrix_ptr,sparse_ptr, &
       boflag,pulseflag,conflag,isplit,time,imc)
    use configptrmod
    use sparseptrmod
    use walkmod
    use mpimod    !! nprocs,myrank
    use fileptrmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: conflag, boflag, pulseflag, isplit,imc
    DATATYPE,intent(in) :: invector(www%firstconfig:www%lastconfig)
    DATATYPE,intent(out) :: outvector(www%firstconfig:www%lastconfig)
    Type(CONFIGPTR),intent(in) :: matrix_ptr
    Type(SPARSEPTR),intent(in) :: sparse_ptr
    real*8,intent(in) :: time
    DATATYPE :: workvector(www%maxconfigsperproc), workvector2(www%maxconfigsperproc),&
         outtemp(www%botconfig:www%topconfig+1)      !! AUTOMATIC
    integer :: iproc,prevproc,nextproc,deltaproc

    if (www%parconsplit.eq.0) then
       OFLWR "ERROR PARCON"; CFLST
    endif     !! now botconfig:topconfig = firstconfig:lastconfig

!! doing circ mult slightly different than e.g. SINCDVR/coreproject.f90 
!!     and ftcore.f90, holding hands in a circle, prevproc and nextproc, 
!!     each chunk gets passed around the circle
    prevproc=mod(nprocs+myrank-2,nprocs)+1
    nextproc=mod(myrank,nprocs)+1

    workvector=0; workvector2=0; outtemp=0

    if (www%topconfig.ge.www%botconfig) then
       outvector(:)=0d0
       workvector(1:www%topconfig-www%botconfig+1)=invector(:)
    endif

    do deltaproc=0,nprocs-1

!! PASSING BACKWARD (plus deltaproc)
       iproc=mod(myrank-1+deltaproc,nprocs)+1

       if (www%topconfig.ge.www%botconfig.and.www%alltopconfigs(iproc).ge.www%allbotconfigs(iproc)) then
          call sparseconfigmult_byproc(iproc,iproc,www,workvector,outtemp,matrix_ptr,sparse_ptr, &
               boflag, 0, pulseflag, conflag,time,0,isplit,isplit,0,imc)
          outvector(:)=outvector(:)+outtemp(www%botconfig:www%topconfig)
       endif

!! PASSING BACKWARD
!! mympisendrecv(sendbuf,recvbuf,dest,source,...)

       call mympisendrecv(workvector,workvector2,prevproc,nextproc,deltaproc,www%maxconfigsperproc)
       workvector(:)=workvector2(:)

    enddo

  end subroutine sparseconfigmultone_circ

#endif

  subroutine sparseconfigmultone(www,invector,outvector,matrix_ptr,sparse_ptr, &
       boflag,pulseflag,conflag,isplit,time,imc)
    use sparse_parameters
    use configptrmod
    use sparseptrmod
    use walkmod
    use fileptrmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: conflag, boflag, pulseflag, isplit,imc
    DATATYPE,intent(in) :: invector(www%firstconfig:www%lastconfig)
    DATATYPE,intent(out) :: outvector(www%firstconfig:www%lastconfig)
    Type(CONFIGPTR),intent(in) :: matrix_ptr
    Type(SPARSEPTR),intent(in) :: sparse_ptr
    real*8,intent(in) :: time

#ifdef MPIFLAG
    if (www%parconsplit.eq.0) then
#endif
       call sparseconfigmultone_noparcon(www,invector,outvector,matrix_ptr,sparse_ptr, &
            boflag,pulseflag,conflag,isplit,time,imc)
#ifdef MPIFLAG
    else
       select case(sparsesummaflag)
       case(0)
          call sparseconfigmultone_gather(www,invector,outvector,matrix_ptr,sparse_ptr, &
               boflag,pulseflag,conflag,isplit,time,imc)
       case(1)
          call sparseconfigmultone_summa(www,invector,outvector,matrix_ptr,sparse_ptr, &
               boflag,pulseflag,conflag,isplit,time,imc)
       case(2)
          call sparseconfigmultone_circ(www,invector,outvector,matrix_ptr,sparse_ptr, &
               boflag,pulseflag,conflag,isplit,time,imc)
       case default
          OFLWR "Error sparsesummaflag ",sparsesummaflag; CFLST
       end select
    endif
#endif

  end subroutine sparseconfigmultone

  subroutine sparseconfigmultxxx_noparcon(www,invector,outvector,matrix_ptr,sparse_ptr, &
       boflag, nucflag, pulseflag, conflag,time,onlytdflag,imc)
    use r_parameters
    use configptrmod
    use sparseptrmod
    use mpimod    !! myrank,nprocs
    use walkmod
    use fileptrmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: conflag,boflag,nucflag,pulseflag,onlytdflag,imc
    DATATYPE,intent(in) :: invector(numr,www%firstconfig:www%lastconfig)
    DATATYPE,intent(out) :: outvector(numr,www%firstconfig:www%lastconfig)
    Type(CONFIGPTR),intent(in) :: matrix_ptr
    Type(SPARSEPTR),intent(in) :: sparse_ptr
    real*8,intent(in) :: time

    if (www%parconsplit.ne.0) then
       OFLWR "ERROR NOPARCON XXX", www%parconsplit; CFLST
    endif

    if (www%topconfig.ge.www%botconfig) then
       call sparseconfigmult_byproc(1,nprocs,www,invector,outvector(:,www%botconfig:www%topconfig),&
            matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,onlytdflag,1,numr,0,imc)
    endif
    call mpiallgather(outvector,www%numconfig*numr,www%configsperproc(:)*numr,www%maxconfigsperproc*numr)

  end subroutine sparseconfigmultxxx_noparcon

#ifdef MPIFLAG

  subroutine sparseconfigmultxxx_gather(www,invector,outvector,matrix_ptr,sparse_ptr, &
       boflag, nucflag, pulseflag, conflag,time,onlytdflag,imc)
    use r_parameters
    use configptrmod
    use sparseptrmod
    use mpimod    !! myrank,nprocs
    use walkmod
    use fileptrmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: conflag,boflag,nucflag,pulseflag,onlytdflag,imc
    DATATYPE,intent(in) :: invector(numr,www%firstconfig:www%lastconfig)
    DATATYPE,intent(out) :: outvector(numr,www%firstconfig:www%lastconfig)
    Type(CONFIGPTR),intent(in) :: matrix_ptr
    Type(SPARSEPTR),intent(in) :: sparse_ptr
    real*8,intent(in) :: time
    DATATYPE, allocatable :: workvector(:,:)

    if (www%parconsplit.eq.0) then
       OFLWR "ERROR PARCON XXX"; CFLST
    endif   !! now botconfig:topconfig = firstconfig:lastconfig

    allocate(workvector(numr,www%numconfig))
    workvector(:,:)=0d0

    if (www%topconfig.ge.www%botconfig) then
       workvector(:,www%botconfig:www%topconfig) = invector(:,:)
    endif
    call mpiallgather(workvector,www%numconfig*numr,www%configsperproc(:)*numr,&
         www%maxconfigsperproc*numr)

    if (www%topconfig.ge.www%botconfig) then
       call sparseconfigmult_byproc(1,nprocs,www,workvector,outvector,&
            matrix_ptr,sparse_ptr, boflag, nucflag, pulseflag, conflag,time,&
            onlytdflag,1,numr,0,imc)
    endif
    deallocate(workvector)

  end subroutine sparseconfigmultxxx_gather

  subroutine sparseconfigmultxxx_summa(www,invector,outvector,matrix_ptr,sparse_ptr, &
       boflag, nucflag, pulseflag, conflag,time,onlytdflag,imc)
    use r_parameters
    use configptrmod
    use sparseptrmod
    use mpimod    !! myrank,nprocs
    use walkmod
    use fileptrmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: conflag,boflag,nucflag,pulseflag,onlytdflag,imc
    DATATYPE,intent(in) :: invector(numr,www%firstconfig:www%lastconfig)
    DATATYPE,intent(out) :: outvector(numr,www%firstconfig:www%lastconfig)
    Type(CONFIGPTR),intent(in) :: matrix_ptr
    Type(SPARSEPTR),intent(in) :: sparse_ptr
    real*8,intent(in) :: time
    DATATYPE :: workvector(numr,www%maxconfigsperproc), &
         outtemp(numr,www%botconfig:www%topconfig+1)    !!AUTOMATIC
    integer :: iproc

    if (www%parconsplit.eq.0) then
       OFLWR "ERROR PARCON XXX"; CFLST
    endif   !! now botconfig:topconfig = firstconfig:lastconfig

    workvector=0; outtemp=0

    if (www%lastconfig.ge.www%firstconfig) then
       outvector=0d0
    endif

    do iproc=1,nprocs

       if (myrank.eq.iproc.and.www%lastconfig.ge.www%firstconfig) then
          workvector(:,1:www%topconfig-www%botconfig+1) = invector(:,:)
       endif
       if (www%alltopconfigs(iproc).ge.www%allbotconfigs(iproc)) then
          call mympibcast(workvector,iproc,&
               (www%alltopconfigs(iproc)-www%allbotconfigs(iproc)+1)*numr)
          if (www%lastconfig.ge.www%firstconfig) then
             call sparseconfigmult_byproc(iproc,iproc,www,workvector,outtemp,matrix_ptr,sparse_ptr, &
                  boflag, nucflag, pulseflag, conflag,time,onlytdflag,1,numr,0,imc)
             outvector(:,:)=outvector(:,:) + outtemp(:,www%botconfig:www%topconfig)
          endif
       endif
    enddo

  end subroutine sparseconfigmultxxx_summa

  subroutine sparseconfigmultxxx_circ(www,invector,outvector,matrix_ptr,sparse_ptr, &
       boflag, nucflag, pulseflag, conflag,time,onlytdflag,imc)
    use r_parameters
    use configptrmod
    use sparseptrmod
    use mpimod    !! myrank,nprocs
    use walkmod
    use fileptrmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: conflag,boflag,nucflag,pulseflag,onlytdflag,imc
    DATATYPE,intent(in) :: invector(numr,www%firstconfig:www%lastconfig)
    DATATYPE,intent(out) :: outvector(numr,www%firstconfig:www%lastconfig)
    Type(CONFIGPTR),intent(in) :: matrix_ptr
    Type(SPARSEPTR),intent(in) :: sparse_ptr
    real*8,intent(in) :: time
    DATATYPE :: workvector(numr,www%maxconfigsperproc), workvector2(numr,www%maxconfigsperproc),&
         outtemp(numr,www%botconfig:www%topconfig+1)      !! AUTOMATIC
    integer :: iproc,deltaproc,prevproc,nextproc

    if (www%parconsplit.eq.0) then
       OFLWR "ERROR PARCON XXX"; CFLST
    endif   !! now botconfig:topconfig = firstconfig:lastconfig

!! doing circ mult slightly different than e.g. SINCDVR/coreproject.f90 
!!     and ftcore.f90, holding hands in a circle, prevproc and nextproc, 
!!     each chunk gets passed around the circle
    prevproc=mod(nprocs+myrank-2,nprocs)+1
    nextproc=mod(myrank,nprocs)+1

    workvector=0; workvector2=0; outtemp=0

    if (www%topconfig.ge.www%botconfig) then
       outvector(:,:)=0d0
       workvector(:,1:www%topconfig-www%botconfig+1)=invector(:,:)
    endif

    do deltaproc=0,nprocs-1

!! PASSING BACKWARD (plus deltaproc)
       iproc=mod(myrank-1+deltaproc,nprocs)+1

       if (www%topconfig.ge.www%botconfig.and.www%alltopconfigs(iproc).ge.www%allbotconfigs(iproc)) then
          call sparseconfigmult_byproc(iproc,iproc,www,workvector,outtemp,matrix_ptr,sparse_ptr, &
               boflag, nucflag, pulseflag, conflag,time,onlytdflag,1,numr,0,imc)
          outvector(:,:)=outvector(:,:) + outtemp(:,www%botconfig:www%topconfig)
       endif

!! PASSING BACKWARD
!! mympisendrecv(sendbuf,recvbuf,dest,source,...)

       call mympisendrecv(workvector,workvector2,prevproc,nextproc,deltaproc,&
            numr*www%maxconfigsperproc)
       workvector(:,:)=workvector2(:,:)

    enddo

  end subroutine sparseconfigmultxxx_circ

#endif

  subroutine sparseconfigmultxxx(www,invector,outvector,matrix_ptr,sparse_ptr, &
       boflag, nucflag, pulseflag, conflag,time,onlytdflag,imc)
    use r_parameters
    use configptrmod
    use sparseptrmod
    use sparse_parameters
    use walkmod
    use fileptrmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: conflag,boflag,nucflag,pulseflag,onlytdflag,imc
    DATATYPE,intent(in) :: invector(numr,www%firstconfig:www%lastconfig)
    DATATYPE,intent(out) :: outvector(numr,www%firstconfig:www%lastconfig)
    Type(CONFIGPTR),intent(in) :: matrix_ptr
    Type(SPARSEPTR),intent(in) :: sparse_ptr
    real*8,intent(in) :: time

#ifdef MPIFLAG
    if (www%parconsplit.eq.0) then
#endif
       call sparseconfigmultxxx_noparcon(www,invector,outvector,matrix_ptr,sparse_ptr, &
            boflag, nucflag, pulseflag, conflag,time,onlytdflag,imc)
#ifdef MPIFLAG
    else
       select case (sparsesummaflag)
       case(0)
          call sparseconfigmultxxx_gather(www,invector,outvector,matrix_ptr,sparse_ptr, &
               boflag, nucflag, pulseflag, conflag,time,onlytdflag,imc)
       case(1)
          call sparseconfigmultxxx_summa(www,invector,outvector,matrix_ptr,sparse_ptr, &
               boflag, nucflag, pulseflag, conflag,time,onlytdflag,imc)
       case(2)
          call sparseconfigmultxxx_circ(www,invector,outvector,matrix_ptr,sparse_ptr, &
               boflag, nucflag, pulseflag, conflag,time,onlytdflag,imc)
       case default
          OFLWR "Error sparsesummaflag ",sparsesummaflag; CFLST
       end select
    endif
#endif

  end subroutine sparseconfigmultxxx


!! All purpose subroutine.

  subroutine sparseconfigmult(www,invector,outvector,matrix_ptr,sparse_ptr, &
       boflag, nucflag, pulseflag, conflag,time,imc)
    use r_parameters
    use configptrmod
    use sparseptrmod
    use walkmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: conflag,boflag,nucflag,pulseflag,imc
    DATATYPE,intent(in) :: invector(www%localnconfig*numr)
    DATATYPE,intent(out) :: outvector(www%localnconfig*numr)
    Type(CONFIGPTR),intent(in) :: matrix_ptr
    Type(SPARSEPTR),intent(in) :: sparse_ptr
    real*8,intent(in) :: time
    
    call sparseconfigmultxxx(www,invector,outvector,matrix_ptr,sparse_ptr,&
         boflag, nucflag, pulseflag, conflag,time,0,imc)

  end subroutine sparseconfigmult

!! All purpose subroutine, just the pulse.

  subroutine sparseconfigpulsemult(www,invector,outvector,matrix_ptr, sparse_ptr,which,imc)
    use r_parameters
    use configptrmod
    use sparseptrmod
    use walkmod
    implicit none
    type(walktype),intent(in) :: www
    integer,intent(in) :: which,imc
    DATATYPE,intent(in) :: invector(www%totadim)
    DATATYPE,intent(out) :: outvector(www%totadim)
    Type(CONFIGPTR),intent(in) :: matrix_ptr
    Type(SPARSEPTR),intent(in) :: sparse_ptr

    call sparseconfigmultxxx(www,invector,outvector,matrix_ptr,sparse_ptr, 0,0,1,0,-1d0,which,imc)

  end subroutine sparseconfigpulsemult


  subroutine parsparseconfigpreconmult(www,invector,outvector,matrix_ptr, sparse_ptr, &
       boflag, nucflag, pulseflag, conflag, inenergy,time)
    use r_parameters
    use configptrmod
    use sparseptrmod
    use walkmod
    use mpimod   !! myrank
    implicit none
    type(walktype),intent(in) :: www
    DATATYPE,intent(in) :: invector(numr,www%botconfig:www%topconfig)
    DATATYPE,intent(out) :: outvector(numr,www%botconfig:www%topconfig)
    DATATYPE,intent(in) :: inenergy
    Type(CONFIGPTR),intent(in) :: matrix_ptr
    Type(SPARSEPTR),intent(in) :: sparse_ptr
    integer,intent(in) ::  conflag,boflag,nucflag,pulseflag
    real*8,intent(in) :: time
    DATATYPE :: workvector(numr,www%botconfig:www%topconfig+1),&
         tempvector(numr,www%botconfig:www%topconfig+1)     !!AUTOMATIC

    if (www%topconfig-www%botconfig+1.eq.0) then
       return
    endif

    workvector(:,:)=1d0; tempvector(:,:)=0d0
    call sparseconfigmult_byproc(myrank,myrank,www,workvector,tempvector,matrix_ptr, sparse_ptr, &
         boflag, nucflag, pulseflag, conflag,time,0,1,numr,1,-1)
  
    outvector(:,:)=invector(:,:)/(tempvector(:,www%botconfig:www%topconfig)-inenergy)

  end subroutine parsparseconfigpreconmult

end module sparsemultmod

