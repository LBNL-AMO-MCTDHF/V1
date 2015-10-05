
!! FOR AUTOCORRELATION.

#include "Definitions.INC"

!! SAVE INITIAL WAVE FUNCTION

module automod
  implicit none
  DATATYPE, allocatable :: overlaps(:,:), orig_spfs(:,:), orig_avectors(:,:,:)
  integer, allocatable :: xcalledflag,  calledflag
end module automod

subroutine autocorrelate_initial()
  use automod
  use parameters
  use xxxmod
  implicit none
  integer :: imc
  
  allocate(overlaps(0:autosize,mcscfnum),  orig_spfs(  spfsize, nspf ),  calledflag,  xcalledflag)
  overlaps(:,:)=0d0
  calledflag=0;  xcalledflag=0
  
  allocate(orig_avectors(numr,firstconfig:lastconfig,mcscfnum))

  orig_spfs(:,:)=RESHAPE(yyy%cmfpsivec(spfstart:spfend,0),(/spfsize,nspf/))
  do imc=1,mcscfnum
     orig_avectors(:,:,imc)=RESHAPE(yyy%cmfpsivec(astart(imc):aend(imc),0),(/numr,localnconfig/))
  enddo
  
end subroutine

subroutine autocorrelate()
  use automod
  use parameters
  use mpimod
  use xxxmod
  implicit none
  integer ::  i,imc,sflag
  integer, save :: lastouttime=0
  real*8 :: thistime

  if (mod(xcalledflag,autosteps).eq.0) then
     do imc=1,mcscfnum
        call autocorrelate_one(yyy%cmfpsivec(astart(imc),0),yyy%cmfpsivec(spfstart,0),orig_spfs(:,:), &
             orig_avectors(:,:,imc), overlaps(calledflag,imc),numr)
     enddo

     if (mod(calledflag,dipmodtime).eq.0.and.calledflag.gt.0) then

        thistime=calledflag*par_timestep*autosteps
        sflag=0
        if (floor(thistime/diptime).gt.lastouttime) then
           lastouttime=floor(thistime/diptime)
           sflag=1
        endif
     
        if (myrank.eq.1) then
           open(881,file="Dat/ovl.dat", status="unknown")
           do i=0,calledflag
              write(881,'(I10,100F20.10)') i,overlaps(i,:)
           enddo
           close(881)
        endif
        
        call autocall(calledflag, overlaps(:,:), sflag)

     endif
     calledflag=calledflag+1
  endif
  xcalledflag=xcalledflag+1
  
end subroutine autocorrelate


function getval(jjj,which)
  use parameters
  implicit none
  
  integer, parameter :: max_numelec=20
  integer, save :: kkk(max_numelec)=0
  integer :: jjj(max_numelec), mm,ll,nn, getval, which,kk

  if (which .gt. numelec) then
     OFLWR "Err which"; CFLST
  endif

  if (numelec.gt.max_numelec) then
     OFLWR "Redim max_numelec in autocorrelate";  CFLST
  endif

  kkk(1)=jjj(1)     ! kkk is the permutation we're building
  do kk=2,which   ! position in permutation
     mm=0           
     do ll=1,jjj(kk)
        mm=mm+1
        do nn=1,kk-1
           if (kkk(nn).eq.mm) then
              mm=mm+1
           endif
        enddo
     enddo
     kkk(kk)=mm
  enddo

  getval=kkk(which)

end function getval

!! CALCULATE THE AUTOCORRELATION FUNCTION.

module autobiomod
  use biorthotypemod
  implicit none
  type(biorthotype),target :: autobiovar
end module

subroutine autocorrelate_one(avector,inspfs,orig_spf,orig_avector,inoverlaps,innr)
  use autobiomod
  use biorthomod
  use configmod
  use parameters
  implicit none

  DATATYPE, allocatable :: mobio(:,:),abio(:,:)
  DATATYPE :: inspfs(  spfsize, nspf ), dot
  DATATYPE :: orig_spf(  spfsize, nspf ), orig_avector(innr,firstconfig:lastconfig), inoverlaps, avector(innr,firstconfig:lastconfig)
  DATATYPE,target :: smo(nspf,nspf)
  integer :: innr

  allocate(mobio(spfsize,nspf),abio(innr,firstconfig:lastconfig))

  if(auto_biortho.eq.0) then
     OFLWR "Permoverlaps disabled.  set auto_biortho=1"; CFLST

!!$     call permoverlaps(innr, numelec, spfsize, inspfs, orig_spf, avector,orig_avector,inoverlaps, 1, &
!!$          autopermthresh, autonormthresh, nspf, nspf, numconfig, numconfig, configlist, ndof, configlist,&
!!$          ndof, 0, parorbsplit)

  else

     abio=orig_avector

     call bioset(autobiovar,smo,innr)
     call biortho(orig_spf,inspfs,mobio,abio,autobiovar)

     inoverlaps=dot(avector,abio,localnconfig*innr)
     if (parconsplit.ne.0) then
        call mympireduceone(inoverlaps)
     endif
  endif

  deallocate(mobio,abio)

end subroutine autocorrelate_one


subroutine autocorrelate_final()
  use automod
  implicit none
  deallocate(orig_spfs,   orig_avectors,  overlaps)
end subroutine autocorrelate_final



!!$  !! MULTI-USE OVERLAP ROUTINE, BEFORE BIORTHO.  MANUAL LOOPS
!!$  !!  THROUGH PERMUTATIONS.
!!$  !!070610 subroutine autocorrelate_one(avector,inspfs,orig_spf,orig_avector,overlap,calledflag)
!!$  !!070610b  subroutine permoverlaps(spfmat, avector,inspfs,orig_spf,orig_avector,inoverlaps, number)
!!$  !! With transformflag=1, returns inavector
!!$  !! With transformflag=0, returns inoverlap given inavector.
!!$  
!!$  
!!$  subroutine permoverlaps(rdimension, both_numelec, spfsize, inspfs,orig_spfs, inavector,orig_avector,inoverlap,&
!!$       printflag, sumthresh, normthresh, in_nspf, orig_nspf, in_numconfig, orig_numconfig, in_configlist, &
!!$       inconfiglistlda, orig_configlist, origconfiglistlda, transformflag,parorbsplit)
!!$  
!!$  !! NO PARAMETERS !!
!!$  
!!$    implicit none
!!$  
!!$    integer, parameter :: max_numelec=20
!!$    integer,intent(in) :: both_numelec, spfsize, inconfiglistlda, origconfiglistlda, rdimension, &
!!$         orig_nspf, orig_numconfig, orig_configlist(origconfiglistlda, orig_numconfig), &
!!$         in_nspf, in_numconfig, in_configlist(inconfiglistlda, in_numconfig),parorbsplit
!!$    DATATYPE, intent(in) :: inspfs(  spfsize, in_nspf ), orig_spfs(  spfsize, orig_nspf ), &
!!$         orig_avector(orig_numconfig,rdimension)
!!$    DATATYPE, intent(out) :: inoverlap
!!$    DATATYPE, intent(inout) :: inavector(in_numconfig,rdimension)  
!!$    DATATYPE :: lsums(0:max_numelec), dot
!!$    integer :: i,j, reorder, phase,ir,  getval,thisconfig(2*both_numelec), thatconfig(2*both_numelec), &
!!$         thootconfig(2*both_numelec), qqq, printflag, inskipped, origskipped, k, l, transformflag
!!$    integer, pointer :: jj1,jj2,jj3,jj4,jj5,jj6,jj7,jj8,jj9,jj10, jj11,jj12, jj13, jj14, jj15, jj16, jj17, &
!!$         jj18, jj19, jj20, kk2,kk3,kk4,kk5,kk6,kk7,kk8,kk9,kk10, kk11,kk12, kk13, kk14, kk15, kk16, kk17, &
!!$         kk18, kk19, kk20
!!$    integer, target :: jjj(max_numelec), kkk(max_numelec)
!!$    real*8 :: nperms, permavg, permavgcount, sumthresh, normthresh, rsum, rmax
!!$    DATATYPE :: spfmat(in_nspf,orig_nspf) 
!!$  
!!$  
!!$    print *, "AUTOCORRELATE_ONE DISABLED....  NOT AGREEING W/BIORTHO IN REPLACE_WITHNAT"; call mpistop()
!!$  
!!$  !  replace_withnat works because leads to krylov termination at order 1 with spfflag=0
!!$  !     that's using biotransform.
!!$  
!!$    do i=1,in_nspf
!!$       do j=1,orig_nspf
!!$          spfmat(i,j) = dot(inspfs(:,i),orig_spfs(:,j),spfsize)
!!$       enddo
!!$    enddo
!!$  
!!$    if (parorbsplit.eq.3) then
!!$       call mympireduce(spfmat,in_nspf*orig_nspf)
!!$    endif
!!$  
!!$    if (transformflag==1) then
!!$       inavector=0.d0
!!$    else
!!$       inoverlap=0.d0
!!$    endif
!!$   
!!$    jj1 => jjj(1)
!!$    jj2 => jjj(2)
!!$    jj3 => jjj(3)
!!$    jj4 => jjj(4)
!!$    jj5 => jjj(5)
!!$    jj6 => jjj(6)
!!$    jj7 => jjj(7)
!!$    jj8 => jjj(8)
!!$    jj9 => jjj(9)
!!$    jj10 => jjj(10)
!!$    jj11 => jjj(11)
!!$    jj12 => jjj(12)
!!$    jj13 => jjj(13)
!!$    jj14 => jjj(14)
!!$    jj15 => jjj(15)
!!$    jj16 => jjj(16)
!!$    jj17 => jjj(17)
!!$    jj18 => jjj(18)
!!$    jj19 => jjj(19)
!!$    jj20 => jjj(20)
!!$  
!!$  
!!$    kk2 => kkk(2)
!!$    kk3 => kkk(3)
!!$    kk4 => kkk(4)
!!$    kk5 => kkk(5)
!!$    kk6 => kkk(6)
!!$    kk7 => kkk(7)
!!$    kk8 => kkk(8)
!!$    kk9 => kkk(9)
!!$    kk10 => kkk(10)
!!$    kk11 => kkk(11)
!!$    kk12 => kkk(12)
!!$    kk13 => kkk(13)
!!$    kk14 => kkk(14)
!!$    kk15 => kkk(15)
!!$    kk16 => kkk(16)
!!$    kk17 => kkk(17)
!!$    kk18 => kkk(18)
!!$    kk19 => kkk(19)
!!$    kk20 => kkk(20)
!!$  
!!$  
!!$    permavg=0
!!$    permavgcount=0
!!$  
!!$    if (both_numelec.gt.max_numelec) then
!!$       write(*, *) "Redim max_numelec in autocorrelate";    call mpistop()
!!$    endif
!!$  
!!$  !!!  jjj is a PERMUTATION.  sum over all permutations.
!!$  !!!
!!$  
!!$  
!!$    inskipped=0
!!$    do i=1,in_numconfig
!!$       thisconfig=in_configlist(1:2*both_numelec,i)
!!$  
!!$       if (transformflag==1) then
!!$          rsum=10.d0
!!$       else
!!$          rsum=0.d0
!!$          do l=1,rdimension
!!$             rmax= abs(inavector(i,l))**2
!!$             rsum=rsum+rmax
!!$          enddo
!!$          rsum=sqrt(rsum)
!!$       endif
!!$  
!!$       if (rsum.lt.normthresh) then
!!$          inskipped=inskipped+1
!!$       else
!!$  
!!$       origskipped=0
!!$       do j=1,orig_numconfig
!!$          thatconfig=orig_configlist(1:2*both_numelec,j)
!!$  
!!$          rsum=0.d0
!!$          do k=1,rdimension
!!$             rmax= abs(orig_avector(j,k))**2
!!$             rsum=rsum+rmax
!!$          enddo
!!$          rsum=sqrt(rsum)
!!$          
!!$          if (rsum.lt.normthresh) then
!!$             origskipped=origskipped+1
!!$          else
!!$  
!!$  
!!$          nperms=0.d0
!!$  
!!$          lsums(0)=1.d0
!!$  
!!$          do jj1=1,both_numelec
!!$             lsums(1)=lsums(0)
!!$             if (1.le.both_numelec) then
!!$                thootconfig(1:2)=thatconfig(jj1*2-1:jj1*2)
!!$                if (thootconfig(2).eq.thisconfig(2)) then
!!$                   lsums(1)=lsums(1)*spfmat(thisconfig(1),thootconfig(1))
!!$                else
!!$                   lsums(1)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(1)) .gt.sumthresh) then
!!$          do jj2=1,max(both_numelec-1,1)
!!$             lsums(2)=lsums(1)
!!$             if (2.le.both_numelec) then 
!!$                kk2=getval(jjj,2)
!!$                thootconfig(3:4)=thatconfig(kk2*2-1:kk2*2)
!!$                if (thootconfig(4).eq.thisconfig(4)) then
!!$                   lsums(2)=lsums(2)*spfmat(thisconfig(3),thootconfig(3))
!!$                else
!!$                   lsums(2)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(2)) .gt.sumthresh) then
!!$          do jj3=1,max(both_numelec-2,1)
!!$             lsums(3)=lsums(2)
!!$             if (3.le.both_numelec) then
!!$                kk3=getval(jjj,3)
!!$                thootconfig(5:6)=thatconfig(kk3*2-1:kk3*2)
!!$                if (thootconfig(6).eq.thisconfig(6)) then
!!$                   lsums(3)=lsums(3)*spfmat(thisconfig(5),thootconfig(5))
!!$                else
!!$                   lsums(3)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(3)) .gt.sumthresh) then
!!$          do jj4=1,max(both_numelec-3,1)
!!$             lsums(4)=lsums(3)
!!$             if (4.le.both_numelec) then
!!$                kk4=getval(jjj,4)
!!$                thootconfig(7:8)=thatconfig(kk4*2-1:kk4*2)
!!$                if (thootconfig(8).eq.thisconfig(8)) then
!!$                   lsums(4)=lsums(4)*spfmat(thisconfig(7),thootconfig(7))
!!$                else
!!$                   lsums(4)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(4)) .gt.sumthresh) then
!!$          do jj5=1,max(both_numelec-4,1)
!!$             lsums(5)=lsums(4)
!!$             if (5.le.both_numelec) then
!!$                kk5=getval(jjj,5)
!!$                thootconfig(9:10)=thatconfig(kk5*2-1:kk5*2)
!!$                if (thootconfig(10).eq.thisconfig(10)) then
!!$                   lsums(5)=lsums(5)*spfmat(thisconfig(9),thootconfig(9))
!!$                else
!!$                   lsums(5)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(5)) .gt.sumthresh) then
!!$          do jj6=1,max(both_numelec-5,1)
!!$             lsums(6)=lsums(5)
!!$             if (6.le.both_numelec) then
!!$                kk6=getval(jjj,6)
!!$                thootconfig(11:12)=thatconfig(kk6*2-1:kk6*2)
!!$                if (thootconfig(12).eq.thisconfig(12)) then
!!$                   lsums(6)=lsums(6)*spfmat(thisconfig(11),thootconfig(11))
!!$                else
!!$                   lsums(6)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(6)) .gt.sumthresh) then
!!$          do jj7=1,max(both_numelec-6,1)
!!$             lsums(7)=lsums(6)
!!$             if (7.le.both_numelec) then
!!$                kk7=getval(jjj,7)
!!$                thootconfig(13:14)=thatconfig(kk7*2-1:kk7*2)
!!$                if (thootconfig(14).eq.thisconfig(14)) then
!!$                   lsums(7)=lsums(7)*spfmat(thisconfig(13),thootconfig(13))
!!$                else
!!$                   lsums(7)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(7)) .gt.sumthresh) then
!!$          do jj8=1,max(both_numelec-7,1)
!!$  
!!$             lsums(8)=lsums(7)
!!$             if (8.le.both_numelec) then
!!$                kk8=getval(jjj,8)
!!$                thootconfig(15:16)=thatconfig(kk8*2-1:kk8*2)
!!$                if (thootconfig(16).eq.thisconfig(16)) then
!!$                   lsums(8)=lsums(8)*spfmat(thisconfig(15),thootconfig(15))
!!$                else
!!$                   lsums(8)=0.d0
!!$                endif
!!$             endif
!!$  
!!$             if (abs(lsums(8)) .gt.sumthresh) then
!!$          do jj9=1,max(both_numelec-8,1)
!!$  
!!$             lsums(9)=lsums(8)
!!$             if (9.le.both_numelec) then
!!$                kk9=getval(jjj,9)
!!$                thootconfig(17:18)=thatconfig(kk9*2-1:kk9*2)
!!$                if (thootconfig(18).eq.thisconfig(18)) then
!!$                   lsums(9)=lsums(9)*spfmat(thisconfig(17),thootconfig(17))
!!$                else
!!$                   lsums(9)=0.d0
!!$                endif
!!$             endif
!!$                if (abs(lsums(9)) .gt.sumthresh) then
!!$          do jj10=1,max(both_numelec-9,1)
!!$             lsums(10)=lsums(9)
!!$             if (10.le.both_numelec) then
!!$                kk10=getval(jjj,10)
!!$                thootconfig(19:20)=thatconfig(kk10*2-1:kk10*2)
!!$                if (thootconfig(20).eq.thisconfig(20)) then
!!$                   lsums(10)=lsums(10)*spfmat(thisconfig(19),thootconfig(19))
!!$                else
!!$                   lsums(10)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(10)) .gt.sumthresh) then
!!$          do jj11=1,max(both_numelec-10,1)
!!$             lsums(11)=lsums(10)
!!$             if (11.le.both_numelec) then
!!$                kk11=getval(jjj,11)
!!$                thootconfig(21:22)=thatconfig(kk11*2-1:kk11*2)
!!$                if (thootconfig(22).eq.thisconfig(22)) then
!!$                   lsums(11)=lsums(11)*spfmat(thisconfig(21),thootconfig(21))
!!$                else
!!$                   lsums(11)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(11)) .gt.sumthresh) then
!!$          do jj12=1,max(both_numelec-11,1)
!!$             lsums(12)=lsums(11)
!!$             if (12.le.both_numelec) then
!!$                kk12=getval(jjj,12)
!!$                thootconfig(23:24)=thatconfig(kk12*2-1:kk12*2)
!!$                if (thootconfig(24).eq.thisconfig(24)) then
!!$                   lsums(12)=lsums(12)*spfmat(thisconfig(23),thootconfig(23))
!!$                else
!!$                   lsums(12)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(12)) .gt.sumthresh) then
!!$          do jj13=1,max(both_numelec-12,1)
!!$             lsums(13)=lsums(12)
!!$             if (13.le.both_numelec) then
!!$                kk13=getval(jjj,13)
!!$                thootconfig(25:26)=thatconfig(kk13*2-1:kk13*2)
!!$                if (thootconfig(26).eq.thisconfig(26)) then
!!$                   lsums(13)=lsums(13)*spfmat(thisconfig(25),thootconfig(25))
!!$                else
!!$                   lsums(13)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(13)) .gt.sumthresh) then
!!$          do jj14=1,max(both_numelec-13,1)
!!$             lsums(14)=lsums(13)
!!$             if (14.le.both_numelec) then
!!$                kk14=getval(jjj,14)
!!$                thootconfig(27:28)=thatconfig(kk14*2-1:kk14*2)
!!$                if (thootconfig(28).eq.thisconfig(28)) then
!!$                   lsums(14)=lsums(14)*spfmat(thisconfig(27),thootconfig(27))
!!$                else
!!$                   lsums(14)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(14)) .gt.sumthresh) then
!!$          do jj15=1,max(both_numelec-14,1)
!!$             lsums(15)=lsums(14)
!!$             if (15.le.both_numelec) then
!!$                kk15=getval(jjj,15)
!!$                thootconfig(29:30)=thatconfig(kk15*2-1:kk15*2)
!!$                if (thootconfig(30).eq.thisconfig(30)) then
!!$                   lsums(15)=lsums(15)*spfmat(thisconfig(29),thootconfig(29))
!!$                else
!!$                   lsums(15)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(15)) .gt.sumthresh) then
!!$  
!!$  
!!$  
!!$          do jj16=1,max(both_numelec-15,1)
!!$             qqq=16
!!$             lsums(qqq)=lsums(qqq-1)
!!$             if (qqq.le.both_numelec) then
!!$                kkk(qqq)=getval(jjj,qqq)
!!$                thootconfig(qqq*2-1:qqq*2)=thatconfig(kkk(qqq)*2-1:kkk(qqq)*2)
!!$                if (thootconfig(qqq*2).eq.thisconfig(qqq*2)) then
!!$                   lsums(qqq)=lsums(qqq)*spfmat(thisconfig(qqq*2-1),thootconfig(qqq*2-1))
!!$                else
!!$                   lsums(qqq)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(qqq)) .gt.sumthresh) then
!!$  
!!$          do jj17=1,max(both_numelec-16,1)
!!$             qqq=17
!!$             lsums(qqq)=lsums(qqq-1)
!!$             if (qqq.le.both_numelec) then
!!$                kkk(qqq)=getval(jjj,qqq)
!!$                thootconfig(qqq*2-1:qqq*2)=thatconfig(kkk(qqq)*2-1:kkk(qqq)*2)
!!$                if (thootconfig(qqq*2).eq.thisconfig(qqq*2)) then
!!$                   lsums(qqq)=lsums(qqq)*spfmat(thisconfig(qqq*2-1),thootconfig(qqq*2-1))
!!$                else
!!$                   lsums(qqq)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(qqq)) .gt.sumthresh) then
!!$  
!!$          do jj18=1,max(both_numelec-17,1)
!!$             qqq=18
!!$             lsums(qqq)=lsums(qqq-1)
!!$             if (qqq.le.both_numelec) then
!!$                kkk(qqq)=getval(jjj,qqq)
!!$                thootconfig(qqq*2-1:qqq*2)=thatconfig(kkk(qqq)*2-1:kkk(qqq)*2)
!!$                if (thootconfig(qqq*2).eq.thisconfig(qqq*2)) then
!!$                   lsums(qqq)=lsums(qqq)*spfmat(thisconfig(qqq*2-1),thootconfig(qqq*2-1))
!!$                else
!!$                   lsums(qqq)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(qqq)) .gt.sumthresh) then
!!$  
!!$          do jj19=1,max(both_numelec-18,1)
!!$             qqq=19
!!$             lsums(qqq)=lsums(qqq-1)
!!$             if (qqq.le.both_numelec) then
!!$                kkk(qqq)=getval(jjj,qqq)
!!$                thootconfig(qqq*2-1:qqq*2)=thatconfig(kkk(qqq)*2-1:kkk(qqq)*2)
!!$                if (thootconfig(qqq*2).eq.thisconfig(qqq*2)) then
!!$                   lsums(qqq)=lsums(qqq)*spfmat(thisconfig(qqq*2-1),thootconfig(qqq*2-1))
!!$                else
!!$                   lsums(qqq)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(qqq)) .gt.sumthresh) then
!!$  
!!$          do jj20=1,max(both_numelec-19,1)
!!$             qqq=20
!!$             lsums(qqq)=lsums(qqq-1)
!!$             if (qqq.le.both_numelec) then
!!$                kkk(qqq)=getval(jjj,qqq)
!!$                thootconfig(qqq*2-1:qqq*2)=thatconfig(kkk(qqq)*2-1:kkk(qqq)*2)
!!$                if (thootconfig(qqq*2).eq.thisconfig(qqq*2)) then
!!$                   lsums(qqq)=lsums(qqq)*spfmat(thisconfig(qqq*2-1),thootconfig(qqq*2-1))
!!$                else
!!$                   lsums(qqq)=0.d0
!!$                endif
!!$             endif
!!$             if (abs(lsums(qqq)) .gt.sumthresh) then
!!$                   phase=reorder(thootconfig)
!!$                   if (transformflag==1) then
!!$                      do ir=1,rdimension 
!!$                         inavector(i,ir)=inavector(i,ir) + lsums(qqq) * orig_avector(j,ir) * phase
!!$                      enddo
!!$                   else
!!$                      do ir=1,rdimension
!!$                         inoverlap=inoverlap+lsums(qqq)*ALLCON(inavector(i,ir)) * orig_avector(j,ir) * phase
!!$                      enddo
!!$                   endif
!!$                   nperms=nperms+1
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$   !!*
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$                  endif
!!$          enddo
!!$          permavg=permavg+nperms  !!/floatfac((numelec-restrictms)/2)/floatfac((numelec+restrictms)/2)
!!$       endif
!!$       enddo
!!$       permavgcount=permavgcount+1.d0
!!$    endif
!!$  enddo
!!$  
!!$  if (printflag==1) then
!!$     write(*,*) "Overlaps done, avg # of perms per config ", permavg/permavgcount
!!$  !   write(*,*) "  Number skipped and proportion due to insufficient A-norm:"
!!$  !   write(*,*) "  Current: ", inskipped, in_numconfig, real(inskipped,8)/real(in_numconfig,8)
!!$  !   write(*,*) "  Original ", origskipped, orig_numconfig, real(origskipped,8)/real(orig_numconfig,8)
!!$  endif
!!$  
!!$  end subroutine permoverlaps

