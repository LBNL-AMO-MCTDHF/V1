
#include "Definitions.INC"

function qbox(notusedint)
  implicit none
  integer :: qbox,notusedint,notusedint2
  notusedint2=0*notusedint !! avoid warn unused
  qbox=1
end function qbox


subroutine transferparams(innumspf,inspfrestrictflag,inspfmvals,inspfugrestrict,inspfugvals,outspfsmallsize,logorbpar)
  use myparams
  implicit none
  integer,intent(in) :: innumspf,inspfrestrictflag,inspfmvals(innumspf), inspfugrestrict,inspfugvals(innumspf)
  integer,intent(out) :: outspfsmallsize
  integer :: ii
  logical, intent(out) :: logorbpar

  logorbpar=.false.

  numspf=innumspf;  spfrestrictflag=inspfrestrictflag; spfugrestrict=inspfugrestrict

  if (spfugrestrict.ne.0) then
     if (spfrestrictflag.eq.0) then
        OFLWR "spfugrestrict requires spfrestrictflag ", spfugrestrict, spfrestrictflag; CFLST
     endif
     do ii=1,numspf
        if (abs(inspfugvals(ii)).ne.1) then
           OFLWR "TWOECHECK UG ",ii,inspfugvals(:); CFLST
        endif
     enddo
!!     outspfsmallsize=(lbig+1)/2*numerad
     outspfsmallsize=(lbig+2)/2*numerad     !! lbig even ok
  else
     if (spfrestrictflag.ne.0) then
        outspfsmallsize=(lbig+1)*numerad
     else
        outspfsmallsize=(lbig+1)*numerad*(2*mbig+1)
     endif
  endif

  allocate(spfmvals(numspf));  spfmvals(:)=inspfmvals(:)
  allocate(spfugvals(numspf));  spfugvals(:)=inspfugvals(:)

end subroutine transferparams



module mycdotmod
contains
recursive function mycdot(one,two,n)
  implicit none
  integer,intent(in) :: n
  DATATYPE,intent(in) :: one(n), two(n)
  DATATYPE :: mycdot, sum
  integer :: i
  sum=0.d0
  do i=1,n
     sum = sum + one(i) * two(i) 
  enddo
  mycdot=sum
end function mycdot
recursive function mydot(one,two,n)
  implicit none
  integer,intent(in) :: n
  DATATYPE,intent(in) :: one(n), two(n)
  DATATYPE :: mydot, sum
  integer :: i
  sum=0.d0
  do i=1,n
     sum = sum + CONJUGATE(one(i)) * two(i) 
  enddo
  mydot=sum
end function mydot
end module

!  DATAECS :: rmatrix(numerad,numerad,mseriesmax+1,lseriesmax+1)
!  real*8 :: ylmvals(0:2*mbig, 1:lbig+1, lseriesmax+1)

module tinvsubmod
contains

subroutine op_tinv(mlow,mhigh,howmany,indensity,outpotential)
  use myparams
  use myprojectmod
  implicit none
  integer, intent(in) :: mlow,mhigh,howmany
  DATATYPE,intent(in) :: indensity(numerad,lbig+1,mlow:mhigh,howmany)
  DATATYPE,intent(out) :: outpotential(numerad,lbig+1,mlow:mhigh,howmany)
  integer :: deltam,lsum,j1,j2,k1,ii
  DATATYPE :: twoeden2(numerad),twoeden3(numerad),mydensity(numerad,lbig+1),&
       mypotential(numerad,lbig+1)    !! AUTOMATIC

  outpotential(:,:,:,:)=0d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(deltam,lsum,j1,j2,k1,ii,twoeden2,twoeden3,mydensity,mypotential)

!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
  do ii=1,howmany
     do deltam=mlow,mhigh
        mydensity(:,:)=indensity(:,:,deltam,ii)
        mypotential(:,:)=0d0
        do lsum=1+abs(deltam),jacobisummax+1
           twoeden2(:)=0d0
           do j1=1,lbig+1
              twoeden2(:)=twoeden2(:) + mydensity(:,j1) * ylmvals(abs(deltam),j1,lsum-abs(deltam))
           enddo
           twoeden3(:)=0d0
           do j2=1,numerad
              if (atomflag==0) then
                 twoeden3(:)=twoeden3(:) + twoeden2(j2) * rmatrix(:,j2,abs(deltam)+1,lsum-abs(deltam)) * 4.d0
              else
                 twoeden3(:)=twoeden3(:) + twoeden2(j2) * rmatrix(:,j2,1,lsum) 
              endif
           enddo

           do k1=1,lbig+1
              mypotential(:,k1) = mypotential(:,k1) - &
                   twoeden3(:)*ylmvals(abs(deltam),k1,lsum-abs(deltam))
           enddo
        enddo
        outpotential(:,:,deltam,ii)=mypotential(:,:)
     enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

end subroutine op_tinv

end module tinvsubmod

!  DATAECS :: rmatrix(numerad,numerad,mseriesmax+1,lseriesmax+1)
!  real*8 :: ylmvals(0:2*mbig, 1:lbig+1, lseriesmax+1)

subroutine call_twoe_matel00(lowspf,highspf,inspfs1,inspfs2,twoematel,twoereduced,xtimingdir,xnotiming) !! ok unused
  use myparams
  use myprojectmod
  use mycdotmod
  use clockmod   !! IN PARENT DIRECTORY
  use tinvsubmod
  implicit none
  integer, intent(in) :: xnotiming,lowspf,highspf
  DATATYPE,intent(out) :: twoereduced(numerad,lbig+1,-2*mbig:2*mbig, numspf,lowspf:highspf),&
       twoematel(numspf,numspf,numspf,lowspf:highspf)
  DATATYPE,intent(in) :: inspfs1(numerad,lbig+1,-mbig:mbig,numspf),inspfs2(numerad,lbig+1,-mbig:mbig,numspf)
  character,intent(in) :: xtimingdir*(*)
  integer :: mvalue2a, mvalue1b, mvalue2b, mvalue1a, itime, jtime, &
       spf1a,spf1b,spf2a,spf2b, deltam,qq,rr,qq2,rr2,times(100),nnnspf
  DATATYPE :: myden(numerad,lbig+1,-2*mbig:2*mbig),myred(numerad,lbig+1,-2*mbig:2*mbig)   !! AUTOMATIC
  DATATYPE,allocatable :: twoeden(:,:,:,:)

  nnnspf=highspf-lowspf+1

  call myclock(itime)

  if (highspf.gt.numspf.or.lowspf.lt.1) then
     print *, "programmer fail xxx",lowspf,highspf,numspf; stop
  endif

  allocate(twoeden(numerad,lbig+1,-2*mbig:2*mbig,numspf))
  twoeden(:,:,:,:)=0d0
  if (nnnspf.gt.0) then
     twoematel(:,:,:,:)=0d0; twoereduced(:,:,:,:,:)=0d0
  endif

  do spf2b=lowspf,highspf
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(spf2a,qq,qq2,rr,rr2,myden,mvalue2a,mvalue2b,deltam)
!$OMP DO SCHEDULE(DYNAMIC)
     do spf2a=1,numspf
        ! integrating over electron 2

        if (spfrestrictflag==1) then
           qq=spfmvals(spf2a);           rr=spfmvals(spf2a)
           qq2=spfmvals(spf2b);           rr2=spfmvals(spf2b)
        else
           qq=-mbig;           rr=mbig
           qq2=-mbig;           rr2=mbig
        endif

        myden(:,:,:)=0d0
        do mvalue2a=qq,rr
           do mvalue2b=qq2,rr2
              deltam=mvalue2b-mvalue2a
              myden(:,:,deltam) = myden(:,:,deltam) + &
!! switching 2-2016       CONJUGATE(inspfs1(:,:,mvalue2a,spf2a)) * inspfs2(:,:,mvalue2b,spf2b)
                   inspfs1(:,:,mvalue2a,spf2a) * CONJUGATE(inspfs2(:,:,mvalue2b,spf2b))
           enddo
        enddo
        twoeden(:,:,:,spf2a)=myden(:,:,:)
     enddo  !! spf2a
!$OMP END DO
!$OMP END PARALLEL

     call op_tinv(-2*mbig,2*mbig,numspf,twoeden,twoereduced(:,:,:,:,spf2b))

  enddo

  call myclock(jtime);           times(1)=times(1)+jtime-itime;    itime=jtime

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(spf1a,spf1b,spf2a,spf2b,qq,qq2,rr,rr2,myden,myred,mvalue1a,mvalue1b,deltam)
!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
  do spf1b=1,numspf
     do spf1a=1,numspf
              
        if (spfrestrictflag==1) then
           qq=spfmvals(spf1a);                    rr=spfmvals(spf1a)
           qq2=spfmvals(spf1b);                    rr2=spfmvals(spf1b)
        else
           qq=-mbig;                 rr=mbig
           qq2=-mbig;                 rr2=mbig
        endif

        myden(:,:,:)=0d0
        do mvalue1a=qq,rr
           do mvalue1b=qq2,rr2

              deltam=mvalue1a-mvalue1b
              myden(:,:,deltam) = myden(:,:,deltam) + &
!! switching 2-2016        CONJUGATE(inspfs1(:,:,mvalue1a,spf1a)) * inspfs2(:,:,mvalue1b,spf1b)
                   inspfs1(:,:,mvalue1a,spf1a) * CONJUGATE(inspfs2(:,:,mvalue1b,spf1b))
           enddo
        enddo

        do spf2b=lowspf,highspf
           do spf2a=1,numspf
              myred(:,:,:)=twoereduced(:,:,:,spf2a,spf2b)

!! same diff it's symmetric         twoematel(spf2a,spf2b,spf1a,spf1b) = ...

              twoematel(spf1a,spf1b,spf2a,spf2b) = mycdot(myden(:,:,:),myred(:,:,:),numerad*(lbig+1)*(4*mbig+1))
              
           enddo
        enddo
     enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

  call myclock(itime);   times(2)=times(2)+itime-jtime

  deallocate(twoeden)

!  if ((myrank.eq.1).and.(notiming.eq.0)) then
!     if (timingflag==1) then
!        if (mod(numcalledhere,mytimingout).eq.0)  then
!           open(853,file="twoe.time.dat",status="unknown", position="append")
!           write(853, '(F13.3, T16, 100I15)') time, times(1:3)/1000
!           close(853)
!        endif
!     endif
!  endif

!  if (numcalledhere.eq.1) then
!     if (checknan2(twoematel,numspf**4)) then
!        call openfile()
!        write(mpifileptr,*) "Error.  Twoematel has NaN's."
!        write(mpifileptr,*) "For the moment, fix this by changing the grid."
!        write(mpifileptr,*) "Or, compile without -fast-math."
!        CFLST
!     endif
!  endif


end subroutine call_twoe_matel00



subroutine hatom_op(howmany,inspfs, outspfs, hatomreduced)
  use myparams
  use myprojectmod
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: inspfs(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE,intent(out) :: outspfs(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE,intent(in) :: hatomreduced(numerad,lbig+1,-2*mbig:2*mbig)
  integer ::  mvalue1b, mvalue1a,deltam,ivect

  outspfs(:,:,:,:)=0d0

  if (numhatoms.eq.0) then
     return
  endif

  if (spfrestrictflag==1) then
     OFLWR "Hey, don't use spfrestrictflag if you have hatoms!";CFLST
  endif
  if (numr.gt.1) then
     print *, "Hatom not for calc with multiple r gridpoints.sorry.";     stop
  endif

  do ivect=1,howmany
     do mvalue1a=-mbig,mbig
        do mvalue1b=-mbig,mbig
           deltam=mvalue1a-mvalue1b
           outspfs(:,:,mvalue1a,ivect)= outspfs(:,:,mvalue1a,ivect)- &
                inspfs(:,:,mvalue1b,ivect) * hatomreduced(:,:,deltam)
        enddo
     enddo
  enddo

end subroutine hatom_op


subroutine op_contact(howmany,notusedin,spfout)
  use myparams
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: notusedin(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE,intent(out) :: spfout(numerad,lbig+1,-mbig:mbig,howmany)
  spfout(:,:,:,:) = 0d0
end subroutine op_contact


subroutine call_frozen_matels_core(infrozens,numfrozen,frozenkediag,frozenpotdiag,frozenreduced,hatomreduced)
  use myparams
  use myprojectmod
  use mycdotmod
  use orbmultsubmod   !! IN PARENT DIRECTORY
  use tinvsubmod
  implicit none
  integer,intent(in) :: numfrozen
  DATATYPE,intent(in) :: infrozens(numerad,lbig+1,-mbig:mbig,numfrozen),&
       hatomreduced(numerad,lbig+1,-2*mbig:2*mbig)
  DATATYPE,intent(out) :: frozenkediag,frozenpotdiag
  DATATYPE,intent(out) :: frozenreduced(numerad,lbig+1,-2*mbig:2*mbig)
  DATATYPE :: direct(numfrozen,numfrozen), exch(numfrozen,numfrozen)
  integer :: mvalue2a, mvalue2b,  spf2a,spf2b, deltam,i,  iispf,ispf, sizespf, bigsize
  DATATYPE, allocatable :: frodensity(:,:,:,:), tempreduced(:,:,:,:), cfrodensity(:,:,:,:), &
       tempmult(:,:), tempmult2(:,:)
  DATATYPE :: myden(numerad,lbig+1,-2*mbig:2*mbig),temppotmatel(numfrozen)

  frozenkediag=0;    frozenpotdiag=0

  if (numfrozen.eq.0) then
     return
  endif

  sizespf=numerad*(lbig+1)*(2*mbig+1)
  bigsize=numerad*(lbig+1)*(4*mbig+1)

  allocate(frodensity(numerad,lbig+1,-2*mbig:2*mbig,numfrozen),&
       tempreduced(numerad,lbig+1,-2*mbig:2*mbig,numfrozen), &
       cfrodensity(numerad,lbig+1,-2*mbig:2*mbig,numfrozen), &
       tempmult(sizespf,numfrozen), tempmult2(sizespf,numfrozen))
  frodensity=0d0;  tempreduced=0; cfrodensity=0; tempmult=0; tempmult2=0

  frodensity(:,:,:,1)=0d0
  do spf2a=1,numfrozen
     myden(:,:,:)=0d0
     do mvalue2a=-mbig,mbig
        do mvalue2b=-mbig,mbig
           deltam=mvalue2b-mvalue2a
           myden(:,:,deltam) = myden(:,:,deltam) + &
                CONJUGATE(infrozens(:,:,mvalue2a,spf2a)) * infrozens(:,:,mvalue2b,spf2a)
        enddo
     enddo
     frodensity(:,:,:,1)=frodensity(:,:,:,1)+myden(:,:,:) * 2
  enddo  !! spf2a

  frozenreduced=0
  call op_tinv(-2*mbig,2*mbig,1,frodensity(:,:,:,1),frozenreduced(:,:,:))

  exch(:,:)=0d0
  do spf2b=1,numfrozen
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(spf2a,myden,mvalue2a,mvalue2b,deltam)
!$OMP DO SCHEDULE(DYNAMIC)
     do spf2a=1,numfrozen
        ! integrating over electron 2

        myden(:,:,:)=0d0
        do mvalue2a=-mbig,mbig
           do mvalue2b=-mbig,mbig
              deltam=mvalue2b-mvalue2a
              myden(:,:,deltam) = myden(:,:,deltam) + &
                   CONJUGATE(infrozens(:,:,mvalue2a,spf2a)) * infrozens(:,:,mvalue2b,spf2b)
           enddo
        enddo
        frodensity(:,:,:,spf2a)=myden(:,:,:)
     enddo  !! spf2a
!$OMP END DO
!$OMP END PARALLEL

     cfrodensity(:,:,:,:)=CONJUGATE(frodensity(:,:,:,:))

     call op_tinv(-2*mbig,2*mbig,numfrozen,frodensity(:,:,:,:),tempreduced(:,:,:,:))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(spf2a)
!$OMP DO SCHEDULE(DYNAMIC)
     do spf2a=1,numfrozen
        exch(spf2a,spf2b)=mycdot(cfrodensity(:,:,:,spf2a),tempreduced(:,:,:,spf2a),bigsize)
     enddo
!$OMP END DO
!$OMP END PARALLEL
  enddo  !! spf2b


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(spf2a,myden,mvalue2a,mvalue2b,deltam)
!$OMP DO SCHEDULE(DYNAMIC)
  do spf2a=1,numfrozen
     myden(:,:,:)=0d0
     do mvalue2a=-mbig,mbig
        do mvalue2b=-mbig,mbig
           deltam=mvalue2b-mvalue2a
           myden(:,:,deltam) = myden(:,:,deltam) + &
                CONJUGATE(infrozens(:,:,mvalue2a,spf2a)) * infrozens(:,:,mvalue2b,spf2a)
        enddo
     enddo
     frodensity(:,:,:,spf2a)=myden(:,:,:)
  enddo  !! spf2a
!$OMP END DO
!$OMP END PARALLEL

  call op_tinv(-2*mbig,2*mbig,numfrozen,frodensity(:,:,:,:),tempreduced(:,:,:,:))

  direct(:,:)=0d0
  call MYGEMM('T','N',numfrozen,numfrozen,bigsize,DATAONE,frodensity(:,:,:,:),bigsize, &
       tempreduced(:,:,:,:),bigsize,DATAZERO,direct(:,:),numfrozen)

  frozenpotdiag=0
  do ispf=1,numfrozen
     do iispf=1,numfrozen
        frozenpotdiag=frozenpotdiag+2*direct(iispf,ispf)-exch(iispf,ispf)
     enddo
  enddo


  temppotmatel=0d0

  call mult_pot(numfrozen,infrozens(:,:,:,:),tempmult(:,:))
  if (numhatoms.gt.0) then
     call hatom_op(numfrozen,infrozens(:,:,:,:),tempmult2(:,:),hatomreduced)
     tempmult(:,:)=tempmult(:,:)+tempmult2(:,:)
  endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(spf2a)
!$OMP DO SCHEDULE(DYNAMIC)
  do spf2a=1,numfrozen
     temppotmatel(spf2a)=mydot(infrozens(:,:,:,spf2a),tempmult(:,spf2a),sizespf)
  enddo
!$OMP END DO
!$OMP END PARALLEL

  do i=1,numfrozen
     frozenpotdiag=frozenpotdiag+2*temppotmatel(i)
  enddo
  if (bornopflag/=1) then
     print *, "nuclear not done for frozen.";     stop
  endif
  temppotmatel=0d0

  call mult_ke(infrozens(:,:,:,:),tempmult(:,:),numfrozen)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(spf2a)
!$OMP DO SCHEDULE(DYNAMIC)
  do spf2a=1,numfrozen
     temppotmatel(spf2a)=mydot(infrozens(:,:,:,spf2a),tempmult(:,spf2a),sizespf)
  enddo
!$OMP END DO
!$OMP END PARALLEL

  frozenkediag=0d0
  do i=1,numfrozen
     frozenkediag=frozenkediag+2*temppotmatel(i)
  enddo

  deallocate(frodensity,cfrodensity,tempreduced,tempmult,tempmult2)

end subroutine call_frozen_matels_core


!! EXCHANGE (direct is op_frozenreduced)

subroutine op_frozen_exchange0(howmany,inspfs,outspfs,infrozens,numfrozen,inspfmvals)
  use myparams
  use tinvsubmod
  implicit none
  integer,intent(in) :: numfrozen,howmany,inspfmvals(howmany)
  DATATYPE,intent(in) :: infrozens(numerad,lbig+1,-mbig:mbig,numfrozen), &
       inspfs(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE,intent(out) :: outspfs(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE, allocatable :: tempreduced(:,:,:,:),twoeden(:,:,:,:)
  DATATYPE :: twoemat2(numerad,lbig+1,-2*mbig:2*mbig), myspf(numerad,lbig+1,-mbig:mbig)   !! AUTOMATIC
  integer :: mvalue2a, mvalue2b, spf2a,spf2b, deltam, qq,rr,qq2,rr2

  outspfs(:,:,:,:)=0d0

  if (numfrozen.eq.0) then
     return
  endif

  allocate(tempreduced(numerad,lbig+1,-2*mbig:2*mbig,howmany),&
       twoeden(numerad,lbig+1,-2*mbig:2*mbig,howmany))
  tempreduced(:,:,:,:)=0d0; twoeden(:,:,:,:)=0d0

  do spf2b=1,numfrozen
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(spf2a,qq,qq2,rr,rr2,twoemat2,mvalue2a,mvalue2b,deltam)
!$OMP DO SCHEDULE(DYNAMIC)
     do spf2a=1,howmany
        ! integrating over electron 2

        if (spfrestrictflag.eq.1) then
           qq=inspfmvals(spf2a)
           rr=inspfmvals(spf2a)
        else
           qq=-mbig
           rr=mbig
        endif

        qq2=-mbig
        rr2=mbig

        twoemat2(:,:,:)=0d0
        do mvalue2a=qq,rr
           do mvalue2b=qq2,rr2

!! deltam is frozen m minus spf m
              deltam=mvalue2b-mvalue2a

!! THIS WAY (CONJG FROZEN, NO CONJG INSPFS)
              twoemat2(:,:,deltam) = twoemat2(:,:,deltam) + &
                   inspfs(:,:,mvalue2a,spf2a) * CONJUGATE(infrozens(:,:,mvalue2b,spf2b))
           enddo
        enddo
        twoeden(:,:,:,spf2a)=twoemat2(:,:,:)
     enddo   !! spf2a
!$OMP END DO
!$OMP END PARALLEL

     call op_tinv(-2*mbig,2*mbig,howmany,twoeden,tempreduced(:,:,:,:))

!! deltam (twoemat2 index) is frozen m (not conjugated above) minus spf m (conjugated above)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(spf2a,myspf,twoemat2,mvalue2a,mvalue2b,qq,qq2,rr,rr2)
!$OMP DO SCHEDULE(DYNAMIC)
     do spf2a=1,howmany

        twoemat2(:,:,:)=tempreduced(:,:,:,spf2a)

        if (spfrestrictflag.eq.1) then
           qq=inspfmvals(spf2a)
           rr=inspfmvals(spf2a)
        else
           qq=-mbig
           rr=mbig
        endif

        qq2=-mbig
        rr2=mbig

        myspf(:,:,:)=0d0
        do mvalue2a=qq,rr   !! spf
           do mvalue2b=qq2,rr2  !! frozen

!! YES MINUS 1 TIMES              
              myspf(:,:,mvalue2a) = myspf(:,:,mvalue2a) -    & 
!!NO                   twoemat2(:,:,mvalue2b-mvalue2a) * CONJUGATE(infrozens(:,:,mvalue2b,spf2b))
!! THIS WAY
                   twoemat2(:,:,mvalue2b-mvalue2a) * (infrozens(:,:,mvalue2b,spf2b))
           enddo
        enddo
        outspfs(:,:,:,spf2a) = outspfs(:,:,:,spf2a) + myspf(:,:,:)
     enddo
!$OMP END DO
!$OMP END PARALLEL

  enddo

  deallocate(tempreduced,twoeden)
  
end subroutine op_frozen_exchange0


subroutine getdensity(density, indenmat, inspfs,howmany)
  use myparams
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: indenmat(howmany,howmany), inspfs(numerad,lbig+1,-mbig:mbig,howmany)
  complex*16,intent(out) :: density(numerad,lbig+1,2*mbig+1)
  integer :: i,j,ii,jj,kk, mm,nn
  real*8 :: phi,pi

  pi=4d0*atan(1d0)
  density=0.d0

  do mm=-mbig,mbig
  do nn=-mbig,mbig
     do i=1,howmany
        do j=1,howmany
           
           do kk =  1,2*mbig+1
              phi = 2*pi*kk/real(2*mbig+1,8)
              do jj=1,lbig+1
                 do ii=1,numerad

!!  integral sum_ij  | phi_j(x) > rho_ji < phi_i (x) |                    

!! 111510 was with transpose denmat                   density(ii,jj,kk) = density(ii,jj,kk) + indenmat(i,j) * inspfs(ii,jj,mm,i)*CONJUGATE(inspfs(ii,jj,nn,j)) * exp((0.d0,1.d0)*(nn-mm)*phi)

!! NOW 111510
                    density(ii,jj,kk) = density(ii,jj,kk) + indenmat(j,i) * inspfs(ii,jj,mm,i)*CONJUGATE(inspfs(ii,jj,nn,j)) * exp((0.d0,1.d0)*(nn-mm)*phi)

                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  enddo

end subroutine getdensity


subroutine op_conjg(in, out)
  use myparams
  implicit none
  DATATYPE,intent(in) :: in(numerad,lbig+1,-mbig:mbig)
  DATATYPE,intent(out) :: out(numerad,lbig+1,-mbig:mbig)
  integer :: i
  do i=-mbig,mbig
     out(:,:,i)=ALLCON(in(:,:,-i))
  enddo
end subroutine op_conjg


subroutine op_reflectz(in, out)
  use myparams
  implicit none
  DATATYPE,intent(in) :: in(numerad,lbig+1,-mbig:mbig)
  DATATYPE,intent(out) :: out(numerad,lbig+1,-mbig:mbig)
  integer :: i
  do i=1,lbig+1
     out(:,i,:)=in(:,lbig+2-i,:)
  enddo
end subroutine op_reflectz

subroutine op_reflecty(in, out)
  use myparams
  implicit none
  DATATYPE,intent(in) :: in(numerad,lbig+1,-mbig:mbig)
  DATATYPE,intent(out) :: out(numerad,lbig+1,-mbig:mbig)
  integer :: i
  do i=-mbig,mbig
     out(:,:,i)=in(:,:,-i)*(-1)**i
  enddo
end subroutine op_reflecty


subroutine op_reflectx(in, out)
  use myparams
  implicit none
  DATATYPE,intent(in) :: in(numerad,lbig+1,-mbig:mbig)
  DATATYPE,intent(out) :: out(numerad,lbig+1,-mbig:mbig)
  integer :: i
  do i=-mbig,mbig
     out(:,:,i)=in(:,:,-i)
  enddo
end subroutine op_reflectx


subroutine mult_zdipole(howmany,in, out, realflag)
  use myparams
  use myprojectmod
  implicit none
  integer,intent(in) :: realflag,howmany
  DATATYPE,intent(in) :: in(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE,intent(out) :: out(numerad,lbig+1,-mbig:mbig,howmany)
  integer :: i,ii

#ifndef CNORMFLAG
  if (realflag.ne.0) then
     do ii=1,howmany
        do i=-mbig,mbig
           out(:,:,i,ii)=in(:,:,i,ii)*real(zdipole(:,:),8)
        enddo
     enddo
  else
#endif
     do ii=1,howmany     
        do i=-mbig,mbig
           out(:,:,i,ii)=in(:,:,i,ii)*zdipole(:,:)
        enddo
     enddo
#ifndef CNORMFLAG
  endif
#endif

end subroutine mult_zdipole


subroutine mult_xdipole(howmany,in, out,realflag)
  use myparams
  use myprojectmod
  implicit none
  integer,intent(in) :: realflag,howmany
  DATATYPE,intent(in) :: in(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE,intent(out) :: out(numerad,lbig+1,-mbig:mbig,howmany)
  integer :: i,ii

  out(:,:,:,:)=0d0

#ifndef CNORMFLAG
  if (realflag.ne.0) then
     do ii=1,howmany
        do i=-mbig+1,mbig
           out(:,:,i,ii)=in(:,:,i-1,ii)*real(xydipole(:,:),8) /2 !!TWOFIX
        enddo
        do i=-mbig,mbig-1
           out(:,:,i,ii)=out(:,:,i,ii)+in(:,:,i+1,ii)*real(xydipole(:,:),8) /2  !!TWOFIX
        enddo
     enddo
  else
#endif
     do ii=1,howmany
        do i=-mbig+1,mbig
           out(:,:,i,ii)=in(:,:,i-1,ii)*xydipole(:,:) /2 !!TWOFIX
        enddo
        do i=-mbig,mbig-1
           out(:,:,i,ii)=out(:,:,i,ii)+in(:,:,i+1,ii)*xydipole(:,:) /2  !!TWOFIX
        enddo
     enddo
#ifndef CNORMFLAG
  endif
#endif

end subroutine mult_xdipole


subroutine mult_ydipole(howmany,in, out, realflag)
  use myparams
  use myprojectmod
  implicit none
  integer,intent(in) :: realflag,howmany
  DATATYPE,intent(in) :: in(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE,intent(out) :: out(numerad,lbig+1,-mbig:mbig,howmany)
  integer :: i,ii

  out(:,:,:,:)=0d0

#ifndef CNORMFLAG
  if (realflag.ne.0) then
     do ii=1,howmany
        do i=-mbig+1,mbig
           out(:,:,i,ii)=in(:,:,i-1,ii)*real(xydipole(:,:),8)*(0d0,1d0)  /2 !!TWOFIX
        enddo
        do i=-mbig,mbig-1
           out(:,:,i,ii)=out(:,:,i,ii)+in(:,:,i+1,ii)*real(xydipole(:,:),8)*(0d0,-1d0)  /2 !!TWOFIX
        enddo
     enddo
  else
#endif
     do ii=1,howmany
        do i=-mbig+1,mbig
           out(:,:,i,ii)=in(:,:,i-1,ii)*xydipole(:,:)*(0d0,1d0)  /2 !!TWOFIX
        enddo
        do i=-mbig,mbig-1
           out(:,:,i,ii)=out(:,:,i,ii)+in(:,:,i+1,ii)*xydipole(:,:)*(0d0,-1d0)  /2 !!TWOFIX
        enddo
     enddo
#ifndef CNORMFLAG
  endif
#endif

end subroutine mult_ydipole



!! OUTPUTS orbitals from firstspf to lastspf (out of numspf)

subroutine mult_reducedpot(firstspf,lastspf,inspfs,outspfs,reducedpot)
  use myparams
  implicit none
  integer,intent(in) :: firstspf,lastspf
  DATATYPE,intent(out) :: outspfs(numerad,lbig+1, -mbig:mbig,firstspf:lastspf)
  DATATYPE, intent(in) :: inspfs(numerad,lbig+1, -mbig:mbig, numspf),&
       reducedpot(numerad,lbig+1,-2*mbig:2*mbig,  numspf,firstspf:lastspf)
  DATATYPE :: mine(numerad,lbig+1)      !! AUTOMATIC
  integer :: ispf,imval,flag,kspf,kmval

  outspfs(:,:,:,:)=0d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ispf,imval,flag,kspf,kmval,mine)
!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
  do ispf=firstspf,lastspf
     do imval=-mbig,mbig
        flag=1
        if ((spfrestrictflag.eq.1).and.(imval.ne.spfmvals(ispf))) then
           flag=0
        endif
        if (flag==1) then
           mine(:,:)=0d0
           do kspf=1,numspf
              do kmval=-mbig,mbig
!! reducedpot is usual notation: <ispf | kspf> so so sum over slow index kspf
                    
                 mine(:,:) = mine(:,:) + & 
!! switching 2-2016         reducedpot(:,:,imval-kmval,ispf,kspf) * inspfs(:,:,kmval,kspf)
                      reducedpot(:,:,kmval-imval,kspf,ispf) * inspfs(:,:,kmval,kspf)
              enddo
           enddo
           outspfs(:,:,imval,ispf) = mine(:,:)
        endif
     enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

end subroutine mult_reducedpot


!! FOR EXPERIMENTAL ADDITION OF ADDITIONAL HYDROGEN ATOMS VIA POISSON SOLVE

subroutine hatomcalc(hatomreduced)
  use myparams
  use myprojectmod
  use tinvsubmod
  implicit none
  DATATYPE :: interpolate
  DATATYPE,intent(out) :: hatomreduced(numerad,lbig+1,-2*mbig:2*mbig)
  DATATYPE ::  hatomden(numerad,lbig+1,-2*mbig:2*mbig)   !! AUTOMATIC
  integer :: mvalue2a, mvalue2b, iatom, deltam, ixi,ieta

  if (numhatoms.eq.0) then
     return
  endif
  if (spfrestrictflag==1) then
     OFLWR "Hey, don't use spfrestrictflag if you have hatoms!";CFLST
  endif
  if (numr.gt.1) then
     OFLWR "Hatom not for calc with multiple r gridpoints.sorry."; CFLST
  endif
  do iatom=1,numhatoms
     if (hlocrealflag.ne.0) then
        OFLWR "HATOM AT r=", hlocreal(1,iatom), " theta=", hlocreal(2,iatom); CFL
     else
        OFLWR "REINSTATE " ; CFLST
   !!print *, "HATOM AT ", radialpoints(hlocs(1,iatom)), thetapoints(hlocs(2,iatom))
     endif
  enddo

  hatomden(:,:,:)=0d0

  do mvalue2a=-mbig,mbig
     do mvalue2b=-mbig,mbig
        deltam=mvalue2b-mvalue2a
        do iatom=1,numhatoms
           if (hlocrealflag.ne.0) then
              do ieta=1,lbig+1
                 do ixi=1,numerad
                    hatomden(ixi,ieta,deltam) = hatomden(ixi,ieta,deltam) + &
                         1.d0/(real(2*mbig+1,8))*hlocs(3,iatom)**(mvalue2a+mvalue2b) * &
                         interpolate(hlocreal(1,iatom),hlocreal(2,iatom),real(rpoints(1),8), abs(deltam),ixi,ieta)
                 enddo
              enddo
           else
              hatomden(hlocs(1,iatom),hlocs(2,iatom),deltam) = &
                   hatomden(hlocs(1,iatom),hlocs(2,iatom),deltam) + &
                   1.d0/(real(2*mbig+1,8))*hlocs(3,iatom)**(mvalue2a+mvalue2b)
           endif
        enddo
     enddo
  enddo

  hatomreduced=0.d0
  call op_tinv(-2*mbig,2*mbig,1,hatomden,hatomreduced)

end subroutine hatomcalc


!! DIRECT ONLY

subroutine op_frozenreduced(howmany,inspfs,outspfs,frozenreduced)
  use myparams
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: inspfs(numerad,lbig+1,-mbig:mbig,howmany), &
       frozenreduced(numerad,lbig+1,-2*mbig:2*mbig)
  DATATYPE,intent(out) :: outspfs(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE :: mine(numerad,lbig+1)
  integer :: kmval,imval,ii

  outspfs(:,:,:,:)=0d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,imval,kmval,mine) 
!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
  do ii=1,howmany
     do imval=-mbig,mbig
        mine(:,:)=0d0
        do kmval=-mbig,mbig
           mine(:,:) = mine(:,:) + & 
                frozenreduced(:,:,imval-kmval) * inspfs(:,:,kmval,ii)
        enddo
        outspfs(:,:,imval,ii)=mine(:,:)
     enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

end subroutine op_frozenreduced


subroutine restrict_spfs(inspfs,howmany,in_spfmvals)
  use myparams
  implicit none
  integer,intent(in) :: howmany,in_spfmvals(howmany)
  DATATYPE,intent(inout) :: inspfs(numerad,lbig+1,-mbig:mbig,howmany)
  call restrict_spfs0(inspfs,howmany,in_spfmvals,1)
end subroutine restrict_spfs


module hermdotmod
contains
  function hermdot(one,two,n)
    implicit none
    integer,intent(in) :: n
    DATATYPE,intent(in) :: one(n), two(n)
    DATATYPE :: hermdot
    hermdot=DOT_PRODUCT(one,two)
  end function hermdot
end module hermdotmod


subroutine restrict_spfs0(inspfs,howmany,in_spfmvals,printflag)
!! using hermdot; want to see if wfn has been chopped.
  use hermdotmod
  use myparams
  implicit none
  integer,intent(in) :: howmany,in_spfmvals(howmany),printflag
  DATATYPE,intent(inout) :: inspfs(numerad,lbig+1,-mbig:mbig,howmany)  
  integer :: ispf
  real*8 :: normsq1, normsq2

  normsq1=real(hermdot(inspfs,inspfs,howmany*edim*(2*mbig+1)),8)  !! ok hermdot
  do ispf=1,howmany
     inspfs(:,:,-mbig:in_spfmvals(ispf)-1,ispf)=0.d0
     inspfs(:,:,in_spfmvals(ispf)+1:mbig,ispf)=0.d0
  enddo
if (printflag.ne.0) then
  normsq2=real(hermdot(inspfs,inspfs,howmany*edim*(2*mbig+1)),8)  !! ok hermdot
  if (abs(normsq1-normsq2)/howmany.gt.1.d-7) then
     OFLWR "   WARNING, in restrict_spfs I lost norm.", normsq1,normsq2;  CFL
  endif
endif
end subroutine restrict_spfs0


subroutine ugrestrict_spfs(inspfs,howmany,in_spfugvals)
  use myparams
  implicit none
  integer,intent(in) :: howmany,in_spfugvals(howmany)
  DATATYPE,intent(inout) :: inspfs(numerad,lbig+1,-mbig:mbig,howmany)  
  call ugrestrict_spfs0(inspfs,howmany,in_spfugvals,1)
end subroutine ugrestrict_spfs

subroutine ugrestrict_spfs0(inspfs,howmany,in_spfugvals,printflag)
!! using hermdot; want to see if wfn has been chopped.
  use hermdotmod
  use myparams
  implicit none
  integer,intent(in) :: howmany,in_spfugvals(howmany),printflag
  DATATYPE,intent(inout) :: inspfs(numerad,lbig+1,-mbig:mbig,howmany)  
  integer :: ispf
  real*8 :: normsq1,normsq2
  DATATYPE :: outspfs(numerad,lbig+1,-mbig:mbig,howmany)   !! AUTOMATIC

  normsq1=real(hermdot(inspfs,inspfs,howmany*edim*(2*mbig+1)),8)  !! ok hermdot
  do ispf=1,howmany
     call ugrestrict(inspfs(:,:,:,ispf),outspfs(:,:,:,ispf), in_spfugvals(ispf))     
  enddo
  inspfs=outspfs
  if (printflag.ne.0) then
     normsq2=real(hermdot(inspfs,inspfs,howmany*edim*(2*mbig+1)),8)  !! ok hermdot
     if (abs(normsq1-normsq2)/howmany/max(normsq1,1d0).gt.1.d-3) then
        OFLWR "   WARNING, in UG restrict_spfs I lost norm.", normsq1,normsq2;  CFL
     endif
  endif

end subroutine ugrestrict_spfs0



function getsmallugvalue(inspf,inmval)
  use myparams
  implicit none
  integer,intent(in) :: inmval
  DATAECS,intent(in) :: inspf(numerad,lbig+1)
  DATAECS :: mytempspf(numerad,lbig+1)   !! AUTOMATIC
  integer :: k,getsmallugvalue

!! g or u: reflection then phi <-> -phi for these m eigenfuncts    
!!  yes I know (-1)^-1 = (-1)^1 who cares its an absolute value

  do k=1,lbig+1
     mytempspf(:,lbig+2-k)=inspf(:,k)*(-1)**abs(inmval)   
  enddo
  getsmallugvalue=nint(real(ecsdot(mytempspf,inspf,edim),4))

contains
  function ecsdot(one,two,n)
    implicit none
    integer,intent(in) :: n
    DATAECS,intent(in) :: one(n), two(n)
    DATAECS :: ecsdot, sum
    integer :: i
    sum=0.d0
    do i=1,n
       sum = sum + one(i) * two(i) 
    enddo
    ecsdot=sum
  end function ecsdot

end function getsmallugvalue


subroutine ugrestrict(inspfs,outspfs,ugval)
  use myparams
  implicit none
  integer,intent(in) :: ugval
  DATATYPE,intent(in) :: inspfs(numerad,lbig+1,-mbig:mbig)
  DATATYPE,intent(out) :: outspfs(numerad,lbig+1,-mbig:mbig)
  integer :: k, inmval

  if (abs(ugval).ne.1) then
     OFLWR " Error, ugrestrictflag when ugval is not 1 or -1 : ", ugval; CFLST
  endif

  outspfs=0d0
  do inmval=-mbig,mbig
     do k=1,lbig+1
        outspfs(:,lbig+2-k,inmval)=inspfs(:,k,inmval)*(-1)**abs(inmval)   !! g or u: reflection then phi <-> -phi for these m eigenfuncts   
     enddo
  enddo
!!  outspfs=0.5d0 * (outspfs + ugval*inspfs)

  outspfs=0.5d0 * (ugval*outspfs + inspfs)

end subroutine ugrestrict
  

subroutine bothcompact_spfs(inspfs,outspfs,howmany,in_spfmvals,in_spfugvals)
  use myparams
  implicit none
  integer,intent(in) :: howmany,in_spfmvals(howmany),in_spfugvals(howmany)
  DATATYPE,intent(in) :: inspfs(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE,intent(out) :: outspfs(numerad,(lbig+2)/2,howmany)  !! lbig even ok
  DATATYPE :: midspfs(numerad,lbig+1,howmany)

  call mcompact_spfs(inspfs,midspfs,howmany,in_spfmvals)
  call ugcompact_spfs(midspfs,outspfs,howmany,in_spfmvals,in_spfugvals)

end subroutine bothcompact_spfs

subroutine bothexpand_spfs(inspfs,outspfs,howmany,in_spfmvals,in_spfugvals)
  use myparams
  implicit none
  integer,intent(in) :: howmany,in_spfmvals(howmany),in_spfugvals(howmany)
  DATATYPE,intent(in) :: inspfs(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE,intent(out) :: outspfs(numerad,(lbig+2)/2,howmany)   !! lbig even ok
  DATATYPE :: midspfs(numerad,lbig+1,howmany)

  call ugexpand_spfs(inspfs,midspfs,howmany,in_spfmvals,in_spfugvals)
  call mexpand_spfs(midspfs,outspfs,howmany,in_spfmvals)

end subroutine bothexpand_spfs


subroutine mcompact_spfs(inspfs,outspfs,howmany,in_spfmvals)
  use myparams
  implicit none
  integer,intent(in) :: howmany,in_spfmvals(howmany)
  DATATYPE,intent(in) :: inspfs(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE,intent(out) :: outspfs(numerad,lbig+1,howmany)
  integer :: ispf
  do ispf=1,howmany
     outspfs(:,:,ispf)=inspfs(:,:,in_spfmvals(ispf),ispf)
  enddo
end subroutine mcompact_spfs

subroutine mexpand_spfs(inspfs,outspfs,howmany,in_spfmvals)
  use myparams
  implicit none
  integer,intent(in) :: howmany,in_spfmvals(howmany)
  DATATYPE,intent(in) :: inspfs(numerad,lbig+1,howmany)
  DATATYPE,intent(out) :: outspfs(numerad,lbig+1,-mbig:mbig,howmany)
  integer :: ispf

  outspfs(:,:,:,:)=0d0

  do ispf=1,howmany
     outspfs(:,:,in_spfmvals(ispf),ispf)=inspfs(:,:,ispf)
  enddo
end subroutine mexpand_spfs

!! g or u: reflection then phi <-> -phi for these m eigenfuncts   

subroutine ugcompact_spfs(inspfs,outspfs,howmany,in_spfmvals,in_spfugvals)
  use myparams
  implicit none
  integer,intent(in) :: howmany,in_spfugvals(howmany),in_spfmvals(howmany)
  DATATYPE,intent(in) :: inspfs(numerad,lbig+1,howmany)
  DATATYPE,intent(out) :: outspfs(numerad,(lbig+2)/2,howmany)   !! lbig even ok
  integer :: ispf,k

!!$  if (mod(lbig+1,2).ne.0) then
!!$     OFLWR "GGG))))iiiiiERROROAAAAUGHGGHG"; CFLST
!!$  endif

  outspfs=0d0
  do ispf=1,howmany
     do k=1,(lbig+2)/2
        outspfs(:,k,ispf)=inspfs(:,lbig+2-k,ispf)*(-1)**abs(in_spfmvals(ispf)) * in_spfugvals(ispf)
     enddo
  enddo
  outspfs(:,:,:)= sqrt(0.5d0) * (outspfs(:,:,:) + inspfs(:,1:(lbig+2)/2,:))

end subroutine ugcompact_spfs


subroutine ugexpand_spfs(inspfs,outspfs,howmany,in_spfmvals,in_spfugvals)
  use myparams
  implicit none
  integer,intent(in) :: howmany,in_spfugvals(howmany),in_spfmvals(howmany)
  DATATYPE,intent(out) :: outspfs(numerad,lbig+1,howmany)
  DATATYPE,intent(in) :: inspfs(numerad,(lbig+2)/2,howmany)    !! lbig even ok
  integer :: ispf,k

!!$  if (mod(lbig+1,2).ne.0) then
!!$     OFLWR "GGG))))iiiiiERROROAAAAUGHGGHG"; CFLST
!!$  endif

  outspfs(:,:,:)=0d0
  outspfs(:,1:(lbig+2)/2,:)=inspfs(:,:,:)

  do ispf=1,howmany
     do k=1,(lbig+2)/2
        outspfs(:,lbig+2-k,ispf)=inspfs(:,k,ispf)*(-1)**abs(in_spfmvals(ispf)) * in_spfugvals(ispf)
     enddo
  enddo
  outspfs(:,:,:)= sqrt(0.5d0) * outspfs(:,:,:) 

end subroutine ugexpand_spfs


module lobstuffmod
contains

function lobatto(numpoints,points2d,n,x)
  implicit none
  integer :: n,i,numpoints
  real*8 :: points2d(numpoints), x,product,lobatto
  product = 1.0d0
  do i=1,numpoints
     if (i/=n) product = product*(x-points2d(i))/(points2d(n)-points2d(i))
  enddo
  lobatto=product
end function lobatto

function  etalobatto(n,x, mvalue, etapoints, etaweights) !! ok unused
  use myparams
  implicit none
  integer :: mvalue,   n
  real*8 :: x,etalobatto, etapoints(lbig+1), etaweights(lbig+1)

  if (x .lt. -1.d0 .or. x .gt. 1.d0) then
     etalobatto=0.0d0;     return
  endif
  etalobatto = lobatto(lbig+1,etapoints,n,x)/sqrt(etaweights(n))  !! /sqrt(thetaweights(n)) put in after
           
! for both xi and eta

  etalobatto=etalobatto*       sqrt(  (x**2 - 1.d0)  /  (etapoints(n)**2 - 1.d0)  )**mod(mvalue+100,2)
  
end function etalobatto


function  etalobattoint(n,x, mvalue, etapoints, etaweights)  !! ok unused
  use myparams
  implicit none
  integer :: mvalue,   n
  real*8 :: x,etalobattoint, etapoints(lbig+1), etaweights(lbig+1)
  if (x .lt. -1.d0 .or. x .gt. 1.d0) then
     etalobattoint=0.0d0;     return
  endif
  etalobattoint = lobatto(lbig+1,etapoints,n,x)
  etalobattoint=etalobattoint*        sqrt(  (x**2 - 1.d0)  /  (etapoints(n)**2 - 1.d0)  )**mod(mvalue+100,2)
end function etalobattoint


function  xilobatto(n,x, mvalue, xinumpoints, firstelpts, secondelpts, &   !! ok unused
     xipoints, xiweights, xielementsizes,xinumelements,xiflag)
  implicit none
  integer :: mvalue, l,el2d, n,num,whichelement,point2d,  xinumpoints, xiflag,xinumelements
  real*8 :: x,x2d,lastboundary,xielementsizes(*), firstelpts(*), secondelpts(*)
  DATAECS :: xilobatto, xiweights(*), xipoints(*)

  xilobatto=0.0d0 
  if (x .le. 1) then
     xilobatto=0.0d0;     return
  endif
  lastboundary = 1d0
  do whichelement=1,xinumelements
     lastboundary = lastboundary + xielementsizes(whichelement)
     if (lastboundary.ge.x) then
        exit
     endif
  enddo
  if (whichelement.eq.xinumelements+1) then
     xilobatto=0d0; return
  endif
  lastboundary = lastboundary - xielementsizes(whichelement)
  num=1
  if (mod(n-1,xinumpoints-1).eq.0) then
     num=2
  endif

  l=mod(n-1,xinumpoints-1) +1;  el2d=(n-l)/(xinumpoints-1)+1
  point2d = n - (el2d-1)*(xinumpoints-1)
  
  if (.not.(el2d.eq.whichelement)) then
     if (num.eq.1) then
        xilobatto=0.0d0
        return
     else
        el2d = el2d -1;        point2d = point2d + xinumpoints-1
        if (.not.(el2d.eq.whichelement)) then
           xilobatto=0.0d0;           return
        endif
     endif
  endif

  x2d = ( (x - lastboundary) / xielementsizes(whichelement) ) * 2.0d0 - 1.0d0

  if (whichelement.eq.1) then
     xilobatto = lobatto(xinumpoints,firstelpts,point2d,x2d)/sqrt(xiweights(n))
  else
     xilobatto = lobatto(xinumpoints,secondelpts,point2d,x2d)/sqrt(xiweights(n))
  endif

! for both xi and eta

  if (xiflag.eq.1) then
     xilobatto=xilobatto* &
          sqrt(  (x**2 - 1.d0)  /  (xipoints(n)**2 - 1.d0)  )**mod(mvalue+100,2)
  endif

end function xilobatto


function  xilobattoint(n,x, mvalue, xinumpoints, firstelpts, secondelpts, xipoints,& 
     xiweights, xielementsizes,xinumelements,xiflag) !! ok unused

  implicit none

  integer :: mvalue, l,el2d, n,num,whichelement,point2d,  xinumpoints, xiflag,xinumelements
  real*8 :: xielementsizes(*), firstelpts(*), secondelpts(*),x,x2d,lastboundary
  DATAECS :: xilobattoint, xiweights(*), xipoints(*)

  xilobattoint=0.0d0 
  if (x .le. 1d0) then
     xilobattoint=0.0d0;     return
  endif
  lastboundary = 1d0
  do whichelement=1,xinumelements
     lastboundary = lastboundary + xielementsizes(whichelement)
     if (lastboundary.ge.x) then
        exit
     endif
  enddo
  if (whichelement.eq.xinumelements+1) then
     xilobattoint=0d0; return
  endif
  lastboundary = lastboundary - xielementsizes(whichelement)
  num=1
  if (mod(n-1,xinumpoints-1).eq.0) then
     num=2
  endif

  l=mod(n-1,xinumpoints-1) +1;  el2d=(n-l)/(xinumpoints-1)+1
  point2d = n - (el2d-1)*(xinumpoints-1)
  
  if (.not.(el2d.eq.whichelement)) then
     if (num.eq.1) then
        xilobattoint=0.0d0
        return
     else
        el2d = el2d -1
        point2d = point2d + xinumpoints-1
        if (.not.(el2d.eq.whichelement)) then
           xilobattoint=0.0d0
           return
        endif
     endif
  endif

  x2d = ( (x - lastboundary) / xielementsizes(whichelement) ) * 2.0d0 - 1.0d0

  if (whichelement.eq.1) then
     xilobattoint = lobatto(xinumpoints,firstelpts,point2d,x2d)
  else
     xilobattoint = lobatto(xinumpoints,secondelpts,point2d,x2d)
  endif

  if (xiflag.eq.1) then
     xilobattoint=xilobattoint* &
          sqrt(  (x**2 - 1.d0)  /  (xipoints(n)**2 - 1.d0)  )**mod(mvalue+100,2)
  endif

end function xilobattoint

end module lobstuffmod


!! all proddrhos are the same. (and asymmetric individually)

subroutine velmultiply(howmany,spfin,spfout, myxtdpot0,myytdpot0,myztdpot)
  use myparams
  use myprojectmod
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: spfin(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE,intent(out) :: spfout(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE,intent(in) :: myxtdpot0,myztdpot,myytdpot0
  integer :: imval, qq, ieta , ixi, i ,ivect
  complex*16 :: csum1,csum2,cfacreal,cfacimag
  real*8 :: myrhotdpotreal,myrhotdpotimag
  DATATYPE :: work(lbig+1),work2(lbig+1)

#ifdef REALGO
OFLWR "Velocity gauge not available for real time propagation"; CFLST
#endif

  spfout(:,:,:,:)=0d0

  if (abs(myztdpot).gt.1.d-10) then
     do ivect=1,howmany
        do imval=-mbig,mbig
           qq=lbig+1
           do ixi=1,numerad
              work2(:)=spfin(ixi,:,imval,ivect)
              call MYGEMV('N',qq,qq,DATAONE,sparseddz_eta(:,:,ixi,abs(imval)+1),qq,work2,1,DATAZERO, work, 1)
              spfout(ixi,:,imval,ivect)= spfout(ixi,:,imval,ivect) + work(1:lbig+1) * myztdpot * (0.d0,1.d0)  
           enddo

           i=2*bandwidth+1
           csum1=myztdpot * (0.d0,1.d0) 
           do ieta=1,lbig+1
              call MYGBMV('N',numerad,numerad,bandwidth,bandwidth,csum1,sparseddz_xi_banded(:,:,ieta,abs(imval)+1),i,&
                   spfin(:,ieta,imval,ivect),1,DATAONE, spfout(:,ieta,imval,ivect), 1)
           enddo
           spfout(:,:,imval,ivect) = spfout(:,:,imval,ivect) - &
                sparseddz_diag(:,:,abs(imval)+1) * spfin(:,:,imval,ivect) * myztdpot * (0.0d0, 1.d0) 
        enddo
     enddo
  endif

  cfacreal=0d0;  myrhotdpotreal=sqrt(real(myytdpot0,8)**2+real(myxtdpot0,8)**2)
  if (myrhotdpotreal.ne.0d0) then
     cfacreal=exp((0d0,-1d0)*atan2(real(myytdpot0,8),real(myxtdpot0,8)))     
  endif
  cfacimag=0d0;  myrhotdpotimag=sqrt(imag((0d0,0d0)+myytdpot0)**2+imag((0d0,0d0)+myxtdpot0)**2)
  if (myrhotdpotimag.ne.0d0) then
     cfacimag=exp((0d0,-1d0)*atan2(imag((0d0,0d0)+myytdpot0),imag((0d0,0d0)+myxtdpot0)))     
  endif

  if (abs(myrhotdpotreal**2+myrhotdpotimag**2).gt.0d0) then

     csum1= (0.d0,1.d0) * ( myrhotdpotreal *cfacreal + (0d0,1d0) * myrhotdpotimag *cfacimag )                     /2 !!TWOFIX
     csum2= (0.d0,1.d0) * ( myrhotdpotreal * CONJG(cfacreal) + (0d0,1d0) * myrhotdpotimag * CONJG(cfacimag))      /2 !!TWOFIX

     do ivect=1,howmany
        do imval = -mbig,mbig
           if (mod(imval,2).eq.0) then  !! even.  no transpose.
              qq=lbig+1
              do ixi=1,numerad
                 work2(:)=spfin(ixi,:,imval,ivect)
                 call MYGEMV('N',qq,qq,DATAONE,sparseddrho_eta(:,:,ixi,1),qq,work2,1,DATAZERO, work, 1)
                 if (imval.gt.-mbig) then
                    spfout(ixi,:,imval-1,ivect)= spfout(ixi,:,imval-1,ivect) + work(1:lbig+1) * csum1
                 endif
                 if (imval .lt. mbig) then
                    spfout(ixi,:,imval+1,ivect)= spfout(ixi,:,imval+1,ivect) + work(1:lbig+1) * csum2
                 endif
              enddo
              i=2*bandwidth+1
              do ieta=1,lbig+1
                 if (imval.gt.-mbig) then
                    call MYGBMV('N',numerad,numerad,bandwidth,bandwidth,csum1,sparseddrho_xi_banded(:,:,ieta,1),i,&
                         spfin(:,ieta,imval,ivect),1,DATAONE, spfout(:,ieta,imval-1,ivect), 1)
                 endif

                 if (imval.lt.mbig) then
                    call MYGBMV('N',numerad,numerad,bandwidth,bandwidth,csum2,sparseddrho_xi_banded(:,:,ieta,1),i,&
                         spfin(:,ieta,imval,ivect),1,DATAONE, spfout(:,ieta,imval+1,ivect), 1)
                 endif
              enddo
              if (imval.gt.-mbig) then                       
                 spfout(:,:,imval-1,ivect) = spfout(:,:,imval-1,ivect) - sparseddrho_diag(:,:,1) * spfin(:,:,imval,ivect) *csum1
              endif
              if (imval.lt.mbig) then                       
                 spfout(:,:,imval+1,ivect) = spfout(:,:,imval+1,ivect) - sparseddrho_diag(:,:,1) * spfin(:,:,imval,ivect) *csum2
              endif

           !! lowering
              if (imval.gt.-mbig) then
                 spfout(:,:,imval-1,ivect)= spfout(:,:,imval-1,ivect) + ddrhopot * (imval-1) * spfin(:,:,imval,ivect) * csum1
              endif
           
           !! raising
              if (imval.lt.mbig) then
                 spfout(:,:,imval+1,ivect)= spfout(:,:,imval+1,ivect) - ddrhopot * (imval+1) * spfin(:,:,imval,ivect) * csum2 
              endif
           
           else  !! even or odd:  odd. transpose.

              qq=lbig+1
              do ixi=1,numerad
                 work2(:)=spfin(ixi,:,imval,ivect)
                 call MYGEMV('T',qq,qq,DATANEGONE,sparseddrho_eta(:,:,ixi,1),qq,work2,1,DATAZERO, work, 1)
                 if (imval.gt.-mbig) then
                    spfout(ixi,:,imval-1,ivect)= spfout(ixi,:,imval-1,ivect) + work(1:lbig+1) * csum1
                 endif
                 if (imval .lt. mbig) then
                    spfout(ixi,:,imval+1,ivect)= spfout(ixi,:,imval+1,ivect) + work(1:lbig+1) * csum2
                 endif
              enddo

              i=2*bandwidth+1
              do ieta=1,lbig+1
                 if (imval.gt.-mbig) then
                    call MYGBMV('T',numerad,numerad,bandwidth,bandwidth,csum1*(-1),sparseddrho_xi_banded(:,:,ieta,1),i,&
                         spfin(:,ieta,imval,ivect),1,DATAONE, spfout(:,ieta,imval-1,ivect), 1)
                 endif

                 if (imval.lt.mbig) then
                    call MYGBMV('T',numerad,numerad,bandwidth,bandwidth,csum2*(-1),sparseddrho_xi_banded(:,:,ieta,1),i,&
                         spfin(:,ieta,imval,ivect),1,DATAONE, spfout(:,ieta,imval+1,ivect), 1)
                 endif
              enddo
              if (imval.gt.-mbig) then                       
                 spfout(:,:,imval-1,ivect) = spfout(:,:,imval-1,ivect) + sparseddrho_diag(:,:,1) * spfin(:,:,imval,ivect) * csum1
              endif
              if (imval.lt.mbig) then                       
                 spfout(:,:,imval+1,ivect) = spfout(:,:,imval+1,ivect) + sparseddrho_diag(:,:,1) * spfin(:,:,imval,ivect) * csum2
              endif

           !! lowering
              if (imval.gt.-mbig) then
                 spfout(:,:,imval-1,ivect)= spfout(:,:,imval-1,ivect) + ddrhopot * (imval) * spfin(:,:,imval,ivect) * csum1
              endif
           
           !! raising
              if (imval.lt.mbig) then
                 spfout(:,:,imval+1,ivect)= spfout(:,:,imval+1,ivect) - ddrhopot  * (imval) * spfin(:,:,imval,ivect) * csum2
              endif
           endif  ! mval even or odd
        enddo     ! imval
     enddo        ! ivect
  endif

end subroutine velmultiply



#ifndef REALGO
#define XXMVXX zgemv 
#define XXBBXX zgbmv 
#else
#define XXMVXX dgemv 
#define XXBBXX dgbmv 
#endif


!! needs factor of 1/r^2 for ham

subroutine mult_ke(in, out,howmany)
  use myparams
  use myprojectmod  
  implicit none
  integer,intent(in) :: howmany
  DATATYPE,intent(in) :: in(numerad,lbig+1,-mbig:mbig,howmany)
  DATATYPE,intent(out) :: out(numerad,lbig+1,-mbig:mbig,howmany)
  integer :: m2val, ixi,ieta, i,jj
  DATATYPE :: work(lbig+1), work2(lbig+1),myin(numerad,lbig+1),myout(numerad,lbig+1)

  out=0.d0

  do jj=1,howmany
     do m2val=-mbig,mbig
        myin(:,:)=in(:,:,m2val,jj)
        myout(:,:)=0d0
        do ixi=1,numerad
           work2=myin(ixi,:)
           call XXMVXX('N',lbig+1,lbig+1,DATAONE,sparseops_eta(:,:,ixi,abs(m2val)+1),lbig+1,work2,1,DATAZERO, work, 1)
           myout(ixi,:)=work(1:lbig+1)
        enddo
        i=2*bandwidth+1
        do ieta=1,lbig+1
           call XXBBXX('N',numerad,numerad,bandwidth,bandwidth,DATAONE,sparseops_xi_banded(:,:,ieta,abs(m2val)+1),i,myin(:,ieta),1,DATAONE, myout(:,ieta), 1)
        enddo
        out(:,:,m2val,jj) = myout(:,:) - sparseops_diag(:,:,abs(m2val)+1) * myin(:,:)
     enddo
  enddo

end subroutine mult_ke



subroutine reinterpolate_orbs_real(rspfs,dims,num)
  use myparams
  implicit none
  integer, intent(in) :: dims(3),num
  real*8 :: rspfs(dims(1),dims(2),dims(3),num)
  OFLWR "Reinterpolate orbs not supported atom/diatom"; CFLST
  rspfs(:,:,:,:)=0
end subroutine reinterpolate_orbs_real


subroutine reinterpolate_orbs_complex(cspfs,dims,num)
  use myparams
  implicit none
  integer, intent(in) :: dims(3),num
  complex*16 :: cspfs(dims(1),dims(2),dims(3),num)
  OFLWR "Reinterpolate orbs not supported atom/diatom"; CFLST
  cspfs(:,:,:,:)=0
end subroutine reinterpolate_orbs_complex

subroutine bigdimsub(localdims,alldims)
  use myparams
  implicit none
  integer, intent(in) :: localdims(3)
  integer, intent(out) :: alldims(3)
  OFLWR "What?  calling bigdimsub for atom or diatom not allowed"; CFLST
  alldims(:)=localdims(:)
end subroutine bigdimsub


subroutine splitscatterv(inbig,outlocal)
  use myparams
  implicit none
  integer :: inbig,outlocal
  OFLWR "What? don't call split gather/scatter routines for atom/diatom.  Programmer fail"; CFLST
  inbig=outlocal
end subroutine splitscatterv

subroutine splitscatterv_complex(inbig,outlocal)
  use myparams
  implicit none
  integer :: inbig,outlocal
  OFLWR "What? don't call split gather/scatter routines for atom/diatom.  Programmer fail"; CFLST
  inbig=outlocal
end subroutine splitscatterv_complex

subroutine splitscatterv_real(inbig,outlocal)
  use myparams
  implicit none
  integer :: inbig,outlocal
  OFLWR "What? don't call split gather/scatter routines for atom/diatom.  Programmer fail"; CFLST
  inbig=outlocal
end subroutine splitscatterv_real

subroutine splitgatherv(inbig,outlocal,ilog)
  use myparams
  implicit none
  logical :: ilog
  integer :: inbig,outlocal
  OFLWR "What? don't call split gather/scatter routines for atom/diatom.  Programmer fail"; CFLST
  inbig=outlocal
  ilog=.true.
end subroutine splitgatherv

subroutine splitgatherv_complex(inbig,outlocal,ilog)
  use myparams
  implicit none
  logical :: ilog
  integer :: inbig,outlocal
  OFLWR "What? don't call split gather/scatter routines for atom/diatom.  Programmer fail"; CFLST
  inbig=outlocal
  ilog=.true.
end subroutine splitgatherv_complex

subroutine splitgatherv_real(inbig,outlocal,ilog)
  use myparams
  implicit none
  logical :: ilog
  integer :: inbig,outlocal
  OFLWR "What? don't call split gather/scatter routines for atom/diatom.  Programmer fail"; CFLST
  inbig=outlocal
  ilog=.true.
end subroutine splitgatherv_real




!!$  
!!$  !! flag=1 means flux, otherwise whole op
!!$  subroutine call_flux_op_twoe00(lowspf,highspf,mobra,moket,V2,flag) 
!!$  !! determine the 2-electron matrix elements in the orbital basis for the flux operator i(H-H^{\dag}) 
!!$  !! input :
!!$  !! mobra - the orbitals that contained in the bra 
!!$  !! moket - the orbitals that contained in the ket 
!!$  !! output : 
!!$  !! V2 - the 2-electron matrix elements corresponding with potential energy (contract with 1/R) 
!!$    use myparams
!!$    use myprojectmod   !! rmatrix,ylmvals
!!$    implicit none
!!$    integer,intent(in) :: flag,lowspf,highspf
!!$    DATATYPE,intent(in) :: mobra(numerad,lbig+1,-mbig:mbig,numspf),moket(numerad,lbig+1,-mbig:mbig,numspf)
!!$    DATATYPE,intent(out) :: V2(numspf,numspf,numspf,lowspf:highspf)
!!$    DATATYPE :: twoemat2(numerad,lbig+1,-2*mbig:2*mbig),twoeden(numerad,lbig+1),&   !! AUTOMATIC
!!$         twoeden2(numerad),twoeden3(numerad)
!!$    integer :: i,a,j,b,mvali,mvala,mvalj,mvalb,ixi,ieta,deltam,lsum,qq,rr,qq2,rr2 
!!$  
!!$    if (highspf-lowspf+1.gt.0) then
!!$       V2=0d0
!!$    endif
!!$  
!!$  !! The bra determinant is <ij|
!!$  !! The ket determinant is |ab>
!!$  !! without the ecs rmatrix and ylm are real, sooooooo V-V^dagger is zero by def...
!!$  !! ALWAYS USING ALLLCON TO MAKE SURE ITS ALWAYS hermitian elements of V-V^\dagger and V\dagger must have transpose 
!!$  
!!$    do b=lowspf,highspf
!!$      do j=1,numspf !! bra vector
!!$  !! integrating over electron 2
!!$        twoemat2=0d0
!!$        if(spfrestrictflag==1) then
!!$          qq=spfmvals(j);        rr=spfmvals(j)
!!$          qq2=spfmvals(b);        rr2=spfmvals(b)
!!$        else
!!$          qq=-mbig;        rr=mbig
!!$          qq2=-mbig;        rr2=mbig
!!$        endif
!!$        do mvalj=qq,rr
!!$          do mvalb=qq2,rr2
!!$            deltam=mvalb-mvalj
!!$  
!!$  !! switching 2-2016     twoeden(:,:) = CONJUGATE(mobra(:,:,mvalj,j)) * moket(:,:,mvalb,b)
!!$  
!!$            twoeden(:,:) = mobra(:,:,mvalj,j) * CONJUGATE(moket(:,:,mvalb,b))
!!$  
!!$            do lsum=1+abs(deltam),jacobisummax+1
!!$  !! this is always real and serves as a delta function, can use for both V and V^\dagger
!!$  !! We do not need conjugation and separate V and V^\dagger densities as Ylm is real
!!$  !! Ylm's here mean it is a delta function for e-2 in eta 
!!$              twoeden2=0d0
!!$              do ieta=1,lbig+1
!!$                twoeden2(:)=twoeden2(:) + twoeden(:,ieta) * ylmvals(abs(deltam),ieta,lsum-abs(deltam)) 
!!$              enddo
!!$  !! we need to be able to do the conjugate here
!!$  !! thankfully the rmatrix is diagonal in xi, so transposing rmatrix only switches e-s which doesn't make any sense
!!$              twoeden3=0d0
!!$              select case(flag)
!!$              case(1)  ! flux imag
!!$  #ifdef ECSFLAG
!!$                 do ixi=1,numerad
!!$                    if(atomflag==0) then
!!$                       twoeden3(:)=twoeden3(:) + twoeden2(ixi) * 4d0 * (rmatrix(:,ixi,abs(deltam)+1,lsum-abs(deltam)) - ALLCON(rmatrix(:,ixi,abs(deltam)+1,lsum-abs(deltam))))
!!$                    else
!!$                       twoeden3(:)=twoeden3(:) + twoeden2(ixi) * (rmatrix(:,ixi,1,lsum) - ALLCON(rmatrix(:,ixi,1,lsum)))
!!$                    endif
!!$                 enddo
!!$  #endif
!!$              case(2)  ! flux real
!!$  #ifdef ECSFLAG
!!$                 do ixi=1,numerad
!!$                    if(atomflag==0) then
!!$                       twoeden3(:)=twoeden3(:) + twoeden2(ixi) * 4d0 * (rmatrix(:,ixi,abs(deltam)+1,lsum-abs(deltam)) + ALLCON(rmatrix(:,ixi,abs(deltam)+1,lsum-abs(deltam))))
!!$                    else
!!$                       twoeden3(:)=twoeden3(:) + twoeden2(ixi) * (rmatrix(:,ixi,1,lsum) + ALLCON(rmatrix(:,ixi,1,lsum)))
!!$                    endif
!!$                enddo
!!$  #endif
!!$                 case default
!!$                 do ixi=1,numerad
!!$                    if(atomflag==0) then
!!$                       twoeden3(:)=twoeden3(:) + twoeden2(ixi) * 4d0 * (rmatrix(:,ixi,abs(deltam)+1,lsum-abs(deltam)))
!!$                    else
!!$                       twoeden3(:)=twoeden3(:) + twoeden2(ixi) * (rmatrix(:,ixi,1,lsum))
!!$                    endif
!!$                 enddo
!!$              end select
!!$  
!!$  !! this is always real and serves as a delta function, can use for both V and V^\dagger
!!$  !! this uses the same Ylm realness as above to cleanly finish the matrix element
!!$              do ieta=1,lbig+1
!!$                twoemat2(:,ieta,deltam)=twoemat2(:,ieta,deltam) - twoeden3(:) * ylmvals(abs(deltam),ieta,lsum-abs(deltam))
!!$              enddo
!!$            enddo
!!$          enddo
!!$        enddo
!!$  !! slap on the other two orbitals, twoemat2 already holds a half transformed (for only e-2) V-V\dagger
!!$        do a=1,numspf !! ket vector
!!$          do i=1,numspf !!  bra vector
!!$            if(spfrestrictflag==1) then
!!$              qq=spfmvals(i);            rr=spfmvals(i)
!!$              qq2=spfmvals(a);            rr2=spfmvals(a)
!!$            else
!!$              qq=-mbig;            rr=mbig
!!$              qq2=-mbig;            rr2=mbig
!!$            endif
!!$            do mvali=qq,rr
!!$              do mvala=qq2,rr2
!!$                deltam=mvali-mvala
!!$                do ieta=1,lbig+1
!!$                  do ixi=1,numerad
!!$  !! saving the integral in Chemist's form
!!$  !! V(i,a,j,b) = (ia|jb)
!!$  !! i and j are in the bra, a and b are in the ket ie <ij|V-V\dagger|ab>
!!$                    !V2(i,a,j,b) = V2(i,a,j,b) + (0d0,1d0) * &
!!$                    !      CONJUGATE(mobra(ixi,ieta,mvali,i)) * moket(ixi,ieta,mvala,a) * twoemat2(ixi,ieta,deltam)
!!$  
!!$                    V2(i,a,j,b) = V2(i,a,j,b) +  &
!!$  !! switching 2-2016        CONJUGATE(mobra(ixi,ieta,mvali,i)) * moket(ixi,ieta,mvala,a) * twoemat2(ixi,ieta,deltam)
!!$                      mobra(ixi,ieta,mvali,i) * CONJUGATE(moket(ixi,ieta,mvala,a)) * twoemat2(ixi,ieta,deltam)
!!$  
!!$                  enddo
!!$                enddo
!!$              enddo
!!$            enddo
!!$          enddo
!!$        enddo
!!$      enddo
!!$    enddo
!!$  
!!$    if (highspf-lowspf+1.gt.0) then
!!$       if(flag.eq.1) then
!!$          V2=(0d0,1d0)*V2       !!! yes because I got the imaginary part above with conjugates...  returns imaginary... mult by i to get -2 x imag part (real valued number) like others
!!$       else if (flag.eq.2) then
!!$          V2=(-1)*V2
!!$       endif
!!$    endif
!!$  
!!$  
!!$  end subroutine call_flux_op_twoe00
!!$  
