
#include "Definitions.INC"



subroutine transferparams(innumspf,inspfrestrictflag,inspfmvals,inspfugrestrict,inspfugvals,outspfsmallsize,outorbparflag) ! ok unused
  use myparams
  implicit none
  integer :: innumspf,inspfrestrictflag,inspfmvals(innumspf), inspfugrestrict,inspfugvals(innumspf), outspfsmallsize
  logical, intent(out) :: outorbparflag
  numspf=innumspf;  
  outspfsmallsize=totpoints
  outorbparflag=orbparflag
end subroutine transferparams

subroutine twoedealloc()
end subroutine twoedealloc

subroutine call_flux_op_twoe() !mobra,moket,V2,flag) 
print *, "DOME flux_op_twoe"; stop
end subroutine call_flux_op_twoe

!!$ subroutine hatom_matel(inspfs1, inspfs2, hatommatel,numberspf)   !!!rmatrix,ylmvals, 

!!$ subroutine hatom_op(inspfs, outspfs)    !! , rmatrix,ylmvals

subroutine call_frozen_matels0() !infrozens,numfrozen,frozenkediag,frozenpotdiag)  !! returns last two.  a little cloogey
print *, "DOME CALL FROZEN MATELS."; stop
end subroutine call_frozen_matels0

subroutine call_frozen_exchange0() !inspfs,outspfs,infrozens,numfrozen)   !! rmatrix ylmvals
print *, "DOME CALL FROZEN EXCHANGE"; stop
end subroutine call_frozen_exchange0

subroutine getdensity() !density, indenmat, inspfs,numspf)
print *, "DOME GETDENSITY"; stop
end subroutine getdensity

subroutine mult_zdipole(in,out)
  use myparams
  use myprojectmod
  implicit none
  DATATYPE :: in(totpoints),out(totpoints)
  out(:)=in(:)*dipoles(:,3)
end subroutine mult_zdipole

subroutine mult_imzdipole(in, out)
  use myparams
  use myprojectmod
  implicit none
  DATATYPE :: in(totpoints),out(totpoints)
#ifdef REALGO
  out(:)=0
#else
  out(:)=in(:)*imag(dipoles(:,3))
#endif
end subroutine mult_imzdipole

subroutine mult_ydipole(in,out)
  use myparams
  use myprojectmod
  implicit none
  DATATYPE :: in(totpoints),out(totpoints)
  out(:)=in(:)*dipoles(:,2)
end subroutine mult_ydipole

subroutine mult_imydipole(in, out)
  use myparams
  use myprojectmod
  implicit none
  DATATYPE :: in(totpoints),out(totpoints)
#ifdef REALGO
  out(:)=0
#else
  out(:)=in(:)*imag(dipoles(:,2))
#endif
end subroutine mult_imydipole

subroutine mult_xdipole(in,out)
  use myparams
  use myprojectmod
  implicit none
  DATATYPE :: in(totpoints),out(totpoints)
  out(:)=in(:)*dipoles(:,1)
end subroutine mult_xdipole

subroutine mult_imxdipole() !in, out)
  use myparams
  use myprojectmod
  implicit none
  DATATYPE :: in(totpoints),out(totpoints)
#ifdef REALGO
  out(:)=0
#else
  out(:)=in(:)*imag(dipoles(:,1))
#endif
end subroutine mult_imxdipole



subroutine hatomcalc()
end subroutine hatomcalc

subroutine op_frozenreduced() !inspfs,outspfs)
print *, "DOME OP FROZENREDUCED"; stop
end subroutine op_frozenreduced

subroutine restrict_spfs() !inspfs,numspf,spfmvals)
end subroutine restrict_spfs

subroutine ugrestrict_spfs() !inspfs,numspf,spfmvals)
end subroutine ugrestrict_spfs

subroutine restrict_spfs0() !inspfs,numspf,spfmvals,printflag)
end subroutine restrict_spfs0

subroutine ugrestrict_spfs0() !inspfs,numspf,spfmvals,printflag)
end subroutine ugrestrict_spfs0

subroutine bothcompact_spfs() !inspfs,outspfs,numspf,spfmvals,spfugvals)
print *, "NOT APPLICABLE BOTHCOMPACT SINCDVR"
end subroutine bothcompact_spfs

subroutine bothexpand_spfs() !inspfs,outspfs,numspf,spfmvals,spfugvals)
print *, "NOT APPLICABLE BOTHEXPAND SINCDVR"
end subroutine bothexpand_spfs

subroutine mcompact_spfs() !inspfs,outspfs,numspf,spfmvals)
print *, "NOT APPLICABLE MCOMPACT SINCDVR"; stop
end subroutine mcompact_spfs

subroutine mexpand_spfs() !inspfs,outspfs,numspf,spfmvals)
print *, "NOT APPLICABLE MEXPAND SINCDVR"; stop
end subroutine mexpand_spfs


subroutine velmultiply(spfin,spfout, myxtdpot0,myytdpot0,myztdpot)
  use myparams
  implicit none
  DATATYPE :: spfin(totpoints),spfout(totpoints),myxtdpot0,myytdpot0,myztdpot,&
       work(totpoints)
  spfout(:)=0d0
  if (abs(myxtdpot0).gt.0d0) then
     call mult_xderiv(spfin,work,1)
     spfout(:)=spfout(:)+work(:)*(0d0,1d0)
  endif
  if (abs(myytdpot0).gt.0d0) then
     call mult_yderiv(spfin,work,1)
     spfout(:)=spfout(:)+work(:)*(0d0,1d0)
  endif
  if (abs(myztdpot).gt.0d0) then
     call mult_zderiv(spfin,work,1)
     spfout(:)=spfout(:)+work(:)*(0d0,1d0)
  endif

end subroutine velmultiply

subroutine imvelmultiply() !spfin,spfout, myxtdpot0,myytdpot0,myztdpot)
print *, "DOME imVELMULTIPLY SINCDVR"; stop
end subroutine imvelmultiply

subroutine mult_imke() !in, out)
print *, "DOME MULT IMKE"; stop
end subroutine mult_imke

subroutine mult_reke() !in, out)
print *, "DOME MULT reKE"; stop
end subroutine mult_reke


!!! from circ_sub_mpi, fttimes:
!!! times(6) = circ math
!!! from myzfft3d_par:
!!! times(1) = zero   times(2)=fourier
!!! from mytranspose times(3) = transpose   times(4) = mpi  times(5) = copy

recursive subroutine call_twoe_matel(inspfs10,inspfs20,twoematel,twoereduced,timingdir,notiming) 
  use myparams
  use myprojectmod
  implicit none
  integer ::  spf1a, spf1b, spf2a, spf2b, ii,jj,&
       itime,jtime,getlen,notiming,qqstart(nprocs),qqend(nprocs),qqblocks(nprocs),&
       kk21,kk22,kk23,  ii21,ii22,ii23, ibox,jproc
  integer, save :: xcount=0, times(10)=0,fttimes(10)=0,qqcount=0
  character :: timingdir*(*)
  DATATYPE :: inspfs10(totpoints,numspf),inspfs20(totpoints,numspf), &
       twoematel(numspf,numspf,numspf,numspf),twoereduced(totpoints,numspf,numspf)
!!$  DATATYPE, allocatable :: reducedwork2d(:,:,:,:),reducedwork1d(:,:,:),reducedwork3d(:,:,:,:,:),&
!!$       reducedhuge(:,:,:,:,:,:,:),reducedtemp(:,:,:)
  DATATYPE :: reducedwork3d(gridsize(1),gridsize(2),gridsize(3),numspf,numspf), &
       reducedhuge(numpoints(1),2,numpoints(2),2,numpoints(3),1,2),&
       reducedtemp(numpoints(1),numpoints(2),numpoints(3))
  integer :: pointsperproc(nprocs),procstart(nprocs),procend(nprocs),firsttime,lasttime
  integer, save :: maxpointsperproc
  DATATYPE :: & !!$ twoeden01(numpoints(1)),twoeden02(numpoints(1),numpoints(2)),&
       twoeden03(numpoints(1),numpoints(2),numpoints(3)),&
       tempden03(numpoints(1),numpoints(2),numpoints(3))
  DATATYPE,allocatable :: twoeden03huge(:,:,:,:,:,:,:)
  DATATYPE :: twoeden03big(numpoints(1),numpoints(2),numpoints(3),nbox(3))
!! TRYING ZDOTU !! BUGGY IN THEPAST
  complex*16 :: ZDOTU
  real*8 :: DDOT
!!  DATATYPE ::  myden(totpoints), myreduced(totpoints)  !! I GET SEGFAULTS THIS WAY
  DATATYPE,allocatable ::  myden(:), myreduced(:)


!! ZEROING TIMES... not cumulative
  times(:)=0; fttimes(:)=0

  call myclock(firsttime); itime=firsttime

  if ((myrank.eq.1).and.(notiming.eq.0)) then
     if (debugflag.eq.10) then
        open(8853, file=timingdir(1:getlen(timingdir)-1)//"/twoematel.abs.time.dat", status="unknown", position="append")
        write(8853,*) "****"
        close(8853)
        call system("date --rfc-3339=ns >>"//timingdir(1:getlen(timingdir)-1)//"/twoematel.abs.time.dat")
     endif
  endif

  if (griddim.ne.3) then
     OFLWR "OOGA DIM"; CFLST
  endif

!! for toepflag 0. if doing toepflag=0, need to reduce.
  maxpointsperproc=0
  if (orbparflag) then
     procstart(:)=1
     procend(:)=gridpoints(3)
  else
     do ii=1,nprocs
        procstart(ii)=(ii-1)*gridpoints(3)/nprocs+1
        procend(ii)=ii*gridpoints(3)/nprocs
     enddo
     do ii=1,nprocs
        pointsperproc(ii)=procend(ii)-procstart(ii)+1
        if (maxpointsperproc.lt.pointsperproc(ii)) then
           maxpointsperproc=pointsperproc(ii)
        endif
     enddo
  endif
  
  if (notwoflag.ne.0) then
     return
  endif

!$OMP PARALLEL
!$OMP MASTER
  twoereduced(:,:,:)=0d0
  reducedwork3d(:,:,:,:,:)=0d0;
!$OMP END MASTER
!$OMP END PARALLEL

  call myclock(jtime); times(1)=times(1)+jtime-itime;  

  do spf2b=1,numspf
  do spf2a=1,numspf

 ! integrating over electron 2
     call myclock(itime)
        twoeden03(:,:,:)=RESHAPE(CONJUGATE(inspfs10(:,spf2a)) * inspfs20(:,spf2b),&
             (/numpoints(1),numpoints(2),numpoints(3)/))
     call myclock(jtime); times(2)=times(2)+jtime-itime;

     if (toepflag.ne.0) then
        do ii=2,griddim
           if (gridpoints(ii).ne.gridpoints(1)) then
              OFLWR "DOME NONCUBE",gridpoints(:); CFLST
           endif
        enddo
#ifdef MPIFLAG
        if (orbparflag) then
           if (numpoints(1).ne.numpoints(2).or.numpoints(1).ne.numpoints(3)*nbox(3).or.nbox(1).gt.1.or.nbox(2).gt.1) then
              OFLWR "WOOTTTFFFF"; CFLST
           endif
           if (.not.localflag) then
              
              call myclock(itime)
              
!!$                 allocate(twoeden03big(numpoints(1),numpoints(2),numpoints(3),nbox(3)))              
              qqblocks(:)=totpoints
              do ii=1,nprocs
                 qqend(ii)=ii*totpoints; qqstart(ii)=(ii-1)*totpoints+1; 
              enddo
#ifdef REALGO
              call mygatherv_real(twoeden03,twoeden03big,totpoints*nprocs,qqstart(myrank),qqend(myrank),qqblocks(:),qqstart(:),.true.)
#else
              call mygatherv_complex(twoeden03,twoeden03big,totpoints*nprocs,qqstart(myrank),qqend(myrank),qqblocks(:),qqstart(:),.true.)
#endif
              allocate(twoeden03huge(numpoints(1),2,numpoints(2),2,numpoints(3),nbox(3),2))  
              twoeden03huge(:,:,:,:,:,:,:)=0d0; 
              twoeden03huge(:,1,:,1,:,:,1)=twoeden03big(:,:,:,:)
!!$                 deallocate(twoeden03big)
              
              call myclock(jtime); times(3)=times(3)+jtime-itime;
              
           else
              call myclock(itime)
              allocate(twoeden03huge(numpoints(1),2,numpoints(2),2,numpoints(3),1,2))  
              twoeden03huge(:,:,:,:,:,:,:)=0d0; 
              do ibox=1,nbox(3)  !! processor sending
                 jproc=(ibox+1)/2
                 if (ibox.eq.myrank.and.jproc.eq.myrank) then
                    twoeden03huge(:,1,:,1,:,1,mod(ibox-1,2)+1)=twoeden03(:,:,:)
                 else if (ibox.eq.myrank) then
                    call mympisend(twoeden03,jproc,999,totpoints)
                 else if (jproc.eq.myrank) then
                    call mympirecv(tempden03,ibox,999,totpoints)
                    twoeden03huge(:,1,:,1,:,1,mod(ibox-1,2)+1)=tempden03(:,:,:)
                 endif
              enddo
              call myclock(jtime); times(3)=times(3)+jtime-itime;
           endif
        else
#endif

           call myclock(itime)
           allocate(twoeden03huge(numpoints(1),2,numpoints(2),2,numpoints(3),1,2))  
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP MASTER
           twoeden03huge(numpoints(1),2,numpoints(2),2,numpoints(3),1,2)=0;twoeden03huge(1,1,1,1,1,1,1)=0
!$OMP END MASTER
!$OMP BARRIER
!$OMP END PARALLEL

           twoeden03huge(:,:,:,:,:,:,:)=0d0; 
           twoeden03huge(:,1,:,1,:,1,1)=twoeden03(:,:,:)
           call myclock(jtime); times(1)=times(1)+jtime-itime;
           
#ifdef MPIFLAG
        endif  !! orbparflag
        if (localflag) then
           call myclock(itime)
!!$              allocate(reducedhuge(numpoints(1),2,numpoints(2),2,numpoints(3),1,2),&
!!$                   reducedtemp(numpoints(1),numpoints(2),numpoints(3)))

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP MASTER
           reducedhuge(numpoints(1),2,numpoints(2),2,numpoints(3),1,2)=0;reducedhuge(1,1,1,1,1,1,1)=0
           reducedtemp(numpoints(1),numpoints(2),numpoints(3))=0;reducedtemp(1,1,1)=0
!$OMP END MASTER
!$OMP BARRIER
!$OMP END PARALLEL

           call myclock(jtime); times(1)=times(1)+jtime-itime; itime=jtime
#ifdef REALGO
           call circ3d_sub_real_mpi(threed_two,twoeden03huge,reducedhuge,gridpoints(3),numpoints(3),fttimes)
#else
           call circ3d_sub_mpi(threed_two,twoeden03huge,reducedhuge,gridpoints(3),numpoints(3),fttimes)
#endif
           call myclock(jtime); times(4)=times(4)+jtime-itime; itime=jtime
           
           do ibox=1,nbox(3)  !! processor receiving
              jproc=(ibox+nbox(3)+1)/2
              if (ibox.eq.myrank.and.jproc.eq.myrank) then
                 reducedwork3d(:,:,:,spf2a,spf2b)=RESHAPE(reducedhuge(:,2,:,2,:,1,mod(ibox+nbox(3)-1,2)+1),&
                      (/gridsize(1),gridsize(2),gridsize(3)/))
              else if (ibox.eq.myrank) then
                 call mympirecv(reducedwork3d(:,:,:,spf2a,spf2b),jproc,999,totpoints)
              else if (jproc.eq.myrank) then
                 reducedtemp(:,:,:)=reducedhuge(:,2,:,2,:,1,mod(ibox+nbox(3)-1,2)+1)
                 call mympisend(reducedtemp,ibox,999,totpoints)
              endif
           enddo
!!$              deallocate(reducedtemp)
           call myclock(jtime); times(5)=times(5)+jtime-itime
        else
#endif
           call myclock(itime)
!!$              allocate(reducedhuge(numpoints(1),2,numpoints(2),2,numpoints(3),nbox(3),2))
           call myclock(jtime); times(1)=times(1)+jtime-itime; itime=jtime
#ifdef REALGO
           call circ3d_sub_real(threed_two,twoeden03huge,reducedhuge,gridpoints(3))
#else
           call circ3d_sub(threed_two,twoeden03huge,reducedhuge,gridpoints(3))
#endif
           reducedwork3d(:,:,:,spf2a,spf2b)=RESHAPE(reducedhuge(:,2,:,2,:,:,2),&
                (/gridsize(1),gridsize(2),gridsize(3)/))
           call myclock(jtime); times(4)=times(4)+jtime-itime
#ifdef MPIFLAG
        endif
#endif
!!$           deallocate(reducedhuge,twoeden03huge)
        deallocate(twoeden03huge)
     else
        call myclock(itime)
        do ii23=procstart(myrank),procend(myrank)
        do kk23=1,numpoints(3)                    
        do kk22=1,numpoints(2)                 
        do kk21=1,numpoints(1)                 
           
        do ii22=1,gridsize(2)
        do ii21=1,gridsize(1)
           reducedwork3d(ii21,ii22,ii23,spf2a,spf2b)=reducedwork3d(ii21,ii22,ii23,spf2a,spf2b)+ &
                twoeden03(kk21,kk22,kk23) * threed_two(ii21-kk21,ii22-kk22,ii23-kk23-gridoffset)
        enddo
        enddo
     
        enddo
        enddo
        enddo
        enddo

        call myclock(jtime); times(4)=times(4)+jtime-itime
     endif !!TOEPLITZ

  enddo
  enddo

!! If orbparflag=.false., toepflag.ne.0, then distribute effort and allgather.
!! If orbparflag=.true., toepflag .ne.0, and MPIFLAG is set, then only have local block; no communication afterwards is needed
!! If orbparflag=.true., toepflag .eq.0, then allocate big grid and reduce (need more memory with toepflag 0)

  if (griddim.eq.3) then
     call myclock(itime)
     if (orbparflag.and.toepflag.eq.0) then
        call mympireduce(reducedwork3d,gridsize(1)*gridsize(2)*gridsize(3)*numspf**2)
     else if (toepflag.eq.0) then    !! no orbparflag
        jj=gridsize(1)*gridsize(2)
        do spf2b=1,numspf
           do spf2a=1,numspf
              call mpiallgather(reducedwork3d(:,:,:,spf2a,spf2b),jj*gridsize(3),jj*pointsperproc(:),jj*maxpointsperproc)
           enddo
        enddo
     endif
     call myclock(jtime); times(5)=times(5)+jtime-itime; 
  endif

  call myclock(itime)
  do spf2b=1,numspf
  do spf2a=1,numspf
        twoereduced(:,spf2a,spf2b) =RESHAPE(reducedwork3d(:,:,gridlow:gridhigh,spf2a,spf2b),(/totpoints/))
  enddo
  enddo

!  call mpibarrier()
!  OFLWR "***  005"; CFL
!  call mpibarrier()

  call myclock(jtime); times(1)=times(1)+jtime-itime; itime=jtime


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! reduction is performed in main MCTDHF routines NOT HERE !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(spf2a,spf2b,myreduced,myden)

!$OMP MASTER
  twoematel(:,:,:,:)=0
!$OMP END MASTER

!  call mpibarrier()
!  OFLWR "***  006"; CFL
!  call mpibarrier()

!! (YONG FYI)
!! I GET SEGFAULTS IF MYDEN AND MYREDUCED ARE AUTOMATIC VARIABLES NOT ALLOCATABLE.

!$OxxxxMP PARALLEL DEFAULT(SHARED) PRIVATE(spf2a,spf2b,myreduced,myden)

  allocate(myden(totpoints), myreduced(totpoints))

  do spf1b=1,numspf
  do spf1a=1,numspf
     myden(:)=CONJUGATE(inspfs10(:,spf1a)) * inspfs20(:,spf1b)
!$OMP DO
  do spf2b=1,numspf
  do spf2a=1,numspf
     myreduced(:)=twoereduced(:,spf2a,spf2b)

!! TRYING ZDOTU !! BUGGY IN THE PAST !!

#ifdef REALGO
     twoematel(spf2a,spf2b,spf1a,spf1b) = DDOT(totpoints,myden,1,myreduced,1)
#else
     twoematel(spf2a,spf2b,spf1a,spf1b) = ZDOTU(totpoints,myden,1,myreduced,1)
#endif
     enddo
     enddo
!$OMP END DO
  enddo
  enddo
  deallocate(myden,myreduced)
!$OMP END PARALLEL

!  call mpibarrier()
!  OFLWR "***  007"; CFL
!  call mpibarrier()


  call myclock(jtime); times(6)=times(6)+jtime-itime;  
  lasttime=jtime
  times(7)=times(7)+lasttime-firsttime

!!$  if (getpot.ne.0) then
!!$     jj=numspf**2;     ii=totpoints
!!$     call MYGEMM('N','N', ii,jj,jj,DATAONE, twoereduced(:,:,:),ii,twoden(:,:,:,:),jj, DATAONE,reducedpot(:,:,:),ii)
!!$  endif

  xcount=xcount+1

  if ((myrank.eq.1).and.(notiming.eq.0)) then
     if (xcount==1) then
        open(853, file=timingdir(1:getlen(timingdir)-1)//"/twoematel.time.dat", status="unknown")
        write(853,'(100A11)')   "etc", "1den", "1mpi", "doit", "2mpi", "dot","all","ft_zero","ft_ft","ft_tr","ft_mpi","ft_copy","ft_circ"
        close(853)
     endif
     open(853, file=timingdir(1:getlen(timingdir)-1)//"/twoematel.time.dat", status="unknown", position="append")
     write(853,'(100I11)')  times(1:7),fttimes(1:6);        close(853)
  endif

!!!fttimes
!!! times(6) = circ math
!!! from myzfft3d_par:
!!! times(1) = zero   times(2)=fourier
!!! from mytranspose times(3) = transpose   times(4) = mpi  times(5) = copy


  if (myrank.eq.1.and.(notiming.eq.0).and.debugflag.eq.10) then
        call system("date --rfc-3339=ns >>"//timingdir(1:getlen(timingdir)-1)//"/twoematel.abs.time.dat")
        open(8853, file=timingdir(1:getlen(timingdir)-1)//"/twoematel.abs.time.dat", status="unknown", position="append")
        write(8853,*) "****"
        close(8853)
  endif

  qqcount=qqcount+1
  if ((notiming.eq.0).and.debugflag.eq.10.and.qqcount.gt.1) then
     OFLWR "DEBUG10STOP"; CFLST
  endif

!contains
!  subroutine myclock(mytime)
!    integer :: values(10),mytime
!    integer, parameter :: fac(5:8)=(/60*60*1000,60*1000,1000,1/)  !! hour,minute,second,millisecond
!    call date_and_time(values=values)
!    mytime=values(8)+values(7)*fac(7)+values(6)*fac(6)+values(5)*fac(5)
!  end subroutine myclock
end subroutine call_twoe_matel



!! NOW ONLY OUTPUTS ONE. CALL IN LOOP. FOR OPENMPI TRY.

subroutine mult_reducedpot(inspfs,outspf,whichspf,reducedpot)
  use myparams
  implicit none

  integer :: ispf,kspf,whichspf
  DATATYPE :: reducedpot(totpoints, numspf,numspf),  outspf(totpoints)
  DATATYPE, intent(in) :: inspfs(totpoints, numspf)

  outspf(:)=0d0

  do ispf=whichspf,whichspf
     do kspf=1,numspf

!! CLASSES?
!        if (she lls(kspf).eq.ish ell) then
           
!! reducedpot is usual notation: <ispf | kspf> so so sum over slow index kspf
                    
        outspf(:) = outspf(:) + &
             reducedpot(:,ispf,kspf) * inspfs(:,kspf)

!        endif
     enddo
  enddo
end subroutine mult_reducedpot


function qbox(idim)
  use myparams
  implicit none
  integer :: qbox,idim
  if (orbparflag.and.idim.eq.3) then
     qbox=myrank
  else
     qbox=1
  endif
end function qbox


subroutine get_pot(outpot)
  use myparams
  implicit none
  integer :: jj,ii,kk,qbox
  DATATYPE :: outpot(totpoints)
  outpot(:)=0d0
  ii=1;  kk=totpoints
  do jj=1,griddim 
     kk=kk/numpoints(jj)
     call get_pot_onedim(outpot,jj,qbox(jj),ii,kk)
     ii=ii*numpoints(jj)
  enddo
end subroutine get_pot


subroutine get_pot_onedim(out,idim,whichbox,nnn,mmm)
  use myparams
  use myprojectmod  
  implicit none
  integer :: mmm,idim,nnn,jj,whichbox,ii
  DATATYPE :: out(nnn,numpoints(idim),mmm)
  do jj=1,mmm
     do ii=1,numpoints(idim)
        out(:,ii,jj)=out(:,ii,jj) + littlepot(idim)%mat(ii,whichbox)
     enddo
  enddo
end subroutine get_pot_onedim


subroutine get_rad(outpot)
  use myparams
  implicit none
  integer :: jj,ii,kk,qbox
  DATATYPE :: outpot(totpoints)
  outpot(:)=0d0
  ii=1;  kk=totpoints
  do jj=1,griddim 
     kk=kk/numpoints(jj)
     call get_rsq_onedim(outpot,jj,qbox(jj),ii,kk)
     ii=ii*numpoints(jj)
  enddo
  outpot(:)=sqrt(outpot)
end subroutine get_rad


subroutine get_rsq_onedim(out,idim,whichbox,nnn,mmm)
  use myparams
  use myprojectmod  
  implicit none
  integer :: mmm,idim,nnn,jj,whichbox,ii
  DATATYPE :: out(nnn,numpoints(idim),mmm)
  do jj=1,mmm
     do ii=1,numpoints(idim)
        out(:,ii,jj)=out(:,ii,jj) + sinepoints(idim)%mat(ii,whichbox)**2
     enddo
  enddo
end subroutine get_rsq_onedim


subroutine get_dipoles()
  use myparams
  use myprojectmod  
  implicit none
  integer :: jj,ii,kk,qbox
  ii=1;  kk=totpoints
  do jj=1,griddim 
     kk=kk/numpoints(jj)
     call get_one_dipole(dipoles(:,jj),jj,qbox(jj),ii,kk)
     ii=ii*numpoints(jj)
  enddo
end subroutine get_dipoles


subroutine get_one_dipole(out,idim,whichbox,nnn,mmm)
  use myparams
  use myprojectmod
  implicit none
  integer :: mmm,idim,nnn,jj,whichbox,ii
  DATATYPE :: out(nnn,numpoints(idim),mmm)
  do jj=1,mmm
     do ii=1,numpoints(idim)
        out(:,ii,jj)=sinepoints(idim)%mat(ii,whichbox)
     enddo
  enddo
end subroutine get_one_dipole


subroutine mult_ke(in, out,howmany,timingdir,notiming)
  use myparams
  implicit none
  integer :: howmany,notiming,ii
  DATATYPE :: in(totpoints,howmany), out(totpoints,howmany)
  character :: timingdir*(*)
  if (toepflag.gt.1) then
     do ii=1,howmany
        call mult_ke_toep(in(:,ii),out(:,ii))    !! slower for grids tried (up to 249)
     enddo
  else
     call mult_ke_old(in,out,howmany,timingdir,notiming)
  endif
end subroutine mult_ke


subroutine mult_ke_toep(in, out)
  use myparams
  use myprojectmod
  implicit none
  DATATYPE :: in(totpoints), out(totpoints),mywork(totpoints)
  integer :: qstart(nprocs),qblocks(nprocs),qend(nprocs),ibox,jproc,nulltimes(10)
  DATATYPE, allocatable :: bigin(:,:,:,:,:,:),bigout(:,:,:,:,:,:)
  DATATYPE, allocatable, save :: bigke(:,:,:), hugeke(:,:,:)
  integer, save :: allocated=0

#ifndef MPIFLAG
  if (orbparflag) then
     OFLWR "Ack, with orbparflag need MPIFLAG set during compilation."; CFLST
  endif
#endif
  if (mod((nprocs-1)*numpoints(3),2).eq.1) then
     OFLWR "WOOOOOD"; CFLST
  endif

  allocate(bigin(gridsize(1),2,gridsize(2),2,gridsize(3),2),&
       bigout(gridsize(1),2,gridsize(2),2,gridsize(3),2))
  bigin(:,:,:,:,:,:)=0d0

#ifdef MPIFLAG
  if (orbparflag) then
     do ibox=1,nbox(3)  !! processor sending
        jproc=(ibox+1)/2
        if (ibox.eq.myrank.and.jproc.eq.myrank) then
           bigin(:,1,:,1,:,mod(ibox-1,2)+1)=RESHAPE(in,(/numpoints(1),numpoints(2),numpoints(3)/))
        else if (ibox.eq.myrank) then
           call mympisend(in,jproc,999,totpoints)
        else if (jproc.eq.myrank) then
           call mympirecv(mywork,ibox,999,totpoints)
           bigin(:,1,:,1,:,mod(ibox-1,2)+1)=RESHAPE(mywork(:),(/numpoints(1),numpoints(2),numpoints(3)/))
        endif
     enddo
  else
#endif
  bigin(:,1,:,1,:,1)=RESHAPE(in,(/numpoints(1),numpoints(2),numpoints(3)/))
#ifdef MPIFLAG
  endif
#endif


  if (allocated.eq.0) then
     allocated=1
     allocate(bigke(0-gridsize(1):gridsize(1)-1,0-gridsize(1):gridsize(1)-1,0-gridsize(1):gridsize(1)-1))
#ifdef MPIFLAG
     if (myrank.eq.1.or.(.not.orbparflag)) then
#endif
        allocate(hugeke(0-gridpoints(1):gridpoints(1)-1,0-gridpoints(1):gridpoints(1)-1,0-gridpoints(1):gridpoints(1)-1))
        hugeke(:,:,:)=0d0
        hugeke(1-gridsize(1):gridsize(1)-1,0,0)=hugeke(1-gridsize(1):gridsize(1)-1,0,0)+&
             kevect(1)%cmat(1-gridsize(1):gridsize(1)-1)
        hugeke(0,1-gridsize(2):gridsize(2)-1,0)=hugeke(0,1-gridsize(2):gridsize(2)-1,0)+&
             kevect(2)%cmat(1-gridsize(2):gridsize(2)-1)
        hugeke(0,0,1-gridsize(3):gridsize(3)-1)=hugeke(0,0,1-gridsize(3):gridsize(3)-1)+&
             kevect(3)%cmat(1-gridsize(3):gridsize(3)-1)
#ifdef MPIFLAG
     else
        allocate(hugeke(1,1,1))
     endif
     if (orbparflag) then
        qblocks(:)=8*totpoints
        do jproc=1, nprocs
           qstart(jproc)=8*totpoints*(jproc-1)+1; qend(jproc)=8*totpoints*jproc
        enddo
#ifdef REALGO
        call myscatterv_real(hugeke,bigke,8*totpoints*nprocs,qstart(myrank),qend(myrank),qblocks(:),qstart(:))
#else
        call myscatterv_complex(hugeke,bigke,8*totpoints*nprocs,qstart(myrank),qend(myrank),qblocks(:),qstart(:))
#endif
     else
#endif
        bigke(:,:,:)=hugeke(:,:,:)
#ifdef MPIFLAG
     endif
#endif
     deallocate(hugeke)
  endif  !! allocated

#ifdef MPIFLAG 
  if (orbparflag) then
#ifdef REALGO
     call circ3d_sub_real_mpi(bigke,bigin,bigout,gridpoints(3),numpoints(3),nulltimes)
#else
     call circ3d_sub_mpi(bigke,bigin,bigout,gridpoints(3),numpoints(3),nulltimes)
#endif
     do ibox=1,nbox(3)  !! processor receiving
        jproc=(ibox+nbox(3)+1)/2
        if (ibox.eq.myrank.and.jproc.eq.myrank) then
           out(:)=RESHAPE(bigout(:,2,:,2,:,mod(ibox+nbox(3)-1,2)+1),(/totpoints/))
        else if (ibox.eq.myrank) then
           call mympirecv(out,jproc,999,totpoints)
        else if (jproc.eq.myrank) then
           mywork(:)=RESHAPE(bigout(:,2,:,2,:,mod(ibox+nbox(3)-1,2)+1),(/totpoints/))
           call mympisend(mywork,ibox,999,totpoints)
        endif
     enddo
  else
#endif
#ifdef REALGO
     call circ3d_sub_real(bigke,bigin,bigout,gridpoints(3))
#else
     call circ3d_sub(bigke,bigin,bigout,gridpoints(3))
#endif
     out(:)=RESHAPE(bigout(:,2,:,2,:,2),(/totpoints/))
#ifdef MPIFLAG
  endif
#endif
end subroutine mult_ke_toep


subroutine mult_ke_old(in, out,howmany,timingdir,notiming)
  use myparams
  implicit none
  integer :: howmany,notiming
  character :: timingdir*(*)
  DATATYPE :: in(totpoints,howmany), out(totpoints,howmany)
  call mult_allpar(in,out,1,howmany,timingdir,notiming)
end subroutine mult_ke_old

subroutine mult_xderiv(in, out,howmany)
  use myparams
  implicit none
  integer :: howmany
  character (len=200) :: timingdir="booga"
  DATATYPE :: in(totpoints,howmany), out(totpoints,howmany)
  call mult_allpar(in,out,2,howmany,timingdir,2)
end subroutine mult_xderiv

subroutine mult_yderiv(in, out,howmany)
  use myparams
  implicit none
  integer :: howmany
  character (len=200) :: timingdir="booga"
  DATATYPE :: in(totpoints,howmany), out(totpoints,howmany)
  call mult_allpar(in,out,3,howmany,timingdir,2)
end subroutine mult_yderiv

subroutine mult_zderiv(in, out,howmany)
  use myparams
  implicit none
  integer :: howmany
  character (len=200) :: timingdir="booga"
  DATATYPE :: in(totpoints,howmany), out(totpoints,howmany)
  call mult_allpar(in,out,4,howmany,timingdir,2)
end subroutine mult_zderiv


recursive subroutine mult_allpar(in, out,inoption,howmany,timingdir,notiming)
  use myparams
  implicit none
  integer :: idim,inoption,option,howmany,notiming
  DATATYPE :: in(totpoints,howmany), out(totpoints,howmany), temp(totpoints,howmany)
  logical :: dodim(3)
  character :: timingdir*(*)

  if (griddim.ne.3) then
     OFLWR "ERWRESTOPPP"; CFLST
  endif

  dodim(:)=.false.
  if (inoption.eq.1) then
     option=1
     dodim(:)=.true.
  else if (inoption.lt.1.or.inoption.gt.4) then
     OFLWR "OWWOWORE WHAT?"; CFLST
     option=99
  else
     option=2
     dodim(inoption-1)=.true.
  endif

  out(:,:)=0d0

  if (.not.orbparflag) then
     if (nbox(1).ne.1.or.nbox(2).ne.1.or.nbox(3).ne.1) then
        OFLWR "OOOFSSxxxF",nbox; CFLST
     endif
     do idim=1,griddim
        if (dodim(idim)) then
           call mult_allone(in,temp,idim,option,1,1,howmany)
           out(:,:)=out(:,:)+temp(:,:)
        endif
     enddo
  else
     if (nbox(1).ne.1.or.nbox(2).ne.1.or.nbox(3).ne.nprocs) then
        OFLWR "OOOFSSF"; CFLST
     endif
     do idim=1,2
        if (dodim(idim)) then
           call mult_allone(in,temp,idim,option,1,1,howmany)
           out(:,:)=out(:,:)+temp(:,:)
        endif
     enddo
     idim=3
     if (dodim(idim)) then
        call mult_alltoall_z(in,temp,option,howmany,timingdir,notiming)
        out(:,:)=out(:,:)+temp(:,:)

!!$        do ibox=1,nprocs
!!$           do jbox=1,nprocs
!!$              call mult_allone(in,temp,idim,option,ibox,jbox,howmany)
!!$              if (jbox.eq.myrank) then
!!$                 out(:,:)=out(:,:)+temp(:,:)
!!$              endif
!!$           enddo
!!$        enddo
        
     endif
  endif
end subroutine mult_allpar



recursive subroutine mult_alltoall_z(in, out,option,howmany,timingdir,notiming)
  use myparams
  use myprojectmod  
  implicit none
  integer :: nnn,option,ii,howmany
  DATATYPE :: in(numpoints(1)*numpoints(2),numpoints(3),howmany),&
       out(numpoints(1)*numpoints(2),numpoints(3),howmany),&
       work(numpoints(1)*numpoints(2),numpoints(3),nprocs,howmany),&
       work2(numpoints(1)*numpoints(2),numpoints(3),howmany,nprocs),&
       work3(numpoints(1)*numpoints(2),numpoints(3),howmany,nprocs)
  integer :: times(10),atime,btime,notiming,getlen
  character :: timingdir*(*)
  integer, save :: xcount=0

  times(:)=0

  if (nprocs.ne.nbox(3)) then
     OFLWR "EEGNOT STOP",nprocs,nbox(3); CFLST
  endif
  if (totpoints.ne.numpoints(1)*numpoints(2)*numpoints(3)) then
     OFLWR "WHAAAAAZZZZ?",totpoints,numpoints(1),numpoints(2),numpoints(3); CFLST
  endif

  call myclock(atime)

  nnn=numpoints(1)*numpoints(2)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,atime,btime)
  select case(option)
  case(1)  !! KE
!$OMP DO SCHEDULE(STATIC)
     do ii=1,howmany
        call MYGEMM('N','T',nnn,gridpoints(3),numpoints(3),DATAONE,in(:,:,ii),nnn,ketot(3)%mat(1,1,1,myrank),gridpoints(3),DATAZERO, work(:,:,:,ii), nnn)
     enddo
!$OMP END DO
  case(2) 
!$OMP DO SCHEDULE(STATIC)
     do ii=1,howmany
        call MYGEMM('N','T',nnn,gridpoints(3),numpoints(3),DATAONE,in(:,:,ii),nnn,fdtot(3)%mat(1,1,1,myrank),gridpoints(3),DATAZERO, work(:,:,:,ii), nnn)
     enddo
!$OMP END DO
  case default 
     OFLWR "WHAAAAT"; CFLST
  end select

!$OMP BARRIER

call myclock(btime); times(1)=times(1)+btime-atime; atime=btime

!$OMP DO SCHEDULE(STATIC)
  do ii=1,nprocs
     work2(:,:,:,ii)=work(:,:,ii,:)
  enddo
!$OMP END DO

!$OMP BARRIER

call myclock(btime); times(2)=times(2)+btime-atime; atime=btime

!$OMP MASTER
  call mympialltoall(work2,work3,totpoints*howmany)
!$OMP END MASTER

!$OMP BARRIER
!$OMP END PARALLEL

call myclock(btime); times(3)=times(3)+btime-atime; atime=btime

  out(:,:,:)=0d0
  do ii=1,nprocs
     out(:,:,:)=out(:,:,:)+work3(:,:,:,ii)
  enddo

  call myclock(btime); times(4)=times(4)+btime-atime

  if (debugflag.eq.42.and.myrank.eq.1) then
     xcount=xcount+1
     if (xcount==1) then
        open(2853, file=timingdir(1:getlen(timingdir)-1)//"/zke.time.dat", status="unknown")
        write(2853,'(100A11)')   "mult", "arrange","mpi","reduce"
        close(2853) 
     endif
     open(2853, file=timingdir(1:getlen(timingdir)-1)//"/zke.time.dat", status="unknown", position="append")
     write(2853,'(100I11)')  times(1:4);        close(2853)
  endif
end subroutine mult_alltoall_z






subroutine mult_allone(in, out,idim,option,ibox,jbox,howmany)
  use myparams
  implicit none
  integer :: mmm,idim,nnn,jdim,option,ibox,jbox,howmany
  DATATYPE :: in(totpoints,howmany), out(totpoints,howmany)
  if (ibox.lt.1.or.jbox.lt.1) then
     OFLWR "ERGSTTTTOP"; CFLST
  endif
  if (ibox.gt.nbox(idim).or.jbox.gt.nbox(idim)) then
     OFLWR "OOGSTOP"; CFLST
  endif
  if ((.not.orbparflag).or.idim.ne.3) then
     if (ibox.ne.1.or.jbox.ne.1) then
        OFLWR "BOOBOBOBXSTOP"; CFLST
     endif
  endif
  if ((.not.orbparflag).or.(idim.ne.3).or.(ibox.eq.myrank)) then
     nnn=1; mmm=1
     do jdim=1,idim-1
        nnn=nnn*numpoints(jdim)
     enddo
     do jdim=idim+1,griddim
        mmm=mmm*numpoints(jdim)
     enddo
     call mult_all0(in,out,ibox,jbox,idim,nnn,mmm*howmany,option)
#ifdef MPIFLAG
     if (orbparflag.and.(idim.eq.3).and.(jbox.ne.ibox)) then
        call mympisend(out,jbox,999,totpoints*howmany)
     endif
#endif
  endif
#ifdef MPIFLAG
  if (orbparflag.and.(idim.eq.3).and.(jbox.ne.ibox).and.(jbox.eq.myrank)) then
     call mympirecv(out,ibox,999,totpoints*howmany)
  endif
#endif

end subroutine mult_allone



subroutine mult_all0(in, out,ibox,jbox,idim,nnn,mmm,option)
  use myparams
  use myprojectmod  
  implicit none
  integer :: mmm,idim,nnn,jj,ibox,jbox,option
  DATATYPE :: in(nnn,numpoints(idim),mmm)
  DATATYPE :: out(nnn,numpoints(idim),mmm)
  select case(option)
  case(1)  !! KE
     do jj=1,mmm
        call MYGEMM('N','T',nnn,numpoints(idim),numpoints(idim),DATAONE,in(:,:,jj),nnn,ketot(idim)%mat(1,jbox,1,ibox),gridpoints(idim),DATAZERO, out(:,:,jj), nnn)
     enddo
  case(2)  !! X Y or Z derivative (real valued antisymmetric)
     do jj=1,mmm
        call MYGEMM('N','T',nnn,numpoints(idim),numpoints(idim),DATAONE,in(:,:,jj),nnn,fdtot(idim)%mat(1,jbox,1,ibox),gridpoints(idim),DATAZERO, out(:,:,jj), nnn)
     enddo
  case default 
     OFLWR "WHAAAAT"; CFLST
  end select
end subroutine mult_all0





!! (do need this -- set to zero)
subroutine hatom_matel(inspfs1, inspfs2, hatommatel,numberspf) 
  use myparams
  implicit none
  integer :: numberspf
  DATATYPE :: inspfs1(totpoints,numberspf), inspfs2(totpoints,numberspf),hatommatel(numberspf,numberspf)
  hatommatel(:,:)=0d0
end subroutine hatom_matel

subroutine hatom_op(inspf, outspf)
  use myparams
  implicit none
  DATATYPE :: inspf(totpoints),outspf(totpoints)
  outspf(:)=0d0
end subroutine hatom_op


!!$subroutine mult_ke_toep_old(in, out)
!!$  use myparams
!!$  implicit none
!!$  integer :: idim
!!$  DATATYPE :: in(totpoints), out(totpoints), temp(totpoints)
!!$  if (orbparflag) then
!!$     OFLWR "MULT_ke_toep not supported for orbparflag yet"; CFLST
!!$  endif
!!$  out(:)=0d0
!!$  do idim=1,griddim
!!$     call mult_ke0toep(in, temp,idim)      !!   1BOX   1,1
!!$     out(:)=out(:)+temp(:)
!!$!not debugged noncube
!!$     call liltranspose(in,numpoints(idim),totpoints/numpoints(idim))
!!$     call liltranspose(out,numpoints(idim),totpoints/numpoints(idim))
!!$  enddo
!!$end subroutine mult_ke_toep_old
!!$
!!$subroutine liltranspose(mat,idim,jdim)
!!$  implicit none
!!$  integer :: idim,jdim
!!$  DATATYPE :: mat(idim*jdim),temp(idim*jdim)
!!$  temp(:)=RESHAPE(TRANSPOSE(RESHAPE(mat,(/idim,jdim/))),(/idim*jdim/))
!!$  mat(:)=temp(:)
!!$end subroutine liltranspose
!!$
!!$subroutine mult_ke0toep(in, out,idim)
!!$  use myparams
!!$  use myprojectmod  
!!$  implicit none
!!$  integer :: mmm,idim
!!$  DATATYPE :: in(numpoints(idim),totpoints/numpoints(idim))
!!$  DATATYPE :: out(numpoints(idim),totpoints/numpoints(idim))
!!$  if (orbparflag) then
!!$     OFLWR "KE TOEP NOT SUPPORTED ORBPARFLAG"; CFLST
!!$  endif
!!$  mmm=totpoints/numpoints(idim)
!!$#ifndef REALGO
!!$     call toeplitz1d_sub(kevect(idim)%cmat(:),in(:,:),out(:,:),numpoints(idim),mmm)
!!$#else
!!$     call toeplitz1d_sub_real(kevect(idim)%cmat(:),in(:,:),out(:,:),numpoints(idim),mmm)
!!$#endif
!!$end subroutine mult_ke0toep

!!! is not fast.  for debug of fft
!!!subroutine mult_ke_fast(in, out)
!!!  use myparams
!!!  implicit none
!!!  integer :: idim
!!!  DATATYPE :: in(totpoints), out(totpoints), temp(totpoints)
!!!
!!!  out(:)=0d0
!!!  do idim=1,griddim
!!!     call mult_ke0(in, temp,1,1,idim,1,totpoints/numpoints(idim))    !!   1BOX   1,1
!!!     out(:)=out(:)+temp(:)
!!!!     call liltranspose(in,totpoints/numpoints(idim),numpoints(idim))   think it's the other way but haven't debugged
!!!!     call liltranspose(out,totpoints/numpoints(idim),numpoints(idim))  with noncube
!!!
!!!     call liltranspose(in,numpoints(idim),totpoints/numpoints(idim))
!!!     call liltranspose(out,numpoints(idim),totpoints/numpoints(idim))
!!!
!!!
!!!  enddo
!!!
!!!end subroutine mult_ke_fast
!!!
