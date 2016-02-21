
!! SUBROUTINES FOR COMPUTING MEAN FIELD !!

!! REDUCEDPOT is usual notation bra, ket   other reduced mats opposite notation for ease of actreducedconjg

#include "Definitions.INC"


subroutine get_allden()
  implicit none
  call get_reducedham();  call getdenmatx();
end subroutine get_allden


subroutine get_reducedham()
  use parameters
  use configmod
  use xxxmod
  implicit none
  integer :: itime,jtime,times(100)

  call system_clock(itime)
  call get_tworeducedx(www,yyy%reducedpottally,yyy%cmfavec(:,:,0),yyy%cmfavec(:,:,0),mcscfnum)
  call system_clock(jtime);  times(4)=times(4)+jtime-itime

  call system_clock(itime)
  call get_reducedproderiv(www,yyy%reducedproderiv,yyy%cmfavec(:,:,0),yyy%cmfavec(:,:,0),mcscfnum)
  call system_clock(jtime);  times(6)=times(6)+jtime-itime

  call system_clock(itime)
  call get_reducedr(www,yyy%reducedinvr,yyy%reducedinvrsq,yyy%reducedr,&
       yyy%cmfavec(:,:,0),yyy%cmfavec(:,:,0),mcscfnum)
  call system_clock(itime);  times(7)=times(7)+itime-jtime

end subroutine get_reducedham


!! IN THE END reducedpottally is bra2,ket2,bra1,ket1 ; 
!! I sum it up by conjugating ket2 bra2 ket1 bra1.  Notice a1 is conjg not a2.
!! reducedpottally is the TWO ELECTRON REDUCED DENSITY MATRIX.

subroutine get_reducedpot()
  use xxxmod
  use configmod
  use opmod
  implicit none
  call get_reducedpot0(www,yyy%reducedpottally(:,:,:,:,0), yyy%reducedpot(:,:,:,0),twoereduced(:,:,:))
end subroutine get_reducedpot


subroutine get_reducedpot0(www,intwoden,outpot,twoereduced)
  use walkmod
  use spfsize_parameters
  implicit none
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: twoereduced(reducedpotsize, www%nspf,www%nspf),&
       intwoden(www%nspf,www%nspf,www%nspf,www%nspf)
  DATATYPE,intent(out) :: outpot(reducedpotsize,www%nspf,www%nspf)
  integer :: ii

! I have bra2,ket2,bra1,ket1. (chemists' notation for contraction)

! e.g. M= -2, -1, 2, 1 all different classes. 
!   no zeroed entries even if constraining, in general, I believe.

  ii=reducedpotsize

!  if(deb ugflag.ne.0) then
!     temptwoden(:,:,:,:)=intwoden(:,:,:,:)
!     do ispf=1,www%nspf
!        do jspf=1,www%nspf
!           if (orbclass(ispf).ne.orbclass(jspf)) then
!              temptwoden(:,ispf,:,jspf)=0d0
!              temptwoden(ispf,:,jspf,:)=0d0
!           endif
!        enddo
!     enddo
!     call MYGEMM('N','N', ii,www%nspf**2,www%nspf**2,(1.0d0,0.d0), &
!         twoereduced,ii,temptwoden,www%nspf**2, (0.d0,0.d0),outpot,ii)
!  else

     call MYGEMM('N','N', ii,www%nspf**2,www%nspf**2,(1.0d0,0.d0), &
          twoereduced,ii,intwoden,www%nspf**2, (0.d0,0.d0),outpot,ii)


end subroutine get_reducedpot0



subroutine get_tworeducedx(www,reducedpottally,avector1,in_avector2,numvects)
  use fileptrmod
  use r_parameters
  use walkmod
  implicit none
  integer,intent(in) :: numvects
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: avector1(numr,www%firstconfig:www%lastconfig,numvects),  &
       in_avector2(numr,www%firstconfig:www%lastconfig,numvects)
  DATATYPE,intent(out) :: reducedpottally(www%nspf,www%nspf,www%nspf,www%nspf)
  DATATYPE,allocatable :: avector2(:,:,:)
  DATATYPE,allocatable :: mytally(:,:,:,:)
  DATATYPE ::  a1(numr,numvects), a2(numr,numvects), dot, csum
  DATAECS :: rvalues(numr)
  integer ::   ispf, jspf, iispf, jjspf ,  config2, config1,dirphase, iwalk,ii,ihop

!! DO SUMMA (parconsplit.ne.0 and sparsesummaflag.eq.2, "circ")

  allocate(avector2(numr,www%numconfig,numvects));    avector2(:,:,:)=0d0
  if (www%lastconfig.ge.www%firstconfig) then
     avector2(:,www%firstconfig:www%lastconfig,:)=in_avector2(:,:,:)
  endif
  if (www%parconsplit.ne.0) then
     do ii=1,numvects
        call mpiallgather(avector2(:,:,ii),www%numconfig*numr,&
             www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
     enddo
  endif

  reducedpottally(:,:,:,:)=0.d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rvalues,config1,a1,a2,iwalk,dirphase,config2,ii,ispf,iispf,jspf,jjspf,ihop,csum,mytally)

!! REDUCTION(+:reducedpottally) caused crashes!  OMP CRITICAL instead below.

  allocate(mytally(www%nspf,www%nspf,www%nspf,www%nspf))
  mytally(:,:,:,:)=0d0

  rvalues(:)=bondpoints(:)

        !! doubly off diagonal walks

!$OMP DO SCHEDULE(DYNAMIC)
  do config1=www%botconfig,www%topconfig

     a1(:,:)=avector1(:,config1,:)
        
     do ihop=1,www%numdoublehops(config1)
        config2=www%doublehop(ihop,config1)

        do ii=1,numvects
           a2(:,ii)=avector2(:,config2,ii)/rvalues(:)
        enddo

        csum=dot(a1,a2,numvects*numr)           !! 1/R factor above

        do iwalk=www%doublehopwalkstart(ihop,config1),www%doublehopwalkend(ihop,config1)
           
           dirphase=www%doublewalkdirphase(iwalk,config1)

           ispf=www%doublewalkdirspf(1,iwalk,config1)   !BRA2 
           jspf=www%doublewalkdirspf(2,iwalk,config1)   !KET2 (walk)
           iispf=www%doublewalkdirspf(3,iwalk,config1)  !BRA1
           jjspf=www%doublewalkdirspf(4,iwalk,config1)  !KET1 (walk)
           
           mytally(ispf,jspf,iispf,jjspf) =  &    
                mytally(ispf,jspf,iispf,jjspf) +  &
                dirphase*csum           !! 1/R factor above
           
           mytally(iispf,jjspf,ispf,jspf) =  &
                mytally(iispf,jjspf,ispf,jspf) + &
                dirphase*csum           !! 1/R factor above

        enddo

     enddo
  enddo   ! config1
!$OMP END DO
!$OMP CRITICAL
  reducedpottally(:,:,:,:)=reducedpottally(:,:,:,:)+mytally(:,:,:,:)
!$OMP END CRITICAL
  deallocate(mytally)
!$OMP END PARALLEL

  deallocate(avector2)

  call mympireduce(reducedpottally(:,:,:,:), www%nspf**4)

end subroutine get_tworeducedx




subroutine get_reducedproderiv(www,reducedproderiv,avector1,in_avector2,numvects)
  use fileptrmod
  use r_parameters
  use opmod   !! rkemod, proderivmod
  use walkmod
  implicit none
  integer,intent(in) :: numvects
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: avector1(numr,www%firstconfig:www%lastconfig,numvects), &
       in_avector2(numr,www%firstconfig:www%lastconfig,numvects)
  DATATYPE,intent(out) :: reducedproderiv(www%nspf,www%nspf)
  DATATYPE,allocatable :: avector2(:,:,:)
  DATATYPE :: a1(numr,numvects), a2(numr,numvects), a2mult(numr,numvects), mypro(numr,numr), &
       myredpro(www%nspf,www%nspf),dot,csum
  integer ::  config1,config2,  ispf,jspf,  dirphase,     iwalk,ii,ihop

!! DO SUMMA (parconsplit.ne.0 and sparsesummaflag.eq.2, "circ")

  allocate(avector2(numr,www%numconfig,numvects))
  avector2(:,:,:)=0d0
  if (www%lastconfig.ge.www%firstconfig) then
     avector2(:,www%firstconfig:www%lastconfig,:)=in_avector2(:,:,:)
  endif

  if (www%parconsplit.ne.0) then
     do ii=1,numvects
        call mpiallgather(avector2(:,:,ii),www%numconfig*numr,&
             www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
     enddo
  endif

  reducedproderiv(:,:)=0.d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(a1,a2,mypro,myredpro,config1,iwalk,config2,dirphase,ispf,jspf,a2mult,ihop,csum) 

!! removing REDUCTION(+:reducedproderiv) due to get_twoereducedx problems OMP CRITICAL instead

  mypro(:,:)=proderivmod(:,:)
  myredpro(:,:)=0d0

!$OMP DO SCHEDULE(DYNAMIC)
  do config1=www%botconfig,www%topconfig

     do ihop=1,www%numsinglehops(config1)
        config2=www%singlehop(ihop,config1)

        a1(:,:)=avector1(:,config1,:)
        a2(:,:)=avector2(:,config2,:)

        call MYGEMM('N','N',numr,numvects,numr,DATAONE,&
             mypro(:,:),numr,a2(:,:),numr,DATAZERO,a2mult(:,:),numr)

        csum=dot(a1,a2mult,numvects*numr)

        do iwalk=www%singlehopwalkstart(ihop,config1),www%singlehopwalkend(ihop,config1)

           dirphase=www%singlewalkdirphase(iwalk,config1)
           ispf=www%singlewalkopspf(1,iwalk,config1)
           jspf=www%singlewalkopspf(2,iwalk,config1)

           myredpro(jspf,ispf)=myredpro(jspf,ispf)+ &
                dirphase*csum

        enddo

! PREVIOUS (jj,ii private)
!        do jj=1,numr
!           do ii=1,numr
!              a1(:)=avector1(ii,config1,:)
!              a2(:)=avector2(jj,config2,:)
!              
!              reducedproderiv(jspf,ispf)=reducedproderiv(jspf,ispf)+ &
!                   dirphase*dot(a1,a2,numvects)*proderivmod(ii,jj)
!           enddo
!        enddo

     enddo
  enddo
!$OMP END DO
!$OMP CRITICAL
  reducedproderiv(:,:)=reducedproderiv(:,:)+myredpro(:,:)
!$OMP END CRITICAL
!$OMP END PARALLEL

  deallocate(avector2)

  call mympireduce(reducedproderiv(:,:), www%nspf**2)

end subroutine get_reducedproderiv


subroutine get_reducedr(www,reducedinvr,reducedinvrsq,reducedr,avector1,in_avector2,numvects)
  use fileptrmod
  use r_parameters
  use walkmod
  implicit none
  integer,intent(in) :: numvects
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: avector1(numr,www%firstconfig:www%lastconfig,numvects), &
       in_avector2(numr,www%firstconfig:www%lastconfig,numvects)
  DATATYPE,intent(out) :: reducedinvr(www%nspf,www%nspf),reducedr(www%nspf,www%nspf), &
       reducedinvrsq(www%nspf,www%nspf)
  DATATYPE,allocatable :: avector2(:,:,:)
  DATATYPE ::  a1(numr,numvects), a2(numr,numvects), a2r(numr,numvects), &
       a2inv(numr,numvects), a2invsq(numr,numvects), dot, rdot,invdot,invsqdot
  integer ::  config1,config2,   ispf,jspf,  dirphase,    iwalk,ii, ihop
  DATAECS ::  invrvalues(numr),invrsqvalues(numr),rvalues(numr)
  DATATYPE :: myinvr(www%nspf,www%nspf),myr(www%nspf,www%nspf),  myinvrsq(www%nspf,www%nspf)

  reducedinvr(:,:)=0.d0;  reducedr(:,:)=0.d0;  reducedinvrsq(:,:)=0.d0

  if (numr.eq.1) then
     return
  endif

!! DO SUMMA (parconsplit.ne.0 and sparsesummaflag.eq.2, "circ")

  allocate(avector2(numr,www%numconfig,numvects))
  avector2(:,:,:)=0d0

  if (www%lastconfig.ge.www%firstconfig) then
     avector2(:,www%firstconfig:www%lastconfig,:) = in_avector2(:,:,:)
  endif
  if (www%parconsplit.ne.0) then
     do ii=1,numvects
        call mpiallgather(avector2(:,:,ii),www%numconfig*numr,&
             www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
     enddo
  endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rvalues,invrvalues,invrsqvalues,config1,a1,a2,a2r,a2inv,a2invsq,iwalk,config2,dirphase,ispf,jspf,myinvr,myr,myinvrsq,ihop,invdot,rdot,invsqdot) 

!! removing REDUCTION(+:reducedinvr,reducedr,reducedinvrsq) due to get_twoereducedx problems CRITICAL instead

  rvalues(:)=bondpoints(:)
  invrvalues(:)=(1.d0/bondpoints(:))
  invrsqvalues(:)=(1.d0/bondpoints(:)**2)

  myinvr(:,:)=0.d0;  myr(:,:)=0.d0;  myinvrsq(:,:)=0.d0

!$OMP DO SCHEDULE(DYNAMIC)
  do config1=www%botconfig,www%topconfig
     a1(:,:)=avector1(:,config1,:)

     do ihop=1,www%numsinglehops(config1)
        config2=www%singlehop(ihop,config1)

        a2(:,:)=avector2(:,config2,:)
        do ii=1,numvects
           a2r(:,ii)=a2(:,ii)*rvalues(:)
           a2inv(:,ii)=a2(:,ii)*invrvalues(:)
           a2invsq(:,ii)=a2(:,ii)*invrsqvalues(:)
        enddo

        invdot=dot(a1,a2inv,numvects*numr)
        rdot=dot(a1,a2r,numvects*numr)
        invsqdot=dot(a1,a2invsq,numvects*numr)

        do iwalk=www%singlehopwalkstart(ihop,config1),www%singlehopwalkend(ihop,config1)

           dirphase=www%singlewalkdirphase(iwalk,config1)
           ispf=www%singlewalkopspf(1,iwalk,config1)
           jspf=www%singlewalkopspf(2,iwalk,config1)

           myinvr(jspf,ispf)=myinvr(jspf,ispf)+         dirphase*invdot
           myr(jspf,ispf)=myr(jspf,ispf)+         dirphase*rdot
           myinvrsq(jspf,ispf)=myinvrsq(jspf,ispf)+     dirphase*invsqdot

        enddo

     enddo
  enddo
!$OMP END DO
!$OMP CRITICAL
  reducedinvr(:,:)=reducedinvr(:,:)+myinvr(:,:)
  reducedinvrsq(:,:)=reducedinvrsq(:,:)+myinvrsq(:,:)
  reducedr(:,:)=reducedr(:,:)+myr(:,:)
!$OMP END CRITICAL
!$OMP END PARALLEL

  deallocate(avector2)
  
  call mympireduce(reducedr(:,:), www%nspf**2)
  call mympireduce(reducedinvr(:,:), www%nspf**2)
  call mympireduce(reducedinvrsq(:,:), www%nspf**2)

end subroutine get_reducedr

