
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
  call get_tworeducedx(www,yyy%reducedpottally,yyy%cmfpsivec(astart(1),0),yyy%cmfpsivec(astart(1),0),mcscfnum)
  call system_clock(jtime);  times(4)=times(4)+jtime-itime

  call system_clock(itime)
  call get_reducedproderiv(www,yyy%reducedproderiv,yyy%cmfpsivec(astart(1),0),yyy%cmfpsivec(astart(1),0),mcscfnum)
  call system_clock(jtime);  times(6)=times(6)+jtime-itime

  call system_clock(itime)
  call get_reducedr(www,yyy%reducedinvr,yyy%reducedinvrsq,yyy%reducedr,yyy%cmfpsivec(astart(1),0),yyy%cmfpsivec(astart(1),0),mcscfnum)
  call system_clock(itime);  times(7)=times(7)+itime-jtime

end subroutine get_reducedham


!! IN THE END reducedpottally is bra2,ket2,bra1,ket1 ; I sum it up by conjugating ket2 bra2 ket1 bra1.  Notice a1 is conjg not a2.

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
  DATATYPE :: TWOEreduced(reducedpotsize, www%nspf,www%nspf)
  DATATYPE :: intwoden(www%nspf,www%nspf,www%nspf,www%nspf),outpot(reducedpotsize,www%nspf,www%nspf)
  integer :: ii

! I have bra2,ket2,bra1,ket1. (chemists' notation for contraction)

! e.g. M= -2, -1, 2, 1 all different classes.  no zeroed entries even if constraining, in general, I believe.

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
!     call MYGEMM('N','N', ii,www%nspf**2,www%nspf**2,(1.0d0,0.d0), twoereduced,ii,temptwoden,www%nspf**2, (0.d0,0.d0),outpot,ii)
!  else

     call MYGEMM('N','N', ii,www%nspf**2,www%nspf**2,(1.0d0,0.d0), twoereduced,ii,intwoden,www%nspf**2, (0.d0,0.d0),outpot,ii)


end subroutine get_reducedpot0



subroutine get_tworeducedx(www,reducedpottally,avector1,in_avector2,numvects)
  use fileptrmod
  use r_parameters
  use walkmod
  implicit none
  integer,intent(in) :: numvects
  type(walktype),intent(in) :: www
  integer ::   ispf, jspf, iispf, jjspf ,  config2, config1,dirphase, iwalk, qq
  DATATYPE,intent(in) :: avector1(numr,www%firstconfig:www%lastconfig,numvects),  in_avector2(numr,www%firstconfig:www%lastconfig,numvects)
  DATATYPE,intent(out) :: reducedpottally(www%nspf,www%nspf,www%nspf,www%nspf)
  DATATYPE :: avector2(numr,www%numconfig,numvects)             !! AUTOMATIC
  DATATYPE ::  a1(numvects), a2(numvects), dot
  DATAECS :: thisrvalue

  avector2(:,:,:)=0d0
  avector2(:,www%firstconfig:www%lastconfig,:)=in_avector2(:,:,:)

!! DO SUMMA
  if (www%parconsplit.ne.0) then
     do qq=1,numvects
        call mpiallgather(avector2(:,:,qq),www%numconfig*numr,www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
     enddo
  endif

  reducedpottally(:,:,:,:)=0.d0

  !! by configurations; following matel.f90

  do qq=1,numr
     thisrvalue=bondpoints(qq)

        !! doubly off diagonal walks
        
     do config1=www%botconfig,www%topconfig

        a1(:)=avector1(qq,config1,:)   !! NO MORE
        
        do iwalk=1,www%numdoublewalks(config1)
           
           dirphase=www%doublewalkdirphase(iwalk,config1)
           config2=www%doublewalk(iwalk,config1)
           a2(:)=avector2(qq,config2,:)
           
           ispf=www%doublewalkdirspf(1,iwalk,config1)   !BRA2 
           jspf=www%doublewalkdirspf(2,iwalk,config1)   !KET2 (walk)
           iispf=www%doublewalkdirspf(3,iwalk,config1)  !BRA1
           jjspf=www%doublewalkdirspf(4,iwalk,config1)  !KET1 (walk)
           
           reducedpottally(ispf,jspf,iispf,jjspf) =  &    
                reducedpottally(ispf,jspf,iispf,jjspf) +  &
                dirphase*dot(a1,a2,numvects)/thisrvalue
           
           reducedpottally(iispf,jjspf,ispf,jspf) =  &
                reducedpottally(iispf,jjspf,ispf,jspf) + &
                dirphase*dot(a1,a2,numvects)/thisrvalue

        enddo
     enddo   ! config1
  enddo  !! qq

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
  DATATYPE,intent(in) :: avector1(numr,www%firstconfig:www%lastconfig,numvects), in_avector2(numr,www%firstconfig:www%lastconfig,numvects)
  DATATYPE,intent(out) :: reducedproderiv(www%nspf,www%nspf)
  DATATYPE :: avector2(numr,www%numconfig,numvects)     !! AUTOMATIC
  DATATYPE :: a1(numvects), a2(numvects), dot
  integer ::  config1,config2,  ispf,jspf,  dirphase,     iwalk,ii,jj

  avector2(:,:,:)=0d0
  avector2(:,www%firstconfig:www%lastconfig,:)=in_avector2(:,:,:)

!! DO SUMMA
  if (www%parconsplit.ne.0) then
     do ii=1,numvects
        call mpiallgather(avector2(:,:,ii),www%numconfig*numr,www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
     enddo
  endif

  reducedproderiv(:,:)=0.d0

  do config1=www%botconfig,www%topconfig

     do iwalk=1,www%numsinglewalks(config1)
        config2=www%singlewalk(iwalk,config1);        dirphase=www%singlewalkdirphase(iwalk,config1)
        ispf=www%singlewalkopspf(1,iwalk,config1);        jspf=www%singlewalkopspf(2,iwalk,config1)

        do jj=1,numr
           do ii=1,numr
              a1(:)=avector1(ii,config1,:)
              a2(:)=avector2(jj,config2,:)
              
              reducedproderiv(jspf,ispf)=reducedproderiv(jspf,ispf)+ &
                   dirphase*dot(a1,a2,numvects)*proderivmod(ii,jj)
           enddo
        enddo
     enddo
  enddo

  call mympireduce(reducedproderiv(:,:), www%nspf**2)

end subroutine get_reducedproderiv


subroutine get_reducedr(www,reducedinvr,reducedinvrsq,reducedr,avector1,in_avector2,numvects)
  use fileptrmod
  use r_parameters
  use walkmod
  implicit none
  integer,intent(in) :: numvects
  type(walktype),intent(in) :: www
  DATATYPE,intent(in) :: avector1(numr,www%firstconfig:www%lastconfig,numvects), in_avector2(numr,www%firstconfig:www%lastconfig,numvects)
  DATATYPE,intent(out) :: reducedinvr(www%nspf,www%nspf),reducedr(www%nspf,www%nspf),  reducedinvrsq(www%nspf,www%nspf)
  DATATYPE :: avector2(numr,www%numconfig,numvects)      !! AUTOMATIC
  DATATYPE ::  a1(numvects), a2(numvects), dot
  integer ::  config1,config2,   ispf,jspf,  dirphase,    iwalk,ii
  DATAECS :: thisrvalue,  csum,csum2

  reducedinvr(:,:)=0.d0;  reducedr(:,:)=0.d0;  reducedinvrsq(:,:)=0.d0
  if (numr.eq.1) then
     return
  endif

!! DO SUMMA

  avector2(:,:,:)=0d0
  avector2(:,www%firstconfig:www%lastconfig,:) = in_avector2(:,:,:)
  if (www%parconsplit.ne.0) then
     do ii=1,numvects
        call mpiallgather(avector2(:,:,ii),www%numconfig*numr,www%configsperproc(:)*numr,www%maxconfigsperproc*numr)
     enddo
  endif

  do ii=1,numr

     thisrvalue=bondpoints(ii);        csum=(1.d0/thisrvalue);        csum2=(1.d0/thisrvalue**2)

     do config1=www%botconfig,www%topconfig

        a1(:)=avector1(ii,config1,:)
        
        do iwalk=1,www%numsinglewalks(config1)
           config2=www%singlewalk(iwalk,config1);              dirphase=www%singlewalkdirphase(iwalk,config1)
           a2(:)=avector2(ii,config2,:)
           
           ispf=www%singlewalkopspf(1,iwalk,config1);              jspf=www%singlewalkopspf(2,iwalk,config1)
           
           
           reducedinvr(jspf,ispf)=reducedinvr(jspf,ispf)+         dirphase*dot(a1,a2,numvects)*csum
           reducedr(jspf,ispf)=reducedr(jspf,ispf)+         dirphase*dot(a1,a2,numvects)*thisrvalue
           reducedinvrsq(jspf,ispf)=reducedinvrsq(jspf,ispf)+     dirphase*dot(a1,a2,numvects)*csum2
        enddo
     enddo
  enddo !! numr
  
  call mympireduce(reducedr(:,:), www%nspf**2);  call mympireduce(reducedinvr(:,:), www%nspf**2)
  call mympireduce(reducedinvrsq(:,:), www%nspf**2)

end subroutine get_reducedr



