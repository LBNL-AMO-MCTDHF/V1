
!! SUBROUTINES FOR COMPUTING MEAN FIELD !!

!! REDUCEDPOT is usual notation bra, ket   other reduced mats opposite notation for ease of actreducedconjg

#include "Definitions.INC"


subroutine get_allden()
  implicit none
  call get_reducedham();  call getdenmatx();
end subroutine get_allden


subroutine get_reducedham()
  use parameters
  use xxxmod
  implicit none
  integer :: itime,jtime,times(100)

  call system_clock(itime)
  call get_tworeducedx(yyy%reducedpottally,yyy%cmfpsivec(astart(1),0),yyy%cmfpsivec(astart(1),0))
  call system_clock(jtime);  times(4)=times(4)+jtime-itime

  call system_clock(itime)
  call get_reducedproderiv(yyy%reducedproderiv,yyy%cmfpsivec(astart(1),0),yyy%cmfpsivec(astart(1),0))
  call system_clock(jtime);  times(6)=times(6)+jtime-itime

  call system_clock(itime)
  call get_reducedr(yyy%reducedinvr,yyy%reducedinvrsq,yyy%reducedr,yyy%cmfpsivec(astart(1),0),yyy%cmfpsivec(astart(1),0))
  call system_clock(itime);  times(7)=times(7)+itime-jtime

end subroutine get_reducedham


!! IN THE END reducedpottally is bra2,ket2,bra1,ket1 ; I sum it up by conjugating ket2 bra2 ket1 bra1.  Notice a1 is conjg not a2.

subroutine get_reducedpot()
  use xxxmod
  use opmod
  implicit none
  call get_reducedpot0(yyy%reducedpottally(:,:,:,:,0), yyy%reducedpot(:,:,:,0),twoereduced(:,:,:))
end subroutine get_reducedpot


subroutine get_reducedpot0(intwoden,outpot,twoereduced)
  use parameters
  implicit none
  DATATYPE :: TWOEreduced(reducedpotsize, nspf,nspf)
  DATATYPE :: intwoden(nspf,nspf,nspf,nspf),outpot(reducedpotsize,nspf,nspf)  !!,temptwoden(nspf,nspf,nspf,nspf)
  integer :: ii

! I have bra2,ket2,bra1,ket1. (chemists' notation for contraction)

! e.g. M= -2, -1, 2, 1 all different classes.  no zeroed entries even if constraining, in general, I believe.

  ii=reducedpotsize

!  if(deb ugflag.ne.0) then
!     temptwoden(:,:,:,:)=intwoden(:,:,:,:)
!     do ispf=1,nspf
!        do jspf=1,nspf
!           if (orbclass(ispf).ne.orbclass(jspf)) then
!              temptwoden(:,ispf,:,jspf)=0d0
!              temptwoden(ispf,:,jspf,:)=0d0
!           endif
!        enddo
!     enddo
!     call MYGEMM('N','N', ii,nspf**2,nspf**2,(1.0d0,0.d0), twoereduced,ii,temptwoden,nspf**2, (0.d0,0.d0),outpot,ii)
!  else

     call MYGEMM('N','N', ii,nspf**2,nspf**2,(1.0d0,0.d0), twoereduced,ii,intwoden,nspf**2, (0.d0,0.d0),outpot,ii)


end subroutine get_reducedpot0



subroutine get_tworeducedx(reducedpottally,avector1,in_avector2)
  use parameters
  use walkmod
  implicit none

  integer ::   ispf, jspf, iispf, jjspf ,  config2, config1,dirphase, iwalk, qq
  DATATYPE,intent(in) :: avector1(numr,firstconfig:lastconfig,mcscfnum),  in_avector2(numr,firstconfig:lastconfig,mcscfnum)
  DATATYPE,intent(out) :: reducedpottally(nspf,nspf,nspf,nspf)
  DATATYPE :: avector2(numr,numconfig,mcscfnum)             !! AUTOMATIC
  DATATYPE ::  a1(mcscfnum), a2(mcscfnum), dot
  DATAECS :: thisrvalue

  avector2(:,:,:)=0d0
  avector2(:,firstconfig:lastconfig,:)=in_avector2(:,:,:)

!! DO SUMMA
  if (parconsplit.ne.0) then
     do qq=1,mcscfnum
        call mpiallgather(avector2(:,:,qq),numconfig*numr,configsperproc(:)*numr,maxconfigsperproc*numr)
     enddo
  endif

  reducedpottally(:,:,:,:)=0.d0

  !! by configurations; following matel.f90

  do qq=1,numr
     thisrvalue=bondpoints(qq)

        !! doubly off diagonal walks
        
     do config1=botconfig,topconfig

        a1(:)=avector1(qq,config1,:)   !! NO MORE
        
        do iwalk=1,numdoublewalks(config1)
           
           dirphase=doublewalkdirphase(iwalk,config1)
           config2=doublewalk(iwalk,config1)
           a2(:)=avector2(qq,config2,:)
           
           ispf=doublewalkdirspf(1,iwalk,config1)   !BRA2 
           jspf=doublewalkdirspf(2,iwalk,config1)   !KET2 (walk)
           iispf=doublewalkdirspf(3,iwalk,config1)  !BRA1
           jjspf=doublewalkdirspf(4,iwalk,config1)  !KET1 (walk)
           
           reducedpottally(ispf,jspf,iispf,jjspf) =  &    
                reducedpottally(ispf,jspf,iispf,jjspf) +  &
                dirphase*dot(a1,a2,mcscfnum)/thisrvalue
           
           reducedpottally(iispf,jjspf,ispf,jspf) =  &
                reducedpottally(iispf,jjspf,ispf,jspf) + &
                dirphase*dot(a1,a2,mcscfnum)/thisrvalue

        enddo
     enddo   ! config1
  enddo  !! qq

  call mympireduce(reducedpottally(:,:,:,:), nspf**4)

end subroutine get_tworeducedx




subroutine get_reducedproderiv(reducedproderiv,avector1,in_avector2)
  use opmod   !! rkemod, proderivmod
  use parameters
  use walkmod
  implicit none

  DATATYPE,intent(in) :: avector1(numr,firstconfig:lastconfig,mcscfnum), in_avector2(numr,firstconfig:lastconfig,mcscfnum)
  DATATYPE,intent(out) :: reducedproderiv(nspf,nspf)
  DATATYPE :: avector2(numr,numconfig,mcscfnum)     !! AUTOMATIC
  DATATYPE :: a1(mcscfnum), a2(mcscfnum), dot
  integer ::  config1,config2,  ispf,jspf,  dirphase,     iwalk,ii,jj

  avector2(:,:,:)=0d0
  avector2(:,firstconfig:lastconfig,:)=in_avector2(:,:,:)

!! DO SUMMA
  if (parconsplit.ne.0) then
     do ii=1,mcscfnum
        call mpiallgather(avector2(:,:,ii),numconfig*numr,configsperproc(:)*numr,maxconfigsperproc*numr)
     enddo
  endif

  reducedproderiv(:,:)=0.d0

  do config1=botconfig,topconfig

     do iwalk=1,numsinglewalks(config1)
        config2=singlewalk(iwalk,config1);        dirphase=singlewalkdirphase(iwalk,config1)
        ispf=singlewalkopspf(1,iwalk,config1);        jspf=singlewalkopspf(2,iwalk,config1)
        
        do jj=1,numr
           do ii=1,numr
              a1(:)=avector1(ii,config1,:)
              a2(:)=avector2(jj,config2,:)
              
              reducedproderiv(jspf,ispf)=reducedproderiv(jspf,ispf)+ &
                   dirphase*dot(a1,a2,mcscfnum)*proderivmod(ii,jj)
           enddo
        enddo
     enddo
  enddo

  call mympireduce(reducedproderiv(:,:), nspf**2)

end subroutine get_reducedproderiv


subroutine get_reducedr(reducedinvr,reducedinvrsq,reducedr,avector1,in_avector2)
  use parameters
  use walkmod
  implicit none

  DATATYPE,intent(in) :: avector1(numr,firstconfig:lastconfig,mcscfnum), in_avector2(numr,firstconfig:lastconfig,mcscfnum)
  DATATYPE,intent(out) :: reducedinvr(nspf,nspf),reducedr(nspf,nspf),  reducedinvrsq(nspf,nspf)
  DATATYPE :: avector2(numr,numconfig,mcscfnum)      !! AUTOMATIC
  DATATYPE ::  a1(mcscfnum), a2(mcscfnum), dot
  integer ::  config1,config2,   ispf,jspf,  dirphase,    iwalk,ii
  DATAECS :: thisrvalue,  csum,csum2

  reducedinvr(:,:)=0.d0;  reducedr(:,:)=0.d0;  reducedinvrsq(:,:)=0.d0
  if (numr.eq.1) then
     return
  endif

  avector2(:,:,:)=0d0
  avector2(:,firstconfig:lastconfig,:) = in_avector2(:,:,:)
  if (parconsplit.ne.0) then
     do ii=1,mcscfnum
        call mpiallgather(avector2(:,:,ii),numconfig*numr,configsperproc(:)*numr,maxconfigsperproc*numr)
     enddo
  endif

  do ii=1,numr

     thisrvalue=bondpoints(ii);        csum=(1.d0/thisrvalue);        csum2=(1.d0/thisrvalue**2)

     do config1=botconfig,topconfig

        a1(:)=avector1(ii,config1,:)
        
        do iwalk=1,numsinglewalks(config1)
           config2=singlewalk(iwalk,config1);              dirphase=singlewalkdirphase(iwalk,config1)
           a2(:)=avector2(ii,config2,:)
           
           ispf=singlewalkopspf(1,iwalk,config1);              jspf=singlewalkopspf(2,iwalk,config1)
           
           
           reducedinvr(jspf,ispf)=reducedinvr(jspf,ispf)+         dirphase*dot(a1,a2,mcscfnum)*csum
           reducedr(jspf,ispf)=reducedr(jspf,ispf)+         dirphase*dot(a1,a2,mcscfnum)*thisrvalue
           reducedinvrsq(jspf,ispf)=reducedinvrsq(jspf,ispf)+     dirphase*dot(a1,a2,mcscfnum)*csum2
        enddo
     enddo
  enddo !! numr
  
  call mympireduce(reducedr(:,:), nspf**2);  call mympireduce(reducedinvr(:,:), nspf**2)
  call mympireduce(reducedinvrsq(:,:), nspf**2)

end subroutine get_reducedr



