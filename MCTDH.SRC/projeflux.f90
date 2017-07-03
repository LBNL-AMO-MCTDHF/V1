
!! ALL MODULES

!! PROJECTED FLUX (partial photoionization ACTIONS 17 and 28)

#include "Definitions.INC"


module projefluxmod !! needed for cation walks and bi-orthonormalization
  implicit none
  type PPTYPE
     integer :: tnumconfig=0, tnumr=0, tnspf=0
     integer :: maxpwalk1=0 
     integer :: eachstate=0
     real*8 :: catfactor=0d0
     integer, allocatable :: tconfiglist(:,:),pphase1(:,:) ,pspf1(:,:,:) ,numpwalk1(:) ,pwalk1(:,:)
     DATATYPE, allocatable :: tmo(:,:),tavec(:,:,:)
  end type PPTYPE
  integer :: catloaded=0
  type(PPTYPE),allocatable :: ppp(:)
  integer :: nstate=0
  integer, allocatable :: alreadystate(:)

  integer ::     NUMANGLES=0, NUMERAD=0

end module projefluxmod


module projloadsubmod
  implicit none
contains

!! here we are doing walks from our BO target state to our wavefunction \Psi(t)
  subroutine projeflux_singlewalks(ifile)
    use fileptrmod
    use basis_parameters !! numpart, num2part
    use aarrmod
    use projefluxmod
    use configmod
    use mpimod
    use configsubmod
    implicit none
    integer, intent(in) :: ifile
    integer :: iconfig,jconfig,idof,iindex,iwalk,flag,dirphase
    integer :: tempconfig(num2part),temporb(2)

!! get the number of single walks from the target N-1 e- state to our regular N e- state
    OFLWR "Getting cation single walks"; CFL

    allocate(ppp(ifile)%numpwalk1(ppp(ifile)%tnumconfig))
    ppp(ifile)%numpwalk1(:)=0

    do iconfig=1,ppp(ifile)%tnumconfig
       iwalk=0
       do iindex=1,2*www%nspf

          if (www%holeflag.eq.0) then
             tempconfig(3:num2part)=ppp(ifile)%tconfiglist(:,iconfig)
             temporb=aarr(iindex)
             tempconfig(1:2)=temporb(:);        
             flag=0
             do idof=2,numpart 
                if(iind(tempconfig(idof*2-1:idof*2)).eq.iindex) flag=1
             enddo
          else   !! holeflag
             flag=1
             do idof=1,numpart+1
                if(iind(ppp(ifile)%tconfiglist(idof*2-1:idof*2,iconfig)).eq.iindex) then
                   flag=0
                   tempconfig(1:(idof-1)*2) = ppp(ifile)%tconfiglist(1:(idof-1)*2,iconfig)
                   tempconfig(idof*2-1:num2part) = ppp(ifile)%tconfiglist(idof*2+1:num2part+2,iconfig)
                   temporb=aarr(iindex)
                   exit
                endif
             enddo
          endif   !! holeflag
          if(flag.eq.0) then
             dirphase=reorder(tempconfig,www%numpart)
             if(allowedconfig0(www,tempconfig,www%dflevel)) then
                jconfig=getconfiguration(tempconfig,www)
                if (jconfig.ge.www%botconfig.and.jconfig.le.www%topconfig) then
                   iwalk=iwalk+1
                endif
             endif
          endif
       enddo
       ppp(ifile)%numpwalk1(iconfig)=iwalk
    enddo

!! figure out the maximum number of single target walks
    ppp(ifile)%maxpwalk1=0
    do iconfig=1,ppp(ifile)%tnumconfig
       if(ppp(ifile)%maxpwalk1.lt.ppp(ifile)%numpwalk1(iconfig)) ppp(ifile)%maxpwalk1=ppp(ifile)%numpwalk1(iconfig)
    enddo
    OFLWR "Max # single walks from cation state on this processor is ",ppp(ifile)%maxpwalk1;CFL
    call mpibarrier()

    allocate(ppp(ifile)%pphase1(ppp(ifile)%maxpwalk1,ppp(ifile)%tnumconfig),&
         ppp(ifile)%pwalk1(ppp(ifile)%maxpwalk1,ppp(ifile)%tnumconfig),&
         ppp(ifile)%pspf1(2,ppp(ifile)%maxpwalk1,ppp(ifile)%tnumconfig))
    if (ppp(ifile)%maxpwalk1.gt.0) then
       ppp(ifile)%pwalk1=0;  ppp(ifile)%pspf1=0;  ppp(ifile)%pphase1=0
    endif

    do iconfig=1,ppp(ifile)%tnumconfig
       iwalk=0
       do iindex=1,2*www%nspf
          dirphase=(-9999);        tempconfig(:)=(-9999);        temporb(:)=(-9999)
          if (www%holeflag.eq.0) then
             flag=0
             do idof=1,numpart-1
                if(iind(ppp(ifile)%tconfiglist(idof*2-1:idof*2,iconfig)).eq.iindex) then
                   flag=1
                   exit
                endif
             enddo
             if (flag.eq.0) then
                temporb=aarr(iindex)
                tempconfig(1:2)=temporb(:)
                tempconfig(3:num2part)=ppp(ifile)%tconfiglist(:,iconfig)
                dirphase=reorder(tempconfig,www%numpart)
             endif
          else    !! holeflag
             flag=1
             do idof=1,numpart+1
                if(iind(ppp(ifile)%tconfiglist(idof*2-1:idof*2,iconfig)).eq.iindex) then
                   flag=0
                   temporb=aarr(iindex)
                   tempconfig(1:(idof-1)*2) = ppp(ifile)%tconfiglist(1:(idof-1)*2,iconfig)
                   tempconfig(idof*2-1:num2part) = ppp(ifile)%tconfiglist(idof*2+1:num2part+2,iconfig)
                   dirphase=(-1)**idof
                   exit
                endif
             enddo
          endif  !! holeflag
          if(flag.eq.0) then
             if(allowedconfig0(www,tempconfig,www%dflevel)) then
                jconfig=getconfiguration(tempconfig,www)
                if (jconfig.ge.www%botconfig.and.jconfig.le.www%topconfig) then
                   iwalk=iwalk+1;          ppp(ifile)%pwalk1(iwalk,iconfig)=jconfig
                   ppp(ifile)%pspf1(:,iwalk,iconfig)=temporb
                   ppp(ifile)%pphase1(iwalk,iconfig)=dirphase
                endif
             endif
          endif
       enddo
!! check and make sure no bad walk error
       if(ppp(ifile)%numpwalk1(iconfig).ne.iwalk) then
          OFLWR "TARGET SINGLE WALK ERROR";CFLST
       endif
    enddo
  
    call mpibarrier()
    OFLWR "DONE getting proj single walks"; CFL

  end subroutine projeflux_singlewalks


  subroutine projeflux_catload(ifile,outnumstate)
    use parameters    !! catavectorfiles and others
    use configmod
    use projefluxmod
    use mpimod
    use mpisubmod
    use utilmod
    use configloadmod
    use loadstuffmod
    implicit none
    integer,intent(in) :: ifile
    integer,intent(out) :: outnumstate
    integer :: i,ir,tnum2part,istate,ierr, cgflag, &
         spfcomplex, acomplex, tdims(3), myiostat, &
         targetms, targetrestrictflag, targetspinproject, targetspinval
    real*8 :: cgfac,aa,bb,cc
    DATATYPE, allocatable :: tmotemp(:,:),readta(:,:,:),tempta(:,:)

!! read in the data from mcscf for our target cation state

    if (myrank.eq.1) then
       open(909,file=catspffiles(ifile),status="unknown",form="unformatted",iostat=myiostat)
       call checkiostat(myiostat,"opening "//catspffiles(ifile))
       open(910,file=catavectorfiles(ifile),status="unknown",form="unformatted",iostat=myiostat)
       call checkiostat(myiostat,"opening "//catavectorfiles(ifile))
       call avector_header_read(910,ppp(ifile)%eachstate,tnum2part,ppp(ifile)%tnumr,ppp(ifile)%tnumconfig,targetrestrictflag,targetms,&
            targetspinproject,targetspinval,acomplex,ierr)
    endif
    call mympiibcastone(ppp(ifile)%eachstate,1); call mympiibcastone(tnum2part,1); 
    call mympiibcastone(ppp(ifile)%tnumr,1); call mympiibcastone(ppp(ifile)%tnumconfig,1); 
    call mympiibcastone(targetrestrictflag,1);  call mympiibcastone(targetms,1);  
    call mympiibcastone(targetspinproject,1);  call mympiibcastone(targetspinval,1); 
    call mympiibcastone(acomplex,1);  call mympiibcastone(ierr,1)

    if (ierr.ne.0) then
       OFLWR "Avector header read error; redone 09/29/2014; recompute vector on disk."; CFLST
    endif

    if (myrank.eq.1) then
       call spf_header_read(909,tdims,ppp(ifile)%tnspf,spfcomplex)
    endif
    call mympiibcast(tdims,1,3);
    call mympiibcastone(ppp(ifile)%tnspf,1);
    call mympiibcastone(spfcomplex,1)

!! have to project on BO wfns.  Otherwise doesn't make sense in prolate.  
!! Not supported anymore, will support
!!   again, mcscf mode w/many r's on file

    if (ppp(ifile)%tnspf.gt.nspf+numfrozen) then
       OFLWR "ERROR, for now can't do more orbs in projection than in calculation",ppp(ifile)%tnspf,nspf+numfrozen; CFLST
    endif
    if ( (holeflag.eq.0.and.tnum2part.ne.num2part-2) .or. (holeflag.ne.0.and.tnum2part.ne.num2part+2) ) then
       OFLWR "Vectors read are not n-1 electron functions"; CFLST
    endif

    cgflag=0
    if (restrictflag.eq.0) then
       OFLWR " WARNING: propagated state is not spin projection restricted.  No CG algebra."
       CFL
    else
       if (targetrestrictflag.eq.0) then
          OFLWR "WARNING: N-1 electron state is not spin projection restricted (restrictflag=0).  No CG algebra."
          CFL
       else
          if (abs(targetms-restrict_ms).ne.1) then
             OFLWR "Targetms should differ from restrictms by 1.  Cation and neutral have no dyson orbital"
             WRFL  "    Targetms=", targetms, " restrictms=",restrict_ms; CFLST
          endif
          if (all_spinproject.eq.0) then
             OFLWR "WARNING: Wave function is not spin restricted (allspinproject=0).  No CG algebra.";CFL
          else
             if (targetspinproject.eq.0) then
                OFLWR "WARNING: N-1 electron state is not spin restricted (allspinproject=0).  No CG algebra."
                CFL
             else
                if (abs(targetspinval-spin_restrictval).gt.1) then
                   OFLWR "ERROR: Spin value of wave function and N-1 electron state differ by more than 1/2.",&
                        targetspinval,spin_restrictval; CFLST
                else
                   OFLWR "Spin adapted wave functions for propagation and projection. doing CG algebra.";CFL
                   cgflag=1
                endif
             endif
          endif
       endif
    endif

!!!!  100912 NEED CG COEFFICIENT RATIO 
              
!! programmed this with restrict_ms, targetms.  
!! Do have option to change spinrestrictval; should use this.
!!
!!    have high spin cation and high spin neutral.
!!
!!    | s s  >  =  sum_m_s'  < s s   | s' s' m_s'  >  | s' m_s' >  x  | 1 +/-1/2 >
!!
!!     neutral                                          cation       outgoing electron
!!     m_s = s
!!
!!    only have |s' s' > m_s'=s'  cation not |s' s'-1> m_s'=s'-1
!!
!!   
!!    if restrict_ms > targetms  then we couple high spin to high spin to get high spin, one CG coef = 1.
!!    
!!        otherwise we couple high spin target down to high spin neutral:
!!
!!         projection is spin down  x high spin target   e.g.   < s' s' |  s s  1/2 -1/2 > = A
!!         also have < s' s' | s s-1 1/2 1/2 > = B
!!  
!!         so mult cross section by ( A^2 + B^2 / A^2 )
!!    
!!   BUT if restrictflag is off, can't do this.  Must realize that you are only projecting one 
!!     component of C.G. sum.  that would be a problem.  But if you were to calculate say the 
!!     3 lowest states for a triplet, you'd get all components.
!!     otherwise cross section is meaningless.
!! all arguments are x 2, integers for half spin.

!!$                 if (restrict_ms.lt.targetms.and.restrictflag.eq.1) then
!!$  !! what we're calculating  targetms-restrict_ms=+/-1
!!$                    aa= doubleclebschsq(targetms,1,targetms,restrict_ms-targetms,restrict_ms) 
!!$  !! only one could be nonzero obviously
!!$                   bb= doubleclebschsq(targetms,1,targetms-2,restrict_ms-targetms+2,restrict_ms)  
!!$                    cc= doubleclebschsq(targetms,1,targetms+2,restrict_ms-targetms-2,restrict_ms)      
!!$                    cgfac = (aa+bb+cc)/aa
!!$                 else
!!$                    cgfac=1
!!$                 endif

    if (cgflag.eq.0) then
       cgfac=1
    else
!! what we're calculating  targetms-restrict_ms=+/-1
       aa= doubleclebschsq(targetspinval,1,targetms,restrict_ms-targetms,spin_restrictval)
!! only one could be nonzero obviously
       bb= doubleclebschsq(targetspinval,1,targetms-2,restrict_ms-targetms+2,spin_restrictval)
       cc= doubleclebschsq(targetspinval,1,targetms+2,restrict_ms-targetms-2,spin_restrictval)      
       cgfac = (aa+bb+cc)/aa
    endif

    ppp(ifile)%catfactor = cgfac * catfacs(ifile)

    allocate(ppp(ifile)%tmo(spfsize,nspf),ppp(ifile)%tavec(numr,ppp(ifile)%tnumconfig,ppp(ifile)%eachstate))
    ppp(ifile)%tmo=0d0;  ppp(ifile)%tavec=0d0; 

    allocate(tmotemp(spfsize,nspf+numfrozen),readta(ppp(ifile)%tnumr,ppp(ifile)%tnumconfig,ppp(ifile)%eachstate), &
         tempta(ppp(ifile)%tnumconfig,ppp(ifile)%eachstate))
    tmotemp=0d0; readta=0d0; tempta=0d0

    OFLWR "Reading", ppp(ifile)%eachstate," Born-Oppenheimer states."; CFL

    if (myrank.eq.1) then
       call simple_load_avectors(910,acomplex,readta(:,:,:),tnum2part,ppp(ifile)%tnumr,ppp(ifile)%tnumconfig,ppp(ifile)%eachstate)
       do ir=1,min(ppp(ifile)%tnumr,numr)
          ppp(ifile)%tavec(ir,:,:)=readta(ir,:,:)
       enddo
!!$
!!$ 07-2016 commenting this out and just normalizing up to tnumr... might want to redo this if it fails
!!$
!!$       do ir=min(ppp(ifile)%tnumr,numr)+1,numr
!!$          call staticvector(tempta(:,:),ppp(ifile)%tnumconfig*ppp(ifile)%eachstate)
!!$          ppp(ifile)%tavec(ir,:,:)=tempta(:,:)
!!$       enddo
!!$
!!$       do ir=1,numr

       do ir=1,min(ppp(ifile)%tnumr,numr)

!! projecting on normalized electronic wfn at each R
!! we go to war with the army we've got

          tempta(:,:)=ppp(ifile)%tavec(ir,:,:)
          do istate=1,ppp(ifile)%eachstate
             tempta(:,istate)=tempta(:,istate)/&
                  sqrt(dot(tempta(:,istate),tempta(:,istate),ppp(ifile)%tnumconfig)) !! no * bondweights(ir)
          enddo
          ppp(ifile)%tavec(ir,:,:)=tempta(:,:)
       enddo
       call spf_read0(909,nspf+numfrozen,spfdims,ppp(ifile)%tnspf,tdims,spfcomplex,spfdimtype,&
            tmotemp(:,:),(/0,0,0/))
       close(909)
       close(910)
    else
       call spf_read0(-42,nspf+numfrozen,spfdims,ppp(ifile)%tnspf,tdims,spfcomplex,spfdimtype,&
            tmotemp(:,:),(/0,0,0/))
    endif
    call mympibcast(ppp(ifile)%tavec(:,:,:),1,ppp(ifile)%tnumconfig*ppp(ifile)%eachstate*numr)

    ppp(ifile)%tmo(:,1:ppp(ifile)%tnspf-numfrozen) = tmotemp(:,numfrozen+1:ppp(ifile)%tnspf)
    ppp(ifile)%tnspf=ppp(ifile)%tnspf-numfrozen

    do i=ppp(ifile)%tnspf+1,nspf
       ppp(ifile)%tmo(:,i)=0d0
       call staticvector(ppp(ifile)%tmo(:,i),spfsize)
       if (parorbsplit.eq.3) then
          call gramschmidt(spfsize,i-1,spfsize,ppp(ifile)%tmo(:,:),ppp(ifile)%tmo(:,i),.true.)
       else
          call gramschmidt(spfsize,i-1,spfsize,ppp(ifile)%tmo(:,:),ppp(ifile)%tmo(:,i),.false.)
       endif
    enddo

    allocate(ppp(ifile)%tconfiglist(tnum2part,ppp(ifile)%tnumconfig))
    ppp(ifile)%tconfiglist=0

    if (myrank.eq.1) then
       open(910,file=catavectorfiles(ifile),status="unknown",form="unformatted",iostat=myiostat)
       call checkiostat(myiostat,"opening "//catavectorfiles(ifile))
       call avector_header_read_simple(910,ppp(ifile)%eachstate,tnum2part,ppp(ifile)%tnumr,ppp(ifile)%tnumconfig,acomplex)
       call get_avectorfile_configlist(910,acomplex,ppp(ifile)%tconfiglist,tnum2part,ppp(ifile)%tnumr,ppp(ifile)%tnumconfig)
       close(910)
    endif

    call mympiibcast(ppp(ifile)%tconfiglist,1,tnum2part*ppp(ifile)%tnumconfig)

!! do the walks from the target state into our final state

    call projeflux_singlewalks(ifile)

    deallocate(tmotemp,readta,tempta)

    outnumstate=ppp(ifile)%eachstate

    call mpibarrier()

  end subroutine projeflux_catload

!! load propagated cation states for strong field ionization if strongcatflag.ne.0

  subroutine projeflux_strongcatload(ifile,curtime)
    use parameters    !! catavectorfiles and others
    use configmod
    use projefluxmod
    use mpimod
    use mpisubmod
    use orbmultsubmod   !! gauge_transform
    use utilmod
    implicit none
    integer,intent(in) :: ifile, curtime
    real*8 :: cattime
    integer :: istate,myiostat,catmolength, catalength, readtime, ispf, ir,&
         prevstate, jfile
    DATATYPE, allocatable :: tmotemp(:,:),readta(:,:,:),tempta(:,:), transxmo(:,:)
    DATATYPE :: csum
    integer, save :: numcatsavesteps(100) = -1

    if (strongcatflag.eq.0) then
       OFLWR "programmer fail strongcatflag.eq.0"; CFLST
    endif

    prevstate=0
    do jfile=1,ifile-1
       prevstate=prevstate+ppp(jfile)%eachstate
    enddo

    if (catenergies(prevstate+ppp(ifile)%eachstate).eq.(0d0,0d0)) then
       OFLWR "For strongcatflag, need catenergies for at least ",&
            prevstate+ppp(ifile)%eachstate," states"; CFLST
    endif

    allocate(readta(ppp(ifile)%tnumr,ppp(ifile)%tnumconfig,ppp(ifile)%eachstate), &
         tempta(ppp(ifile)%tnumconfig,ppp(ifile)%eachstate))

    if (myrank.eq.1) then
       if (parorbsplit.eq.3) then
          allocate(tmotemp(spfsize*nprocs,ppp(ifile)%tnspf))
       else
          allocate(tmotemp(spfsize,ppp(ifile)%tnspf))
       endif
    else
       allocate(tmotemp(1,ppp(ifile)%tnspf))
    endif

    tmotemp=0d0; readta=0d0; tempta=0d0

!!    OFLWR "Reading", ppp(ifile)%eachstate," Born-Oppenheimer states."; CFL

    if (myrank.eq.1) then
       inquire (iolength=catmolength) tmotemp(:,:)
       inquire (iolength=catalength) readta(:,:,:)
    endif
    call mympiibcastone(catmolength,1); call mympiibcastone(catalength,1)

!!    readtime = floor(real(numpropsteps,8)/fluxinterval) - curtime*fluxskipmult
!!    if (readtime.lt.0) then
!!       OFLWR "programmer fail readtime",readtime,curtime, fluxskipmult; CFLST
!!    endif

    if (numcatsavesteps(ifile).eq.(-1)) then
       if (myrank.eq.1) then
          open(10011,file=strongcatspffiles(ifile),status="old",form="unformatted",&
               access="direct",recl=catmolength,iostat=myiostat)
          call checkiostat(myiostat,"opening "//strongcatspffiles(ifile))
          numcatsavesteps(ifile) = -1
          myiostat=0
          do while (myiostat.eq.0)
             numcatsavesteps(ifile)=numcatsavesteps(ifile)+1
             read(10011,rec=numcatsavesteps(ifile)+1,iostat=myiostat) csum
          enddo
          close(10011)
          numcatsavesteps(ifile)=numcatsavesteps(ifile)-1
       endif
       call mympiibcastone(numcatsavesteps(ifile),1)
       OFLWR "     ...numcatsavesteps = ",numcatsavesteps(ifile), ifile; CFL
!!       OFLWR "     ... whereas numpropsteps = ",numpropsteps; CFL
    endif

    readtime = numcatsavesteps(ifile) - curtime*fluxskipmult
    if (readtime.lt.0) then
       readtime=0
    endif

    cattime = (numcatsavesteps(ifile) - readtime) * fluxinterval * par_timestep

!! ORBITALS

    ppp(ifile)%tmo=0d0

    if(myrank.eq.1) then
       open(10011,file=strongcatspffiles(ifile),status="old",form="unformatted",&
            access="direct",recl=catmolength,iostat=myiostat)
       call checkiostat(myiostat,"opening "//strongcatspffiles(ifile))
       read(10011,rec=readtime+1,iostat=myiostat) tmotemp(:,:) 
       call checkiostat(myiostat,"reading "//strongcatspffiles(ifile))
       close(10011)
    endif
    if (parorbsplit.ne.3) then
       if (myrank.eq.1) then
          ppp(ifile)%tmo(:,1:ppp(ifile)%tnspf)=tmotemp(:,:)
       endif
       call mympibcast(ppp(ifile)%tmo(:,:),1,spfsize*ppp(ifile)%tnspf)
    else
       do ispf=1,ppp(ifile)%tnspf
          call splitscatterv(tmotemp(:,ispf),ppp(ifile)%tmo(:,ispf))
       enddo
    endif

    do ispf=ppp(ifile)%tnspf+1,nspf
       ppp(ifile)%tmo(:,ispf)=0d0
       call staticvector(ppp(ifile)%tmo(:,ispf),spfsize)
       if (parorbsplit.eq.3) then
          call gramschmidt(spfsize,ispf-1,spfsize,ppp(ifile)%tmo(:,:),ppp(ifile)%tmo(:,ispf),.true.)
       else
          call gramschmidt(spfsize,ispf-1,spfsize,ppp(ifile)%tmo(:,:),ppp(ifile)%tmo(:,ispf),.false.)
       endif
    enddo

!! gaugefluxflag available to attempt gauge-invariant calculation with photoionization
!! during the pulse.  Best performance seems to be calculate wfn in length gauge, transform
!! to velocity gauge for flux: velflag=0, gaugefluxflag=2 (Volkov phase is only a function of 
!! time in the velocity gauge, no preferred origin in the velocity gauge)

!! gaugefluxflag=1, transform velocity to length; gaugefluxflag=2, transform length to velocity
    if ((gaugefluxflag.eq.1.and.velflag.ne.0).or.(gaugefluxflag.eq.2.and.velflag.eq.0)) then
       allocate(transxmo(spfsize,nspf))
       transxmo=0d0
       call gauge_transform(velflag,curtime*FluxInterval*FluxSkipMult*par_timestep,&
            nspf,ppp(ifile)%tmo(:,:),transxmo(:,:))
       ppp(ifile)%tmo(:,:) = transxmo(:,:)
       deallocate(transxmo)
    endif

!! A-VECTOR
    if(myrank.eq.1) then
       open(10022,file=strongcatavectorfiles(ifile),status="old",form="unformatted",&
            access="direct",recl=catalength,iostat=myiostat)
       call checkiostat(myiostat,"opening "//strongcatavectorfiles(ifile))
       read(10022,rec=readtime+1,iostat=myiostat) readta(:,:,:)
       call checkiostat(myiostat,"reading "//strongcatavectorfiles(ifile))
       close(10022)

       ppp(ifile)%tavec=0d0 
       do ir=1,min(ppp(ifile)%tnumr,numr)
          ppp(ifile)%tavec(ir,:,:)=readta(ir,:,:)
       enddo

!! catenergies here

       do istate=1,ppp(ifile)%eachstate
          ppp(ifile)%tavec(:,:,istate) = ppp(ifile)%tavec(:,:,istate) &
               * exp((0d0, 1d0)*catenergies(prevstate+istate)*cattime)
       enddo

!!$
!!$ 07-2016 commenting this out and just normalizing up to tnumr... might want to redo this if it fails
!!$
!!$       do ir=min(ppp(ifile)%tnumr,numr)+1,numr
!!$          call staticvector(tempta(:,:),ppp(ifile)%tnumconfig*ppp(ifile)%eachstate)
!!$          ppp(ifile)%tavec(ir,:,:)=tempta(:,:)
!!$       enddo
!!$
!!$       do ir=1,numr

       do ir=1,min(ppp(ifile)%tnumr,numr)

!! projecting on normalized electronic wfn at each R
!! we go to war with the army we've got

          tempta(:,:)=ppp(ifile)%tavec(ir,:,:)
          do istate=1,ppp(ifile)%eachstate
             tempta(:,istate)=tempta(:,istate)/&
                  sqrt(dot(tempta(:,istate),tempta(:,istate),ppp(ifile)%tnumconfig)) !! no * bondweights(ir)
          enddo
          ppp(ifile)%tavec(ir,:,:)=tempta(:,:)
       enddo
    endif

    deallocate(tmotemp,readta,tempta)

    call mympibcast(ppp(ifile)%tavec(:,:,:),1,ppp(ifile)%tnumconfig*ppp(ifile)%eachstate*numr)

  end subroutine projeflux_strongcatload

end module projloadsubmod


module projutilsubmod
  implicit none
contains

  subroutine only_one_angle(howmany,iangle,inout,NUMANGLES)
    use spfsize_parameters 
    use fileptrmod
    integer,intent(in) :: howmany,iangle,NUMANGLES
    DATATYPE,intent(inout) :: inout(spfdims(1),spfdims(2),spfdims(3),howmany)
    if (NUMANGLES.ne.spfdims(2)) then
       OFLWR "WHAAAAAAATTTT?? onlyoneangle",spfdims(2),NUMANGLES; CFLST
    endif
    inout(:,1:iangle-1,:,:)=0
    inout(:,iangle+1:NUMANGLES,:,:)=0
  end subroutine only_one_angle

  function myhermdots(bra,ket,totsize,NUMERAD,NUMANGLES)
    use fileptrmod
    use dotmod
    integer,intent(in) :: totsize,NUMERAD,NUMANGLES
    DATATYPE,intent(in) :: bra(NUMERAD,NUMANGLES,totsize/NUMERAD/(NUMANGLES)),&
         ket(NUMERAD,NUMANGLES,totsize/NUMERAD/(NUMANGLES))
    DATATYPE :: myhermdots(NUMANGLES)
    DATATYPE :: tempbra(NUMERAD,totsize/NUMERAD/(NUMANGLES)), &
         tempket(NUMERAD,totsize/NUMERAD/(NUMANGLES)),&            !! AUTOMATIC
         temparray(NUMANGLES)
    integer :: ibb,il

    ibb=(totsize/NUMERAD/(NUMANGLES)) * NUMERAD * (NUMANGLES)
    if (ibb.ne.totsize) then
       OFLWR "MYHERMDOTS ERROR ",ibb,totsize,NUMERAD,NUMANGLES; CFLST
    endif

    temparray=0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(il,tempbra,tempket)
!$OMP DO SCHEDULE(DYNAMIC)
    do il=1,NUMANGLES
       tempbra(:,:)=bra(:,il,:);    tempket(:,:)=ket(:,il,:)
       temparray(il)=hermdot(tempbra,tempket,totsize/(NUMANGLES))
    enddo
!$OMP END DO
!$OMP END PARALLEL
    myhermdots(:) = temparray(:)

  end function myhermdots

!! get the contraction of the flux operator (iH-H^\dag) with our current set of orbitals 
!! F*inspfs=outspfs

  subroutine projeflux_op_onee(inspfs,outspfs)
    use parameters   !! nucfluxopt
    use opmod
    use fluxutilmod
    implicit none
    DATATYPE,intent(in) :: inspfs(spfsize,numr)
    DATATYPE,intent(out) :: outspfs(spfsize,numr)
    DATATYPE :: workspfs(spfsize,numr), mymat(numr,numr)    !! AUTOMATIC
    integer :: r

    workspfs=0; outspfs=0; mymat=0

    if (nucfluxopt.ne.2) then

!! terms from antihermitian part of B-O hamiltonian and hermitian part of Y operator
!! which are nonzero in the external region for electrons and so characterize the
!! electronic flux

!! terms from antihermitian of B-O hamiltonian 

       call mult_imke(numr,inspfs(:,:),workspfs(:,:))
       do r=1,numr
          outspfs(:,r) = outspfs(:,r) + workspfs(:,r) * real(1d0 / bondpoints(r)**2,8)
       enddo
    
       if(FluxOpType.eq.0.or.FluxOpType.eq.2) then
          call mult_impot(numr,inspfs(:,:),workspfs(:,:))
       else if(FluxOpType.eq.1) then
          call mult_imhalfniumpot(numr,inspfs(:,:),workspfs(:,:))
       endif
       do r=1,numr
          outspfs(:,r) = outspfs(:,r) + workspfs(:,r) * real(1d0 / bondpoints(r),8)
       enddo

       if (nonuc_checkflag.eq.0) then

!! hermitian part of antisymmetric operator is imaginary valued
          call op_imyderiv(numr,inspfs,workspfs)

!! antihermitian part of antisymmetric operator is real valued
          mymat(:,:)=real(proderivmod(:,:),8)

!! transpose verified 2nd argument 'N' not 'T' (with 'T' results are way off)
          call MYGEMM('N','N',spfsize,numr,numr,DATAONE,workspfs,spfsize,mymat,numr,&
               DATAONE,outspfs,spfsize)

       endif  !! nonuc_checkflag

    endif   !! nucfluxopt.ne.2

    if (nonuc_checkflag.eq.0.and.nucfluxopt.ne.0) then

!! terms from imaginary part of inverse bond length, hermitian part of Y operator, and
!! antihermitian part of nuclear ke second derivative operator, all of which are 
!! nonzero only in external region for nuclei and so characterize the nuclear flux

       call mult_reke(numr,inspfs(:,:),workspfs(:,:))
       do r=1,numr
          outspfs(:,r) = outspfs(:,r) + workspfs(:,r) * imag((0d0,0d0)+ 1d0 / bondpoints(r)**2)
       enddo
    
       if(FluxOpType.eq.0.or.FluxOpType.eq.2) then
          call mult_repot(numr,inspfs(:,:),workspfs(:,:))
       else if(FluxOpType.eq.1) then
          call mult_rehalfniumpot(numr,inspfs(:,:),workspfs(:,:))
       endif
       do r=1,numr
          outspfs(:,r) = outspfs(:,r) + workspfs(:,r) * imag((0d0,0d0)+ 1d0 / bondpoints(r))
       enddo

!! real part of antisymmetric operator is antihermitian
       call op_reyderiv(numr,inspfs,workspfs)

!! imag part is hermitian
       mymat(:,:)=imag(proderivmod(:,:)+(0d0,0d0))

!! transpose verified 2nd argument 'N' not 'T' (with 'T' results are way off)
       call MYGEMM('N','N',spfsize,numr,numr,DATAONE,workspfs,spfsize,mymat,numr,&
            DATAONE,outspfs,spfsize)

       mymat(:,:)=imag(rkemod(:,:)+(0d0,0d0))
       call MYGEMM('N','N',spfsize,numr,numr,DATAONE,inspfs,spfsize,mymat,numr,&
            DATAONE,outspfs,spfsize)

    endif

!! scale correctly
    outspfs(:,:)=outspfs(:,:)*(-2d0)

  end subroutine projeflux_op_onee

end module projutilsubmod


module projcoresubmod
  implicit none
contains

  subroutine projeflux_load()
    use parameters    !! numcatfiles
    use projefluxmod
    use projloadsubmod
    implicit none
    integer :: ifile,eachstate

    NUMANGLES=0; NUMERAD=0

    if (angularflag.ne.0) then    !! FULLY DIFFERENTIAL PHOTOIONIZATION (function of angle and energy)
       if (spfdimtype(2).ne.2) then
          OFLWR "For angularflag, fully differential photoionizatoin, must have spfdimtype(2).eq.2"
          WRFL "     (atom/diatom only)"; CFLST
       endif

       NUMANGLES=spfdims(2)
       NUMERAD=spfdims(1)

    endif

    if (catloaded.ne.0) then
       OFLWR "error, cations loaded",catloaded; CFLST
    endif
    catloaded=1

    if (nstate.ne.0) then
       OFLWR "error, cation wave functions already loaded? ", nstate; CFLST
    endif

    allocate(ppp(numcatfiles), alreadystate(numcatfiles))
    alreadystate=0

    nstate=0
    do ifile=1,numcatfiles
       alreadystate(ifile)=nstate
       call projeflux_catload(ifile,eachstate)
       nstate=nstate+eachstate
    enddo

  end subroutine projeflux_load


!! construct the one electron functions
  subroutine projeflux_doproj(ifile,cata,neuta,mo,infac,projwfn)
    use parameters     !! nspf, projfluxfile etc.
    use projefluxmod
    use mpimod
    use mpisubmod
    implicit none
    integer,intent(in) :: ifile
    DATATYPE,intent(in) :: cata(numr,ppp(ifile)%tnumconfig),&
         neuta(numr,first_config:last_config),mo(spfsize,nspf)
    real*8, intent(in) :: infac
    DATATYPE,intent(out) :: projwfn(spfsize,2,numr)
    DATATYPE :: projcoefs(numr,nspf,2)             !! AUTOMATIC
    integer :: jconfig,iwalk,iconfig,ispf,ispin,iphase,ir

!! make the single electron wfn

    projcoefs(:,:,:)=0d0
    do jconfig=1,ppp(ifile)%tnumconfig
       do iwalk=1,ppp(ifile)%numpwalk1(jconfig)
          iconfig=ppp(ifile)%pwalk1(iwalk,jconfig);      ispf=ppp(ifile)%pspf1(1,iwalk,jconfig)
          ispin=ppp(ifile)%pspf1(2,iwalk,jconfig);      iphase=ppp(ifile)%pphase1(iwalk,jconfig)

          projcoefs(:,ispf,ispin)=projcoefs(:,ispf,ispin) + &
               CONJUGATE(cata(:,jconfig)) * neuta(:,iconfig) * iphase
       enddo
    enddo

    call mympireduce(projcoefs,numr*nspf*2)

    projwfn(:,:,:)=0d0
    do ir=1,numr
       do ispin=1,2
          do ispf=1,nspf

!! sqrt(infac) to recover infac with bra and ket factor.  infac positive

             projwfn(:,ispin,ir) = projwfn(:,ispin,ir) + mo(:,ispf) * projcoefs(ir,ispf,ispin) * sqrt(infac)
          enddo
       enddo
    enddo

  end subroutine projeflux_doproj


  subroutine projeflux_allproj(ifile,catavec,neutavec,neutmo,projwfn)
    use parameters !! mcscfnum
    use projefluxmod
    implicit none
    integer,intent(in) :: ifile
    DATATYPE, intent(in) :: catavec(numr,ppp(ifile)%tnumconfig,ppp(ifile)%eachstate),&
         neutmo(spfsize,nspf),neutavec(numr,first_config:last_config,mcscfnum)
    DATATYPE,intent(out) :: projwfn(spfsize,2,numr,ppp(ifile)%eachstate,mcscfnum)
    integer :: imc,istate
    DATATYPE :: nullvector(numr)

    nullvector=0

    do imc=1,mcscfnum
       do istate=1,ppp(ifile)%eachstate

   if (tot_adim.gt.0) then
      call projeflux_doproj(ifile,catavec(:,:,istate),neutavec(:,:,imc),neutmo(:,:),ppp(ifile)%catfactor,projwfn(:,:,:,istate,imc))
   else
      call projeflux_doproj(ifile,catavec(:,:,istate),nullvector(:),neutmo(:,:),ppp(ifile)%catfactor,projwfn(:,:,:,istate,imc))
   endif

      enddo
   enddo

  end subroutine projeflux_allproj


  subroutine projeflux_saveproj(ifile,nt,tau,projwfn)
    use parameters !! mcscfnum
    use projefluxmod
    use mpimod
    implicit none
    integer,intent(in) :: ifile,nt,tau
    DATATYPE,intent(in) :: projwfn(spfsize,2,numr,ppp(ifile)%eachstate,mcscfnum)
    DATATYPE,allocatable :: bigprojwfn(:,:)
    integer :: offset,myiostat,mylength,ir,imc,istate,ispin

    if (myrank.eq.1) then
       if (parorbsplit.eq.3) then
          allocate(bigprojwfn(spfsize*nprocs,2));        bigprojwfn(:,:)=0d0
       else
          allocate(bigprojwfn(spfsize,2))
       endif
    else
       allocate(bigprojwfn(1,2));   bigprojwfn(:,:)=0d0
    endif
   
    if (myrank.eq.1) then
       inquire (iolength=mylength) bigprojwfn
       open(1003,file=projfluxfile,status="unknown",form="unformatted",&
            access="direct",recl=mylength,iostat=myiostat)
       call checkiostat(myiostat,"opening projfluxfile")
    endif

    do ir=1,numr
       do imc=1,mcscfnum
          do istate=1,ppp(ifile)%eachstate
             
             if (parorbsplit.eq.3) then
                do ispin=1,2
                   call splitgatherv(projwfn(:,ispin,ir,istate,imc),bigprojwfn(:,ispin),.false.)
                enddo
             else
                bigprojwfn(:,:)=projwfn(:,:,ir,istate,imc)
             endif

             offset=(istate-1+alreadystate(ifile))*mcscfnum*(nt+1)*numr + (imc-1)*(nt+1)*numr + tau*numr + ir

             if (myrank.eq.1) then
                write(1003,rec=offset,iostat=myiostat) bigprojwfn(:,:) 
                call checkiostat(myiostat,"writing projfluxfile")
             endif
          enddo
       enddo
    enddo
    if (myrank.eq.1) then
       close(1003)
    endif

    deallocate(bigprojwfn)

    call mpibarrier()

  end subroutine projeflux_saveproj



  subroutine projeflux_projectdisk(nt,dt)
    use parameters    !! catavectorfiles and others
    use configmod
    use projefluxmod   !! ppp(ifile)%eachstate
    use mpimod
    use mpisubmod
    implicit none
!! necessary working variables
    integer,intent(in) :: nt
    real*8,intent(in) :: dt
    DATATYPE, allocatable ::  mobio(:,:),abio(:,:,:),mymo(:,:),myavec(:,:,:), &
         readmo(:,:),readavec(:,:,:)
    DATATYPE :: nullvector(numr)
!! mcscf specific read variables
!    DATATYPE,target :: smo(nspf,nspf)
    integer :: tau, i, ifile,  myiostat
    DATATYPE,allocatable :: projwfn(:,:,:,:,:)

    nullvector(:)=0

!! allocate all necessary extra memory and io params to do this looping business

    OFLWR "Allocating arrays for projected flux"; CFL

    allocate(mobio(spfsize,nspf),abio(numr,first_config:last_config,mcscfnum),&
         mymo(spfsize,nspf),  myavec(numr,first_config:last_config,mcscfnum))
    mobio=0; mymo=0
    if (tot_adim.gt.0) then
       abio=0; myavec=0
    endif
    if (myrank.eq.1) then
       if (parorbsplit.eq.3) then
          allocate(readmo(spfsize*nprocs,nspf))
       else
          allocate(readmo(spfsize,nspf))
       endif
       allocate(readavec(numr,num_config,mcscfnum))
    else
       allocate(readmo(1,nspf),readavec(1,numr,mcscfnum))
    endif
    readmo=0; readavec=0

!! IFILE LOOP

    do ifile=1,numcatfiles

       allocate(projwfn(spfsize,2,numr,ppp(ifile)%eachstate,mcscfnum))
       projwfn=0

       call mpibarrier()
       OFLWR "    ...go time loop", ifile; CFL

       if (myrank.eq.1) then
          inquire (iolength=i) readmo
          open(1001,file=fluxmofile,status="unknown",form="unformatted",&
               access="direct",recl=i,iostat=myiostat)
          call checkiostat(myiostat,"opening "//fluxmofile)
          inquire (iolength=i) readavec
          open(1002,file=fluxafile,status="unknown",form="unformatted",&
               access="direct",recl=i,iostat=myiostat)
          call checkiostat(myiostat,"opening "//fluxafile)
       endif

!! do the loop over all time, ALL TIME and now in parallel-o-vision

       do tau=0,nt    
!! read in this time's wavefucntion

          if (myrank.eq.1) then
             read(1001,rec=FluxSkipMult*tau+1,iostat=myiostat) readmo(:,:)
             call checkiostat(myiostat,"reading "//fluxmofile)
             read(1002,rec=FluxSkipMult*tau+1,iostat=myiostat) readavec(:,:,:)
             call checkiostat(myiostat,"reading "//fluxafile)
          endif
          if (parorbsplit.ne.3) then
             if (myrank.eq.1) then
                mymo(:,:)=readmo(:,:)
             endif
             call mympibcast(mymo(:,:),1,totspfdim)
          else
             do i=1,nspf
                call splitscatterv(readmo(:,i),mymo(:,i))
             enddo
          endif
          if (par_consplit.eq.0) then
             if (myrank.eq.1) then
                myavec(:,:,:)=readavec(:,:,:)
             endif
             call mympibcast(myavec(:,:,:),1,numr*num_config*mcscfnum)
          else
             do i=1,mcscfnum
                if (tot_adim.gt.0) then
                   call myscatterv(readavec(:,:,i),myavec(:,:,i),configs_perproc(:)*numr)
                else
                   call myscatterv(readavec(:,:,i),nullvector(:),configs_perproc(:)*numr)
                endif
             enddo
          endif

          call projeflux_projectone(ifile, mymo, tau, dt, myavec, projwfn)

          call projeflux_saveproj(ifile,nt,tau,projwfn)


!! write out times cause we're bored
!    if(mod(tau,100).eq.0.or.tau.eq.nt) then
!      call openfile
!      write(mpifileptr,'(A52,F10.4)') " Timing Statistics: producing 1e- proj wfn as of T= ",tau*dt
!      write(mpifileptr,'(100A10)') "Times: ", "All", "IO", "Walks", "Biorth", "1e Bld", "1e Op", &
!          "Fluxeval", "FT g!tau"
!      write(mpifileptr,'(A10,100I10)') " ", times(1:8)/100
!    endif

       enddo   !! DO TAU=0,NT

       deallocate(projwfn)

    enddo  !! DO IFILE

!! clean up
    if (myrank.eq.1) then
       close(1001);  close(1002)
    endif
    call mpibarrier()
    deallocate(readmo,readavec)
    deallocate(mobio,abio,mymo,myavec)

  end subroutine projeflux_projectdisk


  subroutine projeflux_projectone(ifile,inspfs,spftime,dt,inavectors,projwfn)
    use parameters    !! catavectorfiles and others
    use biorthomod
    use biorthotypemod
    use configmod
    use projefluxmod
    use mpimod
    use mpisubmod
    use orbmultsubmod   !! gauge_transform
    use projloadsubmod
    implicit none
    integer,intent(in) :: ifile, spftime
    real*8,intent(in) :: dt
    DATATYPE,intent(in) :: inspfs(spfsize,nspf), inavectors(numr,first_config:last_config,mcscfnum)
    DATATYPE, intent(out) :: projwfn(spfsize,2,numr,ppp(ifile)%eachstate,mcscfnum)
    DATATYPE, allocatable ::  mobio(:,:),myavec(:,:,:), catspfs(:,:), neutspfs(:,:)
    type(biorthotype),target :: projbiovar(numcatfiles)
    DATATYPE,target :: smo(nspf,nspf)
    DATATYPE :: nullvector(numr)
    integer :: imc

    nullvector(:)=0; smo=0d0

    allocate(mobio(spfsize,nspf),catspfs(spfsize,nspf),neutspfs(spfsize,nspf))
    mobio=0; catspfs=0; neutspfs=0

!! gaugefluxflag available to attempt gauge-invariant calculation with photoionization
!! during the pulse.  Best performance seems to be calculate wfn in length gauge, transform
!! to velocity gauge for flux: velflag=0, gaugefluxflag=2 (Volkov phase is only a function of 
!! time in the velocity gauge, no preferred origin in the velocity gauge)

!! gaugefluxflag=1, transform velocity to length; gaugefluxflag=2, transform length to velocity
    if ((gaugefluxflag.eq.1.and.velflag.ne.0).or.(gaugefluxflag.eq.2.and.velflag.eq.0)) then
       call gauge_transform(velflag,spftime*dt,nspf,inspfs(:,:), neutspfs(:,:))
    else
       neutspfs(:,:)=inspfs(:,:)
    endif

    if (strongcatflag.ne.0) then

!! The best treatment with strongcatflag explicitly couples cation 
!!   states when the field is on, propagating backwards from the end 
!!   of the pulse.

!! do cation backwards, timefac=(0d0,1d0), pulse must be exactly flipped

       call projeflux_strongcatload(ifile,spftime)
       catspfs(:,:)=ppp(ifile)%tmo(:,:)

    else

!! without strongcatflag...
!!   transform cation to velocity gauge if flux is calculated in velocity gauge.
!!   notice gaugefluxflag.ne.1 below
!! The cation states are bound states and so are located at small radius.
!!   Assume length gauge cation states are unperturbed.  So transform cation states
!!   to velocity gauge if flux is calculated in velocity gauge.

       if ((gaugefluxflag.ne.1.and.velflag.ne.0).or.(gaugefluxflag.eq.2.and.velflag.eq.0)) then
          call gauge_transform(0,spftime*dt,nspf,ppp(ifile)%tmo(:,:),catspfs(:,:))
       else
          catspfs(:,:)=ppp(ifile)%tmo(:,:)
       endif

    endif

!!$    REINSTATE VARIABLE WHICHSIDEPROJ (from version 0) HERE 
!!$
!!$  want to transform cation orbitals and a-vectors not neutral orbitals and
!!$  a-vectors so that cation orbitals and neutral orbitals are biorthogonal.
!!$  
!!$  Whichsideproj=0 (default) transform neutral   
!!$  Whichsideproj=1 (reinstating) transform cation
!!$
!!$  It is better to transform cation (whichsideproj=1) because then cation orbitals can
!!$  obey whatever orthogonality condition.  Whichsideproj=1 allows cation eigenfunctions
!!$  to be calculated with c-norm (e.g. cmctdhf_atom) for calculating photoionization 
!!$  using usual herm-norm propagation calculation (chmctdhf_atom).

!!$  Do extra cation walks?  Probably not best to use catww, mvalues etc might not line up.
!!$  Save for later.

!!$    if (whichsideproj.eq.0) then     !! WHICHSIDEPROJ.EQ.0 : transform neutral

    allocate(myavec(numr,first_config:last_config,mcscfnum))
    if (tot_adim.gt.0) then
       myavec(:,:,:)=inavectors(:,:,:)
    endif

    call bioset(projbiovar(ifile),smo,numr,bioww)

    if (tot_adim.gt.0) then
       call biortho(neutspfs(:,:),catspfs(:,:),mobio(:,:),myavec(:,:,1),projbiovar(ifile))
       do imc=2,mcscfnum
          call biotransform(neutspfs(:,:),catspfs(:,:),myavec(:,:,imc),projbiovar(ifile))
       enddo
    else
       call biortho(neutspfs(:,:),catspfs(:,:),mobio(:,:),nullvector(:),projbiovar(ifile))
       do imc=2,mcscfnum
          call biotransform(neutspfs(:,:),catspfs(:,:),nullvector(:),projbiovar(ifile))
       enddo
    endif

    call projeflux_allproj(ifile,ppp(ifile)%tavec(:,:,:),myavec,mobio,projwfn)

!!$    else    !! WHICHSIDEPROJ.NE.0: transform cation

!!$    endif  !! WHICHSIDEPROJ

    deallocate(mobio,catspfs,neutspfs,myavec)

  end subroutine projeflux_projectone


  subroutine projeflux_catdealloc()
    use parameters !! numcatfiles
    use projefluxmod
    implicit none
    integer :: ifile
    do ifile=1,numcatfiles
       deallocate(ppp(ifile)%tavec,ppp(ifile)%tmo)
       deallocate(ppp(ifile)%tconfiglist,ppp(ifile)%numpwalk1,ppp(ifile)%pwalk1,ppp(ifile)%pspf1,ppp(ifile)%pphase1)
       ppp(ifile)%tnumconfig=0
       ppp(ifile)%maxpwalk1=0 
       ppp(ifile)%eachstate=0
       ppp(ifile)%catfactor=0d0
    enddo
    catloaded=0
  end subroutine projeflux_catdealloc

!! do the double time integral piece    08-2015 now looping over istate and imc here

  subroutine projeflux_double_time_int(mem,totstate,nt,dt,NUMANGLES,NUMERAD)
    use parameters    !! par_timestep and others
    use mpimod
    use mpisubmod
    use pulsesubmod
    use projutilsubmod
    use utilmod
    implicit none
    integer,intent(in) :: mem,totstate,nt,NUMANGLES,NUMERAD
    real*8, intent(in) :: dt
    integer :: i,k,tlen,istate,curtime,tau,ir ,imc,ioffset, il
    integer :: BatchSize,NBat,ketreadsize,brareadsize,ketbat,brabat,kettime,bratime,&
         bratop,getlen,myiostat,oldtime
    real*8 :: MemTot,MemVal,wfi,estep,myfac,windowfunct,MemNum, tentsum
    DATATYPE, allocatable :: bramo(:,:,:,:),ketmo(:,:,:,:),gtau(:,:,:),gtau_ad(:,:,:,:),ketop(:,:,:,:),&
         read_bramo(:,:,:), read_ketmo(:,:,:), deweighted_bramo(:,:,:), ketmo_ad(:,:,:,:,:),&
         ketop_ad(:,:,:,:,:), gtaudiag(:,:,:), gtaudiag_ad(:,:,:,:), gtausum(:,:,:), gtausum_ad(:,:,:,:),&
         tot_gtausum(:,:), tot_gtausum_ad(:,:,:), gtaunow(:,:), gtaunow_ad(:,:,:), gtausave(:,:), &
         gtausave_ad(:,:,:)
    complex*16, allocatable :: ftgtau(:),pulseft(:,:), total(:),ftgtau_ad(:,:),total_ad(:,:),&
         ftgtausum_ad(:)
    complex*16 :: ftgtausum, csum
    real*8, allocatable :: pulseftsq(:), tentfunction(:)
    DATATYPE :: pots1(3)
    character (len=4) :: xstate0,xmc0
    character (len=3) :: xstate1,xmc1

    if (ceground.eq.(0d0,0d0)) then
       OFLWR "Eground is ZERO.  Are you sure?  If want zero just make it small."
       WRFL  "  Otherwise need eground: initial state energy."; CFLST
    endif

!! determine if we should do batching or not
!! 250,000 words/MB, real*8 2words/#, complex*16 4words/#

#ifdef REALGO
    MemVal = 1.25d5
#else
    MemVal = 6.25d4
#endif

!! rank 1 now only holds one entire orbital vector for parorbsplit.eq.3.
!! read_bramo and read_ketmo not dimensioned with BatchSize any more.
!! MemNum based on size of arrays bramo and ketmo instead.

!!$  if (parorbsplit.eq.3) then
!!$     MemNum=4*spfsize*numr*nprocs    !! 4 for bra and ket and both spun orbitals
!!$  else

    MemNum=4*spfsize*numr

!!$  endif

    call openfile()
    write(mpifileptr,'(A30,F9.3,A3)') " Guess at necessary memory is ",&
         (nt+1)*MemNum/MemVal," MB"
    if(mem.eq.0) then
       write(mpifileptr,*) "g(tau) will be computed with all of psi in core"
       BatchSize=nt+1
    else
       MemTot=real(mem,8)    
       write(mpifileptr,*) "g(tau) will be computed with all psi being read in batches"
       write(mpifileptr,'(A33,F9.3,A3)') " Desired amount of memory to use ",MemTot," MB"
       
       BatchSize=floor(MemTot * MemVal / MemNum) 
       if(BatchSize.lt.1) then
          write(mpifileptr,*) "Tiny amount of memory or huge wavefunction, Batchsize is 1" 
          BatchSize=1
       else if(BatchSize.ge.nt+1) then
          write(mpifileptr,*) "Hooray, there is enough memory, doing it all in core" 
          BatchSize=nt+1
       else
          write(mpifileptr,*) "Batchsize is ",BatchSize,"/",(nt+1)
       endif
    endif
    call closefile()

    allocate(gtau(0:nt,totstate,mcscfnum),ketmo(spfsize,numr,2,BatchSize),ketop(spfsize,numr,2,BatchSize),&
         bramo(spfsize,numr,2,BatchSize),gtaudiag(0:nt,totstate,mcscfnum),gtausum(0:nt,totstate,mcscfnum),&
         tot_gtausum(0:nt,mcscfnum),gtaunow(totstate,mcscfnum),gtausave(totstate,mcscfnum))
    gtau=0; ketmo=0; ketop=0; bramo=0; gtaudiag=0; gtausum=0; tot_gtausum=0; gtaunow=0; gtausave=0

    if (angularflag.ne.0) then
       allocate(gtau_ad(0:nt,totstate,mcscfnum,NUMANGLES),deweighted_bramo(spfsize,numr,2),&
            ketmo_ad(spfsize,numr,2,batchsize,NUMANGLES), ketop_ad(spfsize,numr,2,batchsize,NUMANGLES),&
            gtaudiag_ad(0:nt,totstate,mcscfnum,NUMANGLES),gtausum_ad(0:nt,totstate,mcscfnum,NUMANGLES),&
            tot_gtausum_ad(0:nt,mcscfnum,NUMANGLES), gtaunow_ad(totstate,mcscfnum,NUMANGLES), &
            gtausave_ad(totstate,mcscfnum,NUMANGLES))
    else
       allocate(gtau_ad(0:0,1,1,1),deweighted_bramo(1,1,2),&
            ketmo_ad(1,1,2,1,1), ketop_ad(1,1,2,1,1),&
            gtaudiag_ad(0:0,1,1,1),gtausum_ad(0:0,1,1,1),&
            tot_gtausum_ad(0:0,1,1),gtaunow_ad(1,1,1),gtausave_ad(1,1,1))
    endif

    gtau_ad=0;   deweighted_bramo=0;  ketmo_ad=0; ketop_ad=0; gtaudiag_ad=0; gtausum_ad=0; 
    tot_gtausum_ad=0; gtaunow_ad=0d0; gtausave_ad=0d0

    if (myrank.eq.1) then
       if (parorbsplit.eq.3) then
          allocate(read_ketmo(spfsize*nprocs,numr,2),&
               read_bramo(spfsize*nprocs,numr,2))
       else
          allocate(read_ketmo(spfsize,numr,2),&
               read_bramo(spfsize,numr,2))
       endif
    else
       allocate(read_ketmo(1,numr,2),read_bramo(1,numr,2))
    endif
    read_ketmo=0; read_bramo=0

    NBat=ceiling(real(nt+1)/real(BatchSize))
    ketreadsize=0;  brareadsize=0

    if (myrank.eq.1) then
       inquire (iolength=tlen) read_ketmo(:,1,:)
    endif
    call mympiibcastone(tlen,1)

    OFLWR "Projected 1e- function record length is ",tlen;  CFL

    gtau(:,:,:)=0d0
    gtau_ad(:,:,:,:)=0d0

!! looping here now 08-2015

    do imc=1,mcscfnum

       do istate=1,totstate

          OFLWR "Computing the CrossSection for state ",istate,"wfn",imc; CFL

          if (myrank.eq.1) then
             open(1003,file=projfluxfile,status="unknown",form="unformatted",&
                  access="direct",recl=tlen,iostat=myiostat)
             call checkiostat(myiostat,"opening "//projfluxfile)
          endif

!! lets do this double time integral double batched monster loop here. 
!! as can be seen in total flux eveything but biortho is fast as hell

          do ketbat=1,NBat
             OFLWR "Reading ket batch ", ketbat, " of ", NBat," for state ",istate; CFL
             ketreadsize=min(BatchSize,nt+1-(ketbat-1)*BatchSize)
     
!! read the orbital |\psi(t)>

             do i=1,ketreadsize !! loop over times in this batch
                if(myrank.eq.1) then
                   do ir=1,numr !! loop over all the r's for this time           

  ioffset=(istate-1)*mcscfnum*(nt+1)*numr + (imc-1)*(nt+1)*numr + ((ketbat-1)*BatchSize+i-1)*numr + ir

                      read(1003,rec=ioffset,iostat=myiostat) read_ketmo(:,ir,:)
                      call checkiostat(myiostat,"reading ketmo")
                   enddo
                endif
                if (parorbsplit.ne.3) then
                   if (myrank.eq.1) then
                      ketmo(:,:,:,i)=read_ketmo(:,:,:)
                   endif
                   call mympibcast(ketmo(:,:,:,i),1,spfsize*numr*2)
                else
                   do k=1,2
                      do ir=1,numr
                         call splitscatterv(read_ketmo(:,ir,k),ketmo(:,ir,k,i))
                      enddo
                   enddo
                endif
             enddo       !! do i=1,ketreadsize

             if (fluxoptype.eq.0) then
                OFLWR "PROGRAM ME FLUXOPTYPE 0 PROJEFLUX"; CFLST
             endif

             do i=1,ketreadsize
                do k=1,2 !! loop over the spins...
                   call projeflux_op_onee(ketmo(:,:,k,i), ketop(:,:,k,i))
                enddo
!! for integral dt
                curtime=(ketbat-1)*BatchSize+i-1
                gtaudiag(curtime,istate,imc) = hermdot(ketmo(:,:,:,i),ketop(:,:,:,i),2*spfsize*numr)

                if (curtime.eq.0) then
                   gtausave(istate,imc) = gtaudiag(curtime,istate,imc)
                endif

                if (flux_subtract.ne.0) then
                   gtaudiag(curtime,istate,imc) = gtaudiag(curtime,istate,imc) - gtausave(istate,imc)
                endif

             enddo

  if (angularflag.ne.0) then
     do il=1,NUMANGLES
        ketmo_ad(:,:,:,1:ketreadsize,il) = ketmo(:,:,:,1:ketreadsize)
        call only_one_angle(2*numr*ketreadsize,il,ketmo_ad(:,:,:,:,il),NUMANGLES)
        do i=1,ketreadsize
           do k=1,2 !! loop over the spins...
              call projeflux_op_onee(ketmo_ad(:,:,k,i,il), ketop_ad(:,:,k,i,il))
           enddo
!! for integral dt
           curtime=(ketbat-1)*BatchSize+i-1
           gtaudiag_ad(curtime,istate,imc,il) = hermdot(ketmo_ad(:,:,:,i,il),ketop_ad(:,:,:,i,il),2*spfsize*numr)

           if (curtime.eq.0) then
              gtausave_ad(istate,imc,il) = gtaudiag_ad(curtime,istate,imc,il)
           endif

           if (flux_subtract.ne.0) then
              gtaudiag_ad(curtime,istate,imc,il) = gtaudiag_ad(curtime,istate,imc,il) - gtausave_ad(istate,imc,il)
           endif

        enddo
     enddo
  endif  !! angularflag

!! begin the bra batch read loop
             do brabat=1,ketbat
                OFLWR "Reading bra batch ", brabat, " of ", ketbat," for state ",istate; CFL
                brareadsize=min(BatchSize,nt+1-(brabat-1)*BatchSize)
                do i=1,brareadsize
                   if(myrank.eq.1) then
                      do ir=1,numr 

   ioffset = (istate-1)*mcscfnum*(nt+1)*numr + (imc-1)*(nt+1)*numr + ((brabat-1)*BatchSize+i-1)*numr + ir

                         read(1003,rec=ioffset,iostat=myiostat) read_bramo(:,ir,:)
                         call checkiostat(myiostat,"reading bramo")
                      enddo
                   endif
                   if (parorbsplit.ne.3) then
                      if (myrank.eq.1) then
                         bramo(:,:,:,i)=read_bramo(:,:,:)
                      endif
                      call mympibcast(bramo(:,:,:,i),1,spfsize*numr*2)
                   else
                      do k=1,2
                         do ir=1,numr
                            call splitscatterv(read_bramo(:,ir,k),bramo(:,ir,k,i))
                         enddo
                      enddo
                   endif
                enddo     !! do i=1,brareadsize

!! loop over the specific ket indices & over all previous 
!! times & evaluate the actual g(tau) expression           

                do kettime=1,ketreadsize

                   curtime=(ketbat-1)*BatchSize+kettime-1 

                   if(brabat.lt.ketbat) then
                      bratop=brareadsize
                   else
                      bratop=kettime
                   endif
                   do bratime=1,bratop
                      oldtime=((brabat-1)*BatchSize+bratime-1)
                      tau=curtime-oldtime

                      gtaunow(istate,imc) = hermdot(bramo(:,:,:,bratime),ketop(:,:,:,kettime),2*spfsize*numr)

!! dt factor here
                      gtau(tau,istate,imc) = gtau(tau,istate,imc) + gtaunow(istate,imc) * dt

   if (angularflag.ne.0) then

      do k=1,2
         do ir=1,numr

!! with /2/pi, we have Mb per steradian column 4 angprojspifile... right?  column 5 is nonsense pmctdhf/chmctdhf
!! weights elecweights(:,:,:,3)=1d0.

            deweighted_bramo(:,ir,k) = bramo(:,ir,k,bratime) / &
                 RESHAPE(elecweights(:,:,:,2),(/spfsize/)) / 2 / pi
         enddo
      enddo

!! symmetric flux operator.  total still equals integrated differential.
!! results (at threshold where there is oscillation and usually some negative values) actually 
!! look worse for H2 with this symmetrized operator, below, instead of what's above

      gtaunow_ad(istate,imc,:) = 0.5d0 * &
           myhermdots(deweighted_bramo(:,:,:),ketop(:,:,:,kettime),2*spfsize*numr,NUMERAD,NUMANGLES)

      do il=1,NUMANGLES
         gtaunow_ad(istate,imc,il) = gtaunow_ad(istate,imc,il) + 0.5d0 * &
              hermdot(deweighted_bramo(:,:,:),ketop_ad(:,:,:,kettime,il),2*spfsize*numr)
      enddo

!! dt factor here
      gtau_ad(tau,istate,imc,:) = gtau_ad(tau,istate,imc,:) + gtaunow_ad(istate,imc,:) * dt

   endif  !! angularflag

                   enddo !! do bratime
                enddo !! do kettime
             enddo   !! do brabat
          enddo   !! do ketbat
        
          if (curtime.ne.nt) then
             OFLWR "DEEG CURTIME NT ERR", curtime, nt; CFLST
          endif

          if (myrank.eq.1) then
             close(1003)
          endif

       enddo   !! do istate
    enddo      !! do imc

    if (parorbsplit.eq.3) then
       call mympireduce(gtau(:,:,:), (nt+1)*totstate*mcscfnum)
       call mympireduce(gtaudiag(:,:,:), (nt+1)*totstate*mcscfnum) 
       if (angularflag.ne.0) then
          call mympireduce(gtau_ad(:,:,:,:), (NUMANGLES)*(nt+1)*totstate*mcscfnum)
          call mympireduce(gtaudiag_ad(:,:,:,:), (NUMANGLES)*(nt+1)*totstate*mcscfnum)
       endif
    endif

!! INTEGRAL DT

    gtausum(:,:,:)=0d0;  gtausum_ad(:,:,:,:)=0d0
    gtausum(0,:,:)=gtaudiag(0,:,:)
    do i=1,nt
       gtausum(i,:,:)=gtausum(i-1,:,:) + gtaudiag(i,:,:) * dt / 2d0   !! 2 looks correct
    enddo
    tot_gtausum=0d0
    do istate=1,totstate
       tot_gtausum(:,:)=tot_gtausum(:,:) + gtausum(:,istate,:)
    enddo

    if (angularflag.ne.0) then
       gtausum_ad(0,:,:,:)=gtaudiag_ad(0,:,:,:)
       do i=1,nt
          gtausum_ad(i,:,:,:)=gtausum_ad(i-1,:,:,:) + gtaudiag_ad(i,:,:,:) * dt / 2d0
       enddo
       tot_gtausum_ad=0d0
       do istate=1,totstate
          tot_gtausum_ad(:,:,:)=tot_gtausum_ad(:,:,:) + gtausum_ad(:,istate,:,:)
       enddo
    endif

    if (myrank.eq.1) then

!! nevermind KVLsum.dat and gtau.dat
!! now integrals dt and domega with projfluxtsumfile, angprojfluxtsumfile in namelist &parinp
!! and additional column in projspifile

!! integral dt:

       do imc=1,mcscfnum
          write(xmc0,'(I4)') imc+1000;   xmc1=xmc0(2:4)

          open(171,file=projfluxtsumfile(1:getlen(projfluxtsumfile))//"_all_"//xmc1//".dat", &
               status="unknown",iostat=myiostat)
          call checkiostat(myiostat,"opening projfluxtsumfile tot")
          write(171,*,iostat=myiostat) "#   ", curtime
          call checkiostat(myiostat,"writing projfluxtsumfile tot")
          do i=0,curtime
             write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat) &
                  i*dt,  tot_gtausum(i,imc)
          enddo
          call checkiostat(myiostat,"writing projfluxtsumfile tot")
          close(171)

          if (angularflag.ne.0) then
             open(171,file=angprojfluxtsumfile(1:getlen(angprojfluxtsumfile))//"_all_"//xmc1//".dat", &
                  status="unknown",iostat=myiostat)
             call checkiostat(myiostat,"opening angprojfluxtsumfile tot")
             write(171,*,iostat=myiostat) "#   ", curtime
             call checkiostat(myiostat,"writing angprojfluxtsumfile tot")
             do i=0,curtime
                write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat) &
                     i*dt,  tot_gtausum_ad(i,imc,:)
             enddo
             call checkiostat(myiostat,"writing angprojfluxtsumfile tot")
             close(171)
          endif

          do istate=1,totstate
             write(xstate0,'(I4)') istate+1000;        xstate1=xstate0(2:4)

             open(171,file=projfluxtsumfile(1:getlen(projfluxtsumfile))//"_"//xstate1//"_"//xmc1//".dat", &
                  status="unknown",iostat=myiostat)
             call checkiostat(myiostat,"opening projfluxtsumfile")
             write(171,*,iostat=myiostat) "#   ", curtime
             call checkiostat(myiostat,"writing projfluxtsumfile")
             do i=0,curtime
                write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat) &
                     i*dt,  gtausum(i,istate,imc), gtaudiag(i,istate,imc)
             enddo
             call checkiostat(myiostat,"writing projfluxtsumfile")
             close(171)

             if (angularflag.ne.0) then
                open(171,file=angprojfluxtsumfile(1:getlen(angprojfluxtsumfile))//"_"//xstate1//"_"//xmc1//".dat", &
                     status="unknown",iostat=myiostat)
                call checkiostat(myiostat,"opening angprojfluxtsumfile")
                write(171,*,iostat=myiostat) "#   ", curtime
                call checkiostat(myiostat,"writing angprojfluxtsumfile")
                do i=0,curtime
                   write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat) &
                        i*dt,  gtausum_ad(i,istate,imc,:), gtaudiag_ad(i,istate,imc,:)
                enddo
                call checkiostat(myiostat,"writing angprojfluxtsumfile")
                close(171)
             endif

          enddo !! istate
       enddo    !! imc
    endif       !! myrank

    call mpibarrier()

    OFLWR "Taking the fourier transform of g(tau) to get cross section at T= ",curtime*dt; CFL

    allocate(ftgtau(-curtime:curtime), pulseft(-curtime:curtime,3),&
         pulseftsq(-curtime:curtime),total(-curtime:curtime),tentfunction(-curtime:curtime))
    ftgtau(:)=0d0; pulseft(:,:)=0d0; pulseftsq(:)=0d0; total(:)=0; tentfunction(:)=0

    do i=-curtime,curtime
       tentfunction(i) = (1+curtime-abs(i)) * windowfunct(abs(i),curtime,17)
    enddo

    tentsum = get_rtentsum(curtime,tentfunction(-curtime:curtime))

    if (angularflag.ne.0) then
       allocate(ftgtau_ad(-curtime:curtime,NUMANGLES), total_ad(-curtime:curtime,NUMANGLES),&
            ftgtausum_ad(NUMANGLES))
    else
       allocate(ftgtau_ad(0:0,1), total_ad(0:0,1), ftgtausum_ad(1))
    endif
    ftgtau_ad=0; total_ad=0; ftgtausum_ad=0

    do i=0,curtime
       call vectdpot(i*dt,0,pots1,-1)  !! LENGTH GAUGE
       if (pulsewindowtoo == 0) then
          pulseft(i,:)=pots1(:)
       else
          pulseft(i,:)=pots1(:) * windowfunct(i,curtime,17)  !! action 17
       endif
    enddo

    do i=1,3
       call zfftf_wrap(2*curtime+1,pulseft(-curtime:curtime,i))
    enddo

    pulseft(:,:)=pulseft(:,:)*dt

    do i=-curtime,curtime
       pulseft(i,:)=pulseft(i,:)*exp((0.d0,1.d0)*(curtime+i)*curtime*2*pi/real(2*curtime+1))
    enddo

    pulseftsq(:) = abs(pulseft(:,1)**2) + abs(pulseft(:,2)**2) + abs(pulseft(:,3)**2)

    estep=2*pi/dt/(2*curtime+1)
  
    do imc=1,mcscfnum

       write(xmc0,'(I4)') imc+1000;   xmc1=xmc0(2:4)

       total(:)=0d0
       total_ad(:,:)=0d0

       do istate=1,totstate

          write(xstate0,'(I4)') istate+1000;    xstate1=xstate0(2:4)

          ftgtau(:)=0d0;
          do i=0,curtime
             ftgtau(i) = ALLCON(gtau(i,istate,imc)) * windowfunct(i,curtime,17) &
                  * exp((0d0,-1d0)*ALLCON(ceground)*dt*i)
          enddo
          do i=1,curtime
             ftgtau(-i) = ALLCON(ftgtau(i))
          enddo

!! subtract tent function   (1+curtime-abs(i))/(curtime+1)^2   for better performance

          if (flux_subtract.ne.0) then
             csum = get_ctentsum(curtime, ftgtau(-curtime:curtime))
             ftgtau(-curtime:curtime) = ftgtau(-curtime:curtime) - &
                  csum/tentsum * tentfunction(-curtime:curtime)
          endif

          if (angularflag.ne.0) then
             do i=0,curtime
                ftgtau_ad(i,:) = ALLCON(gtau_ad(i,istate,imc,:))   * windowfunct(i,curtime,17) * &
                     exp((0d0,-1d0)*ALLCON(ceground)*dt*i)
             enddo
             do i=1,curtime
                ftgtau_ad(-i,:) = ALLCON(ftgtau_ad(i,:))
             enddo
             if (flux_subtract.ne.0) then
                do il=1,NUMANGLES
                   csum = get_ctentsum(curtime, ftgtau_ad(-curtime:curtime,il))
                   ftgtau_ad(-curtime:curtime,il) = ftgtau_ad(-curtime:curtime,il) - &
                        csum/tentsum * tentfunction(-curtime:curtime)
                enddo
             endif
          endif

!! gtaufile

          if(myrank.eq.1) then

             open(777,file="Dat/gtau_"//xstate1//"_"//xmc1//".dat",&
                  status="replace",action="readwrite",position="rewind",iostat=myiostat)
             call checkiostat(myiostat,"opening gtau file")
             write(777,*,iostat=myiostat)
             call checkiostat(myiostat,"writing gtau file")
             do i=0,curtime
                write(777,'(F18.12, T22, 400E20.8)',iostat=myiostat)  i*dt, ftgtau(i), &
                     ALLCON(gtau(i,istate,imc)) * exp((0d0,-1d0)*ALLCON(ceground)*dt*i)
             enddo
             call checkiostat(myiostat,"writing gtau file")
             close(777)

          endif

          call zfftf_wrap_diff(2*curtime+1,ftgtau(-curtime:curtime),ftdiff)
        
          ftgtau(:)=ftgtau(:)*dt

          do i=-curtime,curtime
             ftgtau(i)=ftgtau(i)*exp((0.d0,1.d0)*(curtime+i)*curtime*2*pi/real(2*curtime+1))
          enddo

          total(:)=total(:)+ftgtau(:)

          if (angularflag.ne.0) then
             do il=1,NUMANGLES
                call zfftf_wrap_diff(2*curtime+1,ftgtau_ad(-curtime:curtime,il),ftdiff)
             enddo

             ftgtau_ad(:,:)=ftgtau_ad(:,:)*dt

             do i=-curtime,curtime
                ftgtau_ad(i,:)=ftgtau_ad(i,:)*exp((0.d0,1.d0)*(curtime+i)*curtime*2*pi/real(2*curtime+1))
             enddo

             total_ad(:,:)=total_ad(:,:)+ftgtau_ad(:,:)
          endif

          if(myrank.eq.1) then

             open(1004,file=projspifile(1:getlen(projspifile))//"_"//xstate1//"_"//xmc1//".dat",&
                  status="replace",action="readwrite",position="rewind",iostat=myiostat)
             call checkiostat(myiostat,"opening proj spi file")
             write(1004,*,iostat=myiostat)
             call checkiostat(myiostat,"writing proj spi file")
             write(1004,*) "# Omega; pulse ft; proj flux t= ",finaltime

!! nevermind KVLsum.dat and gtau.dat
!! now integrals dt and domega with projfluxtsumfile, angprojfluxtsumfile in namelist &parinp
!! and additional column in projspifile

!! integral domega in ftgtausum

             ftgtausum=0d0
             do i=-curtime,curtime
                wfi=(i+curtime)*estep

                ftgtausum=ftgtausum + ftgtau(i) * estep / PI / 4d0   !! 4 looks correct

!! LENGTH GAUGE WAS FT'ed multiply by wfi don't divide
!! NEVERMIND FACTOR OF 1/3
!!              myfac = 5.291772108d0**2 / 3d0 * 2d0 * PI / 1.37036d2 * wfi 

!! WITH THIS FACTOR, NOW THE QUANTUM MECHANICAL PHOTOIONIZATION CROSS SECTION IN 
!! MEGABARNS (10^-18 cm^2) IS IN COLUMN 3 REAL PART
                myfac = 5.291772108d0**2 * 2d0 * PI / 1.37036d2 * wfi 

                write(1004,'(F18.12, T22, 400E20.8)',iostat=myiostat)  wfi,  &
                     pulseftsq(i), ftgtau(i)/pulseftsq(i) * myfac, ftgtau(i) / 4 / PI, ftgtausum

             enddo
             call checkiostat(myiostat,"writing proj spi file")
             close(1004)

             if (angularflag.ne.0) then
                open(1004,file=angprojspifile(1:getlen(angprojspifile))//"_"//xstate1//"_"//xmc1//".dat",&
                     status="replace",action="readwrite",position="rewind",iostat=myiostat)
                call checkiostat(myiostat,"opening angproj spi file")
                write(1004,*,iostat=myiostat)
                call checkiostat(myiostat,"writing angproj spi file")
                write(1004,*) "# Omega; pulse ft; differential proj flux t= ",finaltime

                ftgtausum_ad(:)=0d0
                do i=-curtime,curtime
                   wfi=(i+curtime)*estep

                   ftgtausum_ad(:)=ftgtausum_ad(:) + ftgtau_ad(i,:) * estep / PI / 4d0   !! 4 looks correct

                   myfac = 5.291772108d0**2 * 2d0 * PI / 1.37036d2 * wfi 

                   do il=1,NUMANGLES
                      write(1004,'(F18.12, I5, 1400E20.8)',iostat=myiostat)  wfi, il, &
                           pulseftsq(i), ftgtau_ad(i,il)/pulseftsq(i) * myfac, &
                           ftgtau_ad(i,il) / 4 / PI, ftgtausum_ad(il)
                   enddo
                   write(1004,*)
                enddo
                call checkiostat(myiostat,"writing angproj spi file")
                close(1004)
             endif
          endif  !! if myrank.eq.1
       enddo  !! do istate

       if(myrank.eq.1) then
          open(1004,file=projspifile(1:getlen(projspifile))//"_all_"//xmc1//".dat",&
               status="replace",action="readwrite",position="rewind",iostat=myiostat)
          call checkiostat(myiostat,"opening total proj spi file")
          write(1004,*,iostat=myiostat)
          call checkiostat(myiostat,"writing total proj spi file")
          write(1004,*) "# Omega; pulse ft; proj flux t= ",finaltime

          ftgtausum=0d0
          do i=-curtime,curtime
             wfi=(i+curtime)*estep

             ftgtausum=ftgtausum + total(i) * estep / PI / 4d0   !! 4 looks correct

             myfac = 5.291772108d0**2 * 2d0 * PI / 1.37036d2 * wfi
             write(1004,'(F18.12, T22, 400E20.8)',iostat=myiostat)  wfi,  pulseftsq(i), &
                  total(i)/pulseftsq(i) * myfac, total(i) / 4 / PI, ftgtausum
          enddo
          call checkiostat(myiostat,"writing total proj spi file")
          close(1004)

          if (angularflag.ne.0) then
             open(1004,file=angprojspifile(1:getlen(angprojspifile))//"_all_"//xmc1//".dat",&
                  status="replace",action="readwrite",position="rewind",iostat=myiostat)
             call checkiostat(myiostat,"opening total ang proj spi file")
             write(1004,*,iostat=myiostat)
             call checkiostat(myiostat,"writing total ang proj spi file")
             write(1004,*) "# Omega; pulse ft; differential proj flux t= ",finaltime

             ftgtausum_ad(:)=0d0
             do i=-curtime,curtime
                wfi=(i+curtime)*estep

                ftgtausum_ad(:)=ftgtausum_ad(:) + total_ad(i,:) * estep / PI / 4d0   !! 4 looks correct

                myfac = 5.291772108d0**2 * 2d0 * PI / 1.37036d2 * wfi

                do il=1,NUMANGLES
                   write(1004,'(F18.12, I5, 400E20.8)',iostat=myiostat)  wfi, il, pulseftsq(i), &
                        total_ad(i,il)/pulseftsq(i) * myfac, &
                        total_ad(i,il) / 4 / PI, ftgtausum_ad(il)
                enddo
                write(1004,*)
             enddo
             call checkiostat(myiostat,"writing total ang proj spi file")
             close(1004)
          endif
       endif

       call mpibarrier()

    enddo  !! do imc

    deallocate(gtau, ketmo, ketop,&
         bramo, gtaudiag, gtausum,&
         tot_gtausum, gtaunow, gtausave)
    deallocate(gtau_ad, deweighted_bramo,&
         ketmo_ad, ketop_ad,&
         gtaudiag_ad, gtausum_ad,&
         tot_gtausum_ad, gtaunow_ad, gtausave_ad)
    deallocate(read_ketmo, read_bramo)
    deallocate(ftgtau, pulseft, pulseftsq, total, tentfunction)
    deallocate(ftgtau_ad, total_ad, ftgtausum_ad)

  end subroutine projeflux_double_time_int

  subroutine projeflux_statedealloc()
    use projefluxmod
    nstate=0
    deallocate(alreadystate,ppp)
  end subroutine projeflux_statedealloc

end module projcoresubmod



module projfluxduringmod
  implicit none
  integer :: allocated=0, curtime=0
  DATATYPE,allocatable :: gtausum(:,:), gtausum_ad(:,:,:), gtausave(:,:), gtausave_ad(:,:,:)
end module projfluxduringmod


module projactionmod
contains

subroutine projeflux_single(mem)
  use parameters    !! fluxinterval, par_timestep, etc.
  use projefluxmod
  use projcoresubmod
  use mpimod   !! myrank
  implicit none
  integer,intent(in) :: mem
  integer :: nt,myiostat
  real*8 :: dt

  OFLWR ;  WRFL   "   *** DOING PROJECTED FLUX. ***    ";  WRFL; CFL

  call projeflux_load()

  nt=floor(real(numpropsteps,8)/fluxinterval/fluxskipmult)
  dt=real(FluxInterval*FluxSkipMult,8)*par_timestep

  call projeflux_projectdisk(nt, dt)

  call projeflux_catdealloc()

  OFLWR "Go double time integral"; CFL


!! do the double time integral
  call projeflux_double_time_int(mem,nstate,nt,dt,NUMANGLES,NUMERAD)

  call mpibarrier()
  OFLWR "Cross Section acquired, cleaning and closing";CFL

!! not necessary
  call projeflux_statedealloc()

  if (myrank.eq.1) then
     open(1003,file=projfluxfile,status="unknown",iostat=myiostat)
     write(1003,*) "Projfluxfile deleted"
     close(1003)
  endif
  call mpibarrier()
  call mpistop()

end subroutine projeflux_single


subroutine projeflux_during(inspfs,inavectors,dt)
  use parameters    !! mcscfnum
  use mpimod   !! myrank
  use projefluxmod
  use projfluxduringmod
  use projutilsubmod
  use projcoresubmod
  use mpisubmod
  implicit none
  real*8,intent(in) :: dt
  DATATYPE,intent(in) :: inspfs(spfsize,nspf), inavectors(tot_adim,mcscfnum)
  DATATYPE,allocatable :: gtaunow(:,:), projwfn(:,:,:,:,:), gtaunow_ad(:,:,:), &
       ketmo(:,:), ketop(:,:),tot_gtausum(:), tot_gtausum_ad(:,:)
  integer :: ifile, k, imc, istate, il, myiostat, getlen
  character (len=4) :: xstate0,xmc0
  character (len=3) :: xstate1,xmc1

  if (allocated.eq.0) then

     call projeflux_load()

     curtime=0
     allocated=1
     allocate(gtausum(nstate,mcscfnum),gtausave(nstate,mcscfnum))
     gtausum=0; gtausave=0
     if (angularflag.ne.0) then
        allocate(gtausum_ad(nstate,mcscfnum,NUMANGLES),gtausave_ad(nstate,mcscfnum,NUMANGLES))
        gtausum_ad=0; gtausave_ad=0
     endif

     if (myrank.eq.1) then
        do imc=1,mcscfnum
           write(xmc0,'(I4)') imc+1000;    xmc1=xmc0(2:4)

           open(171,file=projfluxtsumfile(1:getlen(projfluxtsumfile))//"_all_"//xmc1//".dat", &
                status="unknown",iostat=myiostat)
           call checkiostat(myiostat,"opening projfluxtsumfile tot")
           write(171,*,iostat=myiostat) "#   ", curtime
           call checkiostat(myiostat,"writing projfluxtsumfile tot")
           close(171)

           if (angularflag.ne.0) then
              open(171,file=angprojfluxtsumfile(1:getlen(angprojfluxtsumfile))//"_all_"//xmc1//".dat", &
                   status="unknown",iostat=myiostat)
              call checkiostat(myiostat,"opening angprojfluxtsumfile tot")
              write(171,*,iostat=myiostat) "#   ", curtime
              call checkiostat(myiostat,"writing angprojfluxtsumfile tot")
              close(171)
           endif

           do istate=1,nstate
              write(xstate0,'(I4)') istate+1000;    xstate1=xstate0(2:4)

              open(171,file=projfluxtsumfile(1:getlen(projfluxtsumfile))//"_"//xstate1//"_"//xmc1//".dat", &
                   status="unknown",iostat=myiostat)
              call checkiostat(myiostat,"opening projfluxtsumfile")
              write(171,*,iostat=myiostat) "#   ", curtime
              call checkiostat(myiostat,"writing projfluxtsumfile")
              close(171)

              if (angularflag.ne.0) then
                 open(171,file=angprojfluxtsumfile(1:getlen(angprojfluxtsumfile))//"_"//xstate1//"_"//xmc1//".dat", &
                      status="unknown",iostat=myiostat)
                 call checkiostat(myiostat,"opening angprojfluxtsumfile")
                 write(171,*,iostat=myiostat) "#   ", curtime
                 call checkiostat(myiostat,"writing angprojfluxtsumfile")
                 close(171)
              endif
           enddo
        enddo
     endif
  else
     curtime=curtime+1
  endif

  allocate(gtaunow(nstate,mcscfnum),tot_gtausum(mcscfnum))
  gtaunow=0; tot_gtausum=0
  if (angularflag.ne.0) then
     allocate(gtaunow_ad(nstate,mcscfnum,NUMANGLES),tot_gtausum_ad(mcscfnum,NUMANGLES))
  else
     allocate(gtaunow_ad(1,1,1),tot_gtausum_ad(1,1))
  endif
  gtaunow_ad=0; tot_gtausum_ad=0

  allocate(ketmo(spfsize,numr),ketop(spfsize,numr))

  do ifile=1,numcatfiles
     allocate(projwfn(spfsize,2,numr,ppp(ifile)%eachstate,mcscfnum))
     call projeflux_projectone(ifile, inspfs, curtime, dt, inavectors, projwfn)

     if (fluxoptype.eq.0) then
        OFLWR "PROGRAM ME FLUXOPTYPE 0 PROJEFLUX"; CFLST
     endif

     do imc=1,mcscfnum
        do istate=1,ppp(ifile)%eachstate
           do k=1,2 !! loop over the spins...
              ketmo(:,:)=projwfn(:,k,:,istate,imc)

              call projeflux_op_onee(ketmo,ketop)

              gtaunow(istate+alreadystate(ifile),imc) = gtaunow(istate+alreadystate(ifile),imc) + hermdot(ketmo,ketop,spfsize*numr)

           enddo
        enddo
     enddo

     if (angularflag.ne.0) then
        do il=1,NUMANGLES
           do imc=1,mcscfnum
              do istate=1,ppp(ifile)%eachstate
                 do k=1,2

                    ketmo(:,:)=projwfn(:,k,:,istate,imc)

                    call only_one_angle(numr,il,ketmo(:,:),NUMANGLES)
                    call projeflux_op_onee(ketmo, ketop)

                    gtaunow_ad(istate+alreadystate(ifile),imc,il) = gtaunow_ad(istate+alreadystate(ifile),imc,il) + hermdot(ketmo,ketop,spfsize*numr)
                 enddo
              enddo
           enddo
        enddo
     endif  !! angularflag

     deallocate(projwfn)

  enddo

  if (parorbsplit.eq.3) then
     call mympireduce(gtaunow(:,:),nstate*mcscfnum)
     if (angularflag.ne.0) then
        call mympireduce(gtaunow_ad(:,:,:),nstate*mcscfnum*NUMANGLES)
     endif
  endif

!! 06-16 subtract matrix element at t=0 (nonzero if there is leakage of
!!  eigenfunction onto ecs grid, eigenvalue with imaginary part)

  if (curtime.eq.0) then
     gtausave(:,:) = gtaunow(:,:)
     if (angularflag.ne.0) then
        gtausave_ad(:,:,:) = gtaunow_ad(:,:,:)
     endif
  endif

  if (flux_subtract.ne.0) then
     gtaunow(:,:)=gtaunow(:,:)-gtausave(:,:)
     if (angularflag.ne.0) then
        gtaunow_ad(:,:,:)=gtaunow_ad(:,:,:)-gtausave_ad(:,:,:)
     endif
  endif

  gtausum(:,:)=gtausum(:,:) + gtaunow(:,:) * dt / 2d0   !! 2 looks correct
  tot_gtausum=0d0
  do istate=1,nstate
     tot_gtausum(:)=tot_gtausum(:)+gtausum(istate,:)
  enddo

  if (angularflag.ne.0) then
     gtausum_ad(:,:,:)=gtausum_ad(:,:,:) + gtaunow_ad(:,:,:) * dt / 2d0   !! 2 looks correct
     tot_gtausum_ad=0d0
     do istate=1,nstate
        tot_gtausum_ad(:,:)=tot_gtausum_ad(:,:)+gtausum_ad(istate,:,:)
     enddo
  endif

  if (myrank.eq.1) then
     do imc=1,mcscfnum
        write(xmc0,'(I4)') imc+1000;   xmc1=xmc0(2:4)

        open(171,file=projfluxtsumfile(1:getlen(projfluxtsumfile))//"_all_"//xmc1//".dat", &
             status="old",position="append",iostat=myiostat)
        call checkiostat(myiostat,"opening projfluxtsumfile tot")
        write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat) &
             curtime*dt,  tot_gtausum(imc)
        call checkiostat(myiostat,"writing projfluxtsumfile tot")
        close(171)

        if (angularflag.ne.0) then
           open(171,file=angprojfluxtsumfile(1:getlen(angprojfluxtsumfile))//"_all_"//xmc1//".dat", &
                status="old",position="append",iostat=myiostat)
           call checkiostat(myiostat,"opening projfluxtsumfile tot")           
           write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat) &
                curtime*dt,  tot_gtausum_ad(imc,:)
           call checkiostat(myiostat,"writing angprojfluxtsumfile tot")
           close(171)
        endif

        do istate=1,nstate
           write(xstate0,'(I4)') istate+1000;      xstate1=xstate0(2:4)

           open(171,file=projfluxtsumfile(1:getlen(projfluxtsumfile))//"_"//xstate1//"_"//xmc1//".dat", &
                status="old",position="append",iostat=myiostat)
           call checkiostat(myiostat,"opening projfluxtsumfile")
           write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat) &
                curtime*dt,  gtausum(istate,imc), gtaunow(istate,imc)
           call checkiostat(myiostat,"writing projfluxtsumfile")
           close(171)

           if (angularflag.ne.0) then
              open(171,file=angprojfluxtsumfile(1:getlen(angprojfluxtsumfile))//"_"//xstate1//"_"//xmc1//".dat", &
                   status="old",position="append",iostat=myiostat)
              call checkiostat(myiostat,"opening projfluxtsumfile")           
              write(171,'(F18.12, T22, 400E20.8)',iostat=myiostat) &
                   curtime*dt,  gtausum_ad(istate,imc,:), gtaunow_ad(istate,imc,:)
              call checkiostat(myiostat,"writing angprojfluxtsumfile")
              close(171)
           endif
        enddo
     enddo
  endif

  call mpibarrier()

!QQQ actions_final?  call projeflux_catdealloc()
!QQQ                 call projeflux_statedealloc()

  deallocate(gtaunow,gtaunow_ad,ketmo,ketop,tot_gtausum,tot_gtausum_ad)
  
end subroutine projeflux_during

end module projactionmod



