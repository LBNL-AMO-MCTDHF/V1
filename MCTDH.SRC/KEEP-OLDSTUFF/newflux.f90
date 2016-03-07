
!! DISSOCIATIE FLUX ROUTINES

#include "Definitions.INC"

module fluxmod
  implicit none
  DATATYPE, allocatable :: dots(:,:),  curves(:,:), hdots(:,:)
  complex*16, allocatable :: imagkemod(:,:)
  integer :: maxnumrefconfig=0
  integer :: maxnumrefspf=0
  integer, allocatable ::  calledflag(:)
  integer, allocatable :: numrefconfig(:)   !! numfluxfiles
  integer, allocatable :: numrefspf(:)   !! numfluxfiles

 !! maxfluxcurves is max # of vectors from single input (mcscf.bin) file
  integer :: maxfluxcurves=0
  integer :: totfluxcurves=0
  integer, allocatable :: refconfiglists(:,:,:)  !! ndof, maxnumrefconfig, numfluxfiles
  DATATYPE, allocatable :: refvects(:,:,:,:) !! maxnumrefconfig, maxfluxcurves, numecspoints, numfluxfiles
  DATAECS, allocatable :: refvals(:,:,:)   !! maxfluxcurves, numecspoints, numfluxfiles
  DATATYPE, allocatable :: refspfs(:,:,:,:) !! spfsize, maxnumrefspf, numecspoints, numfluxfiles

end module fluxmod

subroutine fluxinit()
  use parameters
  use fluxmod
  use h2projectmod
  use opmod
  use commonmod
  implicit none
  integer :: i, myiostat, dummyint(100)
  integer :: mynuccharge1, mynuccharge2, myndof, mymrestrictflag, myrestrictflag, mynumr, mymrestrictval, myrestrictms, myspinrestrictval, myspfrestrictflag 

  if (bondgridpoints.ne.numr+2) then
     call openfile();     write(mpifileptr,*) "Bondgridpoints err ", bondgridpoints, numr;     call closefile();     call mpistop()
  endif
  if (bondcelement.gt.bondnumelements) then
     call openfile();     write(mpifileptr,*) "No flux, ECS:", bondcelement, bondnumelements;     call closefile();     call mpistop()
  endif
  allocate(imagkemod(numecspoints,numecspoints))
  allocate(numrefconfig(numfluxfiles), numrefspf(numfluxfiles))
  allocate(calledflag(numprop))
  calledflag=0

  imagkemod(:,:) = imag(rkemod(firstecspoint:numr,firstecspoint:numr) + (0.d0,0.d0))

  maxfluxcurves=0;  totfluxcurves=0;  maxnumrefconfig=0;  maxnumrefspf=0

  do i=1,numfluxfiles
     if (numfluxcurves(i).gt.maxfluxcurves) then
        maxfluxcurves=numfluxcurves(i)
     endif
     totfluxcurves=totfluxcurves + numfluxcurves(i)

     open(9954,file=fluxfilenames(i), form="unformatted", status="old", iostat=myiostat)
     if (myiostat/=0) then
        call openfile();        write(mpifileptr,*) "Error on read fluxfile = ", fluxfilenames(i), "iostat= ", myiostat;        call closefile();        call mpistop()
     endif

!! dummy reads for now.
     read(9954) mynuccharge1, mynuccharge2, myndof
!! 052412: myallspinproject -> numr (need numr not allspinproject)
     read(9954) mymrestrictflag, myspfrestrictflag, myrestrictflag, mynumr
     read(9954) mymrestrictval, myrestrictms, myspinrestrictval
     read(9954) numrefconfig(i), numrefspf(i)
     read(9954) dummyint(1:nspf)
     close(9954)
     if (numrefconfig(i).gt.maxnumrefconfig) then
        maxnumrefconfig=numrefconfig(i)
     endif
     if (numrefspf(i).gt.maxnumrefspf) then
        maxnumrefspf=numrefspf(i)
     endif
  enddo

  call openfile()
  write(mpifileptr,*) " Read info from top of fluxfiles. "
  write(mpifileptr,*) " Totfluxcurves is ", totfluxcurves
  write(mpifileptr,*) " Maxnumrefconfig, maxnumrefspf: ", maxnumrefconfig, maxnumrefspf
  call closefile()

  allocate(  refconfiglists(ndof,maxnumrefconfig,numfluxfiles) )
  allocate(  refvects(maxnumrefconfig, maxfluxcurves, numecspoints, numfluxfiles) )
  allocate(  refvals(maxfluxcurves,numecspoints,numfluxfiles) )
  allocate(  refspfs(spfsize, maxnumrefspf, numecspoints, numfluxfiles) )
  allocate(dots(numecspoints,totfluxcurves), curves(numecspoints,totfluxcurves), hdots(numecspoints,totfluxcurves))

  call fluxread()

end subroutine fluxinit


subroutine flux( avectorin, spfin,  iprop, thistime )
  use configmod
  use fluxmod
  use commonmod
  use parameters
  implicit none

  integer ::  jj, ifile, icurve,  isplit,  iprop
  integer ::  numcalled=0
  real*8 :: thistime
  integer, parameter :: printflag=0
  DATATYPE :: avectorin(numconfig,numr), spfin(spfsize,nspf)
  DATATYPE :: permovl(1, maxfluxcurves) !! AUTOMATIC

  if (iprop.gt.100) then
     call openfile();     write(mpifileptr,*) "REDIM LC!!!";     call closefile();     call mpistop()
  endif

  numcalled=numcalled+1
  icurve=0

  do ifile=1,numfluxfiles
     if (printflag.eq.1) then
        call openfile();        write(mpifileptr,*)  "go flux file ", ifile;        call closefile()
     endif
     do isplit=firstecspoint,numr
        do jj=1,numfluxcurves(ifile)
           call permoverlaps(1,numelec, spfsize, spfin(:,:),refspfs(:,:,isplit-firstecspoint+1,ifile), avectorin(:,isplit), &
                refvects(:,jj,isplit-firstecspoint+1,ifile), permovl(1,jj), 0, fluxpermthresh, fluxnormthresh, nspf, numrefspf(ifile), numconfig, numrefconfig(ifile), configlist, ndof, refconfiglists(:,:,ifile),ndof, 0)
        enddo

        !! permovl(1,jj) contains the overlap of current wfn with reference vector jj  < current | reference(jj) >
        do jj=1,numfluxcurves(ifile)
           dots(isplit-firstecspoint+1, jj+icurve)=CONJUGATE(permovl(1,jj)+(0.d0,0.d0))
           curves(isplit-firstecspoint+1, jj+icurve)=refvals(jj,isplit-firstecspoint+1, ifile)
           call zgemv('N', numecspoints, numecspoints, (1.d0, 0.d0), imagkemod, numecspoints, dots(:,jj+icurve), 1, (0.d0, 0.d0), hdots(:,jj+icurve), 1)
        enddo
     enddo  !! isplit

     icurve=icurve+numfluxcurves(ifile)

  enddo  !! ifile

  if (icurve.ne.totfluxcurves) then
     call openfile()
     write(mpifileptr,*) "Error totflux ", icurve, totfluxcurves
     call closefile();     call mpistop()
  endif
  if (numcalled.eq.1) then
     open(841,file="Flux.Bin", status="unknown", form="unformatted")
     write(841) numecspoints, fluxtimestep, totfluxcurves
     close(841)
  endif
  open(841,file="Flux.Bin", status="old", form="unformatted", position="append")
  write(841) numcalled, thistime
  write(841) dots, hdots
  close(841)

end subroutine flux
  

!! called with action=14, post processing

module fluxftmod
  implicit none

!!030711  complex*16, allocatable :: alldots(:,:,:), allhdots(:,:,:)
  DATATYPE, allocatable :: alldots(:,:,:), allhdots(:,:,:)

!! total number of data points across all flux files (Flux.Bins)
  integer :: totnumflux=0, eachnumflux(11)=0
  real*8 :: firstfluxtime(11)=0.d0
  real*8 :: lastfluxtime(11)=0.d0
  integer :: totfluxcurves=0
  complex*16, allocatable :: dots(:,:),  hdots(:,:), ftarray(:,:), ftout(:,:)
end module


subroutine fluxft()
  use parameters
  use fluxftmod
  use mpimod
  implicit none

  integer :: i, myiostat, flag, numdata, icurve, itime, jtime, numecspoints, timediff, jbin, jcurve
  integer :: innumecspoints, intotfluxcurves, numcalled
  real*8 :: influxtimestep, thistime, sum
  DATATYPE :: dot
  real*8, allocatable :: fluxsum(:,:)
  complex*16, allocatable :: fftrans(:,:), fftranssum(:)
  complex*16 :: tdpotft
  real*8 ::  myft, ffsum, estep
  character (len=2) :: ext(20) = (/ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20' /)

  if (numfluxreadfiles.gt.11) then
     call openfile()
     write(mpifileptr,*) "Redim fluxft"
     call closefile();     call mpistop()
  endif

!! get size of files

  do i=1,numfluxreadfiles
     open(841,file=fluxreadfilenames(i), status="old", form="unformatted", iostat=myiostat)     
     if (myiostat.ne.0) then
        call openfile()
        write(mpifileptr,*) "Error, flux file : ", fluxreadfilenames(i)
        write(mpifileptr,*) "IOSTAT= ", myiostat
        call closefile()
        call mpistop()
     endif
     if (i.eq.1) then
        read(841) numecspoints, fluxtimestep, totfluxcurves
     else
        read(841) innumecspoints, influxtimestep, intotfluxcurves
        if ((innumecspoints.ne.numecspoints).or.(influxtimestep.ne.fluxtimestep).or.(intotfluxcurves.ne.totfluxcurves)) then
           call openfile()
           write(mpifileptr,*) "Error, flux files do not agree: file number ",i
           write(mpifileptr,*) "Totfluxcurves : ", totfluxcurves, intotfluxcurves
           write(mpifileptr,*) "Fluxtimestep :  ", fluxtimestep, influxtimestep
           write(mpifileptr,*) "Numecspoints :  ", numecspoints, innumecspoints
           call closefile()
           call mpistop()
        endif
     endif
     close(841)
  enddo
  
  allocate(dots(numecspoints,totfluxcurves), hdots(numecspoints,totfluxcurves))

  totnumflux=0
  flag=0
  do i=1,numfluxreadfiles

     eachnumflux(i)=0

     open(841,file=fluxreadfilenames(i), status="old", form="unformatted", iostat=myiostat)     

     myiostat=0
     read(841,iostat=myiostat) innumecspoints, influxtimestep, intotfluxcurves

     do while (myiostat.eq.0) 
        if (myiostat==0) then
           read(841,iostat=myiostat) numcalled, thistime
           if (flag.eq.0) then
              flag=1
              firstfluxtime(i)=thistime
           endif
           if (myiostat==0) then
              read(841,iostat=myiostat) dots, hdots
              if (myiostat==0) then
                 totnumflux=totnumflux+1
                 eachnumflux(i)=eachnumflux(i)+1
              endif
           endif
        endif
        lastfluxtime(i)=thistime
     enddo
     close(841)
     call openfile()
     write(mpifileptr,*) "Total number of records on ", fluxreadfilenames(i), " is ", eachnumflux(i)
     write(mpifileptr,*) "First, last time is ", firstfluxtime(i), lastfluxtime(i)
     call closefile()
  enddo

!! points at zero will be duplicated

  totnumflux=totnumflux-numfluxreadfiles+1

  call openfile()
  write(mpifileptr,*) "Now allocating alldots etc."
  call closefile()
  
  allocate(alldots(numecspoints,totnumflux,totfluxcurves), allhdots(numecspoints,totnumflux,totfluxcurves), ftarray(-totnumflux+1:totnumflux-1,totfluxcurves), ftout(-totnumflux+1:totnumflux-1,totfluxcurves))

  ftarray=0.d0

  allocate(fluxsum(totnumflux,totfluxcurves))

  call openfile()
  write(mpifileptr,*) "Allocated.  Now reading."
  call closefile()

  jtime=0
  do i=1,numfluxreadfiles

!! points at 0 will be duplicated
     if (i.gt.1) then
        jtime=jtime-1
     endif
     call openfile()
     write(mpifileptr,*) "File: ", fluxreadfilenames(i)
     call closefile()
     open(841,file=fluxreadfilenames(i), status="old", form="unformatted", iostat=myiostat)     
     read(841,iostat=myiostat) innumecspoints, influxtimestep, intotfluxcurves
     do itime=1,eachnumflux(i)
        jtime=jtime+1
        if (jtime.gt.totnumflux) then
           print *, "Totnumflux error."
           stop
        endif
        read(841,iostat=myiostat) numcalled, thistime
        read(841,iostat=myiostat) alldots(:,jtime,:), allhdots(:,jtime,:)
     enddo
  enddo

  call openfile()
  write(mpifileptr,*) "Done reading.  Getting dots."
  call closefile()

  fluxsum=0.d0

  do itime=1, totnumflux
     jcurve=0
     jbin=1
     sum=0.d0

     do icurve=1,totfluxcurves

        sum=sum + dot(alldots(:,itime,icurve), allhdots(:,itime,icurve),numecspoints) * fluxtimestep
        jcurve=jcurve+1
        if (jcurve.eq.nucfluxcondense(jbin)) then
           if (itime.eq.1) then
              fluxsum(itime,jbin) = sum
           else
              fluxsum(itime,jbin) = fluxsum(itime-1,jbin)  + sum
           endif
           jbin=jbin+1
           jcurve=0
           sum=0.d0
        endif
     enddo
  enddo

  open(177,file="Fluxsum.Dat", status="unknown")
  do itime=1,totnumflux
     write(177,'(F10.5,100F18.12)') (itime-1)*fluxtimestep, fluxsum(itime,1:jbin-1)
  enddo
  close(177)

  do timediff = (-totnumflux+1), totnumflux-1
     do itime=max(1,1-timediff), min(totnumflux,totnumflux-timediff)
        jtime=itime+timediff
        do icurve=1,totfluxcurves
           ftarray(timediff,icurve) = ftarray(timediff,icurve) + &
                dot(alldots(:,itime,icurve), allhdots(:,jtime,icurve),numecspoints) * fluxtimestep
        enddo
     enddo
  enddo

  call openfile()
  write(mpifileptr,*) "Done dots. Calling FT."
  call closefile()

  numdata=2*totnumflux-1
  allocate(fftrans(numdata,totfluxcurves), fftranssum(numdata))

  Estep=2*pi/fluxtimestep/numdata   *   0.5d0 

  write(*, *) 

  do icurve=1,totfluxcurves
     call fluxcall(numdata,ftarray(:,icurve),fluxtimestep, eground, fftrans(:,icurve))  !! returns fftrans
  enddo

  jcurve=0
  jbin=1
  fftranssum=0.d0

  do icurve=1,totfluxcurves
     if (jcurve.eq.0) print *, "Binning curve ", icurve, " through..."
     jcurve=jcurve+1
     fftranssum=fftranssum+fftrans(:,icurve)
     if (jcurve.eq.nucfluxcondense(jbin)) then
        print *, "   .... curve ", icurve
        if (jbin.gt.20) then
           call openfile()
           write(mpifileptr,*) "REDIM FLUX!!!"
           call closefile();           call mpistop()
        endif
        ffsum=0.d0
        do i=1,numdata
           ffsum=ffsum+abs(fftranssum(i))
        enddo
        ffsum=ffsum*Estep
        write(*, *) "FLUX SUM, bin ", jbin," : ", ffsum/pi/2.d0

        open(171,file="Flux.Dat"//ext(jbin)//"",status="unknown")
        write(171,*) "#   ", numdata
        do i=1,numdata
           if (tdflag==0) then
              if (velflag==0) then
                 myft=((i-1)*Estep)
              else
                 myft=1.d0/((i-1)*Estep)
              endif
           else
              myft=4.d0/abs(tdpotft((i-1)*Estep)**2)/((i-1)*Estep)
           endif
           write(171,'(F18.12, T22, 400E20.8)')  (i-1)*Estep, abs(fftranssum(i)), fftranssum(i), abs(fftranssum(i)) * myft, myft
        enddo
        close(171)
        jcurve=0
        jbin=jbin+1
        fftranssum=0.d0
     endif
  enddo

end subroutine fluxft


!! reads MCSCF vectors at start of run

subroutine fluxread()
  use parameters
  use fluxmod
  use commonmod
  implicit none

  integer ::  inmcnumr, innumerad, inlbig, inmbig,  myiostat, ifile, inmcscfnum, i,j, ipoint, dummyint(100)
  integer :: mynuccharge1, mynuccharge2, myndof, mymrestrictflag, myrestrictflag, mynumr, mymrestrictval, myrestrictms, myspinrestrictval, myspfrestrictflag
  real*8 :: inmcstart, inmcfinish, rsum
  DATATYPE, allocatable :: tempconfig(:,:)
  DATAECS :: tempvals(numecspoints)  !! AUTOMATIC

  call openfile()
  write(mpifileptr,*) " READING MCSCF VECTORS FOR FLUX!! "
  call closefile()

  do ifile=1,numfluxfiles

     open(9954,file=fluxfilenames(ifile), form="unformatted", status="old", iostat=myiostat)
     call checkiostat(myiostat)

!! dummy reads for now.
     read(9954) mynuccharge1, mynuccharge2, myndof
!! 052412: myallspinproject -> numr (need numr not allspinproject)
     read(9954) mymrestrictflag, myspfrestrictflag, myrestrictflag, mynumr
     read(9954) mymrestrictval, myrestrictms, myspinrestrictval

     read(9954) i,j   !! numrefconfig, numrefspf already read
     read(9954) dummyint(1:nspf)
     read(9954) inmcscfnum, inmcnumr, inmcstart, inmcfinish, innumerad, inlbig, inmbig

     call openfile()
     write(mpifileptr,*) " Number of mcscf curves on file= ", inmcscfnum
     call closefile()
     if (inmcscfnum.lt.numfluxcurves(ifile)) then
        call openfile()
        write(mpifileptr,*) " Error, not enough mcscf curves; numfluxcurves is ", numfluxcurves(ifile)
        call closefile()
        call mpistop()
     endif
     if ((inmcnumr.ne.numecspoints) .or. (abs(inmcstart-bondpoints(firstecspoint)).gt.1.d-7) .or. (abs(inmcfinish-bondpoints(numr)).gt.1.d-7) .or. (inlbig.ne.lbig) .or. (inmbig .ne. mbig)) then
        call openfile()
        write(mpifileptr,*) "MCSCF vectors on file do not agree with current calculation.  Data on file, current:"
        write(mpifileptr,*) "numecspoints:    ", inmcnumr, numecspoints
        write(mpifileptr,*) "first ecs point: ", inmcstart, bondpoints(firstecspoint)
        write(mpifileptr,*) "last ecs point:  ", inmcfinish, bondpoints(numr)
        write(mpifileptr,*) "lbig:            ", inlbig, lbig
        write(mpifileptr,*) "mbig:            ", inmbig, mbig
        call closefile()
        call mpistop()
     endif
     if (innumerad.gt.numerad) then
        call openfile()
        write(mpifileptr,*) "Error, numerad on file is bigger:         ", innumerad, numerad
        call closefile()
        call mpistop()
     endif
     if (innumerad.lt.numerad) then
        call openfile()
        write(mpifileptr,*) "WARNING, numerad on file is smaller:         ", innumerad, numerad
        call closefile()
     endif
     
!!$  for reference- 
!!$  allocate(  refconfiglists(ndof,maxnumrefconfig,numfluxfiles) )
!!$  allocate(  refvects(maxnumrefconfig, maxfluxcurves, numecspoints, numfluxfiles) )
!!$  allocate(  refvals(maxfluxcurves,numecspoints,numfluxfiles) )
!!$  allocate(  refspfs(spfsize, maxnumrefspf, numecspoints, numfluxfiles) )
 
     read(9954)  refconfiglists(:,1:numrefconfig(ifile),ifile)

     refspfs(:, 1:numrefspf(ifile), :, ifile) = 0.d0

     allocate(tempconfig(numrefconfig(ifile),numecspoints))

     do ipoint=1,numecspoints
        read(9954)  rsum  ! dummy
        read(9954)  dummyint(1:inmcscfnum)

!!THIS ONE.  spf_reads assume spf variable is spfsize big.  integer variables
!!           are size on file.  DONE!!xx

#ifdef REALGO
        call one_spf_real_read(   9954,innumerad,lbig+1,mbig,numrefspf(ifile),refspfs(:,:,ipoint,ifile),maxnumrefspf,myiostat)
#else
        call one_spf_complex_read(9954,innumerad,lbig+1,mbig,numrefspf(ifile),refspfs(:,:,ipoint,ifile),maxnumrefspf,myiostat)
#endif
        call checkiostat(myiostat)

!!!        read(9954)  refspfs(1:innu merad, :, :, 1:numrefspf(ifile), ipoint, ifile) 

        do i=1,inmcscfnum
           read(9954)  tempvals(ipoint)
           read(9954)  tempconfig(:,ipoint)
           do j=1,numfluxcurves(ifile)
              if (whichfluxcurves(j,ifile).eq.i) then
                 refvects(1:numrefconfig(ifile), j, ipoint, ifile) = tempconfig(:,ipoint)
                 refvals(j, ipoint, ifile) = tempvals(ipoint)
              endif
           enddo
        enddo
     enddo
     close(9954)
     deallocate(tempconfig)
  end do  ! ifile

  call openfile()
  write(mpifileptr,*) " DONE FLUX READ! "
  call closefile()


!!$  print *, "READ CHECK, num mcscf= ", mcscfnum
!!$  do ir=1,numecspoints
!!$     print *, "RVAL ", ir
!!$     do i=1,mcscfnum
!!$        do j=1,mcscfnum
!!$           print *, i,j,dot(mcscfvects(:,i,ir),mcscfvects(:,j,ir), numconfig)
!!$        enddo
!!$     enddo
!!$  enddo
!!$  stop

end subroutine fluxread
