

!!!  EXTRA COORINATE SPECIFIC ROUTINES INCLUDING DIATOMSTATS FOR EXPECTATION VALUES (DiatomStats.Dat)


#include "Definitions.INC"

!! boxes subroutine getmyparams(inmpifileptr,inpfile,spfdims,outnumsh ells,outdim)

subroutine getmyparams(inmpifileptr,inpfile,spfdims,spfdimtype,reducedpotsize,outnumr,outnucrepulsion,nonuc_checkflag)
  use myparams
  implicit none

  integer :: spfdims(3),spfdimtype(3),inmpifileptr,nonuc_checkflag,myiostat, reducedpotsize, outnumr,&
       nargs, idim ,i,j,len,getlen
  character :: inpfile*(*)
  real*8 :: outnucrepulsion, rsq(griddim)
  character (len=200) :: buffer
  character (len=200) :: nullbuff="                                                                                "
  NAMELIST /sincparinp/        numpoints,spacing,griddim,notwoflag,coulflag,nuccharges,orblanthresh, &
       numcenters,centershift,orblanorder,toepflag,nonucrepflag,debugflag, &
       toothnbig, toothnsmall, orbparflag,num_skip_orbs,orb_skip,orblancheckmod
#ifdef PGFFLAG
  integer :: myiargc
  nargs=myiargc()
#else
  nargs=iargc()
#endif
  mpifileptr=inmpifileptr

  call getmyranknprocs(myrank,nprocs)

  nonuc_checkflag=1
  outnumr=1

  if (griddim.ne.3) then
     OFLWR "Griddim .ne. 3 not debugged.  checkme",griddim; CFLST
  endif

  open(971,file=inpfile, status="old", iostat=myiostat)

  if (myiostat/=0) then

     OFLWR "No Input.Inp found, not reading sincparams. iostat=",myiostat;CFL
  else
     read(971,nml=sincparinp)
     close(971)
  endif


  call openfile()
  do i=1,nargs
     buffer=nullbuff
     call getarg(i,buffer)
     len=getlen(buffer)
     if (buffer(1:2) .eq. 'N=') then
        read(buffer(3:len),*) numpoints(1)
        write(mpifileptr, *) "numpoints set to ", numpoints(1)
     endif
     if (buffer(1:5) .eq. 'Toep=') then
        read(buffer(6:len),*) toepflag
        write(mpifileptr, *) "toepflag set to ", toepflag
     endif
     if (buffer(1:6) .eq. 'SDebug') then
        if (.not.(buffer(1:7) .eq. 'SDebug=')) then
           write(mpifileptr,*) "Please set SDebug= not SDebug as command line option"; CFLST
        endif
        read(buffer(8:len),*) debugflag
        write(mpifileptr, *) "Debugflag for sinc dvr set to ",debugflag," by command line option"
     endif
  enddo
  call closefile()


  numpoints(:)=numpoints(1); nbox(:)=1
  if (orbparflag) then
     if (nprocs*(numpoints(3)/nprocs).ne.numpoints(3)) then
        OFLWR "ACK, orbparflag is set, but numpoints is not divisible by nprocs",numpoints(3),nprocs; CFLST
     endif
     numpoints(3)=numpoints(3)/nprocs
     nbox(3)=nprocs
  endif
  
  numpoints(griddim+1:)=1

!! NBOX HARDWIRE 1 BUT FOR PAR.

  totpoints=1
  do idim=1,griddim
     totpoints=totpoints*numpoints(idim)
     gridpoints(idim)=nbox(idim)*numpoints(idim)
  enddo

  maxgridpoints=0
  do idim=1,griddim
     maxgridpoints=max(gridpoints(idim),maxgridpoints)
  enddo

  if (griddim.gt.3) then
     OFLWR "redim spfdims"; CFLST
  endif

  spfdims(:)=1
  spfdims(1:griddim)=numpoints(1:griddim)

  spfdimtype(:)=1   !! [ -inf...inf  variable range ]

  reducedpotsize=1
  do idim=1,griddim
     reducedpotsize=reducedpotsize*numpoints(idim)
  enddo

  sumcharge=0d0
  do i=1,numcenters
     sumcharge=sumcharge+nuccharges(i)
  enddo

  nucrepulsion=0d0

  if (griddim.eq.3) then
     if (nonucrepflag.eq.0) then
        do i=1,numcenters
           do j=i+1,numcenters
              rsq(:)=((centershift(:,i)-centershift(:,j))*spacing/2 )**2 
              nucrepulsion=nucrepulsion +  nuccharges(i)*nuccharges(j)/sqrt(rsq(1)+rsq(2)+rsq(3))
           enddo
        enddo
     endif
     
  endif

!  print *, "NUCREPULSION", nucrepulsion; stop

  outnucrepulsion=nucrepulsion


!! If orbparflag=.false., toepflag.ne.0, then distribute effort and allgather.
!! If orbparflag=.true., toepflag .ne.0, and MPIFLAG is set, then only have local block; no communication afterwards is needed
!! If orbparflag=.true., toepflag .eq.0, then allocate big grid and reduce (need more memory with toepflag 0)

  gridlow=1 ; gridhigh=numpoints(3); localflag=.false.; gridoffset=0
  gridsize(1:griddim)=numpoints(1:griddim)

  if (orbparflag) then
     localflag=.true.
     gridoffset=(myrank-1)*numpoints(3)
#ifdef MPIFLAG
     if (toepflag.eq.0.or.debugflag.eq.9886) then
#endif
        gridlow=(myrank-1)*numpoints(3)+1; gridhigh=myrank*numpoints(3)
        gridsize(1:griddim)=gridpoints(1:griddim)
        localflag=.false.
#ifdef MPIFLAG
     endif
#endif

  else
     do j=1,3
        if (nbox(j).ne.1) then
           OFLWR "FDSsdffdfds"; CFLST
        endif
     enddo
  endif

  if (toepflag.eq.0.and.orbparflag) then
     OFLWR "Toepflag=0 and orbparflag=.TRUE. not allowed"; CFLST
  endif

end subroutine getmyparams


subroutine printmyopts()
  use myparams !! mpifileptr
  implicit none
  integer :: ii

!!no  OFLWR
  WRFL "********  SINC DVR PARAMS ******** " 
  WRFL "griddim,spacing",griddim,spacing
  WRFL "numpoints",numpoints(1:griddim)
  WRFL "numcenters",numcenters
  do ii=1,numcenters
     WRFL "nuccharge",nuccharges(ii)
     WRFL "centershift",centershift(:,ii)*0.5d0
  enddo
  WRFL "orblanorder,orblanthresh",orblanorder,orblanthresh
  WRFL "notwoflag,coulflag",notwoflag,coulflag
  WRFL "toepflag,nonucrepflag",toepflag,nonucrepflag
  WRFL "toothnbig, toothnsmall",toothnbig, toothnsmall
  WRFL "orbparflag",orbparflag
  WRFL "  -----  "
  WRFL "NBOX ", nbox(1:griddim)
  WRFL "totpoints",totpoints
  WRFL "**********************************"
!!no  CFL
  
end subroutine printmyopts

subroutine checkmyopts()
end subroutine checkmyopts



function cylindricalvalue() !radpoint, thetapoint,nullrvalue,mvalue, invector)
  DATATYPE :: cylindricalvalue
  cylindricalvalue=0d0
print *, "DOME CYLINDRICAL VALUE SINC (what?)"; stop
end function cylindricalvalue

subroutine get_maxsparse() !nx,ny,nz,xvals,yvals,zvals, maxsparse,povsparse)
  print *, "DOME get_maxsparse SINC"; STOP
end subroutine get_maxsparse
subroutine get_sphericalsparse() !nx,ny,nz,xvals,yvals,zvals, maxsparse,sparsetransmat,sparsestart,sparseend,sparsexyz,povsparse)
  print *, "DOME get_sphericalsparse SINC"; STOP
end subroutine get_sphericalsparse
