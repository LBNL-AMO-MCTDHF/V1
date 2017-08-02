

!!!  EXTRA COORINATE SPECIFIC ROUTINES INCLUDING DIATOMSTATS FOR EXPECTATION VALUES (DiatomStats.Dat)


#include "Definitions.INC"

!! boxes subroutine getmyparams(inmpifileptr,inpfile,spfdims,outnumsh ells,outdim)

subroutine getmyparams(inmpifileptr,inpfile,spfdims,spfdimtype,reducedpotsize,outnumr,&
     outnucrepulsion,nonuc_checkflag)
  use myparams
  use pfileptrmod
  use pmpimod
  use onedfunmod
  implicit none

  integer,intent(in) :: inmpifileptr
  character,intent(in) :: inpfile*(*)
  integer,intent(out) :: spfdims(3),spfdimtype(3),nonuc_checkflag,reducedpotsize, outnumr
  real*8,intent(out) :: outnucrepulsion
  integer :: nargs, i,j, len,getlen,myiostat
  character (len=SLN) :: buffer
  character (len=SLN) :: nullbuff
  real*8, parameter :: pi = 3.14159265358979323844d0
  real*8 :: zz, ss , xfac
  DATATYPE ::  ns(1)
  DATATYPE, parameter :: czero(1)=0
  NAMELIST /sinc1dparinp/        numpoints,spacing,twostrength,nuccharges,orblanthresh, &
       numcenters,centershift,orblanorder,nonucrepflag,debugflag, &
       orbparflag,num_skip_orbs,orb_skip,orblancheckmod,zke_paropt,&
       capflag,capstrength,capstart,cappower,fft_ct_paropt,fft_batchdim,&
       fft_circbatchdim,maxcap,mincap,capmode, &
       scalingflag,scalingdistance,smoothness,scalingtheta,scalingstretch,&
       ivoflag, loadedocc, orbtargetflag,orbtarget,&
       toepflag,softness,twotype,harmstrength, twomode, nucstrength, &
       eigmode, coulmode, sechmode, nucmode, combinesech, nucrangefac, &
       softnesstwoe, softnessnuc

#ifdef PGFFLAG
  integer :: myiargc
  nargs=myiargc()
#else
  nargs=iargc()
#endif

  do i=1,SLN
     nullbuff(i:i)=" "
  enddo

!!! MPI INIT

  call getworldcommgroup(PROJ_COMM_WORLD,PROJ_GROUP_WORLD)
  mpifileptr=inmpifileptr
  call getmyranknprocs(myrank,nprocs)

  nonuc_checkflag=1
  outnumr=1

  open(971,file=inpfile, status="old", iostat=myiostat)

  if (myiostat/=0) then
     OFLWR "No Input.Inp found, not reading sincparams. iostat=",myiostat;CFL
  else
     read(971,nml=sinc1dparinp)
     close(971)

!! enforce defaults that depend on other variables
     softnessnuc=softness
     softnesstwoe=softness

     open(971,file=inpfile, status="old", iostat=myiostat)
     read(971,nml=sinc1dparinp)
     close(971)

  endif

  call openfile()
  do i=1,nargs
     buffer=nullbuff
     call getarg(i,buffer)
     len=getlen(buffer)
     if (buffer(1:2) .eq. 'N=') then
        read(buffer(3:len),*) numpoints
        write(mpifileptr, *) "numpoints set to ", numpoints
     endif
     if (buffer(1:9) .eq. 'Batchdim=') then
        read(buffer(10:len),*) fft_batchdim
        write(mpifileptr, *) "fft_batchdim set to ", fft_batchdim
     endif
     if (buffer(1:10) .eq. 'Circbatch=') then
        read(buffer(11:len),*) fft_circbatchdim
        write(mpifileptr, *) "fft_circbatchdim set to ", fft_circbatchdim
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


  if (myrank.eq.1.and.debugflag.ne.0) then
     call ftset(1,mpifileptr)
  else
     call ftset(0,-99)
  endif

!! matches main
  if (nprocs.eq.1) then
     orbparflag=.false.
  endif

  if (orbparflag) then
     if (mod(numpoints,nprocs).ne.0) then
        OFLWR "ACK, numpoints is not divisible by nprocs",numpoints,nprocs; CFLST
     endif
  endif  

!! NBOX HARDWIRE 1 BUT FOR PAR.

  nbox=1
  qbox=1
  if (orbparflag) then
     numpoints=numpoints/nprocs
     nbox=nprocs
     qbox=myrank
  endif

  totpoints=numpoints
 
  gridpoints=nbox*numpoints

  spfdims(:)=1
  spfdims(3)=numpoints

  spfdimtype(:)=1   !! [ -inf...inf  variable range ]

  reducedpotsize=numpoints

  sumcharge=0d0
  nucrepulsion=0
  do i=1,numcenters
     sumcharge=sumcharge+nuccharges(i)
     if (nucmode.eq.0.or.twomode.ne.0) then  !! physically motivated internuclear repulsion
        !!                                   !! always for coulomb or linear, nucmode sets choice for sech

        if (twomode.ne.0.or.combinesech.eq.0) then  !! coulomb or old version sech
        
           !! straight sum of two-particle interactions just like electrons
           
           xfac=1
           if (twomode.eq.(-1)) then
              xfac = softness/softnessnuc
           endif

           do j=i+1,numcenters        
              ns(:) =  nuccharges(i) * &
                   xfac*onedfun(xfac*spacing*DATAONE*(centershift(i:i)-centershift(j:j))/2,1,1d0,nuccharges(j))
              nucrepulsion = nucrepulsion + ns(1)
           enddo
        
        else   !! combinesech.ne.0 ::

           !! not sure what I am doing here, in defining the internuclear potential.
           !! electron-electron is pairwise.  electron-nuclear has enhancement when nuclei approach.
           !! adding pairwise nuclear-nuclear results in overly attractive and deep BO PES.
           !! so add nuclear potential that is like electron-nuclear potential.
           !! must cancel self-interation to maintain asymptotes.
           
           !! should multiply the interaction by the nuclear charge, not include it in charge product
           !! in sechsq.
           
           !! Nuclear potential is (-1) x electronic potential minus self-interaction
           !!   same as sum of 2-particle interactions for coulomb
           
           !! no
           !!   -2.25000 h2  or -4.46 heh  or -1 he2           united no two elec
           !        ns(:) = 0.5d0 * (elecpot(spacing*DATAONE*centershift(i:i)/2, 1, nuccharges(i), 0) &
           !             - onedfun(czero, 1, nuccharges(i), nuccharges(i)))
           
           !  -3.25000.. h2 united atom ground state without two-electron
           !        ns(:) = 0.5d0 * nuccharges(i) * elecpot( spacing*DATAONE*centershift(i:i)/2, 1, 1d0, i )

           ! -2.5 h2 -6.71 HeH -13 He2 without two-electron, not repulsive enough even with 2-e
           ! ns(:) = nuccharges(i) * elecpot( spacing*DATAONE*centershift(i:i)/2, 1, 1d0, i )

           ! try it again after changing skipcenter thing
           ! -1.5 h2  -5.0 he2    without     -0.905 h2  -1.19 he2     with two-e
           !  again strong binding potential for e.g. He2 away from R=0
           !  ns(:) = nuccharges(i) * elecpot( spacing*DATAONE*centershift(i:i)/2, 1, 1d0, i )
           
           !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!

           ! -2.25 h2     -11 He2    without         -1.655 h2  -7.19 He2 with two-e
           ! strong binding potential for all species away from R=0    ..bad       
           !           ns(:) = 0.5d0 * nuccharges(i) * ( elecpot(spacing*DATAONE*centershift(i:i)/2, 1, 1d0, 0 ) &
           !                - onedfun(czero, 1, 1d0, nuccharges(i)) )

           !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!

           !! -0.5000.. h2   or -0.71  heh or +3 he2   united atom without two electron
           !!                          for he2, orbitals are -8, -5.5, -2, -0.5, hmm  
           !!                          for he2, would be E=-2 without orbital promotion
           !!                          that's the same as the energy of one helium
           !!                          likewise for h2, -0.5, that's the energy of one hydrogen
           !!                          I like it.. but does not bind anything
           
           ns(:) = nuccharges(i) * ( elecpot(spacing*DATAONE*centershift(i:i)/2, 1, 1d0, 0 ) &
                - onedfun(czero, 1, 1d0, nuccharges(i)) )
           
           if (real(nucstrength*ns(1)) < -1d-12) then
              OFLWR "EEEEROROROROR ", ns(1); CFLST
           endif
        
           nucrepulsion = nucrepulsion + ns(1)
        
        endif  !! version of physical potential combinesech

     else  !! nucmode.eq.0.or.twomode.ne.0
        
        !! here, nucmode.ne.0 :: ad hoc internuclear repulsion 
        !! nuclear repulsion designed to make asymptote equal to R=0 with twostrength=0
        !! (for bosons, hypothetically)
        
        do j=i+1,numcenters        
           ns(:) = real( &
                onedfun(spacing*DATAONE*(centershift(i:i)-centershift(j:j))/2/nucrangefac,1,1d0,1d0) / &
                onedfun(czero,1,1d0,1d0),8)
           
           if (combinesech.ne.0) then !! new version sech-squared, sum then square sech
              
              !! one-elec potentials have been scaled to make united atom limit correct
              nucrepulsion = nucrepulsion + 0.5d0 * &
                   ns(1)*( (nuccharges(i)+nuccharges(j))**3 - nuccharges(i)**3 - nuccharges(j)**3 )
              
           else  !! old version sech-squared, sum of pair potentials

              if (sechmode.ne.0) then  !! sechmode 0 required here
                 OFLWR "Error, combinesech.eq.0.and.sechmode.ne.0.and.nucmode.ne.0 not allowed"
                 WRFL  "   can't sum pairwise potentials with different length scales and"
                 WRFL  "   perform ad hoc correction to get R=0 asymptote."; CFLST
              endif
              
              ss = 1d0/softness
              zz = nuccharges(i)*(nuccharges(i)+ss) + nuccharges(j)*(nuccharges(j)+ss)
              zz = 0.25d0 * ( 2*ss**2 + 4*zz - 2*ss*sqrt(ss**2+4*zz) )
              !              print * ,"EFF NUC CHARGE",sqrt(zz)
              nucrepulsion = nucrepulsion + 0.5d0 * &
                   ns(1)*( (nuccharges(i)+nuccharges(j)) * zz &
                   - nuccharges(i)**3 - nuccharges(j)**3 )
              
           endif

        enddo   !! j=1,numcenters

     endif   !! physically motivated internuclear repulsion or not, nucmode

  enddo  !! i=1,numcenters

  outnucrepulsion = nucrepulsion * nucstrength
  
  return
  
end subroutine getmyparams


subroutine printmyopts()
  use myparams
  use pfileptrmod
  implicit none
  integer :: ii

!!no  OFLWR
  WRFL "********  SINC DVR PARAMS ******** " 
  WRFL "spacing",spacing
  WRFL "numpoints",numpoints
  WRFL "********  HAMILTONIAN PARAMS ******** "
  WRFL "numcenters",numcenters
  do ii=1,numcenters
     WRFL "  nuccharge",nuccharges(ii)
     WRFL "  radius   ",centershift(ii)*0.5d0
  enddo
  WRFL "twostrength", twostrength
  WRFL "nucstrength", nucstrength
  if (twotype.eq.0) then
     WRFL "constant two-body interaction, twotype==0"
     WRFL "  SOFTNESS for one-body interactions", softness
  else
     WRFL "potential two-body interaction, twotype.ne.0"
     WRFL "  SOFTNESS ", softness
     if (twomode.eq.(-1)) then
        WRFL "  for coulomb, internuclear softness ", softnessnuc
        WRFL "  for coulomb, interelectronic softness ", softnesstwoe
     endif
  endif
  if (twomode.eq.0) then
     WRFL "sech^2 potential, twomode==0"
     if(sechmode.eq.0) then
        WRFL "   3d hydrogenic eigvals"
     else
        WRFL "   half-integer even parity"
     endif
  elseif (twomode.eq.1) then
     WRFL "soft coulomb potential, twomode==1"
     if (coulmode.eq.0) then
        WRFL "   coulmode == 0 : centrifugal with hydrogenic eigvals"
     elseif (coulmode.eq.1) then
        WRFL "   coulmode == 1 : centrifugal with half integer even parity"
     endif
  else
     WRFL "linear potential"
  endif
  WRFL "********  OTHER PARAMS ********** "
  WRFL "orblanorder,orblanthresh",orblanorder,orblanthresh
  WRFL "nonucrepflag",nonucrepflag
  WRFL "orbparflag",orbparflag
  WRFL "zke_paropt",zke_paropt
  WRFL "   --> fft_ct_paropt",fft_ct_paropt
  WRFL "fft_batchdim",fft_batchdim
  WRFL "fft_circbatchdim",fft_circbatchdim
  WRFL "  -----  "
  WRFL "NBOX ", nbox
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
