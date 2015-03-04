
!! BORN OPPENHEIMER CURVE OUTPUT (PLOTS) FOR NONADIABATIC CALC, AND FOR NUCLEAR FLUX CALCULATION (ACT=14)

#include "Definitions.INC"

!!!!   FOR proper projection on CURVES:  SHOULD REALLY do a true c norm by calculating the c-norm overlap of orbitals.
!!!         ---> true c-norm overlap could be used to verify whether state is accurate.

module curvemod
  use lantypemod
  implicit none
  DATATYPE, allocatable, target :: lanczosdots(:,:,:,:),  lanczoscurves(:,:,:,:)
  Type(lantype), save :: lancurveptr
  integer, allocatable :: korders(:)
  real*8, allocatable :: nnorms(:)    !! fine: starts hermitian krylov recursion.
end module curvemod


subroutine curvealloc()
  use curvemod
  use parameters
  implicit none
  integer :: mmaxord
  if (curvespinproject==1) then
     if (numcurves.gt.spinrank) then
        call openfile();        write(mpifileptr,*) "Resetting numcurves to spinrank= ", spinrank;        call closefile()
        numcurves=spinrank
     endif
  else
     if (numcurves.gt.numconfig) then
        call openfile();        write(mpifileptr,*) "Resetting numcurves to numconfig= ", numconfig;        call closefile()
        numcurves=numconfig
     endif
  endif
  allocate( lanczosdots(numcurves,numr, 2, numprop),  lanczoscurves(numcurves,numr, 2, numprop) )
  if (curvespinproject==0) then
     mmaxord=numconfig
  else
     mmaxord=spinrank
  endif
  call lanallocate(curveorder,numconfig,lancurveptr,0, mmaxord,1)
  allocate(    nnorms(numr) , korders(numr))
end subroutine

subroutine curvedealloc()
  use curvemod
  use parameters
  implicit none
  deallocate( lanczosdots,  lanczoscurves)
  deallocate(  nnorms , korders ) 
  call landeallocate(lancurveptr,0)
end subroutine





subroutine lancurves( matrix_ptr, avectorin, spfin,  iprop, thistime)
  use parameters
  use mpimod
  use curvemod
  use configptrmod
  implicit none

#ifndef REALGO
#ifdef ECSFLAG
  integer, parameter :: cflag=1
#else
  integer, parameter :: cflag=0
#endif
#else
  integer, parameter :: cflag=0
#endif
  DATAECS :: oldvalue, newvalue
  DATATYPE :: hermdot 
  integer ::  flag, startorder,endorder,  isplit,iflag,  iprop,  xxcall=0, kdoneflag, convergedflag  
  Type(CONFIGPTR) :: matrix_ptr
  DATATYPE :: avectorin(numconfig,numr), spfin(spfsize,nspf)
  DATATYPE, pointer :: dots(:,:),  curves(:,:)
  real*8 :: thistime, nextran  !! BE CAREFUL WITH NEXTRAN AND MPI !!   ALL PROCS MUST BE SYNCED !!
  integer, parameter :: printflag=0, checkstep=10
  integer :: conflag, boflag, pulseflag, nucflag,  hamway, i  

!! TO BE SAFE.  DON'T MESS.

  call rand_init(-999957.899d0)
  if (iprop.gt.100) then
     call openfile();     write(mpifileptr,*) "REDIM LC!!!"
     call closefile();     call mpistop()
  endif
  if (curvespinproject==1) then
     if ((spinwalkflag/=1)) then
        call openfile();        write(mpifileptr,*) "Need spinwalkflag=1 for spin projection";        call closefile();        call mpistop()
     endif
  endif
  if (numprop.gt.20) then
     call openfile();     write(mpifileptr,*) "REDIM";     call closefile();     call mpistop()
  endif
  if ((sparseconfigflag/=1)) then
     call openfile();     write(mpifileptr,*) "Must call lancurves with sparseconfigflag=1"
     call closefile();     call mpistop()
  endif

  if (xxcall==0) then
     call system("mkdir LanCurves") 
  endif
  xxcall=xxcall+1


  if (printflag.eq.1) then
     call openfile()
     write(mpifileptr,*) "Lan curves start" , myrank
     call closefile()
  endif

  do hamway=0,2,2   !  Plot H_0 eigencurves or those of H_0 + dipole(t)
     dots => lanczosdots(:,:,hamway/2+1,iprop)
     curves => lanczoscurves(:,:,hamway/2+1,iprop)
     korders=0
     do isplit=1,numr
        iflag=0
        if (mod(isplit-1,nprocs)+1/=myrank) then
           iflag=1
        endif
        if (iflag==0) then

           kdoneflag=0
           convergedflag=0
           oldvalue=1.d+10
           newvalue=1.d+8
           lancurveptr%initvector(:,1)=avectorin(:,isplit)

           !! add small amount of static!!   If you change this, must account for zero vector in call to lansetup which now would produce a (caught) error.

           if (nolanstaticflag/=1) then
              do i=1, numconfig
                lancurveptr%initvector(i,1)=lancurveptr%initvector(i,1) + 1.d-6/real(numconfig*numr,8) * nextran()
              enddo
           endif
           if (curvespinproject==1) then
              call configspin_projectone(lancurveptr%initvector(:,1))
           endif

           nnorms(isplit)=real(sqrt(hermdot(lancurveptr%initvector(:,1), &    !! yes, good, herm normalize the first krylov vect
                lancurveptr%initvector(:,1),numconfig)) ,8)
           lancurveptr%initvector(:,1)=lancurveptr%initvector(:,1)/nnorms(isplit)

           flag=0
           startorder=1
           endorder=checkstep

           do while (flag.eq.0)
              if (endorder.ge.lancurveptr%maxlanorder) then
                 endorder=lancurveptr%maxlanorder
                 flag=1
              endif
              conflag=0
              nucflag=0
              boflag=1
              if (hamway==0) then
                 pulseflag=0
              else
                 pulseflag=1
              endif
              call lansetup( startorder, endorder, matrix_ptr,cflag, printflag,1,0,pulseflag,0, korders(isplit), kdoneflag, 0, 0, lancurveptr, curvespinproject, 1 ,isplit)   !! hamway:  2= LIPs   0=H_0 only

              if (convergedflag/=1) then
                 if (kdoneflag==1) then
                    call openfile()
                    write(mpifileptr,*) "In curves, Krylov space linearly dep, assume convergence at order ", korders(isplit), endorder, "  hamway is  ", hamway, " for isplit= ",isplit
                    call closefile()
                    convergedflag=1
                 endif

!!080810 BAK TO 0612 (not using for flux... normalization not a huge deal)
                call CONFIGEIG(lancurveptr%lanham(:,:),korders(isplit),lancurveptr%maxlanorder,lancurveptr%rightlaneigvects,lancurveptr%lanvalues)

                 if (printflag.eq.1) then
                    call openfile()
                    write(mpifileptr, *) "Lan values, split ", isplit, " order ", korders(isplit)," : "

                    write(mpifileptr,'(10F20.12)') lancurveptr%lanvalues(1:min(korders(isplit),numcurves))
                    call closefile()
                 endif
                 if (korders(isplit).ge.numcurves) then
                    oldvalue=newvalue
                    newvalue=lancurveptr%lanvalues(numcurves)
                 endif
                 if (abs(newvalue-oldvalue).lt.lanthresh) then
                    convergedflag=1
                 endif
              endif !! convergedlfag
              if (convergedflag==1) then
                 flag=1
              endif
              startorder=startorder+checkstep
              endorder=endorder+checkstep

           enddo   !! DO WHILE FLAG

           curves(1:min(korders(isplit), numcurves),isplit)=lancurveptr%lanvalues(1:min(korders(isplit), numcurves))

!! dot of rightlaneigvect with initvec

           dots(1:min(korders(isplit),numcurves),isplit)=CONJUGATE((0.d0,0.d0)+lancurveptr%rightlaneigvects(1,1:min(korders(isplit),numcurves)))*nnorms(isplit)

        endif !! iflag
     enddo  !! isplit

     if (printflag.eq.1) then
        call openfile();     write(mpifileptr,*) "Lan curves done; bcast" , myrank;     call closefile()
     endif
        do isplit=1,numr
           call mympiibcast(korders(isplit),    mod(isplit-1,nprocs)+1     , 1)
           call mympibcast(curves(:,isplit),    mod(isplit-1,nprocs)+1     , min(numcurves,korders(isplit)))
           call mympibcast(dots(:,isplit),    mod(isplit-1,nprocs)+1     , min(numcurves,korders(isplit)))
        enddo
  if (printflag.eq.1) then
     OFLWR "Lan curves done with bcast" , myrank; CFL
  endif

  end do  !! HAMWAY

!!$     if ((iflag /= 0).and.(curveorder.ne.1)) then
!!$        call openfile()
!!$        write(mpifileptr, *) "WARNING, lancurves not converged ", oldvalue, newvalue, numcurves, curveorder
!!$        call closefile()
!!$     endif

end subroutine lancurves


subroutine lancurveplot( iprop, thistime)
  use parameters
  use mpimod
  use curvemod
  use configptrmod
  use commonmod
  implicit none

  integer ::   isplit,iflag,  iprop, iii, nn, getlen2, ilen
  DATATYPE, pointer :: dots(:,:),  curves(:,:)
  real*8 :: thistime, rsum
  real*8, save :: mmax, mmin
  integer :: calledflag(20)=0, hamway, i
  character (len=40) :: filename, filename2, sysline
  character(len=256) :: cwdline  !! cwd
  character(len=256) :: nnbb="                                                                                                  &
&                                                                                                                                                              "

  if (myrank.ne.1) then
     return
  endif

  calledflag(iprop)=calledflag(iprop)+1

  do hamway=0,2,2

     dots => lanczosdots(:,:,hamway/2+1,iprop)
     curves => lanczoscurves(:,:,hamway/2+1,iprop)

     if (hamway.eq.0) then
        iflag=0
        if (calledflag(iprop).lt.10) then
           write(filename,'(A10,I1,A1,I1,A5)') "LanCurves/",calledflag(iprop),"_",iprop,".dat0"
           write(filename2,'(A10,I1,A1,I1,A6)') "LanCurves/",calledflag(iprop),"_",iprop,"W.dat0"
           iflag=1
        else if (calledflag(iprop).lt.100) then
           write(filename,'(A10,I2,A1,I1,A5)') "LanCurves/",calledflag(iprop),"_",iprop,".dat0"
           write(filename2,'(A10,I2,A1,I1,A6)') "LanCurves/",calledflag(iprop),"_",iprop,"W.dat0"
           iflag=1
        else if (calledflag(iprop).lt.1000) then
           write(filename,'(A10,I3,A1,I1,A5)') "LanCurves/",calledflag(iprop),"_",iprop,".dat0"
           write(filename2,'(A10,I3,A1,I1,A6)') "LanCurves/",calledflag(iprop),"_",iprop,"W.dat0"
           iflag=1
        else if (calledflag(iprop).lt.10000) then
           write(filename,'(A10,I4,A1,I1,A5)') "LanCurves/",calledflag(iprop),"_",iprop,".dat0"
           write(filename2,'(A10,I4,A1,I1,A6)') "LanCurves/",calledflag(iprop),"_",iprop,"W.dat0"
           iflag=1
        endif
     else
        iflag=0
        if (calledflag(iprop).lt.10) then
           write(filename,'(A10,I1,A1,I1,A4)') "LanCurves/",calledflag(iprop),"_",iprop,".dat"
           write(filename2,'(A10,I1,A1,I1,A5)') "LanCurves/",calledflag(iprop),"_",iprop,"W.dat"
           iflag=1
        else if (calledflag(iprop).lt.100) then
           write(filename,'(A10,I2,A1,I1,A4)') "LanCurves/",calledflag(iprop),"_",iprop,".dat"
           write(filename2,'(A10,I2,A1,I1,A5)') "LanCurves/",calledflag(iprop),"_",iprop,"W.dat"
           iflag=1
        else if (calledflag(iprop).lt.1000) then
           write(filename,'(A10,I3,A1,I1,A4)') "LanCurves/",calledflag(iprop),"_",iprop,".dat"
           write(filename2,'(A10,I3,A1,I1,A5)') "LanCurves/",calledflag(iprop),"_",iprop,"W.dat"
           iflag=1
        else if (calledflag(iprop).lt.10000) then
           write(filename,'(A10,I4,A1,I1,A4)') "LanCurves/",calledflag(iprop),"_",iprop,".dat"
           write(filename2,'(A10,I4,A1,I1,A5)') "LanCurves/",calledflag(iprop),"_",iprop,"W.dat"
           iflag=1
        endif
     endif
     call openfile()
     write(mpifileptr,*) " --> writing lancurves, numcurves= ", numcurves
     call closefile()
     
     if (iflag==1) then
        open(888,file=filename, status="unknown")
        open(889,file=filename2, status="unknown")
        open(890,file="LanCurves/Points.dat", status="unknown")
        open(891,file="LanCurves/Params.dat", status="unknown")
        write(891,*) numcurves, numr, calledflag(iprop)
        close(891)
        do isplit=1,numr
           write(888,'(100F12.6)') real(curves(1:numcurves, isplit))
           write(889,'(100F12.6)') abs(dots(1:numcurves, isplit)**2/(bondweights(isplit)))
           write(890,'(100F12.6)') real(bondpoints(isplit),8)
        enddo
        write(888,*)
        write(889,*)
        write(890,*)
        
        close(888)
        close(889)
        close(890)
     endif
  enddo

  if (calledflag(iprop)==1) then  ! does this with last hamway curves, whatever
     mmax=-1.d+10
     mmin=1.d+10
     do i=1,numr
        if (real(curves(1,i)).lt.mmin) then
           mmin=real(curves(1,i))
        endif
        if (real(curves(1,i)).gt.mmax) then
           mmax=real(curves(1,i))
        endif
     enddo
     mmax=mmin+1.d0
     mmin=mmin-0.2d0
  endif

#ifndef NOMATHEMATICA

  open(809, file="MC-LinesQQQ.txt", status="old", iostat=iii)
  if ( iii /=0 ) then
     if (calledflag(iprop).lt.2) then
        call openfile() 
        write(mpifileptr,*) 
        write(mpifileptr,*) "  **  I did not find MC-LinesQQQ.txt in the directory, thus MC-Lan.txt will not be made"

        write(mpifileptr,*) 
        call closefile()
     endif
  else
     close(809)
     call system("cp MC-LinesQQQ.txt MC-Lan.txt")
     call system('perl -pi -e "s/KKK/Lan/" MC-Lan.txt')

     cwdline=nnbb
     write(cwdline,*) 'perl -pi -e "s/MMM/'
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A37)') "Dressed State Born Oppenheimer Curves"
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A13)') '/" MC-Lan.txt'
     call system(cwdline)

     cwdline=nnbb
     write(cwdline,*) 'perl -pi -e "s/LLL/'
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A23)') "Born Oppenheimer Curves"
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A13)') '/" MC-Lan.txt'
     call system(cwdline)

     cwdline=nnbb
     write(cwdline,*) 'perl -pi -e "s/RRR/'
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A14,I1,A4)') "MathLanCurves_",iprop,".gif"
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A13)') '/" MC-Lan.txt'
     call system(cwdline)

     cwdline=nnbb
     write(cwdline,*) 'perl -pi -e "s/RXR/'
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A15,I1,A4)') "MathLanCurves0_",iprop,".gif"
     ilen=getlen2(cwdline)+1
     write(cwdline(ilen:),'(A13)') '/" MC-Lan.txt'
     call system(cwdline)

     do nn=1,2
        cwdline=nnbb
        write(cwdline,*) 'perl -pi -e "s/XYZ/'
        ilen=getlen2(cwdline)+1
        write(cwdline(ilen:),'(F12.6)') mmax
        ilen=getlen2(cwdline)+1
        write(cwdline(ilen:),'(A13)') '/" MC-Lan.txt'
        call system(cwdline)

        cwdline=nnbb
        write(cwdline,*) 'perl -pi -e "s/ZYX/'
        ilen=getlen2(cwdline)+1
        write(cwdline(ilen:),'(F12.6)') mmin
        ilen=getlen2(cwdline)+1
        write(cwdline(ilen:),'(A13)') '/" MC-Lan.txt'
        call system(cwdline)

        cwdline=nnbb
        write(cwdline,*) 'perl -pi -e "s/AAA/'
        ilen=getlen2(cwdline)+1
        write(cwdline(ilen:),'(F9.6)') real(bondpoints(1),8)
        ilen=getlen2(cwdline)+1
        write(cwdline(ilen:),'(A13)') '/" MC-Lan.txt'
        call system(cwdline)

        cwdline=nnbb
        write(cwdline,*) 'perl -pi -e "s/BBB/'
        ilen=getlen2(cwdline)+1
        write(cwdline(ilen:),'(F9.6)') real(bondpoints(numr),8)
        ilen=getlen2(cwdline)+1
        write(cwdline(ilen:),'(A13)') '/" MC-Lan.txt'
        call system(cwdline)
     enddo

     cwdline=nnbb
     write(cwdline,*) 'perl -pi -e "s/SSS/'
     ilen=getlen2(cwdline)+2
     write(cwdline(ilen:),'(I2)') iprop
     ilen=getlen2(cwdline)+2
     write(cwdline(ilen:),*) '/" MC-Lan.txt'
     call system(cwdline)

     cwdline=nnbb
     write(cwdline,*) 'perl -pi -e "s/TTT/'
     ilen=getlen2(cwdline)+2
     if (numcurves.lt.10) then
        write(cwdline(ilen:),'(I1)') numcurves
     else 
        write(cwdline(ilen:),'(I2)') numcurves
     endif

     ilen=getlen2(cwdline)+2
     write(cwdline(ilen:),*) '/" MC-Lan.txt'
     call system(cwdline)

     cwdline=nnbb
     write(cwdline,*) 'perl -pi -e "s/UUU/'
     ilen=getlen2(cwdline)+2
     if (numr.lt.10) then
        write(cwdline(ilen:),'(I1)') numr
     else if (numr.lt.100) then
        write(cwdline(ilen:),'(I2)') numr
     else if (numr.lt.1000) then
        write(cwdline(ilen:),'(I3)') numr
     else
        write(cwdline(ilen:),'(I4)') numr
     endif

     ilen=getlen2(cwdline)+2
     write(cwdline(ilen:),*) '/" MC-Lan.txt'
     call system(cwdline)

     cwdline=nnbb
     write(cwdline,*) 'perl -pi -e "s/VVV/'
     ilen=getlen2(cwdline)+2
     if (calledflag(iprop).lt.10) then
        write(cwdline(ilen:),'(I1)') calledflag(iprop)
     else if (calledflag(iprop).lt.100) then
        write(cwdline(ilen:),'(I2)') calledflag(iprop)
     else if (calledflag(iprop).lt.1000) then
        write(cwdline(ilen:),'(I3)') calledflag(iprop)
     else if (calledflag(iprop).lt.10000) then
        write(cwdline(ilen:),'(I4)') calledflag(iprop)
     else
        write(cwdline(ilen:),'(I5)') calledflag(iprop)
     endif

     ilen=getlen2(cwdline)+2
     write(cwdline(ilen:),*) '/" MC-Lan.txt'
     call system(cwdline)

     cwdline=nnbb
     write(cwdline,*) 'perl -pi -e "s/XXX/'
     ilen=getlen2(cwdline)+2
     rsum=plotmodulus*par_timestep*0.02418884d0

     if (rsum.lt.10.d0) then
        write(cwdline(ilen:),'(F7.5)') rsum
     else if (rsum.lt.100.d0) then
        write(cwdline(ilen:),'(F8.5)') rsum
     else if (rsum.lt.1000.d0) then
        write(cwdline(ilen:),'(F9.5)') rsum
     else if (rsum.lt.10000.d0) then
        write(cwdline(ilen:),'(F10.5)') rsum
     else if (rsum.lt.100000.d0) then
        write(cwdline(ilen:),'(F11.5)') rsum
     else if (rsum.lt.1000000.d0) then
        write(cwdline(ilen:),'(F12.5)') rsum
     else 
        write(cwdline(ilen:),'(F13.5)') rsum
     endif

     ilen=getlen2(cwdline)+2
     write(cwdline(ilen:),*) '/" MC-Lan.txt'
     call system(cwdline)

     if (mod(calledflag(iprop),10).eq.1) then
        call openfile()
        !        write(mpifileptr,*) "Calling Mathematica!"
        write(mpifileptr,*) "Not calling Mathematica - do it yourself:"
        write(mpifileptr,*) '     math -noprompt -run "<<MC-Lan.txt"'
        call closefile()
        !        call system('math -noprompt -run "<<MC-Lan.txt"')
     endif

     write(sysline,'(A24,I1)') "mv MC-Lan.txt MC-Lan.txt",iprop
     call system(sysline(1:25))
  endif
     
#endif

end subroutine lancurveplot



subroutine nonlancurves(matrix_ptr, avectorin, spfin,  iprop, thistime)
  use parameters
  use mpimod
  use curvemod
  use configptrmod
  use commonmod
  implicit none

  Type(CONFIGPTR) :: matrix_ptr

  integer ::    isplit,iflag,  iprop,  xxcall=0, i, j, hamway
  DATATYPE :: avectorin(numconfig,numr), spfin(spfsize,nspf), dot
  DATATYPE, pointer :: dots(:,:),  curves(:,:)
  real*8 :: thistime
  integer, parameter :: printflag=0
  DATATYPE, allocatable :: smallconfigvects3(:,:), smallconfigvects(:,:),smallconfigvals(:)

  if ((sparseconfigflag/=0)) then
     call openfile();     write(mpifileptr,*) "Must call nonlancurves with sparseconfigflag=0";     call closefile();     call mpistop()
  endif
  if (xxcall==0) then
     call system("mkdir LanCurves") 
  endif

  xxcall=xxcall+1

  if (printflag.eq.1) then
     OFLWR "Non-Lan curves start" , myrank; CFL
  endif

  do hamway=0,2,2  !  Plot H_0 eigencurves or those of H_0 + dipole(t)
     
     dots => lanczosdots(:,:,hamway/2+1,iprop)
     curves => lanczoscurves(:,:,hamway/2+1,iprop)
     dots=0.d0
     curves=0.d0

     allocate(smallconfigvects3(numconfig,numconfig), smallconfigvects(numconfig,numconfig), smallconfigvals(numconfig))

     do isplit=1,numr
        iflag=0
        if (mod(isplit-1,nprocs)+1/=myrank) then
           iflag=1
        endif

        if (iflag==0) then

#ifndef NEWWALKS           
           smallconfigvects3(:,:) = 1.d0/bondpoints(isplit) * matrix_ptr%cpot(:,:,1) + 1.d0/bondpoints(isplit)**2 * matrix_ptr%cop(:,:,1)

           if (constraintflag/=0) then
              smallconfigvects3(:,:) = smallconfigvects3(:,:) + matrix_ptr%ccon(:,:,1)
           endif
           if (tdflag/=0) then
              smallconfigvects3(:,:) = smallconfigvects3(:,:) + matrix_ptr%cpulse(:,:,1)/bondpoints(isplit)
           endif
           if (curvespinproject==1) then
              smallconfigvects3=transpose(smallconfigvects3)
              do i=1,numconfig
                 call configspin_projectone(smallconfigvects3(:,i))
              enddo
              smallconfigvects3=transpose(smallconfigvects3)
              do i=1,numconfig
                 call configspin_projectone(smallconfigvects3(:,i))
              enddo
           endif

           call CONFIGEIG(smallconfigvects3, &
                numconfig,numconfig, smallconfigvects,  smallconfigvals)

           if (curvespinproject/=0) then
              j=0;              i=0
              do while (j.lt.numcurves)
                 i=i+1
                 if (i.gt.numconfig) then
                    call openfile()
                    write(mpifileptr,*) "SPIN CURVE ERROR? have only found ", j ," curves out of ", numcurves
                    call closefile()
                    j=numcurves
                 else
                    if (abs(smallconfigvals(i)).gt.1.d-6) then
                       j=j+1
                       curves(j,isplit)=smallconfigvals(i)
                       dots(j,isplit)=dot(smallconfigvects(:,i), avectorin(:,isplit),numconfig)
                    endif
                 endif
              enddo
           else
              curves(1:numcurves,isplit)=smallconfigvals(1:numcurves)
              do i=1,numcurves
                 dots(i,isplit)=dot(smallconfigvects(:,i), avectorin(:,isplit),numconfig)
              enddo
           endif
#else
           OFLWR "PROGME CURVES"; CFLST
#endif

        endif !! iflag
     enddo  !! isplit

     deallocate(smallconfigvects,smallconfigvects3,smallconfigvals)

     do isplit=1,numr
        call mympibcast(curves(1:numcurves,isplit),    mod(isplit-1,nprocs)+1     , numcurves)
        call mympibcast(dots(1:numcurves,isplit),    mod(isplit-1,nprocs)+1     , numcurves)
     enddo

  end do  !! HAMWAY

end subroutine nonlancurves


