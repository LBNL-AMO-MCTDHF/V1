
!! ALL MODULES

!! CONFIGURATION SUBROUTINES.  SLATER DETERMINANTS NOT SPIN EIGFUNCTS.
!!  SEE WALKS.F90 FOR MORE; SPIN.F90 FOR SPIN PROJECTION.

#include "Definitions.INC"


module configsubmod
contains

subroutine printconfig(thisconfig,www)
  use walkmod
  use fileptrmod
  implicit none

  type(walktype),intent(in) :: www
  integer,intent(in) :: thisconfig(www%num2part)
  character (len=4) :: mslabels(2) =["a ","b "]
  integer :: i

  write(mpifileptr,'(100(I3,A2))') (thisconfig((i-1)*2+1), &
       mslabels(thisconfig(i*2)), i=1,www%numpart)

end subroutine printconfig


!! ORBITAL ORDERING:  RETURNS SPINORBITAL INDEX GIVEN ORBITAL INDEX AND SPIN
!! (SPIN IS 1 (ALPHA) OR 2 (BETA))

function iind(twoarr)
  implicit none
  integer :: iind
  integer,intent(in) :: twoarr(2)

!!$  if (orderflag==1) then
!!$     iind=twoarr(1) + (twoarr(2)-1)*nspf 
!!$  else

     iind=(twoarr(1)-1)*2 + twoarr(2)

!!$  endif
end function


!! RETURN CONFIGURATION INDEX GIVEN ORBITAL OCCUPANCIES

function getconfiguration(thisconfig,www)
  use fileptrmod
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: thisconfig(www%num2part)
  integer :: getconfiguration,  j,flag,k, dir,newdir, step,aa,bb, ii,kk,jj,flag1,flag2

  getconfiguration=-1

  dir=1;  j=1;  step=max(1,www%numconfig/4);  flag=0

  do while (flag.eq.0)
     flag=1
     do k=1,www%numpart

!!$SP        aa=www%configlist(2*k-1,www%configorder(j));
        aa=www%configlist(2*k-1,j);

        bb=thisconfig(2*k-1)

        if (aa .ne. bb) then
           flag=0
           if (aa.lt.bb) then
              newdir=1
           else
              newdir=-1
           endif
           if (newdir.ne.dir) then
              step=max(1,step/2)
           endif
           dir=newdir;           j=j+dir*step
           if (j.le.1) then
              j=1;              dir=1
           endif
           if (j.ge.www%numconfig) then
              j=www%numconfig;              dir=-1
           endif
           exit
        endif
     enddo
  enddo


     ii=j

!!$   DON'T HAVE MAXSPINSETSIZE YET.
!!$   otherwise would be do j=max(1,ii-maxspinsetsize),min(numconfig,ii+maxspinsetsize)

     kk=(-1); flag1=0; flag2=0
     do while (flag1.eq.0.or.flag2.eq.0)
        kk=kk+1
        do jj=0,1
           j=ii+kk*(-1)**jj
           if (jj.eq.0.and.j.gt.www%numconfig) then
              flag1=1
              cycle
           endif
           if (jj.eq.1.and.j.lt.1) then
              flag2=1
              cycle
           endif

           flag=1
           do k=1,www%num2part
!!$SP          if (thisconfig(k).ne.www%configlist(k,www%configorder(j))) then
              if (thisconfig(k).ne.www%configlist(k,j)) then
                 flag=0
                 exit
              endif
           enddo
           if (flag.eq.1) then
!!$SP              getconfiguration=www%configorder(j)
              getconfiguration=j
              return
           endif
        enddo
     enddo
     OFLWR
     call printconfig(thisconfig,www)
     WRFL "NEWGETCONFIG NEWCONFIG ERROR"; CFLST

end function getconfiguration



!! RETURN CONFIGURATION INDEX GIVEN ORBITAL OCCUPANCIES

!! PUTS AN UN ORDERED CONFIGURATION INTO PROPER (INDEX INCREASING) ORDER
!! AND RETURNS SIGN OF PERMUTATION

function reorder(thisconfig,numpart)
  implicit none
  integer, intent(in) :: numpart
  integer, intent(inout) :: thisconfig(1:2*numpart)
  integer :: reorder, phase, flag, jdof, temporb(2)

  phase=1;  flag=0
  do while (flag==0)
     flag=1
     do jdof=1,numpart-1
        if (iind(thisconfig((jdof-1)*2+1:jdof*2)) .gt. &
             iind(thisconfig(jdof*2+1:(jdof+1)*2))) then
           phase=phase*(-1)
            flag=0
           temporb = thisconfig((jdof-1)*2+1:jdof*2)
           thisconfig((jdof-1)*2+1:jdof*2) = thisconfig(jdof*2+1:(jdof+1)*2)
           thisconfig(jdof*2+1:(jdof+1)*2) = temporb
        endif
     enddo
  enddo
  reorder=phase
end function


!! INDICATES WHETHER A GIVEN CONFIG IS IN OUR CONFIGURATION LIST


function allowedconfig0(www,thisconfig,in_df)
  use basis_parameters
  use walkmod
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: thisconfig(www%num2part), in_df
  integer :: i, isum, j, numwithinshell(numshells), numinshell(numshells), ishell, ii, k
  logical :: allowedconfig0, tempflag

  do j=1,www%numpart
     do i=j+1,www%numpart
        if ( (thisconfig(2*i-1)==thisconfig(2*j-1)).and.(thisconfig(2*i)==thisconfig(2*j)) ) then
           allowedconfig0=.false.;           return
        endif
     enddo
  enddo
  do j=1,www%numpart
     i=thisconfig(j*2-1)
     if ((i.lt.1).or.(i.gt.www%nspf)) then
        allowedconfig0=.false.;        return
     endif
  enddo
  if (restrictflag==1) then    ! by m_s
     isum=0
     do i=2,www%num2part,2
        isum=isum+(thisconfig(i)*2-3)
     enddo
     if (isum /= www%restrictms) then
        allowedconfig0=.false.
        return
     endif
  end if
  if (spfrestrictflag==1) then
     if (www%m_restrictflag==1.or.www%m_restrictmin.gt.-99999.or.www%m_restrictmax.lt.99999) then
        isum=getmval(www,thisconfig)
        if ((www%m_restrictflag==1.and.isum /= www%m_restrictval).or.&
             isum.lt.www%m_restrictmin.or.isum.gt.www%m_restrictmax) then
           allowedconfig0=.false.;        return
        endif
     endif
  end if
  if ((spfugrestrict==1).and.(www%ug_restrictflag==1)) then
     isum=getugval(www,thisconfig)
     if (isum /= www%ug_restrictval) then
        allowedconfig0=.false.
        return
     endif
  end if

  
  numwithinshell(:)=0;      numinshell(:)=0
  do ishell=1,numshells
     do i=allshelltop(ishell-1)+1,allshelltop(ishell)  !! spatial orbital
        do ii=1,2 ! spin 

           k=iind((/ i,ii /))  ! spin orbital
           tempflag=.false.
           do j=1,www%numpart
              if (iind(thisconfig((j*2)-1:j*2)) == k) then  ! got this spin orbital in the configuration
                 tempflag=.true.
                 exit
              endif
           enddo
           if (tempflag) then
!counts how many spin orbitals of shell ishell are included in the configuration
              numinshell(ishell)=numinshell(ishell)+1   
!cumulative
              numwithinshell(ishell:numshells)=numwithinshell(ishell:numshells)+1
           endif
        enddo
     enddo
  enddo

  if (www%holeflag.ne.0) then
     numwithinshell(:)=2*allshelltop(1:numshells)-numwithinshell(:)
     numinshell(:)=2*(allshelltop(1:numshells)-allshelltop(0:numshells-1))-numinshell(:)
  endif
  
do ishell=1,numshells
     if (numwithinshell(ishell).lt.2*allshelltop(ishell)-numexcite(ishell)+in_df) then
        allowedconfig0=.false.
        return
     endif
     if (numinshell(ishell).gt.maxocc(ishell)-in_df) then
        allowedconfig0=.false.
        return
     endif
     if (numinshell(ishell).lt.minocc(ishell)+in_df) then
        allowedconfig0=.false.
        return
     endif
  enddo
  allowedconfig0=.true.

end function allowedconfig0


!! RETURNS M-VALUE OF CONFIGURATION

function getmval(www,thisconfig)
  use walkmod
  use basis_parameters
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: thisconfig(www%num2part)
  integer :: i, isum, getmval
  if ((spfrestrictflag==0)) then
     getmval=0;     return
  endif
  isum=0
  do i=1,www%num2part-1,2
     isum=isum+(spfmvals(thisconfig(i)))
  enddo
  getmval=isum
end function getmval

!! RETURNS UGVAL OF CONFIGURATION (GOES FROM ORIGINAL UG VALUES)
!!  IF UGVALS ARE ALL SPECIFIED (NONZERO)

function getugval(www,thisconfig)
  use walkmod
  use basis_parameters
  implicit none
  type(walktype),intent(in) :: www
  integer,intent(in) :: thisconfig(www%num2part)
  integer :: i, isum, getugval
  isum=1
  do i=1,www%num2part-1,2
     isum=isum*(spfugvals(thisconfig(i)))
  enddo
  getugval=isum

end function getugval


!!$  subroutine re_order_configlist(configlist,configorder,num2part,numconfig,numspinblocks,spinblockstart,spinblockend)
!!$    use fileptrmod
!!$    use sparse_parameters
!!$    implicit none
!!$    integer,intent(in) :: numconfig, num2part, numspinblocks
!!$    integer,intent(inout) :: spinblockstart(numspinblocks),spinblockend(numspinblocks),&
!!$         configlist(num2part,numconfig),configorder(numconfig)
!!$    integer :: ii,jj,kk,iconfig
!!$    integer,allocatable :: newstart(:),newend(:),newlist(:,:),neworder(:)
!!$  
!!$    allocate(newstart(numspinblocks),newend(numspinblocks),newlist(num2part,numconfig), neworder(numconfig))
!!$    newstart=0; newend=0; newlist=0; neworder=0
!!$  
!!$    iconfig=0
!!$    do ii=1,numspinblocks
!!$       jj=mod((ii-1)*sparseprime,numspinblocks)+1
!!$       newstart(ii)=iconfig+1
!!$       newend(ii)=newstart(ii) + (spinblockend(jj)-spinblockstart(jj))
!!$  
!!$  
!!$       newlist(:,newstart(ii):newend(ii)) = configlist(:,spinblockstart(jj):spinblockend(jj))
!!$  
!!$  !!     iconfig=iconfig+(spinblockend(jj)-spinblockstart(jj))+1
!!$  
!!$       do kk=spinblockstart(jj),spinblockend(jj)
!!$          iconfig=iconfig+1
!!$          neworder(kk)=configorder(iconfig)
!!$       enddo
!!$    enddo
!!$    if (iconfig.ne.numconfig) then
!!$       OFLWR "CONFIGCHEK ERROR ",iconfig,numconfig; CFLST
!!$    endif
!!$  
!!$    configlist(:,:)=newlist(:,:)
!!$    spinblockstart(:)=newstart(:)
!!$    spinblockend(:)=newend(:)
!!$    configorder(:)=neworder(:)
!!$  
!!$    deallocate(newstart,newend,newlist,neworder)
!!$  
!!$  end subroutine re_order_configlist



!! GETS CONFIGURATION LIST (SLATER DETERMINANTS NOT SPIN EIGFUNCTS)
!!  AT BEGINNING. 

subroutine fast_newconfiglist(www,domflags,dorestrictflags)
  use output_parameters
  use fileptrmod
  use basis_parameters
  use sparse_parameters
  use walkmod
  use mpimod
  use ham_parameters  !! tdflag and offaxispulseflag for configtypes
  use mpisubmod
  implicit none
  type(walktype) :: www
  logical,intent(in) :: domflags,dorestrictflags
  logical :: alreadycounted
  integer, parameter :: max_numpart=80
  integer, pointer :: &
       ii1,ii2,ii3,ii4,ii5,ii6,ii7,ii8,ii9,ii10, &
       ii11,ii12, ii13, ii14, ii15, ii16, ii17, ii18, ii19, ii20, &
       ii21,ii22, ii23, ii24, ii25, ii26, ii27, ii28, ii29, ii30, &
       ii31,ii32, ii33, ii34, ii35, ii36, ii37, ii38, ii39, ii40, &
       ii41,ii42, ii43, ii44, ii45, ii46, ii47, ii48, ii49, ii50, &
       ii51,ii52, ii53, ii54, ii55, ii56, ii57, ii58, ii59, ii60, &
       ii61,ii62, ii63, ii64, ii65, ii66, ii67, ii68, ii69, ii70, &
       ii71,ii72, ii73, ii74, ii75, ii76, ii77, ii78, ii79, ii80
  integer, pointer :: &
       jj1,jj2,jj3,jj4,jj5,jj6,jj7,jj8,jj9,jj10, &
       jj11,jj12, jj13, jj14, jj15, jj16, jj17, jj18, jj19, jj20, &
       jj21,jj22, jj23, jj24, jj25, jj26, jj27, jj28, jj29, jj30, &
       jj31,jj32, jj33, jj34, jj35, jj36, jj37, jj38, jj39, jj40, &
       jj41,jj42, jj43, jj44, jj45, jj46, jj47, jj48, jj49, jj50, &
       jj51,jj52, jj53, jj54, jj55, jj56, jj57, jj58, jj59, jj60, &
       jj61,jj62, jj63, jj64, jj65, jj66, jj67, jj68, jj69, jj70, &
       jj71,jj72, jj73, jj74, jj75, jj76, jj77, jj78, jj79, jj80
  integer, target :: iii(max_numpart)  !! no, it is set.
  integer, target :: jjj(max_numpart)  !! no, it is set.

  integer,allocatable :: bigspinblockstart(:),bigspinblockend(:)
  integer :: i, idof, ii , lowerarr(max_numpart),upperarr(max_numpart),  thisconfig(www%num2part),&
       nullint,kk,iconfig,mm, single(max_numpart),ishell,jj,maxssize=0,sss,nss,ssflag,numdoubly,&
       mynumexcite(numshells),isums(numshells),jshell,jsums(numshells),&
       ppp,numspinblocks,mytot(numshells)

  if (www%numpart.gt.max_numpart) then
     OFLWR "Resize get_newconfiglist"; CFLST
  endif

  if (dorestrictflags) then
     www%m_restrictflag=mrestrictflag
     www%m_restrictval=mrestrictval
     www%m_restrictmin=mrestrictmin
     www%m_restrictmax=mrestrictmax
     www%ug_restrictflag=ugrestrictflag
     www%ug_restrictval=ugrestrictval
  endif

!!$ orderflag 0 hardwire
!!$  if (orderflag.eq.1) then
!!$     OFLWR "orderflag 1 not supported for fastconfig (would be trivial, a simplification, no sort I think)"; CFLST
!!$  endif

  allocate(www%configsperproc(nprocs),www%alltopconfigs(nprocs),www%allbotconfigs(nprocs))
  www%configsperproc=0; www%alltopconfigs=0; www%allbotconfigs=0

!! 10-2015 NOW INTERNAL LOOP

  do ppp=0,1
     if (ppp.eq.0) then
        alreadycounted=.false.
     else
        alreadycounted=.true.
     endif

  if (alreadycounted) then
     OFLWR "Go fast_newconfiglist.  Allocating...";CFL
     call waitawhile()
     call mpibarrier()
     deallocate(bigspinblockstart,bigspinblockend)
     allocate(www%configlist(www%num2part,www%numconfig),  www%configtypes(www%numconfig), &
          bigspinblockstart(numspinblocks+2*nprocs),bigspinblockend(numspinblocks+2*nprocs))
     bigspinblockstart=0; bigspinblockend=0; 
     www%configlist(:,:)=0;  www%configtypes(:)=0
     call waitawhile()
     call mpibarrier()
     OFLWR "   Allocated.  getting configurations."; CFL
  else
     allocate(bigspinblockstart(1),bigspinblockend(1))  !! avoid warn bounds
     OFLWR "Go fast_newconfiglist";CFL
  endif

  ii1 => iii(1) ; ii2 => iii(2) ; ii3 => iii(3) ; ii4 => iii(4) ; ii5 => iii(5) ; ii6 => iii(6) ;
  ii7 => iii(7) ; ii8 => iii(8) ; ii9 => iii(9) ; ii10 => iii(10) ; ii11 => iii(11) ; ii12 => iii(12) ;
  ii13 => iii(13) ; ii14 => iii(14) ; ii15 => iii(15) ; ii16 => iii(16) ; ii17 => iii(17); ii18 => iii(18)
  ii19 => iii(19);ii20 => iii(20)

  ii21 => iii(21) ; ii22 => iii(22) ; ii23 => iii(23) ; ii24 => iii(24) ; ii25 => iii(25) ; ii26 => iii(26)
  ii27 => iii(27) ; ii28 => iii(28) ; ii29 => iii(29) ; ii30 => iii(30) ; ii31 => iii(31) ; ii32 => iii(32)
  ii33 => iii(33) ; ii34 => iii(34) ; ii35 => iii(35) ; ii36 => iii(36) ; ii37 => iii(37) ; ii38 => iii(38)
  ii39 => iii(39);ii40 => iii(40)

  ii41 => iii(41) ; ii42 => iii(42) ; ii43 => iii(43) ; ii44 => iii(44) ; ii45 => iii(45) ; ii46 => iii(46)
  ii47 => iii(47) ; ii48 => iii(48) ; ii49 => iii(49) ; ii50 => iii(50) ; ii51 => iii(51) ; ii52 => iii(52)
  ii53 => iii(53) ; ii54 => iii(54) ; ii55 => iii(55) ; ii56 => iii(56) ; ii57 => iii(57) ; ii58 => iii(58)
  ii59 => iii(59);ii60 => iii(60)

  ii61 => iii(61) ; ii62 => iii(62) ; ii63 => iii(63) ; ii64 => iii(64) ; ii65 => iii(65) ; ii66 => iii(66)
  ii67 => iii(67) ; ii68 => iii(68) ; ii69 => iii(69) ; ii70 => iii(70) ; ii71 => iii(71) ; ii72 => iii(72)
  ii73 => iii(73) ; ii74 => iii(74) ; ii75 => iii(75) ; ii76 => iii(76) ; ii77 => iii(77) ; ii78 => iii(78)
  ii79 => iii(79);  ii80 => iii(80)


  jj1 => jjj(1) ; jj2 => jjj(2) ; jj3 => jjj(3) ; jj4 => jjj(4) ; jj5 => jjj(5) ; jj6 => jjj(6) ;
  jj7 => jjj(7) ; jj8 => jjj(8) ; jj9 => jjj(9) ; jj10 => jjj(10) ; jj11 => jjj(11) ; jj12 => jjj(12) ;
  jj13 => jjj(13) ; jj14 => jjj(14) ; jj15 => jjj(15) ; jj16 => jjj(16) ; jj17 => jjj(17); jj18 => jjj(18)
  jj19 => jjj(19);jj20 => jjj(20)

  jj21 => jjj(21) ; jj22 => jjj(22) ; jj23 => jjj(23) ; jj24 => jjj(24) ; jj25 => jjj(25) ; jj26 => jjj(26)
  jj27 => jjj(27) ; jj28 => jjj(28) ; jj29 => jjj(29) ; jj30 => jjj(30) ; jj31 => jjj(31) ; jj32 => jjj(32)
  jj33 => jjj(33) ; jj34 => jjj(34) ; jj35 => jjj(35) ; jj36 => jjj(36) ; jj37 => jjj(37) ; jj38 => jjj(38)
  jj39 => jjj(39);jj40 => jjj(40)

  jj41 => jjj(41) ; jj42 => jjj(42) ; jj43 => jjj(43) ; jj44 => jjj(44) ; jj45 => jjj(45) ; jj46 => jjj(46)
  jj47 => jjj(47) ; jj48 => jjj(48) ; jj49 => jjj(49) ; jj50 => jjj(50) ; jj51 => jjj(51) ; jj52 => jjj(52)
  jj53 => jjj(53) ; jj54 => jjj(54) ; jj55 => jjj(55) ; jj56 => jjj(56) ; jj57 => jjj(57) ; jj58 => jjj(58)
  jj59 => jjj(59);jj60 => jjj(60)

  jj61 => jjj(61) ; jj62 => jjj(62) ; jj63 => jjj(63) ; jj64 => jjj(64) ; jj65 => jjj(65) ; jj66 => jjj(66)
  jj67 => jjj(67) ; jj68 => jjj(68) ; jj69 => jjj(69) ; jj70 => jjj(70) ; jj71 => jjj(71) ; jj72 => jjj(72)
  jj73 => jjj(73) ; jj74 => jjj(74) ; jj75 => jjj(75) ; jj76 => jjj(76) ; jj77 => jjj(77) ; jj78 => jjj(78)
  jj79 => jjj(79);  jj80 => jjj(80)


  lowerarr(:)=1000000
  upperarr(:)=1000000

  numdoubly=www%numpart-abs(www%restrictms)

  if (mod(numdoubly,2).ne.0) then
     OFLWR "Looks like restrictms not compatible with numelec - need both even or both odd ",&
          www%restrictms,www%numelec,www%numpart; CFLST
  endif

  do ii=1,numdoubly
     lowerarr(ii)=(ii+1)/2
  enddo
  do ii=1,abs(www%restrictms)
     lowerarr(ii+numdoubly)=numdoubly/2+ii
  enddo

  do ii=1,www%numpart
     upperarr(ii)=www%nspf+1-lowerarr(www%numpart+1-ii)
  enddo

  do ii=www%numpart+3,max_numpart
     upperarr(ii)=1000000+ii-www%numpart+2
     lowerarr(ii)=1000000+ii-www%numpart+2
  enddo

!!#          min        max
!!numpart     ..        nspf
!!numpart+1  1000000    1000000
!!numpart+2  1000000    1000000
!!numpart+3  1000001    1000001  etc.
!!numpart+4  1000002    1000002
!!numpart+4  1000003    1000003

!! numexcite(:) is CUMULATIVE

  mytot(1:numshells)=2*(allshelltop(1:numshells)-allshelltop(0:numshells-1))

!! isums: cumulative maximum excitations given maxocc for higher shells

  isums(:)=0
  do ishell=1,numshells

!! count maximum number of excitations in ishell given maxocc for shells greater than ishell
     do jshell=ishell+1,numshells
        isums(ishell) = isums(ishell) + max(0,min(maxocc(jshell),mytot(jshell)))
     enddo
  enddo

!!$ prev      isums(1:numshells)=max(0,2*allshelltop(1:numshells)-max(0,www%numelec-isums(1:numshells)))

!! minimum number of electrons in shells 1 through ishell given maxocc for higher shells
!!    is then max(0,numelec-isum)
!! maximum number of holes in shells 1 through ishell given maxocc for higher shells is then 
!!     2*allshelltop(ishell)-max(0,numelec-isum) 

  isums(1:numshells)=max(0,2*allshelltop(1:numshells)-max(0,www%numelec-isums(1:numshells)))

!! --> isums(ishell) is now the maximum number of holes in shell 1 through ishell <--
!! --> given maximum electron occupancy restriction on other shells <--

!! count maximum number of excitations in shells 1 through ishell, given minocc for those shells

!! min number of electrons
  jsums=0
  do ishell=1,numshells
     do jshell=1,ishell
        jsums(ishell)=jsums(ishell)+max(0,min(minocc(jshell),mytot(jshell)))
     enddo
  enddo

!! jsums is now minimum number of electrons

!! maximum number of holes in shells 1 through ishell, given minocc for those shells, is 
!!    2*allshelltop(ishell)-jsums(ishell)
!! set jsums to this:

  jsums(:)=max(0,2*allshelltop(1:numshells)-jsums(:))

!! --> jsums(ishell) is now the maximum number of holes in shell 1 through ishell <--
!! --> given minimum occupancy restriction on shells 1 through ishell <--

!! now the final maximum number of holes
!!   for shells 1 through ishell is mynumexcite(ishell) accounting both for
!!   minocc restriction on shells 1 through ishell and for maxocc for higher shells is mynumexcite

  mynumexcite(1:numshells)=min(min(numexcite(1:numshells),isums(:)), jsums(:))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  OFLWR "Excitations per shell"
  write(mpifileptr,'(3I10)') (mynumexcite(ishell),isums(ishell),jsums(ishell),ishell=1,numshells)
  CFL

  kk=0  !! electron

  do ishell=1,numshells

!! maximum number of cumulative excitations mm

     if (www%holeflag.eq.0) then

!! mm is minimum number of electrons in shells 1 through ishell
!! kk is the old minimum

        mm=2*allshelltop(ishell)-min(mynumexcite(ishell),2*www%nspf-www%numpart)

        do jj=kk+1,mm   !! which electron. the jjth electron can start in the qth orbital,
                           !! which is in ishell, where q = allshelltop(ishell)-(mm-jj)/2.  Right?

           upperarr(jj)=min(upperarr(jj),allshelltop(ishell)-(mm-jj)/2)
        enddo

     else
        
!! mm is minimum number of holes in shells 1 through ishell
!! kk is the old minimum

        mm=mynumexcite(ishell)

        do jj=kk+1,mm   !! which electron. the jjth electron can start in the qth orbital,
                           !! which is in ishell, where q = allshelltop(ishell)-(mm-jj)/2.  Right?

           lowerarr(jj)=max(lowerarr(jj),allshelltop(ishell-1)-(mm-jj)/2)
        enddo

     endif

     kk=max(0,mm)

  enddo  ! ishell


  if (.not.alreadycounted) then
     if (www%holeflag.eq.0) then
        OFLWR "Orbital ranges for each electron ::"
     else
        OFLWR "Orbital ranges for each hole ::"
     endif
     do ii=1,www%numpart
        WRFL lowerarr(ii),upperarr(ii)
     enddo
     WRFL; CFL
  endif

  iconfig=0
  nss=0
  maxssize=0
  ssflag=1

  do ii1=  lowerarr(1),upperarr(1)
  do ii2=      max(lowerarr(2), ii1)         ,upperarr(2)
  do ii3=  max(max(lowerarr(3) ,ii2),ii1 +1 ),upperarr(3)
  do ii4=  max(max(lowerarr(4) ,ii3 ),ii2 +1 ),upperarr(4)
  do ii5=  max(max(lowerarr(5) ,ii4 ),ii3 +1 ),upperarr(5)
  do ii6=  max(max(lowerarr(6) ,ii5 ),ii4 +1 ),upperarr(6)
  do ii7=  max(max(lowerarr(7) ,ii6 ),ii5 +1 ),upperarr(7)
  do ii8=  max(max(lowerarr(8) ,ii7 ),ii6 +1 ),upperarr(8)
  do ii9=  max(max(lowerarr(9) ,ii8 ),ii7 +1 ),upperarr(9)

  do ii10= max(max(lowerarr(10),ii9 ),ii8 +1 ),upperarr(10)
  do ii11= max(max(lowerarr(11),ii10),ii9 +1 ),upperarr(11)
  do ii12= max(max(lowerarr(12),ii11),ii10+1 ),upperarr(12)
  do ii13= max(max(lowerarr(13),ii12),ii11+1 ),upperarr(13)
  do ii14= max(max(lowerarr(14),ii13),ii12+1 ),upperarr(14)
  do ii15= max(max(lowerarr(15),ii14),ii13+1 ),upperarr(15)
  do ii16= max(max(lowerarr(16),ii15),ii14+1 ),upperarr(16)
  do ii17= max(max(lowerarr(17),ii16),ii15+1 ),upperarr(17)
  do ii18= max(max(lowerarr(18),ii17),ii16+1 ),upperarr(18)
  do ii19= max(max(lowerarr(19),ii18),ii17+1 ),upperarr(19)

  do ii20= max(max(lowerarr(20),ii19),ii18+1),upperarr(20)
  do ii21= max(max(lowerarr(21),ii20),ii19+1),upperarr(21)
  do ii22= max(max(lowerarr(22),ii21),ii20+1),upperarr(22)
  do ii23= max(max(lowerarr(23),ii22),ii21+1),upperarr(23)
  do ii24= max(max(lowerarr(24),ii23),ii22+1),upperarr(24)
  do ii25= max(max(lowerarr(25),ii24),ii23+1),upperarr(25)
  do ii26= max(max(lowerarr(26),ii25),ii24+1),upperarr(26)
  do ii27= max(max(lowerarr(27),ii26),ii25+1),upperarr(27)
  do ii28= max(max(lowerarr(28),ii27),ii26+1),upperarr(28)
  do ii29= max(max(lowerarr(29),ii28),ii27+1),upperarr(29)

  do ii30= max(max(lowerarr(30),ii29),ii28+1),upperarr(30)
  do ii31= max(max(lowerarr(31),ii30),ii29+1),upperarr(31)
  do ii32= max(max(lowerarr(32),ii31),ii30+1),upperarr(32)
  do ii33= max(max(lowerarr(33),ii32),ii31+1),upperarr(33)
  do ii34= max(max(lowerarr(34),ii33),ii32+1),upperarr(34)
  do ii35= max(max(lowerarr(35),ii34),ii33+1),upperarr(35)
  do ii36= max(max(lowerarr(36),ii35),ii34+1),upperarr(36)
  do ii37= max(max(lowerarr(37),ii36),ii35+1),upperarr(37)
  do ii38= max(max(lowerarr(38),ii37),ii36+1),upperarr(38)
  do ii39= max(max(lowerarr(39),ii38),ii37+1),upperarr(39)

  do ii40= max(max(lowerarr(40),ii39),ii38+1),upperarr(40)
  do ii41= max(max(lowerarr(41),ii40),ii39+1),upperarr(41)
  do ii42= max(max(lowerarr(42),ii41),ii40+1),upperarr(42)
  do ii43= max(max(lowerarr(43),ii42),ii41+1),upperarr(43)
  do ii44= max(max(lowerarr(44),ii43),ii42+1),upperarr(44)
  do ii45= max(max(lowerarr(45),ii44),ii43+1),upperarr(45)
  do ii46= max(max(lowerarr(46),ii45),ii44+1),upperarr(46)
  do ii47= max(max(lowerarr(47),ii46),ii45+1),upperarr(47)
  do ii48= max(max(lowerarr(48),ii47),ii46+1),upperarr(48)
  do ii49= max(max(lowerarr(49),ii48),ii47+1),upperarr(49)

  do ii50= max(max(lowerarr(50),ii49),ii48+1),upperarr(50)
  do ii51= max(max(lowerarr(51),ii50),ii49+1),upperarr(51)
  do ii52= max(max(lowerarr(52),ii51),ii50+1),upperarr(52)
  do ii53= max(max(lowerarr(53),ii52),ii51+1),upperarr(53)
  do ii54= max(max(lowerarr(54),ii53),ii52+1),upperarr(54)
  do ii55= max(max(lowerarr(55),ii54),ii53+1),upperarr(55)
  do ii56= max(max(lowerarr(56),ii55),ii54+1),upperarr(56)
  do ii57= max(max(lowerarr(57),ii56),ii55+1),upperarr(57)
  do ii58= max(max(lowerarr(58),ii57),ii56+1),upperarr(58)
  do ii59= max(max(lowerarr(59),ii58),ii57+1),upperarr(59)

  do ii60= max(max(lowerarr(60),ii59),ii58+1),upperarr(60)
  do ii61= max(max(lowerarr(61),ii60),ii59+1),upperarr(61)
  do ii62= max(max(lowerarr(62),ii61),ii60+1),upperarr(62)
  do ii63= max(max(lowerarr(63),ii62),ii61+1),upperarr(63)
  do ii64= max(max(lowerarr(64),ii63),ii62+1),upperarr(64)
  do ii65= max(max(lowerarr(65),ii64),ii63+1),upperarr(65)
  do ii66= max(max(lowerarr(66),ii65),ii64+1),upperarr(66)
  do ii67= max(max(lowerarr(67),ii66),ii65+1),upperarr(67)
  do ii68= max(max(lowerarr(68),ii67),ii66+1),upperarr(68)
  do ii69= max(max(lowerarr(69),ii68),ii67+1),upperarr(69)

  do ii70= max(max(lowerarr(70),ii69),ii68+1),upperarr(70)
  do ii71= max(max(lowerarr(71),ii70),ii69+1),upperarr(71)
  do ii72= max(max(lowerarr(72),ii71),ii70+1),upperarr(72)
  do ii73= max(max(lowerarr(73),ii72),ii71+1),upperarr(73)
  do ii74= max(max(lowerarr(74),ii73),ii72+1),upperarr(74)
  do ii75= max(max(lowerarr(75),ii74),ii73+1),upperarr(75)
  do ii76= max(max(lowerarr(76),ii75),ii74+1),upperarr(76)
  do ii77= max(max(lowerarr(77),ii76),ii75+1),upperarr(77)
  do ii78= max(max(lowerarr(78),ii77),ii76+1),upperarr(78)
  do ii79= max(max(lowerarr(79),ii78),ii77+1),upperarr(79)

  do ii80= max(max(lowerarr(80),ii79),ii78+1),upperarr(80)

  if (okexcite(iii)) then
     sss=0
     ssflag=0

     do idof=1,www%numpart
        thisconfig(idof*2-1) = iii(idof)
     enddo

     single(:)=1
     single(www%numpart+1:)=0

     do idof=2,www%numpart
        if (iii(idof).eq.iii(idof-1)) then
           thisconfig((idof-1)*2)=1
           thisconfig(idof*2)=2
           single(idof-1)=0
           single(idof)=0
        endif
     enddo

     do jj1=0,single(1)
     do jj2=0,single(2)
     do jj3=0,single(3)
     do jj4=0,single(4)
     do jj5=0,single(5)
     do jj6=0,single(6)
     do jj7=0,single(7)
     do jj8=0,single(8)
     do jj9=0,single(9)
     do jj10=0,single(10)

     do jj11=0,single(11)
     do jj12=0,single(12)
     do jj13=0,single(13)
     do jj14=0,single(14)
     do jj15=0,single(15)
     do jj16=0,single(16)
     do jj17=0,single(17)
     do jj18=0,single(18)
     do jj19=0,single(19)
     do jj20=0,single(20)

     do jj21=0,single(21)
     do jj22=0,single(22)
     do jj23=0,single(23)
     do jj24=0,single(24)
     do jj25=0,single(25)
     do jj26=0,single(26)
     do jj27=0,single(27)
     do jj28=0,single(28)
     do jj29=0,single(29)
     do jj30=0,single(30)

     do jj31=0,single(31)
     do jj32=0,single(32)
     do jj33=0,single(33)
     do jj34=0,single(34)
     do jj35=0,single(35)
     do jj36=0,single(36)
     do jj37=0,single(37)
     do jj38=0,single(38)
     do jj39=0,single(39)
     do jj40=0,single(40)

     do jj41=0,single(41)
     do jj42=0,single(42)
     do jj43=0,single(43)
     do jj44=0,single(44)
     do jj45=0,single(45)
     do jj46=0,single(46)
     do jj47=0,single(47)
     do jj48=0,single(48)
     do jj49=0,single(49)
     do jj50=0,single(50)

     do jj51=0,single(51)
     do jj52=0,single(52)
     do jj53=0,single(53)
     do jj54=0,single(54)
     do jj55=0,single(55)
     do jj56=0,single(56)
     do jj57=0,single(57)
     do jj58=0,single(58)
     do jj59=0,single(59)
     do jj60=0,single(60)

     do jj61=0,single(61)
     do jj62=0,single(62)
     do jj63=0,single(63)
     do jj64=0,single(64)
     do jj65=0,single(65)
     do jj66=0,single(66)
     do jj67=0,single(67)
     do jj68=0,single(68)
     do jj69=0,single(69)
     do jj70=0,single(70)

     do jj71=0,single(71)
     do jj72=0,single(72)
     do jj73=0,single(73)
     do jj74=0,single(74)
     do jj75=0,single(75)
     do jj76=0,single(76)
     do jj77=0,single(77)
     do jj78=0,single(78)
     do jj79=0,single(79)
     do jj80=0,single(80)

        do idof=1,www%numpart
           if (single(idof).eq.1) then
              thisconfig(idof*2)=jjj(idof)+1
           endif
        enddo

        if (allowedconfig0(www,thisconfig,www%dflevel)) then
           sss=sss+1
           iconfig=iconfig+1
           if (ssflag.eq.0) then 
              nss=nss+1
              ssflag=1
              if (alreadycounted) then 
                 bigspinblockstart(nss)=iconfig
              endif
           endif
           if (alreadycounted) then
              bigspinblockend(nss)=iconfig
              www%configlist(:,iconfig)=thisconfig; 
              nullint=reorder(www%configlist(:,iconfig),www%numpart)
           endif
        endif


  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo
  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo
  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo
  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo
  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo
  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo
  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo
  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo

  if (sss.gt.maxssize) then
     maxssize=sss
  endif
endif

  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo
  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo
  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo
  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo
  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo
  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo
  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo
  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo;  enddo

  if (.not.alreadycounted) then
     www%numconfig=iconfig
     OFLWR; WRFL "FASTNEWCONFIG: NUMBER OF CONFIGURATIONS ",www%numconfig; WRFL; CFL
     if (www%numconfig.le.0) then
        OFLWR "NO configs."; CFLST
     endif

     numspinblocks=nss
     OFLWR "NUMSPINBLOCKS, MAXSPINBLOCKSIZE FASTCONFIG",nss,maxssize;CFL
  endif

  if (www%numconfig.eq.0) then
     OFLWR "No configs!! "; CFLST
  endif

  enddo !! ppp


  if (iconfig/=www%numconfig) then
     OFLWR "Configlist err",iconfig,www%numconfig; CFLST
  endif

!!$SP  do iconfig=1,www%numconfig
!!$SP     www%configorder(iconfig)=iconfig
!!$SP  enddo

  if (sparseprime.ne.1) then
     OFLWR "Reordering configlist with sparseprime is NO LONGER SUPPORTED"
     CFLST
  endif

!!$SP       OFLWR "Reordering configlist with sparseprime=",sparseprime; CFL
!!$SP  
!!$SP       call re_order_configlist(www%configlist(:,:),www%configorder(:),&
!!$SP             www%num2part,www%numconfig,numspinblocks,bigspinblockstart(:),&
!!$SP             bigspinblockend(:))
!!$SP    endif

  if (nss.ne.numspinblocks) then
     OFLWR "NUMSPINBLOCKS ERR",nss,numspinblocks; CFLST
  endif

  www%alltopconfigs(:)=0
  jj=1
  do ii=1,nprocs-1

!!     do while (bigspinblockend(jj).lt.www%numconfig*ii/nprocs)

     do while (bigspinblockend(jj).lt.floor(1d0*www%numconfig*ii/nprocs))
        www%alltopconfigs(ii:)=bigspinblockend(jj)
        jj=jj+1
     enddo
  enddo
  www%alltopconfigs(nprocs)=www%numconfig
  www%allbotconfigs(1)=1
  do ii=2,nprocs
     www%allbotconfigs(ii)=www%alltopconfigs(ii-1)+1
  enddo

  OFLWR; WRFL "BOTWALKS /TOPWALKS",www%numconfig
  do ii=1,nprocs
     WRFL www%allbotconfigs(ii),www%alltopconfigs(ii),www%alltopconfigs(ii)-www%allbotconfigs(ii)+1
  enddo
  WRFL; CFL
  do ii=1,nprocs     
     if (www%allbotconfigs(ii).gt.www%alltopconfigs(ii)+1) then   
        OFLWR "ERROR, NUMBER OF CONFIGS PROCESSOR ",ii," IS LESS THAN ZERO", www%alltopconfigs(ii)-www%allbotconfigs(ii)+1 ;CFLST
     endif
  enddo
  
  www%botconfig=www%allbotconfigs(myrank)
  www%topconfig=www%alltopconfigs(myrank)

  www%sparseconfigflag=sparseconfigflag

  if (sparseconfigflag.eq.0) then
     www%configstart=1
     www%configend=www%numconfig
     www%parconsplit=0
     www%startrank=1
     www%endrank=nprocs
  else
     www%configstart=www%botconfig
     www%configend=www%topconfig
     www%startrank=myrank
     www%endrank=myrank
  endif

  www%configsperproc(:)=www%alltopconfigs(:) - www%allbotconfigs(:) + 1

  if (www%parconsplit.eq.0) then
     www%firstconfig=1
     www%lastconfig=www%numconfig
  else
     www%firstconfig=www%botconfig
     www%lastconfig=www%topconfig
  endif

  allocate(          www%configmvals(www%firstconfig:www%lastconfig), &
       www%configugvals(www%firstconfig:www%lastconfig))
  www%configmvals(:)=0; www%configugvals(:)=0; 

!!$SP  allocate(www%configorder(www%firstconfig:www%lastconfig)); www%configorder=0

  www%localnconfig=www%lastconfig-www%firstconfig+1

  ii=0
  www%maxconfigsperproc=0
  do i=1,nprocs
     ii=ii+www%configsperproc(i)
     if (www%configsperproc(i).gt.www%maxconfigsperproc) then
        www%maxconfigsperproc=www%configsperproc(i)
     endif
  enddo
  if (ii.ne.www%numconfig) then
     OFLWR "EERRROROR", ii, www%numconfig; CFLST
  endif
  do i=1,nprocs
     if (www%configsperproc(i).lt.0) then
        OFLWR "Configs per proc lt 0 for proc ",i,"nprocs must be greater than numconfig???  Can't do."; CFLST
     endif
  enddo

  if (iprintconfiglist.ne.0) then
     OFLWR "CONFIGLIST"
     do ii=1,www%numconfig
!        write(mpifileptr,'(A12,I12,A4)',advance='no') "  Config ", ii," is "
        write(mpifileptr,'(A12,I12,A4)',advance='no') "  Config ", 0," is "
        call printconfig(www%configlist(:,ii),www)
     enddo
     WRFL; CFL
  endif

  deallocate(bigspinblockstart,bigspinblockend)

  if (spfrestrictflag.ne.0) then
     do ii=www%firstconfig,www%lastconfig
        www%configmvals(ii)=getmval(www,www%configlist(:,ii))
     enddo
     if (spfugrestrict.ne.0) then
        do ii=www%firstconfig,www%lastconfig
           www%configugvals(ii)=getugval(www,www%configlist(:,ii))
        enddo
     endif
  endif

  if (domflags.and.(tdflag.eq.0.or.offaxispulseflag.eq.0)) then
     if (spfrestrictflag.ne.0) then
        if (spfugrestrict.ne.0.and.tdflag.eq.0) then
           www%configtypes(www%firstconfig:www%lastconfig) = &
                www%configmvals(:)*www%configugvals(:)
        else
           www%configtypes(www%firstconfig:www%lastconfig)=www%configmvals(:)
        endif
     endif
#ifdef MPIFLAG
     if (www%parconsplit.ne.0) then
        call mpiallgather_i(www%configtypes,www%numconfig,www%configsperproc(:),&
             www%maxconfigsperproc)
     endif
#endif
  endif

  OFLWR "     ...Done fast_newconfiglist"; CFL
  return
  return
  return

contains
  function okexcite(kkk)
    implicit none
    logical :: okexcite
    integer,intent(in) :: kkk(www%numpart)
    integer :: numwithinshell(numshells),ii,ishell,numinshell(numshells)

    numwithinshell(:)=0; numinshell(:)=0
    do ii=1,www%numpart
       do ishell=1,numshells
          if (kkk(ii).gt.allshelltop(ishell-1).and.kkk(ii).le.allshelltop(ishell)) then
             numinshell(ishell)=numinshell(ishell)+1
             numwithinshell(ishell:numshells)=numwithinshell(ishell:numshells)+1
             exit
          endif
       enddo
    enddo

    if (www%holeflag.ne.0) then
       numwithinshell(:)=2*allshelltop(1:numshells)-numwithinshell(:)
       numinshell(:)=2*(allshelltop(1:numshells)-allshelltop(0:numshells-1))-numinshell(:)
    endif

    do ishell=1,numshells
       if (numwithinshell(ishell).lt.2*allshelltop(ishell)-numexcite(ishell)) then
          okexcite=.false.
          return
       endif
       if ( numinshell(ishell).lt.minocc(ishell) .or. numinshell(ishell).gt.maxocc(ishell)) then
          okexcite=.false.
          return
       endif
    enddo

    okexcite=.true.
  end function okexcite

  function quickaarr(xind)
    implicit none
    integer, dimension(2) :: quickaarr
    integer, save :: temp(2)
    integer,intent(in) :: xind
    integer :: ind,q

    ind=xind-1
!orderflag 0
    temp(2)           = mod(ind,2)+1
!!       q=(ind-mod(ind,2))/2  no need
    q=ind/2                   !! rounds down.
    temp(1)           = q+1
    quickaarr=temp
  end function quickaarr

end subroutine fast_newconfiglist




subroutine set_newconfiglist(wwin,wwout,domflags)
  use fileptrmod
  use walkmod
  use mpimod
  use basis_parameters
  use ham_parameters  !! tdflag and offaxispulseflag for configtypes
  use mpisubmod
  implicit none
  type(walktype),intent(in) :: wwin
  logical,intent(in) :: domflags
  type(walktype),intent(inout) :: wwout
  integer :: iconfig,jconfig

  wwout%m_restrictflag=wwin%m_restrictflag
  wwout%m_restrictval=wwin%m_restrictval
  wwout%m_restrictmin=wwin%m_restrictmin
  wwout%m_restrictmax=wwin%m_restrictmax
  wwout%ug_restrictflag=wwin%ug_restrictflag
  wwout%ug_restrictval=wwin%ug_restrictval

  wwout%dflevel=wwin%dfrestrictflag

  allocate(wwout%configsperproc(nprocs),wwout%alltopconfigs(nprocs),wwout%allbotconfigs(nprocs))
  wwout%configsperproc(:)=wwin%dfconfsperproc(:)
  wwout%allbotconfigs(:)=wwin%allbotdfconfigs(:)
  wwout%alltopconfigs(:)=wwin%alltopdfconfigs(:)

  wwout%numconfig=wwin%numdfconfigs

  wwout%maxconfigsperproc=wwin%maxdfconfsperproc

  allocate(wwout%configlist(wwout%num2part,wwout%numconfig),wwout%configtypes(wwout%numconfig))
  wwout%configlist=0; wwout%configtypes=0;

  do iconfig=1,wwout%numconfig
     wwout%configlist(:,iconfig)=wwin%configlist(:,wwin%ddd%dfincludedconfigs(iconfig))
  enddo

  jconfig=0
  do iconfig=1,wwin%numconfig
!!$SP     if (wwin%ddd%dfincludedmask(wwin%configorder(iconfig)).ne.0) then
!!$SP        jconfig=jconfig+1
!!$SP        wwout%configorder(jconfig)=wwin%ddd%dfincludedindex(wwin%configorder(iconfig))
!!$SP     endif
     if (wwin%ddd%dfincludedmask(iconfig).ne.0) then
        jconfig=jconfig+1
     endif
  enddo
  if (jconfig.ne.wwout%numconfig) then
     OFLWR "CHECKFAIL JCONFIG", jconfig,wwout%numconfig; CFLST
  endif

  wwout%botconfig=wwin%botdfconfig
  wwout%topconfig=wwin%topdfconfig
  
  wwout%sparseconfigflag=wwin%sparseconfigflag

  if (wwout%sparseconfigflag.eq.0) then
     wwout%configstart=1
     wwout%configend=wwin%numdfconfigs
  else
     wwout%configstart=wwin%botdfconfig
     wwout%configend=wwin%topdfconfig
  endif

  wwout%startrank=wwin%startrank
  wwout%endrank=wwin%endrank

  wwout%parconsplit=wwin%parconsplit

  if (wwout%parconsplit.eq.0) then
     wwout%firstconfig=1
     wwout%lastconfig=wwout%numconfig
  else
     wwout%firstconfig=wwout%botconfig
     wwout%lastconfig=wwout%topconfig
  endif

  wwout%localnconfig=(wwout%lastconfig-wwout%firstconfig+1)

  allocate(  wwout%configmvals(wwout%firstconfig:wwout%lastconfig),&
       wwout%configugvals(wwout%firstconfig:wwout%lastconfig))
  wwout%configmvals=0; wwout%configugvals=0; 
  
!!$SP  allocate(      wwout%configorder(wwout%firstconfig:wwout%lastconfig));   
!!$SP  wwout%configorder=0

  do iconfig=wwout%firstconfig,wwout%lastconfig
     wwout%configmvals(iconfig)=getmval(wwout,wwout%configlist(:,iconfig))
     wwout%configugvals(iconfig)=getugval(wwout,wwout%configlist(:,iconfig))
  enddo
  if (domflags.and.(tdflag.eq.0.or.offaxispulseflag.eq.0)) then
     if (spfrestrictflag.ne.0) then
        if (spfugrestrict.ne.0.and.tdflag.eq.0) then
           wwout%configtypes(wwout%firstconfig:wwout%lastconfig) = &
                wwout%configmvals(:)*wwout%configugvals(:)
        else
           wwout%configtypes(wwout%firstconfig:wwout%lastconfig)=wwout%configmvals(:)
        endif
     endif
#ifdef MPIFLAG
     if (wwout%parconsplit.ne.0) then
        call mpiallgather_i(wwout%configtypes,wwout%numconfig,wwout%configsperproc(:),&
             wwout%maxconfigsperproc)
     endif
#endif
  endif

end subroutine set_newconfiglist

end module
