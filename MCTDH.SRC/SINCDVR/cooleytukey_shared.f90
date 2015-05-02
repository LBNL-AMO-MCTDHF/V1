
#define MAXPRIMES 7

!!$
!!$Apache License
!!$                           Version 2.0, January 2004
!!$                        http://www.apache.org/licenses/
!!$
!!$   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION
!!$
!!$   1. Definitions.
!!$
!!$      "License" shall mean the terms and conditions for use, reproduction,
!!$      and distribution as defined by Sections 1 through 9 of this document.
!!$
!!$      "Licensor" shall mean the copyright owner or entity authorized by
!!$      the copyright owner that is granting the License.
!!$
!!$      "Legal Entity" shall mean the union of the acting entity and all
!!$      other entities that control, are controlled by, or are under common
!!$      control with that entity. For the purposes of this definition,
!!$      "control" means (i) the power, direct or indirect, to cause the
!!$      direction or management of such entity, whether by contract or
!!$      otherwise, or (ii) ownership of fifty percent (50%) or more of the
!!$      outstanding shares, or (iii) beneficial ownership of such entity.
!!$
!!$      "You" (or "Your") shall mean an individual or Legal Entity
!!$      exercising permissions granted by this License.
!!$
!!$      "Source" form shall mean the preferred form for making modifications,
!!$      including but not limited to software source code, documentation
!!$      source, and configuration files.
!!$
!!$      "Object" form shall mean any form resulting from mechanical
!!$      transformation or translation of a Source form, including but
!!$      not limited to compiled object code, generated documentation,
!!$      and conversions to other media types.
!!$
!!$      "Work" shall mean the work of authorship, whether in Source or
!!$      Object form, made available under the License, as indicated by a
!!$      copyright notice that is included in or attached to the work
!!$      (an example is provided in the Appendix below).
!!$
!!$      "Derivative Works" shall mean any work, whether in Source or Object
!!$      form, that is based on (or derived from) the Work and for which the
!!$      editorial revisions, annotations, elaborations, or other modifications
!!$      represent, as a whole, an original work of authorship. For the purposes
!!$      of this License, Derivative Works shall not include works that remain
!!$      separable from, or merely link (or bind by name) to the interfaces of,
!!$      the Work and Derivative Works thereof.
!!$
!!$      "Contribution" shall mean any work of authorship, including
!!$      the original version of the Work and any modifications or additions
!!$      to that Work or Derivative Works thereof, that is intentionally
!!$      submitted to Licensor for inclusion in the Work by the copyright owner
!!$      or by an individual or Legal Entity authorized to submit on behalf of
!!$      the copyright owner. For the purposes of this definition, "submitted"
!!$      means any form of electronic, verbal, or written communication sent
!!$      to the Licensor or its representatives, including but not limited to
!!$      communication on electronic mailing lists, source code control systems,
!!$      and issue tracking systems that are managed by, or on behalf of, the
!!$      Licensor for the purpose of discussing and improving the Work, but
!!$      excluding communication that is conspicuously marked or otherwise
!!$      designated in writing by the copyright owner as "Not a Contribution."
!!$
!!$      "Contributor" shall mean Licensor and any individual or Legal Entity
!!$      on behalf of whom a Contribution has been received by Licensor and
!!$      subsequently incorporated within the Work.
!!$
!!$   2. Grant of Copyright License. Subject to the terms and conditions of
!!$      this License, each Contributor hereby grants to You a perpetual,
!!$      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
!!$      copyright license to reproduce, prepare Derivative Works of,
!!$      publicly display, publicly perform, sublicense, and distribute the
!!$      Work and such Derivative Works in Source or Object form.
!!$
!!$   3. Grant of Patent License. Subject to the terms and conditions of
!!$      this License, each Contributor hereby grants to You a perpetual,
!!$      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
!!$      (except as stated in this section) patent license to make, have made,
!!$      use, offer to sell, sell, import, and otherwise transfer the Work,
!!$      where such license applies only to those patent claims licensable
!!$      by such Contributor that are necessarily infringed by their
!!$      Contribution(s) alone or by combination of their Contribution(s)
!!$      with the Work to which such Contribution(s) was submitted. If You
!!$      institute patent litigation against any entity (including a
!!$      cross-claim or counterclaim in a lawsuit) alleging that the Work
!!$      or a Contribution incorporated within the Work constitutes direct
!!$      or contributory patent infringement, then any patent licenses
!!$      granted to You under this License for that Work shall terminate
!!$      as of the date such litigation is filed.
!!$
!!$   4. Redistribution. You may reproduce and distribute copies of the
!!$      Work or Derivative Works thereof in any medium, with or without
!!$      modifications, and in Source or Object form, provided that You
!!$      meet the following conditions:
!!$
!!$      (a) You must give any other recipients of the Work or
!!$          Derivative Works a copy of this License; and
!!$
!!$      (b) You must cause any modified files to carry prominent notices
!!$          stating that You changed the files; and
!!$
!!$      (c) You must retain, in the Source form of any Derivative Works
!!$          that You distribute, all copyright, patent, trademark, and
!!$          attribution notices from the Source form of the Work,
!!$          excluding those notices that do not pertain to any part of
!!$          the Derivative Works; and
!!$
!!$      (d) If the Work includes a "NOTICE" text file as part of its
!!$          distribution, then any Derivative Works that You distribute must
!!$          include a readable copy of the attribution notices contained
!!$          within such NOTICE file, excluding those notices that do not
!!$          pertain to any part of the Derivative Works, in at least one
!!$          of the following places: within a NOTICE text file distributed
!!$          as part of the Derivative Works; within the Source form or
!!$          documentation, if provided along with the Derivative Works; or,
!!$          within a display generated by the Derivative Works, if and
!!$          wherever such third-party notices normally appear. The contents
!!$          of the NOTICE file are for informational purposes only and
!!$          do not modify the License. You may add Your own attribution
!!$          notices within Derivative Works that You distribute, alongside
!!$          or as an addendum to the NOTICE text from the Work, provided
!!$          that such additional attribution notices cannot be construed
!!$          as modifying the License.
!!$
!!$      You may add Your own copyright statement to Your modifications and
!!$      may provide additional or different license terms and conditions
!!$      for use, reproduction, or distribution of Your modifications, or
!!$      for any such Derivative Works as a whole, provided Your use,
!!$      reproduction, and distribution of the Work otherwise complies with
!!$      the conditions stated in this License.
!!$
!!$   5. Submission of Contributions. Unless You explicitly state otherwise,
!!$      any Contribution intentionally submitted for inclusion in the Work
!!$      by You to the Licensor shall be under the terms and conditions of
!!$      this License, without any additional terms or conditions.
!!$      Notwithstanding the above, nothing herein shall supersede or modify
!!$      the terms of any separate license agreement you may have executed
!!$      with Licensor regarding such Contributions.
!!$
!!$   6. Trademarks. This License does not grant permission to use the trade
!!$      names, trademarks, service marks, or product names of the Licensor,
!!$      except as required for reasonable and customary use in describing the
!!$      origin of the Work and reproducing the content of the NOTICE file.
!!$
!!$   7. Disclaimer of Warranty. Unless required by applicable law or
!!$      agreed to in writing, Licensor provides the Work (and each
!!$      Contributor provides its Contributions) on an "AS IS" BASIS,
!!$      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
!!$      implied, including, without limitation, any warranties or conditions
!!$      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
!!$      PARTICULAR PURPOSE. You are solely responsible for determining the
!!$      appropriateness of using or redistributing the Work and assume any
!!$      risks associated with Your exercise of permissions under this License.
!!$
!!$   8. Limitation of Liability. In no event and under no legal theory,
!!$      whether in tort (including negligence), contract, or otherwise,
!!$      unless required by applicable law (such as deliberate and grossly
!!$      negligent acts) or agreed to in writing, shall any Contributor be
!!$      liable to You for damages, including any direct, indirect, special,
!!$      incidental, or consequential damages of any character arising as a
!!$      result of this License or out of the use or inability to use the
!!$      Work (including but not limited to damages for loss of goodwill,
!!$      work stoppage, computer failure or malfunction, or any and all
!!$      other commercial damages or losses), even if such Contributor
!!$      has been advised of the possibility of such damages.
!!$
!!$   9. Accepting Warranty or Additional Liability. While redistributing
!!$      the Work or Derivative Works thereof, You may choose to offer,
!!$      and charge a fee for, acceptance of support, warranty, indemnity,
!!$      or other liability obligations and/or rights consistent with this
!!$      License. However, in accepting such obligations, You may act only
!!$      on Your own behalf and on Your sole responsibility, not on behalf
!!$      of any other Contributor, and only if You agree to indemnify,
!!$      defend, and hold each Contributor harmless for any liability
!!$      incurred by, or claims asserted against, such Contributor by reason
!!$      of your accepting any such warranty or additional liability.
!!$
!!$   END OF TERMS AND CONDITIONS
!!$
!!$   Copyright 2015 the regents of the University of California
!!$
!!$   Licensed under the Apache License, Version 2.0 (the "License");
!!$   you may not use this file except in compliance with the License.
!!$   You may obtain a copy of the License at
!!$
!!$       http://www.apache.org/licenses/LICENSE-2.0
!!$
!!$   Unless required by applicable law or agreed to in writing, software
!!$   distributed under the License is distributed on an "AS IS" BASIS,
!!$   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!$   See the License for the specific language governing permissions and
!!$   limitations under the License.


module ct_fileptrmod
  implicit none
  integer :: mpifileptr = (-798)
end module ct_fileptrmod

module ct_mpimod
  implicit none
  integer :: myrank = -1
  integer :: nprocs = -1
  integer :: CT_COMM_WORLD = -1
  integer :: CT_GROUP_WORLD = -1
end module ct_mpimod


module ct_options
  integer :: ct_dimensionality=3
  integer :: ct_paropt=1
end module ct_options


module ct_primesetmod
  implicit none
  integer, allocatable :: CT_COMM_EACH(:,:),CT_GROUP_EACH(:,:),CT_PROCSET(:,:,:)
  integer :: CT_MYLOC(MAXPRIMES),CT_MYRANK(MAXPRIMES),CT_MYPOSITION(MAXPRIMES)
  integer :: ct_called=0
  integer :: ct_numprimes = (-1)
  integer :: ct_maxprime = (-1)
  integer :: ct_minprime = (-1)
  integer :: ct_pf(MAXPRIMES)=1,ct_maxposition(MAXPRIMES)=1

end module ct_primesetmod


subroutine ctdim(in_ctdim)
  use ct_fileptrmod
  use ct_mpimod
  use ct_options
  implicit none
  integer, intent(in) :: in_ctdim
  if (in_ctdim.ne.1.and.in_ctdim.ne.3) then 
     write(mpifileptr,*) "CTDIM NOT SUPPORTED", in_ctdim; call mpistop()
  endif
  ct_dimensionality=in_ctdim
end subroutine ctdim


subroutine ct_init(in_ctparopt,in_mpifileptr)
  use ct_fileptrmod
  use ct_mpimod
  use ct_options
  implicit none
  integer, intent(in) :: in_ctparopt,in_mpifileptr
  ct_paropt=in_ctparopt
  call getmyranknprocs(myrank,nprocs)
  call getworldcommgroup(CT_COMM_WORLD,CT_GROUP_WORLD)
  if (myrank.eq.1) then
     mpifileptr=6
  else
!!$     mpifileptr=in_mpifileptr
!!$     open(mpifileptr,file="/dev/null", status="unknown")
     mpifileptr=4910
     open(mpifileptr,file="/dev/null", status="unknown")
  endif

 call ct_getprimeset()

end subroutine ct_init


subroutine twiddlemult_mpi(blocksize,in,out,dim1,howmany,rdd)
  use ct_fileptrmod
  use ct_primesetmod
  implicit none
  integer, intent(in) :: blocksize,dim1,howmany,rdd
  complex*16, intent(in) :: in(blocksize,dim1,howmany)
  complex*16, intent(out) :: out(blocksize,dim1,howmany)
  complex*16 :: twiddle1(dim1,ct_maxposition(rdd)),tt1(dim1)
  integer :: ii,n1

  call gettwiddlesmall(twiddle1(:,:),dim1*ct_maxposition(rdd),ct_pf(rdd))

  tt1(:)=twiddle1(:,ct_myposition(rdd))**(ct_myrank(rdd)-1)
  do ii=1,howmany
     do n1=1,dim1
        out(:,n1,ii) = in(:,n1,ii) * tt1(n1)
     enddo
  enddo

end subroutine twiddlemult_mpi


subroutine myzfft1d_slowindex_mpi(in,out,totsize,rdd)
  use ct_fileptrmod
  use ct_options
  use ct_primesetmod
  implicit none
  integer, intent(in) :: totsize,rdd
  complex*16, intent(in) :: in(totsize)
  complex*16, intent(out) :: out(totsize)
  complex*16 :: fouriermatrix(ct_pf(rdd),ct_pf(rdd)),twiddle(ct_pf(rdd))
  integer :: ii

  call gettwiddlesmall(twiddle,ct_pf(rdd),1)
  do ii=1,ct_pf(rdd)
     fouriermatrix(:,ii)=twiddle(:)**(ii-1)
  enddo
  select case (ct_paropt)
  case(0)
  call simple_circ(in,out,fouriermatrix,totsize,rdd)
  case(1)
  call simple_summa(in,out,fouriermatrix,totsize,rdd)
  case default
     write(mpifileptr,*) "ct_paropt not recognized",ct_paropt; call mpistop()
  end select

end subroutine myzfft1d_slowindex_mpi


subroutine simple_circ(in, out,mat,howmany,rdd)
  use ct_fileptrmod
  use ct_mpimod
  use ct_primesetmod
  implicit none
  integer, intent(in) :: howmany,rdd
  complex*16, intent(in) :: in(howmany), mat(ct_pf(rdd),ct_pf(rdd))
  complex*16, intent(out) :: out(howmany)
  integer :: thisfileptr
#ifdef MPIFLAG
  complex*16 :: work2(howmany),work(howmany)
  integer :: ibox,jbox,deltabox,nnn
#endif

  thisfileptr=6

  if (rdd.lt.1.or.rdd.gt.ct_numprimes) then
     write(*,*) "recursion depth error circ",rdd,ct_numprimes; call mpistop()
  endif

#ifndef MPIFLAG
  out(:)=in(:)
  return
#else

!!$  if (ct_pf(rdd).eq.1) then
!!$     write(thisfileptr,*) "localnumprocs=1 will work, but please edit calling subroutine or "
!!$     write(thisfileptr,*) "  main input file and do not call cooleytukey with nprocs=1."
!!$     call mpistop()
!!$  endif

  nnn=1
  out(:)=0

  do deltabox=0,ct_pf(rdd)-1
     ibox=mod(ct_pf(rdd)+CT_MYRANK(rdd)-1+deltabox,ct_pf(rdd))+1
     jbox=mod(ct_pf(rdd)+CT_MYRANK(rdd)-1-deltabox,ct_pf(rdd))+1

     work(:)=in(:)*mat(ibox,CT_MYRANK(rdd))

     if (deltabox.ne.0) then
        call mympisendrecv_complex_local(work(:),work2(:),ibox,jbox,deltabox,howmany,CT_COMM_EACH(CT_MYLOC(rdd),rdd))
        out(:)=out(:)+work2(:)
     else
        out(:)=out(:)+work(:)
     endif
  enddo
#endif

end subroutine simple_circ


subroutine simple_summa(in, out,mat,howmany,rdd)
  use ct_fileptrmod
  use ct_mpimod
  use ct_primesetmod
  implicit none
  integer, intent(in) :: howmany,rdd
  complex*16, intent(in) :: in(howmany), mat(ct_pf(rdd),ct_pf(rdd))
  complex*16, intent(out) :: out(howmany)
  integer :: thisfileptr
#ifdef MPIFLAG
  complex*16 :: work(howmany)
  integer :: ibox,nnn
#endif

  thisfileptr=6

  if (rdd.lt.1.or.rdd.gt.ct_numprimes) then
     write(*,*) "recursion depth error circ",rdd,ct_numprimes; call mpistop()
  endif

#ifndef MPIFLAG
  out(:)=in(:)
  return
#else

!!$  if (ct_pf(rdd).eq.1) then
!!$     write(thisfileptr,*) "localnumprocs=1 will work, but please edit calling subroutine or "
!!$     write(thisfileptr,*) "  main program input file; do not call cooleytukey with nprocs=1."
!!$     call mpistop()
!!$  endif

  nnn=1
  out(:)=0d0

  do ibox=1,ct_pf(rdd)
     if (CT_MYRANK(rdd).eq.ibox) then
        work(:)=in(:)
     endif
     call mympicomplexbcast_local(work(:),ibox,howmany,CT_COMM_EACH(CT_MYLOC(rdd),rdd))
     out(:)=out(:)+work(:)*mat(CT_MYRANK(rdd),ibox)
  enddo
#endif

end subroutine simple_summa


  

subroutine getprimefactor(dim,myfactor)
  implicit none
  integer, intent(in) :: dim
  integer, intent(out) :: myfactor
  integer :: iprime
  integer, parameter :: numprimes=31
  integer, parameter :: primelist(numprimes)=&
       (/  2,  3,  5,  7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,&
          73, 79, 83, 89, 97,101,103,107,109,113,127 /)  ! no need to go remotely this high

  myfactor=dim
  do iprime=1,numprimes
     if (mod(dim,primelist(iprime)).eq.0) then
        myfactor=primelist(iprime)
        return
     endif
  enddo
end subroutine getprimefactor


subroutine getallprimefactors(dim,factormax,numfactors,allfactors)
  implicit none
  integer, intent(in) :: dim,factormax
  integer, intent(out) :: allfactors(factormax),numfactors
  integer :: thisdim,flag

  allfactors(:)=1
  numfactors=1
  thisdim=dim
  flag=0

798 if (numfactors.eq.factormax) then
     allfactors(factormax)=thisdim
     return
  else
     call getprimefactor(thisdim,allfactors(numfactors))
     if (allfactors(numfactors).eq.thisdim) then
        return
     endif
     thisdim=thisdim/allfactors(numfactors)
     numfactors=numfactors+1
  endif
go to 798
end subroutine getallprimefactors


subroutine gettwiddlesmall(twiddlefacs,dim1,dim2)
  implicit none
  integer, intent(in) :: dim1,dim2
  complex*16, intent(out) :: twiddlefacs(dim1)
  complex*16 :: phi
  integer :: k1, itwiddle(dim1)
  real*8, parameter :: pi=3.14159265358979323846264338327950d0
  phi=exp((0d0,-2d0) * pi / (dim1*dim2))
  do k1=1,dim1
     itwiddle(k1)=(k1-1)
  enddo
  twiddlefacs(:)=phi**itwiddle(:)
end subroutine gettwiddlesmall


!! PRIMESET SUBROUTINES


subroutine ct_getprimeset()
  use ct_fileptrmod
  use ct_mpimod
  use ct_primesetmod
  use ct_options
  implicit none
  integer :: ii

  if (ct_called.ne.0) then
     write(mpifileptr,*) "ONLY CALL CT_GETPRIMESET ONCE (programmer fail)"; call mpistop()
  endif
  ct_called=1

  call getallprimefactors(nprocs,MAXPRIMES,ct_numprimes,ct_pf)

  ct_maxprime=1; ct_minprime=32767
  do ii=1,ct_numprimes
     if (ct_pf(ii).gt.ct_maxprime) then
        ct_maxprime=ct_pf(ii)
     endif
     if (ct_pf(ii).lt.ct_minprime) then
        ct_minprime=ct_pf(ii)
     endif
  enddo
  
  write(mpifileptr,*)
  write(mpifileptr,*) "Go CT_INIT"
  write(mpifileptr,*) "   CT_MAXPRIME IS ",ct_maxprime
  write(mpifileptr,*) "    CT_PRIMEFACTORS ARE"
  write(mpifileptr,*) "  ",ct_pf(1:ct_numprimes)
  allocate(CT_COMM_EACH(nprocs/ct_minprime,ct_numprimes),CT_GROUP_EACH(nprocs/ct_minprime,ct_numprimes))
  CT_COMM_EACH(:,:)=(-42); CT_GROUP_EACH(:,:)=(-42)
  allocate(CT_PROCSET(ct_maxprime, nprocs/ct_minprime,ct_numprimes))
  CT_PROCSET(:,:,:)=1

  write(mpifileptr,*) "Calling ct_construct..."
  call ct_construct()
  write(mpifileptr,*) "   ....Called ct_construct."
  write(mpifileptr,*) 

end subroutine ct_getprimeset




subroutine ct_construct()
  use ct_fileptrmod
  use ct_mpimod
  use ct_primesetmod
  implicit none
#ifdef MPIFLAG
  integer :: thisfileptr,procshift(nprocs),ierr,iprime,&
       allprocs0(nprocs), proc_check, ii, &
       allprocs(ct_pf(7),ct_pf(6),ct_pf(5),ct_pf(4),ct_pf(3),ct_pf(2),ct_pf(1)),&
       qqtop(7),icomm
  integer, target :: qq(7),pp0(7),pp1(7)
  integer, pointer :: qq1,qq2,qq3,qq4,qq5,qq6,qq7

  if (ct_numprimes.gt.7) then
     print *, "REDIM CONSTRUCT",ct_numprimes,7; call mpistop()
  endif

  proc_check=ct_pf(1)*ct_pf(2)*ct_pf(3)*ct_pf(4)*ct_pf(5)*ct_pf(6)*ct_pf(7)
  if (proc_check.ne.nprocs) then
     write(mpifileptr,*) "Proc check programmer error ", proc_check,nprocs,&
          ct_pf(1:7); call mpistop()
  endif

  ct_maxposition(:)=1
  ct_maxposition(1:ct_numprimes)=nprocs
  do ii=1,ct_numprimes
     ct_maxposition(ii:ct_numprimes)=ct_maxposition(ii:ct_numprimes)/ct_pf(ii)
  enddo

  do ii=1,nprocs
     allprocs0(ii)=ii
  enddo

  allprocs(:,:,:,:,:,:,:)=RESHAPE(allprocs0,(/ct_pf(7),ct_pf(6),ct_pf(5),ct_pf(4),ct_pf(3),ct_pf(2),ct_pf(1)/))

  qq1=>qq(1); qq2=>qq(2); qq3=>qq(3); qq4=>qq(4); qq5=>qq(5); 
  qq6=>qq(6); qq7=>qq(7); 

  thisfileptr=6

  CT_MYLOC = (-99)
  CT_MYRANK = (-99)
  CT_MYPOSITION(:) = 1

  do iprime=1,ct_numprimes
     CT_MYPOSITION(iprime)=mod(myrank-1,ct_maxposition(iprime))+1
  enddo

  do iprime=1,ct_numprimes
     qqtop(:)=1
     qqtop(1:ct_numprimes)=ct_pf(1:ct_numprimes)
     qqtop(iprime)=1
     icomm=0

     do qq7=1,qqtop(7)
     do qq6=1,qqtop(6)
     do qq5=1,qqtop(5)
     do qq4=1,qqtop(4)
     do qq3=1,qqtop(3)
     do qq2=1,qqtop(2)
     do qq1=1,qqtop(1)

        pp0(1:7)=qq(1:7)
        pp1(1:7)=qq(1:7)
        pp0(iprime)=1
        pp1(iprime)=ct_pf(iprime)

        icomm=icomm+1
        if (icomm.gt.nprocs/ct_minprime) then
           write(mpifileptr,*) "Error construct",icomm,nprocs,ct_minprime; call mpistop()
        endif

        CT_PROCSET(1:ct_pf(iprime),icomm,iprime)=RESHAPE( &
             allprocs(pp0(7):pp1(7), pp0(6):pp1(6), pp0(5):pp1(5), pp0(4):pp1(4), &
             pp0(3):pp1(3), pp0(2):pp1(2), pp0(1):pp1(1)),(/ct_pf(iprime)/))

        do ii=1,ct_pf(iprime)
           if (CT_PROCSET(ii,icomm,iprime).eq.myrank) then
              if (CT_MYLOC(iprime).gt.0) then
                 print *, "ERROR MYLOC"; call mpistop()
              endif
              CT_MYLOC(iprime)=icomm
              CT_MYRANK(iprime)=ii
           endif
        enddo

        procshift(1:ct_pf(iprime))=CT_PROCSET(1:ct_pf(iprime),icomm,iprime) - 1

        call mpi_group_incl(CT_GROUP_WORLD,ct_pf(iprime),procshift,CT_GROUP_EACH(icomm,iprime),ierr)
        if (ierr.ne.0) then
           write(thisfileptr,*) "Error group incl CT",icomm,iprime,ierr; call mpistop()
        endif
        call mpi_comm_create(CT_COMM_WORLD, CT_GROUP_EACH(icomm,iprime), CT_COMM_EACH(icomm,iprime),ierr)
        if (ierr.ne.0) then
           write(thisfileptr,*) "Error comm create CT",icomm,iprime,ierr; call mpistop()
        endif

!!$ make destroy subroutine; keepme
!!$          call mpi_comm_free(CT_COMM_EACH(icomm,iprime),ierr)
!!$          if (ierr.ne.0) then
!!$             write(thisfileptr,*) "Error comm destroy CT",icomm,iprime,ierr; call mpistop()
!!$          endif
!!$          call mpi_group_free(CT_GROUP_EACH(icomm,iprime),ierr)
!!$          if (ierr.ne.0) then
!!$             write(thisfileptr,*) "Error group destroy CT",icomm,iprime,ierr; call mpistop()
!!$          endif

     enddo
     enddo
     enddo
     enddo
     enddo
     enddo
     enddo
     if (CT_MYLOC(iprime).le.0) then
        print *, "MYLOC ERROR",myrank,CT_MYLOC(iprime),iprime; call mpistop()
     endif
  enddo
#endif

end subroutine ct_construct



