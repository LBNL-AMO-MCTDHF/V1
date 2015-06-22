
#define MAXPRIMES 10

#include "Definitions.INC"

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



module ct_options
  integer :: ct_paropt=1
end module ct_options


module ct_primesetmod
  implicit none

  integer, allocatable :: CT_COMM_EACH(:,:,:),CT_GROUP_EACH(:,:,:)
  integer :: CT_MYLOC(MAXPRIMES,3),CT_MYRANK(MAXPRIMES,3),CT_MYPOSITION(MAXPRIMES,3)

  integer :: ct_called=0
  integer :: ct_numprimes = (-1)
  integer :: ct_maxprime = (-1)
  integer :: ct_minprime = (-1)
  integer :: ct_pf(MAXPRIMES)=1,ct_maxposition(MAXPRIMES)=1

end module ct_primesetmod



subroutine ct_init(in_ctparopt)
  use pmpimod
  use ct_options
  implicit none
  integer, intent(in) :: in_ctparopt

  ct_paropt=in_ctparopt

  call ct_getprimeset()
  
end subroutine ct_init


subroutine twiddlemult_mpi(blocksize,in,out,dim1,howmany,rdd,oplevel)
  use ct_primesetmod
  implicit none
  integer, intent(in) :: blocksize,dim1,howmany,rdd,oplevel
  complex*16, intent(in) :: in(blocksize,dim1,howmany)
  complex*16, intent(out) :: out(blocksize,dim1,howmany)
  complex*16 :: twiddle1(dim1,ct_maxposition(rdd)),tt1(dim1)
  integer :: ii,n1

  call gettwiddlesmall(twiddle1(:,:),dim1*ct_maxposition(rdd),ct_pf(rdd))

  tt1(:)=twiddle1(:,ct_myposition(rdd,oplevel))**(ct_myrank(rdd,oplevel)-1)
  do ii=1,howmany
     do n1=1,dim1
        out(:,n1,ii) = in(:,n1,ii) * tt1(n1)
     enddo
  enddo

end subroutine twiddlemult_mpi


subroutine myzfft1d_slowindex_mpi(in,out,totsize,rdd,oplevel)
  use pfileptrmod
  use ct_options
  use ct_primesetmod
  implicit none
  integer, intent(in) :: totsize,rdd,oplevel
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
  call simple_circ(in,out,fouriermatrix,totsize,rdd,oplevel)
  case(1)
  call simple_summa(in,out,fouriermatrix,totsize,rdd,oplevel)
  case default
     OFLWR "ct_paropt not recognized",ct_paropt; CFLST
  end select

end subroutine myzfft1d_slowindex_mpi


subroutine simple_circ(in, out,mat,howmany,rdd,oplevel)
  use pmpimod
  use ct_primesetmod
  implicit none
  integer, intent(in) :: howmany,rdd,oplevel
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
  out(1)=mat(1,1) !! avoid warn unused
#else

  nnn=1
  out(:)=0

  do deltabox=0,ct_pf(rdd)-1
     ibox=mod(ct_pf(rdd)+CT_MYRANK(rdd,oplevel)-1+deltabox,ct_pf(rdd))+1
     jbox=mod(ct_pf(rdd)+CT_MYRANK(rdd,oplevel)-1-deltabox,ct_pf(rdd))+1

     work(:)=in(:)*mat(ibox,CT_MYRANK(rdd,oplevel))

     if (deltabox.ne.0) then
        call mympisendrecv_complex_local(work(:),work2(:),ibox,jbox,deltabox,howmany,CT_COMM_EACH(CT_MYLOC(rdd,oplevel),rdd,oplevel))
        out(:)=out(:)+work2(:)
     else
        out(:)=out(:)+work(:)
     endif
  enddo
#endif

end subroutine simple_circ


subroutine simple_summa(in, out,mat,howmany,rdd,oplevel)
  use pmpimod
  use ct_primesetmod
  implicit none
  integer, intent(in) :: howmany,rdd,oplevel
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
  out(1)=mat(1,1) !! avoid warn unused
#else

  nnn=1
  out(:)=0d0

  do ibox=1,ct_pf(rdd)
     if (CT_MYRANK(rdd,oplevel).eq.ibox) then
        work(:)=in(:)
     endif
     call mympicomplexbcast_local(work(:),ibox,howmany,CT_COMM_EACH(CT_MYLOC(rdd,oplevel),rdd,oplevel))
     out(:)=out(:)+work(:)*mat(CT_MYRANK(rdd,oplevel),ibox)
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
  use pmpimod
  use pfileptrmod
  use ct_primesetmod
  use ct_options
  implicit none
  integer :: ii, allprocs0(procsplit(1),procsplit(2),procsplit(3)),jj,kk,iproc

  if (ct_called.ne.0) then
     OFLWR "ONLY CALL CT_GETPRIMESET ONCE (programmer fail)"; CFLST
  endif
  ct_called=1

  ct_maxposition(:)=1

  select case(orbparlevel)
  case(3)
     call getallprimefactors(nprocs,MAXPRIMES,ct_numprimes,ct_pf)
     ct_maxposition(1:ct_numprimes)=nprocs
  case(2)
     call getallprimefactors(sqnprocs,MAXPRIMES,ct_numprimes,ct_pf)
     ct_maxposition(1:ct_numprimes)=sqnprocs
  case(1)
     call getallprimefactors(cbnprocs,MAXPRIMES,ct_numprimes,ct_pf)
     ct_maxposition(1:ct_numprimes)=cbnprocs
  end select


  do ii=1,ct_numprimes
     ct_maxposition(ii:ct_numprimes)=ct_maxposition(ii:ct_numprimes)/ct_pf(ii)
  enddo



  ct_maxprime=1; ct_minprime=32767
  do ii=1,ct_numprimes
     if (ct_pf(ii).gt.ct_maxprime) then
        ct_maxprime=ct_pf(ii)
     endif
     if (ct_pf(ii).lt.ct_minprime) then
        ct_minprime=ct_pf(ii)
     endif
  enddo
  
  OFLWR
  WRFL "Go CT_INIT"
  WRFL "   CT_MAXPRIME IS ",ct_maxprime
  WRFL "    CT_PRIMEFACTORS ARE"
  WRFL "  ",ct_pf(1:ct_numprimes)
  CFL

  allocate(CT_COMM_EACH(procsplit(3)/ct_minprime,ct_numprimes,orbparlevel:3),&
       CT_GROUP_EACH(procsplit(3)/ct_minprime,ct_numprimes,orbparlevel:3))
  CT_COMM_EACH(:,:,:)=(-42); 
  CT_GROUP_EACH(:,:,:)=(-42)
 
  OFLWR "Calling ct_construct..."; CFL

  iproc=0
  do ii=1,procsplit(3)
     do jj=1,procsplit(2)
        do kk=1,procsplit(1)
           iproc=iproc+1
           allprocs0(kk,jj,ii)=iproc
        enddo
     enddo
  enddo

  call ct_construct(allprocs0)

!  do ii=orbparlevel,3
!     do jj=1,nbox(2)
!        do kk=1,nbox(1)
!           call ct_construct(allprocs0(:,jj,kk)) !!,CTBOX(kk,jj,ii))
!        enddo
!     enddo
!  enddo

  OFLWR "   ....Called ct_construct."; CFL


end subroutine ct_getprimeset


subroutine ct_construct(allprocs0)
  use pmpimod
  use pfileptrmod
  use ct_primesetmod
  implicit none
  integer :: allprocs0(procsplit(1),procsplit(2),procsplit(3))
#ifdef MPIFLAG
  integer :: thisfileptr,procshift(ct_maxprime),ierr,iprime,&
       proc_check, ii, icomm, ilevel, &
       xxtop(1:ct_numprimes),yytop(1:ct_numprimes),xx,yy
  integer :: procset(ct_maxprime)

  proc_check=1; xxtop(:)=1; yytop(:)=1
  do ii=1,ct_numprimes
     proc_check=proc_check*ct_pf(ii)
     yytop(1:ii-1)=yytop(1:ii-1)*ct_pf(ii)
     xxtop(ii+1: )=xxtop(ii+1:)*ct_pf(ii)
  enddo
  if (proc_check.ne.procsplit(3)) then
     OFLWR "Proc check programmer error ", proc_check,procsplit(3),&
          ct_pf(1:ct_numprimes); CFLST
  endif

  thisfileptr=6

  CT_MYLOC = (-99)
  CT_MYRANK = (-99)
  CT_MYPOSITION = 1

  do ilevel=orbparlevel,3

     do iprime=1,ct_numprimes
        CT_MYPOSITION(iprime,ilevel)=mod(boxrank(ilevel)-1,ct_maxposition(iprime))+1
     enddo

     do iprime=1,ct_numprimes

        icomm=0

        do xx=1, xxtop(iprime)    !!  EITHER LOOP ORDER APPEARS OK
        do yy=1, yytop(iprime)    !!  CHECKME (with qqtop etc was fast index outer, slow inner)

           icomm=icomm+1
           if (icomm.gt.procsplit(ilevel)/ct_minprime) then
              OFLWR  "Error construct",icomm,procsplit(ilevel),ct_minprime; CFLST
           endif
           
           select case(ilevel)
           case(3)
              call ct_cast_assign(allprocs0(boxrank(1),boxrank(2),:),yytop(iprime),ct_pf(iprime),xxtop(iprime),yy,xx,&
                   procset(1:ct_pf(iprime)))
           case(2)
              call ct_cast_assign(allprocs0(boxrank(1),:,boxrank(3)),yytop(iprime),ct_pf(iprime),xxtop(iprime),yy,xx,&
                   procset(1:ct_pf(iprime)))
           case(1)
              call ct_cast_assign(allprocs0(:,boxrank(2),boxrank(3)),yytop(iprime),ct_pf(iprime),xxtop(iprime),yy,xx,&
                   procset(1:ct_pf(iprime)))
           end select

           do ii=1,ct_pf(iprime)
              if (procset(ii).eq.myrank) then
                 if (CT_MYLOC(iprime,ilevel).gt.0) then
                    print *, "ERROR MYLOC"; call mpistop()
                 endif
                 CT_MYLOC(iprime,ilevel)=icomm
                 CT_MYRANK(iprime,ilevel)=ii
              endif
           enddo
           
           procshift(1:ct_pf(iprime))=procset(1:ct_pf(iprime)) - 1

           call mpi_group_incl(PROJ_GROUP_WORLD,ct_pf(iprime),procshift,CT_GROUP_EACH(icomm,iprime,ilevel),ierr)

           if (ierr.ne.0) then
              write(thisfileptr,*) "Error group incl CT",icomm,iprime,ierr; call mpistop()
           endif

           call mpi_comm_create(PROJ_COMM_WORLD, CT_GROUP_EACH(icomm,iprime,ilevel), CT_COMM_EACH(icomm,iprime,ilevel),ierr)

           if (ierr.ne.0) then
              write(thisfileptr,*) "Error comm create CT",icomm,iprime,ierr; call mpistop()
           endif
        enddo
        enddo
        if (CT_MYLOC(iprime,ilevel).le.0) then
           print *, "MYLOC ERROR",myrank,boxrank(ilevel),CT_MYLOC(iprime,ilevel),iprime; call mpistop()
        endif
     enddo
  enddo
#endif

contains
  subroutine ct_cast_assign(threemat,first,middle,last,x1,x3,middle_out)
    implicit none
    integer, intent(in) :: first,middle,last,x1,x3, threemat(first,middle,last)
    integer, intent(out) :: middle_out(middle)
    middle_out(:)=threemat(x1,:,x3)
  end subroutine ct_cast_assign
  
end subroutine ct_construct



