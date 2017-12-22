

!! FOR SAVING/PLOTTING ACTIONS 2-12

#include "Definitions.INC"

subroutine save_density( thistime, inspfs, indenmat, iprop, denfilename)
  use parameters
  use mpimod
  implicit none
  DATATYPE,intent(in) :: inspfs(spfdims(1),spfdims(2),spfdims(3),nspf), indenmat(nspf,nspf)
  real*8,intent(in) :: thistime
  character,intent(in) :: denfilename*(*)
  integer,intent(in) :: iprop
  integer :: i,j, imval
  character (len=headersize) :: header
  complex*16 :: cmdensity(spfdims(1),spfdims(2),spfdims(3)),mtrans(spfdims(3),spfdims(3)) 
  DATATYPE :: mdensity(spfdims(1),spfdims(2),spfdims(3))                       
  complex*16 :: density(spfdims(1),spfdims(2),spfdims(3))                         !! AUTOMATIC

!!$  This was for transforming to phi values from m.
!!$    But menu asks for m-projection.  Perhaps later, re-implement this
!!$
!!$  if (spfdimtype(3).eq.1) then  !! assume fourier basis (-mbig:mbig)
!!$     if (mod(spfdims(3),2).ne.1) then
!!$        OFLWR "FOURRR ERROR"; CFLST
!!$     endif
!!$     do i=1,spfdims(3)
!!$        do j=1,spfdims(3)
!!$           jj=j-(spfdims(3)+1)/2
!!$           mtrans(i,j) = exp((0.d0,1.d0)*jj*2*pi*i/real(spfdims(3)))
!!$        enddo
!!$     enddo
!!$  else

     mtrans(:,:)=0d0
     do i=1,spfdims(3)
        mtrans(i,i)=1d0
     enddo

!!$  endif

  call getdensity(density, indenmat, inspfs,nspf)
  call zgemm('N', 'N', spfdims(1)*spfdims(2),spfdims(3),spfdims(3), (1.d0,0.d0), &
       density, spfdims(1)*spfdims(2), mtrans,spfdims(3), (0.d0,0.d0),&
       cmdensity, spfdims(1)*spfdims(2))

!! why square root?  is this a mistake?
!! no not a mistake; storing DVR coefficients of density
!! another sqrt(weights) is divided in cylindricalvalue.

  do imval=1,spfdims(3)
     mdensity(:,:,imval)=  cmdensity(:,:,imval) / &          !! ok conversion
          sqrt(elecweights(:,:,imval,1)* &                   !! ok conversion
          elecweights(:,:,imval,2)*elecweights(:,:,imval,3)) !! ok conversion
  enddo
  call write_den_header(header,thistime,iprop)
  call save_orbvector(mdensity, spfsize, denfile, denfilename, header)

end subroutine save_density


subroutine save_natorb( thistime, inspfs, indenvects, indenvals, iprop)
  use parameters
  implicit none
  DATATYPE,intent(in) :: inspfs(spfsize,nspf), indenvects(nspf,nspf)
  CNORMTYPE,intent(in) :: indenvals(nspf)
  integer,intent(in) :: iprop
  real*8,intent(in) :: thistime
  DATATYPE :: tempspfs(spfsize,nspf)  !! AUTOMATIC
  integer :: i,j,ispf
  character (len=headersize) :: header
  DATATYPE :: csum

  tempspfs=0.d0
  do j=1,nspf  ! which natorb
     do i=1,nspf  ! which original

!!        tempspfs(:,j)=tempspfs(:,j)+ inspfs(:,i)*CONJUGATE(indenvects(i,j))

        tempspfs(:,j)=tempspfs(:,j)+ inspfs(:,i)*indenvects(i,j)

     enddo
  enddo
  do ispf=1,nspf
     csum=indenvals(ispf)
     call write_nat_header(header,thistime,iprop,ispf,csum)
     call save_orbvector(tempspfs(:,ispf), spfsize, natorbfile, natplotbin, header)
  enddo
end subroutine save_natorb

subroutine save_spf( thistime, inspfs, iprop)
  use parameters
  implicit none
  DATATYPE,intent(in) :: inspfs(spfsize,nspf)
  real*8,intent(in) :: thistime
  integer,intent(in) :: iprop
  character (len=headersize) :: header
  integer :: ispf

  do ispf=1,nspf
     call write_spf_header(header,thistime,iprop,ispf)
     call save_orbvector(inspfs(:,ispf), spfsize, spfplotfile, spfplotbin, header)
  enddo
end subroutine save_spf



subroutine save_rnatorb( thistime, inrdenvects, inrdenvals, iprop)
  use parameters
  implicit none
  DATATYPE,intent(in) :: inrdenvects(numr,numr)
  CNORMTYPE,intent(in) :: inrdenvals(numr)
  real*8,intent(in) :: thistime
  integer :: iorb,iprop
  character (len=headersize) :: header
  DATATYPE :: csum, tempvects(numr,numr)

  tempvects=inrdenvects

!  do iorb=1,numr
!     call fixphase0(tempvects(:,iorb),numr)
!  enddo

  OFLWR " Saving rnatorb !";CFL

!! noo way too many  do iorb=1,numr
  do iorb=1,min(numr,min(num_config,10))
     csum=inrdenvals(iorb)
     call write_nat_header(header,thistime,iprop,iorb,csum)
     call save_orbvector(tempvects(:,iorb), numr, rnatorbfile, rnatplotbin, header)
  enddo
end subroutine save_rnatorb


subroutine write_den_header(header,thistime,readprop)
  implicit none
  character,intent(out) :: header*(*)
  real*8,intent(in) :: thistime
  integer,intent(in) :: readprop
  write(header,*) thistime,readprop
end subroutine

subroutine write_nat_header(header,thistime,readprop,ispf,mydenval)
  implicit none
  character,intent(out) :: header*(*)
  real*8,intent(in) :: thistime
  integer,intent(in) :: readprop, ispf
  DATATYPE,intent(in) :: mydenval
  complex*16 :: csum
  csum=mydenval
  write(header,*) thistime,readprop, ispf, csum
end subroutine

subroutine write_spf_header(header,thistime,readprop,ispf)
  implicit none
  character,intent(out) :: header*(*)
  real*8,intent(in) :: thistime
  integer,intent(in) :: readprop, ispf
  write(header,*) thistime,readprop,ispf
end subroutine write_spf_header

subroutine read_den_header(header,thistime,readprop)
  implicit none
  character,intent(in) :: header*(*)
  real*8,intent(out) :: thistime
  integer,intent(out) :: readprop
  read(header,*) thistime,readprop
end subroutine

subroutine read_nat_header(header,thistime,readprop,ispf,mydenval)
  implicit none
  character,intent(in) :: header*(*)
  real*8,intent(out) :: thistime
  integer,intent(out) :: readprop, ispf
  DATATYPE,intent(out) :: mydenval
  complex*16 :: csum
  read(header,*) thistime,readprop, ispf, csum
  mydenval=csum !! ok imp conv mctdh
end subroutine read_nat_header

subroutine read_spf_header(header,thistime,readprop,ispf)
  implicit none
  character,intent(in) :: header*(*)
  real*8,intent(out) :: thistime
  integer,intent(out) :: readprop, ispf
  read(header,*) thistime,readprop,ispf
end subroutine read_spf_header

subroutine save_natorb_final()
  use parameters
  implicit none
  call close_orbvector(natorbfile)
end subroutine save_natorb_final

subroutine save_rnatorb_final()
  use parameters
  implicit none
  call close_orbvector(rnatorbfile)
end subroutine save_rnatorb_final

subroutine save_spf_final()
  use parameters
  implicit none
  call close_orbvector(spfplotfile)
end subroutine save_spf_final

subroutine save_density_final()
  use parameters
  implicit none
  call close_orbvector(denfile)
end subroutine save_density_final

subroutine save_natorb_initial()
  use parameters
  use xxxmod
  use mpimod
  implicit none
  if (myrank.eq.1) then
     open(natorbfile,file=natplotbin, status="replace", form="unformatted")
     close(natorbfile)
  endif
  call save_natorb( -1.d0, yyy%cmfspfs(:,0), yyy%denvects, yyy%denvals , 1)
end subroutine save_natorb_initial

subroutine save_density_initial(denfilename)
  use parameters
  use xxxmod
  use mpimod
  implicit none
  character,intent(in) :: denfilename*(*)
  if (myrank.eq.1) then
     open(denfile,file=denfilename, status="replace", form="unformatted")
     close(denfile)
  endif
  call save_density( -1.d0, yyy%cmfspfs(:,0), yyy%denmat(:,:,0) , 1, denfilename)
end subroutine save_density_initial

subroutine save_rnatorb_initial()
  use parameters
  use xxxmod
  use mpimod
  implicit none
  OFLWR "REDO SAVE RNAT - NEED TO CALC"; CFLST
!  open(rnatorbfile,file=rnatplotbin, status="replace", form="unformatted")
!  close(rnatorbfile)
!  call save_rnatorb(-1.d0, yyy%rdenvects(:,:), yyy%rdenvals(:), 1)
end subroutine save_rnatorb_initial

subroutine save_spf_initial()
  use parameters
  use xxxmod
  use mpimod
  implicit none
  if (myrank.eq.1) then
     open(spfplotfile,file=spfplotbin, status="replace", form="unformatted")
     close(spfplotfile)
  endif
  call save_spf( -1.d0, yyy%cmfspfs(:,0), 1)
end subroutine save_spf_initial

