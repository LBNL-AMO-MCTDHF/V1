

!! ketot is matrix elements of -1/2 d^2/dx^2 in sinc function basis

subroutine sineDVR(ketot,fdtot, points,gridpoints,spacing)
  implicit none
  integer,intent(in) :: gridpoints
  real*8,intent(in) :: spacing
  real*8,intent(out) :: ketot(1-gridpoints:gridpoints-1),points(gridpoints),&
       fdtot(1-gridpoints:gridpoints-1)
  integer :: i
  real*8 ::  mypi

  myPi   = 4d0*atan(1d0)

  do i=1,gridpoints
     points(i) = i*spacing
  enddo

  points(:)=points(:)-(points(1)+points(gridpoints))/2d0

  do i=1-gridpoints,gridpoints-1
     if (i==0) then
        ketot(i) = myPi**2/3.d0
        fdtot(i) = 0d0
     else
        ketot(i) = 2.d0/real(i,8)**2 * (-1)**(i)
        fdtot(i) = 1d0/real(i,8) * (-1)**(i)
     endif
  enddo

  ketot=ketot/2.d0/spacing**2 
  fdtot=fdtot/spacing

end subroutine sineDVR



  





