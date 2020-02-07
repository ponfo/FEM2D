module FunctionOnSurfacesMOD
  use tools
  implicit none
 contains
  function funcOnSurfaces(i, x, y)
  implicit none
  integer(ikind), intent(in) :: i
  real(rkind), intent(in) :: x
  real(rkind), intent(in) :: y
  real(rkind) :: vecFunX(0)
  real(rkind) :: vecFunY(0)
  real(rkind), dimension(2) :: funcOnSurfaces
  funcOnSurfaces(1) = vecFunX(i)
  funcOnSurfaces(2) = vecFunY(i)
  end function funcOnSurfaces
end module FunctionOnSurfacesMOD
