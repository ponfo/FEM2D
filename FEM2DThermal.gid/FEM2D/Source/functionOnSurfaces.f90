module FunctionOnSurfacesMOD
  use tools
  implicit none
contains
  real(rkind) function funcOnSurfaces(i, x, y)
  implicit none
  integer(ikind), intent(in) :: i
  real(rkind), intent(in) :: x
  real(rkind), intent(in) :: y
  real(rkind) :: vecFunSource(1)
  vecFunSource(1)=   30
  funcOnSurfaces = vecfunSource(i)
  end function funcOnSurfaces
end module FunctionOnSurfacesMOD
