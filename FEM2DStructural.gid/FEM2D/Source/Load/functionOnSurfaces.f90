module FunctionOnSurfacesMOD
  use tools
  implicit none
contains
  real(rkind) function funcOnSurfaces(i, x, y)
  implicit none
  integer(ikind), intent(in) :: i
  real(rkind), intent(in) :: x
  real(rkind), intent(in) :: y
  real(rkind) :: vecFunSource(2)
  vecFunSource(1)=   20
  vecFunSource(2)=   10
  funcOnSurfaces = vecfunSource(i)
  end function funcOnSurfaces
end module FunctionOnSurfacesMOD
