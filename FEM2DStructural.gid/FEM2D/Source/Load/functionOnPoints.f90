module FunctionOnPointsMOD
  use tools
  implicit none
contains
  real(rkind) function funcOnPoints(i, x, y)
  implicit none
  integer(ikind), intent(in) :: i
  real(rkind), intent(in) :: x
  real(rkind), intent(in) :: y
  real(rkind) :: vecFunSource(1)
  vecFunSource(1)=   10
  funcOnPoints = vecfunSource(i)
  end function funcOnPoints
end module FunctionOnPointsMOD
