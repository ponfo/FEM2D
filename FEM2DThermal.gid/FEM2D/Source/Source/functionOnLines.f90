module FunctionOnLinesMOD
  use tools
  implicit none
contains
  real(rkind) function funcOnLines(i, x, y)
  implicit none
  integer(ikind), intent(in) :: i
  real(rkind), intent(in) :: x
  real(rkind), intent(in) :: y
  real(rkind) :: vecFunSource(0)
  funcOnLines = vecfunSource(i)
  end function funcOnLines
end module FunctionOnLinesMOD
