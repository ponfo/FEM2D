module FunctionOnLinesMOD
  use tools
  implicit none
 contains
  function funcOnLines(i, x, y)
  implicit none
  integer(ikind), intent(in) :: i
  real(rkind), intent(in) :: x
  real(rkind), intent(in) :: y
  real(rkind) :: vecFunX(0)
  real(rkind) :: vecFunY(0)
  real(rkind), dimension(2) :: funcOnLines
  funcOnLines(1) = vecFunX(i)
  funcOnLines(2) = vecFunY(i)
  end function funcOnLines
end module FunctionOnLinesMOD
