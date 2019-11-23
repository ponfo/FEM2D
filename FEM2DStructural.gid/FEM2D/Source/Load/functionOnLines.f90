module FunctionOnLinesMOD
  use tools
  implicit none
 contains
  function funcOnLines(i, x, y)
  implicit none
  integer(ikind), intent(in) :: i
  real(rkind), intent(in) :: x
  real(rkind), intent(in) :: y
  real(rkind) :: vecFunX(1)
  real(rkind) :: vecFunY(1)
  real(rkind), dimension(2) :: funcOnLines
  vecFunX(1)=       100
  vecFunX(1)=       0.0
  funcOnLines(1) = vecFunX(i)
  funcOnLines(2) = vecFunY(i)
  end function funcOnLines
end module FunctionOnLinesMOD
