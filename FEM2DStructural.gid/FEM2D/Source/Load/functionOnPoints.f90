module FunctionOnPointsMOD
  use tools
  implicit none
 contains
  function funcOnPoints(i, x, y)
  implicit none
  integer(ikind), intent(in) :: i
  real(rkind), intent(in) :: x
  real(rkind), intent(in) :: y
  real(rkind) :: vecFunX(1)
  real(rkind) :: vecFunY(1)
  real(rkind), dimension(2) :: funcOnPoints
  vecFunX(1)=       100
  vecFunY(1)=       100
  funcOnPoints(1) = vecFunX(i)
  funcOnPoints(2) = vecFunY(i)
  end function funcOnPoints
end module FunctionOnPointsMOD
