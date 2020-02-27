module FunctionOnPointsMOD
  use tools
  implicit none
 contains
  function funcOnPoints(i, x, y)
  implicit none
  integer(ikind), intent(in) :: i
  real(rkind), intent(in) :: x
  real(rkind), intent(in) :: y
  real(rkind) :: vecFunX(5)
  real(rkind) :: vecFunY(5)
  real(rkind), dimension(2) :: funcOnPoints
  vecFunX(1)=       0.0
  vecFunY(1)=     -2000
  vecFunX(2)=       0.0
  vecFunY(2)=     -2000
  vecFunX(3)=       0.0
  vecFunY(3)=     -2000
  vecFunX(4)=       0.0
  vecFunY(4)=     -2000
  vecFunX(5)=       0.0
  vecFunY(5)=     -2000
  funcOnPoints(1) = vecFunX(i)
  funcOnPoints(2) = vecFunY(i)
  end function funcOnPoints
end module FunctionOnPointsMOD
