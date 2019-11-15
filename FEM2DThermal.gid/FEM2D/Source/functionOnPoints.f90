module FunctionOnPointsMOD
  use tools
contains
  real(rkind) function funcOnPoints(i, x, y)
    implicit none
    integer(ikind), intent(in) :: i
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(1) :: vecFunSource
    vecFunSource(1) = 0.d0
    funcOnPoints = vecFunSource(i)
  end function funcOnPoints
end module FunctionOnPointsMOD
