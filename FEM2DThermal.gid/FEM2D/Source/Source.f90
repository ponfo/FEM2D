module SourceMOD
  use tools
  implicit none
  private
  public :: sourceTYPE
  type sourceTYPE
    contains
      procedure :: funcOnPoints
      procedure :: funcOnLines
      procedure :: funcOnSurfaces
  end type sourceTYPE
contains
  function funcOnPoints(this, i, x, y)
  implicit none
  integer(ikind), intent(in) :: i
  real(rkind), intent(in) :: x, y
  real(rkind) :: funcOnPoints
  class(sourceTYPE) :: this
  real(rkind) :: vecFunSource(0)
  funcOnPoints = vecfunSource(i)
  end function funcOnPoints
  function funcOnLines(this, i, x, y)
  implicit none
  integer(ikind), intent(in) :: i
  real(rkind), intent(in) :: x, y
  real(rkind) :: funcOnLines
  class(sourceTYPE) :: this
  real(rkind) :: vecFunSource(0)
  funcOnLines = vecfunSource(i)
  end function funcOnLines
  function funcOnSurfaces(this, i, x, y)
  implicit none
  integer(ikind), intent(in) :: i
  real(rkind), intent(in) :: x, y
  real(rkind) :: funcOnSurfaces
  class(sourceTYPE) :: this
  real(rkind) :: vecFunSource(0)
  funcOnSurfaces = vecfunSource(i)
  end function funcOnSurfaces
end module SourceMOD
