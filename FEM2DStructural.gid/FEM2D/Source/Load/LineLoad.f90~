module LineSourceMOD
  use tools
  use FunctionOnLinesMOD
  use PointMOD
  implicit none
  private
  public :: LineSourceTYPE, lineSource
  type LineSourceTYPE
     integer(ikind) :: iPoint
     integer(ikind) :: iSource
   contains
     procedure :: init
     procedure :: apply
  end type LineSourceTYPE
  interface lineSource
     procedure :: constructor
  end interface lineSource
contains
  type(LineSourceTYPE) function constructor(iPoint, iSource)
    implicit none
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iSource
    call constructor%init(iPoint, iSource)
  end function constructor
  subroutine init(this, iPoint, iSource)
    implicit none
    class(LineSourceTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iSource
    this%iPoint = iPoint
    this%iSource = iSource
  end subroutine init

  subroutine apply(this, point, rhs)
    implicit none
    class(LineSourceTYPE), intent(inout) :: this
    type(PointTYPE), dimension(:), intent(inout) :: point
    real(rkind), dimension(:), intent(inout) :: rhs
    real(rkind) :: val
    val = funcOnLines(this%iSource, point(this%iPoint)%getx(), point(this%iPoint)%gety())
    rhs(this%iPoint) = rhs(this%iPoint) + val
  end subroutine apply
    
end module LineSourceMOD
