module PointSourceMOD
  use tools
  use FunctionOnPointsMOD
  use PointMOD
  implicit none
  private
  public :: PointSourceTYPE, pointSource
  type PointSourceTYPE
     integer(ikind) :: iPoint
     integer(ikind) :: iSource
   contains
     procedure :: init
     procedure :: apply
  end type PointSourceTYPE
  
  interface pointSource
     procedure :: constructor
  end interface pointSource
  
contains
  
  type(PointSourceTYPE) function constructor(iPoint, iSource)
    implicit none
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iSource
    call constructor%init(iPoint, iSource)
  end function constructor
  
  subroutine init(this, iPoint, iSource)
    implicit none
    class(PointSourceTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iSource
    this%iPoint = iPoint
    this%iSource = iSource
  end subroutine init

  subroutine apply(this, point, rhs)
    implicit none
    class(PointSourceTYPE), intent(inout) :: this
    type(PointTYPE), dimension(:), intent(inout) :: point
    real(rkind), dimension(:), intent(inout) :: rhs
    real(rkind) :: val
    val = funcOnPoints(this%iSource, point(this%iPoint)%getx(), point(this%iPoint)%gety())
    rhs(this%iPoint) = rhs(this%iPoint) + val
  end subroutine apply

end module PointSourceMOD
