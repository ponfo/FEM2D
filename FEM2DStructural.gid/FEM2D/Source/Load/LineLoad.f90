module LineLoadMOD
  use tools
  use FunctionOnLinesMOD
  use PointMOD
  implicit none
  private
  public :: LineLoadTYPE, lineLoad
  type LineLoadTYPE
     integer(ikind) :: iPoint
     integer(ikind) :: iLoad
   contains
     procedure :: init
     procedure :: apply
  end type LineLoadTYPE
  interface lineLoad
     procedure :: constructor
  end interface lineLoad
contains
  type(LineLoadTYPE) function constructor(iPoint, iLoad)
    implicit none
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iLoad
    call constructor%init(iPoint, iLoad)
  end function constructor
  subroutine init(this, iPoint, iLoad)
    implicit none
    class(LineLoadTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iLoad
    this%iPoint = iPoint
    this%iLoad = iLoad
  end subroutine init

  subroutine apply(this, point, rhs)
    implicit none
    class(LineLoadTYPE), intent(inout) :: this
    type(PointTYPE), dimension(:), intent(inout) :: point
    real(rkind), dimension(:), intent(inout) :: rhs
    real(rkind) :: val
    val = funcOnLines(this%iLoad, point(this%iPoint)%getx(), point(this%iPoint)%gety())
    rhs(this%iPoint) = rhs(this%iPoint) + val
  end subroutine apply
    
end module LineLoadMOD
