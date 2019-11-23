module PointLoadMOD
  use tools
  use FunctionOnPointsMOD
  use PointMOD
  implicit none
  private
  public :: PointLoadTYPE, pointLoad
  type PointLoadTYPE
     integer(ikind) :: iPoint
     integer(ikind) :: iLoad
   contains
     procedure :: init
     procedure :: apply
  end type PointLoadTYPE
  
  interface pointLoad
     procedure :: constructor
  end interface pointLoad
  
contains
  
  type(PointLoadTYPE) function constructor(iPoint, iLoad)
    implicit none
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iLoad
    call constructor%init(iPoint, iLoad)
  end function constructor
  
  subroutine init(this, iPoint, iLoad)
    implicit none
    class(PointLoadTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iLoad
    this%iPoint = iPoint
    this%iLoad = iLoad
  end subroutine init

  subroutine apply(this, point, rhs)
    implicit none
    class(PointLoadTYPE), intent(inout) :: this
    type(PointTYPE), dimension(:), intent(inout) :: point
    real(rkind), dimension(:), intent(inout) :: rhs
    real(rkind), dimension(2) :: val
    val = funcOnPoints(this%iLoad, point(this%iPoint)%getx(), point(this%iPoint)%gety())
    rhs(2*this%iPoint-1) = rhs(2*this%iPoint-1) + val(1)
    rhs(2*this%iPoint) = rhs(2*this%iPoint) + val(2)
  end subroutine apply

end module PointLoadMOD
