module PointPtrMOD
  use tools
  use PointMOD
  implicit none
  private
  public :: PointPtrTYPE
  type :: PointPtrTYPE
     type(PointTYPE), pointer :: ptr
   contains
     procedure, public :: allocate
     procedure, public :: setID
     procedure, public :: getID
     procedure, public :: setX
     procedure, public :: getX
     procedure, public :: setY
     procedure, public :: getY
  end type PointPtrTYPE
contains
  subroutine allocate(this, point)
    implicit none
    class(PointPtrTYPE), intent(inout) :: this
    type(PointTYPE), target, intent(in) :: point
    this%ptr => point
  end subroutine allocate
  subroutine setID(this, id)
    implicit none
    class(PointPtrTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    call this%ptr%setID(id)
  end subroutine setID
  integer(ikind) function getID(this)
    implicit none
    class(PointPtrTYPE), intent(inout) :: this
    getID = this%ptr%getID()
  end function getID
  subroutine setX(this, x)
    implicit none
    class(PointPtrTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    call this%ptr%setX(x)
  end subroutine setX
  real(rkind) function getX(this)
    implicit none
    class(PointPtrTYPE), intent(inout) :: this
    getX = this%ptr%getX()
  end function getX
  subroutine setY(this, y)
    implicit none
    class(PointPtrTYPE), intent(inout) :: this
    real(rkind), intent(in) :: y
    call this%ptr%setY(y)
  end subroutine setY
  real(rkind) function getY(this)
    implicit none
    class(PointPtrTYPE), intent(inout) :: this
    getY = this%ptr%getY()
  end function getY
end module PointPtrMOD

