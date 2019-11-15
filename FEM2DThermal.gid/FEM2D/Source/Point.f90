module PointMOD
  use tools
  implicit none
  private
  public :: PointTYPE, point
  type PointTYPE
     private
     integer(ikind) :: id
     real(rkind)    :: x
     real(rkind)    :: y
   contains
     procedure, public :: init
     procedure, public :: setID
     procedure, public :: getID
     procedure, public :: setX
     procedure, public :: getX
     procedure, public :: setY
     procedure, public :: getY
  end type PointTYPE
  interface point
     procedure constructor
  end interface point
contains
  type(PointTYPE) function constructor(id, x, y)
    implicit none
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    call constructor%init(id, x, y)
  end function constructor
  subroutine init(this, id, x, y)
    implicit none
    class(PointTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    this%id = id
    this%x = x
    this%y = y
  end subroutine init
  subroutine setID(this, id)
    implicit none
    class(PointTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    this%id = id
  end subroutine setID
  integer(ikind) function getID(this)
    implicit none
    class(PointTYPE), intent(inout) :: this
    getID = this%id
  end function getID
  subroutine setX(this, x)
    implicit none
    class(PointTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    this%x = x
  end subroutine setX
  real(rkind) function getX(this)
    implicit none
    class(PointTYPE), intent(inout) :: this
    getX = this%x
  end function getX
  subroutine setY(this, y)
    implicit none
    class(PointTYPE), intent(inout) :: this
    real(rkind), intent(in) :: y
    this%y = y
  end subroutine setY
  real(rkind) function getY(this)
    implicit none
    class(PointTYPE), intent(inout) :: this
    getY = this%y
  end function getY
end module PointMOD
