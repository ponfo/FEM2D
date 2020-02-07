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
     real(rkind)    :: z
   contains
     procedure, public :: init
     procedure, public :: setID
     procedure, public :: getID
     procedure, public :: setX
     procedure, public :: getX
     procedure, public :: setY
     procedure, public :: getY
     procedure, public :: setZ
     procedure, public :: getZ
  end type PointTYPE
  interface point
     procedure constructor
  end interface point
contains
  type(PointTYPE) function constructor(id, x, y, z)
    implicit none
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), intent(in) :: z
    call constructor%init(id, x, y, z)
  end function constructor
  subroutine init(this, id, x, y, z)
    implicit none
    class(PointTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), intent(in) :: z
    this%id = id
    this%x = x
    this%y = y
    this%z = z
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
  subroutine setZ(this, z)
    implicit none
    class(PointTYPE), intent(inout) :: this
    real(rkind), intent(in) :: z
    this%z = z
  end subroutine setZ
  real(rkind) function getZ(this)
    implicit none
    class(PointTYPE), intent(inout) :: this
    getZ = this%z
  end function getZ
end module PointMOD
