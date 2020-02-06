module DomainMOD
  use tools
  use DebuggerMOD

  use PointMOD
  use PointPtrMOD
  implicit none
  private
  public :: DomainTYPE
  type :: DomainTYPE
     integer(ikind)                                :: nPoint
     integer(ikind)                                :: nLine
     integer(ikind)                                :: nTriang
     integer(ikind)                                :: nQuad
     integer(ikind)                                :: nElem
     type(PointTYPE)   , dimension(:), allocatable :: point
   contains
     procedure, public :: getnPoint
     procedure, public :: getnLine
     procedure, public :: getnTriang
     procedure, public :: getnQuad
     procedure, public :: getnElem
     procedure, public :: getPoint
  end type DomainTYPE

contains

  integer(ikind) function getnPoint(this)
    implicit none
    class(DomainTYPE), intent(inout) :: this
    getnPoint = this%nPoint
  end function getnPoint

  integer(ikind) function getnLine(this)
    implicit none
    class(DomainTYPE), intent(inout) :: this
    getnLine = this%nLine
  end function getnLine
  
  integer(ikind) function getnTriang(this)
    implicit none
    class(DomainTYPE), intent(inout) :: this
    getnTriang = this%nTriang
  end function getnTriang
  
  integer(ikind) function getnQuad(this)
    implicit none
    class(DomainTYPE), intent(inout) :: this
    getnQuad = this%nQuad
  end function getnQuad
  
  integer(ikind) function getnElem(this)
    implicit none
    class(DomainTYPE), intent(inout) :: this
    getnElem = this%nElem
  end function getnElem

  type(PointPtrTYPE) function getPoint(this, idx)
    implicit none
    class(DomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: idx
    call getPoint%allocate(this%point(idx))
  end function getPoint

end module DomainMOD
