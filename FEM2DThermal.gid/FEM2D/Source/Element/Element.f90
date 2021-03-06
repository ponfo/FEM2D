module ElementMOD
  use tools
  
  use PointMOD
  use PointPtrMOD
  
  use IntegratorMOD
  use IntegratorPtrMOD
  
  use MaterialPtrMOD

  use GeometryPtrMOD
  implicit none
  private
  public :: ElementTYPE
  type, abstract :: ElementTYPE
     integer(ikind)                                :: id
     integer(ikind)                                :: nPoint
     integer(ikind)                                :: nDof
     type(PointPtrTYPE), dimension(:), allocatable :: point
     type(IntegratorPtrTYPE)                       :: integrator
     type(MaterialPtrTYPE)                         :: material
     type(GeometryPtrTYPE)                         :: geometry
   contains
     procedure, public  :: getID
     procedure, public  :: getnPoint
     procedure, public  :: getnDof
     procedure, public  :: getPointID
     procedure, public  :: getIntegrator
     procedure, public  :: getMaterial
     procedure, public  :: getGeometry
     procedure, public  :: setID
     procedure, public  :: setnPoint
     procedure, public  :: setnDof
     procedure, public  :: setPoint
     procedure, public  :: setIntegrator
     procedure, public  :: getOnePoint
     procedure, public  :: getAllPoints
  end type ElementTYPE
  
contains

  integer(ikind) function getID(this)
    implicit none
    class(ElementTYPE), intent(inout) :: this
    getID = this%id
  end function getID
  
  integer(ikind) function getnPoint(this)
    implicit none
    class(ElementTYPE), intent(inout) :: this
    getnPoint = this%nPoint
  end function getnPoint

  integer(ikind) function getnDof(this)
    implicit none
    class(ElementTYPE), intent(inout) :: this
    getnDof = this%nDof
  end function getnDof

  type(PointPtrTYPE) function getOnePoint(this, i)
    implicit none
    class(ElementTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getOnePoint = this%point(i)
  end function getOnePoint
  
  function getAllPoints(this)
    implicit none
    class(ElementTYPE), intent(inout) :: this
    type(PointPtrTYPE), dimension(this%nPoint) :: getAllPoints
    getAllPoints(1:this%nPoint) = this%point(1:this%nPoint)
  end function getAllPoints
  
  integer(ikind) function getPointID(this, i)
    implicit none
    class(ElementTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getPointID = this%point(i)%getID()
  end function getPointID

  type(IntegratorPtrTYPE) function getIntegrator(this)
    implicit none
    class(ElementTYPE), intent(inout) :: this
    getIntegrator = this%integrator
  end function getIntegrator

  type(MaterialPtrTYPE) function getMaterial(this)
    implicit none
    class(ElementTYPE), intent(inout) :: this
    getMaterial = this%material
  end function getMaterial

  type(GeometryPtrTYPE) function getGeometry(this)
    implicit none
    class(ElementTYPE), intent(inout) :: this
    getGeometry = this%geometry
  end function getGeometry

  subroutine setID(this, id)
    implicit none
    class(ElementTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    this%id = id
  end subroutine setID

  subroutine setnPoint(this, nPoint)
    implicit none
    class(ElementTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nPoint
    this%nPoint = nPoint
  end subroutine setnPoint

  subroutine setnDof(this, nDof)
    implicit none
    class(ElementTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nDof
    this%nDof = nDof
  end subroutine setnDof

  subroutine setPoint(this, i, point)
    implicit none
    class(ElementTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    class(PointTYPE), target, intent(in) :: point
    this%point(i)%ptr => point
  end subroutine setPoint

  subroutine setIntegrator(this, integrator)
    implicit none
    class(ElementTYPE), intent(inout) :: this
    type(IntegratorTYPE), target, intent(in) :: integrator
    this%integrator%ptr => integrator
  end subroutine setIntegrator

end module ElementMOD
