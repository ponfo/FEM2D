module Element2DPtrMOD
  use tools
  use PointMOD
  use PointPtrMOD
  use MaterialPtrMOD
  use GeometryPtrMOD
  use IntegratorMOD
  use IntegratorPtrMOD
  use ElementMOD
  use Element2DMOD
  implicit none
  private
  public :: Element2DPtrTYPE
  type :: Element2DPtrTYPE
     class(Element2DTYPE), pointer :: ptr
   contains
     procedure, public  :: getID
     procedure, public  :: getnPoint
     procedure, public  :: getnDof
     procedure, public  :: getPointID
     procedure, public  :: getMaterial
     procedure, public  :: getGeometry
     procedure, public  :: getIntegrator
     procedure, public  :: getArea
     procedure, public  :: setID
     procedure, public  :: setnPoint
     procedure, public  :: setnDof
     procedure, public  :: setPoint
     procedure, public  :: setIntegrator
     procedure, public  :: setArea
     procedure, public  :: jacobian
     procedure, public  :: jacobianDet
     procedure, public  :: shapeFunc
     procedure, public  :: dShapeFunc
     procedure, public  :: getStiffness
     generic  , public  :: getPoint => getOnePoint, getAllPoints
     procedure, private :: getOnePoint
     procedure, private :: getAllPoints
  end type Element2DPtrTYPE

contains

  integer(ikind) function getID(this)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    getID = this%ptr%id
  end function getID

  integer(ikind) function getnPoint(this)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    getnPoint = this%ptr%nPoint
  end function getnPoint

  integer(ikind) function getnDof(this)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    getnDof = this%ptr%nDof
  end function getnDof
  
  type(PointPtrTYPE) function getOnePoint(this, i)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getOnePoint = this%ptr%point(i)
  end function getOnePoint
  
  function getAllPoints(this)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    type(PointPtrTYPE), dimension(this%ptr%nPoint) :: getAllPoints
    getAllPoints(1:this%ptr%nPoint) = this%ptr%point(1:this%ptr%nPoint)
  end function getAllPoints
  
  integer(ikind) function getPointID(this, i)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getPointID = this%ptr%point(i)%getID()
  end function getPointID
  
  function getMaterial(this)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    type(MaterialPtrTYPE) :: getMaterial
    getMaterial = this%ptr%material
  end function getMaterial

  function getGeometry(this)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    type(GeometryPtrTYPE) :: getGeometry
    getGeometry = this%ptr%geometry
  end function getGeometry

  function getIntegrator(this)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    type(IntegratorPtrTYPE) :: getIntegrator
    getIntegrator = this%ptr%integrator
  end function getIntegrator

  subroutine setID(this, id)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    this%ptr%id = id
  end subroutine setID

  subroutine setnPoint(this, nPoint)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nPoint
    this%ptr%nPoint = nPoint
  end subroutine setnPoint

  subroutine setnDof(this, nDof)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nDof
    this%ptr%nDof = nDof
  end subroutine setnDof

  subroutine setPoint(this, i, point)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    class(PointTYPE), target, intent(in) :: point
    this%ptr%point(i)%ptr => point
  end subroutine setPoint

  subroutine setIntegrator(this, integrator)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    type(IntegratorTYPE), target, intent(in) :: integrator
    this%ptr%integrator%ptr => integrator
  end subroutine setIntegrator

  subroutine setArea(this)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    call this%ptr%setArea()
  end subroutine setArea

  real(rkind) function getArea(this)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    getArea = this%ptr%area
  end function getArea

  function jacobian(this, x, y)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(2,2) :: jacobian
    jacobian = this%ptr%jacobian(x, y)
  end function jacobian

  real(rkind) function jacobianDet(this, jacobian)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    real(rkind), dimension(2,2) :: jacobian
    jacobianDet = this%ptr%jacobianDet(jacobian)
  end function jacobianDet

  function shapeFunc(this, x, y)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(this%ptr%nPoint*this%ptr%nDof) :: shapeFunc
    shapeFunc = this%ptr%shapeFunc(x, y)
  end function shapeFunc

  function dShapeFunc(this, x, y)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(2, this%ptr%nPoint*this%ptr%nDof) :: dShapeFunc
    dShapeFunc = this%ptr%dShapeFunc(x,y)
  end function dShapeFunc

  function getStiffness(this)
    implicit none
    class(Element2DPtrTYPE), intent(inout) :: this
    real(rkind), dimension(this%ptr%nPoint*this%ptr%nDof,this%ptr%nPoint*this%ptr%nDof) :: getStiffness
    getStiffness = this%ptr%getStiffness()
  end function getStiffness

end module Element2DPtrMOD
