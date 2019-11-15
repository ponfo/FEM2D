module ThLinTriangElementMOD
  use tools
  use PointMOD
  use PointPtrMOD
  use MaterialPtrMOD
  use PointPtrMOD
  use IntegratorMOD
  use Element2DMOD
  use TriangElementMOD
  implicit none
  private
  public :: ThLinTriangElementTYPE, thermalLinearTriangElement
  type, extends(TriangElementTYPE) :: ThLinTriangElementTYPE
     private
   contains
     procedure, public :: init
     procedure, public :: shapeFunc
     procedure, public :: dShapeFunc
  end type ThLinTriangElementTYPE
  
  interface thermalLinearTriangElement
     procedure constructor
  end interface thermalLinearTriangElement
  
contains
  
  type(ThLinTriangElementTYPE) function constructor(material, point, integrator)
    implicit none
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call constructor%init(material, point, integrator)
  end function constructor

  subroutine init(this, material, point, integrator)
    implicit none
    class(ThLinTriangElementTYPE), intent(inout) :: this
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call this%setnPoint(3)
    call this%setnDof(1)
    allocate(this%point(this%nPoint))
    this%material = material
    this%point(1:this%nPoint) = point(1:this%nPoint)
    call this%setIntegrator(integrator)
    call this%setArea()
  end subroutine init
  
  function shapeFunc(this, x, y)
    implicit none
    class(ThLinTriangElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(this%nPoint*this%nDof) :: shapeFunc
    shapeFunc(1) = 1-x-y
    shapeFunc(2) = x
    shapefunc(3) = y
  end function shapeFunc
  
  function dShapeFunc(this, x, y)
    implicit none
    class(ThLinTriangElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(2, this%nPoint*this%nDof) :: dShapeFunc
    dShapeFunc(1,1) = -1
    dShapeFunc(1,2) = 1
    dShapeFunc(1,3) = 0
    dShapeFunc(2,1) = -1
    dShapeFunc(2,2) = 0
    dShapeFunc(2,3) = 1
  end function dShapeFunc
  
end module ThLinTriangElementMOD
