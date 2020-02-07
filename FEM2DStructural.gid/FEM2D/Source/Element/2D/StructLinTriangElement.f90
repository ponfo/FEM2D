module StructLinTriangElementMOD
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
  public :: StructLinTriangElementTYPE, structLinearTriangElement
  type, extends(TriangElementTYPE) :: StructLinTriangElementTYPE
     private
   contains
     procedure, public :: init
     procedure, public :: shapeFunc
     procedure, public :: dShapeFunc
  end type StructLinTriangElementTYPE
  
  interface structLinearTriangElement
     procedure constructor
  end interface structLinearTriangElement
  
contains
  
  type(StructLinTriangElementTYPE) function constructor(id, material, point, integrator)
    implicit none
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call constructor%init(id, material, point, integrator)
  end function constructor

  subroutine init(this, id, material, point, integrator)
    implicit none
    class(StructLinTriangElementTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call this%setID(id)
    call this%setnPoint(3)
    call this%setnDof(2)
    allocate(this%point(this%nPoint))
    this%material = material
    this%point(1:this%nPoint) = point(1:this%nPoint)
    call this%setIntegrator(integrator)
    call this%setArea()
  end subroutine init
  
  function shapeFunc(this, x, y)
    implicit none
    class(StructLinTriangElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(this%nPoint*this%nDof) :: shapeFunc
    shapeFunc(1) = 1-x-y
    shapeFunc(2) = 1-x-y
    shapeFunc(3) = x
    shapeFunc(4) = x
    shapefunc(5) = y
    shapefunc(6) = y
  end function shapeFunc
  
  function dShapeFunc(this, x, y)
    implicit none
    class(StructLinTriangElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(2, this%nPoint*this%nDof) :: dShapeFunc
    dShapeFunc(1,1) = -1
    dShapeFunc(1,2) = -1
    dShapeFunc(1,3) = 1
    dShapeFunc(1,4) = 1
    dShapeFunc(1,5) = 0
    dShapeFunc(1,6) = 0
    dShapeFunc(2,1) = -1
    dShapeFunc(2,2) = -1
    dShapeFunc(2,3) = 0
    dShapeFunc(2,4) = 0
    dShapeFunc(2,5) = 1
    dShapeFunc(2,6) = 1
  end function dShapeFunc
  
end module StructLinTriangElementMOD
