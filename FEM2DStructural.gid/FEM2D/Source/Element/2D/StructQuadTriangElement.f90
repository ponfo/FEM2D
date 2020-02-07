module StructQuadTriangElementMOD
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
  public :: StructQuadTriangElementTYPE, structQuadTriangElement
  type, extends(TriangElementTYPE) :: StructQuadTriangElementTYPE
     private
   contains
     procedure, public :: init
     procedure, public :: shapeFunc
     procedure, public :: dShapeFunc
  end type StructQuadTriangElementTYPE
  
  interface structQuadTriangElement
     procedure constructor
  end interface structQuadTriangElement
  
contains
  
  type(StructQuadTriangElementTYPE) function constructor(id, material, point, integrator)
    implicit none
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call constructor%init(id, material, point, integrator)
  end function constructor

  subroutine init(this, id, material, point, integrator)
    implicit none
    class(StructQuadTriangElementTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call this%setID(id)
    call this%setnPoint(6)
    call this%setnDof(2)
    allocate(this%point(this%nPoint))
    this%material = material
    this%point(1:this%nPoint) = point(1:this%nPoint)
    call this%setIntegrator(integrator)
    call this%setArea()
  end subroutine init
  
  function shapeFunc(this, x, y)
    implicit none
    class(StructQuadTriangElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(this%nPoint*this%nDof) :: shapeFunc
    !Esquinas
    shapeFunc(1) = (1-x-y)*(1-2*x-2*y)
    shapeFunc(2) = (1-x-y)*(1-2*x-2*y)
    shapeFunc(3) = x*(2*x-1)
    shapeFunc(4) = x*(2*x-1)
    shapeFunc(5) = y*(2*y-1)
    shapeFunc(6) = y*(2*y-1)
    !Lados
    shapeFunc(7) = 4*x*(1-x-y)
    shapeFunc(8) = 4*x*(1-x-y)
    shapeFunc(9) = 4*x*y
    shapeFunc(10) = 4*x*y
    shapeFunc(11) = 4*y*(1-x-y)
    shapeFunc(12) = 4*y*(1-x-y)
  end function shapeFunc
  
  function dShapeFunc(this, x, y)
    implicit none
    class(StructQuadTriangElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(2, this%nPoint*this%nDof) :: dShapeFunc
    !Esquinas
    dShapeFunc(1,1) = 4*x+4*y-3
    dShapeFunc(1,2) = 4*x+4*y-3
    dShapeFunc(2,1) = 4*y+4*x-3
    dShapeFunc(2,2) = 4*y+4*x-3
    dShapeFunc(1,3) = 4*x-1
    dShapeFunc(1,4) = 4*x-1
    dShapeFunc(2,3) = 0.d0
    dShapeFunc(2,4) = 0.d0
    dShapeFunc(1,5) = 0.d0
    dShapeFunc(1,6) = 0.d0
    dShapeFunc(2,5) = 4*y-1
    dShapeFunc(2,6) = 4*y-1
    !Lados
    dShapeFunc(1,7) = -8*x-4*y+4
    dShapeFunc(1,8) = -8*x-4*y+4
    dShapeFunc(2,7) = -4*x
    dShapeFunc(2,8) = -4*x
    dShapeFunc(1,9) = 4*y
    dShapeFunc(1,10) = 4*y
    dShapeFunc(2,9) = 4*x
    dShapeFunc(2,10) = 4*x
    dShapeFunc(1,11) = -4*y
    dShapeFunc(1,12) = -4*y
    dShapeFunc(2,11) = -8*y-4*x+4
    dShapeFunc(2,12) = -8*y-4*x+4
  end function dShapeFunc
  
end module StructQuadTriangElementMOD
