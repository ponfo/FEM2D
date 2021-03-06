module ThermalQuadTriangElementMOD
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
  public :: ThermalQuadTriangElementTYPE, thermalQuadTriangElement
  type, extends(TriangElementTYPE) :: ThermalQuadTriangElementTYPE
     private
   contains
     procedure, public :: init
     procedure, public :: shapeFunc
     procedure, public :: dShapeFunc
  end type ThermalQuadTriangElementTYPE
  
  interface thermalQuadTriangElement
     procedure constructor
  end interface thermalQuadTriangElement
  
contains
  
  type(ThermalQuadTriangElementTYPE) function constructor(id, material, point, integrator)
    implicit none
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call constructor%init(id, material, point, integrator)
  end function constructor

  subroutine init(this, id, material, point, integrator)
    implicit none
    class(ThermalQuadTriangElementTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call this%setID(id)
    call this%setnPoint(6)
    call this%setnDof(1)
    allocate(this%point(this%nPoint))
    this%material = material
    this%point(1:this%nPoint) = point(1:this%nPoint)
    call this%setIntegrator(integrator)
    call this%setArea()
  end subroutine init
  
  function shapeFunc(this, x, y)
    implicit none
    class(ThermalQuadTriangElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(this%nPoint*this%nDof) :: shapeFunc
    !Esquinas
    shapeFunc(1) = (1-x-y)*(1-2*x-2*y)
    shapeFunc(2) = x*(2*x-1)
    shapeFunc(3) = y*(2*y-1)
    !Lados
    shapeFunc(4) = 4*x*(1-x-y)
    shapeFunc(5) = 4*x*y
    shapeFunc(6) = 4*y*(1-x-y)
  end function shapeFunc
  
  function dShapeFunc(this, x, y)
    implicit none
    class(ThermalQuadTriangElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(2, this%nPoint*this%nDof) :: dShapeFunc
    !Esquinas
    dShapeFunc(1,1) = 4*x+4*y-3
    dShapeFunc(2,1) = 4*y+4*x-3
    dShapeFunc(1,2) = 4*x-1
    dShapeFunc(2,2) = 0.d0
    dShapeFunc(1,3) = 0.d0
    dShapeFunc(2,3) = 4*y-1
    !Lados
    dShapeFunc(1,4) = -8*x-4*y+4
    dShapeFunc(2,4) = -4*x
    dShapeFunc(1,5) = 4*y
    dShapeFunc(2,5) = 4*x
    dShapeFunc(1,6) = -4*y
    dShapeFunc(2,6) = -8*y-4*x+4
  end function dShapeFunc
  
end module ThermalQuadTriangElementMOD
