module StructLinQuadElementMOD
  use tools
  use PointMOD
  use PointPtrMOD
  use MaterialPtrMOD
  use IntegratorMOD
  use Element2DMOD
  use QuadElementMOD
  implicit none
  private
  public :: StructLinQuadElementTYPE, structLinearQuadElement
  type, extends(QuadElementTYPE) :: StructLinQuadElementTYPE
     private
   contains
     procedure, public :: init
     procedure, public :: shapeFunc
     procedure, public :: dShapeFunc
  end type StructLinQuadElementTYPE
  
  interface structLinearQuadElement
     procedure constructor
  end interface structLinearQuadElement
  
contains
  
  type(StructLinQuadElementTYPE) function constructor(id, material, point, integrator)
    implicit none
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call constructor%init(id, material, point, integrator)
  end function constructor

  subroutine init(this, id, material, point, integrator)
    implicit none
    class(StructLinQuadElementTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call this%setID(id)
    call this%setnPoint(4)
    call this%setnDof(2)
    allocate(this%point(this%nPoint))
    this%material = material
    this%point(1:this%nPoint) = point(1:this%nPoint)
    call this%setIntegrator(integrator)
    call this%setArea()
  end subroutine init
  
  function shapeFunc(this, x, y)
    implicit none
    class(StructLinQuadElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(this%nPoint*this%nDof) :: shapeFunc
    shapeFunc(1) = (1.d0/4.d0)*(1-x)*(1-y)
    shapeFunc(2) = (1.d0/4.d0)*(1-x)*(1-y)
    shapeFunc(3) = (1.d0/4.d0)*(1+x)*(1-y)
    shapeFunc(4) = (1.d0/4.d0)*(1+x)*(1-y)
    shapeFunc(5) = (1.d0/4.d0)*(1+x)*(1+y)
    shapeFunc(6) = (1.d0/4.d0)*(1+x)*(1+y)
    shapeFunc(7) = (1.d0/4.d0)*(1-x)*(1+y)
    shapeFunc(8) = (1.d0/4.d0)*(1-x)*(1+y)
  end function shapeFunc
  
  function dShapeFunc(this, x, y)
    implicit none
    class(StructLinQuadElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(2, this%nPoint*this%nDof) :: dShapeFunc
    dShapeFunc(1,1) = -(1.d0/4.d0)*(1-y)
    dShapeFunc(1,2) = -(1.d0/4.d0)*(1-y)
    dShapeFunc(2,1) = -(1.d0/4.d0)*(1-x)
    dShapeFunc(2,2) = -(1.d0/4.d0)*(1-x)
    dShapeFunc(1,3) = (1.d0/4.d0)*(1-y)
    dShapeFunc(1,4) = (1.d0/4.d0)*(1-y)
    dShapeFunc(2,3) = -(1.d0/4.d0)*(1+x)
    dShapeFunc(2,4) = -(1.d0/4.d0)*(1+x)
    dShapeFunc(1,5) = (1.d0/4.d0)*(1+y)
    dShapeFunc(1,6) = (1.d0/4.d0)*(1+y)
    dShapeFunc(2,5) = (1.d0/4.d0)*(1+x)
    dShapeFunc(2,6) = (1.d0/4.d0)*(1+x)
    dShapeFunc(1,7) = -(1.d0/4.d0)*(1+y)
    dShapeFunc(1,8) = -(1.d0/4.d0)*(1+y)
    dShapeFunc(2,7) = (1.d0/4.d0)*(1-x)
    dShapeFunc(2,8) = (1.d0/4.d0)*(1-x)
  end function dShapeFunc
  
end module StructLinQuadElementMOD
