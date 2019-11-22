module StructQuadLineElementMOD
  use tools
  use Element1DMOD
  use StructElement1DMOD
  use MaterialPtrMOD
  use PointPtrMOD
  use IntegratorMOD
  implicit none
  private
  public :: StructQuadLineElementTYPE, structQuadLineElement
  type, extends(StructElement1DTYPE) :: StructQuadLineElementTYPE
   contains
     procedure, public :: init
     procedure, public :: shapeFunc
     procedure, public :: dShapeFunc
  end type StructQuadLineElementTYPE
  
  interface structQuadLineElement
     procedure constructor
  end interface structQuadLineElement

contains
  type(StructQuadLineElementTYPE) function constructor(id, material, point, integrator)
    implicit none
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call constructor%init(id, material, point, integrator)
  end function constructor

  subroutine init(this, id, material, point, integrator)
    implicit none
    class(StructQuadLineElementTYPE), intent(inout) :: this
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
    call this%setLenght()
  end subroutine init

  function shapeFunc(this, u)
    implicit none
    class(StructQuadLineElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: u
    real(rkind), dimension(this%nPoint*this%nDof) :: shapeFunc
    shapeFunc(1) = 0.5*u*(u-1)
    shapeFunc(2) = 0.5*u*(u-1)
    shapeFunc(3) = (1+u)*(1-u)
    shapeFunc(4) = (1+u)*(1-u)
    shapeFunc(5) = 0.5*u*(1+u)
    shapeFunc(6) = 0.5*u*(1+u)
  end function shapeFunc

  function dShapeFunc(this, u)
    implicit none
    class(StructQuadLineElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: u
    real(rkind), dimension(this%nPoint*this%nDof) :: dShapeFunc
    dShapeFunc(1) = u-0.5
    dShapeFunc(2) = u-0.5
    dShapeFunc(3) = -2*u
    dShapeFunc(4) = -2*u
    dShapeFunc(5) = u+0.5
    dShapeFunc(6) = u+0.5
  end function dShapeFunc

end module StructQuadLineElementMOD
