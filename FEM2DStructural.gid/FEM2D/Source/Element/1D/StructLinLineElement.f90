module StructLinLineElementMOD
  use tools
  use Element1DMOD
  use StructElement1DMOD
  use MaterialPtrMOD
  use PointPtrMOD
  use IntegratorMOD
  implicit none
  private
  public :: StructLinLineElementTYPE, structLinearLineElement
  type, extends(StructElement1DTYPE) :: StructLinLineElementTYPE
   contains
     procedure, public :: init
     procedure, public :: shapeFunc
     procedure, public :: dShapeFunc
  end type StructLinLineElementTYPE
  
  interface structLinearLineElement
     procedure constructor
  end interface structLinearLineElement

contains
  type(StructLinLineElementTYPE) function constructor(id, material, point, integrator)
    implicit none
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call constructor%init(id, material, point, integrator)
  end function constructor

  subroutine init(this, id, material, point, integrator)
    implicit none
    class(StructLinLineElementTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call this%setID(id)
    call this%setnPoint(2)
    call this%setnDof(2)
    allocate(this%point(this%nPoint))
    this%material = material
    this%point(1:this%nPoint) = point(1:this%nPoint)
    call this%setIntegrator(integrator)
    call this%setLenght()
  end subroutine init

  function shapeFunc(this, u)
    implicit none
    class(StructLinLineElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: u
    real(rkind), dimension(this%nPoint*this%nDof) :: shapeFunc
    shapeFunc(1) = .5 - .5*u
    shapeFunc(2) = .5 - .5*u
    shapeFunc(3) = .5 + .5*u
    shapeFunc(4) = .5 + .5*u
  end function shapeFunc

  function dShapeFunc(this, u)
    implicit none
    class(StructLinLineElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: u
    real(rkind), dimension(this%nPoint*this%nDof) :: dShapeFunc
    dShapeFunc(1) = -.5
    dShapeFunc(2) = -.5
    dShapeFunc(3) = .5
    dShapeFunc(4) = .5
  end function dShapeFunc

end module StructLinLineElementMOD
