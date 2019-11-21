module ThQuadLineElementMOD
  use tools
  use Element1DMOD
  use ThElement1DMOD
  use MaterialPtrMOD
  use PointPtrMOD
  use IntegratorMOD
  implicit none
  private
  public :: ThQuadLineElementTYPE, thermalQuadLineElement
  type, extends(ThElement1DTYPE) :: ThQuadLineElementTYPE
   contains
     procedure, public :: init
     procedure, public :: shapeFunc
     procedure, public :: dShapeFunc
  end type ThQuadLineElementTYPE
  
  interface thermalQuadLineElement
     procedure constructor
  end interface thermalQuadLineElement

contains
  type(ThQuadLineElementTYPE) function constructor(material, point, integrator)
    implicit none
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call constructor%init( material, point, integrator)
  end function constructor

  subroutine init(this, material, point, integrator)
    implicit none
    class(ThQuadLineElementTYPE), intent(inout) :: this
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call this%setnPoint(3)
    call this%setnDof(1)
    allocate(this%point(this%nPoint))
    this%material = material
    this%point(1:this%nPoint) = point(1:this%nPoint)
    call this%setIntegrator(integrator)
    call this%setLenght()
  end subroutine init

  function shapeFunc(this, u)
    implicit none
    class(ThQuadLineElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: u
    real(rkind), dimension(this%nPoint*this%nDof) :: shapeFunc
    shapeFunc(1) = 0.5*u*(u-1)
    shapeFunc(2) = (1+u)*(1-u)
    shapeFunc(3) = 0.5*u*(1+u)
  end function shapeFunc

  function dShapeFunc(this, u)
    implicit none
    class(ThQuadLineElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: u
    real(rkind), dimension(this%nPoint*this%nDof) :: dShapeFunc
    dShapeFunc(1) = u-0.5
    dShapeFunc(2) = -2*u
    dShapeFunc(3) = u+0.5
  end function dShapeFunc

end module ThQuadLineElementMOD
