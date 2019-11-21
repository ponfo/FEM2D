module ThLinLineElementMOD
  use tools
  use Element1DMOD
  use ThElement1DMOD
  use MaterialPtrMOD
  use PointPtrMOD
  use IntegratorMOD
  implicit none
  private
  public :: ThLinLineElementTYPE, thermalLinearLineElement
  type, extends(ThElement1DTYPE) :: ThLinLineElementTYPE
   contains
     procedure, public :: init
     procedure, public :: shapeFunc
     procedure, public :: dShapeFunc
  end type ThLinLineElementTYPE
  
  interface thermalLinearLineElement
     procedure constructor
  end interface thermalLinearLineElement

contains
  type(ThLinLineElementTYPE) function constructor(material, point, integrator)
    implicit none
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call constructor%init(material, point, integrator)
  end function constructor

  subroutine init(this, material, point, integrator)
    implicit none
    class(ThLinLineElementTYPE), intent(inout) :: this
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call this%setnPoint(2)
    call this%setnDof(1)
    allocate(this%point(this%nPoint))
    this%material = material
    this%point(1:this%nPoint) = point(1:this%nPoint)
    call this%setIntegrator(integrator)
    call this%setLenght()
  end subroutine init

  function shapeFunc(this, u)
    implicit none
    class(ThLinLineElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: u
    real(rkind), dimension(this%nPoint*this%nDof) :: shapeFunc
    shapeFunc(1) = .5 - .5*u
    shapeFunc(2) = .5 + .5*u
  end function shapeFunc

  function dShapeFunc(this, u)
    implicit none
    class(ThLinLineElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: u
    real(rkind), dimension(this%nPoint*this%nDof) :: dShapeFunc
    dShapeFunc(1) = -.5
    dShapeFunc(2) = .5
  end function dShapeFunc
    

end module ThLinLineElementMOD

    