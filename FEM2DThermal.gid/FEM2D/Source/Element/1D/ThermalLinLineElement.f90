module ThermalLinLineElementMOD
  use tools
  use Element1DMOD
  use ThermalElement1DMOD
  use MaterialPtrMOD
  use PointPtrMOD
  use IntegratorMOD
  implicit none
  private
  public :: ThermalLinLineElementTYPE, thermalLinearLineElement
  type, extends(ThermalElement1DTYPE) :: ThermalLinLineElementTYPE
   contains
     procedure, public :: init
     procedure, public :: shapeFunc
     procedure, public :: dShapeFunc
  end type ThermalLinLineElementTYPE
  
  interface thermalLinearLineElement
     procedure constructor
  end interface thermalLinearLineElement

contains
  type(ThermalLinLineElementTYPE) function constructor(id, material, point, integrator)
    implicit none
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call constructor%init(id, material, point, integrator)
  end function constructor

  subroutine init(this, id, material, point, integrator)
    implicit none
    class(ThermalLinLineElementTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call this%setID(id)
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
    class(ThermalLinLineElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: u
    real(rkind), dimension(this%nPoint*this%nDof) :: shapeFunc
    shapeFunc(1) = .5 - .5*u
    shapeFunc(2) = .5 + .5*u
  end function shapeFunc

  function dShapeFunc(this, u)
    implicit none
    class(ThermalLinLineElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: u
    real(rkind), dimension(this%nPoint*this%nDof) :: dShapeFunc
    dShapeFunc(1) = -.5
    dShapeFunc(2) = .5
  end function dShapeFunc
    

end module ThermalLinLineElementMOD

    
