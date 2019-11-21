module ThQuadQuadElementMOD
  use tools
  use PointMOD
  use PointPtrMOD
  use MaterialPtrMOD
  use IntegratorMOD
  use Element2DMOD
  use QuadElementMOD
  implicit none
  private
  public :: ThQuadQuadElementTYPE, thermalQuadQuadElement
  type, extends(QuadElementTYPE) :: ThQuadQuadElementTYPE
     private
   contains
     procedure, public :: init
     procedure, public :: shapeFunc
     procedure, public :: dShapeFunc
  end type ThQuadQuadElementTYPE
  
  interface thermalQuadQuadElement
     procedure constructor
  end interface thermalQuadQuadElement
  
contains
  
  type(ThQuadQuadElementTYPE) function constructor(material, point, integrator)
    implicit none
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call constructor%init(material, point, integrator)
  end function constructor

  subroutine init(this, material, point, integrator)
    implicit none
    class(ThQuadQuadElementTYPE), intent(inout) :: this
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call this%setnPoint(8)
    call this%setnDof(1)
    allocate(this%point(this%nPoint))
    this%material = material
    this%point(1:this%nPoint) = point(1:this%nPoint)
    call this%setIntegrator(integrator)
    call this%setArea()
  end subroutine init
  
  function shapeFunc(this, x, y)
    implicit none
    class(ThQuadQuadElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(this%nPoint*this%nDof) :: shapeFunc
    !Esquinas
    shapeFunc(1) = (1.d0/4.d0)*(1-y)*(1-x)*(-1-x-y)
    shapeFunc(2) = (1.d0/4.d0)*(1-y)*(1+x)*(-1-y+x)
    shapeFunc(3) = (1.d0/4.d0)*(1+y)*(1+x)*(-1+x+y)
    shapeFunc(4) = (1.d0/4.d0)*(1+y)*(1-x)*(-1-x+y)
    !Lados
    shapefunc(5) = (1.d0/2.d0)*(1-y)*(1-x*x)
    shapeFunc(6) = (1.d0/2.d0)*(1+x)*(1-y*y)
    shapeFunc(7) = (1.d0/2.d0)*(1+y)*(1-x*x)
    shapefunc(8) = (1.d0/2.d0)*(1-x)*(1-y*y)
  end function shapeFunc
  
  function dShapeFunc(this, x, y)
    implicit none
    class(ThQuadQuadElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(2, this%nPoint*this%nDof) :: dShapeFunc
    !Esquinas
    dShapeFunc(1,1) = x/4.d0+y/4.d0-(x*y)/4.d0-(y/4.d0-1.d0/4.d0)*(x+y+1)-1.d0/4.d0
    dShapeFunc(2,1) = x/4.d0+y/4.d0-(x*y)/4.d0-(x/4.d0-1.d0/4.d0)*(x+y+1)-1.d0/4.d0
    dShapeFunc(1,2) = (y/4.d0-1.d0/4.d0)*(y-x+1)-(y/4.d0-1.d0/4.d0)*(x+1)
    dShapeFunc(2,2) = (y/4.d0-1.d0/4.d0)*(x+1)+((x+1)*(y-x+1))/4.d0
    dShapeFunc(1,3) = (y/4.d0+1.d0/4.d0)*(x+1)+(y/4.d0+1.d0/4.d0)*(x+y-1)
    dShapeFunc(2,3) = (y/4.d0+1.d0/4.d0)*(x+1)+((x+1)*(x+y-1))/4.d0
    dShapeFunc(1,4) = (y/4.d0+1.d0/4.d0)*(x-y+1)+(y/4.d0+1/4.d0)*(x-1)
    dShapeFunc(2,4) = ((x-1)*(x-y+1))/4.d0-(y/4.d0+1.d0/4.d0)*(x-1)
    !Lados
    dShapeFunc(1,5) = 2*x*(y/2.d0-1.d0/2.d0)
    dShapeFunc(2,5) = (x**2)/2.d0-1.d0/2.d0
    dShapeFunc(1,6) = 1.d0/2.d0-(y**2)/2.d0
    dShapeFunc(2,6) = -2*y*(x/2.d0+1.d0/2.d0)
    dShapeFunc(1,7) = -2*x*(y/2.d0+1.d0/2.d0)
    dShapeFunc(2,7) = 1.d0/2.d0-(x**2)/2.d0
    dShapeFunc(1,8) = (y**2)/2.d0-1.d0/2.d0
    dShapeFunc(2,8) = 2*y*(x/2.d0-1.d0/2.d0)
  end function dShapeFunc
  
end module ThQuadQuadElementMOD