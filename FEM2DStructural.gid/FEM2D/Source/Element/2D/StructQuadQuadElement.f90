module StructQuadQuadElementMOD
  use tools
  use PointMOD
  use PointPtrMOD
  use MaterialPtrMOD
  use IntegratorMOD
  use Element2DMOD
  use QuadElementMOD
  implicit none
  private
  public :: StructQuadQuadElementTYPE, structQuadQuadElement
  type, extends(QuadElementTYPE) :: StructQuadQuadElementTYPE
     private
   contains
     procedure, public :: init
     procedure, public :: shapeFunc
     procedure, public :: dShapeFunc
  end type StructQuadQuadElementTYPE
  
  interface structQuadQuadElement
     procedure constructor
  end interface structQuadQuadElement
  
contains
  
  type(StructQuadQuadElementTYPE) function constructor(id, material, point, integrator)
    implicit none
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call constructor%init(id, material, point, integrator)
  end function constructor

  subroutine init(this, id, material, point, integrator)
    implicit none
    class(StructQuadQuadElementTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(in) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    type(IntegratorTYPE), intent(in) :: integrator
    call this%setID(id)
    call this%setnPoint(8)
    call this%setnDof(2)
    allocate(this%point(this%nPoint))
    this%material = material
    this%point(1:this%nPoint) = point(1:this%nPoint)
    call this%setIntegrator(integrator)
    call this%setArea()
  end subroutine init
  
  function shapeFunc(this, x, y)
    implicit none
    class(StructQuadQuadElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(this%nPoint*this%nDof) :: shapeFunc
    !Esquinas
    shapeFunc(1) = (1.d0/4.d0)*(1-y)*(1-x)*(-1-x-y)
    shapeFunc(2) = (1.d0/4.d0)*(1-y)*(1-x)*(-1-x-y)
    shapeFunc(3) = (1.d0/4.d0)*(1-y)*(1+x)*(-1-y+x)
    shapeFunc(4) = (1.d0/4.d0)*(1-y)*(1+x)*(-1-y+x)
    shapeFunc(5) = (1.d0/4.d0)*(1+y)*(1+x)*(-1+x+y)
    shapeFunc(6) = (1.d0/4.d0)*(1+y)*(1+x)*(-1+x+y)
    shapeFunc(7) = (1.d0/4.d0)*(1+y)*(1-x)*(-1-x+y)
    shapeFunc(8) = (1.d0/4.d0)*(1+y)*(1-x)*(-1-x+y)
    !Lados
    shapefunc(9) = (1.d0/2.d0)*(1-y)*(1-x*x)
    shapefunc(10) = (1.d0/2.d0)*(1-y)*(1-x*x)
    shapeFunc(11) = (1.d0/2.d0)*(1+x)*(1-y*y)
    shapeFunc(12) = (1.d0/2.d0)*(1+x)*(1-y*y)
    shapeFunc(13) = (1.d0/2.d0)*(1+y)*(1-x*x)
    shapeFunc(14) = (1.d0/2.d0)*(1+y)*(1-x*x)
    shapefunc(15) = (1.d0/2.d0)*(1-x)*(1-y*y)
    shapefunc(16) = (1.d0/2.d0)*(1-x)*(1-y*y)
  end function shapeFunc
  
  function dShapeFunc(this, x, y)
    implicit none
    class(StructQuadQuadElementTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(2, this%nPoint*this%nDof) :: dShapeFunc
    !Esquinas
    dShapeFunc(1,1) = x/4.d0+y/4.d0-(x*y)/4.d0-(y/4.d0-1.d0/4.d0)*(x+y+1)-1.d0/4.d0
    dShapeFunc(1,2) = x/4.d0+y/4.d0-(x*y)/4.d0-(y/4.d0-1.d0/4.d0)*(x+y+1)-1.d0/4.d0
    dShapeFunc(2,1) = x/4.d0+y/4.d0-(x*y)/4.d0-(x/4.d0-1.d0/4.d0)*(x+y+1)-1.d0/4.d0
    dShapeFunc(2,2) = x/4.d0+y/4.d0-(x*y)/4.d0-(x/4.d0-1.d0/4.d0)*(x+y+1)-1.d0/4.d0
    dShapeFunc(1,3) = (y/4.d0-1.d0/4.d0)*(y-x+1)-(y/4.d0-1.d0/4.d0)*(x+1)
    dShapeFunc(1,4) = (y/4.d0-1.d0/4.d0)*(y-x+1)-(y/4.d0-1.d0/4.d0)*(x+1)
    dShapeFunc(2,3) = (y/4.d0-1.d0/4.d0)*(x+1)+((x+1)*(y-x+1))/4.d0
    dShapeFunc(2,4) = (y/4.d0-1.d0/4.d0)*(x+1)+((x+1)*(y-x+1))/4.d0
    dShapeFunc(1,5) = (y/4.d0+1.d0/4.d0)*(x+1)+(y/4.d0+1.d0/4.d0)*(x+y-1)
    dShapeFunc(1,6) = (y/4.d0+1.d0/4.d0)*(x+1)+(y/4.d0+1.d0/4.d0)*(x+y-1)
    dShapeFunc(2,5) = (y/4.d0+1.d0/4.d0)*(x+1)+((x+1)*(x+y-1))/4.d0
    dShapeFunc(2,6) = (y/4.d0+1.d0/4.d0)*(x+1)+((x+1)*(x+y-1))/4.d0
    dShapeFunc(1,7) = (y/4.d0+1.d0/4.d0)*(x-y+1)+(y/4.d0+1/4.d0)*(x-1)
    dShapeFunc(1,8) = (y/4.d0+1.d0/4.d0)*(x-y+1)+(y/4.d0+1/4.d0)*(x-1)
    dShapeFunc(2,7) = ((x-1)*(x-y+1))/4.d0-(y/4.d0+1.d0/4.d0)*(x-1)
    dShapeFunc(2,8) = ((x-1)*(x-y+1))/4.d0-(y/4.d0+1.d0/4.d0)*(x-1)
    !Lados
    dShapeFunc(1,9) = 2*x*(y/2.d0-1.d0/2.d0)
    dShapeFunc(1,10) = 2*x*(y/2.d0-1.d0/2.d0)
    dShapeFunc(2,9) = (x**2)/2.d0-1.d0/2.d0
    dShapeFunc(2,10) = (x**2)/2.d0-1.d0/2.d0
    dShapeFunc(1,11) = 1.d0/2.d0-(y**2)/2.d0
    dShapeFunc(1,12) = 1.d0/2.d0-(y**2)/2.d0
    dShapeFunc(2,11) = -2*y*(x/2.d0+1.d0/2.d0)
    dShapeFunc(2,12) = -2*y*(x/2.d0+1.d0/2.d0)
    dShapeFunc(1,13) = -2*x*(y/2.d0+1.d0/2.d0)
    dShapeFunc(1,14) = -2*x*(y/2.d0+1.d0/2.d0)
    dShapeFunc(2,13) = 1.d0/2.d0-(x**2)/2.d0
    dShapeFunc(2,14) = 1.d0/2.d0-(x**2)/2.d0
    dShapeFunc(1,15) = (y**2)/2.d0-1.d0/2.d0
    dShapeFunc(1,16) = (y**2)/2.d0-1.d0/2.d0
    dShapeFunc(2,15) = 2*y*(x/2.d0-1.d0/2.d0)
    dShapeFunc(2,16) = 2*y*(x/2.d0-1.d0/2.d0)
  end function dShapeFunc
  
end module StructQuadQuadElementMOD
