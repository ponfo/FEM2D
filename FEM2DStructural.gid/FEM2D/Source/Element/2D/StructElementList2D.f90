module StructElementList2DMOD
  use tools
  use DebuggerMOD

  use Element2DMOD
  use Element2DPtrMOD
  use StructLinTriangElementMOD
  use StructLinQuadElementMOD
  use StructQuadTriangElementMOD
  use StructQuadQuadElementMOD

  use PointPtrMOD

  use MaterialPtrMOD

  use IntegratorMOD
  use IntegratorPtrMOD

  implicit none
  private
  public :: StructElementList2DTYPE, structElementList2D
  type :: StructElementList2DTYPE
     private
     integer(ikind)                                    :: gaussOrder
     integer(ikind)                                    :: nElem
     type(Element2DPtrTYPE), dimension(:), allocatable :: element
     type(IntegratorTYPE)                              :: triangIntegrator
     type(IntegratorTYPE)                              :: quadIntegrator
   contains
     procedure, public  :: init
     procedure, public  :: addElement
     procedure, public  :: getGaussOrder
     procedure, public  :: getnElem
     procedure, public  :: getElement
     procedure, public  :: getTriangIntegrator
     procedure, public  :: getQuadIntegrator
     
     procedure, private :: valueTriangle
     procedure, private :: valueQuad
  end type StructElementList2DTYPE

  interface structElementList2D
     procedure constructor
  end interface structElementList2D

  procedure(addElementLin), pointer :: addElementPtr => null()

  integer(ikind)                   , save                                    :: iTriang
  integer(ikind)                   , save                                    :: iQuad
  integer(ikind)                   , save                                    :: iElem
  type(StructLinTriangElementTYPE) , save, dimension(:), allocatable, target :: linTriangElem
  type(StructLinQuadElementTYPE)   , save, dimension(:), allocatable, target :: linQuadElem
  type(StructQuadTriangElementTYPE), save, dimension(:), allocatable, target :: quadTriangElem
  type(StructQuadQuadElementTYPE)  , save, dimension(:), allocatable, target :: quadQuadElem
  
contains

  type(StructElementList2DTYPE) function constructor(isQuadratic, nTriangElem, nQuadElem, nGauss)
    implicit none
    integer(ikind), intent(in) :: isQuadratic
    integer(ikind), intent(in) :: nTriangElem
    integer(ikind), intent(in) :: nQuadElem
    integer(ikind), intent(in) :: nGauss
    call constructor%init(isQuadratic, nTriangElem, nQuadElem, nGauss)
  end function constructor
  
  subroutine init(this, isQuadratic, nTriangElem, nQuadElem, nGauss)
    implicit none
    class(StructElementList2DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: isQuadratic
    integer(ikind), intent(in) :: nTriangElem
    integer(ikind), intent(in) :: nQuadElem
    integer(ikind), intent(in) :: nGauss
    call debugLog('    Initiating ElementList2D')
    iTriang = 0
    iQuad = 0
    iElem = 0
    this%gaussOrder = nGauss
    this%nElem = nTriangElem + nQuadElem
    if(isQuadratic == 0) then
       addElementPtr => addElementLin
       if(nTriangElem > 0) then
          allocate(linTriangElem(nTriangElem))
          call debugLog('      Allocated Linear Triangs: ', size(linTriangElem))
          call linTriangElem(1)%setnPoint(3)
          call linTriangElem(1)%setnDof(2)
          this%triangIntegrator = integrator(nGauss, 'triangle')
          call this%valueTriangle(linTriangElem(1))
       end if
       if(nQuadElem > 0) then
          allocate(linQuadElem(nQuadElem))
          call debugLog('      Allocated Linear Quads: ', size(linQuadElem))
          call linQuadElem(1)%setnPoint(4)
          call linQuadElem(1)%setnDof(2)
          this%quadIntegrator = integrator(nGauss, 'quad')
          call this%valueQuad(linQuadElem(1))
       end if
    else if(isQuadratic == 1) then
       addElementPtr => addElementQuad
       if(nTriangElem > 0) then
          allocate(quadTriangElem(nTriangElem))
          call debugLog('      Allocated Quadratic Triangs: ', size(quadTriangElem))
          call quadTriangElem(1)%setnPoint(6)
          call quadTriangelem(1)%setnDof(2)
          this%triangIntegrator = integrator(nGauss, 'triangle')
          call this%valueTriangle(quadTriangElem(1))
       end if
       if(nQuadElem > 0) then
          allocate(quadQuadElem(nQuadElem))
          call debugLog('      Allocated Quadratic Quads: ', size(quadQuadElem))
          call quadQuadElem(1)%setnPoint(8)
          call quadQuadElem(1)%setnDof(2)
          this%quadIntegrator = integrator(nGauss, 'quad')
          call this%valueQuad(quadQuadElem(1))
       end if
    end if
    allocate(this%element(nTriangElem+nQuadElem))
    call debugLog('      Allocated element pointers: ', size(this%element))
  end subroutine init

  subroutine valueTriangle(this, elementInput)
    implicit none
    class(StructElementList2DTYPE), intent(inout) :: this
    class(Element2DTYPE), intent(inout) :: elementInput
    integer(ikind) :: i
    integer(ikind) :: n
    n = elementInput%getnPoint()*elementInput%getnDof()
    allocate(this%triangIntegrator%shapeFunc(this%triangIntegrator%integTerms, n))
    allocate(this%triangIntegrator%dShapeFunc(this%triangIntegrator%integTerms, 2, n))
    do i = 1, this%triangIntegrator%integTerms
       this%triangIntegrator%shapeFunc(i, 1:n) = &
            elementInput%shapefunc(this%triangIntegrator%gPoint(i,1), this%triangIntegrator%gPoint(i,2))
       this%triangIntegrator%dShapeFunc(i, 1:2, 1:n) = &
            elementInput%dShapeFunc(this%triangIntegrator%gPoint(i,1), this%triangIntegrator%gPoint(i,2))
    end do
  end subroutine valueTriangle

  subroutine valueQuad(this, elementInput)
    implicit none
    class(StructElementList2DTYPE), intent(inout) :: this
    class(Element2DTYPE), intent(inout) :: elementInput
    integer(ikind) :: i
    integer(ikind) :: n
    n = elementInput%getnPoint()*elementInput%getnDof()
    allocate(this%quadIntegrator%shapeFunc(this%quadIntegrator%integTerms, n))
    allocate(this%quadIntegrator%dShapeFunc(this%quadIntegrator%integTerms, 2, n))
    do i = 1, this%quadIntegrator%integTerms
       this%quadIntegrator%shapeFunc(i, 1:n) = &
            elementInput%shapefunc(this%quadIntegrator%gPoint(i,1), this%quadIntegrator%gPoint(i,2))
       this%quadIntegrator%dShapeFunc(i, 1:2, 1:n) = &
            elementInput%dShapeFunc(this%quadIntegrator%gPoint(i,1), this%quadIntegrator%gPoint(i,2))
    end do
  end subroutine valueQuad
  
  subroutine addElement(this, id, material, point)
    implicit none
    class(StructElementList2DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(inout) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    call addElementPtr(this, id, material, point)
  end subroutine addElement
  
  subroutine addElementLin(this, id, material, point)
    class(StructElementList2DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(inout) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    integer(ikind) :: i
    iElem = iElem + 1
    if(size(point) == 3) then
       iTriang = iTriang + 1
       linTriangElem(iTriang) = &
            structLinearTriangElement(id, material, point, this%triangIntegrator)
       this%element(iElem)%ptr => linTriangElem(iTriang)
    else if(size(point) == 4) then
       iQuad = iQuad + 1
       linQuadElem(iQuad) = &
            structLinearQuadElement(id, material, point, this%quadIntegrator)
       this%element(iElem)%ptr => linQuadElem(iQuad)
    end if
  end subroutine addElementLin

  subroutine addElementQuad(this, id, material, point)
    class(StructElementList2DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(inout) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    integer(ikind) :: i
    iElem = iElem + 1
    if(size(point) == 6) then
       iTriang = iTriang + 1
       quadTriangElem(iTriang) = &
            structQuadTriangElement(id, material, point, this%triangIntegrator)
       this%element(iElem)%ptr => quadTriangElem(iTriang)
    else if(size(point) == 8) then
       iQuad = iQuad + 1
       quadQuadElem(iQuad) = &
            structQuadQuadElement(id, material, point, this%quadIntegrator)
       this%element(iElem)%ptr => quadQuadElem(iQuad)
    end if
  end subroutine addElementQuad
    
  integer(ikind) function getGaussOrder(this)
    implicit none
    class(StructElementList2DTYPE), intent(inout) :: this
    getGaussOrder = this%gaussOrder
  end function getGaussOrder

  integer(ikind) function getnElem(this)
    implicit none
    class(StructElementList2DTYPE), intent(inout) :: this
    getnElem = this%nElem
  end function getnElem

  type(Element2DPtrTYPE) function getElement(this, i)
    implicit none
    class(StructElementList2DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getElement = this%element(i)
  end function getElement

  type(IntegratorPtrTYPE) function getTriangIntegrator(this)
    implicit none
    class(StructElementList2DTYPE), intent(inout) :: this
    call getTriangIntegrator%allocate(this%triangIntegrator)
  end function getTriangIntegrator

  type(IntegratorPtrTYPE) function getQuadIntegrator(this)
    implicit none
    class(StructElementList2DTYPE), intent(inout) :: this
    call getQuadIntegrator%allocate(this%quadIntegrator)
  end function getQuadIntegrator

end module StructElementList2DMOD
