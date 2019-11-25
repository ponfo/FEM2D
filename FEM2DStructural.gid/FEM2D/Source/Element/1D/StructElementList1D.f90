module StructElementList1DMOD
  use tools
  use DebuggerMOD

  use Element1DMOD
  use Element1DPtrMOD
  use StructLinLineElementMOD
  use StructQuadLineElementMOD

  use PointPtrMOD

  use MaterialPtrMOD

  use IntegratorMOD
  use IntegratorPtrMOD

  implicit none
  private
  public :: StructElementList1DTYPE, structElementList1D
  type :: StructElementList1DTYPE
     private
     integer(ikind)                                    :: gaussOrder
     integer(ikind)                                    :: nElem
     type(Element1DPtrTYPE), dimension(:), allocatable :: element
     type(IntegratorTYPE)                              :: integrator
   contains
     procedure, public  :: init
     procedure, public  :: addElement
     procedure, public  :: getGaussOrder
     procedure, public  :: getnElem
     procedure, public  :: getElement
     procedure, public  :: getIntegrator
     procedure, private :: valueElement
  end type StructElementList1DTYPE

  interface structElementList1D
     procedure constructor
  end interface structElementList1D

  procedure(addElementLin), pointer :: addElementPtr => null()

  integer(ikind)                 , save                                    :: iElem
  type(StructLinLineElementTYPE) , save, dimension(:), allocatable, target :: linLineElem
  type(StructQuadLineElementTYPE), save, dimension(:), allocatable, target :: quadLineElem
  
contains

  type(StructElementList1DTYPE) function constructor(isQuadratic, nElem, nGauss)
    implicit none
    integer(ikind), intent(in) :: isQuadratic
    integer(ikind), intent(in) :: nElem
    integer(ikind), intent(in) :: nGauss
    call constructor%init(isQuadratic, nElem, nGauss)
  end function constructor
  
  subroutine init(this, isQuadratic, nElem, nGauss)
    implicit none
    class(StructElementList1DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: isQuadratic
    integer(ikind), intent(in) :: nElem
    integer(ikind), intent(in) :: nGauss
    call debugLog('    Initiating ElementList1D')
    iElem = 0
    this%gaussOrder = nGauss
    this%nElem = nElem
    this%integrator = integrator(nGauss, 'line')
    if(isQuadratic == 0) then
       addElementPtr => addElementLin
       allocate(linLineElem(nElem))
       call debugLog('      Allocated linear line elements: ', size(linLineElem))
       call linLineElem(1)%setnPoint(2)
       call linLineElem(1)%setnDof(2)
       call this%valueElement(linLineElem(1))
    else if(isQuadratic == 1) then
       addElementPtr => addElementQuad
       allocate(quadLineElem(nElem))
       call debugLog('      Allocated quad line elements: ', size(quadLineElem))
       call quadLineElem(1)%setnPoint(3)
       call quadLineElem(1)%setnDof(2)
       call this%valueElement(quadLineElem(1))
    end if
    allocate(this%element(nElem))
    call debugLog('      Allocated element1D pointers: ', size(this%element))
  end subroutine init

  subroutine valueElement(this, elementInput)
    implicit none
    class(StructElementList1DTYPE), intent(inout) :: this
    class(Element1DTYPE), intent(inout) :: elementInput
    integer(ikind) :: i
    integer(ikind) :: n
    n = elementInput%getnPoint()*elementInput%getnDof()
    allocate(this%integrator%shapeFunc(this%integrator%integTerms, n))
    allocate(this%integrator%dShapeFunc(this%integrator%integTerms, 1, n))
    do i = 1, this%integrator%integTerms
       this%integrator%shapeFunc(i, 1:n) = &
            elementInput%shapefunc(this%integrator%gPoint(i,1))
       this%integrator%dShapeFunc(i, 1, 1:n) = &
            elementInput%dShapeFunc(this%integrator%gPoint(i,1))
    end do
  end subroutine valueElement
  
  subroutine addElement(this, id, material, point)
    implicit none
    class(StructElementList1DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(inout) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    call addElementPtr(this, id, material, point)
  end subroutine addElement
  
  subroutine addElementLin(this, id, material, point)
    class(StructElementList1DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(inout) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    integer(ikind) :: i
    iElem = iElem + 1
    linLineElem(iElem) = &
         structLinearLineElement(id, material, point, this%integrator)
    this%element(iElem)%ptr => linLineElem(iElem)
  end subroutine addElementLin

  subroutine addElementQuad(this, id, material, point)
    class(StructElementList1DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    type(MaterialPtrTYPE), intent(inout) :: material
    type(PointPtrTYPE), dimension(:), intent(in) :: point
    integer(ikind) :: i
    iElem = iElem + 1
    quadLineElem(iElem) = &
         structQuadLineElement(id, material, point, this%integrator)
    this%element(iElem)%ptr => quadLineElem(iElem)
  end subroutine addElementQuad
    
  integer(ikind) function getGaussOrder(this)
    implicit none
    class(StructElementList1DTYPE), intent(inout) :: this
    getGaussOrder = this%gaussOrder
  end function getGaussOrder

  integer(ikind) function getnElem(this)
    implicit none
    class(StructElementList1DTYPE), intent(inout) :: this
    getnElem = this%nElem
  end function getnElem

  type(Element1DPtrTYPE) function getElement(this, i)
    implicit none
    class(StructElementList1DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getElement = this%element(i)
  end function getElement

  type(IntegratorPtrTYPE) function getIntegrator(this)
    implicit none
    class(StructElementList1DTYPE), intent(inout) :: this
    call getIntegrator%allocate(this%integrator)
  end function getIntegrator

end module StructElementList1DMOD
