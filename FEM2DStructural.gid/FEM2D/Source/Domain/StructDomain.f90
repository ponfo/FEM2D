module StructDomainMOD
  use tools
  use DebuggerMOD

  use DomainMOD

  use PointMOD
  use PointPtrMOD

  use SparseKit

  use StructMaterialMOD
  use MaterialPtrMOD

  use LoadMOD

  use StructBoundaryCondition1DMOD

  use Element1DPtrMOD
  use Element2DPtrMOD

  use StructElementList1DMOD
  use StructElementList2DMOD

  implicit none
  private
  public :: StructDomainTYPE, structDomain
  type, extends(DomainTYPE) :: StructDomainTYPE
     type(StructMaterialTYPE), dimension(:), allocatable :: material
     type(LoadTYPE)                                      :: load
     type(StructElementList1DTYPE)                       :: elementList1D
     type(StructElementList2DTYPE)                       :: elementList2D
     type(StructBoundaryCondition1DTYPE)                 :: bc1D
   contains
     procedure, public :: init

     procedure, public :: addPoint
     procedure, public :: addElement
     procedure, public :: addMaterial
     procedure, public :: addPointLoad
     procedure, public :: addLineLoad
     procedure, public :: addSurfaceLoad
     procedure, public :: addFixDisplacementX
     procedure, public :: addFixDisplacementY
     !procedure, public :: setTemperatureLoad

     procedure, public :: applyLoad
     procedure, public :: applyBC1D

     procedure, private :: addElement1D
     procedure, private :: addElement2D
  end type StructDomainTYPE

  interface structDomain
     procedure :: constructor
  end interface structDomain

  integer(ikind), save :: iElem
  integer(ikind), save :: iMaterial
  integer(ikind), save :: iPoint

contains

  type(StructDomainTYPE) function constructor(nPoint, isQuadratic            &
       , nLine, nTriang, nQuad, nGauss, nMaterial, nPointLoad                &
       , nLineLoad, nSurfaceLoad, nFixDisplacementX, nFixDisplacementY       )
    implicit none
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: isQuadratic
    integer(ikind), intent(in) :: nLine
    integer(ikind), intent(in) :: nTriang
    integer(ikind), intent(in) :: nQuad
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: nMaterial
    integer(ikind), intent(in) :: nPointLoad
    integer(ikind), intent(in) :: nLineLoad
    integer(ikind), intent(in) :: nSurfaceLoad
    integer(ikind), intent(in) :: nFixDisplacementX
    integer(ikind), intent(in) :: nFixDisplacementY
    call constructor%init(nPoint, isQuadratic                                &
       , nLine, nTriang, nQuad, nGauss, nMaterial, nPointLoad                &
       , nLineLoad, nSurfaceLoad, nFixDisplacementX, nFixDisplacementY       )
  end function constructor

  subroutine init(this, nPoint, isQuadratic                                  &
       , nLine, nTriang, nQuad, nGauss, nMaterial, nPointLoad                &
       , nLineLoad, nSurfaceLoad, nFixDisplacementX, nFixDisplacementY       )
    implicit none
    class(StructDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: isQuadratic
    integer(ikind), intent(in) :: nLine
    integer(ikind), intent(in) :: nTriang
    integer(ikind), intent(in) :: nQuad
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: nMaterial
    integer(ikind), intent(in) :: nPointLoad
    integer(ikind), intent(in) :: nLineLoad
    integer(ikind), intent(in) :: nSurfaceLoad
    integer(ikind), intent(in) :: nFixDisplacementX
    integer(ikind), intent(in) :: nFixDisplacementY
    this%nDof = 2
    this%nPoint = nPoint
    this%nLine = nLine
    this%nTriang = nTriang
    this%nQuad = nQuad
    allocate(this%point(this%nPoint*this%nDof))
    call debugLog('    Allocated points: ', size(this%point))
    allocate(this%material(nMaterial))
    call debugLog('    Allocated materials: ', size(this%material))
    iElem = 0
    iMaterial = 0
    iPoint = 0
    if(nLine > 0) this%elementList1D = &
         structElementList1D(isQuadratic, nLine, nGauss)
    if(nTriang > 0 .or. nQuad > 0) this%elementList2D = &
         structElementList2D(isQuadratic, nTriang, nQuad, nGauss)
    this%bc1D = structBoundaryCondition1D(nFixDisplacementX, nFixDisplacementY)
    this%load = load(nPointLoad, nLineLoad, nSurfaceLoad, nGauss, isQuadratic)
  end subroutine init

  subroutine addPoint(this, x, y, z)
    implicit none
    class(StructDomainTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), intent(in) :: z
    iPoint = iPoint + 1
    this%point(iPoint) = point(iPoint, x, y, z)
  end subroutine addPoint

  subroutine addElement(this, type, nPoint, matID, pointList)
    implicit none
    class(StructDomainTYPE), intent(inout) :: this
    character(*), intent(in) :: type
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: matID
    integer(ikind), dimension(nPoint), intent(in) :: pointList
    iElem = iElem + 1
    if(trim(type) == 'Linear' .or. &
       trim(type) == 'linear' .or. &
       trim(type) == 'LINEAR') then
       call this%addElement1D(iElem, nPoint, matID, pointList)
    else if(trim(type) == 'Triangle' .or. &
            trim(type) == 'triangle' .or. &
            trim(type) == 'TRIANGLE') then
       call this%addElement2D(iElem, nPoint, matID, pointList)
    else if(trim(type) == 'Quadrilateral' .or. &
            trim(type) == 'quadrilateral' .or. &
            trim(type) == 'QUADRILATERAL' .or. &
            trim(type) == 'quad')          then
       call this%addElement2D(iElem, nPoint, matID, pointList)
    end if
  end subroutine addElement

  subroutine addElement1D(this, iElem, nPoint, matID, pointList)
    implicit none
    class(StructDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iElem
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: matID
    integer(ikind), dimension(nPoint), intent(in) :: pointList
    integer(ikind) :: i
    type(PointPtrTYPE), dimension(nPoint) :: point
    type(MaterialPtrTYPE) :: material
    do i = 1, nPoint
       call point(i)%allocate(this%point(pointList(i)))
    end do
    call material%allocate(this%material(matID))
    call this%elementList1D%addElement(iElem, material, point)
  end subroutine addElement1D

  subroutine addElement2D(this, iElem, nPoint, matID, pointList)
    implicit none
    class(StructDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iElem
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: matID
    integer(ikind), dimension(nPoint), intent(in) :: pointList
    integer(ikind) :: i
    type(PointPtrTYPE), dimension(nPoint) :: point
    type(MaterialPtrTYPE) :: material
    do i = 1, nPoint
       call point(i)%allocate(this%point(pointList(i)))
    end do
    call material%allocate(this%material(matID))
    call this%elementList2D%addElement(iElem, material, point)
  end subroutine addElement2D

  subroutine addMaterial(this, young, poissonCoef, thermalCoef, area, thickness)
    implicit none
    class(StructDomainTYPE), intent(inout) :: this
    real(rkind), intent(in) :: young
    real(rkind), intent(in) :: poissonCoef
    real(rkind), intent(in) :: thermalCoef
    real(rkind), intent(in) :: area
    real(rkind), intent(in) :: thickness
    iMaterial = iMaterial + 1
    this%material(iMaterial) = structMaterial(young, poissonCoef, thermalCoef, area, thickness)
  end subroutine addMaterial

  subroutine addPointLoad(this, iPoint, iLoad)
    implicit none
    class(StructDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iLoad
    call this%load%addPointLoad(iPoint, iLoad)
  end subroutine addPointLoad

  subroutine addLineLoad(this, pointID, iLoad)
    implicit none
    class(StructDomainTYPE), intent(inout) :: this
    integer(ikind), dimension(:), intent(in) :: pointID
    integer(ikind), intent(in) :: iLoad
    call this%load%addLineLoad(pointID, iLoad)
  end subroutine addLineLoad

  subroutine addSurfaceLoad(this, iElem, iLoad)
    implicit none
    class(StructDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iElem
    integer(ikind), intent(in) :: iLoad
    call this%load%addSurfaceLoad(iElem, iLoad)
  end subroutine addSurfaceLoad

  subroutine addFixDisplacementX(this, id, value)
    implicit none
    class(StructDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    call this%bc1D%addFixDisplacementX(id, value)
  end subroutine addFixDisplacementX

  subroutine addFixDisplacementY(this, id, value)
    implicit none
    class(StructDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    call this%bc1D%addFixDisplacementY(id, value)
  end subroutine addFixDisplacementY

  subroutine applyLoad(this, rhs)
    implicit none
    class(StructDomainTYPE), intent(inout) :: this
    real(rkind), dimension(:), intent(inout) :: rhs
    call this%load%apply(this%elementList2D, this%point, rhs)
  end subroutine applyLoad

  subroutine applyBC1D(this, stiffness, rhs)
    implicit none
    class(StructDomainTYPE), intent(inout) :: this
    class(Sparse), intent(inout) :: stiffness
    real(rkind), dimension(:), intent(inout) :: rhs
    call this%bc1D%apply(stiffness, rhs)
  end subroutine applyBC1D

end module StructDomainMOD
  
    
