Module ThDomainMOD
  use tools
  use DebuggerMOD

  use DomainMOD

  use PointMOD
  use PointPtrMOD

  use SparseKit
  
  use ThMaterialMOD
  use MaterialPtrMOD

  use SourceMOD

  use ThBoundaryCondition1DMOD
  use ThBoundaryCondition2DMOD

  use Element1DPtrMOD
  use Element2DPtrMOD
  
  use ThElementList1DMOD
  use ThElementList2DMOD

  implicit none
  private
  public :: ThDomainTYPE, thermalDomain
  type, extends(DomainTYPE) :: ThDomainTYPE
     type(ThMaterialTYPE), dimension(:), allocatable :: material
     type(SourceTYPE)                                :: source
     type(ThElementList1DTYPE)                       :: elementList1D
     type(ThElementList2DTYPE)                       :: elementList2D
     type(ThBoundaryCondition1DTYPE)                 :: bc1D
     type(ThBoundaryCondition2DTYPE)                 :: bc2D
   contains
     procedure, public :: init
     
     procedure, public :: addPoint
     procedure, public :: addElement
     procedure, public :: addMaterial
     procedure, public :: addPointSource
     procedure, public :: addLineSource
     procedure, public :: addSurfaceSource
     procedure, public :: addDirichletPoint
     procedure, public :: addNormalFluxPoint
     procedure, public :: addNormalFluxLine
     procedure, public :: addConvectionPoint
     procedure, public :: addConvectionLine

     procedure, public :: applySource
     procedure, public :: applyBC1D
     procedure, public :: applyBC2D

     procedure, private :: addElement1D
     procedure, private :: addElement2D
  end type ThDomainTYPE

  interface thermalDomain
     procedure :: constructor
  end interface thermalDomain

  integer(ikind), save :: iMaterial
  integer(ikind) :: iPoint

contains

  type(ThDomainTYPE) function constructor(nPoint, isQuadratic           &
       , nLine, nTriang, nQuad, nGauss, nMaterial, nPointSource         &
       , nLineSource, nSurfaceSource, nDirichletPoint, nNormalFluxPoint &
       , nNormalFluxLine, nConvectionPoint, nConvectionLine             )
    implicit none
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: isQuadratic
    integer(ikind), intent(in) :: nLine
    integer(ikind), intent(in) :: nTriang
    integer(ikind), intent(in) :: nQuad
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: nMaterial
    integer(ikind), intent(in) :: nPointSource
    integer(ikind), intent(in) :: nLineSource
    integer(ikind), intent(in) :: nSurfaceSource
    integer(ikind), intent(in) :: nDirichletPoint
    integer(ikind), intent(in) :: nNormalFluxPoint
    integer(ikind), intent(in) :: nNormalFluxLine
    integer(ikind), intent(in) :: nConvectionPoint
    integer(ikind), intent(in) :: nConvectionLine
    call constructor%init(nPoint, isQuadratic                           &
       , nLine, nTriang, nQuad, nGauss, nMaterial, nPointSource         &
       , nLineSource, nSurfaceSource, nDirichletPoint, nNormalFluxPoint &
       , nNormalFluxLine, nConvectionPoint, nConvectionLine             )
  end function constructor

  subroutine init(this, nPoint, isQuadratic                             &
       , nLine, nTriang, nQuad, nGauss, nMaterial, nPointSource         &
       , nLineSource, nSurfaceSource, nDirichletPoint, nNormalFluxPoint &
       , nNormalFluxLine, nConvectionPoint, nConvectionLine             )
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: isQuadratic
    integer(ikind), intent(in) :: nLine
    integer(ikind), intent(in) :: nTriang
    integer(ikind), intent(in) :: nQuad
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: nMaterial
    integer(ikind), intent(in) :: nPointSource
    integer(ikind), intent(in) :: nLineSource
    integer(ikind), intent(in) :: nSurfaceSource
    integer(ikind), intent(in) :: nDirichletPoint
    integer(ikind), intent(in) :: nNormalFluxPoint
    integer(ikind), intent(in) :: nNormalFluxLine
    integer(ikind), intent(in) :: nConvectionPoint
    integer(ikind), intent(in) :: nConvectionLine
    call debugLog('    Initiating Thermal Domain')
    this%nPoint = nPoint
    this%nLine = nLine
    this%nTriang = nTriang
    this%nQuad = nQuad
    allocate(this%point(nPoint))
    call debugLog('    Allocated points: ', size(this%point))
    allocate(this%material(nMaterial))
    call debugLog('    Allocated materials: ', size(this%material))
    iMaterial = 0
    iPoint = 0
    if(nLine > 0) this%elementList1D = &
         thermalElementList1D(isQuadratic, nLine, nGauss)
    if(nTriang > 0 .or. nQuad > 0) this%elementList2D = &
         thermalElementList2D(isQuadratic, nTriang, nQuad, nGauss)
    this%bc1D = thermalBoundaryCondition1D(nDirichletPoint, nNormalFluxPoint, nConvectionPoint)
    this%bc2D = thermalBoundaryCondition2D(nNormalFluxLine, nConvectionLine, nGauss, isQuadratic)
    this%source = source(nPointSource, nLineSource, nSurfaceSource)
  end subroutine init

  subroutine addPoint(this, x, y)
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    iPoint = iPoint + 1
    this%point(iPoint) = point(iPoint, x, y)
  end subroutine addPoint

  subroutine addElement(this, type, nPoint, matID, pointList)
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    character(*), intent(in) :: type
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: matID
    integer(ikind), dimension(nPoint), intent(in) :: pointList
    if(trim(type) == 'Linear' .or. &
       trim(type) == 'linear' .or. &
       trim(type) == 'LINEAR') then
       call this%addElement1D(nPoint, matID, pointList)
    else if(trim(type) == 'Triangle' .or. &
            trim(type) == 'triangle' .or. &
            trim(type) == 'TRIANGLE') then
       call this%addElement2D(nPoint, matID, pointList)
    else if(trim(type) == 'Quadrilateral' .or. &
            trim(type) == 'quadrilateral' .or. &
            trim(type) == 'QUADRILATERAL' .or. &
            trim(type) == 'quad')          then
       call this%addElement2D(nPoint, matID, pointList)
    end if
  end subroutine addElement

  subroutine addElement1D(this, nPoint, matID, pointList)
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: matID
    integer(ikind), dimension(nPoint), intent(in) :: pointList
    integer(ikind) :: i
    type(PointPtrTYPE), dimension(nPoint) :: point
    type(MaterialPtrTYPE) :: material
    do i = 1, nPoint
       call point(i)%allocate(this%point(i))
    end do
    call material%allocate(this%material(matID))
    call this%elementList1D%addElement(material, point)
  end subroutine addElement1D

  subroutine addElement2D(this, nPoint, matID, pointList)
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: matID
    integer(ikind), dimension(nPoint), intent(in) :: pointList
    integer(ikind) :: i
    type(PointPtrTYPE), dimension(nPoint) :: point
    type(MaterialPtrTYPE) :: material
    do i = 1, nPoint
       call point(i)%allocate(this%point(i))
    end do
    call material%allocate(this%material(matID))
    call this%elementList2D%addElement(material, point)
  end subroutine addElement2D

  subroutine addMaterial(this, kx, ky)
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    real(rkind), intent(in) :: kx
    real(rkind), intent(in) :: ky
    iMaterial = iMaterial + 1
    this%material(iMaterial) = thermalMaterial(kx, ky)
  end subroutine addMaterial

  subroutine addPointSource(this, iPoint, iSource)
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iSource
    call this%source%addPointSource(iPoint, iSource)
  end subroutine addPointSource

  subroutine addLineSource(this, iPoint, iSource)
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iSource
    call this%source%addLineSource(iPoint, iSource)
  end subroutine addLineSource

  subroutine addSurfaceSource(this, iElem, iSource)
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iElem
    integer(ikind), intent(in) :: iSource
    call this%source%addSurfaceSource(iElem, iSource)
  end subroutine addSurfaceSource

  subroutine addDirichletPoint(this, id, value)
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    call this%bc1D%addDirichletPoint(id, value)
  end subroutine addDirichletPoint

  subroutine addNormalFluxPoint(this, id, value)
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    call this%bc1D%addNormalFluxPoint(id, value)
  end subroutine addNormalFluxPoint

  subroutine addNormalFluxLine(this, elemID, pointID, value)
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: elemID
    integer(ikind), dimension(:), intent(in) :: pointID
    real(rkind), intent(in) :: value
    call this%bc2D%addNormalFluxLine(elemID, pointID, value)
  end subroutine addNormalFluxLine

  subroutine addConvectionPoint(this, id, coef, temp)
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: coef
    real(rkind), intent(in) :: temp
    call this%bc1D%addConvectionPoint(id, coef, temp)
  end subroutine addConvectionPoint
  
  subroutine addConvectionLine(this, elemID, pointID, coef, temp)
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: elemID
    integer(ikind), dimension(:), intent(in) :: pointID
    real(rkind), intent(in) :: coef
    real(rkind), intent(in) :: temp
    call this%bc2D%addConvectionLine(elemID, pointID, coef, temp)
  end subroutine addConvectionLine

  subroutine applySource(this, rhs)
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    real(rkind), dimension(:), intent(inout) :: rhs
    call this%source%apply(this%elementList2D, this%point, rhs)
  end subroutine applySource

  subroutine applyBC1D(this, stiffness, rhs)
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    class(Sparse), intent(inout) :: stiffness
    real(rkind), dimension(:), intent(inout) :: rhs
    call this%bc1D%apply(stiffness, rhs)
  end subroutine applyBC1D

  subroutine applyBC2D(this, stiffness, rhs)
    implicit none
    class(ThDomainTYPE), intent(inout) :: this
    class(Sparse), intent(inout) :: stiffness
    real(rkind), dimension(:), intent(inout) :: rhs
    call this%bc2D%apply(this%elementList2D, stiffness, rhs)
  end subroutine applyBC2D
  
end Module ThDomainMOD
    

  
