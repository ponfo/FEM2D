module IODataMOD
  use tools
  use DebuggerMOD

  use ThermalProblemMOD

  use DomainMOD

  use PointMOD
  use PointPtrMOD

  use SparseKit
  
  use ThermalMaterialMOD
  use MaterialPtrMOD

  use SourceMOD

  use ThermalBoundaryCondition1DMOD
  use ThermalBoundaryCondition2DMOD

  use Element1DPtrMOD
  use Element2DPtrMOD
  
  use ThermalElementList1DMOD
  use ThermalElementList2DMOD
  
  use HeatFluxMOD

  implicit none
  private
  public :: IODataTYPE
  type :: IODataTYPE
     type(ThermalProblemTYPE) :: problem
     type(HeatFluxTYPE)       :: heatFlux
   contains
     procedure, public  :: initThermalProblem
     
     procedure, public  :: addPoint
     procedure, public  :: addMaterial
     procedure, public  :: addElement
     
     procedure, public  :: addPointSource
     procedure, public  :: addLineSource
     procedure, public  :: addSurfaceSource
     procedure, public  :: addDirichletPoint
     procedure, public  :: addNormalFluxPoint
     procedure, public  :: addNormalFluxLine
     procedure, public  :: addConvectionPoint
     procedure, public  :: addConvectionLine

     procedure, public  :: setUp
     procedure, public  :: postProcess

     procedure, private :: addElement1D
     procedure, private :: addElement2D
  end type IODataTYPE

  integer(ikind), save :: iElem
  integer(ikind), save :: iMaterial
  integer(ikind), save :: iPoint

contains

  subroutine initThermalProblem(this, nPoint, isQuadratic               &
       , nLine, nTriang, nQuad, nGauss, nMaterial, nPointSource         &
       , nLineSource, nSurfaceSource, nDirichletPoint, nNormalFluxPoint &
       , nNormalFluxLine, nConvectionPoint, nConvectionLine             )
    implicit none
    class(IODataTYPE), intent(inout) :: this
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
    call debugLog('  Initiating Thermalermal Problem')
    call debugLog('    Initiating Thermalermal Domain')
    this%problem%domain%nPoint = nPoint
    this%problem%domain%nLine = nLine
    this%problem%domain%nTriang = nTriang
    this%problem%domain%nQuad = nQuad
    allocate(this%problem%domain%point(nPoint))
    call debugLog('    Allocated points: ', size(this%problem%domain%point))
    allocate(this%problem%domain%material(nMaterial))
    call debugLog('    Allocated materials: ', size(this%problem%domain%material))
    iElem = 0
    iMaterial = 0
    iPoint = 0
    if(nLine > 0) this%problem%domain%elementList1D = &
         thermalElementList1D(isQuadratic, nLine, nGauss)
    if(nTriang > 0 .or. nQuad > 0) this%problem%domain%elementList2D = &
         thermalElementList2D(isQuadratic, nTriang, nQuad, nGauss)
    this%problem%domain%bc1D = &
         thermalBoundaryCondition1D(nDirichletPoint, nNormalFluxPoint, nConvectionPoint)
    this%problem%domain%bc2D = &
         thermalBoundaryCondition2D(nNormalFluxLine, nConvectionLine, nGauss, isQuadratic)
    this%problem%domain%source = &
         source(nPointSource, nLineSource, nSurfaceSource, nGauss, isQuadratic) 
    if(isQuadratic == 0) then
       this%problem%stiffness = &
            sparse(nnz = nLine*4+nTriang*9+nQuad*16, rows = nPoint)
    else if(isQuadratic == 1) then
       this%problem%stiffness = &
            sparse(nnz = nLine*9+nTriang*36+nQuad*64, rows = nPoint)
    end if
    call debugLog('    Allocated Stiffness')
    call debugLog('      Estimated nnz: ', this%problem%stiffness%getnnz())
    call debugLog('      Matrix order: ', this%problem%stiffness%getn())
    allocate(this%problem%rhs(nPoint))
    call debugLog('    Allocated RHS: ', size(this%problem%rhs))
    this%problem%rhs = 0.d0
    allocate(this%problem%dof(nPoint))
    call debugLog('    Allocated DOF: ', size(this%problem%dof))
    this%problem%dof = 0.d0
  end subroutine initThermalProblem

  subroutine addPoint(this, x, y, z)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), intent(in) :: z
    iPoint = iPoint + 1
    this%problem%domain%point(iPoint) = point(iPoint, x, y, z)
  end subroutine addPoint
  
  subroutine addMaterial(this, kx, ky)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    real(rkind), intent(in) :: kx
    real(rkind), intent(in) :: ky
    iMaterial = iMaterial + 1
    this%problem%domain%material(iMaterial) = thermalMaterial(kx, ky)
  end subroutine addMaterial

  subroutine addElement(this, type, nPoint, matID, pointList)
    implicit none
    class(IODataTYPE), intent(inout) :: this
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
    class(IODataTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iElem
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: matID
    integer(ikind), dimension(nPoint), intent(in) :: pointList
    integer(ikind) :: i
    type(PointPtrTYPE), dimension(nPoint) :: point
    type(MaterialPtrTYPE) :: material
    do i = 1, nPoint
       call point(i)%allocate(this%problem%domain%point(pointList(i)))
    end do
    call material%allocate(this%problem%domain%material(matID))
    call this%problem%domain%elementList1D%addElement(iElem, material, point)
  end subroutine addElement1D

  subroutine addElement2D(this, iElem, nPoint, matID, pointList)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iElem
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: matID
    integer(ikind), dimension(nPoint), intent(in) :: pointList
    integer(ikind) :: i
    type(PointPtrTYPE), dimension(nPoint) :: point
    type(MaterialPtrTYPE) :: material
    do i = 1, nPoint
       call point(i)%allocate(this%problem%domain%point(pointList(i)))
    end do
    call material%allocate(this%problem%domain%material(matID))
    call this%problem%domain%elementList2D%addElement(iElem, material, point)
  end subroutine addElement2D
  
  subroutine addPointSource(this, iPoint, iSource)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iSource
    call this%problem%domain%source%addPointSource(iPoint, iSource)
  end subroutine addPointSource

  subroutine addLineSource(this, pointID, iSource)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    integer(ikind), dimension(:), intent(in) :: pointID
    integer(ikind), intent(in) :: iSource
    call this%problem%domain%source%addLineSource(pointID, iSource)
  end subroutine addLineSource

  subroutine addSurfaceSource(this, iElem, iSource)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iElem
    integer(ikind), intent(in) :: iSource
    call this%problem%domain%source%addSurfaceSource(iElem, iSource)
  end subroutine addSurfaceSource
    
  subroutine addDirichletPoint(this, id, value)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    call this%problem%domain%bc1D%addDirichletPoint(id, value)
  end subroutine addDirichletPoint

  subroutine addNormalFluxPoint(this, id, value)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    call this%problem%domain%bc1D%addNormalFluxPoint(id, value)
  end subroutine addNormalFluxPoint

  subroutine addNormalFluxLine(this, elemID, pointID, value)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: elemID
    integer(ikind), dimension(:), intent(in) :: pointID
    real(rkind), intent(in) :: value
    call this%problem%domain%bc2D%addNormalFluxLine(elemID, pointID, value)
  end subroutine addNormalFluxLine

  subroutine addConvectionPoint(this, id, coef, temp)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: coef
    real(rkind), intent(in) :: temp
    call this%problem%domain%bc1D%addConvectionPoint(id, coef, temp)
  end subroutine addConvectionPoint
  
  subroutine addConvectionLine(this, elemID, pointID, coef, temp)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: elemID
    integer(ikind), dimension(:), intent(in) :: pointID
    real(rkind), intent(in) :: coef
    real(rkind), intent(in) :: temp
    call this%problem%domain%bc2D%addConvectionLine(elemID, pointID, coef, temp)
  end subroutine addConvectionLine

  subroutine setUp(this)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    call this%problem%setUp()
  end subroutine setUp

  subroutine postProcess(this)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    call this%heatFlux%calculateFlux(this%problem)
  end subroutine postProcess

end module IODataMOD
