module IODataMOD
  use tools
  use DebuggerMOD

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
  
  use DomainMOD
  use StructDomainMOD

  use ProblemMOD
  use StructProblemMOD

  use NormalStressMOD
  use ShearStressMOD
  use StrainMOD

  implicit none
  private
  public :: IODataTYPE
  type IODataTYPE
     type(StructProblemTYPE) :: problem
     type(NormalStressTYPE)  :: normalStress
     type(ShearStressTYPE)   :: shearStress
     type(StrainTYPE)        :: strain
   contains
     procedure, public  :: initStructProblem
     
     procedure, public  :: addPoint
     procedure, public  :: addMaterial
     procedure, public  :: addElement
     
     procedure, public  :: addPointLoad
     procedure, public  :: addLineLoad
     procedure, public  :: addSurfaceLoad
     procedure, public  :: addPressure
     procedure, public  :: addFixDisplacementX
     procedure, public  :: addFixDisplacementY
     procedure, public  :: setTemperatureLoad

     procedure, public  :: setUp
     procedure, public  :: postProcess

     procedure, private :: addElement1D
     procedure, private :: addElement2D
  end type IODataTYPE

  integer(ikind), save :: iElem
  integer(ikind), save :: iMaterial
  integer(ikind), save :: iPoint

contains

  subroutine initStructProblem(this, nPoint, isQuadratic                          &
       , nLine, nTriang, nQuad, nGauss, nMaterial, nPointLoad                     &
       , nLineLoad, nSurfaceLoad, nPressure, nFixDisplacementX, nFixDisplacementY )
    implicit none
    class(IODataTYPE), intent(inout) :: this
    integer(ikind)   , intent(in)    :: nPoint
    integer(ikind)   , intent(in)    :: isQuadratic
    integer(ikind)   , intent(in)    :: nLine
    integer(ikind)   , intent(in)    :: nTriang
    integer(ikind)   , intent(in)    :: nQuad
    integer(ikind)   , intent(in)    :: nGauss
    integer(ikind)   , intent(in)    :: nMaterial
    integer(ikind)   , intent(in)    :: nPointLoad
    integer(ikind)   , intent(in)    :: nLineLoad
    integer(ikind)   , intent(in)    :: nSurfaceLoad
    integer(ikind)   , intent(in)    :: nPressure
    integer(ikind)   , intent(in)    :: nFixDisplacementX
    integer(ikind)   , intent(in)    :: nFixDisplacementY
    integer(ikind)                   :: nDof
    call debugLog('  Initiating Structural Problem')
    nDof = 2
    call debugLog('    Initiating Structural Domain')
    this%problem%domain%nDof = nDof
    this%problem%domain%nPoint = nPoint
    this%problem%domain%nLine = nLine
    this%problem%domain%nTriang = nTriang
    this%problem%domain%nQuad = nQuad
    allocate(this%problem%domain%point(this%problem%domain%nPoint*this%problem%domain%nDof))
    call debugLog('    Allocated points: ', size(this%problem%domain%point))
    allocate(this%problem%domain%material(nMaterial))
    call debugLog('    Allocated materials: ', size(this%problem%domain%material))
    iElem = 0
    iMaterial = 0
    iPoint = 0
    if(nLine > 0) this%problem%domain%elementList1D = &
         structElementList1D(isQuadratic, nLine, nGauss)
    if(nTriang > 0 .or. nQuad > 0) this%problem%domain%elementList2D = &
         structElementList2D(isQuadratic, nTriang, nQuad, nGauss)
    this%problem%domain%bc1D = &
         structBoundaryCondition1D(nFixDisplacementX, nFixDisplacementY)
    this%problem%domain%load = &
         load(nPointLoad, nLineLoad, nSurfaceLoad, nPressure, nGauss, isQuadratic)
    
    if(isQuadratic == 0) then
       this%problem%stiffness =                                                    &
            sparse(nnz = (nLine*(2*nDof)**2+nTriang*(3*nDof)**2+nQuad*(4*nDof)**2) &
            , rows = nPoint*nDof                                                   )
    else if(isQuadratic == 1) then
       this%problem%stiffness =                                                    &
            sparse(nnz = (nLine*(3*nDof)**2+nTriang*(6*nDof)**2+nQuad*(8*nDof)**2) &
            , rows = nPoint*nDof                                                   )
    end if
    call debugLog('    Allocated Stiffness')
    call debugLog('      Estimated nnz: ', this%problem%stiffness%getnnz())
    call debugLog('      Matrix order: ', this%problem%stiffness%getn())
    allocate(this%problem%rhs(nPoint*nDof))
    call debugLog('    Allocated RHS: ', size(this%problem%rhs))
    this%problem%rhs = 0.d0
    allocate(this%problem%dof(nPoint*nDof))
    call debugLog('    Allocated DOF: ', size(this%problem%dof))
    this%problem%dof = 0.d0
  end subroutine initStructProblem

  subroutine addPoint(this, x, y, z)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), intent(in) :: z
    iPoint = iPoint + 1
    this%problem%domain%point(iPoint) = point(iPoint, x, y, z)
  end subroutine addPoint

  subroutine addMaterial(this, young, poissonCoef, thermalCoef, area, thickness)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    real(rkind), intent(in) :: young
    real(rkind), intent(in) :: poissonCoef
    real(rkind), intent(in) :: thermalCoef
    real(rkind), intent(in) :: area
    real(rkind), intent(in) :: thickness
    iMaterial = iMaterial + 1
    this%problem%domain%material(iMaterial) = &
         structMaterial(young, poissonCoef, thermalCoef, area, thickness)
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

  subroutine addPointLoad(this, iPoint, iLoad)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iLoad
    call this%problem%domain%load%addPointLoad(iPoint, iLoad)
  end subroutine addPointLoad

  subroutine addLineLoad(this, pointID, iLoad)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    integer(ikind), dimension(:), intent(in) :: pointID
    integer(ikind), intent(in) :: iLoad
    call this%problem%domain%load%addLineLoad(pointID, iLoad)
  end subroutine addLineLoad

  subroutine addSurfaceLoad(this, iElem, iLoad)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iElem
    integer(ikind), intent(in) :: iLoad
    call this%problem%domain%load%addSurfaceLoad(iElem, iLoad)
  end subroutine addSurfaceLoad

  subroutine addPressure(this, elemID, pointID, value)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: elemID
    integer(ikind), dimension(:), intent(in) :: pointID
    real(rkind), intent(in) :: value
    call this%problem%domain%load%addPressure(elemID, pointID, value)
  end subroutine addPressure

  subroutine addFixDisplacementX(this, id, value)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    call this%problem%domain%bc1D%addFixDisplacementX(id, value)
  end subroutine addFixDisplacementX

  subroutine addFixDisplacementY(this, id, value)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    call this%problem%domain%bc1D%addFixDisplacementY(id, value)
  end subroutine addFixDisplacementY

  subroutine setTemperatureLoad(this, stableTemp, temperature)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    real(rkind), intent(in) :: stableTemp
    real(rkind), dimension(:), intent(in) :: temperature
    call this%problem%domain%load%setTemperatureLoad(stableTemp, temperature)
  end subroutine setTemperatureLoad

  subroutine setUp(this)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    call this%problem%setUp()
  end subroutine setUp

  subroutine postProcess(this)
    implicit none
    class(IODataTYPE), intent(inout) :: this
    call this%normalStress%calculateNormalStress(this%problem)
    call this%shearStress%calculateShearStress(this%problem)
    call this%strain%calculateStrain(this%problem)
  end subroutine postProcess

end module IODataMOD
