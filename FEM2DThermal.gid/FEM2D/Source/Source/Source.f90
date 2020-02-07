module SourceMOD
  use tools
  use DebuggerMOD

  use PointMOD

  use PointSourceMOD
  use LineSourceMOD
  use SurfaceSourceMOD

  use IntegratorMOD
  use IntegratorPtrMOD

  use ThermalElementList2DMOD
  implicit none
  private
  public :: SourceTYPE, source
  type SourceTYPE
     type(PointSourceTYPE)  , dimension(:), allocatable :: pointSource
     type(LineSourceTYPE)   , dimension(:), allocatable :: lineSource
     type(SurfaceSourceTYPE), dimension(:), allocatable :: surfaceSource
     type(IntegratorTYPE)                               :: integrator1D
   contains
     procedure :: init

     procedure :: addPointSource
     procedure :: getnPointSource
     procedure :: getPointSource
     
     procedure :: addLineSource
     procedure :: getnLineSource
     procedure :: getLineSource

     procedure :: addSurfaceSource
     procedure :: getnSurfaceSource
     procedure :: getSurfaceSource

     procedure :: apply

     procedure, private :: valueIntegrator1DLinear
     procedure, private :: valueIntegrator1DQuadratic
  end type SourceTYPE

  interface source
     procedure constructor
  end interface source

  integer(ikind), save :: iPointSource
  integer(ikind), save :: iLineSource
  integer(ikind), save :: iSurfaceSource

contains

  type(SourceTYPE) function constructor(nPointSource, nLineSource, nSurfaceSource, nGauss, isQuadratic)
    implicit none
    integer(ikind), intent(in) :: nPointSource
    integer(ikind), intent(in) :: nLineSource
    integer(ikind), intent(in) :: nSurfaceSource
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: isQuadratic
    call constructor%init(nPointSource, nLineSource, nSurfaceSource, nGauss, isQuadratic)
  end function constructor

  subroutine init(this, nPointSource, nLineSource, nSurfaceSource, nGauss, isQuadratic)
    implicit none
    class(SourceTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nPointSource
    integer(ikind), intent(in) :: nLineSource
    integer(ikind), intent(in) :: nSurfaceSource
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: isQuadratic
    call debugLog('      Initiating Source')
    allocate(this%pointSource(nPointSource))
    if(isQuadratic == 0) then
       allocate(this%lineSource(nLineSource-1))
    else if(isQuadratic == 1) then
       allocate(this%lineSource((nLineSource-1)/2))
    end if
    allocate(this%surfaceSource(nSurfaceSource))
    call debugLog('        Allocated pointSource: ', size(this%pointSource))
    call debugLog('        Allocated lineSource: ', size(this%lineSource))
    call debugLog('        Allocated surfaceSource: ', size(this%surfaceSource))
    iPointSource = 0
    iLineSource = 0
    iSurfaceSource = 0
    call debugLog('        Valueing Integrator1D for source integrations')
    this%integrator1D = integrator(nGauss, 'line')
    if(isQuadratic == 0) then
       call this%valueIntegrator1DLinear()
    else if(isQuadratic == 1) then
       call this%valueIntegrator1DQuadratic()
    end if
  end subroutine init

  subroutine valueIntegrator1DLinear(this)
    implicit none
    class(SourceTYPE), intent(inout) :: this
    integer(ikind) :: i
    allocate(this%integrator1D%shapeFunc(2, this%integrator1D%integTerms))
    allocate(this%integrator1D%dShapeFunc(2, 1, this%integrator1D%integTerms))
    do i = 1, this%integrator1D%integTerms
       this%integrator1D%shapeFunc(1,i) = 0.5*(1-this%integrator1D%gPoint(i,1))
       this%integrator1D%shapeFunc(2,i) = 0.5*(1+this%integrator1D%gPoint(i,1))
       this%integrator1D%dShapeFunc(1,1,i) = -0.5
       this%integrator1D%dShapeFunc(2,1,i) = 0.5
    end do
  end subroutine valueIntegrator1DLinear

  subroutine valueIntegrator1DQuadratic(this)
    implicit none
    class(SourceTYPE), intent(inout) :: this
    integer(ikind) :: i
    real(rkind) :: u
    allocate(this%integrator1D%shapeFunc(3, this%integrator1D%integTerms))
    allocate(this%integrator1D%dShapeFunc(3, 1, this%integrator1D%integTerms))
    do i = 1, this%integrator1D%integTerms
       u = this%integrator1D%gPoint(i,1)
       this%integrator1D%shapeFunc(1,i) = 0.5*u*(u-1)
       this%integrator1D%shapeFunc(2,i) = (1+u)*(1-u)
       this%integrator1D%shapeFunc(3,i) = 0.5*u*(1+u)
       this%integrator1D%dShapeFunc(1,1,i) = u-0.5
       this%integrator1D%dShapeFunc(2,1,i) = -2*u
       this%integrator1D%dShapeFunc(3,1,i) = u+0.5
    end do
  end subroutine valueIntegrator1DQuadratic
    
  subroutine addPointSource(this, iPoint, iSource)
    implicit none
    class(SourceTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iSource
    iPointSource = iPointSource + 1
    this%pointSource(iPointSource) = pointSource(iPoint, iSource)
  end subroutine addPointSource
  
  type(PointSourceTYPE) function getPointSource(this, i)
    implicit none
    class(SourceTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getPointSource = this%pointSource(i)
  end function getPointSource
  
  integer(ikind) function getnPointSource(this)
    implicit none
    class(SourceTYPE), intent(inout) :: this
    getnPointSource = size(this%pointSource)
  end function getnPointSource
  
  subroutine addLineSource(this, pointID, iSource)
    implicit none
    class(SourceTYPE), intent(inout) :: this
    integer(ikind), dimension(:), intent(in) :: pointID
    integer(ikind), intent(in) :: iSource
    iLineSource = iLineSource + 1
    this%lineSource(iLineSource) = lineSource(pointID, iSource, this%integrator1D)
  end subroutine addLineSource
  
  type(LineSourceTYPE) function getLineSource(this, i)
    implicit none
    class(SourceTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getLineSource = this%lineSource(i)
  end function getLineSource
  
  integer(ikind) function getnLineSource(this)
    implicit none
    class(SourceTYPE), intent(inout) :: this
    getnLineSource = size(this%lineSource)
  end function getnLineSource

  subroutine addSurfaceSource(this, iElem, iSource)
    implicit none
    class(SourceTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iElem
    integer(ikind), intent(in) :: iSource
    iSurfaceSource = iSurfaceSource + 1
    this%surfaceSource(iSurfaceSource) = surfaceSource(iElem, iSource)
  end subroutine addSurfaceSource

  type(SurfaceSourceTYPE) function getSurfaceSource(this, i)
    implicit none
    class(sourceTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getSurfaceSource = this%surfaceSource(i)
  end function getSurfaceSource

  integer(ikind) function getnSurfaceSource(this)
    implicit none
    class(SourceTYPE), intent(inout) :: this
    getnSurfaceSource = size(this%surfaceSource)
  end function getnSurfaceSource

  subroutine apply(this, elementList, point, rhs)
    implicit none
    class(SourceTYPE), intent(inout) :: this
    type(ThermalElementList2DTYPE), intent(inout) :: elementList
    type(PointTYPE), dimension(:), intent(inout) :: point
    real(rkind), dimension(:), intent(inout) :: rhs
    integer(ikind) :: i
    do i = 1, this%getnPointSource()
       call this%pointSource(i)%apply(point, rhs)
    end do
    do i = 1, this%getnLineSource()
       call this%lineSource(i)%apply(point, rhs)
    end do
    do i = 1, this%getnSurfaceSource()
       call this%surfaceSource(i)%apply(elementList, rhs)
    end do
  end subroutine apply

end module SourceMOD

  
    
    
    
    
     
     
  
