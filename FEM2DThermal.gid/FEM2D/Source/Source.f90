module SourceMOD
  use tools
  use DebuggerMOD

  use PointSourceMOD
  use LineSourceMOD
  use SurfaceSourceMOD

  use IntegratorMOD
  use IntegratorPtrMOD

  use PointMOD
  use ThElementList2DMOD
  implicit none
  private
  public :: SourceTYPE, source
  type SourceTYPE
     type(PointSourceTYPE)  , dimension(:), allocatable :: pointSource
     type(LineSourceTYPE)   , dimension(:), allocatable :: lineSource
     type(surfaceSourceTYPE), dimension(:), allocatable :: surfaceSource
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
  end type SourceTYPE

  interface source
     procedure constructor
  end interface source

  integer(ikind), save :: iPointSource
  integer(ikind), save :: iLineSource
  integer(ikind), save :: iSurfaceSource

contains

  type(SourceTYPE) function constructor(nPointSource, nLineSource, nSurfaceSource)
    implicit none
    integer(ikind), intent(in) :: nPointSource
    integer(ikind), intent(in) :: nLineSource
    integer(ikind), intent(in) :: nSurfaceSource
    call constructor%init(nPointSource, nLineSource, nSurfaceSource)
  end function constructor

  subroutine init(this, nPointSource, nLineSource, nSurfaceSource)
    implicit none
    class(SourceTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nPointSource
    integer(ikind), intent(in) :: nLineSource
    integer(ikind), intent(in) :: nSurfaceSource
    call debugLog('      Initiating Source')
    allocate(this%pointSource(nPointSource))
    allocate(this%lineSource(nLineSource))
    allocate(this%surfaceSource(nSurfaceSource))
    call debugLog('        Allocated pointSource: ', size(this%pointSource))
    call debugLog('        Allocated lineSource: ', size(this%lineSource))
    call debugLog('        Allocated surfaceSource: ', size(this%surfaceSource))
    iPointSource = 0
    iLineSource = 0
    iSurfaceSource = 0
  end subroutine init
    
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
  
  subroutine addLineSource(this, iPoint, iSource)
    implicit none
    class(SourceTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iSource
    iLineSource = iLineSource + 1
    this%lineSource(iLineSource) = lineSource(iPoint, iSource)
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
    type(ThElementList2DTYPE), intent(inout) :: elementList
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

  
    
    
    
    
     
     
  
