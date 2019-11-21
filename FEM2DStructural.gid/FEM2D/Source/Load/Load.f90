module LoadMOD
  use tools
  use DebuggerMOD

  use PointMOD

  use PointLoadMOD
  use LineLoadMOD
  use SurfaceLoadMOD

  use IntegratorMOD
  use IntegratorPtrMOD

  use ThElementList2DMOD
  implicit none
  private
  public :: LoadTYPE, source
  type LoadTYPE
     type(PointLoadTYPE)  , dimension(:), allocatable :: pointLoad
     type(LineLoadTYPE)   , dimension(:), allocatable :: lineLoad
     type(surfaceLoadTYPE), dimension(:), allocatable :: surfaceLoad
   contains
     procedure :: init

     procedure :: addPointLoad
     procedure :: getnPointLoad
     procedure :: getPointLoad
     
     procedure :: addLineLoad
     procedure :: getnLineLoad
     procedure :: getLineLoad

     procedure :: addSurfaceLoad
     procedure :: getnSurfaceLoad
     procedure :: getSurfaceLoad

     procedure :: apply
  end type LoadTYPE

  interface source
     procedure constructor
  end interface source

  integer(ikind), save :: iPointLoad
  integer(ikind), save :: iLineLoad
  integer(ikind), save :: iSurfaceLoad

contains

  type(LoadTYPE) function constructor(nPointLoad, nLineLoad, nSurfaceLoad)
    implicit none
    integer(ikind), intent(in) :: nPointLoad
    integer(ikind), intent(in) :: nLineLoad
    integer(ikind), intent(in) :: nSurfaceLoad
    call constructor%init(nPointLoad, nLineLoad, nSurfaceLoad)
  end function constructor

  subroutine init(this, nPointLoad, nLineLoad, nSurfaceLoad)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nPointLoad
    integer(ikind), intent(in) :: nLineLoad
    integer(ikind), intent(in) :: nSurfaceLoad
    call debugLog('      Initiating Load')
    allocate(this%pointLoad(nPointLoad))
    allocate(this%lineLoad(nLineLoad))
    allocate(this%surfaceLoad(nSurfaceLoad))
    call debugLog('        Allocated pointLoad: ', size(this%pointLoad))
    call debugLog('        Allocated lineLoad: ', size(this%lineLoad))
    call debugLog('        Allocated surfaceLoad: ', size(this%surfaceLoad))
    iPointLoad = 0
    iLineLoad = 0
    iSurfaceLoad = 0
  end subroutine init
    
  subroutine addPointLoad(this, iPoint, iLoad)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iLoad
    iPointLoad = iPointLoad + 1
    this%pointLoad(iPointLoad) = pointLoad(iPoint, iLoad)
  end subroutine addPointLoad
  
  type(PointLoadTYPE) function getPointLoad(this, i)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getPointLoad = this%pointLoad(i)
  end function getPointLoad
  
  integer(ikind) function getnPointLoad(this)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    getnPointLoad = size(this%pointLoad)
  end function getnPointLoad
  
  subroutine addLineLoad(this, iPoint, iLoad)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iLoad
    iLineLoad = iLineLoad + 1
    this%lineLoad(iLineLoad) = lineLoad(iPoint, iLoad)
  end subroutine addLineLoad
  
  type(LineLoadTYPE) function getLineLoad(this, i)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getLineLoad = this%lineLoad(i)
  end function getLineLoad
  
  integer(ikind) function getnLineLoad(this)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    getnLineLoad = size(this%lineLoad)
  end function getnLineLoad

  subroutine addSurfaceLoad(this, iElem, iLoad)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iElem
    integer(ikind), intent(in) :: iLoad
    iSurfaceLoad = iSurfaceLoad + 1
    this%surfaceLoad(iSurfaceLoad) = surfaceLoad(iElem, iLoad)
  end subroutine addSurfaceLoad

  type(SurfaceLoadTYPE) function getSurfaceLoad(this, i)
    implicit none
    class(sourceTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getSurfaceLoad = this%surfaceLoad(i)
  end function getSurfaceLoad

  integer(ikind) function getnSurfaceLoad(this)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    getnSurfaceLoad = size(this%surfaceLoad)
  end function getnSurfaceLoad

  subroutine apply(this, elementList, point, rhs)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    type(ThElementList2DTYPE), intent(inout) :: elementList
    type(PointTYPE), dimension(:), intent(inout) :: point
    real(rkind), dimension(:), intent(inout) :: rhs
    integer(ikind) :: i
    do i = 1, this%getnPointLoad()
       call this%pointLoad(i)%apply(point, rhs)
    end do
    do i = 1, this%getnLineLoad()
       call this%lineLoad(i)%apply(point, rhs)
    end do
    do i = 1, this%getnSurfaceLoad()
       call this%surfaceLoad(i)%apply(elementList, rhs)
    end do
  end subroutine apply

end module LoadMOD

  
    
    
    
    
     
     
  
