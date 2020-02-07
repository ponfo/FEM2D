module LoadMOD
  use tools
  use DebuggerMOD

  use PointMOD

  use PointLoadMOD
  use LineLoadMOD
  use SurfaceLoadMOD
  use TemperatureLoadMOD
  use PressureMOD

  use IntegratorMOD
  use IntegratorPtrMOD

  use StructelementList1DMOD
  use StructElementList2DMOD
  implicit none
  private
  public :: LoadTYPE, load
  type LoadTYPE
     type(PointLoadTYPE)      , dimension(:), allocatable :: pointLoad
     type(LineLoadTYPE)       , dimension(:), allocatable :: lineLoad
     type(SurfaceLoadTYPE)    , dimension(:), allocatable :: surfaceLoad
     type(PressureTYPE)       , dimension(:), allocatable :: pressure
     type(TemperatureLoadTYPE)                            :: tempLoad
     type(IntegratorTYPE)                                 :: integrator1D
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

     procedure :: addPressure
     procedure :: getnPressure
     procedure :: getPressure

     procedure :: setTemperatureLoad

     procedure :: apply

     procedure, private :: valueIntegrator1DLinear
     procedure, private :: valueIntegrator1DQuadratic
     procedure :: initPrintPressure
     procedure :: endPrintPressure
  end type LoadTYPE

  interface load
     procedure constructor
  end interface load

  integer(ikind), save :: iPointLoad
  integer(ikind), save :: iLineLoad
  integer(ikind), save :: iSurfaceLoad
  integer(ikind), save :: iPressure
  logical       , save :: isThereTempLoad

contains

  type(LoadTYPE) function constructor(nPointLoad, nLineLoad, nSurfaceLoad, nPressure, nGauss, isQuadratic)
    implicit none
    integer(ikind), intent(in) :: nPointLoad
    integer(ikind), intent(in) :: nLineLoad
    integer(ikind), intent(in) :: nSurfaceLoad
    integer(ikind), intent(in) :: nPressure
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: isQuadratic
    call constructor%init(nPointLoad, nLineLoad, nSurfaceLoad, nPressure, nGauss, isQuadratic)
  end function constructor

  subroutine init(this, nPointLoad, nLineLoad, nSurfaceLoad, nPressure, nGauss, isQuadratic)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nPointLoad
    integer(ikind), intent(in) :: nLineLoad
    integer(ikind), intent(in) :: nSurfaceLoad
    integer(ikind), intent(in) :: nPressure
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: isQuadratic
    call debugLog('      Initiating Load')
    allocate(this%pointLoad(nPointLoad))
    if(isQuadratic == 0) then
       allocate(this%lineLoad(nLineLoad-1))
    else if(isQuadratic == 1) then
       allocate(this%lineLoad((nLineLoad-1)/2))
    end if
    allocate(this%surfaceLoad(nSurfaceLoad))
    allocate(this%pressure(nPressure))
    call debugLog('        Allocated pointLoad: ', size(this%pointLoad))
    call debugLog('        Allocated lineLoad: ', size(this%lineLoad))
    call debugLog('        Allocated surfaceLoad: ', size(this%surfaceLoad))
    call debugLog('        Allocated pressure: ', size(this%pressure))
    iPointLoad = 0
    iLineLoad = 0
    iSurfaceLoad = 0
    iPressure = 0
    isThereTempLoad = .false.
    call debugLog('        Valueing Integrator1D for load integrations')
    this%integrator1D = integrator(nGauss, 'line')
    if(isQuadratic == 0) then
       call this%valueIntegrator1DLinear()
    else if(isQuadratic == 1) then
       call this%valueIntegrator1DQuadratic()
    end if
  end subroutine init

  subroutine valueIntegrator1DLinear(this)
    implicit none
    class(LoadTYPE), intent(inout) :: this
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
    class(LoadTYPE), intent(inout) :: this
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
  
  subroutine addLineLoad(this, pointID, iLoad)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    integer(ikind), dimension(:), intent(in) :: pointID
    integer(ikind), intent(in) :: iLoad
    iLineLoad = iLineLoad + 1
    this%lineLoad(iLineLoad) = lineLoad(pointID, iLoad, this%integrator1D)
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
    class(LoadTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getSurfaceLoad = this%surfaceLoad(i)
  end function getSurfaceLoad

  integer(ikind) function getnSurfaceLoad(this)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    getnSurfaceLoad = size(this%surfaceLoad)
  end function getnSurfaceLoad

  subroutine addPressure(this, elemID, pointID, value)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: elemID
    integer(ikind), dimension(:), intent(in) :: pointID
    real(rkind), intent(in) :: value
    iPressure = iPressure + 1
    this%pressure(iPressure) = pressure(elemID, pointID, value, this%integrator1D)
  end subroutine addPressure

  type(PressureTYPE) function getPressure(this, i)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getPressure = this%pressure(i)
  end function getPressure

  integer(ikind) function getnPressure(this)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    getnPressure = size(this%pressure)
  end function getnPressure

  subroutine setTemperatureLoad(this, stableTemp, temperature)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    real(rkind), intent(in) :: stableTemp
    real(rkind), dimension(:), intent(in) :: temperature
    this%tempLoad = temperatureLoad(stableTemp, temperature)
    isThereTempLoad = .true.
  end subroutine setTemperatureLoad

  subroutine apply(this, elementList1D, elementList2D, point, rhs)
    implicit none
    class(LoadTYPE), intent(inout) :: this
    type(StructElementList1DTYPE), intent(inout) :: elementList1D
    type(StructElementList2DTYPE), intent(inout) :: elementList2D
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
       call this%surfaceLoad(i)%apply(elementList2D, rhs)
    end do
    call this%initPrintPressure()
    !print*, this%getnPressure()
    do i = 1,  this%getnPressure()
       call this%pressure(i)%apply(elementList2D, point, rhs)
    end do
    call this%endPrintPressure()
    if(isThereTempLoad) call this%tempLoad%apply(elementList1D, elementList2D, rhs)
  end subroutine apply

  subroutine initPrintPressure(this)
    class(LoadTYPE), intent(inout) :: this
    character(100) :: projectName
    integer(ikind) :: i
    open(1, file = 'projectData.dat', status = 'old')
    read(1, '(*(A))') projectName
    close(1)
    open(92, file = 'Pressure.flavia.res')
    write(92,'(A)') 'GiD Post Result File 2.0'
    write(92,'(3A)') 'Result "Pressure" "', trim(projectName), '" 1 Vector OnNodes'
    write(92,'(A)') 'Values'
  end subroutine initPrintPressure
  subroutine endPrintPressure(this)
    class(LoadTYPE), intent(inout) :: this
    write(92,'(A)') 'End Values'
    close(92)
  end subroutine endPrintPressure
  
end module LoadMOD
