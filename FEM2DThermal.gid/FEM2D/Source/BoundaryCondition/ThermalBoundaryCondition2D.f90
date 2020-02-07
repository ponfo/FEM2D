module ThermalBoundaryCondition2DMOD
  use tools
  use DebuggerMOD
  use NormalFluxLineMOD
  use ConvectionLineMOD
  use IntegratorMOD
  use ThermalElementList2DMOD
  use SparseKit
  implicit none
  private
  public :: ThermalBoundaryCondition2DTYPE, thermalBoundaryCondition2D
  type ThermalBoundaryCondition2DTYPE
     type(NormalFluxLineTYPE), dimension(:), allocatable :: normalFluxLine
     type(ConvectionLineTYPE), dimension(:), allocatable :: convectionLine
     type(IntegratorTYPE)                                :: integrator1D
   contains
     procedure :: init
     
     procedure :: addNormalFluxLine
     procedure :: getNormalFluxLine
     procedure :: getnNormalFlux
     
     procedure :: addConvectionLine
     procedure :: getConvectionLine
     procedure :: getnConvection

     procedure :: apply

     procedure, private :: valueIntegrator1DLinear
     procedure, private :: valueIntegrator1DQuadratic
  end type ThermalBoundaryCondition2DTYPE
  
  interface thermalBoundaryCondition2D
     procedure :: constructor
  end interface thermalBoundaryCondition2D
  
  integer(ikind), save :: iDirichlet
  integer(ikind), save :: iNormalFlux
  integer(ikind), save :: iConvection
  
contains
  
  type(ThermalBoundaryCondition2DTYPE) function constructor(nNormalFlux, nConvection, nGauss, isQuadratic)
    implicit none
    integer(ikind), intent(in) :: nNormalFlux
    integer(ikind), intent(in) :: nConvection
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: isQuadratic
    call constructor%init(nNormalFlux, nConvection, nGauss, isQuadratic)
  end function constructor
  
  subroutine init(this, nNormalFlux, nConvection, nGauss, isQuadratic)
    implicit none
    class(ThermalBoundaryCondition2DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nNormalFlux
    integer(ikind), intent(in) :: nConvection
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: isQuadratic
    call debugLog('      Initiating boundary conditions 2D')
    allocate(this%normalFluxLine(nNormalFlux))
    allocate(this%convectionLine(nConvection))
    call debugLog('        Allocated normalFluxLines: ', size(this%normalFluxLine))
    call debugLog('        Allocated convectionLines: ', size(this%convectionLine))
    iNormalFlux = 0
    iConvection = 0
    call debugLog('        Valueing Integrator1D for boundary integrations')
    this%integrator1D = integrator(nGauss, 'line')
    if(isQuadratic == 0) then
       call this%valueIntegrator1DLinear()
    else if(isQuadratic == 1) then
       call this%valueIntegrator1DQuadratic()
    end if
  end subroutine init

  subroutine valueIntegrator1DLinear(this)
    implicit none
    class(ThermalBoundaryCondition2DTYPE), intent(inout) :: this
    integer(ikind) :: i
    allocate(this%integrator1D%shapeFunc(2, this%integrator1D%integTerms))
    do i = 1, this%integrator1D%integTerms
       this%integrator1D%shapeFunc(1,i) = 0.5*(1-this%integrator1D%gPoint(i,1))
       this%integrator1D%shapeFunc(2,i) = 0.5*(1+this%integrator1D%gPoint(i,1))
    end do
  end subroutine valueIntegrator1DLinear

  subroutine valueIntegrator1DQuadratic(this)
    implicit none
    class(ThermalBoundaryCondition2DTYPE), intent(inout) :: this
    integer(ikind) :: i
    real(rkind) :: u
    allocate(this%integrator1D%shapeFunc(3, this%integrator1D%integTerms))
    do i = 1, this%integrator1D%integTerms
       u = this%integrator1D%gPoint(i,1)
       this%integrator1D%shapeFunc(1,i) = 0.5*u*(u-1)
       this%integrator1D%shapeFunc(3,i) = (1+u)*(1-u)
       this%integrator1D%shapeFunc(2,i) = 0.5*u*(1+u)
    end do
  end subroutine valueIntegrator1DQuadratic
  
  subroutine addNormalFluxLine(this, elemID, pointID, value)
    implicit none
    class(ThermalBoundaryCondition2DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: elemID
    integer(ikind), dimension(:), intent(in) :: pointID
    real(rkind), intent(in) :: value
    iNormalFlux = iNormalFlux + 1
    this%normalFluxLine(iNormalFlux) = normalFluxLine(elemID, pointID, value, this%integrator1D)
  end subroutine addNormalFluxLine
  
  type(NormalFluxLineTYPE) function getNormalFluxLine(this, i)
    implicit none
    class(ThermalBoundaryCondition2DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getNormalFluxLine = this%normalFluxLine(i)
  end function getNormalFluxLine
  
  integer(ikind) function getnNormalFlux(this)
    implicit none
    class(ThermalBoundaryCondition2DTYPE), intent(inout) :: this
    getnNormalFlux = size(this%normalFluxLine)
  end function getnNormalFlux
  
  subroutine addConvectionLine(this, elemID, pointID, coef, temp)
    implicit none
    class(ThermalBoundaryCondition2DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: elemID
    integer(ikind), dimension(:), intent(in) :: pointID
    real(rkind), intent(in) :: coef
    real(rkind), intent(in) :: temp
    iConvection = iConvection + 1
    this%convectionLine(iConvection) = &
         convectionLine(elemID, pointID, coef, temp, this%integrator1D)
  end subroutine addConvectionLine
  
  type(ConvectionLineTYPE) function getConvectionLine(this, i)
    implicit none
    class(ThermalBoundaryCondition2DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getConvectionLine = this%convectionLine(i)
  end function getConvectionLine
  
  integer(ikind) function getnConvection(this)
    implicit none
    class(ThermalBoundaryCondition2DTYPE), intent(inout) :: this
    getnConvection = size(this%convectionLine)
  end function getnConvection

  subroutine apply(this, elementList, stiffness, rhs)
    implicit none
    class(ThermalBoundaryCondition2DTYPE), intent(inout) :: this
    type(ThermalElementList2DTYPE), intent(inout) :: elementList
    class(Sparse), intent(inout) :: stiffness
    real(rkind), dimension(:), intent(inout) :: rhs
    integer(ikind) :: i
    do i = 1, this%getnNormalFlux()
       call this%normalFluxLine(i)%apply(elementList, rhs)
    end do
    do i = 1, this%getnConvection()
       call this%convectionLine(i)%apply(elementList, stiffness, rhs)
    end do
  end subroutine apply
  
end module ThermalBoundaryCondition2DMOD
