module ThermalBoundaryCondition1DMOD
  use tools
  use DebuggerMOD
  use DirichletPointMOD
  use NormalFluxPointMOD
  use ConvectionPointMOD
  use SparseKit
  implicit none
  private
  public :: ThermalBoundaryCondition1DTYPE, thermalBoundaryCondition1D
  type ThermalBoundaryCondition1DTYPE
     type(DirichletPointTYPE) , dimension(:), allocatable :: dirichletPoint
     type(NormalFluxPointTYPE), dimension(:), allocatable :: normalFluxPoint
     type(ConvectionPointTYPE), dimension(:), allocatable :: convectionPoint
   contains
     procedure :: init

     procedure :: addDirichletPoint
     procedure :: getDirichletPoint
     procedure :: getnDirichlet

     procedure :: addNormalFluxPoint
     procedure :: getNormalFluxPoint
     procedure :: getnNormalFlux

     procedure :: addConvectionPoint
     procedure :: getConvectionPoint
     procedure :: getnConvection

     procedure :: apply
  end type ThermalBoundaryCondition1DTYPE

  interface thermalBoundaryCondition1D
     procedure :: constructor
  end interface thermalBoundaryCondition1D

  integer(ikind), save :: iDirichlet
  integer(ikind), save :: iNormalFlux
  integer(ikind), save :: iConvection

contains

  type(ThermalBoundaryCondition1DTYPE) function constructor(nDirichlet, nNormalFlux, nConvection)
    implicit none
    integer(ikind), intent(in) :: nDirichlet
    integer(ikind), intent(in) :: nNormalFlux
    integer(ikind), intent(in) :: nConvection
    call constructor%init(nDirichlet, nNormalFlux, nConvection)
  end function constructor

  subroutine init(this, nDirichlet, nNormalFlux, nConvection)
    implicit none
    class(ThermalBoundaryCondition1DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nDirichlet
    integer(ikind), intent(in) :: nNormalFlux
    integer(ikind), intent(in) :: nConvection
    call debugLog('      Initiating boundary conditions 1D')
    allocate(this%dirichletPoint(nDirichlet))
    allocate(this%normalFluxPoint(nNormalFlux))
    allocate(this%convectionPoint(nConvection))
    call debugLog('        Allocated dirichletPoints: ', size(this%dirichletPoint))
    call debugLog('        Allocated normalFluxPoint: ', size(this%normalFluxPoint))
    call debugLog('        Allocated convectionPoint: ', size(this%convectionPoint))
    iDirichlet = 0
    iNormalFlux = 0
    iConvection = 0
  end subroutine init

  subroutine addDirichletPoint(this, id, value)
    implicit none
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    class(ThermalBoundaryCondition1DTYPE), intent(inout) :: this
    iDirichlet = iDirichlet + 1
    this%dirichletPoint(iDirichlet) = dirichletPoint(id, value)
  end subroutine addDirichletPoint
  
  type(DirichletPointTYPE) function getDirichletPoint(this, i)
    implicit none
    class(ThermalBoundaryCondition1DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getDirichletPoint = this%dirichletPoint(i)
  end function getDirichletPoint
  
  integer(ikind) function getnDirichlet(this)
    implicit none
    class(ThermalBoundaryCondition1DTYPE), intent(inout) :: this
    getnDirichlet = size(this%dirichletPoint)
  end function getnDirichlet

  subroutine addNormalFluxPoint(this, pointID, value)
    implicit none
    class(ThermalBoundaryCondition1DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: pointID
    real(rkind), intent(in) :: value
    iNormalFlux = iNormalFlux + 1
    this%normalFluxPoint(iNormalFlux) = normalFluxPoint(pointID, value)
  end subroutine addNormalFluxPoint

  type(NormalFluxPointTYPE) function getNormalFluxPoint(this, i)
    implicit none
    class(ThermalBoundaryCondition1DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getNormalFluxPoint = this%normalFluxPoint(i)
  end function getNormalFluxPoint

  integer(ikind) function getnNormalFlux(this)
    implicit none
    class(ThermalBoundarycondition1DTYPE), intent(inout) :: this
    getnNormalFlux = size(this%normalFluxPoint)
  end function getnNormalFlux

  subroutine addConvectionPoint(this, pointID, coef, temp)
    implicit none
    class(ThermalBoundaryCondition1DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: pointID
    real(rkind), intent(in) :: coef
    real(rkind), intent(in) :: temp
    iConvection = iConvection + 1
    this%convectionPoint(iConvection) = convectionPoint(pointID, coef, temp)
  end subroutine addConvectionPoint

  type(ConvectionPointTYPE) function getConvectionPoint(this, i)
    implicit none
    class(ThermalBoundaryCondition1DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: i
    getConvectionPoint = this%convectionPoint(i)
  end function getConvectionPoint

  integer(ikind) function getnConvection(this)
    implicit none
    class(ThermalBoundarycondition1DTYPE), intent(inout) :: this
    getnConvection = size(this%convectionPoint)
  end function getnConvection

  subroutine apply(this, stiffness, rhs)
    implicit none
    class(ThermalBoundaryCondition1DTYPE), intent(inout) :: this
    class(Sparse), intent(inout) :: stiffness
    real(rkind), dimension(:), intent(inout) :: rhs
    integer(ikind) :: i
    do i = 1, this%getnNormalFlux()
       call this%normalFluxPoint(i)%apply(rhs)
    end do
    do i = 1, this%getnConvection()
       call this%convectionPoint(i)%apply(stiffness, rhs)
    end do
    do i = 1, this%getnDirichlet()
       call this%dirichletPoint(i)%apply(stiffness, rhs)
    end do
  end subroutine apply

end module ThermalBoundaryCondition1DMOD


  
    
