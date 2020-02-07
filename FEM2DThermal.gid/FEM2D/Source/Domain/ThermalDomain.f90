Module ThermalDomainMOD
  use tools
  use DebuggerMOD

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

  implicit none
  private
  public :: ThermalDomainTYPE
  type, extends(DomainTYPE) :: ThermalDomainTYPE
     type(ThermalMaterialTYPE), dimension(:), allocatable :: material
     type(SourceTYPE)                                     :: source
     type(ThermalElementList1DTYPE)                       :: elementList1D
     type(ThermalElementList2DTYPE)                       :: elementList2D
     type(ThermalBoundaryCondition1DTYPE)                 :: bc1D
     type(ThermalBoundaryCondition2DTYPE)                 :: bc2D
   contains
     procedure, public :: applySource
     procedure, public :: applyBC1D
     procedure, public :: applyBC2D
  end type ThermalDomainTYPE

contains

  subroutine applySource(this, rhs)
    implicit none
    class(ThermalDomainTYPE), intent(inout) :: this
    real(rkind), dimension(:), intent(inout) :: rhs
    call this%source%apply(this%elementList2D, this%point, rhs)
  end subroutine applySource

  subroutine applyBC1D(this, stiffness, rhs)
    implicit none
    class(ThermalDomainTYPE), intent(inout) :: this
    class(Sparse), intent(inout) :: stiffness
    real(rkind), dimension(:), intent(inout) :: rhs
    call this%bc1D%apply(stiffness, rhs)
  end subroutine applyBC1D

  subroutine applyBC2D(this, stiffness, rhs)
    implicit none
    class(ThermalDomainTYPE), intent(inout) :: this
    class(Sparse), intent(inout) :: stiffness
    real(rkind), dimension(:), intent(inout) :: rhs
    call this%bc2D%apply(this%elementList2D, stiffness, rhs)
  end subroutine applyBC2D
  
end Module ThermalDomainMOD
    

  
