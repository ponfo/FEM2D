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
  public :: ThDomainTYPE
  type, extends(DomainTYPE) :: ThDomainTYPE
     type(ThMaterialTYPE), dimension(:), allocatable :: material
     type(SourceTYPE)                                :: source
     type(ThElementList1DTYPE)                       :: elementList1D
     type(ThElementList2DTYPE)                       :: elementList2D
     type(ThBoundaryCondition1DTYPE)                 :: bc1D
     type(ThBoundaryCondition2DTYPE)                 :: bc2D
   contains
     procedure, public :: applySource
     procedure, public :: applyBC1D
     procedure, public :: applyBC2D
  end type ThDomainTYPE

contains

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
    

  
