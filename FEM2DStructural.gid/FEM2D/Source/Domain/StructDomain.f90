module StructDomainMOD
  use tools
  use DebuggerMOD

  use DomainMOD

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

  implicit none
  private
  public :: StructDomainTYPE
  type, extends(DomainTYPE) :: StructDomainTYPE
     type(StructMaterialTYPE), dimension(:), allocatable :: material
     type(LoadTYPE)                                      :: load
     type(StructElementList1DTYPE)                       :: elementList1D
     type(StructElementList2DTYPE)                       :: elementList2D
     type(StructBoundaryCondition1DTYPE)                 :: bc1D
   contains
     procedure, public :: applyLoad
     procedure, public :: applyBC1D
  end type StructDomainTYPE

contains

  subroutine applyLoad(this, rhs)
    implicit none
    class(StructDomainTYPE), intent(inout) :: this
    real(rkind), dimension(:), intent(inout) :: rhs
    call this%load%apply(this%elementList1D, this%elementList2D, this%point, rhs)
  end subroutine applyLoad

  subroutine applyBC1D(this, stiffness, rhs)
    implicit none
    class(StructDomainTYPE), intent(inout) :: this
    class(Sparse), intent(inout) :: stiffness
    real(rkind), dimension(:), intent(inout) :: rhs
    call this%bc1D%apply(stiffness, rhs)
  end subroutine applyBC1D

end module StructDomainMOD
  
    
