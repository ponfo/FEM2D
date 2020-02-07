module StructBoundaryCondition1DMOD
  use tools
  use DebuggerMOD
  use FixDisplacementMOD
  use SparseKit
  implicit none
  private
  public :: StructBoundaryCondition1DTYPE, structBoundaryCondition1D
  type StructBoundaryCondition1DTYPE
     type(FixDisplacementTYPE), dimension(:), allocatable :: fixDisplacementX
     type(FixDisplacementTYPE), dimension(:), allocatable :: fixDisplacementY
   contains
     procedure :: init

     procedure :: addFixDisplacementX
     procedure :: getnFixDisplacementX

     procedure :: addFixDisplacementY
     procedure :: getnFixDisplacementY

     procedure :: apply
  end type StructBoundaryCondition1DTYPE

  interface StructBoundaryCondition1D
     procedure :: constructor
  end interface StructBoundaryCondition1D

  integer(ikind), save :: iFixDisplacementX
  integer(ikind), save :: iFixDisplacementY

contains

  type(StructBoundaryCondition1DTYPE) function constructor(nFixDisplacementX, nFixDisplacementY)
    implicit none
    integer(ikind), intent(in) :: nFixDisplacementX
    integer(ikind), intent(in) :: nFixDisplacementY
    call constructor%init(nFixDisplacementX, nFixDisplacementY)
  end function constructor

  subroutine init(this, nFixDisplacementX, nFixDisplacementY)
    implicit none
    class(StructBoundaryCondition1DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nFixDisplacementX
    integer(ikind), intent(in) :: nFixDisplacementY
    call debugLog('      Initiating boundary conditions 1D')
    allocate(this%fixDisplacementX(nFixDisplacementX))
    allocate(this%fixDisplacementY(nFixDisplacementY))
    call debugLog('        Allocated Fix Displacement in X: ', size(this%fixDisplacementX))
    call debugLog('        Allocated Fix Displacement in Y: ', size(this%fixDisplacementY))
    iFixDisplacementX = 0
    iFixDisplacementY = 0
  end subroutine init

  subroutine addFixDisplacementX(this, id, value)
    implicit none
    class(StructBoundaryCondition1DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    iFixDisplacementX = iFixDisplacementX + 1
    this%fixDisplacementX(iFixDisplacementX) = fixDisplacement(2*id-1, value)
  end subroutine addFixDisplacementX

  integer(ikind) function getnFixDisplacementX(this)
    implicit none
    class(StructBoundaryCondition1DTYPE), intent(inout) :: this
    getnFixDisplacementX = size(this%fixDisplacementX)
  end function getnFixDisplacementX

  subroutine addFixDisplacementY(this, id, value)
    implicit none
    class(StructBoundaryCondition1DTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    iFixDisplacementY = iFixDisplacementY + 1
    this%fixDisplacementY(iFixDisplacementY) = fixDisplacement(2*id, value)
  end subroutine addFixDisplacementY

  integer(ikind) function getnFixDisplacementY(this)
    implicit none
    class(StructBoundaryCondition1DTYPE), intent(inout) :: this
    getnFixDisplacementY = size(this%fixDisplacementY)
  end function getnFixDisplacementY

  subroutine apply(this, stiffness, rhs)
    implicit none
    class(StructBoundaryCondition1DTYPE), intent(inout) :: this
    class(Sparse), intent(inout) :: stiffness
    real(rkind), dimension(:), intent(inout) :: rhs
    integer(ikind) :: i
    do i = 1, this%getnFixDisplacementX()
       call this%fixDisplacementX(i)%apply(stiffness, rhs)
    end do
    do i = 1, this%getnFixDisplacementY()
       call this%fixDisplacementY(i)%apply(stiffness, rhs)
    end do
  end subroutine apply
end module StructBoundaryCondition1DMOD
    
    
    
    
    
