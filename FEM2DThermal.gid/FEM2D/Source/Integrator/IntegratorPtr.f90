module IntegratorPtrMOD
  use tools
  use DebuggerMOD
  use IntegratorMOD
  implicit none
  private
  public :: IntegratorPtrTYPE
  type :: IntegratorPtrTYPE
     class(IntegratorTYPE), pointer :: ptr
   contains
     procedure, public :: allocate
  end type IntegratorPtrTYPE

contains

  subroutine allocate(this, integrator)
    implicit none
    class(IntegratorPtrTYPE), intent(inout) :: this
    type(IntegratorTYPE), target, intent(in) :: integrator
    this%ptr => integrator
  end subroutine allocate
  
end module IntegratorPtrMOD
