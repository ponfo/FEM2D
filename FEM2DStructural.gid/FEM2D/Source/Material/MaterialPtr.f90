module MaterialPtrMOD
  use tools
  use MaterialMOD
  use StructMaterialMOD
  implicit none
  private
  public :: MaterialPtrTYPE
  type MaterialPtrTYPE
     class(StructMaterialTYPE), pointer :: ptr
   contains
     procedure, public :: allocate
  end type MaterialPtrTYPE

contains

  subroutine allocate(this, material)
    implicit none
    class(MaterialPtrTYPE) :: this
    type(StructMaterialTYPE), target, intent(in) :: material
    this%ptr => material
  end subroutine allocate
  
end module MaterialPtrMOD
