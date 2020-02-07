module MaterialPtrMOD
  use tools
  use MaterialMOD
  use ThermalMaterialMOD
  implicit none
  private
  public :: MaterialPtrTYPE
  type MaterialPtrTYPE
     class(ThermalMaterialTYPE), pointer :: ptr
   contains
     procedure, public :: allocate
  end type MaterialPtrTYPE

contains

  subroutine allocate(this, material)
    implicit none
    class(MaterialPtrTYPE) :: this
    type(ThermalMaterialTYPE), target, intent(in) :: material
    this%ptr => material
  end subroutine allocate
  
end module MaterialPtrMOD
