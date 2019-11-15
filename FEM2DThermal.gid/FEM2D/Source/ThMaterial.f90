module ThMaterialMOD
  use tools
  use MaterialMOD
  implicit none
  private
  public :: ThMaterialTYPE, thermalMaterial
  type, extends(MaterialTYPE) :: ThMaterialTYPE
     real(rkind), dimension(2) :: conductivity
   contains
     procedure :: init
  end type ThMaterialTYPE
   interface thermalMaterial
     procedure :: constructor
  end interface thermalMaterial
contains
  type(ThMaterialTYPE) function constructor(kx,ky)
    implicit none
    real(rkind), intent(in) :: kx, ky
    call constructor%init(kx, ky)
  end function constructor
  subroutine init(this, kx, ky)
    implicit none
    real(rkind), intent(in) :: kx, ky
    class(ThMaterialTYPE), intent(inout) :: this
    this%conductivity = (/kx, ky/)
  end subroutine init
end module ThMaterialMOD
