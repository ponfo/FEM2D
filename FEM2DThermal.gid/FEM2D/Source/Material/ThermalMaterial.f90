module ThermalMaterialMOD
  use tools
  use MaterialMOD
  implicit none
  private
  public :: ThermalMaterialTYPE, thermalMaterial
  type, extends(MaterialTYPE) :: ThermalMaterialTYPE
     real(rkind), dimension(2) :: conductivity
   contains
     procedure :: init
  end type ThermalMaterialTYPE
   interface thermalMaterial
     procedure :: constructor
  end interface thermalMaterial
contains
  type(ThermalMaterialTYPE) function constructor(kx,ky)
    implicit none
    real(rkind), intent(in) :: kx, ky
    call constructor%init(kx, ky)
  end function constructor
  subroutine init(this, kx, ky)
    implicit none
    real(rkind), intent(in) :: kx, ky
    class(ThermalMaterialTYPE), intent(inout) :: this
    this%conductivity = (/kx, ky/)
  end subroutine init
end module ThermalMaterialMOD
