module StructMaterialMOD
  use tools
  use MaterialMOD
  implicit none
  private
  public :: StructMaterialTYPE, structMaterial
  type, extends(MaterialTYPE) :: StructMaterialTYPE
     real(rkind) :: young
     real(rkind) :: poissonCoef
     real(rkind) :: thermalCoef
     real(rkind) :: area
     real(rkind) :: thickness
   contains
     procedure :: init
  end type StructMaterialTYPE

  interface structMaterial
     procedure :: constructor
  end interface structMaterial

contains

  type(StructMaterialTYPE) function constructor(young, poissonCoef, thermalCoef, area, thickness)
    implicit none
    real(rkind), intent(in) :: young
    real(rkind), intent(in) :: poissonCoef
    real(rkind), intent(in) :: thermalCoef
    real(rkind), intent(in) :: area
    real(rkind), intent(in) :: thickness
    call constructor%init(young, poissonCoef, thermalCoef, area, thickness)
  end function constructor

  subroutine init(this, young, poissonCoef, thermalCoef, area, thickness)
    implicit none
    class(StructMaterialTYPE), intent(inout) :: this
    real(rkind), intent(in) :: young
    real(rkind), intent(in) :: poissonCoef
    real(rkind), intent(in) :: thermalCoef
    real(rkind), intent(in) :: area
    real(rkind), intent(in) :: thickness
    this%young       = young
    this%poissonCoef = poissonCoef
    this%thermalCoef = thermalCoef
    this%area        = area
    this%thickness   = thickness
  end subroutine init

end module StructMaterialMOD
