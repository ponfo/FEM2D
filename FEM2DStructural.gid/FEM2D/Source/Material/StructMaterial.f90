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
     real(rkind) :: d11, d12, d21, d22, d33
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
    real(rkind) :: factor
    this%young       = young
    this%poissonCoef = poissonCoef
    this%thermalCoef = thermalCoef
    this%area        = area
    this%thickness   = thickness
!!$    !Deformación plana:
!!$    factor = young/((1+poissonCoef)*(1-2*poissonCoef))
!!$    this%d11 = factor*(1-poissonCoef)
!!$    this%d12 = factor*poissonCoef
!!$    this%d21 = factor*poissonCoef
!!$    this%d22 = factor*(1-poissonCoef)
!!$    this%d33 = factor*(1-2*poissonCoef)/2.d0
    !Tensión plana
    factor = young/(1-poissonCoef**2)
    this%d11 = factor
    this%d12 = factor*poissonCoef
    this%d21 = factor*poissonCoef
    this%d22 = factor
    this%d33 = factor*(1-poissonCoef)/2.d0
  end subroutine init

end module StructMaterialMOD
