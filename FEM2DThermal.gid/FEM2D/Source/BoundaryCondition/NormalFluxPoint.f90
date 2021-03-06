module NormalFluxPointMOD
  use tools
  implicit none
  private
  public :: NormalFluxPointTYPE, normalFluxPoint
  type NormalFluxPointTYPE
     integer(ikind) :: pointID
     real(rkind)    :: value
   contains
     procedure :: init
     procedure :: apply
  end type NormalFluxPointTYPE

  interface normalFluxPoint
     procedure :: constructor
  end interface normalFluxPoint

contains

  type(NormalFluxPointTYPE) function constructor(pointID, value)
    implicit none
    integer(ikind), intent(in) :: pointID
    real(rkind), intent(in) :: value
    call constructor%init(pointID, value)
  end function constructor

  subroutine init(this, pointID, value)
    implicit none
    class(NormalFluxPointTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: pointID
    real(rkind), intent(in) :: value
    this%pointID = pointID
    this%value = value
  end subroutine init

  subroutine apply(this, rhs)
    implicit none
    class(NormalFluxPointTYPE), intent(inout) :: this
    real(rkind), dimension(:), intent(inout) :: rhs
    rhs(this%pointID) = rhs(this%pointID) - this%value
  end subroutine apply

end module NormalFluxPointMOD
