module ConvectionPointMOD
  use tools
  use SparseKit
  implicit none
  private
  public :: ConvectionPointTYPE, convectionPoint
  type ConvectionPointTYPE
     integer(ikind) :: pointID
     real(rkind)    :: coef
     real(rkind)    :: temp
   contains
     procedure :: init
     procedure :: apply
  end type ConvectionPointTYPE

  interface convectionPoint
     procedure :: constructor
  end interface convectionPoint

contains

  type(ConvectionPointTYPE) function constructor(pointID, coef, temp)
    implicit none
    integer(ikind), intent(in) :: pointID
    real(rkind), intent(in) :: coef
    real(rkind), intent(in) :: temp
    call constructor%init(pointID, coef, temp)
  end function constructor

  subroutine init(this, pointID, coef, temp)
    implicit none
    class(ConvectionPointTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: pointID
    real(rkind), intent(in) :: coef
    real(rkind), intent(in) :: temp
    this%pointID = pointID
    this%coef = coef
    this%temp = temp
  end subroutine init

  subroutine apply(this, stiffness, rhs)
    implicit none
    Class(ConvectionPointTYPE), intent(inout) :: this
    type(Sparse), intent(inout) :: stiffness
    real(rkind), dimension(:), intent(inout) :: rhs
    integer(ikind) :: index
    index = stiffness%AI(this%pointID)
    do while(index < stiffness%AI(this%pointID+1))
       if(stiffness%AJ(index) == this%pointID) then
          stiffness%A(index) = stiffness%A(index) + this%coef
          exit
       end if
       index = index + 1
    end do
    rhs(this%pointID) = rhs(this%pointID) + this%coef*this%temp
  end subroutine apply
  
end module ConvectionPointMOD

    
