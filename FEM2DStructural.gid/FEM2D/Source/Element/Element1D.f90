module Element1DMOD
  use tools
  use ElementMOD
  use IntegratorPtrMOD
  implicit none
  private
  public :: Element1DTYPE
  type, abstract, extends(ElementTYPE) :: Element1DTYPE
     real(rkind) :: lenght
   contains
     procedure                                 :: setLenght
     procedure                                 :: getLenght
     procedure                                 :: jacobian
     procedure(getStiffnessInterf)  , deferred :: getStiffness
     procedure(shapeFuncInterf)     , deferred :: shapeFunc
     procedure(dShapeFuncInterf)    , deferred :: dShapeFunc
  end type Element1DTYPE

  abstract interface
     function shapeFuncInterf(this, u)
       use tools
       import Element1DTYPE
       class(Element1DTYPE), intent(inout) :: this
       real(rkind), intent(in) :: u
       real(rkind), dimension(this%nPoint*this%nDof) :: shapeFuncInterf
     end function shapeFuncInterf
  end interface
  
  abstract interface
     function dShapeFuncInterf(this, u)
       use tools
       import Element1DTYPE
       class(Element1DTYPE), intent(inout) :: this
       real(rkind), intent(in) :: u
       real(rkind), dimension(this%nPoint*this%nDof) :: dShapeFuncInterf
     end function dShapeFuncInterf
  end interface

  abstract interface
     function getStiffnessInterf(this)
       use tools
       use IntegratorPtrMOD
       import Element1DTYPE
       class(Element1DTYPE), intent(inout) :: this
       real(rkind), dimension(this%nPoint*this%nDof,this%nPoint*this%nDof) :: getStiffnessInterf
     end function getStiffnessInterf
  end interface
  
contains

  subroutine setLenght(this)
    implicit none
    class(Element1DTYPE), intent(inout) :: this
    this%lenght = sqrt((this%point(1)%getx()-this%point(2)%getx())**2 &
         + (this%point(1)%gety()-this%point(2)%gety())**2             )
  end subroutine setLenght

  real(rkind) function getLenght(this)
    implicit none
    class(Element1DTYPE), intent(inout) :: this
    getLenght = this%lenght
  end function getLenght

  real(rkind) function jacobian(this, u)
    implicit none
    class(Element1DTYPE), intent(inout) :: this
    real(rkind), intent(in) :: u
    integer(ikind) :: i
    real(rkind) :: term1, term2
    real(rkind), dimension(this%nPoint*this%nDof) :: dsf
    dsf = this%dShapeFunc(u)
    term1 = 0
    term2 = 0
    do i = 1, this%getnPoint()
       term1 = term1 + dsf(i)*this%point(i)%getx()
       term2 = term2 + dsf(i)*this%point(i)%gety()
    end do
    jacobian = sqrt(term1**2+term2**2)
  end function jacobian
  
end module Element1DMOD
