module Element2DMOD
  use tools
  use ElementMOD
  use IntegratorPtrMOD
  implicit none
  private
  public :: Element2DTYPE
  type, abstract, extends(ElementTYPE) :: Element2DTYPE
     real(rkind) :: area
   contains
     procedure                                 :: getArea
     procedure(getStiffnessInterf)  , deferred :: getStiffness
     procedure(setAreaInterf)       , deferred :: setArea
     procedure(shapeFuncInterf)     , deferred :: shapeFunc
     procedure(dShapeFuncInterf)    , deferred :: dShapeFunc
     procedure                                 :: jacobian
     procedure                                 :: jacobianDet
  end type Element2DTYPE
  
  abstract interface
     subroutine setAreaInterf(this)
       use tools
       import Element2DTYPE
       class(Element2DTYPE), intent(inout) :: this
     end subroutine setAreaInterf
  end interface
  
  abstract interface
     function shapeFuncInterf(this, x, y)
       use tools
       import Element2DTYPE
       class(Element2DTYPE), intent(inout) :: this
       real(rkind), intent(in) :: x
       real(rkind), intent(in) :: y
       real(rkind), dimension(this%nPoint*this%nDof) :: shapeFuncInterf
     end function shapeFuncInterf
  end interface
  
  abstract interface
     function dShapeFuncInterf(this, x, y)
       use tools
       import Element2DTYPE
       class(Element2DTYPE), intent(inout) :: this
       real(rkind), intent(in) :: x
       real(rkind), intent(in) :: y
       real(rkind), dimension(2, this%nPoint*this%nDof) :: dShapeFuncInterf
     end function dShapeFuncInterf
  end interface

  abstract Interface
     function getStiffnessInterf(this)
       use tools
       use IntegratorPtrMOD
       import Element2DTYPE
       class(Element2DTYPE), intent(inout) :: this
       real(rkind), dimension(this%nPoint*this%nDof,this%nPoint*this%nDof) :: getStiffnessInterf
     end function getStiffnessInterf
  end Interface
  
contains
  
  real(rkind) function getArea(this)
    implicit none
    class(Element2DTYPE), intent(inout) :: this
    getArea = this%area
  end function getArea

  function jacobian(this, x, y)
    implicit none
    class(Element2DTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), dimension(2,2) :: jacobian
    integer(ikind) :: i
    real(rkind), dimension(2,this%nPoint*this%nDof) :: dsf
    jacobian = 0.d0
    dsf = this%dShapeFunc(x,y)
    do i = 1, this%nPoint*this%nDof, this%nDof
       jacobian(1,1) = jacobian(1,1) + dsf(1,i)*this%point(i)%getx() !dx/d(xi)
       jacobian(1,2) = jacobian(1,2) + dsf(1,i)*this%point(i)%gety() !dy/d(xi)
       jacobian(2,1) = jacobian(2,1) + dsf(2,i)*this%point(i)%getx() !dx/d(eta)
       jacobian(2,2) = jacobian(2,2) + dsf(2,i)*this%point(i)%gety() !dy/d(eta)
    end do
  end function jacobian

  real(rkind) function jacobianDet(this, jacobian)
    implicit none
    class(Element2DTYPE), intent(inout) :: this
    real(rkind), dimension(2,2), intent(in) :: jacobian
    jacobianDet = jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1)
  end function jacobianDet
  
end module Element2DMOD
