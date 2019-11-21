module QuadElementMOD
  use tools
  use PointMOD
  use PointPtrMOD
  use Element2DMOD
  use ThElement2DMOD
  implicit none
  private
  public :: QuadElementTYPE
  type, extends(ThElement2DTYPE), abstract :: QuadElementTYPE
   contains
     procedure, public :: setArea
  end type QuadElementTYPE
contains
  subroutine setArea(this)
    implicit none
    class(QuadElementTYPE), intent(inout) :: this
    this%area = this%point(1)%getx()*this%point(2)%gety() &
         - this%point(1)%gety()*this%point(2)%getx()     &
         + this%point(2)%getx()*this%point(3)%gety()     &
         - this%point(2)%gety()*this%point(3)%getx()     &
         + this%point(3)%getx()*this%point(4)%gety()     &
         - this%point(3)%gety()*this%point(4)%getx()     &
         + this%point(4)%getx()*this%point(1)%gety()     &
         - this%point(4)%gety()*this%point(1)%getx()
    this%area = this%area/2.d0
  end subroutine setArea
end module QuadElementMOD
