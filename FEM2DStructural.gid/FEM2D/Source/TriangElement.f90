module TriangElementMOD
  use tools
  use PointMOD
  use PointPtrMOD
  use Element2DMOD
  use ThElement2DMOD
  implicit none
  private
  public :: TriangElementTYPE
  type, extends(ThElement2DTYPE), abstract :: TriangElementTYPE
     private
   contains
     procedure, public :: setArea
  end type TriangElementTYPE
contains
  subroutine setArea(this)
    implicit none
    class(TriangElementTYPE), intent(inout) :: this
    this%area = this%point(1)%getx()                   &
         * (this%point(2)%gety()-this%point(3)%gety()) &
         + this%point(2)%getx()                        &
         * (this%point(3)%gety()-this%point(1)%gety()) &
         + this%point(3)%getx()                        &
         * (this%point(1)%gety()-this%point(2)%gety())
    this%area = this%area/2.d0
  end subroutine setArea
end module TriangElementMOD
