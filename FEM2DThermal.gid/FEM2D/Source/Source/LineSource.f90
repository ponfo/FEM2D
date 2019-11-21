module LineSourceMOD
  use tools
  use FunctionOnLinesMOD
  use PointMOD
  use IntegratorMOD
  use IntegratorPtrMOD
  implicit none
  private
  public :: LineSourceTYPE, lineSource
  type LineSourceTYPE
     integer(ikind), dimension(:), allocatable :: pointID
     integer(ikind)                            :: iSource
     type(IntegratorPtrTYPE)                   :: integrator1D
   contains
     procedure :: init
     procedure :: apply
  end type LineSourceTYPE
  
  interface lineSource
     procedure :: constructor
  end interface lineSource
  
contains
  
  type(LineSourceTYPE) function constructor(pointID, iSource, integrator1D)
    implicit none
    integer(ikind), dimension(:), intent(in) :: pointID
    integer(ikind), intent(in) :: iSource
    type(IntegratorTYPE), target, intent(in) :: integrator1D
    call constructor%init(pointID, iSource, integrator1D)
  end function constructor
  
  subroutine init(this, pointID, iSource, integrator1D)
    implicit none
    class(LineSourceTYPE), intent(inout) :: this
    integer(ikind), dimension(:), intent(in) :: pointID
    integer(ikind), intent(in) :: iSource
    type(IntegratorTYPE), target, intent(in) :: integrator1D
    allocate(this%pointID(size(pointID)))
    this%pointID = pointID
    this%iSource = iSource
    this%integrator1D%ptr => integrator1D
  end subroutine init

  subroutine apply(this, point, rhs)
    implicit none
    class(LineSourceTYPE), intent(inout) :: this
    type(PointTYPE), dimension(:), intent(inout) :: point
    real(rkind), dimension(:), intent(inout) :: rhs
    integer(ikind) :: i, j
    real(rkind) :: x, y
    real(rkind), dimension(:), allocatable :: integral
    real(rkind), dimension(:), allocatable :: jacobian
    real(rkind), dimension(:), allocatable :: valuedSource
    if(.not.allocated(integral)) allocate(integral(size(this%pointID)))
    if(.not.allocated(jacobian)) allocate(jacobian(this%integrator1D%ptr%integTerms))
    if(.not.allocated(valuedSource)) allocate(valuedSource(this%integrator1D%ptr%integTerms))
    integral = 0
    jacobian = 0
    do i = 1, this%integrator1D%ptr%integTerms
       x = 0
       y = 0
       do j = 1, size(this%pointID)
          x = x + this%integrator1D%ptr%dShapeFunc(j,1,i)*point(this%pointID(j))%getx()
          y = y + this%integrator1D%ptr%dShapeFunc(j,1,i)*point(this%pointID(j))%gety()
       end do
       jacobian(i) = sqrt(x**2+y**2)
       valuedSource(i) = funcOnLines(this%iSource, x, y)
       do j = 1, size(this%pointID)
          integral(j) = integral(j) + this%integrator1D%ptr%weight(i) &
               *this%integrator1D%ptr%shapeFunc(j,i)*valuedSource(i)*jacobian(i)
       end do
    end do
    do i = 1, size(this%pointID)
       rhs(point(this%pointID(i))%getID()) = rhs(point(this%pointID(i))%getID()) + integral(i)
    end do
  end subroutine apply
    
end module LineSourceMOD
