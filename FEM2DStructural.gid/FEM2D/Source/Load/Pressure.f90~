module PressureMOD
  use tools

  implicit none
  private
  public :: PressureTYPE, pressure
  type PressureTYPE
     integer(ikind)                            :: elemID
     integer(ikind), dimension(:), allocatable :: pointID
     real(rkind)                               :: value
     type(IntegratorPtrTYPE)                   :: integrator1D
   contains
     procedure, public :: init

     procedure, public :: apply
  end type PressureTYPE

  interface pressure
     procedure :: constructor
  end interface pressure

contains

  type(PressureTYPE) function constructor(elemID, pointID, value, integrator1D)
    implicit none
    integer(ikind), intent(in) :: elemID
    integer(ikind), dimension(:), intent(in) :: pointID
    real(rkind), intent(in) :: value
    type(IntegratorTYPE), target, intent(in) :: integrator1D
    call constructor%init(elemID, pointID, value, integrator1D)
  end function constructor

  subroutine init(this, elemID, pointID, value, integrator1D)
    implicit none
    class(PressureTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: elemID
    integer(ikind), dimension(:), intent(in) :: pointID
    real(rkind), intent(in) :: value
    type(IntegratorTYPE), target, intent(in) :: integrator1D
    this%elemID = elemID
    this%pointID = pointID
    this%value = value
    this%integrator1D%ptr => integrator1D
  end subroutine init

  subroutine apply(this, elementList, rhs)
    implicit none
    class(PressureTYPE), intent(inout) :: this
    type(StructElementList2DTYPE), intent(inout) :: elementList
    real(rkind), dimension(:), intent(inout) :: rhs
    type(Element2DPtrTYPE) :: element
    integer(ikind) :: i, j
    real(rkind) :: x, y, dx, dy
    real(rkind), dimension(:,:), allocatable :: integral
    real(rkind), dimension(:), allocatable :: jacobian
    real(rkind), dimension(:), allocatable :: valuex
    real(rkind), dimension(:), allocatable :: valuey
    allocate(integral(2,size(this%pointID)))
    allocate(jacobian(this%integrator1D%ptr%integTerms))
    allocate(valuex(this%integrator1D%ptr%integTerms))
    allocate(valuey(this%integrator1D%ptr%integTerms))
    integral = 0
    jacobian = 0
    integral = 0
    element = elementList%getElement(this%elemID)
    do i = 1, this%integrator1D%ptr%integTerms
       x = 0
       y = 0
       dx = 0
       dy = 0
       do j = 1, size(pointID)
          x = x + this%integrator1D%ptr%shapeFunc(j,i)*point(this%pointID(j))%getx()
          y = y + this%integrator1D%ptr%shapeFunc(j,i)*point(this%pointID(j))%gety()
          dx = dx + this%integrator1D%ptr%dShapeFunc(j,1,i)*point(this%pointID(j))%getx()
          dy = dy + this%integrator1D%ptr%dShapeFunc(j,1,i)*point(this%pointID(j))%gety()
       end do
       jacobian(i) = sqrt(dx**2+dy**2)
       valuex(i) = value*(-dy)*element%ptr%material%ptr%thickness
       valuey(i) = value*dx*element%ptr%material%ptr%thickness
       do j = 1, size(this%pointID)
          integral(1,j) = integral(1,j) + this%integrator1D%ptr%weight(i) &
               *this%integrator1D%ptr%shapeFunc(j,i)*valuex(i)*jacobian(i)
          integral(2,j) = integral(2,j) + this%integrator1D%ptr%weight(i) &
               *this%integrator1D%ptr%shapeFunc(j,i)*valuey(i)*jacobian(i)
       end do
    end do
    do i = 1, size(this%pointID)
       rhs(2*this%pointID(i)-1) = rhs(2*this%pointID(i)-1) + integral(1,i)
       rhs(2*this%pointID(i)) = rhs(2*this%pointID(i)) + integral(2,i)
    end do
    deallocate(integral)
    deallocate(jacobian)
    deallocate(valuex)
    deallocate(valuey)
  end subroutine apply

  subroutine printPressure(this, elementList, valuex, valuey)
    implicit none
    class(PressureTYPE), intent(inout) :: this
    type(StructElementList2DTYPE), intent(inout) :: elementList
    real(rkind), dimension(:) :: valuex
    real(rkind), dimension(:) :: valuey
    write(92,'(I0,2E16.8)') this%elemID, integral(1,1), integral(2,1)
    do i = 2, size(this%pointID)
       write(92,'(2E16.8)') integral(1,i), integral(2,i)
    end do
  end subroutine printPressure
    
    

end module PressureMOD
