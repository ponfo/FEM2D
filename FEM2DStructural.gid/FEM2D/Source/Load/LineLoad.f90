module LineLoadMOD
  use tools
  use FunctionOnLinesMOD
  use PointMOD
  use IntegratorMOD
  use IntegratorPtrMOD
  implicit none
  private
  public :: LineLoadTYPE, lineLoad
  type LineLoadTYPE
     integer(ikind), dimension(:), allocatable :: pointID
     integer(ikind)                            :: iLoad
     type(IntegratorPtrTYPE)                   :: integrator1D
   contains
     procedure :: init
     procedure :: apply
  end type LineLoadTYPE
  
  interface lineLoad
     procedure :: constructor
  end interface lineLoad
  
contains
  
  type(LineLoadTYPE) function constructor(pointID, iLoad, integrator1D)
    implicit none
    integer(ikind), dimension(:), intent(in) :: pointID
    integer(ikind), intent(in) :: iLoad
    type(IntegratorTYPE), target, intent(in) :: integrator1D
    call constructor%init(pointID, iLoad, integrator1D)
  end function constructor
  
  subroutine init(this, pointID, iLoad, integrator1D)
    implicit none
    class(LineLoadTYPE), intent(inout) :: this
    integer(ikind), dimension(:), intent(in) :: pointID
    integer(ikind), intent(in) :: iLoad
    type(IntegratorTYPE), target, intent(in) :: integrator1D
    allocate(this%pointID(size(pointID)))
    this%pointID = pointID
    this%iLoad = iLoad
    this%integrator1D%ptr => integrator1D
  end subroutine init

  subroutine apply(this, point, rhs)
    implicit none
    class(LineLoadTYPE), intent(inout) :: this
    type(PointTYPE), dimension(:), intent(inout) :: point
    real(rkind), dimension(:), intent(inout) :: rhs
    integer(ikind) :: i, j
    real(rkind) :: x, y
    real(rkind), dimension(:,:), allocatable :: integral
    real(rkind), dimension(:), allocatable :: jacobian
    real(rkind), dimension(:,:), allocatable :: valuedLoad
    allocate(integral(2,size(this%pointID)))
    allocate(jacobian(this%integrator1D%ptr%integTerms))
    allocate(valuedLoad(2,this%integrator1D%ptr%integTerms))
    integral = 0
    jacobian = 0
    valuedLoad = 0
    do i = 1, this%integrator1D%ptr%integTerms
       x = 0
       y = 0
       do j = 1, size(this%pointID)
          x = x + this%integrator1D%ptr%dShapeFunc(j,1,i)*point(this%pointID(j))%getx()
          y = y + this%integrator1D%ptr%dShapeFunc(j,1,i)*point(this%pointID(j))%gety()
       end do
       jacobian(i) = sqrt(x**2+y**2)
       valuedLoad(1:2,i) = funcOnLines(this%iLoad, x, y)
       do j = 1, size(this%pointID)
          integral(1,j) = integral(1,j) + this%integrator1D%ptr%weight(i) &
               *this%integrator1D%ptr%shapeFunc(j,i)*valuedLoad(1,i)*jacobian(i)
          integral(2,j) = integral(2,j) + this%integrator1D%ptr%weight(i) &
               *this%integrator1D%ptr%shapeFunc(j,i)*valuedLoad(2,i)*jacobian(i)
       end do
    end do
    do i = 1, size(this%pointID)
       rhs(2*this%pointID(i)-1) = rhs(2*this%pointID(i)-1) + integral(1,i)
       rhs(2*this%pointID(i)) = rhs(2*this%pointID(i)) + integral(2,i)
    end do
    deallocate(integral)
    deallocate(jacobian)
    deallocate(valuedLoad)
  end subroutine apply
    
end module LineLoadMOD
