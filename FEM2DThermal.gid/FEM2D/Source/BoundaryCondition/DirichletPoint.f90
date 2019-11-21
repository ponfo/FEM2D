module DirichletPointMOD
  use tools
  use SparseKit
  implicit none
  private
  public :: DirichletPointTYPE, dirichletPoint
  type DirichletPointTYPE
     integer(ikind) :: id
     real(rkind) :: value
   contains
     procedure :: init
     procedure :: apply
  end type DirichletPointTYPE
  
  interface dirichletPoint
     procedure :: constructor
  end interface dirichletPoint
  
contains
  
  type(DirichletPointTYPE) function constructor(id, value)
    implicit none
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    call constructor%init(id, value)
  end function constructor
  
  subroutine init(this, id, value)
    implicit none
    class(DirichletPointTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    this%id = id
    this%value = value
  end subroutine init

  subroutine apply(this, stiffness, rhs)
    implicit none
    class(DirichletPointTYPE), intent(inout) :: this
    class(Sparse), intent(inout) :: stiffness
    real(rkind), dimension(:), intent(inout) :: rhs
    integer(ikind) :: i, index
    index = stiffness%AI(this%id)
    do while(index < stiffness%AI(this%id+1))
       stiffness%A(index) = 0.d0
       if(stiffness%AJ(index) == this%id) then
          stiffness%A(index) = 1.d0
       end if
       index = index + 1
    end do
    rhs(this%id) = this%value
  end subroutine apply
    
end module DirichletPointMOD
