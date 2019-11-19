module ProblemMOD
  use tools

  use SparseKit
  implicit none
  private
  public :: ProblemTYPE
  type :: ProblemTYPE
     type(Sparse)                           :: stiffness
     real(rkind), dimension(:), allocatable :: rhs
     real(rkind), dimension(:), allocatable :: dof
  end type ProblemTYPE
end module ProblemMOD
