module SolverMOD
  use tools
  use DebuggerMOD
  use ProblemMOD
  use ThProblemMOD
  use SparseKit
  implicit none
  private
  public :: staticSolver
  interface staticSolver
     procedure :: staticSolver
  end interface staticSolver
contains
  subroutine staticSolver(problem)
    implicit none
    class(ThProblemTYPE), intent(inout) :: problem
    real(rkind) :: start, finish
    call debugLog('Solving linear system')
    print*, 'Solving linear system'
    call cpu_time(start)
    problem%dof = CGOMP(problem%stiffness, problem%rhs)
    !problem%dof = gmres(problem%stiffness, problem%rhs)
    call cpu_time(finish)
    call debugLog('Done solving')
    print*, 'Solver time = ', (finish-start)/4
  end subroutine staticSolver
end module SolverMOD
