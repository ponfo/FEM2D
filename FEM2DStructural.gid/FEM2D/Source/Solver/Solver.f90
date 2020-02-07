module SolverMOD
  use tools
  use DebuggerMOD
  use IODataMOD
  use SparseKit
  implicit none
  private
  public :: staticSolver
  interface staticSolver
     procedure :: staticSolver
  end interface staticSolver
contains
  subroutine staticSolver(io)
    implicit none
    class(IODataTYPE), intent(inout) :: io
    real(rkind) :: start, finish
    call debugLog('Solving linear system')
    print*, 'Solving linear system'
    call cpu_time(start)
    !io%problem%dof = CGOMP(io%problem%stiffness, io%problem%rhs)
    io%problem%dof = bicGrad(io%problem%stiffness, io%problem%rhs)
    call cpu_time(finish)
    call debugLog('Done solving')
    print*, 'Solver time = ', (finish-start)/4
  end subroutine staticSolver
end module SolverMOD
