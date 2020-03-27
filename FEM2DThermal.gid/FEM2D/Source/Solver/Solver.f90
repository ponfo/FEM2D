
module solvermod
    use tools
    use debuggermod
    use iodatamod
    use sparsekit
    implicit none
    private
    public :: staticsolver
    interface staticsolver
       procedure :: staticsolv
    end interface staticsolver
  contains
    subroutine staticsolv(io)
      implicit none
      class(iodatatype), intent(inout) :: io
      real(rkind) :: start, finish, pt(64)
      integer(ikind) :: maxfct, mnum, mtype, phase, idum(1)
      integer(ikind) :: nrhs, iparm(64), msglvl, error
      call debuglog('solving linear system')
      print'(a)', 'solving linear system'
      print'(/,a)', 'pardiso time'
      call cpu_time(start)
      pt       = 0
      nrhs     = 1
      maxfct   = 1
      mnum     = 1
      iparm    = 0
      iparm(1) = 1 ! user defined iparms
      iparm(2) = 2 ! 3 The parallel (OpenMP) version of the nested dissection algorithm.
                   ! 2 The nested dissection algorithm from the METIS package.
                   ! 0 The minimum degree algorithm.
      iparm(4) = 61! LU-preconditioned CGS iteration with a stopping criterion of
                   ! 1.0E-6 for nonsymmetric matrices.
      iparm(24)= 1 ! two-level factorization algorithm.
      iparm(60)= 1 ! in-core mode or out-core mode
      msglvl   = 0 ! non-print statistical information.
      mtype    = 1 ! real and structurally symmetric 
      phase    = 13! analisys, numerical factorization, solve, iterative refinement
      call pardisomkl(pt, maxfct, mnum, mtype, phase, io%problem%stiffness,&
           idum, nrhs, iparm, msglvl, io%problem%rhs, io%problem%dof, error)
      call cpu_time(finish)
      print'(a,e14.7,/)', 'solver time = ', (finish-start)
      call debuglog('done solving')
    end subroutine staticsolv
end module solvermod
