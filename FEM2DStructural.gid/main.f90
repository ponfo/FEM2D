program FEM2DStaticStruct
  use tools
  use DataInputMOD
  use StructProblemMOD
  use SolverMOD
  use DataOutputMOD
  implicit none
  type(StructProblemTYPE) :: problem
  print'(A)', 'Initiating FEM2DStructural'
  call initFEM2D(problem)
  call problem%setUp()
  call staticSolver(problem)
  call printResults(resultName = 'Displacement'                 &
       , step         = 1                                       &
       , graphType    = 'Vector'                                &
       , locationName = 'onNodes'                               &
       , resultNumber = problem%domain%nPoint*problem%domain%nDof &
       , component1   = problem%dof            )
  call finishProgram() 
end program FEM2DStaticStruct
