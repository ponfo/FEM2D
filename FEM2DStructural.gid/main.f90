program FEM2DStaticStruct
  use tools
  use DataInputMOD
  use StructProblemMOD
  use SolverMOD
  use PostProcessMOD
  use DataOutputMOD
  implicit none
  type(StructProblemTYPE) :: problem
  type(StressTYPE) :: stress
  print'(A)', 'Initiating FEM2DStructural'
  call initFEM2D(problem)
  call problem%setUp()
  call staticSolver(problem)
  call stress%calculateStress(problem)
  call printResults(resultName = 'Displacement'                 &
       , step         = 1                                       &
       , graphType    = 'Vector'                                &
       , locationName = 'onNodes'                               &
       , resultNumber = problem%domain%nPoint                    &
       , component1   = problem%dof                              )
    call printResults(resultName = 'StressOnTriangs'        &
       , type         = 'Triangle'                      &
       , step         = 1                               &
       , graphType    = 'Vector'                        &
       , locationName = 'onGaussPoints'                 &
       , gaussPoints  = stress%triangGPoint            &
       , resultNumber = size(Stress%triangElemID)      &
       , elemID       = stress%triangElemID            &
       , component1   = stress%triangSx                &
       , component2   = stress%triangSy                )
  call printResults(resultName = 'StressOnQuads'          &
       , type         = 'Quadrilateral'                 &
       , step         = 1                               &
       , graphType    = 'Vector'                        &
       , locationName = 'onGaussPoints'                 &
       , gaussPoints  = stress%quadGPoint              &
       , resultNumber = size(stress%quadElemID)        &
       , elemID       = stress%quadElemID              &
       , component1   = stress%quadSx                  &
       , component2   = stress%quadSy                  )
  call finishProgram() 
end program FEM2DStaticStruct
