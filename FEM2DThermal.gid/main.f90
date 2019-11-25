program FEM2DStaticThermal
  use tools
  use DataInputMOD
  use ThProblemMOD
  use SolverMOD
  use HeatFluxMOD
  use DataOutputMOD
  implicit none
  type(ThProblemTYPE) :: problem
  type(HeatFluxTYPE) :: heatFlux
  real(rkind), dimension(:), allocatable :: zeros
  print'(A)', 'Initiating FEM2DThermal'
  call initFEM2D(problem)
  call problem%setUp()
  call staticSolver(problem)
  call heatFlux%calculateFlux(problem)
  call printResults(resultName = 'Temperature' &
       , step         = 1                      &
       , graphType    = 'Scalar'               &
       , locationName = 'onNodes'              &
       , resultNumber = problem%domain%nPoint   &
       , component1   = problem%dof            )
  allocate(zeros(size(heatFlux%lineQ)))
  zeros = 0.d0
  call printResults(resultName = 'FluxOnLines'          &
       , type         = 'Linear'                        &
       , step         = 1                               &
       , graphType    = 'Vector'                        &
       , locationName = 'onGaussPoints'                 &
       , gaussPoints  = heatFlux%lineGPoint              &
       , resultNumber = size(heatFlux%lineElemID)        &
       , elemID       = heatFlux%lineElemID              &
       , component1   = heatFlux%lineQ                   &
       , component2   = zeros                            )
  call printResults(resultName = 'FluxOnTriangs'        &
       , type         = 'Triangle'                      &
       , step         = 1                               &
       , graphType    = 'Vector'                        &
       , locationName = 'onGaussPoints'                 &
       , gaussPoints  = heatFlux%triangGPoint            &
       , resultNumber = size(heatFlux%triangElemID)      &
       , elemID       = heatFlux%triangElemID            &
       , component1   = heatFlux%triangQx                &
       , component2   = heatFlux%triangQy                )
  call printResults(resultName = 'FluxOnQuads'          &
       , type         = 'Quadrilateral'                 &
       , step         = 1                               &
       , graphType    = 'Vector'                        &
       , locationName = 'onGaussPoints'                 &
       , gaussPoints  = heatFlux%quadGPoint              &
       , resultNumber = size(heatFlux%quadElemID)        &
       , elemID       = heatFlux%quadElemID              &
       , component1   = heatFlux%quadQx                  &
       , component2   = heatFlux%quadQy                  )
  call printForCoupling(problem%domain%nPoint, problem%dof)
  call finishProgram() 
end program FEM2DStaticThermal
