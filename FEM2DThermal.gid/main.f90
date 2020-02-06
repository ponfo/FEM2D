program FEM2DStaticThermal
  use tools
  use DataInputMOD
  use IODataMOD
  use DataOutputMOD
  implicit none
  type(IODataTYPE) :: io
  real(rkind), dimension(:), allocatable :: zeros
  print'(A)', 'Initiating FEM2DThermal'
  call initFEM2D(io)
  call io%setUp()
  call io%solveStatic()
  call io%postProcess()
  call printResults(resultName = 'Temperature'    &
       , step         = 1                         &
       , graphType    = 'Scalar'                  &
       , locationName = 'onNodes'                 &
       , resultNumber = io%problem%domain%nPoint   &
       , component1   = io%problem%dof             )
  allocate(zeros(size(io%heatFlux%lineQ)))
  zeros = 0.d0
  call printResults(resultName = 'FluxOnLines'             &
       , type         = 'Linear'                           &
       , step         = 1                                  &
       , graphType    = 'Vector'                           &
       , locationName = 'onGaussPoints'                    &
       , gaussPoints  = io%heatFlux%lineGPoint              &
       , resultNumber = size(io%heatFlux%lineElemID)        &
       , elemID       = io%heatFlux%lineElemID              &
       , component1   = io%heatFlux%lineQ                   &
       , component2   = zeros                               )
  call printResults(resultName = 'FluxOnTriangs'           &
       , type         = 'Triangle'                         &
       , step         = 1                                  &
       , graphType    = 'Vector'                           &
       , locationName = 'onGaussPoints'                    &
       , gaussPoints  = io%heatFlux%triangGPoint            &
       , resultNumber = size(io%heatFlux%triangElemID)      &
       , elemID       = io%heatFlux%triangElemID            &
       , component1   = io%heatFlux%triangQx                &
       , component2   = io%heatFlux%triangQy                )
  call printResults(resultName = 'FluxOnQuads'             &
       , type         = 'Quadrilateral'                    &
       , step         = 1                                  &
       , graphType    = 'Vector'                           &
       , locationName = 'onGaussPoints'                    &
       , gaussPoints  = io%heatFlux%quadGPoint              &
       , resultNumber = size(io%heatFlux%quadElemID)        &
       , elemID       = io%heatFlux%quadElemID              &
       , component1   = io%heatFlux%quadQx                  &
       , component2   = io%heatFlux%quadQy                  )
  call printForCoupling(io%problem%domain%nPoint, io%problem%dof)
  call finishProgram() 
end program FEM2DStaticThermal
