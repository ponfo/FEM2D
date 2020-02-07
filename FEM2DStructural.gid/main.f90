program FEM2DStaticStruct
  use tools
  use DataInputMOD
  use IODataMOD
  use SolverMOD
  use DataOutputMOD
  implicit none
  type(IODataTYPE) :: io
  print'(A)', 'Initiating FEM2DStructural'
  call initFEM2D(io)
  call io%setUp()
  call staticSolver(io)
  call io%postProcess()
  call printResults(resultName = 'Displacement'                 &
       , step         = 1                                       &
       , graphType    = 'Vector'                                &
       , locationName = 'onNodes'                               &
       , resultNumber = io%problem%domain%nPoint                 &
       , component1   = io%problem%dof                           )
  call printResults(resultName = 'NormalStressOnTriangs'        &
       , type         = 'Triangle'                              &
       , step         = 1                                       &
       , graphType    = 'Vector'                                &
       , locationName = 'onGaussPoints'                         &
       , gaussPoints  = io%normalStress%triangGPoint             &
       , resultNumber = size(io%normalStress%triangElemID)       &
       , elemID       = io%normalStress%triangElemID             &
       , component1   = io%normalStress%triangNSx                &
       , component2   = io%normalStress%triangNSy                )
  call printResults(resultName = 'NormalStressOnQuads'          &
       , type         = 'Quadrilateral'                         &
       , step         = 1                                       &
       , graphType    = 'Vector'                                &
       , locationName = 'onGaussPoints'                         &
       , gaussPoints  = io%normalStress%quadGPoint               &
       , resultNumber = size(io%normalStress%quadElemID)            &
       , elemID       = io%normalStress%quadElemID               &
       , component1   = io%normalStress%quadNSx                  &
       , component2   = io%normalStress%quadNSy                   )
  call printResults(resultName = 'ShearStressOnTriangs'         &
       , type         = 'Triangle'                              &
       , step         = 1                                       &
       , graphType    = 'Scalar'                                &
       , locationName = 'onGaussPoints'                         &
       , gaussPoints  = io%shearStress%triangGPoint              &
       , resultNumber = size(io%shearStress%triangElemID)        &
       , elemID       = io%shearStress%triangElemID              &
       , component1   = io%shearStress%triangShS                 )
  call printResults(resultName = 'ShearStressOnQuads'           &
       , type         = 'Quadrilateral'                         &
       , step         = 1                                       &
       , graphType    = 'Scalar'                                &
       , locationName = 'onGaussPoints'                         &
       , gaussPoints  = io%shearStress%quadGPoint                  &
       , resultNumber = size(io%normalStress%quadElemID)           &
       , elemID       = io%shearStress%quadElemID                  &
       , component1   = io%shearStress%quadShS                     )
  call printResults(resultName = 'StrainOnTriangs'              &
       , type         = 'Triangle'                              &
       , step         = 1                                       &
       , graphType    = 'Vector'                                &
       , locationName = 'onGaussPoints'                         &
       , gaussPoints  = io%strain%triangGPoint                   &
       , resultNumber = size(io%strain%triangElemID)             &
       , elemID       = io%strain%triangElemID                   &
       , component1   = io%strain%triangEpx                      &
       , component2   = io%strain%triangEpy                      )
  call printResults(resultName = 'StrainOnQuads'                &
       , type         = 'Quadrilateral'                         &
       , step         = 1                                       &
       , graphType    = 'Vector'                                &
       , locationName = 'onGaussPoints'                         &
       , gaussPoints  = io%strain%quadGPoint                     &
       , resultNumber = size(io%strain%quadElemID)               &
       , elemID       = io%strain%quadElemID                     &
       , component1   = io%strain%quadEpx                        &
       , component2   = io%strain%quadEpy                        )
  call finishProgram() 
end program FEM2DStaticStruct
