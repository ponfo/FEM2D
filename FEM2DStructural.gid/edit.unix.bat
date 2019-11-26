#!/bin/bash

emacs   FEM2D/Source/Lib/Debugger.f90                                \
	FEM2D/Source/Lib/utilities.f90                               \
	FEM2D/Source/Lib/quicksort.f90                               \
	FEM2D/Source/Lib/SparseKit.f90                               \
	FEM2D/Source/Point/Point.f90                                 \
	FEM2D/Source/Point/PointPtr.f90                              \
                                                                     \
	FEM2D/Source/Material/Material.f90                           \
	FEM2D/Source/Material/MaterialPtr.f90                        \
	FEM2D/Source/Material/StructMaterial.f90                     \
                                                                     \
	FEM2D/Source/Geometry/Geometry.f90                           \
	FEM2D/Source/Geometry/GeometryPtr.f90                        \
                                                                     \
	FEM2D/Source/Load/FunctionOnPoints.f90                       \
	FEM2D/Source/Load/FunctionOnLines.f90                        \
	FEM2D/Source/Load/FunctionOnSurfaces.f90                     \
	FEM2D/Source/Load/PointLoad.f90                              \
	FEM2D/Source/Load/LineLoad.f90                               \
	FEM2D/Source/Load/SurfaceLoad.f90                            \
	FEM2D/Source/Load/Pressure.f90                               \
	FEM2D/Source/Load/TemperatureLoad.f90                        \
	FEM2D/Source/Load/Load.f90                                   \
                                                                     \
	FEM2D/Source/BoundaryCondition/FixDisplacement.f90           \
	FEM2D/Source/BoundaryCondition/StructBoundaryCondition1D.f90 \
                                                                     \
	FEM2D/Source/Integrator/Integrator.f90                       \
	FEM2D/Source/Integrator/IntegratorPtr.f90                    \
                                                                     \
	FEM2D/Source/Element/Element.f90                             \
                                                                     \
	FEM2D/Source/Element/1D/Element1D.f90                        \
	FEM2D/Source/Element/1D/StructElement1D.f90                  \
	FEM2D/Source/Element/1D/Element1DPtr.f90                     \
	FEM2D/Source/Element/1D/StructLinLineElement.f90             \
	FEM2D/Source/Element/1D/StructQuadLineElement.f90            \
	FEM2D/Source/Element/1D/StructElementList1D.f90              \
                                                                     \
	FEM2D/Source/Element/2D/Element2D.f90                        \
	FEM2D/Source/Element/2D/StructElement2D.f90                  \
	FEM2D/Source/Element/2D/Element2DPtr.f90                     \
	FEM2D/Source/Element/2D/TriangElement.f90                    \
	FEM2D/Source/Element/2D/StructLinTriangElement.f90           \
	FEM2D/Source/Element/2D/StructQuadTriangElement.f90          \
	FEM2D/Source/Element/2D/QuadElement.f90                      \
	FEM2D/Source/Element/2D/StructLinQuadElement.f90             \
	FEM2D/Source/Element/2D/StructQuadQuadElement.f90            \
	FEM2D/Source/Element/2D/StructElementList2D.f90              \
                                                                     \
	FEM2D/Source/Domain/Domain.f90                               \
	FEM2D/Source/Domain/StructDomain.f90                         \
                                                                     \
	FEM2D/Source/Problem/Problem.f90                             \
	FEM2D/Source/Problem/StructProblem.f90                       \
                                                                     \
	FEM2D/Source/Solver/Solver.f90                               \
                                                                     \
	FEM2D/Source/PostProcess/NormalStress.f90                    \
	FEM2D/Source/PostProcess/ShearStress.f90                     \
	FEM2D/Source/PostProcess/Strain.f90                          \
	                                                             \
	FEM2D/Source/DataIO/DataInput.f90                            \
	FEM2D/Source/DataIO/DataOutput.f90                           &




	
	
