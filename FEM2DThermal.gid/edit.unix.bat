#!/bin/bash

emacs   FEM2D/Source/Lib/Debugger.f90                        \
	FEM2D/Source/Lib/utilities.f90                           \
	FEM2D/Source/Lib/quicksort.f90                           \
	FEM2D/Source/Lib/SparseKit.f90                           \
                                                             \
	FEM2D/Source/Point/Point.f90                             \
	FEM2D/Source/Point/PointPtr.f90                          \
                                                             \
	FEM2D/Source/Material/Material.f90                       \
	FEM2D/Source/Material/MaterialPtr.f90                    \
	FEM2D/Source/Material/ThMaterial.f90                     \
                                                             \
	FEM2D/Source/Geometry/Geometry.f90                       \
	FEM2D/Source/Geometry/GeometryPtr.f90                    \
                                                             \
	FEM2D/Source/Source/FunctionOnPoints.f90                 \
	FEM2D/Source/Source/FunctionOnLines.f90                  \
	FEM2D/Source/Source/FunctionOnSurfaces.f90               \
	FEM2D/Source/Source/PointSource.f90                      \
	FEM2D/Source/Source/LineSource.f90                       \
	FEM2D/Source/Source/SurfaceSource.f90                    \
	FEM2D/Source/Source/Source.f90                           \
                                                             \
	FEM2D/Source/BoundaryCondition/DirichletPoint.f90        \
	FEM2D/Source/BoundaryCondition/NormalFluxPoint.f90       \
	FEM2D/Source/BoundaryCondition/ConvectionPoint.f90       \
	FEM2D/Source/BoundaryCondition/ThBoundaryCondition1D.f90 \
	FEM2D/Source/BoundaryCondition/NormalFluxLine.f90        \
	FEM2D/Source/BoundaryCondition/ConvectionLine.f90        \
	FEM2D/Source/BoundaryCondition/ThBoundaryCondition2D.f90 \
                                                             \
	FEM2D/Source/Integrator/Integrator.f90                   \
	FEM2D/Source/Integrator/IntegratorPtr.f90                \
                                                             \
	FEM2D/Source/Element/Element.f90                         \
                                                             \
	FEM2D/Source/Element/1D/Element1D.f90                    \
	FEM2D/Source/Element/1D/ThElement1D.f90                  \
	FEM2D/Source/Element/1D/Element1DPtr.f90                 \
	FEM2D/Source/Element/1D/ThLinLineElement.f90             \
	FEM2D/Source/Element/1D/ThQuadLineElement.f90            \
	FEM2D/Source/Element/1D/ThElementList1D.f90              \
                                                             \
	FEM2D/Source/Element/2D/Element2D.f90                    \
	FEM2D/Source/Element/2D/ThElement2D.f90                  \
	FEM2D/Source/Element/2D/Element2DPtr.f90                 \
	FEM2D/Source/Element/2D/TriangElement.f90                \
	FEM2D/Source/Element/2D/ThLinTriangElement.f90           \
	FEM2D/Source/Element/2D/ThQuadTriangElement.f90          \
	FEM2D/Source/Element/2D/QuadElement.f90                  \
	FEM2D/Source/Element/2D/ThLinQuadElement.f90             \
	FEM2D/Source/Element/2D/ThQuadQuadElement.f90            \
	FEM2D/Source/Element/2D/ThElementList2D.f90              \
                                                             \
	FEM2D/Source/Domain/Domain.f90                           \
	FEM2D/Source/Domain/ThDomain.f90                         \
                                                             \
	FEM2D/Source/Problem/Problem.f90                         \
	FEM2D/Source/Problem/ThProblem.f90                       \
                                                             \
	FEM2D/Source/Solver/Solver.f90                           \
	                                                         \
	FEM2D/Source/PostProcess/HeatFlux.f90                    \
	                                                         \
	FEM2D/Source/DataIO/IOData.f90                           \
	FEM2D/Source/DataIO/DataInput.f90                        \
	FEM2D/Source/DataIO/DataOutput.f90                       &




	
	
