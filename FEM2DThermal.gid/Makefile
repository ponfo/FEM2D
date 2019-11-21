# Makefile

COMPILER	:=  ifort
VPATH		:=  $(PWD)/FEM2D/Source
BINDIR		:=  $(PWD)/FEM2D/Bin
OBJECTDIR	:=  $(PWD)/FEM2D/Objects
FFLAGS		:=  -Ofast -free -static -module $(OBJECTDIR)
FFLAGSDebug 	:=  -O0 -fpp -check bounds -traceback -warn nounused -module $(OBJECTDIR)

OBJECTS := $(BINDIR)/Debugger.o                 \
	$(BINDIR)/utilities.o                   \
	$(BINDIR)/quicksort.o                   \
	$(BINDIR)/SparseKit.o                   \
                                                \
	$(BINDIR)/Point.o                       \
	$(BINDIR)/PointPtr.o                    \
                                                \
	$(BINDIR)/Integrator.o                  \
	$(BINDIR)/IntegratorPtr.o               \
                                                \
	$(BINDIR)/Material.o                    \
	$(BINDIR)/ThMaterial.o                  \
	$(BINDIR)/MaterialPtr.o                 \
                                                \
	$(BINDIR)/Geometry.o                    \
	$(BINDIR)/GeometryPtr.o                 \
                                                \
	$(BINDIR)/Element.o                     \
                                                \
	$(BINDIR)/Element1D.o                   \
	$(BINDIR)/ThElement1D.o                 \
	$(BINDIR)/Element1DPtr.o                \
	$(BINDIR)/ThLinLineElement.o            \
	$(BINDIR)/ThQuadLineElement.o           \
	$(BINDIR)/ThElementList1D.o             \
                                                \
	$(BINDIR)/Element2D.o                   \
	$(BINDIR)/ThElement2D.o                 \
	$(BINDIR)/Element2DPtr.o                \
	$(BINDIR)/TriangElement.o               \
	$(BINDIR)/ThLinTriangElement.o          \
	$(BINDIR)/ThQuadTriangElement.o         \
	$(BINDIR)/QuadElement.o                 \
	$(BINDIR)/ThLinQuadElement.o            \
	$(BINDIR)/ThQuadQuadElement.o           \
	$(BINDIR)/ThElementList2D.o             \
                                                \
	$(BINDIR)/functionOnPoints.o            \
	$(BINDIR)/functionOnLines.o             \
	$(BINDIR)/functionOnSurfaces.o          \
	$(BINDIR)/PointSource.o                 \
	$(BINDIR)/LineSource.o                  \
	$(BINDIR)/SurfaceSource.o               \
	$(BINDIR)/Source.o                      \
                                                \
	$(BINDIR)/DirichletPoint.o              \
	$(BINDIR)/NormalFluxPoint.o             \
	$(BINDIR)/ConvectionPoint.o             \
	$(BINDIR)/ThBoundaryCondition1D.o       \
	$(BINDIR)/NormalFluxLine.o              \
	$(BINDIR)/ConvectionLine.o              \
	$(BINDIR)/ThBoundaryCondition2D.o       \
                                                \
	$(BINDIR)/Domain.o                      \
	$(BINDIR)/ThDomain.o                    \
                                                \
	$(BINDIR)/Problem.o                     \
	$(BINDIR)/ThProblem.o                   \
                                                \
	$(BINDIR)/Solver.o                      \
	                                        \
	$(BINDIR)/HeatFlux.o                    \
	                                        \
	$(BINDIR)/DataInput.o                   \
	$(BINDIR)/DataOutput.o                  \
	$(BINDIR)/main.o

main: $(OBJECTS)
	$(COMPILER) $(FFLAGS) $(OBJECTS) -o main

precomp:
	$(COMPILER) $(VPATH)/Lib/Debugger.f90 \
	$(VPATH)/Lib/utilities.f90 precomp.f90 -o precomp; rm -f *.mod

debug: $(OBJECTS)
	$(COMPILER) $(FFLAGSDebug) $(OBJECTS) -o main

$(BINDIR)/main.o : main.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@ 

$(BINDIR)/Debugger.o : $(VPATH)/Lib/Debugger.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/utilities.o : $(VPATH)/Lib/utilities.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/quicksort.o : $(VPATH)/Lib/quicksort.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/SparseKit.o : $(VPATH)/Lib/SparseKit.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Point.o : $(VPATH)/Point/Point.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/PointPtr.o : $(VPATH)/Point/PointPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Integrator.o : $(VPATH)/Integrator/Integrator.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/IntegratorPtr.o : $(VPATH)/Integrator/IntegratorPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Material.o : $(VPATH)/Material/Material.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThMaterial.o : $(VPATH)/Material/ThMaterial.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/MaterialPtr.o : $(VPATH)/Material/MaterialPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Geometry.o : $(VPATH)/Geometry/Geometry.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/GeometryPtr.o : $(VPATH)/Geometry/GeometryPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Element.o : $(VPATH)/Element/Element.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Element1D.o : $(VPATH)/Element/1D/Element1D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThElement1D.o : $(VPATH)/Element/1D/ThElement1D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Element1DPtr.o : $(VPATH)/Element/1D/Element1DPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThLinLineElement.o : $(VPATH)/Element/1D/ThLinLineElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThQuadLineElement.o : $(VPATH)/Element/1D/ThQuadLineElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThElementList1D.o : $(VPATH)/Element/1D/ThElementList1D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Element2D.o : $(VPATH)/Element/2D/Element2D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThElement2D.o : $(VPATH)/Element/2D/ThElement2D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Element2DPtr.o : $(VPATH)/Element/2D/Element2DPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/TriangElement.o : $(VPATH)/Element/2D/TriangElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThLinTriangElement.o : $(VPATH)/Element/2D/ThLinTriangElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThQuadTriangElement.o : $(VPATH)/Element/2D/ThQuadTriangElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/QuadElement.o : $(VPATH)/Element/2D/QuadElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThLinQuadElement.o : $(VPATH)/Element/2D/ThLinQuadElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThQuadQuadElement.o : $(VPATH)/Element/2D/ThQuadQuadElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThElementList2D.o : $(VPATH)/Element/2D/ThElementList2D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/functionOnPoints.o : $(VPATH)/Source/functionOnPoints.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/functionOnLines.o : $(VPATH)/Source/functionOnLines.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/functionOnSurfaces.o : $(VPATH)/Source/functionOnSurfaces.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/PointSource.o : $(VPATH)/Source/PointSource.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/LineSource.o : $(VPATH)/Source/LineSource.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/SurfaceSource.o : $(VPATH)/Source/SurfaceSource.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Source.o : $(VPATH)/Source/Source.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/DirichletPoint.o : $(VPATH)/BoundaryCondition/DirichletPoint.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/NormalFluxPoint.o : $(VPATH)/BoundaryCondition/NormalFluxPoint.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ConvectionPoint.o : $(VPATH)/BoundaryCondition/ConvectionPoint.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThBoundaryCondition1D.o : $(VPATH)/BoundaryCondition/ThBoundaryCondition1D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/NormalFluxLine.o : $(VPATH)/BoundaryCondition/NormalFluxLine.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ConvectionLine.o : $(VPATH)/BoundaryCondition/ConvectionLine.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThBoundaryCondition2D.o : $(VPATH)/BoundaryCondition/ThBoundaryCondition2D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Domain.o : $(VPATH)/Domain/Domain.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThDomain.o : $(VPATH)/Domain/ThDomain.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Problem.o : $(VPATH)/Problem/Problem.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThProblem.o : $(VPATH)/Problem/ThProblem.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Solver.o : $(VPATH)/Solver/Solver.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/HeatFlux.o : $(VPATH)/PostProcess/HeatFlux.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/DataInput.o : $(VPATH)/DataIO/DataInput.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/DataOutput.o : $(VPATH)/DataIO/DataOutput.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

clean:
	rm -f $(BINDIR)/*.o main $(OBJECTDIR)/*.mod
