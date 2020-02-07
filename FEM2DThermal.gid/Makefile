# Makefile

COMPILER	:=  ifort
VPATH		:=  $(PWD)/FEM2D/Source
BINDIR		:=  $(PWD)/FEM2D/Bin
OBJECTDIR	:=  $(PWD)/FEM2D/Objects
FFLAGS		:=  -Ofast -free -static -check bounds -traceback -module $(OBJECTDIR)
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
	$(BINDIR)/ThermalMaterial.o             \
	$(BINDIR)/MaterialPtr.o                 \
                                                \
	$(BINDIR)/Geometry.o                    \
	$(BINDIR)/GeometryPtr.o                 \
                                                \
	$(BINDIR)/Element.o                     \
                                                \
	$(BINDIR)/Element1D.o                   \
	$(BINDIR)/ThermalElement1D.o            \
	$(BINDIR)/Element1DPtr.o                \
	$(BINDIR)/ThermalLinLineElement.o       \
	$(BINDIR)/ThermalQuadLineElement.o      \
	$(BINDIR)/ThermalElementList1D.o        \
                                                \
	$(BINDIR)/Element2D.o                   \
	$(BINDIR)/ThermalElement2D.o            \
	$(BINDIR)/Element2DPtr.o                \
	$(BINDIR)/TriangElement.o               \
	$(BINDIR)/ThermalLinTriangElement.o     \
	$(BINDIR)/ThermalQuadTriangElement.o    \
	$(BINDIR)/QuadElement.o                 \
	$(BINDIR)/ThermalLinQuadElement.o       \
	$(BINDIR)/ThermalQuadQuadElement.o      \
	$(BINDIR)/ThermalElementList2D.o        \
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
	$(BINDIR)/ThermalBoundaryCondition1D.o  \
	$(BINDIR)/NormalFluxLine.o              \
	$(BINDIR)/ConvectionLine.o              \
	$(BINDIR)/ThermalBoundaryCondition2D.o  \
                                                \
	$(BINDIR)/Domain.o                      \
	$(BINDIR)/ThermalDomain.o               \
                                                \
	$(BINDIR)/Problem.o                     \
	$(BINDIR)/ThermalProblem.o              \
	                                        \
	$(BINDIR)/HeatFlux.o                    \
	                                        \
	$(BINDIR)/IOData.o                      \
	$(BINDIR)/DataInput.o                   \
	$(BINDIR)/DataOutput.o                  \
                                                \
	$(BINDIR)/Solver.o                      \
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
$(BINDIR)/ThermalMaterial.o : $(VPATH)/Material/ThermalMaterial.f90
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
$(BINDIR)/ThermalElement1D.o : $(VPATH)/Element/1D/ThermalElement1D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Element1DPtr.o : $(VPATH)/Element/1D/Element1DPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThermalLinLineElement.o : $(VPATH)/Element/1D/ThermalLinLineElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThermalQuadLineElement.o : $(VPATH)/Element/1D/ThermalQuadLineElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThermalElementList1D.o : $(VPATH)/Element/1D/ThermalElementList1D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Element2D.o : $(VPATH)/Element/2D/Element2D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThermalElement2D.o : $(VPATH)/Element/2D/ThermalElement2D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Element2DPtr.o : $(VPATH)/Element/2D/Element2DPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/TriangElement.o : $(VPATH)/Element/2D/TriangElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThermalLinTriangElement.o : $(VPATH)/Element/2D/ThermalLinTriangElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThermalQuadTriangElement.o : $(VPATH)/Element/2D/ThermalQuadTriangElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/QuadElement.o : $(VPATH)/Element/2D/QuadElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThermalLinQuadElement.o : $(VPATH)/Element/2D/ThermalLinQuadElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThermalQuadQuadElement.o : $(VPATH)/Element/2D/ThermalQuadQuadElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThermalElementList2D.o : $(VPATH)/Element/2D/ThermalElementList2D.f90
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
$(BINDIR)/ThermalBoundaryCondition1D.o : $(VPATH)/BoundaryCondition/ThermalBoundaryCondition1D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/NormalFluxLine.o : $(VPATH)/BoundaryCondition/NormalFluxLine.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ConvectionLine.o : $(VPATH)/BoundaryCondition/ConvectionLine.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThermalBoundaryCondition2D.o : $(VPATH)/BoundaryCondition/ThermalBoundaryCondition2D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Domain.o : $(VPATH)/Domain/Domain.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThermalDomain.o : $(VPATH)/Domain/ThermalDomain.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Problem.o : $(VPATH)/Problem/Problem.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThermalProblem.o : $(VPATH)/Problem/ThermalProblem.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Solver.o : $(VPATH)/Solver/Solver.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/HeatFlux.o : $(VPATH)/PostProcess/HeatFlux.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/IOData.o : $(VPATH)/DataIO/IOData.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/DataInput.o : $(VPATH)/DataIO/DataInput.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/DataOutput.o : $(VPATH)/DataIO/DataOutput.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

clean:
	rm -f $(BINDIR)/*.o main $(OBJECTDIR)/*.mod
