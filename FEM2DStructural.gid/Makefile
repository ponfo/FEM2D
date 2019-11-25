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
	$(BINDIR)/StructMaterial.o              \
	$(BINDIR)/MaterialPtr.o                 \
                                                \
	$(BINDIR)/Geometry.o                    \
	$(BINDIR)/GeometryPtr.o                 \
                                                \
	$(BINDIR)/Element.o                     \
                                                \
	$(BINDIR)/Element1D.o                   \
	$(BINDIR)/StructElement1D.o             \
	$(BINDIR)/Element1DPtr.o                \
	$(BINDIR)/StructLinLineElement.o        \
	$(BINDIR)/StructQuadLineElement.o       \
	$(BINDIR)/StructElementList1D.o         \
                                                \
	$(BINDIR)/Element2D.o                   \
	$(BINDIR)/StructElement2D.o             \
	$(BINDIR)/Element2DPtr.o                \
	$(BINDIR)/TriangElement.o               \
	$(BINDIR)/StructLinTriangElement.o      \
	$(BINDIR)/StructQuadTriangElement.o     \
	$(BINDIR)/QuadElement.o                 \
	$(BINDIR)/StructLinQuadElement.o        \
	$(BINDIR)/StructQuadQuadElement.o       \
	$(BINDIR)/StructElementList2D.o         \
                                                \
	$(BINDIR)/functionOnPoints.o            \
	$(BINDIR)/functionOnLines.o             \
	$(BINDIR)/functionOnSurfaces.o          \
	$(BINDIR)/PointLoad.o                   \
	$(BINDIR)/LineLoad.o                    \
	$(BINDIR)/SurfaceLoad.o	                \
	$(BINDIR)/TemperatureLoad.o             \
	$(BINDIR)/Load.o                        \
                                                \
	$(BINDIR)/FixDisplacement.o             \
	$(BINDIR)/StructBoundaryCondition1D.o   \
                                                \
	$(BINDIR)/Domain.o                      \
	$(BINDIR)/StructDomain.o                \
                                                \
	$(BINDIR)/Problem.o                     \
	$(BINDIR)/StructProblem.o               \
                                                \
	$(BINDIR)/Solver.o                      \
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
$(BINDIR)/StructMaterial.o : $(VPATH)/Material/StructMaterial.f90
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
$(BINDIR)/StructElement1D.o : $(VPATH)/Element/1D/StructElement1D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Element1DPtr.o : $(VPATH)/Element/1D/Element1DPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/StructLinLineElement.o : $(VPATH)/Element/1D/StructLinLineElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/StructQuadLineElement.o : $(VPATH)/Element/1D/StructQuadLineElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/StructElementList1D.o : $(VPATH)/Element/1D/StructElementList1D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Element2D.o : $(VPATH)/Element/2D/Element2D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/StructElement2D.o : $(VPATH)/Element/2D/StructElement2D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Element2DPtr.o : $(VPATH)/Element/2D/Element2DPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/TriangElement.o : $(VPATH)/Element/2D/TriangElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/StructLinTriangElement.o : $(VPATH)/Element/2D/StructLinTriangElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/StructQuadTriangElement.o : $(VPATH)/Element/2D/StructQuadTriangElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/QuadElement.o : $(VPATH)/Element/2D/QuadElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/StructLinQuadElement.o : $(VPATH)/Element/2D/StructLinQuadElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/StructQuadQuadElement.o : $(VPATH)/Element/2D/StructQuadQuadElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/StructElementList2D.o : $(VPATH)/Element/2D/StructElementList2D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/functionOnPoints.o : $(VPATH)/Load/functionOnPoints.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/functionOnLines.o : $(VPATH)/Load/functionOnLines.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/functionOnSurfaces.o : $(VPATH)/Load/functionOnSurfaces.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/PointLoad.o : $(VPATH)/Load/PointLoad.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/LineLoad.o : $(VPATH)/Load/LineLoad.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/SurfaceLoad.o : $(VPATH)/Load/SurfaceLoad.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/TemperatureLoad.o : $(VPATH)/Load/TemperatureLoad.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Load.o : $(VPATH)/Load/Load.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/FixDisplacement.o : $(VPATH)/BoundaryCondition/FixDisplacement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/StructBoundaryCondition1D.o : $(VPATH)/BoundaryCondition/StructBoundaryCondition1D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Domain.o : $(VPATH)/Domain/Domain.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/StructDomain.o : $(VPATH)/Domain/StructDomain.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Problem.o : $(VPATH)/Problem/Problem.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/StructProblem.o : $(VPATH)/Problem/StructProblem.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Solver.o : $(VPATH)/Solver/Solver.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/DataInput.o : $(VPATH)/DataIO/DataInput.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/DataOutput.o : $(VPATH)/DataIO/DataOutput.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

clean:
	rm -f $(BINDIR)/*.o main $(OBJECTDIR)/*.mod
