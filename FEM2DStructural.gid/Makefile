# Makefile

COMPILER	:=  ifort
VPATH		:=  $(PWD)/FEM2D/Source
BINDIR		:=  $(PWD)/FEM2D/Bin
OBJECTDIR	:=  $(PWD)/FEM2D/Objects
FFLAGS		:=  -Ofast -check bounds -traceback -module $(OBJECTDIR)
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
	$(COMPILER) $(VPATH)/Debugger.f90 \
	$(VPATH)/utilities.f90 precomp.f90 -o precomp

debug: $(OBJECTS)
	$(COMPILER) $(FFLAGS) $(OBJECTS) -o main

$(BINDIR)/main.o : main.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@ 

$(BINDIR)/%.o : $(VPATH)/%.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@ 
clean:
	rm -f $(BINDIR)/*.o main $(OBJECTDIR)/*.mod

cleanPrecomp:
	rm -f *.mod
