module DataInputMOD
  use tools
  use DebuggerMOD
  use IODataMOD
  implicit none
  private
  public :: initFEM2D
  integer(ikind), parameter    :: projectData = 1
  integer(ikind), parameter    :: project = 2
  integer(ikind), parameter    :: functions = 4
  integer(ikind), dimension(8) :: date_time
  integer(ikind)               :: nElem
  integer(ikind)               :: nLinearElem
  integer(ikind)               :: nTriangElem
  integer(ikind)               :: nRectElem
  integer(ikind)               :: nPoint
  integer(ikind)               :: iPoint
  integer(ikind)               :: nNormalFluxOP
  integer(ikind)               :: nNormalFluxOL
  integer(ikind)               :: nDirichlet
  integer(ikind)               :: nMaterial
  integer(ikind)               :: nGauss
  integer(ikind)               :: nConvectionOP
  integer(ikind)               :: nConvectionOL
  integer(ikind)               :: isQuadratic
  integer(ikind)               :: nSourceOP
  integer(ikind)               :: nSourceOL
  integer(ikind)               :: nSourceOS
  integer(ikind)               :: nPointSource
  integer(ikind)               :: nLineSource
  integer(ikind)               :: nSurfaceSource
  character(100)               :: projectName
  character(100)               :: path
  character(100)               :: aux
  logical       , parameter    :: verbose = .false.
  logical                      :: isMaterialAsigned = .true.
  interface initFEM2D
     procedure :: initFEM2D
  end interface initFEM2D
contains
  subroutine initFEM2D(io)
    implicit none
    type(IODataTYPE), intent(inout) :: io
    call initLog(.true., 'log.dat')
    call debugLog('  Reading project data')
    call readProjectData
    call debugLog('  Reading mesh data')
    call initMesh(io)
    if(isMaterialAsigned) then
       call debugLog('  Reading materials properties')
       call initMaterials(io)
       call debugLog('  Reading elements')
       call initElements(io)
    else
       call debugLog('  Auto asigning properties')
       call autoAsignMaterial(io)
       call debugLog('  Reading elements')
       call initElementsDefaultMat(io)
    end if
    call debugLog('  Reading point and line Sources')
    call readPointLineSurfaceSources(io)
    call debugLog('  Reading Boundary Conditions')
    call readBoundaryConditions(io)
    call debugLog('End loading data')
  end subroutine initFEM2D
  
  subroutine readProjectData
    implicit none
    open(projectData, file = 'projectData.dat')
    read(projectData, '(*(A))') projectName
    read(projectData, '(*(A))') path
    close(projectData)
    call debugLog('    Project name: ', trim(projectName))
    call debugLog('    Path: ', trim(path))
  end subroutine readProjectData
  
  subroutine initMesh(io)
    implicit none
    type(IODataTYPE), intent(inout) :: io
    integer(ikind) :: i
    real(rkind)    :: x, y, z
    open(project, file = trim(projectName)//'.dat')
    do i = 1, 9
       read(project,*)
    end do
    read(project,*)  aux, nElem
    read(project,*)  aux, nPoint
    read(project,*)  aux, isQuadratic
    read(project,*)  aux, nLinearElem
    read(project,*)  aux, nTriangElem
    read(project,*)  aux, nRectElem
    read(project,*)  aux, nMaterial
    call checknMaterial(nMaterial)
    read(project,*)  aux, nGauss    
    read(project,*)  aux, nDirichlet
    read(project,*)  aux, nNormalFluxOP
    read(project,*)  aux, nConvectionOP
    read(project,*)  aux, nNormalFluxOL
    read(project,*)  aux, nConvectionOL
    read(project,*)  aux, nPointSource
    read(project,*)  aux, nSourceOP
    read(project,*)  aux, nLineSource
    read(project,*)  aux, nSourceOL
    read(project,*)  aux, nSurfaceSource
    read(project,*)  aux, nSourceOS
    call debugLog('    Number of Elements.............................: ', nElem)
    call debugLog('    Are Elements Quadratic.........................: ', isQuadratic)
    call debugLog('    Number of linear elements......................: ', nLinearElem)
    call debugLog('    Number of Triangular elements..................: ', nTriangElem)
    call debugLog('    Number of Rectangular elements.................: ', nRectElem)
    call debugLog('    Number of Nodes................................: ', nPoint)
    call debugLog('    Number of Dirichlet conditions.................: ', nDirichlet)
    call debugLog('    Number of NormalFluxOnPoints conditions........: ', nNormalFluxOP)
    call debugLog('    Number of NormalFluxOnLines conditions.........: ', nNormalFluxOL)    
    call debugLog('    Number of ConvectionOnPoints conditions........: ', nConvectionOP)
    call debugLog('    Number of ConvectionOnLines conditions.........: ', nConvectionOL)    
    call debugLog('    Number of Sources on points....................: ', nSourceOP)
    call debugLog('    Number of points with pointSource..............: ', nPointSource)
    call debugLog('    Number of Sources on lines.....................: ', nSourceOL)
    call debugLog('    Number of points with lineSource...............: ', nLineSource)
    call debugLog('    Number of Sources on surfaces..................: ', nSourceOS)
    call debugLog('    Number of Surfaces with surfaceSource..........: ', nSurfaceSource)
    call debugLog('    Number of Materials............................: ', nMaterial)
    call debugLog('    Gauss cuadrature order.........................: ', nGauss)
    call io%initThermalProblem(nPoint, isQuadratic, nLinearElem, nTriangElem, nRectElem    &
         , nGauss, nMaterial, nPointSource, nLineSource, nSurfaceSource, nDirichlet        &
         , nNormalFluxOP, nNormalFluxOL, nConvectionOP, nConvectionOL                      )
    do i = 1, 6
       read(project,*)
    end do
    if(verbose) print'(A)', 'Coordinates:'
    if(verbose) print'(1X,A)', 'NODES      X          Y          Z'
    do i = 1, nPoint
       read(project,*) iPoint, x, y, z
       if(verbose) print'(3X,I0,3X,E10.3,3X,E10.3,3X,E10.3)', iPoint, x, y, z
       call io%addPoint(x, y, z)
    end do
  end subroutine initMesh
  
  subroutine initMaterials(io)
    implicit none
    type(IODataTYPE), intent(inout) :: io
    integer(ikind) :: i, iMat
    real(rkind) :: kx, ky
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(A)', 'Material         Kx            Ky    '
    do i = 1, nMaterial
       read(project,*) iMat, kx, ky
       call io%addMaterial(kx, ky)
       if(verbose) print'(4X,I0,7X,2(E10.3,3X))', iMat, kx, ky 
    end do
  end subroutine initMaterials
  
  subroutine initElements(io)
    type(IODataTYPE), intent(inout) :: io
    integer(ikind) :: i, j, iElem, iMat, iNode, Conectivities(8)
    character(len=13) :: type
    Conectivities = 0
    do i = 1, 28+nSourceOP+nSourceOL+nSourceOS
       read(project,*)
    end do
    if(verbose) print'(A)', 'Element  |      Type      |  material index  |  nNodes  |  connectivities'
    do i = 1, nElem
       read(project,*) iElem, type, iMat, iNode, (Conectivities(j),j=1,iNode)
       if(verbose) print'(I5,A15,I18,I14,5X,*(I5,X))', iElem, type, iMat, iNode, (Conectivities(j),j=1,iNode)
       call io%addElement(type, iNode, iMat, Conectivities)
    end do
  end subroutine initElements
  
  subroutine readPointLineSurfaceSources(io)
    implicit none
    type(IODataTYPE), intent(inout) :: io
    integer(ikind)                     :: i
    integer(ikind), dimension(:), allocatable :: iNode, iElem, iSource
    allocate(iNode(max(nPointSource, nLineSource)))
    allocate(iElem(nSurfaceSource))
    allocate(iSource(max(nPointSource, nLineSource, nSurfaceSource)))
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'pointSources'
    if(verbose) print'(A)', 'Node    Source'
    do i = 1, nPointSource
       read(project,*) iNode(i), iSource(i)
       if(verbose) print'(I0,5X,I0)', iNode(i), iSource(i)
       call io%addPointSource(iNode(i), iSource(i))
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'lineSources'
    if(verbose) print'(A)', 'Node   Source'
    do i = 1, nLineSource
       read(project,*) iNode(i), iSource(i)
       if(verbose) print'(I0,5X,I0)', iNode(i), iSource(i)
    end do
    if(isQuadratic == 0) then
       do i = 1, nLineSource-1
          call io%addLineSource(iNode(i:i+1), iSource(i))
       end do
    else if(isQuadratic == 1) then
       do i = 1, nLineSource-2, 2
          call io%addLineSource(iNode(i:i+2), iSource(i))
       end do
    end if
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'surfaceSources'
    if(verbose) print'(A)', 'Element   Source'
    do i = 1, nSurfaceSource
       read(project,*) iElem(i), iSource(i)
       if(verbose) print'(I0,5X,I0)', iElem(i), iSource(i)
       call io%addSurfaceSource(iElem(i), iSource(i))
    end do
  end subroutine readPointLineSurfaceSources

  subroutine readBoundaryConditions(io)
    implicit none
    type(IODataTYPE), intent(inout) :: io
    integer(ikind)                   :: i, j, id, elemID, nPointID, iPoint
    integer(ikind), dimension(:), allocatable :: pointID
    real(rkind)                      :: value
    real(rkind)                      :: coef, temp
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Dirichlet conditions'
    if(verbose) print'(A)', 'Node    Value'
    do i = 1, nDirichlet
       read(Project,*) id, value
       if(verbose) print'(I0,5X,E10.3)', id, value
       call io%addDirichletPoint(id, value)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(isQuadratic == 0) then
       nPointID = 2
    else
       nPointID = 3
    end if
    allocate(pointID(nPointID))
    if(verbose) print'(/,A)', 'Normal Flux On Points conditions'
    if(verbose) print'(A)', 'Nodes     Value'
    do i = 1, nNormalFluxOP
       read(Project,*) iPoint, value
       if(verbose) print*, iPoint, value
       call io%addNormalFluxPoint(iPoint, value)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Normal Flux On Lines conditions'
    if(verbose) print'(A)', 'Elem    Nodes     Value'
    do i = 1, nNormalFluxOL
       read(Project,*) elemID, (pointID(j),j=1,nPointID), value
       if(verbose) print*, elemID, (pointID(j),j=1,nPointID), value
       call io%addNormalFluxLine(elemID, pointID, value)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Convection On Points conditions'
    if(verbose) print'(A)', 'Nodes     Coef     Temp'
    do i = 1, nConvectionOP
       read(Project,*) iPoint, coef, temp
       if(verbose) print*, iPoint, coef, temp
       call io%addConvectionPoint(iPoint, coef, temp)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Convection On Lines conditions'
    if(verbose) print'(A)', 'Elem    Nodes     Coef     Temp'
    do i = 1, nConvectionOL
       read(Project,*) elemID, (pointID(j),j=1,nPointID), coef, temp
       if(verbose) print*, elemID, (pointID(j),j=1,nPointID), coef, temp
       call io%addConvectionLine(elemID, pointID, coef, temp)
    end do
    close(project)
  end subroutine readBoundaryConditions

  subroutine checknMaterial(nMaterial)
    implicit none
    integer(ikind), intent(inout) :: nMaterial
    if(nMaterial == 0) then
       print*, '*************************************************************'
       print*, '**  Material not asigned, auto asigning values equal to 1  **'
       print*, '*************************************************************'
       nMaterial = 1
       isMaterialAsigned = .false.
    end if
  end subroutine checknMaterial

  subroutine autoAsignMaterial(io)
    implicit none
    type(IODataTYPE), intent(inout) :: io
    integer(ikind) :: i, iMat
    real(rkind) :: kx, ky
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(A)', 'Material         Kx            Ky    '
    do i = 1, nMaterial
       iMat = 1
       kx = 1
       ky = 1
       call io%addMaterial(kx, ky)
       if(verbose) print'(4X,I0,7X,2(E10.3,3X),A)', iMat, kx, ky, '*AUTO ASIGNED*'
    end do
  end subroutine autoAsignMaterial

  subroutine initElementsDefaultMat(io)
    implicit none
    type(IODataTYPE), intent(inout) :: io
    integer(ikind) :: i, j, iElem, iMat, iSource, iNode, Conectivities(8), auxInt
    character(len=13) :: type
    iMat = 1
    do i = 1, 28+nSourceOP+nSourceOL+nSourceOS
       read(project,*)
    end do
    if(verbose) print'(A)', 'Element  |      Type      |  material index  |  nNodes  |  connectivities'
    do i = 1, nElem
       read(project,*) iElem, type, iMat, iNode, (Conectivities(j),j=1,iNode)
       iMat = 1
       if(verbose) print'(I5,A15,I18,I14,5X,*(I5,X))', iElem, type, iMat, iNode, (Conectivities(j),j=1,iNode)
       call io%addElement(type, iNode, iMat, Conectivities)
    end do
  end subroutine initElementsDefaultMat
end module DataInputMOD


