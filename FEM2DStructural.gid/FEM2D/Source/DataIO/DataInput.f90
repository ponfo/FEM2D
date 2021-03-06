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
  integer(ikind)               :: nDirichletX
  integer(ikind)               :: nDirichletY
  integer(ikind)               :: nPressureLine
  integer(ikind)               :: nMaterial
  integer(ikind)               :: nGauss
  integer(ikind)               :: isQuadratic
  integer(ikind)               :: nLoadOP
  integer(ikind)               :: nLoadOL
  integer(ikind)               :: nLoadOS
  integer(ikind)               :: nPointLoad
  integer(ikind)               :: nLineLoad
  integer(ikind)               :: nSurfaceLoad
  integer(ikind)               :: isThereThermalCoupling
  real(rkind)                  :: stableTemp
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
    call debugLog('  Reading point and line Loads')
    call readPointLineSurfaceLoads(io)
    call debugLog('  Reading Boundary Conditions')
    call readBoundaryConditions(io)
    call debugLog('  Checking if there is thermal-structural coupling')
    call checkThermalCoupling(io)
    call debugLog('End loading data')
    close(project)
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
    read(project,*)  aux, nDirichletX
    read(project,*)  aux, nDirichletY
    read(project,*)  aux, nPressureLine
    read(project,*)  aux, nPointLoad
    read(project,*)  aux, nLoadOP
    read(project,*)  aux, nLineLoad
    read(project,*)  aux, nLoadOL
    read(project,*)  aux, nSurfaceLoad
    read(project,*)  aux, nLoadOS
    call debugLog('    Number of Elements.............................: ', nElem)
    call debugLog('    Are Elements Quadratic.........................: ', isQuadratic)
    call debugLog('    Number of linear elements......................: ', nLinearElem)
    call debugLog('    Number of Triangular elements..................: ', nTriangElem)
    call debugLog('    Number of Rectangular elements.................: ', nRectElem)
    call debugLog('    Number of Nodes................................: ', nPoint)
    call debugLog('    Number of Dirichlet X conditions...............: ', nDirichletX)    
    call debugLog('    Number of Dirichlet Y conditions...............: ', nDirichletY)
    call debugLog('    Number of Lines With Pressure..................: ', nPressureLine)
    call debugLog('    Number of Loads on points......................: ', nLoadOP)
    call debugLog('    Number of points with pointLoad................: ', nPointLoad)
    call debugLog('    Number of Loads on lines.......................: ', nLoadOL)
    call debugLog('    Number of points with lineLoad.................: ', nLineLoad)
    call debugLog('    Number of Loads on surfaces....................: ', nLoadOS)
    call debugLog('    Number of Surfaces with surfaceLoad............: ', nSurfaceLoad)
    call debugLog('    Number of Materials............................: ', nMaterial)
    call debugLog('    Gauss cuadrature order.........................: ', nGauss)
    call io%initStructProblem(nPoint, isQuadratic, nLinearElem, nTriangElem, nRectElem        &
         , nGauss, nMaterial, nPointLoad, nLineLoad, nSurfaceLoad, nPressureLine, nDirichletX &
         , nDirichletY)
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
    real(rkind) :: alpha, E, nu, A, t
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(A)', 'Material         alpha        E        nu       A       t    '
    do i = 1, nMaterial
       read(project,*) iMat, alpha, E, nu, A, t
       call io%addMaterial(E, nu, alpha, A, t)
       if(verbose) print'(4X,I0,7X,2(E10.3,3X))', iMat, alpha, E, nu, A, t 
    end do
  end subroutine initMaterials
  
  subroutine initElements(io)
    type(IODataTYPE), intent(inout) :: io
    integer(ikind) :: i, j, iElem, iMat, iNode, Conectivities(8)
    character(len=13) :: type
    Conectivities = 0
    do i = 1, 28+nLoadOP+nLoadOL+nLoadOS
       read(project,*)
    end do
    if(verbose) print'(A)', 'Element  |      Type      |  material index  |  nNodes  |  connectivities'
    do i = 1, nElem
       read(project,*) iElem, type, iMat, iNode, (Conectivities(j),j=1,iNode)
       if(verbose) print'(I5,A15,I18,I14,5X,*(I5,X))', iElem, type, iMat, iNode, (Conectivities(j),j=1,iNode)
       call io%addElement(type, iNode, iMat, Conectivities)
    end do
  end subroutine initElements
  
  subroutine readPointLineSurfaceLoads(io)
    implicit none
    type(IODataTYPE), intent(inout) :: io
    integer(ikind)                     :: i
    integer(ikind), dimension(:), allocatable :: iNode, iElem, iLoad
    allocate(iNode(max(nPointLoad, nLineLoad)))
    allocate(iElem(nSurfaceLoad))
    allocate(iLoad(max(nPointLoad, nLineLoad, nSurfaceLoad)))
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'pointLoads'
    if(verbose) print'(A)', 'Node    Load'
    do i = 1, nPointLoad
       read(project,*) iNode(i), iLoad(i)
       if(verbose) print'(I0,5X,I0)', iNode(i), iLoad(i)
       call io%addPointLoad(iNode(i), iLoad(i))
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'lineLoads'
    if(verbose) print'(A)', 'Node   Load'
    do i = 1, nLineLoad
       read(project,*) iNode(i), iLoad(i)
       if(verbose) print'(I0,5X,I0)', iNode(i), iLoad(i)
    end do
    if(isQuadratic == 0) then
       do i = 1, nLineLoad-1
          call io%addLineLoad(iNode(i:i+1), iLoad(i))
       end do
    else if(isQuadratic == 1) then
       do i = 1, nLineLoad-2, 2
          call io%addLineLoad(iNode(i:i+2), iLoad(i))
       end do
    end if
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'surfaceLoads'
    if(verbose) print'(A)', 'Element   Load'
    do i = 1, nSurfaceLoad
       read(project,*) iElem(i), iLoad(i)
       if(verbose) print'(I0,5X,I0)', iElem(i), iLoad(i)
       call io%addSurfaceLoad(iElem(i), iLoad(i))
    end do
  end subroutine readPointLineSurfaceLoads
  
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
    if(verbose) print'(/,A)', 'Dirichlet X conditions'
    if(verbose) print'(A)', 'Node    Value'
    do i = 1, nDirichletX
       read(Project,*) id, value
       if(verbose) print'(I0,5X,E10.3)', id, value
       call io%addFixDisplacementX(id, value)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Dirichlet Y conditions'
    if(verbose) print'(A)', 'Node    Value'
    do i = 1, nDirichletY
       read(Project,*) id, value
       if(verbose) print'(I0,5X,E10.3)', id, value
       call io%addFixDisplacementY(id, value)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(isQuadratic == 0) then
       nPointID = 2
    else if(isQuadratic == 1) then
       nPointID = 3
    end if
    allocate(pointID(nPointID))
    if(verbose) print'(/,A)', 'Pressure On Lines conditions'
    if(verbose) print'(A)', 'Elem    Nodes     Value'
    do i = 1, nPressureLine
       read(Project,*) elemID, (pointID(j),j=1,nPointID), value
       if(verbose) print*, elemID, (pointID(j),j=1,nPointID), value
       call io%addPressure(elemID, pointID, value)
    end do
  end subroutine readBoundaryConditions

  subroutine checkThermalCoupling(io)
    implicit none
    type(IODataTYPE), intent(inout) :: io
    read(project,*)
    read(project,*)
    read(project,*)
    read(project,*) aux, isThereThermalCoupling
    if(isThereThermalCoupling == 1) then
       read(project,*) aux, stableTemp
       call setTemperatureLoad(io)
    end if
  end subroutine checkThermalCoupling

  subroutine setTemperatureLoad(io)
    implicit none
    type(IODataTYPE), intent(inout) :: io
    integer(ikind) :: i, nPoint, auxInt
    real(rkind), dimension(:), allocatable :: temp
    open(8, file = 'Temperatures.dat', status = 'old')
    read(8,*)
    read(8,*)
    read(8,*) nPoint
    allocate(temp(nPoint))
    read(8,*)
    do i = 1, nPoint
       read(8,*) auxInt, temp(i)
    end do
    close(8)
    call io%setTemperatureLoad(stableTemp, temp)
  end subroutine setTemperatureLoad
  
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
    real(rkind) :: alpha, E, nu, A, t
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(A)', 'Material         alpha        E        nu       A       t    '
    do i = 1, nMaterial
       iMat = 1
       alpha = 0.0000016
       E = 210000000000
       nu = 0.3
       A = 1
       t = 1
       call io%addMaterial(alpha, E, nu, A, t)
       if(verbose) print'(4X,I0,7X,2(E10.3,3X),A)', iMat, alpha, E, nu, A, t, '*AUTO ASIGNED*'
    end do
  end subroutine autoAsignMaterial
  
  subroutine initElementsDefaultMat(io)
    implicit none
    type(IODataTYPE), intent(inout) :: io
    integer(ikind) :: i, j, iElem, iMat, iLoad, iNode, Conectivities(8), auxInt
    character(len=13) :: type
    iMat = 1
    do i = 1, 28+nLoadOP+nLoadOL+nLoadOS
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


