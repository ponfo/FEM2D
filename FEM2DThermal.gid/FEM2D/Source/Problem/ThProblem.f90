module ThProblemMOD
  use tools
  use DebuggerMOD

  use PointMOD
  use PointPtrMOD

  use Element1DPtrMOD
  use Element2DPtrMOD

  use SparseKit

  use DomainMOD
  use ThDomainMOD

  use ProblemMOD
  implicit none
  private
  public :: ThProblemTYPE, thermalProblem
  type, extends(ProblemTYPE) :: ThProblemTYPE
     type(ThDomainTYPE) :: domain
   contains
     procedure, public  :: init

     procedure, public  :: addPoint
     procedure, public  :: addElement
     procedure, public  :: addMaterial
     procedure, public  :: addPointSource
     procedure, public  :: addLineSource
     procedure, public  :: addSurfaceSource
     procedure, public  :: addDirichletPoint
     procedure, public  :: addNormalFluxPoint
     procedure, public  :: addNormalFluxLine
     procedure, public  :: addConvectionPoint
     procedure, public  :: addConvectionLine

     procedure, public  :: setUp => assembleSystem

     procedure, private :: assembleStiffness
  end type ThProblemTYPE

  interface thermalProblem
     procedure constructor
  end interface thermalProblem

contains
  
  type(ThProblemTYPE) function constructor(nPoint, isQuadratic           &
       , nLine, nTriang, nQuad, nGauss, nMaterial, nPointSource         &
       , nLineSource, nSurfaceSource, nDirichletPoint, nNormalFluxPoint &
       , nNormalFluxLine, nConvectionPoint, nConvectionLine             )
    implicit none
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: isQuadratic
    integer(ikind), intent(in) :: nLine
    integer(ikind), intent(in) :: nTriang
    integer(ikind), intent(in) :: nQuad
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: nMaterial
    integer(ikind), intent(in) :: nPointSource
    integer(ikind), intent(in) :: nLineSource
    integer(ikind), intent(in) :: nSurfaceSource
    integer(ikind), intent(in) :: nDirichletPoint
    integer(ikind), intent(in) :: nNormalFluxPoint
    integer(ikind), intent(in) :: nNormalFluxLine
    integer(ikind), intent(in) :: nConvectionPoint
    integer(ikind), intent(in) :: nConvectionLine
    call constructor%init(nPoint, isQuadratic                           &
       , nLine, nTriang, nQuad, nGauss, nMaterial, nPointSource         &
       , nLineSource, nSurfaceSource, nDirichletPoint, nNormalFluxPoint &
       , nNormalFluxLine, nConvectionPoint, nConvectionLine             )
  end function constructor

  subroutine init(this, nPoint, isQuadratic                             &
       , nLine, nTriang, nQuad, nGauss, nMaterial, nPointSource         &
       , nLineSource, nSurfaceSource, nDirichletPoint, nNormalFluxPoint &
       , nNormalFluxLine, nConvectionPoint, nConvectionLine             )
    implicit none
    class(ThProblemTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: isQuadratic
    integer(ikind), intent(in) :: nLine
    integer(ikind), intent(in) :: nTriang
    integer(ikind), intent(in) :: nQuad
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: nMaterial
    integer(ikind), intent(in) :: nPointSource
    integer(ikind), intent(in) :: nLineSource
    integer(ikind), intent(in) :: nSurfaceSource
    integer(ikind), intent(in) :: nDirichletPoint
    integer(ikind), intent(in) :: nNormalFluxPoint
    integer(ikind), intent(in) :: nNormalFluxLine
    integer(ikind), intent(in) :: nConvectionPoint
    integer(ikind), intent(in) :: nConvectionLine
    call debugLog('  Initiating Thermal Problem')
    this%domain = thermalDomain(nPoint, isQuadratic                     &
       , nLine, nTriang, nQuad, nGauss, nMaterial, nPointSource         &
       , nLineSource, nSurfaceSource, nDirichletPoint, nNormalFluxPoint &
       , nNormalFluxLine, nConvectionPoint, nConvectionLine             )
    if(isQuadratic == 0) then
       this%stiffness = &
            sparse(nnz = nLine*4+nTriang*9+nQuad*16, rows = nPoint)
    else if(isQuadratic == 1) then
       this%stiffness = &
            sparse(nnz = nLine*9+nTriang*36+nQuad*64, rows = nPoint)
    end if
    call debugLog('    Allocated Stiffness')
    call debugLog('      Estimated nnz: ', this%stiffness%getnnz())
    call debugLog('      Matrix order: ', this%stiffness%getn())
    allocate(this%rhs(nPoint))
    call debugLog('    Allocated RHS: ', size(this%rhs))
    allocate(this%dof(nPoint))
    call debugLog('    Allocated DOF: ', size(this%dof))
    this%dof = 0.d0
  end subroutine init

  subroutine addPoint(this, x, y, z)
    implicit none
    class(ThProblemTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), intent(in) :: z
    call this%domain%addPoint(x, y, z)
  end subroutine addPoint

  subroutine addElement(this, type, nPoint, matID, pointList)
    implicit none
    class(ThProblemTYPE), intent(inout) :: this
    character(*), intent(in) :: type
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: matID
    integer(ikind), dimension(nPoint), intent(in) :: pointList
    call this%domain%addElement(type, nPoint, matID, pointList)
  end subroutine addElement

  subroutine addMaterial(this, kx, ky)
    implicit none
    class(ThProblemTYPE), intent(inout) :: this
    real(rkind), intent(in) :: kx
    real(rkind), intent(in) :: ky
    call this%domain%addMaterial(kx, ky)
  end subroutine addMaterial

  subroutine addPointSource(this, iPoint, iSource)
    implicit none
    class(ThProblemTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iSource
    call this%domain%addPointSource(iPoint, iSource)
  end subroutine addPointSource

  subroutine addLineSource(this, iPoint, iSource)
    implicit none
    class(ThProblemTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iSource
    call this%domain%addLineSource(iPoint, iSource)
  end subroutine addLineSource

  subroutine addSurfaceSource(this, iElem, iSource)
    implicit none
    class(ThProblemTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iElem
    integer(ikind), intent(in) :: iSource
    call this%domain%addSurfaceSource(iElem, iSource)
  end subroutine addSurfaceSource
    
  subroutine addDirichletPoint(this, id, value)
    implicit none
    class(ThProblemTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    call this%domain%addDirichletPoint(id, value)
  end subroutine addDirichletPoint

  subroutine addNormalFluxPoint(this, id, value)
    implicit none
    class(ThProblemTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    call this%domain%addNormalFluxPoint(id, value)
  end subroutine addNormalFluxPoint

  subroutine addNormalFluxLine(this, elemID, pointID, value)
    implicit none
    class(ThProblemTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: elemID
    integer(ikind), dimension(:), intent(in) :: pointID
    real(rkind), intent(in) :: value
    call this%domain%addNormalFluxLine(elemID, pointID, value)
  end subroutine addNormalFluxLine

  subroutine addConvectionPoint(this, id, coef, temp)
    implicit none
    class(ThProblemTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: coef
    real(rkind), intent(in) :: temp
    call this%domain%addConvectionPoint(id, coef, temp)
  end subroutine addConvectionPoint
  
  subroutine addConvectionLine(this, elemID, pointID, coef, temp)
    implicit none
    class(ThProblemTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: elemID
    integer(ikind), dimension(:), intent(in) :: pointID
    real(rkind), intent(in) :: coef
    real(rkind), intent(in) :: temp
    call this%domain%addConvectionLine(elemID, pointID, coef, temp)
  end subroutine addConvectionLine

  subroutine assembleSystem(this)
    implicit none
    class(ThProblemTYPE), intent(inout) :: this
    call debugLog('  Assembling stiffness matrix and right hand side vector')
    print'(A)', 'Assembling stiffness matrix and right hand side vector'
    call this%assembleStiffness()
    call this%domain%applySource(this%rhs)
    call this%domain%applyBC2D(this%stiffness, this%rhs)
    call this%domain%applyBC1D(this%stiffness, this%rhs)
  end subroutine assembleSystem

  subroutine assembleStiffness(this)
    implicit none
    class(ThProblemTYPE), intent(inout) :: this
    integer(ikind) :: i, j, iElem, nElem1D, nElem2D, nPoint
    real(rkind), dimension(:,:), allocatable :: localStiffness
    type(Element1DPtrTYPE) :: element1D
    type(Element2DPtrTYPE) :: element2D
    nElem1D = this%domain%nLine
    do iElem = 1, nElem1D
       element1D = this%domain%elementList1D%getElement(iElem)
       nPoint = element1D%getnPoint()
       allocate(localStiffness(nPoint,nPoint))
       localStiffness = element1D%getStiffness()
       do i = 1, nPoint
          do j = 1, nPoint
             call this%stiffness%append(val = localStiffness(i,j) &
                  , row = element1D%getPointID(i)                 &
                  , col = element1D%getPointID(j)                 )
          end do
       end do
       deallocate(localStiffness)
    end do
    nElem2D = this%domain%nTriang + this%domain%nQuad
    do iElem = 1, nElem2D
       element2D = this%domain%elementList2D%getElement(iElem)
       nPoint = element2D%getnPoint()
       allocate(localStiffness(nPoint,nPoint))
       localStiffness = element2D%getStiffness()
       do i = 1, nPoint
          do j = 1, nPoint
             call this%stiffness%append(val = localStiffness(i,j) &
                  , row = element2D%getPointID(i)                 &
                  , col = element2D%getPointID(j)                 )
          end do
       end do
       deallocate(localStiffness)
    end do
    call this%stiffness%makeCRS()
  end subroutine assembleStiffness
  
end module ThProblemMOD

    