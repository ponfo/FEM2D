module StructProblemMOD
  use tools
  use DebuggerMOD

  use PointMOD
  use PointPtrMOD

  use Element1DPtrMOD
  use Element2DPtrMOD

  use SparseKit

  use DomainMOD
  use StructDomainMOD

  use ProblemMOD
  implicit none
  private
  public :: StructProblemTYPE, structProblem
  type, extends(ProblemTYPE) :: StructProblemTYPE
     type(StructDomainTYPE) :: domain
   contains
     procedure, public :: init

     procedure, public :: addPoint
     procedure, public :: addElement
     generic  , public :: addMaterial => addMaterialWithThermalCoef, addMaterialWithoutThermalCoef
     procedure, public :: addPointLoad
     procedure, public :: addLineLoad
     procedure, public :: addSurfaceLoad
     procedure, public :: addFixDisplacementX
     procedure, public :: addFixDisplacementY
     procedure, public :: setTemperatureLoad

     procedure, public :: setUp => assembleSystem

     procedure, private :: assembleStiffness
     procedure, private :: addMaterialWithThermalCoef
     procedure, private :: addMaterialWithoutThermalCoef
  end type StructProblemTYPE

  interface structProblem
     procedure constructor
  end interface structProblem

contains

  type(StructProblemTYPE) function constructor(nPoint, isQuadratic           &
       , nLine, nTriang, nQuad, nGauss, nMaterial, nPointLoad                &
       , nLineLoad, nSurfaceLoad, nFixDisplacementX, nFixDisplacementY       )
    implicit none
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: isQuadratic
    integer(ikind), intent(in) :: nLine
    integer(ikind), intent(in) :: nTriang
    integer(ikind), intent(in) :: nQuad
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: nMaterial
    integer(ikind), intent(in) :: nPointLoad
    integer(ikind), intent(in) :: nLineLoad
    integer(ikind), intent(in) :: nSurfaceLoad
    integer(ikind), intent(in) :: nFixDisplacementX
    integer(ikind), intent(in) :: nFixDisplacementY
    call constructor%init(nPoint, isQuadratic                                &
       , nLine, nTriang, nQuad, nGauss, nMaterial, nPointLoad                &
       , nLineLoad, nSurfaceLoad, nFixDisplacementX, nFixDisplacementY       )
  end function constructor

  subroutine init(this, nPoint, isQuadratic                                  &
       , nLine, nTriang, nQuad, nGauss, nMaterial, nPointLoad                &
       , nLineLoad, nSurfaceLoad, nFixDisplacementX, nFixDisplacementY       )
    implicit none
    class(StructProblemTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: isQuadratic
    integer(ikind), intent(in) :: nLine
    integer(ikind), intent(in) :: nTriang
    integer(ikind), intent(in) :: nQuad
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: nMaterial
    integer(ikind), intent(in) :: nPointLoad
    integer(ikind), intent(in) :: nLineLoad
    integer(ikind), intent(in) :: nSurfaceLoad
    integer(ikind), intent(in) :: nFixDisplacementX
    integer(ikind), intent(in) :: nFixDisplacementY
    integer(ikind)             :: nDof
    call debugLog('  Initiating Structural Problem')
    nDof = 2
    this%domain = structProblem(nPoint, isQuadratic                          &
       , nLine, nTriang, nQuad, nGauss, nMaterial, nPointLoad                &
       , nLineLoad, nSurfaceLoad, nFixDisplacementX, nFixDisplacementY       )
    if(isQuadratic == 0) then
       this%stiffness = &
            sparse(nnz = (nLine*4+nTriang*9+nQuad*16)*nDof, rows = nPoint*nDof)
    else if(isQuadratic == 1) then
       this%stiffness = &
            sparse(nnz = (nLine*9+nTriang*36+nQuad*64)*nDof, rows = nPoint*nDof)
    end if
    call debugLog('    Allocated Stiffness')
    call debugLog('      Estimated nnz: ', this%stiffness%getnnz())
    call debugLog('      Matrix order: ', this%stiffness%getn())
    allocate(this%rhs(nPoint*nDof))
    call debugLog('    Allocated RHS: ', size(this%rhs))
    allocate(this%dof(nPoint*nDof))
    call debugLog('    Allocated DOF: ', size(this%dof))
    this%dof = 0.d0
  end subroutine init

  subroutine addPoint(this, x, y, z)
    implicit none
    class(StructProblemTYPE), intent(inout) :: this
    real(rkind), intent(in) :: x
    real(rkind), intent(in) :: y
    real(rkind), intent(in) :: z
    call this%domain%addPoint(x, y, z)
  end subroutine addPoint

  subroutine addElement(this, type, nPoint, matID, pointList)
    implicit none
    class(StructProblemTYPE), intent(inout) :: this
    character(*), intent(in) :: type
    integer(ikind), intent(in) :: nPoint
    integer(ikind), intent(in) :: matID
    integer(ikind), dimension(nPoint), intent(in) :: pointList
    call this%domain%addElement(type, nPoint, matID, pointList)
  end subroutine addElement

  subroutine addMaterialWithoutThermalCoef(this, young, poissonCoef)
    implicit none
    class(StructProblemTYPE), intent(inout) :: this
    real(rkind), intent(in) :: young
    real(rkind), intent(in) :: poissonCoef
    call this%domain%addMaterial(this, young, poissonCoef)
  end subroutine addMaterialWithoutThermalCoef

  subroutine addMaterialWithThermalCoef(this, young, poissonCoef, thermalCoef)
    implicit none
    class(StructProblemTYPE), intent(inout) :: this
    real(rkind), intent(in) :: young
    real(rkind), intent(in) :: poissonCoef
    real(rkind), intent(in) :: thermalCoef
    call this%domain%addMaterial(this, young, poissonCoef, thermalCoef)
  end subroutine addMaterialWithThermalCoef

  subroutine addPointLoad(this, iPoint, iLoad)
    implicit none
    class(StructProblemTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iLoad
    call this%domain%addPointLoad(iPoint, iLoad)
  end subroutine addPointLoad

  subroutine addLineLoad(this, iPoint, iLoad)
    implicit none
    class(StructProblemTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iPoint
    integer(ikind), intent(in) :: iLoad
    call this%domain%addLineLoad(iPoint, iLoad)
  end subroutine addLineLoad

  subroutine addSurfaceLoad(this, iElem, iLoad)
    implicit none
    class(StructProblemTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iElem
    integer(ikind), intent(in) :: iLoad
    call this%domain%addSurfaceLoad(iElem, iLoad)
  end subroutine addSurfaceLoad

  subroutine addFixDisplacementX(this, id, value)
    implicit none
    class(StructProblemTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    call this%domain%addFixDisplacementX(id, value)
  end subroutine addFixDisplacementX

  subroutine addFixDisplacementY(this, id, value)
    implicit none
    class(StructProblemTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: id
    real(rkind), intent(in) :: value
    call this%domain%addFixDisplacementY(id, value)
  end subroutine addFixDisplacementY
    
  subroutine setTemperatureLoad(this, thermalDof)
    implicit none
    class(StructProblemTYPE), intent(inout) :: this
    real(rkind), dimension(this%domain%nPoint) :: thermalDof
    call this%domain%setTemperatureLoad(thermalDof)
  end subroutine setTemperatureLoad

  subroutine assembleSystem(this)
    implicit none
    class(StructProblemTYPE), intent(inout) :: this
    call debugLog('  Assembling stiffness matrix and right hand side vector')
    print'(A)', 'Assembling stiffness matrix and right hand side vector'
    call this%assembleStiffness()
    call this%domain%applySource(this%rhs)
    call this%domain%applyBC2D(this%stiffness, this%rhs)
    call this%domain%applyBC1D(this%stiffness, this%rhs)
  end subroutine assembleSystem

  subroutine assembleStiffness(this)
    implicit none
    class(StructProblemTYPE), intent(inout) :: this
    integer(ikind) :: i, j, iElem, nElem1D, nElem2D, nPoint, nDof
    real(rkind), dimension(:,:), allocatable :: localStiffness
    type(Element1DPtrTYPE) :: element1D
    type(Element2DPtrTYPE) :: element2D
    nDof = this%domain%nDof
    nElem1D = this%domain%nLine
    do iElem = 1, nElem1D
       element1D = this%domain%elementList1D%getElement(iElem)
       nPoint = element1D%getnPoint()
       allocate(localStiffness(nPoint*nDof,nPoint*nDof))
       localStiffness = element1D%getStiffness()
       do i = 1, nPoint
          do j = 1, nPoint
             do ii = 1, nDof
                do jj = 1, nDof
                   call this%stiffness%append(                                            &
                        val = localStiffness(i*nDof-(nDof-ii),j*nDof-(nDof-jj))          &
                        , row = element1D%getPointID(i)*nDof-(nDof-ii)                    &
                        , col = element1D%getPointID(j)*nDof-(nDof-jj)                    )
                end do
             end do
          end do
       end do
       deallocate(localStiffness)
    end do
    nElem2D = this%domain%nTriang + this%domain%nQuad
    do iElem = 1, nElem2D
       element2D = this%domain%elementList2D%getElement(iElem)
       nPoint = element2D%getnPoint()
       allocate(localStiffness(nPoint*nDof,nPoint*nDof))
       localStiffness = element2D%getStiffness()
       do i = 1, nPoint
          do j = 1, nPoint
             do ii = 1, nDof
                do jj = 1, nDof
                   call this%stiffness%append(                                            &
                        val = localStiffness(i*nDof-(nDof-ii),j*nDof-(nDof-jj))          &
                        , row = element2D%getPointID(i)*nDof-(nDof-ii)                    &
                        , col = element2D%getPointID(j)*nDof-(nDof-jj)                    )
                end do
             end do
          end do
       end do
       deallocate(localStiffness)
    end do
    call this%stiffness%makeCRS()
  end subroutine assembleStiffness
  
end module StructProblemMOD
