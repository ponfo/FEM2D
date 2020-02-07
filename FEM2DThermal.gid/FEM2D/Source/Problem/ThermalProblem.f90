module ThermalProblemMOD
  use tools
  use DebuggerMOD

  use PointMOD
  use PointPtrMOD

  use Element1DPtrMOD
  use Element2DPtrMOD

  use SparseKit

  use DomainMOD
  use ThermalDomainMOD

  use ProblemMOD
  implicit none
  private
  public :: ThermalProblemTYPE
  type, extends(ProblemTYPE) :: ThermalProblemTYPE
     type(ThermalDomainTYPE) :: domain
   contains
     procedure, public  :: setUp => assembleSystem

     procedure, private :: assembleStiffness
  end type ThermalProblemTYPE
  
contains

  subroutine assembleSystem(this)
    implicit none
    class(ThermalProblemTYPE), intent(inout) :: this
    call debugLog('  Assembling stiffness matrix and right hand side vector')
    print'(A)', 'Assembling stiffness matrix and right hand side vector'
    call this%assembleStiffness()
    call this%domain%applySource(this%rhs)
    call this%domain%applyBC2D(this%stiffness, this%rhs)
    call this%domain%applyBC1D(this%stiffness, this%rhs)
  end subroutine assembleSystem

  subroutine assembleStiffness(this)
    implicit none
    class(ThermalProblemTYPE), intent(inout) :: this
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
  
end module ThermalProblemMOD

    
