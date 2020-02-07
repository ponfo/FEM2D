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
  public :: StructProblemTYPE
  type, extends(ProblemTYPE) :: StructProblemTYPE
     type(StructDomainTYPE) :: domain
   contains
     procedure, public :: setUp => assembleSystem
     
     procedure, private :: assembleStiffness
  end type StructProblemTYPE

contains

  subroutine assembleSystem(this)
    implicit none
    class(StructProblemTYPE), intent(inout) :: this
    integer(ikind) :: i
    call debugLog('  Assembling stiffness matrix and right hand side vector')
    print'(A)', 'Assembling stiffness matrix and right hand side vector'
    call this%assembleStiffness()
    call this%domain%applyLoad(this%rhs)
    call this%domain%applyBC1D(this%stiffness, this%rhs)
!!$    print*, 'stiffness'
!!$    call this%stiffness%printNonZeros()
!!$    print*, 'rhs'
!!$    do i = 1, size(this%rhs)
!!$       print'(A,I0,A,E16.8)', 'rhs(', i, ') = ', this%rhs(i)
!!$    end do
  end subroutine assembleSystem

  subroutine assembleStiffness(this)
    implicit none
    class(StructProblemTYPE), intent(inout) :: this
    integer(ikind) :: i, j, ii, jj, iElem, nElem1D, nElem2D, nPoint, nDof
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
                        , row = element2D%getPointID(i)*nDof-(nDof-ii)                   &
                        , col = element2D%getPointID(j)*nDof-(nDof-jj)                   )
                end do
             end do
          end do
       end do
       deallocate(localStiffness)
    end do
    call this%stiffness%makeCRS()
  end subroutine assembleStiffness
  
end module StructProblemMOD
