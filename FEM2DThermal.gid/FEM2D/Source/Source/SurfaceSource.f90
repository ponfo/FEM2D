module SurfaceSourceMOD
  use tools
  use FunctionOnSurfacesMOD
  use IntegratorMOD
  use IntegratorPtrMOD
  use ThElementList2DMOD
  use Element2DPtrMOD
  use PointPtrMOD
  implicit none
  private
  public :: SurfaceSourceTYPE, surfaceSource
  type SurfaceSourceTYPE
     integer(ikind) :: iElem
     integer(ikind) :: iSource
   contains
     procedure, public  :: init
     procedure, public  :: apply
     procedure, private :: setupIntegration
     procedure, private :: getValuedSource
  end type SurfaceSourceTYPE

  interface surfaceSource
     procedure :: constructor
  end interface surfaceSource

contains

  type(SurfaceSourceTYPE) function constructor(iElem, iSource)
    implicit none
    integer(ikind), intent(in) :: iElem
    integer(ikind), intent(in) :: iSource
    call constructor%init(iElem, iSource)
  end function constructor

  subroutine init(this, iElem, iSource)
    implicit none
    class(SurfaceSourceTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iElem
    integer(ikind), intent(in) :: iSource
    this%iElem = iElem
    this%iSource = iSource
  end subroutine init

  subroutine apply(this, elementList, rhs)
    implicit none
    class(SurfaceSourceTYPE), intent(inout) :: this
    type(ThElementList2DTYPE), intent(inout) :: elementList
    real(rkind), dimension(:), intent(inout) :: rhs
    type(Element2DPtrTYPE) :: element
    type(IntegratorPtrTYPE) :: integrator
    real(rkind), dimension(:), allocatable :: valuedSource
    real(rkind), dimension(:), allocatable :: jacobianDet
    integer(ikind) :: i, j, nPoint, pointID
    real(rkind) :: val
    element = elementList%getElement(this%iElem)
    integrator = element%getIntegrator()
    allocate(valuedSource(integrator%ptr%integTerms))
    allocate(jacobianDet(integrator%ptr%integTerms))
    call this%setupIntegration(element, integrator, valuedSource, jacobianDet)
    nPoint = element%getnPoint()
    do i = 1, nPoint
       val = 0
       do j = 1, integrator%ptr%integTerms
          val = val + integrator%ptr%weight(j)*integrator%ptr%shapeFunc(j,i) &
               *valuedSource(j)*jacobianDet(j)
       end do
       pointID = element%getPointID(i)
       rhs(pointID) = rhs(pointID) + val
    end do
    deallocate(valuedSource)
    deallocate(jacobianDet)
  end subroutine apply

  subroutine setupIntegration(this, element, integrator, valuedSource, jacobianDet)
    implicit none
    class(SurfaceSourceTYPE), intent(inout) :: this
    type(Element2DPtrTYPE), intent(inout) :: element
    type(IntegratorPtrTYPE), intent(in) :: integrator
    real(rkind), dimension(integrator%ptr%integTerms), intent(out) :: valuedSource
    real(rkind), dimension(integrator%ptr%integTerms), intent(out) :: jacobianDet
    integer(ikind) :: i
    real(rkind), dimension(2,2) :: jacobian
    valuedSource = this%getValuedSource(element, integrator)
    do i = 1, integrator%ptr%integTerms
       jacobian = element%jacobian(integrator%ptr%gPoint(i,1),integrator%ptr%gPoint(i,2))
       jacobianDet(i) = element%jacobianDet(jacobian)
    end do
  end subroutine setupIntegration

  function getValuedSource(this, element, integrator)
    implicit none
    class(SurfaceSourceTYPE), intent(inout) :: this
    type(Element2DPtrTYPE), intent(inout) :: element
    type(IntegratorPtrTYPE), intent(in) :: integrator
    real(rkind), dimension(integrator%ptr%integTerms) :: getValuedSource
    integer(ikind) :: i, j, nPoint
    real(rkind) :: x, y
    type(PointPtrTYPE), dimension(:), allocatable :: point
    nPoint = element%getnPoint()
    do i = 1, integrator%ptr%integTerms
       point = element%getPoint()
       x = 0
       y = 0
       do j = 1, nPoint
          x = x + integrator%ptr%shapeFunc(i,j)*point(j)%getx()
          y = y + integrator%ptr%shapeFunc(i,j)*point(j)%gety()
       end do
       getValuedSource(i) = funcOnSurfaces(this%iSource, x, y)
    end do
  end function getValuedSource

end module SurfaceSourceMOD
  
