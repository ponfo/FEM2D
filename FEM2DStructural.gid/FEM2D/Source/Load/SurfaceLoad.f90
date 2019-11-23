module SurfaceLoadMOD
  use tools
  use FunctionOnSurfacesMOD
  use IntegratorMOD
  use IntegratorPtrMOD
  use StructElementList2DMOD
  use Element2DPtrMOD
  use PointPtrMOD
  implicit none
  private
  public :: SurfaceLoadTYPE, surfaceLoad
  type SurfaceLoadTYPE
     integer(ikind) :: iElem
     integer(ikind) :: iLoad
   contains
     procedure, public  :: init
     procedure, public  :: apply
     procedure, private :: setupIntegration
     procedure, private :: getValuedLoad
  end type SurfaceLoadTYPE

  interface surfaceLoad
     procedure :: constructor
  end interface surfaceLoad

contains

  type(SurfaceLoadTYPE) function constructor(iElem, iLoad)
    implicit none
    integer(ikind), intent(in) :: iElem
    integer(ikind), intent(in) :: iLoad
    call constructor%init(iElem, iLoad)
  end function constructor

  subroutine init(this, iElem, iLoad)
    implicit none
    class(SurfaceLoadTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: iElem
    integer(ikind), intent(in) :: iLoad
    this%iElem = iElem
    this%iLoad = iLoad
  end subroutine init

  subroutine apply(this, elementList, rhs)
    implicit none
    class(SurfaceLoadTYPE), intent(inout) :: this
    type(StructElementList2DTYPE), intent(inout) :: elementList
    real(rkind), dimension(:), intent(inout) :: rhs
    type(Element2DPtrTYPE) :: element
    type(IntegratorPtrTYPE) :: integrator
    real(rkind), dimension(:,:), allocatable :: valuedLoad
    real(rkind), dimension(:), allocatable :: jacobianDet
    integer(ikind) :: i, j, nPoint, pointID
    real(rkind), dimension(2) :: val
    element = elementList%getElement(this%iElem)
    integrator = element%getIntegrator()
    allocate(valuedLoad(2, integrator%ptr%integTerms))
    allocate(jacobianDet(integrator%ptr%integTerms))
    call this%setupIntegration(element, integrator, valuedLoad, jacobianDet)
    nPoint = element%getnPoint()
    do i = 1, nPoint
       val = 0
       do j = 1, integrator%ptr%integTerms
          val(1) = val(1) + integrator%ptr%weight(j)*integrator%ptr%shapeFunc(j,2*i-1) &
               *valuedLoad(1,j)*jacobianDet(j)
          val(2) = val(2) + integrator%ptr%weight(j)*integrator%ptr%shapeFunc(j,2*i)   &
               *valuedLoad(2,j)*jacobianDet(j)
       end do
       pointID = element%getPointID(i)
       rhs(2*pointID-1) = rhs(2*pointID-1) + val(1)
       rhs(2*pointID) = rhs(2*pointID) + val(2)
    end do
    deallocate(valuedLoad)
    deallocate(jacobianDet)
  end subroutine apply

  subroutine setupIntegration(this, element, integrator, valuedLoad, jacobianDet)
    implicit none
    class(SurfaceLoadTYPE), intent(inout) :: this
    type(Element2DPtrTYPE), intent(inout) :: element
    type(IntegratorPtrTYPE), intent(in) :: integrator
    real(rkind), dimension(2, integrator%ptr%integTerms), intent(out) :: valuedLoad
    real(rkind), dimension(integrator%ptr%integTerms), intent(out) :: jacobianDet
    integer(ikind) :: i
    real(rkind), dimension(2,2) :: jacobian
    valuedLoad = this%getValuedLoad(element, integrator)
    do i = 1, integrator%ptr%integTerms
       jacobian = element%jacobian(integrator%ptr%gPoint(i,1),integrator%ptr%gPoint(i,2))
       jacobianDet(i) = element%jacobianDet(jacobian)
    end do
  end subroutine setupIntegration

  function getValuedLoad(this, element, integrator)
    implicit none
    class(SurfaceLoadTYPE), intent(inout) :: this
    type(Element2DPtrTYPE), intent(inout) :: element
    type(IntegratorPtrTYPE), intent(in) :: integrator
    real(rkind), dimension(2, integrator%ptr%integTerms) :: getValuedLoad
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
       getValuedLoad(1:2,i) = funcOnSurfaces(this%iLoad, x, y)
    end do
  end function getValuedLoad

end module SurfaceLoadMOD
  
