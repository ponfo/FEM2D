module ThElement1DMOD
  use tools
  use Element1DMOD
  use IntegratorPtrMOD
  implicit none
  private
  public :: ThElement1DTYPE
  type, abstract, extends(Element1DTYPE) :: ThElement1DTYPE
   contains
     procedure, public :: getStiffness
  end type ThElement1DTYPE

contains

  function getStiffness(this)
    implicit none
    class(ThElement1DTYPE), intent(inout) :: this
    real(rkind), dimension(this%nPoint*this%nDof,this%nPoint*this%nDof) :: getStiffness
    integer(ikind) :: i, j, k, nPoint
    real(rkind), dimension(:), allocatable :: jacobian
    type(IntegratorPtrTYPE) :: integrator
    integrator = this%getIntegrator()
    allocate(jacobian(integrator%ptr%integTerms))
    do i = 1, integrator%ptr%integTerms
       jacobian(i) = this%jacobian(integrator%ptr%gPoint(i,1))
    end do
    nPoint = this%getnPoint()
    do i = 1, nPoint
       do j = 1, nPoint
          getStiffness(i,j) = 0.d0
          do k = 1, integrator%ptr%integTerms
             getStiffness(i,j) = getStiffness(i,j)   &
                  + integrator%ptr%weight(k)          &
                  * this%material%ptr%conductivity(1)  &
                  * integrator%ptr%dShapeFunc(k,1,i)  &
                  * integrator%ptr%dShapeFunc(k,1,j)  &
                  / jacobian(k)
          end do
       end do
    end do
  end function getStiffness

end module ThElement1DMOD
  
