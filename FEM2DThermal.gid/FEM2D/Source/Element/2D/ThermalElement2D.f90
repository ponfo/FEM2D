module ThermalElement2DMOD
  use tools
  use Element2DMOD
  use IntegratorPtrMOD
  implicit none
  private
  public :: ThermalElement2DTYPE
  type, abstract, extends(Element2DTYPE) :: ThermalElement2DTYPE
   contains
     procedure, public :: getStiffness
  end type ThermalElement2DTYPE

contains

  function getStiffness(this)
    implicit none
    class(ThermalElement2DTYPE), intent(inout) :: this
    real(rkind), dimension(this%nPoint*this%nDof,this%nPoint*this%nDof) :: getStiffness
    integer(ikind) :: i, j, k, nPoint
    real(rkind) :: bi, bj, ci, cj
    real(rkind), dimension(:,:,:), allocatable :: jacobian
    real(rkind), dimension(:), allocatable :: jacobianDet
    type(IntegratorPtrTYPE) :: integrator
    integrator = this%getIntegrator()
    allocate(jacobian(integrator%ptr%integTerms,2,2))
    allocate(jacobianDet(integrator%ptr%integTerms))
    do i = 1, integrator%ptr%integTerms
       jacobian(i,1:2,1:2) = this%jacobian(integrator%ptr%gPoint(i,1),integrator%ptr%gPoint(i,2))
       jacobianDet(i) = this%jacobianDet(jacobian(i,1:2,1:2))
    end do
    nPoint = this%getnPoint()
    do i = 1, nPoint
       do j = 1, nPoint
          getStiffness(i,j) = 0.d0
          do k = 1, integrator%ptr%integTerms
             bi = jacobian(k,2,2)*integrator%ptr%dShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%ptr%dShapeFunc(k,2,i)
             bj = jacobian(k,2,2)*integrator%ptr%dShapeFunc(k,1,j) &
                  - jacobian(k,1,2)*integrator%ptr%dShapeFunc(k,2,j)
             ci = jacobian(k,1,1)*integrator%ptr%dShapeFunc(k,2,i) &
                  -jacobian(k,2,1)*integrator%ptr%dShapeFunc(k,1,i)
             cj = jacobian(k,1,1)*integrator%ptr%dShapeFunc(k,2,j) &
                  -jacobian(k,2,1)*integrator%ptr%dShapeFunc(k,1,j)
             getStiffness(i,j) = getStiffness(i,j)          &
                  + integrator%ptr%weight(k)                 &
                  *(this%material%ptr%conductivity(1)*bi*bj  &
                  + this%material%ptr%conductivity(2)*ci*cj) &
                  / jacobianDet(k)
          end do
       end do
    end do
  end function getStiffness

end module ThermalElement2DMOD
