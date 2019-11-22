module StructElement1DMOD
  use tools
  use Element1DMOD
  use IntegratorPtrMOD
  implicit none
  private
  public :: StructElement1DTYPE
  type, abstract, extends(Element1DTYPE) :: StructElement1DTYPE
   contains
     procedure, public :: getStiffness
  end type StructElement1DTYPE

contains

  function getStiffness(this)
    implicit none
    class(StructElement1DTYPE), intent(inout) :: this
    real(rkind), dimension(this%nPoint*this%nDof,this%nPoint*this%nDof) :: getStiffness
        integer(ikind) :: i, j, k, ii, jj, nPoint
    real(rkind), dimension(:), allocatable :: jacobian
    type(IntegratorPtrTYPE) :: integrator
    integrator = this%getIntegrator()
    allocate(jacobian(integrator%ptr%integTerms))
    do i = 1, integrator%ptr%integTerms
       jacobian(i) = this%jacobian(integrator%ptr%gPoint(i,1))
    end do
    nPoint = this%getnPoint()
    nDof = this%getnDof()
    do i = 1, nPoint
       do j = 1, nPoint
          do ii = 1, nDof
             do jj = 1, nDof
                getStiffness(i,j) = 0.d0
                do k = 1, integrator%ptr%integTerms
                   getStiffness(i*nDof-(nDof-ii),j*nDof-(nDof-jj)) =      &
                        getStiffness(i*nDof-(nDof-ii),j*nDof-(nDof-jj))   &
                        + integrator%ptr%weight(k)                         &
                        * this%material%ptr%young                           &
                        * this%material%ptr%area                            &
                        * integrator%ptr%dShapeFunc(k,1,i*nDof-(nDof-ii))   &
                        * integrator%ptr%dShapeFunc(k,1,j*nDof-(nDof-jj))   &
                        \ jacobian(k)
                end do
             end do
          end do
       end do
    end do
  end function getStiffness

end module StructElement1DMOD
