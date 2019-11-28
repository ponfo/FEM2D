module StructElement2DMOD
  use tools
  use Element2DMOD
  use IntegratorPtrMOD
  implicit none
  private
  public :: StructElement2DTYPE
  type, abstract, extends(Element2DTYPE) :: StructElement2DTYPE
   contains
     procedure, public :: getStiffness
  end type StructElement2DTYPE

contains

  function getStiffness(this)
    implicit none
    class(StructElement2DTYPE), intent(inout) :: this
    real(rkind), dimension(this%nPoint*this%nDof,this%nPoint*this%nDof) :: getStiffness
    integer(ikind) :: i, j, k, nPoint, nDof, ii, jj
    real(rkind) :: bi, bj, ci, cj
    real(rkind), dimension(:,:,:), allocatable :: jacobian
    real(rkind), dimension(:), allocatable :: jacobianDet
    real(rkind), dimension(2,2) :: Kij
    type(IntegratorPtrTYPE) :: integrator
    integrator = this%getIntegrator()
    allocate(jacobian(integrator%ptr%integTerms,2,2))
    allocate(jacobianDet(integrator%ptr%integTerms))
    do i = 1, integrator%ptr%integTerms
       jacobian(i,1:2,1:2) = this%jacobian(integrator%ptr%gPoint(i,1),integrator%ptr%gPoint(i,2))
       jacobianDet(i) = this%jacobianDet(jacobian(i,1:2,1:2))
    end do
    nPoint = this%getnPoint()
    nDof = this%getnDof()
    do i = 1, nPoint
       do j = 1, nPoint
          ii = nDof*i-1
          jj = nDof*j-1
          getStiffness(ii,jj)     = 0.d0
          getStiffness(ii+1,jj)   = 0.d0
          getStiffness(ii,jj+1)   = 0.d0
          getStiffness(ii+1,jj+1) = 0.d0
          do k = 1, integrator%ptr%integTerms
             bi = jacobian(k,2,2)*integrator%ptr%dShapeFunc(k,1,ii) &
                  - jacobian(k,1,2)*integrator%ptr%dShapeFunc(k,2,ii)
             bj = jacobian(k,2,2)*integrator%ptr%dShapeFunc(k,1,jj) &
                  - jacobian(k,1,2)*integrator%ptr%dShapeFunc(k,2,jj)
             ci = jacobian(k,1,1)*integrator%ptr%dShapeFunc(k,2,ii) &
                  -jacobian(k,2,1)*integrator%ptr%dShapeFunc(k,1,ii)
             cj = jacobian(k,1,1)*integrator%ptr%dShapeFunc(k,2,jj) &
                  -jacobian(k,2,1)*integrator%ptr%dShapeFunc(k,1,jj)
             Kij(1,1) = bi*bj*this%material%ptr%d11 + ci*cj*this%material%ptr%d33
             Kij(1,2) = bi*cj*this%material%ptr%d12 + bj*ci*this%material%ptr%d33
             Kij(2,1) = ci*bj*this%material%ptr%d21 + bi*cj*this%material%ptr%d33
             Kij(2,2) = bi*bj*this%material%ptr%d33 + ci*cj*this%material%ptr%d22                   
             getStiffness(ii,jj) = getStiffness(ii,jj)                    &
                  + integrator%ptr%weight(k)*Kij(1,1)/jacobianDet(k)
             getStiffness(ii,jj+1) = getStiffness(ii,jj+1)                &
                  + integrator%ptr%weight(k)*Kij(1,2)/jacobianDet(k)
             getStiffness(ii+1,jj) = getStiffness(ii+1,jj)                &
                  + integrator%ptr%weight(k)*Kij(2,1)/jacobianDet(k)
             getStiffness(ii+1,jj+1) = getStiffness(ii+1,jj+1)            &
                  + integrator%ptr%weight(k)*Kij(2,2)/jacobianDet(k)
          end do
!!$          print*, 'i -> ', i
!!$          print*, 'j -> ', j
!!$          print*, 'ii -> ', ii
!!$          print*, 'jj -> ', jj
!!$          print*, 'K(i,j) = ', getStiffness(ii,jj)
!!$          print*, 'K(i,j+1) = ', getStiffness(ii,jj+1)
!!$          print*, 'K(i+1,j) = ', getStiffness(ii+1,jj)
!!$          print*, 'K(i+1,j+1) = ', getStiffness(ii+1,jj+1)
       end do
    end do
    getStiffness = getStiffness * this%material%ptr%thickness
    print*, 'stiffness for element -> ', this%id
    do i = 1, size(getStiffness,1)
       do j = 1, size(getStiffness,2)
          print'(A,I0,A,I0,A,E16.8)', 'stiffness(', i, ',', j, ') = ', getStiffness(i,j)
       end do
    end do
  end function getStiffness

end module StructElement2DMOD
