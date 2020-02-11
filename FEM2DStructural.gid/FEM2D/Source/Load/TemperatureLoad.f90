module TemperatureLoadMOD!En este momento el thermal y el structural tienen que tener la misma malla
  use tools
  use DebuggerMOD

  use StructElementList1DMOD
  use StructElementList2DMOD

  use Element1DPtrMOD
  use Element2DPtrMOD

  use MaterialMOD
  use MaterialPtrMOD
  use StructMaterialMOD

  use IntegratorPtrMOD
  implicit none
  private
  public :: TemperatureLoadTYPE, temperatureLoad
  type :: TemperatureLoadTYPE
     real(rkind)                            :: stableTemp
     real(rkind), dimension(:), allocatable :: temperature
     real(rkind), dimension(:), allocatable :: strain
   contains
     procedure, public :: init

     procedure, public :: apply

     procedure, private :: getTemperatureStrain
  end type TemperatureLoadTYPE

  interface temperatureLoad
     procedure :: constructor
  end interface temperatureLoad

contains

  type(TemperatureLoadTYPE) function constructor(stableTemp, temperature)
    implicit none
    real(rkind), intent(in) :: stableTemp
    real(rkind), dimension(:), intent(in) :: temperature
    call constructor%init(stableTemp, temperature)
  end function constructor

  subroutine init(this, stableTemp, temperature)
    implicit none
    class(TemperatureLoadTYPE), intent(inout) :: this
    real(rkind), intent(in) :: stableTemp
    real(rkind), dimension(:), intent(in) :: temperature
    allocate(this%temperature(size(temperature)))
    allocate(this%strain(size(temperature)))
    this%stableTemp = stableTemp
    this%temperature = temperature
  end subroutine init

  subroutine apply(this, elementList1D, elementList2D, rhs)
    implicit none
    class(TemperatureLoadTYPE), intent(inout) :: this
    type(StructElementList1DTYPE), intent(inout) :: elementList1D
    type(StructElementList2DTYPE), intent(inout) :: elementList2D
    real(rkind), dimension(:), intent(inout) :: rhs
    integer(ikind) :: iElem, nElem, iPoint, nPoint, nDof, k, pointID
    real(rkind) :: strain, int, int1, int2, bi, ci, d11, d12, d21, d22, thickness
    real(rkind), dimension(:,:,:), allocatable :: jacobian
    real(rkind), dimension(:), allocatable :: jacobianDet
    type(IntegratorPtrTYPE) :: integrator
    type(Element1DPtrTYPE) :: element1D
    type(Element2DPtrTYPE) :: element2D
    call this%getTemperatureStrain(elementList1D, elementList2D)

    nElem = elementList1D%getnElem()
    do iElem = 1, nElem
       element1D = elementList1D%getElement(iElem)
       integrator = element1D%getIntegrator()
       nPoint = element1D%getnPoint()
       allocate(jacobianDet(integrator%ptr%integTerms))
       do k = 1, integrator%ptr%integTerms
          jacobianDet(k) = element1D%jacobian(integrator%ptr%gPoint(k,1))
       end do
       do iPoint = 1, nPoint
          strain = this%strain(element1D%getPointID(iPoint))
          int = 0.d0
          do k = 1, integrator%ptr%integTerms
             int = int + element1D%ptr%material%ptr%young*integrator%ptr%dShapeFunc(k,1,iPoint)
          end do
          rhs(pointID) = rhs(pointID) + int
       end do
       deallocate(jacobianDet)
    end do
    if(allocated(jacobianDet)) deallocate(jacobianDet)
    nElem = elementList2D%getnElem()
    do iElem = 1, nElem
       element2D = elementList2D%getElement(iElem)
       integrator = element2D%getIntegrator()
       nPoint = element2D%getnPoint()
       nDof = element2D%getnDof()
       allocate(jacobian(integrator%ptr%integTerms,2,2))
       allocate(jacobianDet(integrator%ptr%integTerms))
       do k = 1, integrator%ptr%integTerms
          jacobian(k,1:2,1:2) = &
               element2D%jacobian(integrator%ptr%gPoint(k,1),integrator%ptr%gPoint(k,2))
          jacobianDet(k) = element2D%jacobianDet(jacobian(k,1:2,1:2))
       end do
       do iPoint = 1, nPoint
          pointID = element2D%getPointID(iPoint)
          strain = this%strain(pointID)
          int1 = 0.d0
          int2 = 0.d0
          do k = 1, integrator%ptr%integTerms
             bi = jacobian(k,2,2)*integrator%ptr%dShapeFunc(k,1,iPoint*nDof-1)  &
                  - jacobian(k,1,2)*integrator%ptr%dShapeFunc(k,2,iPoint*nDof-1)
             ci = jacobian(k,1,1)*integrator%ptr%dShapeFunc(k,2,iPoint*nDof-1)  &
                  -jacobian(k,2,1)*integrator%ptr%dShapeFunc(k,1,iPoint*nDof-1)
             d11 = element2D%ptr%material%ptr%d11
             d12 = element2D%ptr%material%ptr%d12
             d21 = element2D%ptr%material%ptr%d21
             d22 = element2D%ptr%material%ptr%d22
             thickness = element2D%ptr%material%ptr%thickness
             int1 = int1 + integrator%ptr%weight(k)*bi*(d11*strain+d12*strain)*thickness
             int2 = int2 + integrator%ptr%weight(k)*ci*(d12*strain+d22*strain)*thickness
          end do
          rhs(pointID*nDof-1) = rhs(pointID*nDof-1) + int1
          rhs(pointID*nDof) = rhs(pointID*nDof) + int2
          print'(A,I0,A,E16.8)', 'added temp term in row ', nDof*pointID-1, ' -> ', int1
          print'(A,I0,A,E16.8)', 'added temp term in row ', nDof*pointID, ' -> ', int2
       end do
       deallocate(jacobian)
       deallocate(jacobianDet)
    end do
    if(allocated(jacobian)) deallocate(jacobian)
    if(allocated(jacobianDet)) deallocate(jacobianDet)
  end subroutine apply
  
  subroutine getTemperatureStrain(this, elementList1D, elementList2D)
    implicit none
    class(TemperatureLoadTYPE), intent(inout) :: this
    type(StructElementList1DTYPE), intent(inout) :: elementList1D
    type(StructElementList2DTYPE), intent(inout) :: elementList2D
    integer(ikind) :: iElem, nElem, iPoint, nPoint
    real(rkind) :: alpha, temp
    type(Element1DPtrTYPE) :: element1D
    type(Element2DPtrTYPE) :: element2D
    nElem = elementList1D%getnElem()
    do iElem = 1, nElem
       element1D = elementList1D%getElement(iElem)
       nPoint = element1D%getnPoint()
       do iPoint = 1, nPoint
          alpha = element1D%ptr%material%ptr%thermalCoef
          temp = this%temperature(element1D%getPointID(iPoint))
          !Deformación plana
          this%strain = (1+element1D%ptr%material%ptr%poissonCoef)*alpha*(temp-this%stableTemp)
          !Tensión plana
          !this%strain(element1D%getPointID(iPoint)) = alpha*(temp-this%stableTemp)
       end do
    end do
    nElem = elementList2D%getnElem()
    do iElem = 1, nElem
       element2D = elementList2D%getElement(iElem)
       nPoint = element2D%getnPoint()
       do iPoint = 1, nPoint
          alpha = element2D%ptr%material%ptr%thermalCoef
          temp = this%temperature(element2D%getPointID(iPoint))
          !Deformación plana
          this%strain = (1+element2D%ptr%material%ptr%poissonCoef)*alpha*(temp-this%stableTemp)
          !Tensión plana
          !this%strain(element2D%getPointID(iPoint)) = alpha*(temp-this%stableTemp)
          !print'(A,I0,A,E16.8)', 'dT on node ', element2D%getPointID(iPoint), ' -> ', (temp-this%stableTemp)
          !print'(A,I0,A,E16.8)', 'strain on node ', element2D%getPointID(iPoint), ' -> ', this%strain(element2D%getPointID(iPoint))
       end do
    end do
  end subroutine getTemperatureStrain

end module TemperatureLoadMOD
