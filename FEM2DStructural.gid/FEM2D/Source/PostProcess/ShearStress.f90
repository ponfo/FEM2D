module ShearStressMOD
  use tools
  use DebuggerMOD

  use PointMOD
  use PointPtrMOD

  use MaterialPtrMOD
  use StructMaterialMOD

  use IntegratorMOD
  use IntegratorPtrMOD

  use Element1DPtrMOD
  use Element2DPtrMOD

  use StructProblemMOD
  implicit none
  private
  public :: ShearStressTYPE
  type ShearStressTYPE
     integer(ikind) , dimension(:)  , allocatable :: lineElemID
     real(rkind)    , dimension(:)  , allocatable :: lineShS
     real(rkind)    , dimension(:)  , allocatable :: lineGPoint
     integer(ikind) , dimension(:)  , allocatable :: triangElemID
     real(rkind)    , dimension(:)  , allocatable :: triangShS
     real(rkind)    , dimension(:,:), allocatable :: triangGPoint
     integer(ikind) , dimension(:)  , allocatable :: quadElemID
     real(rkind)    , dimension(:)  , allocatable :: quadShS
     real(rkind)    , dimension(:,:), allocatable :: quadGPoint
   contains
     procedure, public  :: calculateShearStress
     procedure, private :: valueGPoints1D
     procedure, private :: valueGPointsTriang
     procedure, private :: valueGPointsQuad
  end type ShearStressTYPE

  procedure(addTriangShearStress), pointer :: addShearStress => null()
  integer(ikind) :: lineCount
  integer(ikind) :: triangCount
  integer(ikind) :: quadCount

contains

  subroutine calculateShearStress(this, problem)
    implicit none
    class(ShearStressTYPE), intent(inout) :: this
    class(StructProblemTYPE), intent(inout) :: problem
    integer(ikind) :: iElem, iGauss, i, triangElemCount, quadElemCount
    integer(ikind) :: nElem, nPoint, nLine, nTriang, nQuad, nDof
    integer(ikind) :: nGauss, nGaussPointLine, nGaussPointTriang, nGaussPointQuad
    real(rkind) :: xi, eta, x, y,  bi, ci, dNidx, dNidy, jacobianDet
    real(rkind) :: k, d33, ShS, jacobian1D
    real(rkind), dimension(:,:), allocatable :: dsf
    real(rkind), dimension(2,2) :: jacobian
    type(Element1DPtrTYPE) :: element1D
    type(Element2DPtrTYPE) :: element2D
    type(IntegratorPtrTYPE) :: integrator
    type(PointPtrTYPE), dimension(:), allocatable :: point
    call debugLog('Calculating shear stress')
    print'(A)', 'Calculating shear stress'
    call this%valueGPoints1D(problem)
    nLine = problem%domain%nLine
    if(nLine > 0) then
       nGauss = problem%domain%elementList1D%getGaussOrder()
    else
       nGauss = problem%domain%elementList2D%getGaussOrder()
    end if
    nGaussPointLine = getnGaussPointLine(nGauss, nLine)
    allocate(this%lineElemID(nLine))
    allocate(this%lineShS(nGaussPointLine))
    !Shear Stress for one dimensional elements:
    lineCount = 0
    nElem = problem%domain%getnLine()
    do iElem = 1, nElem
       element1D = problem%domain%elementList1D%getElement(iElem)
       this%lineElemID(iElem) = element1D%getID()
       nPoint = element1D%getnPoint()
       integrator = element1D%getIntegrator()
       do iGauss = 1, integrator%ptr%integTerms
          xi = integrator%ptr%gPoint(iGauss,1)
          jacobian1D = element1D%jacobian(xi)
          ShS = 0.d0
          do i = 1, nPoint
             ShS = ShS + integrator%ptr%dShapeFunc(iGauss,1,i) &
                  *problem%dof(element1D%getPointID(i))/jacobian1D
          end do
          d33 = element2D%ptr%material%ptr%d33
          call addLineShearStress(this, ShS)
       end do
    end do
    !Shear Stress for bidimensional elements:
    nTriang = problem%domain%nTriang
    nQuad = problem%domain%nQuad
    if(nTriang == 0 .and. nQuad == 0) return
    call this%valueGPointsTriang(problem)
    call this%valueGPointsQuad(problem)
    nGaussPointTriang = getnGaussPointTriang(nGauss, nTriang)
    nGaussPointQuad = getnGaussPointQuad(nGauss, nQuad)
    allocate(this%triangElemID(nTriang))
    allocate(this%triangShS(nGaussPointTriang))
    allocate(this%quadElemID(nQuad))
    allocate(this%quadShS(nGaussPointQuad))
    triangCount = 0
    quadCount = 0
    triangElemCount = 0
    quadElemCount = 0
    nElem = problem%domain%getnTriang() + problem%domain%getnQuad()
    do iElem = 1, nElem
       element2D = problem%domain%elementList2D%getElement(iElem)
       nPoint = element2D%getnPoint()
       nDof = element2D%getnDof()
       integrator = element2D%getIntegrator()
       if(nPoint == 3 .or. nPoint == 6) then
          triangElemCount = triangElemCount + 1
          addShearStress => addTriangShearStress
          this%triangElemID(triangElemCount) = element2D%getID()
       else if(nPoint == 4 .or. nPoint == 8) then
          quadElemCount = quadElemCount + 1
          addShearStress => addQuadShearStress
          this%quadElemID(quadElemCount) = element2D%getID()
       else
          print'(A)', '** Shear Stress ERROR1 **'
       end if
       do iGauss = 1, integrator%ptr%integTerms
          xi = integrator%ptr%gPoint(iGauss,1)
          eta = integrator%ptr%gPoint(iGauss,2)
          jacobian = element2D%jacobian(xi, eta)
          jacobianDet = element2D%jacobianDet(jacobian)
          allocate(dsf(2,nPoint*nDof))
          dsf = integrator%ptr%dShapeFunc(iGauss,:,:)
          ShS = 0.d0
          do i = 1, nPoint
             bi = jacobian(2,2)*dsf(1,nDof*i-1) - jacobian(1,2)*dsf(2,nDof*i-1)
             ci = jacobian(1,1)*dsf(2,nDof*i-1) - jacobian(2,1)*dsf(1,nDof*i-1)
             dNidx = bi/jacobianDet
             dNidy = ci/jacobianDet
             ShS = ShS + dNidx*problem%dof(element2D%getPointID(i)*nDof) &
             + dNidy*problem%dof(element2D%getPointID(i)*nDof-1)
          end do
          d33 = element2D%ptr%material%ptr%d33
          call addShearStress(this, d33*ShS)
          deallocate(dsf)
       end do
    end do
  end subroutine calculateShearStress

  subroutine valueGPoints1D(this, problem)
    implicit none
    class(ShearStressTYPE), intent(inout) :: this
    class(StructProblemTYPE), intent(inout) :: problem
    type(IntegratorPtrTYPE) :: integrator
    integrator = problem%domain%elementList1D%getIntegrator()
    if(allocated(integrator%ptr%gPoint)) then
       allocate(this%lineGPoint(size(integrator%ptr%gPoint,1)))
       this%lineGPoint(:) = integrator%ptr%gPoint(:,1)
    end if
  end subroutine valueGPoints1D

  subroutine valueGPointsTriang(this, problem)
    implicit none
    class(ShearStressTYPE), intent(inout) :: this
    class(StructProblemTYPE), intent(inout) :: problem
    type(IntegratorPtrTYPE) :: integrator
    integrator = problem%domain%elementList2D%getTriangIntegrator()
    if(allocated(integrator%ptr%gPoint)) then
       allocate(this%triangGPoint(size(integrator%ptr%gPoint,1),size(integrator%ptr%gPoint,2)))
       this%triangGPoint = integrator%ptr%gPoint
    end if
  end subroutine valueGPointsTriang

  subroutine valueGPointsQuad(this, problem)
    implicit none
    class(ShearStressTYPE), intent(inout) :: this
    class(structProblemTYPE), intent(inout) :: problem
    type(IntegratorPtrTYPE) :: integrator
    integrator = problem%domain%elementList2D%getQuadIntegrator()
    if(allocated(integrator%ptr%gPoint)) then
       allocate(this%quadGPoint(size(integrator%ptr%gPoint,1),size(integrator%ptr%gPoint,2)))
       this%quadGPoint = integrator%ptr%gPoint
    end if
  end subroutine valueGPointsQuad

  subroutine addLineShearStress(this, ShS)
    implicit none
    class(ShearStressTYPE), intent(inout) :: this
    real(rkind), intent(in) :: ShS
    lineCount = lineCount + 1
    this%lineShS(lineCount) = ShS
  end subroutine addLineShearStress

  subroutine addTriangShearStress(this, ShS)
    implicit none
    class(ShearStressTYPE), intent(inout) :: this
    real(rkind), intent(in) :: ShS
    triangCount = triangCount + 1
    this%triangShS(triangCount) = ShS
  end subroutine addTriangShearStress

  subroutine addQuadShearStress(this, ShS)
    implicit none
    class(ShearStressTYPE), intent(inout) :: this
    real(rkind), intent(in) :: ShS
    quadCount = quadCount + 1
    this%quadShS(quadCount) = ShS
  end subroutine addQuadShearStress

  integer(ikind) function getnGaussPointLine(nGauss, nLine)
    implicit none
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: nLine
    getnGaussPointLine = nLine*nGauss
  end function getnGaussPointLine

  integer(ikind) function getnGaussPointTriang(nGauss, nTriang)
    implicit none
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: nTriang
    if(nGauss == 1) then
       getnGaussPointTriang = nTriang
    else if(nGauss == 2) then
       getnGaussPointTriang = 3*nTriang
    else if(nGauss == 3) then
       getnGaussPointTriang = 4*nTriang
    else
       print*, '** Input Gauss Order not supported! **'
    end if
  end function getnGaussPointTriang
  
  integer(ikind) function getnGaussPointQuad(nGauss, nQuad)
    implicit none
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: nQuad
    if(nGauss == 1) then
       getnGaussPointQuad = nQuad
    else if(nGauss == 2) then
       getnGaussPointQuad = 4*nQuad
    else if(nGauss == 3) then
       getnGaussPointQuad = 9*nQuad
    else
       print*, '** Input Gauss Order not supported! **'
    end if
  end function getnGaussPointQuad

end module ShearStressMOD
