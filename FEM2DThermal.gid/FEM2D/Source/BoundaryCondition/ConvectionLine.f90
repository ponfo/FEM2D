module ConvectionLineMOD
  use tools
  use SparseKit
  use PointPtrMOD
  use IntegratorMOD
  use IntegratorPtrMOD
  use ThermalElementList2DMOD
  use Element2DPtrMOD
  implicit none
  private
  public :: ConvectionLineTYPE, convectionLine
  type ConvectionLineTYPE
     integer(ikind)                            :: elemID
     integer(ikind), dimension(:), allocatable :: pointID
     real(rkind)                               :: coef
     real(rkind)                               :: temp
     type(IntegratorPtrTYPE)                   :: integrator1D
   contains
     procedure :: init
     procedure :: getElemID
     procedure :: getPointID
     procedure :: getCoef
     procedure :: getTemp
     procedure :: getnPoint
     procedure :: apply
     procedure, private :: setupIntegration
  end type ConvectionLineTYPE
  
  interface convectionLine
     procedure :: constructor
  end interface convectionLine
  
contains
  
  type(ConvectionLineTYPE) function constructor(elemID, pointID, coef, temp, integrator1D)
    implicit none
    integer(ikind), intent(in) :: elemID
    integer(ikind), dimension(:), intent(in) :: pointID
    real(rkind), intent(in) :: coef
    real(rkind), intent(in) :: temp
    type(IntegratorTYPE), target, intent(in) :: integrator1D
    call constructor%init(elemID, pointID, coef, temp, integrator1D)
  end function constructor
  
  subroutine init(this, elemID, pointID, coef, temp, integrator1D)
    implicit none
    class(ConvectionLineTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: elemID
    integer(ikind), dimension(:), intent(in) :: pointID
    real(rkind), intent(in) :: coef
    real(rkind), intent(in) :: temp
    type(IntegratorTYPE), target, intent(in) :: integrator1D
    this%elemID = elemID
    this%pointID = pointID
    this%coef = coef
    this%temp = temp
    this%integrator1D%ptr => integrator1D
  end subroutine init
  
  integer(ikind) function getElemID(this)
    implicit none
    class(ConvectionLineTYPE), intent(inout) :: this
    getElemID = this%elemID
  end function getElemID
  
  function getPointID(this)
    implicit none
    class(ConvectionLineTYPE), intent(inout) :: this
    integer(ikind), dimension(size(this%pointID)) :: getPointID
    getPointID = this%pointID
  end function getPointID
  
  real(rkind) function getCoef(this)
    implicit none
    class(ConvectionLineTYPE), intent(inout) :: this
    getCoef = this%coef
  end function getCoef
  
  real(rkind) function getTemp(this)
    implicit none
    class(ConvectionLineTYPE), intent(inout) :: this
    getTemp = this%temp
  end function getTemp
  
  integer(ikind) function getnPoint(this)
    implicit none
    class(ConvectionLineTYPE), intent(inout) :: this
    getnPoint = size(this%pointID)
  end function getnPoint

  subroutine apply(this, elementList, stiffness, rhs)
    implicit none
    class(ConvectionLineTYPE), intent(inout) :: this
    type(ThermalElementList2DTYPE), intent(inout) :: elementList
    class(Sparse), intent(inout) :: stiffness
    real(rkind), dimension(:), intent(inout) :: rhs
    type(Element2DPtrTYPE) :: element
    integer(ikind) :: i, j, k, id, idi, idj, index
    real(rkind) :: magicFactor, int
    real(rkind), dimension(this%integrator1D%ptr%integTerms) :: jacobianDet
    element = elementList%getElement(this%elemID)
    call this%setupIntegration(element, magicFactor, jacobianDet)
    do i = 1, this%getnPoint()
       do j = 1, this%getnPoint()
          int = 0
          do k = 1, this%integrator1D%ptr%integTerms
             int = int + this%integrator1D%ptr%weight(k)*this%integrator1D%ptr%shapeFunc(i,k) &
                  *this%integrator1D%ptr%shapeFunc(j,k)*magicFactor*this%coef*jacobianDet(k)
          end do
          idi = element%getPointID(this%pointID(i))
          idj = element%getPointID(this%pointID(j))
          index = stiffness%AI(idi)
          do while(index < stiffness%AI(idi+1))
             if(stiffness%AJ(index) == idj) then
                stiffness%A(index) = stiffness%A(index) + int
             end if
             index = index + 1
          end do
       end do
    end do
    do i = 1, this%getnPoint()
       int = 0
       do j = 1, this%integrator1D%ptr%integTerms
          int = int + this%integrator1D%ptr%weight(j)*this%integrator1D%ptr%shapeFunc(i,j) &
               *this%coef*this%temp*magicFactor*jacobianDet(j)
       end do
       id = element%getPointID(this%pointID(i))
       rhs(id) = rhs(id) + int
    end do
  end subroutine apply

  subroutine setupIntegration(this, element, magicFactor, jacobianDet)
    implicit none
    class(ConvectionLineTYPE), intent(inout) :: this
    type(Element2DPtrTYPE), intent(inout) :: element
    real(rkind), intent(out) :: magicFactor
    real(rkind), dimension(this%integrator1D%ptr%integTerms), intent(out) :: jacobianDet
    integer(ikind) :: i, j, nPoint
    real(rkind), dimension(this%integrator1D%ptr%integTerms,2,2) :: jacobian
    real(rkind), dimension(this%integrator1D%ptr%integTerms) :: u, v
    real(rkind), dimension(this%integrator1D%ptr%integTerms,2,element%ptr%nPoint) :: dsf
    type(PointPtrTYPE), dimension(size(this%pointID)) :: point
    do i = 1, this%getnPoint()
       point(i) = element%getPoint(this%pointID(i))
    end do
    jacobian = 0
    nPoint = element%getnPoint()
    if(nPoint == 3) then
       magicFactor = 1.d0/2.d0
       u = transformGPoint1D(this%integrator1D,0.d0,1.d0)
       v = 0.d0
       do j = 1, this%integrator1D%ptr%integTerms
          dsf(j,:,:) = element%dShapeFunc(u(j),v(j))
       end do
       do j = 1, this%integrator1D%ptr%integTerms
          jacobian(j,1,1) = dsf(j,1,1)*point(1)%getx() &
               + dsf(j,1,2)*point(2)%getx()
          jacobian(j,1,2) = dsf(j,1,1)*point(1)%gety() &
               + dsf(j,1,2)*point(2)%gety()
          jacobianDet(j) = sqrt(jacobian(j,1,1)**2+jacobian(j,1,2)**2)
       end do
    else if(nPoint == 6) then
       magicFactor = 1.d0/2.d0
       u = transformGPoint1D(this%integrator1D,0.d0,1.d0)
       v = 0.d0
       do j = 1, this%integrator1D%ptr%integTerms
          dsf(j,:,:) = element%dShapeFunc(u(j),v(j))
       end do
       do j = 1, this%integrator1D%ptr%integTerms
          jacobian(j,1,1) = dsf(j,1,1)*point(1)%getx() &
               + dsf(j,1,2)*point(2)%getx()            &
               + dsf(j,1,4)*point(3)%getx()
          jacobian(j,1,2) = dsf(j,1,1)*point(1)%gety() &
               + dsf(j,1,2)*point(2)%gety()            &
               + dsf(j,1,4)*point(3)%gety()
          jacobianDet(j) = sqrt(jacobian(j,1,1)**2+jacobian(j,1,2)**2)
       end do
    else if(nPoint == 4) then
       magicFactor = 1.d0
       u = this%integrator1D%ptr%gPoint(:,1)
       v = -1
       do j = 1, this%integrator1D%ptr%integTerms
          dsf(j,:,:) = element%dShapeFunc(u(j),v(j))
       end do
       do j = 1, this%integrator1D%ptr%integTerms
          jacobian(j,1,1) = dsf(j,1,1)*point(1)%getx() &
               + dsf(j,1,2)*point(2)%getx()
          jacobian(j,1,2) = dsf(j,1,1)*point(1)%gety() &
               + dsf(j,1,2)*point(2)%gety()
          jacobianDet(j) = sqrt(jacobian(j,1,1)**2+jacobian(j,1,2)**2)
       end do
    else if(nPoint == 8) then
       magicFactor = 1.d0
       u = this%integrator1D%ptr%gPoint(:,1)
       v = -1
       do j = 1, this%integrator1D%ptr%integTerms
          dsf(j,:,:) = element%dShapeFunc(u(j),v(j))
       end do
       do j = 1, this%integrator1D%ptr%integTerms
          jacobian(j,1,1) = dsf(j,1,1)*point(1)%getx() &
               + dsf(j,1,2)*point(2)%getx()            &
               + dsf(j,1,5)*point(3)%getx()
          jacobian(j,1,2) = dsf(j,1,1)*point(1)%gety() &
               + dsf(j,1,2)*point(2)%gety()            &
               + dsf(j,1,5)*point(3)%gety()
          jacobianDet(j) = sqrt(jacobian(j,1,1)**2+jacobian(j,1,2)**2)
       end do
    else
       print'(A)', '** convectionLineInt ERROR1 **'
    end if
  end subroutine setupIntegration

  function transformGPoint1D(integrator1D, lower, upper) result(tGPoint)
    implicit none
    class(IntegratorPtrTYPE), intent(inout) :: integrator1D
    real(rkind), intent(in) :: lower
    real(rkind), intent(in) :: upper
    real(rkind), dimension(size(integrator1D%ptr%gPoint,1)) :: tGPoint
    integer(ikind) :: i
    do i = 1, size(integrator1D%ptr%shapeFunc,2)
       tGPoint(i) = (upper+lower)/2.d0 + ((upper-lower)/2.d0)*integrator1D%ptr%gPoint(i,1)
    end do
  end function transformGPoint1D
  
end module ConvectionLineMOD
