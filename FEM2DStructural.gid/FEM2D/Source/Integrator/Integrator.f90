module IntegratorMOD
  use tools
  use DebuggerMOD
  implicit none
  private
  public :: IntegratorTYPE, integrator
  type :: IntegratorTYPE
     integer(ikind)                             :: gaussOrder
     integer(ikind)                             :: integTerms
     real(rkind), dimension(:)    , allocatable :: weight
     real(rkind), dimension(:,:)  , allocatable :: gPoint
     real(rkind), dimension(:,:)  , allocatable :: shapeFunc
     real(rkind), dimension(:,:,:), allocatable :: dShapeFunc
   contains
     procedure, public :: init
     procedure, public :: valueGPoints
     procedure, private :: getG1D
     procedure, private :: getGTriangle
     procedure, private :: getGSquare
  end type IntegratorTYPE

  interface integrator
     procedure constructor
  end interface integrator
  
contains

  type(IntegratorTYPE) function constructor(gaussOrder, type)
    implicit none
    integer(ikind), intent(in) :: gaussOrder
    character(*), intent(in) :: type
    call constructor%init(gaussOrder, type)
  end function constructor

  subroutine init(this, gaussOrder, type)
    implicit none
    class(IntegratorTYPE), intent(inout) :: this
    integer(ikind), intent(in) :: gaussOrder
    character(*), intent(in) :: type
    this%gaussOrder = gaussOrder
    call this%valueGPoints(type)
  end subroutine init
  
  subroutine valueGPoints(this, type)
    implicit none
    class(IntegratorTYPE), intent(inout) :: this
    character(*), intent(in) :: type
    if(trim(type) == 'Line' .or. &
       trim(type) == 'line' .or. &
       trim(type) == 'LINE') then
       call this%getG1D()
    else if(trim(type) == 'Triangle' .or. &
            trim(type) == 'triangle' .or. &
            trim(type) == 'triang'   .or. & 
            trim(type) == 'TRIANGLE') then
       call this%getGTriangle()
    else if(trim(type) == 'Quadrilateral' .or. &
            trim(type) == 'quadrilateral' .or. &
            trim(type) == 'QUADRILATERAL' .or. &
            trim(type) == 'quad')          then
       call this%getGSquare()
    end if
  end subroutine valueGPoints

  subroutine getG1D(this)
    implicit none
    class(IntegratorTYPE), intent(inout) :: this
    real(rkind), parameter :: pi = dacos(-1.d0)
    real(rkind) :: f, df, dx, r
    integer(ikind) :: i, j, iter, k, counter
    real(rkind), dimension(:), allocatable :: p0, p1, tmp
    real(rkind), dimension(:,:), allocatable :: G
    this%integTerms = this%gaussOrder
    allocate(G(2,this%gaussOrder))
    p0 = [1.d0]
    p1 = [1.d0, 0.d0]
    do k = 2, this%gaussOrder
       tmp = ((2*k-1)*[p1,0.d0]-(k-1)*[0.d0, 0.d0,p0])/k
       p0 = p1; p1 = tmp
    end do
    do i = 1, this%gaussOrder
       r = cos(pi*(i-0.25)/(this%gaussOrder+0.5))
       do iter = 1, 10
          f = p1(1); df = 0.
          do k = 2, size(p1)
             df = f + r*df
             f  = p1(k) + r * f
          end do
          dx =  f / df
          r = r - dx
          if (abs(dx)<10*epsilon(dx)) exit
       end do
       G(1,i) = r
       G(2,i) = 2/((1-r**2)*df**2)
    end do
    allocate(this%weight(this%integTerms))
    allocate(this%gPoint(this%integTerms,1))
    this%weight(:) = G(2,:)
    this%gPoint(:,1) = G(1,:)
    call debugLog('        Allocated weights: ', size(this%weight))
    call debugLog('        Allocated gPoints: ', size(this%gPoint))
    deallocate(G)
  end subroutine getG1D

    subroutine getGTriangle(this)
    implicit none
    class(IntegratorTYPE), intent(inout) :: this
    if(this%gaussOrder == 1) then
       this%integTerms  = 1
       allocate(this%weight(this%integTerms))
       allocate(this%gPoint(this%integTerms,2))
       this%weight(1)   = 1.d0/2.d0
       this%gPoint(1,1) = 1.d0/3.d0
       this%gPoint(1,2) = 1.d0/3.d0
    else if(this%gaussOrder == 2) then
       this%integTerms  = 3
       allocate(this%weight(this%integTerms))
       allocate(this%gPoint(this%integTerms,2))
       this%weight(1)   = 1.d0/6.d0
       this%weight(2)   = 1.d0/6.d0
       this%weight(3)   = 1.d0/6.d0
       this%gPoint(1,1) = 1.d0/2.d0
       this%gPoint(1,2) = 0.d0
       this%gPoint(2,1) = 1.d0/2.d0
       this%gPoint(2,2) = 1.d0/2.d0
       this%gPoint(3,1) = 0.d0
       this%gPoint(3,2) = 1.d0/2.d0
    else if(this%gaussOrder == 3) then
       this%integTerms  = 4
       allocate(this%weight(this%integTerms))
       allocate(this%gPoint(this%integTerms,2))
       this%weight(1)   = -27.d0/96.d0
       this%weight(2)   = 25.d0/96.d0
       this%weight(3)   = 25.d0/96.d0
       this%weight(4)   = 25.d0/96.d0
       this%gPoint(1,1) = 1.d0/3.d0
       this%gPoint(1,2) = 1.d0/3.d0
       this%gPoint(2,1) = 1.d0/5.d0
       this%gPoint(2,2) = 1.d0/5.d0
       this%gPoint(3,1) = 3.d0/5.d0
       this%gPoint(3,2) = 1.d0/5.d0
       this%gPoint(4,1) = 1.d0/5.d0
       this%gPoint(4,2) = 3.d0/5.d0
    else
       print'(A)', '** Input Gauss Order not supported for triangular elements! **'
    end if
    call debugLog('        Allocated weights: ', size(this%weight))
    call debugLog('        Allocated gPoints: ', size(this%gPoint))
  end subroutine getGTriangle

  subroutine getGSquare(this)
    implicit none
    class(IntegratorTYPE), intent(inout) :: this
    real(rkind), parameter :: pi = dacos(-1.d0)
    real(rkind) :: f, df, dx, r
    integer(ikind) :: i, j, iter, k, counter
    real(rkind), dimension(:), allocatable :: p0, p1, tmp
    real(rkind), dimension(:,:), allocatable :: G
    this%integTerms = this%gaussOrder**2
    allocate(G(2,this%gaussOrder))
    p0 = [1.d0]
    p1 = [1.d0, 0.d0]
    do k = 2, this%gaussOrder
       tmp = ((2*k-1)*[p1,0.d0]-(k-1)*[0.d0, 0.d0,p0])/k
       p0 = p1; p1 = tmp
    end do
    do i = 1, this%gaussOrder
       r = cos(pi*(i-0.25)/(this%gaussOrder+0.5))
       do iter = 1, 10
          f = p1(1); df = 0.
          do k = 2, size(p1)
             df = f + r*df
             f  = p1(k) + r * f
          end do
          dx =  f / df
          r = r - dx
          if (abs(dx)<10*epsilon(dx)) exit
       end do
       G(1,i) = r
       G(2,i) = 2/((1-r**2)*df**2)
    end do
    allocate(this%weight(this%integTerms))
    allocate(this%gPoint(this%integTerms,2))
    counter = 0
    do i = 1, this%gaussOrder
       do j = 1, this%gaussOrder
          counter = counter + 1
          this%weight(counter) = G(2,i)*G(2,j)
          this%gPoint(counter,1) = G(1,i)
          this%gPoint(counter,2) = G(1,j)
       end do
    end do
    call debugLog('        Allocated weights: ', size(this%weight))
    call debugLog('        Allocated gPoints: ', size(this%gPoint))
    deallocate(G)
  end subroutine getGSquare
end module IntegratorMOD
