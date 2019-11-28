!*************************************************************
!           Instituto Universitario Aeronautico
!                Dpto. Mecanica Aeronautica
!*************************************************************
! Filename      : SparseKit.f90
! Version       : 0.9
! Date          : 14-08-2019
! Programmer(s) : F. Airaudo(fairaudo574@alumnos.iua.edu.ar)
!                 M. Zuñiga(mzuniga433@alumnos.iua.edu.ar)
!*************************************************************
! Description (brief):
!                     Module for all algebraic and functional
!                     Operations with Sparse Matrices. That is,
!                     For matrices where the number of nonzeros
!                     highly exceeds the number of zeros, this
!                     tools allow you to declare a derived
!                     data type to work with them. 
!************************************************************* 
! Dependecies:
!             Use tools        - Filename: Utilities.f90
!             Use QuickSortMod - Filename: quicksort.f90
!*************************************************************
! Public procedures:
!              Type(Sparse)
!              Derived type Procedures 
!                   Subroutine init
!                   Subroutine update
!                   Subroutine append
!                   Subroutine makeCRS
!                   Function get
!                   Function getn
!                   Function getNonZeros
!                   Subroutine printValue
!                   Subroutine printNonZeros
!                   Subroutine printAll
!                   Subroutine deleteRowAndCol
!                   Subroutine free
!              Module Procedures:
!                   Function transpose
!                   Function inverse
!                   Function norm
!                   Function gmres
!                   Function jacobiEigen
!                   Function trace
!                   Function inverseGMRESD
!                   Function id
!                   Subroutine sparse_sparse_prod ->  Operator:
!                   Subroutine coef_sparse_prod   ->    
!                   Subroutine sparse_vect_prod   ->     (*)
!                   Subroutine sparse_sparse_add  ->  Operator:
!                                                        (+)
!                   Subroutine sparse_sparse_sub  ->  Operator:
!                                                        (-)
!*************************************************************
module SparseKit
  !***********************************************
  !*                 EXTERNS                     *
  !***********************************************
  use tools
  use quickSortMod
  implicit none
  private
  public :: Sparse, operator(*), operator(+), operator(-), transpose&
       , norm, gmres, inverse, jacobiEigen, trace, inverseGMRESD, id&
       , lcholesky, det, biCGOMP, bicGrad, CGOMP
  type Triplet
     real(rkind), dimension(:), allocatable :: A
     integer(ikind), dimension(:), allocatable :: row
     integer(ikind), dimension(:), allocatable :: col
  end type Triplet
  type Sparse
     !private
     real(rkind), dimension(:), allocatable :: A
     integer(ikind), dimension(:), allocatable :: AI
     integer(ikind), dimension(:), allocatable :: AJ
     integer(ikind), dimension(:), allocatable :: rowCounter
     integer(ikind) :: n
     integer(ikind) :: nnz
     integer(ikind) :: counter
     type(triplet) :: triplet
   contains
     procedure, public :: init
     procedure, public :: update
     procedure, public :: append
     procedure, public :: makeCRS

     procedure, public :: get
     procedure, public :: getnnz
     procedure, public :: getn

     procedure, public :: printValue
     procedure, public :: printNonZeros
     procedure, public :: printAll

     procedure, public :: deleteRowAndCol

     procedure, public :: free

     procedure, private :: handleDuplicates
  end type Sparse

  interface operator(*)
     module procedure sparse_sparse_prod
     module procedure sparse_vect_prod
     module procedure coef_sparse_prod
  end interface operator(*)
  interface operator(+)
     module procedure sparse_sparse_add
  end interface operator(+)
  interface operator(-)
     module procedure sparse_sparse_sub
  end interface operator(-)
  interface sparse
     procedure constructor
  end interface sparse
  interface transpose
     module procedure transpose
  end interface transpose
  interface norm
     module procedure norm
  end interface norm
  interface trace
     module procedure trace
  end interface trace
  interface det
     module procedure det
  end interface det
  interface lcholesky
     module procedure lcholesky
  end interface lcholesky
  interface gmres
     module procedure gmres
  end interface gmres
  interface inverse
     module procedure inverse
  end interface inverse
  interface jacobiEigen
     module procedure jacobiEigen
  end interface jacobiEigen
  interface biCGOMP
     module procedure biCGOMP
  end interface biCGOMP
  interface bicGrad
     module procedure bicGrad
  end interface bicGrad
  interface CGOMP
     module procedure CGOMP
  end interface CGOMP
  !***********************************************
  !*          LOCAL PRIVATE VARIABLES            *
  !***********************************************
  logical :: isCRSDone = .false.
  real(rkind), dimension(:), allocatable :: valueVector
  real(rkind), dimension(:), allocatable :: auxA
  integer(ikind), dimension(:), allocatable :: auxAJ
  integer(ikind), dimension(:), allocatable :: rowVector
  integer(ikind) :: repeats, l

contains
  !***************************************************
  ! Constructor:
  !     Initializes the object with and estimate of
  !     the matrix's non zeros and the number of rows
  !     Won't work if given an inferior nnz
  !  
  ! Parameters:
  !     Input, nnz(int), rows(int)
  !     Output, -
  !***************************************************
  type(Sparse) function constructor(nnz, rows)
    implicit none
    integer(ikind), intent(in) :: nnz
    integer(ikind), intent(in) :: rows
    call constructor%init(nnz, rows)
  end function constructor
  subroutine init(this, nnz, rows)
    implicit none
    class(Sparse), intent(inout) :: this
    integer(ikind), intent(in) :: nnz
    integer(ikind), intent(in) :: rows
    this%nnz = nnz
    this%n = rows
    allocate(this%triplet%A(this%nnz))
    allocate(this%triplet%row(this%nnz))
    allocate(this%triplet%col(this%nnz))
    this%triplet%A = 0
    this%triplet%row = 0
    this%triplet%col = 0
    this%counter = 0
  end subroutine init
  !***************************************************
  ! update:
  !     updates basic values
  !  
  ! Parameters:
  !     Input, nnz(int), rows(int)
  !     Output, -
  !***************************************************
  subroutine update(this, nnz, rows)
    implicit none
    class(sparse), intent(inout) :: this
    integer(ikind), intent(in) :: nnz
    integer(ikind), intent(in) :: rows
    isCRSDone = .false.
    if(allocated(this%triplet%A)) then
       deallocate(this%triplet%A, this%triplet%row, this%triplet%col)
    end if
    if(allocated(this%A)) then
       deallocate(this%A, this%AI, this%AJ)
    end if
    this%nnz = nnz
    this%n = rows
    allocate(this%triplet%A(this%nnz))
    allocate(this%triplet%row(this%nnz))
    allocate(this%triplet%col(this%nnz))
    this%triplet%A = 0
    this%triplet%row = 0
    this%triplet%col = 0
    this%counter = 0
  end subroutine update
  !***************************************************
  ! append:
  !     takes values one by one and appends it to a
  !     triplet format
  !  
  ! Parameters:
  !     Input, value(realrkind), row(int), col(int)
  !     Output, -
  !***************************************************
  subroutine append(this, val, row, col)
    implicit none
    class(Sparse), intent(inout) :: this
    real(rkind), intent(In) :: val
    integer(ikind), intent(In) :: row, col
    this%counter = this%counter + 1
    this%triplet%A(this%counter) = val
    this%triplet%row(this%counter) = row
    this%triplet%col(this%counter) = col
  end subroutine append
  !***************************************************
  ! makeCRS:
  !     Once all values have been appended to the
  !     triplet call this routine to make the CRS
  !  
  ! Parameters:
  !     Input, -
  !     Output, CRS is usable
  !***************************************************
  subroutine makeCRS(this, sortRows)
    implicit none
    class(Sparse), intent(inout) :: this
    logical, intent(in), optional :: sortRows
    integer(ikind) :: i
    isCRSDone = .true.
    allocate(this%rowCounter(this%n))
    !This%Counter entries in each row, including duplicates
    this%rowCounter = 0
    do i = 1, this%counter
       this%rowCounter(this%triplet%row(i)) = this%rowCounter(this%triplet%row(i)) + 1
    end do
    !Allocate auxA and auxAJ with nnz
    allocate(auxA(this%counter), auxAJ(this%counter))
    !Order A and AJ
    if(.not.present(sortRows) .or. sortRows) then
       call quicksort(this%triplet%row, this%triplet%col, this%triplet%A, 1, this%counter)
    end if
    !sum up duplicates
    call this%handleDuplicates()
    !Allocate A and AJ with nnz
    allocate(this%A(this%counter), this%AJ(this%counter))
    do i = 1, this%counter
       this%A(i) = auxA(i)
       this%AJ(i) = auxAJ(i)
    end do
    !Counstruct row pointers
    allocate(this%AI(this%n+1))
    this%AI = 1
    do i = 2, this%n+1
       this%AI(i) = this%AI(i-1) + this%rowCounter(i-1)
    end do
    deallocate(this%triplet%A, this%triplet%row, this%triplet%col)
    deallocate(auxA, auxAJ)
    deallocate(this%rowCounter)
    this%nnz = size(this%AJ)
  end subroutine makeCRS
  !***************************************************
  ! handleDuplicates:
  !     Sums up every duplicate on the triplet, also
  !     orders values in each row
  !  
  ! Parameters:
  !     Input, -
  !     Output, -
  !***************************************************
  subroutine handleDuplicates(this)
    implicit none
    class(Sparse), intent(inout) :: this
    logical :: mask
    integer(ikind) :: i, j, k
!!$    allocate(rowVector(maxval(this%rowCounter)))
!!$    allocate(valueVector(maxval(this%rowCounter)))
    this%counter = 1
    repeats = 0
    do i = 1, this%n
       rowVector = this%triplet%col(this%counter+repeats:this%counter+repeats+this%rowCounter(i)-1)
       valueVector = this%triplet%A(this%counter+repeats:this%counter+repeats+this%rowCounter(i)-1)
       call quicksort(rowVector, valueVector, 1, this%rowCounter(i))
       j = 0
       do while(j < this%rowCounter(i))
          j = j + 1
          auxA(this%counter) = valueVector(j)
          auxAJ(this%counter) = rowVector(j)
          k = j
          do while(k < this%rowCounter(i))
             k = k + 1
             mask = rowVector(j).eq.rowVector(k)
             if(mask) then
                !add values from k to j
                auxA(this%counter) = auxA(this%counter) + valueVector(k)
                !move k to the back
                call swap(rowVector(k), rowVector(this%rowCounter(i)))
                call swap(valueVector(k), valueVector(this%rowCounter(i)))
                this%rowCounter(i) = this%rowCounter(i) - 1
                repeats = repeats + 1
                k = k - 1
             end if
          end do
          this%counter = this%counter + 1
       end do
    end do
    this%counter = this%counter - 1
!!$    deallocate(rowVector)
!!$    deallocate(valueVector)
  end subroutine HandleDuplicates
  !***************************************************
  ! get:
  !     Gives sparse matrix's value from row i
  !     and col j
  !  
  ! Parameters:
  !     Input, i(int), j(int)
  !     Output, get(i,j)(realrkind)
  !***************************************************
  real(rkind) function get(this, i, j)
    implicit none
    class(Sparse), intent(inout) :: this
    integer(ikind), intent(in) :: i
    integer(ikind), intent(in) :: j
    integer(ikind) :: k
    k = this%AI(i)
    do while(k < this%AI(i+1))
       if(this%AJ(k) == j) then
          get = this%A(k)
          return
       end if
       k = k + 1
    end do
    get = 0.d0
  end function get
  !***************************************************
  ! getnnz:
  !     given ammount of non zeros
  !  
  ! Parameters:
  !     Input, -
  !     Output, getnnz()(integer(ikind))
  !***************************************************
  integer(ikind) function getnnz(this)
    implicit none
    class(Sparse), intent(inout) :: this
    getnnz = this%nnz
  end function getnnz
  !***************************************************
  ! getn:
  !     get order of matrix
  !  
  ! Parameters:
  !     Input, -
  !     Output, getn()(integer(ikind))
  !***************************************************
  integer(ikind) function getn(this)
    implicit none
    class(Sparse), intent(inout) :: this
    getn = this%n
  end function getn
  !***************************************************
  ! printValue:
  !     prints a single value either on console or a
  !     given filename
  !  
  ! Parameters:
  !     Input, i(int), j(int), filename(char)[opt]
  !     Output, value printed
  !***************************************************
  subroutine printValue(this, i, j, filename)
    implicit none
    integer(ikind), parameter :: fileunit = 90
    class(Sparse), intent(inout) :: this
    integer(ikind), intent(in) :: i
    integer(ikind), intent(in) :: j
    character(*), intent(in), optional :: filename
    if(present(filename)) then
       open(fileunit, file = trim(filename), access = 'append')
       write(fileunit,'(A,I0,A,I0,A,E14.7)') 'Matriz value at row ', i, ' column ', j, ' is ', this%get(i,j)
       close(fileunit)
    else
       write(*,'(A,I0,A,I0,A,E14.7)') 'Matriz value at row ', i, ' column ', j, ' is ', this%get(i,j)
    end if
  end subroutine printValue
  !***************************************************
  ! printNonZeros:
  !     Prints all non zeros in a list like format,
  !     either on console or on a given filename
  !  
  ! Parameters:
  !     Input, filename(char)[opt]
  !     Output, non zeros printed
  !***************************************************
  subroutine printNonZeros(this, filename)
    implicit none
    integer(ikind), parameter :: fileunit = 91
    class(Sparse), intent(inout) :: this
    character(*), intent(in), optional :: filename
    integer(ikind) :: i, j
    if(present(filename)) then
       open(fileunit, file = trim(filename), access = 'append')
       write(fileunit,'(A,I0,A)') 'Printing ', this%nnz, ' non zeros'
       write(fileunit,'(A,I0,A,I0)') 'Matrix dimension: ', this%n, 'x', this%n
       write(fileunit,'(A)') '------------------------------------'
       write(fileunit,'(A8,4X,2A12)') 'value', 'row', 'column'
       write(fileunit,'(A)') '------------------------------------'
       do i = 1, this%n
          do j = this%AI(i), this%AI(i+1)-1
             write(fileunit,'(E12.5,I12,I12)')  this%A(j), i, this%AJ(j)
          end do
       end do
       close(fileunit)
       return
    end if
    write(*,'(A,I0,A)') 'Printing ', this%nnz, ' non zeros'
    write(*,'(A,I0,A,I0)') 'Matrix dimension: ', this%n, 'x', this%n
    write(*,'(A)') '------------------------------------'
    write(*,'(A8,4X,2A12)') 'value', 'row', 'column'
    write(*,'(A)') '------------------------------------'
    do i = 1, size(this%AI)-1
       do j = this%AI(i), this%AI(i+1)-1
          write(*,'(E12.5,I12,I12)')  this%A(j), i, this%AJ(j)
       end do
    end do
  end subroutine printNonZeros
  !***************************************************
  ! printAll:
  !     Prints the whole matrix, zeros included, on
  !     console or on a filename if given.
  !  
  ! Parameters:
  !     Input, filename(char)[opt]
  !     Output, whole matrix printed
  !***************************************************
  subroutine printAll(this, filename)
    implicit none
    integer(ikind), parameter :: fileunit = 92
    class(Sparse), intent(inout) :: this
    character(*), intent(in), optional :: filename
    integer(ikind) :: i, j
    if(present(filename)) then
       open(fileunit, file = trim(filename), access = 'append')
       write(fileunit, '(/,A,I0,A,I0)') 'Printing Sparse Matrix, size: ', size(this%AI)-1, 'x', size(this%AI)-1 
       do i = 1, this%n
          write(fileunit,'(*(E14.7,2X))') (this%get(i,j),j=1,this%n)
       end do
       close(fileunit)
       return
    end if
    write(*, '(A,I0,A,I0)') 'Printing Sparse Matrix, size: ', this%n, 'x', this%n
    do i = 1, this%n
       write(*,'(*(E14.7,2X))') (this%get(i,j),j=1,this%n)
    end do
  end subroutine printAll
  !***************************************************
  ! deleteRowAndCol:
  !     Deletes a given row and column.
  !  
  ! Parameters:
  !     Input, row(int), col(int)
  !     Output, -
  !***************************************************
  subroutine deleteRowAndCol(this, row, col)
    Implicit none
    class(Sparse), intent(inout) :: this
    integer(ikind), intent(in) :: row
    integer(ikind), intent(in) :: col
    integer(ikind) :: i, j, k
    integer(ikind) :: rowSize
    integer(ikind), dimension(:), allocatable :: AI
    if(isCRSDone) then
       allocate(AI(size(this%AI)))
       do i = 1, size(AI)
          AI(i) = this%AI(i)
       end do
       deallocate(this%AI)
       do i = size(AI)-1, 1, -1
          if(i == row) then
             rowSize = AI(i+1)-AI(i)
             k = AI(i)
             do while(k < this%nnz-rowSize+1)
                this%AJ(k) = this%AJ(k+rowSize)
                this%A(k) = this%A(k+rowSize)
                k = k + 1
             end do
             this%nnz = this%nnz - rowSize
             do k = row, size(AI)-1
                AI(k) = AI(k+1) - rowSize
             end do
          else 
             do j = AI(i), AI(i+1)-1
                if(this%AJ(j) == col) then
                   k = j
                   do while(k < this%nnz)
                      this%A(k) = this%A(k+1)
                      this%AJ(k) = this%AJ(k+1)
                      k = k + 1
                   end do
                   do k = i+1, size(AI)
                      AI(k) = AI(k)-1
                   end do
                   this%nnz = this%nnz - 1
                end if
             end do
          end if
       end do
       do i = 1, this%nnz
          if(this%AJ(i) > col) then
             this%AJ(i) = this%AJ(i)-1
          end if
       end do
       allocate(this%AI(size(AI)-1))
       do i = 1, size(this%AI)
          this%AI(i) = AI(i)
       end do
    else
       i = 1
       do while(i < this%counter)
          if(this%triplet%row(i) == row .or. this%triplet%col(i) == col) then
             this%triplet%row(i) = this%triplet%row(this%counter)
             this%triplet%col(i) = this%triplet%col(this%counter)
             this%triplet%A(i) = this%triplet%A(this%counter)
             this%counter = this%counter - 1
          else
             i = i + 1
          end if
       end do
       if(this%triplet%row(i) == row .or. this%triplet%col(i) == col) then
          this%counter = this%counter - 1
       end if
       do i = 1, this%counter
          if(this%triplet%row(i) > row) then
             this%triplet%row(i) = this%triplet%row(i) - 1
          end if
          if(this%triplet%col(i) > col) then
             this%triplet%col(i) = this%triplet%col(i) - 1
          end if
       end do
    end if
    this%n = this%n - 1
  end subroutine deleteRowAndCol
  !***************************************************
  ! free:
  !     Clears memory space taken by the sparse matrix
  !  
  ! Parameters:
  !     Input, -
  !     Output, -
  !***************************************************
  subroutine free(this)
    implicit none
    class(Sparse), intent(inout) :: this
    if(allocated(this%A)) deallocate(this%A)
    if(allocated(this%AJ)) deallocate(this%AJ)
    if(allocated(this%AI)) deallocate(this%AI)
  end subroutine free



  !***********************************************
  !*              MODULE PROCEDURES              *
  !***********************************************

  !***************************************************
  ! sparse_sparse_prod:
  !     performs the product between two given sparse
  !     matrices. Operator: (*).
  !  
  ! Parameters:
  !     Input, a(Sparse), b(Sparse)
  !     Output, c(Sparse)
  !***************************************************
  function sparse_sparse_prod(a, b) result(c)
    implicit none
    class(Sparse), intent(in) :: a
    class(Sparse), intent(in) :: b
    type(Sparse) :: c
    type(Sparse) :: bTranspose
    real(rkind) :: Cij
    integer(ikind) :: i, j, k
    integer(ikind) :: nnz
    logical :: nnzFound
    if(a%n /= b%n) then
       print'(A)', '** diferent sizes in input sparse matrices! **'
       return
    end if
    bTranspose = transpose(b)
    !find ammount of nnz
    nnz = 0
    do i = 1, a%n
       do j = 1, bTranspose%n
          nnzFound = .false.
          do k = a%AI(i), a%AI(i+1)-1
             do l = bTranspose%AI(j), bTranspose%AI(j+1)-1
                if(a%AJ(k) == bTranspose%AJ(l)) nnzFound = .true.
             end do
          end do
          if(nnzFound) nnz = nnz + 1
       end do
    end do
    c = sparse(nnz = nnz, rows = a%n)
    do i = 1, a%n
       do j = 1, bTranspose%n
          Cij = 0
          do k = a%AI(i), a%AI(i+1)-1
             do l = bTranspose%AI(j), bTranspose%AI(j+1)-1
                if(a%AJ(k) == bTranspose%AJ(l)) then
                   Cij = Cij + a%A(k)*bTranspose%A(l)
                end if
             end do
          end do
          if(abs(Cij) > 1d-12) call c%append(Cij, i, j)
       end do
    end do
    call bTranspose%free()
    call c%makeCRS(.false.)
  end function sparse_sparse_prod
  !***************************************************
  ! sparse_vect_prod(*):
  !      Performs the product between a sparse matrix
  !      and a condensed real vector.
  !      Operator: (*).
  !  
  ! Parameters:
  !     Input, mat(Sparse), vect(realrkind)
  !     Output, res(realrkind)
  !***************************************************
  function sparse_vect_prod(mat, vect) result(res)
    implicit none
    class(Sparse), intent(in) :: mat
    real(rkind), dimension(:), intent(in) :: vect
    real(rkind), dimension(size(vect)) :: res
    real(rkind) :: val
    integer(ikind) :: i, k
    !$OMP PARALLEL DO PRIVATE(i, k, val)
    do i = 1, mat%n
       val = 0.d0
       do k = mat%AI(i), mat%AI(i+1)-1
          val = val + mat%A(k)*vect(mat%AJ(k))
       end do
       res(i) = val
    end do
    !$OMP END PARALLEL DO 
  end function sparse_vect_prod
  !***************************************************
  ! coef_sparse_prod(*):
  !      Performs the product between a sparse matrix
  !      and a coeficient real vector.
  !      Operator: (*).
  !  
  ! Parameters:
  !     Input, mat(Sparse), coef(realrkind)
  !     Output, res(realrkind)
  !***************************************************
  function coef_sparse_prod(coef ,mat) result(res)
    implicit none
    type(Sparse), intent(in) :: mat
    type(Sparse) :: res
    real(rkind), intent(in) :: coef
    real(rkind) :: c
    integer(ikind) :: i
    c = coef
    res = mat
    do i = 1, res%nnz
       res%A(i) = res%A(i)*c
    end do
  end function coef_sparse_prod
  !***************************************************
  ! sparse_sparse_add:
  !     performs the addition of sparse a plus
  !     sparse b. Operator: (+).
  !  
  ! Parameters:
  !     Input, a(Sparse), b(Sparse)
  !     Output, c(Sparse)
  !***************************************************
  function sparse_sparse_add(a, b) result(c)
    implicit none
    class(Sparse), intent(in) :: a
    class(Sparse), intent(in) :: b
    type(Sparse) :: c
    integer(ikind) :: counter, rowSize, i, k
    c = sparse(nnz = a%nnz+b%nnz, rows = a%n)
    counter = 1
    do i = 1, a%n
       rowSize = a%AI(i+1) - a%AI(i)
       do k = counter, counter+rowSize-1
          call c%append(a%A(k), i, a%AJ(k))
       end do
       counter = counter + rowSize
    end do
    counter = 1
    do i = 1, b%n
       rowSize = b%AI(i+1) - b%AI(i)
       do k = counter, counter+rowSize-1
          call c%append(b%A(k), i, b%AJ(k))
       end do
       counter = counter + rowSize
    end do
    call c%makeCRS
  end function sparse_sparse_add
  !***************************************************
  ! sparse_sparse_sub:
  !     performs the subtraction of sparse a and
  !     sparse b. Operator: (-).
  !  
  ! Parameters:
  !     Input, a(Sparse), b(Sparse)
  !     Output, c(Sparse)
  !***************************************************
  function sparse_sparse_sub(a, b) result(c)
    implicit none
    class(Sparse), intent(in) :: a
    class(Sparse), intent(in) :: b
    type(Sparse) :: c
    integer(ikind) :: counter, rowSize, i, k
    c = sparse(nnz = a%nnz+b%nnz, rows = a%n)
    counter = 1
    do i = 1, a%n
       rowSize = a%AI(i+1) - a%AI(i)
       do k = counter, counter+rowSize-1
          call c%append(a%A(k), i, a%AJ(k))
       end do
       counter = counter + rowSize
    end do
    counter = 1
    do i = 1, b%n
       rowSize = b%AI(i+1) - b%AI(i)
       do k = counter, counter+rowSize-1
          call c%append(-1.d0*b%A(k), i, b%AJ(k))
       end do
       counter = counter + rowSize
    end do
    call c%makeCRS
  end function sparse_sparse_sub
  !***************************************************
  ! transpose:
  !     Obtains the transpose of sparse matrix a
  !  
  ! Parameters:
  !     Input, a(Sparse)
  !     Output, b(Sparse)
  !***************************************************
  function transpose(a) result(b)
    implicit none
    class(Sparse), intent(in) :: a
    type(Sparse) :: b
    integer(ikind), dimension(a%n) :: colCounter
    integer(ikind) :: i, j, k, counter
    b = sparse(nnz = a%nnz, rows = a%n)
    colCounter = 0
    do i = 1, b%nnz
       colCounter(a%AJ(i)) = colCounter(a%AJ(i)) + 1
    end do
    allocate(b%AI(b%n+1))
    allocate(b%A(b%nnz))
    allocate(b%AJ(b%nnz))
    b%AI = 1
    do i = 2, b%n+1
       b%AI(i) = b%AI(i-1) + colCounter(i-1)
    end do
    counter = 1
    colCounter = 0
    do i = 1, b%n
       do k = a%AI(i), a%AI(i+1)-1
          j = a%AJ(k)
          b%AJ(b%AI(j)+colCounter(j)) = i
          b%A(b%AI(j)+colCounter(j)) = a%A(counter)
          colCounter(j) = colCounter(j) + 1
          counter = counter + 1
       end do
    end do
  end function transpose
  !***************************************************
  ! norm:
  !     Computes the frobenius-norm of a sparse matrix
  !     https://en.wikipedia.org/wiki/Matrix_norm#Frobenius_norm
  !  
  ! Parameters:
  !     Input, a(Sparse)
  !     Output, norm(Realrkind)
  !***************************************************
  real(rkind) function norm(a)
    implicit none
    class(Sparse), intent(in) :: a
    integer(ikind) :: i, j, counter
    norm = 0
    counter = 1
    do while(counter <= a%nnz)
       norm = norm + a%A(counter)*a%A(counter)
       counter = counter + 1
    end do
    norm = sqrt(norm)
  end function norm
  !***************************************************
  ! trace:
  !     Computes the sum of the elements of the
  !     diagonal of a sparse matrix
  !  
  ! Parameters:
  !     Input, a(Sparse)
  !     Output, trace(Realrkind)
  !***************************************************
  real(rkind) function trace(a)
    implicit none
    class(Sparse), intent(in) :: a
    integer(ikind) :: i, j
    trace = 0
    do i = 1, a%n
       do j = a%AI(i), a%AI(i+1)-1
          if(a%AJ(j) == i) then
             trace = trace + a%A(j)
             exit
          end if
       end do
    end do
  end function trace
  !***************************************************
  ! det:
  !     Computes the determinant
  !  
  ! Parameters:
  !     Input, a(Sparse)
  !     Output, det(Realrkind)
  !***************************************************
  real(rkind) function det(a)
    implicit none
    class(Sparse), intent(inout) :: a
    type(Sparse) :: m
    integer(ikind) :: i, j
    det = 1.
    m = lcholesky(a)
    do i = 1, m%n
       do j = m%AI(i), m%AI(i+1)-1
          if(m%AJ(j) == i) then
             det = det * m%A(j)
             exit
          end if
       end do
    end do
    det = det**2
  end function det
  !***************************************************
  ! id:
  !     Computes the identity matrix of order n
  !  
  ! Parameters:
  !     Input, n(integer(ikind))
  !     Output, id(Sparse)
  !***************************************************
  function id(n) result(a)
    implicit none
    type(Sparse) :: a
    integer(ikind), intent (in) :: n
    integer(ikind) :: i
    a = sparse(nnz = n, rows = n)
    do i = 1, n
       call a%append(1.d0, i, i)
    end do
    call a%makeCRS
  end function id
  !***************************************************
  ! lcholesky:
  !     Computes Incomplete Cholesky factorization
  !  
  ! Parameters:
  !     Input, a(Sparse)
  !     Output, L(Sparse)
  !***************************************************
  function lcholesky(a) result(L)
    implicit none
    class(Sparse), intent(inout) :: a
    type(Sparse) :: L
    integer(ikind) :: i, j, k
    real(rkind) :: adder1, adder2, m( a%n, a%n)
    L = sparse(nnz = (a%n**2+a%n)/2, rows = a%n)
    m = 0.
    do i = 1, a%n
       do j = 1, a%n
          adder1 = 0.
          do k = 1, j-1
             adder1 = adder1 + m(j,k)**2
          end do
          m(j,j) = sqrt(a%get(j,j)-adder1)
          if( i > j) then
             adder2 = 0.
             do k = 1 , j-1
                adder2 = adder2 + ( m(i,k) * m(j,k))
             end do
             m(i,j) = (1./m(j,j))*(a%get(i,j)-adder2)
          end if
       end do
    end do
    do i = 1, a%n
       call L%append( m(i,i), i, i)
       do j = 1, a%n
          if( i > j) then
             call L%append( m(i,j), i, j)
          end if
       end do
    end do
    call L%makeCRS
    return
  end function lcholesky
  !***************************************************
  ! gmres: 
  !   (generalized minimal residual method)
  !   Approximates the solution of a nonsymmetric system 
  !   of linear equations by the vector in a Krylov 
  !   subspace with minimal residual.
  !  
  ! Parameters:
  !     Input, A(Sparse)
  !     Input, rhs(realrkind)
  !     Output, x(realrkind)
  !     Input, ITR_MAX(ikind), the maximum number of (outer) 
  !                       iterations to take.
  !     Input, MR(ikind), the maximum number of (inner) iterations 
  !                       to take.  MR must be less than N.
  !     Input, TOL_ABS(rkind), an absolute tolerance applied to the
  !                       current residual.
  !     Input, TOL_REL(rkind), a relative tolerance comparing the
  !                       current residual to the initial residual.
  !***************************************************
  function gmres(A, rhs, itrMaxIn, mrIn, tolAbsIn, tolRelIn) result(x)
    implicit none
    class(Sparse), intent(in) :: A
    Integer(ikind), optional, value, intent(in) :: itrMaxIn
    Integer(ikind), optional, value, intent(in) :: mrIn
    Integer(ikind) :: itr_max
    Integer(ikind) :: mr
    Integer :: i, j, k, itr, itr_used
    Integer :: iw(A%n), jj, jrow, jw, k_copy, ua(A%n)
    Real(rkind), intent(in)  :: rhs(A%n)
    Real(rkind), parameter :: delta = 1.0D-03
    Real(rkind), optional, value, intent(in) :: tolAbsIn
    Real(rkind), optional, value, intent(in) :: tolRelIn
    Real(rkind) :: tol_abs, tol_rel
    Real(rkind) :: av, htmp, tl, mu
    Real(rkind) :: g1, g2, rho, rho_tol
    Real(rkind) :: x(A%n), r(A%n), l(A%AI(A%n+1)+1), w(A%n)
    real(rkind), dimension(:), allocatable :: c, g, s, y
    real(rkind), dimension(:,:), allocatable :: h, v
    logical, parameter :: verbose = .true.
    itr_max = 1000
    mr = 100
    tol_abs = 1d-15
    tol_rel = 1d-15
    if(present(itrMaxIn))then
       itr_max = itrMaxIn
    end if
    if(present(mrIn))then
       mr = mrIn
    end if
    if(present(tolAbsIn))then
       tol_abs = tolAbsIn
    end if
    if(present(tolRelIn))then
       tol_rel = tolRelIn
    end if

    allocate(h(mr+1,mr), v(A%n,mr+1))
    allocate(c(mr+1), g(mr+1), s(mr+1), y(mr+1))

    itr_used = 0

    call rearrange_cr ( A%n, A%nnz, A%AI, A%AJ, A%A )

    call diagonal_pointer_cr ( A%n, A%nnz, A%AI, A%AJ, ua )

    call ilu_cr (  A%n, A%nnz, A%AI, A%AJ, A%A, ua, l )
!!$    if ( verbose ) then
!!$       write ( *, '(a)' ) ' '
!!$       write ( *, '(a)' ) 'PMGMRES_ILU_CR'
!!$       write ( *, '(a,i4)' ) '  Number of unknowns = ', A%n
!!$    end if

    do itr = 1, itr_max
       call ax_cr (  A%n, A%nnz, A%AI, A%AJ, A%A, x, r )
       r(1:A%n) = rhs(1:A%n) - r(1:A%n)

       call lus_cr ( A%n, A%nnz, A%AI, A%AJ, l, ua, r, r )

       rho = sqrt ( dot_product ( r, r ) )

!!$       if ( verbose ) then
!!$          write ( *, '(a,i4,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
!!$       end if

       if ( itr == 1 ) then
          rho_tol = rho * tol_rel
       end if

       v(1:A%n,1) = r(1:A%n) / rho

       g(1) = rho
       g(2:mr+1) = 0.0D+00

       h(1:mr+1,1:mr) = 0.0D+00

       do k = 1, mr

          k_copy = k

          call ax_cr ( A%n, A%nnz, A%AI, A%AJ, A%A, v(1:A%n,k), v(1:A%n,k+1) ) 

          call lus_cr ( A%n, A%nnz, A%AI, A%AJ, l, ua, v(1:A%n,k+1), v(1:A%n,k+1) )

          av = sqrt ( dot_product ( v(1:A%n,k+1), v(1:A%n,k+1) ) )

          do j = 1, k
             h(j,k) = dot_product ( v(1:A%n,k+1), v(1:A%n,j) )
             v(1:A%n,k+1) = v(1:A%n,k+1) - v(1:A%n,j) * h(j,k)
          end do

          h(k+1,k) = sqrt ( dot_product ( v(1:A%n,k+1), v(1:A%n,k+1) ) )

          if ( ( av + delta * h(k+1,k)) == av ) then
             do j = 1, k
                htmp = dot_product ( v(1:A%n,k+1), v(1:A%n,j) )
                h(j,k) = h(j,k) + htmp
                v(1:A%n,k+1) = v(1:A%n,k+1) - htmp * v(1:A%n,j)
             end do
             h(k+1,k) = sqrt ( dot_product ( v(1:A%n,k+1), v(1:A%n,k+1) ) )
          end if

          if ( h(k+1,k) /= 0.0D+00 ) then
             v(1:A%n,k+1) = v(1:A%n,k+1) / h(k+1,k)
          end if

          if ( 1 < k ) then
             y(1:k+1) = h(1:k+1,k)
             do j = 1, k - 1
                call mult_givens ( c(j), s(j), j, y )
             end do
             h(1:k+1,k) = y(1:k+1)
          end if

          mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

          c(k) = h(k,k) / mu
          s(k) = -h(k+1,k) / mu
          h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
          h(k+1,k) = 0.0D+00
          call mult_givens ( c(k), s(k), k, g )

          rho = abs ( g(k+1) )

          itr_used = itr_used + 1

!!$          if ( verbose ) then
!!$             write ( *, '(a,i4,a,g14.6)' ) '  K = ', k, '  Residual = ', rho
!!$          end if

          if ( rho <= rho_tol .and. rho <= tol_abs ) then
             exit
          end if

       end do

       k = k_copy - 1

       y(k+1) = g(k+1) / h(k+1,k+1)

       do i = k, 1, -1
          y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
       end do

       do i = 1, A%n
          x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
       end do

       if ( rho <= rho_tol .and. rho <= tol_abs ) then
          exit
       end if

    end do

    if ( verbose ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'PMGMRES_ILU_CR:'
       write ( *, '(a,i6)' ) '  Iterations = ', itr_used
       write ( *, '(a,g14.6)' ) '  Final residual = ', rho
    end if

    return    

  end function gmres
  !***************************************************
  ! Inverse:
  !    Obtains the inverse of sparse matrix A
  !  
  ! Parameters:
  !     Input, A(Sparse)
  !     Output, B(Sparse)
  !***************************************************
  function inverse(A) result(B)
    implicit none
    class(Sparse), intent(in) :: A
    type(Sparse) :: B
    real(rkind) :: y(A%n), x(A%n)
    integer(ikind) :: i, j
    B = sparse(nnz = A%n**2, rows = A%n)
    do j = 1,A%n
       y = 0.
       y(j) = 1.
       x = CGOMP(A,y)
       do i = 1, A%n
          if(abs(x(i)).gt.1d-5)then
             call B%append(x(i), i, j)
          end if
       end do
    end do
    call B%makeCRS
    return
  end function inverse
  !***************************************************
  ! inverseGMRESD:
  !    Obtains the inverse of sparse matrix A
  !    (Global Minimal Residual descent algorithm)
  ! Parameters:
  !     Input, A(Sparse)
  !     Output, B(Sparse)
  !***************************************************
  function inverseGMRESD(A) result(M)
    implicit none
    class(Sparse), intent(in) :: A
    type(Sparse) :: B
    type(Sparse) :: C
    type(Sparse) :: G
    type(Sparse) :: M
    real(rkind) :: alpha
    alpha = 1.
    M = transpose(A)
    do while (abs(alpha) > 1e-30)
       C = A * M
       G = Id(A%n) - C
       B = A * G
       alpha = trace(transpose(G)*B)/(norm(B))**2
       M = M + alpha * G
    end do
    return
  end function inverseGMRESD
  !***************************************************
  ! jacobiEigen:
  !    Obtains the eigenvalues and eigenvectors of a
  !    real symmetric matrix A, using Rutishauser's
  !    modfications of the classical Jacobi rotation
  !    method with threshold pivoting.
  !
  ! Parameters:
  !     Input, A(Sparse)
  !     Output, Eigenvec(realrkind), the matrix of eigenvectors
  !     Output, Eigenval(realrkind), the eigenvalues.
  !***************************************************
  subroutine jacobiEigen(A_input, Eigenval, Eigenvec) 
    implicit none
    class(Sparse), intent(in) :: A_input
    real(rkind), intent(out) :: Eigenval(A_input%n),Eigenvec(A_input%n,A_input%n)
    integer(ikind), parameter :: it_max = 1000
    integer(ikind) :: i, j, k, l, m, p, q, it_num, rot_num
    real(rkind) :: a(A_input%n,A_input%n), bw(A_input%n),  w(A_input%n), zw(A_input%n) 
    real(rkind) :: h, s, t, tau, term, termp, termq, theta, c, g, thresh, gapq
    a = 0.
    Do i = 1, A_input%n
       Do j = A_input%AI(i), A_input%AI(i+1)-1
          a(i, A_input%AJ(j))  = A_input%A(j)
       End Do
    End Do
    Eigenvec = 0.
    do i = 1, A_input%n
       Eigenvec(i,i) = 1.
    end do
    do i = 1, A_input%n
       Eigenval(i) = a(i,i)
    end do
    bw(1:A_input%n) = Eigenval(1:A_input%n)
    zw(1:A_input%n) = 0.
    it_num = 0
    rot_num = 0
    do while ( it_num < it_max )
       it_num = it_num + 1
       !  The convergence threshold is based on the size of the elements in
       !  the strict upper triangle of the matrix.
       thresh = 0.
       do j = 1, A_input%n
          do i = 1, j - 1
             thresh = thresh + a(i,j) ** 2
          end do
       end do
       thresh = sqrt ( thresh ) / (4.*A_input%n)
       if ( thresh == 0. ) then
          exit 
       end if
       do p = 1, A_input%n
          do q = p + 1, A_input%n
             gapq = 10. * abs ( a(p,q) )
             termp = gapq + abs ( Eigenval(p) )
             termq = gapq + abs ( Eigenval(q) )
             !  Annihilate tiny offdiagonal elements.
             if ( 4 < it_num .and. &
                  termp == abs ( Eigenval(p) ) .and. &
                  termq == abs ( Eigenval(q) ) ) then
                a(p,q) = 0.
                !  Otherwise, apply a rotation.
             else if ( thresh <= abs ( a(p,q) ) ) then
                h = Eigenval(q) - Eigenval(p)
                term = abs ( h ) + gapq
                if ( term == abs ( h ) ) then
                   t = a(p,q) / h
                else
                   theta = 0.5 * h / a(p,q)
                   t = 1. / ( abs ( theta ) + sqrt ( 1. + theta * theta ) )
                   if ( theta < 0. ) then 
                      t = - t
                   end if
                end if
                c = 1. / sqrt ( 1. + t * t )
                s = t * c
                tau = s / ( 1. + c )
                h = t * a(p,q)
                !  Accumulate corrections to diagonal elements.
                zw(p) = zw(p) - h                  
                zw(q) = zw(q) + h
                Eigenval(p) = Eigenval(p) - h
                Eigenval(q) = Eigenval(q) + h
                a(p,q) = 0.
                !  Rotate, using information from the upper triangle of A only.
                do j = 1, p - 1
                   g = a(j,p)
                   h = a(j,q)
                   a(j,p) = g - s * ( h + g * tau )
                   a(j,q) = h + s * ( g - h * tau )
                end do
                do j = p + 1, q - 1
                   g = a(p,j)
                   h = a(j,q)
                   a(p,j) = g - s * ( h + g * tau )
                   a(j,q) = h + s * ( g - h * tau )
                end do
                do j = q + 1, A_input%n
                   g = a(p,j)
                   h = a(q,j)
                   a(p,j) = g - s * ( h + g * tau )
                   a(q,j) = h + s * ( g - h * tau )
                end do
                !  Accumulate information in the eigenvector matrix.
                do j = 1, A_input%n
                   g = Eigenvec(j,p)
                   h = Eigenvec(j,q)
                   Eigenvec(j,p) = g - s * ( h + g * tau )
                   Eigenvec(j,q) = h + s * ( g - h * tau )
                end do
                rot_num = rot_num + 1
             end if
          end do
       end do
       bw(1:A_input%n) = bw(1:A_input%n) + zw(1:A_input%n)
       Eigenval(1:A_input%n) = bw(1:A_input%n)
       zw(1:A_input%n) = 0.
    end do
    !  Restore upper triangle of input matrix.
    do j = 1, A_input%n
       do i = 1, j - 1
          a(i,j) = a(j,i)
       end do
    end do
    !  Ascending sort the eigenvalues and eigenvectors.
    do k = 1, A_input%n - 1
       m = k
       do l = k + 1, A_input%n
          if ( Eigenval(l) < Eigenval(m) ) then
             m = l
          end if
       end do
       if ( m /= k ) then
          t    = Eigenval(m)
          Eigenval(m) = Eigenval(k)
          Eigenval(k) = t
          w(1:A_input%n)   = Eigenvec(1:A_input%n,m)
          Eigenvec(1:A_input%n,m) = Eigenvec(1:A_input%n,k)
          Eigenvec(1:A_input%n,k) = w(1:A_input%n)
       end if
    end do
  end subroutine jacobiEigen
  subroutine rearrange_cr ( n, nz_num, ia, ja, a )

    !*****************************************************************************80
    !
    !! REARRANGE_CR sorts a sparse compressed row matrix.
    !
    !  Discussion:
    !
    !    This routine guarantees that the entries in the CR matrix
    !    are properly sorted.
    !
    !    After the sorting, the entries of the matrix are rearranged in such
    !    a way that the entries of each column are listed in ascending order
    !    of their column values.
    !
    !    The matrix A is assumed to be stored in compressed row format.  Only
    !    the nonzero entries of A are stored.  The vector JA stores the
    !    column index of the nonzero value.  The nonzero values are sorted
    !    by row, and the compressed row vector IA then has the property that
    !    the entries in A and JA that correspond to row I occur in indices
    !    IA(I) through IA(I+1)-1.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 July 2007
    !
    !  Author:
    !
    !    Original C version by Lili Ju.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    !    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    !    Charles Romine, Henk van der Vorst,
    !    Templates for the Solution of Linear Systems:
    !    Building Blocks for Iterative Methods,
    !    SIAM, 1994.
    !    ISBN: 0898714710,
    !    LC: QA297.8.T45.
    !
    !    Tim Kelley,
    !    Iterative Methods for Linear and Nonlinear Equations,
    !    SIAM, 2004,
    !    ISBN: 0898713528,
    !    LC: QA297.8.K45.
    !
    !    Yousef Saad,
    !    Iterative Methods for Sparse Linear Systems,
    !    Second Edition,
    !    SIAM, 2003,
    !    ISBN: 0898715342,
    !    LC: QA188.S17.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the order of the system.
    !
    !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
    !
    !    Input, integer ( kind = 4 ) IA(N+1), the compressed row indices.
    !
    !    Input/output, integer ( kind = 4 ) JA(NZ_NUM), the column indices.
    !    On output, these may have been rearranged by the sorting.
    !
    !    Input/output, real ( kind = 8 ) A(NZ_NUM), the matrix values.  On output,
    !    the matrix values may have been moved somewhat because of the sorting.
    !
    implicit none

    integer ( kind = 4 ) n
    integer ( kind = 4 ) nz_num

    real ( kind = 8 ) a(nz_num)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ia(n+1)
    integer ( kind = 4 ) i4temp
    integer ( kind = 4 ) ja(nz_num)
    integer ( kind = 4 ) k
    integer ( kind = 4 ) l
    real ( kind = 8 ) r8temp

    do i = 1, n

       do k = ia(i), ia(i+1) - 2
          do l = k + 1, ia(i+1) - 1

             if ( ja(l) < ja(k) ) then
                i4temp = ja(l)
                ja(l)  = ja(k)
                ja(k)  = i4temp

                r8temp = a(l)
                a(l)   = a(k)
                a(k)   = r8temp
             end if

          end do
       end do

    end do

    return
  end subroutine rearrange_cr
  subroutine diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

    !*****************************************************************************80
    !
    !! DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.
    !
    !  Discussion:
    !
    !    The matrix A is assumed to be stored in compressed row format.  Only
    !    the nonzero entries of A are stored.  The vector JA stores the
    !    column index of the nonzero value.  The nonzero values are sorted
    !    by row, and the compressed row vector IA then has the property that
    !    the entries in A and JA that correspond to row I occur in indices
    !    IA[I] through IA[I+1]-1.
    !
    !    The array UA can be used to locate the diagonal elements of the matrix.
    !
    !    It is assumed that every row of the matrix includes a diagonal element,
    !    and that the elements of each row have been ascending sorted.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 July 2007
    !
    !  Author:
    !
    !    Original C version by Lili Ju.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the order of the system.
    !
    !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
    !
    !    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
    !    indices of the matrix values.  The row vector has been compressed.
    !    On output, the order of the entries of JA may have changed because of
    !    the sorting.
    !
    !    Output, integer ( kind = 4 ) UA(N), the index of the diagonal element
    !    of each row.
    !
    implicit none

    integer ( kind = 4 ) n
    integer ( kind = 4 ) nz_num

    integer ( kind = 4 ) i
    integer ( kind = 4 ) ia(n+1)
    integer ( kind = 4 ) k
    integer ( kind = 4 ) ja(nz_num)
    integer ( kind = 4 ) ua(n)

!!$    ua(1:n) = -1

    do i = 1, n
       do k = ia(i), ia(i+1) - 1
          if ( ja(k) == i ) then
             ua(i) = k
          end if
       end do
    end do

    return
  end subroutine diagonal_pointer_cr
  subroutine ilu_cr ( n, nz_num, ia, ja, a, ua, l )

    !*****************************************************************************80
    !
    !! ILU_CR computes the incomplete LU factorization of a matrix.
    !
    !  Discussion:
    !
    !    The matrix A is assumed to be stored in compressed row format.  Only
    !    the nonzero entries of A are stored.  The vector JA stores the
    !    column index of the nonzero value.  The nonzero values are sorted
    !    by row, and the compressed row vector IA then has the property that
    !    the entries in A and JA that correspond to row I occur in indices
    !    IA(I) through IA(I+1)-1.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    27 July 2007
    !
    !  Author:
    !
    !    Original C version by Lili Ju.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the order of the system.
    !
    !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
    !
    !    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
    !    indices of the matrix values.  The row vector has been compressed.
    !
    !    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
    !
    !    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
    !    of each row.
    !
    !    Output, real ( kind = 8 ) L(NZ_NUM), the ILU factorization of A.
    !
    implicit none

    integer ( kind = 4 ) n
    integer ( kind = 4 ) nz_num

    real ( kind = 8 ) a(nz_num)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ia(n+1)
    integer ( kind = 4 ) iw(n)
    integer ( kind = 4 ) j
    integer ( kind = 4 ) ja(nz_num)
    integer ( kind = 4 ) jj
    integer ( kind = 4 ) jrow
    integer ( kind = 4 ) jw
    integer ( kind = 4 ) k
    real ( kind = 8 ) l(nz_num)
    real ( kind = 8 ) tl
    integer ( kind = 4 ) ua(n)
    !
    !  Copy A.
    !
    l(1:nz_num) = a(1:nz_num)

    do i = 1, n
       !
       !  IW points to the nonzero entries in row I.
       !
       iw(1:n) = -1

       do k = ia(i), ia(i+1) - 1
          iw(ja(k)) = k
       end do

       do j = ia(i), ia(i+1) - 1
          jrow = ja(j)
          if ( i <= jrow ) then
             exit
          end if
          tl = l(j) * l(ua(jrow))
          l(j) = tl
          do jj = ua(jrow) + 1, ia(jrow+1) - 1
             jw = iw(ja(jj))
             if ( jw /= -1 ) then
                l(jw) = l(jw) - tl * l(jj)
             end if
          end do
       end do

       ua(i) = j

       if ( jrow /= i ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'ILU_CR - Fatal error!'
          write ( *, '(a)' ) '  JROW ~= I'
          write ( *, '(a,i8)' ) '  JROW = ', jrow
          write ( *, '(a,i8)' ) '  I    = ', i
          stop
       end if

       if ( l(j) == 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'ILU_CR - Fatal error!'
          write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
          write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
          stop
       end if

       l(j) = 1.0D+00 / l(j)

    end do

    l(ua(1:n)) = 1.0D+00 / l(ua(1:n))

    return
  end subroutine ilu_cr
  subroutine ax_cr ( n, nz_num, ia, ja, a, x, w )

    !*****************************************************************************80
    !
    !! AX_CR computes A*x for a matrix stored in sparse compressed row form.
    !
    !  Discussion:
    !
    !    The Sparse Compressed Row storage format is used.
    !
    !    The matrix A is assumed to be sparse.  To save on storage, only
    !    the nonzero entries of A are stored.  The vector JA stores the
    !    column index of the nonzero value.  The nonzero values are sorted
    !    by row, and the compressed row vector IA then has the property that
    !    the entries in A and JA that correspond to row I occur in indices
    !    IA[I] through IA[I+1]-1.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 July 2007
    !
    !  Author:
    !
    !    Original C version by Lili Ju.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    !    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    !    Charles Romine, Henk van der Vorst,
    !    Templates for the Solution of Linear Systems:
    !    Building Blocks for Iterative Methods,
    !    SIAM, 1994.
    !    ISBN: 0898714710,
    !    LC: QA297.8.T45.
    !
    !    Tim Kelley,
    !    Iterative Methods for Linear and Nonlinear Equations,
    !    SIAM, 2004,
    !    ISBN: 0898713528,
    !    LC: QA297.8.K45.
    !
    !    Yousef Saad,
    !    Iterative Methods for Sparse Linear Systems,
    !    Second Edition,
    !    SIAM, 2003,
    !    ISBN: 0898715342,
    !    LC: QA188.S17.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the order of the system.
    !
    !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
    !
    !    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
    !    indices of the matrix values.  The row vector has been compressed.
    !
    !    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
    !
    !    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
    !
    !    Output, real ( kind = 8 ) W(N), the value of A*X.
    !
    implicit none

    integer ( kind = 4 ) n
    integer ( kind = 4 ) nz_num

    real ( kind = 8 ) a(nz_num)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ia(n+1)
    integer ( kind = 4 ) ja(nz_num)
    integer ( kind = 4 ) k
    integer ( kind = 4 ) k1
    integer ( kind = 4 ) k2
    real ( kind = 8 ) w(n)
    real ( kind = 8 ) x(n)

    w(1:n) = 0.0D+00

    do i = 1, n
       k1 = ia(i)
       k2 = ia(i+1) - 1
       w(i) = w(i) + dot_product ( a(k1:k2), x(ja(k1:k2)) )
    end do

    return
  end subroutine ax_cr
  subroutine lus_cr ( n, nz_num, ia, ja, l, ua, r, z )

    !*****************************************************************************80
    !
    !! LUS_CR applies the incomplete LU preconditioner.
    !
    !  Discussion:
    !
    !    The linear system M * Z = R is solved for Z.  M is the incomplete
    !    LU preconditioner matrix, and R is a vector supplied by the user.
    !    So essentially, we're solving L * U * Z = R.
    !
    !    The matrix A is assumed to be stored in compressed row format.  Only
    !    the nonzero entries of A are stored.  The vector JA stores the
    !    column index of the nonzero value.  The nonzero values are sorted
    !    by row, and the compressed row vector IA then has the property that
    !    the entries in A and JA that correspond to row I occur in indices
    !    IA(I) through IA(I+1)-1.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 July 2007
    !
    !  Author:
    !
    !    Original C version by Lili Ju.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the order of the system.
    !
    !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
    !
    !    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
    !    indices of the matrix values.  The row vector has been compressed.
    !
    !    Input, real ( kind = 8 ) L(NZ_NUM), the matrix values.
    !
    !    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
    !    of each row.
    !
    !    Input, real ( kind = 8 ) R(N), the right hand side.
    !
    !    Output, real ( kind = 8 ) Z(N), the solution of the system M * Z = R.
    !
    implicit none

    integer ( kind = 4 ) n
    integer ( kind = 4 ) nz_num

    integer ( kind = 4 ) i
    integer ( kind = 4 ) ia(n+1)
    integer ( kind = 4 ) j
    integer ( kind = 4 ) ja(nz_num)
    real ( kind = 8 ) l(nz_num)
    real ( kind = 8 ) r(n)
    integer ( kind = 4 ) ua(n)
    real ( kind = 8 ) w(n)
    real ( kind = 8 ) z(n)
    !
    !  Copy R in.
    !
    w(1:n) = r(1:n)
    !
    !  Solve L * w = w where L is unit lower triangular.
    !
    do i = 2, n
       do j = ia(i), ua(i) - 1
          w(i) = w(i) - l(j) * w(ja(j))
       end do
    end do
    !
    !  Solve U * w = w, where U is upper triangular.
    !
    do i = n, 1, -1
       do j = ua(i) + 1, ia(i+1) - 1
          w(i) = w(i) - l(j) * w(ja(j))
       end do
       w(i) = w(i) / l(ua(i))
    end do
    !
    !  Copy Z out.
    !
    z(1:n) = w(1:n)

    return
  end subroutine lus_cr
  subroutine mult_givens ( c, s, k, g )

    !*****************************************************************************80
    !
    !! MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
    !
    !  Discussion:
    !
    !    In order to make it easier to compare this code with the Original C,
    !    the vector indexing is 0-based.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    08 August 2006
    !
    !  Author:
    !
    !    Original C version by Lili Ju.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    !    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    !    Charles Romine, Henk van der Vorst,
    !    Templates for the Solution of Linear Systems:
    !    Building Blocks for Iterative Methods,
    !    SIAM, 1994.
    !    ISBN: 0898714710,
    !    LC: QA297.8.T45.
    !
    !    Tim Kelley,
    !    Iterative Methods for Linear and Nonlinear Equations,
    !    SIAM, 2004,
    !    ISBN: 0898713528,
    !    LC: QA297.8.K45.
    !
    !    Yousef Saad,
    !    Iterative Methods for Sparse Linear Systems,
    !    Second Edition,
    !    SIAM, 2003,
    !    ISBN: 0898715342,
    !    LC: QA188.S17.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) C, S, the cosine and sine of a Givens
    !    rotation.
    !
    !    Input, integer ( kind = 4 ) K, indicates the location of the first
    !    vector entry.
    !
    !    Input/output, real ( kind = 8 ) G(1:K+1), the vector to be modified.
    !    On output, the Givens rotation has been applied to entries G(K) and G(K+1).
    !
    implicit none

    integer ( kind = 4 ) k

    real ( kind = 8 ) c
    real ( kind = 8 ) g(1:k+1)
    real ( kind = 8 ) g1
    real ( kind = 8 ) g2
    real ( kind = 8 ) s

    g1 = c * g(k) - s * g(k+1)
    g2 = s * g(k) + c * g(k+1)

    g(k)   = g1
    g(k+1) = g2

    return
  end subroutine mult_givens

  function biCGOMP(A, rhs) result(x)
    implicit none
    class(Sparse), intent(in) :: A
    real(rkind), dimension(:), intent(in) :: rhs
    real(rkind), dimension(:), allocatable :: x
    integer(ikind) :: i
    integer(ikind), dimension(1) :: auxVect
    integer(ikind), dimension(A%n) :: diagPtr
    real(rkind), dimension(A%n) :: diag 

    allocate(x(A%n))
    x = 0.d0

    call diagonal_pointer_cr(A%n, A%nnz, A%AI, A%AJ, diagPtr)

    do i = 1, A%n
       diag(i) = A%A(diagPtr(i))
    end do

    call biCG(A%A, A%AJ, A%AI, diag, x, rhs, rhs, auxVect, A%n, 0)


  end function biCGOMP

  subroutine biCG(spMtx, spIdx, spRowptr, diagMtx, x, b, x_fix, x_fixIdx, npoin, nfix)
    !--------------------------------------------------------------------------------------------------
    !					 Calcula x dado A.x = b, mediante gradiente biconjugado
    !--------------------------------------------------------------------------------------------------
    implicit none
    integer(ikind) :: npoin, nfix, k
    integer(ikind) :: spRowptr(npoin + 1), spIdx(spRowptr(npoin + 1)-1)
    integer(ikind) :: x_fixIdx(nfix)
    real(rkind) :: spMtx(spRowptr(npoin + 1)-1), x(npoin), b(npoin), diagMtx(npoin)
    real(rkind)	:: x_fix(nfix)
    real(rkind) :: tol, err_old, err_new, alfa, beta, py
    real(rkind), allocatable, dimension(:) :: y, p, r, z

    if(.not.allocated(y)) allocate(y(npoin))
    if(.not.allocated(p)) allocate(p(npoin))
    if(.not.allocated(r)) allocate(r(npoin))
    if(.not.allocated(z)) allocate(z(npoin))

    k = 0
    tol = 1.d-15

    call copy1(nfix, 1.d0, x_fix, x_fixIdx, npoin, x)
    call SpMV(spMtx, spIdx, spRowptr, x, y, npoin, spRowptr(npoin+1)-1)
    call copy2(npoin, 1.d30, x, nfix, x_fixIdx, y)
    call vecsum(npoin, -1.d0, y, b, r)
    call assign2(npoin, r, nfix, x_fixIdx, 0.d0)

    if(vecdot(npoin, r, r) < tol) return

    call vecdiv(npoin, r, diagMtx, p)
    err_new = vecdot(npoin, r, p)
    call SpMV(spMtx, spIdx, spRowptr, p, y, npoin, spRowptr(npoin+1)-1)
    call copy2(npoin, 1.d30, p, nfix, x_fixIdx, y)
    py = vecdot(npoin, p, y)
    !DEBUGGGGG
    alfa = err_new/py

    call vecsum(npoin, alfa, p, x, x)
    err_old = err_new

    do while(dabs(err_old) > tol .and. k < 1000)
       k = k + 1
       call vecsum(npoin, -alfa, y, r, r)
       call vecdiv(npoin, r, diagMtx, z)
       err_new = vecdot(npoin, r, z)
       beta = err_new/err_old
       call vecsum(npoin, beta, p, z, p)
       call SpMV(spMtx, spIdx, spRowptr, p, y, npoin, spRowptr(npoin+1)-1)
       call copy2(npoin, 1.d30, p, nfix, x_fixIdx, y)
       py = vecdot(npoin, p, y)
       !DEBUGGGGG
       alfa = err_new/py
       call vecsum(npoin, alfa, p, x, x)
       err_old = err_new
    end do
  end subroutine biCG

  subroutine vecdiv(n, x, y, z)
    !---------------------------------
    !		Calcula z = x/y
    !---------------------------------
    implicit none
    integer(ikind) :: n, i
    real(rkind) :: x(n), y(n), z(n)

    !$OMP PARALLEL DO PRIVATE(i)
    do i = 1, n
       z(i) = x(i)/y(i)	
    end do
    !$OMP END PARALLEL DO
  end subroutine vecdiv

  subroutine vecsum(n, alfa, x, y, z)
    !-----------------------------------
    !		Calcula z = alfa*x + y
    !-----------------------------------
    implicit none
    integer(ikind) :: n, i
    real(rkind) :: x(n), y(n), z(n), alfa

    !$OMP PARALLEL DO PRIVATE(i)
    do i = 1, n
       z(i) = alfa*x(i) + y(i)
    end do
    !$OMP END PARALLEL DO
  end subroutine vecsum

  subroutine assign(n, y, scal)
    !-----------------------------
    !		y(:) = scal
    !-----------------------------
    implicit none
    integer(ikind) :: n, i
    real(rkind) :: scal, y(n)

    !$OMP PARALLEL DO PRIVATE(i)
    do i = 1, n
       y(i) = scal
    end do
    !$OMP END PARALLEL DO
  end subroutine assign

  subroutine assign2(n, y, m, idx, scal)
    implicit none
    integer(ikind) :: n, m, i
    real(rkind) :: y(n), scal
    integer(ikind) :: idx(m)

    !$OMP PARALLEL DO PRIVATE(i)
    do i = 1, m
       y(idx(i)) = scal
    end do
    !$OMP END PARALLEL DO
  end subroutine assign2

  subroutine copy1(m, alfa, x, idx, n, y)
    !-----------------------------------------
    !			y(idx(:)) = alfa*x
    !-----------------------------------------
    implicit none
    integer(ikind) :: n, m, i
    real(rkind) :: x(m), y(n), alfa
    integer(ikind) :: idx(m)
    !$OMP PARALLEL DO PRIVATE(i)
    do i=1, m
       y(idx(i)) = alfa*x(i)
    end do
    !$OMP END PARALLEL DO
  end subroutine copy1

  subroutine copy2(n, alfa, x, m, idx, y)
    !-----------------------------------------
    !			y(idx(:)) = alfa*x(idx(:))
    !-----------------------------------------
    implicit none
    integer(ikind) :: n, m, i
    real(rkind) :: x(n), y(n), alfa
    integer(ikind) :: idx(m)

    !$OMP PARALLEL DO PRIVATE(i)
    do i=1, m
       y(idx(i)) = alfa*x(idx(i))
    end do
    !$OMP END PARALLEL DO
  end subroutine copy2

  real(rkind) function vecdot(n, x, y)
    !-----------------------------------
    !	Devuelve producto punto <x,y>
    !-----------------------------------
    implicit none
    integer(ikind) :: n, i
    real(rkind) :: x(n), y(n), res
    res = 0.d0

    !$OMP PARALLEL DO REDUCTION(+:res)
    do i = 1, n
       res = res + x(i)*y(i)
    end do
    !$OMP END PARALLEL DO

    vecdot = res
  end function vecdot

  subroutine SpMV(spMtx, spIdx, spRowptr, v, y, npoin, npos)
    !---------------------------------------------------------
    !			Multiplicacion sparse y = M.v
    !---------------------------------------------------------
    implicit none
    integer(ikind) :: npoin, npos
    integer(ikind) :: i, j
    real(rkind) :: spMtx(npos), v(npoin), y(npoin), dot
    integer(ikind) :: spIdx(npos), spRowptr(npoin+1)

    !$OMP PARALLEL DO PRIVATE(i, j, dot)
    do i = 1, npoin
       dot = 0.d0
       do j = spRowptr(i), spRowptr(i+1)-1
          dot = dot + spMtx(j)*v(spIdx(j))
       end do
       y(i) = dot
    end do
    !$OMP END PARALLEL DO 
  end subroutine SpMV
  function bicGrad(A,rhs) result(dof)
    implicit none
    class(sparse), intent(in) :: A
    real(rkind), intent(in) :: rhs(A%n)
    real(rkind)  :: dof(A%n)
    integer(ikind) :: diag(A%n), AI(A%nnz)
    integer(ikind) :: i, j, k
    !
    !  Solve the linear system using the conjugate gradient method.
    !
    do i = 1, size(A%AI)-1
       do j = A%AI(i), A%AI(i+1)-1
          AI(j) = i
       end do
    end do
    !
    !  Index the diagonal elements for use by the CG solver.
    !
    call diag_index (A%nnz, AI, A%AJ,A%n, diag )
    do k = 1, A%nnz
       if ( AI(k) < 1 .or. A%n < AI(k) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FEM2D_POISSON_CG - Fatal error!'
          write ( *, '(a)' ) '  Illegal IA(K)'
          stop
       end if
       if ( A%AJ(k) < 1 .or. A%n < A%AJ(k) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FEM2D_POISSON_CG - Fatal error!'
          write ( *, '(a)' ) '  Illegal JA(K)'
          stop
       end if
    end do
    call solve_cg (A%n, diag, A%nnz, AI, A%AJ, A%A, rhs, dof )
  end function bicGrad
  subroutine diag_index ( m, ia, ja, n, diag )
    !
    !  DIAG_INDEX determines where the diagonal matrix entries are stored.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 January 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, the number of adjacencies.
    !
    !    Input, integer ( kind = 4 ) IA(M), JA(M), the row and column indices 
    !    of adjacencies.
    !
    !    Input, integer ( kind = 4 ) N, the number of nodes.
    !
    !    Output, integer ( kind = 4 ) DIAG(N), contains for each index 1 <= I <= N, 
    !    the unique index J such that IA[J] = JA[J] = I.
    !
    implicit none
    integer ( kind = 4 ) m
    integer ( kind = 4 ) n
    integer ( kind = 4 ) diag(n)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ia(m)
    integer ( kind = 4 ) j
    integer ( kind = 4 ) ja(m)
    diag(1:n) = -1
    do j = 1, m
       if ( ia(j) == ja(j) ) then
          i = ia(j)
          if ( diag(i) /= -1 ) then
             write ( *, '(a)' ) ' '
             write ( *, '(a)' ) 'DIAG_INDEX - Fatal error!'
             write ( *, '(a)' ) '  Multiple occurrences of diagonal pairs.'
             write ( *, '(a,i4,a,i4,a,i4,a)' ) &
                  '  IA(', j, ') = JA(', j, ') = ', ia(j), ' and'
             write ( *, '(a,i4,a,i4,a,i4)'   ) &
                  '  IA(', diag(i), ') = JA(', diag(i), ') = ', ia(j)
             stop
          end if
          diag(i) = j
       end if
    end do
    do i = 1,  n
       if ( diag(i) == -1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DIAG_INDEX - Fatal error!'
          write ( *, '(a,i4,a)' ) '  DIAG(', i, ') = -1.'
          stop
       end if
    end do
    return
  end subroutine diag_index
  subroutine cg_rc ( n, b, x, r, z, p, q, job )
    !
    !  CG_RC is a reverse communication conjugate gradient routine.
    !
    !  Discussion:
    !
    !    This routine seeks a solution of the linear system A*x=b
    !    where b is a given right hand side vector, A is an n by n
    !    symmetric positive definite matrix, and x is an unknown vector
    !    to be determined.
    !
    !    Under the assumptions that the matrix A is large and sparse,
    !    the conjugate gradient method may provide a solution when
    !    a direct approach would be impractical because of excessive
    !    requirements of storage or even of time.
    !
    !    The conjugate gradient method presented here does not require the 
    !    user to store the matrix A in a particular way.  Instead, it only 
    !    supposes that the user has a way of calculating
    !      y = alpha * A * x + b * y
    !    and of solving the preconditioned linear system
    !      M * x = b
    !    where M is some preconditioning matrix, which might be merely
    !    the identity matrix, or a diagonal matrix containing the
    !    diagonal entries of A.
    !
    !    This routine was extracted from the "templates" package.
    !    There, it was not intended for direct access by a user;
    !    instead, a higher routine called "cg()" was called once by
    !    the user.  The cg() routine then made repeated calls to 
    !    cgrevcom() before returning the result to the user.
    !
    !    The reverse communication feature of cgrevcom() makes it, by itself,
    !    a very powerful function.  It allows the user to handle issues of
    !    storage and implementation that would otherwise have to be
    !    mediated in a fixed way by the function argument list.  Therefore,
    !    this version of cgrecom() has been extracted from the templates
    !    library and documented as a stand-alone procedure.
    !
    !    The user sets the value of JOB to 1 before the first call,
    !    indicating the beginning of the computation, and to the value of
    !    2 thereafter, indicating a continuation call.  
    !    The output value of JOB is set by cgrevcom(), which
    !    will return with an output value of JOB that requests a particular
    !    new action from the user.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    12 January 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    !    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    !    Charles Romine, Henk van der Vorst,
    !    Templates for the Solution of Linear Systems:
    !    Building Blocks for Iterative Methods,
    !    SIAM, 1994,
    !    ISBN: 0898714710,
    !    LC: QA297.8.T45.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the dimension of the matrix.
    !
    !    Input, real ( kind = 8 ) B(N), the right hand side vector.
    !
    !    Input/output, real ( kind = 8 ) X(N).  On first call, the user 
    !    should store an initial guess for the solution in X.  On return with
    !    JOB = 4, X contains the latest solution estimate.
    !
    !    Input/output, real ( kind = 8 ) R(N), Z(N), P(N), Q(N),
    !    information used by the program during the calculation.  The user
    !    does not need to initialize these vectors.  However, specific
    !    return values of JOB may require the user to carry out some computation
    !    using data in some of these vectors.
    !
    !    Input/output, integer ( kind = 4 ) JOB, communicates the task to be done.
    !    The user needs to set the input value of JOB to 1, before the first call,
    !    and then to 2 for every subsequent call for the given problem.
    !    The output value of JOB indicates the requested user action.  
    !    * JOB = 1, compute Q = A * P;
    !    * JOB = 2: solve M*Z=R, where M is the preconditioning matrix;
    !    * JOB = 3: compute R = R - A * X;
    !    * JOB = 4: check the residual R for convergence.  
    !               If satisfactory, terminate the iteration.
    !               If too many iterations were taken, terminate the iteration.
    !
    implicit none
    integer ( kind = 4 ) n
    real ( kind = 8 ) alpha
    real ( kind = 8 ) b(n)
    real ( kind = 8 ) beta
    integer ( kind = 4 ) iter
    integer ( kind = 4 ) job
    real ( kind = 8 ) p(n)
    real ( kind = 8 ) pdotq
    real ( kind = 8 ) q(n)
    real ( kind = 8 ) r(n)
    real ( kind = 8 ) rho
    real ( kind = 8 ) rho_old
    integer ( kind = 4 ) rlbl
    real ( kind = 8 ) x(n)
    real ( kind = 8 ) z(n)
    !
    !  Some local variables must be preserved between calls.
    !
    save iter
    save rho
    save rho_old
    save rlbl
    !
    !  Initialization.
    !  Ask the user to compute the initial residual.
    !
    if ( job == 1 ) then
       r(1:n) = b(1:n)
       job = 3
       rlbl = 2
       !
       !  Begin first conjugate gradient loop.
       !  Ask the user for a preconditioner solve.
       !
    else if ( rlbl == 2 ) then
       iter = 1
       job = 2
       rlbl = 3
       !
       !  Compute the direction.
       !  Ask the user to compute ALPHA.
       !  Save A*P to Q.
       !
    else if ( rlbl == 3 ) then
       rho = dot_product ( r, z )
       if ( 1 < iter ) then
          beta = rho / rho_old
          z(1:n) = z(1:n) + beta * p(1:n)
       end if
       p(1:n) = z(1:n)
       job = 1
       rlbl = 4
       !
       !  Compute current solution vector.
       !  Ask the user to check the stopping criterion.
       !
    else if ( rlbl == 4 ) then
       pdotq = dot_product ( p, q )
       alpha = rho / pdotq
       x(1:n) = x(1:n) + alpha * p(1:n)
       r(1:n) = r(1:n) - alpha * q(1:n)
       job = 4
       rlbl = 5
       !
       !  Begin the next step.
       !  Ask for a preconditioner solve.
       !
    else if ( rlbl == 5 ) then
       rho_old = rho
       iter = iter + 1
       job = 2
       rlbl = 3
    end if
    return
  end subroutine cg_rc
  subroutine solve_cg ( n, diag, nz_num, ia, ja, a, b, x )
    !
    !  SOLVE_CG solves a linear system using the conjugate gradient method.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 January 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of nodes.
    !
    !    Input, integer ( kind = 4 ) DIAG(N), contains for each index 1 <= I <= N, 
    !    the unique index J such that IA(J) = JA(J) = I.
    !
    !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries.
    !
    !    Input, integer ( kind = 4 ) IA(NZ_NUM), JA(NZ_NUM), the row and column
    !    indices of the nonzero entries.
    !
    !    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero entries of the matrix.
    !
    !    Input, real ( kind = 8 ) B(N), the right hand side.
    !
    !    Output, real ( kind = 8 ) X(N), the solution of the linear system.
    !
    implicit none
    integer ( kind = 4 ) n
    integer ( kind = 4 ) nz_num
    real ( kind = 8 ) a(nz_num)
    real ( kind = 8 ) aii
    real ( kind = 8 ) b(n)
    real ( kind = 8 ) bnrm2
    integer ( kind = 4 ) diag(n)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ia(nz_num)
    integer ( kind = 4 ) it
    integer ( kind = 4 ) it_max
    integer ( kind = 4 ) j
    integer ( kind = 4 ) ja(nz_num)
    integer ( kind = 4 ) job
    integer ( kind = 4 ) k
    real ( kind = 8 ) p(n)
    real ( kind = 8 ) q(n)
    real ( kind = 8 ) r(n)
    real ( kind = 8 ) rnrm2
    real ( kind = 8 ) tol
    real ( kind = 8 ) x(n)
    real ( kind = 8 ) z(n)
    it = 0
    it_max = 10000
    tol = 1.0D-30
    bnrm2 = sqrt ( sum ( b(1:n)**2 ) )
    x(1:n) = b(1:n) / a(diag(1:n))
!!$    write ( *, '(a)' ) ''
!!$    write ( *, '(a)' ) '  Step        Residual'
!!$    write ( *, '(a)' ) ''
    job = 1
    do
       call cg_rc ( n, b, x, r, z, p, q, job )
       !
       !  Compute q = A * p.
       !
       if ( job == 1 ) then
          q(1:n) = 0.0D+00
          do k = 1, nz_num
             i = ia(k)
             j = ja(k)
             q(i) = q(i) + a(k) * p(j)
          end do
          !
          !  Solve M * z = r.
          !
       else if ( job == 2 ) then
          z(1:n) = r(1:n) / a(diag(1:n))
          !
          !  Compute r = r - A * x.
          !
       else if ( job == 3 ) then
          do k = 1, nz_num
             i = ia(k)
             j = ja(k)
             r(i) = r(i) - a(k) * x(j)
          end do
          !
          !  Stopping test.
          !
       else if ( job == 4 ) then
          rnrm2 = sqrt ( sum ( r(1:n)**2 ) )
          if ( bnrm2 == 0.0D+00 ) then
             if ( rnrm2 <= tol ) then
                exit
             end if
          else
             if ( rnrm2 <= tol * bnrm2 ) then
                exit
             end if
          end if
          it = it + 1
!!$          write ( *, '(2x,i4,2x,g14.6)' ) it, rnrm2
          if ( it_max <= it ) then
!!$             write ( *, '(a)' ) ''
             write ( *, '(a)' ) '  Iteration limit exceeded.'
             write ( *, '(a)' ) '  Terminating early.'
             exit
          end if
       end if
       job = 2
    end do
!!$    write ( *, '(a)' ) ''
    write ( *, '(a,i6)' ) '  Number of iterations was ', it
    write ( *, '(a,g14.6)' ) '  Estimated error is ', rnrm2
    return
  end subroutine solve_cg
  function CGOMP(M,b) result(x)
    implicit none
    class(sparse), intent(in) :: M
    real(rkind), intent(in) :: b(M%n)
    integer(ikind) :: i
    integer(ikind) :: j
    integer(ikind) :: it
    integer(ikind) :: it_max
    real(rkind) :: p(M%n)
    real(rkind) :: y(M%n)
    real(rkind) :: q(M%n)
    real(rkind) :: r(M%n)
    real(rkind) :: x(M%n)
    real(rkind) :: z(M%n)
    real(rkind) :: diag(M%n)
    real(rkind) :: bb
    real(rkind) :: rr
    real(rkind) :: tol
    real(rkind) :: alpha
    real(rkind) :: beta
    real(rkind) :: pq
    real(rkind) :: rho
    real(rkind) :: rho_old
     x = 0.d0
    do i = 1, size(M%AI)-1
       do j = M%AI(i), M%AI(i+1)-1
          if ( i == M%AJ(j) ) then
             diag(i) = M%A(j)
          end if
       end do
    end do
    it_max = 10000
    tol = 1.0D-15
    it = 0
    call vecdiv(M%n,b,diag,x)
    bb = sqrt (dot_product(b,b))
    do i = 1, M%n
       r(i) = b(i)                                
    end do
    y = M * x
    call vecsum(M%n,1.d0,r,-y,r)
    do
       it = it + 1
       call vecdiv(M%n,r,diag,z)
       rho = dot_product(r,z)
       if (it>1)then
          beta = rho / rho_old
          call vecsum(M%n,beta,p,z,z)
       end if
       do i = 1, M%n
          p(i) = z(i)
       end do
       call assign(M%n,q,0.d0)
       y = M * p
       call vecsum(M%n,1.d0,y,q,q)
       alpha = rho / dot_product(p,q)
       call vecsum(M%n,alpha,p,x,x)
       call vecsum(M%n,-alpha,q,r,r)
       rr = dot_product(r,r)
       if ( bb == 0.d0 ) then
          if ( rr <= tol ) then
             exit
          end if
       else
          if ( rr <= tol * bb ) then
             exit
          end if
       end if
       if ( it_max <= it ) then
          exit
       end if
       rho_old = rho
    end do
    write ( *,'(A,I5)' ) '  Number of iterations: ', it
    write ( *,'(A,E14.7)' ) '  Estimated error: ', rr
  end function CGOMP
end module SparseKit
