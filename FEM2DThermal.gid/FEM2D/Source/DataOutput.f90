module DataOutputMOD
  use tools
  implicit none
  private
  public :: printResults, finishProgram
  interface printResults
     procedure :: printResultsVec1
     procedure :: printResultsVec2
  end interface printResults
  interface finishProgram
     procedure :: finishProgram
  end interface finishProgram
  integer(ikind), parameter    :: projectData = 1
  integer(ikind), parameter    :: results = 3
  integer(ikind), dimension(8) :: date_time
  character(100)               :: projectName, path
contains
  subroutine init()
    implicit none
    open(projectData, file = 'projectData.dat')
    read(projectData, '(A)') projectName
    read(projectData, '(A)') path
    close(projectData)
  end subroutine init
  subroutine printResultsVec1(resultName, step, graphType, locationName, resultNumber&
       , component1)
    implicit none
    integer(ikind)                 :: iPoint
    integer(ikind), intent(in)     :: step, resultNumber
    real(rkind), intent(in), dimension(resultNumber) :: component1
    character(*), intent(in) :: resultName, graphType, locationName
    call init()
    open(results, file = trim(projectName)//'.flavia.res')
    write(results,'(A)')      'GiD Post Result File 1.0'
    write(results,*)   'Result "',trim(resultName),'" "',trim(projectName)&
         ,'" ',step,' ',trim(graphType),' ',trim(locationName)
    write(results,*)   'Values'
    Do iPoint = 1, resultNumber
       Write(results,*) iPoint, component1(iPoint)
    End Do
    write(results,*)   'End Values'
  end subroutine printResultsVec1
  subroutine printResultsVec2(resultName, type, step, graphType, locationName, gaussPoints &
       , resultNumber, elemID, component1, component2)
    implicit none
    character(*), intent(in) :: resultName
    character(*), intent(in) :: type
    integer(ikind), intent(in) :: step
    character(*), intent(in) :: graphType
    character(*), intent(in) :: locationName
    real(rkind), dimension(:,:), intent(in) :: gaussPoints
    integer(ikind), intent(in) :: resultNumber
    integer(ikind), dimension(resultNumber), intent(in) :: elemID
    real(rkind), dimension(:), intent(in) :: component1
    real(rkind), dimension(:), intent(in) :: component2
    real(rkind) :: prom(2)
    integer(ikind) :: i, j, k, count, numberGP
    if(resultNumber == 0) return
    call init()
    write(results,'(/,3A)') 'GaussPoints "Points'//trim(resultName), '" ElemType ', trim(type)
    write(results,'(A,I0)') 'Number of GaussPoints: ', size(gaussPoints,1)
    write(results,'(A)') 'Natural Coordinates: Given'
    do i = 1, size(gaussPoints,1)
       write(results,'(F26.16,2X,F26.16)') gaussPoints(i,1), gaussPoints(i,2)
    end do
    write(results,'(A)') 'End gausspoints'
    write(results,'(5A,I0,6A)') 'Result "', trim(resultName), '" "', trim(projectName), '" ', step &
         , ' ', trim(graphType), ' ', trim(locationName), ' "Points'//trim(resultName) , '"'
    write(results,'(A)') 'Values'
    count = 0
    do i = 1, resultNumber
       count = count + 1
       write(results,'(I0,2X,F26.16,2X,F26.16)') elemID(i), component1(count), component2(count)
       do j = 2, size(gaussPoints,1)
          count = count + 1
          write(results,'(6X,F26.16,2X,F26.16)') component1(count), component2(count)
       end do
    end do
    write(results,'(A)') 'End Values'
  end subroutine printResultsVec2

  subroutine finishProgram()
    implicit none
    close(results)
    !call free()
    call date_and_time(VALUES=date_time)
    print'(A)','::::::::::::::: Finish FEM2D :::::::::::::::'
    print'(A,I0,A,I0,A,I0)', 'Date: ', date_time(3), "/", date_time(2), "/", date_time(1)
    print'(A,I0,A,I0,A,I0)', 'Hour: ', date_time(5), ":", date_time(6), ":", date_time(7)
  end subroutine finishProgram
!!$  subroutine free()
!!$    integer(ikind), parameter :: function=1, projectData=2
!!$    character(50)             :: projectName, path
!!$    open(projectData, file = 'projectData.dat')
!!$    read(projectData,'(A50)') projectName
!!$    read(projectData,'(A50)') path
!!$    close(projectData)
!!$    open(function, file = trim(path)//'/FEM2D/Source/Source.f90')
!!$    write(function,'(A)')       'module SourceMOD'
!!$    write(function,'(A)')       '  use tools'
!!$    write(function,'(A)')       '  implicit none'
!!$    write(function,'(A)')       '  private'
!!$    write(function,'(A)')       '  public :: sourceTYPE'
!!$    write(function,'(A)')       '  type sourceTYPE'
!!$    write(function,'(A)')       '   contains'
!!$    write(function,'(A)')       '     procedure :: func'
!!$    write(function,'(A)')       '  end type sourceTYPE'
!!$    write(function,'(A)')       'contains'
!!$    write(function,'(A)')       '  function func( this, i, x, y)' 
!!$    write(function,'(A)')       '    implicit none'
!!$    write(function,'(A)')       '    class(sourceTYPE) :: this'
!!$    write(function,'(A)')       '    integer(ikind), intent(in) :: i'
!!$    write(function,'(A)')       '    real(rkind), intent(in) :: x, y'
!!$    write(function,'(A)')       '    real(rkind) :: func'
!!$    write(function,'(A)')       '  end function func'
!!$    write(function,'(A)')       'end module SourceMOD'
!!$    close(function)
!!$  end subroutine free
end module DataOutputMOD
