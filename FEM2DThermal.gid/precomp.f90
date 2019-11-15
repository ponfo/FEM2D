program precomp
  use tools
  implicit none
  integer(ikind), parameter :: functionData=1 , function=2, projectData=3
  integer(ikind)            :: i, iSource, nmat, numberSource , nPoints ,iPoint, date_time(8)
  integer(ikind)            :: nSourceOS, nSourceOL, nSourceOP
  character(100)            :: sourceFun, projectName, path
  character(40)             :: aux
  open(projectData, file = 'projectData.dat')
  read(projectData,'(*(A))') projectName
  read(projectData,'(*(A))') path
  close(projectData)
  print'(/,A)','::::::::::::::: Starting FEM2D :::::::::::::::'
  call date_and_time(VALUES=date_time)
  print'(A,X,I0,A,I0,A,I0)','Date:',date_time(3),"/",date_time(2),"/", date_time(1)
  print'(A,X,I0,A,I0,A,I0)','Hour:',date_time(5),":",date_time(6),":",date_time(7)
  print'(/,A,/)','Writing functions...'
  open(functionData, file = trim(projectName)//'.dat')
  do i = 1, 10
     read(functionData,*)
  end do
  read(functionData,*)  aux, nPoints
  read(functionData,*)  
  read(functionData,*)  
  read(functionData,*)  
  read(functionData,*)
  read(functionData,*)  aux, nMat
  do i = 1, 7
     read(functionData,*)
  end do
  read(functionData,*)  aux, nSourceOP
  read(functionData,*)
  read(functionData,*)  aux, nSourceOL
  read(functionData,*)
  read(functionData,*)  aux, nSourceOS
  
  do iPoint = 1, 20+nPoints+nMat
     read(functionData,*)
  end do

  open(function, file = trim(path)//'/FEM2D/Source/functionOnPoints.f90')
  write(function,'(A)')       'module FunctionOnPointsMOD'
  write(function,'(A)')       '  use tools'
  write(function,'(A)')       '  implicit none'
  write(function,'(A)')       'contains'
  write(function,'(A)')       '  real(rkind) function funcOnPoints(this, i, x, y)' 
  write(function,'(A)')       '  implicit none'
  write(function,'(A)')       '  integer(ikind), intent(in) :: i'
  write(function,'(A)')       '  real(rkind), intent(in) :: x'
  write(function,'(A)')       '  real(rkind), intent(in) :: y'
  write(function,'(A,I0,A)')  '  real(rkind) :: vecFunSource(',nSourceOP,')'
  
  do iSource = 1, nSourceOP
     read(functionData,'(I10,X,A)') numberSource, sourceFun
     write(function,'(A,I0,A)')'  vecFunSource(',numberSource,')='//trim(sourceFun)
     print'(A,I0,A)','pointFunc(',numberSource,')='//trim(sourceFun)
  end do
  
  write(function,'(A)')       '  funcOnPoints = vecfunSource(i)'
  write(function,'(A)')       '  end function funcOnPoints'
  write(function,'(A)')       'end module FunctionOnPointsMOD'
  close(function)

  do iPoint = 1, 7
     read(functionData,*)
  end do

  open(function, file = trim(path)//'/FEM2D/Source/functionOnLines.f90')
  write(function,'(A)')       'module FunctionOnLinesMOD'
  write(function,'(A)')       '  use tools'
  write(function,'(A)')       '  implicit none'
  write(function,'(A)')       'contains'
  write(function,'(A)')       '  real(rkind) function funcOnLines(this, i, x, y)' 
  write(function,'(A)')       '  implicit none'
  write(function,'(A)')       '  integer(ikind), intent(in) :: i'
  write(function,'(A)')       '  real(rkind), intent(in) :: x'
  write(function,'(A)')       '  real(rkind), intent(in) :: y'
  write(function,'(A,I0,A)')  '  real(rkind) :: vecFunSource(',nSourceOL,')'
  
  do iSource = 1, nSourceOL
     read(functionData,'(I10,X,A)') numberSource, sourceFun
     write(function,'(A,I0,A)')'  vecFunSource(',numberSource,')='//trim(sourceFun)
     print'(A,I0,A)','lineFunc(',numberSource,')='//trim(sourceFun)
  end do
  
  write(function,'(A)')       '  funcOnLines = vecfunSource(i)'
  write(function,'(A)')       '  end function funcOnLines'
  write(function,'(A)')       'end module FunctionOnLinesMOD'
  close(function)
  
  do iPoint = 1, 7
     read(functionData,*)
  end do

  open(function, file = trim(path)//'/FEM2D/Source/functionOnSurfaces.f90')
  write(function,'(A)')       'module FunctionOnSurfacesMOD'
  write(function,'(A)')       '  use tools'
  write(function,'(A)')       '  implicit none'
  write(function,'(A)')       'contains'
  write(function,'(A)')       '  real(rkind) function funcOnSurfaces(this, i, x, y)' 
  write(function,'(A)')       '  implicit none'
  write(function,'(A)')       '  integer(ikind), intent(in) :: i'
  write(function,'(A)')       '  real(rkind), intent(in) :: x'
  write(function,'(A)')       '  real(rkind), intent(in) :: y'
  write(function,'(A,I0,A)')  '  real(rkind) :: vecFunSource(',nSourceOS,')'

  do iSource = 1, nSourceOS
     read(functionData,'(I10,X,A)') numberSource, sourceFun
     write(function,'(A,I0,A)')'  vecFunSource(',numberSource,')='//trim(sourceFun)
     print'(A,I0,A)','surfaceFunc(',numberSource,')='//trim(sourceFun)
  end do
  
  write(function,'(A)')       '  funcOnSurfaces = vecfunSource(i)'
  write(function,'(A)')       '  end function funcOnSurfaces'
  write(function,'(A)')       'end module FunctionOnSurfacesMOD'
  close(function)

  close(functionData)
  print'(A)', 'Compiling functions...'
  call execute_command_line('cd;cd '//trim(path)//'/ ;make')
  print'(A)', 'Finish compiling functions'
end program precomp
