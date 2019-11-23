program precomp
  use tools
  implicit none
  integer(ikind), parameter :: functionData=1 , function=2, projectData=3
  integer(ikind)            :: i, iLoad, nmat, numberLoad , nPoints ,iPoint, date_time(8)
  integer(ikind)            :: nLoadOS, nLoadOL, nLoadOP
  character(100)            :: loadFunX, loadFunY, projectName, path
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
  do i = 1, 4
     read(functionData,*)
  end do
  read(functionData,*)  aux, nLoadOP
  read(functionData,*)
  read(functionData,*)  aux, nLoadOL
  read(functionData,*)
  read(functionData,*)  aux, nLoadOS
  
  do iPoint = 1, 20+nPoints+nMat
     read(functionData,*)
  end do

  open(function, file = trim(path)//'/FEM2D/Source/Load/functionOnPoints.f90')
  write(function,'(A)')       'module FunctionOnPointsMOD'
  write(function,'(A)')       '  use tools'
  write(function,'(A)')       '  implicit none'
  write(function,'(A)')       ' contains'
  write(function,'(A)')       '  function funcOnPoints(i, x, y)' 
  write(function,'(A)')       '  implicit none'
  write(function,'(A)')       '  integer(ikind), intent(in) :: i'
  write(function,'(A)')       '  real(rkind), intent(in) :: x'
  write(function,'(A)')       '  real(rkind), intent(in) :: y'
  write(function,'(A,I0,A)')  '  real(rkind) :: vecFunX(',nLoadOP,')'
  write(function,'(A,I0,A)')  '  real(rkind) :: vecFunY(',nLoadOP,')'
  write(function,'(A)')       '  real(rkind), dimension(2) :: funcOnPoints'
  
  do iLoad = 1, nLoadOP
     read(functionData,'(I5,X,A10,X,A10)') numberLoad, loadFunX, loadFunY
     write(function,'(A,I0,A)')'  vecFunX(',numberLoad,')='//trim(loadFunX)
     write(function,'(A,I0,A)')'  vecFunX(',numberLoad,')='//trim(loadFunY)
     write(*,'(A,I0,A)')'  vecFunX(',numberLoad,')='//trim(loadFunX)
     write(*,'(A,I0,A)')'  vecFunY(',numberLoad,')='//trim(loadFunY)
  end do
  
  write(function,'(A)')       '  funcOnPoints(1) = vecFunX(i)'
  write(function,'(A)')       '  funcOnPoints(2) = vecFunY(i)'
  write(function,'(A)')       '  end function funcOnPoints'
  write(function,'(A)')       'end module FunctionOnPointsMOD'
  close(function)

  do iPoint = 1, 7
     read(functionData,*)
  end do

  open(function, file = trim(path)//'/FEM2D/Source/Load/functionOnLines.f90')
  write(function,'(A)')       'module FunctionOnLinesMOD'
  write(function,'(A)')       '  use tools'
  write(function,'(A)')       '  implicit none'
  write(function,'(A)')       ' contains'
  write(function,'(A)')       '  function funcOnLines(i, x, y)' 
  write(function,'(A)')       '  implicit none'
  write(function,'(A)')       '  integer(ikind), intent(in) :: i'
  write(function,'(A)')       '  real(rkind), intent(in) :: x'
  write(function,'(A)')       '  real(rkind), intent(in) :: y'
  write(function,'(A,I0,A)')  '  real(rkind) :: vecFunX(',nLoadOL,')'
  write(function,'(A,I0,A)')  '  real(rkind) :: vecFunY(',nLoadOL,')'
  write(function,'(A)')       '  real(rkind), dimension(2) :: funcOnLines'
  
  do iLoad = 1, nLoadOL
     read(functionData,'(I5,X,A10,X,A10)') numberLoad, loadFunX, loadFunY
     write(function,'(A,I0,A)')'  vecFunX(',numberLoad,')='//trim(loadFunX)
     write(function,'(A,I0,A)')'  vecFunX(',numberLoad,')='//trim(loadFunY)
     write(*,'(A,I0,A)')'  vecFunX(',numberLoad,')='//trim(loadFunX)
     write(*,'(A,I0,A)')'  vecFunY(',numberLoad,')='//trim(loadFunY)
  end do
  
  write(function,'(A)')       '  funcOnLines(1) = vecFunX(i)'
  write(function,'(A)')       '  funcOnLines(2) = vecFunY(i)'
  write(function,'(A)')       '  end function funcOnLines'
  write(function,'(A)')       'end module FunctionOnLinesMOD'
  close(function)
  
  do iPoint = 1, 7
     read(functionData,*)
  end do

  open(function, file = trim(path)//'/FEM2D/Source/Load/functionOnSurfaces.f90')
  write(function,'(A)')       'module FunctionOnSurfacesMOD'
  write(function,'(A)')       '  use tools'
  write(function,'(A)')       '  implicit none'
  write(function,'(A)')       ' contains'
  write(function,'(A)')       '  function funcOnSurfaces(i, x, y)' 
  write(function,'(A)')       '  implicit none'
  write(function,'(A)')       '  integer(ikind), intent(in) :: i'
  write(function,'(A)')       '  real(rkind), intent(in) :: x'
  write(function,'(A)')       '  real(rkind), intent(in) :: y'
  write(function,'(A,I0,A)')  '  real(rkind) :: vecFunX(',nLoadOS,')'
  write(function,'(A,I0,A)')  '  real(rkind) :: vecFunY(',nLoadOS,')'
  write(function,'(A)')       '  real(rkind), dimension(2) :: funcOnSurfaces'
  
  do iLoad = 1, nLoadOS
     read(functionData,'(I5,X,A10,X,A10)') numberLoad, loadFunX, loadFunY
     write(function,'(A,I0,A)')'  vecFunX(',numberLoad,')='//trim(loadFunX)
     write(function,'(A,I0,A)')'  vecFunY(',numberLoad,')='//trim(loadFunY)
     write(*,'(A,I0,A)')'  vecFunX(',numberLoad,')='//trim(loadFunX)
     write(*,'(A,I0,A)')'  vecFunX(',numberLoad,')='//trim(loadFunY)
  end do
  
  write(function,'(A)')       '  funcOnSurfaces(1) = vecFunX(i)'
  write(function,'(A)')       '  funcOnSurfaces(2) = vecFunY(i)'
  write(function,'(A)')       '  end function funcOnSurfaces'
  write(function,'(A)')       'end module FunctionOnSurfacesMOD'
  close(function)

  close(functionData)
  print'(A)', 'Compiling functions...'
  call execute_command_line('cd;cd '//trim(path)//'/ ;make')
  print'(A)', 'Finish compiling functions'
end program precomp
