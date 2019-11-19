module GeometryPtrMOD
  use tools
  use GeometryMOD
  implicit none
  private
  public :: GeometryPtrTYPE
  type GeometryPtrTYPE
     class(GeometryTYPE), pointer :: ptr
  end type GeometryPtrTYPE
end module GeometryPtrMOD
