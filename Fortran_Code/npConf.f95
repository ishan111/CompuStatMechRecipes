module npConf
  use  system_props
  implicit none

  real,dimension(:,:),allocatable   :: npConf
  integer :: npNo
  real,allocatable,dimension(:) :: nSysLen
end
