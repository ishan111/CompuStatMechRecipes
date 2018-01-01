program HelloRand
  use random
  use new_config
  implicit none
  integer :: kill,a
  integer :: npNo = 3
  integer,allocatable,dimension(:,:) :: npConf

  print *,"hello",rand()
  a=randz(5,9)
  do kill=1,100
    print *, rand(),"   ",randz(5,9),"     ",rand_norm()
  end do
end program HelloRand
