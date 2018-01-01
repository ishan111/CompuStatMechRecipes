program HelloRand
  use random
  implicit none
  integer :: kill,a
  print *,"hello",rand()
  a=randz(5,9)
  do kill=1,100
    print *, rand(),"   ",randz(5,9),"     ",rand_norm()
  end do
end program HelloRand
