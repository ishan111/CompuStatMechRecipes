program HelloRand
  implicit none
  integer :: kill,a
  print *,"hello",rand()
  a=randz(5,9)
  do kill=1,100
    print *, rand(),"   ", randz(5,9)
  enddo
contains
  integer function randz(zMin,zMax)
    integer :: zMin, zMax, z
    real :: n,randFract
    randFract = rand()
    n=zMax-zMin
    randz=zMin + nint(randFract*n)
  end function randz

end program HelloRand
