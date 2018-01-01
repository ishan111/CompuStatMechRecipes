module random
contains
  integer function randz(zMin,zMax)
    ! creates random Integers
    integer :: zMin, zMax, z
    real :: n,randFract
    randFract = rand()
    n = zMax - zMin
    randz = zMin + nint(randFract*n)
  end function randz
  real function rand_norm()
    real :: uRand1, uRand2 , nRand1, nRand2
    real,parameter :: PI=4.D0*DATAN(1.D0)
    uRand1 = rand()
    uRand2 = rand()
    nRand1 = ((-2*log(uRand1))**.5)*cos(2*PI*uRand2)
    nRand2 = ((-2*log(uRand1))**.5)*sin(2*PI*uRand2)
    rand_norm = nRand1
  end function rand_norm
end module random
