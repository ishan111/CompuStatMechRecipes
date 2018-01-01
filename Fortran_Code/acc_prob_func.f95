module acc_prob_func
  use system_props

contains
  subroutine acc_prob_disp(pSel,dist,ndist,accProb)
    integer,intent(in) :: pSel
    real,intent(in) ,dimension(:) :: dist,ndist
    real,intent(out) :: delEnergy,accProb
    delEnergy = del_energy(ndist)-del_energy(dist)
    accProb = -delEnergy/K*T
  end subroutine acc_prob_disp

  subroutine acc_prob_ins(pSel,ndist,accProb)
    integer,intent(in) :: pSel
    real,intent(in) ,dimension(:) :: ndist
    real,intent(out) :: delEnergy,accProb
    delEnergy = del_energy(ndist)
    accProb = -delEnergy/K*T
  end subroutine acc_prob_ins

  subroutine acc_prob_rem(pSel,dist,accProb)
    integer,intent(in) :: pSel
    real,intent(in) ,dimension(:) :: dist
    real,intent(out) :: delEnergy,accProb
    delEnergy = -del_energy(dist)
    accProb = -delEnergy/K*T
  end subroutine acc_prob_rem

  real function del_energy(dist,pNo)
    integer,intent(in) :: pSel
    real,intent(in) ,dimension(:) :: dist
    real ::delEnergy = 0
    do i=1,pNo
      if (i/=pSel) then
        delEnergy = delEnergy + part_energy(dist(i))
      end if
    enddo
    del_energy = delEnergy
  end function del_energy

  real function part_energy(dist)
    real :: dist
    real :: partEnergy
    partEnergy = 4*eps((sig/dist)**12 - (sig/dist)**6)
    part_energy = partEnergy
  end function part_energy



  subroutine DistCalc(pConf,pNo, pSel,dist, sysLen, dimes)
    real,dimension(:,:),intent(in) :: pConf
    integer,intent(in) :: pNo
    integer,intent(in) :: pSel
    real,intent(out) ,dimension(:) :: dist
    real,intent(in) ,dimension(:) :: sysLen
    integer, intent(in) :: dimes
    integer :: i,j

    real,intent(out) ,dimension(:,dimes) :: dist


    do i=1,pNo
      dist(i)=0
      do j=1,dimes
        distx(i,j) = pConf(i,j)-pConf(pSel,j)
        distx(i,j)= min(abs(distx(i,j)),abs(sysLen(j)-distx(i,j)))
        dist(i)=dist(i)+distx(i,j)**2
      end do
      dist(i) = dist(i)**0.5
    endDo

  end subroutine DistCalc
end module acc_prob_func
