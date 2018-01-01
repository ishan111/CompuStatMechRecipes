module new_config
  use system_props
  use npConf
  use random
  use sims_props
contains
  !to remove particle
  subroutine remove_part( pSelect )
    !integer , intent(inout) ::  npNo,pSelect, dimes
    integer , intent(in) :: pSelect
    !real, dimension(:,:), allocatable ,intent(inout) ::  npConf
    real, dimension(:,:) , allocatable ::  temp
    integer :: i,j
    allocate(temp(npNo-1,dimes))

    do i=1,npNo
      do j=1:dimes
        if (i<pSelect) then
          temp(i,j) =  npConf(i,j)
        elseif (i>pSelect) then
          temp(i-1,j)=npConf(i,j)
        endif
      end do
    end do
    npConf = intent
    npNo = npNo - 1
  end subroutine remove_part

  subroutine insert_part( pSelect)
    integer , intent(in) :: pSelect
    real, dimension(:,:) , allocatable ::  temp
    integer :: i,j
    allocate(temp(npNo+1,dimes))

    do i=1,npNo+1
      do j=1:dimes
        if (i<pSelect) then
          temp(i,j) =  npConf(i,j)
        else
          temp(i-1,j)=rand()*nSysLen(j)
        endif
      end do
    end do
    npConf = intent
    npNo = npNo + 1
  end subroutine insert_part

  subroutine disp_part(pSelect)
    integer , intent(in) :: pSelect
    real,allocatable,dimension(:) :: disp_vector
    real : disp_len = 0
    allocate(disp_vector(dimes))
    integer :: j
    do j=1,dimes
      disp_vector(j) = rand_norm()
      disp_len = disp_len + disp_vector(j)**2
    end do
    disp_len =  max_disp_len * (disp_len**0.5)
    do j=1,dimes
      disp_vector(j) = disp_vector(j) * disp_len
      pConf(pSelect,j) = pConf(pSelect,j) + disp_vector(j)
      pConf(pSelect,j) = modulo(pConf(pSelect,j),nSysLen(j))
    end do


end subroutine disp_part
end module new_config
