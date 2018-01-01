module copyConfig
  use pConf
  use npConf
  use system_props
contains
  subroutine copyToAttemptBuffer
    integer :: i,j
    npNo = pNo
    if (allocated(npNo)) then
      deallocate(npNo)
    endif
    allocate(npConf(npNo,dimes))
    do i=1,pNo
      do j=1,dimes
        npConf(i,j)=pConf(i,j)
      end do
    end do

  end subroutine copyToAttemptBuffer

  subroutine updateConfig
    integer :: i,j
    pNo = npNo
    if (allocated(pNo)) then
      deallocate(pNo)
    endif
    allocate(pConf(pNo,dimes))
    do i=1,pNo
      do j=1,dimes
        pConf(i,j)=npConf(i,j)
      end do
    end do
  end subroutine updateConfig
end module copyConfig
