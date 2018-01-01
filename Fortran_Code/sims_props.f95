module sims_props
  integer :: noTrials, trialNo
  integer,allocatable, dimension(:) :: nAtt,nAcc
  !move thing
  integer :: mDisp = 1
  integer :: mRem = 2
  integer :: mIns = 3
  integer :: mTot = 4

  real :: max_disp_len
  real,allocatable,dimension(:) :: moveProbs
end module sims_props
