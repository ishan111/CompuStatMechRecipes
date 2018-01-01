program mc_simulation
  use new_config
  use sims_props
  use system_props
  use Initialize
  use random
  use sample
  use tmmc
  implicit none
  integer :: mSel,pSelect
  call Initialize()
  do trialNo = 1,noTrials
    mSel = select_move()
    call copyToAttemptBuffer()
    select case(mSel)
    case(mDisp)

      pSelect = randz(1,npNo)
      call disp_part(pSelect)
      call acc_prob_disp()

      call sample_disp()
    case(mRem)

      pSelect = randz(1,npNo)
      call disp_part(pSelect)
      call acc_prob_rem()
      call sample_rem()
    case(mIns)

      pSelect = NnpNo + 1
      call disp_part_ins(pSelect)
      call acc_prob()
      call sample_ins()
    end select
  end do
end program mc_simulation
