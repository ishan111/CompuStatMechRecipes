module sample
  use results
  use copyConfig
  use npConf
  use pConf
  use tmmc
contains
  subroutine sample_disp(accProb)

    real :: accProb
    if () then
      call tm_update()
    endif
    if () then
      call tm_bias()
    endif
    attRand = log(rand())
    if (attRand<=accProb) then
      call updateConfig()
      call results()

    else
      call results()

    end if
  end subroutine sample_disp
end module sample
