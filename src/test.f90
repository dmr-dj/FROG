

program test

  use main_lib_FROG, only: INITIALIZE_VAMP, STEPFWD_VAMP

  implicit none

  logical :: well_done = .FALSE.
  integer :: s, steps=10


  well_done = INITIALIZE_VAMP()
  if (well_done) then
    WRITE(*,*) "FROG INITIALIZATION COMPLETE"
  endif

  do s=1, steps

    well_done = STEPFWD_VAMP()

  enddo

  if (well_done) then
    WRITE(*,*) "FROG INTEGRATION COMPLETED"
  endif

 end program test
