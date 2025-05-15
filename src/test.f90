

program test

  use main_lib_VAMPER, only: INITIALIZE_VAMP, STEPFWD_VAMP

  implicit none

  logical :: well_done = .FALSE.


  well_done = INITIALIZE_VAMP()
  if (well_done) then
    WRITE(*,*) "VAMPER INITIALIZATION COMPLETE"
  endif

  well_done = STEPFWD_VAMP()
  if (well_done) then
    WRITE(*,*) "VAMPER INTEGRATION COMPLETED"
  endif

 end program test
