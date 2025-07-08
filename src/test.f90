

program test

  use main_lib_VAMPER, only: INITIALIZE_VAMP, STEPFWD_VAMP

  implicit none

  logical :: well_done = .FALSE.
  integer :: s, steps=10


  well_done = INITIALIZE_VAMP()
  if (well_done) then
    WRITE(*,*) "VAMPER INITIALIZATION COMPLETE"
  endif

  do s=1, steps

    well_done = STEPFWD_VAMP()

  enddo

  if (well_done) then
    WRITE(*,*) "VAMPER INTEGRATION COMPLETED"
  endif

 end program test
