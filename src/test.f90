

program test

  use main_lib_FROG, only: INITIALIZE_FROG, STEPFWD_FROG

  implicit none

  logical :: well_done = .FALSE.
  integer :: s, steps=5


  well_done = INITIALIZE_FROG()
  if (well_done) then
    WRITE(*,*) "FROG INITIALIZATION COMPLETE"
  endif

  do s=1, steps

    well_done = STEPFWD_FROG()

  enddo

  if (well_done) then
    WRITE(*,*) "FROG INTEGRATION COMPLETED"
  endif

 end program test
