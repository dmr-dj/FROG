

program test

  use main_lib_FROG, only: INITIALIZE_FROG, INITIALIZE_FROGVARS, STEPFWD_FROG

  implicit none

  logical :: well_done = .FALSE.
  integer :: s, steps=5


  well_done = INITIALIZE_FROG()
  if (well_done) then
    WRITE(*,*) "FROG GRID INITIALIZATION DONE"
  endif

  well_done = INITIALIZE_FROGVARS()
  if (well_done) then
    WRITE(*,*) "FROG VARS INITIALIZATION DONE"
  endif

  do s=1, steps

    well_done = STEPFWD_FROG()

  enddo

  if (well_done) then
    WRITE(*,*) "FROG INTEGRATION COMPLETED"
  endif

 end program test
