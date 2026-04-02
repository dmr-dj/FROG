

program test

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

  use main_lib_FROG, only: INITIALIZE_FROG, INITIALIZE_FROGVARS, STEPFWD_FROG, WRITE_FROGRESTART

  implicit none

  logical :: well_done = .FALSE.
  integer :: s, steps=5

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

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

  well_done = WRITE_FROGRESTART()


 end program test
