

program test

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

  use main_lib_FROG, only: INITIALIZE_FROG, INITIALIZE_FROGVARS, STEPFWD_FROG, WRITE_FROGRESTART
  use parameter_mod, only: TotTime   !dmr number of simulated years (from /Param/ namelist)

  implicit none

  logical :: well_done = .FALSE.
  integer :: s, steps
!dmr&clo --- The offline driver advances one coupling step (= one year,
!dmr&clo     nb_coupling_steps = YearType) per STEPFWD_FROG call, so the number
!dmr&clo     of calls is exactly TotTime (years). Previously steps was hard-coded
!dmr&clo     to 5, which silently ignored TotTime from the namelist: a run
!dmr&clo     configured for e.g. 17 years still only did 5. TotTime is available
!dmr&clo     here because INITIALIZE_FROG() below runs read_namelist first.
!dmr&clo     (In the coupled build the host drives the loop; this file is the
!dmr&clo     offline driver only.)

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

  steps = TotTime   !dmr [offline] run length in years, from the namelist
  WRITE(*,*) "FROG offline run: ", steps, " year(s)"

  do s=1, steps

    well_done = STEPFWD_FROG()

  enddo

  if (well_done) then
    WRITE(*,*) "FROG INTEGRATION COMPLETED"
  endif

  well_done = WRITE_FROGRESTART(-1)


 end program test
