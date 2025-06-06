!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2024 Didier M. Roche (a.k.a. dmr)

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#include "constant.h"

    MODULE timer_mod

     USE parameter_mod, only: YearType

     IMPLICIT NONE

     PUBLIC

     TYPE :: cell_time
       INTEGER :: current_step
       LOGICAL :: end_year
       LOGICAL :: end_month
       LOGICAL :: is_daily
     END TYPE cell_time

     integer, parameter :: Year360 = 360, Year365 = 365, Month_30days=30, nb_months = 12
     integer, dimension(nb_months), parameter :: list_endmonth365 = (/31,59,90,120,151,181,212,243,273,304,334,365/)

     CONTAINS

     function update_time_cell(cell_time_var) result(success)

       type(cell_time), intent(inout) :: cell_time_var

       integer :: remainder_year

       logical :: success

       cell_time_var%current_step = cell_time_var%current_step + 1

       remainder_year = MOD(cell_time_var%current_step,YearType)

       if (remainder_year.EQ.0) then
         cell_time_var%end_year = .true.
       else
         cell_time_var%end_year = .false.
       endif

       if (cell_time_var%is_daily) then
         if (YearType.EQ.Year360) then
           if (MOD(remainder_year,Month_30days).EQ.0) then
             cell_time_var%end_month = .true.
           else
             cell_time_var%end_month = .false.
           endif
         elseif (YearType.EQ.Year365) then
           if (FINDLOC(list_endmonth365,remainder_year,1).NE.0) then
             cell_time_var%end_month = .true.
           else
             cell_time_var%end_month = .false.
           endif
         else
           WRITE(*,*) "Unknown YearType", Yeartype
           STOP
         endif
       else
         cell_time_var%end_month = .true.
       endif

!~       ! CHECK
!~         WRITE(*,*) "Cell_time updated::"
!~         WRITE(*,*) cell_time_var%current_step, "step"
!~         WRITE(*,*) cell_time_var%end_year    , "end_year"
!~         WRITE(*,*) cell_time_var%end_month   , "end_month"
!~         WRITE(*,*) cell_time_var%is_daily    , "is_daily"

!~       ! CHECK


     end function update_time_cell


     function init_time_cell(step_init,log_year,log_month,log_daily) result(initialized_time_cell)

       integer, intent(in) :: step_init
       logical, intent(in) :: log_year,log_month,log_daily

       type(cell_time) :: initialized_time_cell

       initialized_time_cell%current_step = step_init
       initialized_time_cell%end_year     = log_year
       initialized_time_cell%end_month    = log_month
       initialized_time_cell%is_daily     = log_daily

     end function init_time_cell

    END MODULE timer_mod

!dmr The End of All Things (op. cit.)
