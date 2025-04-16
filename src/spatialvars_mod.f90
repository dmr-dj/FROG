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


    MODULE spatialvars_mod

    IMPLICIT NONE

    PRIVATE


            ! SPATIAL GLOBAL VARIABLES

     real, dimension(:,:),allocatable, PUBLIC  ::    Temp      & !dmr [SPAT_VAR], soil temperature over the vertical // prognostic
                                                    ,Kp        & !dmr [CNTST]     heat conductivity constant over the depth, current value is 2
                                                    ,n         & !dmr [SPAT_VAR], porosity on the vertical
                                                    ,Cp        & !dmr [SPAT_VAR]  specific heat capacity
                                                    ,pori      & !dmr [???  TBC]
                                                    ,porf        !dmr [???  TBC]


     PUBLIC:: allocate_spatialvars

     CONTAINS


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Allocation of two dimensional variables (VERTCL, SPAT_VAR)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     SUBROUTINE allocate_spatialvars ! VERTCL, SPAT_VAR

       use parameter_mod, only: gridNoMax, z_num

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       allocate(Temp(1:z_num,1:gridNoMax)) !dmr SPAT_VAR
       allocate(Kp(1:z_num-1,1:gridNoMax)) !dmr SPAT_VAR
       allocate(n(1:z_num,1:gridNoMax))    !dmr SPAT_VAR
       allocate(Cp(1:z_num,1:gridNoMax))   !dmr SPAT_VAR
       allocate(pori(1:z_num,1:gridNoMax)) !dmr SPAT_VAR
       allocate(porf(1:z_num,1:gridNoMax)) !dmr SPAT_VAR


     END SUBROUTINE allocate_spatialvars


    END MODULE spatialvars_mod
