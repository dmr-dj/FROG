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


      MODULE main_lib_VAMPER


      IMPLICIT NONE

      PRIVATE

      PUBLIC :: INITIALIZE_VAMP

        ! [Initialization of GRID]
        !
        !             -> INITIALIZE VAMP

        ! [defining vertical domain]
        !
        !

     contains


     function INITIALIZE_VAMP() result(is_a_success)

       use grids_more,      only: INIT_maskGRID, nb_unmaskedp
       use parameter_mod,   only: read_namelist, set_numbergridpoints, t_disc, z_disc
       use spatialvars_mod, only: spatialvars_allocate, spatialvars_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical :: is_a_success

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! Function to initialize the domain (SPATIAL / Unique point)
        !   -> sets nb_unmaskedp the number of unmasked points, i.e. nb of computation points
        call INIT_maskGRID

        call set_numbergridpoints(nb_unmaskedp)

        ! Read the namelist to define most global constants
        call read_namelist

        ! Time initialization
        call t_disc

        ! Vertical discretization
        call z_disc

        ! Allocation of main variables
        call spatialvars_allocate

        ! Initialization of spatial variables (2D : z_num,gridpoints)

        call spatialvars_init

     end function INITIALIZE_VAMP






      END MODULE main_lib_VAMPER
