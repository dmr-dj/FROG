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

#include "constant.h"

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: INITIALIZE_VAMP

     contains


     function INITIALIZE_VAMP() result(is_a_success)

       use grids_more,      only: INIT_maskGRID, nb_unmaskedp, forcing_timelength
       use parameter_mod,   only: read_namelist, set_numbergridpoints, set_numberforcingsteps, t_disc, z_disc
       use spatialvars_mod, only: spatialvars_allocate, spatialvars_init
#if ( CARBON == 1 )
       use carbon,          only: carbon_first_init
#endif
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
        call set_numberforcingsteps(forcing_timelength)

        ! Read the namelist to define most global constants
        call read_namelist

        ! Time initialization
        call t_disc

        ! Vertical discretization
        call z_disc

        ! Allocation of base 1-D variables in Carbon and init constants
#if ( CARBON == 1 )
        call carbon_first_init
#endif

        ! Allocation of main variables
        call spatialvars_allocate

        ! Initialization of spatial variables (2D : z_num,gridpoints)
        ! [NOTA] spatialvars_init needs a T_init and a GeoHFlux ...
        !        need to read them before if to be spatialized
        !        For now [2025-04-16], fixed to constants in parameter_mod
        call spatialvars_init

        is_a_success = .TRUE.

     end function INITIALIZE_VAMP

     function STEPFWD_VAMP() result(is_a_success)

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


        is_a_success = .TRUE.

        ! UPDATE_CLIMATE_FORCING
        ! DO_VAMPER_STEP


     end function STEPFWD_VAMP


    END MODULE main_lib_VAMPER
