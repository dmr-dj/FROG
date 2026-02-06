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


    MODULE main_lib_FROG

#include "constant.h"

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: INITIALIZE_VAMP, STEPFWD_VAMP, GET_COUPLING_STEP
#if ( CARBON > 0 )
      PUBLIC :: INITIALIZE_CARBON_STOCK
#endif

      INTEGER :: nb_coupling_steps

      TYPE cpl_fields
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: TempForc
#if ( CARBON == 1 )
        !REAL, DIMENSION(:,:),   ALLOCATABLE :: B3_vegForc
        REAL, DIMENSION(:,:),   ALLOCATABLE :: B4_vegForc
        REAL, DIMENSION(:,:),   ALLOCATABLE :: Fv_vegForc
        REAL, DIMENSION(:,:),   ALLOCATABLE :: r_leaf_vegForc
        REAL, DIMENSION(:,:),   ALLOCATABLE :: fracgr_vegForc
        REAL, DIMENSION(:,:),   ALLOCATABLE :: darea_vegForc
#endif

!#if ( SNOW_EFFECT == 1 )
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: dsnow_thick
!#endif

      END TYPE cpl_fields

      PUBLIC :: cpl_fields

     contains


     function SET_COUPLING_STEP() result(is_a_success)

     use parameter_mod, only: YearType

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical :: is_a_success

       nb_coupling_steps = YearType ! for now, implement a first version of the coupling on a yearly basis
       is_a_success = .true.

     end function SET_COUPLING_STEP


     function GET_COUPLING_STEP() result(the_coupling_step)

       integer :: the_coupling_step

       the_coupling_step = nb_coupling_steps

     end function GET_COUPLING_STEP



     function INITIALIZE_VAMP() result(is_a_success)

       use grids_more,      only: INIT_maskGRID, nb_unmaskedp, forcing_timelength, INIT_netCDF_output
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

        ! Coupling timestep (a.k.a. internal timestep)
        is_a_success = SET_COUPLING_STEP()

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

        call INIT_netCDF_output

        is_a_success = .TRUE.

     end function INITIALIZE_VAMP

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!--------------------------
  function INITIALIZE_CARBON_STOCK(coupled_fields) result(is_a_success)

#if (OFFLINE_RUN == 1)
       USE spatialvars_mod, ONLY: UPDATE_climate_forcing
#else
       USE spatialvars_mod, ONLY: SET_coupled_climate_forcing
       USE grids_more, ONLY: flatten_it_3D, flatten_it
#endif
       use spatialvars_mod, only: spatialvars_init_carbon


    logical :: is_a_success
    REAL, DIMENSION(:,:), ALLOCATABLE :: temperature_forcing_nextsteps, snowthickness_forcing_nextsteps

#if ( OFFLINE_RUN == 0 )
       type(cpl_fields), intent(in), optional :: coupled_fields
#else
       logical         , intent(in), optional :: coupled_fields ! dummy unused
#endif

#if (OFFLINE_RUN == 1)
        ! UPDATE_CLIMATE_FORCING
    CALL UPDATE_climate_forcing(nb_coupling_steps,temperature_forcing_nextsteps                                                 &
#if ( SNOW_EFFECT == 1 )
                                  , snowthickness_forcing_nextsteps                                                             &
#endif
                                   )
#else
        ! FORCING is coming from the coupled component
        if (PRESENT(coupled_fields)) then

    CALL SET_coupled_climate_forcing(nb_coupling_steps, temperature_forcing_nextsteps,                                          &
                   coupled_temp_set = flatten_it_3D(coupled_fields%TempForc,UBOUND(coupled_fields%TempForc,dim=3))              &
#if ( CARBON == 1 )
                 , b4_content = flatten_it(TRANSPOSE(coupled_fields%B4_vegForc(:,:)))                                           &
                 , Fv_content = flatten_it(TRANSPOSE(coupled_fields%Fv_vegForc(:,:)))                                           &
                 , r_leaf_content = flatten_it(TRANSPOSE(coupled_fields%r_leaf_vegForc(:,:)))                                           &
                 , fracgr_content = flatten_it(TRANSPOSE(coupled_fields%fracgr_vegForc(:,:)))                                   &
                 , darea_content = flatten_it(TRANSPOSE(coupled_fields%darea_vegForc(:,:)))                                     &
#endif
#if ( SNOW_EFFECT == 1 )
                 , snowthick_forc_nxt = snowthickness_forcing_nextsteps                                                         &
                 , coupled_dsnow_set  = flatten_it_3D(coupled_fields%dsnow_thick,UBOUND(coupled_fields%dsnow_thick,dim=3))      &
#endif
                                          )

!dmr [TODO] UPDATED NEED FOR THE COUPLED CASE RE. SNOW THICKNESS !!!

        else
          WRITE(*,*) "[ABORT] :: we are in coupled setup, need a forcing input fields"
        endif

#endif
        call spatialvars_init_carbon

        is_a_success = .TRUE.


  end function INITIALIZE_CARBON_STOCK
!--------------------------

     function STEPFWD_VAMP(coupled_fields) result(is_a_success)

       USE spatialvars_mod, ONLY: DO_spatialvars_step
#if (OFFLINE_RUN == 1)
       USE spatialvars_mod, ONLY: UPDATE_climate_forcing
#else
       USE spatialvars_mod, ONLY: SET_coupled_climate_forcing
       USE grids_more, ONLY: flatten_it_3D, flatten_it
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical :: is_a_success

       REAL, DIMENSION(:,:), ALLOCATABLE :: temperature_forcing_nextsteps, snowthickness_forcing_nextsteps

#if ( OFFLINE_RUN == 0 )
       type(cpl_fields), intent(in), optional :: coupled_fields
#else
       logical         , intent(in), optional :: coupled_fields ! dummy unused
#endif

!~        REAL, DIMENSION(:,:,:), INTENT(in), OPTIONAL :: coupled_temp_set     ! will be (spat_coord1, spat_coord2, 1:stepstoDo)
!~                                                                             ! must be (1:gridNoMax,1:stepstoDO)
!~        REAL, DIMENSION(:,:), INTENT(in), OPTIONAL :: coupled_b3, coupled_b4 ! will be (spat_coord1, spat_coord2)
!~                                                                             ! must be (1:gridNoMax)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if (OFFLINE_RUN == 1)
        ! UPDATE_CLIMATE_FORCING
        CALL UPDATE_climate_forcing(nb_coupling_steps,temperature_forcing_nextsteps &
#if ( SNOW_EFFECT == 1 )
                                  , snowthickness_forcing_nextsteps                 &
#endif
                                   )
#else
        ! FORCING is coming from the coupled component
        if (PRESENT(coupled_fields)) then

          CALL SET_coupled_climate_forcing(nb_coupling_steps, temperature_forcing_nextsteps,                                    &
                   coupled_temp_set = flatten_it_3D(coupled_fields%TempForc,UBOUND(coupled_fields%TempForc,dim=3))              &
#if ( CARBON == 1 )
                 !, b3_content = flatten_it(TRANSPOSE(coupled_fields%B3_vegForc(:,:)))                                           &
                 , b4_content = flatten_it(TRANSPOSE(coupled_fields%B4_vegForc(:,:)))                                           &
                 , Fv_content = flatten_it(TRANSPOSE(coupled_fields%Fv_vegForc(:,:)))                                           &
                 , r_leaf_content = flatten_it(TRANSPOSE(coupled_fields%r_leaf_vegForc(:,:)))                                           &
                 , fracgr_content = flatten_it(TRANSPOSE(coupled_fields%fracgr_vegForc(:,:)))                                   &
                 , darea_content = flatten_it(TRANSPOSE(coupled_fields%darea_vegForc(:,:)))                                     &
#endif
#if ( SNOW_EFFECT == 1 )
                 , snowthick_forc_nxt = snowthickness_forcing_nextsteps                                                         &
                 , coupled_dsnow_set  = flatten_it_3D(coupled_fields%dsnow_thick,UBOUND(coupled_fields%dsnow_thick,dim=3))      &
#endif
                                          )

!dmr [TODO] UPDATED NEED FOR THE COUPLED CASE RE. SNOW THICKNESS !!!

        else
          WRITE(*,*) "[ABORT] :: we are in coupled setup, need a forcing input fields"
        endif

#endif
        ! DO_FROG_STEP
        CALL DO_spatialvars_step(nb_coupling_steps, temperature_forcing_nextsteps                                               &
#if ( SNOW_EFFECT == 1 )
                               , snowthickness_forcing_nextsteps                                                                &
#endif
                                )

        is_a_success = .TRUE.

     end function STEPFWD_VAMP


    END MODULE main_lib_FROG
