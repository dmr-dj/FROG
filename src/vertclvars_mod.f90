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

    MODULE vertclvars_mod

     IMPLICIT NONE

     PRIVATE

     PUBLIC :: vertclvars_init, DO_vertclvars_step

     CONTAINS



     SUBROUTINE DO_vertclvars_step(nb_steps_toDO, Kp, T_bottom, Temp, T_air, n, Per_depth,                                       &
                                                                      ! Temp is the one column temperature of soil
                                                                      ! T_air is the surface soil temperature forcing for the timesteps
                                                                      ! T_air should be T_air(nb_steps_toDO) exactly
                                                                      ! n is porosity in the vertical
                                                                      ! Per_depth is the diagnosed "permafrost" or freezing depth (in meters)
                                   ALT, altmax_lastyear, compteur_time_step, organic_indd, deepSOM_a, deepSOM_s, deepSOM_p       &
                                   , deepSOM, fc,  b4_lok, Fv_lok, r_leaf_lok, fracgr_lok, darea_lok                             &
                                   , alpha_a_lok, alpha_s_lok, alpha_p_lok, mu_soil_rev_lok, beta_a_lok, beta_s_lok, beta_p_lok  &
                                   , deepSOM_tot                                                                                 &
                                   , snowlayer_thick_forcing, Temp_snow_col                                                      &
                                   , snowlayer_depth, snowlayer_nb                                                               &
                                   , Tmean_col, Tmmin_col, Tmmax_col)





        USE parameter_mod,     ONLY: z_num ! organic_ind
        USE parameter_mod,     ONLY: D, dt, dz
        use Fonction_temp,     ONLY: diagnose_frost_Depth
        use Fonction_implicit, ONLY: Implicit_T
        use timer_mod,         ONLY: cell_time, update_time_cell

#if ( CARBON == 1 )
        USE parameter_mod, ONLY: bio_diff_k_const, diff_k_const, bioturbation_depth, min_cryoturb_alt, max_cryoturb_alt, zf_soil
        USE parameter_mod, ONLY: YearType
        USE Carbon,        ONLY: carbon_main, ncarb  !carbon_redistribute, decomposition, cryoturbation
#endif

#if ( SNOW_EFFECT == 1 )
        USE simple_snow,   ONLY: get_snow_profile, depth_layer
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        INTEGER                         , INTENT(in)        :: nb_steps_toDO
        REAL, DIMENSION(1:nb_steps_toDO), INTENT(in)        :: T_air
        REAL, DIMENSION(1:z_num)        , INTENT(inout)     :: Kp
        REAL, DIMENSION(1:z_num)        , INTENT(IN)        :: n
        REAL, DIMENSION(1:z_num)        , INTENT(INOUT)     :: Temp
        REAL                            , INTENT(in)        :: T_bottom
        REAL, DIMENSION(3)              , INTENT(out)       :: Per_depth
        REAL, DIMENSION(1:z_num)        , INTENT(INOUT)     :: Tmean_col
        REAL, DIMENSION(1:z_num)        , INTENT(OUT)       :: Tmmin_col
        REAL, DIMENSION(1:z_num)        , INTENT(OUT)       :: Tmmax_col


!~         INTEGER,DIMENSION(1:z_num), INTENT(inout), OPTIONAL :: Temp_positive                    ! Where temp is once positive over one year
        REAL                      , INTENT(inout)           :: ALT, altmax_lastyear             ! Active Layer Thickness
!~         REAL                      , INTENT(in),    OPTIONAL :: clay
        TYPE(cell_time)           , INTENT(inout)           :: compteur_time_step
        INTEGER                   , INTENT(in)              :: organic_indd                     ! index of the end of the organic layer
!~         LOGICAL                   , INTENT(inout), OPTIONAL :: end_year

#if ( CARBON == 1 )
        REAL   ,DIMENSION(1:z_num), INTENT(inout),OPTIONAL  :: deepSOM_a, deepSOM_s, deepSOM_p
        REAL   ,DIMENSION(1:z_num), INTENT(inout), OPTIONAL :: deepSOM
        REAL, DIMENSION(ncarb,ncarb), intent(inout),OPTIONAL:: fc !! flux fractions within carbon pools

        REAL, OPTIONAL,                   INTENT(in)        ::  b4_lok !b3_lok,
        REAL, OPTIONAL,                   INTENT(in)        ::  Fv_lok 
        REAL, OPTIONAL,                   INTENT(in)        ::  r_leaf_lok 
        REAL, OPTIONAL,                   INTENT(in)        ::  fracgr_lok 
        REAL, OPTIONAL,                   INTENT(in)        ::  darea_lok 
        REAL, DIMENSION(1:z_num), OPTIONAL,            INTENT(inout)        ::  alpha_a_lok 
        REAL, DIMENSION(1:z_num), OPTIONAL,            INTENT(inout)        ::  alpha_s_lok 
        REAL, DIMENSION(1:z_num), OPTIONAL,            INTENT(inout)        ::  alpha_p_lok 
        REAL, DIMENSION(1:z_num), OPTIONAL,            INTENT(inout)        ::  beta_a_lok 
        REAL, DIMENSION(1:z_num), OPTIONAL,            INTENT(inout)        ::  beta_s_lok 
        REAL, DIMENSION(1:z_num), OPTIONAL,            INTENT(inout)        ::  beta_p_lok 
        REAL,                     OPTIONAL,            INTENT(inout)        ::  mu_soil_rev_lok 
        REAL, OPTIONAL,            INTENT(inout)     ::  deepSOM_tot 
#endif

          ! SNOW VARIABLES
        REAL, DIMENSION(1:nb_steps_toDO), OPTIONAL, INTENT(in)        :: snowlayer_thick_forcing !! a time series of the snowlayer thickness [m]
        REAL, DIMENSION(:)              , OPTIONAL, INTENT(INOUT)     :: Temp_snow_col   !! Snow temperature, whole column (max_nb_snow_layers)
        INTEGER                         , OPTIONAL, INTENT(INOUT)     :: snowlayer_nb    !! Nb of active snow layers
        REAL                            , OPTIONAL, INTENT(INOUT)     :: snowlayer_depth !! Thickness of snow layers



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        INTEGER                  :: ll
        REAL                     :: T_soil
        REAL, DIMENSION(1:z_num) :: T_old
        LOGICAL                  :: success
        LOGICAL                  :: end_year
        REAL                     :: altmax_thisyear

#if ( SNOW_EFFECT == 1 )
        integer                  :: nb_snowlayers = 0, nb_layers_WC
        real, dimension(:), allocatable :: rho_snow, dz_snowlayers

        real, dimension(:), allocatable :: T_old_WC, dz_WC, Temp_WC, Kp_WC
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Computes the forward step(s) in time for one column
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        altmax_thisyear=0.0

        do ll=1,nb_steps_toDO !boucle temporelle / nombre de pas de temps à faire en avant


         success = update_time_cell(compteur_time_step)
         end_year = compteur_time_step%end_year


         T_soil = T_air(ll)
         T_old(1:z_num) = Temp(1:z_num)

#if ( SNOW_EFFECT == 1 )

         ! [SNOW]
         !    1. From snow thickness, compute number of snow layers and vertical density
         !    2. Construct the merged array for Temperature containing the snow and the normal soil

         if (snowlayer_thick_forcing(ll).GT.depth_layer(1)) then ! there is snow on the ground

            call get_snow_profile(snowlayer_thick_forcing(ll),nb_snowlayers,rho_snow, dz_snowlayers)


            !dmr --- Create merged variables to provide WholeColumn (WC) data

            nb_layers_WC = nb_snowlayers+z_num

            allocate(T_old_WC(nb_layers_WC))
            allocate(dz_WC(nb_layers_WC))
            allocate(Temp_WC(nb_layers_WC))
            allocate(Kp_WC(nb_layers_WC))

            T_old_WC(1:nb_snowlayers) = Temp_snow_col(1:nb_snowlayers) ! do not take into account the fact that the nb_snow layers may have varied ... data presence?
            dz_WC(1:nb_snowlayers)    = dz_snowlayers(1:nb_snowlayers)
            T_old_WC(nb_snowlayers+1:nb_layers_WC) = T_old(:)
            dz_WC(nb_snowlayers+1:nb_layers_WC) = dz(:)


!~ PROGNOSTIC VARIABLES FROM SNOW THAT NEED UPDATING
!~ snowlayer_thick_forcing, Temp_snow_col     &
!~                                    , snowlayer_depth, snowlayer_nb
!~          WRITE(*,*) "Column CALL", Temp_snow_col
!~             WRITE(*,*) "SNOW CASE nb_snowlayers = ", nb_snowlayers, snowlayer_thick_forcing(ll)

         ! would need a forcing here in terms of Snow
         ! This need to provide z_snow, rho_snow
         ! i.e. the number of snow layers and their respective density
         ! Temperature of the snow layers should be provided as well (through an extended T_variable ...)

!~          real, dimension(1:z_max)           , intent(in) :: T_old    !dmr Previous time step soild temperature [C]
!~          real, dimension(1:z_max)           , intent(in) :: dz       !dmr layer thickness in the soil          [m]
!~          real, dimension(1:z_max)           , intent(out):: Timp     !dmr placeholder for new temperature      [C]
!~          real, dimension(1:z_max)           , intent(out):: Kp       !dmr placeholder for new Kp per layer     [?]


            !! WRITE(*,*) "organic_ind === ", organic_indd
            call Implicit_T(nb_layers_WC,nb_snowlayers,T_old_WC,T_soil,T_bottom,dt,dz_WC,n,organic_indd,Temp_WC,Kp_WC,&
                             rho_snow=rho_snow(1:nb_snowlayers))

            !dmr Then update Temp with the lower part
            Temp(1:z_num) = Temp_WC(nb_snowlayers+1:nb_layers_WC)
            Kp(1:z_num)   = Kp_WC(nb_snowlayers+1:nb_layers_WC)
            Temp_snow_col(1:nb_snowlayers) = Temp_WC(1:nb_snowlayers)
            snowlayer_depth = SUM(dz_snowlayers(1:nb_snowlayers),DIM=1)
            snowlayer_nb    = nb_snowlayers

!~             WRITE(*,*) "STATS WITH SNOW"
!~             WRITE(*,*) Temp_snow_col(1:nb_snowlayers)
!~             WRITE(*,*) snowlayer_depth, snowlayer_nb
!~             WRITE(*,*) Kp(1:nb_snowlayers)

            deallocate(T_old_WC)
            deallocate(dz_WC)
            deallocate(Temp_WC)
            deallocate(Kp_WC)

         else ! No Snow, proceeding normally

            !-------------- Numerical difference routine when there is no snow --------!
            call Implicit_T(z_num, 0, T_old,T_soil,T_bottom,dt,dz,n,organic_indd,Temp,Kp)


         endif


#else

       !-------------- Numerical difference routine when there is no snow --------!

         call Implicit_T(z_num, 0, T_old,T_soil,T_bottom,dt,dz,n,organic_indd,Temp,Kp)



#endif



                       !> Returns Per_depth as depth of the 0C isotherm
                       !>    in meters, not per cell, also the current ALT as last term
         Per_depth = diagnose_frost_Depth(Temp,D)
         altmax_thisyear = MAX(altmax_thisyear,Per_depth(3)) ! the maximum active layer thickness so far this year

#if ( CARBON == 1 )
       ! nb and mbv carbon cycle call
       ! at the end of each year computes the actve layer thickness, needed for redistribution
!~          call compute_alt(Temp, Temp_positive, ALT, compteur_time_step, end_year, altmax_lastyear, D)
       !write(*,*) 'ALT', ALT
       call carbon_main (Temp, ALT, deepSOM_a, deepSOM_s, deepSOM_p, max_cryoturb_alt,                &
                          min_cryoturb_alt, diff_k_const, bio_diff_k_const, bioturbation_depth,       &
                          deepSOM, fc, Fv_lok, r_leaf_lok, fracgr_lok, darea_lok, deepSOM_tot,        &
                          alpha_a_lok, alpha_s_lok, alpha_p_lok,                                      &
                          mu_soil_rev_lok, beta_a_lok, beta_s_lok, beta_p_lok)
#endif

         if (end_year) then
           altmax_lastyear = ALT
           ALT = altmax_thisyear
           altmax_thisyear = 0.0

           !dmr [TODO] send the stored yearly temp result to the routine that prepares output

         endif

        !dmr OUTPUT section
         
         Tmean_col(1:z_num) = Tmean_col(1:z_num) + Temp(1:z_num) / nb_steps_toDo 
         where (Temp(1:z_num) .LT. Tmmin_col(1:z_num))
            Tmmin_col = Temp
         endwhere
         where (Temp(1:z_num) .GT. Tmmax_col(1:z_num))
            Tmmax_col = Temp
         endwhere

        enddo ! boucle temporelle

     END SUBROUTINE DO_vertclvars_step


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Computes the 1-D characteristics, inputs/outputs ar for one column or one point
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     SUBROUTINE vertclvars_init(Gfx_loc, Tinit_loc, Kp_loc, Cp_loc, organic_ind_loc, porosity_profvertcl, temperature_profvertcl  &
                              , T_bottom)

       use parameter_mod, only: PorosityType, D, Bool_Organic, organic_depth, z_num !organic_ind,
       use parameter_mod, only: dz, T_freeze
       use Fonction_temp, only: AppHeatCapacity, ThermalConductivity
       use simple_snow,   only: init_snow_profile

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real,                       intent(in)  :: Gfx_loc, Tinit_loc                                   ! Geothermal flux and mean temperature locally


       real, dimension(1:z_num),   intent(out) :: porosity_profvertcl, temperature_profvertcl   ! vertical 1-D porosity and temp profile
       real, dimension(1:z_num-1), intent(out) :: Kp_loc                                        ! Kp first constant, then computed
       real, dimension(1:z_num-1), intent(out) :: Cp_loc                                        ! Kp first constant, then computed
       integer,                    intent(out) :: organic_ind_loc                               ! indx of the bottom of organic layer
       real,                       intent(out) :: T_bottom                                      ! Computation of bottom temperature (last level deep)


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer                  :: kk, ll
       real, dimension(z_num-1) :: h_n, h_pori, h_porf
       real, dimension(z_num)   :: porf,pori


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr    [2024-06-28] CALCULATION OF POROSITY in the vertical
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       !dmr
       !dmr intent(in)                PorosityType == 1 or 2 [TO_BE_CLARIFIED]
       !dmr intent(in)                Bool_Organic == 1 use organic layer, else not
       !dmr intent(in)                organic_depth = depth of organic layer (I guess), in meters
       !dmr intent(in) (z_num)        D depth of the layer considered in meters
       !dmr intent(out) (allocatable) n = porosity of each layer in the vertical
       !dmr intent(out)               organic_ind, index in vertical of the end of organic layer (1:organic_ind)
       call Porosity_init(PorosityType, D, Bool_Organic, organic_depth, porosity_profvertcl, organic_ind_loc)  !CALCULATION OF POROSITY


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr    [NOTA] For now it seems that Kp is constant, are there reasons to have it spatially or vertically variable?
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       do kk=1,UBOUND(Kp_loc,DIM=1) !dmr Kp is size z_num-1
         Kp_loc(kk)=2
       end do

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr    Computes the initial vertical 1-D profile of temperature in the soil from T_init and Geoheatflow Gfx
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       call GeoHeatFlow(Gfx_loc, Kp_loc, dz, Tinit_loc, temperature_profvertcl)

       T_bottom = temperature_profvertcl(z_num)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr    Then computes Cp and porf, pori and thermal conductivity of soil Kp
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       do ll = 1,2 !dmr WhatIs ll <- Repeating the steps two times?

                                          !Calculation of heat capacity of soil, Cp, porf & pori are intent(out)
          call AppHeatCapacity(z_num,temperature_profvertcl,T_freeze,porosity_profvertcl, organic_ind_loc, Cp_loc, porf, pori)

!~           do kk=1,z_num-1

!~              h_pori(kk) = (pori(kk) + pori(kk+1))/2
!~              h_porf(kk) = (porf(kk) + porf(kk+1))/2
!~              h_n(kk) = (porosity_profvertcl(kk) + porosity_profvertcl(kk+1))/2
!~              call ThermalConductivity(kk,h_n(kk),h_pori(kk),h_porf(kk), organic_ind_loc, temperature_profvertcl(kk), Kp_loc(kk))
!~           end do

                                          !Calculation of the thermal condutivity of soil, Kp is the only intent(out)
          call ThermalConductivity(h_n,h_pori,h_porf, organic_ind_loc, temperature_profvertcl, Kp_loc)

       end do

#if ( SNOW_EFFECT == 1 )
       call init_snow_profile()
#endif


     END SUBROUTINE vertclvars_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr    Computes an initial guess of the 1-D vertical temperature profile from a T0 and knowing Kp(Z) and Gfx at the bottom
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     SUBROUTINE GeoHeatFlow(Gfx, Kp, dz, T0, T)

       use parameter_mod, only: z_num

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real, intent(in)                             :: Gfx, T0 ! Gfx = Earth's geothermal heat flux, T0 is modified top temperature (?)
       real, dimension(z_num)  , intent(in)         :: dz ! vertical geometry (thickness of layer)
       real, dimension(z_num-1), intent(in)         :: Kp ! dmr [2024-06-28] : Kp is z_num-1 only and constant to 2
       real, dimension(z_num)  , intent(out)        :: T  ! new vertical temperature profile

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer :: kk ! index

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       T(1) = T0

       do kk = 2, z_num
          T(kk) = T(kk-1) + (((Gfx/1000.0)/Kp(kk-1))*dz(kk-1))
       end do

     END SUBROUTINE GeoHeatFlow

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr    Computes the vertical 1-D porosity profile
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     SUBROUTINE Porosity_init(Porosity_Type, Depth_layer, Bool_Organic, organic_depth, n, organic_ind)

       use parameter_mod, only : n_organic, Porosity_soil, z_num

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer, intent(in)            :: Porosity_Type, Bool_Organic
       real, intent(in)               :: organic_depth
       real, dimension(:), intent(in) :: Depth_layer
       integer, intent(out)           :: organic_ind
       real, dimension(:), intent(out):: n                              ! this is the porosity calculated

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real, dimension(:), allocatable :: DepthCalcOrg
       real                            :: pente
       real                            :: origine
       integer                         :: pas_z

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       allocate(DepthCalcOrg(1:z_num))

       if (Bool_organic == 1) then

         do pas_z = 1,z_num
          DepthCalcOrg(pas_z) = abs(Depth_layer(pas_z) - organic_depth)
         end do

         call indice_minimum(DepthCalcOrg, z_num, organic_ind)

         do pas_z = 1,organic_ind
            n(pas_z) = n_organic
         end do

       else ! Bool_organic != 1

         n(1) = Porosity_soil
         organic_ind = 1

       end if


       if (Porosity_Type == 1) then

          pente = (Porosity_soil -0.15)/(Depth_layer(organic_ind) - Depth_layer(z_num))
          origine = Porosity_soil - pente*Depth_layer(organic_ind)

          do pas_z = organic_ind, z_num
            n(pas_z) = pente * Depth_layer(pas_z) + origine
          end do

       elseif(Porosity_Type == 2) then

          do pas_z = organic_ind, z_num
            n(pas_z) = Porosity_soil
          end do

       else ! Porosity_Type != 1,2

          do pas_z = organic_ind, z_num
            n(pas_z) = Porosity_soil*exp(-0.000395*Depth_layer(pas_z))
          end do

       end if

     END SUBROUTINE Porosity_init


     SUBROUTINE indice_minimum(tab, taille, ind_min)

       implicit none

       integer, intent(in) :: taille
       real, dimension(:), intent(in) :: tab
       integer, intent(out) :: ind_min

       integer :: kk

       ind_min = 1

       do kk = 2,taille,1

          if(tab(kk)<tab(ind_min)) then
             ind_min = kk
          end if

       end do

     END SUBROUTINE indice_minimum


    END MODULE vertclvars_mod
