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


    MODULE vertclvars_mod

    IMPLICIT NONE

    PRIVATE


     PUBLIC :: vertclvars_init


     CONTAINS


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Computes the 1-D characteristics, inputs/outputs ar for one column or one point
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     SUBROUTINE vertclvars_init(Gfx_loc, Tinit_loc, Kp_loc, organic_ind_loc, porosity_profvertcl, temperature_profvertcl)

       use parameter_mod, only: PorosityType, D, Bool_Organic, organic_depth, organic_ind, z_num
       use parameter_mod, only: dz

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real,                       intent(in)  :: Gfx_loc, Tinit_loc                                   ! Geothermal flux and mean temperature locally


       real, dimension(1:z_num),   intent(out) :: porosity_profvertcl, temperature_profvertcl   ! vertical 1-D porosity and temp profile
       real, dimension(1:z_num-1), intent(out) :: Kp_loc                                        ! Kp constant at present
       integer,                    intent(out) :: organic_ind_loc                               ! indx of the bottom of organic layer

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer                         :: kk

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


  subroutine indice_minimum(tab, taille, ind_min)

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



  end subroutine indice_minimum


    END MODULE vertclvars_mod
