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

#include "constant.h"

    IMPLICIT NONE

    PRIVATE


            ! SPATIAL GLOBAL VARIABLES

     real, dimension(:,:),allocatable, PUBLIC :: Temp      & !dmr [SPAT_VAR], soil temperature over the vertical // prognostic
                                                ,Kp        & !dmr [CNTST]     heat conductivity constant over the depth, current value is 2
                                                ,n         & !dmr [SPAT_VAR], porosity on the vertical
                                                ,Cp        & !dmr [SPAT_VAR]  specific heat capacity
                                                ,pori      & !dmr [???  TBC]
                                                ,porf        !dmr [???  TBC]

     real, dimension(:), allocatable, PUBLIC :: GeoHFlux   &
                                               ,Tinit_SV


     integer, dimension(:), allocatable, PUBLIC :: orgalayer_indx


#if ( CARBON == 1 )
     real,dimension(:,:)  , allocatable, PUBLIC :: deepSOM_a & !dmr [TBD]
                                                 , deepSOM_s & !dmr [TBD]
                                                 , deepSOM_p   !dmr [TBD]
     real, dimension(:)   , allocatable, PUBLIC :: clay_SV
     real,dimension(:,:,:), allocatable, PUBLIC :: fc_SV       !dmr [TBD]
#endif



     PUBLIC:: spatialvars_allocate, spatialvars_init

     CONTAINS


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Allocation of two dimensional variables (VERTCL, SPAT_VAR)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     SUBROUTINE spatialvars_allocate ! VERTCL, SPAT_VAR

       use parameter_mod, only: gridNoMax, z_num
       use carbon       , only: ncarb

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

       allocate(GeoHFlux(1:gridNoMax))
       allocate(Tinit_SV(1:gridNoMax))


       allocate(orgalayer_indx(1:gridNoMax))

#if ( CARBON == 1 )
                        !nb and mbv Carbon cycle
       allocate(deepSOM_a(1:z_num,1:gridNoMax))
       allocate(deepSOM_s(1:z_num,1:gridNoMax))
       allocate(deepSOM_p(1:z_num,1:gridNoMax))
       allocate(fc_SV(1:ncarb,1:ncarb,1:gridNoMax))
       allocate(clay_SV(1:gridNoMax))
#endif

     END SUBROUTINE spatialvars_allocate


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Allocation of two dimensional variables (VERTCL, SPAT_VAR)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     SUBROUTINE spatialvars_init ! VERTCL, SPAT_VAR

       use parameter_mod,  only: gridNoMax, z_num
       use vertclvars_mod, only: vertclvars_init

           ! Temporary addendum [2025-04-16]
       use parameter_mod,  only: Gfx, T_init

#if ( CARBON == 1 )
       use carbon        , only: carbon_init
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       integer :: gridp
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

            !dmr [NOTA] For now, dummy init of spatial variables based on constants
        GeoHFlux(:) = Gfx
        Tinit_SV(:) = T_init

            !dmr Initialization of all columns, one by one
        do gridp = 1, gridNoMax
          call vertclvars_init(GeoHFlux(gridp), Tinit_SV(gridp), Kp(:,gridp),Cp(:,gridp), orgalayer_indx(gridp), n(:,gridp) &
                             , Temp(:,gridp))

#if ( CARBON == 1 )
          call carbon_init(deepSOM_a(:,gridp), deepSOM_s(:,gridp), deepSOM_p(:,gridp), fc_SV(:,:,gridp), clay_SV(gridp))
#endif
        enddo


     END SUBROUTINE spatialvars_init



    END MODULE spatialvars_mod
