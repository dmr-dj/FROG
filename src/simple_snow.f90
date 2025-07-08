!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!
!   Copyright 2025 FRATRES-E (https://github.com/FRATRES-E)
!     FRamework for fAst TRansient Earth-system Studies and Education

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

      module simple_snow

       !! version: v1.0
       !! display: public private protected
       !! proc_internals: true


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [Simple Snow]

!!     @author  Didier M. Roche  (dmr)
!!     @date Creation date: July, 03rd, 2024

!!     @brief This module is intended to provide the base computation for simple snow calculations for FROG

!>
!>     DESCRIPTION : Here add the long_description of the module ...
!>        - Subroutine get_snow_profile :
!>        - Formula:
!> $$ \frac{d\lambda}{dt} , \frac{d\phi}{dt},  \frac{dz}{dt} $$
!>     @reference References: papers or other documents to be cited... [Site the website if possible](https://iloveclim.eu)
!>     @date Last modification: 2025-07-08
!>     @author Last modified by : dmr

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       implicit none
       private

       ! dmr -- PUBLIC variables, functions, subroutines.
       public :: get_snow_profile, init_snow_profile

       integer             , parameter, public :: max_nb_snow_layers = 8
       real, dimension(0:8), parameter, public :: depth_layer = (/0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 10.0/) ! in [m]
       real, dimension(0:8), parameter         :: thick_layer = (/0.0, 0.1, 0.15, 0.25, 0.5, 1.0, 2.0, 2.0, 4.0/) ! in [m]
       real, dimension(0:8)                    :: rho_layer, midpt_layer

      contains


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      SUBROUTINE init_snow_profile

!!      AUTHOR : Didier M. Roche a.k.a. dmr
!!      DESCRIPTION: Initialize the vertical profiles from the given fixed data and the density function
!!      REFERENCES:

       integer :: i

       rho_layer(:) = 0
       midpt_layer(:) = 0

!~        WRITE(*,*) "CHECK RHO_SNOW", midpt_layer(0), rho_layer(0)
       do i=1,max_nb_snow_layers
         midpt_layer(i) = depth_layer(i-1) + thick_layer(i)/2.
         rho_layer(i) = get_rhowsnow(midpt_layer(i))
       enddo

      END SUBROUTINE init_snow_profile

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      SUBROUTINE get_snow_profile(snowlayer_thick,nb_snowlayers,rho_snow, dz_snowlayers)


!!      AUTHOR : Didier M. Roche a.k.a. dmr
!!      DESCRIPTION: Subroutine is providing the number of snow layers to be included and their vertical densities
!!      REFERENCES:

       real, intent(in)         :: snowlayer_thick                   !! Snow thickness [m]
       integer, intent(out)     :: nb_snowlayers                     !! Nb of snow_layers to use
       real, dimension(:), allocatable, intent(out) :: rho_snow      !! Density of snow layers computed from Sakahashi, 1965
       real, dimension(:), allocatable, intent(out) :: dz_snowlayers !! Actual thickness of the snow layers [m]

       ! Local variables
       logical :: add_layer=.false.
       ! Begining of the subroutine

           !dmr find what the thickness means in terms of nb of layers given the discretization above
       nb_snowlayers = get_indxby_low(depth_layer,snowlayer_thick)-1 !dmr need to remove one since index of the array starts at zero ...

       if (((snowlayer_thick-depth_layer(nb_snowlayers)).GT.depth_layer(1)).and.(nb_snowlayers.LT.max_nb_snow_layers)) then !dmr there is quite a bit of snow below last layer, add
          add_layer = .true.
          nb_snowlayers = nb_snowlayers + 1
       endif

       allocate(rho_snow(1:nb_snowlayers))
       allocate(dz_snowlayers(1:nb_snowlayers))

       rho_snow(1:nb_snowlayers) = rho_layer(1:nb_snowlayers)

       !dmr There is one last incomplete layer
       if (add_layer) then
         dz_snowlayers(1:nb_snowlayers-1) = thick_layer(1:nb_snowlayers-1)
         dz_snowlayers(nb_snowlayers) = snowlayer_thick-depth_layer(nb_snowlayers-1)
       else
         dz_snowlayers(1:nb_snowlayers) = thick_layer(1:nb_snowlayers)
       endif

!~        WRITE(*,*) "Finally check discretization", SUM(dz_snowlayers(1:nb_snowlayers),DIM=1), " | ", snowlayer_thick, "//" &
!~                 , nb_snowlayers


      END SUBROUTINE get_snow_profile

      function get_indxby_low(array,value) result(index_low)

        real, dimension(:), intent(in) :: array
        real,               intent(in) :: value

        integer :: index_low

        index_low = minloc(abs(array(:)-value),DIM=1) !dmr returns closest, but might be above (arrays assumed sorted)
        if (array(index_low).GT.value) then
          index_low = index_low - 1
        endif

        return
      end function get_indxby_low
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function get_rhowsnow(depth_inlayer) result(densityofsnow)


!!      AUTHOR : Didier M. Roche a.k.a. dmr
!!      DESCRIPTION: Subroutine is for computing the density of snow knowing the depth in the snow layer
!!      REFERENCES: Sakahashi, 1965, https://repository.kulib.kyoto-u.ac.jp/dspace/bitstream/2433/178487/1/sgk00500_065.pdf


       real, intent(in)    :: depth_inlayer     !! [m]
       real                :: densityofsnow     !! [kg.m-3]

       real, parameter     :: aa = 500.0 !! density of snow max       [kg.m-3]
       real, parameter     :: bb = 150.0 !! density of snow min       [kg.m-3]
       real, parameter     :: cc = 0.6   !! exponential decrease coef [m-1]

         densityofsnow = aa - (aa-bb)*exp((-1.0)*cc*depth_inlayer)

      end function get_rhowsnow


end module simple_snow

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
