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

      module Fonction_implicit


        implicit none

        private
        public :: Implicit_T

        contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine Implicit_T(z_max,z_snow, T_old,Tu,Tb,dt,dz,n,org_ind,Timp,Kp, rho_snow)


         use parameter_mod, only : z_num, rho_ice, Gfx, T_freeze,rho_snow_freeze,s_l_max
         use Fonction_temp, only : AppHeatCapacity, ThermalConductivity, AppHeatCapacitySnow, ThermalConductivitySnow

         integer                            , intent(in) :: z_max    !dmr maximum number of vertical layers (all included, snow + soil)
         integer,                             intent(in) :: z_snow   !dmr current number of snow layers        [1]

         real, dimension(1:z_max)           , intent(in) :: T_old    !dmr Previous time step soild temperature [C]
         real                               , intent(in) :: Tu       !dmr Temperature forcing at the surface   [C]
         real                               , intent(in) :: Tb       !dmr Temperature fixed at the bottom      [C]
         real                               , intent(in) :: dt       !dmr timestep duration                    [s]
         real, dimension(1:z_max)           , intent(in) :: dz       !dmr layer thickness in the soil          [m]
         real, dimension(:)                 , intent(in) :: n        !dmr porosity in each soil layer          [1]
         integer                            , intent(in) :: org_ind  !dmr organic index ...
         real, dimension(1:z_max)           , intent(out):: Timp     !dmr placeholder for new temperature      [C]
         real, dimension(1:z_max)           , intent(out):: Kp       !dmr placeholder for new Kp per layer     [?]
         real, dimension(1:z_snow), optional, intent(in) :: rho_snow !dmr density of snow per snow layer       [?]

         real, dimension(1:z_max) :: pori, porf, Cp_temp
         real, dimension(1:z_max,1:z_max) :: MM
         real, dimension(1:z_max) :: Knows
         real :: m_Gfx, A, B, C, Z1
         integer :: kk, ll, ll_soil
         real, dimension(1:z_max) :: T_last, Kp_m
         real, dimension(1:z_max) :: T_iter
         real, dimension(1:z_max) :: DD
         real, dimension(1:z_max-1) :: DL, DU
         integer :: info_dgesv
         integer :: z_eff

#if ( SNOW_EFFECT == 1 )
         ! dmr locally computed
         real, dimension(:), allocatable :: Cp_s, Kp_s
#endif


         m_Gfx = gfx/1000.0
         z_eff = z_max-z_num+1


         T_last(1:z_max) = T_old(1:z_max)

#if ( SNOW_EFFECT == 1 )
         if (PRESENT(rho_snow)) then
           allocate(Cp_s(z_snow))
           allocate(Kp_s(z_snow))
         endif
#endif

         do kk=1,5 ! dmr --- why doing this 5 times ??

           MM(1:z_max,1:z_max) = 0
           DD(1:z_max)=0
           DL(1:z_max-1)=0
           DU(1:z_max-1)=0

           if (kk==1) then
            T_iter(1:z_max) = T_last(1:z_max)
           else
            T_iter(1:z_max) = 0.5*(T_iter(1:z_max)+T_last(1:z_max))
           end if


#if ( SNOW_EFFECT == 1 )
           if (PRESENT(rho_snow)) then

              !dmr [NOTA] Proposing a new structure where there could be two layer types: snow and soil in that order
              !dmr        snow: 1->snow
              !dmr        soil: snow->z_max, with the constraint that it entails z_num layers from (z_max-z_num+1):z_max

             call AppHeatCapacitySnow(rho_snow,Cp_s)
             call ThermalConductivitySnow(rho_snow,Kp_s)

             Cp_temp(1:z_snow) = Cp_s(:)
             Kp_m(1:z_snow)    = Kp_s(:)
           endif
#endif

             !dmr Given Temperature and porosity (n), this computes a new Cp value and porf, pori on the vertical
           call AppHeatCapacity(z_num,T_iter(z_eff:z_max),T_freeze,n,org_ind,Cp_temp(z_eff:z_max)      &
                               ,porf(z_eff:z_max),pori(z_eff:z_max))

!~            do ll=z_eff,z_max-1
!~              ll_soil = ll-z_eff+1
!~               !dmr Given the number of layer, porosity, porosities, Temperature, computes the Kp of the layer
!~              call ThermalConductivity(ll,n(ll_soil),pori(ll_soil),porf(ll_soil),org_ind,T_iter(ll),Kp(ll_soil))
!~              Kp(z_max) = 2
!~            end do
           call ThermalConductivity(n,pori(z_eff:z_max),porf(z_eff:z_max),org_ind,T_iter(z_eff:z_max)   &
                              ,Kp_m(z_eff:z_max))


           do ll=1+1,z_max-1

             Z1 = T_last(ll)

             A=(dt/((dz(ll-1)+dz(ll))*0.5*dz(ll)) * (Kp_m(ll-1)/Cp_temp(ll)))
             B=(dt/((dz(ll+1)+dz(ll))*0.5*dz(ll)) * (Kp_m(ll)/Cp_temp(ll)))

             C= 1+A+B

             MM(ll,ll-1) = -A
             MM(ll,ll) = C
             MM(ll,ll+1) = -B
             DL(ll-1) = -A
             DU(ll) = -B
             DD(ll) = C
             Knows(ll) = Z1

           end do

           A=(dt/((dz(z_max-1)+dz(z_max))*0.5*dz(z_max-1)) * (Kp_m(z_max-1)/Cp_temp(z_max-1)))
           C=1.0+A
           Knows(1) = Tu
           Knows(z_max)=Tb
           MM(z_max,z_max)=1
           MM(1,1)=1
           DD(1) = 1
           DD(z_max) = 1+A
           DL(z_max-1) = -A

!~ dmr [INFO] LAPACK CALL
!~ subroutine sgtsv  (  integer  n,
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!~       integer  nrhs,
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!~       real, dimension( * )    dl,
!>          DL is REAL array, dimension (N-1)
!>          On entry, DL must contain the (n-1) sub-diagonal elements of
!>          A.
!>
!>          On exit, DL is overwritten by the (n-2) elements of the
!>          second super-diagonal of the upper triangular matrix U from
!>          the LU factorization of A, in DL(1), ..., DL(n-2).
!~       real, dimension( * )    d,
!>          D is REAL array, dimension (N)
!>          On entry, D must contain the diagonal elements of A.
!>
!>          On exit, D is overwritten by the n diagonal elements of U.
!~       real, dimension( * )    du,
!>          DU is REAL array, dimension (N-1)
!>          On entry, DU must contain the (n-1) super-diagonal elements
!>          of A.
!>
!>          On exit, DU is overwritten by the (n-1) elements of the first
!>          super-diagonal of U.
!~       real, dimension( ldb, * )  b,
!>          B is REAL array, dimension (LDB,NRHS)
!>          On entry, the N by NRHS matrix of right hand side matrix B.
!>          On exit, if INFO = 0, the N by NRHS solution matrix X.
!~       integer  ldb,
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!~       integer  info )
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, U(i,i) is exactly zero, and the solution
!>               has not been computed.  The factorization has not been
!>               completed unless i = N.
!~ dmr [INFO] LAPACK CALL

           call sgtsv(z_num,1,DL,DD,DU,Knows,z_num,info_dgesv)

           T_iter(1:z_max) = Knows(1:z_max)

         end do !! on kk 1,5


         Timp(1:z_max) = T_iter(1:z_max)

        end subroutine Implicit_T

      end module Fonction_implicit
