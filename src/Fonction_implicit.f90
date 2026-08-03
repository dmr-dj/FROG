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

         !dmr [NOTA] This routine calls sgtsv, the SINGLE precision LAPACK
         !dmr        tridiagonal solver. FROG must therefore be built in single
         !dmr        precision. If the working precision is ever changed (via
         !dmr        -fdefault-real-8, or by introducing an explicit kind
         !dmr        parameter), switch to dgtsv accordingly: passing 8-byte
         !dmr        reals to sgtsv corrupts memory silently.

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
         real, dimension(1:z_max) :: Knows
         real :: m_Gfx, A, B, C, Z1
         integer :: kk, ll
         real, dimension(1:z_max) :: T_last, Kp_m
         real, dimension(1:z_max) :: T_iter
!dmr&clo --- G3: fixed-point convergence control. The loop iterates on the
!dmr&clo     Cp(T) non-linearity (latent heat near freezing). Previously it ran
!dmr&clo     a fixed 5 times ("why doing this 5 times ??"); now it stops once
!dmr&clo     the solution stops moving. max_iter is kept at 5 so the behaviour
!dmr&clo     is IDENTICAL to before when convergence is not reached earlier --
!dmr&clo     the criterion can only save iterations, never change a converged
!dmr&clo     result. tol is on temperature [same units as T, i.e. degC/K].
         real, dimension(1:z_max) :: T_prev_sol
         real, parameter :: tol = 1.0e-4
         integer, parameter :: max_iter = 10
         real :: dT_max
         real, dimension(1:z_max) :: DD
         real, dimension(1:z_max-1) :: DL, DU
         integer :: info_lapack
         integer :: z_eff

#if ( SNOW_EFFECT == 1 )
         ! dmr locally computed
         real, dimension(:), allocatable :: Cp_s, Kp_s
#endif


         m_Gfx = gfx/1000.0
         z_eff = z_max-z_num+1

         !dmr Guard: sgtsv expects 4-byte reals. Fail loudly rather than
         !dmr corrupting memory if the build precision ever changes.
         if (kind(Knows) /= kind(1.0e0)) then
           WRITE(*,*) "[ABORT] Implicit_T: sgtsv requires single precision reals, got kind =", kind(Knows)
           STOP
         endif


         T_last(1:z_max) = T_old(1:z_max)

#if ( SNOW_EFFECT == 1 )
         if (PRESENT(rho_snow)) then
           allocate(Cp_s(z_snow))
           allocate(Kp_s(z_snow))
         endif
#endif

!dmr&clo --- G3: iterate to convergence (was a fixed 5 passes)
         T_prev_sol(1:z_max) = T_last(1:z_max)   !dmr&clo seed for the 1st delta

         do kk=1,max_iter

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

             DL(ll-1) = -A
             DU(ll) = -B
             DD(ll) = C
             Knows(ll) = Z1

           end do

           A=(dt/((dz(z_max-1)+dz(z_max))*0.5*dz(z_max-1)) * (Kp_m(z_max-1)/Cp_temp(z_max-1)))
           C=1.0+A
           Knows(1) = Tu
           Knows(z_max)=Tb
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

           !dmr The system spans z_max rows (snow + soil), not z_num. Using
           !dmr z_num here left the deepest nb_snowlayers rows unsolved
           !dmr whenever snow was present, detaching the bottom boundary
           !dmr condition (geothermal heat flux).
           call sgtsv(z_max,1,DL,DD,DU,Knows,z_max,info_lapack)

           if (info_lapack /= 0) then
             WRITE(*,*) "[ABORT] Implicit_T: sgtsv failed, info =", info_lapack
             WRITE(*,*) "        z_max =", z_max, " z_snow =", z_snow
             STOP
           endif

!dmr&clo --- G3: convergence on the raw solution. Knows is the freshly solved
!dmr&clo     temperature; T_iter below is the damped average used to seed the
!dmr&clo     next pass, so the criterion must compare Knows to the previous
!dmr&clo     pass's Knows, NOT T_iter to T_last.
           dT_max = MAXVAL(ABS(Knows(1:z_max) - T_prev_sol(1:z_max)))
           T_prev_sol(1:z_max) = Knows(1:z_max)

           T_iter(1:z_max) = Knows(1:z_max)

           if (dT_max < tol) exit    !dmr&clo Cp(T) has stopped moving

         end do !! fixed-point loop on Cp(T), up to max_iter passes
!!         WRITE(*,*) "kk convergence == ", kk, dT_max


         Timp(1:z_max) = T_iter(1:z_max)

        end subroutine Implicit_T

      end module Fonction_implicit
