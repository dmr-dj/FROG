      module Fonction_implicit


        use parameter_mod, only : z_num, rho_ice, Gfx, T_freeze,rho_snow_freeze,s_l_max
        use Fonction_temp, only : AppHeatCapacity, ThermalConductivity

        implicit none

        private
        public :: Implicit_T

        contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine Implicit_T(T_old,Tu,Tb,dt,dz,n,org_ind,Timp,Kp,z_max)

         integer, intent(in) :: z_max
         integer, intent(in) :: org_ind
         real, intent(in) :: dt,Tu,Tb
         real, dimension(:), intent(in) :: T_old, n, dz
         real, dimension(1:z_max), intent(out) :: Timp, Kp

         real, dimension(1:z_max) :: pori, porf, Cp_temp
         real, dimension(1:z_max,1:z_max) :: MM
         real, dimension(1:z_max) :: Knows
         real :: m_Gfx, A, B, C, Z1
         integer :: kk, ll
         real, dimension(1:z_max) :: T_last, Kp_s
         real, dimension(1:z_max) :: T_iter
         real, dimension(1:z_max) :: DD
         real, dimension(1:z_max-1) :: DL, DU
         integer :: info_dgesv

         m_Gfx = gfx/1000.0

         T_last(1:z_max) = T_old(1:z_max)

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

      !dmr [FUTURE] Only concerned with the soil part for this section

             !dmr Given Temperature and porosity (n), this computes a new Cp value and porf, pori on the vertical
           call AppHeatCapacity(z_num,T_iter,T_freeze,n,org_ind,Cp_temp,porf,pori)

           do ll=1,z_num-1
              !dmr Given the number of layer, thickness, porosity, Temperature, computes the Kp of the layer
             call ThermalConductivity(ll,n(ll),pori(ll),porf(ll),org_ind,T_iter(ll),Kp(ll))
             Kp(z_num) = 2
           end do

           Kp_s(1:z_num-1) = (Kp(1:z_num-1)+Kp(2:z_num))*0.5

      !dmr [FUTURE] Only concerned with the soil part for this section

!~            do ll=2,z_num-1
           do ll=1+1,z_max-1

             Z1 = T_last(ll)

             A=(dt/((dz(ll-1)+dz(ll))*0.5*dz(ll)) * (Kp_s(ll-1)/Cp_temp(ll)))
             B=(dt/((dz(ll+1)+dz(ll))*0.5*dz(ll)) * (Kp_s(ll)/Cp_temp(ll)))

             C= 1+A+B

             MM(ll,ll-1) = -A
             MM(ll,ll) = C
             MM(ll,ll+1) = -B
             DL(ll-1) = -A
             DU(ll) = -B
             DD(ll) = C
             Knows(ll) = Z1

           end do

           A=(dt/((dz(z_max-1)+dz(z_max))*0.5*dz(z_max-1)) * (Kp_s(z_max-1)/Cp_temp(z_max-1)))
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

         end do


         Timp(1:z_max) = T_iter(1:z_max)

        end subroutine Implicit_T

      end module Fonction_implicit
