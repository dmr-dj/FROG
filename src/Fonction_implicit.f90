module Fonction_implicit


  use parameter_mod, only : z_num, rho_ice, Gfx, T_freeze,rho_snow_freeze,s_l_max
  use Fonction_temp, only : AppHeatCapacity, ThermalConductivity

  implicit none

  private
  public :: Implicit_T

  contains


    subroutine Implicit_T(T_old,Tu,Tb,dt,dz,n,org_ind,Timp,Kp)

      integer, intent(in) :: org_ind
      real, intent(in) :: dt,Tu,Tb
      real, dimension(:), intent(in) :: T_old, n, dz
      real, dimension(z_num), intent(out) :: Timp, Kp
!dmr [UNUSED]      real, dimension(z_num), intent(out) :: Cp

      real, dimension(z_num) :: pori, porf, Cp_temp
      real, dimension(z_num,z_num) :: MM
      real, dimension(z_num) :: Knows
      real :: m_Gfx, A, B, C, Z1
      integer :: kk, ll
      real, dimension(z_num) :: T_last, Kp_s
!dmr [UNUSED]      real, dimension(z_num) :: T_last, T_new, Kp_s
      real, dimension(z_num) :: T_iter
      real, dimension(z_num) :: DD
      real, dimension(z_num-1) :: DL, DU
!dmr [UNUSED]      integer, dimension(z_num) :: IPIV
      integer :: info_dgesv

      m_Gfx = gfx/1000.0

      T_last(1:z_num) = T_old(1:z_num)

      do kk=1,5

         MM(1:z_num,1:z_num) = 0
         DD(1:z_num)=0
         DL(1:z_num-1)=0
         DU(1:z_num-1)=0

         if (kk==1) then

            T_iter(1:z_num) = T_last(1:z_num)

         else

            T_iter(1:z_num) = 0.5*(T_iter(1:z_num)+T_last(1:z_num))

         end if

         call AppHeatCapacity(z_num,T_iter,T_freeze,n,org_ind,Cp_temp,porf,pori)

         do ll=1,z_num-1

            call ThermalConductivity(ll,n(ll),pori(ll),porf(ll),org_ind,T_iter(ll),Kp(ll))
            Kp(z_num) = 2
            !Kp(ll) = 2
            !Cp_temp(ll) = Cp_temp(ll)

         end do


         Kp_s(1:z_num-1) = (Kp(1:z_num-1)+Kp(2:z_num))*0.5
         !Kp_s(1:z_num-1) = Kp(1:z_num-1)


         do ll=2,z_num-1

            Z1 = T_last(ll)

            A=(dt/((dz(ll-1)+dz(ll))*0.5*dz(ll)) * (Kp_s(ll-1)/Cp_temp(ll)))
            B=(dt/((dz(ll+1)+dz(ll))*0.5*dz(ll)) * (Kp_s(ll)/Cp_temp(ll)))

            !A=(dt/((dz(ll-1))*dz(ll)) * (Kp_s(ll-1)/Cp_temp(ll)))
            !B=(dt/((dz(ll+1))*dz(ll)) * (Kp_s(ll)/Cp_temp(ll)))
            C= 1+A+B
            !write(*,*) 1-A-B , ll

            MM(ll,ll-1) = -A
            MM(ll,ll) = C
            MM(ll,ll+1) = -B
            DL(ll-1) = -A
            DU(ll) = -B
            DD(ll) = C
            Knows(ll) = Z1

         end do

         A=(dt/((dz(z_num-1)+dz(z_num))*0.5*dz(z_num-1)) * (Kp_s(z_num-1)/Cp_temp(z_num-1)))
         C=1.0+A
         Knows(1) = Tu
         Knows(z_num)=Tb
         MM(z_num,z_num)=1
         MM(1,1)=1
         DD(1) = 1
         DD(z_num) = 1+A
         DL(z_num-1) = -A






         !MM(z_num,z_num-1)=-A
         !MM(z_num,z_num)=C
         !MM(1,1) = 1

         !DD(z_num) = C
         !DD(1) = 1

         !Knows(z_num) = T_last(z_num) + ((dt/Cp_temp(z_num))*m_Gfx/dz(z_num-1))

         !call sgesv(z_num,1,MM,z_num,IPIV,Knows,z_num,info_dgesv)
         call sgtsv(z_num,1,DL,DD,DU,Knows,z_num,info_dgesv)
         T_iter(1:z_num) = Knows(1:z_num)

      end do


      Timp(1:z_num) = T_iter(1:z_num)
      !write(*,*) dz(z_num),dz(z_num-1)

    end subroutine Implicit_T


end module Fonction_implicit
