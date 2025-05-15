module vamper_step_mod


!~   use parameter_mod, only : z_num,TotTime,nb_mon_per_year,YearType,Depth,GridType,PorosityType,T_init,Bool_glacial
!~   use parameter_mod, only : Bool_Organic,organic_depth,Gfx, T_freeze, EQ_Tr, EQ1_EQ2, Bool_delta,t_fin, alpha
!~   use parameter_mod, only : Forcage_Month_day,Bool_Swe_Snw,Bool_Model_Snow,Bool_Bessi,s_l_max ! Bool_layer_temp,
!~   use Fonction_temp, only : AppHeatCapacity, ThermalConductivity, Permafrost_Depth

!~   use Model_snow, only : snw_average_swe, snw_proc, snw_average_snw, snw_average_snw_tot
!~   use Fonction_implicit, only : Implicit_snow, Implicit_T

!~   Implicit none

!~   private !dmr making sure local variables are local, hence private

!~   public ::  Vamper_step ! Vamper_init, Lecture_forcing


!~ #include "constant.h"

!~ contains


!~   subroutine Vamper_step(T_air,swe_f_t,Temp,Tb,Cp,Kp,n,organic_ind,dim_temp,dim_swe,z_num,dz,dt, &  !nb_lines, t_step,
!~                          porf,pori,D, t_num                                           &    ! ,rho_snow_t,snw_dp_t,T_snw_t, spy,
!~ #if ( CARBON == 0 )
!~                           )
!~ #else
!~                   ,end_year , deepSOM_a, deepSOM_s, deepSOM_p, clay                         &   !!  fc,
!~                   ,diff_k_const, bio_diff_k_const, min_cryoturb_alt, max_cryoturb_alt, zf_soil &
!~                   ,bioturbation_depth, ALT, altmax_lastyear)

!~ #endif

!~     use parameter_mod, only: namerun

!~ #if ( CARBON == 1 )
!~     use Carbon, only : carbon_init, compute_alt, carbon_redistribute, decomposition, cryoturbation
!~ #endif

!~     integer, intent(inout) ::  organic_ind, dim_swe, dim_temp ! nb_lines, , t_step
!~     real, intent(in) :: dt, Tb
!~     integer, intent(in) :: z_num
!~     integer, intent(in) :: t_num   !! spy,

!~     real,dimension(z_num),intent(inout) :: dz,n,porf,pori
!~     real,dimension(z_num),intent(inout) :: Kp,Cp,D
!~     real,dimension(dim_swe),intent(in) :: swe_f_t
!~     real,dimension(dim_temp),intent(in) :: T_air
!~     real, dimension(z_num),intent(inout) :: Temp


!~     integer :: kk, ll
!~     real :: T_soil, snw_tot, swe_f,Per_depth

!~     real, dimension(z_num-1) ::  h_n, h_pori, h_porf
!~     real, dimension(z_num) :: T_old


!~ #if ( CARBON == 1 )
!~ !nb and mbv Carbon cycle
!~     real                                  :: ALT, altmax_lastyear
!~     real, dimension(z_num)                :: Temp_positive
!~     real, dimension(z_num), intent(inout) :: deepSOM_a, deepSOM_s, deepSOM_p!, fc
!~     real                                  :: clay
!~     real, intent(in)                      :: diff_k_const, bio_diff_k_const
!~     real, intent(in)                      :: min_cryoturb_alt, max_cryoturb_alt
!~     real, dimension(z_num) , intent(in)   :: zf_soil
!~     real, intent(in)                      :: bioturbation_depth
!~ #endif


!~ !dmr ---
!~ !dmr [NOTA] ---> This bit looks very much like an init, not being part of a step ...
!~ !dmr ---


!~     write(*,*) "[PRINC] organic_ind: ", organic_ind

!~     do kk=1,z_num

!~        Cp(kk) = 1

!~     end do

!~     do ll = 1,2

!~        call AppHeatCapacity(z_num,Temp,T_freeze,n, organic_ind, Cp, porf, pori)         !Calculation of heat capacity of soil

!~        do kk=1,z_num-1

!~           h_pori(kk) = (pori(kk) + pori(kk+1))/2
!~           h_porf(kk) = (porf(kk) + porf(kk+1))/2
!~           h_n(kk) = (n(kk) + n(kk+1))/2

!~           call ThermalConductivity(kk,h_n(kk),h_pori(kk),h_porf(kk), organic_ind, Temp(kk), Kp(kk))    !Calculation of thermal condutivity of soil

!~        end do

!~     end do


!~   end subroutine Vamper_step


end module vamper_step_mod
