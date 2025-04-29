module vamper_step_mod


  use parameter_mod, only : z_num,TotTime,nb_mon_per_year,YearType,Depth,GridType,PorosityType,T_init,Bool_glacial
  use parameter_mod, only : Bool_Organic,organic_depth,Gfx, T_freeze, EQ_Tr, EQ1_EQ2, Bool_delta,t_fin, alpha
  use parameter_mod, only : Forcage_Month_day,Bool_Swe_Snw,Bool_Model_Snow,Bool_Bessi,s_l_max ! Bool_layer_temp,
  use Fonction_temp, only : AppHeatCapacity, ThermalConductivity, Permafrost_Depth

  use Model_snow, only : snw_average_swe, snw_proc, snw_average_snw, snw_average_snw_tot
  use Fonction_implicit, only : Implicit_snow, Implicit_T

  Implicit none

  private !dmr making sure local variables are local, hence private

  public ::  Vamper_step ! Vamper_init, Lecture_forcing


#include "constant.h"

contains


  subroutine Vamper_step(T_air,swe_f_t,Temp,Tb,Cp,Kp,n,organic_ind,dim_temp,dim_swe,z_num,dz,dt, &  !nb_lines, t_step,
                         porf,pori,D, t_num                                           &    ! ,rho_snow_t,snw_dp_t,T_snw_t, spy,
#if ( CARBON == 0 )
                          )
#else
                  ,end_year , deepSOM_a, deepSOM_s, deepSOM_p, clay                         &   !!  fc,
                  ,diff_k_const, bio_diff_k_const, min_cryoturb_alt, max_cryoturb_alt, zf_soil &
                  ,bioturbation_depth, ALT, altmax_lastyear)

#endif

    use parameter_mod, only: namerun

#if ( CARBON == 1 )
    use Carbon, only : carbon_init, compute_alt, carbon_redistribute, decomposition, cryoturbation
#endif

    integer, intent(inout) ::  organic_ind, dim_swe, dim_temp ! nb_lines, , t_step
    real, intent(in) :: dt, Tb
    integer, intent(in) :: z_num
    integer, intent(in) :: t_num   !! spy,

    real,dimension(z_num),intent(inout) :: dz,n,porf,pori
    real,dimension(z_num),intent(inout) :: Kp,Cp,D
    real,dimension(dim_swe),intent(in) :: swe_f_t
    real,dimension(dim_temp),intent(in) :: T_air
    real, dimension(z_num),intent(inout) :: Temp
#if ( CARBON == 1 )
    integer                              :: compteur_time_step !nb of day or month of the year
    integer                              :: end_year !=0 if not end of year, =1 if end of year
#endif

    integer :: kk, ll
    real :: T_soil, snw_tot, swe_f,Per_depth

    real, dimension(z_num-1) ::  h_n, h_pori, h_porf
    real, dimension(z_num) :: T_old


#if ( CARBON == 1 )
!nb and mbv Carbon cycle
    real                                  :: ALT, altmax_lastyear
    real, dimension(z_num)                :: Temp_positive
    real, dimension(z_num), intent(inout) :: deepSOM_a, deepSOM_s, deepSOM_p!, fc
    real                                  :: clay
    real, intent(in)                      :: diff_k_const, bio_diff_k_const
    real, intent(in)                      :: min_cryoturb_alt, max_cryoturb_alt
    real, dimension(z_num) , intent(in)   :: zf_soil
    real, intent(in)                      :: bioturbation_depth
#endif


!dmr ---
!dmr [NOTA] ---> This bit looks very much like an init, not being part of a step ...
!dmr ---


    write(*,*) "[PRINC] organic_ind: ", organic_ind

    do kk=1,z_num

       Cp(kk) = 1

    end do

    do ll = 1,2

       call AppHeatCapacity(z_num,Temp,T_freeze,n, organic_ind, Cp, porf, pori)         !Calculation of heat capacity of soil

       do kk=1,z_num-1

          h_pori(kk) = (pori(kk) + pori(kk+1))/2
          h_porf(kk) = (porf(kk) + porf(kk+1))/2
          h_n(kk) = (n(kk) + n(kk+1))/2

          call ThermalConductivity(kk,h_n(kk),h_pori(kk),h_porf(kk), organic_ind, Temp(kk), Kp(kk))    !Calculation of thermal condutivity of soil

       end do

    end do

#if ( CARBON == 1 )
!nb and mbv
    end_year=0
    compteur_time_step=0
#endif

!dmr ---
!dmr [NOTA] <--- This bit looks very much like an init, not being part of a step ...
!dmr ---

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
    do ll=1,t_num !boucle temporelle
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( CARBON == 1 )
       !nb and mbv for carbon cycle, is it the end of the year?
       compteur_time_step=compteur_time_step+1
#if ( DAILY == 1 )
       if (compteur_time_step.eq.YearType) then !nomber of days per year
#else
       if (compteur_time_step.eq.nb_mon_per_year) then !nomber of months per year
#endif
           end_year=1
           compteur_time_step=0
       else
           end_year=0
       endif
#endif

![DELSNOW]        snw_old = snw_tot

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr    This is the section that updates the climate forcing, ill-placed
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       T_soil = T_air(mod(ll,dim_temp)+1)
       swe_f  = swe_f_t(mod(ll,dim_swe)+1)
       snw_tot = swe_f_t(mod(ll,dim_swe)+1)

       T_old(1:z_num) = Temp(1:z_num)

       !-------------- Numerical difference routine when there is snow or not --------!

       call Implicit_T(T_old,T_soil,Tb,dt,dz,n,organic_ind,Temp,Kp)

                       ! Returns Per_depth as depth of the 0C isotherm
                       !    in meters, not per cell
       call Permafrost_Depth(Temp,D,Per_depth)


#if ( CARBON == 1 )
       ! nb and mbv carbon cycle call
       ! at the end of each year computes the actve layer thickness, needed for redistribution
       call compute_alt(Temp, Temp_positive, ALT, compteur_time_step, end_year, altmax_lastyear, D)
       !write(*,*) 'ALT', ALT
       ! redistribute carbon from biosphere model
       call carbon_redistribute(Temp, deepSOM_a, deepSOM_s, deepSOM_p, dz, ALT)
       ! computes the decomposition in permafrost as a function of temperature (later : humidity and soil type?)
       !! ICI verifier D pour zi_soil?
       call decomposition(Temp, D, dt, deepSOM_a, deepSOM_s, deepSOM_p, clay)
       ! cryoturbation et bioturbation
       !! ICI chercher zi_soil et zf_soil dans VAMPER ?? D ???
       call cryoturbation(Temp, deepSOM_a, deepSOM_s, deepSOM_p, altmax_lastyear, max_cryoturb_alt, &
            min_cryoturb_alt, D, zf_soil, diff_k_const, bio_diff_k_const, dt,         &
            bioturbation_depth)
#endif



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
    end do ! on the time loop
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

  end subroutine Vamper_step


end module vamper_step_mod
