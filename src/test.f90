

program test_fonctions

#include "constant.h"

  use parameter_mod, only: z_num, gridNoMax, dt, t_num, D, dz

  !dmr [2024-06-28] Functions used in the main
  use vamper_step_mod, only :  Vamper_step ! Vamper_init,Lecture_forcing,

  use Fonction_temp, only : AppHeatCapacity, ThermalConductivity
!~   use Fonction_init, only : GeoHeatFlow ! Porosity_init, , Glacial_index


  ! use Para_fonctions, only : z_disc ! t_disc,

  use Model_snow, only : snw_average_swe, snw_proc, snw_average_snw, snw_average_snw_tot

#if ( CARBON == 1 )
  use parameter_mod, only: alt, altmax_lastyear, bio_diff_k_const, diff_k_const, max_cryoturb_alt, min_cryoturb_alt, zf_soil
  use spatialvars_mod, only: deepSOM_a, deepSOM_s, deepSOM_p,  clay_SV ! fc_SV,
  use Carbon, only : carbon_first_init, carbon_init
#endif

  use main_lib_VAMPER, only: INITIALIZE_VAMP

  use spatialvars_mod, only: Temp, Kp, Cp, n, porf, pori

  implicit none

  integer :: kk, ll,organic_ind,dim_temp,dim_swe,t_step ! ,t_deb ! ,t_num, nb_lines, spy,
!~   real,dimension(:),allocatable :: time_gi, glacial_ind ! dmr glacial indexes, to be checked with Amaury
  real :: Tb ! ,dt moved to parameter_mod

  real, dimension(:),allocatable::     &
                                       ! FORCING VARIABLES
                                        T_air     & !dmr [SPAT_VAR], surface air temperature // forcing
                                       ,swe_f_t   & !dmr [SPAT_VAR]  SWE forcing, allocated to (1:dim_swe) and values read in external text file (unit_nb_2)
                                       ,snow_dp_t & !dmr [SPAT_VAR] snow depth over time forcing ???
                                       ,rho_snow_t& !dmr [SPAT_VAR] density of snow over time forcing ???
                                       ,T_snw_t     !dmr [SPAT_VAR] temperature of snow over time forcing ???

#if ( CARBON == 1 )

!~     real                              :: clay
    REAL                              :: bioturbation_depth
    integer                           :: end_year
#endif

  integer :: gridNo = 1 ! index for spatial loops
  logical :: well_done = .FALSE.


  well_done = INITIALIZE_VAMP()
  if (well_done) then
    WRITE(*,*) "VAMPER INITIALIZATION COMPLETE"
  endif

  READ(*,*)

!~   t_deb = 0
  kk=1
  ll=1

  do gridNo = 1, gridNoMax



!~   !dmr [2024-06-28] [ADDING COMMENTS]
!~   !dmr This subroutine will be reading the forcing from external files. There is a suite of options as to how the forcing is done.
!~   !dmr Basically, it needs to fill in:
!~   !dmr
!~   !dmr

!~   !dmr intent(inout) (z_num)     Temp         Temperature of the soil for each layer, read in external file unit_nb_3, one timestep (initial?)

!~   !dmr intent(out)               dim_swe      Length of the SWE forcing, calculated from the length of the input text file
!~   !dmr intent(out) (allocatable) swe_f_t      SWE forcing, allocated to (1:dim_swe) and values read in external text file (unit_nb_2)

!~   !dmr intent(out)               dim_temp     Length of the temp forcing, calculated from the length of the input text file
!~   !dmr intent(out) (allocatable) T_air        temperature of the air forcing, allocated to (1:dim_temp) and read in external file (unit_nb_1)


!~   !dmr Those three are allocated to (1:dim_temp) and if BESSI, read from external files.

!~   !dmr intent(out) (allocatable) snow_dp_t     if BESSI, unit_nb_4 else set to zero
!~   !dmr intent(out) (allocatable) rho_snow_t    if BESSI, unit_nb_5 else set to zero
!~   !dmr intent(out) (allocatable) T_snw_t       if BESSI, unit_nb_6 else ... commented read, nothing done [NOTA: UNINITIALIZED]

!~   call Lecture_forcing(z_num,T_air,swe_f_t,snow_dp_t,rho_snow_t,T_snw_t,Temp(:,gridNo),dim_temp,dim_swe)

  write(*,*) "[MAIN] D: ", D
  write(*,*) "[MAIN] t_num:", t_num
  write(*,*) "[MAIN] forcing: ", dim_temp, dim_swe
  write(*,*) "[Prof]", dz
!~   write(*,*) "[MAIN] 1|Temp: ",Temp

  t_step = dim_temp

  READ(*,*)

  !do kk = 1,800


  !dmr [2024-06-28] [ADDING COMMENTS]
  !dmr Main stepping of the VAMPER model
  !dmr

  !dmr intent(inout) (z_num)     dz
  !dmr intent(inout) (z_num)     n
  !dmr intent(inout) (z_num)     porf
  !dmr intent(inout) (z_num)     pori
  !dmr intent(inout) (z_num)     Kp            / Conductivité thermique, dépend de la température
  !dmr intent(inout) (z_num)     Cp
  !dmr intent(inout) (z_num)     D
  !dmr intent(inout) (z_num)     Temp          / Amaury dixit: température à l'initialisation puis température du sol calculée

  !dmr intent(in)    (dim_swe)   swe_f_t       / SWE forcing data
  !dmr intent(in)    (dim_temp)  T_air         / Air temperature forcing data

  !dmr intent(in)    (dim_temp)  snw_dp_t      / SNOW depth forcing
  !dmr intent(in)    (dim_temp)  rho_snow_t    / DENSITY of snow forcing
  !dmr intent(in)    (dim_temp)  T_snw_t       / TEMP of snow forcing
  !dmr intent(in)    (nb_lines)  glacial_ind   / glacial index modifier

  call Vamper_step(T_air,swe_f_t,Temp(:,gridNo),Tb,Cp(:,gridNo),Kp(:,gridNo),n(:,gridNo),organic_ind &
                  ,dim_temp,dim_swe,z_num,dz,dt                         & !nb_lines, ,t_step
                  ,porf(:,gridNo),pori(:,gridNo), D, t_num    &    ! rho_snow_t,snow_dp_t,T_snw_t, spy,
#if ( CARBON == 0 )
                   )
#else

                  ,end_year , deepSOM_a, deepSOM_s, deepSOM_p, clay_SV(gridNo)                &   !! fc_SV,
                  ,diff_k_const, bio_diff_k_const, min_cryoturb_alt, max_cryoturb_alt, zf_soil       &
                  ,bioturbation_depth, ALT, altmax_lastyear)
#endif

  write(*,*) "[MAIN] 2|Temp: ",Temp(:,gridNo)

  deallocate(T_air)
  deallocate(swe_f_t)
  deallocate(snow_dp_t)
  deallocate(rho_snow_t)
  deallocate(T_snw_t)

  enddo ! loop on gridNo

  write(*,*) "ok"

 end program test_fonctions
