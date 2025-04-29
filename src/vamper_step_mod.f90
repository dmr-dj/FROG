module vamper_step_mod


  use parameter_mod, only : z_num,TotTime,nb_mon_per_year,YearType,Depth,GridType,PorosityType,T_init,Bool_glacial
  use parameter_mod, only : Bool_Organic,organic_depth,Gfx, T_freeze, EQ_Tr, EQ1_EQ2, Bool_delta,t_fin, alpha
  use parameter_mod, only : Forcage_Month_day,Bool_Swe_Snw,Bool_Model_Snow,Bool_Bessi,s_l_max ! Bool_layer_temp,
  use Fonction_temp, only : AppHeatCapacity, ThermalConductivity, Permafrost_Depth

  use Model_snow, only : snw_average_swe, snw_proc, snw_average_snw, snw_average_snw_tot
  use Fonction_implicit, only : Implicit_snow, Implicit_T

  Implicit none

  private !dmr making sure local variables are local, hence private

  public :: Lecture_forcing, Vamper_step ! Vamper_init,

  integer :: u_n_ml
  integer :: layer_temp23,layer_temp53,layer_temp93,layer_temp143,layer_temp250,layer_temp350,layer_temp550,layer_temp900
  integer :: unit_nb_1,unit_nb_2,unit_nb_3,unit_nb_4,unit_nb_5,unit_nb_6

#include "constant.h"

contains

  subroutine Lecture_forcing(z_num,T_air,swe_f_t,snw_dp_t,rho_snow_t,T_snw,Temp,dim_temp,dim_swe)

    use parameter_mod, only: Tempsolinit, Tempairmonth, Tempairday, Tempsnowmonth ,Tempsnowday

    integer, intent(in) :: z_num
    real, dimension(z_num),intent(inout) :: Temp
    integer, intent(out) :: dim_temp,dim_swe
    real, dimension(:), allocatable, intent(out):: T_air, swe_f_t, snw_dp_t,rho_snow_t,T_snw
    integer :: kk,ii,ll
    real :: ligne

    ll = 0
    dim_swe = 0
    dim_temp = 0

   ! write(*,*) "[Condition] boucle:," ,Forcage_Month_day


    open(newunit=unit_nb_3,file=Tempsolinit,status="old",action='read')
        !WRITE(*,*) "READ "//Tempsolinit

#if ( DAILY == 1 )
    open(newunit=unit_nb_1,file=Tempairday,status="old",action='read')
          open(newunit=unit_nb_2,file=Tempsnowday,status="old",action='read')
          !WRITE(*,*) "READ "//Tempairday
          !WRITE(*,*) "READ "//Tempsnowday
#else
    open(newunit=unit_nb_1,file=Tempairmonth,status="old",action='read')
          open(newunit=unit_nb_2,file=Tempsnowmonth,status="old",action='read')
         ! WRITE(*,*) "READ "//Tempairmonth
          !WRITE(*,*) "READ "//Tempsnowmonth
#endif


    !stop ici

![DELSNOW]    if (Bool_Bessi==1)then
![DELSNOW]
![DELSNOW]       open(newunit=unit_nb_4,file="Donnee/Snow_dp_1998.txt",status="old",action='read')
![DELSNOW]       open(newunit=unit_nb_5,file="Donnee/Rho_snow_1998.txt",status="old",action='read')
![DELSNOW]
![DELSNOW]    end if

    !open(newunit=unit_nb_6,file="Donnee/T_snw.txt",status="old",action='read')

   ! if(EQ1_EQ2==2)then

       do kk=1,z_num
          read(unit_nb_3,*) Temp(kk)
       end do

       !do kk=1,z_num
       !   Temp(kk) = Temp(kk)
       !end do

       close(unit_nb_3)

   ! end if


    do

       read(unit_nb_2,*,iostat=ii) ligne

       if(ii/=0)exit

       dim_swe = dim_swe + 1

    end do

    allocate(swe_f_t(1:dim_swe))
    rewind(unit_nb_2)

    do kk = 1,dim_swe

       read(unit_nb_2,*) swe_f_t(kk)

    end do


    do

       read(unit_nb_1,*,iostat=ii) ligne

       if(ii/=0)exit

       dim_temp = dim_temp + 1

    end do

    allocate(T_air(1:dim_temp))
    rewind(unit_nb_1)

    do kk = 1,dim_temp

       read(unit_nb_1,*) T_air(kk)

    end do


    if (Bool_Bessi==1)then
       allocate(rho_snow_t(1:dim_temp))
       allocate(snw_dp_t(1:dim_temp))
       allocate(T_snw(1:dim_temp))

       do kk = 1,dim_temp

          read(unit_nb_5,*) rho_snow_t(kk)
          read(unit_nb_4,*) snw_dp_t(kk)
          read(unit_nb_6,*) T_snw(kk)

       end do

    else

       allocate(rho_snow_t(1:dim_temp))
       allocate(snw_dp_t(1:dim_temp))
       allocate(T_snw(1:dim_temp))

       do kk = 1,dim_temp

          rho_snow_t(kk) = 0
          snw_dp_t(kk) = 0
          !read(unit_nb_6,*) T_snw(kk)

       end do

    end if


    close(unit_nb_1)
    close(unit_nb_2)


  end subroutine Lecture_forcing


  subroutine Vamper_step(T_air,swe_f_t,Temp,Tb,Cp,Kp,n,organic_ind,glacial_ind,nb_lines,dim_temp,dim_swe,z_num,dz,dt,t_step, &
                         porf,pori,t_deb,rho_snow_t,snw_dp_t,T_snw_t,D, spy, t_num                                           &
#if ( CARBON == 0 )
                          )
#else
                  ,end_year , deepSOM_a, deepSOM_s, deepSOM_p, fc, clay                         &
                  ,diff_k_const, bio_diff_k_const, min_cryoturb_alt, max_cryoturb_alt, zf_soil &
                  ,bioturbation_depth, ALT, altmax_lastyear)

#endif

    use parameter_mod, only: namerun

#if ( CARBON == 1 )
    use Carbon, only : carbon_init, compute_alt, carbon_redistribute, decomposition, cryoturbation
#endif

    integer, intent(inout) ::  organic_ind, nb_lines, dim_swe, dim_temp, t_step,t_deb
    real, intent(in) :: dt, Tb
    integer, intent(in) :: z_num
    integer, intent(in) :: spy, t_num

    real,dimension(z_num),intent(inout) :: dz,n,porf,pori
    real,dimension(z_num),intent(inout) :: Kp,Cp,D
    real,dimension(dim_swe),intent(in) :: swe_f_t
    real,dimension(dim_temp),intent(in) :: T_air,snw_dp_t,rho_snow_t,T_snw_t
    real,dimension(nb_lines),intent(in) :: glacial_ind
    real, dimension(z_num),intent(inout) :: Temp
#if ( CARBON == 1 )
    integer                              :: compteur_time_step !nb of day or month of the year
    integer                              :: end_year !=0 if not end of year, =1 if end of year
#endif

    integer :: kk, ll, indice_tab, s_l_t,ind_snw
    real :: T_soil, T_glacial, swe_tot,snw_tot,rho_snow,swe_f,frac_snw,k_s,Cp_snow,snw_old,dz_snow,Per_depth
    real, dimension(s_l_max) :: Tsnw
    real, dimension(z_num-1) ::  h_n, h_pori, h_porf
    real, dimension(z_num) :: T_old
    real, dimension(:),allocatable :: T_layer23,T_layer53,T_layer93,T_layer143,T_layer250,T_layer350,T_layer550,T_layer900

! dmr&mbv --- Added for cleaner output at given fixed levels
    integer, dimension(10)            :: indx_min
    real   , dimension(10), parameter :: fixed_levs = [ 0.0, 0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 500.0, 750.0, 1000.0 ]
    integer :: zzz , index_claque , tempmens , ALLTT

#if ( CARBON == 1 )
!nb and mbv Carbon cycle
    real                                  :: ALT, altmax_lastyear
    real, dimension(z_num)                :: Temp_positive
    real, dimension(z_num), intent(inout) :: deepSOM_a, deepSOM_s, deepSOM_p, fc
    real                                  :: clay
    real, intent(in)                      :: diff_k_const, bio_diff_k_const
    real, intent(in)                      :: min_cryoturb_alt, max_cryoturb_alt
    real, dimension(z_num) , intent(in)   :: zf_soil
    real, intent(in)                      :: bioturbation_depth
#endif

!dmr [TBRMD] Useless in the context of more than one column ... they would overwrite each other !
!~      if (Bool_layer_temp==1)then

!~         open(newunit=u_n_ml,file="Results/"//trim(namerun)//".txt", status="replace",action='write')
!~         !open(newunit=u_n_ml,file="Results/toto.txt", status="replace",action='write')

!~        DO zzz = LBOUND(indx_min,DIM=1), UBOUND(indx_min,DIM=1)
!~           indx_min(zzz) = minloc(abs(D-fixed_levs(zzz)), DIM=1)
!~        ENDDO

!~     end if

!dmr ---
!dmr [NOTA] This bit looks very much like an init, not being part of a step ... --->
!dmr ---

![DELSNOW]     swe_tot = 0
![DELSNOW]     snw_tot = 0
![DELSNOW]     snw_old = 0
![DELSNOW]     ind_snw = 0
![DELSNOW]     s_l_t = 1

    write(*,*) "[PRINC] organic_ind: ", organic_ind

![DELSNOW]     do kk =1,s_l_max
![DELSNOW]
![DELSNOW]        Tsnw(kk) = -4
![DELSNOW]
![DELSNOW]     end do

!dmr [CHECK] The following code seems to be modifying the porosity where there is snow ... even if snow not used !
!dmr         Commenting it out, but will change the results possibly
![DELSNOW]     dz_snow= 0.01
![DELSNOW]     D(1) = dz_snow*s_l_max
![DELSNOW]     do kk = 1,s_l_max-1
![DELSNOW]        D(kk+1) = D(kk) - dz_snow
![DELSNOW]     end do

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

![DELSNOW]       if (Bool_Bessi==1)then
![DELSNOW]          rho_snow = rho_snow_t(mod(ll,dim_temp)+1)
![DELSNOW]          snw_tot = snw_dp_t(mod(ll,dim_temp)+1)
![DELSNOW]       end if

!~        if (Bool_glacial==1)then

!~           indice_tab = nb_lines-floor((-(ll/real(spy))+t_deb)/100.0)
!~           T_glacial=alpha*(glacial_ind(indice_tab-1)+(100-mod((-(ll/real(spy))+t_deb),100.0))*(glacial_ind(indice_tab)- &
!~ glacial_ind(indice_tab-1))/100.0)
!~           T_soil = (T_glacial+T_soil)
!~        end if

![DELSNOW]       if (snw_tot > 0.000001 .or. swe_f > 0.000001  ) then

![DELSNOW]          if (Bool_Bessi==0) then
![DELSNOW]             call snw_average_swe(swe_f, swe_tot, snw_tot, rho_snow)
![DELSNOW]          end if

![DELSNOW]          if (abs(swe_tot - swe_f )< 0.000001 .and. Bool_Bessi==0) then

![DELSNOW]             Tsnw = T_soil
![DELSNOW]             K_s = 0.07
![DELSNOW]             frac_snw = 1

![DELSNOW]          end if

![DELSNOW]          if (abs(snw_old )< 0.000001) then

![DELSNOW]             Tsnw(1) = T_soil
![DELSNOW]             !Tsnw(2) = (
![DELSNOW]             K_s = 0.07
![DELSNOW]             frac_snw = 1

![DELSNOW]          end if
![DELSNOW]       end if

![dmr] This states that in the case that Tsoil is positive, snow instantaneously melts????
![DELSNOW]       if (T_soil>0.0) then
![DELSNOW]
![DELSNOW]          snw_tot = 0
![DELSNOW]
![DELSNOW]       end if

       T_old(1:z_num) = Temp(1:z_num)

       !-------------- Numerical difference routine when there is snow or not --------!

![DELSNOW]       if (snw_tot > 0.00001) then

![DELSNOW]          s_l_t = 1

![DELSNOW]          !call Implicit_snow(snw_tot,rho_snow,Tsnw,T_old,T_soil,dt,dz,n,organic_ind,Temp,Cp,Kp,Cp_snow,s_l_t)


![DELSNOW]          !T_soil = T_snw_t(mod(ll,dim_temp)+1)

![DELSNOW]          call Implicit_T(T_old,T_soil,Tb,dt,dz,n,organic_ind,Temp,Kp)

![DELSNOW]          if (Bool_Bessi==0) then
![DELSNOW]             call snw_proc(Tsnw(1), snw_tot, swe_tot, frac_snw, Cp_snow, rho_snow, dt)
![DELSNOW]          end if

![DELSNOW]       else

![DELSNOW]           swe_tot = 0.0
![DELSNOW]           snw_tot = 0.0

![DELSNOW]           do kk =1,s_l_max
![DELSNOW]
![DELSNOW]              Tsnw(kk) = 0
![DELSNOW]
![DELSNOW]           end do

          call Implicit_T(T_old,T_soil,Tb,dt,dz,n,organic_ind,Temp,Kp)

![DELSNOW]       end if

!~        if (Bool_layer_temp==1)then
!~          if (ll>t_num-7300)then
!~             write(u_n_ml,'(10F12.3)') ( Temp(indx_min(zzz)) &
!~                  , zzz=LBOUND(indx_min,DIM=1),UBOUND(indx_min,DIM=1))
!~           end if
!~        end if

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
       write(ALLTT,*) ALT
#endif

       !MBV test write temp par mois
      ! if (mod(ll,30).eq.0) then
         write(tempmens,*) (Temp(index_claque),                                &
                          index_claque=LBOUND(Temp,dim=1), UBOUND(Temp, dim=1))
       !endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
    end do ! on the time loop
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

    t_deb = t_deb - floor(t_step/365.0)

!~     if (Bool_glacial.eq.1) then !dmr indice_tab is undefined if Bool_glacial is not 1
!~     !  write(*,*) "[PRINC] indice_tab,t_num,organic_ind: ", indice_tab,t_num,organic_ind
!~     endif

    close(u_n_ml)
    close(ALLTT)

    close(tempmens)


  end subroutine Vamper_step


  subroutine vamp_step_dt


  end subroutine vamp_step_dt


end module vamper_step_mod
