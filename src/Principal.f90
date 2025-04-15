module Principal


  use parameter_mod, only : z_num,TotTime,nb_mon_per_year,YearType,Depth,GridType,PorosityType,T_init,Bool_glacial 
  use parameter_mod, only : Bool_Organic,organic_depth,Gfx, T_freeze, EQ_Tr, EQ1_EQ2, Bool_delta,t_fin, alpha
  use parameter_mod, only : Bool_layer_temp,Forcage_Month_day,Bool_Swe_Snw,Bool_Model_Snow,Bool_Bessi,s_l_max
  use Fonction_temp, only : AppHeatCapacity, ThermalConductivity, Permafrost_Depth
  use Fonction_init, only : Porosity_init, GeoHeatFlow, Glacial_index
  use Para_fonctions, only : t_disc, z_disc
  use Model_snow, only : snw_average_swe, snw_proc, snw_average_snw, snw_average_snw_tot
  use Fonction_implicit, only : Implicit_snow, Implicit_T
 
  Implicit none
  
  private !dmr making sure local variables are local, hence private
  
  public :: Vamper_init, Lecture_forcing, Vamper_step
  
!dmr&mbv  integer :: u_n_23,u_n_53,u_n_93,u_n_143,u_n_250,u_n_350,u_n_550,u_n_900 [TOREMOVE]
  integer :: u_n_ml
  integer :: layer_temp23,layer_temp53,layer_temp93,layer_temp143,layer_temp250,layer_temp350,layer_temp550,layer_temp900
  integer :: unit_nb_1,unit_nb_2,unit_nb_3,unit_nb_4,unit_nb_5,unit_nb_6

# include "constant.h"

contains

  subroutine Vamper_init(dz,D,Temp,time_gi,glacial_ind,nb_lines,Kp,Cp,n,organic_ind,Tb)

    use parameter_mod, only: z_num

!~     integer, intent(in) :: z_num                                     [TBRMD]
    
    real, dimension(z_num),intent(in)          :: dz,D ! geometry of the layers
    real, dimension(z_num),intent(out)         :: Cp   ! Cp constant that is not constant but computed
    real, dimension(z_num-1),intent(out)       :: Kp   ! Kp constant ? [TO_BE_CLARIFIED]
    
    real, dimension(:),allocatable,intent(out) :: time_gi,glacial_ind
    real, dimension(:),intent(out)             :: n,Temp
    
    real, intent(out) :: Tb
    integer, intent(out) :: nb_lines
    integer, intent(out) :: organic_ind ! depth of the organic layer

    real, dimension(z_num-1) ::  h_n, h_pori, h_porf
    real, dimension(z_num) :: porf,pori
    integer :: kk,ll
    
    !dmr [2024-06-28] CALCULATION OF POROSITY in the vertical
    !dmr   
    !dmr [NOTA]  The organic layer will be spatial in space, so "n" should be spatial in space in the end, also organic_ind and organic_depth
    !dmr         Input could be a map of organic_depth
    !dmr 
    !dmr intent(in)                PorosityType == 1 or 2 [TO_BE_CLARIFIED]
    !dmr intent(in)                Bool_Organic == 1 use organic layer, else not
    !dmr intent(in)                organic_depth = depth of organic layer (I guess), in meters
    !dmr intent(in) (z_num)        D depth of the layer considered in meters
    !dmr intent(out) (allocatable) n = porosity of each layer in the vertical
    !dmr intent(out)               organic_ind, index in vertical of the end of organic layer (1:organic_ind)
    call Porosity_init(PorosityType, D, Bool_Organic, organic_depth, n, organic_ind )  !CALCULATION OF POROSITY

    !dmr [NOTA] The code below seems to be used to increase the value of n (Porosity ?) if Depth is greater than 1.4 meter ...
    !write(*,*) n

    !do kk=1,z_num
       !if(D(kk)>1.4)then
       !   n(kk) = n(kk) + 0.5
       !end if
    !end do


    if (Bool_glacial == 1)then

       call Glacial_index(time_gi,glacial_ind,nb_lines)    ! GLACIAL INDEX FOR THE FORCING OF TEMPERATURE

    else

       allocate(glacial_ind(1:1))
       glacial_ind(1) = 0

    end if
                                
    !dmr [NOTA] For now it seems that Kp is constant, are there reasons to have it spatially or vertically variable?
    do kk=1,UBOUND(Kp,DIM=1) !dmr Kp is size z_num-1   
       Kp(kk)=2
    end do
          
    !dmr [2024-06-28] CALCULATION OF HeatFlow
    !dmr   

!dmr [TBRMD] already allocated    allocate(Temp(z_num))
          
    call GeoHeatFlow(Gfx, Kp, dz, T_init, Temp)            

    Tb = Temp(z_num)                         ! Lower boundary condition 
    
    write(*,*) "[PRinc ] Tb, Temp" , Tb ,Temp

    do ll = 1,2 !dmr WhatIs ll ? 
       
       call AppHeatCapacity(z_num,Temp,T_freeze,n, organic_ind, Cp, porf, pori)         !Calculation of heat capacity of soil

       do kk=1,z_num-1
          
          h_pori(kk) = (pori(kk) + pori(kk+1))/2
          h_porf(kk) = (porf(kk) + porf(kk+1))/2
          h_n(kk) = (n(kk) + n(kk+1))/2
          call ThermalConductivity(kk,h_n(kk),h_pori(kk),h_porf(kk), organic_ind, Temp(kk), Kp(kk))    !Calculation of thermal condutivity of soil
          
       end do
       
    end do


  end subroutine Vamper_init


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
    
    if (Bool_Bessi==1)then

       open(newunit=unit_nb_4,file="Donnee/Snow_dp_1998.txt",status="old",action='read')
       open(newunit=unit_nb_5,file="Donnee/Rho_snow_1998.txt",status="old",action='read')

    end if

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

     if (Bool_layer_temp==1)then

        open(newunit=u_n_ml,file="Results/"//trim(namerun)//".txt", status="replace",action='write')
        !open(newunit=u_n_ml,file="Results/toto.txt", status="replace",action='write')

       DO zzz = LBOUND(indx_min,DIM=1), UBOUND(indx_min,DIM=1)
          indx_min(zzz) = minloc(abs(D-fixed_levs(zzz)), DIM=1)
       ENDDO

    end if


    swe_tot = 0
    snw_tot = 0
    snw_old = 0
    ind_snw = 0
    s_l_t = 1
    !write(*,*) "[PRINC] organic_ind: ", organic_ind
    
    do kk =1,s_l_max

       Tsnw(kk) = -4

    end do

    dz_snow= 0.01
    D(1) = dz_snow*s_l_max
    do kk = 1,s_l_max-1
       D(kk+1) = D(kk) - dz_snow
    end do
    
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

! [TOREMOVE]
! dmr&mbv --- already defined in t_disc
!# if Daily == 0
!
!       spy = 12
!       t_num = t_step*12
!
!#else 
!
!       spy = 365
!       ! t_num = t_step 
!       !t_num = 2
!       t_num = 365*30
!
!#endif
! [TOREMOVE]
     ! ICI 
   ! open(newunit=tempmens, file="Results/testmoismens.txt", status="replace")
   ! open(newunit=ALLTT, file="Results/testaltmens.txt", status="replace")
    
    !write(*,*) n

#if ( CARBON == 1 )
!nb and mbv
    end_year=0
    compteur_time_step=0
#endif

!boucle temporelle
    do ll=1,t_num

       ! write(*,*) "Time step:: ", ll," / ", t_num

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


       snw_old = snw_tot

       T_soil = T_air(mod(ll,dim_temp)+1)
       swe_f  = swe_f_t(mod(ll,dim_swe)+1)
       snw_tot = swe_f_t(mod(ll,dim_swe)+1)
       !snw_tot = 0

       if (Bool_Bessi==1)then
          rho_snow = rho_snow_t(mod(ll,dim_temp)+1)
          snw_tot = snw_dp_t(mod(ll,dim_temp)+1)
       end if
       !write(*,*) snw_tot,T_soil,rho_snow
       if (Bool_glacial==1)then

          indice_tab = nb_lines-floor((-(ll/real(spy))+t_deb)/100.0)
          T_glacial=alpha*(glacial_ind(indice_tab-1)+(100-mod((-(ll/real(spy))+t_deb),100.0))*(glacial_ind(indice_tab)- &
glacial_ind(indice_tab-1))/100.0)
          T_soil = (T_glacial+T_soil)
          !write(*,*) indice_tab,T_soil,glacial_ind(indice_tab-1),glacial_ind(indice_tab),T_glacial,mod((-(ll/spy)+t_deb),100.0)

       end if

       if (snw_tot > 0.000001 .or. swe_f > 0.000001  ) then

          if (Bool_Bessi==0) then
             call snw_average_swe(swe_f, swe_tot, snw_tot, rho_snow)
          end if
          
          if (abs(swe_tot - swe_f )< 0.000001 .and. Bool_Bessi==0) then

             Tsnw = T_soil
             K_s = 0.07
             frac_snw = 1

          end if
          
          if (abs(snw_old )< 0.000001) then

             Tsnw(1) = T_soil
             !Tsnw(2) = (
             K_s = 0.07
             frac_snw = 1

          end if
       end if
       
       if (T_soil>0.0) then

          snw_tot = 0

       end if

       T_old(1:z_num) = Temp(1:z_num)
       !write(*,*) Temp

       !-------------- Numerical difference routine when there is snow or not --------!

       if (snw_tot > 0.00001) then
          
          s_l_t = 1
          
          !call Implicit_snow(snw_tot,rho_snow,Tsnw,T_old,T_soil,dt,dz,n,organic_ind,Temp,Cp,Kp,Cp_snow,s_l_t)
       

          !T_soil = T_snw_t(mod(ll,dim_temp)+1)

          call Implicit_T(T_old,T_soil,Tb,dt,dz,n,organic_ind,Temp,Kp) 
!dmr [UNUSED]          call Implicit_T(T_old,T_soil,Tb,dt,dz,n,organic_ind,Temp,Cp,Kp) 


          if (Bool_Bessi==0) then
             call snw_proc(Tsnw(1), snw_tot, swe_tot, frac_snw, Cp_snow, rho_snow, dt)
!dmr [UNUSED]             call snw_proc(Tsnw(1),Temp(1), snw_tot, swe_tot, frac_snw, Cp_snow, rho_snow, dt)
          end if
          
       else

          swe_tot = 0.0
          snw_tot = 0.0
          
          do kk =1,s_l_max

             Tsnw(kk) = 0

          end do

          !T_soil = T_snw_t(mod(ll,dim_temp)+1)

          call Implicit_T(T_old,T_soil,Tb,dt,dz,n,organic_ind,Temp,Kp) 
!dmr [UNUSED]          call Implicit_T(T_old,T_soil,Tb,dt,dz,n,organic_ind,Temp,Cp,Kp)           
          

       end if
       
       !write(*,*) Kp(5)*dt/(Cp(5)*dz(5)*dz(5)) 

       if (Bool_layer_temp==1)then
         if (ll>t_num-7300)then
            write(u_n_ml,'(10F12.3)') ( Temp(indx_min(zzz)) &
                 , zzz=LBOUND(indx_min,DIM=1),UBOUND(indx_min,DIM=1))
          end if
       end if

       call Permafrost_Depth(Temp,D,Per_depth)
       ! write(*,*) "[PRINC] Per_depth: ", Per_depth
       
       !write(*,*) Kp

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
    end do

    t_deb = t_deb - floor(t_step/365.0)

    if (Bool_glacial.eq.1) then !dmr indice_tab is undefined if Bool_glacial is not 1
    !  write(*,*) "[PRINC] indice_tab,t_num,organic_ind: ", indice_tab,t_num,organic_ind
    endif
   
    close(u_n_ml)
    close(ALLTT) 

    close(tempmens)

  
  end subroutine Vamper_step

  
  
end module Principal
