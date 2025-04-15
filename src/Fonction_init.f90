

module Fonction_init

  use parameter_mod, only : n_organic, Porosity_soil, z_num

  implicit none

  private
  public :: Porosity_init, GeoHeatFlow, Glacial_index

  contains

  !Routine paramettrant la porosité du sol en fonction de la profondeur!

  subroutine Porosity_init(Porosity_Type, Depth_layer, Bool_Organic, organic_depth, n, organic_ind)

    ! Entrées et sorties !

    integer, intent(in) :: Porosity_Type, Bool_Organic
    real, intent(in) ::  organic_depth
    real, dimension(:), intent(in) :: Depth_layer
    integer, intent(out) :: organic_ind
    real, dimension(:), intent(out) :: n

    ! Variables locales !

    real, dimension(:), allocatable :: DepthCalcOrg
    real :: pente
    real :: origine
    integer :: pas_z

!dmr [TBRMD] already allocated    allocate(n(1:z_num))
    allocate(DepthCalcOrg(1:z_num))

    if (Bool_organic == 1) then


       do pas_z = 1,z_num

          DepthCalcOrg(pas_z) = abs(Depth_layer(pas_z) - organic_depth)

       end do

       call indice_minimum(DepthCalcOrg, z_num, organic_ind)

       do pas_z = 1,organic_ind

          n(pas_z) = n_organic

       end do

    else

       n(1) = Porosity_soil
       organic_ind = 1

    end if



    if (Porosity_Type == 1) then


       pente = (Porosity_soil -0.15)/(Depth_layer(organic_ind) - Depth_layer(z_num))
       origine = Porosity_soil - pente*Depth_layer(organic_ind)

       do pas_z = organic_ind, z_num

          n(pas_z) = pente * Depth_layer(pas_z) + origine

       end do

    elseif(Porosity_Type == 2) then

       do pas_z = organic_ind, z_num

          n(pas_z) = Porosity_soil

       end do

    else

       do pas_z = organic_ind, z_num

          n(pas_z) = Porosity_soil*exp(-0.000395*Depth_layer(pas_z))

       end do

    end if

  end subroutine Porosity_init




  subroutine GeoHeatFlow(Gfx, Kp, dz, T0, T)

!~     integer, intent(in) :: z_num
    real, intent(in)                             :: Gfx, T0 ! Gfx = Earth's geothermal heat flux, T0 is modified top temperature (?)
    real, dimension(z_num)  , intent(in)         :: dz ! vertical geometry (thickness of layer)
    real, dimension(z_num-1), intent(in)         :: Kp ! dmr [2024-06-28] : Kp is z_num-1 only ...
    real, dimension(z_num)  , intent(out)        :: T  ! new vertical temperature profile

    integer :: kk ! index

    T(1) = T0

    do kk = 2, z_num
       T(kk) = T(kk-1) + (((Gfx/1000.0)/Kp(kk-1))*dz(kk-1))
    end do



  end subroutine GeoHeatFlow


  subroutine Glacial_index(time_gi, glacial_ind, nb_lines)

    real,dimension(:),allocatable,intent(out) :: time_gi, glacial_ind
    integer,intent(out) :: nb_lines
    integer :: unit_number4,ll
    real :: time_temp, niveau_eau,glacial_ind_temp,kk
    character(len=8) :: char ! [NOTA : BAD NAME] (char is a reserved word in FORTRAN)


    open(newunit=unit_number4,file="Donnee/Glacial_index.txt",status="old",action='read')

    read(unit_number4,*) char, nb_lines

    write(*,*) "[FONC_INIT] nb_lines: ", nb_lines

    allocate(time_gi(1:nb_lines))
    allocate(glacial_ind(1:nb_lines))

    do ll=1,nb_lines-1

       read(unit_number4,*) time_temp,niveau_eau,glacial_ind_temp,kk

       glacial_ind(ll) = glacial_ind_temp
       time_gi(ll) = time_temp

    end do


  end subroutine Glacial_index


  subroutine indice_minimum(tab, taille, ind_min)

    implicit none

    integer, intent(in) :: taille
    real, dimension(:), intent(in) :: tab
    integer, intent(out) :: ind_min

    integer :: kk

    ind_min = 1

    do kk = 2,taille,1

       if(tab(kk)<tab(ind_min)) then

          ind_min = kk

       end if

    end do



  end subroutine indice_minimum


end module Fonction_init



