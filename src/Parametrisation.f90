module Parametrisation


  implicit none

  public !dmr here can be fully public, these are only parameter constants


  integer, parameter  :: str_len = 256


  character(len=str_len) :: namerun
  integer :: TotTime         !temps total en année
  integer :: nb_day_per_month           !
  integer :: nb_mon_per_year
  integer :: t_fin
  integer :: YearType         !nombre de jour par an
  integer :: z_num            !nombre de couches étudiée
  integer :: GridType           !(1) Log-generated, (2) Linear-generated
  integer :: PorosityType       !(1) linéaire, (autre) exponentiellement décroissante en fonction de la profondeur
  integer :: Bool_Snow          ! forçage en neige ou non (1 ou 0)
  integer :: Bool_Organic       ! prise en compte de la couche organique ou non (1 ou 0)
  integer :: EQ_Tr             ! Equilibrum run (0) or Transient run (1) -> using different forcing Temperature and snow
  integer :: EQ1_EQ2           ! EQ1(1) initial temperature calculated with the Geothermal heat flux. EQ2 initial temperature read in a file .txt
  integer :: Bool_delta        !
  integer :: Bool_glacial          ! Using glacial index to modify air temperature
  integer :: Bool_layer_temp       ! Creation of .txt with the temperature of the soil at different layer
  integer :: Forcage_Month_day     ! (1) Daily or (0) monthly forcing
  integer :: Bool_Swe_Snw          ! (1) Snow forcing, (0) Swe forcing
  integer :: Bool_Model_Snow       ! (1) Usinsg snow model to find snow_depth, (0) Forcing with snow_depth
  integer :: Bool_Bessi
  integer :: Bool_geometric

  real :: Depth              !profondeur de la modélisation
  real :: T_init                !température initiale a la surface
  real :: T_freeze             !température où l'eau est considérée comme gelée
  real :: freezing_range          ! Temperature at wich the snow start to melt
  real :: Gfx                  ! flux géothermique de la terre (a modifier peut être)
  real :: Porosity_soil       ! porosité du sol
  real :: organic_depth     ! profondeur de la couche organique
  real :: n_organic           ! porosité de la couche organique
  real :: n_soil_bot         ! porosité en bas de la couche de sol
  real :: alpha

  !       DENSITÉ DE DIFFÉRENTES MATIÈRES (en kg/m³) !

  real :: rho_snow_freeze            ! densité de la neige
  real :: rho_water          ! densité de l'eau
  real :: rho_ice             ! densité de la glace
  real :: rho_organic        ! densité de la matière organique
  real :: rho_soil           ! densité du sol
  real :: rho_snow_fresh            ! densité de la neige fraiche

  !      Capacité thermique massique  (en J/(K*kg))    !

  real :: C_water          ! Capacité thermique massique de l'eau
  real :: C_ice            ! Capacité thermique massique de la glace
  real :: C_organic        ! Capacité thermique massique de la matière organique
  real :: C_dry_soil        ! Capacité thermique massique du sol


  !      Conductivité thermique (en W/(m*K))          !

  real :: K_other_minerals    ! Conductivité thermique des autres minéraux
  real :: K_quartz            ! conductivité thermique du quartz
  real :: K_organic        ! conductivité thermique de la matière organique
  real :: K_ice            ! conductivité thermique de la glace
  real :: K_fluids         ! conductivité thermique des fluides

  real :: q_quartz            ! pourcentage de quartz dans le sol

  real :: gravity          ! accéleration gravitationnelle
  real :: Latent_heat    ! en J/kg

  integer :: s_l_max      ! nombre de couche de neige (marche que avec 1)

  !dmr spatialisation
  integer :: gridNoMax    !dmr spatial index, no assumption of the spatial arrangement

  character(len=str_len) :: Tempsolinit
  character(len=str_len) :: Tempairmonth
  character(len=str_len) :: Tempairday
  character(len=str_len) :: Tempsnowmonth
  character(len=str_len) :: Tempsnowday

contains

  subroutine lecture_namelist

    implicit none

!déclaration des variables

    INTEGER :: rc,fu, stderr, fo
    CHARACTER(len= 19)  :: file_path ="vamper_namelist.nml", file_path_chck = "vamper_namelist.chk"
    CHARACTER(len=1), PARAMETER :: tab=char(9)
    CHARACTER(len=3), PARAMETER :: spa=""//tab//tab//tab

    NAMELIST /Param/ namerun,TotTime,nb_day_per_month,nb_mon_per_year, t_fin,YearType,Bool_glacial,alpha,PorosityType, &
                     Bool_Organic,Porosity_soil,organic_depth,n_organic,n_soil_bot, q_quartz,Gfx,Bool_Snow,     &
                     Bool_Swe_Snw,Bool_Model_Snow,Bool_Bessi,s_l_max,z_num,GridType,gridNoMax,Bool_layer_temp,  &
                     Depth,T_init,Bool_delta,Bool_geometric, EQ_Tr, EQ1_EQ2

    NAMELIST /Physique/ rho_snow_freeze,rho_water,rho_ice,rho_organic,rho_soil,rho_snow_fresh,C_water,C_ice, &
            C_organic,C_dry_soil,K_other_minerals,K_quartz,K_organic,K_ice,K_fluids,T_freeze,freezing_range,&
            gravity,Latent_heat

    NAMELIST /Tempdata/ Tempsolinit, Tempairmonth, Tempairday, Tempsnowmonth, Tempsnowday

    INQUIRE (file=file_path, iostat=rc)

    IF (rc /= 0) THEN
                WRITE (stderr, '("Error: input file ", a, " does not exist")') file_path
                STOP
    ENDIF

    ! Open and read Namelist file.
    OPEN (action='read', file=file_path, iostat=rc, newunit=fu)

    IF (rc /= 0) WRITE (stderr, '("Error: invalid Namelist format")')
    READ (nml=Param, iostat=rc, unit=fu)

    IF (rc /= 0) WRITE (stderr, '("Error: invalid Namelist format")')
    READ (nml=Physique, iostat=rc, unit=fu)

    IF (rc /= 0) WRITE (stderr, '("Error: invalid Namelist format")')
    READ (nml=Tempdata, iostat=rc, unit=fu)

    IF (rc /= 0) WRITE (stderr, '("Error: invalid Namelist format")')

    CLOSE (fu)

    ! Write out a file to check the parameters that have been read

    OPEN (action='write', file=file_path_chck, iostat=rc, newunit=fo)

    write(fo,*) "=== PARAMETERS CHECK ==="

    write(fo,*) "namerun"//spa, namerun
    write(fo,*) "TotTime"//spa, TotTime         !temps total en année
    write(fo,*) "nb_day_per_month"//spa, nb_day_per_month           !nombre de jour entre chaque pas de temps
    write(fo,*) "t_fin"//spa, t_fin
    write(fo,*) "YearType"//spa, YearType         !nombre de jour par an
    write(fo,*) "z_num"//spa, z_num            !nombre de couches étudiée
    write(fo,*) "GridType"//spa, GridType           !(1) Log-generated, (2) Linear-generated
    write(fo,*) "PorosityType"//spa, PorosityType       !(1) linéaire, (autre) exponentiellement décroissante en fonction de la profondeur
    write(fo,*) "Bool_Snow"//spa, Bool_Snow          ! forçage en neige ou non (1 ou 0)
    write(fo,*) "Bool_Organic"//spa, Bool_Organic       ! prise en compte de la couche organique ou non (1 ou 0)
    write(fo,*) "EQ_Tr"//spa, EQ_Tr             ! Equilibrum run (0) or Transient run (1) -> using different forcing Temperature and snow
    write(fo,*) "EQ1_EQ2"//spa, EQ1_EQ2           ! EQ1(1) initial temperature calculated with the Geothermal heat flux. EQ2 initial temperature read in a file .txt
    write(fo,*) "Bool_delta"//spa, Bool_delta        !
    write(fo,*) "Bool_glacial"//spa, Bool_glacial          ! Using glacial index to modify air temperature
    write(fo,*) "Bool_layer_temp"//spa, Bool_layer_temp       ! Creation of .txt with the temperature of the soil at different layer
    write(fo,*) "Bool_Swe_Snw"//spa, Bool_Swe_Snw          ! (1) Snow forcing, (0) Swe forcing
    write(fo,*) "Bool_Model_Snow"//spa, Bool_Model_Snow       ! (1) Usinsg snow model to find snow_depth, (0) Forcing with snow_depth
    write(fo,*) "Bool_Bessi"//spa, Bool_Bessi
    write(fo,*) "Bool_geometric"//spa, Bool_geometric
    write(fo,*) "Depth"//spa, Depth              !profondeur de la modélisation
    write(fo,*) "T_init"//spa, T_init                !température initiale a la surface
    write(fo,*) "T_freeze"//spa, T_freeze             !température où l'eau est considérée comme gelée
    write(fo,*) "freezing_range"//spa, freezing_range          ! Temperature at wich the snow start to melt
    write(fo,*) "Gfx"//spa, Gfx                  ! flux géothermique de la terre (a modifier peut être)
    write(fo,*) "Porosity_soil"//spa, Porosity_soil       ! porosité du sol
    write(fo,*) "organic_depth"//spa, organic_depth     ! profondeur de la couche organique
    write(fo,*) "n_organic"//spa, n_organic           ! porosité de la couche organique
    write(fo,*) "n_soil_bot"//spa, n_soil_bot         ! porosité en bas de la couche de sol
    write(fo,*) "alpha"//spa, alpha
    write(fo,*) "rho_snow_freeze"//spa, rho_snow_freeze            ! densité de la neige
    write(fo,*) "rho_water"//spa, rho_water          ! densité de l'eau
    write(fo,*) "rho_ice"//spa, rho_ice             ! densité de la glace
    write(fo,*) "rho_organic"//spa, rho_organic        ! densité de la matière organique
    write(fo,*) "rho_soil"//spa, rho_soil           ! densité du sol
    write(fo,*) "rho_snow_fresh"//spa, rho_snow_fresh            ! densité de la neige fraiche
    write(fo,*) "C_water"//spa, C_water          ! Capacité thermique massique de l'eau
    write(fo,*) "C_ice"//spa, C_ice            ! Capacité thermique massique de la glace
    write(fo,*) "C_organic"//spa, C_organic        ! Capacité thermique massique de la matière organique
    write(fo,*) "C_dry_soil"//spa, C_dry_soil        ! Capacité thermique massique du sol
    write(fo,*) "K_other_minerals"//spa, K_other_minerals    ! Conductivité thermique des autres minéraux
    write(fo,*) "K_quartz"//spa, K_quartz            ! conductivité thermique du quartz
    write(fo,*) "K_organic"//spa, K_organic        ! conductivité thermique de la matière organique
    write(fo,*) "K_ice"//spa, K_ice            ! conductivité thermique de la glace
    write(fo,*) "K_fluids"//spa, K_fluids         ! conductivité thermique des fluides
    write(fo,*) "q_quartz"//spa, q_quartz            ! pourcentage de quartz dans le sol
    write(fo,*) "gravity"//spa, gravity          ! accéleration gravitationnelle
    write(fo,*) "Latent_heat"//spa, Latent_heat    ! en J/kg
    write(fo,*) "s_l_max"//spa, s_l_max      ! nombre de couche de neige (marche que avec 1)
    write(fo,*) "gridNoMax"//spa, gridNoMax    !dmr spatial index, no assumption of the spatial arrangement
    write(fo,*) "Tempsolinit"//spa, Tempsolinit
    write(fo,*) "Tempairmonth"//spa, Tempairmonth
    write(fo,*) "Tempairday"//spa, Tempairday
    write(fo,*) "Tempsnowmonth"//spa, Tempsnowmonth
    write(fo,*) "Tempsnowday"//spa, Tempsnowday

    write(fo,*) "=== ************* ==="

    close(fo)

  end subroutine lecture_namelist


end module Parametrisation




















