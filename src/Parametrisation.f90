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
    CHARACTER(len=30) :: to_print

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

    write(to_print,'(a30)') "namerun"
    write(fo,*) adjustl(to_print), trim(namerun)
    write(to_print,'(a30)') "TotTime"
    write(fo,*) adjustl(to_print), TotTime                    !temps total en année
    write(to_print,'(a30)') "nb_day_per_month"
    write(fo,*) adjustl(to_print), nb_day_per_month           !nombre de jour entre chaque pas de temps
    write(to_print,'(a30)') "t_fin"
    write(fo,*) adjustl(to_print), t_fin
    write(to_print,'(a30)') "YearType"
    write(fo,*) adjustl(to_print), YearType         !nombre de jour par an
    write(to_print,'(a30)') "z_num"
    write(fo,*) adjustl(to_print), z_num            !nombre de couches étudiée
    write(to_print,'(a30)') "GridType"
    write(fo,*) adjustl(to_print), GridType           !(1) Log-generated, (2) Linear-generated
    write(to_print,'(a30)') "PorosityType"
    write(fo,*) adjustl(to_print), PorosityType       !(1) linéaire, (autre) exponentiellement décroissante en fonction de la profondeur
    write(to_print,'(a30)') "Bool_Snow"
    write(fo,*) adjustl(to_print), Bool_Snow          ! forçage en neige ou non (1 ou 0)
    write(to_print,'(a30)') "Bool_Organic"
    write(fo,*) adjustl(to_print), Bool_Organic       ! prise en compte de la couche organique ou non (1 ou 0)
    write(to_print,'(a30)') "EQ_Tr"
    write(fo,*) adjustl(to_print), EQ_Tr             ! Equilibrum run (0) or Transient run (1) -> using different forcing Temperature and snow
    write(to_print,'(a30)') "EQ1_EQ2"
    write(fo,*) adjustl(to_print), EQ1_EQ2           ! EQ1(1) initial temperature calculated with the Geothermal heat flux. EQ2 initial temperature read in a file .txt
    write(to_print,'(a30)') "Bool_delta"
    write(fo,*) adjustl(to_print), Bool_delta        !
    write(to_print,'(a30)') "Bool_glacial"
    write(fo,*) adjustl(to_print), Bool_glacial          ! Using glacial index to modify air temperature
    write(to_print,'(a30)') "Bool_layer_temp"
    write(fo,*) adjustl(to_print), Bool_layer_temp       ! Creation of .txt with the temperature of the soil at different layer
    write(to_print,'(a30)') "Bool_Swe_Snw"
    write(fo,*) adjustl(to_print), Bool_Swe_Snw          ! (1) Snow forcing, (0) Swe forcing
    write(to_print,'(a30)') "Bool_Model_Snow"
    write(fo,*) adjustl(to_print), Bool_Model_Snow       ! (1) Usinsg snow model to find snow_depth, (0) Forcing with snow_depth
    write(to_print,'(a30)') "Bool_Bessi"
    write(fo,*) adjustl(to_print), Bool_Bessi
    write(to_print,'(a30)') "Bool_geometric"
    write(fo,*) adjustl(to_print), Bool_geometric
    write(to_print,'(a30)') "Depth"
    write(fo,*) adjustl(to_print), Depth              !profondeur de la modélisation
    write(to_print,'(a30)') "T_init"
    write(fo,*) adjustl(to_print), T_init                !température initiale a la surface
    write(to_print,'(a30)') "T_freeze"
    write(fo,*) adjustl(to_print), T_freeze             !température où l'eau est considérée comme gelée
    write(to_print,'(a30)') "freezing_range"
    write(fo,*) adjustl(to_print), freezing_range          ! Temperature at wich the snow start to melt
    write(to_print,'(a30)') "Gfx"
    write(fo,*) adjustl(to_print), Gfx                  ! flux géothermique de la terre (a modifier peut être)
    write(to_print,'(a30)') "Porosity_soil"
    write(fo,*) adjustl(to_print), Porosity_soil       ! porosité du sol
    write(to_print,'(a30)') "organic_depth"
    write(fo,*) adjustl(to_print), organic_depth     ! profondeur de la couche organique
    write(to_print,'(a30)') "n_organic"
    write(fo,*) adjustl(to_print), n_organic           ! porosité de la couche organique
    write(to_print,'(a30)') "n_soil_bot"
    write(fo,*) adjustl(to_print), n_soil_bot         ! porosité en bas de la couche de sol
    write(to_print,'(a30)') "alpha"
    write(fo,*) adjustl(to_print), alpha
    write(to_print,'(a30)') "rho_snow_freeze"
    write(fo,*) adjustl(to_print), rho_snow_freeze            ! densité de la neige
    write(to_print,'(a30)') "rho_water"
    write(fo,*) adjustl(to_print), rho_water          ! densité de l'eau
    write(to_print,'(a30)') "rho_ice"
    write(fo,*) adjustl(to_print), rho_ice             ! densité de la glace
    write(to_print,'(a30)') "rho_organic"
    write(fo,*) adjustl(to_print), rho_organic        ! densité de la matière organique
    write(to_print,'(a30)') "rho_soil"
    write(fo,*) adjustl(to_print), rho_soil           ! densité du sol
    write(to_print,'(a30)') "rho_snow_fresh"
    write(fo,*) adjustl(to_print), rho_snow_fresh            ! densité de la neige fraiche
    write(to_print,'(a30)') "C_water"
    write(fo,*) adjustl(to_print), C_water          ! Capacité thermique massique de l'eau
    write(to_print,'(a30)') "C_ice"
    write(fo,*) adjustl(to_print), C_ice            ! Capacité thermique massique de la glace
    write(to_print,'(a30)') "C_organic"
    write(fo,*) adjustl(to_print), C_organic        ! Capacité thermique massique de la matière organique
    write(to_print,'(a30)') "C_dry_soil"
    write(fo,*) adjustl(to_print), C_dry_soil        ! Capacité thermique massique du sol
    write(to_print,'(a30)') "K_other_minerals"
    write(fo,*) adjustl(to_print), K_other_minerals    ! Conductivité thermique des autres minéraux
    write(to_print,'(a30)') "K_quartz"
    write(fo,*) adjustl(to_print), K_quartz            ! conductivité thermique du quartz
    write(to_print,'(a30)') "K_organic"
    write(fo,*) adjustl(to_print), K_organic        ! conductivité thermique de la matière organique
    write(to_print,'(a30)') "K_ice"
    write(fo,*) adjustl(to_print), K_ice            ! conductivité thermique de la glace
    write(to_print,'(a30)') "K_fluids"
    write(fo,*) adjustl(to_print), K_fluids         ! conductivité thermique des fluides
    write(to_print,'(a30)') "q_quartz"
    write(fo,*) adjustl(to_print), q_quartz            ! pourcentage de quartz dans le sol
    write(to_print,'(a30)') "gravity"
    write(fo,*) adjustl(to_print), gravity          ! accéleration gravitationnelle
    write(to_print,'(a30)') "Latent_heat"
    write(fo,*) adjustl(to_print), Latent_heat    ! en J/kg
    write(to_print,'(a30)') "s_l_max"
    write(fo,*) adjustl(to_print), s_l_max      ! nombre de couche de neige (marche que avec 1)
    write(to_print,'(a30)') "gridNoMax"
    write(fo,*) adjustl(to_print), gridNoMax    !dmr spatial index, no assumption of the spatial arrangement
    write(to_print,'(a30)') "Tempsolinit"
    write(fo,*) adjustl(to_print), trim(Tempsolinit)
    write(to_print,'(a30)') "Tempairmonth"
    write(fo,*) adjustl(to_print), trim(Tempairmonth)
    write(to_print,'(a30)') "Tempairday"
    write(fo,*) adjustl(to_print), trim(Tempairday)
    write(to_print,'(a30)') "Tempsnowmonth"
    write(fo,*) adjustl(to_print), trim(Tempsnowmonth)
    write(to_print,'(a30)') "Tempsnowday"
    write(fo,*) adjustl(to_print), trim(Tempsnowday)

    write(fo,*) "=== ************* ==="

    close(fo)

  end subroutine lecture_namelist


end module Parametrisation




















