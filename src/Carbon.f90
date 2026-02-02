#include "constant.h"
module Carbon

  !use blbla, only: var

    IMPLICIT NONE

    INTEGER, PARAMETER :: iactive = 1      !! Index for active carbon pool (unitless)
    INTEGER, PARAMETER :: islow = 2        !! Index for slow carbon pool (unitless)
    INTEGER, PARAMETER :: ipassive = 3     !! Index for passive carbon pool (unitless)
    INTEGER, PARAMETER :: ncarb = 3        !! Number of soil carbon pools (unitless)
    REAL :: f_a, f_s, f_p
    INTEGER :: c_perm_fich

#if (CARBON > 0)

contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
  subroutine carbon_main (Temp, altmax_lastyear, deepSOM_a, deepSOM_s, deepSOM_p, max_cryoturb_alt,          &
                          min_cryoturb_alt, diff_k_const, bio_diff_k_const, bioturbation_depth,      &
                          deepSOM, fc, Fv_l, fracgr, darea, deepSOM_tot,                             &
                          alpha_a, alpha_s, alpha_p, mu_soil_rev, beta_a, beta_s, beta_p) !b3_l, 

      use parameter_mod,  only : z_num, dt ! time step in seconds
      use parameter_mod, only: zf_soil
      use parameter_mod, only: D, dz


      real, dimension(z_num), intent(in)                  :: Temp
      !integer, dimension(z_num), intent(inout)            :: Temp_positive
      !integer, intent(in)                                 :: compteur_time_step
      !logical, intent(in)                                 :: end_year
      !REAL, intent(in)                                    :: clay
      REAL, intent(inout)                                 :: altmax_lastyear
      real, dimension(z_num), intent(inout)                :: deepSOM_a, deepSOM_s, deepSOM_p
      !real, dimension(z_num), intent(inout)               :: dz ! epaisseur de chaque couche
      !real, dimension(z_num), intent(in)                  :: D
      real, intent(in)                                    ::  max_cryoturb_alt, min_cryoturb_alt
      !REAL, DIMENSION(0:z_num),  INTENT(in)               :: zf_soil        !! depths of full levels (m)
      REAL, INTENT(in)                                    :: diff_k_const
      REAL, INTENT(in)                                    :: bio_diff_k_const
      real, intent(in)                                    :: bioturbation_depth
      REAL, DIMENSION(ncarb,ncarb), intent(in)            :: fc                         !! flux fractions within carbon pools
      REAL, dimension(z_num) , INTENT(out)                :: deepSOM
      REAL, INTENT(out)                                   :: deepSOM_tot
      REAL, intent(in)                                    :: Fv_l !b3_l,
      REAL, intent(in)                                    :: fracgr ! fraction of land in cell
      REAL, intent(in)                                    :: darea ! fraction of land in cell
      REAL                                  :: mu_soil_rev
      real, dimension(z_num)                :: alpha_a, alpha_s, alpha_p, beta_a, beta_s, beta_p


      !call compute_alt(Temp, Temp_positive, altmax_lastyear, compteur_time_step, end_year, altmax_lastyear, D)
      ! write(*,*) 'altmax_lastyear', altmax_lastyear

       ! redistribute carbon from biosphere model
      call carbon_redistribute(Temp, deepSOM_a, deepSOM_s, deepSOM_p, altmax_lastyear, Fv_l/360.) !Attention Fv in kg/m2/YEAR -> per day

       ! computes the decomposition in permafrost as a function of temperature (later : humidity and soil type?)
      call decomposition(Temp,  deepSOM_a, deepSOM_s, deepSOM_p, fc)

       ! cryoturbation et bioturbation
      call bio_cryoturbation(Temp, deepSOM_a, deepSOM_s, deepSOM_p, altmax_lastyear, max_cryoturb_alt, &
             min_cryoturb_alt, diff_k_const, bio_diff_k_const, bioturbation_depth,     &
             alpha_a, alpha_s, alpha_p, mu_soil_rev, beta_a, beta_s, beta_p)


      deepSOM(:) = deepSOM_a(:) + deepSOM_p(:) + deepSOM_s(:)
      deepSOM_tot = sum(deepSOM(:)*dz(:))*darea*fracgr


  end subroutine carbon_main

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!--------------------------
  subroutine carbon_first_init
  !initialise all the variables before the horizontal and time loops
  ! initialize carbon timestep

    use parameter_mod, only: YearType, z_num, one_day
    use parameter_mod, only: diff_k_const, bio_diff_k_const, cryoturbation_diff_k_in, bioturbation_diff_k_in
!~     use parameter_mod, only: ALT, altmax_lastyear
    use parameter_mod, only: zf_soil ! depth at full levels = lower layer-interface (m)
    use parameter_mod, only: zi_soil ! depth of intermediate levels (m)
    use parameter_mod, only: D, dz ! D =depth of lower level(m), dz=thickness (m)

    integer :: i

    diff_k_const=cryoturbation_diff_k_in/(one_day*YearType) !a verifier  si temps ok ?
    bio_diff_k_const=bioturbation_diff_k_in/(one_day*YearType)


    allocate(zf_soil(0:z_num))
    allocate(zi_soil(1:z_num))

    ! Define the soil layers
    zf_soil(0) = 0.0
    do i = 1, z_num
      zf_soil(i)=D(i)
      zi_soil(i)=zf_soil(i-1)+dz(i)/2.
    enddo
    !write(*,*) 'zf_soil', zf_soil(:)
    !write(*,*) 'zi_soil', zi_soil(:)

! dmr&mv --- moved to the spatial module (a.k.a. SpaceX)
!~     !Active layer depth
!~     ALT=0.0
!~     altmax_lastyear=0.0

    call open_carbon_output

  end subroutine carbon_first_init
!--------------------------


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  subroutine carbon_init(deepSOM_a, deepSOM_s, deepSOM_p,  fc, ALT, b4_spat, Temp) !b3_spat,
  subroutine carbon_init(deepSOM_a, deepSOM_s, deepSOM_p, fc, altmax_lastyear, b4_spat, Temp, deepSOM, deepSOM_tot, fracgr, darea, &
          alpha_a, alpha_s, alpha_p, mu_soil_rev, beta_a, beta_s, beta_p)

  ! initialize the carbon variables to 0
  !initialize carbon stocks with b4

    use parameter_mod, only: z_num
    use parameter_mod, only: dz

    !! carbon pools: indices

    real, dimension(z_num), intent(in)                  :: Temp
    !real, dimension(z_num), intent(inout)               :: dz ! epaisseur de chaque couche
    real,dimension(z_num), intent(inout) :: deepSOM_a, deepSOM_s, deepSOM_p
    REAL, DIMENSION(ncarb,ncarb), intent(out)  :: fc                         !! flux fractions within carbon pools
    REAL, DIMENSION(ncarb)                     :: fr                         !! fraction of decomposed carbon that goes into the atmosphere
    REAL, intent(out)                          :: altmax_lastyear ! , altmax_lastyear b3_spat,
    REAL, intent(out)                          :: b4_spat

    REAL, PARAMETER :: clay=1.0

    REAL, dimension(z_num) , INTENT(out)                :: deepSOM
    REAL, INTENT(out)                                   :: deepSOM_tot
    REAL, intent(in)                                    :: fracgr ! fraction of land in cell
    REAL, intent(in)                                    :: darea ! fraction of land in cell
    REAL                                  :: mu_soil_rev
    real, dimension(z_num)                :: alpha_a, alpha_s, alpha_p, beta_a, beta_s, beta_p


    deepSOM_a(:) = 0.0
    deepSOM_s(:) = 0.0
    deepSOM_p(:) = 0.0

    !altmax_lastyear = 0.0
    altmax_lastyear = 30 ! initialised to 30 m

    ! Only in mode uncoupled, fixed to a specific value
    !b3_spat = 1.
    !b4_spat = 1.

    ! calculate soil organic matter flux fractions
    fc(iactive,iactive) = 0.0
    fc(iactive,ipassive) = 0.004
    fc(iactive,islow) = 1.0 - (.85-.68*clay) - fc(iactive,ipassive)
    !
    fc(islow,islow) = .0
    fc(islow,iactive) = .42
    fc(islow,ipassive) = .03
    !
    fc(ipassive,ipassive) = .0
    fc(ipassive,iactive) = .45
    fc(ipassive,islow) = .0
    !
    fr(:) = 1.0-fc(:,iactive)-fc(:,islow)-fc(:,ipassive) ! pour verifier, faire un print, doit etre 0
    !firstcall = .FALSE.

! init to zero
      alpha_a(:)=0.0
      alpha_s(:)=0.0
      alpha_p(:)=0.0
      mu_soil_rev=0.0
      beta_a(:)=0.0
      beta_s(:)=0.0
      beta_p(:)=0.0



    !Redistribution of carbon from veget with Fv at beginning of year

    !parameters for initialisation of carbon pools
    f_a=0.01
    f_s=0.14
    f_p=0.85

    !write(*,*) 'b4_spat dans carbon_init', b4_spat
    call carbon_redistribute(Temp, deepSOM_a, deepSOM_s, deepSOM_p, altmax_lastyear, b4_spat) !b4 in kgC/m2/day

    deepSOM(:) = deepSOM_a(:) + deepSOM_p(:) + deepSOM_s(:)
    deepSOM_tot = sum(deepSOM(:)*dz(:))*darea*fracgr
    !write(*,*) 'deepSOM_tot in carbon_init in GtC', deepSOM_tot*1e-15
    !write(*,*) 'deepSOM_a, s, p in carbon_init', sum(deepSOM_a(:)*dz(:))*darea*fracgr, sum(deepSOM_s(:)*dz(:))*darea*fracgr,  sum(deepSOM_p(:)*dz(:))*darea*fracgr

    !modify values for redistribution of carbon
    f_a=1/3.
    f_s=1/3.
    f_p=1/3.


  end subroutine carbon_init
!--------------------------

!--------------------------
  subroutine carbon_redistribute(Temp, deepSOM_a, deepSOM_s, deepSOM_p, altmax_lastyear, b4) ! Temp,
  ! initialize the carbon in the active layer by redistributing carbon
  ! from litter and soil from vecode

    use parameter_mod, only: z_num, one_day, dt   ! z_num nombre d'horizons (couches)
    use parameter_mod, only: zf_soil
    use parameter_mod, only: dz

    real, dimension(z_num),intent(in)              :: Temp
    real,                      intent(in)          :: altmax_lastyear
    real,dimension(z_num),intent(inout)            :: deepSOM_a, deepSOM_s, deepSOM_p !soil organic matter (g /m^3 )
    real                , intent(in)               :: b4 ! Carbon in litter and soil, should be in g/m2 b3
    integer                                        :: k, il
    real                                           :: intdep !depth of integration
    real                                           :: z_lit !e folding depth
    real                                           :: som_input_TS
    real                                           :: dsom_litter
    real, dimension(z_num)                         :: som_profile, dsom_litter_z !  double precision?
    real                                           :: diff, verif_som ! double precision ?
    real                                           ::totcarbon

    !write(*,*) 'f_a in carbon_redistribute', f_a

!add input zfsoil
!define som_profile vecteur
    !totcarbon = SUM ((deepSOM_a (:) + deepSOM_s(:) + deepSOM_p(:) )*dz(:))
    !write (*,*) "debut routine redistribution" , totcarbon

!~     b3=1.! en g/m²  !!!! (du coup on divise som input par dz )
!~     b4=1 ! en g /m²
    !som_input_TS=b3+b4 ! matiere organique qui arrive dans le sol (litiere + sol) ! attention  a l unite
    som_input_TS=b4*1000*dt/one_day ! b4 de ilovecim en kg/m2/jour-> en g/m2/jour , matiere organique qui arrive dans le sol (sol) !
!add up the soil carbon from all veg pools, and change units from (gC/(m**2 of ground)/day) to gC/m^2 per timestep


    z_lit=0.2 ! arbitrary, peut aussi dependre de la prof des racines
    intdep = max(altmax_lastyear, z_lit) !depth of integration pour l instant = ALT, peut dependre de
!    la profondeur des racines... et ne peut etre plus petit que z_lit

    ! initialise valeur carbon dans couche active

    ! si couche active plus petite que les 2 premieres couches on repartit dans les 2 premieres couches
    !write(*,*) 'profondeur deux premieres couches', (dz(1)+dz(2))
    if (altmax_lastyear.le.(dz(1)+dz(2))) then
   ! if (altmax_lastyear.le.(dz(1))) then
        deepSOM_a(1)=deepSOM_a(1)+ f_a * (som_input_TS/2)/ dz(1)
        deepSOM_a(2)=deepSOM_a(2)+ f_a * (som_input_TS/2) / dz(2)
        deepSOM_s(1)=deepSOM_s(1)+ f_s * (som_input_TS/2)/ dz(1)
        deepSOM_s(2)=deepSOM_s(2)+ f_s * (som_input_TS/2) / dz(2)
        deepSOM_p(1)=deepSOM_p(1)+ f_p * (som_input_TS/2)/ dz(1)
        deepSOM_p(2)=deepSOM_p(2)+ f_p * (som_input_TS/2) / dz(2)
        !write(*,*) 'case 1'
    else

      dsom_litter=0.0
      dsom_litter_z=0.0

      som_profile(:)=0.0
      do il = 1, z_num ! sur la verticale
        som_profile(il) = 1.0 / ( 1.0 - EXP( -intdep / z_lit ) ) * &
                        ( EXP(-zf_soil(il-1)/z_lit)  - &
                        EXP( -zf_soil(il)/z_lit ) )
      enddo
     !The above calculation should result in a som_profile of 0.9999999 which is
                   ! pretty good but not good enough to keep the mass balance
                   ! closed. Put the
                   ! residual in the top layer (just because that layer should
                   ! always be present)
      som_profile(1) = 1.0 - SUM(som_profile(2:z_num))
      !endif


      ! Use the profile to redistribute the som_input
      DO il = 1, z_num ! ngrnd
                      ! Divide by the tickness of the layer, not the tickness of the whole
                      ! profile (which is the case when using z_lit). The unit of som_input_Ts
                      ! is gC/m^2 per timestep. We want dsom_litter_z to be in gC/m^3 per
                      ! timestep so divide by the layer depth over which the som input is to
                      ! be distributed.
       dsom_litter_z(il) = som_input_TS * &
                           som_profile(il)/ dz(il) !z_thickness(il)
                      ! The original code already used zf_soil. This may have resulted in an
                      ! inconsitency. Given that zi_soil has now been changed to zf_soil. The
                      ! original code is consistent with the new approach that uses zf_soil
       dsom_litter = dsom_litter + &
                           dsom_litter_z(il) * dz(il) ! z_thickness(il)
      END DO

     !write(*,*) 'dsom_litter, som_input_TS', dsom_litter, som_input_TS
     diff=dsom_litter-som_input_TS
     dsom_litter_z(1)=dsom_litter_z(1)-diff/dz(1)
     verif_som=sum(dsom_litter_z(:)* dz(:))
     !write(*,*) 'dsom_litter, som_input_TS after, diff, dsom_litter_z', verif_som, som_input_TS, diff, dsom_litter_z

     !deepSOM_a(:)=deepSOM_a(:)+1/3.*dsom_litter_z(:) !* dz (:)
     !deepSOM_s(:)=deepSOM_s(:)+1/3.*dsom_litter_z(:) !* dz (:)
     !deepSOM_p(:)=deepSOM_p(:)+1/3.*dsom_litter_z(:) !* dz (:)
     !deepSOM_a(:)=deepSOM_a(:)+0.01*dsom_litter_z(:) !* dz (:)
     !deepSOM_s(:)=deepSOM_s(:)+0.09*dsom_litter_z(:) !* dz (:)
     !deepSOM_p(:)=deepSOM_p(:)+0.9*dsom_litter_z(:) !* dz (:)
     deepSOM_a(:)=deepSOM_a(:)+f_a*dsom_litter_z(:) !* dz (:)
     deepSOM_s(:)=deepSOM_s(:)+f_s*dsom_litter_z(:) !* dz (:)
     deepSOM_p(:)=deepSOM_p(:)+f_p*dsom_litter_z(:) !* dz (:)


     endif

         !totcarbon = SUM ((deepSOM_a (:) + deepSOM_s(:) + deepSOM_p(:) ) *dz(:))
         !write (*,*) "deepSOM fin routine redistribution" , totcarbon


   !write(*,*) 'deepSOM_a dans redistribute', sum(deepSOM_a(:)*dz(:))
   !write(*,*) 'deepSOM_s dans redistribute', sum(deepSOM_s(:)*dz(:))
   !write(*,*) 'deepSOM_p dans redistribute', sum(deepSOM_p(:)*dz(:))


  end subroutine carbon_redistribute
!--------------------------


!--------------------------
  subroutine decomposition(Temp, deepSOM_a, deepSOM_s, deepSOM_p, fc) !h_pori?, clay
  !This routine calculates soil organic matter decomposition
  !decomposition depending on temperature, humidity, soil texture

    use parameter_mod, only: z_num, dt
    use parameter_mod, only: ZeroCelsius=>tK_zero_C
    use parameter_mod, only: dz, zi_soil

!~     !! carbon pools: indices
!~     INTEGER, PARAMETER :: iactive = 1      !! Index for active carbon pool (unitless)
!~     INTEGER, PARAMETER :: islow = 2        !! Index for slow carbon pool (unitless)
!~     INTEGER, PARAMETER :: ipassive = 3     !! Index for passive carbon pool (unitless)
!~     !INTEGER, PARAMETER :: isurface = 4     !! Index for passive carbon pool (unitless)
!~     !INTEGER, PARAMETER :: ncarb = 4        !! Number of soil carbon pools (unitless)
!~     INTEGER, PARAMETER :: ncarb = 3        !! Number of soil carbon pools (unitless)

    real, dimension(z_num),intent(in) :: Temp !, h_pori !, moist_in =humidity
    real, dimension(z_num)           :: temp_kelvin
    real, parameter                             :: q10 = 2.0

    logical, parameter                          :: limit_decomp_moisture = .false. ! to be changed to true if humidity in model
    real, dimension(z_num)           :: moistfunc_result
    REAL                                        :: stomate_tau = 4.699E6
    REAL                                        :: depth_modifier = 1.E6             !! e-folding depth of turnover rates,following Koven et al.,2013,Biogeosciences. A very large value means no depth modification
!    REAL, DIMENSION(:), INTENT(in)           :: clay            !! clay content
    REAL                                        :: fbact_a, fbact_s, fbact_p
    REAL, DIMENSION(z_num)                      :: fbact                      !! turnover constant (day)
    REAL, dimension(z_num)                      :: tempfunc_result
    REAL, DIMENSION(ncarb,ncarb), intent(in)    :: fc                         !! flux fractions within carbon pools
    REAL                                        :: dC
    REAL                                        :: fpassive = 1617.45                !! convertiing factor to go from active pool turnover to passive pool turnover from Guimberteau et al 2018 GMD
    REAL                                        :: fslow = 37.0                      !! convertiing factor to go from active pool turnover to slow pool turnover from Guimberteau et al 2018 GMD
    INTEGER                                     :: ij, il, iv, ip
    REAL                                        :: temp_local
    REAL, DIMENSION(z_num), intent(inout)       :: deepSOM_a, deepSOM_s, deepSOM_p !soil organic matter (g /m^3 )
    REAL, DIMENSION(ncarb, ncarb)               :: somflux
    REAL :: totalcarbon
    !real, dimension(z_num), intent(in)          :: dz


    totalcarbon = SUM ((deepSOM_a (:)  + deepSOM_s (:) + deepSOM_p (:))*dz(:))
    !write (*,*) " decomposition start ", totalcarbon
! Donne le taux de decomposition fbact (temps de résidence)
    temp_kelvin(:)=Temp(:)+ZeroCelsius
! cutoff respiration when T < -1C
       WHERE (temp_kelvin(:) .GT. ZeroCelsius ) ! normal
          tempfunc_result(:) = EXP(log(q10) * ( temp_kelvin(:) - (ZeroCelsius+30.) ) / 10. )
       ELSEWHERE (temp_kelvin(:) .GT. ZeroCelsius - 1. )  ! linear dropoff to zero
          tempfunc_result(:) = (temp_kelvin(:) - (ZeroCelsius - 1.)) * &
               EXP( log(q10) * ( ZeroCelsius - (ZeroCelsius+30.) ) / 10. )
       ELSEWHERE  ! zero
          tempfunc_result(:) = EPSILON(0.)
       ENDWHERE

       tempfunc_result(:) = MAX(MIN( 1.0, tempfunc_result(:) ), EPSILON(0.))

    IF ( limit_decomp_moisture ) THEN
       ! stomate moisture control function
!       moistfunc_result(:,:) = -1.1 * moist_in(:,:) * moist_in(:,:) + 2.4 * moist_in(:,:) - 0.29
!       moistfunc_result(:,:) = max( 0.25, min( 1.0, moistfunc_result(:,:) ) )
    ELSE
       !moistfunc_result(:) = 1.0
       moistfunc_result(:) = 0.25 !ici parameter to be modified
    ENDIF

    DO ij = 1, z_num
      fbact(ij) = stomate_tau/(moistfunc_result(ij) * tempfunc_result(ij)) / EXP(-zi_soil(ij)/depth_modifier)
    !write(*,*)'stomate_tau, moistfunc_result(ij), tempfunc_result(ij), zi_soil(ij), depth_modifier' &
    !         , stomate_tau, moistfunc_result(ij), tempfunc_result(ij), zi_soil(ij), depth_modifier

    ENDDO



    !write(*,*) 'fbact', fbact(:)



    !DO ip = 1, gridNoMax ! sur tous les points de grille

!! a voir          IF (  veget_mask_2d(ip,iv) ) THEN
             !
        DO il = 1, z_num! pour chaque niveau vertical !z_num=nombre de couches de sol
            !
            ! 1 function that gives soil organic matter residence time as a function of
            !     soil temperature (in seconds)
            !tprof= Temp+ZeroCelsius! temperature en Kelvin
            !temp = tprof(ip,il,iv) - ZeroCelsius
            temp_local=Temp(il)
            fbact_a = fbact(il) ! a =actif, puis modifie pour slow et passif
            fbact_a = MAX(fbact_a,dt)
                !
            IF ( fbact_a/HUGE(1.) .GT. .1 ) THEN
               fbact_s = fbact_a
               fbact_p = fbact_a
            ELSE
               fbact_s = fbact_a * fslow
               fbact_p = fbact_a * fpassive
            ENDIF


            ! 2 oxic decomposition: carbon and oxygen consumption
            !
            ! 2.1 active
            !
            !cas oxic (pour nous)
            dC = deepSOM_a(il) * dt/fbact_a ! difference (deltacarbon)
            !write(*,*) 'dC active', dC

            !dC = dC * ( 1. - som_turn_iactive_clay_frac * clay(ip) ) ! modification avec la texture du sol, a voir ?

            ! flux vers les autres reservoirs
            somflux(iactive,ipassive) = fc(iactive,ipassive) * dC ! d'ou ca part, ou ca va
            somflux(iactive,islow) = fc(iactive,islow) * dC

            deepSOM_a(il) = deepSOM_a(il) - dC
            !dO2 = wO2/wC * dC(icarbon)*fr(ip,iactive,iv) / totporO2_soil(ip,il,iv)
            !O2_soil(ip,il,iv) = MAX( O2_soil(ip,il,iv) - dO2, zero)
            ! keep delta C * fr in memory (generates energy)

            !deltaSOM1_a(ip,il,iv,iele) = dC(iele)*fr(ip,iactive,iv) !!this line!!! A l'air pas utilise, a verifier


            ! 2.2 slow
            !

            dC = deepSOM_s(il) * dt/fbact_s

            ! flux vers les autres reservoirs
            somflux(islow,iactive) = fc(islow,iactive) * dC
            somflux(islow,ipassive) = fc(islow,ipassive) * dC

            !
            deepSOM_s(il) = deepSOM_s(il) - dC
            !dO2 = wO2/wC * dC(iele)*fr(ip,islow,iv) / totporO2_soil(ip,il,iv)
            !O2_soil(ip,il,iv) = MAX( O2_soil(ip,il,iv) - dO2, zero)
            ! keep delta C * fr in memory (generates energy)
            !deltaSOM1_s(il, ip) = dC*fr(ip,islow) !!this line!!!

            !
            ! 2.3 passive
            !

            dC = deepSOM_p(il) * dt/fbact_p

            ! flux vers les autres reservoirs

            somflux(ipassive,iactive) = fc(ipassive,iactive) * dC
            somflux(ipassive,islow) = fc(ipassive,islow) * dC

            deepSOM_p(il) = deepSOM_p(il) - dC

            ! keep delta C * fr in memory (generates energy)

            !deltaSOM1_p(il, ip) = dC*fr(ip,ipassive) !!this line!!!



            ! 4 add fluxes between reservoirs

            deepSOM_a(il)=deepSOM_a(il)+somflux(islow,iactive)+somflux(ipassive,iactive)
            deepSOM_s(il)=deepSOM_s(il)+somflux(iactive,islow)+somflux(ipassive,islow)
            deepSOM_p(il)=deepSOM_p(il)+somflux(iactive,ipassive)+somflux(islow,ipassive)

        ENDDO !z_num

        !ENDIF ! veget_mask_2d

   ! ENDDO ! End loop over GridNoMax
       totalcarbon = SUM ((deepSOM_a (:)  + deepSOM_s (:) + deepSOM_p (:))*dz(:))
       !write (*,*) " decomposition fin ", totalcarbon

  end subroutine decomposition
!--------------------------


!--------------------------
  subroutine bio_cryoturbation(Temp, deepSOM_a, deepSOM_s, deepSOM_p, altmax_lastyear, max_cryoturb_alt, &
   min_cryoturb_alt, diff_k_const, bio_diff_k_const, bioturbation_depth,                 &
   alpha_a, alpha_s, alpha_p, mu_soil_rev, beta_a, beta_s, beta_p)

  ! Vertical mixing by cryoturbation with a diffusion scheme

    use parameter_mod, only: z_num, YearType ! nombre de jours par an
    use parameter_mod, only: zf_soil, zi_soil, dz
    use parameter_mod, only: dt
    use parameter_mod, only: min_stomate

    real, dimension(z_num), intent(in)           :: Temp !, h_pori
    REAL, DIMENSION(z_num), intent(inout) :: deepSOM_a, deepSOM_s, deepSOM_p !soil organic matter (g /m^3 )
    real, intent(in)                             :: altmax_lastyear
    real, intent(in)                      :: max_cryoturb_alt, min_cryoturb_alt !altmax_lastyear
    real                                  :: cryoturbation_depth
    logical                                      :: cryoturb_location, bioturb_location
    REAL, DIMENSION(z_num)                :: diff_k               !! Diffusion constant (m^2/s)
    REAL, INTENT(in)                      :: diff_k_const
    REAL, INTENT(in)                      :: bio_diff_k_const
    real, dimension(z_num)                :: xc_cryoturb, xd_cryoturb
    real, dimension(z_num)                :: alpha_a, alpha_s, alpha_p, beta_a, beta_s, beta_p
    REAL                                  :: mu_soil_rev
    real                                  :: xe_a, xe_s, xe_p
    real, intent(in)                      :: bioturbation_depth
    integer                                      :: il
    real                                  :: totalcarbon1, totalcarbon2
    real                                  :: altSOM_a_old , altSOM_s_old, altSOM_p_old !soil organic matter (g /m^2 )
    real                                  :: altSOM_a , altSOM_s, altSOM_p !soil organic matter (g /m^2 )
    real                                  :: surfC_totake_a, surfC_totake_s, surfC_totake_p
    integer                                      :: igrnd
    real                                  :: pool_start, pool_end

! Test carbon tot

       totalcarbon1 = SUM ((deepSOM_a (:)  + deepSOM_s (:) + deepSOM_p (:))*dz(:))
       !write (*,*) " tot carbon start biocryo" , totalcarbon1

! Calcul des coefficients

      cryoturb_location = .false.
      bioturb_location = .false.

! Si on a du gel - degel en surface : on fait de la cryoturbation (max_cryoturb_alt =3m)
       cryoturb_location =  ( altmax_lastyear .LT. max_cryoturb_alt ) &
!In the former vertical discretization scheme the first level was at 0.016 cm; now it's only 0.00048 so we set an equivalent threshold directly as a fixed depth of 1 cm,
            .AND. ( altmax_lastyear .GE. min_cryoturb_alt )
!       IF (use_fixed_cryoturbation_depth) THEN
!          cryoturbation_depth(:,:) = fixed_cryoturbation_depth(:,:) ! garder ce cas pour faire des tests en particulier comparaisons avec donnes sur site
!       ELSE
          cryoturbation_depth = altmax_lastyear ! ce cas a utiliser par defaut
!       ENDIF

! partout ou pas de cryoturbation - > bioturbation avec coefficient plus petit
       bioturb_location = ( altmax_lastyear .GE. max_cryoturb_alt )

! si on est dans une zone de cryoturbation
       IF ( cryoturb_location ) THEN
        ! cas 4 dans orchidee
           DO il = 1, z_num ! linear dropoff to zero between alt and 3*alt
               IF ( zi_soil(il) .LE. cryoturbation_depth ) THEN
                   diff_k(il) = diff_k_const ! dans la couche active : constante
               ELSE ! en dessous couche active
                   diff_k(il) = diff_k_const*(1.-MAX(MIN((zi_soil(il)-cryoturbation_depth)/ &
                                 (2.*cryoturbation_depth),1.),0.))
               ENDIF
               IF ( zf_soil(il) .GT. max_cryoturb_alt ) THEN ! tout en dessous mis a zero
                   diff_k(il) = 0.0
               ENDIF
           END DO

! sinon bioturbation
       ELSE IF ( bioturb_location ) THEN
           DO il = 1, z_num
               IF ( zi_soil(il) .LE. bioturbation_depth ) THEN ! bioturbation_depth=2m
                   diff_k(il) = bio_diff_k_const
               ELSE
                   diff_k(il) = 0.0
               ENDIF
           END DO
       ELSE
           diff_k(:) = 0.0
       END IF

! lien avec la surface, mu dans annexe these Frederique Hourdin
       mu_soil_rev=diff_k(1)*dt/(zf_soil(1)-zf_soil(0))/(zi_soil(2)-zi_soil(1))

       IF ( cryoturb_location .OR. bioturb_location ) then
          DO il = 1,z_num-1 ! on calcule couche par couche de la surface vers le fond
             xc_cryoturb(il) = (zf_soil(il)-zf_soil(il-1))  / dt
             xd_cryoturb(il) = diff_k(il) / (zi_soil(il+1)-zi_soil(il))
          ENDDO

          !bottom
          xc_cryoturb(z_num) = (zf_soil(z_num)-zf_soil(z_num-1))  / dt

          xe_a = xc_cryoturb(z_num)+xd_cryoturb(z_num-1)
          xe_s = xc_cryoturb(z_num)+xd_cryoturb(z_num-1)
          xe_p = xc_cryoturb(z_num)+xd_cryoturb(z_num-1)
          alpha_a(z_num-1) = xd_cryoturb(z_num-1) / xe_a
          alpha_s(z_num-1) = xd_cryoturb(z_num-1) / xe_s
          alpha_p(z_num-1) = xd_cryoturb(z_num-1) / xe_p
          beta_a(z_num-1) = xc_cryoturb(z_num)*deepSOM_a(z_num) / xe_a
          beta_s(z_num-1) = xc_cryoturb(z_num)*deepSOM_s(z_num) / xe_s
          beta_p(z_num-1) = xc_cryoturb(z_num)*deepSOM_p(z_num) / xe_p

          !other levels
          DO il = z_num-2,1,-1
                xe_a = xc_cryoturb(il+1) + (1.-alpha_a(il+1))*xd_cryoturb(il+1) + xd_cryoturb(il)
                xe_s = xc_cryoturb(il+1) + (1.-alpha_s(il+1))*xd_cryoturb(il+1) + xd_cryoturb(il)
                xe_p = xc_cryoturb(il+1) + (1.-alpha_p(il+1))*xd_cryoturb(il+1) + xd_cryoturb(il)
                alpha_a(il) = xd_cryoturb(il) / xe_a
                alpha_s(il) = xd_cryoturb(il) / xe_s
                alpha_p(il) = xd_cryoturb(il) / xe_p
                beta_a(il) = (xc_cryoturb(il+1)*deepSOM_a(il+1) + &
                     xd_cryoturb(il+1)*beta_a(il+1)) / xe_a
                beta_s(il) = (xc_cryoturb(il+1)*deepSOM_s(il+1) + &
                     xd_cryoturb(il+1)*beta_s(il+1)) / xe_s
                beta_p(il) = (xc_cryoturb(il+1)*deepSOM_p(il+1) + &
                     xd_cryoturb(il+1)*beta_p(il+1)) / xe_p
          ENDDO

       ENDIF

! verification pour conservation masse de carbone (et azote)
          pool_start = 0.0
          DO igrnd = 1, z_num !ngrnd
             pool_start = pool_start + &
                  (deepSOM_a(igrnd) + deepSOM_s(igrnd) + deepSOM_p(igrnd)) * &
                  (zf_soil(igrnd)-zf_soil(igrnd-1))
                  !write(*,*) 'diff zf , dz', zf_soil(igrnd)-zf_soil(igrnd-1), dz(igrnd)
          END DO


! Calcul de la diffusion

         ! 1. calculate the total soil organic matter in the active layer
          altSOM_a_old = 0.0
          altSOM_s_old = 0.0
          altSOM_p_old = 0.0
          altSOM_a = 0.0
          altSOM_s = 0.0
          altSOM_p = 0.0

       IF ( cryoturb_location .OR. bioturb_location )THEN
                   ! 1. calculate the total soil organic matter
                   DO il = 1, z_num
                      altSOM_a_old = altSOM_a_old + deepSOM_a(il)*(zf_soil(il)-zf_soil(il-1))
                      altSOM_s_old = altSOM_s_old + deepSOM_s(il)*(zf_soil(il)-zf_soil(il-1))
                      altSOM_p_old = altSOM_p_old + deepSOM_p(il)*(zf_soil(il)-zf_soil(il-1))
                   ENDDO

                   ! 2. diffuse the soil organic matter
                   deepSOM_a(1) = (deepSOM_a(1)+mu_soil_rev*beta_a(1)) / &
                        (1.+mu_soil_rev*(1.-alpha_a(1)))
                   !write(*,*) 'a ', deepSOM_a(1), mu_soil_rev, beta_a(1), 1.+mu_soil_rev*(1.-alpha_a(1))
                   deepSOM_s(1) = (deepSOM_s(1)+mu_soil_rev*beta_s(1)) / &
                        (1.+mu_soil_rev*(1.-alpha_s(1)))
                   deepSOM_p(1) = (deepSOM_p(1)+mu_soil_rev*beta_p(1)) / &
                        (1.+mu_soil_rev*(1.-alpha_p(1)))

                   DO il = 2, z_num
                      deepSOM_a(il) = alpha_a(il-1)*deepSOM_a(il-1) + beta_a(il-1)
                      deepSOM_s(il) = alpha_s(il-1)*deepSOM_s(il-1) + beta_s(il-1)
                      deepSOM_p(il) = alpha_p(il-1)*deepSOM_p(il-1) + beta_p(il-1)
                   ENDDO

                   !! 3. recalculate the total soil organic matter
                   DO il = 1, z_num
                      altSOM_a = altSOM_a + deepSOM_a(il)*(zf_soil(il)-zf_soil(il-1))
                      altSOM_s = altSOM_s + deepSOM_s(il)*(zf_soil(il)-zf_soil(il-1))
                      altSOM_p = altSOM_p + deepSOM_p(il)*(zf_soil(il)-zf_soil(il-1))
                   ENDDO

                   IF ( altSOM_a_old > min_stomate .AND. &
                        (ABS(altSOM_a-altSOM_a_old)/altSOM_a_old.GT.min_stomate) ) THEN
                      WRITE (*,*) 'DZ warn: cryoturbate: total C not conserved, A ', &
                           'diff=',altSOM_a,altSOM_a_old,altSOM_a-altSOM_a_old,                  &
                           (altSOM_a-altSOM_a_old)/altSOM_a_old
                      !CALL ipslerr_p (3,'cryoturbate','','','')
                      deepSOM_a(1)=deepSOM_a(1)-(altSOM_a-altSOM_a_old)/(zf_soil(1)-zf_soil(0))
                   ENDIF

                   IF ( altSOM_s_old > min_stomate .AND. &
                        (ABS(altSOM_s-altSOM_s_old)/altSOM_s_old.GT.min_stomate) ) THEN
                      WRITE (*,*) 'DZ warn: cryoturbate: total C not conserved, S ', &
                           'diff=',altSOM_s,altSOM_s_old,altSOM_s-altSOM_s_old,     &
                           (altSOM_s-altSOM_s_old)/altSOM_s_old
                      !CALL ipslerr_p (3,'cryoturbate','','','')
                      deepSOM_s(1)=deepSOM_s(1)-(altSOM_s-altSOM_s_old)/(zf_soil(1)-zf_soil(0))
                   ENDIF

                   IF ( altSOM_p_old > min_stomate .AND. &
                        (ABS(altSOM_p-altSOM_p_old)/altSOM_p_old.GT.min_stomate) ) THEN
                      WRITE (*,*) 'DZ warn: cryoturbate: total C not conserved, P ', &
                           'diff=',altSOM_p,altSOM_p_old,altSOM_p-altSOM_p_old,     &
                           (altSOM_p-altSOM_p_old)/altSOM_p_old
                      !CALL ipslerr_p (3,'cryoturbate','','','')
                      deepSOM_p(1)=deepSOM_p(1)-(altSOM_p-altSOM_p_old)/(zf_soil(1)-zf_soil(0))
                   ENDIF


                   ! dans orchidee pas utilise

                   ! 4. subtract the organic matter in the top layer(s) so that the total organic matter content of the active layer is conserved.
                   ! for now remove this correction term...
!                   surfC_totake_a(ip,iv) = (altC_a(ip,iv)-altC_a_old(ip,iv))/(zf_soil(altmax_ind(ip,iv))-zf_soil(0))
!                   surfC_totake_s(ip,iv) = (altC_s(ip,iv)-altC_s_old(ip,iv))/(zf_soil(altmax_ind(ip,iv))-zf_soil(0))
!                   surfC_totake_p(ip,iv) = (altC_p(ip,iv)-altC_p_old(ip,iv))/(zf_soil(altmax_ind(ip,iv))-zf_soil(0))
!                   deepC_a(ip,1:altmax_ind(ip,iv),iv) = deepC_a(ip,1:altmax_ind(ip,iv),iv) - surfC_totake_a(ip,iv)
!                   deepC_s(ip,1:altmax_ind(ip,iv),iv) = deepC_s(ip,1:altmax_ind(ip,iv),iv) - surfC_totake_s(ip,iv)
!                   deepC_p(ip,1:altmax_ind(ip,iv),iv) = deepC_p(ip,1:altmax_ind(ip,iv),iv) - surfC_totake_p(ip,iv)

!                   ! if negative values appear, we don't subtract the delta-C
!                   from top layers
!                   IF (ANY(deepC_a(ip,1:altmax_ind(ip,iv),iv) .LT. zero) ) THEN
!                      deepC_a(ip,1:altmax_ind(ip,iv),iv)=deepC_a(ip,1:altmax_ind(ip,iv),iv)+surfC_totake_a(ip,iv)
!                      IF (altC_a(ip,iv) .GT. zero) THEN
!                         deepC_a(ip,:,iv)=deepC_a(ip,:,iv)*altC_a_old(ip,iv)/altC_a(ip,iv)
!                      ENDIF
!                   ENDIF
!                   IF (ANY(deepC_s(ip,1:altmax_ind(ip,iv),iv) .LT. zero) ) THEN
!                      deepC_s(ip,1:altmax_ind(ip,iv),iv)=deepC_s(ip,1:altmax_ind(ip,iv),iv)+surfC_totake_s(ip,iv)
!                      IF (altC_s(ip,iv) .GT. zero) THEN
!                         deepC_s(ip,:,iv)=deepC_s(ip,:,iv)*altC_s_old(ip,iv)/altC_s(ip,iv)
!                      ENDIF
!                   ENDIF
!                   IF (ANY(deepC_p(ip,1:altmax_ind(ip,iv),iv) .LT. zero) ) THEN
!                      deepC_p(ip,1:altmax_ind(ip,iv),iv)=deepC_p(ip,1:altmax_ind(ip,iv),iv)+surfC_totake_p(ip,iv)
!                      IF (altC_p(ip,iv) .GT. zero) THEN
!                         deepC_p(ip,:,iv)=deepC_p(ip,:,iv)*altC_p_old(ip,iv)/altC_p(ip,iv)
!                      ENDIF
!                   ENDIF


       ENDIF

       !! 4.1.1 Calculate components of the mass balance
       pool_end = 0.0
          DO igrnd = 1, z_num !ngrnd
             pool_end = pool_end + &
                  (deepSOM_a(igrnd) + deepSOM_s(igrnd) + deepSOM_p(igrnd)) * &
                  (zf_soil(igrnd)-zf_soil(igrnd-1))
          END DO

        !write(*,*) 'pool_start, pool_end', pool_start, pool_end

        totalcarbon2 = SUM ((deepSOM_a (:)  + deepSOM_s (:) + deepSOM_p (:)) *dz(:))
        !write (*,*) " tot carbon fin biocryo" , totalcarbon2
        !if (abs(totalcarbon2-totalcarbon1).gt.min_stomate) then
        !  write (*,*) " tot carbon debut fin biocryo" , totalcarbon1, totalcarbon2
        !endif

  end subroutine bio_cryoturbation
!--------------------------

  subroutine update_orgalayer_indx(deepSOM, orgalayer_indx)
  ! Update of organic layer depth index depending on the max depth of carbon
  ! (deepSOM)

    use parameter_mod, only: z_num

    REAL, dimension(z_num) , INTENT(in)                :: deepSOM
    INTEGER, INTENT(out)                               :: orgalayer_indx
    integer igrnd

      orgalayer_indx=0
      do igrnd = 1, z_num !ngrnd for each vertical level
         if (deepSOM(igrnd).gt. 1e-12) then !if some carbon in this level
             !write(*,*) 'deepSOM, index', deepSOM(igrnd), igrnd
             orgalayer_indx=igrnd ! keep index
         endif
      enddo

  end subroutine update_orgalayer_indx
!--------------------------

  subroutine open_carbon_output()
  ! Open file for carbon output

        OPEN(newunit=c_perm_fich, file='outputdata/carbon/C_permafrost.txt',status='unknown')

        WRITE(c_perm_fich,'(1A20)') "DeepSOM_tot (GtC)"

  end subroutine open_carbon_output
!--------------------------

  subroutine write_carbon_output(deepSOM_tot)
  ! Write total carbon content in output file

    REAL           :: deepSOM_tot

    !write (c_perm_fich,'(1f16.2)') deepSOM_tot
    write (c_perm_fich,*) deepSOM_tot*1e-15
    write(*,*) 'depSOMtot write', deepSOM_tot*1e-15

  end subroutine write_carbon_output
!--------------------------

  subroutine close_carbon_output()
  ! Close file for carbon output

     CLOSE(c_perm_fich)

  end subroutine close_carbon_output
!--------------------------



#endif

end module Carbon



