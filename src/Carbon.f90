module Carbon

  !use blbla, only: var


  Implicit none


contains


!--------------------------
  subroutine carbon_first_init(dt, max_cryoturb_alt, min_cryoturb_alt, diff_k_const, bio_diff_k_const, &
             bioturbation_depth, zi_soil, dz, zf_soil, ALT, altmax_lastyear)
  !initialise all the variables before the horizontal and time loops
  ! initialize carbon timestep

    use parameter_mod, only: nb_day_per_month, YearType, z_num

    REAL, intent(out)                 :: dt !time step in seconds
    REAL, intent(out)                 :: max_cryoturb_alt
    REAL, intent(out)                 :: min_cryoturb_alt
    REAL                              :: cryoturbation_diff_k_in !! input time constant of cryoturbation (m^2/y)
    REAL                              :: bioturbation_diff_k_in  !! input time constant of bioturbation (m^2/y)
    REAL, intent(out)                 :: diff_k_const
    REAL, intent(out)                 :: bio_diff_k_const
    REAL                              :: one_day
    REAL, intent(out)                 :: bioturbation_depth
    REAL, dimension (z_num), intent(in)  :: zi_soil
    REAL, dimension (z_num), intent(in)  :: dz
    REAL, dimension (0:z_num), intent(out) :: zf_soil
    REAL, intent(inout)               :: ALT
    REAL, intent(inout)               :: altmax_lastyear

    one_day=24*60*60
    max_cryoturb_alt= 3 !m
    min_cryoturb_alt=0.01 !1cm
    cryoturbation_diff_k_in = .001
    diff_k_const=cryoturbation_diff_k_in/(one_day*YearType) !a verifier  si temps ok ?
    bioturbation_diff_k_in = 0.0001
    bio_diff_k_const=bioturbation_diff_k_in/(one_day*YearType)
    bioturbation_depth=2 !m

    ! Define the soil layers
    zf_soil(:) = 0.0
    zf_soil(1:z_num)=zi_soil(:)+dz(:)/2.
    zf_soil(0) = 0.

    !Active layer depth
    ALT=0.0
    altmax_lastyear=0.0
    

  end subroutine carbon_first_init
!--------------------------

!--------------------------
  subroutine carbon_init(deepSOM_a, deepSOM_s, deepSOM_p,  fc, clay)
  ! initialize the carbon variables to 0

    use parameter_mod, only: z_num

    !! carbon pools: indices
    INTEGER, PARAMETER :: iactive = 1      !! Index for active carbon pool (unitless)
    INTEGER, PARAMETER :: islow = 2        !! Index for slow carbon pool (unitless)
    INTEGER, PARAMETER :: ipassive = 3     !! Index for passive carbon pool (unitless)
    INTEGER, PARAMETER :: ncarb = 3        !! Number of soil carbon pools (unitless)

    real,dimension(z_num), intent(inout) :: deepSOM_a, deepSOM_s, deepSOM_p
    REAL, DIMENSION(ncarb,ncarb), intent(out)  :: fc                         !! flux fractions within carbon pools
    REAL, DIMENSION(ncarb)                     :: fr                         !! fraction of decomposed carbon that goes into the atmosphere
    REAL, intent(out)                          :: clay

    deepSOM_a(:) = 0.0
    deepSOM_s(:) = 0.0
    deepSOM_p(:) = 0.0

    clay=1.0

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

  end subroutine carbon_init
!--------------------------

!--------------------------
  subroutine compute_alt(Temp, Temp_positive, ALT, compteur_time_step, end_year, altmax_lastyear, D)
  ! computes alt (active layer thickness) each year

    use parameter_mod, only: z_num, YearType, nb_day_per_month

    real, dimension(z_num), intent(in)                  ::   Temp
    real, dimension(z_num), intent(inout)               :: Temp_positive
    real,  intent(inout)                                :: ALT ! active layer thickness
    real,  intent(out)                                  :: altmax_lastyear ! alt last year
    integer, intent(in)                                 :: compteur_time_step
    integer, intent(in)                                 :: end_year
    integer                                             :: k, i, ALT_index
    real, dimension(z_num), intent(in)                  :: D

   !initialisation au debut de l annee
   !write(*,*) 'compteur_time_step', compteur_time_step
   if (compteur_time_step.eq.1) then 
      Temp_positive(:) = 0.0
   endif

   !au cours de l annee passe a 1 la ou temp > 0
   !where(Temp(:).gt.0.0)
   !   Temp_positive(:)=1.0
   !end where  

   if (Temp(1).gt.0.0) then ! si en surface temperature positive
    Temp_positive(1)=1.0
    do i=2, z_num
     if ((Temp(i).gt.0.0).and.(Temp_positive(i-1).gt.0.0)) then
       Temp_positive(i)=1.0
     endif
    enddo
   endif
   !write(*,*) 'Temp', Temp
   !write(*,*) 'Temp_positive', Temp_positive

   ! fin  de l'annee
   !write(*,*) 'end_year', end_year
   if (end_year.eq.1) then
   ALT_index=1
     altmax_lastyear=ALT
     if (Temp_positive(1).gt.0.0) then ! si en surface on a une temperature positive
      do k=1, z_num
       if (Temp_positive(k).gt.0.0) then ! tant que la temperature est positive on descend 
          ALT_index=k
       else
          exit ! A VERIFIER
       endif
      enddo
     endif
   ALT=D(ALT_index)
   !write(*,*) ' dans compute_alt : indice, ALT=', ALT_index , ALT
   !write(*,*) 'ialtmax_lastyear', altmax_lastyear
   endif

   !write(*,*) 'compteur_time_step', compteur_time_step
   !write(*,*) ' dans compute_alt : indice, ALT=', ALT_index , ALT
   !write(*,*) 'ialtmax_lastyear', altmax_lastyear

  end subroutine compute_alt
!--------------------------


!--------------------------
  subroutine carbon_redistribute(Temp, deepSOM_a, deepSOM_s, deepSOM_p, dz, ALT)
  ! initialize the carbon in the active layer by redistributing carbon 
  ! from litter and soil from vecode

    use parameter_mod, only: z_num   ! nombre d'horizons (couches)

    real, dimension(z_num),intent(in)              :: Temp
    real,                      intent(in)          :: ALT
    real, dimension(z_num), intent(inout)          :: dz ! epaisseur de chaque couche
    real                                           :: b3, b4 ! Carbon in litter and soil, for now fixed
    real,dimension(z_num),intent(inout)            :: deepSOM_a, deepSOM_s, deepSOM_p
    integer                                        :: k
    real                                           :: intdep !depth of integration
    real                                           :: z_lit !e folding depth 

!add input zfsoil
!define som_profile vecteur

    b3=50.! GtC ?
    b4=100.

    z_lit=0.2 ! arbitrary, peut aussi dependre de la prof des racines
    intdep = max(ALT, z_lit) !depth of integration pour l instant = ALT, peut dependre de
!    la profondeur des racines... et ne peut etre plus petit que z_lit

    ! initialise valeur carbon dans couche active

    ! si couche active plus petite que les 2 premires couches on repartit dans les 2 premieres couches
      !if (ALT.le.2) then
        deepSOM_a(1)=b3/2+b4/2
        deepSOM_a(2)=b3/2+b4/2
        !mettre les autres aussi
      !else

    !som_profile(:)=0.0
    !do il = 1, z_num ! sur la verticale
!        som_profile(il) = 1.0 / ( 1.0 - EXP( -intdep / z_lit ) ) * &
!                        ( EXP(-zf_soil(il-1)/z_lit)  - &
!                        EXP( -zf_soil(il)/z_lit ) )
    !enddo
    ! The above calculation should result in a som_profile of 0.9999999 which is
                   ! pretty good but not good enough to keep the mass balance
                   ! closed. Put the 
                   ! residual in the top layer (just because that layer should
                   ! always be present)
    !som_profile(1) = 1.0 - SUM(som_profile(2:z_num))
      !endif

   !write(*,*) 'deepSOM_a dans redistribute', deepSOM_a(:)


  end subroutine carbon_redistribute
!--------------------------


!--------------------------
  subroutine decomposition(Temp, zi_soil, dt, deepSOM_a, deepSOM_s, deepSOM_p, clay) !h_pori?
  !This routine calculates soil organic matter decomposition
  !decomposition depending on temperature, humidity, soil texture

    use parameter_mod, only: z_num

    !! carbon pools: indices
    INTEGER, PARAMETER :: iactive = 1      !! Index for active carbon pool (unitless)
    INTEGER, PARAMETER :: islow = 2        !! Index for slow carbon pool (unitless)
    INTEGER, PARAMETER :: ipassive = 3     !! Index for passive carbon pool (unitless)
    !INTEGER, PARAMETER :: isurface = 4     !! Index for passive carbon pool (unitless)
    !INTEGER, PARAMETER :: ncarb = 4        !! Number of soil carbon pools (unitless)
    INTEGER, PARAMETER :: ncarb = 3        !! Number of soil carbon pools (unitless)

    real, dimension(z_num),intent(in) :: Temp !, h_pori !, moist_in =humidity
    real, dimension(z_num)           :: temp_kelvin
    real, parameter                             :: q10 = 2.0
    real, parameter                             :: ZeroCelsius = 273.15
    logical, parameter                          :: limit_decomp_moisture = .false. ! to be changed to true if humidity in model
    real, dimension(z_num)           :: moistfunc_result
    REAL, DIMENSION(z_num), INTENT(in)          :: zi_soil ! profondeur de la couche = D dans VAMPER
    REAL                                        :: stomate_tau = 4.699E6
    REAL                                        :: depth_modifier = 1.E6             !! e-folding depth of turnover rates,following Koven et al.,2013,Biogeosciences. A very large value means no depth modification
!    REAL, DIMENSION(:), INTENT(in)           :: clay            !! clay content
    REAL , intent(in)                                       :: clay            !! clay content
    REAL, INTENT(in)                            :: dt          !! time step in seconds
    REAL                                        :: fbact_a, fbact_s, fbact_p
    REAL, DIMENSION(z_num)                      :: fbact                      !! turnover constant (day)
    REAL, dimension(z_num)                      :: tempfunc_result
    REAL, DIMENSION(ncarb,ncarb)                :: fc                         !! flux fractions within carbon pools
    REAL                                        :: dC
    REAL                                        :: fpassive = 1617.45                !! convertiing factor to go from active pool turnover to passive pool turnover from Guimberteau et al 2018 GMD
    REAL                                        :: fslow = 37.0                      !! convertiing factor to go from active pool turnover to slow pool turnover from Guimberteau et al 2018 GMD
    INTEGER                                     :: ij, il, iv, ip
    REAL                                        :: temp_local
    REAL, DIMENSION(z_num), intent(inout)       :: deepSOM_a, deepSOM_s, deepSOM_p
    REAL, DIMENSION(ncarb, ncarb)               :: somflux


! Donne le taux de decomposition fbact (temps de résidence)
    temp_kelvin(:)=Temp(:)+ZeroCelsius
! cutoff respiration when T < -1C
       WHERE (temp_kelvin(:) .GT. ZeroCelsius ) ! normal
          tempfunc_result(:) = EXP(log(q10) * ( temp_kelvin(:) - (ZeroCelsius+30.) ) / 10. )
       ELSEWHERE (temp_kelvin(:) .GT. ZeroCelsius - 1. )  ! linear dropoff to zero
          tempfunc_result(:) = (temp_kelvin(:) - (ZeroCelsius - 1.)) * &
               EXP( log(q10) * ( ZeroCelsius - (ZeroCelsius+30.) ) / 10. )
       ELSEWHERE  ! zero
          tempfunc_result(:) = 0.0
       ENDWHERE

       tempfunc_result(:) = MAX(MIN( 1.0, tempfunc_result(:) ), 0.0)

    IF ( limit_decomp_moisture ) THEN
       ! stomate moisture control function
!       moistfunc_result(:,:) = -1.1 * moist_in(:,:) * moist_in(:,:) + 2.4 * moist_in(:,:) - 0.29
!       moistfunc_result(:,:) = max( 0.25, min( 1.0, moistfunc_result(:,:) ) )
    ELSE
       moistfunc_result(:) = 1.0
    ENDIF

    DO ij = 1, z_num
      fbact(ij) = stomate_tau/(moistfunc_result(ij) * tempfunc_result(ij)) / EXP(-zi_soil(ij)/depth_modifier)
    ENDDO



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
             
!    ENDDO ! End loop over GridNoMax
            


  end subroutine decomposition
!--------------------------


!--------------------------
  subroutine cryoturbation(Temp, deepSOM_a, deepSOM_s, deepSOM_p, altmax_lastyear,max_cryoturb_alt,&
   min_cryoturb_alt, zi_soil, zf_soil, diff_k_const, bio_diff_k_const, dt, &
       bioturbation_depth)
  ! Vertical mixing by cryoturbation with a diffusion scheme

    use parameter_mod, only: z_num, YearType ! nombre de jours par an

    real, dimension(z_num), intent(in)    ::Temp !, h_pori
    REAL, DIMENSION(z_num), intent(inout) :: deepSOM_a, deepSOM_s, deepSOM_p
    real, intent(in)                      :: altmax_lastyear, max_cryoturb_alt, min_cryoturb_alt
    real                                  :: cryoturbation_depth
    logical                               :: cryoturb_location, bioturb_location
    REAL, DIMENSION(z_num), INTENT(in)    :: zi_soil ! !! depths of intermediate levels (m)! profondeur de la couche = D dans VAMPER
    REAL, DIMENSION(0:z_num),  INTENT(in) :: zf_soil        !! depths of full levels (m)
    REAL, DIMENSION(z_num)	          :: diff_k               !! Diffusion constant (m^2/s)
    REAL, INTENT(in)                         :: diff_k_const
    REAL, INTENT(in)                          :: bio_diff_k_const
    REAL                           :: mu_soil_rev
    REAL, INTENT(in)                          :: dt
    real, dimension(z_num) :: xc_cryoturb, xd_cryoturb 
    real, dimension(z_num) :: alpha_a, alpha_s, alpha_p, beta_a, beta_s, beta_p
    real :: xe_a, xe_s, xe_p
    real, intent(in) :: bioturbation_depth
    integer :: il


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
        !
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


! Calcul de la diffusion
 
       IF ( cryoturb_location .OR. bioturb_location )THEN 
                   ! 1. calculate the total soil organic matter 
                   !DO il = 1, ngrnd
                   !   altSOM_a_old(ip,iv,iele) = altSOM_a_old(ip,iv,iele) + deepSOM_a(ip,il,iv,iele)*(zf_soil(il)-zf_soil(il-1))
                   !   altSOM_s_old(ip,iv,iele) = altSOM_s_old(ip,iv,iele) + deepSOM_s(ip,il,iv,iele)*(zf_soil(il)-zf_soil(il-1))
                   !   altSOM_p_old(ip,iv,iele) = altSOM_p_old(ip,iv,iele) + deepSOM_p(ip,il,iv,iele)*(zf_soil(il)-zf_soil(il-1))
                   !ENDDO
                   
                   ! 2. diffuse the soil organic matter                 
                   deepSOM_a(1) = (deepSOM_a(1)+mu_soil_rev*beta_a(1)) / &
                        (1.+mu_soil_rev*(1.-alpha_a(1)))
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
                   !DO il = 1, ngrnd
                   !   altSOM_a(ip,iv,iele) = altSOM_a(ip,iv,iele) + deepSOM_a(ip,il,iv,iele)*(zf_soil(il)-zf_soil(il-1))
                   !   altSOM_s(ip,iv,iele) = altSOM_s(ip,iv,iele) + deepSOM_s(ip,il,iv,iele)*(zf_soil(il)-zf_soil(il-1))
                   !   altSOM_p(ip,iv,iele) = altSOM_p(ip,iv,iele) + deepSOM_p(ip,il,iv,iele)*(zf_soil(il)-zf_soil(il-1))
                   !ENDDO

                   
                 
                   ! 4. subtract the organic matter in the top layer(s) so that the total organic matter content of the active layer is conserved.             
                   ! for now remove this correction term...
!                   surfC_totake_a(ip,iv) = (altC_a(ip,iv)-altC_a_old(ip,iv))/(zf_soil(altmax_ind(ip,iv))-zf_soil(0))
!                   surfC_totake_s(ip,iv) = (altC_s(ip,iv)-altC_s_old(ip,iv))/(zf_soil(altmax_ind(ip,iv))-zf_soil(0))
!                   surfC_totake_p(ip,iv) = (altC_p(ip,iv)-altC_p_old(ip,iv))/(zf_soil(altmax_ind(ip,iv))-zf_soil(0))
!                   deepC_a(ip,1:altmax_ind(ip,iv),iv) = deepC_a(ip,1:altmax_ind(ip,iv),iv) - surfC_totake_a(ip,iv)
!                   deepC_s(ip,1:altmax_ind(ip,iv),iv) = deepC_s(ip,1:altmax_ind(ip,iv),iv) - surfC_totake_s(ip,iv)
!                   deepC_p(ip,1:altmax_ind(ip,iv),iv) = deepC_p(ip,1:altmax_ind(ip,iv),iv) - surfC_totake_p(ip,iv)
!
!                   ! if negative values appear, we don't subtract the delta-C from top layers
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



  end subroutine cryoturbation
!--------------------------

end module Carbon



