

module Fonction_temp

  use parameter_mod, only : C_ice,C_dry_soil,C_organic,C_water,freezing_range,Latent_Heat,rho_ice,z_num,K_fluids
  use parameter_mod, only :rho_organic,rho_soil,rho_water,K_other_minerals,K_organic,K_quartz,q_quartz,Bool_geometric,K_ice

  implicit none

  private
  public:: AppHeatCapacity, ThermalConductivity, diagnose_frost_Depth

  contains



  subroutine AppHeatCapacity(z_num, T, Tf, n, org_ind, Cp, porf, pori)
    !dmr Given Temperature and porosity (n), this computes a new Cp value and porf, pori on the vertical
    implicit none
    integer, intent(in) :: z_num, org_ind
    real, intent(in) :: Tf
    real, dimension(z_num), intent(in) :: T,n
    real, dimension(z_num), intent(out) :: Cp
    real, dimension(z_num), intent(out) ::  porf, pori
    integer :: kk
    real :: dTheta, theta, a, Csoil


    if (org_ind > 1) then

       Csoil = (1.0 - n(1)) * rho_organic * C_organic
    else
       Csoil = (1.0 - n(1)) * rho_soil * C_dry_soil
    end if
    do kk = 1, z_num

       if (kk <= 0) then
          Csoil=((1 - n(kk))* rho_organic * C_organic)
       else
          Csoil= 1.E6
       end if

       if (T(kk) < Tf) then

          a = - (((T(kk) - Tf) / freezing_range) ** 2.0)
          theta = exp(a)
          a = -2.0 / (freezing_range * freezing_range)
          dTheta = a * (T(kk) - Tf) * theta
          porf(kk) = n(kk) * theta
          pori(kk) = n(kk) - porf(kk)
          Cp(kk) = Csoil + (pori(kk) * C_ice * rho_ice) + (porf(kk) * C_water *rho_water)+(n(kk)*rho_water*Latent_Heat*dTheta)
       else
          porf(kk) = n(kk)
          pori(kk) = 0.0
          Cp(kk) = Csoil + (porf(kk) * C_water * rho_water)


       end if

    end do

  end subroutine AppHeatCapacity


  subroutine ThermalConductivity(layer , h_n, h_pori, h_porf, org_ind, Temp, Ther_cond )

    integer, intent(in) :: layer, org_ind
    real, intent(in) :: h_n, h_pori, h_porf, Temp
    real, intent(out) :: Ther_cond

    real :: Ksoil, Kice, Kfluids

    if (org_ind > layer) then

       Ksoil = K_organic

    else

       Ksoil = (K_quartz ** q_quartz) * (K_other_minerals ** (1.0-q_quartz))


    end if

    Kfluids = 0.1145 + 0.0016318 * (273.15 + Temp)

    Kice = 0.4865 + 488.19/(273.15 +Temp)

    if (Bool_geometric == 1) then

       Ther_cond = (Ksoil**(1-h_n) * Kice**(h_pori) * Kfluids**(h_porf))

    else

       Ther_cond = ((Ksoil**0.5)*(1-h_n) + (K_ice**0.5)*(h_pori) + (K_fluids**0.5)*(h_porf))**2

    end if

  end subroutine ThermalConductivity


  FUNCTION diagnose_frost_Depth(Temp,Depth_vals) RESULT(freeze_temp_max_min)

!~      use grids_more, only: undefined_value

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     REAL, DIMENSION(z_num), INTENT(in) :: Temp                ! Vertical temperature profile (should be °C or K ???)
     REAL, DIMENSION(z_num), INTENT(in) :: Depth_vals          ! Vertical depth values (m)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     REAL, DIMENSION(3)                 :: freeze_temp_max_min ! maximum and minimum depth with freezing conditions + active_layer_depth

     INTEGER :: indx_max, indx_min

     REAL, parameter           :: zero = 273.15
     LOGICAL, DIMENSION(z_num) :: mask_depth

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Initialization
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     indx_min = z_num
     indx_max = 1

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

    WHERE(Temp.LT.zero)
      mask_depth = .TRUE.
    ELSEWHERE
      mask_depth = .FALSE.
    ENDWHERE


    IF (ANY(mask_depth,dim=1)) then
                               ! There is at least one index where the temperature is below zero

      indx_max = FINDLOC(mask_depth,.TRUE.,DIM=1,BACK=.TRUE.)      ! find the deepest    index that matches freeze conditions
      indx_min = FINDLOC(mask_depth,.TRUE., DIM=1)                 ! find the shallowest index that matches freeze conditions

      IF (ALL(mask_depth,dim=1)) then ! Frozen everywhere
        freeze_temp_max_min(1) = Depth_vals(z_num)
        freeze_temp_max_min(2) = Depth_vals(1)
        freeze_temp_max_min(2) = 0.0

      ELSE                            ! not frozen everywhere

        if (indx_min.LT.indx_max) then ! There is an actual frozen layer

          if ((indx_min.GT.1).and.(indx_max.lt.z_num)) then ! general case: surface not frozen and not frozen to the bottom in the form F F F T T T T T T T F F F F
            freeze_temp_max_min(1) = interpol_value(Temp(indx_max:indx_max+1),Depth_vals(indx_max:indx_max+1))
            freeze_temp_max_min(2) = interpol_value(Temp(indx_min-1:indx_min),Depth_vals(indx_min-1:indx_min))
            freeze_temp_max_min(3) = freeze_temp_max_min(2)
          else if (mask_depth(1)) then ! surface is frozen
            freeze_temp_max_min(2) = Depth_vals(1)
            freeze_temp_max_min(3) = 0.0 ! no active layer
            if (mask_depth(z_num)) then ! bottom is frozen
              freeze_temp_max_min(1) = Depth_vals(z_num)
            else
              freeze_temp_max_min(1) = interpol_value(Temp(indx_max:indx_max+1),Depth_vals(indx_max:indx_max+1))
            endif

          else ! surface not frozen

            freeze_temp_max_min(2) = interpol_value(Temp(indx_min-1:indx_min),Depth_vals(indx_min-1:indx_min))
            freeze_temp_max_min(3) = freeze_temp_max_min(2)
            if (mask_depth(z_num)) then ! bottom is frozen
              freeze_temp_max_min(1) = Depth_vals(z_num)
            else
              freeze_temp_max_min(1) = interpol_value(Temp(indx_max:indx_max+1),Depth_vals(indx_max:indx_max+1))
            endif

          endif

        endif ! on not frozen everywhere

      ENDIF ! on frozen / not frozen everywhere

    ELSE ! Always above zero, nothing to do, no active layer per se
        freeze_temp_max_min(1:2) = Depth_vals(1)
        freeze_temp_max_min(3)   = 0.0
    ENDIF

!~     WRITE(*,*) "ALT et al.", freeze_temp_max_min, Temp(1), mask_depth(1)

  END FUNCTION diagnose_frost_Depth

  FUNCTION interpol_value(temp,Depth) result(depth_interpol)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

    REAL, DIMENSION(2), INTENT(in) :: temp, Depth

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

    REAL :: depth_interpol

    REAL    :: alpha, dist

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


    alpha = ABS(Temp(1))+ABS(Temp(2)) ! Temperature distance (°C) around the zero
    dist  = ABS(Depth(1)-Depth(2))     ! distance (m) between the two grid points

    if (Temp(1).GT.0.0) then
      depth_interpol = Depth(1) + Temp(1) * dist/alpha
    else
      depth_interpol = Depth(2) + Temp(2) * dist/alpha
    endif

  END FUNCTION interpol_value


end module Fonction_temp

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
