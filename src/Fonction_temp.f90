

module Fonction_temp

  use parameter_mod, only : C_ice,C_dry_soil,C_organic,C_water,freezing_range,Latent_Heat,rho_ice,z_num,K_fluids
  use parameter_mod, only :rho_organic,rho_soil,rho_water,K_other_minerals,K_organic,K_quartz,q_quartz,Bool_geometric,K_ice

  implicit none

  private
  public:: AppHeatCapacity, ThermalConductivity, Permafrost_Depth

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


  subroutine Permafrost_Depth(Temp,D,Per_depth)
     real, dimension(z_num), intent(in) :: Temp,D
     real, intent(out) :: Per_depth
     integer :: kk, z_m, z_p
     real :: alpha

     z_p = 0
     z_m = 0

     do kk=0,z_num-1

        if (Temp(z_num-kk)<0) then

           z_m=z_num-kk
           z_p=z_num-kk+1
           exit

        end if

     end do

     if ((z_p.GT.0).AND.(z_m.GT.0)) then

       alpha = Temp(z_p)-Temp(z_m)

       Per_depth = D(z_p) - Temp(z_p)/alpha

     else
       Per_depth = -9999.
     endif

  end subroutine Permafrost_Depth


end module Fonction_temp
