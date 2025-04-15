

module Model_snow

  use parameter_mod, only : Latent_heat, rho_water, rho_snow_freeze, gravity, rho_snow_fresh

  implicit none

  private
  
  public :: snw_average_swe, snw_average_snw, snw_average_snw_tot, snw_proc
  
  contains

    
    subroutine snw_average_swe(swe_f, swe_tot, snw_tot, rho_snow)
      
      real, intent(in) :: swe_f
      real, intent(inout) :: swe_tot, snw_tot, rho_snow
      real :: snw_f
      
      
      if (swe_f > 0.000001) then

         snw_f = swe_f * (rho_water/rho_snow_freeze)
         snw_tot = snw_tot + snw_f
         swe_tot = swe_tot + swe_f

         if (snw_tot > 0.000001) then
            
            rho_snow = swe_tot*rho_water/snw_tot

         else
            
            rho_snow = rho_snow_freeze

         end if

      end if


    end subroutine snw_average_swe


    subroutine snw_average_snw(snw_f, swe_tot, snw_tot, rho_snow,swe_f)
      
      real, intent(in) :: snw_f
      real, intent(inout) :: swe_tot, snw_tot, rho_snow, swe_f

      
      
      if (snw_f > 0.000001) then

         swe_f = snw_f * (rho_snow_freeze/rho_water)
         snw_tot = snw_tot + snw_f
         swe_tot = swe_tot + swe_f

         if (snw_tot > 0.000001) then
            
            rho_snow = swe_tot*rho_water/snw_tot

         else
            
            rho_snow = rho_snow_freeze

         end if

      end if


    end subroutine snw_average_snw


    subroutine snw_average_snw_tot(rho_snow)
      
      real, intent(out) :: rho_snow
            
      rho_snow = 10.0*36*2.3


    end subroutine snw_average_snw_tot




!~     subroutine snw_average_snw_tot(snw_tot,snw_tot_old, rho_snow,swe_tot,dt)
      
!~       real, intent(in) :: snw_tot, snw_tot_old,dt
!~       real, intent(inout) :: rho_snow,swe_tot
!~       real :: swe_f

      
!~       rho_snow = 10.0*36*2.3


!~     end subroutine snw_average_snw_tot
 




    subroutine snw_proc(Tsnw, snw_tot, swe_tot, frac_snw, Cp_snow, rho_snow, dt)
!dmr [UNUSED]    subroutine snw_proc(Tsnw,T_soil, snw_tot, swe_tot, frac_snw, Cp_snow, rho_snow, dt)    

      real, intent(inout) :: Tsnw, snw_tot, swe_tot, frac_snw, rho_snow
!dmr [UNUSED]      real, intent(inout) :: Tsnw,T_soil, snw_tot, swe_tot, frac_snw, rho_snow      
      real, intent(in) :: Cp_snow, dt
      real :: frac_snw_old, rho_snow_old, swe_melt, H_1, H_2, N, var, rho_snow_new

      swe_melt = 0.0
      frac_snw_old = frac_snw
      rho_snow_old = rho_snow

      if ( Tsnw >= 0.0 ) then
         
         H_1 = (Tsnw)*Cp_snow*snw_tot
         H_2 = rho_water * Latent_heat * swe_tot
         !write(*,*) H_1-H_2
         if (H_1<=H_2) then
            
            frac_snw = 1-(H_1/H_2)
            !write(*,*) snw_tot
            snw_tot = frac_snw * snw_tot
            write(*,*) "[MOD_SNW] frac_snw, Tsnw, snw_tot: ", frac_snw, Tsnw, snw_tot
            Tsnw = 0.0
            !write(*,*) "H1<H2"
         else
            
            frac_snw = 0.0
            snw_tot = 0.0
            swe_tot = 0.0
            Tsnw = 0.0
         end if
        
      else
         
         frac_snw = frac_snw_old
         
      end if
         
      N = 0.5*snw_tot*rho_snow*5.0*0.00000001*rho_snow*gravity
      var = 14.643 - (4000.0/min(Tsnw+273.16,273.16))-0.02*rho_snow
      rho_snow_new = rho_snow + N*dt*exp(var)
      rho_snow = rho_snow_new

      snw_tot = (rho_snow_old/rho_snow_new)*snw_tot
      
    end subroutine snw_proc

   
end module Model_snow
