module Fonction_implicit


  use parameter_mod, only : z_num, rho_ice, Gfx, T_freeze,rho_snow_freeze,s_l_max
  use Fonction_temp, only : AppHeatCapacity, ThermalConductivity

  implicit none

  private
  public :: Implicit_T, Implicit_snow

  contains

    subroutine Implicit_snow(snw_dp,Tsnw,T_old,Tu,dt,dz,n,org_ind,Timp,Cp,Kp,Cp_snow,s_l)
!dmr [UNUSED]    subroutine Implicit_snow(snw_dp,rho_snow,Tsnw,T_old,Tu,dt,dz,n,org_ind,Timp,Cp,Kp,Cp_snow,s_l)
          
      integer, intent(in) :: org_ind,s_l
      real, intent(in) ::  snw_dp, dt,Tu
!dmr [UNUSED]      real, intent(in) ::  snw_dp, rho_snow, dt,Tu      
      real, dimension(z_num), intent(inout) :: T_old, n, dz
      real, dimension(s_l_max), intent(inout) :: Tsnw
      real, intent(out) :: Cp_snow
      real, dimension(z_num), intent(out) :: Timp, Cp
      real, dimension(z_num), intent(out) :: Kp
    
      ! Variables locales !
     
      real, dimension(z_num+s_l) ::  Kp_s, T_last, n_s, dz_s, K_s
      real, dimension(:), allocatable :: Cp_s, porf, pori
      real, dimension(z_num) ::Cp_temp, porf_temp, pori_temp !dmr [UNUSED], T_new
      real, dimension(z_num+s_l) :: T_iter
      real, dimension(z_num+s_l,z_num+s_l) :: MM
      real, dimension(z_num+s_l) :: DD
      real, dimension(z_num+s_l-1) :: DL, DU
      real, dimension(z_num+s_l) :: Knows
!dmr [UNUSED]      real, dimension(s_l) :: D
      real :: m_Gfx, A, B, C, Z1
      integer :: kk, ll, z_s
!dmr [UNUSED]      integer :: kk, ll, z_s,ind_snw      

!dmr [UNUSED]      integer, dimension(z_num+s_l) :: IPIV
      integer :: info_dgesv


      !write(*,*) s_l
    
      allocate(porf(1:z_num+s_l))
      allocate(pori(1:z_num+s_l))
      allocate(Cp_s(1:z_num+s_l))

!      write(*,*) Tsnw
!      do ll=1,s_l
!         Tsnw(1:s_l_max-s_l) = Tu
!      end do
!      write(*,*) Tsnw 
!      do ll = s_l,s_l_max
!         T_last(s_l-ll+1) = Tsnw(ll)
!      end do
        
      
      

      z_s = z_num+s_l
      m_Gfx = Gfx/1000.0


      T_last(s_l+1:z_s) = T_old(1:z_num)
      dz_s(s_l+1:z_s) = dz(1:z_num)
     
      n_s(s_l+1:z_s) = n(1:z_num)
      
      do kk = 1,s_l
         Cp_s(kk) = (2.108*1000000.0)*rho_snow_freeze/rho_ice
         K_s(kk) = 2.9*(rho_snow_freeze**2)*0.000001
         !Cp_s(kk) = 0.42*(10**6)
         !K_s(kk) = 0.011
         dz_s(kk) = snw_dp
         T_last(kk) = Tsnw(kk)
      end do
      
      write(*,*) "[FONC_IMP] K_s, Cp_s: ", K_s(1), Cp_s(1)
      
      !write(*,*)"T_last=", T_last
     
      
      do kk=1,5 
         
         
         MM(1:z_s,1:z_s)=0
         DD(1:z_s)=0
         DL(1:z_s-1)=0
         DU(1:z_s-1)=0
         
         if (kk==1) then
            
            T_iter(1:z_s) = T_last(1:z_s)
            
         else

            T_iter(1:z_s)=0.5*(T_iter(1:z_s)+T_last(1:z_s))

         end if

         
         call AppHeatCapacity(z_num,T_iter(s_l+1:z_s),T_freeze,n_s(s_l+1:z_s),org_ind+s_l,Cp_temp,porf_temp,pori_temp)
         
         
         do ll=s_l+1,z_s
            
            call ThermalConductivity(ll,n_s(ll),pori_temp(ll-1),porf_temp(ll-1),org_ind+s_l,T_iter(ll),K_s(ll))
            Cp_s(ll) = Cp_temp(ll-1)

         end do
         
         !write(*,*) Cp_s
         !write(*,*) Cp_s
         !write(*,*) dz_s
         !write(*,*) K_s

         Kp_s(1:z_s-1) = (K_s(1:z_s-1)+K_s(2:z_s))*0.5
         !Kp_s(1:z_s-1) = K_s(1:z_s-1)

         Knows(1) = Tu
         
         do ll=2,z_s-1
            
            Z1 = T_last(ll)
            
            A=(dt/((dz_s(ll-1)+dz_s(ll))*0.5*dz_s(ll)) * (Kp_s(ll-1)/Cp_s(ll-1)))
            B=(dt/((dz_s(ll+1)+dz_s(ll))*0.5*dz_s(ll)) * (Kp_s(ll)/Cp_s(ll)))

            A=(dt/((dz_s(ll-1))*dz_s(ll)) * (Kp_s(ll-1)/Cp_s(ll-1)))
            B=(dt/((dz_s(ll+1))*dz_s(ll)) * (Kp_s(ll)/Cp_s(ll)))
            C= 1+A+B
            !write(*,*) 1-A-B, ll

            MM(ll,ll-1) = -A
            MM(ll,ll) = C
            MM(ll,ll+1) = -B
            DL(ll-1) = -A
            DU(ll) = -B
            DD(ll) = C

            Knows(ll) = Z1

         end do
           
         A=(dt/((dz_s(z_s-1)+dz_s(z_s))*0.5*dz_s(z_s)) * (Kp_s(z_s-1)/Cp_s(z_s)))
         C=1.0+A
         MM(z_s,z_s-1)=-A
         MM(z_s,z_s)=C
         MM(1,1) = 1
         DL(z_s-1) = -A
         DD(z_s) = C
         DD(1) = 1

         Knows(z_s) = T_last(z_s) + ((dt/Cp_s(z_s))*m_Gfx/dz_s(z_s-1))
         
         

         !call sgesv(z_s,1,MM,z_s,IPIV,Knows,z_s,info_dgesv)
         call sgtsv(z_s,1,DL,DD,DU,Knows,z_s,info_dgesv)
         
         T_iter(1:z_s) = Knows(1:z_s)
         !write(*,*) T_iter
         !write(*,*) info_dgesv
         
      end do

      !write(*,*) Kp_s
      !write(*,*) Cp_s

      Timp(1:z_num) = T_iter(s_l+1:z_s)
      Cp(1:z_num) = Cp_temp(1:z_num)
      Kp(1:z_num) = K_s(s_l+1:z_s)
      Cp_snow = Cp_s(1)
      Tsnw(1:s_l) = T_iter(1:s_l)

      deallocate(porf)
      deallocate(pori)
      deallocate(Cp_s)


    end subroutine Implicit_snow


    subroutine Implicit_T(T_old,Tu,Tb,dt,dz,n,org_ind,Timp,Kp)
!dmr [UNUSED]   subroutine Implicit_T(T_old,Tu,Tb,dt,dz,n,org_ind,Timp,Cp,Kp)    
      
      integer, intent(in) :: org_ind
      real, intent(in) :: dt,Tu,Tb
      real, dimension(:), intent(in) :: T_old, n, dz
      real, dimension(z_num), intent(out) :: Timp, Kp
!dmr [UNUSED]      real, dimension(z_num), intent(out) :: Cp
      
      real, dimension(z_num) :: pori, porf, Cp_temp
      real, dimension(z_num,z_num) :: MM
      real, dimension(z_num) :: Knows
      real :: m_Gfx, A, B, C, Z1
      integer :: kk, ll
      real, dimension(z_num) :: T_last, Kp_s
!dmr [UNUSED]      real, dimension(z_num) :: T_last, T_new, Kp_s      
      real, dimension(z_num) :: T_iter
      real, dimension(z_num) :: DD
      real, dimension(z_num-1) :: DL, DU
!dmr [UNUSED]      integer, dimension(z_num) :: IPIV
      integer :: info_dgesv
      
      m_Gfx = gfx/1000.0

      T_last(1:z_num) = T_old(1:z_num)
 
      do kk=1,5 
         
         MM(1:z_num,1:z_num) = 0
         DD(1:z_num)=0
         DL(1:z_num-1)=0
         DU(1:z_num-1)=0

         if (kk==1) then
            
            T_iter(1:z_num) = T_last(1:z_num)
            
         else

            T_iter(1:z_num) = 0.5*(T_iter(1:z_num)+T_last(1:z_num))

         end if

         call AppHeatCapacity(z_num,T_iter,T_freeze,n,org_ind,Cp_temp,porf,pori)

         do ll=1,z_num-1

            call ThermalConductivity(ll,n(ll),pori(ll),porf(ll),org_ind,T_iter(ll),Kp(ll))
            Kp(z_num) = 2
            !Kp(ll) = 2
            !Cp_temp(ll) = Cp_temp(ll) 

         end do


         Kp_s(1:z_num-1) = (Kp(1:z_num-1)+Kp(2:z_num))*0.5
         !Kp_s(1:z_num-1) = Kp(1:z_num-1)


         do ll=2,z_num-1
            
            Z1 = T_last(ll)
            
            A=(dt/((dz(ll-1)+dz(ll))*0.5*dz(ll)) * (Kp_s(ll-1)/Cp_temp(ll)))
            B=(dt/((dz(ll+1)+dz(ll))*0.5*dz(ll)) * (Kp_s(ll)/Cp_temp(ll)))

            !A=(dt/((dz(ll-1))*dz(ll)) * (Kp_s(ll-1)/Cp_temp(ll)))
            !B=(dt/((dz(ll+1))*dz(ll)) * (Kp_s(ll)/Cp_temp(ll)))
            C= 1+A+B
            !write(*,*) 1-A-B , ll

            MM(ll,ll-1) = -A
            MM(ll,ll) = C
            MM(ll,ll+1) = -B
            DL(ll-1) = -A
            DU(ll) = -B
            DD(ll) = C
            Knows(ll) = Z1

         end do
         
         A=(dt/((dz(z_num-1)+dz(z_num))*0.5*dz(z_num-1)) * (Kp_s(z_num-1)/Cp_temp(z_num-1)))
         C=1.0+A
         Knows(1) = Tu
         Knows(z_num)=Tb 
         MM(z_num,z_num)=1
         MM(1,1)=1
         DD(1) = 1
         DD(z_num) = 1+A
         DL(z_num-1) = -A

         


         
         
         !MM(z_num,z_num-1)=-A
         !MM(z_num,z_num)=C
         !MM(1,1) = 1
         
         !DD(z_num) = C
         !DD(1) = 1

         !Knows(z_num) = T_last(z_num) + ((dt/Cp_temp(z_num))*m_Gfx/dz(z_num-1))

         !call sgesv(z_num,1,MM,z_num,IPIV,Knows,z_num,info_dgesv) 
         call sgtsv(z_num,1,DL,DD,DU,Knows,z_num,info_dgesv) 
         T_iter(1:z_num) = Knows(1:z_num)
         
      end do
      
   
      Timp(1:z_num) = T_iter(1:z_num)
      !write(*,*) dz(z_num),dz(z_num-1)

    end subroutine Implicit_T


end module Fonction_implicit
