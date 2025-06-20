     MODULE glacial_index 

     IMPLICIT NONE
     
     PRIVATE
    
    
     character(len=256) :: filename = "Data_woelbroeck.csv"

     real, dimension(:), allocatable :: index_glagla      ! rescaled d18O
     real, dimension(:), allocatable :: dates_glagla      ! dates of rescaled d18O


     PUBLIC :: read_glacial_index, get_glacial_index

     CONTAINS
    
     SUBROUTINE read_glacial_index 


    ! test avec l'équation de Niu et al,. 2019 avec I pour O18 et temp 

    !Variables pour I (O18 et O16 )
    real :: OPD     ! o18 present day 
    real :: OMax   ! o18 le plus vieux  

    !Variables pour scale de température 
     
    real      :: t_inter  
    real      :: O18inter ! le O18 qui est calculé entre deux pas de temps 
    real      :: date_min, date_max, pas_index 
    real      :: indiceO18 !permet de creer un index pour aller chercher la valeur de O18
    integer   :: t_deb ! date de début de la simu 
    real      :: t_glacial ! coef de modification de la température 
    
    real      :: valeur_lue, date_lue
    
        ! To store data from file
     real, dimension(:), allocatable :: O18      ! o18 à l'instant t 
     real, dimension(:), allocatable :: datesO18  ! dates des variables 
    
    !Variables pour lire et calculer 
    
    integer :: unit_nb, ii, nb_lines = 0 , ll
    
    
    !Ici on lit les dataset en 1er O18 2eme Temperature  
    
    open(newunit= unit_nb, file=filename , status = "old" , action ="read" ) ! Importer data pour chemin 
      
      do
          read(unit_nb,*,iostat = ii) date_lue , valeur_lue
          if(ii/=0) exit
          nb_lines = nb_lines + 1
      end do

      allocate(datesO18(1:nb_lines))
      allocate(O18(1:nb_lines))
      allocate(index_glagla(1:nb_lines))
      allocate(dates_glagla(1:nb_lines))
      
      rewind(unit_nb)

      do ll=1,nb_lines
         read(unit_nb,*) date_lue , valeur_lue
         datesO18(ll) = date_lue 
         O18(ll) = valeur_lue
      end do 
      
      close(unit_nb) 
    
    
            
    !Ici on calcul I 
    OPD = O18(1) 
    OMax = MAXVAL(O18, dim=1)
        
        
    dates_glagla(:) = datesO18(:)
        
    ! Ici test calcul l'interpolation de O18  

    do ll=1,nb_lines 


      index_glagla(ll) = ((O18(ll) - OPD ) / (OMax - OPD))
      
    

        
    end do 
    
!~     do ll=1,nb_lines 
!~       WRITE(*,*) "Data_WALBROECK_SCALED",  dates_glagla(ll), index_glagla(ll)
!~     enddo
 
    END SUBROUTINE read_glacial_index

    FUNCTION get_glacial_index(date_kBP,nb_pasdetemps) result(g_indx)
    
    INTEGER, INTENT(in) :: date_kBP, nb_pasdetemps ! nb_pasdetemps doit etre en annees
    real :: pas_index, t_inter, A_inter, B_inter
    
    REAL, DIMENSION(:), ALLOCATABLE :: g_indx
    
    INTEGER :: indx_start_date, i, indx_end_date
    INTEGER, DIMENSION(2) :: indx_interpol
   
   
    ALLOCATE(g_indx(nb_pasdetemps))
   
    if ((date_kBP.GT.dates_glagla(UBOUND(dates_glagla,DIM=1))).OR.                    &
       ((date_kBP-nb_pasdetemps/1000.0).LT.dates_glagla(LBOUND(dates_glagla,DIM=1)))) &
       then
         WRITE(*,*) "INDEX GLAGLA OUT OF BOUNDS !!!"
         STOP
    endif
    
    indx_start_date = MINLOC(ABS(dates_glagla(:)-REAL(date_kBP)), DIM=1)
    indx_end_date = MINLOC(ABS(dates_glagla(:)                          &
                           -REAL(date_kBP-nb_pasdetemps/1000.0)), DIM=1)

    WRITE(*,*) "indx_start_date", indx_start_date, dates_glagla, date_kBP
   
    if ((dates_glagla(indx_start_date)-date_kBP).GE. 0.0) then ! data_kBP is more recent than dates_glagla
      indx_interpol(1) =  indx_start_date           ! indx right
      indx_interpol(2) = indx_start_date - 1        ! indx left                 T(2)----i-----T(1)
    else ! data_kBP is older than dates_glagla
      indx_interpol(1) =  indx_start_date + 1
      indx_interpol(2) = indx_start_date
    endif
      
    do i=1 , nb_pasdetemps
    
       t_inter = date_kBP - i/1000. - 1 ! nous sommes en BP !!!! 
       
       if (t_inter.LE.dates_glagla(indx_interpol(2))) then
         indx_interpol(:) = indx_interpol(:) - 1
       endif
       
       pas_index = dates_glagla(indx_interpol(1)) - dates_glagla(indx_interpol(2))
       
       
!~        g_indx(i) = index_glagla (indx_interpol(1)) + (( t_inter - dates_glagla(indx_interpol(1)))/ pas_index) &
!~                 * (index_glagla(indx_interpol(2)) -index_glagla(indx_interpol(1)))

       A_inter = ((index_glagla(indx_interpol(1))-index_glagla(indx_interpol(2)))/(pas_index))
       B_inter = index_glagla(indx_interpol(1)) - A_inter * dates_glagla(indx_interpol(1))
       g_indx(i) = t_inter * A_inter + B_inter
                                          
       write (*,*) t_inter, g_indx(i) 

    end do    
    
    
    END FUNCTION get_glacial_index



    END MODULE glacial_index 
