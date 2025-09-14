!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2024 Didier M. Roche (a.k.a. dmr)

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

    MODULE spatialvars_mod

#include "constant.h"

    USE timer_mod, only: cell_time

    IMPLICIT NONE

    PRIVATE


            ! SPATIAL GLOBAL VARIABLES

     real, dimension(:,:),allocatable         :: Temp      & !dmr [SPAT_VAR], soil temperature over the vertical // prognostic
                                                ,Kp        & !dmr [CNTST]     heat conductivity constant over the depth, current value is 2
                                                ,n         & !dmr [SPAT_VAR], porosity on the vertical
                                                ,Cp        & !dmr [SPAT_VAR]  specific heat capacity
                                                ,pori      & !dmr [???  TBC]
                                                ,porf        !dmr [???  TBC]

     real, dimension(:), allocatable          ::  GeoHFlux    &
                                                , Tinit_SV    &
                                                , T_bottom_SV

     real, dimension(:,:), allocatable :: freeze_depth_SV

     integer, dimension(:), allocatable :: orgalayer_indx

            ! CLIMATE FORCING VARIABLES
     REAL, DIMENSION(:,:), ALLOCATABLE :: forcing_surface_temp ! two dimensions / spatial and time in that order
#if ( SNOW_EFFECT == 1 )
     REAL, DIMENSION(:,:), ALLOCATABLE :: forcing_snow_thick
#endif

            ! could add swe_f_t, snw_dp_t,rho_snow_t,T_snw
     REAL, DIMENSION(:,:), ALLOCATABLE :: restart_temperature  ! VERT/SPAT

     real, dimension(:)   ,  allocatable :: ALT_SV, altmax_ly_SV


#if ( CARBON == 1 )
     real,dimension(:,:)  ,  allocatable         :: deepSOM_a & !dmr [TBD]
                                                  , deepSOM_s & !dmr [TBD]
                                                  , deepSOM_p & !dmr [TBD]
                                                  , deepSOM
     real,dimension(:,:,:),  allocatable         :: fc_SV       !dmr [TBD]

     ! Timer variable carbon only
     logical, dimension(:), allocatable          :: end_year_SV          !=0 if not end of year, =1 if end of year
     real, dimension(:)   , allocatable          :: b3_SV, b4_SV
#endif

#if ( SNOW_EFFECT == 1 )
     real, dimension(:,:) ,allocatable           :: Temp_snow           !dmr [VERTCL, SPAT_VAR] temperature in snow layers               [C]
     real, dimension(:)   ,allocatable           :: depth_snow_layer    !dmr [SPAT_VAR]         depth of snow in the snow layers         [m]
     integer, dimension(:),allocatable           :: nb_snow_layer       !dmr [SPAT_VAR]         number of snow layers from discretization [1]
#endif



!   Main Timer variable

           !> number of timesteps, can be months or days. Spatial variable for parallel execution
     TYPE(cell_time), dimension(:),  allocatable :: compteur_tstep_SV



     PUBLIC:: spatialvars_allocate, spatialvars_init, UPDATE_climate_forcing, DO_spatialvars_step, SET_coupled_climate_forcing

     CONTAINS







!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Allocation of two dimensional variables (VERTCL, SPAT_VAR)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     SUBROUTINE spatialvars_allocate ! VERTCL, SPAT_VAR

       use parameter_mod, only: gridNoMax, timFNoMax, z_num
#if ( CARBON == 1 )
       use carbon       , only: ncarb
#endif
#if ( SNOW_EFFECT == 1)
       use simple_snow, only: max_nb_snow_layers
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       allocate(Temp(1:z_num,1:gridNoMax)) !dmr SPAT_VAR
       allocate(Kp(1:z_num-1,1:gridNoMax)) !dmr SPAT_VAR
       allocate(n(1:z_num,1:gridNoMax))    !dmr SPAT_VAR
       allocate(Cp(1:z_num,1:gridNoMax))   !dmr SPAT_VAR
       allocate(pori(1:z_num,1:gridNoMax)) !dmr SPAT_VAR
       allocate(porf(1:z_num,1:gridNoMax)) !dmr SPAT_VAR

       allocate(GeoHFlux(1:gridNoMax))
       allocate(Tinit_SV(1:gridNoMax))
       allocate(T_bottom_SV(1:gridNoMax))
       allocate(freeze_depth_SV(3,1:gridNoMax))

       allocate(orgalayer_indx(1:gridNoMax))

       allocate(ALT_SV(1:gridNoMax))

#if ( CARBON == 1 )
                        !nb and mbv Carbon cycle
       allocate(deepSOM_a(1:z_num,1:gridNoMax))
       allocate(deepSOM_s(1:z_num,1:gridNoMax))
       allocate(deepSOM_p(1:z_num,1:gridNoMax))
       allocate(deepSOM(1:z_num,1:gridNoMax))
!~        allocate(temp_oncepositive(1:z_num,1:gridNoMax))
       allocate(fc_SV(1:ncarb,1:ncarb,1:gridNoMax))
!~        allocate(clay_SV(1:gridNoMax))
       allocate(altmax_ly_SV(1:gridNoMax))
       allocate(compteur_tstep_SV(1:gridNoMax))
       allocate(end_year_SV(1:gridNoMax))
       allocate(b3_SV(1:gridNoMax))
       allocate(b4_SV(1:gridNoMax))

#endif

#if ( SNOW_EFFECT == 1)
       allocate(Temp_snow(1:max_nb_snow_layers,1:gridNoMax))    !dmr [VERTCL, SPAT_VAR] temperature in snow layers  [C]
       allocate(depth_snow_layer(1:gridNoMax))     !dmr [SPAT_VAR]         depth of snow in the snow layers         [m]
       allocate(nb_snow_layer(1:gridNoMax))        !dmr [SPAT_VAR]         nuber of snow layers from discretization [1]
#endif

       allocate(forcing_surface_temp(1:gridNoMax,1:timFNoMax))
       allocate(restart_temperature(1:z_num,1:gridNoMax)) ! contains the restart temperature at init

#if ( SNOW_EFFECT == 1)
       allocate(forcing_snow_thick(1:gridNoMax,1:timFNoMax))
#endif


     END SUBROUTINE spatialvars_allocate


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Allocation of two dimensional variables (VERTCL, SPAT_VAR)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     SUBROUTINE spatialvars_init ! VERTCL, SPAT_VAR

       use parameter_mod,  only: gridNoMax, z_num ! , timFNoMax
       use vertclvars_mod, only: vertclvars_init

           ! Temporary addendum [2025-04-16]
       use parameter_mod,  only: Gfx, T_init
       use parameter_mod,  only: forc_tas_file, name_tas_variable

       use parameter_mod,  only: GHF_spatial_file, GHF_variable_name, Tinit_spatial_file, Tinit_variable_name

!~        use parameter_mod,  only: Forcage_Month_day

       use timer_mod,      only: init_time_cell

#if ( CARBON == 1 )
       use carbon        , only: carbon_init
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       integer :: gridp
       logical :: logic_month_day
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        if (DAILY.EQ.1) then
            logic_month_day = .TRUE.
        else
            logic_month_day = .FALSE.
        endif

#if ( OFFLINE_RUN == 1 )
        call get_clim_forcing(forc_tas_file, name_tas_variable,forcing_surface_temp)
        call fix_Kelvin_or_Celsius(forcing_surface_temp)
!dmr --- Would need to repeat the work for get_clim_forcing to get the snow thickness ...
!dmr [TBD] another day ... 2025-07-08
#if ( SNOW_EFFECT == 1 )
        forcing_snow_thick(:,:) = 1.5
#endif

#endif



#if (SP_GHF == 1)
        GeoHFlux(1:gridNoMax) = get_Spatial_2Dforcing(GHF_spatial_file,GHF_variable_name)
#else
        GeoHFlux(:) = Gfx
#endif

#if (SP_Tinit == 1)
        Tinit_SV(:) = get_Spatial_2Dforcing(Tinit_spatial_file,Tinit_variable_name)
#else
        !Tinit_SV(:) = T_init
        Tinit_SV(:) = SUM(forcing_surface_temp(:,:),DIM=2)/UBOUND(forcing_surface_temp,DIM=2)
#endif

            !dmr Initialization of all columns, one by one
        do gridp = 1, gridNoMax
          call vertclvars_init(GeoHFlux(gridp), Tinit_SV(gridp), Kp(:,gridp),Cp(:,gridp), orgalayer_indx(gridp), n(:,gridp) &
                             , Temp(:,gridp), T_bottom_SV(gridp))

#if ( CARBON == 1 )
            !dmr Initialization of all columns, one by one
          call carbon_init(deepSOM_a(:,gridp), deepSOM_s(:,gridp), deepSOM_p(:,gridp), fc_SV(:,:,gridp),ALT_SV(gridp)       &
                          ,b3_SV(gridp), b4_SV(gridp))
#endif

          compteur_tstep_SV(gridp) = init_time_cell(0,.FALSE.,.FALSE.,logic_month_day)


        enddo

#if ( SNOW_EFFECT == 1 )
        Temp_snow(:,:) = 0.0 !dmr simplest possible init without data ... To be fixed!
        depth_snow_layer(:) = 3.0
        nb_snow_layer(:) = 0

#endif


     END SUBROUTINE spatialvars_init

     SUBROUTINE fix_Kelvin_or_Celsius(Temp_field)

      use parameter_mod, only: tK_zero_C

      REAL, DIMENSION(:,:), INTENT(inout) :: Temp_field

      if (MINVAL(Temp_field).GT.100.0) then
        Temp_field(:,:) = Temp_field(:,:) - tK_zero_C
      endif

     END SUBROUTINE fix_Kelvin_or_Celsius


     FUNCTION get_Spatial_2Dforcing(nc_file_to_read,name_surf_variable) result(ReadInit_SV)

        use netcdf
        use parameter_mod, only: gridNoMax

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        CHARACTER(len=*), INTENT(in):: nc_file_to_read, name_surf_variable

        REAL, DIMENSION(1:gridNoMax) :: ReadInit_SV

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        INTEGER :: ncid, ret_stat, nDims, nvars, nGlobalAtts, unlimdimid
        INTEGER :: d,v,a, i, j

        CHARACTER(len=NF90_MAX_NAME), DIMENSION(:),   ALLOCATABLE :: dimNAMES
        INTEGER                     , DIMENSION(:),   ALLOCATABLE :: dimLEN

        CHARACTER(len=NF90_MAX_NAME), DIMENSION(:),   ALLOCATABLE :: varNAMES
        INTEGER                     , DIMENSION(:),   ALLOCATABLE :: varXTYPE, varNDIMS, varNATTS
        INTEGER                     , DIMENSION(:,:), ALLOCATABLE :: varDIMIDS

        CHARACTER(len=NF90_MAX_NAME), DIMENSION(:),   ALLOCATABLE :: attNAMES

        LOGICAL                     , DIMENSION(:),   ALLOCATABLE :: varisDIM

        INTEGER ::  dim1, dim2, varunmasked ! maskVarID=0,

        REAL                      , DIMENSION(:,:),   ALLOCATABLE :: varDATA
        REAL                                                      :: var_undef = 0.0

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!      dmr GETTING THE VARIABLE NEEDED

      ret_stat = nf90_open(nc_file_to_read, nf90_nowrite, ncid)
      ret_stat = nf90_inquire(ncid, nDims, nVars, nGlobalAtts, unlimdimid)

      ! I look forward to get a grid with one spatial variable and two dimensions for now
      ! Hence I should get nDims = 3, nVars = 4 at least (could be more if more variables), unlimited without time, hence unlimdimid == -1

!~       if ( (nDims.EQ.2).AND.(nVars.GE.3).AND.(unlimdimid.EQ.-1) ) then ! valid file, no time expected
      if ( (nDims.GE.2).AND.(nVars.GE.3).AND.(unlimdimid.EQ.-1) ) then ! valid file, no time expected

       ALLOCATE(dimNAMES(nDims))
       ALLOCATE(dimLEN(nDims))

       DO d=1,nDims
         ret_stat = nf90_inquire_dimension(ncid, d, dimNAMES(d), dimLEN(d))
!~          WRITE(*,*) "Found dimensions ::", dimNAMES(d), d
       ENDDO

       ! really only need the length of the unlimited (time) variable

       ret_stat = nf90_inq_varid(ncid,name_surf_variable,v)

       ALLOCATE(varNAMES(nVars))
       ALLOCATE(varXTYPE(nVars))
       ALLOCATE(varNDIMS(nVars))
       ALLOCATE(varDIMIDS(nVars,nDims))
       ALLOCATE(varNATTS(nVars))
       ALLOCATE(varisDIM(nVars))

       ret_stat = nf90_inquire_variable(ncid, v, varNAMES(v), varXTYPE(v), varNDIMS(v), varDIMIDS(v,:), varNATTS(v))

      else
         WRITE(*,*) "netCDF file for VAR does not match expectations", name_surf_variable
         STOP
      endif

      if (v.NE.0) then ! I have found my variable to read, I have enough to define it

      if (varNDIMS(v).NE.2) then
         WRITE(*,*) "Current version only support 2D var file for GHF"
         STOP
      else
         ! Get the actual values of the thing

         dim1 = dimLEN(varDIMIDS(v,1))
         dim2 = dimLEN(varDIMIDS(v,2))

         ALLOCATE(varDATA(dim1,dim2))

         ret_stat = NF90_GET_VAR(ncid,v,varDATA)

!~          WRITE(*,*) "DIMS", dim1, dim2, dim3 ! lon lat time is getting out ... correct? (somehow netCDF reads backward)

         ALLOCATE(attNAMES(varNATTS(v)))

         do a=1,varNATTS(v)
           ret_stat = NF90_INQ_ATTNAME(ncid, v, a, attNAMES(a))
           if ((INDEX(attNAMES(a), "missing_value").GT.0) .OR. (INDEX(attNAMES(a), "_FillValue").GT.0) ) then
               ret_stat = NF90_GET_ATT(ncid, v, attNAMES(a), var_undef)
           endif
         enddo

!~          ! Get the missing value of the thing if exists ...
!~          WRITE(*,*) "Value for undef :: ", var_undef
      endif

      ! variable is read in ... now need to check where I have an actual value (not masked)


      varunmasked = 0

      DO j=1, dim2
      DO i=1, dim1
         if (varDATA(i,j).NE.var_undef) then
             ! count it in !
            varunmasked = varunmasked + 1
            ReadInit_SV(varunmasked) = varDATA(i,j)
         endif
      ENDDO
      ENDDO

      WRITE(*,*) "READ Variable from netCDF file: ", name_surf_variable, MINVAL(ReadInit_SV), MAXVAL(ReadInit_SV)

      else
        WRITE(*,*) "var id is zero", v
      endif


      DEALLOCATE(dimNAMES)
      DEALLOCATE(dimLEN)
      DEALLOCATE(varNAMES)
      DEALLOCATE(varXTYPE)
      DEALLOCATE(varNDIMS)
      DEALLOCATE(varDIMIDS)
      DEALLOCATE(varNATTS)
      DEALLOCATE(varisDIM)
      DEALLOCATE(varDATA)



     END FUNCTION get_Spatial_2Dforcing



     SUBROUTINE get_clim_forcing(forc_surf_file, name_surf_variable, forcing_surface_var)

        use netcdf
        use parameter_mod, only: str_len

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        CHARACTER(len=str_len), intent(in) :: forc_surf_file, name_surf_variable

        REAL, DIMENSION(:,:),   intent(out):: forcing_surface_var ! two dimensions / spatial and time in that order

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        INTEGER :: ncid, ret_stat, nDims, nvars, nGlobalAtts, unlimdimid
        INTEGER :: d,v,a, i, j

        CHARACTER(len=NF90_MAX_NAME), DIMENSION(:),   ALLOCATABLE :: dimNAMES
        INTEGER                     , DIMENSION(:),   ALLOCATABLE :: dimLEN

        CHARACTER(len=NF90_MAX_NAME), DIMENSION(:),   ALLOCATABLE :: varNAMES
        INTEGER                     , DIMENSION(:),   ALLOCATABLE :: varXTYPE, varNDIMS, varNATTS
        INTEGER                     , DIMENSION(:,:), ALLOCATABLE :: varDIMIDS

        CHARACTER(len=NF90_MAX_NAME), DIMENSION(:),   ALLOCATABLE :: attNAMES

        LOGICAL                     , DIMENSION(:),   ALLOCATABLE :: varisDIM

        INTEGER :: dim1, dim2, dim3, varunmasked ! maskVarID=0,

        REAL                      , DIMENSION(:,:,:), ALLOCATABLE :: varDATA
        REAL                                                      :: var_undef = 0.0

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!      dmr GETTING THE VARIABLE NEEDED

      ret_stat = nf90_open(forc_surf_file, nf90_nowrite, ncid)
      ret_stat = nf90_inquire(ncid, nDims, nVars, nGlobalAtts, unlimdimid)

      ! I look forward to get a grid with one spatial variable and two dimensions for now
      ! Hence I should get nDims = 3, nVars = 4 at least (could be more if more variables), unlimited with time, hence unlimdimid != -1

      if ( (nDims.EQ.3).AND.(nVars.GE.4).AND.(unlimdimid.NE.-1) ) then ! valid file

       ALLOCATE(dimNAMES(nDims))
       ALLOCATE(dimLEN(nDims))

       DO d=1,nDims
         ret_stat = nf90_inquire_dimension(ncid, d, dimNAMES(d), dimLEN(d))
!~          WRITE(*,*) "Found dimensions ::", dimNAMES(d), d
       ENDDO

       ! really only need the length of the unlimited (time) variable

       ret_stat = nf90_inq_varid(ncid,name_surf_variable,v)

       ALLOCATE(varNAMES(nVars))
       ALLOCATE(varXTYPE(nVars))
       ALLOCATE(varNDIMS(nVars))
       ALLOCATE(varDIMIDS(nVars,nDims))
       ALLOCATE(varNATTS(nVars))
       ALLOCATE(varisDIM(nVars))

       ret_stat = nf90_inquire_variable(ncid, v, varNAMES(v), varXTYPE(v), varNDIMS(v), varDIMIDS(v,:), varNATTS(v))

      else
         WRITE(*,*) "netCDF file for VAR does not match expectations", name_surf_variable
         STOP
      endif

      if (v.NE.0) then ! I have found my variable to read, I have enough to define it

      if (varNDIMS(v).NE.3) then
         WRITE(*,*) "Current version only support 3D var file"
         STOP
      else
         ! Get the actual values of the thing

         dim1 = dimLEN(varDIMIDS(v,1))
         dim2 = dimLEN(varDIMIDS(v,2))
         dim3 = dimLEN(varDIMIDS(v,3))

         ALLOCATE(varDATA(dim1,dim2,dim3))

         ret_stat = NF90_GET_VAR(ncid,v,varDATA)

!~          WRITE(*,*) "DIMS", dim1, dim2, dim3 ! lon lat time is getting out ... correct? (somehow netCDF reads backward)


         ALLOCATE(attNAMES(varNATTS(v)))

         do a=1,varNATTS(v)
           ret_stat = NF90_INQ_ATTNAME(ncid, v, a, attNAMES(a))
           if ((INDEX(attNAMES(a), "missing_value").GT.0) .OR. (INDEX(attNAMES(a), "_FillValue").GT.0) ) then
               ret_stat = NF90_GET_ATT(ncid, v, attNAMES(a), var_undef)
           endif
         enddo

!~          ! Get the missing value of the thing if exists ...
!~          WRITE(*,*) "Value for undef :: ", var_undef
      endif

      ! variable is read in ... now need to check where I have an actual value (not masked)


      varunmasked = 0

      DO j=1, dim2
      DO i=1, dim1
         if (varDATA(i,j,1).NE.var_undef) then
             ! count it in !
            varunmasked = varunmasked + 1
            forcing_surface_var(varunmasked,:) = varDATA(i,j,:)
         endif
      ENDDO
      ENDDO

      else
        WRITE(*,*) "var id is zero", v
      endif


      DEALLOCATE(dimNAMES)
      DEALLOCATE(dimLEN)
      DEALLOCATE(varNAMES)
      DEALLOCATE(varXTYPE)
      DEALLOCATE(varNDIMS)
      DEALLOCATE(varDIMIDS)
      DEALLOCATE(varNATTS)
      DEALLOCATE(varisDIM)
      DEALLOCATE(varDATA)

     END SUBROUTINE get_clim_forcing


     SUBROUTINE DO_spatialvars_step(stepstoDO,forcage_temperature_surface, forcage_epaisseurneige)

        use parameter_mod,  only: gridNoMax
        use vertclvars_mod, only: DO_vertclvars_step
        use grids_more,     only: WRITE_netCDF_output, indx_var_temp_ig, indx_var_palt, indx_var_plt

#if (CARBON == 1 )
        use grids_more,     only: indx_var_carb
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        INTEGER, INTENT(in) :: stepstoDO
        REAL, DIMENSION(1:gridNoMax,1:stepstoDO), INTENT(in) :: forcage_temperature_surface
        REAL, DIMENSION(1:gridNoMax,1:stepstoDO), INTENT(in), OPTIONAL :: forcage_epaisseurneige


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       integer :: gridp
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       ! This is where the parallelization could find place ...
!$omp parallel
!$omp do
       do gridp = 1, gridNoMax

!~          WRITE(*,*) "INTEGRATING ... ", gridp, "/", gridNoMax



         CALL DO_vertclvars_step(stepstoDO,Kp(:,gridp),T_bottom_SV(gridp),Temp(:,gridp), forcage_temperature_surface(gridp,:) &
                               , n(:,gridp),freeze_depth_SV(:,gridp), ALT_SV(gridp), altmax_ly_SV(gridp)                      &
                               , compteur_tstep_SV(gridp)                                                                     &
#if ( CARBON > 0 )
            ! CARBON ONLY VARIABLES
                               , deepSOM_a = deepSOM_a(:,gridp),deepSOM_s = deepSOM_s(:,gridp), deepSOM_p = deepSOM_p(:,gridp)&
                               , deepSOM = deepSOM(:,gridp), fc = fc_SV(:,:,gridp), b3_lok=b3_SV(gridp), b4_lok=b4_SV(gridp)  &
#endif
#if ( SNOW_EFFECT == 1 )
            ! SNOW ONLY VARIABLES
                               , snowlayer_thick_forcing = forcage_epaisseurneige,  Temp_snow_col=Temp_snow(:,gridp)          &
                               , snowlayer_depth = depth_snow_layer(gridp), snowlayer_nb = nb_snow_layer(gridp)               &
#endif
                               )

       enddo
!$omp end do
!$omp end parallel

       ! WRITE OUTPUT

       CALL WRITE_netCDF_output(Temp, indx_var_temp_ig)
       CALL WRITE_netCDF_output(ALT_SV, indx_var_palt)
       CALL WRITE_netCDF_output(freeze_depth_SV(1,:)-freeze_depth_SV(2,:), indx_var_plt)
#if ( CARBON == 1 )
       CALL WRITE_netCDF_output(deepSOM, indx_var_carb)
#endif
     END SUBROUTINE DO_spatialvars_step



     SUBROUTINE UPDATE_climate_forcing(stepstoDO,temperature_forcing_next, snowthickness_forcing_next)

       use parameter_mod,  only: gridNoMax, timFNoMax
       use timer_mod, only: cell_time


       INTEGER, INTENT(in) :: stepstoDO
       REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: temperature_forcing_next ! will be (1:gridNoMax,1:stepstoDO)
       REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(out), OPTIONAL :: snowthickness_forcing_next

! Need to assume that there, the grid is worked upon in its entirety (no parallelization) ...


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       INTEGER :: start_step, end_step, interim_nb, interim_end
       TYPE(cell_time) :: current_step

       current_step = compteur_tstep_SV(1)

       start_step = mod(current_step%current_step + 1,timFNoMax)

       if (current_step%current_step.EQ.0) then ! first ever time step
           start_step = 1 ! for the forcing index
       endif

       if (.NOT.ALLOCATED(temperature_forcing_next)) then
           allocate(temperature_forcing_next(1:gridNoMax,1:stepstoDO))
       endif

       if (PRESENT(snowthickness_forcing_next)) then
           if (.NOT.ALLOCATED(snowthickness_forcing_next)) then
               allocate(snowthickness_forcing_next(1:gridNoMax,1:stepstoDO))
           endif
       endif

       end_step = start_step + stepstoDO - 1

      WRITE(*,*) "Updating forcing ", current_step, start_step, end_step, timFNoMax

       if (end_step.GT.timFNoMax) then ! forcing serie requested is too long ... split !!
        WRITE(*,*) "endstep>timFNoMax "

        interim_nb = timFNoMax - start_step + 1
        temperature_forcing_next(:,1:interim_nb) = forcing_surface_temp(:,start_step:timFNoMax)
        ! then loop on the forcing
        interim_end = stepstoDO - interim_nb + 1
        temperature_forcing_next(:,interim_nb:stepstoDO) = forcing_surface_temp(:,1:interim_end)

#if ( SNOW_EFFECT == 1 )
        if (PRESENT(snowthickness_forcing_next)) then
          snowthickness_forcing_next(:,1:interim_nb) = forcing_snow_thick(:,start_step:timFNoMax)
          snowthickness_forcing_next(:,interim_nb:stepstoDO) = forcing_snow_thick(:,1:interim_end)
        endif
#endif

       else ! enough data already
         temperature_forcing_next(:,:) = forcing_surface_temp(:,start_step:end_step)
#if ( SNOW_EFFECT == 1 )
         if (PRESENT(snowthickness_forcing_next)) then
            snowthickness_forcing_next(:,:) = forcing_snow_thick(:,start_step:end_step)
         endif
#endif
       endif

     END SUBROUTINE UPDATE_climate_forcing

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


     SUBROUTINE SET_coupled_climate_forcing(stepstoDO,temperature_forcing_next, coupled_temp_set, b3_content, b4_content &
                                          , snowthick_forc_nxt, coupled_dsnow_set )

       use parameter_mod,  only: gridNoMax

       INTEGER, INTENT(in) :: stepstoDO
       REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(out)           :: temperature_forcing_next   ! will be (1:gridNoMax,1:stepstoDO)
       REAL, DIMENSION(:,:),              INTENT(in)            :: coupled_temp_set           ! must be (1:gridNoMax,1:stepstoDO)
       REAL, DIMENSION(:)               , INTENT(in),  OPTIONAL :: b3_content, b4_content     ! will be (1:gridNoMax)
       REAL, DIMENSION(:,:),              INTENT(in),  OPTIONAL :: coupled_dsnow_set          ! must be (1:gridNoMax,1:stepstoDO)
       REAL, DIMENSION(:,:),ALLOCATABLE , INTENT(out), OPTIONAL :: snowthick_forc_nxt ! must be (1:gridNoMax,1:stepstoDO)

! Need to assume that there, the grid is worked upon in its entirety (no parallelization) ...


       if (.NOT.ALLOCATED(temperature_forcing_next)) then
           allocate(temperature_forcing_next(1:gridNoMax,1:stepstoDO))
       endif

       temperature_forcing_next(:,:) = coupled_temp_set(:,:)
#if ( SNOW_EFFECT == 1 )
       if (.NOT.ALLOCATED(snowthick_forc_nxt)) then
         allocate(snowthick_forc_nxt(1:gridNoMax,1:stepstoDO))
       endif
       if (PRESENT(coupled_dsnow_set)) then
         snowthick_forc_nxt(:,:) = coupled_dsnow_set(:,:)
       else
         WRITE(*,*) "[ABORT] Missing forcing for dsnow in coupled mode"
       endif
#endif
       call fix_Kelvin_or_Celsius(temperature_forcing_next)

#if ( CARBON == 1 )
       if ((PRESENT(b3_content)).AND.(PRESENT(b4_content))) then
         b3_SV = b3_content
         b4_SV = b4_content
       else
         WRITE(*,*) "[ABORT] Missing forcing for b3, b4 in coupled mode"
       endif
#endif


     END SUBROUTINE SET_coupled_climate_forcing


    END MODULE spatialvars_mod
