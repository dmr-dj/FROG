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


      MODULE grids_more



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use, intrinsic :: iso_fortran_env, only: stdin=>input_unit, stdout=>output_unit, stderr=>error_unit

      IMPLICIT NONE

      PRIVATE

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   History
! dmr           Change from 0.0.0: Created a first version that can test read a topographic/mask grid assumed regular so far
! dmr           Change from 0.1.0: Modified so that the temperature forcing can be used as a grid definition file
! dmr           Change from 0.2.0: Modified so as to get all the information from a time forcing file directly
! dmr           Change from 0.3.0: added the netCDF initialization
! dmr           Change from 0.4.0: added the netCDF writing, expanded variables
! dmr           Change from 0.5.0: streamlined the netCDF writing for the variables part
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      CHARACTER(LEN=5), PARAMETER, PUBLIC :: version_mod ="0.6.0"

      integer, parameter  :: str_len =256


      INTEGER                     , DIMENSION(:,:), ALLOCATABLE :: two_to_oneDims
      INTEGER                     , DIMENSION(:)  , ALLOCATABLE :: one_to_twoDims_i, one_to_twoDims_j


                ! THESE ARE THE TWO MAIN VARIABLES READ IN THE NetCDF files
      INTEGER                                     , PUBLIC      :: nb_unmaskedp
      INTEGER                                     , PUBLIC      :: forcing_timelength


! --- temporary file names that will need to be filled in

      CHARACTER(len=str_len) :: mask_file ! = "tas_ewembi_1979-2016-r128x64-maskocean.nc4"! "tas_ewembi_1979-2016-r128x64-maskocean.nc4" ! "mask_ocean_r128x64.nc"


! --- spatial forcing needed [NOTA] to be implemented if needed in the future ...
!      REAL, DIMENSION(:), ALLOCATABLE :: BC_Kp, BC_Cp, BC_n, BC_pori, BC_porf


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! --- for output generation / netCDF
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      INTEGER         :: spat_dim2, spat_dim1
      REAL, PUBLIC    :: undefined_value
      CHARACTER(LEN=str_len) :: typology_file  ! ="file_typology-r128x64.nc"
      CHARACTER(LEN=str_len) :: netCDFout_file ! ="VAMPER-output.nc"



      INTEGER, PARAMETER :: nb_out_vars = 4, nb_dim_vars = 3

      CHARACTER(LEN=str_len), DIMENSION(nb_dim_vars), PARAMETER:: output_dim_names=[CHARACTER(len=str_len) :: "lat", "lon", "lev"]
      CHARACTER(LEN=str_len), DIMENSION(nb_out_vars), PARAMETER::                                     &
                              output_var_names=[CHARACTER(len=str_len) :: "temp_ig", "palt", "plt", "carb"],  &
                              output_unt_names=[CHARACTER(len=str_len) :: "K", "m", "m","g"],             &
                              output_std_names=[CHARACTER(len=str_len) :: "temperature_in_ground",    &
                                 "permafrost_active_layer_thickness", "permafrost_layer_thickness",""],  &
                              output_lng_names=[CHARACTER(len=str_len) ::                             &
                                 "solid_earth_subsurface_temperature", "", "",""],                       &
                              output_dms_names=[CHARACTER(len=str_len) :: "lev lon lat time",          &
                                 "lon lat time", "lon lat time", "lev lon lat time"]




      INTEGER, DIMENSION(0:nb_dim_vars) :: output_dim_len, output_dim_dimid
      INTEGER :: current_time_record

      INTEGER, DIMENSION(nb_out_vars) :: output_var_dimid
      INTEGER, PARAMETER              :: indx_var_temp_ig = 1, indx_var_palt=2, indx_var_plt=3

      INTEGER, PARAMETER              :: indx_var_carb = 4

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      INTERFACE WRITE_netCDF_output
        MODULE PROCEDURE WRITE_netCDF_output3D, WRITE_netCDF_output2D
      END INTERFACE WRITE_netCDF_output

      PUBLIC :: INIT_maskGRID, INIT_netCDF_output, indx_var_temp_ig, indx_var_palt, indx_var_plt, WRITE_netCDF_output

      PUBLIC :: indx_var_carb


      CONTAINS

! ---

     function flatten_it(spatial_array) result(flattened_array)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real, dimension(:,:), intent(in) :: spatial_array

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real, dimension(nb_unmaskedp) :: flattened_array

       integer :: i,j

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       do j=1, UBOUND(two_to_oneDims,DIM=2)
         do i=1, UBOUND(two_to_oneDims,DIM=1)
           flattened_array(two_to_oneDims(i,j)) = spatial_array(i,j)
         enddo
       enddo

     end function flatten_it

     function un_flatten_it(flattened_array) result(spatial_array)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real, dimension(:), intent(in) :: flattened_array

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real, dimension(spat_dim1,spat_dim2) :: spatial_array

       integer :: k

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       spatial_array(:,:) = undefined_value

       do k=1,UBOUND(flattened_array,dim=1)
         spatial_array(one_to_twoDims_i(k),one_to_twoDims_j(k)) = flattened_array(k)
       enddo

      end function un_flatten_it


      SUBROUTINE INIT_maskGRID()

        use netcdf

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer :: ncid, ret_stat, nDims, nvars, nGlobalAtts, unlimdimid
        integer :: d,v,a, i, j

        CHARACTER(len=NF90_MAX_NAME), DIMENSION(:),   ALLOCATABLE :: dimNAMES
        INTEGER                     , DIMENSION(:),   ALLOCATABLE :: dimLEN

        CHARACTER(len=NF90_MAX_NAME), DIMENSION(:),   ALLOCATABLE :: varNAMES
        INTEGER                     , DIMENSION(:),   ALLOCATABLE :: varXTYPE, varNDIMS, varNATTS
        INTEGER                     , DIMENSION(:,:), ALLOCATABLE :: varDIMIDS

        CHARACTER(len=NF90_MAX_NAME), DIMENSION(:),   ALLOCATABLE :: attNAMES

        LOGICAL                     , DIMENSION(:),   ALLOCATABLE :: varisDIM

        INTEGER                                                   :: maskVarID=0, dim1=-99, dim2=-99, dim3=-99
        INTEGER                     , DIMENSION(:),   ALLOCATABLE :: start_rd, count_rd

        REAL                        , DIMENSION(:,:), ALLOCATABLE :: maskdata
        REAL                                                      :: mask_undef = 0.0

! dmr   For namelist reading ...

        CHARACTER(len=str_len), PARAMETER                         :: file_path ="vamper_inputsGrid.nml"
        INTEGER                                                   :: rc,fu
        NAMELIST /inputsGrid/ mask_file, typology_file, netCDFout_file

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       THINGS TO DO ONCE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Read the appropriate namelist
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      INQUIRE (file=file_path, iostat=rc)
      IF (rc /= 0) THEN
         WRITE (stderr, '("Error: input file ", a, " does not exist")') file_path
         STOP
      ENDIF

      ! Open and read Namelist file.
      OPEN (action='read', file=file_path, iostat=rc, newunit=fu)
      IF (rc /= 0) WRITE (stderr, '("Error: Cannot open namelist file")')

      READ (nml=inputsGrid, iostat=rc, unit=fu)
      IF (rc /= 0) WRITE (stderr, '("Error: invalid Namelist format")')

      CLOSE (fu)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      ret_stat = nf90_open(mask_file, nf90_nowrite, ncid)

      ret_stat = nf90_inquire(ncid, nDims, nVars, nGlobalAtts, unlimdimid)

      ! I look forward to get a grid with one spatial variable and two dimensions for now, plus a time dimension
      ! Hence I should get nDims = 3, nVars = 3,  time, hence unlimdimid != -1


      WRITE(*,*) "info file", nDims, nVars, unlimdimid

      if ( (nDims.EQ.3).AND.(nVars.EQ.4).AND.(unlimdimid.NE.-1) ) then ! valid file: lat,lon,time + 1 spatial dimension
!~       if ( (nDims.EQ.2).AND.(nVars.EQ.3).AND.(unlimdimid.EQ.-1) ) then ! valid file: lat,lon,time + 1 spatial dimension

       ALLOCATE(dimNAMES(nDims))
       ALLOCATE(dimLEN(nDims))

       DO d=1,nDims
         ret_stat = nf90_inquire_dimension(ncid, d, dimNAMES(d), dimLEN(d))
       ENDDO

       ALLOCATE(varNAMES(nVars))
       ALLOCATE(varXTYPE(nVars))
       ALLOCATE(varNDIMS(nVars))
       ALLOCATE(varDIMIDS(nVars,nDims))
       ALLOCATE(varNATTS(nVars))
       ALLOCATE(varisDIM(nVars))

       varisDIM(:) = .FALSE.

       DO v=1,nVars
         ret_stat = nf90_inquire_variable(ncid, v, varNAMES(v), varXTYPE(v), varNDIMS(v), varDIMIDS(v,:), varNATTS(v))

         DO d=1,nDims
           if (TRIM(varNAMES(v)).EQ.TRIM(dimNAMES(d))) then
             varisDIM(v) = .TRUE.
           endif
         ENDDO

         if (.NOT. varisDIM(v)) then
            maskVarID = v
         endif

       ENDDO

      else
         WRITE(*,*) "netCDF file for mask does not match expectations"
         WRITE(*,*) "nVars ==", nVars
         WRITE(*,*) "nDims ==", nDims
         STOP
      endif

      if (maskVarID.NE.0) then ! I have found my variable to read, I have enough to define it


       v = maskVarID
       dim1 = dimLEN(varDIMIDS(v,1))
       dim2 = dimLEN(varDIMIDS(v,2))


       if (varNDIMS(maskVarID).NE.2) then
         if (varNDIMS(maskVarID).GT.3) then
            WRITE(*,*) "Current version only support 2D/3D masks"
            STOP
         else
!~             WRITE(*,*) "Provided with a 3D variable :", varNAMES(maskVarID)
!~             WRITE(*,*) "Time length is: ", dimLEN(unlimdimid), unlimdimid
            dim3 = dimLEN(varDIMIDS(v,3))

         endif
       endif


       ALLOCATE(start_rd(nDims))
       ALLOCATE(count_rd(nDims))

       if (dim3.gt.0) then
         start_rd(3) = 1
         count_rd(3) = 1
       endif

       start_rd(1:2) = 1

       count_rd(1) = dim1
       count_rd(2) = dim2


       ALLOCATE(maskdata(dim1,dim2))

       ret_stat = NF90_GET_VAR(ncid,v,maskdata,start_rd,count_rd)

         ALLOCATE(attNAMES(varNATTS(v)))

         do a=1,varNATTS(v)
           ret_stat = NF90_INQ_ATTNAME(ncid, v, a, attNAMES(a))
           if ((INDEX(attNAMES(a), "missing_value").GT.0) .OR. (INDEX(attNAMES(a), "_FillValue").GT.0) ) then
               ret_stat = NF90_GET_ATT(ncid, v, attNAMES(a), mask_undef)
           endif
         enddo

         ! Get the missing value of the thing if exists ...

      ! variable is read in ... now need to check where I have an actual value (not masked)

      nb_unmaskedp = 0

      DO j=1, dim2
      DO i=1, dim1
         if (maskdata(i,j).NE.mask_undef) then
             ! count it in !
            nb_unmaskedp = nb_unmaskedp + 1
         endif
      ENDDO
      ENDDO

!~       WRITE(*,*) "Count of unmasked datapoints", nb_unmaskedp

      ALLOCATE(two_to_oneDims(dim1,dim2))
      ALLOCATE(one_to_twoDims_i(nb_unmaskedp))
      ALLOCATE(one_to_twoDims_j(nb_unmaskedp))

      nb_unmaskedp = 0

      two_to_oneDims(:,:) = -1
      one_to_twoDims_i(:) = -1
      one_to_twoDims_j(:) = -1

      DO j=1, dim2
      DO i=1, dim1
         if (maskdata(i,j).NE.mask_undef) then
             ! count it in !
             nb_unmaskedp = nb_unmaskedp + 1
             two_to_oneDims(i,j) = nb_unmaskedp
             one_to_twoDims_i(nb_unmaskedp) = i
             one_to_twoDims_j(nb_unmaskedp) = j
         endif
      ENDDO
      ENDDO

      undefined_value = mask_undef
      spat_dim1 = dim1
      spat_dim2 = dim2

      ! So far I have created a list of coordinates and a list of points, need still to have the lat/lon values for all these points somewhere ? Or not ?

      endif

      forcing_timelength = dimLEN(unlimdimid)

      WRITE(*,*) "Currently, VAMPER seems to be off for a run with", nb_unmaskedp, " datapoints."
      WRITE(*,*) "Forcing file provides data for ", dimLEN(unlimdimid), "time steps"

      END SUBROUTINE INIT_maskGRID

! ---

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      SUBROUTINE INIT_netCDF_output

       USE netcdf
       USE parameter_mod, only: z_num, depth_levels=>D

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       CHARACTER(len=str_len) :: command_to_copy ! ="cp -f "//TRIM(typology_file)//" "//TRIM(netCDFout_file)

       INTEGER :: ncid, nDims, nvars, nGlobalAtts, unlimdimid, d, n, depth_varid

       CHARACTER(len=NF90_MAX_NAME), DIMENSION(:),   ALLOCATABLE :: dimNAMES
       INTEGER                     , DIMENSION(:),   ALLOCATABLE :: dimLEN

       LOGICAL, DIMENSION(nb_dim_vars) :: dim_exists_file

       INTEGER :: c, indx_var

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       ! Copy the typology file into the new output file
       command_to_copy="cp -f "//TRIM(typology_file)//" "//TRIM(netCDFout_file)
       WRITE(*,*) "COMMAND // ", TRIM(command_to_copy)
       call execute_command_line(TRIM(command_to_copy))

       ! Now need to check and define the dimension if needed (levels for sure)

       dim_exists_file(:) = .FALSE.

       call handle_err(                                                                                &
            nf90_open(path = TRIM(netCDFout_file), mode = NF90_WRITE, ncid = ncid)                     & ! open existing netCDF dataset
                      , __LINE__)

       ! Check that lat and lon already exist in the file ...

       call handle_err(                                                                                &
           nf90_inquire(ncid, nDims, nVars, nGlobalAtts, unlimdimid)                                   & ! inquire existing dimensions in the file
                      , __LINE__)

       output_dim_dimid(0) = unlimdimid

       ALLOCATE(dimNAMES(nDims))
       ALLOCATE(dimLEN(nDims))

       DO d=1,nDims

       call handle_err(                                                                                &
            nf90_inquire_dimension(ncid, d, dimNAMES(d), dimLEN(d))                                    &
                      , __LINE__)
       DO n=1, nb_dim_vars
          if (TRIM(dimNAMES(d)) == TRIM(output_dim_names(n))) then
             dim_exists_file(n) = .TRUE.
             call handle_err(                                                                          &
                  nf90_inq_dimid(ncid, TRIM(dimNAMES(d)), output_dim_dimid(n))                         &
                      , __LINE__)
             output_dim_len(n) = dimLEN(d)
          endif
       ENDDO

       ENDDO

       if (.NOT. ALL(dim_exists_file)) then                                                       ! at least one dimension does not exist

       DO n=1, nb_dim_vars

         if (.NOT.dim_exists_file(n)) then

         if (n.eq.3) then ! lev
           output_dim_len(n) = z_num
         endif

         call handle_err(                                                                              &
              nf90_redef(ncid=ncid)                                                                    & ! put it into define mode
                      , __LINE__)

         call handle_err(                                                                              &
              nf90_def_dim(ncid, output_dim_names(n), output_dim_len(n), output_dim_dimid(n))          & ! define additional dimensions
                      , __LINE__)

         endif

       ENDDO

       endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       ! define the depth variable (in meters below the surface of the ground, positive downwards)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       call handle_err(                                                                                &
            nf90_def_var(ncid, "depth", NF90_FLOAT, output_dim_dimid(3), depth_varid)                  &  ! define additional variables
                      , __LINE__)

       call handle_err(                                                                                &
            nf90_put_att(ncid, depth_varid, "units", "meters")                                         &
                      , __LINE__)

       call handle_err(                                                                                &
            nf90_put_att(ncid, depth_varid, "standard_name", "depth")                                  &
                      , __LINE__)

       call handle_err(                                                                                &
            nf90_put_att(ncid, depth_varid, "long_name", "depth_below_ground_surface")                 &
                      , __LINE__)


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       ! define the variables iteratively
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        DO indx_var=1, nb_out_vars

           c = count_words_string(output_dms_names(indx_var))

           SELECT CASE (c)
             CASE (4)
               call handle_err(                                                                                           &
                  nf90_def_var(ncid, output_var_names(indx_var), NF90_FLOAT, output_dim_dimid(3:0:-1)                     &
                            , output_var_dimid(indx_var))                                                                 &
                            , __LINE__)
             CASE (3)
               call handle_err(                                                                                           &
                  nf90_def_var(ncid, output_var_names(indx_var), NF90_FLOAT, output_dim_dimid(2:0:-1)                     &
                            , output_var_dimid(indx_var))                                                                 &
                            , __LINE__)
             CASE DEFAULT
               WRITE(*,*) "NetCDF init is not setup for a variable with ", c, " axes"

           END SELECT

           call handle_err(                                                                                               &
              nf90_put_att(ncid, output_var_dimid(indx_var), "units", output_unt_names(indx_var))                         &
                        , __LINE__)

           call handle_err(                                                                                               &
              nf90_put_att(ncid, output_var_dimid(indx_var), "standard_name", output_std_names(indx_var))                 &
                        , __LINE__)


           if (output_lng_names(indx_var) /= "") then
             call handle_err(                                                                                             &
                nf90_put_att(ncid, output_var_dimid(indx_var), "long_name", output_lng_names(indx_var))                   &
                          , __LINE__)
           endif

           call handle_err(                                                                                               &
              nf90_put_att(ncid, output_var_dimid(indx_var), "_FillValue", undefined_value)                               &
                        , __LINE__)

        ENDDO

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       ! Finalize definitions
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       call handle_err(                                                                                &
            nf90_enddef(ncid)                                                                          & ! check definitions, leave define mode
                      , __LINE__)

       ! Write the depth(z_num) variable in the output file
       call handle_err(                                                                                &
            nf90_put_var(ncid, depth_varid, depth_levels)                                              &  ! provide new variable values
                      , __LINE__)


       call handle_err(                                                                                & ! close netcdf dataset
            nf90_close(ncid)                                                                           &
                      , __LINE__)


       current_time_record = 0

      END SUBROUTINE INIT_netCDF_output



      SUBROUTINE WRITE_netCDF_output3D(var_to_write, indx_var)

       USE netcdf
       USE parameter_mod, only: z_num

       INTEGER :: ncid, z
       REAL, DIMENSION(:,:), INTENT(in) :: var_to_write ! variables are 1:z_num, 1:gridNoMax
       INTEGER, INTENT(in) :: indx_var
       REAL, DIMENSION(z_num,spat_dim1,spat_dim2) :: nc_vartowrite ! lev, lon, lat, time

       current_time_record = current_time_record + 1

       do z=1,z_num
         nc_vartowrite(z,:,:) = un_flatten_it(var_to_write(z,:))
       enddo

       !WRITE(*,*) "VARIABLE TO NC WRITE: ", MINVAL(nc_vartowrite), MAXVAL(nc_vartowrite)
       !WRITE(*,*) "dims :: ", z_num, spat_dim1, spat_dim2

       call handle_err(                                                                                &
            nf90_open(path = TRIM(netCDFout_file), mode = NF90_WRITE, ncid = ncid)                     & ! open existing netCDF dataset
                      , __LINE__)

       call handle_err(                                                                                &
            nf90_put_var(ncid, output_var_dimid(indx_var) ,                                            &
                         nc_vartowrite(1:z_num,1:spat_dim1,1:spat_dim2),                               &
                         start=(/1,1,1,current_time_record/), count=(/z_num, spat_dim1, spat_dim2,1/)) &  ! provide new variable values
                      , __LINE__)


       call handle_err(                                                                                & ! close netcdf dataset
            nf90_close(ncid)                                                                           &
                      , __LINE__)


      END SUBROUTINE WRITE_netCDF_output3D

      SUBROUTINE WRITE_netCDF_output2D(var_to_write, indx_var)

       USE netcdf
!~        USE parameter_mod, only: z_num

       INTEGER :: ncid
       REAL, DIMENSION(:), INTENT(in) :: var_to_write ! variables are 1:gridNoMax
       INTEGER, INTENT(in) :: indx_var
       REAL, DIMENSION(spat_dim1,spat_dim2) :: nc_vartowrite ! lon, lat, time

!~        current_time_record = current_time_record + 1

!~        do z=1,z_num
         nc_vartowrite(:,:) = un_flatten_it(var_to_write(:))
!~        enddo

       !WRITE(*,*) "VARIABLE TO NC WRITE: ", MINVAL(nc_vartowrite), MAXVAL(nc_vartowrite)
       !WRITE(*,*) "dims :: ", spat_dim1, spat_dim2

       call handle_err(                                                                                &
            nf90_open(path = TRIM(netCDFout_file), mode = NF90_WRITE, ncid = ncid)                     & ! open existing netCDF dataset
                      , __LINE__)

       call handle_err(                                                                                &
            nf90_put_var(ncid, output_var_dimid(indx_var) ,                                            &
                         nc_vartowrite(1:spat_dim1,1:spat_dim2),                                       &
                         start=(/1,1,current_time_record/), count=(/spat_dim1, spat_dim2,1/))          & ! provide new variable values
                      , __LINE__)

       call handle_err(                                                                                & ! close netcdf dataset
            nf90_close(ncid)                                                                           &
                      , __LINE__)


      END SUBROUTINE WRITE_netCDF_output2D


      SUBROUTINE handle_err(nf90_code, line_location)

        USE netcdf, ONLY: nf90_noerr, nf90_strerror

        INTEGER, INTENT(IN) :: nf90_code
        INTEGER, INTENT(in) :: line_location

        if (nf90_code /= nf90_noerr) then
           WRITE(*,*) "Error in netCDF operation", line_location, nf90_strerror(nf90_code)
           STOP
        endif

      END SUBROUTINE handle_err

! dmr adapted from: https://www.reddit.com/r/fortran/comments/yzhkw5/reading_a_file_into_individual_words/

      FUNCTION count_words_string(input_string) result(word_count)

       CHARACTER(LEN=*), INTENT(in) :: input_string
       INTEGER :: word_count

       CHARACTER(LEN=32), DIMENSION(10) :: words

       INTEGER :: i, ios, counter

        counter = 1
        i = 0
        do

          read(input_string, *, iostat=ios) words(counter:counter+i)

          if (ios /= 0) then
            words(counter+i) = ''
            counter = counter + i
            exit
          else
            i = i + 1
          endif

        enddo

        words(counter) = ''
        word_count = counter - 1


      END FUNCTION count_words_string

      END MODULE grids_more


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
