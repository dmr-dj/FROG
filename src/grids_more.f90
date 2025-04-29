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

      IMPLICIT NONE

      PRIVATE

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   History
! dmr           Change from 0.0.0: Created a first version that can test read a topographic/mask grid assumed regular so far
! dmr           Change from 0.1.0: Modified so that the temperature forcing can be used as a grid definition file
! dmr           Change from 0.2.0: Modified so as to get all the information from a time forcing file directly
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      CHARACTER(LEN=5), PARAMETER, PUBLIC :: version_mod ="0.3.0"

      integer, parameter  :: str_len =256


      INTEGER                     , DIMENSION(:,:), ALLOCATABLE :: two_to_oneDims
      INTEGER                     , DIMENSION(:)  , ALLOCATABLE :: one_to_twoDims_i, one_to_twoDims_j


                ! THESE ARE THE TWO MAIN VARIABLES READ IN THE NetCDF files
      INTEGER                                     , PUBLIC      :: nb_unmaskedp
      INTEGER                                     , PUBLIC      :: forcing_timelength


! --- temporary file names that will need to be filled in

      CHARACTER(len=str_len) :: mask_file = "tas_ewembi_1979-2016-r128x64-maskocean.nc4"! "tas_ewembi_1979-2016-r128x64-maskocean.nc4" ! "mask_ocean_r128x64.nc"

! --- spatial forcing needed
      REAL, DIMENSION(:), ALLOCATABLE :: BC_Kp, BC_Cp, BC_n, BC_pori, BC_porf

      PUBLIC :: INIT_maskGRID

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

       integer :: i,j,k

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       do j=1, UBOUND(two_to_oneDims,DIM=2)
         do i=1, UBOUND(two_to_oneDims,DIM=1)
           flattened_array(two_to_oneDims(i,j)) = spatial_array(i,j)
         enddo
       enddo

     end function flatten_it

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

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       THINGS TO DO ONCE
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

      ! So far I have created a list of coordinates and a list of points, need still to have the lat/lon values for all these points somewhere ? Or not ?

      endif

      forcing_timelength = dimLEN(unlimdimid)

      WRITE(*,*) "Currently, VAMPER seems to be off for a run with", nb_unmaskedp, " datapoints."
      WRITE(*,*) "Forcing file provides data for ", dimLEN(unlimdimid), "time steps"

      END SUBROUTINE INIT_maskGRID

! ---

      END MODULE grids_more


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
