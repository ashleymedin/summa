! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module read_icond_module
USE nr_type
USE netcdf
USE globalData,only: ixHRUfile_min,ixHRUfile_max
USE globalData,only: nTimeDelay      ! number of hours in the time delay histogram
USE globalData,only: nSpecBand       ! number of spectral bands
USE globalData,only: int8Missing     ! missing integer
USE globalData,only:verySmaller      ! a smaller number used as an additive constant to check if substantial difference among real numbers

! access domain types
USE globalData,only:upland             ! domain type for upland areas
USE globalData,only:glacCln1           ! first domain type for glacier clean areas
USE globalData,only:glacCln2           ! second domain type for glacier clean areas
USE globalData,only:glacDbr            ! domain type for glacier debris areas
USE globalData,only:wetland            ! domain type for wetland areas
USE globalData,only:glacieret          ! domain type for glaciers considered too small for flow

implicit none
private
public::read_icond
public::read_icondGlac
public::read_icond_nlayers

contains

 ! ************************************************************************************************
 ! public subroutine read_icond_nlayers: read model initial conditions file for number of snow/soil layers
 ! ************************************************************************************************
 subroutine read_icond_nlayers(iconFile,nGRU,nDOM,indx_meta,err,message)
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nr_type
 USE var_lookup,only:iLookINDEX                        ! variable lookup structure
 USE globalData,only:gru_struc                         ! gru-hru mapping structures
 USE netcdf_util_module,only:nc_file_close             ! close netcdf file
 USE netcdf_util_module,only:nc_file_open              ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                ! netcdf error handling
 USE data_types,only:var_info                          ! metadata
 implicit none

 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)  ,intent(in)   :: iconFile            ! name of input (restart) file
 integer(i4b)  ,intent(in)   :: nGRU                ! total # of GRUs in run space
 integer(i4b)  ,intent(inout):: nDOM                ! max number of domains in any HRU
 type(var_info),intent(in)   :: indx_meta(:)        ! metadata
 integer(i4b)  ,intent(out)  :: err                 ! error code
 character(*)  ,intent(out)  :: message             ! returned error message
 ! locals
 integer(i4b)                :: i,j                 ! loop indices
 integer(i4b)                :: ncID                ! netcdf file id
 integer(i4b)                :: dimID               ! netcdf file dimension id
 integer(i4b)                :: varID               ! netcdf variable id
 integer(i4b)                :: fileGRU             ! number of GRUs in netcdf file
 integer(i4b)                :: fileHRU             ! number of HRUs in netcdf file
 integer(i4b)                :: fileDOM             ! number of domains in netcdf file
 integer(i4b)                :: snowID, soilID      ! netcdf variable ids
 integer(i4b)                :: glceID, lakeID      ! netcdf variable ids
 integer(i4b)                :: iGRU, iHRU, iDOM    ! loop indexes
 integer(i4b)                :: iHRU_global         ! index of HRU in the netcdf file
 logical(lgt)                :: no_glceData         ! flag that no ice data in icond
 logical(lgt)                :: no_lakeData         ! flag that no lake data in icond
 logical(lgt)                :: no_dom              ! flag that no domain variable in file
 integer(i4b),allocatable    :: snowData1(:)        ! number of snow layers in all HRUs
 integer(i4b),allocatable    :: soilData1(:)        ! number of soil layers in all HRUs
 integer(i4b),allocatable    :: glceData1(:)        ! number of glacier ice layers in all HRUs
 integer(i4b),allocatable    :: lakeData1(:)        ! number of lake layers in all HRUs
 integer(i4b),allocatable    :: snowData2(:,:)      ! number of snow layers in all HRUs
 integer(i4b),allocatable    :: soilData2(:,:)      ! number of soil layers in all HRUs
 integer(i4b),allocatable    :: glceData2(:,:)      ! number of glacier ice layers in all HRUs
 integer(i4b),allocatable    :: lakeData2(:,:)      ! number of lake layers in all HRUs
 integer(i4b),allocatable    :: dom_type(:,:)       ! read domain type in from netcdf file
 character(len=256)          :: cmessage            ! downstream error message
 integer(i8b),allocatable    :: gru_id(:)           ! GRU id
 integer(i8b),allocatable    :: hru_id(:)           ! HRU id
 integer(i4b),allocatable    :: gruid_to_index(:)   ! mapping from gru_id to index in gru_struc
 integer(i4b),allocatable    :: hrunc_to_index(:,:) ! mapping from hru_nc to index in gru_struc
 

 ! --------------------------------------------------------------------------------------------------------
 ! initialize error message
 err=0
 message = 'read_icond_nlayers/'
 no_glceData = .false.
 no_lakeData = .false.
 no_dom = .false.

 ! open netcdf file
 call nc_file_open(iconFile,nf90_nowrite,ncID,err,cmessage);
 if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

 ! get number of HRUs in file
 err = nf90_inq_dimid(ncID,"hru",dimId);               if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncID,dimId,len=fileHRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru dimension/'//trim(nf90_strerror(err)); return; end if

! get number of domains with type in file, if present
 err = nf90_inq_dimid(ncID,"dom",dimId)               
 if(err/=nf90_noerr)then
  fileDOM = 1 ! backwards compatible, just upland domain
  no_dom = .true.
  allocate(dom_type(1,fileHRU))
  dom_type = upland
  err=nf90_noerr    ! reset this err
 else
  err = nf90_inquire_dimension(ncID,dimId,len=fileDOM); if(err/=nf90_noerr)then; message=trim(message)//'problem reading dom dimension/'//trim(nf90_strerror(err)); return; end if
  ! read dom_type from netcdf file
  allocate(dom_type(fileDOM,fileHRU))
  err = nf90_inq_varid(ncID,"domType",varID);  if (err/=0) then; message=trim(message)//'problem finding domType'; return; end if
  err = nf90_get_var(ncID,varID,dom_type);     if (err/=0) then; message=trim(message)//'problem reading domType'; return; end if
 end if
 nDOM = fileDOM

 ! read hru_id from netcdf file
 allocate(hru_id(fileHRU))
 err = nf90_inq_varid(ncID,"hruId",varID);     if (err/=0) then; message=trim(message)//'problem finding hruId'; return; end if
 err = nf90_get_var(ncID,varID,hru_id);        if (err/=0) then; message=trim(message)//'problem reading hruId'; return; end if

 ! check if the file has the GRU dimension
 err = nf90_inq_dimid(ncID,"gru",dimID);    
 if(err/=nf90_noerr)then         
   write(*,*) 'WARNING: GRU is not in the initial conditions file ... assuming HRUs in attribute order'
   fileGRU = size(gru_struc(:)%gru_id)
   err=nf90_noerr    ! reset this err
   allocate(gru_id(fileHRU))
   gru_id = gru_struc(:)%gru_id
 else
   err = nf90_inquire_dimension(ncID,dimID,len=fileGRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading gru dimension/'//trim(nf90_strerror(err)); return; end if
   ! read gru_id from netcdf file
   allocate(gru_id(fileGRU))
   err = nf90_inq_varid(ncID,"gruId",varID);   if (err/=0) then; message=trim(message)//'problem finding gruId'; return; end if
   err = nf90_get_var(ncID,varID,gru_id);      if (err/=0) then; message=trim(message)//'problem reading gruId'; return; end if
 end if

 ! Allocate the mapping arrays
 allocate(gruid_to_index(fileGRU), hrunc_to_index(fileGRU,maxval(gru_struc(:)%hruCount)))

 ! Populate the mapping arrays
 do i = 1, fileGRU
   gruid_to_index(i) = -1  ! Initialize with an invalid index
   do iGRU = 1, nGRU
     if (gru_struc(iGRU)%gru_id == gru_id(i)) then
       gruid_to_index(i) = iGRU
       do j = 1, gru_struc(iGRU)%hruCount
         hrunc_to_index(i,j) = -1 
         do iHRU = 1, fileHRU ! assumes HRUs are in GRU order
           if (gru_struc(iGRU)%hruInfo(j)%hru_id == hru_id(iHRU)) then
             hrunc_to_index(i,j) = iHRU
             exit
           endif
         end do
       end do ! HRU id loop
       exit
     endif
   end do
 end do

 ! allocate storage for reading from file (allocate entire file size, even when doing subdomain run)
 allocate(snowData1(fileHRU),snowData2(fileDOM,fileHRU))
 allocate(soilData1(fileHRU),soilData2(fileDOM,fileHRU))
 allocate(glceData1(fileHRU),glceData2(fileDOM,fileHRU))
 allocate(lakeData1(fileHRU),lakeData2(fileDOM,fileHRU))
 snowData1 = 0
 soilData1 = 0
 glceData1 = 0
 lakeData1 = 0
 snowData2 = 0
 soilData2 = 0
 glceData2 = 0
 lakeData2 = 0

 ! count domains and set domain type
 ! NOTE: dom_type 0 will be no domain
 do i = 1,fileGRU
   iGRU = gruid_to_index(i)
   do j = 1,gru_struc(iGRU)%hruCount
     iHRU = hrunc_to_index(i,j)- minval(hrunc_to_index(i,:)) + 1 ! assumes HRUs are in GRU order
     iHRU_global = hrunc_to_index(i,j)
     gru_struc(iGRU)%hruInfo(iHRU)%domCount = 1                                              ! upland domain always present, for changing size glaciers and lakes
     if (any(dom_type(1:fileDOM,iHRU_global)==glacCln1)) &
       gru_struc(iGRU)%hruInfo(iHRU)%domCount = gru_struc(iGRU)%hruInfo(iHRU)%domCount + 1   ! glacier clean domain 1 possible
     if (any(dom_type(1:fileDOM,iHRU_global)==glacCln2)) &
       gru_struc(iGRU)%hruInfo(iHRU)%domCount = gru_struc(iGRU)%hruInfo(iHRU)%domCount + 1   ! glacier clean domain 2 possible
      if (any(dom_type(1:fileDOM,iHRU_global)==glacDbr)) &
       gru_struc(iGRU)%hruInfo(iHRU)%domCount = gru_struc(iGRU)%hruInfo(iHRU)%domCount + 1   ! glacier debris domain possible
     if (any(dom_type(1:fileDOM,iHRU_global)==wetland)) &
       gru_struc(iGRU)%hruInfo(iHRU)%domCount = gru_struc(iGRU)%hruInfo(iHRU)%domCount + 1   ! wetland domain possible
     if (any(dom_type(1:fileDOM,iHRU_global)==glacieret)) &
       gru_struc(iGRU)%hruInfo(iHRU)%domCount = gru_struc(iGRU)%hruInfo(iHRU)%domCount + 1   ! glacieret domain possible
     allocate(gru_struc(iGRU)%hruInfo(iHRU)%domInfo(gru_struc(iGRU)%hruInfo(iHRU)%domCount)) ! allocate third level of gru to hru map
     gru_struc(iGRU)%hruInfo(iHRU)%domInfo(:)%dom_type = dom_type(1:gru_struc(iGRU)%hruInfo(iHRU)%domCount,iHRU_global)
   enddo
 enddo

 ! get netcdf ids for the variables holding number of layers in each domain or hru
 err = nf90_inq_varid(ncID,trim(indx_meta(iLookINDEX%nSnow)%varName),snowID); call netcdf_err(err,message)
 err = nf90_inq_varid(ncID,trim(indx_meta(iLookINDEX%nLake)%varName),lakeID)
 if(err/=nf90_noerr ) no_lakeData = .true.
 err = nf90_inq_varid(ncID,trim(indx_meta(iLookINDEX%nSoil)%varName),soilID); call netcdf_err(err,message)
 err = nf90_inq_varid(ncID,trim(indx_meta(iLookINDEX%nGlce)%varName), glceID)
 if(err/=nf90_noerr) no_glceData = .true.

 ! get nLayer data (reads entire state file)
 if(no_dom)then
   err = nf90_get_var(ncID,snowID,snowData1); call netcdf_err(err,message)
   err = nf90_get_var(ncID,soilID,soilData1); call netcdf_err(err,message)
   if (.not. no_glceData) err = nf90_get_var(ncID,glceID,glceData1); call netcdf_err(err,message)
   if (.not. no_lakeData) err = nf90_get_var(ncID,lakeID,lakeData1); call netcdf_err(err,message)
 else
   err = nf90_get_var(ncID,snowID,snowData2); call netcdf_err(err,message)
   err = nf90_get_var(ncID,soilID,soilData2); call netcdf_err(err,message)
   if (.not. no_glceData) err = nf90_get_var(ncID,glceID,glceData2); call netcdf_err(err,message)
   if (.not. no_lakeData) err = nf90_get_var(ncID,lakeID,lakeData2); call netcdf_err(err,message)
 endif
 ixHRUfile_min=huge(1)
 ixHRUfile_max=0
 ! find the min and max hru indices in the state file
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   iHRU_global = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
   if(iHRU_global < ixHRUfile_min) ixHRUfile_min = iHRU_global
   if(iHRU_global > ixHRUfile_max) ixHRUfile_max = iHRU_global
  end do
 end do

 ! loop over grus in current run to update snow/soil layer information
 do i = 1,fileGRU
  iGRU = gruid_to_index(i)
  do j = 1,gru_struc(iGRU)%hruCount
   iHRU = hrunc_to_index(i,j)- minval(hrunc_to_index(i,:)) + 1 ! assumes HRUs are in GRU order
   iHRU_global = hrunc_to_index(i,j)
   do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount
    if(no_dom)then
      gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSnow = snowData1(iHRU_global)
      gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nLake = lakeData1(iHRU_global)
      gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSoil = soilData1(iHRU_global)
      gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nGlce = glceData1(iHRU_global)
    else
     gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSnow = snowData2(iDOM,iHRU_global)
     gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nLake = lakeData2(iDOM,iHRU_global)
     gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSoil = soilData2(iDOM,iHRU_global)
     gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nGlce = glceData2(iDOM,iHRU_global)
    endif
   end do
  end do
 end do

 ! close file
 call nc_file_close(ncID,err,cmessage)
 if(err/=0)then;message=trim(message)//trim(cmessage);return;end if

 ! cleanup
 deallocate(snowData1,lakeData1,soilData1,glceData1,snowData2,lakeData2,soilData2,glceData2,dom_type)
 deallocate(gru_id,hru_id,gruid_to_index,hrunc_to_index)

 end subroutine read_icond_nlayers


 ! ************************************************************************************************
 ! public subroutine read_icond: read model initial conditions
 ! ************************************************************************************************
 subroutine read_icond(iconFile,                      & ! intent(in):    name of initial conditions file
                       nGRU,                          & ! intent(in):    number of GRUs
                       mparData,                      & ! intent(in):    model parameters
                       progData,                      & ! intent(inout): model prognostic variables
                       bvarData,                      & ! intent(inout): model basin (GRU) variables
                       indxData,                      & ! intent(inout): model indices
                       no_icond_enth,                 & ! intent(out):   flag that enthalpy variables are not in the file
                       err,message)                     ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nr_type
 USE var_lookup,only:iLookVarType                       ! variable lookup structure
 USE var_lookup,only:iLookPROG                          ! variable lookup structure
 USE var_lookup,only:iLookPARAM                         ! variable lookup structure
 USE var_lookup,only:iLookBVAR                          ! variable lookup structure
 USE var_lookup,only:iLookINDEX                         ! variable lookup structure
 USE globalData,only:prog_meta                          ! metadata for prognostic variables
 USE globalData,only:bvar_meta                          ! metadata for basin (GRU) variables
 USE globalData,only:gru_struc                          ! gru-hru mapping structures
 USE globalData,only:startGRU                           ! index of first gru for parallel runs
 USE globalData,only:iname_soil,iname_snow,iname_glce,iname_lake ! named variables to describe the type of layer
 USE netcdf_util_module,only:nc_file_open               ! open netcdf file
 USE netcdf_util_module,only:nc_file_close              ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                 ! netcdf error handling
 USE data_types,only:gru_hru_dom_doubleVec              ! full double precision structure
 USE data_types,only:gru_hru_dom_intVec                 ! full integer structure
 USE data_types,only:gru_doubleVec                      ! gru-length double precision structure (basin variables)
 USE get_ixName_module,only:get_varTypeName             ! to access type strings for error messages
 USE updateState_module,only:updateSoil                  ! update soil states

 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)           ,intent(in)        :: iconFile                 ! name of netcdf file containing the initial conditions
 integer(i4b)           ,intent(in)        :: nGRU                     ! number of grouped response units in simulation domain
 type(gru_hru_dom_doubleVec),intent(in)    :: mparData                 ! model parameters
 type(gru_hru_dom_doubleVec),intent(inout) :: progData                 ! model prognostic variables
 type(gru_doubleVec)    ,intent(inout)     :: bvarData                 ! model basin (GRU) variables
 type(gru_hru_dom_intVec),intent(inout)    :: indxData                 ! model indices
 logical                ,intent(out)       :: no_icond_enth            ! flag that enthalpy variables are not in the file
 integer(i4b)           ,intent(out)       :: err                      ! error code
 character(*)           ,intent(out)       :: message                  ! returned error message
 ! locals
 character(len=256)                        :: cmessage                 ! downstream error message
 integer(i4b)                              :: fileHRU                  ! number of HRUs in file
 integer(i4b)                              :: fileGRU                  ! number of GRUs in file
 integer(i4b)                              :: fileDOM                  ! number of domains in netcdf file
 integer(i4b)                              :: iVar,i,j                 ! loop indices
 integer(i4b),dimension(1)                 :: nrdx                     ! intermediate array of loop indices
 integer(i4b),dimension(5)                 :: ngdx                     ! intermediate array of loop indices
 integer(i4b)                              :: iGRU                     ! loop index
 integer(i4b)                              :: iHRU                     ! loop index
 integer(i4b)                              :: iDOM                     ! loop index
 integer(i4b)                              :: dimID                    ! varible dimension ids
 integer(i4b)                              :: ncVarID                  ! variable ID in netcdf file
 character(256)                            :: dimName                  ! not used except as a placeholder in call to inq_dim function
 integer(i4b)                              :: dimLen                   ! data dimensions
 integer(i4b)                              :: ncID                     ! netcdf file ID
 integer(i4b)                              :: iHRU_global              ! index of HRU in the netcdf file
 real(rkind),allocatable                   :: varData(:)              ! variable data storage
 real(rkind),allocatable                   :: varData2(:,:)            ! variable data storage
 real(rkind),allocatable                   :: varData3(:,:,:)          ! variable data storage
 integer(i4b)                              :: nSnow,nLake,nSoil,nGlce,nToto !# layers
 integer(i4b)                              :: noThetaChange            ! number of layers with no change in total water content (bottom layers)
 integer(i4b)                              :: nTDH                     ! number of points in time-delay 
 integer(i4b)                              :: nGlacier                 ! number of glaciers in basin (attribute files
 integer(i4b)                              :: fileglac                 ! max number of glaciers in any GRU
 integer(i4b)                              :: iLayer,jLayer            ! layer indices
 logical(lgt)                              :: has_glacier              ! flag for glacier presence in at least one GRU
 logical(lgt)                              :: has_wetland              ! flag for wetland/lake presence in at least one GRU
 logical(lgt)                              :: no_dom                   ! flag for no domain variable in file
 ! currently only writing restart for progressive variables with these dimensions
 character(len=32),parameter               :: scalDimName   ='scalarv' ! dimension name for scalar data
 character(len=32),parameter               :: midSoilDimName='midSoil' ! dimension name for soil-only layers
 character(len=32),parameter               :: midTotoDimName='midToto' ! dimension name for layered varaiables
 character(len=32),parameter               :: ifcTotoDimName='ifcToto' ! dimension name for layered varaiables
 character(len=32),parameter               :: tdhDimName    ='tdh'     ! dimension name for time-delay basin variables
 character(len=32),parameter               :: nglDimName    ='grid'    ! dimension name for grid variables (glaciers)
 integer(i8b),allocatable                  :: gru_id(:)                ! GRU id
 integer(i8b),allocatable                  :: hru_id(:)                ! HRU id
 integer(i8b),allocatable                  :: glac_id(:,:)             ! glac id
 integer(i4b),allocatable                  :: gruid_to_index(:)        ! mapping from gru_id to index in gru_struc
 integer(i4b),allocatable                  :: hrunc_to_index(:,:)      ! mapping from hru_nc to index in gru_struc

 ! --------------------------------------------------------------------------------------------------------
 ! Start procedure here
 err=0; message="read_icond/"

 ! --------------------------------------------------------------------------------------------------------
 ! (1) read the file
 ! --------------------------------------------------------------------------------------------------------
 ! open netcdf file
 call nc_file_open(iconFile,nf90_nowrite,ncID,err,cmessage)
 if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
 no_dom = .false.

 ! get number of HRUs in file
 err = nf90_inq_dimid(ncID,"hru",dimID);               if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncID,dimID,len=fileHRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru dimension/'//trim(nf90_strerror(err)); return; end if

 ! get max number of DOMs any HRU in file, if present
 err = nf90_inq_dimid(ncID,"dom",dimId)               
 if(err/=nf90_noerr)then
  fileDOM = 1
  no_dom = .true.
  err=nf90_noerr    ! reset this err
 else
  err = nf90_inquire_dimension(ncID,dimId,len=fileDOM); if(err/=nf90_noerr)then; message=trim(message)//'problem reading dom dimension/'//trim(nf90_strerror(err)); return; end if
 end if

 ! read hru_id from netcdf file
 allocate(hru_id(fileHRU))
 err = nf90_inq_varid(ncID,"hruId",ncVarID);     if (err/=0) then; message=trim(message)//'problem finding hruId'; return; end if
 err = nf90_get_var(ncID,ncVarID,hru_id);        if (err/=0) then; message=trim(message)//'problem reading hruId'; return; end if

 ! check if the file has the GRU dimension
 err = nf90_inq_dimid(ncID,"gru",dimID);    
 if(err/=nf90_noerr)then         
   fileGRU = size(gru_struc(:)%gru_id)
   err=nf90_noerr    ! reset this err
   allocate(gru_id(fileHRU))
   gru_id = gru_struc(:)%gru_id
 else
   err = nf90_inquire_dimension(ncID,dimID,len=fileGRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading gru dimension/'//trim(nf90_strerror(err)); return; end if
   ! read gru_id from netcdf file
   allocate(gru_id(fileGRU))
   err = nf90_inq_varid(ncID,"gruId",ncVarID);   if (err/=0) then; message=trim(message)//'problem finding gruId'; return; end if
   err = nf90_get_var(ncID,ncVarID,gru_id);      if (err/=0) then; message=trim(message)//'problem reading gruId'; return; end if
 end if

 ! Allocate the mapping arrays
 allocate(gruid_to_index(fileGRU), hrunc_to_index(fileGRU,maxval(gru_struc(:)%hruCount)))

 ! Populate the mapping arrays
 do i = 1, fileGRU
   gruid_to_index(i) = -1  ! Initialize with an invalid index
   do iGRU = 1, nGRU
     if (gru_struc(iGRU)%gru_id == gru_id(i)) then
       gruid_to_index(i) = iGRU
       do j = 1, gru_struc(iGRU)%hruCount
         hrunc_to_index(i,j) = -1 
         do iHRU = 1, fileHRU ! assumes HRUs are in GRU order
           if (gru_struc(iGRU)%hruInfo(j)%hru_id == hru_id(iHRU)) then
             hrunc_to_index(i,j) = iHRU
             exit
           endif
         end do
       end do ! HRU id loop
       exit
     endif
   end do
 end do

 ! loop through prognostic variables
 no_icond_enth=.false.
 do iVar = 1,size(prog_meta)
  ! skip variables that are computed later
  if(prog_meta(iVar)%varName=='scalarCanopyWat'           .or. &
     prog_meta(iVar)%varName=='spectralSnowAlbedoDiffuse' .or. &
     prog_meta(iVar)%varName=='scalarSurfaceTemp'         .or. &
     prog_meta(iVar)%varName=='mLayerVolFracWat'          .or. &
     prog_meta(iVar)%varName=='mLayerHeight'                   ) cycle
  ! get variable id
  err = nf90_inq_varid(ncID,trim(prog_meta(iVar)%varName),ncVarID)
  if(err/=nf90_noerr)then
   if(prog_meta(iVar)%varName=='DOMarea'              .or. &
      prog_meta(iVar)%varName=='DOMelev'              .or. &
      prog_meta(iVar)%varName=='scalarAblFrac'        .or. &
      prog_meta(iVar)%varName=='glacMass4AreaChange'       )then; err=nf90_noerr; cycle; endif ! backwards compatible, may be missing, correct in check_icond
   if(prog_meta(iVar)%varName=='scalarCanairEnthalpy' .or. &
      prog_meta(iVar)%varName=='scalarCanopyEnthalpy' .or. &  
      prog_meta(iVar)%varName=='mLayerEnthalpy'            )then; err=nf90_noerr; no_icond_enth=.true.; cycle; endif ! skip enthalpy variables if not in file
   call netcdf_err(err,message)
   message=trim(message)//': problem with getting variable id, var='//trim(prog_meta(iVar)%varName)
   return
  endif

  ! get variable dimension IDs
  select case (prog_meta(iVar)%varType)
   case (iLookVarType%scalarv); err = nf90_inq_dimid(ncID,trim(scalDimName)   ,dimID); call netcdf_err(err,message)
   case (iLookVarType%midSoil); err = nf90_inq_dimid(ncID,trim(midSoilDimName),dimID); call netcdf_err(err,message)
   case (iLookVarType%midToto); err = nf90_inq_dimid(ncID,trim(midTotoDimName),dimID); call netcdf_err(err,message)
   case (iLookVarType%ifcToto); err = nf90_inq_dimid(ncID,trim(ifcTotoDimName),dimID); call netcdf_err(err,message)
   case default
    message=trim(message)//"unexpectedVariableType[name='"//trim(prog_meta(iVar)%varName)//"';type='"//trim(get_varTypeName(prog_meta(iVar)%varType))//"']"
    err=20; return
  end select

  ! check errors
  if(err/=0)then
   message=trim(message)//': problem with dimension ids, var='//trim(prog_meta(iVar)%varName)
   return
  endif

  ! get the dimension length
  err = nf90_inquire_dimension(ncID,dimID,dimName,dimLen); call netcdf_err(err,message)
  if(err/=0)then; message=trim(message)//': problem getting the dimension length'; return; endif

  ! initialize the variable data
  allocate(varData3(fileDOM,fileHRU,dimLen),stat=err)
  allocate(varData2(fileHRU,dimLen),stat=err)
  if(err/=0)then; message=trim(message)//'problem allocating HRU variable data'; return; endif

  ! get data
  if(no_dom)then
   err = nf90_get_var(ncID,ncVarID,varData2); call netcdf_err(err,message)
  else
   err = nf90_get_var(ncID,ncVarID,varData3); call netcdf_err(err,message)
  endif
  if(err/=0)then; message=trim(message)//': problem getting the data for variable '//trim(prog_meta(iVar)%varName); return; endif

  ! store data in prognostics structure
  ! loop through GRUs
  has_glacier = .false.
  has_wetland = .false.
  do i = 1,fileGRU
   iGRU = gruid_to_index(i)
   do j = 1,gru_struc(iGRU)%hruCount
    iHRU = hrunc_to_index(i,j)- minval(hrunc_to_index(i,:)) + 1 ! assumes HRUs are in GRU order
    iHRU_global = hrunc_to_index(i,j)
    do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount
     ! get the number of layers
     nSnow = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSnow
     nLake = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nLake
     nSoil = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSoil
     nGlce = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nGlce
     nToto = nSnow + nLake + nSoil + nGlce
     if(nGlce>0) has_glacier = .true.
     if(nLake>0) has_wetland = .true.

     ! put the data into data structures and check that none of the values are set to nf90_fill_double
     if(no_dom)then
      select case (prog_meta(iVar)%varType)
       case (iLookVarType%scalarv)
        progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1) = varData2(iHRU_global,1)
        if(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1) - nf90_fill_double) < epsilon(varData2))then; err=20; endif
       case (iLookVarType%midSoil)
        progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nSoil) = varData2(iHRU_global,1:nSoil)
        if(any(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nSoil) - nf90_fill_double) < epsilon(varData2)))then; err=20; endif
       case (iLookVarType%midToto)
        progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nToto) = varData2(iHRU_global,1:nToto)
        if(any(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nToto) - nf90_fill_double) < epsilon(varData2)))then; err=20; endif
       case (iLookVarType%ifcToto)
        progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(0:nToto) = varData2(iHRU_global,1:nToto+1)
        if(any(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(0:nToto) - nf90_fill_double) < epsilon(varData2)))then; err=20; endif
       case default
        message=trim(message)//"unexpectedVariableType[name='"//trim(prog_meta(iVar)%varName)//"';type='"//trim(get_varTypeName(prog_meta(iVar)%varType))//"']"
        err=20; return
      end select
     else
      select case (prog_meta(iVar)%varType)
       case (iLookVarType%scalarv)
        progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1)       = varData3(iDOM,iHRU_global,1)
        if(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1) - nf90_fill_double) < epsilon(varData3))then; err=20; endif
       case (iLookVarType%midSoil)
        progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nSoil) = varData3(iDOM,iHRU_global,1:nSoil)
        if(any(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nSoil) - nf90_fill_double) < epsilon(varData3)))then; err=20; endif
       case (iLookVarType%midToto)
        progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nToto) = varData3(iDOM,iHRU_global,1:nToto)
        if(any(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nToto) - nf90_fill_double) < epsilon(varData3)))then; err=20; endif
       case (iLookVarType%ifcToto)
        progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(0:nToto) = varData3(iDOM,iHRU_global,1:nToto+1)
        if(any(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(0:nToto) - nf90_fill_double) < epsilon(varData3)))then; err=20; endif
       case default
        message=trim(message)//"unexpectedVariableType[name='"//trim(prog_meta(iVar)%varName)//"';type='"//trim(get_varTypeName(prog_meta(iVar)%varType))//"']"
        err=20; return
      end select
     endif

     if(err==20)then; message=trim(message)//"data set to the fill value (name='"//trim(prog_meta(iVar)%varName)//"')"; return; endif

     if(iVar==size(prog_meta))then ! last variable in the loop, so we can correct prognostic variables
      ! make sure snow albedo is not negative
      if(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1) < 0._rkind)then
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMax)%dat(1)
      endif

      ! initialize the spectral albedo
      progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1:nSpecBand) = progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1)

      ! make sure canopy ice + liq is positive, otherwise add liquid water to canopy and make total water consistent later
      if( (progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarCanopyLiq)%dat(1) + progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarCanopyIce)%dat(1)) < verySmaller*10._rkind)then
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarCanopyLiq)%dat(1) = verySmaller*10._rkind
       print*, 'WARNING: Canopy water is zero ... setting canopy liquid water to a tiny value.'
      endif
     endif ! (if last variable in the loop)

    end do ! iDOM
   end do ! iHRU
  end do ! iGRU

  ! deallocate storage vector for next variable
  deallocate(varData2,varData3,stat=err)
  if(err/=0)then; message=trim(message)//'problem deallocating HRU variable data'; return; endif

 end do ! end looping through prognostic variables (iVar)

 ! --------------------------------------------------------------------------------------------------------
 ! (2) set number of layers
 ! --------------------------------------------------------------------------------------------------------
 do i = 1,fileGRU
  iGRU = gruid_to_index(i)
  do j = 1,gru_struc(iGRU)%hruCount
   iHRU = hrunc_to_index(i,j)- minval(hrunc_to_index(i,:)) + 1 ! assumes HRUs are in GRU order
   do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount

    ! save the number of layers
    nSnow = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSnow
    nLake = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nLake
    nSoil = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSoil
    nGlce = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nGlce
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nSnow)%dat(1)   = nSnow
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nLake)%dat(1)   = nLake
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nSoil)%dat(1)   = nSoil
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nGlce)%dat(1)   = nGlce
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nLayers)%dat(1) = nSnow + nLake + nSoil + nGlce

    ! define layers that will not have a change in total water content
    noThetaChange = 0
    if(nGlce>0)then
      noThetaChange = nGlce - 1 ! This is a hard-coded value saying only the top glacier layer has change in theta can be increasing in the future to allow more glacier layers to have a change
      ! need at least one glacier top layer with a theta change
      if(noThetaChange>=nGlce)then; err=20; message=trim(message)//'number of glacier ice layers without a change in total water content is not less than the number of glacier ice layers'; return; endif
    endif
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%noThetaChange)%dat(1) = noThetaChange

    ! set layer type
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%layerType)%dat(1:nSnow) = iname_snow
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%layerType)%dat((nSnow+1):(nSnow+nLake)) = iname_lake
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%layerType)%dat((nSnow+nLake+1):(nSnow+nLake+nSoil)) = iname_soil
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%layerType)%dat((nSnow+nLake+nSoil+1):(nSnow+nLake+nSoil+nGlce)) = iname_glce
   end do
  end do
 end do

 ! --------------------------------------------------------------------------------------------------------
 ! (3) update soil layers (diagnostic variables)
 ! --------------------------------------------------------------------------------------------------------
 ! loop through GRUs and HRUs
 do i = 1,fileGRU
  iGRU = gruid_to_index(i)
  do j = 1,gru_struc(iGRU)%hruCount
   iHRU = hrunc_to_index(i,j)- minval(hrunc_to_index(i,:)) + 1 ! assumes HRUs are in GRU order
   do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount

    ! loop through soil layers
    do iLayer = 1,indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nSoil)%dat(1)

     ! get layer in the total vector
     jLayer = iLayer+indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nSnow)%dat(1)+indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nLake)%dat(1)

     ! update soil layers
     call updateSoil(&
                    ! input
                    progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerTemp          )%dat(jLayer),& ! intent(in): temperature vector (K)
                    progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerMatricHead    )%dat(iLayer),& ! intent(in): matric head (m)
                    mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%vGn_alpha          )%dat(iLayer),& ! intent(in): van Genutchen "alpha" parameter
                    mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%vGn_n              )%dat(iLayer),& ! intent(in): van Genutchen "n" parameter
                    mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%theta_sat          )%dat(iLayer),& ! intent(in): soil porosity (-)
                    mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%theta_res          )%dat(iLayer),& ! intent(in): soil residual volumetric water content (-)
                    1._rkind - 1._rkind/mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%vGn_n)%dat(iLayer),& ! intent(in): van Genutchen "m" parameter (-)
                    ! output
                    progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerVolFracWat    )%dat(jLayer),& ! intent(out): volumetric fraction of total water (-)
                    progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerVolFracLiq    )%dat(jLayer),& ! intent(out): volumetric fraction of liquid water (-)
                    progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerVolFracIce    )%dat(jLayer),& ! intent(out): volumetric fraction of ice (-)
                    err,cmessage)                                                                            ! intent(out): error control
     if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

    end do  ! looping through soil layers
   end do  ! looping through DOMs
  end do  ! looping through HRUs
 end do  ! looping through GRUs

 ! --------------------------------------------------------------------------------------------------------
 ! (2) now get the basin variable(s)
 ! --------------------------------------------------------------------------------------------------------

 ! get dimension of time delay histogram (TDH) from initial conditions file
 err = nf90_inq_dimid(ncID,"tdh",dimID)
 if(err/=nf90_noerr)then
  write(*,*) 'WARNING: routingRunoffFuture is not in the initial conditions file ... using zeros'  ! previously created in var_derive.f90
  err=nf90_noerr    ! reset this err

 else
  ! the state file *does* have the basin variable(s), so process them
  err = nf90_inquire_dimension(ncID,dimID,len=nTDH)
  if(err/=nf90_noerr)then; message=trim(message)//'problem reading tdh dimension from initial condition file/'//trim(nf90_strerror(err)); return; end if

  ! check vs hardwired value set in globalData.f90
  if(nTDH /= nTimeDelay)then
   write(*,*) 'tdh=',nTDH,' nTimeDelay=',nTimeDelay
   message=trim(message)//': state file time delay dimension tdh does not match summa expectation of nTimeDelay set in globalData()'
   return
  endif

  ! loop through specific basin variables (currently 1 but loop provided to enable inclusion of others)
  nrdx = (/iLookBVAR%routingRunoffFuture/)   ! array of desired variable indices
  do i = 1,size(nrdx)
   iVar = nrdx(i)

   ! get tdh dimension Id in file (should be 'tdh')
   err = nf90_inq_dimid(ncID,trim(tdhDimName), dimID)
   if(err/=0)then; message=trim(message)//': problem with dimension ids for tdh vars'; return; endif

   ! get the tdh dimension length (dimName and dimLen are outputs of this call)
   err = nf90_inquire_dimension(ncID,dimID,dimName,dimLen); call netcdf_err(err,message)
   if(err/=0)then; message=trim(message)//': problem getting the dimension length for tdh vars'; return; endif

   ! get tdh-based variable id
   err = nf90_inq_varid(ncID,trim(bvar_meta(iVar)%varName),ncVarID); call netcdf_err(err,message)
   if(err/=0)then; message=trim(message)//': problem with getting basin variable id, var='//trim(bvar_meta(iVar)%varName); return; endif

   ! initialize the tdh variable data
   allocate(varData2(fileGRU,dimLen))

   ! get data
   err = nf90_get_var(ncID,ncVarID,varData2); call netcdf_err(err,message)
   if(err/=0)then; message=trim(message)//': problem getting the data'; return; endif

   ! store data in basin var (bvar) structure
   do j = 1,fileGRU
    iGRU = gruid_to_index(j)
    ! put the data into data structures
    bvarData%gru(iGRU)%var(iVar)%dat(1:nTDH) = varData2((j+startGRU-1),1:nTDH)
    ! check whether the first values is set to nf90_fill_double
    if(any(abs(bvarData%gru(iGRU)%var(iVar)%dat(1:nTDH) - nf90_fill_double) < epsilon(varData2)))then; err=20; endif
    if(err==20)then; message=trim(message)//"data set to the fill value (name='"//trim(bvar_meta(iVar)%varName)//"')"; return; endif
   end do ! end iGRU loop

   ! deallocate temporary data array for next variable
   deallocate(varData2)

  end do ! end looping through basin variables
 endif  ! end if case for tdh variables being in init. cond. file

 if (has_glacier)then

  ! get variables that are specific to glaciers but not by glacier (by GRU)
  iVar = iLookBVAR%basin__GlacierStorage
  err = nf90_inq_varid(ncID,trim(bvar_meta(iVar)%varName),ncVarID)
  if(err/=0)then; message=trim(message)//': problem with getting basin variable id, var='//trim(bvar_meta(iVar)%varName); return; endif

  ! initialize the glac variable data
  allocate(varData(fileGRU))

  ! get data
  err = nf90_get_var(ncID,ncVarID,varData); call netcdf_err(err,message)
  if(err/=0)then; message=trim(message)//': problem getting the data'; return; endif

  ! store data in basin var (bvar) structure
  do j = 1,fileGRU
    iGRU = gruid_to_index(j)
    ! put the data into data structures
    bvarData%gru(iGRU)%var(iVar)%dat(1) = varData(j+startGRU-1)
    ! check whether the first values is set to nf90_fill_double
    if(abs(bvarData%gru(iGRU)%var(iVar)%dat(1) - nf90_fill_double) < epsilon(varData))then; err=20; endif
    if(err==20)then; message=trim(message)//"data set to the fill value (name='"//trim(bvar_meta(iVar)%varName)//"')"; return; endif
   end do ! end iGRU loop

   ! deallocate temporary data array for next variable
   deallocate(varData)

  ! get dimension of basin glac variables from initial conditions file (glaciers)
  err = nf90_inq_dimid(ncID,"glac",dimID) ! max number of glaciers in any GRU

  if(err/=nf90_noerr)then
   write(*,*) 'WARNING: glac is not in the initial conditions file ... assuming number of glaciers is the same as in attribute file and zero area for all glaciers'
    do iGRU = 1,nGRU
      nGlacier = gru_struc(iGRU)%nGlacier ! get dimension of basin glacier variables from attribute file, per GRU
      bvarData%gru(iGRU)%var(iLookBVAR%glacAblArea)%dat(1:nGlacier) = 0._rkind
      bvarData%gru(iGRU)%var(iLookBVAR%glacAccArea)%dat(1:nGlacier) = 0._rkind
      if (.not. allocated(gru_struc(iGRU)%glacInfo)) then
        allocate(gru_struc(iGRU)%glacInfo(gru_struc(iGRU)%nGlacier))
      endif
      gru_struc(iGRU)%glacInfo(1:nGlacier)%glac_id = int8Missing  ! set grid information (id) later
    end do
    err=nf90_noerr    ! reset this err
  
  else
   ! the state file *does* have the basin variable(s), so process them
   err = nf90_inquire_dimension(ncID,dimID,len=fileglac);
   if(err/=nf90_noerr)then; message=trim(message)//'problem reading glac dimension from initial condition file/'//trim(nf90_strerror(err)); return; end if

   ! read glac_id from netcdf file
   allocate(glac_id(filegru,fileglac))
   err = nf90_inq_varid(ncID,"glacId",ncVarID);   if (err/=0) then; message=trim(message)//'problem finding glacId'; return; end if
   err = nf90_get_var(ncID,ncVarID,glac_id);      if (err/=0) then; message=trim(message)//'problem reading glacId'; return; end if

   do iGRU = 1,nGRU
     nGlacier = gru_struc(iGRU)%nGlacier ! get dimension of basin glacier variables from attribute file, per GRU
     if (.not. allocated(gru_struc(iGRU)%glacInfo)) then
       allocate(gru_struc(iGRU)%glacInfo(gru_struc(iGRU)%nGlacier))
     endif
     gru_struc(iGRU)%glacInfo(1:nGlacier)%glac_id = glac_id(iGRU,1:nGlacier)  ! set grid information (id)
   end do

   ! loop through specific basin variables
   ngdx = (/iLookBVAR%glacAblArea,iLookBVAR%glacAccArea, iLookBVAR%glacIceRunoffFuture,iLookBVAR%glacSnowRunoffFuture,iLookBVAR%glacFirnRunoffFuture/)   ! array of desired variable indices
   do i = 1,size(ngdx)
    iVar = ngdx(i)

    ! get glac-based variable id
    err = nf90_inq_varid(ncID,trim(bvar_meta(iVar)%varName),ncVarID)
    if(err/=0)then
      if (iVar == iLookBVAR%glacIceRunoffFuture)then ! either all glacier runoff variables are in the file or none
        write(*,*) 'WARNING: glac(Ice,Snow,Firn)RunoffFuture is not in the initial conditions file ... using zeros'  ! previously created in var_derive.f90
        err=nf90_noerr    ! reset this err
        exit ! exit the loop, don't need to check the other glacier runoff variables
      else
        message=trim(message)//': problem with getting basin variable id, var='//trim(bvar_meta(iVar)%varName); return
      endif
    endif

    ! initialize the glac variable data
    allocate(varData2(fileGRU,fileglac))

    ! get data
    err = nf90_get_var(ncID,ncVarID,varData2); call netcdf_err(err,message)
    if(err/=0)then; message=trim(message)//': problem getting the data'; return; endif

    ! store data in basin var (bvar) structure
    do j = 1,fileGRU
     iGRU = gruid_to_index(j)
     nGlacier = gru_struc(iGRU)%nGlacier ! get dimension of basin glacier variables from attribute file, per GRU
     ! put the data into data structures
     bvarData%gru(iGRU)%var(iVar)%dat(1:nGlacier) = varData2(j+startGRU-1,1:nGlacier)
     ! check whether the first values is set to nf90_fill_double
     if(any(abs(bvarData%gru(iGRU)%var(iVar)%dat(1:nGlacier) - nf90_fill_double) < epsilon(varData2)))then; err=20; endif
     if(err==20)then; message=trim(message)//"data set to the fill value (name='"//trim(bvar_meta(iVar)%varName)//"')"; return; endif
    end do ! end iGRU loop

    deallocate(varData2) ! deallocate temporary data array for next variable
   end do ! end looping through basin variables
  endif  ! end if case for tdh variables being in init. cond. file
  deallocate(glac_id)
 endif  ! end if has glacier
 deallocate(hru_id,gru_id,gruid_to_index,hrunc_to_index)
 call nc_file_close(ncID,err,cmessage)
 if(err/=0)then;message=trim(message)//trim(cmessage);return;end if

 end subroutine read_icond

 ! ************************************************************************************************
 ! public subroutine read_icond: read model initial conditions
 ! ************************************************************************************************
 subroutine read_icondGlac(iconGlacFile,                  & ! intent(in):    name of glacier initial conditions file (surface topography)
                           nGRU,                          & ! intent(in):    number of GRUs
                           bvarData,                      & ! intent(inout): basin variable data structure
                           gridData,                      & ! intent(inout): grid data structure
                           err,message)                     ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nr_type
 USE var_lookup,only:iLookVarType                       ! variable lookup structure
 USE var_lookup,only:iLookGRID                          ! variable glacier grid 
 USE var_lookup,only:iLookBVAR                          ! variable basin variables
 USE globalData,only:grid_meta                          ! metadata for grid variables
 USE globalData,only:gru_struc                          ! gru-hru mapping structures
 USE globalData,only:startGRU                           ! index of first gru for parallel runs
 USE globalData,only:maxGlaciers                        ! maximum number of glaciers in any GRU
 USE netcdf_util_module,only:nc_file_open               ! open netcdf file
 USE netcdf_util_module,only:nc_file_close              ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                 ! netcdf error handling
 USE data_types,only:gru_grid_double                    ! full double precision structure
 USE data_types,only:gru_doubleVec                      ! gru-length double precision structure (basin variables)
 USE get_ixName_module,only:get_varTypeName             ! to access type strings for error messages

 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)           ,intent(in)        :: iconGlacFile             ! name of netcdf file containing the glacier initial conditions
 integer(i4b)           ,intent(in)        :: nGRU                     ! number of grouped response units in simulation domain
 type(gru_doubleVec)    ,intent(in)        :: bvarData                 ! basin variable data structure
 type(gru_grid_double)  ,intent(inout)     :: gridData                 ! grid data structure
 integer(i4b)           ,intent(out)       :: err                      ! error code
 character(*)           ,intent(out)       :: message                  ! returned error message
 ! locals
 character(len=256)                        :: cmessage                 ! downstream error message
 integer(i4b)                              :: fileGRU                  ! number of GRUs in file
 integer(i4b)                              :: filexgrid                ! number of xgrid points in file
 integer(i4b)                              :: fileygrid                ! number of ygrid points in file
 integer(i4b)                              :: ny,nx                    ! number of grid points for a glacier
 integer(i4b)                              :: filegrid                 ! max number of glacier grids in any GRU
 integer(i4b)                              :: iVar,i,j,k               ! loop indices
 integer(i4b),dimension(2)                 :: ngdx                     ! intermediate array of loop indices
 integer(i4b)                              :: iGRU,iGrid               ! loop index
 integer(i4b)                              :: dimID                    ! varible dimension ids
 integer(i4b)                              :: ncVarID                  ! variable ID in netcdf file
 integer(i4b)                              :: ncID                     ! netcdf file ID
 real(rkind),allocatable                   :: varData(:,:,:,:)         ! variable data storage
 integer(i8b),allocatable                  :: gru_id(:)                ! GRU id
 integer(i8b),allocatable                  :: grid_id(:,:)             ! Glacier id
 integer(i4b),allocatable                  :: gruid_to_index(:)        ! mapping from gru_id to index in gridData
 integer(i4b),allocatable                  :: glacid_to_index(:,:)     ! mapping from glac_id to index in gridData
 integer(i4b),allocatable                  :: gridid_to_index(:,:)     ! mapping from grid_id to index in gridData
 real(rkind),parameter                     :: areaTol=1.e-2_rkind      ! tolerance to address precision issues in glacier area summation
 real(rkind),parameter                     :: thick4area=0.1_rkind ! an arbitrary small threshold for glacier thickness to be considered as glacier area
 real(rkind)                               :: area,area_grid           ! glacier area for a single glacier (m2)
 ! --------------------------------------------------------------------------------------------------------
 ! Start procedure here
 err=0; message="read_icondGlac/"

 ! --------------------------------------------------------------------------------------------------------
 ! (1) read the file
 ! --------------------------------------------------------------------------------------------------------
 ! open netcdf file
 call nc_file_open(iconGlacFile,nf90_nowrite,ncID,err,cmessage)
 if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

 ! check if the file has the GRU dimension
 err = nf90_inq_dimid(ncID,"gru",dimID);               if(err/=nf90_noerr)then; message=trim(message)//'problem finding gru dimension/'//trim(nf90_strerror(err)); return; end if  
 err = nf90_inquire_dimension(ncID,dimID,len=fileGRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading gru dimension/'//trim(nf90_strerror(err)); return; end if

 ! read gru_id from netcdf file
 allocate(gru_id(fileGRU))
 err = nf90_inq_varid(ncID,"gruId",ncVarID);   if (err/=0) then; message=trim(message)//'problem finding gruId'; return; end if
 err = nf90_get_var(ncID,ncVarID,gru_id);      if (err/=0) then; message=trim(message)//'problem reading gruId'; return; end if
    
 ! get dimension of basin glacier variables from initial conditions file
 err = nf90_inq_dimid(ncID,"grid",dimID) ! max number of grids in any GRU
 err = nf90_inquire_dimension(ncID,dimID,len=filegrid);
 if(err/=nf90_noerr)then; message=trim(message)//'problem reading grid dimension from initial condition file/'//trim(nf90_strerror(err)); return; end if

 ! read grid_id from netcdf file
 allocate(grid_id(fileGRU,filegrid))
 err = nf90_inq_varid(ncID,"gridId",ncVarID); if (err/=0) then; message=trim(message)//'problem finding gridId'; return; end if
 err = nf90_get_var(ncID,ncVarID,grid_id);    if (err/=0) then; message=trim(message)//'problem reading gridId'; return; end if

 ! get ygrid dimension of whole file
 err = nf90_inq_dimid(ncID,"ygrid",dimId);                   if(err/=nf90_noerr)then; message=trim(message)//'problem finding ygrid dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncID, dimId, len = fileygrid); if(err/=nf90_noerr)then; message=trim(message)//'problem reading ygrid dimension/'//trim(nf90_strerror(err)); return; end if
 ! get xgrid dimension of whole file
 err = nf90_inq_dimid(ncID,"xgrid",dimId);                   if(err/=nf90_noerr)then; message=trim(message)//'problem finding xgrid dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncID, dimId, len = filexgrid); if(err/=nf90_noerr)then; message=trim(message)//'problem reading xgrid dimension/'//trim(nf90_strerror(err)); return; end if

 ! Allocate the mapping array
 allocate(gruid_to_index(fileGRU),glacid_to_index(fileGRU,maxGlaciers),gridid_to_index(fileGRU,filegrid))

 ! Populate the mapping arrays
 do i = 1, fileGRU
   gruid_to_index(i) = -1  ! Initialize with an invalid index
   do iGRU = 1, nGRU
     if (gru_struc(iGRU)%gru_id == gru_id(i)) then
       gruid_to_index(i) = iGRU
       do j = 1, gru_struc(iGRU)%nGlacier
         glacid_to_index(i,j) = -1
         if (gru_struc(iGRU)%glacInfo(j)%glac_id==int8Missing) then
           gru_struc(iGRU)%glacInfo(j)%glac_id = grid_id(i,j) ! if there is no glac_id in the intial conditions file, use grid_id
           glacid_to_index(i,j) = j
         else
           do iGrid = 1, size(gru_struc(iGRU)%gridInfo(:)%grid_id)
             if (gru_struc(iGRU)%glacInfo(j)%glac_id == grid_id(i,iGrid)) then
               glacid_to_index(i,j) = iGrid
               exit
             endif
           end do
         endif
         if (glacid_to_index(i,j) == -1) then
           message=trim(message)//'glacier glacId needs to match a gridId'; err=20; return
         endif
       end do ! glac id loop
       do j = 1, size(gru_struc(iGRU)%gridInfo(:)%grid_id)
         gridid_to_index(i,j) = -1
         do iGrid = 1, size(gru_struc(iGRU)%gridInfo(:)%grid_id)
           if (gru_struc(iGRU)%gridInfo(j)%grid_id == grid_id(i,iGrid)) then
             gridid_to_index(i,j) = iGrid
             exit
           endif
         end do  ! grid id loop
       end do ! grid id loop
       exit
     endif
   end do ! gru id loop
 end do ! fileGRU loop

 ! loop through specific basin grid variables
 ngdx = (/iLookGRID%surface_elev, iLookGRID%debris_thick/)   ! array of desired variable indices
 do i = 1,size(ngdx)
  iVar = ngdx(i)

  ! get grid-based variable id
  err = nf90_inq_varid(ncID,trim(grid_meta(iVar)%varName),ncVarID)
  if(err/=0)then; message=trim(message)//': problem with getting grid variable id, var= '//trim(grid_meta(iVar)%varName); return; endif

  ! initialize the grid variable data
  allocate(varData(fileGRU,filegrid,filexgrid,fileygrid),stat=err)
  if(err/=0)then; print*, 'err= ',err; message=trim(message)//'problem allocating GRU variable data'; return; endif

  ! get data
  err = nf90_get_var(ncID,ncVarID,varData); call netcdf_err(err,message)
  if(err/=0)then; message=trim(message)//': problem getting the data'; return; endif

  ! store data in grid structure
  do k = 1, fileGRU
    iGRU = gruid_to_index(k)
    do j = 1, size(gru_struc(iGRU)%gridInfo(:)%grid_id)
      iGrid = gridid_to_index(k,j)
      ny = gru_struc(iGRU)%gridInfo(iGrid)%ny
      nx = gru_struc(iGRU)%gridInfo(iGrid)%nx
      ! put the data into data structures
      gridData%gru(iGRU)%grid(iGrid)%var(iVar)%dat2(1:nx,1:ny) = varData(k+startGRU-1,j,1:nx,1:ny) 
    enddo
    ! correct for negative glacier thickness and note area differences with grid approximation
    if (iVar==iLookGRID%surface_elev) then
      do j = 1, gru_struc(iGRU)%nGlacier
        iGrid = glacid_to_index(k,j)   
        ny = gru_struc(iGRU)%gridInfo(iGrid)%ny
        nx = gru_struc(iGRU)%gridInfo(iGrid)%nx
        associate(& 
          area_abl     => bvarData%gru(k)%var(iLookBVAR%glacAblArea)%dat(j)                          ,&
          area_acc     => bvarData%gru(k)%var(iLookBVAR%glacAccArea)%dat(j)                          ,&
          bed_elev     => gridData%gru(iGRU)%grid(iGrid)%var(iLookGRID%bed_elev)%dat2(1:nx,1:ny)     ,&
          surface_elev => gridData%gru(iGRU)%grid(iGrid)%var(iLookGRID%surface_elev)%dat2(1:nx,1:ny) ,&
          dx           => gru_struc(iGRU)%gridInfo(iGrid)%dx                                         ,&
          dy           => gru_struc(iGRU)%gridInfo(iGrid)%dy                                          &
          )
          area = area_abl + area_acc
          ! correct if glacier thickness is negative
          surface_elev = merge(bed_elev, surface_elev, surface_elev - bed_elev <= thick4area)
          area_grid = sum(merge(0._rkind, dx*dy, surface_elev - bed_elev <= 0._rkind))
          ! check if area is significantly different from grid approximation
          ! this is okay for the first year as it is considered spin up, subsequent years will use the grid data to compute the area
          if (area_grid>0._rkind) then 
            if (abs(area/area_grid-1._rkind) >areaTol) then
              write(*,*) 'WARNING: Area of glacier ',iGrid,' in GRU ',iGRU,' starts at ', area/area_grid, ' times the grid approximation but will be calculated from the grid data after a year.'
              write(*,*) 'If this is not expected, check the glacier grid data in the initial grid conditions file (perhaps the grid should be refined or the glacier area should be adjusted).'
            endif      
          endif
        end associate
      end do
    endif
    ! check whether the first values is set to nf90_fill_double
    if(any(abs(gridData%gru(iGRU)%grid(iGrid)%var(iVar)%dat2(1:nx,1:ny) - nf90_fill_double) < epsilon(varData)))then; err=20; endif
    if(err==20)then; message=trim(message)//"data set to the fill value (name='"//trim(grid_meta(iVar)%varName)//"')"; return; endif
  end do ! end iGRU loop

  ! deallocate temporary data array for next variable
  deallocate(varData, stat=err)
  if(err/=0)then; message=trim(message)//'problem deallocating GRU variable data'; return; endif
 enddo ! end looping through basin variables
 deallocate(gru_id,grid_id,gruid_to_index,glacid_to_index,gridid_to_index)
 call nc_file_close(ncID,err,cmessage)
 if(err/=0)then;message=trim(message)//trim(cmessage);return;end if

 end subroutine read_icondGlac


end module read_icond_module
