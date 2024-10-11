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
USE nrtype
USE netcdf
USE globalData,only: ixHRUfile_min,ixHRUfile_max
USE globalData,only: nTimeDelay   ! number of hours in the time delay histogram
USE globalData,only: nSpecBand    ! number of spectral bands

! access domain types
USE globalData,only:upland             ! domain type for upland areas
USE globalData,only:glacAcc            ! domain type for glacier accumulation areas
USE globalData,only:glacAbl            ! domain type for glacier ablation areas
USE globalData,only:wetland            ! domain type for wetland areas

implicit none
private
public::read_icond
public::read_icond_nlayers


contains

 ! ************************************************************************************************
 ! public subroutine read_icond_nlayers: read model initial conditions file for number of snow/soil layers
 ! ************************************************************************************************
 subroutine read_icond_nlayers(iconFile,nGRU,nHRU,nDOM,indx_meta,err,message)
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nrtype
 USE var_lookup,only:iLookINDEX                        ! variable lookup structure
 USE globalData,only:gru_struc                         ! gru-hru mapping structures
 USE globalData,only:startGRU                          ! index of first gru for parallel runs
 USE netcdf_util_module,only:nc_file_close             ! close netcdf file
 USE netcdf_util_module,only:nc_file_open              ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                ! netcdf error handling
 USE data_types,only:var_info                          ! metadata
 implicit none

 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)  ,intent(in)   :: iconFile           ! name of input (restart) file
 integer(i4b)  ,intent(in)   :: nGRU               ! total # of GRUs in run space
 integer(i4b)  ,intent(in)   :: nHRU               ! total # of HRUs in run space
 integer(i4b)  ,intent(inout):: nDOM               ! total # of domains in run space
 type(var_info),intent(in)   :: indx_meta(:)       ! metadata
 integer(i4b)  ,intent(out)  :: err                ! error code
 character(*)  ,intent(out)  :: message            ! returned error message
 ! locals
 integer(i4b)                :: ncID               ! netcdf file id
 integer(i4b)                :: ixFile             ! index in file
 integer(i4b)                :: dimID              ! netcdf file dimension id
 integer(i4b)                :: varID              ! netcdf variable id
 integer(i4b)                :: fileHRU            ! number of HRUs in netcdf file
 integer(i4b)                :: fileDOM            ! number of domains in netcdf file
 integer(i4b)                :: snowID, soilID     ! netcdf variable ids
 integer(i4b)                :: iceID, lakeID      ! netcdf variable ids
 integer(i4b)                :: iGRU, iHRU, iDOM   ! loop indexes
 integer(i4b)                :: iHRU_global        ! index of HRU in the netcdf file
 logical(lgt)                :: no_iceData         ! flag that no ice data in icond
 logical(lgt)                :: no_lakeData        ! flag that no lake data in icond
 integer(i4b),allocatable    :: snowData(:,:)      ! number of snow layers in all HRUs
 integer(i4b),allocatable    :: soilData(:,:)      ! number of soil layers in all HRUs
 integer(i4b),allocatable    :: iceData(:,:)       ! number of ice layers in all HRUs
 integer(i4b),allocatable    :: lakeData(:,:)      ! number of lake layers in all HRUs
 integer(i8b),allocatable    :: dom_type(:,:)      ! read domain type in from netcdf file
 character(len=256)          :: cmessage           ! downstream error message

 ! --------------------------------------------------------------------------------------------------------
 ! initialize error message
 err=0
 message = 'read_icond_nlayers/'
 no_iceData = .false.
 no_lakeData = .false.

 ! open netcdf file
 call nc_file_open(iconFile,nf90_nowrite,ncID,err,cmessage);
 if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

 ! get number of HRUs in file (the GRU variable(s), if present, are processed at the end)
 err = nf90_inq_dimid(ncID,"hru",dimId);               if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncID,dimId,len=fileHRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru dimension/'//trim(nf90_strerror(err)); return; end if

! get number of domains with type in file, if present
 err = nf90_inq_dimid(ncID,"dom",dimId)               
 if(err/=nf90_noerr)then
  fileDOM = 1 ! backwards compatible, just upland domain
  allocate(dom_type(1,fileHRU))
  dom_type = upland
 else
  err = nf90_inquire_dimension(ncID,dimId,len=fileDOM); if(err/=nf90_noerr)then; message=trim(message)//'problem reading dom dimension/'//trim(nf90_strerror(err)); return; end if
  ! read dom_type from netcdf file
  allocate(dom_type(fileDOM,fileHRU))
  err = nf90_inq_varid(ncID,"domType",varID);  if (err/=0) then; message=trim(message)//'problem finding domType'; return; end if
  err = nf90_get_var(ncID,varID,dom_type);     if (err/=0) then; message=trim(message)//'problem reading domType'; return; end if

 end if

 ! allocate storage for reading from file (allocate entire file size, even when doing subdomain run)
 allocate(snowData(fileDOM,fileHRU))
 allocate(soilData(fileDOM,fileHRU))
 allocate( iceData(fileDOM,fileHRU))
 allocate(lakeData(fileDOM,fileHRU))
 snowData = 0
 soilData = 0
 iceData  = 0
 lakeData = 0

 do iHRU = 1,gru_struc(iGRU)%hruCount
  gru_struc(iGRU)%hruInfo(iHRU)%domCount = 1                                              ! upland domain always present, for changing size glaciers and lakes
  if (any(dom_type(1:fileDOM,gru_struc(iGRU)%hruInfo(iHRU)%hru_nc)==glacAcc)) &
    gru_struc(iGRU)%hruInfo(iHRU)%domCount = gru_struc(iGRU)%hruInfo(iHRU)%domCount + 2   ! accumulation and ablation domains possible
  if (any(dom_type(1:fileDOM,gru_struc(iGRU)%hruInfo(iHRU)%hru_nc)==wetland)) &
    gru_struc(iGRU)%hruInfo(iHRU)%domCount = gru_struc(iGRU)%hruInfo(iHRU)%domCount + 1   ! wetland domain possible
  allocate(gru_struc(iGRU)%hruInfo(iHRU)%domInfo(gru_struc(iGRU)%hruInfo(iHRU)%domCount)) ! allocate third level of gru to hru map
  gru_struc(iGRU)%hruInfo(iHRU)%domInfo(:)%dom_type = dom_type(:,gru_struc(iGRU)%hruInfo(iHRU)%hru_nc)
enddo

 ! get netcdf ids for the variables holding number of layers in each domain or hru
 err = nf90_inq_varid(ncID,trim(indx_meta(iLookINDEX%nSnow)%varName),snowID); call netcdf_err(err,message)
 err = nf90_inq_varid(ncID,trim(indx_meta(iLookINDEX%nLake)%varName),lakeID)
 if(err/=nf90_noerr ) no_lakeData = .true.
 err = nf90_inq_varid(ncID,trim(indx_meta(iLookINDEX%nSoil)%varName),soilID); call netcdf_err(err,message)
 err = nf90_inq_varid(ncID,trim(indx_meta(iLookINDEX%nIce)%varName),iceID)
 if(err/=nf90_noerr) no_iceData = .true.

 ! get nSnow and nSoil data (reads entire state file)
 err = nf90_get_var(ncID,snowID,snowData); call netcdf_err(err,message)
 err = nf90_get_var(ncID,soilID,soilData); call netcdf_err(err,message)
 if (.not. no_iceData)  err = nf90_get_var(ncID,iceID,iceData);   call netcdf_err(err,message)
 if (.not. no_lakeData) err = nf90_get_var(ncID,lakeID,lakeData); call netcdf_err(err,message)
 ixHRUfile_min=huge(1)
 ixHRUfile_max=0
 ! find the min and max hru indices in the state file
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   if(gru_struc(iGRU)%hruInfo(iHRU)%hru_nc < ixHRUfile_min) ixHRUfile_min = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
   if(gru_struc(iGRU)%hruInfo(iHRU)%hru_nc > ixHRUfile_max) ixHRUfile_max = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
  end do
 end do

 ! loop over grus in current run to update snow/soil layer information
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   iHRU_global = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
   do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount
    ixFile = iHRU_global
    gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSnow = snowData(iDOM,ixFile)
    gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nLake = lakeData(iDOM,ixFile)
    gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSoil = soilData(iDOM,ixFile)
    gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nIce  =  iceData(iDOM,ixFile)
   end do
  end do
 end do

 ! close file
 call nc_file_close(ncID,err,cmessage)
 if(err/=0)then;message=trim(message)//trim(cmessage);return;end if

 ! cleanup
 deallocate(snowData,lakeData,soilData,iceData,dom_type)

 end subroutine read_icond_nlayers


 ! ************************************************************************************************
 ! public subroutine read_icond: read model initial conditions
 ! ************************************************************************************************
 subroutine read_icond(iconFile,                      & ! intent(in):    name of initial conditions file
                       nGRU,                          & ! intent(in):    number of GRUs
                       nHRU,                          & ! intent(in):    number of HRUs
                       nDOM,                          & ! intent(in):    number of domains
                       mparData,                      & ! intent(in):    model parameters
                       progData,                      & ! intent(inout): model prognostic variables
                       bvarData,                      & ! intent(inout): model basin (GRU) variables
                       indxData,                      & ! intent(inout): model indices
                       no_icond_enth,                 & ! intent(out):   flag that enthalpy variables are not in the file
                       err,message)                     ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nrtype
 USE var_lookup,only:iLookVarType                       ! variable lookup structure
 USE var_lookup,only:iLookPROG                          ! variable lookup structure
 USE var_lookup,only:iLookPARAM                         ! variable lookup structure
 USE var_lookup,only:iLookBVAR                          ! variable lookup structure
 USE var_lookup,only:iLookINDEX                         ! variable lookup structure
 USE globalData,only:prog_meta                          ! metadata for prognostic variables
 USE globalData,only:bvar_meta                          ! metadata for basin (GRU) variables
 USE globalData,only:gru_struc                          ! gru-hru mapping structures
 USE globalData,only:startGRU                           ! index of first gru for parallel runs
 USE globalData,only:iname_soil,iname_snow,iname_ice,iname_lake ! named variables to describe the type of layer
 USE netcdf_util_module,only:nc_file_open               ! open netcdf file
 USE netcdf_util_module,only:nc_file_close              ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                 ! netcdf error handling
 USE data_types,only:gru_hru_dom_doubleVec              ! full double precision structure
 USE data_types,only:gru_hru_dom_intVec                 ! full integer structure
 USE data_types,only:gru_doubleVec                      ! gru-length double precision structure (basin variables)
 USE data_types,only:var_dlength                        ! double precision structure for a single HRU
 USE data_types,only:var_info                           ! metadata
 USE get_ixName_module,only:get_varTypeName             ! to access type strings for error messages
 USE updatState_module,only:updateSoil                  ! update soil states

 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)           ,intent(in)        :: iconFile                 ! name of netcdf file containing the initial conditions
 integer(i4b)           ,intent(in)        :: nGRU                     ! number of grouped response units in simulation domain
 integer(i4b)           ,intent(in)        :: nHRU                     ! number of hydrological response units in simulation domain
 integer(i4b)           ,intent(in)        :: nDOM                     ! number of domains in simulation domain
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
 integer(i4b)                              :: iVar, i                  ! loop indices
 integer(i4b),dimension(1)                 :: nrdx                      ! intermediate array of loop indices
 integer(i4b),dimension(5)                 :: ngdx                     ! intermediate array of loop indices
 integer(i4b)                              :: iGRU                     ! loop index
 integer(i4b)                              :: iHRU                     ! loop index
 integer(i4b)                              :: iDOM                     ! loop index
 integer(i4b)                              :: dimID                    ! varible dimension ids
 integer(i4b)                              :: ncVarID                  ! variable ID in netcdf file
 character(256)                            :: dimName                  ! not used except as a placeholder in call to inq_dim function
 integer(i4b)                              :: dimLen                   ! data dimensions
 integer(i4b)                              :: ncID                     ! netcdf file ID
 integer(i4b)                              :: ixFile                   ! index in file
 integer(i4b)                              :: iHRU_global              ! index of HRU in the netcdf file
 real(rkind),allocatable                   :: varData2(:,:)            ! variable data storage
 real(rkind),allocatable                   :: varData3(:,:,:)          ! variable data storage
 integer(i4b)                              :: nSnow,nLake,nSoil,nIce,nToto !# layers
 integer(i4b)                              :: nTDH                     ! number of points in time-delay 
 integer(i4b)                              :: nGlacier                 ! number of glaciers in basin (attribute files
 integer(i4b)                              :: nGlacier_max             ! max number of glaciers in any GRU
 integer(i4b)                              :: nWetland                 ! number of wetlands in basin
 integer(i4b)                              :: iLayer,jLayer            ! layer indices
 logical(lgt)                              :: has_glacier              ! flag for glacier presence in at least one GRU
 logical(lgt)                              :: has_wetland              ! flag for wetland/lake presence in at least one GRU
 ! currently only writing restart for progressive variables with these dimensions
 character(len=32),parameter               :: scalDimName   ='scalarv' ! dimension name for scalar data
 character(len=32),parameter               :: midSoilDimName='midSoil' ! dimension name for soil-only layers
 character(len=32),parameter               :: midTotoDimName='midToto' ! dimension name for layered varaiables
 character(len=32),parameter               :: ifcTotoDimName='ifcToto' ! dimension name for layered varaiables
 character(len=32),parameter               :: tdhDimName    ='tdh'     ! dimension name for time-delay basin variables
 character(len=32),parameter               :: nglDimName    ='ngl'     ! dimension name for glacier variables

 ! --------------------------------------------------------------------------------------------------------
 ! Start procedure here
 err=0; message="read_icond/"

 ! --------------------------------------------------------------------------------------------------------
 ! (1) read the file
 ! --------------------------------------------------------------------------------------------------------
 ! open netcdf file
 call nc_file_open(iconFile,nf90_nowrite,ncID,err,cmessage)
 if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

 ! get number of HRUs in file
 err = nf90_inq_dimid(ncID,"hru",dimID);               if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncID,dimID,len=fileHRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru dimension/'//trim(nf90_strerror(err)); return; end if

 ! get max number of DOMs any HRU in file, if present
 err = nf90_inq_dimid(ncID,"dom",dimId)               
 if(err/=nf90_noerr)then
  fileDOM = nDOM
 else
  err = nf90_inquire_dimension(ncID,dimId,len=fileDOM); if(err/=nf90_noerr)then; message=trim(message)//'problem reading dom dimension/'//trim(nf90_strerror(err)); return; end if
 end if

 ! loop through prognostic variables
 no_icond_enth=.false.
 do iVar = 1,size(prog_meta)

  ! skip variables that are computed later
  if(prog_meta(iVar)%varName=='scalarCanopyWat'           .or. &
     prog_meta(iVar)%varName=='spectralSnowAlbedoDiffuse' .or. &
     prog_meta(iVar)%varName=='scalarSurfaceTemp'         .or. &
     prog_meta(iVar)%varName=='mLayerVolFracWat'          .or. &
     prog_meta(iVar)%varName=='mLayerHeight'                   ) err=nf90_noerr; cycle

  ! get variable id
  err = nf90_inq_varid(ncID,trim(prog_meta(iVar)%varName),ncVarID)
  if(err/=nf90_noerr)then
   if(prog_meta(iVar)%varName=='DOMarea'              .or. &
      prog_meta(iVar)%varName=='DOMelev'                   ) err=nf90_noerr; cycle ! backwards compatible, may be missing, correct in check_icond
   if(prog_meta(iVar)%varName=='scalarCanairEnthalpy' .or. &
      prog_meta(iVar)%varName=='scalarCanopyEnthalpy' .or. &  
      prog_meta(iVar)%varName=='mLayerEnthalpy'            ) err=nf90_noerr; no_icond_enth=.true.; cycle ! skip enthalpy variables if not in file
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
  if(err/=0)then; message=trim(message)//'problem allocating HRU variable data'; return; endif

  ! get data
  err = nf90_get_var(ncID,ncVarID,varData3); call netcdf_err(err,message)
  if(err/=0)then; message=trim(message)//': problem getting the data for variable '//trim(prog_meta(iVar)%varName); return; endif

  ! store data in prognostics structure
  ! loop through GRUs
  has_glacier = .false.
  has_wetland = .false.
  do iGRU = 1,nGRU
   do iHRU = 1,gru_struc(iGRU)%hruCount
    iHRU_global = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
    do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount
     ! get the number of layers
     nSnow = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSnow
     nLake = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nLake
     nSoil = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSoil
     nIce  = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nIce
     nToto = nSnow + nLake + nSoil + nIce
     if(nIce>0) has_glacier = .true.
     if(nLake>0) has_wetland = .true.

     ixFile = iHRU_global
     ! put the data into data structures and check that none of the values are set to nf90_fill_double
     select case (prog_meta(iVar)%varType)
      case (iLookVarType%scalarv)
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1)       = varData3(iDOM,ixFile,1)
       if(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1) - nf90_fill_double) < epsilon(varData3))then; err=20; endif
      case (iLookVarType%midSoil)
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nSoil) = varData3(iDOM,ixFile,1:nSoil)
       if(any(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nSoil) - nf90_fill_double) < epsilon(varData3)))then; err=20; endif   
      case (iLookVarType%midToto)
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nToto) = varData3(iDOM,ixFile,1:nToto)
       if(any(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nToto) - nf90_fill_double) < epsilon(varData3)))then; err=20; endif
      case (iLookVarType%ifcToto)
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(0:nToto) = varData3(iDOM,ixFile,1:nToto+1)
       if(any(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(0:nToto) - nf90_fill_double) < epsilon(varData3)))then; err=20; endif
      case default
       message=trim(message)//"unexpectedVariableType[name='"//trim(prog_meta(iVar)%varName)//"';type='"//trim(get_varTypeName(prog_meta(iVar)%varType))//"']"
       err=20; return
     end select

     if(err==20)then; message=trim(message)//"data set to the fill value (name='"//trim(prog_meta(iVar)%varName)//"')"; return; endif

     ! make sure snow albedo is not negative
     if(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1) < 0._rkind)then
      progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMax)%dat(1)
     endif

     ! make sure canopy ice + liq is positive, otherwise add liquid water to canopy and make total water consistent later
     if( (progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarCanopyLiq)%dat(1) + progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarCanopyIce)%dat(1)) < 0.0001_rkind)then
      progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarCanopyLiq)%dat(1) = 0.0001_rkind
     endif

     ! initialize the spectral albedo
     progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1:nSpecBand) = progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1)

    end do ! iDOM
   end do ! iHRU
  end do ! iGRU

  ! deallocate storage vector for next variable
  deallocate(varData3, stat=err)
  if(err/=0)then; message=trim(message)//'problem deallocating HRU variable data'; return; endif

 end do ! end looping through prognostic variables (iVar)

 ! --------------------------------------------------------------------------------------------------------
 ! (2) set number of layers
 ! --------------------------------------------------------------------------------------------------------
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount

    ! save the number of layers
    nSnow = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSnow
    nLake = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nLake
    nSoil = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSoil
    nIce  = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nIce
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nSnow)%dat(1)   = nSnow
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nLake)%dat(1)   = nLake
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nSoil)%dat(1)   = nSoil
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nIce)%dat(1)    = nIce
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nLayers)%dat(1) = nSnow + nLake + nSoil + nIce

    ! set layer type
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%layerType)%dat(1:nSnow) = iname_snow
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%layerType)%dat((nSnow+1):(nSnow+nLake)) = iname_lake
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%layerType)%dat((nSnow+nLake+1):(nSnow+nLake+nSoil)) = iname_soil
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%layerType)%dat((nSnow+nLake+nSoil+1):(nSnow+nLake+nSoil+nIce)) = iname_ice

   end do
  end do
 end do

 ! --------------------------------------------------------------------------------------------------------
 ! (3) update soil layers (diagnostic variables)
 ! --------------------------------------------------------------------------------------------------------
 ! loop through GRUs and HRUs
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
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
  ! check if the file has the GRU dimension
  err = nf90_inq_dimid(ncID,"gru",dimID);    
  if(err/=nf90_noerr)then         
    write(*,*) 'WARNING: GRU is not in the initial conditions file ... assuming 1 GRU'
    fileGRU = 1
    err=nf90_noerr    ! reset this err
  else
    err = nf90_inquire_dimension(ncID,dimID,len=fileGRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading gru dimension/'//trim(nf90_strerror(err)); return; end if
  end if

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
   allocate(varData2(fileGRU,dimLen),stat=err)
   if(err/=0)then; print*, 'err= ',err; message=trim(message)//'problem allocating GRU variable data'; return; endif

   ! get data
   err = nf90_get_var(ncID,ncVarID,varData2); call netcdf_err(err,message)
   if(err/=0)then; message=trim(message)//': problem getting the data'; return; endif

   ! store data in basin var (bvar) structure
   do iGRU = 1,nGRU

    ! put the data into data structures
    bvarData%gru(iGRU)%var(iVar)%dat(1:nTDH) = varData2((iGRU+startGRU-1),1:nTDH)
    ! check whether the first values is set to nf90_fill_double
    if(any(abs(bvarData%gru(iGRU)%var(iVar)%dat(1:nTDH) - nf90_fill_double) < epsilon(varData2)))then; err=20; endif
    if(err==20)then; message=trim(message)//"data set to the fill value (name='"//trim(bvar_meta(iVar)%varName)//"')"; return; endif

   end do ! end iGRU loop

   ! deallocate temporary data array for next variable
   deallocate(varData2, stat=err)
   if(err/=0)then; message=trim(message)//'problem deallocating GRU variable data'; return; endif

  end do ! end looping through basin variables
 endif  ! end if case for tdh variables being in init. cond. file

 if (has_glacier)then
  ! get dimension of basin glacier variables from initial conditions file
  err = nf90_inq_dimid(ncID,"ngl",dimID) ! max number of glaciers in any GRU

  if(err/=nf90_noerr)then
   write(*,*) 'WARNING: ngl is not in the initial conditions file ... assuming same as attribute file and zero area for all glaciers'
    do iGRU = 1,nGRU
      nGlacier = gru_struc(iGRU)%nGlacier ! get dimension of basin glacier variables from attribute file, per GRU
      bvarData%gru(iGRU)%var(iLookBVAR%glacAblArea)%dat(1:nGlacier) = 0._rkind
      bvarData%gru(iGRU)%var(iLookBVAR%glacAccArea)%dat(1:nGlacier) = 0._rkind
    end do
    err=nf90_noerr    ! reset this err
  
  else
   ! the state file *does* have the basin variable(s), so process them
   err = nf90_inquire_dimension(ncID,dimID,len=nGlacier_max);
   if(err/=nf90_noerr)then; message=trim(message)//'problem reading ngl dimension from initial condition file/'//trim(nf90_strerror(err)); return; end if

   ! loop through specific basin variables
   ngdx = (/iLookBVAR%glacAblArea,iLookBVAR%glacAccArea,iLookBVAR%glacIceRunoffFuture,iLookBVAR%glacSnowRunoffFuture,iLookBVAR%glacFirnRunoffFuture/)   ! array of desired variable indices
   do i = 1,size(ngdx)
    iVar = ngdx(i)

    ! get ngl-based variable id
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

    ! initialize the ngl variable data
    allocate(varData2(fileGRU,nGlacier_max),stat=err)
    if(err/=0)then; print*, 'err= ',err; message=trim(message)//'problem allocating GRU variable data'; return; endif

    ! get data
    err = nf90_get_var(ncID,ncVarID,varData2); call netcdf_err(err,message)
    if(err/=0)then; message=trim(message)//': problem getting the data'; return; endif

    ! store data in basin var (bvar) structure
    do iGRU = 1,nGRU
     nGlacier = gru_struc(iGRU)%nGlacier ! get dimension of basin glacier variables from attribute file, per GRU
     ! put the data into data structures
     bvarData%gru(iGRU)%var(iVar)%dat(1:nGlacier) = varData2(iGRU+startGRU-1,1:nGlacier)
     ! check whether the first values is set to nf90_fill_double
     if(any(abs(bvarData%gru(iGRU)%var(iVar)%dat(1:nGlacier) - nf90_fill_double) < epsilon(varData2)))then; err=20; endif
     if(err==20)then; message=trim(message)//"data set to the fill value (name='"//trim(bvar_meta(iVar)%varName)//"')"; return; endif
    end do ! end iGRU loop

    ! deallocate temporary data array for next variable
    deallocate(varData2, stat=err)
    if(err/=0)then; message=trim(message)//'problem deallocating GRU variable data'; return; endif
   enddo ! end looping through basin variables
  endif  ! end if case for variables being in init. cond. file
 endif ! end if has glacier

 end subroutine read_icond

end module read_icond_module
