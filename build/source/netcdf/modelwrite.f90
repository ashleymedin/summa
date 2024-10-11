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

module modelwrite_module

! NetCDF types
USE netcdf
USE netcdf_util_module,only:netcdf_err  ! netcdf error handling function

! top-level data types
USE nrtype

! missing values
USE globalData,only: integerMissing, realMissing
! output constraints
USE globalData,only: maxSnowLayers      ! maximum number of snow layers
USE globalData,only: maxSoilLayers      ! maximum number of soil layers
USE globalData,only: maxIceLayers       ! maximum number of ice layers
USE globalData,only: maxLakeLayers      ! maximum number of lake layers
USE globalData,only: maxGlaciers        ! maximum number of glaciers
USE globalData,only: nTimeDelay         ! maximum number of time delay routing vectors
USE globalData,only: nSpecBand          ! maximum number of spectral bands
! provide access to global data
USE globalData,only: gru_struc          ! gru->hru mapping structure

! provide access to the derived types to define the data structures
USE data_types,only:&
                    ! final data vectors
                    dlength,             & ! var%dat
                    ilength,             & ! var%dat
                    ! no spatial dimension
                    var_i,               & ! x%var(:)            (i4b)
                    var_i8,              & ! x%var(:)            (i8b)
                    var_d,               & ! x%var(:)            (dp)
                    var_ilength,         & ! x%var(:)%dat        (i4b)
                    var_dlength,         & ! x%var(:)%dat        (dp)
                    ! no variable dimension
                    hru_i,               & ! x%hru(:)            (i4b)
                    hru_d,               & ! x%hru(:)            (dp)
                    ! gru dimension
                    gru_int,             & ! x%gru(:)%var(:)     (i4b)
                    gru_double,          & ! x%gru(:)%var(:)     (dp)
                    gru_intVec,          & ! x%gru(:)%var(:)%dat (i4b)
                    gru_doubleVec,       & ! x%gru(:)%var(:)%dat (dp)
                    ! gru+hru dimension
                    gru_hru_int,         & ! x%gru(:)%hru(:)%var(:)     (i4b)
                    gru_hru_int8,        & ! x%gru(:)%hru(:)%var(:)     (i8b)
                    gru_hru_double,      & ! x%gru(:)%hru(:)%var(:)     (dp)
                    gru_hru_intVec,      & ! x%gru(:)%hru(:)%var(:)%dat (i4b)
                    gru_hru_doubleVec,   & ! x%gru(:)%hru(:)%var(:)%dat (dp)
                    ! gru+hru+dom dimension
                    gru_hru_dom_int,     & ! x%gru(:)%hru(:)%dom(:)%var(:)     (i4b)
                    gru_hru_dom_int8,    & ! x%gru(:)%hru(:)%dom(:)%var(:)     integer(8)
                    gru_hru_dom_double,  & ! x%gru(:)%hru(:)%dom(:)%var(:)     (dp)
                    gru_hru_dom_intVec,  & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat (i4b)
                    gru_hru_dom_doubleVec  ! x%gru(:)%hru(:)%dom(:)%var(:)%dat (dp)


! vector lengths
USE var_lookup, only: maxvarFreq ! number of output frequencies
USE var_lookup, only: maxvarStat ! number of statistics
   

implicit none
private
public::writeParm
public::writeData
public::writeBasin
public::writeTime
public::writeRestart


contains

 ! **********************************************************************************************************
 ! public subroutine writeParm: write model parameters
 ! **********************************************************************************************************
 subroutine writeParm(iDOM,iSpatial,struct,meta,err,message)
 USE globalData,only:ncid                        ! netcdf file ids
 USE data_types,only:var_info                    ! metadata info
 USE var_lookup,only:iLookSTAT                   ! index in statistics vector
 USE var_lookup,only:iLookFREQ                   ! index in vector of model output frequencies
 implicit none

 ! declare input variables
 integer(i4b)  ,intent(in)   :: iDOM             ! domain index
 integer(i4b)  ,intent(in)   :: iSpatial         ! hydrologic response unit
 class(*)      ,intent(in)   :: struct           ! data structure
 type(var_info),intent(in)   :: meta(:)          ! metadata structure
 integer(i4b)  ,intent(out)  :: err              ! error code
 character(*)  ,intent(out)  :: message          ! error message
 ! local variables
 integer(i4b)                :: iVar             ! loop through variables

 ! initialize error control
 err=0;message="writeParm/"

 ! loop through local column model parameters
 do iVar = 1,size(meta)

  ! check that the variable is desired
  if (meta(iVar)%statIndex(iLookFREQ%timestep)==integerMissing) cycle

  ! initialize message
  message=trim(message)//trim(meta(iVar)%varName)//'/'

  select type (struct)
   class is (var_i)
     if (iDOM==0) err = nf90_put_var(ncid(iLookFREQ%timestep),meta(iVar)%ncVarID(iLookFREQ%timestep),(/struct%var(iVar)/),start=(/iSpatial/),count=(/1/))
     if (iDOM>0)  err = nf90_put_var(ncid(iLookFREQ%timestep),meta(iVar)%ncVarID(iLookFREQ%timestep),(/struct%var(iVar)/),start=(/iDOM,iSpatial/),count=(/1,1/))
    class is (var_i8)
     if (iDOM==0) err = nf90_put_var(ncid(iLookFREQ%timestep),meta(iVar)%ncVarID(iLookFREQ%timestep),(/struct%var(iVar)/),start=(/iSpatial/),count=(/1/))
     if (iDOM>0)  err = nf90_put_var(ncid(iLookFREQ%timestep),meta(iVar)%ncVarID(iLookFREQ%timestep),(/struct%var(iVar)/),start=(/iDOM,iSpatial/),count=(/1,1/))
    class is (var_d)
     if (iDOM==0) err = nf90_put_var(ncid(iLookFREQ%timestep),meta(iVar)%ncVarID(iLookFREQ%timestep),(/struct%var(iVar)/),start=(/iSpatial/),count=(/1/))
     if (iDOM>0)  err = nf90_put_var(ncid(iLookFREQ%timestep),meta(iVar)%ncVarID(iLookFREQ%timestep),(/struct%var(iVar)/),start=(/iDOM,iSpatial/),count=(/1,1/))
    class is (var_dlength)
     if (iDOM==0) err = nf90_put_var(ncid(iLookFREQ%timestep),meta(iVar)%ncVarID(iLookFREQ%timestep),(/struct%var(iVar)%dat/),start=(/iSpatial,1/),count=(/1,size(struct%var(iVar)%dat)/))
     if (iDOM>0)  err = nf90_put_var(ncid(iLookFREQ%timestep),meta(iVar)%ncVarID(iLookFREQ%timestep),(/struct%var(iVar)%dat/),start=(/iDOM,iSpatial,1/),count=(/1,size(struct%var(iVar)%dat),1/))
   class default; err=20; message=trim(message)//'unknown variable type'; return
  end select
  call netcdf_err(err,message); if (err/=0) return

  ! re-initialize message
  message="writeParm/"
 end do  ! looping through local column model parameters

 end subroutine writeParm

 ! **************************************************************************************
 ! public subroutine writeData: write model time-dependent data
 ! **************************************************************************************
 subroutine writeData(finalizeStats,outputTimestep,maxDOM,nUNITrun,maxLayers,meta,stat,dat,map,indx,err,message)
 USE data_types,only:var_info                       ! metadata type
 USE var_lookup,only:maxVarStat                     ! index into stats structure
 USE var_lookup,only:iLookVarType                   ! index into type structure
 USE var_lookup,only:iLookINDEX                     ! index into index structure
 USE var_lookup,only:iLookSTAT                      ! index into stat structure
 USE globalData,only:outFreq,ncid                   ! output file information
 USE get_ixName_module,only:get_varTypeName         ! to access type strings for error messages
 USE get_ixName_module,only:get_statName            ! to access type strings for error messages
 implicit none
 ! declare dummy variables
 logical(lgt)  ,intent(in)            :: finalizeStats(:)  ! flags to finalize statistics
 integer(i4b)  ,intent(in)            :: outputTimestep(:) ! output time step
 integer(i4b)  ,intent(in)            :: maxDOM            ! maximum number of domains in any HRU
 integer(i4b)  ,intent(in)            :: nUNITrun          ! number of HRUs in the run space (for var)
 integer(i4b)  ,intent(in)            :: maxLayers         ! maximum number of layers
 type(var_info),intent(in)            :: meta(:)           ! meta data
 class(*)      ,intent(in)            :: stat              ! stats data
 class(*)      ,intent(in)            :: dat               ! timestep data
 integer(i4b)  ,intent(in)            :: map(:)            ! map into stats child struct
 type(gru_hru_dom_intVec) ,intent(in) :: indx              ! index data
 integer(i4b)  ,intent(out)           :: err               ! error code
 character(*)  ,intent(out)           :: message           ! error message
 ! local variables
 integer(i4b)                         :: iGRU              ! grouped response unit counter
 integer(i4b)                         :: iHRU              ! hydrologic response unit counter
 integer(i4b)                         :: iDOM              ! domain counter
 integer(i4b)                         :: iVar              ! variable index
 integer(i4b)                         :: iStat             ! statistics index
 integer(i4b)                         :: iFreq             ! frequency index
 integer(i4b)                         :: ncVarID           ! used only for time
 integer(i4b)                         :: nSnow             ! number of snow layers
 integer(i4b)                         :: nLake             ! number of lake layers
 integer(i4b)                         :: nSoil             ! number of soil layers
 integer(i4b)                         :: nIce              ! number of glacier ice layers
 integer(i4b)                         :: nLayers           ! total number of layers
 ! output arrays
 integer(i4b)                         :: datLength         ! length of each data vector
 integer(i4b)                         :: maxLength         ! maximum length of each data vector
 real(rkind)                          :: realVec(nUNITrun) ! real vector for all HRUs in the run space
 real(rkind)                          :: realVecDom(maxDOM,nUNITrun)  ! real array for all HRUs and DOMs in the run space
 real(rkind)                          :: realArray(maxDOM,nUNITrun,maxLayers+1)  ! real array for all HRUs and DOMs in the run space
 integer(i4b)                         :: intArray(maxDOM,nUNITrun,maxLayers+1)   ! integer array for all HRUs and DOMs in the run space
 integer(i4b)                         :: dataType          ! type of data
 integer(i4b),parameter               :: ixInteger=1001    ! named variable for integer
 integer(i4b),parameter               :: ixReal=1002       ! named variable for real
 ! initialize error control
 err=0;message="writeData/"

 ! loop through output frequencies
 do iFreq=1,maxvarFreq

  ! skip frequencies that are not needed
  if(.not.outFreq(iFreq)) cycle

  ! check that we have finalized statistics for a given frequency
  if(.not.finalizeStats(iFreq)) cycle

  ! loop through model variables
  do iVar = 1,size(meta)

   ! handle time first
   if (meta(iVar)%varName=='time')then
    ! get variable index
    err = nf90_inq_varid(ncid(iFreq),trim(meta(iVar)%varName),ncVarID)
    call netcdf_err(err,message); if (err/=0) return
    ! define HRUs and GRUs (only write once)
    iGRU=1; iHRU=1
    ! data bound write
    select type(dat) ! forcStruc
     class is (gru_hru_double)   ! x%gru(:)%hru(:)%var(:)
      err = nf90_put_var(ncid(iFreq),ncVarID,(/dat%gru(iGRU)%hru(iHRU)%var(iVar)/),start=(/outputTimestep(iFreq)/),count=(/1/))
      call netcdf_err(err,message); if (err/=0) return
      cycle ! move onto the next variable
     class default; err=20; message=trim(message)//'time variable must be of type gru_hru_double (forcing data structure)'; return
    end select
   end if  ! id time

   ! define the statistics index
   iStat = meta(iVar)%statIndex(iFreq)

   ! check that the variable is desired
   if (iStat==integerMissing.or.trim(meta(iVar)%varName)=='unknown') cycle

   ! stats output: only scalar variable type
   if(meta(iVar)%varType==iLookVarType%scalarv) then
    select type(stat)
     class is (gru_hru_doubleVec)

      ! loop through HRUs and GRUs, and place data in the single vector
      do iGRU=1,size(gru_struc)
       do iHRU=1,gru_struc(iGRU)%hruCount
         realVec(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix) = stat%gru(iGRU)%hru(iHRU)%var(map(iVar))%dat(iFreq)
       end do
      end do

      ! write data
      err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realVec(1:nUNITrun),start=(/1,outputTimestep(iFreq)/),count=(/nUNITrun,1/))

     class is (gru_hru_dom_doubleVec)
      realVecDom(:,:) = realMissing;    dataType=ixReal ! initialize the data array

      ! loop through DOMs, HRUs, and GRUs, and place data in the single vector
      do iGRU=1,size(gru_struc)
       do iHRU=1,gru_struc(iGRU)%hruCount
        do iDOM=1,gru_struc(iGRU)%hruInfo(iHRU)%domCount
         realVecDom(iDOM,gru_struc(iGRU)%hruInfo(iHRU)%hru_ix) = stat%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(map(iVar))%dat(iFreq)
        end do
       end do
      end do

      ! write data
      err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realVecDom(1:maxDOM,1:nUNITrun),start=(/1,1,outputTimestep(iFreq)/),count=(/maxDOM,nUNITrun,1/))

     class default; err=20; message=trim(message)//'stats must be scalarv and of type gru_hru_dom_doubleVec'; return
    end select  ! stat

   ! non-scalar variables: regular data structures
   else

    ! initialize the data vectors
    select type (dat)
     class is (gru_hru_dom_doubleVec); realArray(:,:,:) = realMissing;    dataType=ixReal
     class is (gru_hru_dom_intVec);     intArray(:,:,:) = integerMissing; dataType=ixInteger
     class default; err=20; message=trim(message)//'data must not be scalarv and either of type gru_hru_dom_doubleVec or gru_hru_dom_intVec'; return
    end select

    ! loop thru GRUs and HRUs
    do iGRU=1,size(gru_struc)
     do iHRU=1,gru_struc(iGRU)%hruCount
      do iDOM=1,gru_struc(iGRU)%hruInfo(iHRU)%domCount

       ! get the model layers
       nSnow   = indx%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nSnow)%dat(1)
       nLake   = indx%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nLake)%dat(1)
       nSoil   = indx%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nSoil)%dat(1)
       nIce    = indx%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nIce)%dat(1)
       nLayers = indx%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nLayers)%dat(1)

       ! get the length of each data vector
       select case (meta(iVar)%varType)
        case(iLookVarType%wLength); datLength = nSpecBand
        case(iLookVarType%midToto); datLength = nLayers
        case(iLookVarType%midSnow); datLength = nSnow
        case(iLookVarType%midLake); datLength = nLake
        case(iLookVarType%midSoil); datLength = nSoil
        case(iLookVarType%midIce ); datLength = nIce
        case(iLookVarType%ifcToto); datLength = nLayers+1
        case(iLookVarType%ifcSnow); datLength = nSnow+1
        case(iLookVarType%ifcLake); datLength = nLake+1
        case(iLookVarType%ifcSoil); datLength = nSoil+1
        case(iLookVarType%ifcIce ); datLength = nIce+1
        case default; cycle
       end select ! vartype
       
       ! get the data vectors
       select type (dat)
        class is (gru_hru_dom_doubleVec); realArray(iDOM,gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,1:datLength) = dat%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(:)
        class is (gru_hru_dom_intVec);     intArray(iDOM,gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,1:datLength) = dat%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(:)
       end select

      end do ! DOM loop
     end do ! HRU loop
    end do ! GRU loop

    ! get the maximum length of each data vector over all domains
    select case (meta(iVar)%varType)
     case(iLookVarType%wLength); maxLength = nSpecBand
     case(iLookVarType%midToto); maxLength = maxLayers
     case(iLookVarType%midSnow); maxLength = maxSnowLayers
     case(iLookVarType%midSoil); maxLength = maxSoilLayers
     case(iLookVarType%midLake); maxLength = maxLakeLayers
     case(iLookVarType%midIce ); maxLength = maxIceLayers
     case(iLookVarType%ifcToto); maxLength = maxLayers+1
     case(iLookVarType%ifcSnow); maxLength = maxSnowLayers+1
     case(iLookVarType%ifcSoil); maxLength = maxSoilLayers+1
     case(iLookVarType%ifcLake); maxLength = maxLakeLayers+1
     case(iLookVarType%ifcIce ); maxLength = maxIceLayers+1
     case default; cycle
    end select ! vartype

    ! write the data vectors
    select case(dataType)
     case(ixReal);    err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realArray(1:maxDOM,1:nUNITrun,1:maxLength),start=(/1,1,1,outputTimestep(iFreq)/),count=(/maxDOM,nUNITrun,maxLength,1/))
     case(ixInteger); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq), intArray(1:maxDOM,1:nUNITrun,1:maxLength),start=(/1,1,1,outputTimestep(iFreq)/),count=(/maxDOM,nUNITrun,maxLength,1/))
     case default; err=20; message=trim(message)//'data must be of type integer or real'; return
    end select ! data type

   end if ! not scalarv

   ! process error code
   if (err/=0) message=trim(message)//trim(meta(iVar)%varName)//'_'//trim(get_statName(iStat))
   call netcdf_err(err,message); if (err/=0) return

  end do ! iVar
 end do ! iFreq

 end subroutine writeData

 ! **************************************************************************************
 ! public subroutine writeBasin: write basin-average variables
 ! **************************************************************************************
 subroutine writeBasin(iGRU,finalizeStats,outputTimestep,meta,stat,dat,map,err,message)
 USE data_types,only:var_info                       ! metadata type
 USE var_lookup,only:maxVarStat                     ! index into stats structure
 USE var_lookup,only:iLookVarType                   ! index into type structure
 USE globalData,only:outFreq,ncid                   ! output file information
 USE get_ixName_module,only:get_varTypeName         ! to access type strings for error messages
 USE get_ixName_module,only:get_statName            ! to access type strings for error messages
 implicit none

 ! declare dummy variables
 integer(i4b)  ,intent(in)     :: iGRU              ! GRU index
 logical(lgt)  ,intent(in)     :: finalizeStats(:)  ! flags to finalize statistics
 integer(i4b)  ,intent(in)     :: outputTimestep(:) ! output time step
 type(var_info),intent(in)     :: meta(:)           ! meta data
 type(dlength) ,intent(in)     :: stat(:)           ! stats data
 type(dlength) ,intent(in)     :: dat(:)            ! timestep data
 integer(i4b)  ,intent(in)     :: map(:)            ! map into stats child struct
 integer(i4b)  ,intent(out)    :: err               ! error code
 character(*)  ,intent(out)    :: message           ! error message
 ! local variables
 integer(i4b)                  :: iVar              ! variable index
 integer(i4b)                  :: iStat             ! statistics index
 integer(i4b)                  :: iFreq             ! frequency index
 ! initialize error control
 err=0;message="f-writeBasin/"

 ! loop through output frequencies
 do iFreq=1,maxvarFreq

  ! skip frequencies that are not needed
  if(.not.outFreq(iFreq)) cycle

  ! check that we have finalized statistics for a given frequency
  if(.not.finalizeStats(iFreq)) cycle

  ! loop through model variables
  do iVar = 1,size(meta)

   ! define the statistics index
   iStat = meta(iVar)%statIndex(iFreq)

   ! check that the variable is desired
   if (iStat==integerMissing.or.trim(meta(iVar)%varName)=='unknown') cycle

   ! stats/data output - select data type
   select case (meta(iVar)%varType)

    case (iLookVarType%scalarv)
     err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),(/stat(map(iVar))%dat(iFreq)/),start=(/iGRU,outputTimestep(iFreq)/),count=(/1,1/))

    case (iLookVarType%routing)
     if (iFreq==1 .and. outputTimestep(iFreq)==1) then
      err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),(/dat(iVar)%dat/),start=(/1/),count=(/nTimeDelay/))
     end if

    case default
     err=40; message=trim(message)//"unknownVariableType[name='"//trim(meta(iVar)%varName)//"';type='"//trim(get_varTypeName(meta(iVar)%varType))//    "']"; return
   end select ! variable type

   ! process error code
   if (err.ne.0) message=trim(message)//trim(meta(iVar)%varName)//'_'//trim(get_statName(iStat))
   call netcdf_err(err,message); if (err/=0) return

  end do ! iVar
 end do ! iFreq

 end subroutine writeBasin

 ! **************************************************************************************
 ! public subroutine writeTime: write current time to all files
 ! **************************************************************************************
 subroutine writeTime(finalizeStats,outputTimestep,meta,dat,err,message)
 USE data_types,only:var_info                       ! metadata type
 USE globalData,only:ncid                           ! output file IDs
 USE var_lookup,only:iLookSTAT                      ! index into stat structure
 implicit none

 ! declare dummy variables
 logical(lgt)  ,intent(in)     :: finalizeStats(:)  ! flags to finalize statistics
 integer(i4b)  ,intent(in)     :: outputTimestep(:) ! output time step
 type(var_info),intent(in)     :: meta(:)           ! meta data
 integer       ,intent(in)     :: dat(:)            ! timestep data
 integer(i4b)  ,intent(out)    :: err               ! error code
 character(*)  ,intent(out)    :: message           ! error message
 ! local variables
 integer(i4b)                  :: iVar              ! variable index
 integer(i4b)                  :: iFreq             ! frequency index
 integer(i4b)                  :: ncVarID           ! used only for time
 ! initialize error control
 err=0;message="f-writeTime/"

 ! loop through output frequencies
 do iFreq=1,maxvarFreq

  ! check that we have finalized statistics for a given frequency
  if(.not.finalizeStats(iFreq)) cycle

  ! loop through model variables
  do iVar = 1,size(meta)

   ! check instantaneous
   if (meta(iVar)%statIndex(iFreq)/=iLookSTAT%inst) cycle

   ! get variable id in file

   err = nf90_inq_varid(ncid(iFreq),trim(meta(iVar)%varName),ncVarID)
   if (err/=0) message=trim(message)//trim(meta(iVar)%varName)
   call netcdf_err(err,message)
   if (err/=0) then; err=20; return; end if

   ! add to file
   err = nf90_put_var(ncid(iFreq),ncVarID,(/dat(iVar)/),start=(/outputTimestep(iFreq)/),count=(/1/))
   if (err/=0) message=trim(message)//trim(meta(iVar)%varName)
   print*, err, trim(meta(iVar)%varName)
   call netcdf_err(err,message)
   if (err/=0) then; err=20; return; end if

  end do ! iVar
 end do ! iFreq

 end subroutine writeTime

 ! *********************************************************************************************************
 ! public subroutine printRestartFile: print a re-start file
 ! *********************************************************************************************************
 subroutine writeRestart(filename,         & ! intent(in): name of restart file
                         nGRU,             & ! intent(in): number of global GRUs
                         nHRU,             & ! intent(in): number of global HRUs
                         maxDOM,           & ! intent(in): max number of domains in any HRU
                         prog_meta,        & ! intent(in): prognostics metadata
                         prog_data,        & ! intent(in): prognostics data
                         bvar_meta,        & ! intent(in): basin (gru) variable metadata
                         bvar_data,        & ! intent(in): basin (gru) variable data
                         maxLayers,        & ! intent(in): maximum number of layers
                         indx_meta,        & ! intent(in): index metadata
                         indx_data,        & ! intent(in): index data
                         err,message)        ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------
 ! access the derived types to define the data structures
 USE data_types,only:var_info               ! metadata
 ! access named variables defining elements in the data structures
 USE var_lookup,only:iLookINDEX             ! named variables for structure elements
 USE var_lookup,only:iLookVarType           ! named variables for structure elements
 USE var_lookup,only:iLookBVAR              ! named variables for structure elements
 ! constants
 USE globalData,only:gru_struc              ! gru-hru mapping structures
 ! external routines
 USE netcdf_util_module,only:nc_file_close  ! close netcdf file
 USE netcdf_util_module,only:nc_file_open   ! open netcdf file
 USE globalData,only:nTimeDelay             ! number of timesteps in the time delay histogram
 USE globalData,only:maxGlaciers            ! maximum number of glaciers
 USE def_output_module,only: write_hru_info ! write HRU information to netcdf file
 
 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! input
 character(len=256),intent(in)          :: filename      ! name of the restart file
 integer(i4b),intent(in)                :: nGRU          ! number of global GRUs
 integer(i4b),intent(in)                :: nHRU          ! number of global HRUs
 integer(i4b),intent(in)                :: maxDOM        ! max number of domains in any HRU
 type(var_info),intent(in)              :: prog_meta(:)  ! prognostic variable metadata
 type(gru_hru_dom_doubleVec),intent(in) :: prog_data     ! prognostic vars
 type(var_info),intent(in)              :: bvar_meta(:)  ! basin variable metadata
 type(gru_doubleVec),intent(in)         :: bvar_data     ! basin variables
 integer(i4b), intent(in)               :: maxLayers     ! maximum number of total layers
 type(var_info),intent(in)              :: indx_meta(:)  ! metadata
 type(gru_hru_dom_intVec),intent(in)    :: indx_data     ! indexing vars
 ! output: error control
 integer(i4b),intent(out)               :: err           ! error code
 character(*),intent(out)               :: message       ! error message
 ! --------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                       :: ncid          ! netcdf file id
 integer(i4b),allocatable           :: ncVarID(:)    ! netcdf variable id
 integer(i4b)                       :: ncSnowID      ! index variable id
 integer(i4b)                       :: ncSoilID      ! index variable id
 integer(i4b)                       :: ncIceID       ! index variable id
 integer(i4b)                       :: ncLakeID      ! index variable id
 integer(i4b)                       :: nSnow         ! number of snow layers
 integer(i4b)                       :: nLake         ! number of lake layers
 integer(i4b)                       :: nSoil         ! number of soil layers
 integer(i4b)                       :: nIce          ! number of glacier ice layers
 integer(i4b)                       :: nLayers       ! number of total layers
 integer(i4b),parameter             :: nScalar=1     ! size of a scalar
 integer(i4b)                       :: nProgVars     ! number of prognostic variables written to state file
 integer(i4b)                       :: hruDimID      ! variable dimension ID
 integer(i4b)                       :: gruDimID      ! variable dimension ID
 integer(i4b)                       :: domDimID      ! variable dimension ID
 integer(i4b)                       :: tdhDimID      ! variable dimension ID
 integer(i4b)                       :: nglDimID      ! variable dimension ID
 integer(i4b)                       :: scalDimID     ! variable dimension ID
 integer(i4b)                       :: specDimID     ! variable dimension ID
 integer(i4b)                       :: midTotoDimID  ! variable dimension ID
 integer(i4b)                       :: ifcTotoDimID  ! variable dimension ID
 integer(i4b)                       :: midSoilDimID  ! variable dimension ID
 integer(i4b)                       :: ifcSoilDimID  ! variable dimension ID
 integer(i4b)                       :: midSnowDimID  ! variable dimension ID
 integer(i4b)                       :: ifcSnowDimID  ! variable dimension ID
 integer(i4b)                       :: midIceDimID   ! variable dimension ID
 integer(i4b)                       :: ifcIceDimID   ! variable dimension ID
 integer(i4b)                       :: midLakeDimID  ! variable dimension ID
 integer(i4b)                       :: ifcLakeDimID  ! variable dimension ID
 character(len=32),parameter        :: hruDimName    ='hru'      ! dimension name for HRUs
 character(len=32),parameter        :: gruDimName    ='gru'      ! dimension name for GRUs
 character(len=32),parameter        :: domDimName    ='dom'      ! dimension name for DOMs
 character(len=32),parameter        :: tdhDimName    ='tdh'      ! dimension name for time-delay basin variables
 character(len=32),parameter        :: nglDimName    ='ngl'      ! dimension name for global variables
 character(len=32),parameter        :: scalDimName   ='scalarv'  ! dimension name for scalar data
 character(len=32),parameter        :: specDimName   ='spectral' ! dimension name for spectral bands
 character(len=32),parameter        :: midTotoDimName='midToto'  ! dimension name for layered varaiables
 character(len=32),parameter        :: ifcTotoDimName='ifcToto'  ! dimension name for layered variables
 character(len=32),parameter        :: midSoilDimName='midSoil'  ! dimension name for soil-only layers
 character(len=32),parameter        :: ifcSoilDimName='ifcSoil'  ! dimension name for soil-only layers
 character(len=32),parameter        :: midSnowDimName='midSnow'  ! dimension name for snow-only layers
 character(len=32),parameter        :: ifcSnowDimName='ifcSnow'  ! dimension name for snow-only layers
 character(len=32),parameter        :: midIceDimName ='midIce'   ! dimension name for glacier ice-only layers
 character(len=32),parameter        :: ifcIceDimName ='ifcIce'   ! dimension name for glacier ice-only layers
 character(len=32),parameter        :: midLakeDimName='midLake'  ! dimension name for lake-only layers
 character(len=32),parameter        :: ifcLakeDimName='ifcLake'  ! dimension name for lake-only layers
 integer(i4b)                       :: cHRU          ! count of HRUs
 integer(i4b)                       :: iDOM          ! index of DOMs
 integer(i4b)                       :: iHRU          ! index of HRUs
 integer(i4b)                       :: iGRU          ! index of GRUs
 integer(i4b)                       :: iVar          ! variable index
 integer(i4b)                       :: nGlacier      ! number of glaciers in GRU
 integer(i4b)                       :: nWetland      ! number of wetlands in GRU
 logical(lgt)                       :: okLength      ! flag to check if the vector length is OK
 character(len=256)                 :: cmessage      ! downstream error message
 ! --------------------------------------------------------------------------------------------------------

 ! initialize error control
 err=0; message='writeRestart/'

 ! size of prognostic variable vector
 nProgVars = size(prog_meta)
 ! include additional basin variable in ID array
 if (maxIceLayers > 0)then
   allocate(ncVarID(nProgVars+6))
 else
   allocate(ncVarID(nProgVars+1))
 end if

 ! create file
 err = nf90_create(trim(filename),NF90_NETCDF4,ncid)
 message='iCreate[create]'; call netcdf_err(err,message); if(err/=0)return

 ! define dimensions
                      err = nf90_def_dim(ncid,trim(gruDimName)    ,nGRU           ,    gruDimID); message='iCreate[gru]'     ; call netcdf_err(err,message); if(err/=0)return
                      err = nf90_def_dim(ncid,trim(hruDimName)    ,nHRU           ,    hruDimID); message='iCreate[hru]'     ; call netcdf_err(err,message); if(err/=0)return
                      err = nf90_def_dim(ncid,trim(domDimName)    ,maxDOM         ,    domDimID); message='iCreate[dom]'     ; call netcdf_err(err,message); if(err/=0)return
                      err = nf90_def_dim(ncid,trim(tdhDimName)    ,nTimeDelay     ,    tdhDimID); message='iCreate[tdh]'     ; call netcdf_err(err,message); if(err/=0)return
                      err = nf90_def_dim(ncid,trim(nglDimName)    ,maxGlaciers    ,    nglDimID); message='iCreate[ngl]'     ; call netcdf_err(err,message); if(err/=0)return
                      err = nf90_def_dim(ncid,trim(scalDimName)   ,nScalar        ,   scalDimID); message='iCreate[scalar]'  ; call netcdf_err(err,message); if(err/=0)return
                      err = nf90_def_dim(ncid,trim(specDimName)   ,nSpecBand      ,   specDimID); message='iCreate[spectral]'; call netcdf_err(err,message); if(err/=0)return
                      err = nf90_def_dim(ncid,trim(midTotoDimName),maxLayers      ,midTotoDimID); message='iCreate[midToto]' ; call netcdf_err(err,message); if(err/=0)return
                      err = nf90_def_dim(ncid,trim(ifcTotoDimName),maxLayers+1    ,ifcTotoDimID); message='iCreate[ifcToto]' ; call netcdf_err(err,message); if(err/=0)return
 if (maxSoilLayers>0) err = nf90_def_dim(ncid,trim(midSoilDimName),maxSoilLayers  ,midSoilDimID); message='iCreate[midSoil]' ; call netcdf_err(err,message); if(err/=0)return
 if (maxSoilLayers>0) err = nf90_def_dim(ncid,trim(ifcSoilDimName),maxSoilLayers+1,ifcSoilDimID); message='iCreate[ifcSoil]' ; call netcdf_err(err,message); if(err/=0)return
 if (maxSnowLayers>0) err = nf90_def_dim(ncid,trim(midSnowDimName),maxSnowLayers  ,midSnowDimID); message='iCreate[midSnow]' ; call netcdf_err(err,message); if(err/=0)return
 if (maxSnowLayers>0) err = nf90_def_dim(ncid,trim(ifcSnowDimName),maxSnowLayers+1,ifcSnowDimID); message='iCreate[ifcSnow]' ; call netcdf_err(err,message); if(err/=0)return
 if (maxIceLayers >0) err = nf90_def_dim(ncid,trim(midIceDimName) ,maxIceLayers   ,midIceDimID);  message='iCreate[midIce]'  ; call netcdf_err(err,message); if(err/=0)return
 if (maxIceLayers >0) err = nf90_def_dim(ncid,trim(ifcIceDimName) ,maxIceLayers+1 ,ifcIceDimID);  message='iCreate[ifcIce]'  ; call netcdf_err(err,message); if(err/=0)return
 if (maxLakeLayers>0) err = nf90_def_dim(ncid,trim(midLakeDimName),maxLakeLayers  ,midLakeDimID); message='iCreate[midLake]' ; call netcdf_err(err,message); if(err/=0)return
 if (maxLakeLayers>0) err = nf90_def_dim(ncid,trim(ifcLakeDimName),maxLakeLayers+1,ifcLakeDimID); message='iCreate[ifcLake]' ; call netcdf_err(err,message); if(err/=0)return
 ! re-initialize error control
 err=0; message='writeRestart/'

 ! define prognostic variables
 do iVar = 1,nProgVars
  if (prog_meta(iVar)%varType==iLookvarType%unknown) cycle

  ! define variable
  select case(prog_meta(iVar)%varType)
   case(iLookvarType%scalarv);                      err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/domDimID,hruDimID,  scalDimID /),ncVarID(iVar))
   case(iLookvarType%wLength);                      err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/domDimID,hruDimID,  specDimID /),ncVarID(iVar))
   case(iLookvarType%midToto);                      err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/domDimID,hruDimID,midTotoDimID/),ncVarID(iVar))
   case(iLookvarType%ifcToto);                      err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/domDimID,hruDimID,ifcTotoDimID/),ncVarID(iVar))
   case(iLookvarType%midSoil); if (maxSoilLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/domDimID,hruDimID,midSoilDimID/),ncVarID(iVar))
   case(iLookvarType%ifcSoil); if (maxSoilLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/domDimID,hruDimID,ifcSoilDimID/),ncVarID(iVar))
   case(iLookvarType%midSnow); if (maxSnowLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/domDimID,hruDimID,midSnowDimID/),ncVarID(iVar))
   case(iLookvarType%ifcSnow); if (maxSnowLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/domDimID,hruDimID,ifcSnowDimID/),ncVarID(iVar))
   case(iLookvarType%midIce);  if (maxIceLayers >0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/domDimID,hruDimID,midIceDimID /),ncVarID(iVar))
   case(iLookvarType%ifcIce);  if (maxIceLayers >0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/domDimID,hruDimID,ifcIceDimID /),ncVarID(iVar))
   case(iLookvarType%midLake); if (maxLakeLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/domDimID,hruDimID,midLakeDimID/),ncVarID(iVar))
   case(iLookvarType%ifcLake); if (maxLakeLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/domDimID,hruDimID,ifcLakeDimID/),ncVarID(iVar))
  end select

  ! check errors
  if(err/=0)then
   message=trim(message)//trim(cmessage)//' [variable '//trim(prog_meta(iVar)%varName)//']'
   return
  end if

  ! add parameter description
  err = nf90_put_att(ncid,ncVarID(iVar),'long_name',trim(prog_meta(iVar)%vardesc))
  call netcdf_err(err,message)

  ! add parameter units
  err = nf90_put_att(ncid,ncVarID(iVar),'units',trim(prog_meta(iVar)%varunit))
  call netcdf_err(err,message)

 end do ! iVar
 
 ! define selected basin variables (derived) -- e.g., hillslope routing, number of glaciers, area of glaciers, etc.
 err = nf90_def_var(ncid, trim(bvar_meta(iLookBVAR%routingRunoffFuture)%varName), nf90_double, (/gruDimID, tdhDimID /), ncVarID(nProgVars+1))
 err = nf90_put_att(ncid,ncVarID(nProgVars+1),'long_name',trim(bvar_meta(iLookBVAR%routingRunoffFuture)%vardesc));   call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncVarID(nProgVars+1),'units'    ,trim(bvar_meta(iLookBVAR%routingRunoffFuture)%varunit));   call netcdf_err(err,message)

 if (maxIceLayers > 0)then
   err = nf90_def_var(ncid, trim(bvar_meta(iLookBVAR%glacAblArea)%varName), nf90_double, (/gruDimID, nglDimID/), ncVarID(nProgVars+2))
   err = nf90_put_att(ncid,ncVarID(nProgVars+2),'long_name',trim(bvar_meta(iLookBVAR%glacAblArea)%vardesc));   call netcdf_err(err,message)
   err = nf90_put_att(ncid,ncVarID(nProgVars+2),'units'    ,trim(bvar_meta(iLookBVAR%glacAblArea)%varunit));   call netcdf_err(err,message)

   err = nf90_def_var(ncid, trim(bvar_meta(iLookBVAR%glacAccArea)%varName), nf90_double, (/gruDimID, nglDimID/), ncVarID(nProgVars+3))
   err = nf90_put_att(ncid,ncVarID(nProgVars+3),'long_name',trim(bvar_meta(iLookBVAR%glacAccArea)%vardesc));   call netcdf_err(err,message)
   err = nf90_put_att(ncid,ncVarID(nProgVars+3),'units'    ,trim(bvar_meta(iLookBVAR%glacAccArea)%varunit));   call netcdf_err(err,message)

   err = nf90_def_var(ncid, trim(bvar_meta(iLookBVAR%glacIceRunoffFuture)%varName), nf90_double, (/gruDimID, nglDimID/), ncVarID(nProgVars+4))
   err = nf90_put_att(ncid,ncVarID(nProgVars+4),'long_name',trim(bvar_meta(iLookBVAR%glacIceRunoffFuture)%vardesc));   call netcdf_err(err,message)
   err = nf90_put_att(ncid,ncVarID(nProgVars+4),'units'    ,trim(bvar_meta(iLookBVAR%glacIceRunoffFuture)%varunit));   call netcdf_err(err,message)

   err = nf90_def_var(ncid, trim(bvar_meta(iLookBVAR%glacSnowRunoffFuture)%varName), nf90_double, (/gruDimID, nglDimID/), ncVarID(nProgVars+5))
   err = nf90_put_att(ncid,ncVarID(nProgVars+5),'long_name',trim(bvar_meta(iLookBVAR%glacSnowRunoffFuture)%vardesc));   call netcdf_err(err,message)
   err = nf90_put_att(ncid,ncVarID(nProgVars+5),'units'    ,trim(bvar_meta(iLookBVAR%glacSnowRunoffFuture)%varunit));   call netcdf_err(err,message)  
   
   err = nf90_def_var(ncid, trim(bvar_meta(iLookBVAR%glacFirnRunoffFuture)%varName), nf90_double, (/gruDimID, nglDimID/), ncVarID(nProgVars+6))
   err = nf90_put_att(ncid,ncVarID(nProgVars+6),'long_name',trim(bvar_meta(iLookBVAR%glacFirnRunoffFuture)%vardesc));   call netcdf_err(err,message)
   err = nf90_put_att(ncid,ncVarID(nProgVars+6),'units'    ,trim(bvar_meta(iLookBVAR%glacFirnRunoffFuture)%varunit));   call netcdf_err(err,message)  
  endif
  
 ! define index variables - snow
 err = nf90_def_var(ncid,trim(indx_meta(iLookINDEX%nSnow)%varName),nf90_int,(/domDimID,hruDimID/),ncSnowID); call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncSnowID,'long_name',trim(indx_meta(iLookINDEX%nSnow)%vardesc));           call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncSnowID,'units'    ,trim(indx_meta(iLookINDEX%nSnow)%varunit));           call netcdf_err(err,message)

! define index variables - lake, 0 if no lake in HRU
 err = nf90_def_var(ncid,trim(indx_meta(iLookINDEX%nLake)%varName),nf90_int,(/domDimID,hruDimID/),ncLakeID); call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncLakeID,'long_name',trim(indx_meta(iLookINDEX%nLake)%vardesc));           call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncLakeID,'units'    ,trim(indx_meta(iLookINDEX%nLake)%varunit));           call netcdf_err(err,message)

 ! define index variables - soil
 err = nf90_def_var(ncid,trim(indx_meta(iLookINDEX%nSoil)%varName),nf90_int,(/domDimID,hruDimID/),ncSoilID); call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncSoilID,'long_name',trim(indx_meta(iLookINDEX%nSoil)%vardesc));           call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncSoilID,'units'    ,trim(indx_meta(iLookINDEX%nSoil)%varunit));           call netcdf_err(err,message)

 ! define index variables - ice, 0 if no glacier in HRU
 err = nf90_def_var(ncid,trim(indx_meta(iLookINDEX%nIce)%varName),nf90_int,(/domDimID,hruDimID/),ncIceID); call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncIceID,'long_name',trim(indx_meta(iLookINDEX%nIce)%vardesc));           call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncIceID,'units'    ,trim(indx_meta(iLookINDEX%nIce)%varunit));           call netcdf_err(err,message)

 ! end definition phase
 err = nf90_enddef(ncid); call netcdf_err(err,message); if (err/=0) return

 ! write variables
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
    cHRU = gru_struc(iGRU)%hruInfo(iHRU)%hru_ix
    do iDOM = 1,gru_struc(iGRU)%hruInfo(iHRU)%domCount
     do iVar = 1,size(prog_meta)

      ! excape if this variable is not used
      if (prog_meta(iVar)%varType==iLookvarType%unknown) cycle

      ! actual number of layers
      nSnow = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSnow
      nLake = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nLake
      nSoil = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSoil
      nIce  = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nIce
      nLayers = nSnow + nLake + nSoil + nIce

      ! check size
      ! NOTE: this may take time that we do not wish to use
      okLength=.true.
      select case (prog_meta(iVar)%varType)
       case(iLookVarType%scalarv);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat) == nScalar  )
       case(iLookVarType%wlength);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat) == nSpecBand)
       case(iLookVarType%midToto);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat) == nLayers  )
       case(iLookVarType%ifcToto);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat) == nLayers+1)
       case(iLookVarType%midSnow); if (nSnow>0) okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat) == nSnow    )
       case(iLookVarType%ifcSnow); if (nSnow>0) okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat) == nSnow+1  )
       case(iLookVarType%midLake); if (nLake>0) okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat) == nLake    )
       case(iLookVarType%ifcLake); if (nLake>0) okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat) == nLake+1  )
       case(iLookVarType%midSoil); if (nSoil>0) okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat) == nSoil    )
       case(iLookVarType%ifcSoil); if (nSoil>0) okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat) == nSoil+1  )
       case(iLookVarType%midIce);  if (nIce>0)  okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat) == nIce     )
       case(iLookVarType%ifcIce);  if (nIce>0)  okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat) == nIce+1   )
       case default; err=20; message=trim(message)//'unknown var type'; return
      end select

      ! error check
      if(.not.okLength)then
       message=trim(message)//'bad vector length for variable '//trim(prog_meta(iVar)%varname)
       err=20; return
      endif

      ! write data
      select case (prog_meta(iVar)%varType)
       case(iLookVarType%scalarv);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat/),start=(/iDOM,cHRU,1/),count=(/1,1,nScalar  /))
       case(iLookVarType%wlength);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat/),start=(/iDOM,cHRU,1/),count=(/1,1,nSpecBand/))
       case(iLookVarType%midToto);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat/),start=(/iDOM,cHRU,1/),count=(/1,1,nLayers  /))
       case(iLookVarType%ifcToto);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat/),start=(/iDOM,cHRU,1/),count=(/1,1,nLayers+1/))
       case(iLookVarType%midSnow); if (nSnow>0) err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat/),start=(/iDOM,cHRU,1/),count=(/1,1,nSnow    /))
       case(iLookVarType%ifcSnow); if (nSnow>0) err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat/),start=(/iDOM,cHRU,1/),count=(/1,1,nSnow+1  /))
       case(iLookVarType%midLake); if (nLake>0) err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat/),start=(/iDOM,cHRU,1/),count=(/1,1,nLake    /))
       case(iLookVarType%ifcLake); if (nLake>0) err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat/),start=(/iDOM,cHRU,1/),count=(/1,1,nLake+1  /))
       case(iLookVarType%midSoil); if (nSoil>0) err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat/),start=(/iDOM,cHRU,1/),count=(/1,1,nSoil    /))
       case(iLookVarType%ifcSoil); if (nSoil>0) err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat/),start=(/iDOM,cHRU,1/),count=(/1,1,nSoil+1  /))
       case(iLookVarType%midIce);  if (nIce>0)  err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat/),start=(/iDOM,cHRU,1/),count=(/1,1,nIce     /))
       case(iLookVarType%ifcIce);  if (nIce>0)  err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat/),start=(/iDOM,cHRU,1/),count=(/1,1,nIce+1   /))
       case default; err=20; message=trim(message)//'unknown var type'; return
      end select

      ! error check
      if (err.ne.0) message=trim(message)//'writing variable:'//trim(prog_meta(iVar)%varName)
      call netcdf_err(err,message); if (err/=0) return
      err=0; message='writeRestart/'

    end do ! iVar loop

    ! write index variables
    err=nf90_put_var(ncid,ncSnowID,(/indx_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nSnow)%dat/),start=(/iDOM,cHRU/),count=(/1,1/))
    err=nf90_put_var(ncid,ncLakeID,(/indx_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nLake)%dat/),start=(/iDOM,cHRU/),count=(/1,1/))
    err=nf90_put_var(ncid,ncSoilID,(/indx_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nSoil)%dat/),start=(/iDOM,cHRU/),count=(/1,1/))
    err=nf90_put_var(ncid,ncIceID, (/indx_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nIce)%dat/), start=(/iDOM,cHRU/),count=(/1,1/))
  
   end do ! iDOM loop
  end do ! iHRU loop
  
  ! write selected basin variables
  err=nf90_put_var(ncid,ncVarID(nProgVars+1),(/bvar_data%gru(iGRU)%var(iLookBVAR%routingRunoffFuture)%dat/), start=(/iGRU,1/),count=(/1,nTimeDelay/))
  if (maxIceLayers > 0)then
    nGlacier = gru_struc(iGRU)%nGlacier
    err=nf90_put_var(ncid,ncVarID(nProgVars+2),(/bvar_data%gru(iGRU)%var(iLookBVAR%glacAblArea)%dat/), start=(/iGRU,1/),count=(/1,nGlacier/))
    err=nf90_put_var(ncid,ncVarID(nProgVars+3),(/bvar_data%gru(iGRU)%var(iLookBVAR%glacAccArea)%dat/), start=(/iGRU,1/),count=(/1,nGlacier/))
    err=nf90_put_var(ncid,ncVarID(nProgVars+4),(/bvar_data%gru(iGRU)%var(iLookBVAR%glacIceRunoffFuture)%dat/),  start=(/iGRU,1/),count=(/1,nGlacier/))
    err=nf90_put_var(ncid,ncVarID(nProgVars+5),(/bvar_data%gru(iGRU)%var(iLookBVAR%glacSnowRunoffFuture)%dat/), start=(/iGRU,1/),count=(/1,nGlacier/))
    err=nf90_put_var(ncid,ncVarID(nProgVars+6),(/bvar_data%gru(iGRU)%var(iLookBVAR%glacFirnRunoffFuture)%dat/), start=(/iGRU,1/),count=(/1,nGlacier/))
  endif
  
 end do  ! iGRU loop

 ! write HRU dimension and ID for file
 call write_hru_info(ncid, err, cmessage); if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

 ! close file
 call nc_file_close(ncid,err,cmessage)
 if(err/=0)then;message=trim(message)//trim(cmessage);return;end if

 ! cleanup
 deallocate(ncVarID)

 end subroutine writeRestart

end module modelwrite_module
