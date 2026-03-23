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
USE nr_type

! missing values
USE globalData,only: integerMissing, realMissing
! output constraints
USE globalData,only:maxSnowLayers       ! maximum number of snow layers
USE globalData,only:maxSoilLayers       ! maximum number of soil layers
USE globalData,only:maxGlceLayers       ! maximum number of glacier ice layers
USE globalData,only:maxLakeLayers       ! maximum number of lake layers
USE globalData,only:maxGlaciers         ! maximum number of glaciers in a GRU
USE globalData,only:maxGrid             ! maximum number of grids in a GRU
USE globalData,only:maxGridX            ! maximum number of grid cells in the x-direction
USE globalData,only:maxGridY            ! maximum number of grid cells in the y-direction
USE globalData,only:nTimeDelay          ! number of timesteps in the time delay histogram
USE globalData,only:nSpecBand           ! maximum number of spectral bands
! provide access to global data
USE globalData,only:nGRUrun             ! number of GRUs in the run
USE globalData,only:nHRUrun             ! number of HRUs in the run
USE globalData,only:maxDOM              ! maximum number of domains in any HRU
USE globalData,only:maxLayers           ! maximum number of layers
USE globalData,only:gru_struc           ! gru->hru mapping structure
USE globalData,only:allowRoutingOutput  ! flag to allow routing variable output

! provide access to the derived types to define the data structures
USE data_types,only:&
                    ! final data vectors
                    dlength,               & ! var%dat
                    ilength,               & ! var%dat
                    ! no spatial dimension
                    var_i,                 & ! x%var(:)                          (i4b)
                    var_i8,                & ! x%var(:)                          (i8b)
                    var_d,                 & ! x%var(:)                          (rkind)
                    var_ilength,           & ! x%var(:)%dat                      (i4b)
                    var_dlength,           & ! x%var(:)%dat                      (rkind)
                    var_dlength2,          & ! x%var(:)%dat2(:,:)                (rkind)
                    ! gru dimension
                    gru_int,               & ! x%gru(:)%var(:)                   (i4b)
                    gru_int8,              & ! x%gru(:)%var(:)                   (i8b)                  
                    gru_double,            & ! x%gru(:)%var(:)                   (rkind)
                    gru_intVec,            & ! x%gru(:)%var(:)%dat               (i4b)
                    gru_doubleVec,         & ! x%gru(:)%var(:)%dat               (rkind)
                    ! gru+hru dimension
                    gru_hru_int,           & ! x%gru(:)%hru(:)%var(:)            (i4b)
                    gru_hru_int8,          & ! x%gru(:)%hru(:)%var(:)            (i8b)
                    gru_hru_double,        & ! x%gru(:)%hru(:)%var(:)            (rkind)
                    gru_hru_intVec,        & ! x%gru(:)%hru(:)%var(:)%dat        (i4b)
                    gru_hru_doubleVec,     & ! x%gru(:)%hru(:)%var(:)%dat        (rkind)
                    ! gru+hru+dom dimension
                    gru_hru_dom_int,       & ! x%gru(:)%hru(:)%dom(:)%var(:)     (i4b)
                    gru_hru_dom_int8,      & ! x%gru(:)%hru(:)%dom(:)%var(:)     (i8b)
                    gru_hru_dom_double,    & ! x%gru(:)%hru(:)%dom(:)%var(:)     (rkind)
                    gru_hru_dom_intVec,    & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat (i4b)
                    gru_hru_dom_doubleVec, & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat (rkind)
                     ! gru+grid dimension
                    gru_grid_double          ! x%gru(:)%grid(:)%var(:)%dat2(:,:) (rkind)
          
! vector lengths
USE var_lookup, only: maxvarFreq ! number of output frequencies
USE var_lookup, only: maxvarStat ! number of statistics

implicit none
private
public::writeParam
public::writeGridParam
public::writeData
public::writeGridData
public::writeTime
public::writeRestart

contains

 ! **********************************************************************************************************
 ! public subroutine writeParam: write model parameters
 ! **********************************************************************************************************
 subroutine writeParam(iDOM,iSpatial,struct,meta,err,message)
 USE globalData,only:ncid                        ! netcdf file ids
 USE data_types,only:var_info                    ! metadata info
 USE var_lookup,only:iLookSTAT                   ! index in statistics vector
 USE var_lookup,only:iLookFREQ                   ! index in vector of model output frequencies
 implicit none

 ! declare input variables
 integer(i4b)  ,intent(in)   :: iDOM             ! domain index (0 if no domain dimension)
 integer(i4b)  ,intent(in)   :: iSpatial         ! HRU index or GRU index
 class(*)      ,intent(in)   :: struct           ! data structure
 type(var_info),intent(in)   :: meta(:)          ! metadata structure
 integer(i4b)  ,intent(out)  :: err              ! error code
 character(*)  ,intent(out)  :: message          ! error message
 ! local variables
 integer(i4b)                :: iVar             ! loop through variables

 ! initialize error control
 err=0;message="writeParam/"

 ! loop through local column model parameters
 do iVar = 1,size(meta)

  ! check that the variable is desired
  if (meta(iVar)%statIndex(iLookFREQ%timestep)==integerMissing) cycle

  ! initialize message
  message=trim(message)//trim(meta(iVar)%varName)//':'

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
   class default; err=20; message=trim(message)//'parameter type must be var_i, var_i8, var_d, or var_dlength'; return
  end select
  call netcdf_err(err,message); if (err/=0) return

  ! re-initialize message
  message="writeParam/"
 end do  ! looping through local column model parameters

 end subroutine writeParam

 ! **********************************************************************************************************
 ! public subroutine writeGridParam: write grid model parameters
 ! **********************************************************************************************************
 subroutine writeGridParam(iSpatial,iGrid,struct,meta,err,message)
 USE globalData,only:ncid                        ! netcdf file ids
 USE data_types,only:var_info                    ! metadata info
 USE var_lookup,only:iLookSTAT                   ! index in statistics vector
 USE var_lookup,only:iLookFREQ                   ! index in vector of model output frequencies
 implicit none

 ! declare input variables
 integer(i4b)  ,intent(in)   :: iSpatial         ! HRU index or GRU index
 integer(i4b)  ,intent(in)   :: iGrid            ! grid index
 class(*)      ,intent(in)   :: struct           ! data structure
 type(var_info),intent(in)   :: meta(:)          ! metadata structure
 integer(i4b)  ,intent(out)  :: err              ! error code
 character(*)  ,intent(out)  :: message          ! error message
 ! local variables
 integer(i4b)                :: nx,ny            ! grid dimensions
 integer(i4b)                :: iVar             ! loop through variables

 ! initialize error control
 err=0;message="writeGridParam/"

 ! loop through local column model parameters
 do iVar = 1,size(meta)

  ! check that the variable is desired
  if (meta(iVar)%statIndex(iLookFREQ%annual)==integerMissing) cycle

  ! only write parameters that are not ids (currently only discludes surface_elev, debris_thick, and cell2hru, but could be others in the future)
  if(trim(meta(iVar)%varName)/='surface_elev' .and. trim(meta(iVar)%varName)/='debris_thick' .and. trim(meta(iVar)%varName)/='cell2hru') then

   ! initialize message
   message=trim(message)//trim(meta(iVar)%varName)//':'

   select type (struct)
    class is (var_dlength2)
     nx = size(struct%var(iVar)%dat2, 1)
     ny = size(struct%var(iVar)%dat2, 2)
     err = nf90_put_var(ncid(iLookFREQ%annual),meta(iVar)%ncVarID(iLookFREQ%annual),(/struct%var(iVar)%dat2/),start=(/iSpatial,iGrid,1,1/),count=(/1,1,nx,ny/))
    class default; err=20; message=trim(message)//'parameter type must be var_dlength2'; return
   end select
   call netcdf_err(err,message); if (err/=0) return

  end if ! if parameter

  ! re-initialize message
  message="writeGridParam/"
 end do  ! looping through local column model parameters

 end subroutine writeGridParam

 ! **************************************************************************************
 ! public subroutine writeData: write model time-dependent data for each HRU
 ! **************************************************************************************
 subroutine writeData(is_bufferedWrite,finalizeStats,outputTimestep,maxWrite,meta,stat,datt,map,indx,err,message)
 USE data_types,only:var_info                       ! metadata type
 USE var_lookup,only:maxvarStat                     ! index into stats structure
 USE var_lookup,only:iLookVarType                   ! index into type structure
 USE var_lookup,only:iLookINDEX                     ! index into index structure
 USE var_lookup,only:iLookFREQ                      ! index into freq structure
 USE globalData,only:outFreq,ncid                   ! output file information
 USE get_ixName_module,only:get_varTypeName         ! to access type strings for error messages
 USE get_ixName_module,only:get_statName            ! to access type strings for error messages
 implicit none
 ! declare dummy variables
 logical(lgt)  ,intent(in)          :: is_bufferedWrite                     ! flag for buffered write
 logical(lgt)  ,intent(in)          :: finalizeStats(:)                     ! flags to finalize statistics
 integer(i4b)  ,intent(in)          :: outputTimestep(:)                    ! output time step
 integer(i4b)  ,intent(in)          :: maxWrite                             ! maximum number of steps written
 type(var_info),intent(in)          :: meta(:)                              ! meta data
 class(*)      ,intent(in)          :: stat                                 ! stats data
 class(*)      ,intent(in)          :: datt(:)                              ! timestep or buffer data
 integer(i4b)  ,intent(in)          :: map(:)                               ! map into stats child struct
 type(gru_hru_dom_intVec),intent(in):: indx                                 ! index data
 integer(i4b)  ,intent(out)         :: err                                  ! error code
 character(*)  ,intent(out)         :: message                              ! error message
 ! local variables
 integer(i4b)                       :: iGRU                                 ! grouped response unit counter
 integer(i4b)                       :: iHRU                                 ! hydrologic response unit counter
 integer(i4b)                       :: iDOM                                 ! domain counter
 integer(i4b)                       :: iVar                                 ! variable index
 integer(i4b)                       :: iStat                                ! statistics index
 integer(i4b)                       :: iFreq                                ! frequency index
 integer(i4b)                       :: iTime                                ! time index
 integer(i4b)                       :: ncVarID                              ! used only for time
 integer(i4b)                       :: nSnow                                ! number of snow layers
 integer(i4b)                       :: nLake                                ! number of lake layers
 integer(i4b)                       :: nSoil                                ! number of soil layers
 integer(i4b)                       :: nGlce                                ! number of glacier ice layers
 integer(i4b)                       :: nGlac                                ! number of glaciers in the GRU
 integer(i4b)                       :: nLayers                              ! total number of layers
 integer(i4b)                       :: ixStart                              ! index of the start of data write
 integer(i4b)                       :: nSpace                               ! number of spatial data elements
 ! output arrays
 integer(i4b)                       :: datLength                            ! length of each data vector
 integer(i4b)                       :: maxLength                            ! maximum length of each data vector
 real(rkind)                        :: timeBuffer(maxWrite)                 ! buffer for all time steps
 real(rkind)                        :: realBuffer(nHRUrun,maxWrite)         ! buffer for all HRUs in the run domain + time steps
 real(rkind)                        :: realBuffer3(maxDOM,nHRUrun,maxWrite) ! buffer for all HRUs and DOMs in the run domain + time steps
 real(rkind),allocatable            :: realArray(:,:)                       ! real array for all HRUs in the run domain
 integer(i4b),allocatable           :: intArray(:,:)                        ! integer array for all HRUs in the run domain
 real(rkind),allocatable            :: realArray3(:,:,:)                    ! real array for all HRUs and DOMs in the run domain
 integer(i4b),allocatable           :: intArray3(:,:,:)                     ! integer array for all HRUs and DOMs in the run domain
 integer(i4b)                       :: dataType                             ! type of data
 integer(i4b),parameter             :: ixInteger=1001                       ! named variable for integer
 integer(i4b),parameter             :: ixReal=1002                          ! named variable for real
 integer(i4b),parameter             :: ixInteger3=1003                      ! named variable for integer array with 3 dimensions (e.g. dom, hru, time)
 integer(i4b),parameter             :: ixReal3=1004                         ! named variable for real array with 3 dimensions (e.g. dom, hru, time)

 ! initialize error control
 err=0;message="writeData/"

 ! allocate real and integer arrays for non-scalar variables to longest possible length
 maxLength = max(nSpecBand,maxLayers+1)
 maxLength = max(maxLength,maxGlaciers)
 if(allowRoutingOutput) maxLength = max(maxLength, nTimeDelay)
 allocate(realArray(nHRUrun,maxLength),intArray(nHRUrun,maxLength),realArray3(maxDOM,nHRUrun,maxLength),intArray3(maxDOM,nHRUrun,maxLength))

 ! loop through output frequencies
 do iFreq=1,maxvarFreq

  ! skip frequencies that are not needed
  if(.not.outFreq(iFreq)) cycle

  ! restrict attention to the timestep data if buffered write
  if(is_bufferedWrite .and. iFreq/=iLookFREQ%timestep) cycle

  ! check that we have finalized statistics for a given frequency
  if(.not.finalizeStats(iFreq)) cycle

  ! get the start index
  if(is_bufferedWrite)then
   ixStart = 1
  else
   ixStart = outputTimestep(iFreq)
  endif

  ! loop through model variables
  do iVar = 1,size(meta)
    
   ! initialize message
   message=trim(message)//trim(meta(iVar)%varName)

   ! ****************************************************************************
   ! *** write time information -- instantaneous
   ! ****************************************************************************

   ! handle time first
   if(meta(iVar)%varName=='time')then
    message=trim(message)//':' ! add statistic (none) to message 

    ! get variable index
    err = nf90_inq_varid(ncid(iFreq),trim(meta(iVar)%varName),ncVarID)
    call netcdf_err(err,message); if (err/=0) return

    ! define HRUs and GRUs (only write once)
    iGRU=1; iHRU=1

    ! data bound array access
    select type(datt) ! forcStruc
     class is (gru_hru_double) ! x%gru(:)%hru(:)%var(:)

      ! put data in time buffer
      do iTime=1,maxWrite
       timeBuffer(iTime) = datt(iTime)%gru(iGRU)%hru(iHRU)%var(iVar)
      end do

     ! check we found the class
     class default; err=20; message=trim(message)//'time variable must be of type gru_hru_double (forcing data structure)'; return
    end select  ! type of data structure

    ! write time
    err = nf90_put_var(ncid(iFreq),ncVarID,(/timeBuffer/),start=(/ixStart/),count=(/maxWrite/))
    call netcdf_err(err,message); if (err/=0) return
    message="writeData/" ! re-initialize message
    cycle ! move onto the next variable

   end if  ! if time

   ! ****************************************************************************
   ! *** write scalar variables
   ! ****************************************************************************

   ! define the statistics index
   iStat = meta(iVar)%statIndex(iFreq)
   message=trim(message)//'_'//trim(get_statName(iStat))//':' ! add statistic to message

   ! check that the variable is desired, currently do not write large variables (unknown and routing) as they are large and slow things down a lot
   if (iStat==integerMissing .or. meta(iVar)%varType==iLookVarType%unknown .or. meta(iVar)%varType==integerMissing)then; message="writeData/"; cycle; endif 
   if (meta(iVar)%varType==iLookVarType%routing .and. .not.allowRoutingOutput)then; message="writeData/"; cycle; endif ! routing variable write can be turned on with the allowRoutingOutput flag

   ! stats output: only scalar variable type
   if(meta(iVar)%varType==iLookVarType%scalarv) then

    ! ----- writing buffered output data ---------------------------------------
    if(is_bufferedWrite)then

     ! loop through time, DOMS, HRUs and GRUs, and place data in the buffer
     do iTime=1,maxWrite
      do iGRU=1,size(gru_struc)

       ! identify data structures
       select type(datt)

        ! *** HRU DOM structures (indices, ...)
        class is (gru_hru_dom_int)
         do iHRU=1,gru_struc(iGRU)%hruCount
          do iDOM=1,gru_struc(iGRU)%hruInfo(iHRU)%domCount
           realBuffer3(iDOM,gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,iTime) = datt(iTime)%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(map(iVar))
          end do  ! dom
         end do  ! hru

        class is (gru_hru_dom_int8)
         do iHRU=1,gru_struc(iGRU)%hruCount
          do iDOM=1,gru_struc(iGRU)%hruInfo(iHRU)%domCount
           realBuffer3(iDOM,gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,iTime) = datt(iTime)%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(map(iVar))
          end do  ! dom
         end do  ! hru
  
        ! *** HRU DOM structures (prognostic, diagnostic, ...)
        class is (gru_hru_dom_double)
        do iHRU=1,gru_struc(iGRU)%hruCount
         do iDOM=1,gru_struc(iGRU)%hruInfo(iHRU)%domCount
          realBuffer3(iDOM,gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,iTime) = datt(iTime)%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(map(iVar))
         end do  ! dom
        end do  ! hru

        class is (gru_hru_int)
         do iHRU=1,gru_struc(iGRU)%hruCount
          realBuffer(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,iTime) = datt(iTime)%gru(iGRU)%hru(iHRU)%var(map(iVar))
         end do  ! hru
          
        class is (gru_hru_int8)
         do iHRU=1,gru_struc(iGRU)%hruCount
          realBuffer(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,iTime) = datt(iTime)%gru(iGRU)%hru(iHRU)%var(map(iVar))
         end do  ! hru

        ! *** HRU structures (forcing, ...)
        class is (gru_hru_double)
         do iHRU=1,gru_struc(iGRU)%hruCount
          realBuffer(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,iTime) = datt(iTime)%gru(iGRU)%hru(iHRU)%var(map(iVar))
         end do  ! hru

        class is (gru_int)
         realBuffer(iGRU,iTime) = datt(iTime)%gru(iGRU)%var(map(iVar))
 
        class is (gru_int8)
         realBuffer(iGRU,iTime) = datt(iTime)%gru(iGRU)%var(map(iVar))

        ! *** GRU structures (basin-average variables, ...)
        class is (gru_double)
         realBuffer(iGRU,iTime) = datt(iTime)%gru(iGRU)%var(map(iVar))

        class default; err=20; message=trim(message)//'scalarv variables must be of type gru_hru_dom_[double or int*], gru_hru_[double or int*], or gru_[double or int*]'; return
       end select  ! time step data structure

      end do  ! gru
     end do  ! time
                  
     ! write data
     select type(datt)
      class is (gru_hru_dom_int);    err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realBuffer3(1:maxDOM,1:nHRUrun,1:maxWrite),start=(/1,1,1/),count=(/maxDOM,nHRUrun,maxWrite/))
      class is (gru_hru_dom_int8);   err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realBuffer3(1:maxDOM,1:nHRUrun,1:maxWrite),start=(/1,1,1/),count=(/maxDOM,nHRUrun,maxWrite/))
      class is (gru_hru_dom_double); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realBuffer3(1:maxDOM,1:nHRUrun,1:maxWrite),start=(/1,1,1/),count=(/maxDOM,nHRUrun,maxWrite/))
      class is (gru_hru_int);    err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realBuffer(1:nHRUrun,1:maxWrite),start=(/1,1/),count=(/nHRUrun,maxWrite/))
      class is (gru_hru_int8);   err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realBuffer(1:nHRUrun,1:maxWrite),start=(/1,1/),count=(/nHRUrun,maxWrite/))
      class is (gru_hru_double); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realBuffer(1:nHRUrun,1:maxWrite),start=(/1,1/),count=(/nHRUrun,maxWrite/))
      class is (gru_int);        err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realBuffer(1:nGRUrun,1:maxWrite),start=(/1,1/),count=(/nGRUrun,maxWrite/))
      class is (gru_int8);       err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realBuffer(1:nGRUrun,1:maxWrite),start=(/1,1/),count=(/nGRUrun,maxWrite/))
      class is (gru_double);     err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realBuffer(1:nGRUrun,1:maxWrite),start=(/1,1/),count=(/nGRUrun,maxWrite/))
     end select
     call netcdf_err(err,message); if (err/=0) return

    ! ----- writing statistics -------------------------------------------------

    ! check that we are not writing buffered data -- writing statistics
    else

     ! check that maxWrite==1
     if(maxWrite/=1)then
      message=trim(message)//'expect maxWrite=1 when not writing buffered output'
      err=20; return
     endif

     ! loop through DOMS, HRUs and GRUs, and place data in the buffer
     do iGRU=1,size(gru_struc)
     
      ! identify data structures
      select type(stat)

       ! *** HRU DOM structures (prognostic, diagnostic, ...)
       class is (gru_hru_dom_doubleVec)
        do iHRU=1,gru_struc(iGRU)%hruCount
         do iDOM=1,gru_struc(iGRU)%hruInfo(iHRU)%domCount
          realBuffer3(iDOM,gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,1) = stat%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(map(iVar))%dat(iFreq)
         end do ! dom
        end do ! hru

       ! *** HRU structures (forcing, prognostic, diagnostic, ...)
       class is (gru_hru_doubleVec)
        do iHRU=1,gru_struc(iGRU)%hruCount
         realBuffer(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,1) = stat%gru(iGRU)%hru(iHRU)%var(map(iVar))%dat(iFreq)
        end do
       
       ! *** GRU structures (basin-average variables, ...)
       class is (gru_doubleVec)
        realBuffer(iGRU,1) = stat%gru(iGRU)%var(map(iVar))%dat(iFreq)

       ! check statistics type
       class default; message=trim(message)//'stats must be scalarv and of type gru_hru_dom_doubleVec, gru_hru_doubleVec, or gru_doubleVec'; err=20; return
      end select  ! stat data structure

     end do  ! gru

     ! write data 
     select type(stat)
      class is (gru_hru_dom_doubleVec); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realBuffer3(1:maxDOM,1:nHRUrun,1),start=(/1,1,outputTimestep(iFreq)/),count=(/maxDOM,nHRUrun,1/))
      class is (gru_hru_doubleVec); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realBuffer(1:nHRUrun,1),start=(/1,outputTimestep(iFreq)/),count=(/nHRUrun,1/))
      class is (gru_doubleVec);     err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realBuffer(1:nGRUrun,1),start=(/1,outputTimestep(iFreq)/),count=(/nGRUrun,1/))
     end select
     call netcdf_err(err,message); if (err/=0) return
     
    endif  ! (if not buffered write -- statistics)

   ! ****************************************************************************
   ! *** write non-scalar variables (regular data structures -- instantaneous)
   ! ****************************************************************************

   ! non-scalar variables: regular data structures
   else

    ! cannot write non-scalar variables in buffered write -- too complicated and slow, so not currently supported
    if(is_bufferedWrite)then
      write(*,*)'WARNING: cannot output non-scalar type data when using the buffered write option (writeFullSeries), skipping variable '//trim(meta(iVar)%varName)
      message="writeData/"; cycle
    endif 

    ! initialize the data vectors
    select type (datt)
     class is (gru_hru_dom_doubleVec); nSpace = nHRUrun; realArray3(:,:,:) = realMissing;   dataType=ixReal3
     class is (gru_hru_dom_intVec);    nSpace = nHRUrun; intArray3(:,:,:) = integerMissing; dataType=ixInteger3
     class is (gru_hru_doubleVec); nSpace = nHRUrun; realArray(:,:) = realMissing;   dataType=ixReal
     class is (gru_hru_intVec);    nSpace = nHRUrun; intArray(:,:) = integerMissing; dataType=ixInteger
     class is (gru_doubleVec);     nSpace = nGRUrun; realArray(:,:) = realMissing;   dataType=ixReal
     class is (gru_intVec);        nSpace = nGRUrun; intArray(:,:) = integerMissing; dataType=ixInteger
     class default; err=20; message=trim(message)//'data is not scalarv so should be either of type gru_hru_dom_[double or int]Vec, gru_hru_[double or int]Vec, or gru_[double or int]Vec'; return
    end select

    ! loop thru GRUs and HRUs
    gruLoop: do iGRU=1,size(gru_struc)
     hruLoop: do iHRU = 1, gru_struc(iGRU)%hruCount
      domLoop: do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount

      ! get the model layers
      nSnow   = indx%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nSnow)%dat(1)
      nLake   = indx%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nLake)%dat(1)
      nSoil   = indx%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nSoil)%dat(1)
      nGlce   = indx%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nGlce)%dat(1)
      nLayers = indx%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nLayers)%dat(1)
      nGlac   = gru_struc(iGRU)%nGlac

      ! get the length of each data vector
      select case (meta(iVar)%varType)
       case(iLookVarType%wLength); datLength = nSpecBand
       case(iLookVarType%midToto); datLength = nLayers
       case(iLookVarType%midSnow); datLength = nSnow
       case(iLookVarType%midLake); datLength = nLake
       case(iLookVarType%midSoil); datLength = nSoil
       case(iLookVarType%midGlce); datLength = nGlce
       case(iLookVarType%ifcToto); datLength = nLayers+1
       case(iLookVarType%ifcSnow); datLength = nSnow+1
       case(iLookVarType%ifcLake); datLength = nLake+1
       case(iLookVarType%ifcSoil); datLength = nSoil+1
       case(iLookVarType%ifcGlce); datLength = nGlce+1
       case(iLookVarType%routing); datLength = nTimeDelay
       case(iLookVarType%glacier); datLength = nGlac
       case default; cycle
      ! case parSoil only in parameters (mpar, not written here) 
      ! case unknown skipped above
      end select ! varType
      
      ! get the data vectors
      select type (datt)
       class is (gru_hru_dom_doubleVec); realArray3(iDOM,gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,1:datLength) = datt(1)%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(:)
       class is (gru_hru_dom_intVec);     intArray3(iDOM,gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,1:datLength) = datt(1)%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(:)
       class is (gru_hru_doubleVec); realArray(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,1:datLength) = datt(1)%gru(iGRU)%hru(iHRU)%var(iVar)%dat(:); if(iDOM==1) exit domLoop ! only need to get the HRU-level data once
       class is (gru_hru_intVec);     intArray(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,1:datLength) = datt(1)%gru(iGRU)%hru(iHRU)%var(iVar)%dat(:); if(iDOM==1) exit domLoop ! only need to get the HRU-level data once
       class is (gru_doubleVec); realArray(iGRU,1:datLength) = datt(1)%gru(iGRU)%var(iVar)%dat(:); if(iHRU==1) exit hruLoop ! only need to get the GRU-level data once
       class is (gru_intVec);     intArray(iGRU,1:datLength) = datt(1)%gru(iGRU)%var(iVar)%dat(:); if(iHRU==1) exit hruLoop ! only need to get the GRU-level data once
      end select

      end do domLoop ! DOM loop
     end do hruLoop ! HRU loop
    end do gruLoop ! GRU loop

    ! get the maximum length of each data vector over all domains
    select case (meta(iVar)%varType)
     case(iLookVarType%wLength); maxLength = nSpecBand
     case(iLookVarType%midToto); maxLength = maxLayers
     case(iLookVarType%midSnow); maxLength = maxSnowLayers
     case(iLookVarType%midSoil); maxLength = maxSoilLayers
     case(iLookVarType%midLake); maxLength = maxLakeLayers
     case(iLookVarType%midGlce); maxLength = maxGlceLayers
     case(iLookVarType%ifcToto); maxLength = maxLayers+1
     case(iLookVarType%ifcSnow); maxLength = maxSnowLayers+1
     case(iLookVarType%ifcSoil); maxLength = maxSoilLayers+1
     case(iLookVarType%ifcLake); maxLength = maxLakeLayers+1
     case(iLookVarType%ifcGlce); maxLength = maxGlceLayers+1
     case(iLookVarType%routing); maxLength = nTimeDelay
     case(iLookVarType%glacier); maxLength = maxGlaciers
     case default; cycle
    end select ! varType

    ! write the data vectors
    if(maxLength==0) cycle ! skip if there is no length
    select case(dataType)
     case(ixReal3);    err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realArray3(1:maxDOM,1:nHRUrun,1:maxLength),start=(/1,1,1,outputTimestep(iFreq)/),count=(/maxDOM,nHRUrun,maxLength,1/))
     case(ixInteger3); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq), intArray3(1:maxDOM,1:nHRUrun,1:maxLength),start=(/1,1,1,outputTimestep(iFreq)/),count=(/maxDOM,nHRUrun,maxLength,1/))
     case(ixReal);    err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realArray(1:nSpace,1:maxLength),start=(/1,1,outputTimestep(iFreq)/),count=(/nSpace,maxLength,1/))
     case(ixInteger); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq), intArray(1:nSpace,1:maxLength),start=(/1,1,outputTimestep(iFreq)/),count=(/nSpace,maxLength,1/))
     case default; err=20; message=trim(message)//'data must be of type integer or real'; return
    end select ! data type

   end if ! not scalarv

   ! process error code
   call netcdf_err(err,message); if (err/=0) return
   message="writeData/" ! re-initialize message

  end do ! iVar
 end do ! iFreq
 deallocate(realArray,intArray,realArray3,intArray3)

 end subroutine writeData

 ! **************************************************************************************
 ! public subroutine writeGridData: write model grid data for each HRU, only writes in non-buffer mode
 ! **************************************************************************************
 subroutine writeGridData(finalizeStats,outputTimestep,meta,datt,err,message)
 USE data_types,only:var_info                       ! metadata type
 USE var_lookup,only:maxVarStat                     ! index into stats structure
 USE var_lookup,only:iLookFREQ                      ! index into freq structure
 USE globalData,only:outFreq,ncid                   ! output file information
 USE get_ixName_module,only:get_varTypeName         ! to access type strings for error messages
 USE get_ixName_module,only:get_statName            ! to access type strings for error messages
 implicit none
 ! declare dummy variables
 logical(lgt)  ,intent(in)   :: finalizeStats(:)       ! flags to finalize statistics
 integer(i4b)  ,intent(in)   :: outputTimestep(:)      ! output time step
 type(var_info),intent(in)   :: meta(:)                ! meta data
 class(*)      ,intent(in)   :: datt                   ! timestep data
 integer(i4b)  ,intent(out)  :: err                    ! error code
 character(*)  ,intent(out)  :: message                ! error message
 ! local variables
 integer(i4b)                :: iGRU                   ! grouped response unit counter
 integer(i4b)                :: iGrid                  ! grid counter
 integer(i4b)                :: iVar                   ! variable index
 integer(i4b)                :: iStat                  ! statistics index
 integer(i4b)                :: iFreq                  ! frequency index
 integer(i4b)                :: nGrid                  ! number of grids in the GRU
 integer(i4b)                :: nx,ny                  ! number of grid cells in x,y directions
 integer(i4b)                :: ixStart                ! index of the start of data write
 ! output arrays
 real(rkind)                 :: realArray4(nGRUrun,maxGrid,maxGridX,maxGridY) ! real array for all GRUs and grids in the run domain
 integer(i4b)                :: dataType               ! type of data
 integer(i4b),parameter      :: ixReal4=1001           ! named variable for real

 ! initialize error control
 err=0;message="writeGridData/"

 ! loop through output frequencies
 ! NOTE: grid variables will only be at annual frequency for now, but leaving in in case this changes in the future
 do iFreq=1,maxvarFreq

  ! skip frequencies that are not needed
  if(.not.outFreq(iFreq)) cycle

  ! check that we have finalized statistics for a given frequency
  if(.not.finalizeStats(iFreq)) cycle

  ! get the start index
  ixStart = outputTimestep(iFreq)

  ! loop through model variables
  do iVar = 1,size(meta)
    
   ! initialize message
   message=trim(message)//trim(meta(iVar)%varName)

   ! ****************************************************************************
   ! *** write grid variables (regular data structures -- instantaneous)
   ! ****************************************************************************

   ! only write variables that change with time (currently only surface_elev and debris_thick, but could be others in the future)
   if(meta(iVar)%varName=='surface_elev' .or. meta(iVar)%varName=='debris_thick') then

    ! define the statistics index
    iStat = meta(iVar)%statIndex(iFreq)
    message=trim(message)//'_'//trim(get_statName(iStat))//':' ! add statistic to message

    ! check that the variable is desired
    if (iStat==integerMissing)then; message="writeGridData/"; cycle; endif 

    ! initialize the data vectors
    select type (datt)
     class is (gru_grid_double); realArray4(:,:,:,:) = realMissing; dataType=ixReal4
     class default; err=20; message=trim(message)//'grid structure data should be of type gru_grid_double'; return
    end select

    ! loop thru GRUs and grids
    do iGRU=1,size(gru_struc)
     nGrid = gru_struc(iGRU)%nGrid
     do iGrid = 1,nGrid ! for now all grids are glaciers
      ! get the length of each data vector
      nx = gru_struc(iGRU)%gridInfo(iGrid)%nx
      ny = gru_struc(iGRU)%gridInfo(iGrid)%ny

      ! get the data vectors
      select type (datt)
       class is (gru_grid_double); realArray4(iGRU,iGrid,1:nx,1:ny) = datt%gru(iGRU)%grid(iGrid)%var(iVar)%dat2(:,:)
      end select

     end do  ! glac loop
    end do  ! GRU loop

    ! write the data vectors
    if(maxGrid==0) cycle ! skip if there is no length
    select case(dataType)
     case(ixReal4); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iFreq),realArray4(1:nGRUrun,1:maxGrid,1:maxGridX,1:maxGridY),start=(/1,1,1,1,outputTimestep(iFreq)/),count=(/nGRUrun,maxGrid,maxGridX,maxGridY,1/))
     case default; err=20; message=trim(message)//'data must be of type real'; return
    end select ! data type

   end if ! not changing with time

   ! process error code
   call netcdf_err(err,message); if (err/=0) return
   message="writeGridData/" ! re-initialize message

  end do ! iVar
 end do ! iFreq

 end subroutine writeGridData

 ! **************************************************************************************
 ! public subroutine writeTime: write current time to all files
 ! **************************************************************************************
 subroutine writeTime(finalizeStats,outputTimestep,meta,datt,err,message)
 USE data_types,only:var_info                       ! metadata type
 USE globalData,only:ncid                           ! output file IDs
 USE var_lookup,only:iLookSTAT                      ! index into stat structure
 implicit none

 ! declare dummy variables
 logical(lgt)  ,intent(in)     :: finalizeStats(:)  ! flags to finalize statistics
 integer(i4b)  ,intent(in)     :: outputTimestep(:) ! output time step
 type(var_info),intent(in)     :: meta(:)           ! meta data
 integer       ,intent(in)     :: datt(:)           ! timestep data
 integer(i4b)  ,intent(out)    :: err               ! error code
 character(*)  ,intent(out)    :: message           ! error message
 ! local variables
 integer(i4b)                  :: iVar              ! variable index
 integer(i4b)                  :: iFreq             ! frequency index
 integer(i4b)                  :: ncVarID           ! used only for time
 ! initialize error control
 err=0;message="writeTime/"

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
   err = nf90_put_var(ncid(iFreq),ncVarID,(/datt(iVar)/),start=(/outputTimestep(iFreq)/),count=(/1/))
   if (err/=0) message=trim(message)//trim(meta(iVar)%varName)
   call netcdf_err(err,message)
   if (err/=0) then; err=20; return; end if

  end do ! iVar
 end do ! iFreq

 end subroutine writeTime

 ! *********************************************************************************************************
 ! public subroutine writeRestart: print a re-start file
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
                         grid_meta,        & ! intent(in): grid metadata
                         grid_data,        & ! intent(in): grid data
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
 USE def_output_module,only: write_id_info ! write HRU information to netcdf file
 
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
 type(var_info),intent(in)              :: grid_meta(:)  ! grid metadata
 type(gru_grid_double),intent(in)       :: grid_data     ! grid data
 ! output: error control
 integer(i4b),intent(out)               :: err           ! error code
 character(*),intent(out)               :: message       ! error message
 ! --------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                       :: ncid          ! netcdf file id
 integer(i4b),allocatable           :: ncVarID(:)    ! netcdf variable id
 integer(i4b),dimension(7)          :: ngdx          ! intermediate array of loop indices for glacier variables
 integer(i4b),dimension(4)          :: nidx          ! intermediate array of loop indices for index variables
 integer(i4b)                       :: nSnow         ! number of snow layers
 integer(i4b)                       :: nLake         ! number of lake layers
 integer(i4b)                       :: nSoil         ! number of soil layers
 integer(i4b)                       :: nGlce         ! number of glacier ice layers
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
 integer(i4b)                       :: midGlceDimID  ! variable dimension ID
 integer(i4b)                       :: ifcGlceDimID  ! variable dimension ID
 integer(i4b)                       :: midLakeDimID  ! variable dimension ID
 integer(i4b)                       :: ifcLakeDimID  ! variable dimension ID
 character(len=32),parameter        :: hruDimName    ='hru'      ! dimension name for HRUs
 character(len=32),parameter        :: gruDimName    ='gru'      ! dimension name for GRUs
 character(len=32),parameter        :: domDimName    ='dom'      ! dimension name for DOMs
 character(len=32),parameter        :: tdhDimName    ='tdh'      ! dimension name for time-delay basin variables
 character(len=32),parameter        :: nglDimName    ='glac'     ! dimension name for glacier
 character(len=32),parameter        :: scalDimName   ='scalarv'  ! dimension name for scalar data
 character(len=32),parameter        :: specDimName   ='spectral' ! dimension name for spectral bands
 character(len=32),parameter        :: midTotoDimName='midToto'  ! dimension name for layered varaiables
 character(len=32),parameter        :: ifcTotoDimName='ifcToto'  ! dimension name for layered variables
 character(len=32),parameter        :: midSoilDimName='midSoil'  ! dimension name for soil-only layers
 character(len=32),parameter        :: ifcSoilDimName='ifcSoil'  ! dimension name for soil-only layers
 character(len=32),parameter        :: midSnowDimName='midSnow'  ! dimension name for snow-only layers
 character(len=32),parameter        :: ifcSnowDimName='ifcSnow'  ! dimension name for snow-only layers
 character(len=32),parameter        :: midGlceDimName='midGlce'  ! dimension name for glacier ice-only layers
 character(len=32),parameter        :: ifcGlceDimName='ifcGlce'  ! dimension name for glacier ice-only layers
 character(len=32),parameter        :: midLakeDimName='midLake'  ! dimension name for lake-only layers
 character(len=32),parameter        :: ifcLakeDimName='ifcLake'  ! dimension name for lake-only layers
 integer(i4b)                       :: cHRU          ! count of HRUs
 integer(i4b)                       :: iDOM          ! index of DOMs
 integer(i4b)                       :: iHRU          ! index of HRUs
 integer(i4b)                       :: iGRU          ! index of GRUs
 integer(i4b)                       :: i             ! loop index
 integer(i4b)                       :: iVar          ! variable index
 integer(i4b)                       :: nGlac         ! number of glaciers in GRU
 logical(lgt)                       :: okLength      ! flag to check if the vector length is OK
 integer(i4b)                       :: size_prog     ! size of prognostic variable vector without index variables
 character(len=256)                 :: cmessage      ! downstream error message
 ! --------------------------------------------------------------------------------------------------------

 ! initialize error control
 err=0; message='writeRestart/'

 ! size of prognostic variable vector
 nProgVars = size(prog_meta)

 ! index variables
 nidx = (/iLookINDEX%nSnow, iLookINDEX%nLake, iLookINDEX%nSoil, iLookINDEX%nGlce/)

 ! include additional basin variable in ID array
 size_prog = nProgVars+1 ! +1 for future runoff variable
 if (maxGlaciers > 0)then
   ngdx = (/iLookBVAR%basin__GlacierStorage,iLookBVAR%updateJulDay,iLookBVAR%glacierAblArea,iLookBVAR%glacierAccArea,iLookBVAR%glacIceRunoffFuture,iLookBVAR%glacSnowRunoffFuture,iLookBVAR%glacFirnRunoffFuture/)
   size_prog =  size_prog+size(ngdx)
 endif
 allocate(ncVarID(size_prog+size(nidx)))

 ! create file
 err = nf90_create(trim(filename),NF90_NETCDF4,ncid)
 message='iCreate[create]'; call netcdf_err(err,message); if(err/=0)return

 ! define dimensions
                     err = nf90_def_dim(ncid,trim(gruDimName)    ,nGRU           ,    gruDimID); message='iCreate[gru]'     ; call netcdf_err(err,message); if(err/=0)return
                     err = nf90_def_dim(ncid,trim(hruDimName)    ,nHRU           ,    hruDimID); message='iCreate[hru]'     ; call netcdf_err(err,message); if(err/=0)return
                     err = nf90_def_dim(ncid,trim(domDimName)    ,maxDOM         ,    domDimID); message='iCreate[dom]'     ; call netcdf_err(err,message); if(err/=0)return
                     err = nf90_def_dim(ncid,trim(tdhDimName)    ,nTimeDelay     ,    tdhDimID); message='iCreate[tdh]'     ; call netcdf_err(err,message); if(err/=0)return
 if( maxGlaciers>0)  err = nf90_def_dim(ncid,trim(nglDimName)    ,maxGlaciers    ,    nglDimID); message='iCreate[glac]'    ; call netcdf_err(err,message); if(err/=0)return
                     err = nf90_def_dim(ncid,trim(scalDimName)   ,nScalar        ,   scalDimID); message='iCreate[scalar]'  ; call netcdf_err(err,message); if(err/=0)return
                     err = nf90_def_dim(ncid,trim(specDimName)   ,nSpecBand      ,   specDimID); message='iCreate[spectral]'; call netcdf_err(err,message); if(err/=0)return
                     err = nf90_def_dim(ncid,trim(midTotoDimName),maxLayers      ,midTotoDimID); message='iCreate[midToto]' ; call netcdf_err(err,message); if(err/=0)return
                     err = nf90_def_dim(ncid,trim(ifcTotoDimName),maxLayers+1    ,ifcTotoDimID); message='iCreate[ifcToto]' ; call netcdf_err(err,message); if(err/=0)return
 if(maxSoilLayers>0) err = nf90_def_dim(ncid,trim(midSoilDimName),maxSoilLayers  ,midSoilDimID); message='iCreate[midSoil]' ; call netcdf_err(err,message); if(err/=0)return
 if(maxSoilLayers>0) err = nf90_def_dim(ncid,trim(ifcSoilDimName),maxSoilLayers+1,ifcSoilDimID); message='iCreate[ifcSoil]' ; call netcdf_err(err,message); if(err/=0)return
 if(maxSnowLayers>0) err = nf90_def_dim(ncid,trim(midSnowDimName),maxSnowLayers  ,midSnowDimID); message='iCreate[midSnow]' ; call netcdf_err(err,message); if(err/=0)return
 if(maxSnowLayers>0) err = nf90_def_dim(ncid,trim(ifcSnowDimName),maxSnowLayers+1,ifcSnowDimID); message='iCreate[ifcSnow]' ; call netcdf_err(err,message); if(err/=0)return
 if(maxGlceLayers>0) err = nf90_def_dim(ncid,trim(midGlceDimName),maxGlceLayers  ,midGlceDimID); message='iCreate[midGlce]' ; call netcdf_err(err,message); if(err/=0)return
 if(maxGlceLayers>0) err = nf90_def_dim(ncid,trim(ifcGlceDimName),maxGlceLayers+1,ifcGlceDimID); message='iCreate[ifcGlce]' ; call netcdf_err(err,message); if(err/=0)return
 if(maxLakeLayers>0) err = nf90_def_dim(ncid,trim(midLakeDimName),maxLakeLayers  ,midLakeDimID); message='iCreate[midLake]' ; call netcdf_err(err,message); if(err/=0)return
 if(maxLakeLayers>0) err = nf90_def_dim(ncid,trim(ifcLakeDimName),maxLakeLayers+1,ifcLakeDimID); message='iCreate[ifcLake]' ; call netcdf_err(err,message); if(err/=0)return
 err=0; message='writeRestart/' ! reset message

 ! define prognostic variables
 do iVar = 1,nProgVars
  ! define variable
  select case(prog_meta(iVar)%varType)
   case(iLookVarType%scalarv);                      err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/domDimID,hruDimID,  scalDimID /),ncVarID(iVar))
   case(iLookVarType%wLength);                      err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/domDimID,hruDimID,  specDimID /),ncVarID(iVar))
   case(iLookVarType%midToto);                      err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/domDimID,hruDimID,midTotoDimID/),ncVarID(iVar))
   case(iLookVarType%ifcToto);                      err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/domDimID,hruDimID,ifcTotoDimID/),ncVarID(iVar))
   case(iLookVarType%midSoil); if (maxSoilLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/domDimID,hruDimID,midSoilDimID/),ncVarID(iVar))
   case(iLookVarType%ifcSoil); if (maxSoilLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/domDimID,hruDimID,ifcSoilDimID/),ncVarID(iVar))
   case(iLookVarType%midSnow); if (maxSnowLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/domDimID,hruDimID,midSnowDimID/),ncVarID(iVar))
   case(iLookVarType%ifcSnow); if (maxSnowLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/domDimID,hruDimID,ifcSnowDimID/),ncVarID(iVar))
   case(iLookVarType%midGlce); if (maxGlceLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/domDimID,hruDimID,midGlceDimID/),ncVarID(iVar))
   case(iLookVarType%ifcGlce); if (maxGlceLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/domDimID,hruDimID,ifcGlceDimID/),ncVarID(iVar))
   case(iLookVarType%midLake); if (maxLakeLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/domDimID,hruDimID,midLakeDimID/),ncVarID(iVar))
   case(iLookVarType%ifcLake); if (maxLakeLayers>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varName),nf90_double,(/domDimID,hruDimID,ifcLakeDimID/),ncVarID(iVar))
   case(iLookVarType%unknown); cycle
  end select
  if(err/=0)then; message=trim(message)//' [variable '//trim(prog_meta(iVar)%varName)//']'; return; end if

  ! add parameter description and units
  err = nf90_put_att(ncid,ncVarID(iVar),'long_name',trim(prog_meta(iVar)%vardesc)); call netcdf_err(err,message)
  err = nf90_put_att(ncid,ncVarID(iVar),'units',trim(prog_meta(iVar)%varunit)); call netcdf_err(err,message)
 end do ! iVar
 
 ! define selected basin variables (derived) -- e.g., hillslope routing, number of glaciers, area of glaciers, etc.
 err = nf90_def_var(ncid, trim(bvar_meta(iLookBVAR%routingRunoffFuture)%varName), nf90_double, (/gruDimID, tdhDimID /), ncVarID(nProgVars+1))
 if(err/=0)then; message=trim(message)//' [variable '//trim(bvar_meta(iLookBVAR%routingRunoffFuture)%varName)//']'; return; end if
 err = nf90_put_att(ncid,ncVarID(nProgVars+1),'long_name',trim(bvar_meta(iLookBVAR%routingRunoffFuture)%vardesc));   call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncVarID(nProgVars+1),'units'    ,trim(bvar_meta(iLookBVAR%routingRunoffFuture)%varunit));   call netcdf_err(err,message)

 if(maxGlaciers > 0)then ! if glaciers are present, include glacier variables
   do i = 1,size(ngdx)
    iVar = ngdx(i)
    select case(bvar_meta(iVar)%varType)
     case(iLookVarType%scalarv); err = nf90_def_var(ncid,trim(bvar_meta(iVar)%varName),nf90_double,(/gruDimID,scalDimID /),ncVarID(nProgVars+1+i))
     case(iLookVarType%glacier); err = nf90_def_var(ncid,trim(bvar_meta(iVar)%varName),nf90_double,(/gruDimID,nglDimID/),ncVarID(nProgVars+1+i))
    end select
    if(err/=0)then; message=trim(message)//' [variable '//trim(bvar_meta(iVar)%varName)//']';return; end if

    ! add parameter description and units
    err = nf90_put_att(ncid,ncVarID(nProgVars+1+i),'long_name',trim(bvar_meta(iVar)%vardesc)); call netcdf_err(err,message)
    err = nf90_put_att(ncid,ncVarID(nProgVars+1+i),'units',trim(bvar_meta(iVar)%varunit)); call netcdf_err(err,message)
   end do ! iVar
 endif ! (if glaciers)
  
 ! define index variables
 do i=1,size(nidx)
   iVar = nidx(i)
   err = nf90_def_var(ncid,trim(indx_meta(iVar)%varName),nf90_int,(/domDimID,hruDimID/),ncVarID(size_prog+i))
   if(err/=0)then; message=trim(message)//' [variable '//trim(indx_meta(iVar)%varName)//']'; return; end if
   err = nf90_put_att(ncid,ncVarID(size_prog+i),'long_name',trim(indx_meta(iVar)%vardesc));   call netcdf_err(err,message)
   err = nf90_put_att(ncid,ncVarID(size_prog+i),'units'    ,trim(indx_meta(iVar)%varunit));   call netcdf_err(err,message)
 end do

 ! end definition phase
 err = nf90_enddef(ncid); call netcdf_err(err,message); if (err/=0) return

 ! write variables
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
    cHRU = gru_struc(iGRU)%hruInfo(iHRU)%hru_ix
    do iDOM = 1,gru_struc(iGRU)%hruInfo(iHRU)%domCount
     do iVar = 1,size(prog_meta)

      ! escape if this variable is not used
      if (prog_meta(iVar)%varType==iLookVarType%unknown) cycle

      ! actual number of layers
      nSnow = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSnow
      nLake = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nLake
      nSoil = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSoil
      nGlce = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nGlce
      nLayers = nSnow + nLake + nSoil + nGlce

      ! check size, NOTE: this may take time that we do not wish to use
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
       case(iLookVarType%midGlce); if (nGlce>0) okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat) == nGlce    )
       case(iLookVarType%ifcGlce); if (nGlce>0) okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat) == nGlce+1  )
       case default; err=20; message=trim(message)//'unknown var type'; return
      end select
      if(.not.okLength)then; message=trim(message)//'bad vector length for variable '//trim(prog_meta(iVar)%varName); err=20; return; endif

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
       case(iLookVarType%midGlce); if (nGlce>0) err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat/),start=(/iDOM,cHRU,1/),count=(/1,1,nGlce    /))
       case(iLookVarType%ifcGlce); if (nGlce>0) err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat/),start=(/iDOM,cHRU,1/),count=(/1,1,nGlce+1  /))
       case default; err=20; message=trim(message)//'unknown var type'; return
      end select
      if (err/=0) message=trim(message)//'writing variable:'//trim(prog_meta(iVar)%varName); call netcdf_err(err,message); if (err/=0) return; err=0; message='writeRestart/'

    end do ! iVar loop

    ! write index variables
    do i=1,size(nidx)
      iVar = nidx(i)
      err=nf90_put_var(ncid,ncVarID(size_prog+i),(/indx_data%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat/),start=(/iDOM,cHRU/),count=(/1,1/))
      if (err/=0) message=trim(message)//'writing variable:'//trim(indx_meta(iVar)%varName); call netcdf_err(err,message); if (err/=0) return; err=0; message='writeRestart/'
    end do
  
   end do ! iDOM loop
  end do ! iHRU loop
  
  ! write selected basin variables
  err=nf90_put_var(ncid,ncVarID(nProgVars+1),(/bvar_data%gru(iGRU)%var(iLookBVAR%routingRunoffFuture)%dat/), start=(/iGRU,1/),count=(/1,nTimeDelay/))
  if (err/=0) message=trim(message)//'writing variable:'//trim(bvar_meta(iLookBVAR%routingRunoffFuture)%varName); call netcdf_err(err,message); if (err/=0) return; err=0; message='writeRestart/'

  if (maxGlaciers > 0)then ! if glaciers are present, include glacier variables
    nGlac = gru_struc(iGRU)%nGlac
    do i=1,size(ngdx)
      iVar = ngdx(i)
      select case(bvar_meta(iVar)%varType)
       case(iLookVarType%scalarv); err=nf90_put_var(ncid,ncVarID(nProgVars+1+i),(/bvar_data%gru(iGRU)%var(iVar)%dat/), start=(/iGRU,1/),count=(/1,nScalar/))
       case(iLookVarType%glacier); err=nf90_put_var(ncid,ncVarID(nProgVars+1+i),(/bvar_data%gru(iGRU)%var(iVar)%dat/), start=(/iGRU,1/),count=(/1,nGlac/))
       case default; err=20; message=trim(message)//'unknown var type'; return
      end select
      if (err/=0) message=trim(message)//'writing variable:'//trim(bvar_meta(iVar)%varName); call netcdf_err(err,message); if (err/=0) return; err=0; message='writeRestart/'
    end do

    ! include grids
    call writeRestartGrid(ncid, nGRU, gruDimID, grid_meta, grid_data, err, cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif
  
 end do  ! iGRU loop

 ! write dimensions and ID for file
 call write_id_info(ncid, gruDimID, hruDimID, domDimID, nglDimID, err, cmessage); if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

 ! close file
 call nc_file_close(ncid,err,cmessage)
 if(err/=0)then;message=trim(message)//trim(cmessage);return;end if

 ! cleanup
 deallocate(ncVarID)

 end subroutine writeRestart

 ! *********************************************************************************************************
 ! private subroutine writeRestartGrid: print a re-start file for grids
 ! *********************************************************************************************************
 subroutine writeRestartGrid(ncid,             & ! intent(in): netcdf file id
                             nGRU,             & ! intent(in): number of global GRUs
                             gruDimID,         & ! intent(in): dimension ID for GRUs
                             grid_meta,        & ! intent(in): grid metadata
                             grid_data,        & ! intent(in): grid data
                             err,message)        ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------
 ! access the derived types to define the data structures
 USE data_types,only:var_info                ! metadata
 ! access named variables defining elements in the data structures
 USE var_lookup,only:iLookGRID               ! named variables for structure elements
 ! constants
 USE globalData,only:gru_struc               ! gru-hru mapping structures
 ! external routines
 USE def_output_module,only: write_gridid_info ! write glacier grid information to netcdf file
 
 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! input
 integer(i4b),intent(in)                :: ncid          ! netcdf file id
 integer(i4b),intent(in)                :: nGRU          ! number of global GRUs
 integer(i4b),intent(in)                :: gruDimID      ! variable dimension ID
 type(var_info),intent(in)              :: grid_meta(:)  ! grid metadata
 type(gru_grid_double),intent(in)       :: grid_data     ! grid data
 ! output: error control
 integer(i4b),intent(out)               :: err           ! error code
 character(*),intent(out)               :: message       ! error message
 ! --------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b),dimension(2)          :: ncVarID            ! netcdf variable id, only one variable currently
 integer(i4b),dimension(2)          :: ngdx               ! intermediate array of loop indices
 integer(i4b)                       :: gridDimID          ! variable dimension ID
 integer(i4b)                       :: xDimID             ! variable dimension ID
 integer(i4b)                       :: yDimID             ! variable dimension ID
 character(len=32),parameter        :: gridDimName='grid' ! dimension name for grid variables
 character(len=32),parameter        :: xDimName='xgrid'   ! dimension name for xgrid
 character(len=32),parameter        :: yDimName='ygrid'   ! dimension name for ygrid
 integer(i4b)                       :: iGRU               ! index of GRUs
 integer(i4b)                       :: iGrid              ! index of grids
 integer(i4b)                       :: i                  ! loop index
 integer(i4b)                       :: nx,ny              ! grid dimensions
 integer(i4b)                       :: iVar               ! variable index
 integer(i4b)                       :: nGrid              ! number of grids in GRU
 character(len=256)                 :: cmessage           ! downstream error message

 ! --------------------------------------------------------------------------------------------------------

 ! initialize error control
 err=0; message='writeRestartGrid/'

 ! grid variables
 ngdx = (/iLookGRID%surface_elev, iLookGRID%debris_thick/) ! array of desired variable indices

 ! define dimensions
 err = nf90_def_dim(ncid,trim(gridDimName)   ,maxGrid     , gridDimID);message='iCreate[grid]'    ; call netcdf_err(err,message); if(err/=0)return
 err = nf90_def_dim(ncid,trim(xDimName)      ,maxGridX    , xDimID);   message='iCreate[xgrid]'   ; call netcdf_err(err,message); if(err/=0)return
 err = nf90_def_dim(ncid,trim(yDimName)      ,maxGridY    , yDimID);   message='iCreate[ygrid]'   ; call netcdf_err(err,message); if(err/=0)return
 ! re-initialize error control
 err=0; message='writeRestartGrid/'

 ! define grid variables
 do i = 1,size(ngdx)
  iVar = ngdx(i)
  err = nf90_def_var(ncid, trim(grid_meta(iVar)%varName), nf90_double, (/gruDimID, gridDimID, xDimID, yDimID/), ncVarID(i))
 
  ! check errors
  if(err/=0)then
    message=trim(message)//' [variable '//trim(grid_meta(iVar)%varName)//']'
    return
  end if

  ! add parameter description
  err = nf90_put_att(ncid,ncVarID(i),'long_name',trim(grid_meta(iVar)%vardesc))
  call netcdf_err(err,message)

  ! add parameter units
  err = nf90_put_att(ncid,ncVarID(i),'units',trim(grid_meta(iVar)%varunit))
  call netcdf_err(err,message)

end do ! iVar

! end definition phase
err = nf90_enddef(ncid); call netcdf_err(err,message); if (err/=0) return

! write variables
do iGRU = 1,nGRU
  nGrid = gru_struc(iGRU)%nGrid
  do iGrid = 1,nGrid
    do i = 1,size(ngdx)
      iVar = ngdx(i)
      nx = gru_struc(iGRU)%gridInfo(iGrid)%nx
      ny = gru_struc(iGRU)%gridInfo(iGrid)%ny
      err=nf90_put_var(ncid,ncVarID(i),(/grid_data%gru(iGRU)%grid(iGrid)%var(iVar)%dat2/), start=(/iGRU,iGrid,1,1/),count=(/1,1,nx,ny/))
      ! error check
      if (err/=0) message=trim(message)//'writing variable:'//trim(grid_meta(iVar)%varName)
      call netcdf_err(err,message); if (err/=0) return
      message = 'writeRestartGrid/'
    end do
  end do

end do  ! iGRU loop

! write grid dimension and ID for file
call write_gridid_info(ncid, gruDimID, gridDimID, err, cmessage); if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

end subroutine writeRestartGrid

end module modelwrite_module
