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

module summa_writeOutput
! used to define/write output files

 ! named variables to define new output files
USE globalData, only: noNewFiles              ! no new output files
USE globalData, only: newFileEveryOct1        ! create a new file on Oct 1 every year (start of the USA water year)

!  model decisions
USE globalData,only:model_decisions           ! model decision structure

! metadata
USE globalData,only:time_meta                 ! metadata on the model time
USE globalData,only:forc_meta                 ! metadata on the model forcing data
USE globalData,only:diag_meta                 ! metadata on the model diagnostic variables
USE globalData,only:prog_meta                 ! metadata on the model prognostic variables
USE globalData,only:flux_meta                 ! metadata on the model fluxes
USE globalData,only:indx_meta                 ! metadata on the model index variables
USE globalData,only:bvar_meta                 ! metadata on basin-average variables
USE globalData,only:grid_meta                 ! metadata on the glacier grid information

! child metadata for stats
USE globalData,only:statForc_meta             ! child metadata for stats
USE globalData,only:statProg_meta             ! child metadata for stats
USE globalData,only:statDiag_meta             ! child metadata for stats
USE globalData,only:statFlux_meta             ! child metadata for stats
USE globalData,only:statIndx_meta             ! child metadata for stats
USE globalData,only:statBvar_meta             ! child metadata for stats

! index of the child data structure
USE globalData,only:forcChild_map             ! index of the child data structure: stats forc
USE globalData,only:progChild_map             ! index of the child data structure: stats prog
USE globalData,only:diagChild_map             ! index of the child data structure: stats diag
USE globalData,only:fluxChild_map             ! index of the child data structure: stats flux
USE globalData,only:indxChild_map             ! index of the child data structure: stats indx
USE globalData,only:bvarChild_map             ! index of the child data structure: stats bvar

! named variables
USE var_lookup,only:maxvarFreq                ! maximum number of output files
USE var_lookup,only:iLookTIME                 ! named variables for time data structure
USE var_lookup,only:iLookDIAG                 ! named variables for local column model diagnostic variables
USE var_lookup,only:iLookPROG                 ! named variables for local column model prognostic variables
USE var_lookup,only:iLookINDEX                ! named variables for local column index variables
USE var_lookup,only:iLookFREQ                 ! named variables for the frequency structure
USE var_lookup,only:iLookGRID                 ! named variables for the glacier grid information
USE var_lookup,only:iLookDECISIONS            ! named variables for elements of the decision structure
USE var_lookup,only:iLookVarType              ! named variables for variable types

! generic variable types
USE nr_type                                   ! variable types, etc.
USE data_types,only:&
                    ! gru dimension
                    gru_double,          &    ! x%gru(:)%var(:)     (dp)
                    ! gru+hru dimension
                    gru_hru_double,      &    ! x%gru(:)%hru(:)%var(:)     (dp)
                    ! gru+hru+dom dimension
                    gru_hru_dom_int,     &    ! x%gru(:)%hru(:)%dom(:)%var(:)     (i4b)
                    gru_hru_dom_double        ! x%gru(:)%hru(:)%dom(:)%var(:)     (dp)
                       
! metadata structure
USE data_types,only:var_info                  ! data type for metadata
USE data_types,only:extended_info             ! data type for extended metadata
USE data_types,only:struct_info               ! summary information on all data structures

! model write options
USE mDecisions_module,only: &
 writePerStep             , &                 ! read forcing data per time step (defualt)
 writeFullSeries                              ! read full forcing series

! safety: set private unless specified otherwise
implicit none
private
public::summa_writeOutputFiles
contains

 ! used to define/write output files
 subroutine summa_writeOutputFiles(modelTimeStep, summa1_struc, err, message)
 ! ---------------------------------------------------------------------------------------
 ! * desired modules
 ! ---------------------------------------------------------------------------------------
 ! data types
 USE summa_type,only:summa1_type_dec                         ! master summa data type
 ! subroutines and functions
 USE time_utils_module,only:elapsedSec                       ! calculate the elapsed time
 USE mDecisions_module,only:mDecisions                       ! module to read model decisions
 USE summa_alarms,only:summa_setWriteAlarms                  ! set alarms to control model output
 USE summa_defineOutput,only:summa_defineOutputFiles         ! define summa output files
 USE modelwrite_module,only:writeRestart,writeRestartGlac    ! module to write model Restart
 USE modelwrite_module,only:writeData                        ! module to write model output
 USE modelwrite_module,only:writeTime                        ! module to write model time
 USE output_stats,only:calcStats                             ! module for compiling output statistics
 ! global data: general
 USE globalData,only:forcingStep                             ! index of current time step in current forcing file
 USE globalData,only:gru_struc                               ! gru-hru mapping structures
 USE globalData,only:structInfo                              ! information on the data structures
 ! global data: time structures
 USE globalData,only:oldTime                                 ! time from the previous time step
 USE globalData,only:finshTime                               ! end time of simulation
 ! global data: decisions for model alarms
 USE globalData,only:ixProgress                              ! define frequency to write progress
 USE globalData,only:ixRestart                               ! define frequency to write restart files
 USE globalData,only:newOutputFile                           ! define option for new output files
 ! buffered write
 USE globalData,only:numtim                                  ! number of time steps
 ! controls on statistics output
 USE globalData,only:statCounter                             ! time counter for stats
 USE globalData,only:resetStats                              ! flags to reset statistics
 USE globalData,only:finalizeStats                           ! flags to finalize statistics
 USE globalData,only:outputTimeStep                          ! timestep in output files
 ! output constraints
 USE globalData,only:maxLayers                               ! maximum number of layers
 USE globalData,only:maxGlaciers                             ! maximum number of glaciers in a GRU
 ! timing variables
 USE globalData,only:startWrite,endWrite                     ! date/time for the start and end of the model writing
 USE globalData,only:elapsedWrite                            ! elapsed time to write data
 ! file information
 USE summaFileManager,only:OUTPUT_PATH,OUTPUT_PREFIX         ! define output file
 USE summaFileManager,only:STATE_PATH                        ! optional path to state output files (defaults to OUTPUT_PATH)
 USE globalData,only:output_fileSuffix                       ! suffix for the output & state files (optional summa argument)
 USE globalData,only:maxDOM                                  ! maximum number of domains in any HRU
 USE globalData,only:nHRUrun                                 ! number of HRU in the run space
 USE globalData,only:nGRUrun                                 ! number of GRU in the run space
 ! ---------------------------------------------------------------------------------------
 ! * variables
 ! ---------------------------------------------------------------------------------------
 implicit none
 ! dummy variables
 integer(i4b),intent(in)               :: modelTimeStep            ! time step index
 type(summa1_type_dec),intent(inout)   :: summa1_struc             ! master summa data structure
 integer(i4b),intent(out)              :: err                      ! error code
 character(*),intent(out)              :: message                  ! error message
 ! local variables
 character(len=256)                    :: timeString                 ! portion of restart file name that contains the write-out time
 character(len=256)                    :: restartFile                ! restart file name
 character(len=256)                    :: restartGlacFile            ! glacier restart file name
 logical(lgt)                          :: printRestart=.false.       ! flag to print a re-start file
 logical(lgt)                          :: printProgress=.false.      ! flag to print simulation progress
 logical(lgt)                          :: defNewOutputFile=.false.   ! flag to define new output files
 logical(lgt)                          :: is_writingOutput=.false.   ! flag to write model output
 logical(lgt)                          :: is_bufferedWrite=.false.   ! flag for buffered write
 integer(i4b)                          :: iTime                      ! indices of time
 integer(i4b)                          :: iGRU,iHRU,iDOM             ! indices of GRUs and HRUs
 integer(i4b)                          :: iStruct                    ! index of model structure
 integer(i4b)                          :: iFreq                      ! index of the output frequency
 integer(i4b)                          :: maxWrite                   ! maximum number of time steps written 
 type(var_info)      , allocatable     :: meta(:)                    ! metadata
 type(extended_info) , allocatable     :: stat_meta(:)               ! statistics metadata (includes only desired variables)
 integer(i4b)        , allocatable     :: child_map(:)               ! index of element in child data structure -- meta(map(ivar)) = stat_meta(ivar)
 class(*)            , allocatable     :: timestepData(:)            ! vector timestep data (unlimited polymorphic structure) 
 class(*)            , allocatable     :: bufferData(:)              ! vector buffer data (unlimited polymorphic structure) 
 class(*)            , allocatable     :: statsData(:)               ! vector stats data (unlimited polymorphic structure)
 ! error control
 integer(i4b)                          :: ierr                       ! error code of downwind routine
 character(LEN=256)                    :: cmessage                   ! error message of downwind routine
 ! ---------------------------------------------------------------------------------------
 ! associate to elements in the data structure
 summaVars: associate(&
  ! statistics structures
  forcStat             => summa1_struc%forcStat    , & ! x%gru(:)%hru(:)%var(:)%dat        -- model forcing data
  progStat             => summa1_struc%progStat    , & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model prognostic (state) variables
  diagStat             => summa1_struc%diagStat    , & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model diagnostic variables
  fluxStat             => summa1_struc%fluxStat    , & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model fluxes
  indxStat             => summa1_struc%indxStat    , & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model indices
  bvarStat             => summa1_struc%bvarStat    , & ! x%gru(:)%var(:)%dat               -- basin-average variabl
  ! primary data structures
  timeStruct           => summa1_struc%timeStruct  , & ! x%var(:)                          -- model time data
  forcStruct           => summa1_struc%forcStruct  , & ! x%gru(:)%hru(:)%var(:)            -- model forcing data
  indxStruct           => summa1_struc%indxStruct  , & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model indices
  progStruct           => summa1_struc%progStruct  , & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model prognostic (state) variables
  diagStruct           => summa1_struc%diagStruct  , & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model diagnostic variables
  fluxStruct           => summa1_struc%fluxStruct  , & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model fluxes
  bvarStruct           => summa1_struc%bvarStruct  , & ! x%gru(:)%var(:)%dat               -- basin-average variables
  gridStruct           => summa1_struc%gridStruct  , & ! x%gru(:)%grid(:)%var(:)%dat2(:,:) -- basin grid parameters and variables
  ! miscellaneous variables
  nGRU                 => summa1_struc%nGRU        , & ! number of grouped response units
  nHRU                 => summa1_struc%nHRU        , & ! number of global hydrologic response units
  nDOM                 => summa1_struc%nDOM          & ! number of global domains (max in any HRU)
 ) ! assignment to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='summa_writeOutputFiles/'

 ! identify the start of the writing
 call date_and_time(values=startWrite)

 ! *****************************************************************************
 ! *** initialize statistics
 ! *****************************************************************************

 ! initialize the statistics flags
 if(modelTimeStep==1)then

  ! initialize time step index
  statCounter(1:maxVarFreq) = 1
  outputTimeStep(1:maxVarFreq) = 1

  ! initialize flags to reset/finalize statistics
  resetStats(:)    = .true.   ! start by resetting statistics
  finalizeStats(:) = .false.  ! do not finalize stats on the first time step

  ! set stats flag for the timestep-level output
  finalizeStats(iLookFREQ%timestep)=.true.

  ! initialize number of domains, hru, and gru in global data
  nGRUrun = nGRU
  nHRUrun = nHRU
  maxDOM  = nDOM

 endif  ! if the first time step

 ! *****************************************************************************
 ! *** initialize/populate data structures for the buffered write
 ! *****************************************************************************

 ! if a buffered write
 if(model_decisions(iLookDECISIONS%write_buff)%iDecision == writeFullSeries)then

  ! initialize data structures for the buffered write
  if(modelTimeStep == 1)then
   call initBufferedWrite(structInfo,numtim,ierr,cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); err=20; return; endif
  endif   ! (if the first time step)

  ! populate data structures
  call popBufferStruct(structInfo,summa1_struc,modelTimeStep,ierr,cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); err=20; return; endif

  ! set the maximum number of data points to write
  maxWrite = numtim

 else   ! (if buffered write)

  ! standard case of write one data value per time step
  maxWrite = 1

 endif  ! (if not buffered write)

 ! *****************************************************************************
 ! *** set alarms for writing data
 ! *****************************************************************************

 ! set alarms to control model output
 call summa_setWriteAlarms(modelTimeStep,                              &   ! time index
                           oldTime%var, timeStruct%var, finshTime%var, &   ! time vectors
                           newOutputFile, defNewOutputFile,            &   ! flag to define new output file
                           ixRestart,     printRestart,                &   ! flag to print the restart file
                           ixProgress,    printProgress,               &   ! flag to print simulation progress
                           resetStats,    finalizeStats,               &   ! flags to reset and finalize stats
                           statCounter,                                &   ! statistics counter
                           err, cmessage)                                  ! error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! identify the need to write output file
 select case(model_decisions(iLookDECISIONS%write_buff)%iDecision)
  case(writePerStep);    is_writingOutput = .true.
  case(writeFullSeries); is_writingOutput = (modelTimeStep == numtim)
  case default
   err=10; message=trim(message)//"unknown option for method used to write model output [option="//trim(model_decisions(iLookDECISIONS%write_buff)%cDecision)//"]"; return
 end select

 ! check if the buffered write
 is_bufferedWrite = (model_decisions(iLookDECISIONS%write_buff)%iDecision == writeFullSeries .and. modelTimeStep == numtim)

 ! print progress
 if(printProgress) write(*,'(i4,1x,5(i2,1x))') timeStruct%var(1:5)

 ! *****************************************************************************
 ! *** define summa output files
 ! *****************************************************************************

 ! check the need to create a new output file
 if(defNewOutputFile .or. modelTimeStep==1)then

  ! define summa output files
  call summa_defineOutputFiles(modelTimeStep, summa1_struc, err, cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! re-initialize the indices for model writing
  outputTimeStep(:)=1

 end if  ! if defining a new file

 ! ****************************************************************************
 ! *** calculate output statistics
 ! ****************************************************************************

 ! loop through GRUs and HRUs
 do iGRU=1,nGRU
  do iHRU=1,gru_struc(iGRU)%hruCount
    do iDOM=1,gru_struc(iGRU)%hruInfo(iHRU)%domCount

     ! calculate output statistics
     do iStruct=1,size(structInfo)
      select case(trim(structInfo(iStruct)%structName))
       case('prog'); call calcStats(progStat%gru(iGRU)%hru(iHRU)%dom(iDOM)%var,progStruct%gru(iGRU)%hru(iHRU)%dom(iDOM)%var,statProg_meta,resetStats,finalizeStats,statCounter,err,cmessage)
       case('diag'); call calcStats(diagStat%gru(iGRU)%hru(iHRU)%dom(iDOM)%var,diagStruct%gru(iGRU)%hru(iHRU)%dom(iDOM)%var,statDiag_meta,resetStats,finalizeStats,statCounter,err,cmessage)
       case('flux'); call calcStats(fluxStat%gru(iGRU)%hru(iHRU)%dom(iDOM)%var,fluxStruct%gru(iGRU)%hru(iHRU)%dom(iDOM)%var,statFlux_meta,resetStats,finalizeStats,statCounter,err,cmessage)
       case('indx'); call calcStats(indxStat%gru(iGRU)%hru(iHRU)%dom(iDOM)%var,indxStruct%gru(iGRU)%hru(iHRU)%dom(iDOM)%var,statIndx_meta,resetStats,finalizeStats,statCounter,err,cmessage)
      end select
      if(err/=0)then; message=trim(message)//trim(cmessage)//'['//trim(structInfo(iStruct)%structName)//']'; return; endif
     end do  ! (looping through structures)

    end do  ! (looping through domains)

    ! calculate output statistics for the HRU-level variables
    call calcStats(forcStat%gru(iGRU)%hru(iHRU)%var,forcStruct%gru(iGRU)%hru(iHRU)%var,statForc_meta,resetStats,finalizeStats,statCounter,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage)//'[forc stats]'; return; endif

  end do  ! (looping through HRUs)

  ! calculate basin stats
  call calcStats(bvarStat%gru(iGRU)%var(:),bvarStruct%gru(iGRU)%var(:),statBvar_meta,resetStats,finalizeStats,statCounter,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage)//'[bvar stats]'; return; endif

 end do  ! (looping through GRUs)

 ! ****************************************************************************
 ! *** write integer time information
 ! ****************************************************************************

 ! NOTE: This is uncommon -- users rarely require integer time variables (iyyy, im, id, ...) because these can be retrieved from julday

 ! NOTE: writing integer time variables is currently restricted to the writePerStep option
 if(model_decisions(iLookDECISIONS%write_buff)%iDecision == writePerStep)then
  call writeTime(finalizeStats,outputTimeStep,time_meta,timeStruct%var,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 else
  if(any(time_meta(:)%varDesire))then
   message=trim(message)//"currently no capabilities to output integer time data when using the buffered write option (writeFullSeries)"
   err=10; return
  endif  ! (if desire integer write)
 endif  ! (if not the writePerStep option)

 ! ****************************************************************************
 ! *** write variables for each HRU
 ! ****************************************************************************

 ! write the model output to the NetCDF file
 if(is_writingOutput)then

  ! loop through data structures
  do iStruct=1,size(structInfo)  ! loop means we can apply error code at the end
   select case(trim(structInfo(iStruct)%structName))

    ! define names of desired data structures
    case('indx','forc','diag','prog','flux','bvar')  ! restrict attention to the variables that we are interested in

    ! get metadata for desired structures
    call get_metadata(trim(structInfo(iStruct)%structName), meta, stat_meta, child_map, ierr, cmessage)
    if(ierr/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! get statistics data structure as a 1-element vector (maxWrite=1)
    call get_statisticVec(trim(structInfo(iStruct)%structName), maxWrite, summa1_struc, statsData, ierr, cmessage)
    if(ierr/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! ----- write buffered data --------------------------------------------------
    if(is_bufferedWrite)then

     ! get saved data for the buffered write
     call get_savedBuffer(trim(structInfo(iStruct)%structName), bufferData, ierr, cmessage) 
     if(ierr/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

     ! will only write the (buffered) timestep data but full routine called for convenience
     call writeData(is_bufferedWrite,finalizeStats,outputTimeStep,maxWrite,statsData(1),bufferData,meta,child_map,indxStruct,ierr,cmessage)
     if(ierr/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! ----- write data and statistics structures ---------------------------------
    else  

     ! check per step write
     if(maxWrite/=1)then; err=20; message=trim(message)//'expect maxWrite=1'; return; endif

     ! get timestep data structure as a 1-element vector (maxWrite=1)
     call get_timestepVec(trim(structInfo(iStruct)%structName), maxWrite, summa1_struc, timestepData, ierr, cmessage)
     if(ierr/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

     ! Passes the full metadata structure and the child map (rather than the stats metadata structure) because
     !  we have the option to write out data of types other than statistics.
     call writeData(is_bufferedWrite,finalizeStats,outputTimeStep,maxWrite,statsData(1),timestepData,meta,child_map,indxStruct,ierr,cmessage)
     if(ierr/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    endif  ! (write data and statistics structures)

    ! ----------------------------------------------------------------------------

    ! just keep going if not interested in a data structure
    case default; cycle
   end select  ! select data structure

  end do  ! (looping through data structures)
 endif  ! (if writing output)

 ! *****************************************************************************
 ! *** write restart file
 ! *****************************************************************************

 ! print a restart file if requested
 if(printRestart)then
  write(timeString,'(i4,3(i2.2))') timeStruct%var(iLookTIME%iyyy),timeStruct%var(iLookTIME%im),timeStruct%var(iLookTIME%id),timeStruct%var(iLookTIME%ih)
  
  if(STATE_PATH == '') then
    restartFile=trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//'_restart_'//trim(timeString)//trim(output_fileSuffix)//'.nc'
  else
    restartFile=trim(STATE_PATH)//trim(OUTPUT_PREFIX)//'_restart_'//trim(timeString)//trim(output_fileSuffix)//'.nc'
  endif

  call writeRestart(restartFile,nGRU,nHRU,nDOM,prog_meta,progStruct,bvar_meta,bvarStruct,maxLayers,indx_meta,indxStruct,err,cmessage)  
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  if(maxGlaciers>0)then
    if(STATE_PATH == '') then
      restartGlacFile=trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//'_restartGrid_'//trim(timeString)//trim(output_fileSuffix)//'.nc'
    else
      restartGlacFile=trim(STATE_PATH)//trim(OUTPUT_PREFIX)//'_restartGrid_'//trim(timeString)//trim(output_fileSuffix)//'.nc'
    endif
    call writeRestartGlac(restartGlacFile,nGRU,grid_meta,gridStruct,err,cmessage) 
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif
 end if

 ! *****************************************************************************
 ! *** update counters
 ! *****************************************************************************

 ! increment output file timestep
 do iFreq = 1,maxvarFreq
  statCounter(iFreq) = statCounter(iFreq)+1
  if(finalizeStats(iFreq)) outputTimeStep(iFreq) = outputTimeStep(iFreq) + 1
 end do

 ! increment forcingStep
 forcingStep=forcingStep+1

 ! if finalized stats, then reset stats on the next time step
 resetStats(:) = finalizeStats(:)

 ! save time vector
 oldTime%var(:) = timeStruct%var(:)

 ! *****************************************************************************
 ! *** finalize
 ! *****************************************************************************

 ! identify the end of the writing
 call date_and_time(values=endWrite)

 ! aggregate the elapsed time for model writing
 elapsedWrite = elapsedWrite + elapsedSec(startWrite, endWrite)

 ! end associate statements
 end associate summaVars

 end subroutine summa_writeOutputFiles

 ! *****************************************************************************
 ! *****************************************************************************
 
 ! private subroutine: initialize data structures for buffered write
 subroutine initBufferedWrite(structInfo,numtim,err,message)
 ! use modules
 USE allocspace_module,only:allocGlobal                      ! module to allocate space
 ! use global data 
 USE globalData,only:fullIndxSave                            ! x(:)%gru(:)%hru(:)%dom(:)%var(:) -- saved output for indices
 USE globalData,only:fullForcSave                            ! x(:)%gru(:)%hru(:)%var(:)        -- saved output for forcing
 USE globalData,only:fullProgSave                            ! x(:)%gru(:)%hru(:)%dom(:)%var(:) -- saved output for prognostic variables
 USE globalData,only:fullDiagSave                            ! x(:)%gru(:)%hru(:)%dom(:)%var(:) -- saved output for diagnostic variables
 USE globalData,only:fullFluxSave                            ! x(:)%gru(:)%hru(:)%dom(:)%var(:) -- saved output for flux variables
 USE globalData,only:fullBvarSave                            ! x(:)%gru(:)%var(:)               -- saved output for basin variables
 implicit none
 ! dummy variables
 type(struct_info)    , intent(in)     :: structInfo(:)      ! information on the data structures
 integer(i4b)         , intent(in)     :: numtim             ! number of model time steps
 integer(i4b)         , intent(out)    :: err                ! error code
 character(*)         , intent(out)    :: message            ! error message
 ! local variables -- temporary data structures
 integer(i4b)                          :: iStruct            ! index of data structure
 type(gru_hru_dom_int),   allocatable  :: tempIndxStruct     ! Indx temp structure: x%gru(:)%hru(:)%dom(:)%var(:) (rkind)
 type(gru_hru_double),    allocatable  :: tempForcStruct     ! Forc temp structure: x%gru(:)%hru(:)%var(:)        (rkind)
 type(gru_hru_dom_double),allocatable  :: tempProgStruct     ! Prog temp structure: x%gru(:)%hru(:)%dom(:)%var(:) (rkind)
 type(gru_hru_dom_double),allocatable  :: tempDiagStruct     ! Diag temp structure: x%gru(:)%hru(:)%dom(:)%var(:) (rkind)
 type(gru_hru_dom_double),allocatable  :: tempFluxStruct     ! Flux temp structure: x%gru(:)%hru(:)%dom(:)%var(:) (rkind)
 type(gru_double),        allocatable  :: tempBvarStruct     ! Bvar temp structure: x%gru(:)%var(:)               (rkind)
 ! error control
 integer(i4b)                          :: ierr               ! error code of downwind routine
 character(LEN=256)                    :: cmessage           ! error message of downwind routine
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='initBufferedWrite/'

 ! allocate space for local data structures
 allocate(tempIndx_struct, tempForc_struct, tempProg_struct, tempDiag_struct, tempFlux_struct, tempBvar_struct, stat=ierr)
 if(ierr/=0)then; err=20; message=trim(message)//'problem allocating temporary data structures'; return; endif

 ! allocate space for temporary data structuress
 do iStruct=1,size(structInfo)  ! loop means we can apply error code at the end
  select case(trim(structInfo(iStruct)%structName))
   case('indx'); call allocGlobal(statIndx_meta%var_info, tempIndx_struct, ierr, cmessage)
   case('forc'); call allocGlobal(statForc_meta%var_info, tempForc_struct, ierr, cmessage)
   case('prog'); call allocGlobal(statProg_meta%var_info, tempProg_struct, ierr, cmessage)
   case('diag'); call allocGlobal(statDiag_meta%var_info, tempDiag_struct, ierr, cmessage)
   case('flux'); call allocGlobal(statFlux_meta%var_info, tempFlux_struct, ierr, cmessage)
   case('bvar'); call allocGlobal(statBvar_meta%var_info, tempBvar_struct, ierr, cmessage)
   case default; cycle ! do not expect additional data structures
  end select
  if(ierr/=0)then; message=trim(message)//trim(cmessage)//'['//trim(structInfo(iStruct)%structName)//']'; err=20; return; endif
 end do  ! (looping through structures)

 ! add a time dimension
 do iStruct=1,size(structInfo)  ! loop means we can apply error code at the end
  select case(trim(structInfo(iStruct)%structName))
   case('indx'); allocate(fullIndxSave(numtim), source=tempIndx_struct, stat=ierr)
   case('forc'); allocate(fullForcSave(numtim), source=tempForc_struct, stat=ierr)
   case('prog'); allocate(fullProgSave(numtim), source=tempProg_struct, stat=ierr)
   case('diag'); allocate(fullDiagSave(numtim), source=tempDiag_struct, stat=ierr)
   case('flux'); allocate(fullFluxSave(numtim), source=tempFlux_struct, stat=ierr)
   case('bvar'); allocate(fullBvarSave(numtim), source=tempBvar_struct, stat=ierr)
   case default; cycle ! do not expect additional data structures
  end select
  if(ierr/=0)then; message=trim(message)//trim(cmessage)//'['//trim(structInfo(iStruct)%structName)//']'; err=20; return; endif
 end do  ! (looping through structues)

 ! deallocate space for data structures
 deallocate(tempIndx_struct, tempForc_struct, tempProg_struct, tempDiag_struct, tempFlux_struct, tempBvar_struct, stat=ierr)
 if(ierr/=0)then; err=20; message=trim(message)//'problem deallocating temporary data structures'; return; endif

 end subroutine initBufferedWrite

 ! *****************************************************************************
 ! *****************************************************************************
 ! private subroutine: populate data structures for buffered write
 subroutine popBufferStruct(structInfo,summaStruct,iTime,err,message)
 ! global data: structures for buffered write
 USE globalData,only:fullIndxSave                         ! x(:)%gru(:)%hru(:)%dom(:)%var(:) -- saved output for indices
 USE globalData,only:fullForcSave                         ! x(:)%gru(:)%hru(:)%var(:)        -- saved output for forcing
 USE globalData,only:fullProgSave                         ! x(:)%gru(:)%hru(:)%dom(:)%var(:) -- saved output for prognostic variables
 USE globalData,only:fullDiagSave                         ! x(:)%gru(:)%hru(:)%dom(:)%var(:) -- saved output for diagnostic variables
 USE globalData,only:fullFluxSave                         ! x(:)%gru(:)%hru(:)%dom(:)%var(:) -- saved output for flux variables
 USE globalData,only:fullBvarSave                         ! x(:)%gru(:)%var(:)               -- saved output for basin variables
 ! global data: structures for GRU-HRU topology
 USE globalData,only:gru_struc                            ! gru-hru mapping structures
 ! derived types
 USE summa_type,only:summa1_type_dec                      ! master summa data type
 implicit none
 ! input variables
 type(struct_info)    , intent(in)    :: structInfo(:)    ! information on the data structures
 type(summa1_type_dec), intent(in)    :: summaStruct      ! master summa data structure
 integer(i4b)         , intent(in)    :: iTime            ! index of time (modelTimeStep)
 ! output variables
 integer(i4b)         , intent(out)   :: err              ! error code
 character(*)         , intent(out)   :: message          ! error message
 ! local variables
 integer(i4b)                         :: iStruct          ! index of data structure
 integer(i4b)                         :: iGRU             ! index of GRU
 integer(i4b)                         :: iHRU             ! index of HRU
 integer(i4b)                         :: iDOM             ! index of domain
 integer(i4b)                         :: iVar             ! index of variable
 integer(i4b)                         :: pVar             ! index of "parent" variable (i.e., index in the data structure)
 type(var_info)       , allocatable   :: meta(:)          ! metadata
 type(extended_info)  , allocatable   :: stat_meta(:)     ! statistics metadata (includes only desired variables)
 integer(i4b)         , allocatable   :: child_map(:)     ! index of element in child data structure -- meta(map(ivar)) = stat_meta(ivar) 
 ! error control
 integer(i4b)                         :: ierr             ! local error code
 character(len=256)                   :: cmessage         ! error message of the downwind routine
 ! associate to elements in the data structure
 ! ----------------------------------------------------------------------------------------------------------------------------
 ! primary data structures
 summaAssociate: associate(&
  indxStruct           => summaStruct%indxStruct  , & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model indices
  forcStruct           => summaStruct%forcStruct  , & ! x%gru(:)%hru(:)%var(:)            -- model forcing data
  progStruct           => summaStruct%progStruct  , & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model prognostic (state) variables
  diagStruct           => summaStruct%diagStruct  , & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model diagnostic variables
  fluxStruct           => summaStruct%fluxStruct  , & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model fluxes
  bvarStruct           => summaStruct%bvarStruct  , & ! x%gru(:)%var(:)%dat               -- basin-average variables
  nGRU                 => summaStruct%nGRU          &
  ) ! assignment to variables in the data structures
 ! -------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='popBufferStruct/'

 ! loop through GRUs and HRUs
 do iGRU=1,nGRU
  do iHRU=1,gru_struc(iGRU)%hruCount
   do iDOM=1,gru_struc(iGRU)%hruInfo(iHRU)%domCount

    ! loop through data structures
    do iStruct=1,size(structInfo)  ! loop means we can apply error code at the end

     ! ----- get metadata for desired structure -----------------------------------------------------------------------------

     ! select structure
     select case(trim(structInfo(iStruct)%structName))

      ! get metadata for desired structures
      case('indx','forc','diag','prog','flux','bvar')  ! restrict attention to the variables that we are interested in
      call get_metadata(trim(structInfo(iStruct)%structName), meta, stat_meta, child_map, ierr, cmessage)
      if(ierr>0)then; err=20; message=trim(message)//trim(cmessage); return; endif

      ! just keep going if not interested in a data structure
      case default; cycle 

     end select  ! select data structure

     ! ----- populate data for desired variables ----------------------------------------------------------------------------

     ! NOTE: use stat_meta because these are the variables desired in the output
     do iVar=1,size(stat_meta) ! skip if size=0
     
      ! don't do anything if var is not requested or not a scalar variable
      if (.not.stat_meta(iVar)%varDesire) cycle                 ! variable not requested
      if (stat_meta(iVar)%varType/=iLookVarType%outstat) cycle  ! not a scalar variable
     
      ! index in parent structure
      pVar = stat_meta(iVar)%ixParent
      if(trim(stat_meta(iVar)%varName) /= trim(meta(pVar)%varName))then
       message=trim(message)//'variable names do not match'
       err=20; return
      endif
     
      ! populate GRU+HRU+DOM structures
      select case(trim(structInfo(iStruct)%structName))
       case('indx'); fullIndxSave(iTime)%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar) = indxStruct%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(pVar)%dat(1)
       case('forc')  ! HRU-only data structure
        if(iDOM==1)  fullForcSave(iTime)%gru(iGRU)%hru(iHRU)%var(iVar) = forcStruct%gru(iGRU)%hru(iHRU)%var(pVar)
       case('diag'); fullProgSave(iTime)%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar) = progStruct%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(pVar)%dat(1)
       case('prog'); fullDiagSave(iTime)%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar) = diagStruct%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(pVar)%dat(1)
       case('flux'); fullFluxSave(iTime)%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar) = fluxStruct%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(pVar)%dat(1)
       case('bvar')  ! GRU-only data structure
        if(iHRU==1 .and. iDOM==1) fullBvarSave(iTime)%gru(iGRU)%var(iVar) = bvarStruct%gru(iGRU)%var(pVar)%dat(1)
       case default; err=20; message=trim(message)//'do not expect any other structures than what is listed'; return
      end select

     end do  ! (looping through variables)

    end do  ! (looping through structures)

    ! deallocate space for metadata structures
    deallocate(meta, stat_meta, child_map, stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating metadata'; err=20; return; endif

   end do  ! (looping through domains)
  end do  ! (looping through HRUs)
 end do  ! (looping through GRUs)

 end associate summaAssociate

 end subroutine popBufferStruct

 ! *****************************************************************************
 ! *****************************************************************************
                                     
 subroutine get_metadata(tag, meta, stat_meta, child_map, err, message)
 ! dummy variables
 character(*)        , intent(in)               :: tag            ! name of data structure
 type(var_info)      , intent(out), allocatable :: meta(:)        ! metadata
 type(extended_info) , intent(out), allocatable :: stat_meta(:)   ! statistics metadata
 integer(i4b)        , intent(out), allocatable :: child_map(:)   ! index of the child data structure
 integer(i4b)        , intent(out)              :: err            ! error code
 character(*)        , intent(out)              :: message        ! error message
 ! local variables
 integer(i4b)                                   :: nVar           ! number of variables in the metadata structure
 integer(i4b)                                   :: nSubset        ! number of variables in the statistics metadata structure
 integer(i4b)                                   :: ierr           ! local error code 
 ! initalize error control
 err=0; message='get_metadata/'

 ! get size of metadata structures for a given structure
 select case(trim(tag))
  case('indx'); nVar = size(indx_meta); nSubset = size(statIndx_meta)
  case('forc'); nVar = size(forc_meta); nSubset = size(statForc_meta)
  case('prog'); nVar = size(prog_meta); nSubset = size(statProg_meta)
  case('diag'); nVar = size(diag_meta); nSubset = size(statDiag_meta)
  case('flux'); nVar = size(flux_meta); nSubset = size(statFlux_meta)
  case('bvar'); nVar = size(bvar_meta); nSubset = size(statBvar_meta)
  case default; err=10; message=trim(message)//'data structure not needed -- should not get here'; return  
 end select

 ! allocate space for the metadata structures
 allocate(meta(nVar), stat_meta(nSubset), child_map(nVar), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating metadata'; err=20; return; endif

 ! save metadata
 select case(trim(tag))
  case('indx'); meta(:) = indx_meta(:); stat_meta(:) = statIndx_meta(:); child_map(:) = indxChild_map(:)
  case('forc'); meta(:) = forc_meta(:); stat_meta(:) = statForc_meta(:); child_map(:) = forcChild_map(:)
  case('prog'); meta(:) = prog_meta(:); stat_meta(:) = statProg_meta(:); child_map(:) = progChild_map(:)
  case('diag'); meta(:) = diag_meta(:); stat_meta(:) = statDiag_meta(:); child_map(:) = diagChild_map(:)
  case('flux'); meta(:) = flux_meta(:); stat_meta(:) = statFlux_meta(:); child_map(:) = fluxChild_map(:)
  case('bvar'); meta(:) = bvar_meta(:); stat_meta(:) = statBvar_meta(:); child_map(:) = bvarChild_map(:)
 end select

 end subroutine get_metadata

 ! *****************************************************************************
 ! *****************************************************************************

 subroutine get_savedBuffer(tag, structVec, err, message)
 ! saved timestep data
 USE globalData,only:fullIndxSave                         ! x(:)%gru(:)%hru(:)%dom(:)%var(:) -- saved output for indices
 USE globalData,only:fullForcSave                         ! x(:)%gru(:)%hru(:)%var(:)        -- saved output for forcing
 USE globalData,only:fullProgSave                         ! x(:)%gru(:)%hru(:)%dom(:)%var(:) -- saved output for prognostic variables
 USE globalData,only:fullDiagSave                         ! x(:)%gru(:)%hru(:)%dom(:)%var(:) -- saved output for diagnostic variables
 USE globalData,only:fullFluxSave                         ! x(:)%gru(:)%hru(:)%dom(:)%var(:) -- saved output for flux variables
 USE globalData,only:fullBvarSave                         ! x(:)%gru(:)%var(:)               -- saved output for basin variables
 ! dummy variables
 character(*) , intent(in)               :: tag            ! name of data structure
 class(*)     , intent(out), allocatable :: structVec(:)   ! vector of any data structure
 integer(i4b) , intent(out)              :: err            ! error code
 character(*) , intent(out)              :: message        ! error message
 ! initialize errror control
 err=0; message='get_savedBuffer/'

 ! allocate data
 select case(trim(tag))
  case('indx'); allocate(structVec, source=fullIndxSave, stat=err)
  case('forc'); allocate(structVec, source=fullForcSave, stat=err)
  case('prog'); allocate(structVec, source=fullProgSave, stat=err)
  case('diag'); allocate(structVec, source=fullDiagSave, stat=err)
  case('flux'); allocate(structVec, source=fullFluxSave, stat=err)
  case('bvar'); allocate(structVec, source=fullBvarSave, stat=err)
  case default; err=20; message=trim(message)//'cannot identify data structure'; return
 end select

 ! check errors
 if(err/=0)then
  message=trim(message)//'problem allocating space for structure '//trim(tag)
  err=20; return
 endif

 end subroutine get_savedBuffer

 ! *****************************************************************************
 ! *****************************************************************************

 subroutine get_timestepVec(tag, nVec, summaStruct, structVec, err, message)
 USE summa_type,only:summa1_type_dec                                ! master summa data type
 ! dummy variables
 character(*)          , intent(in)               :: tag            ! name of data structure
 integer(i4b)          , intent(in)               :: nVec           ! number of vector elements
 type(summa1_type_dec) , intent(in)               :: summaStruct    ! master summa data structure
 class(*)              , intent(out), allocatable :: structVec(:)   ! vector of any data structure
 integer(i4b)          , intent(out)              :: err            ! error code
 character(*)          , intent(out)              :: message        ! error message
 ! initialize errror control
 err=0; message='get_timestepVec/'

 ! allocate data 
 select case(trim(tag))
  case('indx'); allocate(structVec(nVec), source=summaStruct%indxStruct, stat=err)
  case('forc'); allocate(structVec(nVec), source=summaStruct%forcStruct, stat=err)
  case('prog'); allocate(structVec(nVec), source=summaStruct%progStruct, stat=err)
  case('diag'); allocate(structVec(nVec), source=summaStruct%diagStruct, stat=err)
  case('flux'); allocate(structVec(nVec), source=summaStruct%fluxStruct, stat=err)
  case('bvar'); allocate(structVec(nVec), source=summaStruct%bvarStruct, stat=err)
  case default; err=20; message=trim(message)//'cannot identify data structure'; return
 end select

 ! check errors
 if(err/=0)then
  message=trim(message)//'problem allocating space for structure '//trim(tag)
  err=20; return
 endif

 end subroutine get_timestepVec

 ! *****************************************************************************
 ! *****************************************************************************

 subroutine get_statisticVec(tag, nVec, summaStruct, structVec, err, message)
 USE summa_type,only:summa1_type_dec                                ! master summa data type
 ! dummy variables
 character(*)          , intent(in)               :: tag            ! name of data structure
 integer(i4b)          , intent(in)               :: nVec           ! number of vector elements
 type(summa1_type_dec) , intent(in)               :: summaStruct    ! master summa data structure
 class(*)              , intent(out), allocatable :: structVec(:)   ! vector of any data structure
 integer(i4b)          , intent(out)              :: err            ! error code
 character(*)          , intent(out)              :: message        ! error message
 ! initialize errror control
 err=0; message='get_statisticVec/'

 ! allocate data
 select case(trim(tag))
  case('indx'); allocate(structVec(nVec), source=summaStruct%indxStat, stat=err)
  case('forc'); allocate(structVec(nVec), source=summaStruct%forcStat, stat=err)
  case('prog'); allocate(structVec(nVec), source=summaStruct%progStat, stat=err)
  case('diag'); allocate(structVec(nVec), source=summaStruct%diagStat, stat=err)
  case('flux'); allocate(structVec(nVec), source=summaStruct%fluxStat, stat=err)
  case('bvar'); allocate(structVec(nVec), source=summaStruct%bvarStat, stat=err)
  case default; err=20; message=trim(message)//'cannot identify data structure'; return
 end select

 ! check errors
 if(err/=0)then
  message=trim(message)//'problem allocating space for structure '//trim(tag)
  err=20; return
 endif

 end subroutine get_statisticVec

 ! *****************************************************************************
 ! *****************************************************************************

end module summa_writeOutput


