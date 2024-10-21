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

module run_oneGRU_module

! numerical recipes data types
USE nrtype

! global constants
USE globalData,only: yes,no             ! .true. and .false.
USE globalData,only: data_step          ! length of data step (s)

! define data types
USE data_types,only:&
                    ! GRU-to-HRU mapping
                    gru2hru_map,       & ! HRU info
                    ! no spatial dimension
                    var_i,             & ! x%var(:)            (i4b)
                    var_d,             & ! x%var(:)            (rkind)
                    var_ilength,       & ! x%var(:)%dat        (i4b)
                    var_dlength,       & ! x%var(:)%dat        (rkind)
                    ! no variable dimension
                    hru_i,             & ! x%hru(:)            (i4b)
                    hru_dom_d,         & ! x%hru(:)%dom(:)     (rkind)
                    ! hru dimension
                    hru_int,           & ! x%hru(:)%var(:)     (i4b)
                    hru_int8,          & ! x%hru(:)%var(:)     (i8b)
                    hru_double,        & ! x%hru(:)%var(:)     (rkind)
                    hru_intVec,        & ! x%hru(:)%var(:)%dat (i4b)
                    !hru+dom dimension
                    hru_dom_intVec,    & ! x%hru(:)%dom(:)%var(:)%dat (i4b)
                    hru_dom_double,    & ! x%hru(:)%dom(:)%var(:)     (rkind)
                    hru_dom_doubleVec, & ! x%hru(:)%dom(:)%var(:)%dat (rkind)
                    ! hru+dom+ z dimension
                    hru_dom_z_vLookup    ! x%hru(:)%z(:)%var(:)%lookup(:)

! provide access to the named variables that describe elements of parameter structures
USE var_lookup,only:iLookTYPE          ! look-up values for classification of veg, soils etc.
USE var_lookup,only:iLookID            ! look-up values for hru and gru IDs
USE var_lookup,only:iLookATTR          ! look-up values for local attributes
USE var_lookup,only:iLookINDEX         ! look-up values for local column index variables
USE var_lookup,only:iLookFLUX          ! look-up values for local column model fluxes
USE var_lookup,only:iLookBPAR          ! look-up values for basin-average model parameters
USE var_lookup,only:iLookBVAR          ! look-up values for basin-average model variables
USE var_lookup,only:iLookTIME          ! look-up values for model time data
USE var_lookup,only:iLookPROG          ! look-up values for model prognostic (state) variables

! provide access to model decisions
USE globalData,only:model_decisions    ! model decision structure
USE var_lookup,only:iLookDECISIONS     ! look-up values for model decisions
USE globalData,only:data_step          ! length of data step (s)

! access domain types
USE globalData,only:upland             ! domain type for upland areas
USE globalData,only:glacAcc            ! domain type for glacier accumulation areas
USE globalData,only:glacAbl            ! domain type for glacier ablation areas
USE globalData,only:wetland            ! domain type for wetland areas

! provide access to the named variables that describe model decisions
USE mDecisions_module,only:&           ! look-up values for the choice of method for the spatial representation of groundwater
 localColumn, &                        ! separate groundwater representation in each local soil column
 singleBasin, &                        ! single groundwater store over the entire basin
 bigBucket                             ! a big bucket (lumped aquifer model)
! -----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public::run_oneGRU
contains

! ************************************************************************************************
! public subroutine run_oneGRU: simulation for a single GRU
! ************************************************************************************************

! simulation for a single GRU
subroutine run_oneGRU(&
                      ! model control
                      gruInfo,            & ! intent(inout): HRU information for given GRU (# HRUs, #snow+soil layers)
                      dt_init,            & ! intent(inout): used to initialize the length of the sub-step for each HRU
                      ixComputeVegFlux,   & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
                      ! data structures (input)
                      timeVec,            & ! intent(in):    model time data
                      typeHRU,            & ! intent(in):    local classification of soil veg etc. for each HRU
                      idHRU,              & ! intent(in):    local classification of hru and gru IDs
                      attrHRU,            & ! intent(in):    local attributes for each HRU
                      lookupHRU,          & ! intent(in):    local lookup tables for each HRU
                      ! data structures (input-output)
                      mparHRU,            & ! intent(in):    local model parameters
                      bparData,           & ! intent(in):    basin model parameters
                      indxHRU,            & ! intent(inout): model indices
                      forcHRU,            & ! intent(inout): model forcing data
                      progHRU,            & ! intent(inout): prognostic variables for a local HRU
                      diagHRU,            & ! intent(inout): diagnostic variables for a local HRU
                      fluxHRU,            & ! intent(inout): model fluxes for a local HRU
                      bvarData,           & ! intent(inout): basin-average variables
                      ! error control
                      err,message)          ! intent(out):   error control

  ! ----- define downstream subroutines -----------------------------------------------------------------------------------
  USE run_oneHRU_module,only:run_oneHRU                       ! module to run for one HRU
  USE time_utils_module,only:compjulday                       ! convert calendar date to julian day
  USE qTimeDelay_module,only:qGlacier                         ! module to route water through glacier (time lapse)
  USE qTimeDelay_module,only:qOverland                        ! module to route water through an "unresolved" river network
  ! ----- define dummy variables ------------------------------------------------------------------------------------------
  implicit none
  ! model control
  type(gru2hru_map)       , intent(inout) :: gruInfo              ! HRU information for given GRU (# HRUs, #snow+soil layers)
  type(hru_dom_d)         , intent(inout) :: dt_init              ! used to initialize the length of the sub-step for each domain
  type(hru_i)             , intent(inout) :: ixComputeVegFlux     ! flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
  ! data structures (input)
  type(var_i)             , intent(in)    :: timeVec              ! x%var(:)                   -- model time data
  type(hru_int)           , intent(in)    :: typeHRU              ! x%hru(:)%var(:)            -- local classification of soil veg etc. for each HRU
  type(hru_int8)          , intent(in)    :: idHRU                ! x%hru(:)%var(:)            -- local classification of hru and gru IDs
  type(hru_double)        , intent(in)    :: attrHRU              ! x%hru(:)%var(:)            -- local attributes for each HRU
  type(hru_dom_z_vLookup) , intent(in)    :: lookupHRU            ! x%hru(:)%dom(:)%z(:)%var(:)%lookup(:) -- lookup values for each HRU
  ! data structures (input-output)
  type(hru_dom_doubleVec) , intent(in)    :: mparHRU              ! x%hru(:)%dom(:)%var(:)%dat -- local (HRU) model parameters
  type(var_d)             , intent(in)    :: bparData             ! x%var                      -- basin-average parameters
  type(hru_dom_intVec)    , intent(inout) :: indxHRU              ! x%hru(:)%dom(:)%var(:)%dat -- model indices
  type(hru_double)        , intent(inout) :: forcHRU              ! x%hru(:)%dom(:)%var(:)     -- model forcing data
  type(hru_dom_doubleVec) , intent(inout) :: progHRU              ! x%hru(:)%dom(:)%var(:)%dat -- model prognostic (state) variables
  type(hru_dom_doubleVec) , intent(inout) :: diagHRU              ! x%hru(:)%dom(:)%var(:)%dat -- model diagnostic variables
  type(hru_dom_doubleVec) , intent(inout) :: fluxHRU              ! x%hru(:)%dom(:)%var(:)%dat -- model fluxes
  type(var_dlength)       , intent(inout) :: bvarData             ! x%var(:)%dat               -- basin-average variables
  ! error control
  integer(i4b)            , intent(out)   :: err                  ! error code
  character(*)            , intent(out)   :: message              ! error message
  ! ----- define local variables ------------------------------------------------------------------------------------------
  ! general local variables
  character(len=256)                  :: cmessage               ! error message
  integer(i4b)                        :: iHRU                   ! HRU index
  integer(i4b)                        :: jHRU,kHRU              ! index of the hydrologic response unit
  integer(i4b)                        :: iDOM                   ! domain index
  real(rkind)                         :: fracDOM                ! fractional area of a given HRU domain in GRU (-)
  integer(i4b)                        :: ndom_glacGRU           ! number of glacier domains in the GRU
  real(rkind), allocatable            :: elev(:)                ! median elevation of the each glacier cell (m)
  real(rkind), allocatable            :: GWE_deltaYr(:)         ! change in glacier water equivalent per year (m) in each glacier cell
  logical(lgt)                        :: computeVegFluxFlag     ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
  logical(lgt)                        :: updateGlacArea         ! flag to update glacier area
  logical(lgt)                        :: updateLakeArea         ! flag to update lake area
  real(rkind)                         :: currentJulDay          ! current julian day
  real(rkind)                         :: updateJulDay           ! julian day to update glacier area
  real(rkind)                         :: remaining_area         ! remaining area to be distributed
  real(rkind)                         :: remaining_elev         ! remaining elevation to be distributed
  logical(lgt)                        :: runHRU                 ! flag to run the HRU
  logical(lgt)                        :: check_updateGlacArea   ! flag to check if glacier area needs to be updated
  real(rkind)                         :: glacIceMelt            ! glacier ice reservoir melt (m3 s-1)
  real(rkind)                         :: glacSnowMelt           ! glacier snow reservoir melt (m3 s-1)
  real(rkind)                         :: glacFirnMelt           ! glacier firn reservoir melt (m3 s-1)

  ! initialize error control
  err=0; write(message, '(A21,I0,A10,I0,A2)' ) 'run_oneGRU (gru_nc = ',gruInfo%gru_nc,', gruId = ',gruInfo%gru_id,')/'

  ! ----- basin initialization --------------------------------------------------------------------------------------------
  ! initialize runoff variables
  bvarData%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)    = 0._rkind  ! surface runoff (m s-1)
  bvarData%var(iLookBVAR%basin__SoilDrainage)%dat(1)     = 0._rkind  ! soil drainage (m s-1)
  bvarData%var(iLookBVAR%basin__ColumnOutflow)%dat(1)    = 0._rkind  ! outflow from all "outlet" HRUs (those with no downstream HRU)
  bvarData%var(iLookBVAR%basin__TotalRunoff)%dat(1)      = 0._rkind  ! total runoff to the channel from all active components (m s-1)

  ! initialize baseflow variables
  bvarData%var(iLookBVAR%basin__AquiferRecharge)%dat(1)  = 0._rkind ! recharge to the aquifer (m s-1)
  bvarData%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  = 0._rkind ! baseflow from the aquifer (m s-1)
  bvarData%var(iLookBVAR%basin__AquiferTranspire)%dat(1) = 0._rkind ! transpiration loss from the aquifer (m s-1)

  ! initialize glacier variables
  glacIceMelt  = 0._rkind ! glacier ice reservoir melt (m3 s-1)
  glacSnowMelt = 0._rkind ! glacier snow reservoir melt (m3 s-1)
  glacFirnMelt = 0._rkind ! glacier firn reservoir melt (m3 s-1)
  bvarData%var(iLookBVAR%basin__GlacierArea)%dat(1) = 0._rkind ! basin glacier area (m2)

  updateGlacArea = .false. ! initialize updateGlacArea flag
  updateLakeArea = .false. ! initialize updateLakeArea flag

  ! initialize total inflow for each layer in a soil column and glacier size allocation
  ndom_glacGRU = 0 ! initialize number of glacier domains in the GRU
  check_updateGlacArea = .true.
  do iHRU=1,gruInfo%hruCount
    do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
      fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%mLayerColumnInflow)%dat(:) = 0._rkind
      if (progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)==0._rkind) cycle ! skip domains with no area
      if (gruInfo%hruInfo(iHRU)%domInfo(iDOM)%dom_type==glacAcc .or. gruInfo%hruInfo(iHRU)%domInfo(iDOM)%dom_type==glacAbl)then
        if (check_updateGlacArea) then
          ! update glacier area every October 1st
          ! compute the julian day at the start of the year
          call compjulday(timeVec%var(iLookTIME%iyyy),           & ! input  = year
                          10, 1, 1, 1, 0._rkind,                 & ! input  = month, day, hour, minute, second
                          updateJulDay,err,cmessage)               ! output = julian day (fraction of day) + error control
          if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

          ! compute the fractional julian day for the current time step
          call compjulday(timeVec%var(iLookTIME%iyyy),           & ! input  = year
                          timeVec%var(iLookTIME%im),             & ! input  = month
                          timeVec%var(iLookTIME%id),             & ! input  = day
                          timeVec%var(iLookTIME%ih),             & ! input  = hour
                          timeVec%var(iLookTIME%imin),0._rkind,  & ! input  = minute/second
                          currentJulDay,err,cmessage)              ! output = julian day (fraction of day) + error control
          if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
          if (updateJulDay == currentJulDay) updateGlacArea = .true. ! update glacier area if a year passed from last update
          check_updateGlacArea = .false. ! only check this once
        end if

        if (updateGlacArea) then ! allocate space for glacier area and GWE_deltaYr
          ndom_glacGRU = ndom_glacGRU + 1
        else
          ndom_glacGRU = 1 ! allocate at some size
        endif
      endif
    end do
  end do
  allocate(elev(ndom_glacGRU), GWE_deltaYr(ndom_glacGRU))

  ! ********** RUN FOR ONE HRU ********************************************************************************************
  ! loop through HRUs
  ndom_glacGRU = 0 ! initialize number of glacier domains in the GRU

  do iHRU=1,gruInfo%hruCount
    
    ! skip HRUs with no area
    runHRU = .false.
    do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
      if (progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)>0._rkind) runHRU = .true.
    end do
    if (.not. runHRU) cycle

    computeVegFluxFlag = (ixComputeVegFlux%hru(iHRU) == yes)  ! initialize the flag to compute the vegetation flux
    ! ----- run the model --------------------------------------------------------------------------------------------------

    ! simulation for a single HRU
    call run_oneHRU(&
                   ! model control
                   gruInfo%hruInfo(iHRU)%hru_nc,   & ! intent(in):    hru count Id
                   gruInfo%hruInfo(iHRU)%hru_id,   & ! intent(in):    hruId
                   dt_init%hru(iHRU),              & ! intent(inout): initial time step
                   computeVegFluxFlag,             & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
                   gruInfo%hruInfo(iHRU)%domCount, & ! intent(in):    total number of domains
                   gruInfo%hruInfo(iHRU)%domInfo,  & ! intent(inout): domain type and layer information
                   ! data structures (input)
                   typeHRU%hru(iHRU),              & ! intent(in):    local classification of soil veg etc. for each HRU
                   attrHRU%hru(iHRU),              & ! intent(in):    local attributes for each HRU
                   lookupHRU%hru(iHRU),            & ! intent(in):    local lookup tables for each HRU
                   bvarData,                       & ! intent(in):    basin-average model variables
                   ! data structures (input-output)
                   mparHRU%hru(iHRU),              & ! intent(in):    model parameters
                   indxHRU%hru(iHRU),              & ! intent(inout): model indices
                   forcHRU%hru(iHRU),              & ! intent(inout): model forcing data
                   progHRU%hru(iHRU),              & ! intent(inout): model prognostic variables for a local HRU
                   diagHRU%hru(iHRU),              & ! intent(inout): model diagnostic variables for a local HRU
                   fluxHRU%hru(iHRU),              & ! intent(inout): model fluxes for a local HRU
                   ! error control
                   err,cmessage)                      ! intent(out):   error control
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! save the flag for computing the vegetation fluxes
    if(computeVegFluxFlag)       ixComputeVegFlux%hru(iHRU) = yes
    if(.not. computeVegFluxFlag) ixComputeVegFlux%hru(iHRU) = no

    ! ----- compute fluxes across HRUs --------------------------------------------------------------------------------------------------
    ! identify lateral connectivity
    ! (Note:  for efficiency, this could this be done as a setup task, not every timestep)
    kHRU = 0
    ! identify the downslope HRU
    dsHRU: do jHRU=1,gruInfo%hruCount
      if(typeHRU%hru(iHRU)%var(iLookTYPE%downHRUindex) == idHRU%hru(jHRU)%var(iLookID%hruId))then
        if(kHRU==0)then  ! check there is a unique match
          kHRU=jHRU
          exit dsHRU
        end if  ! (check there is a unique match)
      end if  ! (if identified a downslope HRU)
    end do dsHRU

    do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
      if (progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)==0._rkind) cycle ! skip domains with no area
      associate(typeDOM => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%dom_type)

      ! identify the area covered by the current domain
      fracDOM = progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)/ bvarData%var(iLookBVAR%basin__totalArea)%dat(1)

      ! if lateral flows are active, add inflow to the downslope HRU (not active on glacier)
      if(kHRU > 0)then  ! if there is a downslope HRU
        fluxHRU%hru(kHRU)%dom(iDOM)%var(iLookFLUX%mLayerColumnInflow)%dat(:) = fluxHRU%hru(kHRU)%dom(iDOM)%var(iLookFLUX%mLayerColumnInflow)%dat(:)  + fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%mLayerColumnOutflow)%dat(:)
      else ! otherwise just increment basin (GRU) column outflow (m3 s-1) with the hru fraction
        bvarData%var(iLookBVAR%basin__ColumnOutflow)%dat(1) = bvarData%var(iLookBVAR%basin__ColumnOutflow)%dat(1) + sum(fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%mLayerColumnOutflow)%dat(:))
      end if

      ! ----- calculate weighted basin (GRU) fluxes --------------------------------------------------------------------------------------
      if (typeDOM==upland)then
         ! increment basin surface runoff (m s-1)
        bvarData%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)  = bvarData%var(iLookBVAR%basin__SurfaceRunoff)%dat(1) + fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%scalarSurfaceRunoff)%dat(1)*fracDOM

        ! increment basin soil drainage (m s-1)
        bvarData%var(iLookBVAR%basin__SoilDrainage)%dat(1)   = bvarData%var(iLookBVAR%basin__SoilDrainage)%dat(1)  + fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%scalarSoilDrainage)%dat(1) *fracDOM

        ! increment aquifer variables -- ONLY if aquifer baseflow is computed individually for each HRU and aquifer is run
        ! NOTE: groundwater computed later for singleBasin
        ! NOTE: no groundwater for glacier
        if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == localColumn .and. model_decisions(iLookDECISIONS%groundwatr)%iDecision == bigBucket) then
          bvarData%var(iLookBVAR%basin__AquiferRecharge)%dat(1)  = bvarData%var(iLookBVAR%basin__AquiferRecharge)%dat(1)  + fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%scalarAquiferRecharge)%dat(1) *fracDOM
          bvarData%var(iLookBVAR%basin__AquiferTranspire)%dat(1) = bvarData%var(iLookBVAR%basin__AquiferTranspire)%dat(1) + fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%scalarAquiferTranspire)%dat(1)*fracDOM
          bvarData%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  = bvarData%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  + fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%scalarAquiferBaseflow)%dat(1) *fracDOM
        end if
      else if (typeDOM==glacAcc .or. typeDOM==glacAbl)then
        if (typeDOM==glacAcc)then ! collect glacier accumulation melt m s-1
          glacFirnMelt = glacFirnMelt + fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%scalarGlacierMelt)%dat(1) *fracDOM
        else if (typeDOM==glacAbl)then ! collect glacier ablation melt m s-1
          if (progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowDepth)%dat(1)>0._rkind)then
            glacSnowMelt = glacSnowMelt + fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%scalarGlacierMelt)%dat(1) *fracDOM
          else
            glacIceMelt  = glacIceMelt  + fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%scalarGlacierMelt)%dat(1) *fracDOM
          endif
        end if
        bvarData%var(iLookBVAR%basin__GlacierArea)%dat(1) = bvarData%var(iLookBVAR%basin__GlacierArea)%dat(1) + progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)
        ! placeholder line, add actual kg/m2 of glacier storage instead of SWE, or is it delSWE + snowMelt + scalarGlacierMelt
        bvarData%var(iLookBVAR%basin__GlacierStorage)%dat(1) = bvarData%var(iLookBVAR%basin__GlacierStorage)%dat(1) + progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSWE)%dat(1) * progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)
        ! if a year passed from last glacier area update, write fluxes to the output file so that the glacier area can be updated
        if (updateGlacArea) then
          ! save the glacier mass balance associated with this elevation
          ndom_glacGRU = ndom_glacGRU + 1
          elev(ndom_glacGRU) = progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMelev)%dat(1)
          ! placeholder line, add actual kg/m2 of glacier storage instead of SWE
          GWE_deltaYr(ndom_glacGRU) = progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSWE)%dat(1) !- progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSWE_yrend)%dat(1) 
        end if
      else if (typeDOM==wetland)then ! collect wetland fluxes
        ! STUB:  wetland fluxes not yet implemented
        print*, 'WARNING:  wetland fluxes not yet implemented'
      end if
      end associate
    end do ! (looping through domains)

    ! averaging more fluxes (and/or states) can be added to this section as desired

  end do  ! (looping through HRUs)
  ! ********** END LOOP THROUGH HRUS **************************************************************************************
  ! lapse glacier fluxes to the basin by routing through each glacier
  call qGlacier(&
                ! input
                data_step,                                          & ! intent(in):    length of data step (s)
                bparData%var(iLookBPAR%glacStor_kIce),              & ! intent(in):    storage coefficient ice reservoir (hours)
                bparData%var(iLookBPAR%glacStor_kFirn),             & ! intent(in):    storage coefficient snow reservoir (hours)
                bparData%var(iLookBPAR%glacStor_kFirn),             & ! intent(in):    storage coefficient firn reservoir (hours)
                glacIceMelt,                                        & ! intent(in):    total melt into ice reservoirs (m s-1)
                glacSnowMelt,                                       & ! intent(in):    total melt into snow reservoirs (m s-1)
                glacFirnMelt,                                       & ! intent(in):    total melt into firn reservoirs (m s-1)
                bvarData%var(iLookBVAR%glacAblArea)%dat,            & ! intent(in):    per glacier ablation area (m2)
                bvarData%var(iLookBVAR%glacAccArea)%dat,            & ! intent(in):    per glacier accumulation area (m2)
                ! output
                bvarData%var(iLookBVAR%glacIceRunoffFuture)%dat,    & ! intent(inout): per glacier ice reservoir runoff in future time steps (m s-1)
                bvarData%var(iLookBVAR%glacSnowRunoffFuture)%dat,   & ! intent(inout): per glacier snow reservoir runoff in future time steps (m s-1)
                bvarData%var(iLookBVAR%glacFirnRunoffFuture)%dat,   & ! intent(inout): per glacier firn reservoir runoff in future time steps (m s-1)
                bvarData%var(iLookBVAR%glacierRoutedRunoff)%dat(1), & ! intent(out):   routed glacier runoff (m s-1)
                err,message)              ! error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
 
  ! perform the routing
  associate(totalArea => bvarData%var(iLookBVAR%basin__totalArea)%dat(1) )
  
  ! compute water balance for the basin aquifer
  if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == singleBasin)then
    message=trim(message)//'multi_driver/bigBucket groundwater code not transferred from old code base yet'
    err=20; return
  end if

  ! calculate total runoff depending on whether aquifer is connected
  if(model_decisions(iLookDECISIONS%groundwatr)%iDecision == bigBucket) then
    ! aquifer
    bvarData%var(iLookBVAR%basin__TotalRunoff)%dat(1) = bvarData%var(iLookBVAR%basin__SurfaceRunoff)%dat(1) + bvarData%var(iLookBVAR%basin__ColumnOutflow)%dat(1)/totalArea + bvarData%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)
  else
    ! no aquifer
    bvarData%var(iLookBVAR%basin__TotalRunoff)%dat(1) = bvarData%var(iLookBVAR%basin__SurfaceRunoff)%dat(1) + bvarData%var(iLookBVAR%basin__ColumnOutflow)%dat(1)/totalArea + bvarData%var(iLookBVAR%basin__SoilDrainage)%dat(1)
  endif
 
  call qOverland(&
                 ! input
                 model_decisions(iLookDECISIONS%subRouting)%iDecision, & ! intent(in):    index for routing method
                 bvarData%var(iLookBVAR%basin__TotalRunoff)%dat(1),    & ! intent(in):    total runoff to the channel from all active components (m s-1)
                 bvarData%var(iLookBVAR%routingFractionFuture)%dat,    & ! intent(in):    fraction of runoff in future time steps (m s-1)
                 bvarData%var(iLookBVAR%routingRunoffFuture)%dat,      & ! intent(inout): runoff in future time steps (m s-1)
                 ! output
                 bvarData%var(iLookBVAR%averageInstantRunoff)%dat(1),  & ! intent(out):   instantaneous runoff (m s-1)
                 bvarData%var(iLookBVAR%averageRoutedRunoff)%dat(1),   & ! intent(out):   routed runoff (m s-1)
                 err,message)                                            ! intent(out):   error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

  ! add glacier runoff to overland runoff
  bvarData%var(iLookBVAR%averageInstantRunoff)%dat(1) = bvarData%var(iLookBVAR%averageInstantRunoff)%dat(1) + bvarData%var(iLookBVAR%glacierRoutedRunoff)%dat(1)
  bvarData%var(iLookBVAR%averageRoutedRunoff)%dat(1) = bvarData%var(iLookBVAR%averageRoutedRunoff)%dat(1) + bvarData%var(iLookBVAR%glacierRoutedRunoff)%dat(1)

  end associate

  ! Need to update the glacier area
  if (updateGlacArea) then
    ! need to save length, bottom topo, and elevation of glaciers from the end of previous update for this GRU in file associated with gruInfo%gru_id
    ! need to associate each glacier with an HRU and domain
    ! for nGlacier  = gru_struc(iGRU)%nGlacier, skip if no area
    !call flow_MUSCL(& 
    !               ! input
    !               gruInfo%gru_id, & ! intent(in): GRU ID
    !               GWE_deltaYr,    & ! intent(in): change in glacier water equivalent per year (m)
    !               elev,           & ! intent(in): median elevation of the glacier domain (m)
    !               ! output ??
    !               elev_ELA,       & ! intent(out): elevation of the equilibrium line altitude (m)
    !               glac_length,    & ! intent(out): length of the glacier (m)
    !               glac_area,      & ! intent(out): area of the glacier (m2)
    !               err,message)      ! intent(out): error control
    !if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! update the glacier area and elevation by HRU, need to do this inside flow_MUSCL somehow, then as adjust domain area and elevation adjust upland area and elevation
    !progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMelev)%dat(1)
    !progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)
    ! and bvarData%var(iLookBVAR%glacAblArea)%dat,
    !     bvarData%var(iLookBVAR%glacAccArea)%dat,

    ! STUB:  glacier area update not yet implemented
    ! call glacier flow model here using the file above ...
    ! get new fracDOM for each HRU and domain, and new glacier area and new elevMedians
  end if
  deallocate(elev, GWE_deltaYr)

  if (updateGlacArea .or. updateLakeArea) then
    do iHRU=1,gruInfo%hruCount
      ! update the HRU area and elevation
      remaining_area = attrHRU%hru(iHRU)%var(iLookATTR%HRUarea)
      remaining_elev = attrHRU%hru(iHRU)%var(iLookATTR%HRUarea)*attrHRU%hru(iHRU)%var(iLookATTR%elevation)
      do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
        if (gruInfo%hruInfo(iHRU)%domInfo(iDOM)%dom_type.ne.upland) then
          remaining_area = remaining_area - progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)
          remaining_elev = remaining_elev - progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1) * progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMelev)%dat(1)
        end if
      end do
      do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
        if (gruInfo%hruInfo(iHRU)%domInfo(iDOM)%dom_type==upland) then
          progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1) = remaining_area
          if(remaining_area>0.0_rkind) then 
            progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMelev)%dat(1) = remaining_elev/remaining_area
          else
            progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMelev)%dat(1) = 0.0_rkind
            progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1) = attrHRU%hru(iHRU)%var(iLookATTR%elevation)
          end if
        end if
      end do
    end do
  end if

 end subroutine run_oneGRU

end module run_oneGRU_module
