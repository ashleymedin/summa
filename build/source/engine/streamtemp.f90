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

module streamtemp_module

! Stream temperature on the coupled river network.
!
! Water is routed by mizuRoute; heat is routed here. Each GRU may hold one stream HRU, whose stream domain
! is the water column of the reach the GRU drains to (nSnow + nLake + nSoil layers and the aquifer, no
! upland domain). The physics of that column live in the usual SUMMA routines (the lake layers of
! coupled_em, with the advective source of snowLakeSoilGlceNrgFlux). This module does the network part,
! following Wanders et al. (2019, WRR, DynWat) after van Beek et al. (2012, WRR): once the land HRUs have
! run and mizuRoute has routed this step's runoff, it walks the reaches from upstream to downstream and
!   * mixes the outflow of the upstream reaches into the temperature of the inflow,
!   * hands the stream domain its reach inflow, lateral inflow and their temperatures, its outflow and
!     the reach depth and velocity, and runs the column (run_oneHRU with streamPass=.true.),
!   * reads back the column temperature as the temperature of the water leaving the reach.
! A reach without a stream HRU just mixes what reaches it, so every reach of the network has an outlet
! temperature and the stream HRUs can be sparse.

! data types
USE nr_type
USE data_types,only:stream_network        ! per-reach hydraulics and temperatures
USE data_types,only:gru2hru_map           ! GRU-to-HRU mapping
USE data_types,only:var_i                 ! x%var(:)            (i4b)
USE data_types,only:var_d                 ! x%var(:)            (rkind)
USE data_types,only:gru_hru_i             ! x%gru(:)%hru(:)     (i4b)
USE data_types,only:gru_hru_dom_d         ! x%gru(:)%hru(:)%dom(:) (rkind)
USE data_types,only:gru_hru_int           ! x%gru(:)%hru(:)%var(:) (i4b)
USE data_types,only:gru_hru_double        ! x%gru(:)%hru(:)%var(:) (rkind)
USE data_types,only:gru_hru_dom_doubleVec ! x%gru(:)%hru(:)%dom(:)%var(:)%dat (rkind)
USE data_types,only:gru_hru_dom_intVec    ! x%gru(:)%hru(:)%dom(:)%var(:)%dat (i4b)
USE data_types,only:gru_hru_dom_z_vLookup ! x%gru(:)%hru(:)%dom(:)%z(:)%var(:)%lookup(:)
USE data_types,only:gru_doubleVec         ! x%gru(:)%var(:)%dat (rkind)

! physical constants
USE multiconst,only:Tfreeze               ! freezing point of pure water (K)

! named variables
USE var_lookup,only:iLookTYPE             ! look-up values for classification of veg, soils etc.
USE var_lookup,only:iLookFLUX             ! look-up values for local column model fluxes
USE var_lookup,only:iLookDIAG             ! look-up values for local column model diagnostic variables
USE var_lookup,only:iLookPROG             ! look-up values for local column model prognostic variables
USE var_lookup,only:iLookBVAR             ! look-up values for basin-average variables

! constants
USE globalData,only:yes,no                ! .true. and .false.
USE globalData,only:stream                ! horizontal domain type for a stream reach

implicit none
private
public :: stream_domain_map
public :: run_streamNetwork

contains

! ************************************************************************************************
! public subroutine stream_domain_map: locate the stream HRU and domain of each GRU
! ************************************************************************************************
! Per GRU: the reach id the stream HRU asks for (attribute streamSegId, 0 = the reach the GRU drains
! to), the index of that HRU within the GRU and of the stream domain within the HRU (0 = none).
! check_icond has already made sure there is at most one stream HRU per GRU.
subroutine stream_domain_map(nGRU,gru_struc,typeStruct,streamSegId,ixStreamHRU,ixStreamDOM,nStream)
  implicit none
  integer(i4b),intent(in)          :: nGRU              ! number of GRUs
  type(gru2hru_map),intent(in)     :: gru_struc(:)      ! gru-hru mapping structure
  type(gru_hru_int),intent(in)     :: typeStruct        ! local classification of soil veg etc. for each HRU
  integer(i4b),intent(out)         :: streamSegId(nGRU) ! reach id of the stream HRU (0 = reach mapped from the GRU id)
  integer(i4b),intent(out)         :: ixStreamHRU(nGRU) ! index of the stream HRU within the GRU (0 = none)
  integer(i4b),intent(out)         :: ixStreamDOM(nGRU) ! index of the stream domain within that HRU (0 = none)
  integer(i4b),intent(out)         :: nStream           ! number of stream HRUs
  integer(i4b)                     :: iGRU,iHRU,iDOM    ! loop indices
  streamSegId(:) = 0; ixStreamHRU(:) = 0; ixStreamDOM(:) = 0; nStream = 0
  do iGRU=1,nGRU
    do iHRU=1,gru_struc(iGRU)%hruCount
      do iDOM=1,gru_struc(iGRU)%hruInfo(iHRU)%domCount
        if(gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%dom_type==stream)then
          ixStreamHRU(iGRU) = iHRU
          ixStreamDOM(iGRU) = iDOM
          streamSegId(iGRU) = typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%streamSegId)
          nStream = nStream + 1
        end if
      end do
    end do
  end do
end subroutine stream_domain_map

! ************************************************************************************************
! public subroutine run_streamNetwork: the network pass over the stream domains
! ************************************************************************************************
subroutine run_streamNetwork(&
                      ! the river network
                      net,                & ! intent(inout): per-reach hydraulics (in) and temperatures (out)
                      ! model control
                      gru_struc,          & ! intent(inout): HRU information for each GRU (# HRUs, #layers)
                      dt_init,            & ! intent(inout): used to initialize the length of the sub-step for each HRU
                      computeVegFlux,     & ! intent(inout): flag to indicate if we are computing fluxes over vegetation
                      ! data structures (input)
                      typeStruct,         & ! intent(in):    local classification of soil veg etc. for each HRU
                      attrStruct,         & ! intent(in):    local attributes for each HRU
                      lookupStruct,       & ! intent(in):    local lookup tables for each HRU
                      ! data structures (input-output)
                      mparStruct,         & ! intent(in):    local model parameters
                      indxStruct,         & ! intent(inout): model indices
                      forcStruct,         & ! intent(inout): model forcing data
                      progStruct,         & ! intent(inout): prognostic variables for a local HRU
                      diagStruct,         & ! intent(inout): diagnostic variables for a local HRU
                      fluxStruct,         & ! intent(inout): model fluxes for a local HRU
                      bvarStruct,         & ! intent(inout): basin-average variables
                      ! error control
                      err,message)          ! intent(out):   error control
  USE run_oneHRU_module,only:run_oneHRU   ! module to run for one HRU
  implicit none
  ! the river network
  type(stream_network),intent(inout)        :: net              ! per-reach hydraulics and temperatures
  ! model control
  type(gru2hru_map),intent(inout)           :: gru_struc(:)     ! HRU information for each GRU
  type(gru_hru_dom_d),intent(inout)         :: dt_init          ! used to initialize the length of the sub-step for each domain
  type(gru_hru_i),intent(inout)             :: computeVegFlux   ! flag to indicate if we are computing fluxes over vegetation
  ! data structures (input)
  type(gru_hru_int),intent(in)              :: typeStruct       ! local classification of soil veg etc. for each HRU
  type(gru_hru_double),intent(in)           :: attrStruct       ! local attributes for each HRU
  type(gru_hru_dom_z_vLookup),intent(in)    :: lookupStruct     ! lookup tables for each HRU
  ! data structures (input-output)
  type(gru_hru_dom_doubleVec),intent(in)    :: mparStruct       ! local model parameters
  type(gru_hru_dom_intVec),intent(inout)    :: indxStruct       ! model indices
  type(gru_hru_double),intent(inout)        :: forcStruct       ! model forcing data
  type(gru_hru_dom_doubleVec),intent(inout) :: progStruct       ! model prognostic (state) variables
  type(gru_hru_dom_doubleVec),intent(inout) :: diagStruct       ! model diagnostic variables
  type(gru_hru_dom_doubleVec),intent(inout) :: fluxStruct       ! model fluxes
  type(gru_doubleVec),intent(inout)         :: bvarStruct       ! basin-average variables
  ! error control
  integer(i4b),intent(out)                  :: err              ! error code
  character(*),intent(out)                  :: message          ! error message
  ! local variables
  character(len=256)                        :: cmessage         ! error message of downwind routine
  integer(i4b)                              :: iSeq,iSeg,iUps   ! loop indices
  integer(i4b)                              :: jSeg             ! index of an upstream reach
  integer(i4b)                              :: iGRU,iHRU,iDOM   ! indices of the stream domain
  real(rkind)                               :: qUpSum,eUpSum    ! discharge and energy flux arriving from upstream
  logical(lgt)                              :: computeVegFluxFlag
  real(rkind),parameter                     :: minFlow=1.e-12_rkind ! flow below which a temperature is not defined (m3 s-1)
  ! ----------------------------------------------------------------------------------------------------------------------
  err=0; message='run_streamNetwork/'

  ! walk the reaches from upstream to downstream, so the outflow temperature of every upstream reach is known
  do iSeq=1,net%nSeg
    iSeg = net%rchOrder(iSeq)

    ! temperature of the water arriving from upstream: flow-weighted over the upstream reaches
    qUpSum = 0._rkind; eUpSum = 0._rkind
    do iUps=1,net%nUps(iSeg)
      jSeg   = net%ixUps(iUps,iSeg)
      qUpSum = qUpSum + net%qOut(jSeg)
      eUpSum = eUpSum + net%qOut(jSeg)*net%tOut(jSeg)
    end do
    if(qUpSum > minFlow)then
      net%tUp(iSeg) = eUpSum/qUpSum
    else
      net%tUp(iSeg) = net%tLat(iSeg)
    end if
    net%qUp(iSeg) = qUpSum ! consistent with the temperature weighting (the routing sees the same upstream discharge)

    if(net%ixDOM(iSeg) > 0)then
      ! ***** a stream domain stands for this reach: run its water column
      iGRU = net%ixGRU(iSeg); iHRU = net%ixHRU(iSeg); iDOM = net%ixDOM(iSeg)
      associate(fluxDOM => fluxStruct%gru(iGRU)%hru(iHRU)%dom(iDOM), diagDOM => diagStruct%gru(iGRU)%hru(iHRU)%dom(iDOM))
        fluxDOM%var(iLookFLUX%scalarStreamInflow)%dat(1)        = net%qUp(iSeg)
        fluxDOM%var(iLookFLUX%scalarStreamInflowTemp)%dat(1)    = net%tUp(iSeg)
        fluxDOM%var(iLookFLUX%scalarStreamLatInflow)%dat(1)     = net%qLat(iSeg)
        fluxDOM%var(iLookFLUX%scalarStreamLatInflowTemp)%dat(1) = net%tLat(iSeg)
        fluxDOM%var(iLookFLUX%scalarStreamOutflow)%dat(1)       = net%qOut(iSeg)
        diagDOM%var(iLookDIAG%scalarStreamDepth)%dat(1)         = net%depth(iSeg)
        diagDOM%var(iLookDIAG%scalarStreamVelocity)%dat(1)      = net%velocity(iSeg)
      end associate
      computeVegFluxFlag = (computeVegFlux%gru(iGRU)%hru(iHRU) == yes)
      call run_oneHRU(&
                     ! model control
                     gru_struc(iGRU)%hruInfo(iHRU)%hru_nc,   & ! intent(in):    hru count Id
                     gru_struc(iGRU)%hruInfo(iHRU)%hru_id,   & ! intent(in):    hruId
                     dt_init%gru(iGRU)%hru(iHRU),            & ! intent(inout): initial time step
                     computeVegFluxFlag,                     & ! intent(inout): flag to indicate if we are computing fluxes over vegetation
                     gru_struc(iGRU)%hruInfo(iHRU)%domCount, & ! intent(in):    total number of domains
                     gru_struc(iGRU)%hruInfo(iHRU)%domInfo,  & ! intent(inout): domain type and layer information
                     ! data structures (input)
                     typeStruct%gru(iGRU)%hru(iHRU),         & ! intent(in):    local classification of soil veg etc. for each HRU
                     attrStruct%gru(iGRU)%hru(iHRU),         & ! intent(in):    local attributes for each HRU
                     lookupStruct%gru(iGRU)%hru(iHRU),       & ! intent(in):    local lookup tables for each HRU
                     bvarStruct%gru(iGRU),                   & ! intent(in):    basin-average model variables
                     ! data structures (input-output)
                     mparStruct%gru(iGRU)%hru(iHRU),         & ! intent(in):    model parameters
                     indxStruct%gru(iGRU)%hru(iHRU),         & ! intent(inout): model indices
                     forcStruct%gru(iGRU)%hru(iHRU),         & ! intent(inout): model forcing data
                     progStruct%gru(iGRU)%hru(iHRU),         & ! intent(inout): model prognostic variables for a local HRU
                     diagStruct%gru(iGRU)%hru(iHRU),         & ! intent(inout): model diagnostic variables for a local HRU
                     fluxStruct%gru(iGRU)%hru(iHRU),         & ! intent(inout): model fluxes for a local HRU
                     ! error control
                     err,cmessage,                           & ! intent(out):   error control
                     streamPass=.true.)                        ! intent(in):    only the stream domain
      if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
      ! the column temperature is the temperature of the water leaving the reach
      net%tOut(iSeg) = diagStruct%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookDIAG%scalarStreamTemp)%dat(1)

    else
      ! ***** no stream domain: the reach only mixes what reaches it (no exchange with the atmosphere or bed)
      if(net%qUp(iSeg) + net%qLat(iSeg) > minFlow)then
        net%tOut(iSeg) = (net%qUp(iSeg)*net%tUp(iSeg) + net%qLat(iSeg)*net%tLat(iSeg))/(net%qUp(iSeg) + net%qLat(iSeg))
      else
        net%tOut(iSeg) = net%tLat(iSeg)
      end if
    end if
    net%tOut(iSeg) = max(net%tOut(iSeg), Tfreeze)

  end do ! (looping through reaches)

end subroutine run_streamNetwork

end module streamtemp_module
