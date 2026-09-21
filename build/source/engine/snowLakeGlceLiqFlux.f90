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

module snowLakeGlceLiqFlux_module

! access modules
USE nr_type                               ! numerical recipes data types
! physical constants
USE multiconst,only:Tfreeze,     &         ! freezing point of pure water (K)
                    iden_ice,iden_water    ! intrinsic density of ice and water (kg m-3)

! access missing values
USE globalData,only:integerMissing         ! missing integer
USE globalData,only:realMissing            ! missing real number

! physical constants
USE globalData,only:maxVolIceContent       ! snow maximum volumetric ice content to store water (-)
USE globalData,only:iceResidWaterFrac      ! residual volumetric liquid water content in ice (-)

! horizontal domain types
USE globalData,only:wetland                ! sub-GRU lake or pothole: mass balance solved here (not yet implemented)
USE globalData,only:stream                 ! reach water column: mass balance solved by the coupled river network

! named variables
USE var_lookup,only:iLookINDEX             ! named variables for structure elements
USE var_lookup,only:iLookPARAM             ! named variables for structure elements
USE var_lookup,only:iLookPROG              ! named variables for structure elements
USE var_lookup,only:iLookDIAG              ! named variables for structure elements
USE var_lookup,only:iLookFLUX              ! named variables for structure elements
USE var_lookup,only:iLookDERIV             ! named variables for structure elements

! data types
USE data_types,only:var_dlength            ! x%var(:)%dat [rkind]
USE data_types,only:var_ilength            ! x%var(:)%dat [i4b]
USE data_types,only:in_type_snowLakeGlceLiqFlux     ! data type for intent(in) arguments
USE data_types,only:io_type_snowLakeGlceLiqFlux     ! data type for intent(inout) arguments
USE data_types,only:out_type_snowLakeGlceLiqFlux    ! data type for intent(out) arguments

! privacy
implicit none
private
public :: snowIceLiqFlux
public :: lakeLiqFlux
contains
! ************************************************************************************************
! public subroutine snowIceLiqFlux: compute liquid water flux through the snowpack, lake ice, and glacier ice layers
! ************************************************************************************************
subroutine snowIceLiqFlux(&
                      ! input: model control, forcing, and model state vector
                      in_snowLakeGlceLiqFlux,           & ! intent(in):    model control, forcing, and model state vector
                      ! input-output: data structures
                      mpar_data,               & ! intent(in):    model parameters
                      indx_data,               & ! intent(in):    model indices
                      prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                      ! input-output: fluxes and derivatives
                      io_snowLakeGlceLiqFlux,           & ! intent(inout): fluxes and derivatives
                      ! output: error control
                      out_snowLakeGlceLiqFlux)            ! intent(out):   error control
  ! ------------------------------------------------------------------------------------------------------------------------------------------
  ! downwind routines
  USE snow_utils_module,only:fracliquid         ! compute the fraction of liquid water (snow)
  implicit none
  ! input: model control, forcing, and model state vector
  type(in_type_snowLakeGlceLiqFlux) :: in_snowLakeGlceLiqFlux     ! model control, forcing, and model state vector
  ! input-output: data structures
  type(var_dlength),intent(in)      :: mpar_data                  ! model parameters
  type(var_ilength),intent(in)      :: indx_data                  ! model indices
  type(var_dlength),intent(in)      :: prog_data                  ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)   :: diag_data                  ! diagnostic variables for a local HRU
  ! input-output: fluxes and derivatives
  type(io_type_snowLakeGlceLiqFlux)          :: io_snowLakeGlceLiqFlux              ! fluxes and derivatives
  ! output: error control
  type(out_type_snowLakeGlceLiqFlux)         :: out_snowLakeGlceLiqFlux             ! error control
  ! ------------------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  integer(i4b)                      :: nLayers,nStart             ! number of snow/glce layers with water movement and starting layer
  integer(i4b)                      :: iLayer                     ! layer index
  integer(i4b)                      :: ixLayerDesired(1)          ! layer desired (scalar solution)
  integer(i4b)                      :: ixTop                      ! top layer in subroutine call
  integer(i4b)                      :: ixBot                      ! bottom layer in subroutine call
  real(rkind)                       :: multResid                  ! multiplier for the residual water content (-)
  real(rkind)                       :: residThrs                  ! ice density threshold to reduce residual liquid water content (kg m-3)
  real(rkind),parameter             :: residScal=10._rkind        ! scaling factor for residual liquid water content reduction factor (kg m-3)
   real(rkind)                      :: maxVolIceContent_use       ! maximum volumetric ice content depending if snow or firn
  real(rkind)                       :: availCap                   ! available storage capacity [0,1] (-)
  real(rkind)                       :: relSaturn                  ! relative saturation [0,1] (-)
  real(rkind)                       :: iLayerLiqFluxSnLaGl(0:in_snowLakeGlceLiqFlux % nLayers)
  real(rkind)                       :: iLayerLiqFluxSnLaGlDeriv(0:in_snowLakeGlceLiqFlux % nLayers)  
  ! ------------------------------------------------------------------------------------------------------------------------------------------
  ! make association of local variables with information in the data structures
  nLayers = in_snowLakeGlceLiqFlux % nLayers ! get number of snow/glce layers to get water fluxes over
  nStart = in_snowLakeGlceLiqFlux % nStart ! get the start index for the layers
  associate(&
    is_glac           => in_snowLakeGlceLiqFlux % is_glac,                       & ! intent(in): flag to denote if processing a glacier domain
    do_snow           => in_snowLakeGlceLiqFlux % do_snow,                       & ! intent(in): flag to denote if snow is present
    ! input: model control
    firstFluxCall           => in_snowLakeGlceLiqFlux % firstFluxCall,           & ! intent(in): the first flux call
    scalarSolution          => in_snowLakeGlceLiqFlux % scalarSolution,          & ! intent(in): flag to denote if implementing the scalar solution
    ! input: forcing for the top layer
    surface_flux            => in_snowLakeGlceLiqFlux % surface_flux,            & ! intent(in): liquid water flux at the surface (m s-1)
    ! input: water flux at the bottom if already computed
    bottom_flux             => in_snowLakeGlceLiqFlux % bottom_flux,             & ! intent(in): liquid water flux at the bottom if already computed (m s-1)
    ! input: model state vector
    mLayerVolFracLiqTrial   => in_snowLakeGlceLiqFlux % mLayerVolFracLiqTrial,   & ! intent(in): trial value of volumetric fraction of liquid water at the current iteration (-)
    ! input: layer indices
    ixLayerState     => indx_data%var(iLookINDEX%ixLayerState)%dat,              & ! intent(in):    list of indices for all model layers
    ixSnowOnlyHyd    => indx_data%var(iLookINDEX%ixSnowOnlyHyd)%dat,             & ! intent(in):    index in the state subset for hydrology state variables in the snow domain
    ixGlceOnlyHyd    => indx_data%var(iLookINDEX%ixGlceOnlyHyd)%dat,             & ! intent(in):    index in the state subset for hydrology state variables in the glacier ice domain
    ! input: snow properties and parameters
    mLayerVolFracIce => prog_data%var(iLookPROG%mLayerVolFracIce)%dat(nStart+1:nStart+nLayers), & ! intent(in):    volumetric ice content at the start of the time step (-)
    Fcapil           => mpar_data%var(iLookPARAM%Fcapil)%dat(1),                                & ! intent(in):    capillary retention as a fraction of the total pore volume (-)
    k_snow           => mpar_data%var(iLookPARAM%k_snow)%dat(1),                                & ! intent(in):    hydraulic conductivity of snow (m s-1)    
    mw_exp           => mpar_data%var(iLookPARAM%mw_exp)%dat(1),                                & ! intent(in):    exponent for meltwater flow (-)
    snowfrz_scale    => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1),                         & ! intent(in):    freezing curve scaling factor for snow (-)
    ! input-output: diagnostic variables -- only computed for the first iteration
    mLayerPoreSpace  => diag_data%var(iLookDIAG%mLayerPoreSpace)%dat(nStart+1:nStart+nLayers),  & ! intent(inout): pore space in each layer (-)
    mLayerThetaResid => diag_data%var(iLookDIAG%mLayerThetaResid)%dat(nStart+1:nStart+nLayers), & ! intent(inout): residual volumetric liquid water content in each layer (-)
    ! input-output: fluxes and derivatives
    iLayerLiqFluxSnLaGl0      => io_snowLakeGlceLiqFlux % iLayerLiqFluxSnLaGl,           & ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
    iLayerLiqFluxSnLaGlDeriv0 => io_snowLakeGlceLiqFlux % iLayerLiqFluxSnLaGlDeriv,      & ! intent(inout): derivative in vertical liquid water flux at layer interfaces (m s-1)
    ! output: error control
    err                    => out_snowLakeGlceLiqFlux % err,                             & ! intent(out):   error code
    message                => out_snowLakeGlceLiqFlux % cmessage                         & ! intent(out):   error message
    ) ! end association of local variables with information in the data structures
    ! ------------------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='snowIceLiqFlux/'

    ! initialize with index 0
    iLayerLiqFluxSnLaGl = iLayerLiqFluxSnLaGl0(nStart:nLayers+nStart)
    iLayerLiqFluxSnLaGlDeriv = iLayerLiqFluxSnLaGlDeriv0(nStart:nLayers+nStart)

    ! check that the input vectors match nLayers
    if (size(mLayerVolFracLiqTrial)/=nLayers .or. size(mLayerVolFracIce)/=nLayers .or. &
        size(iLayerLiqFluxSnLaGl)/=nLayers+1 .or. size(iLayerLiqFluxSnLaGlDeriv)/=nLayers+1) then
      err=20; message=trim(message)//'size mismatch of input/output vectors'; return
    end if

    ! check the meltwater exponent is >=1
    if (mw_exp<1._rkind) then; err=20; message=trim(message)//'meltwater exponent < 1'; return; end if

    ! get the inputs for the snow layers
    if(is_glac) then  ! snow can be firn
      residThrs = 550._rkind ! firn density threshold to reduce residual liquid water content (kg m-3), maybe should be 800?
      maxVolIceContent_use = min(maxVolIceContent+0.15,0.85_rkind) ! firn maximum volumetric ice content to store water (-)
    else ! snow
      residThrs = 550._rkind ! snow density threshold to reduce residual liquid water
      maxVolIceContent_use = maxVolIceContent ! snow maximum volumetric ice content to store water (-)
    end if

    ! get the indices for the layers
    ixTop = integerMissing
    if (scalarSolution) then
      if (do_snow)then
        ixLayerDesired = pack(ixLayerState, ixSnowOnlyHyd/=integerMissing)
      else
        ixLayerDesired = pack(ixLayerState, ixGlceOnlyHyd/=integerMissing)
      end if
      ixTop = ixLayerDesired(1)
      ixBot = ixLayerDesired(1)
    else
      ixTop = 1
      ixBot = nLayers
    end if

    ! compute properties fixed over the time step
    if (firstFluxCall) then
      ! loop through snow/firn layers
      mLayerPoreSpace  = 1._rkind - mLayerVolFracIce ! compute the pore space (-)
      if(do_snow)then
        do iLayer=1,nLayers ! loop through snow layers
          multResid = 1._rkind/(1._rkind + exp((mLayerVolFracIce(iLayer)*iden_ice - residThrs)/residScal)) ! compute the reduction in liquid water holding capacity at high snow/ice density (-)
          mLayerThetaResid(iLayer) = Fcapil*mLayerPoreSpace(iLayer)*multResid ! compute the residual volumetric liquid water content (-)
        end do  ! end looping through snow/firn layers
      else ! glacier ice
        mLayerThetaResid = iceResidWaterFrac ! 3% reasonable for temperate ice, or non-black ice on lake
      end if  ! end if snow or ice
    end if  ! end if the first flux call
     
    ! compute fluxes
    if (do_snow) then
      if(ixTop==1)then ! compute the liquid flux at the upper boundary (m s-1) if computing top layer as snow
        iLayerLiqFluxSnLaGl(0)      = surface_flux
        iLayerLiqFluxSnLaGlDeriv(0) = 0._rkind ! computed inside computJacob*
      endif
      do iLayer=ixTop,ixBot  ! loop through snow layers
        if (mLayerVolFracLiqTrial(iLayer) > mLayerThetaResid(iLayer)) then ! check that flow occurs
          ! compute the relative saturation (-)
          availCap  = mLayerPoreSpace(iLayer) - mLayerThetaResid(iLayer) ! available capacity
          relSaturn = (mLayerVolFracLiqTrial(iLayer) - mLayerThetaResid(iLayer)) / availCap ! relative saturation
          iLayerLiqFluxSnLaGl(iLayer)      = k_snow*relSaturn**mw_exp
          iLayerLiqFluxSnLaGlDeriv(iLayer) = ( (k_snow*mw_exp)/availCap ) * relSaturn**(mw_exp - 1._rkind)
          if (mLayerVolFracIce(iLayer) > maxVolIceContent_use) then ! NOTE: use start-of-step ice content, to avoid convergence problems
            ! ** allow liquid water to pass through under very high ice density
            iLayerLiqFluxSnLaGl(iLayer) = iLayerLiqFluxSnLaGl(iLayer) + iLayerLiqFluxSnLaGl(iLayer-1)
          end if
        else  ! flow does not occur
          iLayerLiqFluxSnLaGl(iLayer)      = 0._rkind
          iLayerLiqFluxSnLaGlDeriv(iLayer) = 0._rkind
        end if  ! storage above residual content
      end do  ! end loop through snow layers
    else ! ice
      if(ixTop==1) ixTop = 0 ! include the 0 index if the top layer is included, since surface flux downwards is 0 (impermeable) 
      do iLayer=ixBot,ixTop,-1 ! loop through glacier ice layers
        ! ** liquid water goes up since ice is impermeable (upwards direction is negative)
        if (iLayer == nLayers) then ! bottom layer (note, this is nGlce-noThetaChange)
          iLayerLiqFluxSnLaGl(iLayer) = 0._rkind ! no liquid water flux at the bottom of the melting ice layers
          iLayerLiqFluxSnLaGlDeriv(iLayer) = 0._rkind
        else  ! not the bottom layer
          availCap  = min(mLayerVolFracLiqTrial(iLayer+1),mLayerThetaResid(iLayer+1)) ! available capacity
          iLayerLiqFluxSnLaGl(iLayer) = -(mLayerVolFracLiqTrial(iLayer+1) - availCap)
          iLayerLiqFluxSnLaGlDeriv(iLayer) = merge(-1._rkind,0._rkind,mLayerVolFracLiqTrial(iLayer+1)>mLayerThetaResid(iLayer+1)) ! after cancelation, derivative is -1
          ! ** liquid water to passes through ice layers immediately
          iLayerLiqFluxSnLaGl(iLayer) = iLayerLiqFluxSnLaGl(iLayer+1) + iLayerLiqFluxSnLaGl(iLayer)
        end if
      end do  ! end loop through ice layers
    end if  ! end if snow or ice
    if(ixBot==nLayers)then
      iLayerLiqFluxSnLaGl(nLayers) = iLayerLiqFluxSnLaGl(nLayers) + bottom_flux   ! set the bottom flux if already computed
      iLayerLiqFluxSnLaGlDeriv(nLayers) = iLayerLiqFluxSnLaGlDeriv(nLayers) ! may be modified computed inside computJacob, currently bottom flux is always 0 so not needed
    end if

    ! save the results with index 0
    iLayerLiqFluxSnLaGl0(nStart:nLayers+nStart) = iLayerLiqFluxSnLaGl
    iLayerLiqFluxSnLaGlDeriv0(nStart:nLayers+nStart) = iLayerLiqFluxSnLaGlDeriv

  end associate ! end association of local variables with information in the data structures

end subroutine snowIceLiqFlux

! **********************************************************************************************************
! public subroutine lakeLiqFlux: liquid water fluxes through the liquid (unfrozen) lake layers
! **********************************************************************************************************
!  * stream: the reach water column. The water mass is routed by the coupled river network, so no water 
!            moves vertically through the liquid lake layers: every interface flux is zero and the liquid
!            depth is prescribed from the reach volume once per data step. What arrives at the top joins the
!            reach flow at once and is passed on so that its heat can be added to the water column;
!            evaporation leaves the same way.
!  * wetland: a sub-GRU lake or pothole whose water balance would be solved here (inflow at the top, spill
!            over the outlet, seepage to the soil through the bottom flux). Not yet implemented.
subroutine lakeLiqFlux(&
                       ! input: model control, forcing, and model state vector
                       in_snowLakeGlceLiqFlux,  & ! intent(in):    model control, forcing, and model state vector
                       domType,                 & ! intent(in):    horizontal domain type (stream or wetland)
                       surfaceFluxTemp,         & ! intent(in):    temperature of the water arriving at the top (K)
                       topLayerTemp,            & ! intent(in):    trial temperature of the top liquid lake layer (K)
                       ! input-output: data structures
                       indx_data,               & ! intent(in):    model indices
                       prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                       diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                       flux_data,               & ! intent(inout): model fluxes for a local HRU
                       ! input-output: fluxes and derivatives
                       io_snowLakeGlceLiqFlux,  & ! intent(inout): interface fluxes and derivatives
                       ! output: error control
                       out_snowLakeGlceLiqFlux)   ! intent(out):   error control
  implicit none
  ! input: model control, forcing, and model state vector
  type(in_type_snowLakeGlceLiqFlux),intent(in)  :: in_snowLakeGlceLiqFlux     ! model control, forcing, and model state vector
  integer(i4b),intent(in)           :: domType                    ! horizontal domain type
  real(rkind),intent(in)            :: surfaceFluxTemp            ! temperature of the water arriving at the top (K)
  real(rkind),intent(in)            :: topLayerTemp               ! trial temperature of the top liquid lake layer (K)
  ! input-output: data structures
  type(var_ilength),intent(in)      :: indx_data                  ! model indices
  type(var_dlength),intent(in)      :: prog_data                  ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)   :: diag_data                  ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout)   :: flux_data                  ! model fluxes for a local HRU
  ! input-output: fluxes and derivatives
  type(io_type_snowLakeGlceLiqFlux),intent(inout) :: io_snowLakeGlceLiqFlux   ! interface fluxes and derivatives
  ! output: error control
  type(out_type_snowLakeGlceLiqFlux),intent(out)  :: out_snowLakeGlceLiqFlux  ! error control
  ! local variables
  integer(i4b)                      :: iLayer                     ! layer index
  ! ------------------------------------------------------------------------------------------------------------------------------------------
  associate(&
    ! input: model control
    nLayers                   => in_snowLakeGlceLiqFlux % nLayers,                        & ! intent(in):  number of liquid lake layers
    nStart                    => in_snowLakeGlceLiqFlux % nStart,                         & ! intent(in):  index of the layer above the liquid lake layers (nSnow + nLakeFrz)
    surface_flux              => in_snowLakeGlceLiqFlux % surface_flux,                   & ! intent(in):  water arriving at the top of the liquid lake layers (m s-1)
    bottom_flux               => in_snowLakeGlceLiqFlux % bottom_flux,                    & ! intent(in):  flux at the bottom of the lake if already computed (m s-1)
    ! input-output: interface fluxes
    iLayerLiqFluxSnLaGl       => io_snowLakeGlceLiqFlux % iLayerLiqFluxSnLaGl,            & ! intent(inout): [dp(0:)] liquid flux at snow, lake, glce layer interfaces (m s-1)
    iLayerLiqFluxSnLaGlDeriv  => io_snowLakeGlceLiqFlux % iLayerLiqFluxSnLaGlDeriv,       & ! intent(inout): [dp(0:)] derivative in the interface flux w.r.t. the layer above
    ! stream domain
    scalarGroundEvaporation   => flux_data%var(iLookFLUX%scalarGroundEvaporation)%dat(1),    & ! intent(in):  [dp] evaporation from the open water surface (kg m-2 s-1)
    scalarGroundSublimation   => flux_data%var(iLookFLUX%scalarGroundSublimation)%dat(1),    & ! intent(in):  [dp] sublimation from the water surface once it froze (kg m-2 s-1)
    scalarStreamSfcInflow     => flux_data%var(iLookFLUX%scalarStreamSfcInflow)%dat(1),      & ! intent(out): [dp] rain plus melt entering the open water column (m s-1)
    scalarStreamRunoff        => flux_data%var(iLookFLUX%scalarStreamRunoff)%dat(1),         & ! intent(out): [dp] net water the stream domain adds to the reach (m s-1)
    scalarStreamSfcInflowTemp => diag_data%var(iLookDIAG%scalarStreamSfcInflowTemp)%dat(1),  & ! intent(out): [dp] temperature of the rain plus melt entering the open water column (K)
    ! output: error control
    err                       => out_snowLakeGlceLiqFlux % err,                           & ! intent(out): error code
    message                   => out_snowLakeGlceLiqFlux % cmessage                       ) ! intent(out): error message
    ! ------------------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='lakeLiqFlux/'

    ! NOTE: the domain types are module variables rather than parameters, so no select case
    if(domType==stream)then
      ! ***** stream: the reach water column, mass handled by the river network
      ! what arrives at the top joins the flow at once rather than entering the column: no interface flux carries it down
      ! (the interface above, index nStart, keeps the drainage the layer above computed, since that layer needs its bottom flux)
      do iLayer=nStart+1,nStart+nLayers
        iLayerLiqFluxSnLaGl(iLayer)      = 0._rkind
        iLayerLiqFluxSnLaGlDeriv(iLayer) = 0._rkind
      end do
      iLayerLiqFluxSnLaGl(nStart+nLayers) = bottom_flux ! zero for now: no seepage through the bed
      ! rain and melt join the flow: through the water column when it is open, straight to the reach when it is ice covered
      if(topLayerTemp > Tfreeze)then
        scalarStreamSfcInflow     = surface_flux
        scalarStreamSfcInflowTemp = surfaceFluxTemp
      else
        scalarStreamSfcInflow     = 0._rkind
        scalarStreamSfcInflowTemp = Tfreeze
      end if
      ! net water the stream domain itself hands to the reach (evaporation is negative water; so is sublimation from a
      ! surface that froze within the step, when nothing lies on the water: with snow or an ice cover on top the layer
      ! above loses it, and the ice branch sets the runoff)
      scalarStreamRunoff = surface_flux - scalarGroundEvaporation/iden_water
      if(nStart==0) scalarStreamRunoff = scalarStreamRunoff - scalarGroundSublimation/iden_water

    else if(domType==wetland)then
      ! ***** wetland: sub-GRU lake or pothole, mass balance solved here
      ! Intended interface: surface_flux is the inflow at the top (rain, melt, and the upland runoff of the same HRU);
      ! the interface fluxes within the lake are zero (well mixed); the bottom flux is the seepage to the soil,
      ! and spill over the outlet leaves from the top layer; the area then follows the stored volume
      ! (updateLakeArea in run_oneGRU).
      err=20; message=trim(message)//'wetland lake mass balance not implemented'; return

    else
      err=20; message=trim(message)//'lake layers are only expected in stream or wetland domains'; return
    end if

  end associate

end subroutine lakeLiqFlux

end module snowLakeGlceLiqFlux_module
