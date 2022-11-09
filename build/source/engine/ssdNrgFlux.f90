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

module ssdNrgFlux_module

! data types
USE nrtype

! data types
USE data_types,only:var_d           ! x%var(:)       (dp)
USE data_types,only:var_dlength     ! x%var(:)%dat   (dp)
USE data_types,only:var_ilength     ! x%var(:)%dat   (i4b)

! physical constants
USE multiconst,only:&
                    sb,          & ! Stefan Boltzman constant      (W m-2 K-4)
                    Em_Sno,      & ! emissivity of snow            (-)
                    LH_fus,      & ! latent heat of fusion         (J kg-1)
                    LH_vap,      & ! latent heat of vaporization   (J kg-1)
                    LH_sub,      & ! latent heat of sublimation    (J kg-1)
                    gravity,     & ! gravitational acceleteration  (m s-2)
                    Tfreeze,     & ! freezing point of pure water  (K)
                    iden_air,    & ! intrinsic density of air      (kg m-3)
                    iden_ice,    & ! intrinsic density of ice      (kg m-3)
                    iden_water,  &  ! intrinsic density of water    (kg m-3)
                    ! specific heat
                    Cp_air,      & ! specific heat of air          (J kg-1 K-1)
                    Cp_water,    & ! specific heat of liquid water (J kg-1 K-1)
                    ! thermal conductivity
                    lambda_air,  & ! thermal conductivity of air   (J s-1 m-1)
                    lambda_ice,  & ! thermal conductivity of ice   (J s-1 m-1)
                    lambda_water   ! thermal conductivity of water (J s-1 m-1)


! missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! named variables for snow and soil
USE globalData,only:iname_snow     ! named variables for snow
USE globalData,only:iname_soil     ! named variables for soil

! named variables
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements

! model decisions
USE globalData,only:model_decisions                         ! model decision structure
USE var_lookup,only:iLookDECISIONS                          ! named variables for elements of the decision structure

! provide access to look-up values for model decisions
USE mDecisions_module,only:      &
 ! look-up values for method used to compute derivative
 numerical,                      & ! numerical solution
 analytical,                     & ! analytical solution
 ! look-up values for choice of boundary conditions for thermodynamics
 prescribedTemp,                 & ! prescribed temperature
 energyFlux,                     & ! energy flux
 zeroFlux,                       & ! zero flux
 ! look-up values for choice of boundary conditions for soil hydrology
 prescribedHead,                 & ! prescribed head
 ! look-up values for choice of thermal conductivity representation for snow
 Yen1965,                        & ! Yen (1965)
 Mellor1977,                     & ! Mellor (1977)
 Jordan1991,                     & ! Jordan (1991)
 Smirnova2000,                   & ! Smirnova et al. (2000)
 ! look-up values for choice of thermal conductivity representation for soil
 funcSoilWet,                    & ! function of soil wetness
 mixConstit,                     & ! mixture of constituents
 hanssonVZJ,                     & ! test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004
 ! look-up values for the form of Richards' equation
 moisture,                       & ! moisture-based form of Richards' equation
 mixdform                          ! mixed form of Richards' equation

! -------------------------------------------------------------------------------------------------
implicit none
private
public::ssdNrgFlux
! global parameters
real(rkind),parameter            :: dx=1.e-10_rkind             ! finite difference increment (K)
contains

 ! ************************************************************************************************
! public subroutine ssdNrgFlux: compute energy fluxes and derivatives at layer interfaces
 ! ************************************************************************************************
subroutine ssdNrgFlux(&
                      ! input: model control
                      scalarSolution,                     & ! intent(in):    flag to indicate the scalar solution
                      deriv_desired,                & ! intent(in): flag indicating if derivatives are desired
                      ! input: fluxes and derivatives at the upper boundary
                      groundNetFlux,                      & ! intent(in):    total flux at the ground surface (W m-2)
                      dGroundNetFlux_dGroundTemp,         & ! intent(in):    derivative in total ground surface flux w.r.t. ground temperature (W m-2 K-1)
                      ! input: liquid water fluxes
                      iLayerLiqFluxSnow,                  & ! intent(in):    liquid flux at the interface of each snow layer (m s-1)
                      iLayerLiqFluxSoil,                  & ! intent(in):    liquid flux at the interface of each soil layer (m s-1)
                      ! input: trial value of model state variables
                      mLayerTempTrial,                    & ! intent(in):    trial temperature at the current iteration (K)
                      mLayerMatricHeadTrial,              & ! intent(in):    trial matric head at the current iteration(m)
                      mLayerVolFracLiqTrial,              & ! intent(in):    trial volumetric fraction of liquid water at the current iteration(-)
                      mLayerVolFracIceTrial,              & ! intent(in):    trial volumetric fraction of ice water at the current iteration(-)
                      ! input: pre-computed derivatives
                      mLayerdTheta_dTk,                   & ! intent(in):    derivative in volumetric liquid water content w.r.t. temperature (K-1)
                      mLayerFracLiqSnow,                  & ! intent(in):    fraction of liquid water (-)
                      ! input-output: data structures
                      mpar_data,                          & ! intent(in):    model parameters
                      indx_data,                          & ! intent(in):    model indices
                      prog_data,                          & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,                          & ! intent(inout):    model diagnostic variables for a local HRU
                      flux_data,                          & ! intent(inout): model fluxes for a local HRU
                      ! output: fluxes and derivatives at all layer interfaces
                      iLayerNrgFlux,                      & ! intent(out):   energy flux at the layer interfaces (W m-2)
                      dFlux_dTempAbove,                   & ! intent(out):   derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                      dFlux_dTempBelow,                   & ! intent(out):   derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                    ! output: error control
                      err,message)                          ! intent(out): error control
  ! utility modules
  USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water
  USE soil_utils_module,only:matricHead     ! compute the matric head based on volumetric water content
  USE soil_utils_module,only:crit_soilT     ! compute critical temperature below which ice exists
  USE snow_utils_module,only:fracliquid     ! compute fraction of liquid water at a given temperature
  USE soil_utils_module,only:dTheta_dPsi    ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
  USE soil_utils_module,only:dPsi_dTheta    ! compute derivative of the soil moisture characteristic w.r.t. theta (m)

  ! constants
  USE multiconst, only: gravity, &                          ! gravitational acceleration (m s-1)
                        Tfreeze, &                          ! freezing point of water (K)
                        iden_water,iden_ice,&      ! intrinsic density of water and ice (kg m-3)
                        LH_fus                              ! latent heat of fusion (J kg-1)
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! input: model control
  logical(lgt),intent(in)            :: scalarSolution              ! flag to denote if implementing the scalar solution
  logical(lgt),intent(in)            :: deriv_desired               ! flag indicating if derivatives are desired
  ! input: fluxes and derivatives at the upper boundary
  real(rkind),intent(in)             :: groundNetFlux               ! net energy flux for the ground surface (W m-2)
  real(rkind),intent(inout)          :: dGroundNetFlux_dGroundTemp  ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  ! input: liquid water fluxes
  real(rkind),intent(in)             :: iLayerLiqFluxSnow(0:)       ! liquid flux at the interface of each snow layer (m s-1)
  real(rkind),intent(in)             :: iLayerLiqFluxSoil(0:)       ! liquid flux at the interface of each soil layer (m s-1)
  ! input: trial model state variables
  real(rkind),intent(in)              :: mLayerTempTrial(:)         ! temperature in each layer at the current iteration (m)
  real(rkind),intent(in)              :: mLayerMatricHeadTrial(:)   ! matric head in each layer at the current iteration (m)
  real(rkind),intent(in)              :: mLayerVolFracLiqTrial(:)   ! volumetric fraction of liquid at the current iteration (-)
  real(rkind),intent(in)              :: mLayerVolFracIceTrial(:)   ! volumetric fraction of ice at the current iteration (-)
  ! input: pre-computed derivatives
  real(rkind),intent(in)              :: mLayerdTheta_dTk(:)        ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in)              :: mLayerFracLiqSnow(:)       ! fraction of liquid water (-)
  ! input-output: data structures
  type(var_dlength),intent(in)        :: mpar_data                  ! model parameters
  type(var_ilength),intent(in)        :: indx_data                  ! state vector geometry
  type(var_dlength),intent(in)        :: prog_data                  ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)        :: diag_data                  ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout)     :: flux_data                  ! model fluxes for a local HRU
  ! output: fluxes and derivatives at all layer interfaces
  real(rkind),intent(out)             :: iLayerNrgFlux(0:)          ! energy flux at the layer interfaces (W m-2)
  real(rkind),intent(out)             :: dFlux_dTempAbove(0:)       ! derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
  real(rkind),intent(out)             :: dFlux_dTempBelow(0:)       ! derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
  ! output: error control
  integer(i4b),intent(out)            :: err                        ! error code
  character(*),intent(out)            :: message                    ! error message
  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  character(LEN=256)               :: cmessage                     ! error message of downwind routine
  integer(i4b)                     :: i,j,iLayer                   ! index of model layers
  integer(i4b)                     :: ixLayerDesired(1)            ! layer desired (scalar solution)
  integer(i4b)                     :: ixTop                        ! top layer in subroutine call
  integer(i4b)                     :: ixBot                        ! bottom layer in subroutine call
  real(rkind)                      :: qFlux                        ! liquid flux at layer interfaces (m s-1)
  real(rkind)                      :: dz                           ! height difference (m)
  ! additional variables to compute numerical derivatives
  integer(i4b)                     :: nFlux                        ! number of flux calculations required (>1 = numerical derivatives with one-sided finite differences)
  integer(i4b)                     :: itry                         ! index of different flux calculations
  integer(i4b),parameter           :: unperturbed=0                ! named variable to identify the case of unperturbed state variables
  integer(i4b),parameter           :: perturbState=1               ! named variable to identify the case where we perturb the state in the current layer
  integer(i4b),parameter           :: perturbStateTempAbove=2      ! named variable to identify the case where we perturb the state layer above
  integer(i4b),parameter           :: perturbStateTempBelow=3      ! named variable to identify the case where we perturb the state layer below
  integer(i4b),parameter           :: perturbStateWatAbove=4       ! named variable to identify the case where we perturb the state layer above
  integer(i4b),parameter           :: perturbStateWatBelow=5       ! named variable to identify the case where we perturb the state layer below
  integer(i4b)                     :: ixPerturb                    ! index of element in 2-element vector to perturb
  integer(i4b)                     :: ixOriginal                   ! index of perturbed element in the original vector
  real(rkind)                      :: scalarThermCFlux               ! thermal conductivity (W m-1 K-1)
  real(rkind)                      :: scalarThermCFlux_dTempAbove    ! thermal conductivity with perturbation to the temperature state above (W m-1 K-1)
  real(rkind)                      :: scalarThermCFlux_dTempBelow    ! thermal conductivity with perturbation to the temperature state below (W m-1 K-1)
  real(rkind)                      :: scalarThermCFlux_dWatAbove     ! thermal conductivity with perturbation to the water state above
  real(rkind)                      :: scalarThermCFlux_dWatBelow     ! thermal conductivity with perturbation to the water state below
  real(rkind)                      :: flux0,flux1,flux2            ! fluxes used to calculate derivatives (W m-2)
  ! compute fluxes and derivatives at layer interfaces
  integer(i4b),dimension(2)        :: mLayer_ind                   ! indices of above and below layers
  integer(i4b),dimension(2)        :: iLayer_ind                   ! indices of above and below interfaces
  real(rkind)                      :: matricFHead                  ! matric head for frozen soil
  real(rkind)                      :: Tcrit                        ! temperature where all water is unfrozen (K)
  real(rkind)                      :: fLiq                         ! fraction of liquid water (-)
  real(rkind),dimension(2)         :: vectorMatricHeadTrial        ! trial value of matric head (m)
  real(rkind),dimension(2)         :: vectorVolFracLiqTrial        ! trial value of volumetric liquid content (-)
  real(rkind),dimension(2)         :: vectorVolFracIceTrial        ! trial value of volumetric ice content (-)
  real(rkind),dimension(2)         :: vectorTempTrial              ! trial value of temperature (K)
  real(rkind),dimension(2)         :: vectordTheta_dPsi            ! derivative in the soil water characteristic w.r.t. psi (m-1)
  real(rkind),dimension(2)         :: vectordPsi_dTheta            ! derivative in the soil water characteristic w.r.t. theta (m)
  real(rkind),dimension(2)         :: vectorFracLiqSnow            ! fraction of liquid water (-)
  real(rkind),dimension(2)         :: vectortheta_sat              ! layer above and below soil porosity (-)
  real(rkind),dimension(2)         :: vectoriden_soil              ! layer above and below density of soil (kg m-3)
  real(rkind),dimension(2)         :: vectorthCond_soil            ! layer above and below thermal conductivity of soil (W m-1 K-1)
  real(rkind),dimension(2)         :: vectorfrac_sand              ! layer above and below fraction of sand (-)
  real(rkind),dimension(2)         :: vectorfrac_clay              ! layer above and below fraction of clay (-)
  ! recompute the perturbed version of iLayerThermalC, this could be the only version and remove the computThermConduct_module
  real(rkind)                      :: dThermalC_dHydStateAbove     ! derivative in the thermal conductivity w.r.t. water state in the layer above
  real(rkind)                      :: dThermalC_dHydStateBelow     ! derivative in the thermal conductivity w.r.t. water state in the layer above
  real(rkind)                      :: dThermalC_dNrgStateAbove     ! derivative in the thermal conductivity w.r.t. energy state in the layer above
  real(rkind)                      :: dThermalC_dNrgStateBelow     ! derivative in the thermal conductivity w.r.t. energy state in the layer above
  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  ! make association of local variables with information in the data structures
  associate(&
    ixDerivMethod           => model_decisions(iLookDECISIONS%fDerivMeth)%iDecision,      & ! intent(in): method used to calculate flux derivatives
    ix_bcUpprTdyn           => model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision,      & ! intent(in): method used to calculate the upper boundary condition for thermodynamics
    ix_bcLowrTdyn           => model_decisions(iLookDECISIONS%bcLowrTdyn)%iDecision,      & ! intent(in): method used to calculate the lower boundary condition for thermodynamics
    ixRichards                 => model_decisions(iLookDECISIONS%f_Richards)%iDecision,   & ! intent(in): index of the form of Richards' equation
    ixThCondSnow            => model_decisions(iLookDECISIONS%thCondSnow)%iDecision,      & ! intent(in): choice of method for thermal conductivity of snow
    ixThCondSoil            => model_decisions(iLookDECISIONS%thCondSoil)%iDecision,      & ! intent(in): choice of method for thermal conductivity of soil
    ! input: model coordinates
    nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1),               & ! intent(in): number of snow layers
    nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1),             & ! intent(in): total number of layers
    layerType               => indx_data%var(iLookINDEX%layerType)%dat,              & ! intent(in): layer type (iname_soil or iname_snow)
    ixLayerState            => indx_data%var(iLookINDEX%ixLayerState)%dat,           & ! intent(in): list of indices for all model layers
    ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat,          & ! intent(in): index in the state subset for energy state variables in the snow+soil domain
    ! input: thermal properties
    mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat,             & ! intent(in): depth of each layer (m)
    mLayerHeight            => prog_data%var(iLookPROG%mLayerHeight)%dat,            & ! intent(in): height at the mid-point of each layer (m)
    upperBoundTemp          => mpar_data%var(iLookPARAM%upperBoundTemp)%dat(1),      & ! intent(in): temperature of the upper boundary (K)
    lowerBoundTemp          => mpar_data%var(iLookPARAM%lowerBoundTemp)%dat(1),      & ! intent(in): temperature of the lower boundary (K)
    iLayerHeight            => prog_data%var(iLookPROG%iLayerHeight)%dat,            & ! intent(in): height at the interface of each layer (m)
    fixedThermalCond_snow   => mpar_data%var(iLookPARAM%fixedThermalCond_snow)%dat(1),    & ! intent(in): temporally constant thermal conductivity of snow (W m-1 K-1)
    iLayerThermalC          => diag_data%var(iLookDIAG%iLayerThermalC)%dat,          & ! intent(inout): thermal conductivity at the interface of each layer (W m-1 K-1)
    ! input: depth varying soil parameters
    iden_soil               => mpar_data%var(iLookPARAM%soil_dens_intr)%dat,              & ! intent(in): intrinsic density of soil (kg m-3)
    thCond_soil             => mpar_data%var(iLookPARAM%thCond_soil)%dat,                 & ! intent(in): thermal conductivity of soil (W m-1 K-1)
    theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat,                   & ! intent(in): soil porosity (-)
    frac_sand               => mpar_data%var(iLookPARAM%frac_sand)%dat,                   & ! intent(in): fraction of sand (-)
    frac_clay               => mpar_data%var(iLookPARAM%frac_clay)%dat,                   & ! intent(in): fraction of clay (-)
    vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat,                  & ! intent(in):  [dp(:)] van Genutchen "m" parameter (-)
    vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)%dat,                       & ! intent(in):  [dp(:)] van Genutchen "n" parameter (-)
    vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)%dat,                   & ! intent(in):  [dp(:)] van Genutchen "alpha" parameter (m-1)
    theta_res               => mpar_data%var(iLookPARAM%theta_res)%dat,                   & ! intent(in):  [dp(:)] soil residual volumetric water content (-)
    ! input: snow parameters
    snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1),            & ! intent(in):  [dp] scaling parameter for the snow freezing curve (K-1)
    ! output: diagnostic fluxes
    iLayerConductiveFlux => flux_data%var(iLookFLUX%iLayerConductiveFlux)%dat,    & ! intent(out): conductive energy flux at layer interfaces at end of time step (W m-2)
    iLayerAdvectiveFlux  => flux_data%var(iLookFLUX%iLayerAdvectiveFlux)%dat      & ! intent(out): advective energy flux at layer interfaces at end of time step (W m-2)
    )  ! association of local variables with information in the data structures
    ! ------------------------------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='ssdNrgFlux/'

    ! set conductive and advective fluxes to missing in the upper boundary
    ! NOTE: advective flux at the upper boundary is included in the ground heat flux
    iLayerConductiveFlux(0) = realMissing
    iLayerAdvectiveFlux(0)  = realMissing

    ! check the need to compute numerical derivatives
    if(ixDerivMethod==numerical)then
      nFlux=5  ! compute the derivatives and cross derivates using one-sided finite differences
    else
      nFlux=0  ! compute analytical derivatives
    end if

    ! get the indices for the snow+soil layers
    if(scalarSolution)then
      ixLayerDesired = pack(ixLayerState, ixSnowSoilNrg/=integerMissing)
      ixTop = ixLayerDesired(1)
      ixBot = ixLayerDesired(1)
    else
  ixTop = 1
      ixBot = nLayers
    endif

    ! -------------------------------------------------------------------------------------------------------------------------
    ! ***** compute the conductive fluxes at layer interfaces *****
    ! -------------------------------------------------------------------------------------------------------------------------
    do iLayer=ixTop,ixBot ! (loop through model layers)

  ! compute fluxes at the lower boundary -- positive downwards
  if(iLayer==nLayers)then
      ! flux depends on the type of lower boundary condition
      select case(ix_bcLowrTdyn) ! (identify the lower boundary condition for thermodynamics
    case(prescribedTemp); iLayerConductiveFlux(nLayers) = -iLayerThermalC(iLayer)*(lowerBoundTemp - mLayerTempTrial(iLayer))/(mLayerDepth(iLayer)*0.5_rkind)
    case(zeroFlux);       iLayerConductiveFlux(nLayers) = 0._rkind
    case default;         err=20; message=trim(message)//'unable to identify lower boundary condition for thermodynamics'; return
      end select  ! (identifying the lower boundary condition for thermodynamics)

  ! compute fluxes within the domain -- positive downwards
  else
        iLayerConductiveFlux(iLayer)  = -iLayerThermalC(iLayer)*(mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer)) / &
                                        (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))

        !write(*,'(a,i4,1x,2(f9.3,1x))') 'iLayer, iLayerConductiveFlux(iLayer), iLayerThermalC(iLayer) = ', iLayer, iLayerConductiveFlux(iLayer), iLayerThermalC(iLayer)
      end if ! (the type of layer)
 end do

    ! -------------------------------------------------------------------------------------------------------------------------
    ! ***** compute the advective fluxes at layer interfaces *****
    ! -------------------------------------------------------------------------------------------------------------------------
 do iLayer=ixTop,ixBot
  ! get the liquid flux at layer interfaces
        select case(layerType(iLayer))
          case(iname_snow); qFlux = iLayerLiqFluxSnow(iLayer)
          case(iname_soil); qFlux = iLayerLiqFluxSoil(iLayer-nSnow)
          case default; err=20; message=trim(message)//'unable to identify layer type'; return
        end select
        ! compute fluxes at the lower boundary -- positive downwards
        if(iLayer==nLayers)then
          iLayerAdvectiveFlux(iLayer) = -Cp_water*iden_water*qFlux*(lowerBoundTemp - mLayerTempTrial(iLayer))
        ! compute fluxes within the domain -- positive downwards
        else
          iLayerAdvectiveFlux(iLayer) = -Cp_water*iden_water*qFlux*(mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer))
        end if
    end do  ! looping through layers

    ! -------------------------------------------------------------------------------------------------------------------------
    ! ***** compute the total fluxes at layer interfaces *****
    ! -------------------------------------------------------------------------------------------------------------------------
    ! NOTE: ignore advective fluxes for now
 iLayerNrgFlux(0)           = groundNetFlux
    iLayerNrgFlux(ixTop:ixBot) = iLayerConductiveFlux(ixTop:ixBot)
 !print*, 'iLayerNrgFlux(0:4) = ', iLayerNrgFlux(0:4)

 ! -------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the derivative in fluxes at layer interfaces w.r.t temperature in the layer above and the layer below *****
 ! -------------------------------------------------------------------------------------------------------------------------

 ! initialize un-used elements
 dFlux_dTempBelow(nLayers) = -huge(lowerBoundTemp)  ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems

 ! ***** the upper boundary
 dFlux_dTempBelow(0) = dGroundNetFlux_dGroundTemp

 ! loop through INTERFACES...
 do iLayer=ixTop,ixBot

  ! ***** the lower boundary
  if(iLayer==nLayers)then  ! (lower boundary)

   ! identify the lower boundary condition
   select case(ix_bcLowrTdyn)

    ! * prescribed temperature at the lower boundary
    case(prescribedTemp)

     dz = mLayerDepth(iLayer)*0.5_rkind
     if(ixDerivMethod==analytical)then    ! ** analytical derivatives
      dFlux_dTempAbove(iLayer) = iLayerThermalC(iLayer)/dz
     else                              ! ** numerical derivatives
      flux0 = -iLayerThermalC(iLayer)*(lowerBoundTemp - (mLayerTempTrial(iLayer)   ))/dz
      flux1 = -iLayerThermalC(iLayer)*(lowerBoundTemp - (mLayerTempTrial(iLayer)+dx))/dz
      dFlux_dTempAbove(iLayer) = (flux1 - flux0)/dx
    end if

     ! * zero flux at the lower boundary
     case(zeroFlux)
      dFlux_dTempAbove(iLayer) = 0._rkind

     case default; err=20; message=trim(message)//'unable to identify lower boundary condition for thermodynamics'; return

   end select  ! (identifying the lower boundary condition for thermodynamics)

  ! ***** internal layers
            else
   dz = (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
   if(ixDerivMethod==analytical)then    ! ** analytical derivatives
    dFlux_dTempAbove(iLayer) =  iLayerThermalC(iLayer)/dz
    dFlux_dTempBelow(iLayer) = -iLayerThermalC(iLayer)/dz
   else                              ! ** numerical derivatives
    flux0 = -iLayerThermalC(iLayer)*( mLayerTempTrial(iLayer+1)     -  mLayerTempTrial(iLayer)    ) / dz
    flux1 = -iLayerThermalC(iLayer)*( mLayerTempTrial(iLayer+1)     - (mLayerTempTrial(iLayer)+dx)) / dz
    flux2 = -iLayerThermalC(iLayer)*((mLayerTempTrial(iLayer+1)+dx) -  mLayerTempTrial(iLayer)    ) / dz
    dFlux_dTempAbove(iLayer) = (flux1 - flux0)/dx
    dFlux_dTempBelow(iLayer) = (flux2 - flux0)/dx
        end if

  end if  ! type of layer (upper, internal, or lower)

 end do  ! (looping through layers)

 ! end association of local variables with information in the data structures
 end associate

 end subroutine ssdNrgFlux

end module ssdNrgFlux_module

