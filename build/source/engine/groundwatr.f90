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

module groundwatr_module

! data types
USE nr_type

! model constants
USE multiconst,only:iden_water   ! density of water (kg m-3)

! derived types to define the data structures
USE data_types,only:&
                    var_d,              & ! data vector (rkind)
                    var_dlength,        & ! data vector with variable length dimension (rkind)
                    in_type_groundwatr, & ! intent(in) arguments for groundwatr call
                    io_type_groundwatr, & ! intent(inout) arguments for groundwatr call
                    out_type_groundwatr   ! intent(out) arguments for groundwatr call

! look-up values for the form of Richards' equation
USE mDecisions_module,only:       &
 moisture,                        & ! moisture-based form of Richards' equation
 mixdform                           ! mixed form of Richards' equation

! named variables defining elements in the data structures
USE var_lookup,only:iLookATTR    ! named variables for structure elements
USE var_lookup,only:iLookPROG    ! named variables for structure elements
USE var_lookup,only:iLookDIAG    ! named variables for structure elements
USE var_lookup,only:iLookFLUX    ! named variables for structure elements
USE var_lookup,only:iLookPARAM   ! named variables for structure elements

! privacy
implicit none
private
public :: groundwatr
contains

! ************************************************************************************************
! public subroutine groundwatr: compute the groundwater sink term in Richards' equation
! ************************************************************************************************
!
! Method
! ------
!
! Here we assume that water available for shallow groundwater flow includes is all water above
! "field capacity" below the depth zCrit, where zCrit is defined as the lowest point in the soil
! profile where the volumetric liquid water content is less than field capacity.
!
! We further assume that transmssivity (m2 s-1) for each layer is defined assuming that the water
! available for saturated flow is located at the bottom of the soil profile. Specifically:
!  trTotal(iLayer) = tran0*(zActive(iLayer)/soilDepth)**zScale_TOPMODEL
!  trSoil(iLayer)  = trTotal(iLayer) - trTotal(iLayer+1)
! where zActive(iLayer) is the effective water table thickness for all layers up to and including
! the current layer (working from the bottom to the top).
!
! The outflow from each layer is then (m3 s-1)
!  mLayerOutflow(iLayer) = trSoil(iLayer)*tan_slope*contourLength
! where contourLength is the width of a hillslope (m) parallel to a stream
!
! ************************************************************************************************
subroutine groundwatr(&
                      ! input: model control, state variables, and diagnostic variables
                      in_groundwatr,                          & ! intent(in): model control, state variables, and diagnostic variables
                      ! input/output: data structures
                      attr_data,                              & ! intent(in):    spatial attributes
                      mpar_data,                              & ! intent(in):    model parameters
                      prog_data,                              & ! intent(in):    model prognostic variables for a local HRU
                      flux_data,                              & ! intent(inout): model fluxes for a local HRU
                      ! input-output: baseflow
                      io_groundwatr,                          & ! intent(inout): index of lowest saturated layer (NOTE: only computed on the first iteration)
                      ! output: baseflow and error control
                      out_groundwatr)                           ! intent(out):   baseflow and error control
  ! ---------------------------------------------------------------------------------------
  implicit none
  ! input: model control, state variables, and diagnostic variables
  type(in_type_groundwatr),intent(in)    :: in_groundwatr     ! model control, state variables, and diagnostic variables
  ! input-output: data structures
  type(var_d),intent(in)                 :: attr_data         ! spatial attributes
  type(var_dlength),intent(in)           :: mpar_data         ! model parameters
  type(var_dlength),intent(in)           :: prog_data         ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)        :: flux_data         ! model fluxes for a local HRU
  ! input-output: baseflow
  type(io_type_groundwatr),intent(inout) :: io_groundwatr     ! index of lowest saturated layer (NOTE: only computed on the first iteration)
  ! output: baseflow and error control
  type(out_type_groundwatr),intent(out)  :: out_groundwatr    ! baseflow and error control
  ! general local variables
  integer(i4b)                           :: iLayer            ! index of soil layer
  character(len=256)                     :: cmessage          ! error message
  ! ***************************************************************************************
  ! associate variables in data structures
  allocate(out_groundwatr % mLayerBaseflow(in_groundwatr%nSoil),out_groundwatr % dBaseflow_dWat(in_groundwatr%nSoil,in_groundwatr%nSoil),&
           out_groundwatr % dBaseflow_dTk(in_groundwatr%nSoil,in_groundwatr%nSoil)) ! allocate intent(out) data structure components
  associate(&
    ! input: model control
    nSnow               => in_groundwatr % nSnow,                              & ! intent(in):    [i4b] number of snow layers
    nSoil               => in_groundwatr % nSoil,                              & ! intent(in):    [i4b] number of soil layers
    nLayers             => in_groundwatr % nLayers,                            & ! intent(in):    [i4b] total number of layers
    getSatDepth         => in_groundwatr % firstFluxCall,                      & ! intent(in):    [lgt] logical flag to compute index of the lowest saturated layer
    ixRichards          => in_groundwatr % ixRichards,                         & ! intent(in):    [i4b] index of the form of Richards' equation
    ! input: diagnostic variables
    mLayerVolFracLiq    => in_groundwatr % mLayerVolFracLiqTrial,              & ! intent(in):    [dp] volumetric fraction of liquid water (-)
    mLayerVolFracIce    => in_groundwatr % mLayerVolFracIceTrial,              & ! intent(in):    [dp] volumetric fraction of ice (-)
    ! input: derivatives
    dVolTot_dPsi0       => in_groundwatr % dVolTot_dPsi0,                      & ! intent(in):    [dp] derivative in total volumetric water content w.r.t. matric head (m-1)
    mLayerdTheta_dPsi   => in_groundwatr % mLayerdTheta_dPsi,                  & ! intent(in):    [dp] derivative in liquid water content w.r.t. matric potential (m-1)
    mLayerdTheta_dTk    => in_groundwatr % mLayerdTheta_dTk,                   & ! intent(in):    [dp] derivative in volumetric liquid water content w.r.t. temperature (K-1)
    ! input: baseflow parameters
    fieldCapacity       => mpar_data%var(iLookPARAM%fieldCapacity)%dat(1),     & ! intent(in):    [dp] field capacity (-)
    theta_sat           => mpar_data%var(iLookPARAM%theta_sat)%dat,            & ! intent(in):    [dp] soil porosity (-)
    ! input-output: baseflow
    ixSaturation        => io_groundwatr % ixSaturation,                       & ! intent(inout): [i4b] index of lowest saturated layer (NOTE: only computed on the first iteration)
    ! output: diagnostic variables
    scalarExfiltration  => flux_data%var(iLookFLUX%scalarExfiltration)%dat(1), & ! intent(out):   [dp]    exfiltration from the soil profile (m s-1)
    mLayerColumnOutflow => flux_data%var(iLookFLUX%mLayerColumnOutflow)%dat,   & ! intent(out):   [dp(:)] column outflow from each soil layer (m3 s-1)
    ! output: baseflow
    mLayerBaseflow      => out_groundwatr % mLayerBaseflow,                    & ! intent(out):   [dp(:)]   baseflow from each soil layer (m s-1)
    dBaseflow_dWat      => out_groundwatr % dBaseflow_dWat,                    & ! intent(out):   [dp(:,:)] derivative in baseflow w.r.t. soil water characteristic
    dBaseflow_dTk       => out_groundwatr % dBaseflow_dTk,                     & ! intent(out):   [dp(:,:)] derivative in baseflow w.r.t. temperature (m s-1 K-1)
    ! output: error control
    err                 => out_groundwatr % err,                               & ! intent(out):   [i4b]       error code
    message             => out_groundwatr % cmessage                           & ! intent(out):   [character] error message
    )  ! end association to variables in data structures

    ! initialize error control
    err=0; message='groundwatr/'

    ! ************************************************************************************************
    ! (1) compute the "active" portion of the soil profile
    ! ************************************************************************************************

   ! get index of the layer closest to surface that is more than field capacity (NOTE: only compute on the first flux call)
    if (getSatDepth) then
      ixSaturation = nSoil+1  ! unsaturated profile when ixSaturation>nSoil
      do iLayer=nSoil,1,-1  ! start at the lowest soil layer and work upwards to the top layer
        if (mLayerVolFracLiq(iLayer) > fieldCapacity) then; ixSaturation = iLayer  ! index of saturated layer -- keeps getting over-written as move upwards
        else; exit; end if                                                         ! only consider saturated layer at the bottom of the soil profile
      end do  ! end looping through soil layers
    end if

    ! check for an early return (no layers are "active")
    if (ixSaturation > nSoil) then
      scalarExfiltration     = 0._rkind   ! exfiltration from the soil profile (m s-1)
      mLayerColumnOutflow(:) = 0._rkind   ! column outflow from each soil layer (m3 s-1)
      mLayerBaseflow(:)      = 0._rkind   ! baseflow from each soil layer (m s-1)
      dBaseflow_dWat(:,:)    = 0._rkind   ! derivative in baseflow w.r.t. soil water characteristic
      dBaseflow_dTk(:,:)     = 0._rkind   ! derivative in baseflow w.r.t. temperature (m s-1 K-1)
      return
    end if  ! if some layers are saturated

    ! ************************************************************************************************
    ! (2) compute the baseflow flux and its derivative w.r.t matric head and temperature
    ! ************************************************************************************************

    ! use private subroutine to compute baseflow (for multiple calls for numerical Jacobian)
    call computBaseflow(&
                          ! input: control and state variables
                          nSnow,                   & ! intent(in):    number of snow layers
                          nSoil,                   & ! intent(in):    number of soil layers
                          nLayers,                 & ! intent(in):    total number of layers
                          ixSaturation,            & ! intent(in):    index of upper-most "saturated" layer
                          ixRichards,              & ! intent(in):    index of the form of Richards' equation
                          mLayerVolFracLiq,        & ! intent(in):    volumetric fraction of liquid water in each soil layer (-)
                          mLayerVolFracIce,        & ! intent(in):    volumetric fraction of ice in each soil layer (-)
                          ! input/output: data structures
                          attr_data,               & ! intent(in):    spatial attributes
                          mpar_data,               & ! intent(in):    model parameters
                          prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                          flux_data,               & ! intent(inout): model fluxes for a local HRU
                          ! derivatives
                          dVolTot_dPsi0,           & ! intent(in):    derivative in total volumetric water content w.r.t. matric head (m-1)
                          mLayerdTheta_dPsi,       & ! intent(in):    derivative in liquid water content w.r.t. matric potential (m-1)
                          mLayerdTheta_dTk,        & ! intent(in):    derivative in volumetric liquid water content w.r.t. temperature (K-1)
                          ! output: fluxes and derivatives
                          mLayerBaseflow,          & ! intent(out):   baseflow flux in each soil layer (m s-1)
                          dBaseflow_dWat,          & ! intent(out):   derivative in baseflow w.r.t. soil water characteristic
                          dBaseflow_dTk,           & ! intent(out):   derivative in baseflow w.r.t. temperature (m s-1 K-1)
                          err, cmessage)             ! intent(out):   error control
   if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
   ! end association to variables in data structures
  end associate

end subroutine groundwatr


! ***********************************************************************************************************************
! * private subroutine computBaseflow:  compute the baseflow flux and its derivative w.r.t. volumetric liquid water content
! ***********************************************************************************************************************
subroutine computBaseflow(&
                          ! input: control and state variables
                          nSnow,                         & ! intent(in):    number of snow layers
                          nSoil,                         & ! intent(in):    number of soil layers
                          nLayers,                       & ! intent(in):    total number of layers
                          ixSaturation,                  & ! intent(in):    index of upper-most "saturated" layer
                          ixRichards,                    & ! intent(in):    index of the form of Richards' equation
                          mLayerVolFracLiq,              & ! intent(in):    volumetric fraction of liquid water in each soil layer (-)
                          mLayerVolFracIce,              & ! intent(in):    volumetric fraction of ice in each soil layer (-)
                          ! input/output: data structures
                          attr_data,                     & ! intent(in):    spatial attributes
                          mpar_data,                     & ! intent(in):    model parameters
                          prog_data,                     & ! intent(in):    model prognostic variables for a local HRU
                          flux_data,                     & ! intent(inout): model fluxes for a local HRU
                          ! derivatives
                          dVolTot_dPsi0,                 & ! intent(in):    derivative in total volumetric water content w.r.t. matric head (m-1)
                          mLayerdTheta_dPsi,             & ! intent(in):    derivative in liquid water content w.r.t. matric potential (m-1)
                          mLayerdTheta_dTk,              & ! intent(in):    derivative in volumetric liquid water content w.r.t. temperature (K-1)
                          ! output: fluxes and derivatives
                          mLayerBaseflow,                & ! intent(out):   baseflow flux in each soil layer (m s-1)
                          dBaseflow_dWat,                & ! intent(out):   derivative in baseflow w.r.t. soil water characteristic
                          dBaseflow_dTk,                 & ! intent(out):   derivative in baseflow w.r.t. temperature (m s-1 K-1)
                          ! error handling
                          err, message)                    ! intent(out):   error control
  implicit none
  ! ---------------------------------------------------------------------------------------
  ! * dummy variables
  ! ---------------------------------------------------------------------------------------
  ! input: control and state variables
  integer(i4b),intent(in)          :: nSnow                   ! number of snow layers
  integer(i4b),intent(in)          :: nSoil                   ! number of soil layers
  integer(i4b),intent(in)          :: nLayers                 ! total number of layers
  integer(i4b),intent(in)          :: ixSaturation            ! index of upper-most "saturated" layer
  integer(i4b),intent(in)          :: ixRichards              ! index of the form of Richards' equation
  real(rkind),intent(in)           :: mLayerVolFracLiq(:)     ! volumetric fraction of liquid water (-)
  real(rkind),intent(in)           :: mLayerVolFracIce(:)     ! volumetric fraction of ice (-)
  ! input/output: data structures
  type(var_d),intent(in)           :: attr_data               ! spatial attributes
  type(var_dlength),intent(in)     :: mpar_data               ! model parameters
  type(var_dlength),intent(in)     :: prog_data               ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)  :: flux_data               ! model fluxes for a local HRU
  ! derivatives
  real(rkind),intent(in)           :: dVolTot_dPsi0(:)        ! derivative in total volumetric water content w.r.t. matric head (m-1)
  real(rkind),intent(in)           :: mLayerdTheta_dPsi(:)    ! derivative in liquid water content w.r.t. matric potential (m-1)
  real(rkind),intent(in)           :: mLayerdTheta_dTk(:)     ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
  ! output: baseflow
  real(rkind),intent(out)          :: mLayerBaseflow(:)       ! baseflow from each soil layer (m s-1)
  real(rkind),intent(out)          :: dBaseflow_dWat(:,:)     ! derivative in baseflow w.r.t. soil water characteristic
  real(rkind),intent(out)          :: dBaseflow_dTk(:,:)      ! derivative in baseflow w.r.t. temperature (m s-1 K-1)
  ! error handling
  integer(i4b), intent(out)        :: err                     ! error code
  character(*), intent(out)        :: message                 ! error message
  ! ---------------------------------------------------------------------------------------
  ! * local variables
  ! ---------------------------------------------------------------------------------------
  ! general local variables
  integer(i4b)                       :: iLayer,jLayer         ! index of model layer
  ! local variables for the exfiltration
  real(rkind)                        :: totalColumnInflow     ! total column inflow (m s-1)
  real(rkind)                        :: totalColumnOutflow    ! total column outflow (m s-1)
  real(rkind)                        :: availStorage          ! available storage (m)
  real(rkind),parameter              :: xMinEval=0.002_rkind  ! minimum value to evaluate the exfiltration function (m)
  real(rkind),parameter              :: xCenter=0.001_rkind   ! center of the exfiltration function (m)
  real(rkind),parameter              :: xWidth=0.0001_rkind   ! width of the exfiltration function (m)
  real(rkind)                        :: expF,logF             ! logistic smoothing function (-)
  ! local variables for the lateral flux among soil columns
  real(rkind)                        :: activePorosity        ! "active" porosity associated with storage above a threshold (-)
  real(rkind)                        :: drainableWater        ! drainable water in each layer (m)
  real(rkind)                        :: tran0                 ! maximum transmissivity (m2 s-1)
  real(rkind),dimension(nSoil)       :: zActive               ! water table thickness associated with storage below and including the given layer (m)
  real(rkind),dimension(nSoil)       :: trTotal               ! total transmissivity associated with total water table depth zActive (m2 s-1)
  real(rkind),dimension(nSoil)       :: trSoil                ! transmissivity of water in a given layer (m2 s-1)
  ! local variables for the derivatives
  real(rkind),dimension(nSoil,nSoil) :: dBaseflow_dVolLiq     ! derivative in baseflow w.r.t. volumetric liquid water content (s-1)
  real(rkind)                        :: qbTotal               ! total baseflow (m s-1)
  real(rkind)                        :: length2area           ! ratio of hillslope width to hillslope area (m m-2)
  real(rkind),dimension(nSoil)       :: depth2capacity        ! ratio of layer depth to total subsurface storage capacity (-)
  real(rkind),dimension(nSoil)       :: dXdS                  ! change in dimensionless flux w.r.t. change in dimensionless storage (-)
  real(rkind),dimension(nSoil)       :: dLogFunc_dWat         ! derivative in the logistic function w.r.t. soil water characteristic
  real(rkind),dimension(nSoil)       :: dExfiltrate_dWat      ! derivative in exfiltration w.r.t. soil water characteristic
  real(rkind),dimension(nSoil)       :: dExfiltrate_dTk       ! derivative in exfiltration w.r.t. temperature (K-1)
  ! ---------------------------------------------------------------------------------------
  ! * association to data in structures
  ! ---------------------------------------------------------------------------------------
  associate(&
    ! input: coordinate variables
    soilDepth               => prog_data%var(iLookPROG%iLayerHeight)%dat(nLayers),       & ! intent(in):  [dp]    total soil depth (m)
    mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat(nSnow+1:nLayers),& ! intent(in):  [dp(:)] depth of each soil layer (m)
    ! input: diagnostic variables
    surfaceHydCond          => flux_data%var(iLookFLUX%mLayerSatHydCondMP)%dat(1),       & ! intent(in):  [dp]    saturated hydraulic conductivity at the surface (m s-1)
    mLayerColumnInflow      => flux_data%var(iLookFLUX%mLayerColumnInflow)%dat,          & ! intent(in):  [dp(:)] inflow into each soil layer (m3/s)
    ! input: local attributes
    HRUarea                 => attr_data%var(iLookATTR%HRUarea),                         & ! intent(in):  [dp]    HRU area (m2)
    tan_slope               => attr_data%var(iLookATTR%tan_slope),                       & ! intent(in):  [dp]    tan water table slope, taken as tan local ground surface slope (-)
    contourLength           => attr_data%var(iLookATTR%contourLength),                   & ! intent(in):  [dp]    length of contour at downslope edge of HRU (m)
    ! input: baseflow parameters
    zScale_TOPMODEL         => mpar_data%var(iLookPARAM%zScale_TOPMODEL)%dat(1),         & ! intent(in):  [dp]    TOPMODEL exponent (-)
    kAnisotropic            => mpar_data%var(iLookPARAM%kAnisotropic)%dat(1),            & ! intent(in):  [dp]    anisotropy factor for lateral hydraulic conductivity (-
    fieldCapacity           => mpar_data%var(iLookPARAM%fieldCapacity)%dat(1),           & ! intent(in):  [dp]    field capacity (-)
    theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat,                  & ! intent(in):  [dp(:)] soil porosity (-)
    ! output: diagnostic variables
    scalarExfiltration      => flux_data%var(iLookFLUX%scalarExfiltration)%dat(1),       & ! intent(out): [dp]    exfiltration from the soil profile (m s-1)
    mLayerColumnOutflow     => flux_data%var(iLookFLUX%mLayerColumnOutflow)%dat          & ! intent(out): [dp(:)] column outflow from each soil layer (m3 s-1)
    )  ! end association to variables in data structures
    ! -----------------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message="computBaseflow/"
    ! ***********************************************************************************************************************
    ! (1) compute the baseflow flux in each soil layer
    ! ***********************************************************************************************************************

    ! compute the maximum transmissivity
    ! NOTE: this can be done as a pre-processing step
    tran0 = kAnisotropic*surfaceHydCond*soilDepth/zScale_TOPMODEL   ! maximum transmissivity (m2 s-1)

    ! compute the water table thickness (m) and transmissivity in each layer (m2 s-1)
    do iLayer=nSoil,ixSaturation,-1  ! loop through "active" soil layers, from lowest to highest
      ! define drainable water in each layer (m)
      activePorosity = theta_sat(iLayer) - fieldCapacity ! "active" porosity (-)
      drainableWater = mLayerDepth(iLayer)*(max(0._rkind,mLayerVolFracLiq(iLayer) - fieldCapacity))/activePorosity
      ! compute layer transmissivity
      if (iLayer==nSoil) then
        zActive(iLayer) = drainableWater                                       ! water table thickness associated with storage in a given layer (m)
        trTotal(iLayer) = tran0*(zActive(iLayer)/soilDepth)**zScale_TOPMODEL   ! total transmissivity for total depth zActive (m2 s-1)
        trSoil(iLayer)  = trTotal(iLayer)                                      ! transmissivity of water in a given layer (m2 s-1)
      else
        zActive(iLayer) = zActive(iLayer+1) + drainableWater
        trTotal(iLayer) = tran0*(zActive(iLayer)/soilDepth)**zScale_TOPMODEL
        trSoil(iLayer)  = trTotal(iLayer) - trTotal(iLayer+1)
      end if
    end do  ! end looping through soil layers

    ! set un-used portions of the vectors to zero
    if (ixSaturation>1) then
      zActive(1:ixSaturation-1) = 0._rkind
      trTotal(1:ixSaturation-1) = 0._rkind
      trSoil(1:ixSaturation-1)  = 0._rkind
    end if

    ! compute the outflow from each layer (m3 s-1)
    mLayerColumnOutflow(1:nSoil) = trSoil(1:nSoil)*tan_slope*contourLength

    ! compute total column inflow and total column outflow (m s-1)
    totalColumnInflow  = sum(mLayerColumnInflow(1:nSoil))/HRUarea
    totalColumnOutflow = sum(mLayerColumnOutflow(1:nSoil))/HRUarea

    ! compute the available storage (m)
    availStorage = sum(mLayerDepth(1:nSoil)*(theta_sat(1:nSoil) - (mLayerVolFracLiq(1:nSoil)+mLayerVolFracIce(1:nSoil))))

    ! compute the smoothing function (-)
    if (availStorage < xMinEval) then
      ! compute the logistic function
      expF = exp((availStorage - xCenter)/xWidth)
      logF = 1._rkind / (1._rkind + expF)
      ! compute the derivative in the logistic function w.r.t. volumetric water content in each soil layer, NOTE dLogFunc_dTemp = 0
      select case (ixRichards)
       case(moisture); dLogFunc_dWat(1:nSoil) = mLayerDepth(1:nSoil)*(expF/xWidth)/(1._rkind + expF)**2_i4b
       case(mixdform); dLogFunc_dWat(1:nSoil) = mLayerDepth(1:nSoil)*(expF/xWidth)/(1._rkind + expF)**2_i4b * dVolTot_dPsi0(1:nSoil)
       case default; err=20; message=trim(message)//'expect ixRichards to be moisture or mixdform'; return
      end select
    else
      logF             = 0._rkind
      dLogFunc_dWat(:) = 0._rkind
    end if

    ! compute the exfiltration (m s-1)
    if (totalColumnInflow > totalColumnOutflow .and. logF > tiny(1._rkind)) then
      scalarExfiltration = logF*(totalColumnInflow - totalColumnOutflow)  ! m s-1
    else
      scalarExfiltration = 0._rkind
    end if

    ! compute the baseflow in each layer (m s-1)
    mLayerBaseflow(1:nSoil) = (mLayerColumnOutflow(1:nSoil) - mLayerColumnInflow(1:nSoil))/HRUarea

    ! compute the total baseflow
    qbTotal = sum(mLayerBaseflow)

    ! add exfiltration to the baseflow flux at the top layer
    mLayerBaseflow(1)      = mLayerBaseflow(1) + scalarExfiltration
    mLayerColumnOutflow(1) = mLayerColumnOutflow(1) + scalarExfiltration*HRUarea

    ! ***********************************************************************************************************************
    ! (2) compute the derivative in the baseflow flux w.r.t. volumetric liquid water content (m s-1)
    ! ***********************************************************************************************************************

    ! initialize the derivative matrix
    dBaseflow_dVolLiq(:,:) = 0._rkind
    dBaseflow_dWat(:,:) = 0._rkind
    dBaseflow_dTk(:,:) = 0._rkind

    ! compute ratio of hillslope width to hillslope area (m m-2)
    length2area = tan_slope*contourLength/HRUarea

    ! compute the ratio of layer depth to maximum water holding capacity (-)
    depth2capacity(1:nSoil) = mLayerDepth(1:nSoil)/(theta_sat(1:nSoil) - fieldCapacity)/soilDepth
    do iLayer=1,nSoil
      if (mLayerVolFracLiq(iLayer) <= fieldCapacity) depth2capacity(iLayer) = 0._rkind
    end do

    ! compute the change in dimensionless flux w.r.t. change in dimensionless storage (-)
    dXdS(1:nSoil) = zScale_TOPMODEL*(zActive(1:nSoil)/soilDepth)**(zScale_TOPMODEL - 1._rkind)

    ! loop through soil layers
    do iLayer=1,nSoil
      ! compute diagonal terms (s-1)
      dBaseflow_dVolLiq(iLayer,iLayer) = tran0*dXdS(iLayer)*depth2capacity(iLayer)*length2area
      select case (ixRichards)
       case(moisture); dBaseflow_dWat(iLayer,iLayer) = dBaseflow_dVolLiq(iLayer,iLayer)*mLayerdTheta_dPsi(iLayer)
       case(mixdform); dBaseflow_dWat(iLayer,iLayer) = dBaseflow_dVolLiq(iLayer,iLayer)
      end select
      dBaseflow_dTk(iLayer,iLayer) = dBaseflow_dVolLiq(iLayer,iLayer)*mLayerdTheta_dTk(iLayer)
      ! compute off-diagonal terms
      do jLayer=iLayer+1,nSoil  ! only dependent on layers below
        dBaseflow_dVolLiq(iLayer,jLayer) = tran0*(dXdS(iLayer) - dXdS(iLayer+1))*depth2capacity(jLayer)*length2area
        select case (ixRichards)
         case(moisture); dBaseflow_dWat(iLayer,jLayer) = dBaseflow_dVolLiq(iLayer,jLayer)*mLayerdTheta_dPsi(jLayer)
         case(mixdform); dBaseflow_dWat(iLayer,jLayer) = dBaseflow_dVolLiq(iLayer,jLayer)
         end select
        dBaseflow_dTk(iLayer,jLayer) = dBaseflow_dVolLiq(iLayer,jLayer)*mLayerdTheta_dTk(jLayer)
      end do  ! end looping through soil layers
    end do  ! end looping through soil layers

    ! compute the derivative in the exfiltration flux and add to the baseflow derivative matrix
    if (totalColumnInflow > totalColumnOutflow .and. logF > tiny(1._rkind)) then
      do iLayer=1,nSoil
        dExfiltrate_dWat(iLayer) = -sum(dBaseflow_dWat(1:nSoil,iLayer))*logF - dLogFunc_dWat(iLayer)*qbTotal
        dExfiltrate_dTk(iLayer) = -sum(dBaseflow_dTk(1:nSoil,iLayer))*logF
      end do  ! end looping through soil layers
      dBaseflow_dWat(1,1:nSoil) = dBaseflow_dWat(1,1:nSoil) + dExfiltrate_dWat(1:nSoil)
      dBaseflow_dTk(1,1:nSoil) = dBaseflow_dTk(1,1:nSoil) + dExfiltrate_dTk(1:nSoil)
    end if

  end associate ! end association to data in structures

end subroutine computBaseflow

end module groundwatr_module
