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

module snLaIcLiqFlx_module

! access modules
USE nrtype                                 ! numerical recipes data types
USE multiconst,only:iden_ice,iden_water    ! intrinsic density of ice and water (kg m-3)

! access missing values
USE globalData,only:integerMissing         ! missing integer
USE globalData,only:realMissing            ! missing real number

! named variables
USE var_lookup,only:iLookINDEX             ! named variables for structure elements
USE var_lookup,only:iLookPARAM             ! named variables for structure elements
USE var_lookup,only:iLookPROG              ! named variables for structure elements
USE var_lookup,only:iLookDIAG              ! named variables for structure elements

! data types
USE data_types,only:var_dlength            ! x%var(:)%dat [rkind]
USE data_types,only:var_ilength            ! x%var(:)%dat [i4b]
USE data_types,only:in_type_snLaIcLiqFlx     ! data type for intent(in) arguments
USE data_types,only:io_type_snLaIcLiqFlx     ! data type for intent(inout) arguments
USE data_types,only:out_type_snLaIcLiqFlx    ! data type for intent(out) arguments

! privacy
implicit none
private
public :: snLaIcLiqFlx
contains
! ************************************************************************************************
! public subroutine snLaIcLiqFlx: compute liquid water flux through the snowpack
! ************************************************************************************************
subroutine snLaIcLiqFlx(&
                      ! input: model control, forcing, and model state vector
                      in_snLaIcLiqFlx,           & ! intent(in):    model control, forcing, and model state vector
                      ! input-output: data structures
                      mpar_data,               & ! intent(in):    model parameters
                      indx_data,               & ! intent(in):    model indices
                      prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                      ! input-output: fluxes and derivatives
                      io_snLaIcLiqFlx,           & ! intent(inout): fluxes and derivatives
                      ! output: error control
                      out_snLaIcLiqFlx)            ! intent(out):   error control
  implicit none
  ! input: model control, forcing, and model state vector
  type(in_type_snLaIcLiqFlx)          :: in_snLaIcLiqFlx              ! model control, forcing, and model state vector
  ! input-output: data structures
  type(var_dlength),intent(in)      :: mpar_data                  ! model parameters
  type(var_ilength),intent(in)      :: indx_data                  ! model indices
  type(var_dlength),intent(in)      :: prog_data                  ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)   :: diag_data                  ! diagnostic variables for a local HRU
  ! input-output: fluxes and derivatives
  type(io_type_snLaIcLiqFlx)          :: io_snLaIcLiqFlx              ! fluxes and derivatives
  ! output: error control
  type(out_type_snLaIcLiqFlx)         :: out_snLaIcLiqFlx             ! error control
  ! ------------------------------  ------------------------------------------------------------------------------------------------------------
  ! local variables
  integer(i4b)                      :: nLayers,nStart             ! number of snow/ice layers and starting layer
  integer(i4b)                      :: iLayer                     ! layer index
  integer(i4b)                      :: ixLayerDesired(1)          ! layer desired (scalar solution)
  integer(i4b)                      :: ixTop                      ! top layer in subroutine call
  integer(i4b)                      :: ixBot                      ! bottom layer in subroutine call
  real(rkind)                       :: multResid                  ! multiplier for the residual water content (-)
  real(rkind)                       :: residThrs                  ! ice density threshold to reduce residual liquid water content (kg m-3)
  real(rkind),parameter             :: residScal=10._rkind        ! scaling factor for residual liquid water content reduction factor (kg m-3)
  real(rkind)                       :: maxVolIceContent           ! maximum volumetric ice content to store water (-)
  real(rkind)                       :: availCap                   ! available storage capacity [0,1] (-)
  real(rkind)                       :: relSaturn                  ! relative saturation [0,1] (-)
  real(rkind)                       :: k_param                    ! hydraulic conductivity parameter (m s-1)
  logical(lgt)                      :: do_snow                    ! flag to denote if snow is present
  real(rkind)                       :: iLayerLiqFluxSnLaIc(0:in_snLaIcLiqFlx % nLayers)
  real(rkind)                       :: iLayerLiqFluxSnLaIcDeriv(0:in_snLaIcLiqFlx % nLayers)  
  ! ------------------------------------------------------------------------------------------------------------------------------------------
  ! make association of local variables with information in the data structures
  do_snow = in_snLaIcLiqFlx % do_snow ! flag to denote if snow is present
  nLayers = in_snLaIcLiqFlx % nLayers ! get number of snow/ice layers
  nStart = in_snLaIcLiqFlx % nStart ! get the start index for the layers

  associate(&
    ! input: model control
    firstFluxCall           => in_snLaIcLiqFlx % firstFluxCall,           & ! intent(in): the first flux call
    scalarSolution          => in_snLaIcLiqFlx % scalarSolution,          & ! intent(in): flag to denote if implementing the scalar solution
    ! input: forcing for the top layer
    surface_flux            => in_snLaIcLiqFlx % surface_flux,            & ! intent(in): liquid water flux at the surface (m s-1)
    ! input: model state vector
    mLayerVolFracLiqTrial   => in_snLaIcLiqFlx % mLayerVolFracLiqTrial,   & ! intent(in): trial value of volumetric fraction of liquid water at the current iteration (-)
    ! input: layer indices
    ixLayerState     => indx_data%var(iLookINDEX%ixLayerState)%dat,             & ! intent(in):    list of indices for all model layers
    ixSnowOnlyHyd    => indx_data%var(iLookINDEX%ixSnowOnlyHyd)%dat,            & ! intent(in):    index in the state subset for hydrology state variables in the snow domain
    ixIceOnlyHyd     => indx_data%var(iLookINDEX%ixIceOnlyHyd)%dat,             & ! intent(in):    index in the state subset for hydrology state variables in the ice domain
    ! input: snow properties and parameters
    mLayerVolFracIce => prog_data%var(iLookPROG%mLayerVolFracIce)%dat(nStart+1:nStart+nLayers), & ! intent(in):    volumetric ice content at the start of the time step (-)
    Fcapil           => mpar_data%var(iLookPARAM%Fcapil)%dat(1),                & ! intent(in):    capillary retention as a fraction of the total pore volume (-)
    k_snow           => mpar_data%var(iLookPARAM%k_snow)%dat(1),                & ! intent(in):    hydraulic conductivity of snow (m s-1)    
    k_ice            => mpar_data%var(iLookPARAM%k_ice)%dat(1),                 & ! intent(in):    hydraulic conductivity of ice (m s-1)
    mw_exp           => mpar_data%var(iLookPARAM%mw_exp)%dat(1),                & ! intent(in):    exponent for meltwater flow (-)
    ! input-output: diagnostic variables -- only computed for the first iteration
    mLayerPoreSpace  => diag_data%var(iLookDIAG%mLayerPoreSpace)%dat(nStart+1:nStart+nLayers),  & ! intent(inout): pore space in each layer (-)
    mLayerThetaResid => diag_data%var(iLookDIAG%mLayerThetaResid)%dat(nStart+1:nStart+nLayers), & ! intent(inout): residual volumetric liquid water content in each layer (-)
    ! input-output: fluxes and derivatives
    iLayerLiqFluxSnLaIc0      => io_snLaIcLiqFlx % iLayerLiqFluxSnLaIc,               & ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
    iLayerLiqFluxSnLaIcDeriv0 => io_snLaIcLiqFlx % iLayerLiqFluxSnLaIcDeriv,          & ! intent(inout): derivative in vertical liquid water flux at layer interfaces (m s-1)
    ! output: error control
    err                    => out_snLaIcLiqFlx % err,                             & ! intent(out):   error code
    message                => out_snLaIcLiqFlx % cmessage                         & ! intent(out):   error message
    ) ! end association of local variables with information in the data structures
    ! ------------------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='snLaIcLiqFlx/'

    ! initialize with index 0
    iLayerLiqFluxSnLaIc = iLayerLiqFluxSnLaIc0
    iLayerLiqFluxSnLaIcDeriv = iLayerLiqFluxSnLaIcDeriv0

    ! check that the input vectors match nLayers
    if (size(mLayerVolFracLiqTrial)/=nLayers .or. size(mLayerVolFracIce)/=nLayers .or. &
        size(iLayerLiqFluxSnLaIc)/=nLayers+1 .or. size(iLayerLiqFluxSnLaIcDeriv)/=nLayers+1) then
      err=20; message=trim(message)//'size mismatch of input/output vectors'; return
    end if

    ! check the meltwater exponent is >=1
    if (mw_exp<1._rkind) then; err=20; message=trim(message)//'meltwater exponent < 1'; return; end if

    ! get the inputs for the layers
    if (do_snow)then
      residThrs = 550._rkind
      maxVolIceContent = 0.7_rkind
      k_param = k_snow
    else ! ice
      residThrs = 800._rkind
      maxVolIceContent = 0.9_rkind
      k_param = k_ice
    end if

    ! get the indices for the layers
    ixTop = integerMissing
    if (scalarSolution) then
      if (do_snow)then
        ixLayerDesired = pack(ixLayerState, ixSnowOnlyHyd/=integerMissing)
      else
        ixLayerDesired = pack(ixLayerState, ixIceOnlyHyd/=integerMissing)
      end if
      ixTop = ixLayerDesired(1)
      ixBot = ixLayerDesired(1)
    else
      ixTop = 1
      ixBot = nLayers
    end if

    ! define the liquid flux at the upper boundary (m s-1)
    iLayerLiqFluxSnLaIc(0)      = surface_flux
    iLayerLiqFluxSnLaIcDeriv(0) = 0._rkind !computed inside computJacob

    ! compute properties fixed over the time step
    if (firstFluxCall) then
      ! loop through snow/ice layers
      do iLayer=1,nLayers ! loop through snow/ice layers
        multResid = 1._rkind/(1._rkind + exp((mLayerVolFracIce(iLayer)*iden_ice - residThrs)/residScal)) ! compute the reduction in liquid water holding capacity at high snow/ice density (-)
        mLayerPoreSpace(iLayer)  = 1._rkind - mLayerVolFracIce(iLayer) ! compute the pore space (-)
        mLayerThetaResid(iLayer) = Fcapil*mLayerPoreSpace(iLayer)*multResid ! compute the residual volumetric liquid water content (-)
      end do  ! end looping through snow/ice layers
    end if  ! end if the first flux call
     
    ! compute fluxes
    do iLayer=ixTop,ixBot  ! loop through snow/ice layers
      if (mLayerVolFracLiqTrial(iLayer) > mLayerThetaResid(iLayer)) then ! check that flow occurs
        ! compute the relative saturation (-)
        availCap  = mLayerPoreSpace(iLayer) - mLayerThetaResid(iLayer)                 ! available capacity
        relSaturn = (mLayerVolFracLiqTrial(iLayer) - mLayerThetaResid(iLayer)) / availCap    ! relative saturation
        iLayerLiqFluxSnLaIc(iLayer)      = k_param*relSaturn**mw_exp
        iLayerLiqFluxSnLaIcDeriv(iLayer) = ( (k_param*mw_exp)/availCap ) * relSaturn**(mw_exp - 1._rkind)
        if (mLayerVolFracIce(iLayer) > maxVolIceContent) then ! NOTE: use start-of-step ice content, to avoid convergence problems
          ! ** allow liquid water to pass through under very high ice density
          iLayerLiqFluxSnLaIc(iLayer) = iLayerLiqFluxSnLaIc(iLayer) + iLayerLiqFluxSnLaIc(iLayer-1) !NOTE: derivative may need to be updated in future.
        end if
      else  ! flow does not occur
        iLayerLiqFluxSnLaIc(iLayer)      = 0._rkind
        iLayerLiqFluxSnLaIcDeriv(iLayer) = 0._rkind
      end if  ! storage above residual content
    end do  ! end loop through snow/ice layers

    ! save the results with index 0
    iLayerLiqFluxSnLaIc0 = iLayerLiqFluxSnLaIc
    iLayerLiqFluxSnLaIcDeriv0 = iLayerLiqFluxSnLaIcDeriv

  end associate ! end association of local variables with information in the data structures

end subroutine snLaIcLiqFlx

end module snLaIcLiqFlx_module
