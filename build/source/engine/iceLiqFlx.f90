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

module iceLiqFlx_module

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
USE data_types,only:in_type_iceLiqFlx     ! data type for intent(in) arguments
USE data_types,only:io_type_iceLiqFlx     ! data type for intent(inout) arguments
USE data_types,only:out_type_iceLiqFlx    ! data type for intent(out) arguments

! privacy
implicit none
private
public :: iceLiqFlx
contains
! ************************************************************************************************
! public subroutine iceLiqFlx: compute liquid water flux through the icepack
! ************************************************************************************************
subroutine iceLiqFlx(&
                      ! input: model control, forcing, and model state vector
                      in_iceLiqFlx,           & ! intent(in):    model control, forcing, and model state vector
                      ! input-output: data structures
                      mpar_data,               & ! intent(in):    model parameters
                      indx_data,               & ! intent(in):    model indices
                      prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                      ! input-output: fluxes and derivatives
                      io_iceLiqFlx,           & ! intent(inout): fluxes and derivatives
                      ! output: error control
                      out_iceLiqFlx)            ! intent(out):   error control
  implicit none
  ! input: model control, forcing, and model state vector
  type(in_type_iceLiqFlx)          :: in_iceLiqFlx              ! model control, forcing, and model state vector
  ! input-output: data structures
  type(var_dlength),intent(in)      :: mpar_data                  ! model parameters
  type(var_ilength),intent(in)      :: indx_data                  ! model indices
  type(var_dlength),intent(in)      :: prog_data                  ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)   :: diag_data                  ! diagnostic variables for a local HRU
  ! input-output: fluxes and derivatives
  type(io_type_iceLiqFlx)          :: io_iceLiqFlx              ! fluxes and derivatives
  ! output: error control
  type(out_type_iceLiqFlx)         :: out_iceLiqFlx             ! error control
  ! ------------------------------  ------------------------------------------------------------------------------------------------------------
  ! local variables
  integer(i4b)                      :: nIce                      ! number of ice layers
  integer(i4b)                      :: i                          ! search index for scalar solution
  integer(i4b)                      :: iLayer                     ! layer index
  integer(i4b)                      :: ixTop                      ! top layer in subroutine call
  integer(i4b)                      :: ixBot                      ! bottom layer in subroutine call
  real(rkind)                       :: multResid                  ! multiplier for the residual water content (-)
  real(rkind),parameter             :: residThrs=550._rkind       ! ice density threshold to reduce residual liquid water content (kg m-3)
  real(rkind),parameter             :: residScal=10._rkind        ! scaling factor for residual liquid water content reduction factor (kg m-3)
  real(rkind),parameter             :: maxVolIceContent=0.99_rkind! maximum volumetric ice content to store water (-)
  real(rkind)                       :: availCap                   ! available storage capacity [0,1] (-)
  real(rkind)                       :: relSaturn                  ! relative saturation [0,1] (-)
  ! ------------------------------------------------------------------------------------------------------------------------------------------
  ! make association of local variables with information in the data structures
  nIce=in_iceLiqFlx % nIce ! get number of ice layers
  associate(&
    ! input: model control
    firstFluxCall           => in_iceLiqFlx % firstFluxCall,           & ! intent(in): the first flux call
    scalarSolution          => in_iceLiqFlx % scalarSolution,          & ! intent(in): flag to denote if implementing the scalar solution
    ! input: forcing for the ice domain
    scalarSoilDrainage => in_iceLiqFlx % scalarSoilDrainage, & ! intent(in): 
    scalarRainPlusMelt => in_iceLiqFlx % scalarRainPlusMelt, & ! intent(in): 
    ! input: model state vector
    mLayerVolFracLiqTrial   => in_iceLiqFlx % mLayerVolFracLiqTrial,   & ! intent(in): trial value of volumetric fraction of liquid water at the current iteration (-)
    ! input: layer indices
    ixLayerState     => indx_data%var(iLookINDEX%ixLayerState)%dat,            & ! intent(in):    list of indices for all model layers
    ixIceOnlyHyd     => indx_data%var(iLookINDEX%ixIceOnlyHyd)%dat,            & ! intent(in):    index in the state subset for hydrology state variables in the ice domain
    ! input: ice properties and parameters
    mLayerVolFracIce => prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1:nIce), & ! intent(in):    volumetric ice content at the start of the time step (-)
    Fcapil           => mpar_data%var(iLookPARAM%Fcapil)%dat(1),               & ! intent(in):    capillary retention as a fraction of the total pore volume (-)
    k_ice            => mpar_data%var(iLookPARAM%k_ice)%dat(1),                & ! intent(in):    hydraulic conductivity of ice (m s-1), 2.1e-6 = approx. 0.05 mm/hr, from Stevens et al 2017
    mw_exp           => mpar_data%var(iLookPARAM%mw_exp)%dat(1),               & ! intent(in):    exponent for meltwater flow (-)
    ! input-output: diagnostic variables -- only computed for the first iteration
    mLayerPoreSpace  => diag_data%var(iLookDIAG%mLayerPoreSpace)%dat,           & ! intent(inout): pore space in each ice layer (-)
    mLayerThetaResid => diag_data%var(iLookDIAG%mLayerThetaResid)%dat,          & ! intent(inout): residual volumetric liquid water content in each ice layer (-)
    ! input-output: fluxes and derivatives
    iLayerLiqFluxIce       => io_iceLiqFlx % iLayerLiqFluxIce,                 & ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
    iLayerLiqFluxIceDeriv  => io_iceLiqFlx % iLayerLiqFluxIceDeriv,            & ! intent(inout): derivative in vertical liquid water flux at layer interfaces (m s-1)
    ! output: error control
    err                    => out_iceLiqFlx % err,                             & ! intent(out):   error code
    message                => out_iceLiqFlx % cmessage                         & ! intent(out):   error message
    ) ! end association of local variables with information in the data structures
    ! ------------------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='iceLiqFlx/'

    ! check that the input vectors match nIce
    if (size(mLayerVolFracLiqTrial)/=nIce .or. size(mLayerVolFracIce)/=nIce .or. &
        size(iLayerLiqFluxIce)/=nIce+1 .or. size(iLayerLiqFluxIceDeriv)/=nIce+1) then
      err=20; message=trim(message)//'size mismatch of input/output vectors'; return
    end if

    ! check the meltwater exponent is >=1
    if (mw_exp<1._rkind) then; err=20; message=trim(message)//'meltwater exponent < 1'; return; end if

    ! get the indices for the ice+soil layers
    ixTop = integerMissing
    if (scalarSolution) then
      do i=1,size(ixIceOnlyHyd)
        if (ixIceOnlyHyd(i) /= integerMissing) then
          ixTop=ixLayerState(i)
          ixBot=ixTop
          exit  ! break out of loop once found
        end if
      end do
      if (ixTop == integerMissing) then
        err=20; message=trim(message)//'Unable to identify ice layer for scalar solution!'; return
      end if
    else
      ixTop = 1
      ixBot = nIce
    end if

    ! define the liquid flux at the upper boundary (m s-1)
    iLayerLiqFluxIce(0)      = scalarSurfaceRunoff + scalarSoilDrainage
    iLayerLiqFluxIceDeriv(0) = ?? 

    ! compute properties fixed over the time step
    if (firstFluxCall) then
      ! loop through ice layers
      do iLayer=1,nIce ! loop through ice layers
        multResid = 1._rkind/(1._rkind + exp((mLayerVolFracIce(iLayer)*iden_ice - residThrs)/residScal)) ! compute the reduction in liquid water holding capacity at high ice density (-)
        mLayerPoreSpace(iLayer)  = 1._rkind - mLayerVolFracIce(iLayer) ! compute the pore space (-)
        mLayerThetaResid(iLayer) = Fcapil*mLayerPoreSpace(iLayer)*multResid ! compute the residual volumetric liquid water content (-)
      end do  ! end looping through ice layers
    end if  ! end if the first flux call
     
    ! compute fluxes
    do iLayer=ixTop,ixBot  ! loop through ice layers
      if (mLayerVolFracLiqTrial(iLayer) > mLayerThetaResid(iLayer)) then ! check that flow occurs
        ! compute the relative saturation (-)
        availCap  = mLayerPoreSpace(iLayer) - mLayerThetaResid(iLayer)                 ! available capacity
        relSaturn = (mLayerVolFracLiqTrial(iLayer) - mLayerThetaResid(iLayer)) / availCap    ! relative saturation
        iLayerLiqFluxIce(iLayer)      = k_ice*relSaturn**mw_exp
        iLayerLiqFluxIceDeriv(iLayer) = ( (k_ice*mw_exp)/availCap ) * relSaturn**(mw_exp - 1._rkind)
        if (mLayerVolFracIce(iLayer) > maxVolIceContent) then ! NOTE: use start-of-step ice content, to avoid convergence problems
          ! ** allow liquid water to pass through under very high ice density
          iLayerLiqFluxIce(iLayer) = iLayerLiqFluxIce(iLayer) + iLayerLiqFluxIce(iLayer-1) !NOTE: derivative may need to be updated in future.
        end if
      else  ! flow does not occur
        iLayerLiqFluxIce(iLayer)      = 0._rkind
        iLayerLiqFluxIceDeriv(iLayer) = 0._rkind
      end if  ! storage above residual content
    end do  ! end loop through ice layers

  end associate ! end association of local variables with information in the data structures

end subroutine iceLiqFlx

end module iceLiqFlx_module
