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

module lakeIceCover_module

! The ice cover of a lake or stream, after van Beek et al. (2012) and Wanders et al. (2019):
! an ice layer at the freezing point on top of a mixed water body. The lake layers are of two
! kinds, both of type iname_lake: the top nLakeFrz are ice (ice density, residual liquid, treated
! like glacier ice: impermeable, melt squeezed upward, sharp freezing curve) and the rest are
! water (full of liquid, mixed, advected by the reach flow). This module moves ice between the
! two once per outer step, next to the snow layer merging and dividing:
!   * freeze-up: ice that formed in the water layers (their energy dropped below the plateau at
!     the freezing point; the water is mixed, so it rises) is moved into the ice cover, creating
!     the ice layer when the cover is thick enough to stand on its own, otherwise thickening it;
!   * breakup: an ice cover thinner than lakeIceMinThick (5 mm in flowing water, Wanders et al.)
!     returns to the water beneath, mass and latent heat conserved.
! The mass of a lake layer sets its depth (no air), so the 9% expansion on freezing happens here,
! where the ice moves, and not inside the solver.

USE nr_type
USE data_types,only:var_ilength,var_dlength,var_info ! data vectors with variable length dimension, and metadata
USE multiconst,only:Tfreeze,iden_ice,iden_water      ! physical constants
USE globalData,only:realMissing,integerMissing       ! missing values
USE globalData,only:iname_lake                       ! named variable for lake layers
USE globalData,only:iceResidWaterFrac                ! residual volumetric liquid water content in ice (-)
USE globalData,only:nLakeIceLayers_poss              ! number of ice cover layers a lake can grow
USE globalData,only:verySmall                        ! a small number
USE globalData,only:icefrz_mult                      ! freezing curve scaling factor multiplier of snow to ice
USE globalData,only:prog_meta,diag_meta,flux_meta,indx_meta ! metadata, to resize the data structures
USE var_lookup,only:iLookPROG,iLookINDEX,iLookPARAM,iLookVarType ! named variables

implicit none
private
public::lakeIceCover

contains

! ************************************************************************************************
! public subroutine lakeIceCover: grow or break up the ice cover of the lake layers
! ************************************************************************************************
subroutine lakeIceCover(mpar_data,indx_data,prog_data,diag_data,flux_data,modifiedLayers,err,message)
  USE var_derive_module,only:calcHeight   ! compute the height at layer interfaces and layer mid-points
  USE snow_utils_module,only:templiquid   ! temperature at which the freezing curve gives a liquid fraction
  implicit none
  type(var_dlength),intent(in)    :: mpar_data       ! model parameters
  type(var_ilength),intent(inout) :: indx_data       ! model indices
  type(var_dlength),intent(inout) :: prog_data       ! model prognostic variables
  type(var_dlength),intent(inout) :: diag_data       ! model diagnostic variables
  type(var_dlength),intent(inout) :: flux_data       ! model fluxes
  logical(lgt),intent(out)        :: modifiedLayers  ! flag to denote that the layers were modified
  integer(i4b),intent(out)        :: err             ! error code
  character(*),intent(out)        :: message         ! error message
  ! local variables
  character(len=256)              :: cmessage        ! error message of downwind routine
  integer(i4b)                    :: nSnow,nLake,nSoil,nGlce,nLayers,nLakeFrz ! layer counts
  integer(i4b)                    :: ixWat           ! index of the top water layer
  integer(i4b)                    :: iLayer          ! layer index
  integer(i4b)                    :: ixIce           ! index of the bottom ice layer (the one in contact with the water)
  real(rkind)                     :: massIce         ! ice mass moving between the water and the cover (kg m-2)
  real(rkind)                     :: massLiq         ! liquid mass moving with it (kg m-2)
  real(rkind)                     :: depthIce        ! thickness the ice takes at ice density (m)
  real(rkind)                     :: massIceCover    ! ice mass of the cover (kg m-2)
  real(rkind)                     :: massLiqCover    ! liquid mass of the cover (kg m-2)
  real(rkind)                     :: massLiqWat      ! liquid mass of the top water layer (kg m-2)
  real(rkind)                     :: massIceWat      ! ice mass of the top water layer (kg m-2)
  real(rkind)                     :: depthNew        ! updated layer depth (m)
  real(rkind)                     :: minMassIce      ! ice mass of a cover at the minimum thickness (kg m-2)
  ! ----------------------------------------------------------------------------------------------------------------
  err=0; message='lakeIceCover/'
  modifiedLayers = .false.

  ! NOTE: the layer vectors are reallocated when a layer is added or removed, so they are referenced through the
  !       data structure rather than through an associate block
  nSnow    = indx_data%var(iLookINDEX%nSnow)%dat(1)
  nLake    = indx_data%var(iLookINDEX%nLake)%dat(1)
  nSoil    = indx_data%var(iLookINDEX%nSoil)%dat(1)
  nGlce    = indx_data%var(iLookINDEX%nGlce)%dat(1)
  nLayers  = indx_data%var(iLookINDEX%nLayers)%dat(1)
  nLakeFrz = indx_data%var(iLookINDEX%nLakeFrz)%dat(1)
  if(nLake==0) return
  if(nLakeFrz >= nLake)then; err=20; message=trim(message)//'the lake has no water layer beneath its ice'; return; end if

  ! minimum ice cover mass is sized for the minimum thickness at ice density, with no liquid
  minMassIce = mpar_data%var(iLookPARAM%lakeIceMinThick)%dat(1)*iden_ice*(1._rkind - iceResidWaterFrac)

  ! ***** breakup: an ice cover too thin to stand returns to the water beneath
  if(nLakeFrz > 0)then
    ixIce = nSnow + nLakeFrz
    ixWat = ixIce + 1
    massIceCover = depth(ixIce)*ice(ixIce)*iden_ice
    if(massIceCover < minMassIce)then
      massLiqCover = depth(ixIce)*liq(ixIce)*iden_water
      massLiqWat   = depth(ixWat)*liq(ixWat)*iden_water
      massIceWat   = depth(ixWat)*ice(ixWat)*iden_ice
      ! the water layer takes the ice (as ice, so the latent heat stays) and the liquid of the cover;
      ! its temperature is the liquid-weighted one (the cover's liquid is at the freezing point)
      depthNew = (massLiqWat + massLiqCover)/iden_water + (massIceWat + massIceCover)/iden_ice
      call setTemp(ixWat, (massLiqWat*temp(ixWat) + massLiqCover*Tfreeze)/max(massLiqWat + massLiqCover, verySmall))
      call setLiq(ixWat, (massLiqWat + massLiqCover)/(iden_water*depthNew))
      call setIce(ixWat, (massIceWat + massIceCover)/(iden_ice*depthNew))
      call setDepth(ixWat, depthNew)
      ! remove the ice layer from all model vectors
      call rmLakeLayer(prog_data,prog_meta,ixIce,nSnow,nLake,nLayers,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
      call rmLakeLayer(diag_data,diag_meta,ixIce,nSnow,nLake,nLayers,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
      call rmLakeLayer(flux_data,flux_meta,ixIce,nSnow,nLake,nLayers,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
      call rmLakeLayer(indx_data,indx_meta,ixIce,nSnow,nLake,nLayers,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
      nLake = nLake - 1; nLakeFrz = nLakeFrz - 1; nLayers = nLayers - 1
      modifiedLayers = .true.
    end if
  end if

  ! ***** freeze-up: ice that formed in the water (well mixed, so in any water layer) rises into the cover
  ixWat = nSnow + nLakeFrz + 1
  massIceWat = 0._rkind
  do iLayer=ixWat,nSnow+nLake
    massIceWat = massIceWat + depth(iLayer)*ice(iLayer)*iden_ice
  end do
  ! a new cover needs to be thick enough to stand (twice the breakup thickness, so it is not undone at once);
  ! an existing cover takes whatever ice there is
  if( (nLakeFrz==0 .and. massIceWat >= 2._rkind*minMassIce) .or. (nLakeFrz>0 .and. massIceWat > 0._rkind) )then
    massIce  = massIceWat
    ! the ice keeps a residual liquid fraction, taken from the top water layer; the layer is sized for that liquid at
    ! the density of ice, so the ice fraction stays at or below one when the residual freezes later (iceReduce keeps
    ! the fractions as it thins)
    massLiq  = min(iceResidWaterFrac*massIce*(iden_water/iden_ice)/(1._rkind - iceResidWaterFrac), 0.5_rkind*depth(ixWat)*liq(ixWat)*iden_water)
    depthIce = (massIce + massLiq)/iden_ice
    ! each water layer keeps its liquid at full liquid density; the top one also gives up the residual that goes with the ice
    do iLayer=ixWat,nSnow+nLake
      massLiqWat = depth(iLayer)*liq(iLayer)*iden_water
      if(iLayer==ixWat) massLiqWat = massLiqWat - massLiq
      if(massLiqWat < verySmall)then; err=20; message=trim(message)//'a lake water layer froze solid; its ice should have moved to the cover earlier'; return; end if
      call setDepth(iLayer, massLiqWat/iden_water)
      call setLiq(iLayer, 1._rkind)
      call setIce(iLayer, 0._rkind)
      call setTemp(iLayer, max(temp(iLayer), Tfreeze)) ! the water is liquid again; its cold content went into the ice
    end do
    if(nLakeFrz==0)then
      ! create the ice layer above the water layer, then set its state
      if(nLakeIceLayers_poss < 1)then; err=20; message=trim(message)//'no ice layer allowed on the lake (nLakeIceLayers_poss)'; return; end if
      call addLakeLayer(prog_data,prog_meta,ixWat,nSnow,nLake,nLayers,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
      call addLakeLayer(diag_data,diag_meta,ixWat,nSnow,nLake,nLayers,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
      call addLakeLayer(flux_data,flux_meta,ixWat,nSnow,nLake,nLayers,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
      call addLakeLayer(indx_data,indx_meta,ixWat,nSnow,nLake,nLayers,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
      nLake = nLake + 1; nLakeFrz = 1; nLayers = nLayers + 1
      ixIce = ixWat  ! the new layer took the water layer's place, the water layer is now beneath it
      call setDepth(ixIce, depthIce)
      call setIce(ixIce, massIce/(iden_ice*depthIce))
      call setLiq(ixIce, massLiq/(iden_water*depthIce))
      call setTemp(ixIce, iceTemp(liq(ixIce),ice(ixIce)))
    else
      ! thicken the bottom ice layer with the new ice
      ixIce = nSnow + nLakeFrz
      massIceCover = depth(ixIce)*ice(ixIce)*iden_ice
      massLiqCover = depth(ixIce)*liq(ixIce)*iden_water
      depthNew = (massIceCover + massIce + massLiqCover + massLiq)/iden_ice  ! liquid sized at the density of ice, as above
      call setIce(ixIce, (massIceCover + massIce)/(iden_ice*depthNew))
      call setLiq(ixIce, (massLiqCover + massLiq)/(iden_water*depthNew))
      call setDepth(ixIce, depthNew)
      call setTemp(ixIce, min(temp(ixIce), iceTemp(liq(ixIce),ice(ixIce))))
    end if
    modifiedLayers = .true.
  end if

  ! save the layer counts and types, and the heights that follow from the new depths
  if(modifiedLayers)then
    indx_data%var(iLookINDEX%layerType)%dat(nSnow+1:nSnow+nLake) = iname_lake
    indx_data%var(iLookINDEX%nLake)%dat(1)    = nLake
    indx_data%var(iLookINDEX%nLakeFrz)%dat(1) = nLakeFrz
    indx_data%var(iLookINDEX%nLayers)%dat(1)  = nLayers
    call calcHeight(indx_data,prog_data,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  end if

contains

  ! accessors to the layer state, valid across the reallocation of the layer vectors
  function depth(i); integer(i4b),intent(in) :: i; real(rkind) :: depth; depth = prog_data%var(iLookPROG%mLayerDepth)%dat(i);      end function depth
  function temp(i);  integer(i4b),intent(in) :: i; real(rkind) :: temp;  temp  = prog_data%var(iLookPROG%mLayerTemp)%dat(i);       end function temp
  function liq(i);   integer(i4b),intent(in) :: i; real(rkind) :: liq;   liq   = prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(i); end function liq
  function ice(i);   integer(i4b),intent(in) :: i; real(rkind) :: ice;   ice   = prog_data%var(iLookPROG%mLayerVolFracIce)%dat(i); end function ice
  subroutine setDepth(i,x); integer(i4b),intent(in) :: i; real(rkind),intent(in) :: x; prog_data%var(iLookPROG%mLayerDepth)%dat(i)      = x; end subroutine setDepth
  subroutine setTemp(i,x);  integer(i4b),intent(in) :: i; real(rkind),intent(in) :: x; prog_data%var(iLookPROG%mLayerTemp)%dat(i)       = x; end subroutine setTemp
  subroutine setLiq(i,x);   integer(i4b),intent(in) :: i; real(rkind),intent(in) :: x; prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(i) = x; end subroutine setLiq
  subroutine setIce(i,x);   integer(i4b),intent(in) :: i; real(rkind),intent(in) :: x; prog_data%var(iLookPROG%mLayerVolFracIce)%dat(i) = x; end subroutine setIce

  ! the temperature at which the ice freezing curve holds the layer's liquid fraction, so the solver sees ice, not water at 0 C
  function iceTemp(volFracLiq,volFracIce)
    real(rkind),intent(in) :: volFracLiq,volFracIce ! volumetric fractions of liquid water and ice (-)
    real(rkind)            :: iceTemp               ! temperature (K)
    real(rkind)            :: fLiq                  ! liquid fraction of the total water, by mass (-)
    fLiq = volFracLiq/(volFracLiq + volFracIce*(iden_ice/iden_water))
    iceTemp = templiquid(max(fLiq,1.e-6_rkind), mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)*icefrz_mult)
  end function iceTemp

end subroutine lakeIceCover

! ************************************************************************************************
! private subroutine addLakeLayer: insert a lake layer at index ixInsert, copying that layer's state
! ************************************************************************************************
! Follows addModelLayer of layerDivide, for the lake dimension: the layer at ixInsert is duplicated,
! the new copy taking index ixInsert and the original moving to ixInsert+1.
subroutine addLakeLayer(dataStruct,metaStruct,ixInsert,nSnow,nLake,nLayers,err,message)
  USE f2008_funcs_module,only:cloneStruc               ! used to "clone" data structures
  implicit none
  class(*),intent(inout)          :: dataStruct        ! data structure
  type(var_info),intent(in)       :: metaStruct(:)     ! metadata structure
  integer(i4b),intent(in)         :: ixInsert          ! index of the layer to duplicate (in the whole column)
  integer(i4b),intent(in)         :: nSnow,nLake,nLayers ! number of snow layers, lake layers, total number of layers
  integer(i4b),intent(out)        :: err               ! error code
  character(*),intent(out)        :: message           ! error message
  integer(i4b)                    :: iVar              ! index of model variable
  integer(i4b)                    :: ix_lower,ix_upper ! bounds of the vector
  integer(i4b)                    :: ix_divide         ! index of the layer to duplicate within the vector
  logical(lgt)                    :: stateVariable     ! .true. if variable is a state variable
  real(rkind),allocatable         :: tempVec_rkind(:)  ! temporary vector (double precision)
  integer(i4b),allocatable        :: tempVec_i4b(:)    ! temporary vector (integer)
  character(LEN=256)              :: cmessage          ! error message of downwind routine
  err=0; message='addLakeLayer/'

  do iVar=1,size(metaStruct)
    select case(metaStruct(iVar)%varType)
      case(iLookVarType%midLake); ix_lower=1; ix_upper=nLake;   ix_divide=ixInsert-nSnow
      case(iLookVarType%midToto); ix_lower=1; ix_upper=nLayers; ix_divide=ixInsert
      case(iLookVarType%ifcLake); ix_lower=0; ix_upper=nLake;   ix_divide=ixInsert-nSnow
      case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nLayers; ix_divide=ixInsert
      case default; cycle
    end select
    select case(trim(metaStruct(iVar)%varName))
      case('mLayerDepth','mLayerTemp','mLayerVolFracIce','mLayerVolFracLiq','layerType'); stateVariable=.true.
      case default; stateVariable=.false.
    end select

    select type(dataStruct)
      type is (var_dlength)
        if(.not.allocated(dataStruct%var(iVar)%dat))then; err=20; message=trim(message)//'data vector is not allocated'; return; end if
        call cloneStruc(tempVec_rkind, ix_lower, source=dataStruct%var(iVar)%dat, err=err, message=cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
        deallocate(dataStruct%var(iVar)%dat); allocate(dataStruct%var(iVar)%dat(ix_lower:ix_upper+1),stat=err)
        if(err/=0)then; err=20; message=trim(message)//'problem reallocating data vector'; return; end if
        if(stateVariable)then
          dataStruct%var(iVar)%dat(ix_lower:ix_divide) = tempVec_rkind(ix_lower:ix_divide)
          dataStruct%var(iVar)%dat(ix_divide+1)        = tempVec_rkind(ix_divide)
          if(ix_upper > ix_divide) dataStruct%var(iVar)%dat(ix_divide+2:ix_upper+1) = tempVec_rkind(ix_divide+1:ix_upper)
        else
          dataStruct%var(iVar)%dat(:) = realMissing
        end if
        deallocate(tempVec_rkind)
      type is (var_ilength)
        if(.not.allocated(dataStruct%var(iVar)%dat))then; err=20; message=trim(message)//'data vector is not allocated'; return; end if
        call cloneStruc(tempVec_i4b, ix_lower, source=dataStruct%var(iVar)%dat, err=err, message=cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
        deallocate(dataStruct%var(iVar)%dat); allocate(dataStruct%var(iVar)%dat(ix_lower:ix_upper+1),stat=err)
        if(err/=0)then; err=20; message=trim(message)//'problem reallocating data vector'; return; end if
        if(stateVariable)then
          dataStruct%var(iVar)%dat(ix_lower:ix_divide) = tempVec_i4b(ix_lower:ix_divide)
          dataStruct%var(iVar)%dat(ix_divide+1)        = tempVec_i4b(ix_divide)
          if(ix_upper > ix_divide) dataStruct%var(iVar)%dat(ix_divide+2:ix_upper+1) = tempVec_i4b(ix_divide+1:ix_upper)
        else
          dataStruct%var(iVar)%dat(:) = integerMissing
        end if
        deallocate(tempVec_i4b)
      class default; err=20; message=trim(message)//'unable to identify the data structure type'; return
    end select
  end do
end subroutine addLakeLayer

! ************************************************************************************************
! private subroutine rmLakeLayer: remove the lake layer at index ixRemove from all model vectors
! ************************************************************************************************
subroutine rmLakeLayer(dataStruct,metaStruct,ixRemove,nSnow,nLake,nLayers,err,message)
  USE f2008_funcs_module,only:cloneStruc               ! used to "clone" data structures
  implicit none
  class(*),intent(inout)          :: dataStruct        ! data structure
  type(var_info),intent(in)       :: metaStruct(:)     ! metadata structure
  integer(i4b),intent(in)         :: ixRemove          ! index of the layer to remove (in the whole column)
  integer(i4b),intent(in)         :: nSnow,nLake,nLayers ! number of snow layers, lake layers, total number of layers
  integer(i4b),intent(out)        :: err               ! error code
  character(*),intent(out)        :: message           ! error message
  integer(i4b)                    :: iVar              ! index of model variable
  integer(i4b)                    :: ix_lower,ix_upper ! bounds of the vector
  integer(i4b)                    :: ix_rm             ! index of the layer to remove within the vector
  real(rkind),allocatable         :: tempVec_rkind(:)  ! temporary vector (double precision)
  integer(i4b),allocatable        :: tempVec_i4b(:)    ! temporary vector (integer)
  character(LEN=256)              :: cmessage          ! error message of downwind routine
  err=0; message='rmLakeLayer/'

  do iVar=1,size(metaStruct)
    select case(metaStruct(iVar)%varType)
      case(iLookVarType%midLake); ix_lower=1; ix_upper=nLake;   ix_rm=ixRemove-nSnow
      case(iLookVarType%midToto); ix_lower=1; ix_upper=nLayers; ix_rm=ixRemove
      case(iLookVarType%ifcLake); ix_lower=0; ix_upper=nLake;   ix_rm=ixRemove-nSnow
      case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nLayers; ix_rm=ixRemove
      case default; cycle
    end select
    ! NOTE: interface vectors drop the interface below the removed layer (the layers above and below share the remaining one)
    select type(dataStruct)
      type is (var_dlength)
        if(.not.allocated(dataStruct%var(iVar)%dat))then; err=20; message=trim(message)//'data vector is not allocated'; return; end if
        call cloneStruc(tempVec_rkind, ix_lower, source=dataStruct%var(iVar)%dat, err=err, message=cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
        deallocate(dataStruct%var(iVar)%dat); allocate(dataStruct%var(iVar)%dat(ix_lower:ix_upper-1),stat=err)
        if(err/=0)then; err=20; message=trim(message)//'problem reallocating data vector'; return; end if
        dataStruct%var(iVar)%dat(ix_lower:ix_rm-1) = tempVec_rkind(ix_lower:ix_rm-1)
        if(ix_upper > ix_rm) dataStruct%var(iVar)%dat(ix_rm:ix_upper-1) = tempVec_rkind(ix_rm+1:ix_upper)
        deallocate(tempVec_rkind)
      type is (var_ilength)
        if(.not.allocated(dataStruct%var(iVar)%dat))then; err=20; message=trim(message)//'data vector is not allocated'; return; end if
        call cloneStruc(tempVec_i4b, ix_lower, source=dataStruct%var(iVar)%dat, err=err, message=cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
        deallocate(dataStruct%var(iVar)%dat); allocate(dataStruct%var(iVar)%dat(ix_lower:ix_upper-1),stat=err)
        if(err/=0)then; err=20; message=trim(message)//'problem reallocating data vector'; return; end if
        dataStruct%var(iVar)%dat(ix_lower:ix_rm-1) = tempVec_i4b(ix_lower:ix_rm-1)
        if(ix_upper > ix_rm) dataStruct%var(iVar)%dat(ix_rm:ix_upper-1) = tempVec_i4b(ix_rm+1:ix_upper)
        deallocate(tempVec_i4b)
      class default; err=20; message=trim(message)//'unable to identify the data structure type'; return
    end select
  end do
end subroutine rmLakeLayer

end module lakeIceCover_module
