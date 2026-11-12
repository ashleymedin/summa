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

module wetlandArea_module
! The area of a wetland from the liquid water it holds, for a wetland with wetlandMaxArea > 0.
! Filled to its outlet the wetland covers its full area at a mean depth of wetlandMinDepth. Below the
! outlet its stage follows area**wetlandShape (0.5 a cone, 1 a paraboloid, large a pan); above it the
! area grows by wetlandFlatArea for each metre of rise. The full area is wetlandMaxArea, or the part
! of the wetlandMask cells of the wetland's HRU that the glacier has left free of ice if that is
! smaller, so a wetland with no area opens when the ice leaves its basin. Sediment it traps fills the basin:
! the full area and depth shrink with the filled fraction f as (1-f)**wetlandInfillShape and (1-f)**(1-wetlandInfillShape).

USE nr_type
USE data_types,only:var_ilength,var_dlength          ! data vectors with variable length dimension
USE data_types,only:grid_double,grid_info            ! basin grids
USE multiconst,only:iden_ice,iden_water,LH_fus       ! physical constants
USE globalData,only:icefrz_mult                      ! freezing curve scaling factor multiplier of snow to ice
USE globalData,only:thick4area                       ! ice thickness that counts as glacier area (m)
USE var_lookup,only:iLookPROG,iLookDIAG,iLookINDEX,iLookPARAM,iLookGRID ! named variables
implicit none
private
public::lakeImposeLiquid
public::wetlandFullArea
public::wetlandInfill
public::wetlandUpdateArea
real(rkind),parameter,public :: wetMinLiquid=0.01_rkind ! water depth a wetland keeps, liquid or ice (m)
real(rkind),parameter        :: minAreaFrac=0.01_rkind  ! smallest area of an open wetland, as a fraction of wetlandMaxArea (-)
real(rkind),parameter        :: areaTol=0.01_rkind      ! relative change of area that resizes the wetland (-)
contains

! ************************************************************************************************
! public subroutine lakeImposeLiquid: set the liquid depth of the water layers of a stream or wetland
! ************************************************************************************************
! The layers keep their temperature and their ice, so the enthalpy per unit volume follows the new fractions.
subroutine lakeImposeLiquid(liqDepth,mpar_data,indx_data,prog_data,diag_data,err,message)
  USE var_derive_module,only:calcHeight                  ! height at layer interfaces and layer mid-point
  USE snowLakeGlceDepth_module,only:lakePrescribeDepth   ! liquid depth of the lake layers
  USE convertEnthalpyTemp_module,only:T2enthTemp_snLaGl  ! convert temperature to enthalpy for snow, lake, and ice
  implicit none
  real(rkind),intent(in)           :: liqDepth   ! liquid depth to impose (m)
  type(var_dlength),intent(in)     :: mpar_data  ! model parameters
  type(var_ilength),intent(in)     :: indx_data  ! state vector geometry
  type(var_dlength),intent(inout)  :: prog_data  ! prognostic variables
  type(var_dlength),intent(inout)  :: diag_data  ! diagnostic variables
  integer(i4b),intent(out)         :: err        ! error code
  character(*),intent(out)         :: message    ! error message
  ! local variables
  character(len=256)               :: cmessage   ! error message of downwind routine
  integer(i4b)                     :: iLayer     ! layer index
  integer(i4b)                     :: nAbove     ! number of layers above the water layers
  integer(i4b)                     :: nLake      ! number of lake layers
  err=0; message="lakeImposeLiquid/"
  nAbove = indx_data%var(iLookINDEX%nSnow)%dat(1) + indx_data%var(iLookINDEX%nLakeFrz)%dat(1)
  nLake  = indx_data%var(iLookINDEX%nLake)%dat(1)
  call lakePrescribeDepth(&
                  nAbove,nLake-indx_data%var(iLookINDEX%nLakeFrz)%dat(1),                 & ! intent(in):    number of layers above the water layers, and of water layers
                  liqDepth,                                                               & ! intent(in):    liquid depth to impose (m)
                  prog_data%var(iLookPROG%mLayerDepth)%dat,                               & ! intent(inout): depth of each layer (m)
                  prog_data%var(iLookPROG%mLayerVolFracLiq)%dat,                          & ! intent(inout): volumetric fraction of liquid water (-)
                  prog_data%var(iLookPROG%mLayerVolFracIce)%dat,                          & ! intent(inout): volumetric fraction of ice (-)
                  err,cmessage)                                                             ! intent(out):   error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
  call calcHeight(indx_data,prog_data,err,cmessage) ! layer heights follow the new depths
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
  do iLayer=nAbove+1,indx_data%var(iLookINDEX%nSnow)%dat(1)+nLake ! the enthalpy per unit volume of the water layers follows the new liquid and ice fractions
    prog_data%var(iLookPROG%mLayerVolFracWat)%dat(iLayer) = prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(iLayer) &
                                                            + prog_data%var(iLookPROG%mLayerVolFracIce)%dat(iLayer)*(iden_ice/iden_water)
    call T2enthTemp_snLaGl(&
                  .false.,                                                                & ! intent(in):  flag that no liquid water in layer
                  mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)*icefrz_mult,             & ! intent(in):  scaling parameter for the lake freezing curve (K-1)
                  prog_data%var(iLookPROG%mLayerTemp)%dat(iLayer),                        & ! intent(in):  layer temperature (K)
                  prog_data%var(iLookPROG%mLayerVolFracWat)%dat(iLayer),                  & ! intent(in):  volumetric total water content (-)
                  diag_data%var(iLookDIAG%mLayerEnthTemp)%dat(iLayer))                      ! intent(out): temperature component of enthalpy (J m-3)
    prog_data%var(iLookPROG%mLayerEnthalpy)%dat(iLayer) = diag_data%var(iLookDIAG%mLayerEnthTemp)%dat(iLayer) &
                                                          - iden_ice*LH_fus*prog_data%var(iLookPROG%mLayerVolFracIce)%dat(iLayer)
  end do
end subroutine lakeImposeLiquid

! ************************************************************************************************
! public function wetlandFullArea: area of a wetland filled to its outlet (m2)
! ************************************************************************************************
! maxArea, or the ice-free part of the wetlandMask cells of HRU ixHRU if the HRU has any and that is smaller
function wetlandFullArea(ixHRU,maxArea,gridInfo,gridData) result(fullArea)
  implicit none
  integer(i4b),intent(in)          :: ixHRU       ! index of the wetland's HRU in its GRU
  real(rkind),intent(in)           :: maxArea     ! wetlandMaxArea (m2)
  type(grid_info),intent(in)       :: gridInfo(:) ! size and spacing of each grid of the GRU
  type(grid_double),intent(in)     :: gridData    ! grids of the GRU
  real(rkind)                      :: fullArea    ! area of the wetland filled to its outlet (m2)
  ! local variables
  integer(i4b)                     :: iGrid       ! grid index
  integer(i4b)                     :: nBasin      ! number of mask cells in the HRU
  real(rkind)                      :: freeArea    ! area of those cells that is free of ice (m2)
  fullArea = maxArea
  if(.not.allocated(gridData%grid)) return
  nBasin = 0; freeArea = 0._rkind
  do iGrid=1,size(gridInfo)
    associate(nx => gridInfo(iGrid)%nx, ny => gridInfo(iGrid)%ny, var => gridData%grid(iGrid)%var)
      associate(basin   => var(iLookGRID%wetlandMask)%dat2(1:nx,1:ny) > 0.5_rkind .and. nint(var(iLookGRID%cell2hru)%dat2(1:nx,1:ny)) == ixHRU, &
                iceThick => var(iLookGRID%surface_elev)%dat2(1:nx,1:ny) - var(iLookGRID%debris_thick)%dat2(1:nx,1:ny) - var(iLookGRID%bed_elev)%dat2(1:nx,1:ny))
        nBasin   = nBasin + count(basin)
        freeArea = freeArea + gridInfo(iGrid)%dx*gridInfo(iGrid)%dy*count(basin .and. .not.(iceThick > thick4area)) ! a cell with no bed has no ice
      end associate
    end associate
  end do
  if(nBasin > 0) fullArea = min(maxArea, freeArea)
end function wetlandFullArea

! ************************************************************************************************
! public subroutine wetlandInfill: full area and depth of a wetland whose basin holds sediment
! ************************************************************************************************
subroutine wetlandInfill(basinArea,mpar_data,prog_data,diag_data)
  implicit none
  real(rkind),intent(in)           :: basinArea  ! area of the basin filled to its outlet with no sediment (m2)
  type(var_dlength),intent(in)     :: mpar_data  ! model parameters
  type(var_dlength),intent(in)     :: prog_data  ! prognostic variables
  type(var_dlength),intent(inout)  :: diag_data  ! diagnostic variables
  real(rkind)                      :: basinVol   ! water the basin holds to its outlet with no sediment (m3)
  real(rkind)                      :: water      ! fraction of that volume not filled with sediment (-)
  associate(depth0 => mpar_data%var(iLookPARAM%wetlandMinDepth)%dat(1), shapeA => mpar_data%var(iLookPARAM%wetlandInfillShape)%dat(1))
    basinVol = basinArea*depth0
    water = 1._rkind
    if(basinVol > 0._rkind) water = max(1._rkind - prog_data%var(iLookPROG%scalarWetlandSediment)%dat(1)/basinVol, 0._rkind)
    diag_data%var(iLookDIAG%scalarWetlandFullArea)%dat(1)  = basinArea*water**shapeA
    diag_data%var(iLookDIAG%scalarWetlandFullDepth)%dat(1) = depth0*water**(1._rkind - shapeA)
  end associate
end subroutine wetlandInfill

! ************************************************************************************************
! public subroutine wetlandUpdateArea: set the area of a wetland from the liquid water it holds
! ************************************************************************************************
! The water volume and temperature are kept, so the liquid depth changes with the area. A wetland with snow or ice
! keeps its area. A wetland with no area opens at its smallest area, with wetMinLiquid of water, once its full area is not zero.
subroutine wetlandUpdateArea(fullArea,availArea,mpar_data,indx_data,prog_data,diag_data,changed,err,message)
  implicit none
  real(rkind),intent(in)           :: fullArea   ! area of the wetland filled to its outlet (m2)
  real(rkind),intent(in)           :: availArea  ! area of the HRU the wetland can take (m2)
  type(var_dlength),intent(in)     :: mpar_data  ! model parameters
  type(var_ilength),intent(in)     :: indx_data  ! state vector geometry
  type(var_dlength),intent(inout)  :: prog_data  ! prognostic variables
  type(var_dlength),intent(inout)  :: diag_data  ! diagnostic variables
  logical(lgt),intent(out)         :: changed    ! flag that the area changed
  integer(i4b),intent(out)         :: err        ! error code
  character(*),intent(out)         :: message    ! error message
  ! local variables
  character(len=256)               :: cmessage   ! error message of downwind routine
  integer(i4b)                     :: nAbove     ! number of layers above the water layers
  integer(i4b)                     :: nBelow     ! last water layer
  real(rkind)                      :: minArea    ! smallest area of the open wetland (m2)
  real(rkind)                      :: volume     ! liquid water in the wetland (m3)
  real(rkind)                      :: volFull    ! liquid water in the wetland filled to its outlet (m3)
  real(rkind)                      :: areaNew    ! area the volume fills (m2)
  real(rkind)                      :: liqDepth   ! liquid depth at that area (m)
  err=0; message="wetlandUpdateArea/"
  changed = .false.
  associate(DOMarea     => prog_data%var(iLookPROG%DOMarea)%dat(1),             &
            mLayerDepth => prog_data%var(iLookPROG%mLayerDepth)%dat,            &
            volFracLiq  => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat,       &
            volFracIce  => prog_data%var(iLookPROG%mLayerVolFracIce)%dat,       &
            maxArea     => mpar_data%var(iLookPARAM%wetlandMaxArea)%dat(1),     &
            fullDepth   => diag_data%var(iLookDIAG%scalarWetlandFullDepth)%dat(1), &
            stagePower  => mpar_data%var(iLookPARAM%wetlandShape)%dat(1),       &
            flatArea    => mpar_data%var(iLookPARAM%wetlandFlatArea)%dat(1))
    nAbove  = indx_data%var(iLookINDEX%nSnow)%dat(1) + indx_data%var(iLookINDEX%nLakeFrz)%dat(1)
    nBelow  = indx_data%var(iLookINDEX%nSnow)%dat(1) + indx_data%var(iLookINDEX%nLake)%dat(1)
    minArea = min(minAreaFrac*maxArea, availArea)
    if(DOMarea <= 0._rkind)then
      if(fullArea <= 0._rkind .or. minArea <= 0._rkind) return ! the basin is still under the ice
      if(indx_data%var(iLookINDEX%nLake)%dat(1) < 1)then; err=20; message=trim(message)//'a wetland that opens needs lake layers in its initial conditions'; return; end if
      areaNew  = minArea
      liqDepth = wetMinLiquid
    else
      if(nAbove > 0 .or. prog_data%var(iLookPROG%scalarSWE)%dat(1) > 0._rkind .or. any(volFracIce(nAbove+1:nBelow) > 0._rkind)) return
      volume  = DOMarea*sum(mLayerDepth(nAbove+1:nBelow)*volFracLiq(nAbove+1:nBelow))
      volFull = fullArea*fullDepth
      if(volume <= volFull)then
        areaNew = fullArea*(volume/max(volFull, tiny(volFull)))**(1._rkind/(1._rkind + stagePower))
      else
        areaNew = sqrt(fullArea**2 + 2._rkind*flatArea*(volume - volFull))
      end if
      areaNew = min(areaNew, volume/wetMinLiquid) ! no shallower than the liquid the wetland keeps
      areaNew = min(max(areaNew, minArea), availArea)
      if(abs(areaNew - DOMarea) <= areaTol*DOMarea) return
      liqDepth = max(volume/areaNew, wetMinLiquid)
    end if
    call lakeImposeLiquid(liqDepth,mpar_data,indx_data,prog_data,diag_data,err,cmessage)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
    DOMarea = areaNew
    changed = .true.
  end associate
end subroutine wetlandUpdateArea

end module wetlandArea_module
