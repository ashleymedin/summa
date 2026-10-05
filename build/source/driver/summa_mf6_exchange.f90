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

module summa_mf6_exchange
  ! ****************************************************************************************
  ! *** The SUMMA side of the SUMMA <-> MODFLOW 6 coupling                                ***
  ! ****************************************************************************************
  !
  ! Per-HRU getters and setters for the handful of SUMMA quantities the MODFLOW 6 coupler
  ! exchanges, plus the HRU geometry it needs to build its cell map.  They take the top-level
  ! SUMMA structure directly, so both routes into the coupling can use them:
  !
  !   summa_bmi.f90         the BMI wrapper, whose get_value/set_value for the coupler's
  !                         variable names delegate here (used by the standalone couplers)
  !   summa_simulation.f90  the calibration driver, which runs SUMMA without the BMI and
  !                         so calls these directly
  !
  ! HRU order is the same everywhere: GRUs in order, HRUs within each GRU in order, which is
  ! BMI grid 0 and the order the MODFLOW cell map is built in.
  !
  ! NOTE on the flat HRU index: i = (iGRU-1)*gru_struc(iGRU)%hruCount + jHRU is the indexing
  !       SUMMA's BMI has always used, and it is kept verbatim here so the two agree.  It is
  !       only a contiguous numbering when every GRU has the same HRU count; changing it would
  !       silently move every coupled run's HRU->cell mapping, so it is left as it is.

  USE nr_type,    only: i4b, rkind
  USE globalData, only: realMissing  ! a domain that never solved a soil column
  USE globalData, only: verySmall    ! a small number
  USE globalData, only: model_decisions  ! model decision structure
  USE var_lookup, only: iLookDECISIONS   ! named variables for elements of the decision structure
  USE mDecisions_module, only: homegrown_SE ! homegrown saturation excess surface runoff
  USE multiconst, only: iden_water     ! intrinsic density of liquid water (kg m-3)
  USE multiconst, only: Tfreeze        ! freezing point of pure water (K)
  USE summa_type, only: summa1_type_dec

  USE globalData, only: gru_struc            ! HRU information for given GRU
  USE globalData, only: mfAquiferBaseflow    ! MODFLOW 6 coupler aquifer-baseflow feedback channel
  USE globalData, only: mfSurfaceDischarge   ! MODFLOW 6 coupler groundwater-discharge-at-surface channel
  USE globalData, only: mfAquiferTranspire   ! MODFLOW 6 coupler groundwater-ET feedback channel

  USE var_lookup, only: iLookATTR            ! named variables for real valued attribute data structure
  USE var_lookup, only: iLookINDEX           ! named variables for local model indices
  USE var_lookup, only: iLookPROG            ! named variables for local prognostic variables
  USE var_lookup, only: iLookPARAM           ! named variables for local model parameters
  USE var_lookup, only: iLookFLUX            ! named variables for local flux variables
  USE var_lookup, only: iLookDIAG            ! named variables for local diagnostic variables
  USE var_lookup, only: iLookBVAR            ! named variables for basin (GRU) variables

  implicit none
  private

  public :: mf6x_hru_count
  public :: mf6x_hru_longitude
  public :: mf6x_hru_latitude
  public :: mf6x_hru_elevation
  public :: mf6x_hru_area
  public :: mf6x_soil_thickness
  public :: mf6x_root_reach
  public :: mf6x_get_drainage
  public :: mf6x_put_lower_bound_head
  public :: mf6x_put_aquifer_storage
  public :: mf6x_put_aquifer_baseflow
  public :: mf6x_put_surface_discharge
  public :: mf6x_get_aquifer_transpire
  public :: mf6x_put_aquifer_transpire
  public :: mf6x_put_transpire_lim_aqfr
  public :: mf6x_put_infil_lim_aqfr
  public :: mf6x_put_aquifer_reject
  public :: mf6x_root_zone_depth
  public :: mf6x_get_drainage_temp
  public :: mf6x_get_base_nrg_flux
  public :: mf6x_put_aquifer_temp

contains

  ! **************************************************************************************************
  ! Number of HRUs in the run domain of this SUMMA instance (BMI grid 0 size).
  ! **************************************************************************************************
  integer(i4b) function mf6x_hru_count() result(nHRU)
    nHRU = sum(gru_struc(:)%hruCount)
  end function mf6x_hru_count

  ! **************************************************************************************************
  ! HRU centroid longitude (degrees east).
  ! **************************************************************************************************
  subroutine mf6x_hru_longitude(summa_struct, x)
    type(summa1_type_dec), intent(in)  :: summa_struct
    double precision,      intent(out) :: x(:)
    integer(i4b) :: iGRU, jHRU
    associate(attrStruct => summa_struct%attrStruct)   ! x%gru(:)%hru(:)%var(:)
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          x(gru_struc(iGRU)%hruInfo(jHRU)%hru_ix) = &
            attrStruct%gru(iGRU)%hru(jHRU)%var(iLookATTR%longitude)
        end do
      end do
    end associate
  end subroutine mf6x_hru_longitude

  ! **************************************************************************************************
  ! HRU centroid latitude (degrees north).
  ! **************************************************************************************************
  subroutine mf6x_hru_latitude(summa_struct, y)
    type(summa1_type_dec), intent(in)  :: summa_struct
    double precision,      intent(out) :: y(:)
    integer(i4b) :: iGRU, jHRU
    associate(attrStruct => summa_struct%attrStruct)
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          y(gru_struc(iGRU)%hruInfo(jHRU)%hru_ix) = &
            attrStruct%gru(iGRU)%hru(jHRU)%var(iLookATTR%latitude)
        end do
      end do
    end associate
  end subroutine mf6x_hru_latitude

  ! **************************************************************************************************
  ! HRU land-surface elevation (m).
  ! **************************************************************************************************
  subroutine mf6x_hru_elevation(summa_struct, z)
    type(summa1_type_dec), intent(in)  :: summa_struct
    double precision,      intent(out) :: z(:)
    integer(i4b) :: iGRU, jHRU
    associate(attrStruct => summa_struct%attrStruct)
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          z(gru_struc(iGRU)%hruInfo(jHRU)%hru_ix) = &
            attrStruct%gru(iGRU)%hru(jHRU)%var(iLookATTR%elevation)
        end do
      end do
    end associate
  end subroutine mf6x_hru_elevation

  ! **************************************************************************************************
  ! HRU plan area (m2) from attributes.nc, for the coupler's area check and coupled budget.
  ! **************************************************************************************************
  subroutine mf6x_hru_area(summa_struct, area)
    type(summa1_type_dec), intent(in)  :: summa_struct
    double precision,      intent(out) :: area(:)
    integer(i4b) :: iGRU, jHRU
    associate(attrStruct => summa_struct%attrStruct)
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          area(gru_struc(iGRU)%hruInfo(jHRU)%hru_ix) = &
            attrStruct%gru(iGRU)%hru(jHRU)%var(iLookATTR%HRUarea)
        end do
      end do
    end associate
  end subroutine mf6x_hru_area

  ! **************************************************************************************************
  ! How far roots reach below the base of the soil column (m): rootingDepth - soil depth, floored at zero.
  ! **************************************************************************************************
  subroutine mf6x_root_reach(summa_struct, reach)
    type(summa1_type_dec), intent(in)  :: summa_struct
    double precision,      intent(out) :: reach(:)
    integer(i4b) :: iGRU, jHRU, iDOM, i, ixDOM, nSnow, nLake, nSoil
    real(rkind)  :: soilDepth
    associate(progStruct => summa_struct%progStruct, &
              indxStruct => summa_struct%indxStruct, &
              mparStruct => summa_struct%mparStruct)
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          i = gru_struc(iGRU)%hruInfo(jHRU)%hru_ix
          ixDOM = 1
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            if (indxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookINDEX%nGlce)%dat(1) == 0) then
              ixDOM = iDOM; exit
            end if
          end do
          nSnow = indxStruct%gru(iGRU)%hru(jHRU)%dom(ixDOM)%var(iLookINDEX%nSnow)%dat(1)
          nLake = indxStruct%gru(iGRU)%hru(jHRU)%dom(ixDOM)%var(iLookINDEX%nLake)%dat(1)
          nSoil = indxStruct%gru(iGRU)%hru(jHRU)%dom(ixDOM)%var(iLookINDEX%nSoil)%dat(1)
          soilDepth = progStruct%gru(iGRU)%hru(jHRU)%dom(ixDOM)%var(iLookPROG%iLayerHeight)%dat(nSnow+nLake+nSoil) &
                    - progStruct%gru(iGRU)%hru(jHRU)%dom(ixDOM)%var(iLookPROG%iLayerHeight)%dat(nSnow+nLake)
          reach(i) = max(mparStruct%gru(iGRU)%hru(jHRU)%dom(ixDOM)%var(iLookPARAM%rootingDepth)%dat(1) - soilDepth, 0._rkind)
        end do
      end do
    end associate
  end subroutine mf6x_root_reach

  ! **************************************************************************************************
  ! Depth below the surface of the zone whose positive pressure closes the infiltrating area (m), as
  ! soilLiqFlux takes it: the base of the deepest root layer under homegrown_SE, else the soil depth.
  ! **************************************************************************************************
  subroutine mf6x_root_zone_depth(summa_struct, depth)
    type(summa1_type_dec), intent(in)  :: summa_struct
    double precision,      intent(out) :: depth(:)
    integer(i4b) :: iGRU, jHRU, iDOM, i, ixDOM, nSnow, nLake, nSoil, nRoots
    real(rkind)  :: rootingDepth
    associate(progStruct => summa_struct%progStruct, &
              indxStruct => summa_struct%indxStruct, &
              mparStruct => summa_struct%mparStruct)
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          i = gru_struc(iGRU)%hruInfo(jHRU)%hru_ix
          ixDOM = 1
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            if (indxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookINDEX%nGlce)%dat(1) == 0) then
              ixDOM = iDOM; exit
            end if
          end do
          nSnow = indxStruct%gru(iGRU)%hru(jHRU)%dom(ixDOM)%var(iLookINDEX%nSnow)%dat(1)
          nLake = indxStruct%gru(iGRU)%hru(jHRU)%dom(ixDOM)%var(iLookINDEX%nLake)%dat(1)
          nSoil = indxStruct%gru(iGRU)%hru(jHRU)%dom(ixDOM)%var(iLookINDEX%nSoil)%dat(1)
          associate(h => progStruct%gru(iGRU)%hru(jHRU)%dom(ixDOM)%var(iLookPROG%iLayerHeight)%dat)
            nRoots = nSoil
            if (model_decisions(iLookDECISIONS%surfRun_SE)%iDecision == homegrown_SE) then
              rootingDepth = mparStruct%gru(iGRU)%hru(jHRU)%dom(ixDOM)%var(iLookPARAM%rootingDepth)%dat(1)
              nRoots = max(1, count(h(nSnow+nLake:nSnow+nLake+nSoil-1) - h(nSnow+nLake) < rootingDepth-verySmall))
            end if
            depth(i) = h(nSnow+nLake+nRoots) - h(nSnow+nLake)
          end associate
        end do
      end do
    end associate
  end subroutine mf6x_root_zone_depth

  ! **************************************************************************************************
  ! Thickness of the SUMMA soil column for each HRU (m), measured from the ground surface
  ! (iLayerHeight = 0) down to the base of the lowest soil layer.  The coupler forms
  ! lowerBoundHead from it, in place of a hard-coded soil thickness.
  ! **************************************************************************************************
  subroutine mf6x_soil_thickness(summa_struct, thickness)
    type(summa1_type_dec), intent(in)  :: summa_struct
    double precision,      intent(out) :: thickness(:)
    integer(i4b) :: iGRU, jHRU, iDOM, i, ixDOM, nSnow, nLake, nSoil
    associate(progStruct => summa_struct%progStruct, &   ! x%gru(:)%hru(:)%dom(:)%var(:)%dat
              indxStruct => summa_struct%indxStruct)     ! x%gru(:)%hru(:)%dom(:)%var(:)%dat
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          i = gru_struc(iGRU)%hruInfo(jHRU)%hru_ix
          ! prefer the first non-glacier domain; fall back to domain 1
          ixDOM = 1
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            if (indxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookINDEX%nGlce)%dat(1) == 0) then
              ixDOM = iDOM; exit
            end if
          end do
          nSnow = indxStruct%gru(iGRU)%hru(jHRU)%dom(ixDOM)%var(iLookINDEX%nSnow)%dat(1)
          nLake = indxStruct%gru(iGRU)%hru(jHRU)%dom(ixDOM)%var(iLookINDEX%nLake)%dat(1)
          nSoil = indxStruct%gru(iGRU)%hru(jHRU)%dom(ixDOM)%var(iLookINDEX%nSoil)%dat(1)
          thickness(i) = progStruct%gru(iGRU)%hru(jHRU)%dom(ixDOM)%var(iLookPROG%iLayerHeight)%dat(nSnow+nLake+nSoil) &
                       - progStruct%gru(iGRU)%hru(jHRU)%dom(ixDOM)%var(iLookPROG%iLayerHeight)%dat(nSnow+nLake)
        end do
      end do
    end associate
  end subroutine mf6x_soil_thickness

  ! **************************************************************************************************
  ! Drainage out the base of the soil column (m s-1, "scalarSoilDrainage"), averaged over the
  ! HRU's domains by area fraction.  This is the recharge the coupler imposes on MODFLOW 6.
  ! **************************************************************************************************
  subroutine mf6x_get_drainage(summa_struct, drainage)
    type(summa1_type_dec), intent(in)  :: summa_struct
    real,                  intent(out) :: drainage(:)
    integer(i4b) :: iGRU, jHRU, iDOM, i
    real         :: fracDOM
    real(rkind)  :: hruArea, drainDOM
    associate(progStruct => summa_struct%progStruct, &
              fluxStruct => summa_struct%fluxStruct)
      drainage = -999.0
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          i = gru_struc(iGRU)%hruInfo(jHRU)%hru_ix
          drainage(i) = 0._rkind
          ! the flux handed back is per unit of THIS HRU, so the domain weights sum to one over the HRU
          hruArea = 0._rkind
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            hruArea = hruArea + progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)
          end do
          if(hruArea <= 0._rkind) cycle
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            drainDOM = fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarSoilDrainage)%dat(1)
            ! a stream reach has no soil column of its own, so it drains nothing to the aquifer
            if(drainDOM <= realMissing) cycle
            fracDOM = real(progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1) / hruArea)
            drainage(i) = drainage(i) + drainDOM * fracDOM
          end do
        end do
      end do
    end associate
  end subroutine mf6x_get_drainage

  ! **************************************************************************************************
  ! Prescribed-head lower boundary condition for soil hydrology ("lowerBoundHead", m), supplied
  ! by the coupled MODFLOW 6 water table.  Glacier domains keep their own boundary condition.
  ! **************************************************************************************************
  subroutine mf6x_put_lower_bound_head(summa_struct, head)
    type(summa1_type_dec), intent(inout) :: summa_struct
    real,                  intent(in)    :: head(:)
    integer(i4b) :: iGRU, jHRU, iDOM, i
    associate(mparStruct => summa_struct%mparStruct, &
              indxStruct => summa_struct%indxStruct)
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          i = gru_struc(iGRU)%hruInfo(jHRU)%hru_ix
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            if (indxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookINDEX%nGlce)%dat(1) == 0) &
              mparStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPARAM%lowerBoundHead)%dat(1) = head(i)
          end do
        end do
      end do
    end associate
  end subroutine mf6x_put_lower_bound_head

  ! **************************************************************************************************
  ! Aquifer storage from the coupled MODFLOW 6 model ("scalarAquiferStorage", m).
  ! **************************************************************************************************
  subroutine mf6x_put_aquifer_storage(summa_struct, storage)
    type(summa1_type_dec), intent(inout) :: summa_struct
    real,                  intent(in)    :: storage(:)
    integer(i4b) :: iGRU, jHRU, iDOM, i
    associate(progStruct => summa_struct%progStruct, &
              indxStruct => summa_struct%indxStruct)
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          i = gru_struc(iGRU)%hruInfo(jHRU)%hru_ix
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            if (indxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookINDEX%nGlce)%dat(1) == 0) &
              progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPROG%scalarAquiferStorage)%dat(1) = storage(i)
          end do
        end do
      end do
    end associate
  end subroutine mf6x_put_aquifer_storage

  ! **************************************************************************************************
  ! Aquifer baseflow from the coupled MODFLOW 6 model ("scalarAquiferBaseflow", m s-1).
  !
  ! SUMMA overwrites every scalar flux during a step, so this cannot be written straight into
  ! fluxStruct; it goes into the mfAquiferBaseflow channel that the physics reads instead.
  ! **************************************************************************************************
  subroutine mf6x_put_aquifer_baseflow(summa_struct, baseflow)
    type(summa1_type_dec), intent(inout) :: summa_struct
    real,                  intent(in)    :: baseflow(:)
    integer(i4b) :: iGRU, jHRU, i
    if (.not. allocated(mfAquiferBaseflow)) then
      allocate(mfAquiferBaseflow(sum(gru_struc(:)%hruCount))); mfAquiferBaseflow = 0._rkind
    end if
    do iGRU = 1, summa_struct%nGRU_local
      do jHRU = 1, gru_struc(iGRU)%hruCount
        i = gru_struc(iGRU)%hruInfo(jHRU)%hru_ix
        mfAquiferBaseflow(i) = baseflow(i)
      end do
    end do
  end subroutine mf6x_put_aquifer_baseflow

  ! **************************************************************************************************
  ! Groundwater discharge at land surface from MODFLOW (m s-1, + = out of aquifer), added to SUMMA's
  ! surface runoff.  Uses the mfSurfaceDischarge channel, since fluxStruct scalars do not survive a step.
  ! **************************************************************************************************
  subroutine mf6x_put_surface_discharge(summa_struct, discharge)
    type(summa1_type_dec), intent(inout) :: summa_struct
    real,                  intent(in)    :: discharge(:)
    integer(i4b) :: iGRU, jHRU, i
    if (.not. allocated(mfSurfaceDischarge)) then
      allocate(mfSurfaceDischarge(sum(gru_struc(:)%hruCount))); mfSurfaceDischarge = 0._rkind
    end if
    do iGRU = 1, summa_struct%nGRU_local
      do jHRU = 1, gru_struc(iGRU)%hruCount
        i = gru_struc(iGRU)%hruInfo(jHRU)%hru_ix
        mfSurfaceDischarge(i) = discharge(i)
      end do
    end do
  end subroutine mf6x_put_surface_discharge

  ! **************************************************************************************************
  ! Groundwater evapotranspiration actually taken by MODFLOW (m s-1, + = out of aquifer).
  ! **************************************************************************************************
  subroutine mf6x_put_aquifer_transpire(summa_struct, transpire)
    type(summa1_type_dec), intent(inout) :: summa_struct
    real,                  intent(in)    :: transpire(:)
    integer(i4b) :: iGRU, jHRU, i
    if (.not. allocated(mfAquiferTranspire)) then
      allocate(mfAquiferTranspire(sum(gru_struc(:)%hruCount))); mfAquiferTranspire = 0._rkind
    end if
    do iGRU = 1, summa_struct%nGRU_local
      do jHRU = 1, gru_struc(iGRU)%hruCount
        i = gru_struc(iGRU)%hruInfo(jHRU)%hru_ix
        mfAquiferTranspire(i) = transpire(i)
      end do
    end do
  end subroutine mf6x_put_aquifer_transpire

  ! **************************************************************************************************
  ! Aquifer transpiration limiting factor (-) from the coupler, evaluated per MODFLOW cell.
  ! Written into diag, where soilResist reads it back as an input.
  ! **************************************************************************************************
  subroutine mf6x_put_transpire_lim_aqfr(summa_struct, limit)
    type(summa1_type_dec), intent(inout) :: summa_struct
    real,                  intent(in)    :: limit(:)
    integer(i4b) :: iGRU, jHRU, iDOM, i
    associate(diagStruct => summa_struct%diagStruct, &
              indxStruct => summa_struct%indxStruct)
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          i = gru_struc(iGRU)%hruInfo(jHRU)%hru_ix
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            if (indxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookINDEX%nGlce)%dat(1) == 0) &
              diagStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookDIAG%scalarTranspireLimAqfr)%dat(1) = limit(i)
          end do
        end do
      end do
    end associate
  end subroutine mf6x_put_transpire_lim_aqfr

  ! **************************************************************************************************
  ! Aquifer control on the infiltrating area (-) from the coupler, evaluated per MODFLOW cell.
  ! Written into diag, where soilLiqFlux uses it in place of the column's own compression closure.
  ! **************************************************************************************************
  subroutine mf6x_put_infil_lim_aqfr(summa_struct, limit)
    type(summa1_type_dec), intent(inout) :: summa_struct
    real,                  intent(in)    :: limit(:)
    integer(i4b) :: iGRU, jHRU, iDOM, i
    associate(diagStruct => summa_struct%diagStruct, &
              indxStruct => summa_struct%indxStruct)
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          i = gru_struc(iGRU)%hruInfo(jHRU)%hru_ix
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            if (indxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookINDEX%nGlce)%dat(1) == 0) &
              diagStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookDIAG%scalarInfilLimAqfr)%dat(1) = limit(i)
          end do
        end do
      end do
    end associate
  end subroutine mf6x_put_infil_lim_aqfr

  ! **************************************************************************************************
  ! Recharge the unsaturated zone below the soil column rejected (m s-1, + = back into the column).
  ! Written into diag, where the presHead lower boundary returns it to the column base.
  ! **************************************************************************************************
  subroutine mf6x_put_aquifer_reject(summa_struct, reject)
    type(summa1_type_dec), intent(inout) :: summa_struct
    real,                  intent(in)    :: reject(:)
    integer(i4b) :: iGRU, jHRU, iDOM, i
    associate(diagStruct => summa_struct%diagStruct, &
              indxStruct => summa_struct%indxStruct)
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          i = gru_struc(iGRU)%hruInfo(jHRU)%hru_ix
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            if (indxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookINDEX%nGlce)%dat(1) == 0) &
              diagStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookDIAG%scalarAquiferReject)%dat(1) = reject(i)
          end do
        end do
      end do
    end associate
  end subroutine mf6x_put_aquifer_reject

  ! **************************************************************************************************
  ! SUMMA's aquifer transpiration demand, per HRU (m s-1, + = out of aquifer): the aquifer's share of
  ! canopy transpiration, scalarAquiferRootFrac * scalarTranspireLimAqfr / scalarTranspireLim.
  ! **************************************************************************************************
  subroutine mf6x_get_aquifer_transpire(summa_struct, demand)
    type(summa1_type_dec), intent(in)  :: summa_struct
    real,                  intent(out) :: demand(:)
    integer(i4b) :: iGRU, jHRU, iDOM, i
    real(rkind)  :: fracDOM, frac, tlim, hruArea
    ! weighted as mf6x_get_drainage weights drainage
    associate(progStruct => summa_struct%progStruct, &
              diagStruct => summa_struct%diagStruct, &
              fluxStruct => summa_struct%fluxStruct)
      demand = 0.0
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          i = gru_struc(iGRU)%hruInfo(jHRU)%hru_ix
          demand(i) = 0._rkind
          hruArea = 0._rkind
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            hruArea = hruArea + progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)
          end do
          if(hruArea <= 0._rkind) cycle
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            tlim = diagStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookDIAG%scalarTranspireLim)%dat(1)
            if (tlim <= 0._rkind) cycle       ! no transpiration at all, so no aquifer share
            frac = diagStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookDIAG%scalarAquiferRootFrac)%dat(1) &
                 * diagStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookDIAG%scalarTranspireLimAqfr)%dat(1) / tlim
            if (frac <= 0._rkind) cycle       ! no roots below the soil column, or water table out of reach
            fracDOM = progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1) / hruArea
            ! negated: scalarCanopyTranspiration is negative for water leaving the canopy
            demand(i) = demand(i) &
                      - frac * fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarCanopyTranspiration)%dat(1) &
                        / iden_water * fracDOM
          end do
          ! a negative demand is condensation, which the aquifer plays no part in
          if (demand(i) < 0._rkind) demand(i) = 0._rkind
        end do
      end do
    end associate
  end subroutine mf6x_get_aquifer_transpire

  ! **************************************************************************************************
  ! Temperature of the soil drainage (K): the lowest soil layer, floored at freezing since the water
  ! leaves as liquid, averaged over the HRU's soil domains by area.  The recharge temperature for GWE.
  ! **************************************************************************************************
  subroutine mf6x_get_drainage_temp(summa_struct, temp)
    type(summa1_type_dec), intent(in)  :: summa_struct
    real,                  intent(out) :: temp(:)
    integer(i4b) :: iGRU, jHRU, iDOM, i, nLyr
    real(rkind)  :: tsum, asum, areaDOM
    associate(progStruct => summa_struct%progStruct, &
              fluxStruct => summa_struct%fluxStruct, &
              indxStruct => summa_struct%indxStruct)
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          i = gru_struc(iGRU)%hruInfo(jHRU)%hru_ix
          tsum = 0._rkind; asum = 0._rkind
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            ! a stream reach drains nothing to the aquifer, as in mf6x_get_drainage
            if(fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarSoilDrainage)%dat(1) <= realMissing) cycle
            if(indxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookINDEX%nSoil)%dat(1) < 1) cycle
            areaDOM = progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)
            if(areaDOM <= 0._rkind) cycle
            nLyr = indxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookINDEX%nLayers)%dat(1)
            tsum = tsum + areaDOM*max(progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPROG%mLayerTemp)%dat(nLyr), Tfreeze)
            asum = asum + areaDOM
          end do
          temp(i) = real(merge(tsum/asum, Tfreeze, asum > 0._rkind))
        end do
      end do
    end associate
  end subroutine mf6x_get_drainage_temp

  ! **************************************************************************************************
  ! Conduction out the base of the soil column (W m-2, + = down into the aquifer,
  ! "scalarLowerBoundNrgFlux"), per unit of the HRU; a stream reach conducts nothing to the aquifer.
  ! **************************************************************************************************
  subroutine mf6x_get_base_nrg_flux(summa_struct, flux)
    type(summa1_type_dec), intent(in)  :: summa_struct
    real,                  intent(out) :: flux(:)
    integer(i4b) :: iGRU, jHRU, iDOM, i
    real(rkind)  :: fsum, asum, areaDOM
    associate(progStruct => summa_struct%progStruct, &
              fluxStruct => summa_struct%fluxStruct)
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          i = gru_struc(iGRU)%hruInfo(jHRU)%hru_ix
          fsum = 0._rkind; asum = 0._rkind
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            areaDOM = progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)
            asum = asum + areaDOM
            if(fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarSoilDrainage)%dat(1) <= realMissing) cycle
            fsum = fsum + areaDOM*fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarLowerBoundNrgFlux)%dat(1)
          end do
          flux(i) = real(merge(fsum/asum, 0._rkind, asum > 0._rkind))
        end do
      end do
    end associate
  end subroutine mf6x_get_base_nrg_flux

  ! **************************************************************************************************
  ! Aquifer temperature from the coupled MODFLOW 6 GWE model ("scalarAquiferTemp", K); a value <= 0
  ! marks an HRU with no GWE cell.
  ! **************************************************************************************************
  subroutine mf6x_put_aquifer_temp(summa_struct, temp)
    type(summa1_type_dec), intent(inout) :: summa_struct
    real,                  intent(in)    :: temp(:)
    integer(i4b) :: iGRU, jHRU, iDOM, i
    associate(progStruct => summa_struct%progStruct, &
              indxStruct => summa_struct%indxStruct)
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          i = gru_struc(iGRU)%hruInfo(jHRU)%hru_ix
          if (temp(i) <= 0.0) cycle ! no GWE cell under this HRU, so it keeps its own
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            if (indxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookINDEX%nGlce)%dat(1) == 0) &
              progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPROG%scalarAquiferTemp)%dat(1) = temp(i)
          end do
        end do
      end do
    end associate
  end subroutine mf6x_put_aquifer_temp

end module summa_mf6_exchange
