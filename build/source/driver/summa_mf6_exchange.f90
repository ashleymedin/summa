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
  USE multiconst, only: iden_water     ! intrinsic density of liquid water (kg m-3)
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
          x((iGRU-1) * gru_struc(iGRU)%hruCount + jHRU) = &
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
          y((iGRU-1) * gru_struc(iGRU)%hruCount + jHRU) = &
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
          z((iGRU-1) * gru_struc(iGRU)%hruCount + jHRU) = &
            attrStruct%gru(iGRU)%hru(jHRU)%var(iLookATTR%elevation)
        end do
      end do
    end associate
  end subroutine mf6x_hru_elevation

  ! **************************************************************************************************
  ! HRU plan area (m2), as given in attributes.nc.
  !
  ! The coupler checks this against the summed plan area of the MODFLOW cells each HRU maps to.
  ! The scatter conserves recharge RATE, not VOLUME (it is a weight-weighted mean of HRU drainage
  ! rates), so volume is conserved only when the two areas agree - which is why the check exists.
  ! **************************************************************************************************
  subroutine mf6x_hru_area(summa_struct, area)
    type(summa1_type_dec), intent(in)  :: summa_struct
    double precision,      intent(out) :: area(:)
    integer(i4b) :: iGRU, jHRU
    associate(attrStruct => summa_struct%attrStruct)
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          area((iGRU-1) * gru_struc(iGRU)%hruCount + jHRU) = &
            attrStruct%gru(iGRU)%hru(jHRU)%var(iLookATTR%HRUarea)
        end do
      end do
    end associate
  end subroutine mf6x_hru_area

  ! **************************************************************************************************
  ! How far roots reach BELOW the base of the soil column, per HRU (m): rootingDepth - soil depth,
  ! floored at zero.  This is the span of soilResist's aquifer transpiration ramp, and the coupler
  ! uses it to decide how tight the HRU-elevation check has to be: an elevation offset smaller than
  ! this span can still move the limiting factor across its whole range.
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
          i = (iGRU-1) * gru_struc(iGRU)%hruCount + jHRU
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
          i = (iGRU-1) * gru_struc(iGRU)%hruCount + jHRU
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
    associate(progStruct => summa_struct%progStruct, &
              fluxStruct => summa_struct%fluxStruct, &
              bvarStruct => summa_struct%bvarStruct)
      drainage = -999.0
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          i = (iGRU-1) * gru_struc(iGRU)%hruCount + jHRU
          drainage(i) = 0._rkind
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            fracDOM = progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1) &
                    / bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)
            drainage(i) = drainage(i) &
                        + fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarSoilDrainage)%dat(1) * fracDOM
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
          i = (iGRU-1) * gru_struc(iGRU)%hruCount + jHRU
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
          i = (iGRU-1) * gru_struc(iGRU)%hruCount + jHRU
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
        i = (iGRU-1) * gru_struc(iGRU)%hruCount + jHRU
        mfAquiferBaseflow(i) = baseflow(i)
      end do
    end do
  end subroutine mf6x_put_aquifer_baseflow

  ! **************************************************************************************************
  ! Groundwater discharge at land surface from the coupled MODFLOW 6 model (m s-1, + = out of aquifer).
  !
  ! This is the water the aquifer cannot hold once the water table reaches land surface, taken from
  ! whichever boundary package carries role 'surface_discharge' (a DRN at DIS/TOP, as in the upstream
  ! Sagehen model).  It is added to SUMMA's surface runoff, not to the aquifer baseflow term, so it
  ! reaches routing as saturation-excess runoff rather than as baseflow.
  !
  ! Same channel reasoning as mf6x_put_aquifer_baseflow above: fluxStruct scalars do not survive a step.
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
        i = (iGRU-1) * gru_struc(iGRU)%hruCount + jHRU
        mfSurfaceDischarge(i) = discharge(i)
      end do
    end do
  end subroutine mf6x_put_surface_discharge

  ! **************************************************************************************************
  ! Groundwater evapotranspiration actually taken by MODFLOW 6 (m s-1, + = out of aquifer).
  !
  ! SUMMA sends a transpiration DEMAND down (see mf6x_get_aquifer_transpire) and MODFLOW decides how
  ! much of it the water table can actually supply.  The two differ whenever the water table drops
  ! below the EVT extraction depth; the gap is reported by the coupled budget diagnostic rather than
  ! fed back into SUMMA's energy balance, which would need the tight (XMI) coupling.
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
        i = (iGRU-1) * gru_struc(iGRU)%hruCount + jHRU
        mfAquiferTranspire(i) = transpire(i)
      end do
    end do
  end subroutine mf6x_put_aquifer_transpire

  ! **************************************************************************************************
  ! SUMMA's aquifer transpiration DEMAND, per HRU (m s-1, + = out of aquifer).
  !
  ! This is scalarAquiferTranspire as bigAquifer computes it: the aquifer's share of the canopy
  ! transpiration the energy balance has already closed, partitioned by root fraction and the aquifer
  ! transpiration limiting factor.  Sending it to MODFLOW's EVT package therefore moves water between
  ! sources without creating any, and costs no energy-balance consistency.
  ! **************************************************************************************************
  subroutine mf6x_get_aquifer_transpire(summa_struct, demand)
    type(summa1_type_dec), intent(in)  :: summa_struct
    real,                  intent(out) :: demand(:)
    integer(i4b) :: iGRU, jHRU, iDOM, i
    real(rkind)  :: fracDOM, frac, tlim
    ! Reproduces bigAquifer's partition (bigAquifer.f90:112-113) from post-step diagnostics:
    !
    !   aquiferTranspireFrac   = scalarAquiferRootFrac * scalarTranspireLimAqfr / scalarTranspireLim
    !   scalarAquiferTranspire = aquiferTranspireFrac * scalarCanopyTranspiration / iden_water
    !
    ! It is computed here rather than in computFlux because in coupled mode there is no aquifer
    ! state: ixAqWat is integerMissing, so computFlux never enters its aquifer block and bigAquifer
    ! is never called.  Nothing is lost by deriving it afterwards - the quantity is a pure sink that
    ! leaves SUMMA for MODFLOW and enters no SUMMA state equation, and the derivatives bigAquifer
    ! would also produce are only ever used in ixAqWat-indexed Jacobian entries that do not exist
    ! in this mode (computJacob.f90:858,861).
    !
    ! This is a DEMAND, not an extraction.  The energy balance behind scalarCanopyTranspiration has
    ! already closed, and scalarTranspireLimAqfr has already told stomatal resistance how much deep
    ! water is available, so sending this to MODFLOW's EVT package moves water between sources
    ! without creating any.  MODFLOW decides what it can actually supply.
    associate(progStruct => summa_struct%progStruct, &
              diagStruct => summa_struct%diagStruct, &
              fluxStruct => summa_struct%fluxStruct, &
              bvarStruct => summa_struct%bvarStruct)
      demand = 0.0
      do iGRU = 1, summa_struct%nGRU_local
        do jHRU = 1, gru_struc(iGRU)%hruCount
          i = (iGRU-1) * gru_struc(iGRU)%hruCount + jHRU
          demand(i) = 0._rkind
          do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
            tlim = diagStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookDIAG%scalarTranspireLim)%dat(1)
            if (tlim <= 0._rkind) cycle       ! no transpiration at all, so no aquifer share
            frac = diagStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookDIAG%scalarAquiferRootFrac)%dat(1) &
                 * diagStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookDIAG%scalarTranspireLimAqfr)%dat(1) / tlim
            if (frac <= 0._rkind) cycle       ! no roots below the soil column, or water table out of reach
            fracDOM = progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1) &
                    / bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)
            ! SIGN: scalarCanopyTranspiration (kg m-2 s-1) is NEGATIVE while water leaves the
            ! canopy - the same convention bigAquifer inherits, which is why computFlux assembles
            ! the aquifer balance as dS/dt = scalarAquiferTranspire + recharge - baseflow
            ! (computFlux.f90:555) with a negative transpire term.  MODFLOW's EVT RATE is a
            ! positive maximum extraction rate, so negate to get a demand positive out of the
            ! aquifer, matching the convention of every other flux this module hands the coupler.
            demand(i) = demand(i) &
                      - frac * fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarCanopyTranspiration)%dat(1) &
                        / iden_water * fracDOM
          end do
          ! after the negation a negative demand means condensation onto the canopy, which the
          ! aquifer plays no part in
          if (demand(i) < 0._rkind) demand(i) = 0._rkind
        end do
      end do
    end associate
  end subroutine mf6x_get_aquifer_transpire

end module summa_mf6_exchange
