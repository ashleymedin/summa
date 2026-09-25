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

program summa_modflow6
  ! Thin BMI coupler: SUMMA land model <-> MODFLOW 6 groundwater model.
  !
  ! This program is the SUMMA side and the time loop only: it drives SUMMA through its BMI and hands
  ! each step's soil drainage to mf6_coupling, which owns everything MODFLOW.  See mf6_coupling.f90
  ! for the exchange, the &coupler namelist and the HRU->cell map, and
  ! utils/test/test_mflow/README.md for the test cases.
  !
  ! Usage:  summa_modflow6.exe <fileManager.txt> <summa_modflow6.config>
  !         or through utils/test/test_mflow/coupler_commands.sh, which resolves paths and cds into
  !         the MODFLOW case directory.
  !
  ! One exchange per SUMMA data step, explicit with a one-step lag.  Steps 1 and 2 are here, steps 3
  ! to 5 are mf6_coupling's mf6_step:
  !   1. (feedback) the previous step's MODFLOW state is written into SUMMA
  !   2. SUMMA advances one step
  !   3. soil drainage is regridded onto the MODFLOW grid and written into the RCH RECHARGE array,
  !      and the aquifer transpiration demand into the EVT RATE array
  !   4. MODFLOW advances one step (prepare/do/finalize_time_step), leading steady-state stress
  !      periods having been solved out first
  !   5. the new head field and the boundary-package flows are read back and aggregated per HRU
  !
  ! With feedback = .true. these come back each step, so SUMMA's water balance and routed streamflow
  ! include the aquifer:
  !     lowerBoundHead          prescribed head at the base of the soil column (enters the solver)
  !     scalarAquiferStorage    Sy * (MODFLOW water table - soil-column base), diagnostic
  !     scalarAquiferBaseflow   role=baseflow package outflow over the HRU footprint
  !     mfSurfaceDischarge      role=surface_discharge outflow, added to SUMMA's surface runoff
  !     scalarAquiferTranspire  role=gw_et extraction, against the demand SUMMA sent
  !     scalarTranspireLimAqfr  aquifer transpiration limiting factor, evaluated per MODFLOW cell
  ! (scalarAquiferRecharge is not exchanged - SUMMA sets it to its own soil drainage.)
  !
  ! Required SUMMA model decisions.  Build with -DUSE_MODFLOW6=ON (sets MODFLOW_ACTIVE); a plain
  ! summa run with either groundwater option is rejected at start-up.
  !
  !   groundwatr = modflow      bcLowrSoiH = presHead
  !   groundwatr = modLatflow   bcLowrSoiH = presHead, hc_profile = exp_prof,
  !                             infRateMax = topmodel_GA (or noInfExc)
  !
  ! modLatflow additionally runs TOPMODEL-style lateral flow through the soil column above the
  ! MODFLOW water table, for hillslopes where water moves downslope through the soil as well as
  ! recharging the aquifer.  It requires exp_prof because the lateral transmissivity is the vertical
  ! integral of the conductivity over the soil column alone, MODFLOW carrying everything below it,
  ! and exp_prof integrates to a finite base rather than assuming a shallow aquifer of its own.  The
  ! lateral flow is reported as basin__ColumnOutflow and added to total runoff alongside the MODFLOW
  ! baseflow.
  !
  ! The MODFLOW 6 model is read from mfsim.nam in the working directory.  What the coupler
  ! requires of it is checked at start-up by mf6_coupling (units, DIS, RCH READASARRAYS)
  ! and on the first step (one MODFLOW time step per SUMMA data step); see utils/test/test_mflow/README.md.

  use nr_type
  use mf6_coupling,    only : mf6_coupler_type
  use summabmi,        only : summa_bmi
  use globalData,      only : numtim               ! number of SUMMA data steps
  use globalData,      only : data_step            ! length of a SUMMA data step (s)
  use globalData,      only : model_decisions      ! SUMMA model decision structure
  use var_lookup,      only : iLookDECISIONS       ! named indices into model_decisions

  ! look-up values for the choice of groundwater parameterization
  USE mDecisions_module,only:       &
   qbaseTopmodel,                   & ! TOPMODEL-ish baseflow parameterization
   modflowCpl,                      & ! MODFLOW coupled groundwater parameterization
   modLatflow,                      & ! as modflowCpl, plus lateral flow in the soil above
   bigBucket,                       & ! a big bucket (lumped aquifer model)
   noExplicit                         ! no explicit groundwater parameterization

  ! look-up values for the choice of boundary conditions for hydrology
  USE mDecisions_module,only:       &
   prescribedHead,                  & ! prescribed head
   funcBottomHead,                  & ! function of matric head in the lower-most layer
   freeDrainage,                    & ! free drainage
   liquidFlux,                      & ! liquid water flux
   zeroFlux                           ! zero flux

  implicit none

  integer, parameter :: BMI_OK = 0

  ! ---- SUMMA side ----
  type(summa_bmi)        :: summa
  type(mf6_coupler_type) :: coupler
  integer                :: istat, nHRU, modelTimeStep
  integer                :: err
  character(len=1024)    :: message
  character(len=1024)    :: file_manager, config_file
  real, allocatable      :: drain_hru(:)     ! per-HRU soil drainage        (m s-1)
  real, allocatable      :: head_hru(:)      ! per-HRU prescribed head      (m, matric head at soil base)
  real, allocatable      :: bflow_hru(:)     ! per-HRU aquifer baseflow     (m s-1, + = out of aquifer)  -> scalarAquiferBaseflow
  real, allocatable      :: surfdis_hru(:)   ! per-HRU groundwater discharge at land surface (m s-1) -> surface runoff
  real, allocatable      :: gwet_dem_hru(:)  ! per-HRU aquifer transpiration DEMAND from SUMMA (m s-1) -> MODFLOW EVT
  real, allocatable      :: gwet_hru(:)      ! per-HRU groundwater ET actually taken by MODFLOW (m s-1)
  real, allocatable      :: gwet_lim_hru(:)  ! per-HRU aquifer transpiration limiting factor (-), cell-wise mean
  real, allocatable      :: stor_hru(:)      ! per-HRU relative aquifer storage (m of water)              -> scalarAquiferStorage
  double precision, allocatable :: hru_x(:), hru_y(:), hru_z(:)  ! HRU centroid lon/lat and surface elevation
  double precision, allocatable :: soil_thk(:)   ! per-HRU SUMMA soil-column thickness (m), read from SUMMA
  double precision, allocatable :: hru_area(:)   ! per-HRU plan area (m2), for the area check and coupled budget
  double precision, allocatable :: root_reach(:) ! per-HRU root reach below the soil column (m)
  integer :: nlay, nrow, ncol

  call initialize_coupler
  call run_coupler
  call finalize_coupler

contains

  ! ==================================================================================
  subroutine initialize_coupler

    ! -- command line: file manager and this case's coupler config, both required --
    if (command_argument_count() < 2) then
      write(*,*) 'usage: summa_modflow6 <fileManager.txt> <summa_modflow6.config>'
      write(*,*) '  the config is required: each case keeps its own beside its settings, so several'
      write(*,*) '  cases can share one MODFLOW model directory'
      error stop 1
    end if
    call get_command_argument(1, file_manager)
    call get_command_argument(2, config_file)

    ! -- initialize SUMMA through its BMI --
    istat = summa%initialize(trim(file_manager))
    if (istat /= BMI_OK) then; write(*,*) 'summa_modflow6: SUMMA initialize failed'; error stop 1; end if

    ! -- the coupled-groundwater decision must be active --
    if (model_decisions(iLookDECISIONS%groundwatr)%iDecision /= modflowCpl .and. &
        model_decisions(iLookDECISIONS%groundwatr)%iDecision /= modLatflow) then
      write(*,*) 'summa_modflow6: SUMMA model decision groundwatr must be "modflow" or "modLatflow" for the coupler'
      error stop 1
    end if
    if (model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision /= prescribedHead) then
      write(*,*) 'summa_modflow6: SUMMA model decision bcLowrSoiH must be "presHead" for the coupler'
      error stop 1
    end if

    ! -- SUMMA HRU count and geometry (BMI grid 0 = HRU points) --
    istat = summa%get_grid_size(0, nHRU)
    ! the feedback buffers are allocated whether or not feedback is on: they are handed to
    ! the coupler either way, and it simply leaves them alone when there is no feedback
    allocate(drain_hru(nHRU), head_hru(nHRU), bflow_hru(nHRU), stor_hru(nHRU), surfdis_hru(nHRU), &
             gwet_dem_hru(nHRU), gwet_hru(nHRU), gwet_lim_hru(nHRU))
    bflow_hru = 0.0; stor_hru = 0.0; surfdis_hru = 0.0; gwet_dem_hru = 0.0; gwet_hru = 0.0; gwet_lim_hru = 0.0
    allocate(hru_x(nHRU), hru_y(nHRU), hru_z(nHRU), soil_thk(nHRU), hru_area(nHRU), root_reach(nHRU))
    istat = summa%get_grid_x(0, hru_x)   ! HRU longitude  (deg or projected x, must match MODFLOW grid CRS)
    istat = summa%get_grid_y(0, hru_y)   ! HRU latitude   (deg or projected y)
    istat = summa%get_grid_z(0, hru_z)   ! HRU surface elevation (m)
    istat = summa%get_soil_thickness(soil_thk)  ! SUMMA soil-column depth per HRU (m)
    istat = summa%get_hru_area(hru_area)        ! HRU plan area (m2)
    istat = summa%get_root_reach(root_reach)   ! how far roots reach below the soil column (m)
    head_hru = 0.0

    ! -- start MODFLOW 6 and build the HRU -> cell map (run_dir '.': MODFLOW reads mfsim.nam
    !    from the working directory, as this program has always done) --
    call coupler%init(trim(config_file), '.', nHRU, hru_x, hru_y, hru_z, soil_thk, &
                      numtim, dble(data_step), err, message, hru_area=hru_area, root_reach=root_reach)
    if (err /= 0) then; write(*,'(a)') 'summa_modflow6: '//trim(message); error stop 1; end if

    call coupler%grid_shape(nlay, nrow, ncol)
    write(*,'(a,i0,a,i0,a,i0,a,i0,a)') 'summa_modflow6: coupling ', nHRU, ' SUMMA HRUs to a ', &
          nlay, ' x ', nrow, ' x ', ncol, ' MODFLOW 6 DIS grid'
  end subroutine initialize_coupler

  ! ==================================================================================
  subroutine run_coupler
    do modelTimeStep = 1, numtim

      ! 1. push last step's MODFLOW state into SUMMA (lagged one step)
      if (coupler%feedback .and. modelTimeStep > 1) then
        istat = summa%set_value('soil_water_sat-zone_top__head', head_hru)
        if (coupler%have_sy)    istat = summa%set_value('aquifer_water__storage_thickness', stor_hru)
        if (coupler%have_bflow) istat = summa%set_value('land_surface_water__baseflow_volume_flux', bflow_hru)
        if (coupler%have_surfdis) istat = summa%set_value('land_surface_water__domain_outflow_volume_flux', surfdis_hru)
        if (coupler%have_gwet)    istat = summa%set_value('land_vegetation_water__aquifer_transpiration_volume_flux', gwet_hru)
        ! the limiting factor the coupler evaluated per cell; soilResist uses it as given
        if (coupler%have_evt)     istat = summa%set_value('land_vegetation_water__aquifer_transpiration_limit', gwet_lim_hru)
      end if

      ! 2. advance SUMMA one data step (reads forcing, runs physics, writes output)
      istat = summa%update()
      if (istat /= BMI_OK) then; write(*,*) 'summa_modflow6: SUMMA update failed at step ', modelTimeStep; error stop 1; end if

      ! 3-5. SUMMA drainage -> MODFLOW recharge, advance MODFLOW, read the new water table back
      istat = summa%get_value('soil_water__drainage_volume_flux', drain_hru)
      ! SUMMA's aquifer transpiration demand goes down with the drainage; what MODFLOW could
      ! actually supply comes back and is applied (lagged) at the next step, like every other feedback
      if (coupler%have_evt) istat = summa%get_value('land_vegetation_water__aquifer_transpiration_volume_flux', gwet_dem_hru)
      call coupler%step(modelTimeStep, dble(data_step), &
                        drain_hru, head_hru, stor_hru, bflow_hru, err, message, &
                        surfdis_hru=surfdis_hru, gwet_demand_hru=gwet_dem_hru, gwet_hru=gwet_hru, &
                        gwet_lim_hru=gwet_lim_hru)
      if (err /= 0) then; write(*,'(a)') 'summa_modflow6: '//trim(message); error stop 1; end if
    end do
  end subroutine run_coupler

  ! ==================================================================================
  subroutine finalize_coupler
    call coupler%finalize(err, message)
    if (err /= 0) write(*,'(a)') 'summa_modflow6: '//trim(message)
    istat = summa%finalize()
    call sleep(2)   ! let HDF5 close cleanly, as in the stock SUMMA driver
    write(*,'(a)') 'summa_modflow6: finished simulation successfully.'
  end subroutine finalize_coupler

end program summa_modflow6
