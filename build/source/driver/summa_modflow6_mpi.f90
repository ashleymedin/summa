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

program summa_modflow6_mpi
  ! ****************************************************************************************
  ! *** Thin BMI coupler: SUMMA land model  <-->  MODFLOW 6 groundwater model, MPI variant ***
  ! ****************************************************************************************
  !
  ! Same coupling as summa_modflow6.f90 (see there for the full exchange description and required
  ! model decisions, and mf6_coupling.f90 for the MODFLOW side itself), with SUMMA's GRUs split
  ! across MPI ranks the way summa_driver_mpi.f90 splits them for the plain (non-coupled) MPI
  ! build.  MODFLOW 6 itself stays a serial singleton on rank 0: this libmf6 build has no
  ! PETSc/MPI support, and the coupled model here represents one shared aquifer grid under
  ! potentially many GRUs' HRUs, not one MODFLOW grid per GRU, so a singleton MODFLOW fed by
  ! parallel SUMMA is the natural decomposition (real MODFLOW-side parallelism would need a
  ! PETSc-enabled libmf6 and coordinating two independent domain decompositions on the same
  ! ranks - out of scope here; see meson.build's -Dparallel=true / src/Solution/ParallelSolution.f90
  ! upstream if that is ever needed for a MODFLOW-solve-bound regional model).
  !
  ! Every rank runs its own SUMMA instance over its local GRU subset.  Each coupled step, ranks
  ! gather their local soil drainage to rank 0 (MPI_Gatherv), rank 0 alone drives MODFLOW 6 through
  ! mf6_coupling and then scatters the resulting head/storage/baseflow feedback back out
  ! (MPI_Scatterv).  Everything MODFLOW-side lives in mf6_coupling, shared with the serial coupler
  ! and with the calibration driver, so this file is only the SUMMA side and the rank bookkeeping.
  !
  ! See utils/test/test_mflow/README.md for the serial coupler's worked examples and two test cases.
  !
  ! Usage:  mpirun -n <nranks> summa_modflow6_mpi.exe <fileManager.txt> <summa_modflow6.config>

  use nr_type
  use mpi
  use mpi_context,     only : set_mpi_context
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
   modLatFlow,                      & ! as modflowCpl, plus lateral flow in the soil above
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

  ! ---- MPI parallel context ----
  integer(i4b)        :: myrank = 0, nRanks = 1, mpi_ierr
  character(len=1024) :: mpi_errmsg

  ! ---- the MODFLOW 6 side, live on rank 0 only ----
  type(mf6_coupler_type) :: coupler
  logical                :: feedback = .true., have_sy = .false., have_bflow = .false.
  integer                :: err
  character(len=1024)    :: message

  ! ---- SUMMA side ----
  ! nHRU/drain_hru/head_hru/... below are GLOBAL, domain-wide, rank-0-only arrays (exactly as in the
  ! serial coupler; mf6_coupling works on them unmodified).  Each rank's own SUMMA instance only ever
  ! sees its local GRU subset, so its BMI grid-0 size is nHRU_local, and its own exchange arrays are
  ! the "_local" ones - gathered into / scattered out of the global ones each step via
  ! MPI_Gatherv/MPI_Scatterv (see gatherv_real/scatterv_real).
  type(summa_bmi)            :: summa
  integer                    :: istat, nHRU, nHRU_local, modelTimeStep
  character(len=1024)        :: file_manager, config_file
  real, allocatable          :: drain_hru(:)     ! per-HRU soil drainage        (m s-1)
  real, allocatable          :: head_hru(:)      ! per-HRU prescribed head      (m, matric head at soil base)
  real, allocatable          :: bflow_hru(:)     ! per-HRU aquifer baseflow     (m s-1, + = out of aquifer)  -> scalarAquiferBaseflow
  real, allocatable          :: stor_hru(:)      ! per-HRU relative aquifer storage (m of water)              -> scalarAquiferStorage
  double precision, allocatable :: hru_x(:), hru_y(:), hru_z(:)  ! HRU centroid lon/lat and surface elevation
  double precision, allocatable :: soil_thk(:)   ! per-HRU SUMMA soil-column thickness (m), read from SUMMA
  double precision, allocatable :: hru_area(:)   ! per-HRU plan area (m2), for the area check and coupled budget

  ! ---- this rank's local slice (every rank allocates these; sized nHRU_local) ----
  real, allocatable          :: drain_hru_local(:), head_hru_local(:), bflow_hru_local(:), stor_hru_local(:)
  double precision, allocatable :: hru_x_local(:), hru_y_local(:), hru_z_local(:), soil_thk_local(:)
  double precision, allocatable :: hru_area_local(:)

  ! ---- rank-0 bookkeeping for MPI_Gatherv/MPI_Scatterv (counts/displs are per-HRU-variable, in
  !      rank order; valid because summa_work_balance's balance_even gives every rank a contiguous
  !      block of GRUs, hence (GRUs are contiguous in file/HRU order) a contiguous block of HRUs ----
  integer, allocatable        :: hru_counts(:), hru_displs(:)  ! size nRanks, rank 0 only

  integer :: nlay, nrow, ncol

  call MPI_Init(mpi_ierr)
  call set_mpi_context(MPI_COMM_WORLD, myrank, nRanks, mpi_ierr, mpi_errmsg)
  if (mpi_ierr /= MPI_SUCCESS) then
    write(*,*) 'summa_modflow6_mpi: '//trim(mpi_errmsg)
    call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
  end if

  call initialize_coupler
  call run_coupler
  call finalize_coupler

  call MPI_Finalize(mpi_ierr)

contains

  ! ==================================================================================
  subroutine initialize_coupler
    integer :: i

    ! -- command line: file manager and this case's coupler config, both required --
    ! (mpirun/mpiexec replicate argv to every rank, so every rank parses these identically)
    if (command_argument_count() < 2) then
      if (myrank == 0) then
        write(*,*) 'usage: mpirun -n <nranks> summa_modflow6_mpi <fileManager.txt> <summa_modflow6.config>'
        write(*,*) '  the config is required: each case keeps its own beside its settings, so several'
        write(*,*) '  cases can share one MODFLOW model directory'
      end if
      call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
    end if
    call get_command_argument(1, file_manager)
    call get_command_argument(2, config_file)

    ! -- initialize SUMMA through its BMI, over this rank's slice of the run domain --
    istat = summa%initialize_mpi(trim(file_manager), MPI_COMM_WORLD, myrank, nRanks)
    if (istat /= BMI_OK) then
      write(*,*) 'summa_modflow6_mpi: SUMMA initialize failed on rank ', myrank
      call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
    end if

    ! -- the coupled-groundwater decision must be active (same model_decisions on every rank) --
    if (model_decisions(iLookDECISIONS%groundwatr)%iDecision /= modflowCpl .and. &
        model_decisions(iLookDECISIONS%groundwatr)%iDecision /= modLatFlow) then
      write(*,*) 'summa_modflow6_mpi: SUMMA model decision groundwatr must be "modflow" or "modLatflow" for the coupler'
      call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
    end if
    if (model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision /= prescribedHead) then
      write(*,*) 'summa_modflow6_mpi: SUMMA model decision bcLowrSoiH must be "presHead" for the coupler'
      call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
    end if

    ! -- this rank's local HRU count and geometry (BMI grid 0 = HRU points, local to this rank) --
    istat = summa%get_grid_size(0, nHRU_local)
    allocate(drain_hru_local(nHRU_local), head_hru_local(nHRU_local))
    allocate(bflow_hru_local(nHRU_local), stor_hru_local(nHRU_local))
    bflow_hru_local = 0.0; stor_hru_local = 0.0
    allocate(hru_x_local(nHRU_local), hru_y_local(nHRU_local), hru_z_local(nHRU_local))
    istat = summa%get_grid_x(0, hru_x_local)   ! HRU longitude  (deg or projected x, must match MODFLOW grid CRS)
    istat = summa%get_grid_y(0, hru_y_local)   ! HRU latitude   (deg or projected y)
    istat = summa%get_grid_z(0, hru_z_local)   ! HRU surface elevation (m)
    allocate(soil_thk_local(nHRU_local))
    istat = summa%get_soil_thickness(soil_thk_local)  ! SUMMA soil-column depth per HRU (m)
    allocate(hru_area_local(nHRU_local))
    istat = summa%get_hru_area(hru_area_local)        ! HRU plan area (m2)
    head_hru_local = 0.0

    ! -- gather local HRU counts to build the rank-order partition (see hru_counts/hru_displs above),
    !    and the global (domain-wide) HRU geometry that rank 0 needs to build the MODFLOW cell map --
    allocate(hru_counts(nRanks), hru_displs(nRanks))
    call MPI_Gather(nHRU_local, 1, MPI_INTEGER, hru_counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
    if (myrank == 0) then
      nHRU = sum(hru_counts)
      hru_displs(1) = 0
      do i = 2, nRanks
        hru_displs(i) = hru_displs(i-1) + hru_counts(i-1)
      end do
      allocate(drain_hru(nHRU), head_hru(nHRU), bflow_hru(nHRU), stor_hru(nHRU))
      allocate(hru_x(nHRU), hru_y(nHRU), hru_z(nHRU), soil_thk(nHRU), hru_area(nHRU))
      head_hru = 0.0; bflow_hru = 0.0; stor_hru = 0.0
    else
      ! placeholders: never dereferenced off rank 0, but must be allocated to legally pass as the
      ! (rank-0-significant) recvbuf/sendbuf argument of the Gatherv/Scatterv calls below
      allocate(drain_hru(1), head_hru(1), bflow_hru(1), stor_hru(1))
      allocate(hru_x(1), hru_y(1), hru_z(1), soil_thk(1), hru_area(1))
    end if
    call gatherv_dp(hru_x_local, hru_x)
    call gatherv_dp(hru_y_local, hru_y)
    call gatherv_dp(hru_z_local, hru_z)
    call gatherv_dp(soil_thk_local, soil_thk)
    call gatherv_dp(hru_area_local, hru_area)

    ! ================================================================================
    ! MODFLOW 6 runs on rank 0 alone; the other ranks' SUMMA instances are driven
    ! entirely from run_coupler's gather/scatter calls
    ! ================================================================================
    if (myrank == 0) then
      call coupler%init(trim(config_file), '.', nHRU, hru_x, hru_y, hru_z, soil_thk, &
                        numtim, dble(data_step), err, message, hru_area=hru_area)
      if (err /= 0) then
        write(*,'(a)') 'summa_modflow6_mpi: '//trim(message)
        call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
      end if
      feedback   = coupler%feedback
      have_sy    = coupler%have_sy
      have_bflow = coupler%have_bflow

      call coupler%grid_shape(nlay, nrow, ncol)
      write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a)') 'summa_modflow6_mpi: coupling ', nHRU, ' SUMMA HRUs across ', &
            nRanks, ' MPI ranks to a ', nlay, ' x ', nrow, ' x ', ncol, ' MODFLOW 6 DIS grid'
    end if

    ! -- every rank needs to know what the feedback actually carries (decided on rank 0 above,
    !    from MODFLOW's own packages) so it knows which values to set on its own SUMMA instance --
    call MPI_Bcast(feedback,   1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_Bcast(have_sy,    1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_Bcast(have_bflow, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpi_ierr)
  end subroutine initialize_coupler

  ! ==================================================================================
  subroutine run_coupler
    do modelTimeStep = 1, numtim

      ! 1. push last step's MODFLOW state into SUMMA (lagged one step): rank 0 scatters its
      !    global feedback arrays out to every rank's local slice, then every rank applies
      !    its own slice to its own SUMMA instance
      if (feedback .and. modelTimeStep > 1) then
        call scatterv_real(head_hru, head_hru_local)
        istat = summa%set_value('soil_water_sat-zone_top__head', head_hru_local)
        if (have_sy) then
          call scatterv_real(stor_hru, stor_hru_local)
          istat = summa%set_value('aquifer_water__storage_thickness', stor_hru_local)
        end if
        if (have_bflow) then
          call scatterv_real(bflow_hru, bflow_hru_local)
          istat = summa%set_value('land_surface_water__baseflow_volume_flux', bflow_hru_local)
        end if
      end if

      ! 2. advance SUMMA one data step on every rank (reads forcing, runs physics, writes output)
      istat = summa%update()
      if (istat /= BMI_OK) then
        write(*,*) 'summa_modflow6_mpi: SUMMA update failed at step ', modelTimeStep, ' on rank ', myrank
        call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
      end if

      ! 3-5. every rank gets its own drainage slice, rank 0 gathers them into the domain-wide
      !      array and drives MODFLOW 6 with it, reading the new water table back
      istat = summa%get_value('soil_water__drainage_volume_flux', drain_hru_local)
      call gatherv_real(drain_hru_local, drain_hru)
      if (myrank == 0) then
        call coupler%step(modelTimeStep, dble(data_step), &
                          drain_hru, head_hru, stor_hru, bflow_hru, err, message)
        if (err /= 0) then
          write(*,'(a)') 'summa_modflow6_mpi: '//trim(message)
          call MPI_Abort(MPI_COMM_WORLD, 1, mpi_ierr)
        end if
      end if
    end do
  end subroutine run_coupler

  ! ==================================================================================
  subroutine finalize_coupler
    if (myrank == 0) then
      call coupler%finalize(err, message)
      if (err /= 0) write(*,'(a)') 'summa_modflow6_mpi: '//trim(message)
    end if
    istat = summa%finalize()
    call sleep(2)   ! let HDF5 close cleanly, as in the stock SUMMA driver
    if (myrank == 0) write(*,'(a)') 'summa_modflow6_mpi: finished simulation successfully.'
  end subroutine finalize_coupler

  ! ==================================================================================
  ! MPI_Gatherv/MPI_Scatterv wrappers between this rank's local HRU slice and rank 0's
  ! domain-wide array, using the rank-order partition built in initialize_coupler.
  ! ==================================================================================
  subroutine gatherv_real(local_arr, global_arr)
    real, intent(in)    :: local_arr(:)     ! this rank's slice, size nHRU_local
    real, intent(inout) :: global_arr(:)    ! rank 0: size nHRU; other ranks: unused placeholder
    call MPI_Gatherv(local_arr, nHRU_local, MPI_REAL, &
                      global_arr, hru_counts, hru_displs, MPI_REAL, &
                      0, MPI_COMM_WORLD, mpi_ierr)
  end subroutine gatherv_real

  subroutine scatterv_real(global_arr, local_arr)
    real, intent(in)    :: global_arr(:)    ! rank 0: size nHRU; other ranks: unused placeholder
    real, intent(inout) :: local_arr(:)     ! this rank's slice, size nHRU_local
    call MPI_Scatterv(global_arr, hru_counts, hru_displs, MPI_REAL, &
                       local_arr, nHRU_local, MPI_REAL, &
                       0, MPI_COMM_WORLD, mpi_ierr)
  end subroutine scatterv_real

  subroutine gatherv_dp(local_arr, global_arr)
    double precision, intent(in)    :: local_arr(:)   ! this rank's slice, size nHRU_local
    double precision, intent(inout) :: global_arr(:)  ! rank 0: size nHRU; other ranks: unused placeholder
    call MPI_Gatherv(local_arr, nHRU_local, MPI_DOUBLE_PRECISION, &
                      global_arr, hru_counts, hru_displs, MPI_DOUBLE_PRECISION, &
                      0, MPI_COMM_WORLD, mpi_ierr)
  end subroutine gatherv_dp

end program summa_modflow6_mpi
