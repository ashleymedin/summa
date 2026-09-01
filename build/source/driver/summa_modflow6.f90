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
  ! ****************************************************************************************
  ! *** Thin BMI coupler: SUMMA land model  <-->  MODFLOW 6 groundwater model          ***
  ! ****************************************************************************************
  !
  ! Each SUMMA data step:
  !   1. (feedback) the MODFLOW 6 water-table head from the previous step is written into
  !      SUMMA as the prescribed-head lower boundary condition of the soil column
  !      (BMI input  "soil_water_sat-zone_top__head", parameter "lowerBoundHead").
  !   2. SUMMA advances one step.
  !   3. the drainage out the base of the SUMMA soil column
  !      (BMI output "soil_water__drainage_volume_flux", flux "scalarSoilDrainage")
  !      is regridded onto the MODFLOW 6 grid and written into the RCH package
  !      RECHARGE array.
  !   4. MODFLOW 6 advances one step (prepare/do/finalize_time_step).
  !   5. the new MODFLOW 6 head field is read back and aggregated per SUMMA HRU,
  !      ready to be applied at step 1 of the next iteration (explicit, one-step lag).
  !
  ! Coupling requires SUMMA to be built with -DUSE_MODFLOW6=ON (sets MODFLOW_ACTIVE) and
  ! the SUMMA model decision  groundwatr = modflow  with  bcLowrSoiH = presHead.
  !
  ! The MODFLOW 6 model (read from mfsim.nam in the working directory) must:
  !   * use length unit metres and TDIS TIME_UNITS SECONDS,
  !   * have exactly one MODFLOW time step per SUMMA data step (delt == data_step),
  !   * contain an RCH package with READASARRAYS,
  !   * be a single GWF model discretised with DIS.
  !
  ! Configuration (Fortran namelist), default file name "summa_modflow6.config":
  !
  !   &coupler
  !     mf6_model_name   = 'MYMODEL'   ! GWF model name, as in mfsim.nam (upper case)
  !     rch_package_name = 'RCHA'      ! RCH package name, as in the GWF name file (upper case)
  !     soil_thickness   = 2.0         ! m, thickness of the SUMMA soil column (used to place
  !                                    !    the soil base relative to the HRU surface elevation)
  !     map_file         = ''          ! optional HRU->cell weight file; if blank a nearest-cell
  !                                    !    map is built from the MODFLOW 6 DIS grid geometry
  !     feedback         = .true.      ! .false. => one-way (SUMMA drainage -> MODFLOW only)
  !   /
  !
  ! Usage:  summa_modflow6 <fileManager.txt> [summa_modflow6.config]

  use, intrinsic :: iso_c_binding
  use nr_type
  use summabmi,        only : summa_bmi
  use globalData,      only : numtim               ! number of SUMMA data steps
  use globalData,      only : data_step            ! length of a SUMMA data step (s)
  use globalData,      only : model_decisions      ! SUMMA model decision structure
  use var_lookup,      only : iLookDECISIONS       ! named indices into model_decisions
  use mDecisions_module, only : modflowCpl, prescribedHead

  implicit none

  ! ------------------------------------------------------------------------------------
  ! MODFLOW 6 shared library (libmf6): bind directly to the exported C entry points.
  ! ------------------------------------------------------------------------------------
  interface
    integer(c_int) function mf6_initialize() bind(C, name="initialize")
      import :: c_int
    end function mf6_initialize
    integer(c_int) function mf6_finalize() bind(C, name="finalize")
      import :: c_int
    end function mf6_finalize
    integer(c_int) function mf6_get_current_time(t) bind(C, name="get_current_time")
      import :: c_int, c_double
      real(c_double), intent(out) :: t
    end function mf6_get_current_time
    integer(c_int) function mf6_get_end_time(t) bind(C, name="get_end_time")
      import :: c_int, c_double
      real(c_double), intent(out) :: t
    end function mf6_get_end_time
    integer(c_int) function mf6_get_time_step(dt) bind(C, name="get_time_step")
      import :: c_int, c_double
      real(c_double), intent(out) :: dt
    end function mf6_get_time_step
    integer(c_int) function mf6_prepare_time_step(dt) bind(C, name="prepare_time_step")
      import :: c_int, c_double
      real(c_double), intent(in) :: dt
    end function mf6_prepare_time_step
    integer(c_int) function mf6_do_time_step() bind(C, name="do_time_step")
      import :: c_int
    end function mf6_do_time_step
    integer(c_int) function mf6_finalize_time_step() bind(C, name="finalize_time_step")
      import :: c_int
    end function mf6_finalize_time_step
    integer(c_int) function mf6_get_value_ptr_double(addr, cptr) bind(C, name="get_value_ptr_double")
      import :: c_int, c_char, c_ptr
      character(kind=c_char), intent(in)    :: addr(*)
      type(c_ptr),           intent(inout) :: cptr
    end function mf6_get_value_ptr_double
    integer(c_int) function mf6_get_value_ptr_int(addr, cptr) bind(C, name="get_value_ptr_int")
      import :: c_int, c_char, c_ptr
      character(kind=c_char), intent(in)    :: addr(*)
      type(c_ptr),           intent(inout) :: cptr
    end function mf6_get_value_ptr_int
    integer(c_int) function mf6_get_var_itemsize(addr, n) bind(C, name="get_var_itemsize")
      import :: c_int, c_char
      character(kind=c_char), intent(in)  :: addr(*)
      integer(c_int),         intent(out) :: n
    end function mf6_get_var_itemsize
    integer(c_int) function mf6_get_var_nbytes(addr, n) bind(C, name="get_var_nbytes")
      import :: c_int, c_char
      character(kind=c_char), intent(in)  :: addr(*)
      integer(c_int),         intent(out) :: n
    end function mf6_get_var_nbytes
  end interface

  integer, parameter :: BMI_OK = 0

  ! ---- configuration ----
  character(len=256) :: mf6_model_name   = ''
  character(len=256) :: rch_package_name = 'RCHA'
  real(rkind)        :: soil_thickness   = 2.0_rkind
  character(len=256) :: map_file         = ''
  logical            :: feedback         = .true.
  namelist /coupler/ mf6_model_name, rch_package_name, soil_thickness, map_file, feedback

  ! ---- SUMMA side ----
  type(summa_bmi)            :: summa
  integer                    :: istat, nHRU, modelTimeStep
  character(len=1024)        :: file_manager, config_file
  real, allocatable          :: drain_hru(:)     ! per-HRU soil drainage        (m s-1)
  real, allocatable          :: head_hru(:)      ! per-HRU prescribed head      (m, matric head at soil base)
  double precision, allocatable :: hru_x(:), hru_y(:), hru_z(:)  ! HRU centroid lon/lat and surface elevation

  ! ---- MODFLOW 6 side ----
  real(c_double), pointer    :: mf6_head(:) => null()      ! GWF dependent variable  <MODEL>/X       (nodes)
  real(c_double), pointer    :: mf6_rch(:)  => null()      ! RCH recharge array      <MODEL>/<RCH>/RECHARGE
  integer(c_int), pointer    :: mf6_mshape(:) => null()    ! DIS grid shape          <MODEL>/DIS/MSHAPE (nlay,nrow,ncol)
  real(c_double), pointer    :: cellx(:) => null(), celly(:) => null()
  real(c_double), pointer    :: xorigin => null(), yorigin => null(), angrot => null()
  integer                    :: nlay, nrow, ncol
  real(c_double)             :: dt_mf6, tend_mf6

  ! ---- HRU -> MODFLOW cell mapping (sparse) ----
  ! for HRU i, its contributing horizontal cells are map_cell(map_ptr(i):map_ptr(i+1)-1)
  ! with weights map_wgt(...) that sum to 1
  integer, allocatable :: map_ptr(:)   ! size nHRU+1
  integer, allocatable :: map_cell(:)  ! horizontal cell index  (irow-1)*ncol + icol
  real,    allocatable :: map_wgt(:)
  real(c_double), allocatable :: cell_area(:) ! horizontal cell plan area (m2), for the recharge volume split

  call initialize_coupler
  call run_coupler
  call finalize_coupler

contains

  ! ==================================================================================
  subroutine initialize_coupler
    integer :: fu, rc

    ! -- command line: file manager (required), config file (optional) --
    if (command_argument_count() < 1) then
      write(*,*) 'usage: summa_modflow6 <fileManager.txt> [summa_modflow6.config]'
      error stop 1
    end if
    call get_command_argument(1, file_manager)
    config_file = 'summa_modflow6.config'
    if (command_argument_count() >= 2) call get_command_argument(2, config_file)

    open(action='read', file=trim(config_file), iostat=rc, newunit=fu)
    if (rc /= 0) then
      write(*,*) 'summa_modflow6: cannot open configuration file '//trim(config_file)
      error stop 1
    end if
    read(nml=coupler, iostat=rc, unit=fu)
    close(fu)
    if (rc /= 0) then
      write(*,*) 'summa_modflow6: error reading &coupler namelist from '//trim(config_file)
      error stop 1
    end if
    if (len_trim(mf6_model_name) == 0) then
      write(*,*) 'summa_modflow6: mf6_model_name must be set in '//trim(config_file)
      error stop 1
    end if
    call to_upper(mf6_model_name)
    call to_upper(rch_package_name)

    ! -- initialize SUMMA through its BMI --
    istat = summa%initialize(trim(file_manager))
    if (istat /= BMI_OK) then; write(*,*) 'summa_modflow6: SUMMA initialize failed'; error stop 1; end if

    ! -- the coupled-groundwater decision must be active --
    if (model_decisions(iLookDECISIONS%groundwatr)%iDecision /= modflowCpl) then
      write(*,*) 'summa_modflow6: SUMMA model decision groundwatr must be "modflow" for the coupler'
      error stop 1
    end if
    if (model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision /= prescribedHead) then
      write(*,*) 'summa_modflow6: SUMMA model decision bcLowrSoiH must be "presHead" for the coupler'
      error stop 1
    end if

    ! -- SUMMA HRU count and geometry (BMI grid 0 = HRU points) --
    istat = summa%get_grid_size(0, nHRU)
    allocate(drain_hru(nHRU), head_hru(nHRU))
    allocate(hru_x(nHRU), hru_y(nHRU), hru_z(nHRU))
    istat = summa%get_grid_x(0, hru_x)   ! HRU longitude  (deg or projected x, must match MODFLOW grid CRS)
    istat = summa%get_grid_y(0, hru_y)   ! HRU latitude   (deg or projected y)
    istat = summa%get_grid_z(0, hru_z)   ! HRU surface elevation (m)
    head_hru = 0.0

    ! -- initialize MODFLOW 6 (reads mfsim.nam from the current working directory) --
    istat = mf6_initialize()
    if (istat /= BMI_OK) then; write(*,*) 'summa_modflow6: MODFLOW 6 initialize failed'; error stop 1; end if

    ! -- DIS grid dimensions, then pointers to the MODFLOW 6 arrays we exchange --
    call mf6_ptr_int_arr(trim(mf6_model_name)//'/DIS/MSHAPE', mf6_mshape, 3)
    nlay = mf6_mshape(1); nrow = mf6_mshape(2); ncol = mf6_mshape(3)
    call mf6_ptr_double(trim(mf6_model_name)//'/X', mf6_head, &
                        mf6_var_count(trim(mf6_model_name)//'/X'))
    call mf6_ptr_double(trim(mf6_model_name)//'/'//trim(rch_package_name)//'/RECHARGE', &
                        mf6_rch, mf6_var_count(trim(mf6_model_name)//'/'//trim(rch_package_name)//'/RECHARGE'))
    call mf6_ptr_double(trim(mf6_model_name)//'/DIS/CELLX', cellx, ncol)
    call mf6_ptr_double(trim(mf6_model_name)//'/DIS/CELLY', celly, nrow)
    call mf6_ptr_scalar(trim(mf6_model_name)//'/DIS/XORIGIN', xorigin)
    call mf6_ptr_scalar(trim(mf6_model_name)//'/DIS/YORIGIN', yorigin)
    call mf6_ptr_scalar(trim(mf6_model_name)//'/DIS/ANGROT',  angrot)

    ! -- coupling time-step consistency --
    ! NB: MODFLOW's delt is only set once the first time step is prepared, so it is
    ! not meaningful yet right after initialize(); the delt == data_step check is
    ! therefore deferred to the first iteration of run_coupler (see below).
    istat = mf6_get_end_time(tend_mf6)
    if (tend_mf6 < numtim*real(data_step, c_double)*(1.0_c_double - 1.0e-6_c_double)) then
      write(*,'(a)') 'summa_modflow6: WARNING - MODFLOW simulation is shorter than the SUMMA simulation'
    end if

    ! -- build the HRU -> MODFLOW horizontal-cell weight map --
    call build_map

    write(*,'(a,i0,a,i0,a,i0,a,i0,a)') 'summa_modflow6: coupling ', nHRU, ' SUMMA HRUs to a ', &
          nlay, ' x ', nrow, ' x ', ncol, ' MODFLOW 6 DIS grid'
  end subroutine initialize_coupler

  ! ==================================================================================
  subroutine run_coupler
    do modelTimeStep = 1, numtim

      ! 1. push last step's MODFLOW head into SUMMA as the soil-column lower BC (lagged)
      if (feedback .and. modelTimeStep > 1) then
        istat = summa%set_value('soil_water_sat-zone_top__head', head_hru)
      end if

      ! 2. advance SUMMA one data step (reads forcing, runs physics, writes output)
      istat = summa%update()
      if (istat /= BMI_OK) then; write(*,*) 'summa_modflow6: SUMMA update failed at step ', modelTimeStep; error stop 1; end if

      ! 3. SUMMA soil drainage -> MODFLOW 6 RCH recharge array
      istat = summa%get_value('soil_water__drainage_volume_flux', drain_hru)
      call scatter_drainage_to_rch

      ! 4. advance MODFLOW 6 one coupling step
      istat = mf6_prepare_time_step(real(data_step, c_double))
      if (modelTimeStep == 1) then
        ! delt is now set: verify one MODFLOW time step == one SUMMA data step
        istat = mf6_get_time_step(dt_mf6)
        if (abs(dt_mf6 - real(data_step, c_double)) > 1.0e-6_c_double*real(data_step, c_double)) then
          write(*,'(a,g0,a,g0,a)') 'summa_modflow6: MODFLOW time step (', dt_mf6, &
                ') must equal the SUMMA data step (', data_step, ')'
          error stop 1
        end if
      end if
      istat = mf6_do_time_step()
      istat = mf6_finalize_time_step()

      ! 5. read the new MODFLOW head, aggregate per HRU for the next iteration
      if (feedback) call gather_head_to_hru
    end do
  end subroutine run_coupler

  ! ==================================================================================
  subroutine finalize_coupler
    istat = mf6_finalize()
    istat = summa%finalize()
    call sleep(2)   ! let HDF5 close cleanly, as in the stock SUMMA driver
    write(*,'(a)') 'summa_modflow6: finished simulation successfully.'
  end subroutine finalize_coupler

  ! ==================================================================================
  ! SUMMA drainage (m s-1, per HRU)  ->  MODFLOW RECHARGE array (m s-1, per top cell).
  ! Each HRU's flux is spread over its mapped cells by weight; a cell that receives
  ! from several HRUs gets the area-weighted mean flux (so recharge volume is conserved
  ! when the HRU areas equal the covered cell areas).
  ! ==================================================================================
  subroutine scatter_drainage_to_rch
    real(c_double), allocatable :: num(:), den(:)
    integer :: i, k, c

    allocate(num(nrow*ncol), den(nrow*ncol))
    num = 0.0_c_double
    den = 0.0_c_double
    do i = 1, nHRU
      do k = map_ptr(i), map_ptr(i+1) - 1
        c = map_cell(k)
        if (c < 1 .or. c > nrow*ncol) cycle
        num(c) = num(c) + real(map_wgt(k), c_double) * cell_area(c) * real(drain_hru(i), c_double)
        den(c) = den(c) + real(map_wgt(k), c_double) * cell_area(c)
      end do
    end do

    ! RECHARGE is indexed by horizontal cell (row-major) for READASARRAYS RCH
    do c = 1, min(size(mf6_rch), nrow*ncol)
      if (den(c) > 0.0_c_double) then
        mf6_rch(c) = num(c) / den(c)
      else
        mf6_rch(c) = 0.0_c_double
      end if
    end do
    deallocate(num, den)
  end subroutine scatter_drainage_to_rch

  ! ==================================================================================
  ! MODFLOW head (m, per node)  ->  SUMMA prescribed lower-BC matric head (m, per HRU).
  ! The head is sampled at the top-most active node of each mapped column and converted
  ! to a matric (pressure) head at the base of the SUMMA soil column:
  !     lowerBoundHead = h_mf6 - (z_surface_HRU - soil_thickness)
  ! positive => saturated / positive pressure head at the soil base.
  ! ==================================================================================
  subroutine gather_head_to_hru
    integer :: i, k, c, node
    real(c_double) :: hsum, wsum, zbase

    do i = 1, nHRU
      hsum = 0.0_c_double
      wsum = 0.0_c_double
      do k = map_ptr(i), map_ptr(i+1) - 1
        c = map_cell(k)                     ! horizontal (row-major) cell index
        node = top_active_node(c)
        if (node > 0) then
          hsum = hsum + real(map_wgt(k), c_double) * mf6_head(node)
          wsum = wsum + real(map_wgt(k), c_double)
        end if
      end do
      if (wsum > 0.0_c_double) then
        zbase = real(hru_z(i), c_double) - real(soil_thickness, c_double)
        head_hru(i) = real(hsum / wsum - zbase)
      end if
    end do
  end subroutine gather_head_to_hru

  ! For a DIS grid the reduced node number equals the full node number when no cells
  ! are removed; walk down the column and return the first node with a usable head.
  integer function top_active_node(choriz) result(node)
    integer, intent(in) :: choriz          ! (irow-1)*ncol + icol
    integer :: ilay, full
    node = -1
    do ilay = 1, nlay
      full = (ilay-1)*nrow*ncol + choriz
      if (full >= 1 .and. full <= size(mf6_head)) then
        if (mf6_head(full) > -1.0e29_c_double) then   ! HNOFLO / HDRY guard
          node = full
          return
        end if
      end if
    end do
  end function top_active_node

  ! ==================================================================================
  ! Build the HRU -> horizontal-cell weight map, either from an external weight file
  ! or, by default, as a nearest-cell assignment using the MODFLOW 6 DIS geometry.
  ! ==================================================================================
  subroutine build_map
    integer :: i, j, c
    real(c_double) :: dx, dy

    allocate(cell_area(nrow*ncol))
    do i = 1, nrow
      do j = 1, ncol
        c = (i-1)*ncol + j
        ! CELLX/CELLY are grid-local centres; DELR = spacing in x, DELC in y.
        ! Approximate the plan area from local centre spacing (uniform grids exact).
        dx = cell_spacing(cellx, j, ncol)
        dy = cell_spacing(celly, i, nrow)
        cell_area(c) = dx * dy
      end do
    end do

    if (len_trim(map_file) > 0) then
      call read_map_file
    else
      call build_nearest_cell_map
    end if
  end subroutine build_map

  real(c_double) function cell_spacing(centres, idx, ncell) result(d)
    real(c_double), intent(in) :: centres(:)
    integer,        intent(in) :: idx, ncell
    if (ncell == 1) then
      d = 1.0_c_double
    else if (idx == 1) then
      d = abs(centres(2) - centres(1))
    else if (idx == ncell) then
      d = abs(centres(ncell) - centres(ncell-1))
    else
      d = 0.5_c_double * abs(centres(idx+1) - centres(idx-1))
    end if
  end function cell_spacing

  subroutine build_nearest_cell_map
    integer :: i, j, k, cbest
    real(c_double) :: gx, gy, ca, sa, lx, ly, dbest, dist

    allocate(map_ptr(nHRU+1), map_cell(nHRU), map_wgt(nHRU))
    ca = cos(-real(angrot, c_double) * acos(-1.0_c_double) / 180.0_c_double)
    sa = sin(-real(angrot, c_double) * acos(-1.0_c_double) / 180.0_c_double)
    do i = 1, nHRU
      ! world -> grid-local frame (translate by origin, rotate by -angrot)
      gx = real(hru_x(i), c_double) - real(xorigin, c_double)
      gy = real(hru_y(i), c_double) - real(yorigin, c_double)
      lx = ca*gx - sa*gy
      ly = sa*gx + ca*gy
      dbest = huge(dbest)
      cbest = 1
      do j = 1, nrow
        do k = 1, ncol
          dist = (lx - cellx(k))**2 + (ly - celly(j))**2
          if (dist < dbest) then
            dbest = dist
            cbest = (j-1)*ncol + k
          end if
        end do
      end do
      map_ptr(i)  = i
      map_cell(i) = cbest
      map_wgt(i)  = 1.0
    end do
    map_ptr(nHRU+1) = nHRU + 1
  end subroutine build_nearest_cell_map

  ! Weight-file format (whitespace separated, "#" comment lines ignored):
  !   one record per SUMMA HRU, in BMI flatten order:
  !     iHRU   nPairs   cell_1 w_1   cell_2 w_2   ...
  !   cell index is the row-major horizontal index (irow-1)*ncol + icol; weights sum to 1.
  subroutine read_map_file
    integer :: fu, rc, i, ih, np, k, total, pos
    integer, allocatable :: cnt(:)
    character(len=4096) :: line

    open(action='read', file=trim(map_file), iostat=rc, newunit=fu)
    if (rc /= 0) then
      write(*,*) 'summa_modflow6: cannot open map_file '//trim(map_file)
      error stop 1
    end if

    allocate(cnt(nHRU))
    cnt = 0
    total = 0
    ! first pass: count pairs
    do
      read(fu,'(a)',iostat=rc) line
      if (rc /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      read(line,*,iostat=rc) ih, np
      if (rc /= 0 .or. ih < 1 .or. ih > nHRU) cycle
      cnt(ih) = np
      total = total + np
    end do
    rewind(fu)

    allocate(map_ptr(nHRU+1), map_cell(total), map_wgt(total))
    map_ptr(1) = 1
    do i = 1, nHRU
      map_ptr(i+1) = map_ptr(i) + cnt(i)
    end do

    ! second pass: fill
    do
      read(fu,'(a)',iostat=rc) line
      if (rc /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      read(line,*,iostat=rc) ih, np
      if (rc /= 0 .or. ih < 1 .or. ih > nHRU) cycle
      pos = map_ptr(ih)
      read(line,*,iostat=rc) ih, np, (map_cell(pos+k-1), map_wgt(pos+k-1), k = 1, np)
      if (rc /= 0) then
        write(*,*) 'summa_modflow6: bad record in map_file for HRU ', ih
        error stop 1
      end if
    end do
    close(fu)
  end subroutine read_map_file

  ! ==================================================================================
  ! small helpers
  ! ==================================================================================
  subroutine mf6_ptr_double(addr, fptr, n)
    character(len=*), intent(in) :: addr
    real(c_double), pointer, intent(out) :: fptr(:)
    integer, intent(in) :: n                     ! expected element count
    type(c_ptr) :: cptr
    integer :: s
    s = mf6_get_value_ptr_double(cstr(addr), cptr)
    if (s /= BMI_OK .or. .not. c_associated(cptr)) then
      write(*,*) 'summa_modflow6: MODFLOW variable not found: '//trim(addr)
      error stop 1
    end if
    call c_f_pointer(cptr, fptr, [n])
  end subroutine mf6_ptr_double

  subroutine mf6_ptr_scalar(addr, fptr)
    character(len=*), intent(in) :: addr
    real(c_double), pointer, intent(out) :: fptr
    type(c_ptr) :: cptr
    integer :: s
    s = mf6_get_value_ptr_double(cstr(addr), cptr)
    if (s /= BMI_OK .or. .not. c_associated(cptr)) then
      write(*,*) 'summa_modflow6: MODFLOW scalar not found: '//trim(addr)
      error stop 1
    end if
    call c_f_pointer(cptr, fptr)
  end subroutine mf6_ptr_scalar

  subroutine mf6_ptr_int_arr(addr, fptr, n)
    character(len=*), intent(in) :: addr
    integer(c_int), pointer, intent(out) :: fptr(:)
    integer, intent(in) :: n
    type(c_ptr) :: cptr
    integer :: s
    s = mf6_get_value_ptr_int(cstr(addr), cptr)
    if (s /= BMI_OK .or. .not. c_associated(cptr)) then
      write(*,*) 'summa_modflow6: MODFLOW variable not found: '//trim(addr)
      error stop 1
    end if
    call c_f_pointer(cptr, fptr, [n])
  end subroutine mf6_ptr_int_arr

  ! number of elements of a MODFLOW variable = nbytes / itemsize
  integer function mf6_var_count(addr) result(n)
    character(len=*), intent(in) :: addr
    integer(c_int) :: nbytes, isize
    integer :: s
    s = mf6_get_var_nbytes(cstr(addr), nbytes)
    if (s /= BMI_OK) then
      write(*,*) 'summa_modflow6: MODFLOW variable not found: '//trim(addr); error stop 1
    end if
    s = mf6_get_var_itemsize(cstr(addr), isize)
    if (s /= BMI_OK .or. isize <= 0) then
      write(*,*) 'summa_modflow6: bad itemsize for MODFLOW variable: '//trim(addr); error stop 1
    end if
    n = nbytes / isize
  end function mf6_var_count

  ! Convert a Fortran string to a null-terminated C character array.
  function cstr(f) result(c)
    character(len=*), intent(in) :: f
    character(kind=c_char), allocatable :: c(:)
    integer :: i
    allocate(c(len_trim(f)+1))
    do i = 1, len_trim(f)
      c(i) = f(i:i)
    end do
    c(len_trim(f)+1) = c_null_char
  end function cstr

  subroutine to_upper(s)
    character(len=*), intent(inout) :: s
    integer :: i, k
    do i = 1, len_trim(s)
      k = iachar(s(i:i))
      if (k >= iachar('a') .and. k <= iachar('z')) s(i:i) = achar(k - 32)
    end do
  end subroutine to_upper

end program summa_modflow6
