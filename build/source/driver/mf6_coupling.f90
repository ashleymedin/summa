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

module mf6_coupling
  ! ****************************************************************************************
  ! *** The MODFLOW 6 side of the SUMMA <-> MODFLOW 6 coupling                           ***
  ! ****************************************************************************************
  !
  ! Everything here talks to libmf6 and to a grid; nothing here knows about SUMMA data
  ! structures.  The interface is per-HRU arrays in and per-HRU arrays out, so every driver
  ! that can produce a per-HRU soil drainage and accept a per-HRU water table can use it:
  !
  !   summa_modflow6.f90        serial coupler, SUMMA driven through its BMI
  !   summa_modflow6_mpi.f90    as above with SUMMA's GRUs split across ranks; MODFLOW 6
  !                             stays a serial singleton on rank 0, which is the only rank
  !                             that calls into this module
  !   summa_simulation.f90      one coupled run per calibration parameter sample
  !
  ! --- the exchange, once per SUMMA data step ---------------------------------------------
  ! The driver owns steps 1 and 2, this module owns steps 3 to 5 (see mf6_step):
  !   1. (feedback) the MODFLOW 6 water-table head from the previous step is written into
  !      SUMMA as the prescribed-head lower boundary condition of the soil column
  !      (parameter "lowerBoundHead").
  !   2. SUMMA advances one step.
  !   3. the drainage out the base of the SUMMA soil column ("scalarSoilDrainage") is
  !      regridded onto the MODFLOW 6 grid and written into the RCH package RECHARGE array.
  !   4. MODFLOW 6 advances one step (prepare/do/finalize_time_step).
  !   5. the new MODFLOW 6 head field is read back and aggregated per SUMMA HRU, ready to be
  !      applied at step 1 of the next iteration (explicit, one-step lag).
  !
  ! With feedback = .true. two further groundwater quantities come back each step, so SUMMA's
  ! water balance and routed streamflow include the aquifer:
  !     scalarAquiferStorage  = Sy * (MODFLOW water table - soil-column base)
  !     scalarAquiferBaseflow = MODFLOW <bflow_package_name> outflow over the HRU footprint
  ! (scalarAquiferRecharge is not exchanged - SUMMA sets it to its own soil drainage.)
  !
  ! --- configuration -----------------------------------------------------------------------
  !
  !   &coupler
  !     mf6_model_name     = 'MYMODEL' ! GWF model name, as in mfsim.nam (upper case)
  !     rch_package_name   = 'RCHA'    ! RCH package name, as in the GWF name file (upper case)
  !     bflow_package_name = 'CHD'     ! head-dependent boundary package (CHD/DRN/RIV/GHB) whose
  !                                    !    simulated flow feeds back per HRU as scalarAquiferBaseflow
  !                                    !    ('' => skip the baseflow feedback)
  !     map_file           = ''        ! optional HRU->cell weight file; if blank a nearest-cell
  !                                    !    map is built from the MODFLOW 6 DIS grid geometry
  !     mf6_epsg           = 0         ! EPSG code of the MODFLOW grid's projected CRS, used only to
  !                                    !    reproject SUMMA's lon/lat HRU centres before the built-in
  !                                    !    nearest-cell map (nHRU>1; see nearest_hru).  Only WGS84 UTM
  !                                    !    is supported (32601-32660 N, 32701-32760 S); 0 = no
  !                                    !    reprojection (fine when nHRU==1, or with an explicit map_file)
  !     feedback           = .true.    ! .false. => one-way (SUMMA drainage -> MODFLOW only)
  !   /
  !
  ! map_file format: one "iHRU  cell  weight" triple per line (whitespace separated; blank
  ! lines and '#' comments ignored).  cell is the row-major horizontal MODFLOW index
  ! (irow-1)*ncol + icol; weights are normalised per HRU, so put e.g. 1.0 on every line to
  ! spread an HRU over its cells.  An HRU may span any number of lines.
  !
  ! --- the run directory ---------------------------------------------------------------------
  ! libmf6 reads mfsim.nam from the process working directory and writes its listing, head and
  ! budget files there, so two MODFLOW instances in one directory overwrite each other's output.
  ! A driver that runs several coupled models at once (the calibration driver: one model instance
  ! per MPI rank) therefore passes each one its own run_dir, and this module enters and leaves it
  ! around every libmf6 call.  run_dir = '.' (the default, and what the standalone couplers pass)
  ! means "the process working directory", and no directory change happens at all.
  !
  ! --- lifetime ------------------------------------------------------------------------------
  ! mf6_init .. repeated mf6_step .. mf6_finalize is one MODFLOW simulation, and a coupler object
  ! may be initialized again afterwards for the next one; the calibration driver does exactly that,
  ! one MODFLOW simulation per parameter sample, so that every sample starts from the same aquifer
  ! initial condition.  libmf6 is process-global, so only one coupler can be live at a time.

  use, intrinsic :: iso_c_binding
  implicit none
  private

  public :: mf6_coupler_type
  public :: mf6_prepare_run_dir

  integer, parameter :: BMI_OK = 0

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

  ! ------------------------------------------------------------------------------------
  ! POSIX working-directory calls, for run_dir (see "the run directory" above).
  ! ------------------------------------------------------------------------------------
  interface
    integer(c_int) function c_chdir(path) bind(C, name="chdir")
      import :: c_int, c_char
      character(kind=c_char), intent(in) :: path(*)
    end function c_chdir
    type(c_ptr) function c_getcwd(buf, bufsize) bind(C, name="getcwd")
      import :: c_ptr, c_char, c_size_t
      character(kind=c_char), intent(out) :: buf(*)
      integer(c_size_t), value            :: bufsize
    end function c_getcwd
  end interface

  ! ------------------------------------------------------------------------------------
  ! One coupled MODFLOW 6 model.
  ! ------------------------------------------------------------------------------------
  type :: mf6_coupler_type
    private

    ! ---- configuration, from the &coupler namelist ----
    character(len=256)  :: mf6_model_name     = ''
    character(len=256)  :: rch_package_name   = 'RCHA'
    character(len=256)  :: bflow_package_name = 'CHD'
    character(len=256)  :: map_file           = ''
    integer             :: mf6_epsg           = 0
    character(len=1024) :: run_dir            = '.'      ! directory holding mfsim.nam ('.' = process cwd)
    character(len=1024) :: saved_dir          = ''       ! cwd to return to, while inside run_dir

    ! ---- what the driver needs to know; set by mf6_init ----
    logical, public :: feedback   = .true.   ! .false. => one-way (SUMMA drainage -> MODFLOW only)
    logical, public :: have_sy    = .false.  ! STO/SY found, so aquifer storage can be fed back
    logical, public :: have_bflow = .false.  ! bflow package found, so baseflow can be fed back

    ! ---- reprojection state for the built-in nearest-cell map (nHRU>1 only; see nearest_hru) ----
    integer :: utm_zone = 0
    logical :: utm_north = .true., use_utm = .false.

    ! ---- SUMMA-side geometry, kept because gather_head_to_hru needs it every step ----
    integer                       :: nHRU = 0
    double precision, allocatable :: hru_x(:), hru_y(:), hru_z(:)  ! HRU centroid lon/lat and surface elevation
    double precision, allocatable :: soil_thk(:)                   ! per-HRU SUMMA soil-column thickness (m)
    real,             allocatable :: sy_hru(:)                     ! per-HRU MODFLOW specific yield (-), map-weighted

    ! ---- MODFLOW 6 side ----
    real(c_double), pointer :: mf6_head(:) => null()      ! GWF dependent variable  <MODEL>/X       (REDUCED nodes)
    real(c_double), pointer :: mf6_rch(:)  => null()      ! RCH recharge array      <MODEL>/<RCH>/RECHARGE  (user cells, layer 1)
    integer(c_int), pointer :: mf6_mshape(:) => null()    ! DIS grid shape          <MODEL>/DIS/MSHAPE (nlay,nrow,ncol)
    integer(c_int), pointer :: mf6_nodered(:) => null()   ! DIS full->reduced node map <MODEL>/DIS/NODEREDUCED (present only when the grid is reduced)
    real(c_double), pointer :: mf6_sy(:) => null()        ! STO specific yield      <MODEL>/STO/SY  (REDUCED nodes)
    real(c_double), pointer :: bnd_simvals(:) => null()   ! boundary pkg simulated flow <MODEL>/<BPKG>/SIMVALS (m3 s-1, + = into aquifer)
    integer(c_int), pointer :: bnd_nodelist(:) => null()  ! boundary pkg cell list  <MODEL>/<BPKG>/NODELIST (REDUCED node numbers)
    logical                 :: grid_reduced = .false.
    real(c_double), pointer :: cellx(:) => null(), celly(:) => null()
    real(c_double), pointer :: xorigin => null(), yorigin => null(), angrot => null()
    integer                 :: nlay = 0, nrow = 0, ncol = 0

    ! ---- HRU -> MODFLOW cell mapping (sparse) ----
    ! for HRU i, its contributing horizontal cells are map_cell(map_ptr(i):map_ptr(i+1)-1)
    ! with weights map_wgt(...) that sum to 1
    integer, allocatable        :: map_ptr(:)   ! size nHRU+1
    integer, allocatable        :: map_cell(:)  ! horizontal cell index  (irow-1)*ncol + icol
    real,    allocatable        :: map_wgt(:)
    real(c_double), allocatable :: cell_area(:) ! horizontal cell plan area (m2), for the recharge volume split

  contains
    procedure, public :: init       => mf6_init
    procedure, public :: step       => mf6_step
    procedure, public :: finalize   => mf6_finalize_coupler
    procedure, public :: grid_shape => mf6_grid_shape
    ! internals
    procedure, private :: read_config           => mf6_read_config
    procedure, private :: enter_run_dir         => mf6_enter_run_dir
    procedure, private :: leave_run_dir         => mf6_leave_run_dir
    procedure, private :: check_model           => mf6_check_model
    procedure, private :: check_hru_elevation   => mf6_check_hru_elevation
    procedure, private :: build_map             => mf6_build_map
    procedure, private :: build_nearest_cell_map=> mf6_build_nearest_cell_map
    procedure, private :: read_map_file         => mf6_read_map_file
    procedure, private :: build_sy_hru          => mf6_build_sy_hru
    procedure, private :: scatter_drainage_to_rch => mf6_scatter_drainage_to_rch
    procedure, private :: gather_head_to_hru    => mf6_gather_head_to_hru
    procedure, private :: gather_aquifer_to_hru => mf6_gather_aquifer_to_hru
    procedure, private :: nearest_hru           => mf6_nearest_hru
    procedure, private :: top_active_node       => mf6_top_active_node
    procedure, private :: to_reduced            => mf6_to_reduced
  end type mf6_coupler_type

contains

  ! ==================================================================================
  ! Read the &coupler namelist, start MODFLOW 6, and build the HRU -> cell map.
  !
  ! The HRU arrays are in the driver's own HRU order; everything handed back later by
  ! mf6_step is in that same order.
  ! ==================================================================================
  subroutine mf6_init(this, config_file, run_dir, nHRU, hru_x, hru_y, hru_z, soil_thk, &
                      nSummaSteps, summa_data_step, err, message)
    class(mf6_coupler_type), intent(inout) :: this
    character(len=*),        intent(in)    :: config_file       ! &coupler namelist file
    character(len=*),        intent(in)    :: run_dir           ! directory holding mfsim.nam ('.' = process cwd)
    integer,                 intent(in)    :: nHRU
    double precision,        intent(in)    :: hru_x(:), hru_y(:), hru_z(:)  ! HRU lon/lat and surface elevation
    double precision,        intent(in)    :: soil_thk(:)       ! SUMMA soil-column depth per HRU (m)
    integer,                 intent(in)    :: nSummaSteps       ! number of SUMMA data steps in the run
    double precision,        intent(in)    :: summa_data_step   ! length of a SUMMA data step (s)
    integer,                 intent(out)   :: err
    character(len=*),        intent(out)   :: message
    integer        :: nred_len, nvar
    integer        :: istat
    real(c_double) :: tend_mf6

    err = 0; message = ''

    ! -- configuration --
    this%run_dir = run_dir
    call this%read_config(config_file, err, message)
    if (err /= 0) return

    ! -- SUMMA-side geometry (kept: gather_head_to_hru needs hru_z and soil_thk every step) --
    this%nHRU = nHRU
    allocate(this%hru_x(nHRU), this%hru_y(nHRU), this%hru_z(nHRU), this%soil_thk(nHRU))
    this%hru_x    = hru_x(1:nHRU)
    this%hru_y    = hru_y(1:nHRU)
    this%hru_z    = hru_z(1:nHRU)
    this%soil_thk = soil_thk(1:nHRU)

    call this%enter_run_dir(err, message)
    if (err /= 0) return

    ! -- initialize MODFLOW 6 (reads mfsim.nam from the working directory) --
    istat = mf6_initialize()
    if (istat /= BMI_OK) then
      message = 'MODFLOW 6 initialize failed'; err = 20
      call this%leave_run_dir(); return
    end if

    ! -- DIS grid dimensions, then pointers to the MODFLOW 6 arrays we exchange --
    call this%check_model(err, message)
    if (err /= 0) then; call this%leave_run_dir(); return; end if

    call mf6_ptr_int_arr(trim(this%mf6_model_name)//'/DIS/MSHAPE', this%mf6_mshape, 3, err, message)
    if (err /= 0) then; call this%leave_run_dir(); return; end if
    this%nlay = this%mf6_mshape(1); this%nrow = this%mf6_mshape(2); this%ncol = this%mf6_mshape(3)

    ! NB: the element count is taken into a local first - passing mf6_var_count(..., err, message)
    ! straight into a call that also passes err/message would alias those arguments
    nvar = mf6_var_count(trim(this%mf6_model_name)//'/X', err, message)
    call mf6_ptr_double(trim(this%mf6_model_name)//'/X', this%mf6_head, nvar, err, message)
    if (err /= 0) then; call this%leave_run_dir(); return; end if
    nvar = mf6_var_count(trim(this%mf6_model_name)//'/'//trim(this%rch_package_name)//'/RECHARGE', err, message)
    call mf6_ptr_double(trim(this%mf6_model_name)//'/'//trim(this%rch_package_name)//'/RECHARGE', &
                        this%mf6_rch, nvar, err, message)
    if (err /= 0) then; call this%leave_run_dir(); return; end if
    call mf6_ptr_double(trim(this%mf6_model_name)//'/DIS/CELLX', this%cellx, this%ncol, err, message)
    if (err /= 0) then; call this%leave_run_dir(); return; end if
    call mf6_ptr_double(trim(this%mf6_model_name)//'/DIS/CELLY', this%celly, this%nrow, err, message)
    if (err /= 0) then; call this%leave_run_dir(); return; end if
    call mf6_ptr_scalar(trim(this%mf6_model_name)//'/DIS/XORIGIN', this%xorigin, err, message)
    if (err /= 0) then; call this%leave_run_dir(); return; end if
    call mf6_ptr_scalar(trim(this%mf6_model_name)//'/DIS/YORIGIN', this%yorigin, err, message)
    if (err /= 0) then; call this%leave_run_dir(); return; end if
    call mf6_ptr_scalar(trim(this%mf6_model_name)//'/DIS/ANGROT',  this%angrot, err, message)
    if (err /= 0) then; call this%leave_run_dir(); return; end if

    ! -- full(user) -> reduced node map: DIS/NODEREDUCED has length nlay*nrow*ncol
    !    when IDOMAIN removes cells, else length 1 (grid not reduced; node == user node)
    nred_len = mf6_var_count(trim(this%mf6_model_name)//'/DIS/NODEREDUCED', err, message)
    if (err /= 0) then; call this%leave_run_dir(); return; end if
    this%grid_reduced = (nred_len == this%nlay*this%nrow*this%ncol)
    if (this%grid_reduced) then
      call mf6_ptr_int_arr(trim(this%mf6_model_name)//'/DIS/NODEREDUCED', this%mf6_nodered, nred_len, err, message)
      if (err /= 0) then; call this%leave_run_dir(); return; end if
    end if

    ! -- optional groundwater-feedback pointers (only needed when feedback = .true.) --
    if (this%feedback) then
      allocate(this%sy_hru(nHRU)); this%sy_hru = 0.0
      this%have_sy = mf6_try_ptr_double(trim(this%mf6_model_name)//'/STO/SY', this%mf6_sy)
      if (len_trim(this%bflow_package_name) > 0) then
        this%have_bflow = &
          mf6_try_ptr_double(trim(this%mf6_model_name)//'/'//trim(this%bflow_package_name)//'/SIMVALS',  this%bnd_simvals) .and. &
          mf6_try_ptr_int   (trim(this%mf6_model_name)//'/'//trim(this%bflow_package_name)//'/NODELIST', this%bnd_nodelist)
        if (.not. this%have_bflow) write(*,'(a)') 'summa_modflow6: WARNING - boundary package "'// &
          trim(this%bflow_package_name)//'" not found; scalarAquiferBaseflow feedback disabled'
      end if
    end if

    ! -- coupling time-step consistency --
    ! NB: MODFLOW's delt is only set once the first time step is prepared, so it is
    ! not meaningful yet right after initialize(); the delt == data_step check is
    ! therefore deferred to the first call of mf6_step (see below).
    istat = mf6_get_end_time(tend_mf6)
    if (tend_mf6 < nSummaSteps*real(summa_data_step, c_double)*(1.0_c_double - 1.0e-6_c_double)) then
      write(*,'(a)') 'summa_modflow6: WARNING - MODFLOW simulation is shorter than the SUMMA simulation'
    end if

    ! -- build the HRU -> MODFLOW horizontal-cell weight map --
    call this%build_map(err, message)
    if (err /= 0) then; call this%leave_run_dir(); return; end if

    ! -- per-HRU specific yield (map-weighted), fixed for the run --
    if (this%feedback .and. this%have_sy) call this%build_sy_hru

    call this%leave_run_dir()
  end subroutine mf6_init

  ! ==================================================================================
  ! One coupling step: steps 3 to 5 of the exchange described at the top of this module.
  ! The driver has already pushed the previous step's feedback into SUMMA and advanced it.
  !
  ! head_hru/stor_hru/bflow_hru are the driver's own lagged feedback arrays: they are read
  ! back here for the next step, and are left untouched when feedback is off.
  ! ==================================================================================
  subroutine mf6_step(this, istep, summa_data_step, drain_hru, head_hru, stor_hru, bflow_hru, err, message)
    class(mf6_coupler_type), intent(inout) :: this
    integer,                 intent(in)    :: istep             ! 1-based coupling step index
    double precision,        intent(in)    :: summa_data_step   ! length of a SUMMA data step (s)
    real,                    intent(in)    :: drain_hru(:)      ! per-HRU soil drainage (m s-1)
    real,                    intent(inout) :: head_hru(:)       ! per-HRU prescribed head (m, matric head at soil base)
    real,                    intent(inout) :: stor_hru(:)       ! per-HRU relative aquifer storage (m)
    real,                    intent(inout) :: bflow_hru(:)      ! per-HRU aquifer baseflow (m s-1, + = out of aquifer)
    integer,                 intent(out)   :: err
    character(len=*),        intent(out)   :: message
    integer        :: istat
    real(c_double) :: dt_mf6

    err = 0; message = ''

    call this%enter_run_dir(err, message)
    if (err /= 0) return

    ! 3. SUMMA soil drainage -> MODFLOW 6 RCH recharge array
    call this%scatter_drainage_to_rch(drain_hru)

    ! 4. advance MODFLOW 6 one coupling step
    istat = mf6_prepare_time_step(real(summa_data_step, c_double))
    if (istep == 1) then
      ! delt is now set: verify one MODFLOW time step == one SUMMA data step
      istat = mf6_get_time_step(dt_mf6)
      if (abs(dt_mf6 - real(summa_data_step, c_double)) > 1.0e-6_c_double*real(summa_data_step, c_double)) then
        write(message,'(a,g0,a,g0,a)') 'MODFLOW time step (', dt_mf6, &
              ') must equal the SUMMA data step (', summa_data_step, ')'
        err = 20
        call this%leave_run_dir(); return
      end if
    end if
    istat = mf6_do_time_step()
    istat = mf6_finalize_time_step()

    ! 5. read the new MODFLOW state, aggregate per HRU for the next iteration
    if (this%feedback) then
      call this%gather_head_to_hru(head_hru)
      call this%gather_aquifer_to_hru(head_hru, stor_hru, bflow_hru)
    end if

    call this%leave_run_dir()
  end subroutine mf6_step

  ! ==================================================================================
  ! Shut MODFLOW 6 down and release this coupler's state, so the object can be
  ! initialized again for the next coupled run (one per calibration sample).
  ! ==================================================================================
  subroutine mf6_finalize_coupler(this, err, message)
    class(mf6_coupler_type), intent(inout) :: this
    integer,                 intent(out)   :: err
    character(len=*),        intent(out)   :: message
    integer :: istat

    err = 0; message = ''

    call this%enter_run_dir(err, message)
    if (err == 0) then
      istat = mf6_finalize()
      call this%leave_run_dir()
      if (istat /= BMI_OK) then; message = 'MODFLOW 6 finalize failed'; err = 20; end if
    end if

    ! release everything owned here; the MODFLOW-side pointers belong to libmf6, which has
    ! just freed them, so drop them rather than leave them dangling
    if (allocated(this%hru_x))     deallocate(this%hru_x)
    if (allocated(this%hru_y))     deallocate(this%hru_y)
    if (allocated(this%hru_z))     deallocate(this%hru_z)
    if (allocated(this%soil_thk))  deallocate(this%soil_thk)
    if (allocated(this%sy_hru))    deallocate(this%sy_hru)
    if (allocated(this%map_ptr))   deallocate(this%map_ptr)
    if (allocated(this%map_cell))  deallocate(this%map_cell)
    if (allocated(this%map_wgt))   deallocate(this%map_wgt)
    if (allocated(this%cell_area)) deallocate(this%cell_area)
    this%mf6_head     => null()
    this%mf6_rch      => null()
    this%mf6_mshape   => null()
    this%mf6_nodered  => null()
    this%mf6_sy       => null()
    this%bnd_simvals  => null()
    this%bnd_nodelist => null()
    this%cellx        => null()
    this%celly        => null()
    this%xorigin      => null()
    this%yorigin      => null()
    this%angrot       => null()
    this%grid_reduced = .false.
    this%have_sy      = .false.
    this%have_bflow   = .false.
    this%nHRU         = 0
  end subroutine mf6_finalize_coupler

  ! ==================================================================================
  ! Give one model instance its own copy of a MODFLOW 6 model directory.
  !
  ! libmf6 writes its listing, head and budget files into the working directory, so several
  ! instances sharing one directory overwrite each other's output.  The calibration driver
  ! runs one instance per MPI rank, so each rank gets its own directory here, populated once
  ! from the master model directory and reused for every parameter sample afterwards.
  !
  ! The contents are copied rather than symlinked on purpose: MODFLOW opens its output files
  ! for writing, and writing through a symlink would land back in the master directory, which
  ! is the collision this is meant to prevent.
  ! ==================================================================================
  subroutine mf6_prepare_run_dir(master_dir, run_dir, err, message)
    character(len=*), intent(in)  :: master_dir   ! directory holding the MODFLOW 6 mfsim.nam
    character(len=*), intent(in)  :: run_dir      ! this instance's own directory, created if absent
    integer,          intent(out) :: err
    character(len=*), intent(out) :: message
    logical :: ready
    integer :: rc

    err = 0; message = ''

    ! already populated by an earlier sample on this rank: nothing to do
    inquire(file=trim(run_dir)//'/mfsim.nam', exist=ready)
    if (ready) return

    inquire(file=trim(master_dir)//'/mfsim.nam', exist=ready)
    if (.not. ready) then
      message = 'no mfsim.nam in the MODFLOW run directory '//trim(master_dir); err = 20; return
    end if

    call execute_command_line('mkdir -p "'//trim(run_dir)//'"', exitstat=rc)
    if (rc /= 0) then
      message = 'cannot create the MODFLOW run directory '//trim(run_dir); err = 20; return
    end if

    ! -L so a master directory that is itself built from symlinks gives real files here
    call execute_command_line('cp -RL "'//trim(master_dir)//'"/. "'//trim(run_dir)//'"/', exitstat=rc)
    if (rc /= 0) then
      message = 'cannot copy the MODFLOW model from '//trim(master_dir)//' into '//trim(run_dir)
      err = 20; return
    end if
  end subroutine mf6_prepare_run_dir

  ! ==================================================================================
  ! The MODFLOW 6 DIS grid dimensions, for the driver's start-up report.
  ! ==================================================================================
  subroutine mf6_grid_shape(this, nlay, nrow, ncol)
    class(mf6_coupler_type), intent(in)  :: this
    integer,                 intent(out) :: nlay, nrow, ncol
    nlay = this%nlay; nrow = this%nrow; ncol = this%ncol
  end subroutine mf6_grid_shape

  ! ==================================================================================
  ! Read the &coupler namelist.  Namelist groups cannot name derived-type components,
  ! so it is read into locals and copied in.
  ! ==================================================================================
  subroutine mf6_read_config(this, config_file, err, message)
    class(mf6_coupler_type), intent(inout) :: this
    character(len=*),        intent(in)    :: config_file
    integer,                 intent(out)   :: err
    character(len=*),        intent(out)   :: message
    integer :: fu, rc
    character(len=256) :: mf6_model_name, rch_package_name, bflow_package_name, map_file
    integer            :: mf6_epsg
    logical            :: feedback
    namelist /coupler/ mf6_model_name, rch_package_name, bflow_package_name, map_file, mf6_epsg, feedback

    err = 0; message = ''

    ! namelist defaults
    mf6_model_name     = ''
    rch_package_name   = 'RCHA'
    bflow_package_name = 'CHD'
    map_file           = ''
    mf6_epsg           = 0
    feedback           = .true.

    open(action='read', file=trim(config_file), iostat=rc, newunit=fu)
    if (rc /= 0) then
      message = 'cannot open configuration file '//trim(config_file); err = 20; return
    end if
    read(nml=coupler, iostat=rc, unit=fu)
    close(fu)
    if (rc /= 0) then
      message = 'error reading &coupler namelist from '//trim(config_file); err = 20; return
    end if
    if (len_trim(mf6_model_name) == 0) then
      message = 'mf6_model_name must be set in '//trim(config_file); err = 20; return
    end if

    this%mf6_model_name     = mf6_model_name
    this%rch_package_name   = rch_package_name
    this%bflow_package_name = bflow_package_name
    this%map_file           = map_file
    this%mf6_epsg           = mf6_epsg
    this%feedback           = feedback

    call to_upper(this%mf6_model_name)
    call to_upper(this%rch_package_name)
    call to_upper(this%bflow_package_name)

    if (this%mf6_epsg /= 0) then
      this%use_utm = utm_zone_from_epsg(this%mf6_epsg, this%utm_zone, this%utm_north)
      if (.not. this%use_utm) write(*,'(a,i0,a)') 'summa_modflow6: WARNING - mf6_epsg=', this%mf6_epsg, &
        ' is not a supported WGS84 UTM code (32601-32660 N, 32701-32760 S); ignoring it.'
    end if

    ! NOTE: map_file is left as written.  It is read inside the run directory (see mf6_build_map,
    !       which runs between enter_run_dir and leave_run_dir), so a relative name resolves
    !       against the MODFLOW model's own directory either way - which is where it lives.
  end subroutine mf6_read_config

  ! ==================================================================================
  ! Enter/leave the MODFLOW run directory.  A run_dir of '.' means "wherever the process
  ! already is", and then nothing happens at all - which is what the standalone couplers
  ! pass, so their behaviour is exactly as it was before this module existed.
  ! ==================================================================================
  subroutine mf6_enter_run_dir(this, err, message)
    class(mf6_coupler_type), intent(inout) :: this
    integer,                 intent(out)   :: err
    character(len=*),        intent(out)   :: message
    character(kind=c_char) :: buf(4096)
    integer :: i

    err = 0; message = ''
    if (trim(this%run_dir) == '.') return

    ! remember where to come back to
    if (.not. c_associated(c_getcwd(buf, int(size(buf), c_size_t)))) then
      message = 'cannot determine the current working directory'; err = 20; return
    end if
    this%saved_dir = ''
    do i = 1, size(buf)
      if (buf(i) == c_null_char) exit
      this%saved_dir(i:i) = buf(i)
    end do

    if (c_chdir(cstr(trim(this%run_dir))) /= 0) then
      message = 'cannot change into the MODFLOW run directory '//trim(this%run_dir); err = 20; return
    end if
  end subroutine mf6_enter_run_dir

  subroutine mf6_leave_run_dir(this)
    class(mf6_coupler_type), intent(inout) :: this
    integer :: rc
    if (trim(this%run_dir) == '.') return
    if (len_trim(this%saved_dir) == 0) return
    rc = c_chdir(cstr(trim(this%saved_dir)))
  end subroutine mf6_leave_run_dir

  ! ==================================================================================
  ! SUMMA drainage (m s-1, per HRU)  ->  MODFLOW RECHARGE array (m s-1, per top cell).
  ! Each HRU's flux is spread over its mapped cells by weight; a cell that receives
  ! from several HRUs gets the area-weighted mean flux (so recharge volume is conserved
  ! when the HRU areas equal the covered cell areas).
  ! ==================================================================================
  subroutine mf6_scatter_drainage_to_rch(this, drain_hru)
    class(mf6_coupler_type), intent(inout) :: this
    real,                    intent(in)    :: drain_hru(:)
    real(c_double), allocatable :: num(:), den(:)
    integer :: i, k, c

    allocate(num(this%nrow*this%ncol), den(this%nrow*this%ncol))
    num = 0.0_c_double
    den = 0.0_c_double
    do i = 1, this%nHRU
      do k = this%map_ptr(i), this%map_ptr(i+1) - 1
        c = this%map_cell(k)
        if (c < 1 .or. c > this%nrow*this%ncol) cycle
        num(c) = num(c) + real(this%map_wgt(k), c_double) * this%cell_area(c) * real(drain_hru(i), c_double)
        den(c) = den(c) + real(this%map_wgt(k), c_double) * this%cell_area(c)
      end do
    end do

    ! RECHARGE is indexed by horizontal cell (row-major) for READASARRAYS RCH
    do c = 1, min(size(this%mf6_rch), this%nrow*this%ncol)
      if (den(c) > 0.0_c_double) then
        this%mf6_rch(c) = num(c) / den(c)
      else
        this%mf6_rch(c) = 0.0_c_double
      end if
    end do
    deallocate(num, den)
  end subroutine mf6_scatter_drainage_to_rch

  ! ==================================================================================
  ! MODFLOW head (m, per node)  ->  SUMMA prescribed lower-BC matric head (m, per HRU).
  ! The head is sampled at the top-most active node of each mapped column and converted
  ! to a matric (pressure) head at the base of the SUMMA soil column:
  !     lowerBoundHead = h_mf6 - (z_surface_HRU - soil_thk(HRU))
  ! where soil_thk is SUMMA's own per-HRU soil-column depth.
  ! positive => saturated / positive pressure head at the soil base.
  ! ==================================================================================
  subroutine mf6_gather_head_to_hru(this, head_hru)
    class(mf6_coupler_type), intent(inout) :: this
    real,                    intent(inout) :: head_hru(:)
    integer :: i, k, c, node
    real(c_double) :: hsum, wsum, zbase

    do i = 1, this%nHRU
      hsum = 0.0_c_double
      wsum = 0.0_c_double
      do k = this%map_ptr(i), this%map_ptr(i+1) - 1
        c = this%map_cell(k)                ! horizontal (row-major) cell index
        node = this%top_active_node(c)
        if (node > 0) then
          hsum = hsum + real(this%map_wgt(k), c_double) * this%mf6_head(node)
          wsum = wsum + real(this%map_wgt(k), c_double)
        end if
      end do
      if (wsum > 0.0_c_double) then
        zbase = real(this%hru_z(i), c_double) - real(this%soil_thk(i), c_double)
        head_hru(i) = real(hsum / wsum - zbase)
      end if
    end do
  end subroutine mf6_gather_head_to_hru

  ! ==================================================================================
  ! MODFLOW 6  ->  SUMMA aquifer bookkeeping (per HRU), only when feedback=.true.:
  !   scalarAquiferStorage  = Sy * (water table - soil-column base) = sy_hru * head_hru (m)
  !   scalarAquiferBaseflow = -(sum of <bflow_package> flow over the HRU's mapped cells)
  !                           divided by the mapped-cell area  (m s-1, + = out of aquifer)
  ! ==================================================================================
  subroutine mf6_gather_aquifer_to_hru(this, head_hru, stor_hru, bflow_hru)
    class(mf6_coupler_type), intent(inout) :: this
    real,                    intent(in)    :: head_hru(:)
    real,                    intent(inout) :: stor_hru(:)
    real,                    intent(inout) :: bflow_hru(:)
    integer :: i, k, c, kb, nr
    real(c_double) :: fnum, aden

    if (this%have_sy) then
      do i = 1, this%nHRU
        stor_hru(i) = this%sy_hru(i) * head_hru(i)
      end do
    end if

    if (this%have_bflow) then
      do i = 1, this%nHRU
        fnum = 0.0_c_double     ! signed boundary flow over the HRU's mapped cells (m3 s-1, + into aquifer)
        aden = 0.0_c_double     ! matching mapped-cell plan area (m2)
        do k = this%map_ptr(i), this%map_ptr(i+1) - 1
          c = this%map_cell(k)
          if (c < 1 .or. c > this%nrow*this%ncol) cycle
          nr = this%to_reduced(c)                  ! boundary NODELIST holds reduced node numbers
          if (nr > 0) then
            do kb = 1, size(this%bnd_simvals)
              if (this%bnd_nodelist(kb) == nr) &
                fnum = fnum + real(this%map_wgt(k), c_double) * this%bnd_simvals(kb)
            end do
          end if
          aden = aden + real(this%map_wgt(k), c_double) * this%cell_area(c)
        end do
        if (aden > 0.0_c_double) bflow_hru(i) = real(-fnum / aden)
      end do
    end if
  end subroutine mf6_gather_aquifer_to_hru

  ! Map-weighted MODFLOW specific yield per HRU (fixed for the run).
  subroutine mf6_build_sy_hru(this)
    class(mf6_coupler_type), intent(inout) :: this
    integer :: i, k, c, nr
    real(c_double) :: wnum, wden
    do i = 1, this%nHRU
      wnum = 0.0_c_double; wden = 0.0_c_double
      do k = this%map_ptr(i), this%map_ptr(i+1) - 1
        c = this%map_cell(k)
        if (c < 1 .or. c > this%nrow*this%ncol) cycle
        nr = this%to_reduced(c)
        if (nr < 1 .or. nr > size(this%mf6_sy)) cycle
        wnum = wnum + real(this%map_wgt(k), c_double) * this%mf6_sy(nr)
        wden = wden + real(this%map_wgt(k), c_double)
      end do
      if (wden > 0.0_c_double) this%sy_hru(i) = real(wnum / wden)
    end do
  end subroutine mf6_build_sy_hru

  ! Walk down horizontal column choriz and return the REDUCED (solution) node number
  ! of the first layer with a usable head, or -1 if inactive/dry everywhere.
  integer function mf6_top_active_node(this, choriz) result(node)
    class(mf6_coupler_type), intent(in) :: this
    integer,                 intent(in) :: choriz    ! (irow-1)*ncol + icol
    integer :: ilay, nr
    node = -1
    do ilay = 1, this%nlay
      nr = this%to_reduced((ilay-1)*this%nrow*this%ncol + choriz)
      if (nr >= 1 .and. nr <= size(this%mf6_head)) then
        if (this%mf6_head(nr) > -1.0e29_c_double) then   ! HNOFLO / HDRY guard
          node = nr
          return
        end if
      end if
    end do
  end function mf6_top_active_node

  ! full (user) node number -> reduced (solution) node number; <= 0 if the cell is
  ! inactive.  Identity when IDOMAIN removes no cells (grid not reduced).
  integer function mf6_to_reduced(this, nfull) result(nred)
    class(mf6_coupler_type), intent(in) :: this
    integer,                 intent(in) :: nfull
    if (this%grid_reduced) then
      if (nfull >= 1 .and. nfull <= size(this%mf6_nodered)) then
        nred = this%mf6_nodered(nfull)
      else
        nred = -1
      end if
    else
      nred = nfull
    end if
  end function mf6_to_reduced

  ! ==================================================================================
  ! Build the HRU -> horizontal-cell weight map, either from an external weight file
  ! or, by default, as a nearest-cell assignment using the MODFLOW 6 DIS geometry.
  ! ==================================================================================
  subroutine mf6_build_map(this, err, message)
    class(mf6_coupler_type), intent(inout) :: this
    integer,                 intent(out)   :: err
    character(len=*),        intent(out)   :: message
    integer :: i, j, c
    real(c_double) :: dx, dy

    err = 0; message = ''

    allocate(this%cell_area(this%nrow*this%ncol))
    do i = 1, this%nrow
      do j = 1, this%ncol
        c = (i-1)*this%ncol + j
        ! CELLX/CELLY are grid-local centres; DELR = spacing in x, DELC in y.
        ! Approximate the plan area from local centre spacing (uniform grids exact).
        dx = cell_spacing(this%cellx, j, this%ncol)
        dy = cell_spacing(this%celly, i, this%nrow)
        this%cell_area(c) = dx * dy
      end do
    end do

    if (len_trim(this%map_file) > 0) then
      call this%read_map_file(err, message)
    else
      call this%build_nearest_cell_map
    end if
    if (err /= 0) return

    call this%check_hru_elevation(err, message)
  end subroutine mf6_build_map

  ! The coupler makes a handful of assumptions about the MODFLOW 6 model it is handed, and every
  ! one of them fails silently rather than loudly if it is wrong: fluxes off by a fixed factor,
  ! a water table interpreted in the wrong units, recharge landing on the wrong cells.  Check
  ! them here, right after mf6_initialize, before any of it can be mistaken for a bad simulation.
  ! utils/test/test_mflow/README.md states the same requirements for whoever builds the MODFLOW model.
  subroutine mf6_check_model(this, err, message)
    class(mf6_coupler_type), intent(inout) :: this
    integer,                 intent(out)   :: err
    character(len=*),        intent(out)   :: message
    integer(c_int), pointer :: itmuni(:) => null(), lenuni(:) => null(), mshape(:) => null()
    real(c_double), pointer :: rch(:) => null()
    integer(c_int) :: nbytes, isize
    integer :: ncpl
    character(len=12), parameter :: tunit(0:5) = [ character(len=12) :: &
        'undefined', 'seconds', 'minutes', 'hours', 'days', 'years' ]
    character(len=12), parameter :: lunit(0:3) = [ character(len=12) :: &
        'undefined', 'feet', 'meters', 'centimeters' ]

    err = 0; message = ''

    ! -- time unit: SUMMA's data_step and every flux exchanged here are per second --
    if (mf6_try_ptr_int('TDIS/ITMUNI', itmuni)) then
      if (itmuni(1) == 0) then
        write(*,'(a)') 'summa_modflow6: WARNING - TDIS has no TIME_UNITS; assuming SECONDS. '// &
              'The MODFLOW delt check below is then the only guard on the time unit.'
      else if (itmuni(1) /= 1) then
        message = 'TDIS TIME_UNITS must be SECONDS, not '//trim(tunit(min(max(itmuni(1),0),5)))
        err = 20; return
      end if
    end if

    ! -- length unit: heads, elevations and cell geometry are all read as metres --
    if (mf6_try_ptr_int(trim(this%mf6_model_name)//'/DIS/LENUNI', lenuni)) then
      if (lenuni(1) == 0) then
        write(*,'(a)') 'summa_modflow6: WARNING - DIS has no LENGTH_UNITS; assuming METERS. '// &
              'HRU elevations are checked against DIS TOP below, which will catch feet.'
      else if (lenuni(1) /= 2) then
        message = 'DIS LENGTH_UNITS must be METERS, not '//trim(lunit(min(max(lenuni(1),0),3)))
        err = 20; return
      end if
    end if

    ! -- a single GWF model: anything else in the simulation is simply not coupled --
    if (mf6_get_var_nbytes(cstr('__INPUT__/SIM/NAM/MTYPE'), nbytes) == BMI_OK .and. &
        mf6_get_var_itemsize(cstr('__INPUT__/SIM/NAM/MTYPE'), isize) == BMI_OK) then
      if (isize > 0 .and. nbytes/isize > 1) &
        write(*,'(a,i0,a)') 'summa_modflow6: WARNING - mfsim.nam lists ', nbytes/isize, &
              ' models; only '//trim(this%mf6_model_name)//' is coupled to SUMMA'
    end if

    ! -- discretised with DIS: MSHAPE is (nlay,nrow,ncol) for DIS, shorter for DISV/DISU, and
    !    the HRU->cell map, the nearest-cell map and the RCH array all index a structured grid
    if (.not. mf6_try_ptr_int(trim(this%mf6_model_name)//'/DIS/MSHAPE', mshape)) then
      message = 'no DIS grid found for GWF model '//trim(this%mf6_model_name)// &
            ' - check mf6_model_name in the config (upper case, as in mfsim.nam); the coupler '// &
            'requires a single GWF model discretised with DIS, not DISV or DISU'
      err = 20; return
    end if
    if (size(mshape) /= 3) then
      write(message,'(a,i0,a)') trim(this%mf6_model_name)//' is not a DIS grid '// &
            '(MSHAPE has ', size(mshape), ' entries, DIS has 3); DISV and DISU are not supported'
      err = 20; return
    end if
    ncpl = mshape(2) * mshape(3)

    ! -- RCH with READASARRAYS: RECHARGE is then one value per horizontal cell, which is what
    !    scatter_drainage_to_rch writes into; a list-based RCH gives one value per list entry
    if (.not. mf6_try_ptr_double(trim(this%mf6_model_name)//'/'//trim(this%rch_package_name)// &
                                 '/RECHARGE', rch)) then
      message = 'no RCH package '//trim(this%rch_package_name)//' in GWF model '// &
            trim(this%mf6_model_name)//' - check rch_package_name in the config (upper case, as in '// &
            'the GWF name file)'
      err = 20; return
    end if
    if (size(rch) /= ncpl) then
      write(message,'(a,i0,a,i0,a)') 'RCH package '//trim(this%rch_package_name)// &
            ' must be declared READASARRAYS: its RECHARGE array has ', size(rch), &
            ' entries, not one per horizontal cell (', ncpl, ')'
      err = 20; return
    end if
  end subroutine mf6_check_model

  ! Each HRU's SUMMA elevation must be the land surface of the MODFLOW cells it maps to, because
  ! gather_head_to_hru forms  lowerBoundHead = h_mf6 - (z_surface_HRU - soil_thickness).  If the two
  ! disagree the soil column is handed a water table that is metres above or below it, and SUMMA
  ! fails to converge with no indication of why, so check it up front rather than leave it to chance.
  subroutine mf6_check_hru_elevation(this, err, message)
    class(mf6_coupler_type), intent(inout) :: this
    integer,                 intent(out)   :: err
    character(len=*),        intent(out)   :: message
    real(c_double), pointer :: mf6_top(:) => null()
    real(c_double) :: zsum, wsum, ztop
    integer :: i, k, c, node, nbad, nvar
    real(c_double), parameter :: z_tol = 10.0_c_double   ! m, generous: catches wrong-band mistakes

    err = 0; message = ''

    nvar = mf6_var_count(trim(this%mf6_model_name)//'/DIS/TOP', err, message)
    call mf6_ptr_double(trim(this%mf6_model_name)//'/DIS/TOP', mf6_top, nvar, err, message)
    if (err /= 0) return
    if (.not. associated(mf6_top)) return      ! nothing to check against
    nbad = 0
    do i = 1, this%nHRU
      zsum = 0.0_c_double; wsum = 0.0_c_double
      do k = this%map_ptr(i), this%map_ptr(i+1) - 1
        c = this%map_cell(k)
        node = this%top_active_node(c)       ! DIS/TOP is in reduced node numbering, as DIS/X is
        if (node < 1 .or. node > size(mf6_top)) cycle
        zsum = zsum + real(this%map_wgt(k), c_double) * mf6_top(node)
        wsum = wsum + real(this%map_wgt(k), c_double)
      end do
      if (wsum <= 0.0_c_double) cycle
      ztop = zsum / wsum
      if (abs(ztop - real(this%hru_z(i), c_double)) > z_tol) then
        nbad = nbad + 1
        write(*,'(a,i0,2(a,f10.2),a)') 'summa_modflow6: HRU ', i, ' elevation ', this%hru_z(i), &
          ' m disagrees with the mean land surface of its MODFLOW cells, ', ztop, ' m'
      end if
    end do
    if (nbad > 0) then
      message = 'set each HRU elevation in attributes.nc to the mean land '// &
        'surface of the cells it maps to, or fix the map_file; lowerBoundHead is formed from the '// &
        'difference, so a mismatch puts the water table off the soil column.'
      err = 20; return
    end if
  end subroutine mf6_check_hru_elevation

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

  ! Assign every ACTIVE MODFLOW cell to its nearest SUMMA HRU centre (a per-cell
  ! "which HRU owns me" / Thiessen assignment) so an HRU naturally ends up covering
  ! however many cells are nearest to it - unlike a per-HRU nearest-cell search,
  ! this needs no per-HRU cell count ahead of time and gives every active cell to
  ! exactly one HRU (matching what a hand-built map_file for a lumped HRU would do).
  !
  ! --- mapping HRUs to cells ---------------------------------------------------------------
  ! HRUs may each own a disjoint set of cells, down to one cell per HRU; they need not overlap
  ! and need not cover the whole grid.
  !
  ! An HRU's elevation in attributes.nc must be the mean land surface of the cells it maps to,
  ! since gather_head_to_hru forms
  !     lowerBoundHead = h_mf6 - (z_surface_HRU - soil_thickness)
  ! and a mismatch hands the soil column a water table metres above or below it.  That is
  ! checked against DIS/TOP in check_hru_elevation at start-up rather than left to chance.
  !
  ! CAVEAT: SUMMA HRU centres (hru_x, hru_y) are typically longitude/latitude
  ! (degrees) while the MODFLOW DIS grid (cellx, celly) is in a projected CRS
  ! (metres); comparing the two directly is only meaningful when nHRU==1 - the
  ! whole grid then trivially belongs to the one HRU, no distance comparison
  ! needed at all - or when hru_x/hru_y already happen to be in the grid's
  ! projected system. For nHRU>1 with mismatched coordinates this default is not
  ! reliable: supply map_file instead (a warning is printed if that looks likely).
  subroutine mf6_build_nearest_cell_map(this)
    class(mf6_coupler_type), intent(inout) :: this
    integer :: i, j, c, ih, nMapped
    integer, allocatable :: cell_count(:), fill(:)

    if (this%nHRU > 1 .and. .not. this%use_utm .and. &
        all(abs(this%hru_x) <= 180.0_c_double) .and. all(abs(this%hru_y) <= 90.0_c_double)) &
      write(*,'(a)') 'summa_modflow6: WARNING - nearest-cell auto-map with nHRU>1 and HRU '// &
        'coordinates that look like lon/lat: the MODFLOW grid is usually projected (metres), '// &
        'so this default cell assignment is unlikely to be meaningful. Set mf6_epsg (WGS84 UTM) '// &
        'or supply map_file instead.'

    allocate(cell_count(this%nHRU)); cell_count = 0
    nMapped = 0
    do i = 1, this%nrow
      do j = 1, this%ncol
        c = (i-1)*this%ncol + j
        if (this%to_reduced(c) <= 0) cycle   ! skip inactive cells (no HRU should collect from them)
        ih = this%nearest_hru(i, j)
        cell_count(ih) = cell_count(ih) + 1
        nMapped = nMapped + 1
      end do
    end do

    allocate(this%map_ptr(this%nHRU+1))
    this%map_ptr(1) = 1
    do ih = 1, this%nHRU
      this%map_ptr(ih+1) = this%map_ptr(ih) + cell_count(ih)
    end do
    allocate(this%map_cell(nMapped), this%map_wgt(nMapped))
    allocate(fill(this%nHRU)); fill = this%map_ptr(1:this%nHRU)

    do i = 1, this%nrow
      do j = 1, this%ncol
        c = (i-1)*this%ncol + j
        if (this%to_reduced(c) <= 0) cycle
        ih = this%nearest_hru(i, j)
        this%map_cell(fill(ih)) = c
        this%map_wgt(fill(ih))  = 1.0
        fill(ih) = fill(ih) + 1
      end do
    end do
  end subroutine mf6_build_nearest_cell_map

  ! Nearest SUMMA HRU centre to MODFLOW row i, col j (grid-local frame).  Trivially
  ! HRU 1 when there is only one HRU, so no coordinate comparison (and hence no CRS
  ! agreement) is needed in the common single-lumped-HRU case.  For nHRU>1, hru_x/
  ! hru_y (assumed WGS84 lon/lat) are reprojected to the grid's UTM zone first when
  ! mf6_epsg identifies one (use_utm); otherwise they are compared as-is.
  integer function mf6_nearest_hru(this, i, j) result(ibest)
    class(mf6_coupler_type), intent(in) :: this
    integer,                 intent(in) :: i, j
    integer :: ih
    real(c_double) :: wx, wy, gx, gy, ca, sa, lx, ly, dbest, dist
    if (this%nHRU == 1) then
      ibest = 1; return
    end if
    ca = cos(-real(this%angrot, c_double) * acos(-1.0_c_double) / 180.0_c_double)
    sa = sin(-real(this%angrot, c_double) * acos(-1.0_c_double) / 180.0_c_double)
    dbest = huge(dbest)
    ibest = 1
    do ih = 1, this%nHRU
      if (this%use_utm) then
        call lonlat_to_utm(this%hru_x(ih), this%hru_y(ih), this%utm_zone, this%utm_north, wx, wy)
      else
        wx = real(this%hru_x(ih), c_double); wy = real(this%hru_y(ih), c_double)
      end if
      ! world -> grid-local frame (translate by origin, rotate by -angrot)
      gx = wx - real(this%xorigin, c_double)
      gy = wy - real(this%yorigin, c_double)
      lx = ca*gx - sa*gy
      ly = sa*gx + ca*gy
      dist = (lx - this%cellx(j))**2 + (ly - this%celly(i))**2
      if (dist < dbest) then
        dbest = dist
        ibest = ih
      end if
    end do
  end function mf6_nearest_hru

  ! EPSG -> (UTM zone, hemisphere) for the WGS84 UTM series (32601-32660 N, 32701-32760 S).
  ! Returns .false. (zone/north undefined) for any other EPSG code.
  logical function utm_zone_from_epsg(epsg, zone, north) result(ok)
    integer, intent(in)  :: epsg
    integer, intent(out) :: zone
    logical, intent(out) :: north
    ok = .true.
    if (epsg >= 32601 .and. epsg <= 32660) then
      zone = epsg - 32600; north = .true.
    else if (epsg >= 32701 .and. epsg <= 32760) then
      zone = epsg - 32700; north = .false.
    else
      zone = 0; north = .true.; ok = .false.
    end if
  end function utm_zone_from_epsg

  ! WGS84 lon/lat (degrees) -> UTM easting/northing (m): the standard closed-form
  ! forward transverse-Mercator series (Snyder/USGS), accurate to well under a
  ! metre within a UTM zone - ample for nearest-cell matching.
  subroutine lonlat_to_utm(lon_deg, lat_deg, zone, north, easting, northing)
    double precision, intent(in) :: lon_deg, lat_deg
    integer,          intent(in) :: zone
    logical,          intent(in) :: north
    real(c_double),   intent(out) :: easting, northing
    real(c_double), parameter :: a  = 6378137.0_c_double              ! WGS84 semi-major axis (m)
    real(c_double), parameter :: f  = 1.0_c_double/298.257223563_c_double
    real(c_double), parameter :: k0 = 0.9996_c_double                 ! UTM scale factor
    real(c_double) :: pi, e2, ep2, lon0, phi, lam, rn, t, cc, am, m

    pi   = acos(-1.0_c_double)
    e2   = f*(2.0_c_double - f)
    ep2  = e2/(1.0_c_double - e2)
    lon0 = (real(zone, c_double) - 1.0_c_double)*6.0_c_double - 180.0_c_double + 3.0_c_double

    phi = real(lat_deg, c_double)*pi/180.0_c_double
    lam = (real(lon_deg, c_double) - lon0)*pi/180.0_c_double

    rn = a/sqrt(1.0_c_double - e2*sin(phi)**2)
    t  = tan(phi)**2
    cc = ep2*cos(phi)**2
    am = lam*cos(phi)

    m = a*( (1.0_c_double - e2/4.0_c_double - 3.0_c_double*e2**2/64.0_c_double - 5.0_c_double*e2**3/256.0_c_double)*phi &
           - (3.0_c_double*e2/8.0_c_double + 3.0_c_double*e2**2/32.0_c_double + 45.0_c_double*e2**3/1024.0_c_double)*sin(2.0_c_double*phi) &
           + (15.0_c_double*e2**2/256.0_c_double + 45.0_c_double*e2**3/1024.0_c_double)*sin(4.0_c_double*phi) &
           - (35.0_c_double*e2**3/3072.0_c_double)*sin(6.0_c_double*phi) )

    easting = k0*rn*( am + (1.0_c_double-t+cc)*am**3/6.0_c_double &
            + (5.0_c_double-18.0_c_double*t+t**2+72.0_c_double*cc-58.0_c_double*ep2)*am**5/120.0_c_double ) &
            + 500000.0_c_double
    northing = k0*( m + rn*tan(phi)*( am**2/2.0_c_double &
             + (5.0_c_double-t+9.0_c_double*cc+4.0_c_double*cc**2)*am**4/24.0_c_double &
             + (61.0_c_double-58.0_c_double*t+t**2+600.0_c_double*cc-330.0_c_double*ep2)*am**6/720.0_c_double ) )
    if (.not. north) northing = northing + 10000000.0_c_double
  end subroutine lonlat_to_utm

  ! Weight-file format: one  "iHRU  cell  weight"  triple per line (whitespace
  ! separated); blank lines and lines beginning with "#" are ignored.  A given HRU
  ! may span any number of lines.  cell is the row-major horizontal index
  ! (irow-1)*ncol + icol.  Weights are normalised per HRU, so they need only be
  ! relative (e.g. put 1.0 on every line to spread an HRU evenly over its cells).
  subroutine mf6_read_map_file(this, err, message)
    class(mf6_coupler_type), intent(inout) :: this
    integer,                 intent(out)   :: err
    character(len=*),        intent(out)   :: message
    integer :: fu, rc, i, ih, cel
    integer, allocatable :: cell_count(:), fill(:)
    real    :: wgt
    real(c_double) :: wsum
    character(len=256) :: line

    err = 0; message = ''

    open(action='read', file=trim(this%map_file), iostat=rc, newunit=fu)
    if (rc /= 0) then
      message = 'cannot open map_file '//trim(this%map_file); err = 20; return
    end if

    allocate(cell_count(this%nHRU)); cell_count = 0

    ! first pass: count triples per HRU
    do
      read(fu,'(a)',iostat=rc) line
      if (rc /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      read(line,*,iostat=rc) ih, cel, wgt
      if (rc /= 0) then
        close(fu)
        message = 'bad line in map_file: '//trim(line); err = 20; return
      end if
      if (ih < 1 .or. ih > this%nHRU) cycle
      cell_count(ih) = cell_count(ih) + 1
    end do
    rewind(fu)

    allocate(this%map_ptr(this%nHRU+1))
    this%map_ptr(1) = 1
    do i = 1, this%nHRU
      this%map_ptr(i+1) = this%map_ptr(i) + cell_count(i)
    end do
    allocate(this%map_cell(this%map_ptr(this%nHRU+1)-1), this%map_wgt(this%map_ptr(this%nHRU+1)-1))
    allocate(fill(this%nHRU)); fill = this%map_ptr(1:this%nHRU)

    ! second pass: fill
    do
      read(fu,'(a)',iostat=rc) line
      if (rc /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      read(line,*,iostat=rc) ih, cel, wgt
      if (rc /= 0 .or. ih < 1 .or. ih > this%nHRU) cycle
      this%map_cell(fill(ih)) = cel
      this%map_wgt(fill(ih))  = wgt
      fill(ih) = fill(ih) + 1
    end do
    close(fu)

    ! normalise weights per HRU (they are used as relative weights)
    do i = 1, this%nHRU
      wsum = sum(real(this%map_wgt(this%map_ptr(i):this%map_ptr(i+1)-1), c_double))
      if (wsum > 0.0_c_double) &
        this%map_wgt(this%map_ptr(i):this%map_ptr(i+1)-1) = &
          real(real(this%map_wgt(this%map_ptr(i):this%map_ptr(i+1)-1), c_double) / wsum)
    end do
  end subroutine mf6_read_map_file

  ! ==================================================================================
  ! small helpers
  ! ==================================================================================
  subroutine mf6_ptr_double(addr, fptr, n, err, message)
    character(len=*), intent(in) :: addr
    real(c_double), pointer, intent(out) :: fptr(:)
    integer, intent(in) :: n                     ! expected element count
    integer, intent(inout) :: err
    character(len=*), intent(inout) :: message
    type(c_ptr) :: cptr
    integer :: s
    if (err /= 0) return
    s = mf6_get_value_ptr_double(cstr(addr), cptr)
    if (s /= BMI_OK .or. .not. c_associated(cptr)) then
      message = 'MODFLOW variable not found: '//trim(addr); err = 20; return
    end if
    call c_f_pointer(cptr, fptr, [n])
  end subroutine mf6_ptr_double

  subroutine mf6_ptr_scalar(addr, fptr, err, message)
    character(len=*), intent(in) :: addr
    real(c_double), pointer, intent(out) :: fptr
    integer, intent(inout) :: err
    character(len=*), intent(inout) :: message
    type(c_ptr) :: cptr
    integer :: s
    if (err /= 0) return
    s = mf6_get_value_ptr_double(cstr(addr), cptr)
    if (s /= BMI_OK .or. .not. c_associated(cptr)) then
      message = 'MODFLOW scalar not found: '//trim(addr); err = 20; return
    end if
    call c_f_pointer(cptr, fptr)
  end subroutine mf6_ptr_scalar

  subroutine mf6_ptr_int_arr(addr, fptr, n, err, message)
    character(len=*), intent(in) :: addr
    integer(c_int), pointer, intent(out) :: fptr(:)
    integer, intent(in) :: n
    integer, intent(inout) :: err
    character(len=*), intent(inout) :: message
    type(c_ptr) :: cptr
    integer :: s
    if (err /= 0) return
    s = mf6_get_value_ptr_int(cstr(addr), cptr)
    if (s /= BMI_OK .or. .not. c_associated(cptr)) then
      message = 'MODFLOW variable not found: '//trim(addr); err = 20; return
    end if
    call c_f_pointer(cptr, fptr, [n])
  end subroutine mf6_ptr_int_arr

  ! Soft variants of mf6_ptr_double / mf6_ptr_int_arr: return .false. instead of
  ! failing when the MODFLOW variable is absent (used for optional feedback vars).
  ! The pointer is sized from the variable's own nbytes/itemsize.
  logical function mf6_try_ptr_double(addr, fptr) result(ok)
    character(len=*), intent(in) :: addr
    real(c_double), pointer, intent(out) :: fptr(:)
    type(c_ptr) :: cptr
    integer(c_int) :: nbytes, isize
    ok = .false.
    if (mf6_get_value_ptr_double(cstr(addr), cptr) /= BMI_OK) return
    if (.not. c_associated(cptr)) return
    if (mf6_get_var_nbytes(cstr(addr), nbytes) /= BMI_OK) return
    if (mf6_get_var_itemsize(cstr(addr), isize) /= BMI_OK .or. isize <= 0) return
    call c_f_pointer(cptr, fptr, [nbytes/isize])
    ok = .true.
  end function mf6_try_ptr_double

  logical function mf6_try_ptr_int(addr, fptr) result(ok)
    character(len=*), intent(in) :: addr
    integer(c_int), pointer, intent(out) :: fptr(:)
    type(c_ptr) :: cptr
    integer(c_int) :: nbytes, isize
    ok = .false.
    if (mf6_get_value_ptr_int(cstr(addr), cptr) /= BMI_OK) return
    if (.not. c_associated(cptr)) return
    if (mf6_get_var_nbytes(cstr(addr), nbytes) /= BMI_OK) return
    if (mf6_get_var_itemsize(cstr(addr), isize) /= BMI_OK .or. isize <= 0) return
    call c_f_pointer(cptr, fptr, [nbytes/isize])
    ok = .true.
  end function mf6_try_ptr_int

  ! number of elements of a MODFLOW variable = nbytes / itemsize
  integer function mf6_var_count(addr, err, message) result(n)
    character(len=*), intent(in) :: addr
    integer, intent(inout) :: err
    character(len=*), intent(inout) :: message
    integer(c_int) :: nbytes, isize
    integer :: s
    n = 0
    if (err /= 0) return
    s = mf6_get_var_nbytes(cstr(addr), nbytes)
    if (s /= BMI_OK) then
      message = 'MODFLOW variable not found: '//trim(addr); err = 20; return
    end if
    s = mf6_get_var_itemsize(cstr(addr), isize)
    if (s /= BMI_OK .or. isize <= 0) then
      message = 'bad itemsize for MODFLOW variable: '//trim(addr); err = 20; return
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

end module mf6_coupling
