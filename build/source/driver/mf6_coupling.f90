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
  !                                    !    ('' => skip the baseflow feedback).  Shorthand for a
  !                                    !    one-entry bnd_package_names table with role 'baseflow';
  !                                    !    ignored when bnd_package_names is given.
  !     bnd_package_names  = 'CHD','DRN'            ! up to MAXBND head-dependent boundary packages
  !     bnd_package_roles  = 'baseflow','surface_discharge'   ! what SUMMA does with each one's flow:
  !                                    !    'baseflow'          -> scalarAquiferBaseflow, reaches routing
  !                                    !    'surface_discharge' -> added to SUMMA's surface runoff
  !                                    !    'gw_et'             -> groundwater ET actually taken
  !                                    !    Packages sharing a role are summed.
  !     evt_package_name   = ''        ! EVT package (READASARRAYS) whose RATE array is overwritten
  !                                    !    each step with SUMMA's aquifer transpiration demand
  !                                    !    ('' => groundwater ET is not driven from SUMMA)
  !     head_restart_read  = ''        ! read the initial head field from this file instead of IC/STRT
  !     head_restart_write = ''        ! write the final head field here, for a later restart.
  !                                    !    Together these are the coupled restart: run the spin-up
  !                                    !    with _write set, then every later run with _read set, so
  !                                    !    the aquifer and SUMMA's soil column start equilibrated to
  !                                    !    the same thing.  Relative paths resolve against the PROCESS
  !                                    !    working directory, not run_dir, so every calibration rank
  !                                    !    reads the one file the shared spin-up wrote.
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
  ! Boundary-package roles.
  !
  ! What SUMMA does with a returned flux is a property of the ROLE, not of the MODFLOW
  ! package that supplied it.  That indirection is deliberate: the upstream Sagehen model
  ! (the MODFLOW 6 expression of the GSFLOW design) returns water to the land surface from
  ! DRN and groundwater ET from UZF, whereas this coupler starts with DRN and EVT.  Keeping
  ! the role separate from the package means swapping in UZF later is a change here and not
  ! a change in SUMMA.
  !
  !   ROLE_BASEFLOW      -> scalarAquiferBaseflow: reaches routing (basin__AquiferBaseflow)
  !   ROLE_SURFACE_DISCH -> groundwater discharge at land surface: added to surface runoff
  !   ROLE_GW_ET         -> groundwater evapotranspiration actually taken by MODFLOW
  ! ------------------------------------------------------------------------------------
  integer, parameter :: ROLE_BASEFLOW      = 1
  integer, parameter :: ROLE_SURFACE_DISCH = 2
  integer, parameter :: ROLE_GW_ET         = 3
  integer, parameter :: MAXBND             = 8     ! max boundary packages in the &coupler table

  character(len=20), parameter :: ROLE_NAME(3) = [ character(len=20) :: &
      'baseflow', 'surface_discharge', 'gw_et' ]

  ! One head-dependent boundary package whose simulated flow feeds back to SUMMA.
  type :: mf6_bnd_type
    character(len=256)      :: name  = ''
    integer                 :: role  = 0
    logical                 :: found = .false.
    real(c_double), pointer :: simvals(:)  => null()   ! <MODEL>/<PKG>/SIMVALS  (m3 s-1, + into aquifer)
    integer(c_int), pointer :: nodelist(:) => null()   ! <MODEL>/<PKG>/NODELIST (REDUCED node numbers)
  end type mf6_bnd_type

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
    character(len=256)  :: bflow_package_name = 'CHD'   ! back-compatible alias for a single role=baseflow entry
    character(len=256)  :: evt_package_name   = ''      ! EVT package driven with SUMMA's aquifer transpiration demand
    character(len=256)  :: map_file           = ''
    character(len=1024) :: head_restart_read  = ''    ! read the initial head field from here, overriding IC/STRT
    character(len=1024) :: head_restart_write = ''    ! write the final head field here, for a later restart
    integer             :: mf6_epsg           = 0
    character(len=1024) :: run_dir            = '.'      ! directory holding mfsim.nam ('.' = process cwd)
    character(len=1024) :: saved_dir          = ''       ! cwd to return to, while inside run_dir

    ! ---- what the driver needs to know; set by mf6_init ----
    logical, public :: feedback     = .true.   ! .false. => one-way (SUMMA drainage -> MODFLOW only)
    logical, public :: have_sy      = .false.  ! STO/SY found, so aquifer storage can be fed back
    logical, public :: have_bflow   = .false.  ! a role=baseflow package was found
    logical, public :: have_surfdis = .false.  ! a role=surface_discharge package was found
    logical, public :: have_gwet    = .false.  ! a role=gw_et package was found
    logical, public :: have_evt     = .false.  ! an EVT package is available to drive with SUMMA demand

    ! ---- reprojection state for the built-in nearest-cell map (nHRU>1 only; see nearest_hru) ----
    integer :: utm_zone = 0
    logical :: utm_north = .true., use_utm = .false.

    ! ---- SUMMA-side geometry, kept because gather_head_to_hru needs it every step ----
    integer                       :: nHRU = 0
    double precision, allocatable :: hru_x(:), hru_y(:), hru_z(:)  ! HRU centroid lon/lat and surface elevation
    double precision, allocatable :: hru_area(:)                   ! per-HRU plan area from attributes.nc (m2)
    double precision, allocatable :: soil_thk(:)                   ! per-HRU SUMMA soil-column thickness (m)
    real,             allocatable :: sy_hru(:)                     ! per-HRU MODFLOW specific yield (-), map-weighted

    integer :: nss_steps = 0                ! leading steady-state MODFLOW steps run before the coupled period

    ! ---- coupled budget diagnostic (see mf6_budget_report) ----
    logical          :: budget = .false.    ! report the per-step coupled budget
    double precision :: bud_sent = 0.d0     ! cumulative volume SUMMA sent as recharge (m3)
    double precision :: bud_taken = 0.d0    ! cumulative volume MODFLOW's RCH array received (m3)
    double precision :: bud_back(3) = 0.d0  ! cumulative volume returned, by role (m3)
    double precision :: bud_etdem = 0.d0    ! cumulative aquifer-transpiration demand sent (m3)

    ! ---- MODFLOW 6 side ----
    real(c_double), pointer :: mf6_head(:) => null()      ! GWF dependent variable  <MODEL>/X       (REDUCED nodes)
    real(c_double), pointer :: mf6_rch(:)  => null()      ! RCH recharge array      <MODEL>/<RCH>/RECHARGE  (user cells, layer 1)
    integer(c_int), pointer :: mf6_mshape(:) => null()    ! DIS grid shape          <MODEL>/DIS/MSHAPE (nlay,nrow,ncol)
    integer(c_int), pointer :: mf6_nodered(:) => null()   ! DIS full->reduced node map <MODEL>/DIS/NODEREDUCED (present only when the grid is reduced)
    real(c_double), pointer :: mf6_sy(:) => null()        ! STO specific yield      <MODEL>/STO/SY  (REDUCED nodes)
    real(c_double), pointer :: mf6_evt(:) => null()       ! EVT max rate array      <MODEL>/<EVT>/RATE (user cells, READASARRAYS)
    logical                 :: grid_reduced = .false.

    ! ---- head-dependent boundary packages fed back to SUMMA, by role ----
    integer              :: nbnd = 0
    type(mf6_bnd_type)   :: bnd(MAXBND)

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
    procedure, private :: check_hru_area        => mf6_check_hru_area
    procedure, private :: is_steady_state       => mf6_is_steady_state
    procedure, private :: head_restart_in       => mf6_head_restart_read
    procedure, private :: head_restart_out      => mf6_head_restart_write
    procedure, private :: resolve_path          => mf6_resolve_path
    procedure, private :: budget_accumulate     => mf6_budget_accumulate
    procedure, private :: budget_report         => mf6_budget_report
    procedure, private :: build_map             => mf6_build_map
    procedure, private :: build_nearest_cell_map=> mf6_build_nearest_cell_map
    procedure, private :: read_map_file         => mf6_read_map_file
    procedure, private :: build_sy_hru          => mf6_build_sy_hru
    procedure, private :: scatter_drainage_to_rch => mf6_scatter_drainage_to_rch
    procedure, private :: scatter_hru_to_array  => mf6_scatter_hru_to_array
    procedure, private :: gather_head_to_hru    => mf6_gather_head_to_hru
    procedure, private :: gather_aquifer_to_hru => mf6_gather_aquifer_to_hru
    procedure, private :: gather_boundary_to_hru=> mf6_gather_boundary_to_hru
    procedure, private :: gather_role_to_hru    => mf6_gather_role_to_hru
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
                      nSummaSteps, summa_data_step, err, message, hru_area, &
                      restart_read, restart_write)
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
    double precision, optional, intent(in) :: hru_area(:)       ! HRU plan area (m2); enables the area and budget checks
    ! head-restart paths set by the caller; these override whatever the &coupler namelist says, so a
    ! driver can run the same config as a spin-up (write) and then as an evaluation run (read)
    character(len=*), optional, intent(in) :: restart_read, restart_write
    integer        :: nred_len, nvar, ib
    integer        :: istat
    real(c_double) :: tend_mf6

    err = 0; message = ''

    ! -- configuration --
    this%run_dir = run_dir
    call this%read_config(config_file, err, message)
    if (err /= 0) return
    if (present(restart_read))  this%head_restart_read  = restart_read
    if (present(restart_write)) this%head_restart_write = restart_write

    ! -- SUMMA-side geometry (kept: gather_head_to_hru needs hru_z and soil_thk every step) --
    this%nHRU = nHRU
    allocate(this%hru_x(nHRU), this%hru_y(nHRU), this%hru_z(nHRU), this%soil_thk(nHRU))
    this%hru_x    = hru_x(1:nHRU)
    this%hru_y    = hru_y(1:nHRU)
    this%hru_z    = hru_z(1:nHRU)
    this%soil_thk = soil_thk(1:nHRU)
    if (present(hru_area)) then
      allocate(this%hru_area(nHRU)); this%hru_area = hru_area(1:nHRU)
    end if

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

      ! every boundary package in the table: SIMVALS + NODELIST, or a warning and that role off
      do ib = 1, this%nbnd
        this%bnd(ib)%found = &
          mf6_try_ptr_double(trim(this%mf6_model_name)//'/'//trim(this%bnd(ib)%name)//'/SIMVALS',  this%bnd(ib)%simvals) .and. &
          mf6_try_ptr_int   (trim(this%mf6_model_name)//'/'//trim(this%bnd(ib)%name)//'/NODELIST', this%bnd(ib)%nodelist)
        if (.not. this%bnd(ib)%found) write(*,'(a)') 'summa_modflow6: WARNING - boundary package "'// &
          trim(this%bnd(ib)%name)//'" (role '//trim(ROLE_NAME(this%bnd(ib)%role))// &
          ') not found; that feedback is disabled'
      end do
      this%have_bflow   = any(this%bnd(1:this%nbnd)%found .and. this%bnd(1:this%nbnd)%role == ROLE_BASEFLOW)
      this%have_surfdis = any(this%bnd(1:this%nbnd)%found .and. this%bnd(1:this%nbnd)%role == ROLE_SURFACE_DISCH)
      this%have_gwet    = any(this%bnd(1:this%nbnd)%found .and. this%bnd(1:this%nbnd)%role == ROLE_GW_ET)
    end if

    ! -- EVT: the one package SUMMA writes into rather than reads from (Phase 3).  READASARRAYS,
    !    so RATE is one value per horizontal cell, exactly like RCH's RECHARGE.
    if (len_trim(this%evt_package_name) > 0) then
      this%have_evt = mf6_try_ptr_double(trim(this%mf6_model_name)//'/'// &
                                         trim(this%evt_package_name)//'/RATE', this%mf6_evt)
      if (this%have_evt) then
        if (size(this%mf6_evt) /= this%nrow*this%ncol) then
          write(message,'(a,i0,a,i0,a)') 'EVT package '//trim(this%evt_package_name)// &
                ' must be declared READASARRAYS: its RATE array has ', size(this%mf6_evt), &
                ' entries, not one per horizontal cell (', this%nrow*this%ncol, ')'
          err = 20; call this%leave_run_dir(); return
        end if
      else
        write(*,'(a)') 'summa_modflow6: WARNING - EVT package "'//trim(this%evt_package_name)// &
          '" not found; groundwater evapotranspiration is not driven from SUMMA'
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

    ! -- coupled restart: overwrite STRT-derived heads with a previously saved field --
    call this%head_restart_in(err, message)
    if (err /= 0) then; call this%leave_run_dir(); return; end if

    ! -- HRU area vs mapped-cell area: the silent failure mode behind Equation (2) --
    call this%check_hru_area

    ! -- per-HRU specific yield (map-weighted), fixed for the run --
    if (this%feedback .and. this%have_sy) call this%build_sy_hru

    ! the coupled budget needs HRU areas to convert per-HRU fluxes to volumes
    this%budget = allocated(this%hru_area)
    this%bud_sent = 0.d0; this%bud_taken = 0.d0; this%bud_back = 0.d0; this%bud_etdem = 0.d0

    call this%leave_run_dir()
  end subroutine mf6_init

  ! ==================================================================================
  ! One coupling step: steps 3 to 5 of the exchange described at the top of this module.
  ! The driver has already pushed the previous step's feedback into SUMMA and advanced it.
  !
  ! head_hru/stor_hru/bflow_hru are the driver's own lagged feedback arrays: they are read
  ! back here for the next step, and are left untouched when feedback is off.
  ! ==================================================================================
  subroutine mf6_step(this, istep, summa_data_step, drain_hru, head_hru, stor_hru, bflow_hru, err, message, &
                      surfdis_hru, gwet_demand_hru, gwet_hru)
    class(mf6_coupler_type), intent(inout) :: this
    integer,                 intent(in)    :: istep             ! 1-based coupling step index
    double precision,        intent(in)    :: summa_data_step   ! length of a SUMMA data step (s)
    real,                    intent(in)    :: drain_hru(:)      ! per-HRU soil drainage (m s-1)
    real,                    intent(inout) :: head_hru(:)       ! per-HRU prescribed head (m, matric head at soil base)
    real,                    intent(inout) :: stor_hru(:)       ! per-HRU relative aquifer storage (m)
    real,                    intent(inout) :: bflow_hru(:)      ! per-HRU aquifer baseflow (m s-1, + = out of aquifer)
    integer,                 intent(out)   :: err
    character(len=*),        intent(out)   :: message
    ! optional extra exchange channels; a driver that does not pass them behaves exactly as before
    real, optional,          intent(inout) :: surfdis_hru(:)      ! per-HRU groundwater discharge at land surface (m s-1, + = out)
    real, optional,          intent(in)    :: gwet_demand_hru(:)  ! per-HRU aquifer transpiration DEMAND from SUMMA (m s-1, + = out)
    real, optional,          intent(inout) :: gwet_hru(:)         ! per-HRU groundwater ET ACTUALLY taken by MODFLOW (m s-1, + = out)
    integer        :: istat
    real(c_double) :: dt_mf6

    err = 0; message = ''

    call this%enter_run_dir(err, message)
    if (err /= 0) return

    ! 3. prepare the MODFLOW time step FIRST.
    !
    ! The scatter must come after this, not before.  prepare_time_step runs each package's *_rp,
    ! and bnd_rp re-reads the package's PERIOD block whenever a new stress period starts
    ! (BoundaryPackage.f90: "if (this%ionper < kper)").  For RCH READASARRAYS that reload
    ! overwrites the whole RECHARGE array with the value in the input file.  Scattering before
    ! prepare therefore threw SUMMA's drainage away on the first step of every stress period, and
    ! MODFLOW solved that step with the model builder's placeholder constant instead - silently,
    ! because it is one step in seventy-two.  Preparing first and scattering after is the correct
    ! order: the reload happens, then SUMMA's values overwrite it, then the solve sees them.
    ! Leading steady-state stress periods are run out here, before the coupled loop proper.
    !
    ! The upstream Sagehen example equilibrates the water table with a steady-state first stress
    ! period and only then goes transient; that was simplified out when this test model was built
    ! (section 6.1), and section 8.5 wants it back, because otherwise the coupled model starts from
    ! whatever STRT happens to contain.  MODFLOW reports the current period's steady-state flag in
    ! <MODEL>/ISS, but only once the step is prepared, so the shape here is prepare-then-test:
    ! a steady-state step is solved and skipped, and the loop exits on the first transient step,
    ! which is the one the coupling actually drives.
    !
    ! A steady-state step is deliberately NOT given SUMMA's drainage.  It is solved with whatever
    ! the RCH package's own PERIOD block supplies - a climatological mean the model builder
    ! provides - because an instantaneous first-hour drainage rate is not a sensible thing to
    ! equilibrate an aquifer against.  That is also why the scatter sits after this loop.
    !
    ! Models with no steady-state period read ISS = 0 on the first prepare and fall straight
    ! through, so nothing changes for them.
    do
      istat = mf6_prepare_time_step(real(summa_data_step, c_double))
      if (istep /= 1) exit                      ! only leading periods; later steps are transient
      if (.not. this%is_steady_state()) exit
      istat = mf6_do_time_step()
      istat = mf6_finalize_time_step()
      this%nss_steps = this%nss_steps + 1
    end do
    if (istep == 1 .and. this%nss_steps > 0) then
      write(*,'(a,i0,a)') 'summa_modflow6: ran ', this%nss_steps, &
            ' steady-state MODFLOW step(s) before the coupled period, to equilibrate the water table'
    end if

    if (istep == 1) then
      ! delt is now set, and this is the first TRANSIENT step: verify one MODFLOW step == one
      ! SUMMA data step.  Deliberately checked here rather than on a steady-state step, whose
      ! perlen has nothing to do with the coupling interval.
      istat = mf6_get_time_step(dt_mf6)
      if (abs(dt_mf6 - real(summa_data_step, c_double)) > 1.0e-6_c_double*real(summa_data_step, c_double)) then
        write(message,'(a,g0,a,g0,a)') 'MODFLOW time step (', dt_mf6, &
              ') must equal the SUMMA data step (', summa_data_step, ')'
        err = 20
        call this%leave_run_dir(); return
      end if
    end if

    ! 4. SUMMA soil drainage -> MODFLOW 6 RCH recharge array
    call this%scatter_drainage_to_rch(drain_hru)

    ! 4b. SUMMA aquifer transpiration demand -> MODFLOW 6 EVT rate array.  MODFLOW applies its
    !     own water-table-depth limiting, so what it actually takes comes back at step 5.
    if (this%have_evt .and. present(gwet_demand_hru)) &
      call this%scatter_hru_to_array(gwet_demand_hru, this%mf6_evt)

    ! 4c. solve
    istat = mf6_do_time_step()
    istat = mf6_finalize_time_step()

    ! 5. read the new MODFLOW state, aggregate per HRU for the next iteration
    if (this%feedback) then
      call this%gather_head_to_hru(head_hru)
      call this%gather_aquifer_to_hru(head_hru, stor_hru, bflow_hru)
      if (this%have_surfdis .and. present(surfdis_hru)) &
        call this%gather_role_to_hru(ROLE_SURFACE_DISCH, surfdis_hru)
      if (this%have_gwet .and. present(gwet_hru)) &
        call this%gather_role_to_hru(ROLE_GW_ET, gwet_hru)
    end if

    ! 6. accumulate the coupled budget, now that both sides of this step are known
    call this%budget_accumulate(summa_data_step, drain_hru, bflow_hru, surfdis_hru, gwet_hru, gwet_demand_hru)

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
    integer :: istat, ib

    err = 0; message = ''

    ! report before tearing anything down, so a failed finalize still leaves the budget visible
    call this%budget_report

    call this%enter_run_dir(err, message)
    if (err == 0) then
      call this%head_restart_out          ! before finalize: libmf6 is about to free X
      istat = mf6_finalize()
      call this%leave_run_dir()
      if (istat /= BMI_OK) then; message = 'MODFLOW 6 finalize failed'; err = 20; end if
    end if

    ! release everything owned here; the MODFLOW-side pointers belong to libmf6, which has
    ! just freed them, so drop them rather than leave them dangling
    if (allocated(this%hru_x))     deallocate(this%hru_x)
    if (allocated(this%hru_y))     deallocate(this%hru_y)
    if (allocated(this%hru_z))     deallocate(this%hru_z)
    if (allocated(this%hru_area))  deallocate(this%hru_area)
    if (allocated(this%soil_thk))  deallocate(this%soil_thk)
    if (allocated(this%sy_hru))    deallocate(this%sy_hru)
    if (allocated(this%map_ptr))   deallocate(this%map_ptr)
    if (allocated(this%map_cell))  deallocate(this%map_cell)
    if (allocated(this%map_wgt))   deallocate(this%map_wgt)
    if (allocated(this%cell_area)) deallocate(this%cell_area)
    this%mf6_head     => null()
    this%mf6_rch      => null()
    this%mf6_evt      => null()
    this%mf6_mshape   => null()
    this%mf6_nodered  => null()
    this%mf6_sy       => null()
    do ib = 1, MAXBND
      this%bnd(ib)%simvals  => null()
      this%bnd(ib)%nodelist => null()
      this%bnd(ib)%found    = .false.
    end do
    this%cellx        => null()
    this%celly        => null()
    this%xorigin      => null()
    this%yorigin      => null()
    this%angrot       => null()
    this%grid_reduced = .false.
    this%have_sy      = .false.
    this%have_bflow   = .false.
    this%have_surfdis = .false.
    this%have_gwet    = .false.
    this%have_evt     = .false.
    this%budget       = .false.
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
    integer :: fu, rc, ib, ir
    character(len=256) :: mf6_model_name, rch_package_name, bflow_package_name, map_file
    character(len=256) :: evt_package_name
    character(len=1024):: head_restart_read, head_restart_write
    character(len=256) :: bnd_package_names(MAXBND)
    character(len=32)  :: bnd_package_roles(MAXBND)
    integer            :: mf6_epsg
    logical            :: feedback
    namelist /coupler/ mf6_model_name, rch_package_name, bflow_package_name, evt_package_name, &
                       bnd_package_names, bnd_package_roles, map_file, mf6_epsg, feedback, &
                       head_restart_read, head_restart_write

    err = 0; message = ''

    ! namelist defaults
    mf6_model_name     = ''
    rch_package_name   = 'RCHA'
    bflow_package_name = 'CHD'
    evt_package_name   = ''
    bnd_package_names  = ''
    bnd_package_roles  = ''
    head_restart_read  = ''
    head_restart_write = ''
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
    this%evt_package_name   = evt_package_name
    this%head_restart_read  = head_restart_read
    this%head_restart_write = head_restart_write
    this%map_file           = map_file
    this%mf6_epsg           = mf6_epsg
    this%feedback           = feedback

    call to_upper(this%mf6_model_name)
    call to_upper(this%rch_package_name)
    call to_upper(this%bflow_package_name)
    call to_upper(this%evt_package_name)

    ! -- build the boundary-package table --
    ! Two ways in, and the explicit table wins.  If bnd_package_names is left blank the single
    ! bflow_package_name entry is used as role=baseflow, which is exactly what every config
    ! written before this table existed means, so they keep working untouched.
    this%nbnd = 0
    if (len_trim(bnd_package_names(1)) > 0) then
      do ib = 1, MAXBND
        if (len_trim(bnd_package_names(ib)) == 0) cycle
        call to_upper(bnd_package_names(ib))
        call to_lower(bnd_package_roles(ib))
        ir = role_code(bnd_package_roles(ib))
        if (ir == 0) then
          message = 'unknown bnd_package_roles entry "'//trim(bnd_package_roles(ib))// &
                '" for package "'//trim(bnd_package_names(ib))//'" in '//trim(config_file)// &
                ' - expected one of: '//trim(ROLE_NAME(1))//', '//trim(ROLE_NAME(2))//', '//trim(ROLE_NAME(3))
          err = 20; return
        end if
        this%nbnd = this%nbnd + 1
        this%bnd(this%nbnd)%name = bnd_package_names(ib)
        this%bnd(this%nbnd)%role = ir
      end do
    else if (len_trim(this%bflow_package_name) > 0) then
      this%nbnd = 1
      this%bnd(1)%name = this%bflow_package_name
      this%bnd(1)%role = ROLE_BASEFLOW
    end if

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
    call this%scatter_hru_to_array(drain_hru, this%mf6_rch)
  end subroutine mf6_scatter_drainage_to_rch

  ! ==================================================================================
  ! The general version: a per-HRU flux (m s-1) onto a READASARRAYS per-horizontal-cell
  ! MODFLOW array.  RCH's RECHARGE and EVT's RATE are both shaped this way, so both go
  ! through here rather than repeating the weighting.  Cells with no mapped HRU get zero.
  ! ==================================================================================
  subroutine mf6_scatter_hru_to_array(this, flux_hru, arr)
    class(mf6_coupler_type), intent(inout) :: this
    real,                    intent(in)    :: flux_hru(:)
    real(c_double), pointer, intent(inout) :: arr(:)
    real(c_double), allocatable :: num(:), den(:)
    integer :: i, k, c

    if (.not. associated(arr)) return

    allocate(num(this%nrow*this%ncol), den(this%nrow*this%ncol))
    num = 0.0_c_double
    den = 0.0_c_double
    do i = 1, this%nHRU
      do k = this%map_ptr(i), this%map_ptr(i+1) - 1
        c = this%map_cell(k)
        if (c < 1 .or. c > this%nrow*this%ncol) cycle
        num(c) = num(c) + real(this%map_wgt(k), c_double) * this%cell_area(c) * real(flux_hru(i), c_double)
        den(c) = den(c) + real(this%map_wgt(k), c_double) * this%cell_area(c)
      end do
    end do

    ! the array is indexed by horizontal cell (row-major) for a READASARRAYS package
    do c = 1, min(size(arr), this%nrow*this%ncol)
      if (den(c) > 0.0_c_double) then
        arr(c) = num(c) / den(c)
      else
        arr(c) = 0.0_c_double
      end if
    end do
    deallocate(num, den)
  end subroutine mf6_scatter_hru_to_array

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
    integer :: i

    if (this%have_sy) then
      do i = 1, this%nHRU
        stor_hru(i) = this%sy_hru(i) * head_hru(i)
      end do
    end if

    if (this%have_bflow) call this%gather_role_to_hru(ROLE_BASEFLOW, bflow_hru)
  end subroutine mf6_gather_aquifer_to_hru

  ! ==================================================================================
  ! Sum every boundary package carrying the given role into one per-HRU flux.
  !
  ! Several packages may share a role (two drain packages both discharging at land
  ! surface, say), so the contributions add.  Packages that were not found are skipped,
  ! which is what makes a missing optional package a warning at init rather than a crash
  ! here.  head_hru is left alone when no package carries the role, so the driver's
  ! lagged array simply keeps its previous value - the same contract as gather_head_to_hru.
  ! ==================================================================================
  subroutine mf6_gather_role_to_hru(this, role, out_hru)
    class(mf6_coupler_type), intent(inout) :: this
    integer,                 intent(in)    :: role
    real,                    intent(inout) :: out_hru(:)
    integer :: ib
    logical :: first

    first = .true.
    do ib = 1, this%nbnd
      if (this%bnd(ib)%role /= role) cycle
      if (.not. this%bnd(ib)%found) cycle
      call this%gather_boundary_to_hru(ib, out_hru, accumulate = .not. first)
      first = .false.
    end do
  end subroutine mf6_gather_role_to_hru

  ! ==================================================================================
  ! One boundary package's simulated flow  ->  per-HRU flux (m s-1, + = out of aquifer).
  !
  !     out_hru(i) = -( sum_k w_ik * Q_k ) / ( sum_k w_ik * A_k )
  !
  ! MODFLOW reports boundary flow as a volumetric rate (m3 s-1) that is positive INTO the
  ! aquifer, so the sign flips: a drain or a constant head taking water out of the aquifer
  ! becomes a positive outward flux per unit area.  NODELIST holds reduced node numbers,
  ! which differ from full-grid numbering wherever IDOMAIN deactivates cells, hence
  ! to_reduced.  This is the gather that used to be inline in gather_aquifer_to_hru; it is
  ! separate now because roles other than baseflow need exactly the same arithmetic.
  ! ==================================================================================
  subroutine mf6_gather_boundary_to_hru(this, ib, out_hru, accumulate)
    class(mf6_coupler_type), intent(inout) :: this
    integer,                 intent(in)    :: ib           ! index into this%bnd
    real,                    intent(inout) :: out_hru(:)
    logical, optional,       intent(in)    :: accumulate   ! .true. => add to out_hru rather than overwrite
    integer :: i, k, c, kb, nr
    real(c_double) :: fnum, aden
    real           :: flux
    logical        :: add

    add = .false.
    if (present(accumulate)) add = accumulate

    do i = 1, this%nHRU
      fnum = 0.0_c_double     ! signed boundary flow over the HRU's mapped cells (m3 s-1, + into aquifer)
      aden = 0.0_c_double     ! matching mapped-cell plan area (m2)
      do k = this%map_ptr(i), this%map_ptr(i+1) - 1
        c = this%map_cell(k)
        if (c < 1 .or. c > this%nrow*this%ncol) cycle
        nr = this%to_reduced(c)                  ! boundary NODELIST holds reduced node numbers
        if (nr > 0) then
          do kb = 1, size(this%bnd(ib)%simvals)
            if (this%bnd(ib)%nodelist(kb) == nr) &
              fnum = fnum + real(this%map_wgt(k), c_double) * this%bnd(ib)%simvals(kb)
          end do
        end if
        aden = aden + real(this%map_wgt(k), c_double) * this%cell_area(c)
      end do
      if (aden > 0.0_c_double) then
        flux = real(-fnum / aden)
        if (add) then
          out_hru(i) = out_hru(i) + flux
        else
          out_hru(i) = flux
        end if
      end if
    end do
  end subroutine mf6_gather_boundary_to_hru

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

  ! Is the stress period MODFLOW has just prepared a steady-state one?
  !
  ! <MODEL>/ISS is the GWF model's steady-state flag (gwf.f90 allocates it; gwf-sto.f90 points at
  ! it and sets it from each period's STO PERIOD block).  It is only meaningful once a step has
  ! been prepared, which is why callers prepare first and ask afterwards.  A model with no STO
  ! package, or one where the flag cannot be read, is treated as transient - the safe default,
  ! since that is exactly how the coupler behaved before steady-state periods were supported.
  logical function mf6_is_steady_state(this) result(ss)
    class(mf6_coupler_type), intent(in) :: this
    integer(c_int), pointer :: iss(:) => null()
    ss = .false.
    if (.not. mf6_try_ptr_int(trim(this%mf6_model_name)//'/ISS', iss)) return
    if (size(iss) < 1) return
    ss = (iss(1) /= 0)
  end function mf6_is_steady_state

  ! ==================================================================================
  ! Coupled restart of the MODFLOW head field.
  !
  ! Without this, MODFLOW always starts from whatever IC/STRT contains.  In calibration that
  ! is actively wrong: the one-year cold-start spin-up leaves SUMMA equilibrated and shared
  ! across parameter samples, but every sample restarts MODFLOW from raw STRT, so the soil
  ! column and the aquifer are equilibrated to different things (section 8.5).  Writing the
  ! head field at the end of the spin-up and reading it back at the start of every sample
  ! makes the two consistent, at one small file per rank.
  !
  ! Writing X straight into the memory manager after initialize() is enough: MODFLOW's
  ! gwf_ad copies x into xold at the start of every non-retry time step (gwf.f90), so the
  ! value written here becomes both the initial head AND the previous-step head, which is
  ! what a restart means.  There is no need to write XOLD separately.
  !
  ! The file carries the grid shape and the reduced node count so that restarting into a
  ! different model fails loudly instead of scrambling the head field.
  !
  ! Paths are resolved against the process working directory, not the MODFLOW run directory:
  ! a calibration gives every rank its own run_dir, and all of them must read the ONE head
  ! file the shared spin-up wrote.
  ! ==================================================================================
  subroutine mf6_head_restart_read(this, err, message)
    class(mf6_coupler_type), intent(inout) :: this
    integer,                 intent(out)   :: err
    character(len=*),        intent(out)   :: message
    character(len=8)  :: magic
    integer           :: fu, rc, nlay, nrow, ncol, n
    logical           :: there
    character(len=1024) :: path

    err = 0; message = ''
    if (len_trim(this%head_restart_read) == 0) return

    path = this%resolve_path(this%head_restart_read)
    inquire(file=trim(path), exist=there)
    if (.not. there) then
      write(*,'(a)') 'summa_modflow6: WARNING - head_restart_read file "'//trim(path)// &
        '" does not exist; starting MODFLOW from IC/STRT instead'
      return
    end if

    open(newunit=fu, file=trim(path), form='unformatted', access='stream', &
         action='read', status='old', iostat=rc)
    if (rc /= 0) then
      message = 'cannot open head restart file '//trim(path); err = 20; return
    end if
    read(fu, iostat=rc) magic, nlay, nrow, ncol, n
    if (rc /= 0 .or. magic /= 'SUMMF6HD') then
      close(fu); message = trim(path)//' is not a SUMMA-MODFLOW head restart file'; err = 20; return
    end if
    if (nlay /= this%nlay .or. nrow /= this%nrow .or. ncol /= this%ncol .or. n /= size(this%mf6_head)) then
      close(fu)
      write(message,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)') &
        'head restart file '//trim(path)//' was written for a different model: ', &
        nlay, 'x', nrow, 'x', ncol, ' with ', n, ' active nodes, but this model is ', &
        this%nlay, 'x', this%nrow, 'x', this%ncol, ' with ', size(this%mf6_head), ' active nodes'
      err = 20; return
    end if
    read(fu, iostat=rc) this%mf6_head
    close(fu)
    if (rc /= 0) then
      message = 'error reading heads from '//trim(path); err = 20; return
    end if
    write(*,'(a,i0,a)') 'summa_modflow6: restarted MODFLOW from ', n, ' heads in '//trim(path)
  end subroutine mf6_head_restart_read

  subroutine mf6_head_restart_write(this)
    class(mf6_coupler_type), intent(inout) :: this
    integer :: fu, rc
    character(len=1024) :: path

    if (len_trim(this%head_restart_write) == 0) return
    if (.not. associated(this%mf6_head)) return

    path = this%resolve_path(this%head_restart_write)
    open(newunit=fu, file=trim(path), form='unformatted', access='stream', &
         action='write', status='replace', iostat=rc)
    if (rc /= 0) then
      write(*,'(a)') 'summa_modflow6: WARNING - cannot write head restart file '//trim(path)
      return
    end if
    write(fu) 'SUMMF6HD', this%nlay, this%nrow, this%ncol, size(this%mf6_head)
    write(fu) this%mf6_head
    close(fu)
    write(*,'(a,i0,a)') 'summa_modflow6: wrote ', size(this%mf6_head), ' MODFLOW heads to '//trim(path)
  end subroutine mf6_head_restart_write

  ! A relative path is taken relative to the process working directory, which while inside the
  ! run directory is remembered in saved_dir.  run_dir = '.' means no directory change happened
  ! at all, so the path is already correct as written.
  function mf6_resolve_path(this, path) result(full)
    class(mf6_coupler_type), intent(in) :: this
    character(len=*),        intent(in) :: path
    character(len=1024) :: full
    if (path(1:1) == '/' .or. len_trim(this%saved_dir) == 0) then
      full = path
    else
      full = trim(this%saved_dir)//'/'//trim(path)
    end if
  end function mf6_resolve_path

  ! ==================================================================================
  ! Each HRU's area in attributes.nc must equal the summed plan area of the MODFLOW cells
  ! it maps to, or recharge VOLUME is not conserved across the interface.
  !
  ! scatter_hru_to_array writes a weight-weighted mean of HRU drainage RATES, because the
  ! cell area is common to every contribution to a given cell and cancels.  MODFLOW then
  ! multiplies that rate by ITS cell area, so the volume that arrives is
  !     sum_c A_c * q      instead of      A_HRU * q,
  ! and the two agree only when the areas do.  The same applies in reverse to the baseflow
  ! gather, which divides by mapped-cell area and is then applied per unit HRU area.
  !
  ! This is a warning rather than a hard stop: a nearest-cell map over a real basin will
  ! rarely match exactly, and the honest response is to say by how much rather than to
  ! refuse to run.  The elevation check above stops because a bad elevation makes SUMMA
  ! fail to converge; a bad area silently mis-scales a flux, which is why it is reported
  ! here AND accumulated by the budget diagnostic every step.
  ! ==================================================================================
  subroutine mf6_check_hru_area(this)
    class(mf6_coupler_type), intent(inout) :: this
    integer :: i, k, c, nbad
    real(c_double) :: amap, rel, worst
    real(c_double), allocatable :: wcell(:)
    integer :: iworst

    if (.not. allocated(this%hru_area)) return

    ! Total weight landing on each cell, so a cell shared between HRUs is apportioned rather than
    ! counted once per owner.  This is the same denominator the scatter uses: cell c receives the
    ! weighted MEAN rate sum_i w_ic q_i / sum_i w_ic, so HRU i's share of that cell's volume is
    ! w_ic / sum_j w_jc.  The 4-HRU Sagehen map has two HRUs deliberately sharing 1700 cells, and
    ! without this apportioning each would appear to own the whole of them.
    allocate(wcell(this%nrow*this%ncol)); wcell = 0.0_c_double
    do i = 1, this%nHRU
      do k = this%map_ptr(i), this%map_ptr(i+1) - 1
        c = this%map_cell(k)
        if (c < 1 .or. c > this%nrow*this%ncol) cycle
        wcell(c) = wcell(c) + real(this%map_wgt(k), c_double)
      end do
    end do

    nbad = 0; worst = 0.0_c_double; iworst = 0
    do i = 1, this%nHRU
      ! effective mapped area: sum_k (w_ik / sum_j w_jk) * A_k.  Recharge VOLUME is conserved for
      ! HRU i exactly when this equals its area in attributes.nc.
      amap = 0.0_c_double
      do k = this%map_ptr(i), this%map_ptr(i+1) - 1
        c = this%map_cell(k)
        if (c < 1 .or. c > this%nrow*this%ncol) cycle
        if (wcell(c) <= 0.0_c_double) cycle
        amap = amap + (real(this%map_wgt(k), c_double) / wcell(c)) * this%cell_area(c)
      end do
      if (this%hru_area(i) <= 0.0d0 .or. amap <= 0.0_c_double) cycle
      rel = abs(amap - this%hru_area(i)) / this%hru_area(i)
      if (rel > 0.01_c_double) then
        nbad = nbad + 1
        if (rel > worst) then; worst = rel; iworst = i; end if
        if (nbad <= 5) write(*,'(a,i0,a,g0,a,g0,a,f0.1,a)') &
          'summa_modflow6: WARNING - HRU ', i, ' area ', this%hru_area(i), &
          ' m2 but its mapped MODFLOW cells total ', amap, ' m2 (', 100.0_c_double*rel, '% out)'
      end if
    end do

    if (nbad > 0) then
      write(*,'(a,i0,a,i0,a,f0.1,a,i0,a)') 'summa_modflow6: WARNING - ', nbad, ' of ', this%nHRU, &
        ' HRUs disagree with their mapped cell area by more than 1% (worst ', &
        100.0_c_double*worst, '% at HRU ', iworst, '). Recharge VOLUME is not conserved for those '// &
        'HRUs: the scatter conserves rate, not volume. Supply a map_file with exact intersection '// &
        'weights, or correct HRUarea in attributes.nc.'
    end if
    deallocate(wcell)
  end subroutine mf6_check_hru_area

  ! ==================================================================================
  ! Per-step coupled budget: what SUMMA sent, what MODFLOW took, what came back.
  !
  ! Section 8.4 of the write-up notes that SUMMA's coupled water balance is an assembly of
  ! quantities from two models at two different times whose closure has never been audited.
  ! This is that audit.  Volumes are accumulated over the run and reported at finalize:
  !
  !   sent   = sum_i A_HRU,i * q_i * dt                (what SUMMA handed over)
  !   taken  = sum_c A_c * R_c * dt                    (what MODFLOW's RCH array actually got)
  !   back   = sum_i A_HRU,i * f_i * dt  per role      (what came back to SUMMA)
  !
  ! sent /= taken is exactly the area mismatch check_hru_area warns about, measured in m3
  ! rather than percent.  The residual sent - taken is NOT the same thing as MODFLOW's own
  ! budget discrepancy, which is dominated by the one-step lag (section 8.2); this isolates
  ! the mapping error from the lag error, which is the point of reporting it separately.
  ! ==================================================================================
  subroutine mf6_budget_accumulate(this, dt, drain_hru, bflow_hru, surfdis_hru, gwet_hru, gwet_demand_hru)
    class(mf6_coupler_type), intent(inout) :: this
    double precision,        intent(in)    :: dt
    real,                    intent(in)    :: drain_hru(:)
    real,                    intent(in)    :: bflow_hru(:)
    real, optional,          intent(in)    :: surfdis_hru(:), gwet_hru(:), gwet_demand_hru(:)
    integer :: i, c

    if (.not. this%budget) return

    ! SUMMA side: per-HRU flux over the HRU's own area
    do i = 1, this%nHRU
      this%bud_sent = this%bud_sent + this%hru_area(i) * dble(drain_hru(i)) * dt
      if (this%have_bflow)   this%bud_back(ROLE_BASEFLOW) = &
        this%bud_back(ROLE_BASEFLOW) + this%hru_area(i) * dble(bflow_hru(i)) * dt
      if (present(surfdis_hru) .and. this%have_surfdis) this%bud_back(ROLE_SURFACE_DISCH) = &
        this%bud_back(ROLE_SURFACE_DISCH) + this%hru_area(i) * dble(surfdis_hru(i)) * dt
      if (present(gwet_hru) .and. this%have_gwet) this%bud_back(ROLE_GW_ET) = &
        this%bud_back(ROLE_GW_ET) + this%hru_area(i) * dble(gwet_hru(i)) * dt
      if (present(gwet_demand_hru)) &
        this%bud_etdem = this%bud_etdem + this%hru_area(i) * dble(gwet_demand_hru(i)) * dt
    end do

    ! MODFLOW side: the RCH array as it now stands, over the cells' own areas
    do c = 1, min(size(this%mf6_rch), this%nrow*this%ncol)
      this%bud_taken = this%bud_taken + this%cell_area(c) * this%mf6_rch(c) * dt
    end do
  end subroutine mf6_budget_accumulate

  subroutine mf6_budget_report(this)
    class(mf6_coupler_type), intent(in) :: this
    double precision :: resid, pct

    if (.not. this%budget) return

    resid = this%bud_sent - this%bud_taken
    pct   = 0.d0
    if (abs(this%bud_sent) > 0.d0) pct = 100.d0 * resid / this%bud_sent

    write(*,'(a)')        'summa_modflow6: coupled water budget over the run (m3)'
    write(*,'(a,g0)')     '  SUMMA drainage sent      : ', this%bud_sent
    write(*,'(a,g0)')     '  MODFLOW recharge received: ', this%bud_taken
    write(*,'(a,g0,a,f0.3,a)') '  mapping residual         : ', resid, '  (', pct, '% of sent)'
    if (this%have_bflow)   write(*,'(a,g0)') '  returned as baseflow     : ', this%bud_back(ROLE_BASEFLOW)
    if (this%have_surfdis) write(*,'(a,g0)') '  returned at land surface : ', this%bud_back(ROLE_SURFACE_DISCH)
    if (this%have_gwet)    write(*,'(a,g0)') '  taken as groundwater ET  : ', this%bud_back(ROLE_GW_ET)
    if (this%have_evt)     write(*,'(a,g0)') '  ET demand sent           : ', this%bud_etdem
    if (abs(pct) > 1.d0) write(*,'(a)') '  NOTE: a non-zero mapping residual is the HRU/cell area mismatch, '// &
      'not the coupling lag; see the HRU area warnings at start-up.'
  end subroutine mf6_budget_report

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

  subroutine to_lower(s)
    character(len=*), intent(inout) :: s
    integer :: i, k
    do i = 1, len_trim(s)
      k = iachar(s(i:i))
      if (k >= iachar('A') .and. k <= iachar('Z')) s(i:i) = achar(k + 32)
    end do
  end subroutine to_lower

  ! Boundary-package role name -> role code; 0 if the name is not one we know.
  integer function role_code(name) result(ir)
    character(len=*), intent(in) :: name
    integer :: i
    ir = 0
    do i = 1, size(ROLE_NAME)
      if (trim(name) == trim(ROLE_NAME(i))) then; ir = i; return; end if
    end do
  end function role_code

end module mf6_coupling
