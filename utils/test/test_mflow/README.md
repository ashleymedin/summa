# SUMMA / MODFLOW 6 coupled test cases

A thin BMI coupler runs SUMMA as the land model and MODFLOW 6 as the groundwater model.
This folder holds the coupler driver script, two MODFLOW 6 models, and five cases.

## Cases

| script | case | HRUs | `groundwatr` | what it is for |
|---|---|---|---|---|
| `./run_sagehen1.sh` | `domain_sagehen1` | 1 | `modflow` | the basic one-way/two-way coupling |
| `./run_sagehen4.sh` | `domain_sagehen4` | 4 | `modLatflow` | lateral flow between HRUs; guards the cascade ordering in `run_oneGRU` |
| `./run_sagehen4_noLatflow.sh` | `domain_sagehen4` | 4 | `modflow` | the same four HRUs without lateral flow |
| `./run_sagehen1_steady.sh` | `domain_sagehen1` | 1 | `modflow` | steady-state first stress period, and a DRN at land surface returning groundwater discharge as surface runoff |
| `./run_sagehen1_deeproot.sh` | `domain_sagehen1` | 1 | `modflow` | groundwater evapotranspiration through an EVT package |

`domain_sagehen4` carries two file managers, `fileManager_latflow.txt` and
`fileManager_noLatflow.txt`, whose decision files differ on the `groundwatr` line alone, so
running the pair isolates what lateral flow contributes. They write `run1_latflow*` and
`run1_noLatflow*`, so neither overwrites the other.

The first three share the MODFLOW 6 model in `ex-gwf-sagehen`; the last two use
`ex-gwf-sagehen-ss`, which adds a steady-state first stress period plus DRN and EVT packages.
All five run the same 72 hourly steps. Each case carries its own `summa_modflow6.config`, so
they can differ in package names, roles, HRU→cell map and feedback.

`run_sagehen1_deeproot.sh` is a mechanism test rather than a realistic Sagehen setup: it sets
`rootingDepth` to 6 m against a 4 m soil column with `rootProfil = doubleExp` and both root
scale factors at their lower bound, because the bundled parameters otherwise place 99.98% of
roots inside the soil column and groundwater ET never switches on. Note that `powerLaw` clamps
the rooting depth to the soil depth, so it can never place roots below the column.

Build the coupler first, with `-DUSE_MODFLOW6=ON` (see `build/cmake/build_mflow.mac.bash`);
it lands in `srcextern/summa/bin/summa_modflow6.exe`.

## How the coupling works

Each SUMMA data step:

1. **(feedback)** the MODFLOW 6 water-table head from the previous step is written into SUMMA
   as the prescribed-head lower boundary of the soil column
   (BMI input `soil_water_sat-zone_top__head`, parameter `lowerBoundHead`).
2. SUMMA advances one step.
3. the drainage out the base of the SUMMA soil column
   (BMI output `soil_water__drainage_volume_flux`, flux `scalarSoilDrainage`)
   is regridded onto the MODFLOW 6 grid and written into the RCH package `RECHARGE` array,
   and the aquifer transpiration demand into the EVT package `RATE` array.
4. MODFLOW 6 advances one step (prepare / do / finalize_time_step). Leading steady-state
   stress periods are solved out first, using the RCH package's own recharge rather than
   SUMMA's.
5. the new head field and each named boundary package's flow are read back and aggregated per
   HRU, ready for step 1 of the next iteration — so the exchange is **explicit, with a
   one-step lag**.

The scatter happens after `prepare_time_step`, not before: `prepare_time_step` reloads a
package's `PERIOD` block at the start of each stress period, which would overwrite it.

With `feedback = .true.` the coupler also writes back, each step, so that SUMMA's water
balance and routed streamflow include the aquifer:

    scalarAquiferStorage    Sy * (MODFLOW water table - soil-column base), diagnostic only
    scalarAquiferBaseflow   role=baseflow package outflow over the HRU footprint
    mfSurfaceDischarge      role=surface_discharge outflow, added to SUMMA's surface runoff
    scalarAquiferTranspire  role=gw_et extraction, against the demand SUMMA sent
    scalarTranspireLimAqfr  aquifer transpiration limiting factor, per MODFLOW cell

`scalarAquiferRecharge` is not exchanged — SUMMA sets it to its own soil drainage.

Boundary packages are named with a **role** rather than individually, so what SUMMA does with
a returned flux is separate from which MODFLOW package supplied it:

    bnd_package_names  = 'CHD', 'DRN', 'EVTA'
    bnd_package_roles  = 'baseflow', 'surface_discharge', 'gw_et'

Packages sharing a role are summed. `bflow_package_name` still works as a one-entry
`baseflow` table. `evt_package_name` names the EVT package whose `RATE` array SUMMA's demand
is written into; its `SURFACE` should sit at the base of the soil column and its `DEPTH` be
large, because the demand has already been reduced by `scalarTranspireLimAqfr` and a second
extinction ramp in MODFLOW would apply the same limit twice.

`head_restart_read` / `head_restart_write` save and reload the MODFLOW head field, so an
aquifer spun up once can start every later run. The calibration driver wires this up
automatically, per rank.

Start-up reports each HRU's area against its effective mapped cell area, and the run ends with
a coupled water budget: what SUMMA sent, what MODFLOW's RCH array received, what came back by
role, and the residual.

## SUMMA model decisions

| decision | `modflow` | `modLatflow` |
|---|---|---|
| `groundwatr` | `modflow` | `modLatflow` |
| `bcLowrSoiH` | `presHead` | `presHead` |
| `hc_profile` | any | **`exp_prof`** |
| `infRateMax` | any valid for `hc_profile` | `topmodel_GA` or `noInfExc` |

`modLatflow` does everything `modflow` does and additionally runs TOPMODEL-style lateral
flow through the soil column above the water table, for hillslopes where water moves
downslope through the soil as well as recharging the aquifer. It requires `exp_prof`
because the lateral transmissivity is the vertical integral of conductivity over the soil
column alone — MODFLOW carries everything below — and `exp_prof` is the profile that
integrates to a finite base rather than assuming a shallow aquifer of its own. The lateral
flow is reported as `basin__ColumnOutflow` and added to total runoff alongside the MODFLOW
baseflow.

Both options need SUMMA built with `-DUSE_MODFLOW6=ON` (which sets `MODFLOW_ACTIVE`) and
must run through `summa_modflow6.exe`; a plain `summa` run is rejected at start-up.

## Requirements on the MODFLOW 6 model

Read from `mfsim.nam` in the working directory, it must:

* use length unit metres and `TDIS TIME_UNITS SECONDS`;
* have exactly **one MODFLOW time step per SUMMA forcing data step** — MODFLOW `delt` must
  equal SUMMA's `data_step`, taken from the forcing file attribute and constant for the run;
* contain an RCH package with `READASARRAYS`;
* be a single GWF model discretised with `DIS`.

All of these are checked by the coupler and reported by name if violated — the units, grid and
RCH checks at start-up, the time-step check on the first step — so you do not have to discover
them from a run that merely looks wrong.  Units left undefined in `mfsim.nam` draw a warning
rather than a stop, since a missing declaration is not proof of the wrong unit.

## Running

    ./coupler_commands.sh -c CONFIG MODFLOW_CASE SUMMA_FILEMANAGER [summa_modflow6.exe]

`-c/--config` is **required**: each case keeps its own config beside its settings. The
executable defaults to `../bin/summa_modflow6.exe`.

### MPI (`summa_modflow6_mpi.exe`)

Built alongside the plain coupler when both `-DUSE_MODFLOW6=ON` and `-DUSE_MPI=ON` are set.
It splits SUMMA's GRUs across MPI ranks the same way `summa_mpi.exe` does for a plain
(non-coupled) run; MODFLOW 6 itself stays a serial singleton on rank 0, since this libmf6
build has no PETSc/MPI support and the coupled model is one shared aquifer grid under
potentially many GRUs, not one MODFLOW grid per GRU (see the header of
`summa_modflow6_mpi.f90` for the full rationale). In `build/CMakeLists.txt` it links both
halves of that: `parallel_utils` for the GRU split (as `summa_mpi.exe` does) and
`bmif`/`MF6_LIB` for driving MODFLOW 6 (as the plain `summa_modflow6.exe` does). Run it
directly (it takes the same file manager/config arguments as the plain coupler, not through
`coupler_commands.sh`):

    cd MODFLOW_CASE && mpirun -n <nranks> /path/to/summa_modflow6_mpi.exe SUMMA_FILEMANAGER CONFIG

`<nranks>` must not exceed the number of GRUs in the run domain - a rank left with zero
GRUs hits a pre-existing SUMMA-core limitation in forcing-file reading, independent of
MODFLOW coupling (the same constraint applies to the plain `summa_mpi.exe`).

### Parameter calibration (`summa_modflow6_opt.exe`)

The calibration driver evaluates each parameter sample as a full coupled SUMMA/MODFLOW 6
simulation. It is the same program as the uncoupled `summa_opt.exe` - with `-DUSE_MODFLOW6=ON`
the calibration executable is simply renamed - and it couples whenever `simulation.use_modflow`
is set in the TOML configuration, so one binary covers both. See
`utils/test/test_calibration/README.md` for how to configure and run it.

## Where the coupling lives

The MODFLOW side is one module, `build/source/driver/mf6_coupling.f90`: the libmf6 bindings,
the HRU-to-cell map, the start-up checks, and the per-step exchange. Its interface is per-HRU
arrays in and per-HRU arrays out, so it knows nothing about SUMMA data structures, and all
three drivers share it - the serial coupler, the MPI coupler, and the calibration driver.
(The two couplers previously carried the same ~600 lines twice, verbatim.)

The SUMMA side is `build/source/driver/summa_mf6_exchange.f90`: the per-HRU getters and
setters for the exchanged quantities, plus the HRU geometry the cell map is built from. It
too has one implementation with two users - `summa_bmi.f90`'s `get_value`/`set_value` for the
coupler variable names delegate to it, and the calibration driver, which runs SUMMA without
the BMI, calls it directly.

## Coupler configuration

A Fortran namelist, conventionally `summa_modflow6.config` inside the case directory:

    &coupler
      mf6_model_name     = 'MYMODEL' ! GWF model name, as in mfsim.nam (upper case)
      rch_package_name   = 'RCHA'    ! RCH package name, as in the GWF name file (upper case)
      bflow_package_name = 'CHD'     ! head-dependent boundary package (CHD/DRN/RIV/GHB) whose
                                     !   simulated flow feeds back per HRU as scalarAquiferBaseflow
                                     !   ('' => skip the baseflow feedback)
      map_file           = ''        ! optional HRU->cell weight file; if blank a nearest-cell
                                     !   map is built from the MODFLOW 6 DIS grid geometry
      mf6_epsg           = 0         ! EPSG code of the MODFLOW grid's projected CRS, used only to
                                     !   reproject SUMMA's lon/lat HRU centres before the built-in
                                     !   nearest-cell map (nHRU>1).  Only WGS84 UTM is supported
                                     !   (32601-32660 N, 32701-32760 S); 0 = no reprojection, which
                                     !   is fine when nHRU==1 or with an explicit map_file
      feedback           = .true.    ! .false. => one-way (SUMMA drainage -> MODFLOW only)
    /

### `map_file` format

One `iHRU  cell  weight` triple per line, whitespace separated; blank lines and `#` comments
ignored. `cell` is the row-major horizontal MODFLOW index `(irow-1)*ncol + icol`. Weights are
normalised per HRU, so `1.0` on every line spreads an HRU evenly over its cells. An HRU may
span any number of lines.

## Mapping HRUs to cells

HRUs may each own a disjoint set of cells — down to one cell per HRU. They need not overlap
and need not cover the whole grid; `domain_sagehen4` uses four disjoint regions.

The one requirement is that an HRU's `elevation` in `attributes.nc` be the mean land surface
of the cells it maps to, because the coupler forms

    lowerBoundHead = h_mf6 - (z_surface_HRU - soil_thickness)

so a mismatch hands the soil column a water table metres above or below it. The coupler
checks this against `DIS/TOP` at start-up and stops with the offending HRU named, rather
than letting SUMMA fail to converge for no visible reason.
