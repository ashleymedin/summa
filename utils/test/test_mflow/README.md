# SUMMA / MODFLOW 6 coupling

A thin BMI coupler runs SUMMA as the land model and MODFLOW 6 as the groundwater model. This folder
holds one test of the coupling and the calibration through it, `sagehen_calibration/`, and documents
the coupler. The coupled cases it was developed on (Sagehen from 1 to 3396 HRUs, Little River, the
twin experiment) and the tools that build them live in `summa_modflow_workflow`.

## `sagehen_calibration/`

A coupled SUMMA / MODFLOW 6 / mizuRoute calibration of lumped Sagehen Creek, against every kind of
target and through every MODFLOW package the coupler drives.

The domain has 9 GRUs, each a land HRU mapped onto all of its cells and a stream HRU for its reach; mizuRoute
routes the 9 reaches. MODFLOW runs `ex-gwf-sagehen-gwe` with:

| package | role |
|---|---|
| RCH, UZF | SUMMA's drainage, through UZF where the water table is below the soil column |
| UZE, GWE (ESL, CND) | the drainage's and the soil base's heat; GWE's water-table temperature is SUMMA's aquifer temperature |
| SFR | baseflow to the streams (role `baseflow`) |
| DRN | water at land surface (role `surface_discharge`) |
| GHBB, DRNB | the base exchange |
| EVTA | SUMMA's aquifer transpiration demand (role `gw_et`), with roots reaching 6 m |
| OBS6 | the head at ten cells, one per decile of depth to water |

The base exchange is a per-cell correction, so on a lumped HRU it tests the path rather than closing
it better. The calibration spins up over water year 2017 and scores October-November 2017.

`run_test.sh` runs three calibrations:

1. A reference trial with every searched parameter pinned to a true value writes the targets
   (`make_targets.py`), as daily means:

   | target | variable | unit scored | metric |
   |---|---|---|---|
   | discharge | `streamflow` | reach 9 | KGE |
   | stream temperature | `T_reach` | reach 9 | RMSE |
   | heads | `modflow_obs`, ten cells pooled | OBS6 | RMSE |
   | storage | `basin__StorageChange`, accumulated | domain | RMSE |
   | baseflow | `basin__AquiferBaseflow` | GRU 5 | KGE |
   | aquifer | `scalarAquiferStorage`, departure from its mean | HRU 5 | RMSE |

2. The same pinned trial, scored against them, must reproduce every one, and its coupled budget must
   show water or heat through each package above.
3. NSGA-II searches the SUMMA and MODFLOW parameters, and each worker's MODFLOW files must hold its last
   trial's multipliers:

   | parameter | truth | searched | scales |
   |---|---|---|---|
   | `K_mult` | 2.0 | 0.2-5, log | `kh1.txt`, `kv1.txt` (arrays) |
   | `Sy_zone3_mult` | 0.7 | 0.5-1.5 | `sy1.txt` in zone 3 of `syzone1.txt` |
   | `rhk_mult` | 3.0 | 0.1-10, log | SFR PACKAGEDATA column 10, streambed K |
   | `vks_mult` | 0.5 | 0.1-10, log | UZF PACKAGEDATA column 8, VKS |
   | `kts_mult` | 1.5 | 0.5-2 | `kts1.txt`, GWE solid conductivity |
   | `k_soil`, `k_macropore` | 3e-6, 1e-2 m s-1 | their `localParamInfo.txt` bounds, log | SUMMA |

```bash
./run_test.sh [population] [generations] [n_ranks]   # defaults: 4, 2, 5
```

It needs `bin/summa_modflow6_opt_sundials_mizuroute.exe`, built with `-DUSE_MPI=ON -DUSE_MIZUROUTE=ON
-DUSE_MODFLOW6=ON`, and a python3 with netCDF4. It took under 2 hours on a 10-core machine shared with another
calibration, most of it the year of spin-up that every run repeats on every rank. The domain is built by
`build_sagehen_calibration_test.py` in `summa_modflow_workflow`; rebuild it there after changing the
model or the builder.

## How the coupling works

Each SUMMA data step:

1. **(feedback)** the MODFLOW 6 water-table head from the previous step is written into SUMMA
   as the prescribed-head lower boundary of the soil column
   (BMI input `soil_water_sat-zone_top__head`, parameter `lowerBoundHead`).
2. SUMMA advances one step.
3. the drainage out the base of the SUMMA soil column
   (BMI output `soil_water__drainage_volume_flux`, flux `scalarSoilDrainage`)
   is regridded onto the MODFLOW 6 grid. Where a cell's water table is below its top, the soil-column
   base, a positive flux enters the UZF package's infiltration, as GSFLOW's gravity drainage enters
   UZF1, held `uzf_hold_days` and then taken at a constant rate over the next so hourly pulses do
   not each start a UZF wave; elsewhere it enters the RCH package `RECHARGE` array. The aquifer transpiration demand goes
   into the EVT package `RATE` array.
4. with the base exchange named, the coupler sets its GHB stage to the last head where that head
   was at or above the top, its DRN at the top where it was below, and both conductances to SUMMA's
   base conductance (BMI output `soil_bottom_surface__hydraulic_conductance`), so MODFLOW solves
   its correction to SUMMA's flux against its own new head.
5. MODFLOW 6 advances one step (prepare / do / finalize_time_step). Leading steady-state
   stress periods are solved out first, using the RCH package's own recharge rather than
   SUMMA's.
6. the new head field, each named boundary package's flow, UZF's rejected infiltration and what the
   base exchange moved are read back and aggregated per HRU, ready for step 1 of the next iteration —
   so SUMMA sees MODFLOW **with a one-step lag**. Rejected infiltration and the base exchange return to the base of the soil column
   (BMI input `soil_water__rejected_recharge_volume_flux`, `scalarAquiferReject`); SUMMA reports its
   drainage net of it, and the coupler adds it back before regridding.

Below the soil-column base SUMMA's lower boundary is free drainage into UZF; at or above it, the
head-dependent `presHead` flux carries water both ways through RCH. Over the first 10 cm below the
base one gives way linearly to the other, so the flux is continuous in the water table. Without the
base exchange, or with a jump there, cell-matched HRUs whose head sits at `DIS/TOP` reverse their
drainage every step: a trace of water past the top lands in confined storage and reads as metres of
head.

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
automatically, per rank. The serial coupler also takes them on its command line, and writing
heads there makes SUMMA write its restart on the last step as well, so a spin-up leaves a
matched pair. A restarted run hands SUMMA the saved water table on step 1 rather than
`lowerBoundHead`, and restart heads replace a leading steady-state solve rather than being
overwritten by it. Boundary fluxes are not saved, so step 1 of a restarted run carries no
baseflow, seepage or groundwater ET back to SUMMA.

Start-up reports each HRU's area against its effective mapped cell area, and the run ends with
a coupled water budget: what SUMMA sent, what MODFLOW's UZF and RCH received, what UZF rejected
back to the soil column, what came back by role, and the residual.

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
* contain a UZF package with one landflag-1 cell per column it covers, and none of `SIMULATE_ET`,
  `LINEAR_GWET`, `SQUARE_GWET`, `SIMULATE_GWSEEP` or `MOVER` (EVT, SUMMA and the coupler supply those);
  with a GWE model, a UZE package for it;
* be a single GWF model discretised with `DIS`.

All of these are checked by the coupler and reported by name if violated — the units, grid and
RCH and UZF checks at start-up, the time-step check on the first step — so you do not have to discover
them from a run that merely looks wrong.  Units left undefined in `mfsim.nam` draw a warning
rather than a stop, since a missing declaration is not proof of the wrong unit.

## Running

The serial coupler runs from the MODFLOW case directory, which holds `mfsim.nam`:

    cd MODFLOW_CASE && summa_modflow6.exe SUMMA_FILEMANAGER CONFIG [TOML]
                       [--head-restart-read HEADS] [--head-restart-write HEADS]

The TOML is needed only with mizuRoute, whose river network it configures. `--head-restart-read` starts
MODFLOW from a saved head field, and SUMMA from whatever restart its file manager names as
`initConditionFile`. `--head-restart-write` writes the head field at the end, and SUMMA's restart with it.
The cases in `summa_modflow_workflow` run it through their `coupler_commands.sh`, which resolves the
paths and runs a spin-up first when its heads do not exist yet.

### MPI (`summa_modflow6_mpi.exe`)

Built alongside the plain coupler when both `-DUSE_MODFLOW6=ON` and `-DUSE_MPI=ON` are set.
It splits SUMMA's GRUs across MPI ranks the same way `summa_mpi.exe` does for a plain
(non-coupled) run; MODFLOW 6 itself stays a serial singleton on rank 0, since this libmf6
build has no PETSc/MPI support and the coupled model is one shared aquifer grid under
potentially many GRUs, not one MODFLOW grid per GRU (see the header of
`summa_modflow6_mpi.f90` for the full rationale). In `build/CMakeLists.txt` it links both
halves of that: `parallel_utils` for the GRU split (as `summa_mpi.exe` does) and
`bmif`/`MF6_LIB` for driving MODFLOW 6 (as the plain `summa_modflow6.exe` does). Run it
directly; it takes the same file manager/config arguments as the plain coupler:

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
      uzf_package_name   = 'UZF'     ! UZF package below the soil column (upper case)
      uze_package_name   = ''        ! GWE UZE package, required with gwe_model_name
      uzf_hold_days      = 1.0       ! days of drainage UZF takes at a time, at a constant rate
      base_ghb_package_name = ''     ! GHB and DRN, one row per top cell, that MODFLOW solves the
      base_drn_package_name = ''     !   base exchange through; set both or neither
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

`iHRU` is SUMMA's internal HRU index, which runs GRU by GRU. With one GRU that is the row
order of `attributes.nc`; with several it is that order only if the file is itself written
GRU by GRU, which `sagehen_calibration/domain` is.

## Mapping HRUs to cells

HRUs may each own a disjoint set of cells — down to one cell per HRU. They need not overlap
and need not cover the whole grid; `sagehen_calibration/domain` gives each GRU's land HRU its own cells.

The one requirement is that an HRU's `elevation` in `attributes.nc` be the mean land surface
of the cells it maps to, because the coupler forms

    lowerBoundHead = h_mf6 - (z_surface_HRU - soil_thickness)

so a mismatch hands the soil column a water table metres above or below it. The coupler
checks this against `DIS/TOP` plus the column depth at start-up and stops with the offending HRU named, rather
than letting SUMMA fail to converge for no visible reason.

`DIS/TOP` must be the base of SUMMA's soil column, as GSFLOW puts its MODFLOW top at the
soil-zone base, so the two models share no pore space. The test's grid has it 4 m below
the native land surface, which is kept as `land1.txt`.
