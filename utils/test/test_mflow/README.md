# SUMMA / MODFLOW 6 coupled test cases

A thin BMI coupler runs SUMMA as the land model and MODFLOW 6 as the groundwater model.
This folder holds the coupler driver script, the shared MODFLOW 6 model, and two cases.

## Cases

| script | case | HRUs | `groundwatr` | what it is for |
|---|---|---|---|---|
| `./run_sagehen1.sh` | `domain_sagehen1` | 1 | `modflow` | the basic one-way/two-way coupling |
| `./run_sagehen4.sh` | `domain_sagehen4` | 4 | `modLatflow` | lateral flow between HRUs; guards the cascade ordering in `run_oneGRU` |
| `./run_sagehen4_noLatflow.sh` | `domain_sagehen4` | 4 | `modflow` | the same four HRUs without lateral flow |

`domain_sagehen4` carries two file managers, `fileManager_latflow.txt` and
`fileManager_noLatflow.txt`, whose decision files differ on the `groundwatr` line alone, so
running the pair isolates what lateral flow contributes. They write `run1_latflow*` and
`run1_noLatflow*`, so neither overwrites the other.

Both share the MODFLOW 6 model in `ex-gwf-sagehen` and run the same 72 hourly steps.
Each case carries its own `summa_modflow6.config`, so they can differ in package names,
HRU→cell map and feedback while pointing at one model directory.

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
   is regridded onto the MODFLOW 6 grid and written into the RCH package `RECHARGE` array.
4. MODFLOW 6 advances one step (prepare / do / finalize_time_step).
5. the new head field is read back and aggregated per HRU, ready for step 1 of the next
   iteration — so the exchange is **explicit, with a one-step lag**.

With `feedback = .true.` the coupler also writes back, each step, so that SUMMA's water
balance and routed streamflow include the aquifer:

    scalarAquiferStorage  = Sy * (MODFLOW water table - soil-column base)
    scalarAquiferBaseflow = MODFLOW <bflow_package_name> outflow over the HRU footprint

`scalarAquiferRecharge` is not exchanged — SUMMA sets it to its own soil drainage.

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
