# test_calibration

Tests and worked examples for SUMMA's parameter calibration (`summa[_sundials]_opt.exe`).

## Requirements

Calibration needs an executable built with **both** `-DUSE_MPI=ON` and `-DUSE_MIZUROUTE=ON`.

MPI is what builds the calibration driver in the first place, and it is how trials are run
concurrently: one rank coordinates, the others evaluate parameter sets.

mizuRoute is needed because the objective function compares *routed streamflow* against gauge
observations, and routed streamflow is produced by mizuRoute. In a build without it the
simulated flow series is never filled, so there is nothing to score. This is a real constraint
of the current design, not a property of these tests.

## Calibration targets

A calibration scores each parameter trial against one or more **targets**. A target names an observed
series, the simulated variable it is compared against, and the metric that scores the two together.

A configuration that names no target is calibrating the one streamflow series in `[observations]`, and
is read as a calibration with exactly one target - which is what every configuration written before
targets existed does, unchanged. To score a trial against more than one thing, list the targets:

```toml
[[calibration.target]]
name      = "discharge"
variable  = "streamflow"
obs_file  = "CAN_05BB001_daily_flow_observations.nc"
vname_obs = "q_obs"
metric    = "kge"
weight    = 1.0

[[calibration.target]]
name      = "discharge_error"
variable  = "streamflow"
obs_file  = "CAN_05BB001_daily_flow_observations.nc"
metric    = "rmse"
weight    = 0.1
```

Only `obs_file` is required. `obs_path`, `vname_obs`, `metric` and `obs_transform` fall back to the
`[observations]` and `[calibration]` settings when a target does not state its own.

Each trial runs the model **once** and scores that one simulation against every target, because the
targets are different views of the same simulation and have to come from the same run to be
comparable. The trials file records every target's value: `objective` gains a `target` dimension, and
`target_name` labels it.

DDS searches on one number, so it collapses the targets into a weighted sum. Targets are oriented
before they are combined - efficiencies (KGE, KGE', NSE) count as they stand, error metrics (MAE,
RMSE) count negatively - so a larger sum is always a better fit whatever mix of metrics is used. With
a single target of weight one that sum is the metric itself, which is exactly what DDS maximized
before targets existed.

`variable` accepts `streamflow` (or `discharge`) for the routed flow mizuRoute produces, or the name
of any variable SUMMA knows - `basin__StorageChange`, `scalarStreamTemp`, `lowerBoundHead`. A name
SUMMA does not know is refused when the calibration starts, rather than after a simulation has been
paid for.

### Which part of the domain a target scores

By default a SUMMA variable is collected as an area-weighted mean over the whole domain, and
`streamflow` is taken at the network outlet. That is right for a single basin and wrong for anything
else: two basins in one domain hold their own storage, and a stream temperature belongs to one reach.
A target can name the unit it scores, by the id the domain or the river network gives it:

```toml
[[calibration.target]]
name     = "gulkana_storage"
variable = "basin__StorageChange"
gru      = 1                  # one GRU: the area-weighted mean of its HRUs

[[calibration.target]]
name     = "upstream_flow"
variable = "streamflow"
reach    = 710289040          # one reach of the mizuRoute network, by its segId
```

`gru`, `hru` and `reach` are the three choices, and a target names at most one. An id the domain or
the network does not hold is refused at start-up. A variable SUMMA holds once per GRU, like
`basin__StorageChange`, can be scored over a `gru` but not an `hru`, and a `reach` carries routed
`streamflow` and `T_reach` only. Two targets naming the same variable over the same unit share one collected series.

### Stream temperature

`T_reach` (or `stream_temperature`) is the water temperature leaving a reach, in degC, and needs a
`reach`. It is scored in degC rather than K because a KGE's bias ratio would otherwise sit near one
for any plausible error. The network carries a temperature only when some GRU has a stream HRU, and a
`T_reach` target in a network without one is refused at start-up.

```toml
[[calibration.target]]
name      = "stream_temperature"
variable  = "T_reach"
reach     = 9
obs_file  = "USGS_10343500_daily_temperature.nc"   # from utils/pre-processing/acquire_stream_temperature.py
vname_obs = "t_obs"
metric    = "rmse"
```

Without mizuRoute, `streamflow` over a `gru` is that GRU's `averageRoutedRunoff` times its area, in
m3 s-1, so a gauge below a one-GRU domain can be scored with no river network. With mizuRoute running,
streamflow comes from the network only, and a `gru` streamflow target is refused: name the `reach`.

### Observations

Every target reads its observations the same way, from a NetCDF file shaped like the bundled streamflow
record: a `time` coordinate with CF units, and the named variable carrying its own units. Each time
stamps the **end** of the span the value covers, on the clock the forcing is on, as SUMMA stamps its own
time steps; the simulation is averaged over the span since the observation before it, so a monthly product, whose months are 28 to 31 days long, aligns like a daily
gauge record. Getting data into that form is a pre-processing job - see `utils/pre-processing/`, whose
acquisition scripts write it. `obs_units` overrides the units the file states.

### Comparing things that are not measured the way the model carries them

Some observations cannot be compared against a model variable as it stands. GRACE reports basin water
storage once a month, as a departure in millimetres from a multi-year mean; SUMMA carries the rate
storage is changing at, every time step. A target says what has to happen for the two to be
comparable:

```toml
[[calibration.target]]
name           = "grace_tws"
variable       = "basin__StorageChange"   # kg m-2 s-1, which is mm of water per second
gru            = 1
obs_file       = "gulkana_grace_tws.nc"
vname_obs      = "tws_obs"
accumulate     = true                     # integrate the rate into the storage itself
baseline_start = "2009-10-01"             # express both sides as departures from their mean
baseline_end   = "2009-12-31"             #   over the same months
metric         = "rmse"
```

- `accumulate` integrates a rate into the quantity the observations report. Integrating kg m-2 s-1
  over seconds leaves kg m-2, which is millimetres of water.
- `baseline_start`/`baseline_end` express **both** series as departures from their own mean over that
  period. Doing it to both sides is what removes the constant of integration an accumulated series
  carries, which is what makes an integrated rate comparable to a storage anomaly at all.

### Glacier mass balance

A glaciological mass balance is measured between the glacier's own seasonal extremes: winter from the
autumn minimum to the spring maximum, summer from that maximum to the next minimum, annual from
minimum to minimum. `balance` compares each observed balance with the same change in the simulated
storage, between the extremes the simulation itself has:

```toml
[[calibration.target]]
name       = "mass_balance"
variable   = "basin__GlacierMassChange"   # kg m-2 s-1 over the glacier area, not the basin's
gru        = 1
obs_file   = "gulkana_mass_balance_seasonal.nc"
vname_obs  = "mb_obs"
accumulate = true                         # integrate into glacier storage, mm w.e.
balance    = "seasonal"                   # seasonal or annual
metric     = "rmse"
```

- `seasonal` scores the winter and summer balances of one file, in date order. At each observed date
  the simulated maximum and minimum within `balance_window` days (default 60) are found, and the one
  that turns inside the window - the maximum in spring, the minimum in autumn - ends the balance; it
  starts at the opposite extreme in the year before. A melt season that runs early or late is scored
  on its mass, not its dates.
- `annual` uses minima only: the minimum near the observed date less the minimum near a year earlier.
- A balance whose search reaches outside the simulation, or whose turning point is not clear, is left
  out, so the first year of a run scores nothing.
- Score the mass balance or GRACE alongside discharge, not both: they measure the same storage, one
  over the glacier and one over a ~300 km mascon, and on Gulkana they pull the parameters apart.
- `basin__GlacierMassChange` is the glacier domains' total mass change over the glacier area, which is
  what a glacier-wide balance is per unit of.

`utils/pre-processing/acquire_glacier_mass_balance.py` writes the USGS Benchmark Glacier glacier-wide
solutions in this form, seasonal or annual.

## Parameters that must stay in order

Some parameters are only valid in order: SUMMA refuses a soil unless
`theta_res <= critSoilWilting <= critSoilTranspire <= fieldCapacity <= theta_sat`. An ordered chain
says so, and keeps adjacent members at least `gap_fraction` of the chain's whole range apart:

```toml
[[parameter_dependencies.ordered]]
parameters   = ["theta_res", "critSoilWilting", "critSoilTranspire", "fieldCapacity", "theta_sat"]
gap_fraction = 0.025
sample_gaps  = true
```

Without `sample_gaps`, only the members in `param_list` are searched, the rest stay at their
defaults, and a trial out of order is drawn again. So calibrating `theta_sat` alone can only move it
above the default `fieldCapacity`.

With `sample_gaps = true`, every member is searched, whether or not it is in `param_list`. The first is
searched as itself, up to the highest value that still leaves room for the rest; each later one as the
fraction of the room between the one before it (plus the gap) and the highest value it can take. Every
trial is in order and within bounds by construction, so none is drawn again or refused. A later member
cannot also be given a transformation. The trials file holds the values, not the fractions.

## Searching for the trade-offs: `algorithm = "nsga2"`

DDS collapses the targets into a weighted sum and returns one parameter set, so when targets pull
against each other the weights have decided the answer before the search starts. NSGA-II (Deb et al.,
2002) searches on every target at once and keeps the **Pareto front**: the trials no other trial beats
on every target together. Choosing among them is left to whoever reads the output.

```toml
[calibration]
algorithm = "nsga2"      # default "dds"
n_samples = 400          # population_size x generations

[calibration.nsga2]
population_size = 20     # members kept each generation (default 100)
# crossover_probability = 0.9, crossover_eta = 20, mutation_eta = 20
# mutation_probability defaults to 1 / (number of calibrated parameters)
```

The first generation is random. Each later one breeds `population_size` offspring from the one before
- a crowded binary tournament picks the parents, simulated binary crossover recombines them and
polynomial mutation perturbs the children - and keeps the best `population_size` of parents and
offspring together. The operators work in the transformed search space, so a `log` parameter is
recombined on a log scale, and a child that breaks an ordered constraint is bred again. `n_samples` has
to be a whole number of generations, at least two. Target weights play no part; each target is
oriented by its metric. Both searches evaluate trials on the same worker pool.

An NSGA-II trials file carries, beside every trial's parameters and objectives:

| Variable | What it holds |
| --- | --- |
| `pareto_front(sample)` | 1 if no other trial is at least as good on every target and better on one |
| `birth_generation(sample)` | the generation that proposed the trial, 1 the random initial population |
| `population(member, generation)` | the sample index of each member kept after that generation |
| `population_rank(member, generation)` | its non-domination front, 1 the non-dominated |
| `population_crowding(member, generation)` | its crowding distance in that front, infinite at the extremes |

`objective:sense` lists `maximize` or `minimize` per target, which is what a reader needs to recompute
dominance.

## Failed trials

A parameter set the model refuses (`paramCheck`) or fails on during the run (a balance check, a solver
failure) no longer ends the calibration. The trial is flagged in `trial_failed(sample)`, its objectives
are written as fill values, DDS never takes it as its best, and NSGA-II ranks it last on every target.
Any other error, such as a missing observation file, still stops the run.

Each failed trial is saved as a repro case in `<output_path>/failed_trials/`, named
`nsga2_g<generation>_s<sample>/` or `dds_s<sample>/`:

| File | What it holds |
| --- | --- |
| `trial.json` | the parameter values, the error with the GRU, time step and HRU it failed on, the worker, the source revision |
| `config.toml` | the calibration's configuration file |
| `initial_state.nc` | the spun-up state every trial starts from |
| `rerun.sh` | reruns the trial with the serial executable of the same build (override with `SUMMA_EXE`) |

`rerun.sh` writes `output/` and `rerun.log` beside itself, so the directory can be copied elsewhere and
debugged on its own. A run from a manifest is built from a template, so its `rerun.sh` only says so.

## Calibrating a MODFLOW 6 coupled case

A build with `-DUSE_MODFLOW6=ON` names its calibration executable `summa_modflow6_opt.exe`.
It is the same driver as `summa_opt.exe`; what changes is that each parameter sample can be
evaluated as a full coupled SUMMA/MODFLOW 6 simulation rather than a plain SUMMA run. Turn
that on per case in the TOML configuration:

```toml
[simulation]
use_modflow = true

[modflow]
config_file = "/abs/path/to/summa_modflow6.config"   # the &coupler namelist, as for the standalone couplers
run_dir     = "/abs/path/to/the/MODFLOW/model"       # the directory holding mfsim.nam
```

The model decisions must be the coupler's, exactly as for `summa_modflow6.exe`:
`groundwatr = modflow` (or `modLatflow`) and `bcLowrSoiH = presHead`. The run is rejected at
start-up otherwise.

Each MPI rank is an independent model instance, so each gets its own MODFLOW run directory,
`modflow_rank####` under the case's output path, populated once by copying `run_dir` and
reused by every sample that rank evaluates. They are copies rather than symlinks because
MODFLOW writes its listing, head and budget files into its working directory, and ranks
sharing one directory would overwrite each other's output. Budget for the disk this needs:
one copy of the MODFLOW model per rank.

Within a rank, every sample runs its own complete MODFLOW simulation, started and finalized
around the SUMMA run, so each trial sees the same aquifer initial condition and trials stay
comparable.

**The aquifer is spun up with the soil column.** The calibration's one-year cold-start spinup
runs coupled and writes its final head field to `modflow_spinup_heads_rank####.bin` under the
case's output path; every sample that rank evaluates then starts MODFLOW from that field in
place of the `STRT` array in the model's own IC package, just as SUMMA starts from the spun-up
restart file. The spinup itself still begins at `STRT`, so it is worth setting to a sensible
water table, but the calibration period no longer depends on it.

SUMMA's own aquifer follows the same rule. With `aquiferIni = emptyStart` the aquifer is
emptied on the cold start only: the spinup's restart file records that a start has already
happened, so each sample keeps the spun-up aquifer rather than emptying it again.

**This needs mizuRoute as well**, for the same reason every other calibration does: the
objective compares routed streamflow against gauge observations. A coupled calibration build
is therefore `-DUSE_MPI=ON -DUSE_MIZUROUTE=ON -DUSE_MODFLOW6=ON`.

`test_calibration_sagehen.sh` below is the bundled end-to-end test of this path.

## `test_calibration_bow.sh`

A short, self-contained calibration on the Bow River at Banff (CAN_05BB001). Everything it needs
is already in the repository, under `utils/test/test_mizuroute/bow_real_data/`: lumped SUMMA
inputs and RDRS forcing for 2002-2017, the mizuRoute topology and lumped-to-HRU remapping, and
daily streamflow and GRACE observations. It calibrates 2004 after a 2003 spinup, inside the GRACE
record, so a storage target can be added to the same run:

```toml
[[calibration.target]]
name           = "grace_tws"
variable       = "basin__StorageChange"
obs_path       = ".../bow_real_data/observations/"
obs_file       = "CAN_05BB001_grace_tws.nc"
vname_obs      = "tws_obs"
accumulate     = true
baseline_start = "2004-01-01"
baseline_end   = "2004-12-31"
metric         = "rmse"
```

```bash
./test_calibration_bow.sh [n_samples] [n_ranks]     # defaults: 4 samples, 2 ranks
```

It writes a TOML configuration pointing at that data, runs a DDS search over five parameters for
one simulated year after a year of spinup, and checks that finite, plausible objective values were
recorded for the trials.

The period and the sample count are deliberately tiny so the test finishes in minutes. It
exercises the machinery; it does **not** produce a calibrated parameter set. A real calibration
uses thousands of samples over several years — see `n_samples` in the generated config.

## `test_calibration_bow_pareto.sh`

The Bow 2004 calibration scored against two targets, daily discharge (KGE) and the GRACE storage
anomaly (RMSE), with the same budget spent twice: by NSGA-II on the two objectives and by DDS on their
weighted sum. It checks the NSGA-II trials file - finite objectives, a `pareto_front` that is exactly
the non-dominated trials, a population drawn only from parents and offspring - then prints both
fronts and their hypervolumes.

```bash
./test_calibration_bow_pareto.sh [population] [generations] [n_ranks]   # defaults: 8, 5, 5
```

It calibrates `k_soil`, `aquiferScaleFactor` and the two routing parameters. `theta_sat` is left out:
sampling it needs the soil ordering constraint, and under that constraint some trials stop SUMMA
with a soil water balance error.

## `test_calibration_sagehen.sh`

A coupled SUMMA / MODFLOW 6 / mizuRoute calibration on Sagehen Creek, scoring discharge (KGE) and
stream temperature (`T_reach` RMSE, degC) at USGS 10343500, the outlet of reach 9, and the water level
in a synthetic well (`lowerBoundHead` RMSE, m). NSGA-II scores water year 2018 after a water-year-2017
spin-up.

Sagehen has no observation well, so the test makes one up. It sits on the lower valley side two cells
from the channel above the gauge (row 41, column 75, in GRU 7), where the water table is about 15 m
down, below 11 m of UZF. The 2017 snowmelt reaches it in October, 0.62 m above its May low, and it falls
0.47 m through dry water year 2018; the valley floor
is held at land surface by the drains. Its record is the head there in a
reference coupled run at the default parameters, written by
`utils/test/test_mflow/tools/make_sagehen_synthetic_well.py` as a saved USGS daily-values response
(`sagehen/observations/SYNTHETIC-SAGEHEN-1_daily_values.json`). `acquire_groundwater_level.py --features`
turns that into the observation file exactly as it would a real well's. The target scores departures
from the water-year mean, since a well's datum is not the model's; on the `lumped` layout it is GRU 7's
land HRU, averaged over 103 cells, that is compared against one cell's head.

To remake the record, run the lumped domain once for water years 2017-2018 from its settled
`strt1.txt`, with the coupler and the MODFLOW OC saving `HEAD FREQUENCY 6`, then

```bash
make_sagehen_synthetic_well.py <run>/sagehen.hds "2016-10-01 00:00" sagehen/observations/SYNTHETIC-SAGEHEN-1_daily_values.json
../../pre-processing/acquire_groundwater_level.py SYNTHETIC-SAGEHEN-1 2016-10-01 2018-09-30 \
    sagehen/observations/SYNTHETIC-SAGEHEN-1_daily_level.nc --model-utc-offset 0 --gauge-utc-offset -8 \
    --features sagehen/observations/SYNTHETIC-SAGEHEN-1_daily_values.json
```

```bash
./test_calibration_sagehen.sh [lumped|grid] [population] [generations] [n_ranks]   # defaults: lumped, 6, 3, 7
```

The domain is built at run time by `utils/test/test_mflow/tools/build_sagehen9_calibration.py` from the
9 D8 subcatchments of `ex-gwf-sagehen`, each with a stream HRU for its reach, and the bundled basin-mean
forcing in `sagehen/forcing/`. `lumped` gives each GRU one land HRU mapped onto all its cells (18 HRUs,
minutes a water year, UZF taking a week's drainage at a time); `grid` is `domain_sagehen9`, one land HRU per cell (3396 HRUs, hours a
trial). The MODFLOW model is copied with a TDIS as long as the forcing. It needs
`bin/summa_modflow6_opt_sundials_mizuroute.exe` and a python3 with netCDF4 and pyproj.

At 4 x 2 on 5 ranks it takes about 20 minutes, and every trial runs. The front runs from
KGE -0.22 at 2.27 degC RMSE to KGE -1.14 at 2.22 degC, over `frozenPrecipMultip`, `tempCritRain`,
`k_macropore` and `routingGammaScale`; the well misses its record by 0.18 to 0.28 m.

## `test_calibration_wolverine.sh`

A glacier calibration routed by mizuRoute, scoring discharge (KGE) and water temperature (RMSE, degC)
at USGS 15236900, and the glacier-wide seasonal mass balance (RMSE), on the bundled `wolverine/`
domain: three glacier land HRUs and a stream HRU for the one reach. It spins up over 2016 and scores
2017–2019, over `frozenPrecipMultip`, `tempCritRain`, `albedoMax`, `glacierWindFactor`,
`glacierTempReduction`, the three glacier storage constants `glacStor_kIce`, `glacStor_kSnow` and
`glacStor_kFirn`, and `streamWidthMultip` (0.5–3), the reach's surface area. Under
`deepTherml = bedrockLyrs` the groundwater reaching the channel leaves at the bedrock base, which no
parameter sets; with `deepTherml = airTempGW` the search would take `C_ATGW` and `gwTempWindow` instead.
`hypFrac` and `hypLag` are left at their defaults: the hyporheic proxy damps the daily swing, which a
daily-mean target cannot see, and halving `hypLag` leaves the RMSE unchanged.

```bash
./test_calibration_wolverine.sh [population] [generations] [n_ranks]   # defaults: 8, 4, 5
```

It needs `bin/summa_sundials_mizuroute_opt.exe`. At the defaults it takes about 10 minutes, and every
trial runs. The front runs from KGE 0.87 at 0.75 degC and 531 mm, through KGE 0.83 at 0.82 degC and
407 mm, to 0.62 degC at KGE 0.35 and 1,421 mm. Glacier runoff enters the reach at freezing; the water
warms from the groundwater of the unglaciated third of the basin, which leaves at the bedrock's
starting temperature, and from the heat friction dissipates down the reach.
`wolverine/README.md` says how the domain and observations were built.

## `multi_case_example/` -- multi-case calibration

`--manifest <file>` calibrates many basins in one job, from a manifest listing the cases and a
per-case configuration template that the manifest expands.

The test dataset for this is built by `make_stub_century.bash`, which assembles the directory
layout the template expects from the one domain bundled in this repository, the Bow at Banff.
Files are symlinked, so a three-case dataset costs about 128 KB:

```bash
cd multi_case_example
./make_stub_century.bash                  # defaults: stub_century/, 3 cases, 1 at a time
mpirun -np 2 ../../../../bin/summa_sundials_mizuroute_opt.exe \
       --manifest stub_century/stub_manifest.toml
```

Arguments are `[root] [n_cases] [cases_per_node]`. Each case group needs at least two ranks --
one coordinates the search, the rest evaluate samples -- and the groups have to divide the node's
ranks evenly. If the rank count does not fit, the driver warns and reduces `cases_per_node` rather
than refusing to start, so `-np` equal to `2 * cases_per_node` gets the cases you asked for and
anything else still runs. The script prints a matching command when it finishes.

Every case is the same basin under a different name, so the calibrated values are not meaningful.
What it exercises is the part with no other coverage: reading the manifest, expanding the template
per case, distributing cases across ranks, and writing per-case output. Two real bugs in the
multi-case path were found this way, both fixed -- calibration structures were not released
between cases, so the second case failed to allocate, and the template's `work_path` lacked a
trailing slash, so output landed beside the work directory rather than in it.

To run it on a larger dataset, point `home_path` in the template at that dataset and list its
cases in the manifest; the layout is the same.

| File | What it is |
| --- | --- |
| `make_stub_century.bash` | builds the stub dataset, and the manifest and template that drive it |
| `manifest_century.toml` | the full manifest, listing ~115 basins; the stub manifest is derived from it |
| `summa_config_template.toml` | the per-case configuration template the manifest expands |
| `summa_config_CAN_05BB001.toml` | a single filled-in case, useful for seeing what the template produces |
| `century_cases.txt` | the basin list |
| `setup_summa_cases.bash` | builds per-case input directories from a larger dataset, if you have one |
| `check_time.bash` | reports each case's forcing time range, to check a calibration period is covered |

`setup_summa_cases.bash` and `check_time.bash` hard-code a data root (`$HOME/data/century/...`)
and an experiment-specific layout, so they need editing before they run anywhere else. They are
kept because they document how the per-case directories were originally produced;
`make_stub_century.bash` is what builds a dataset you can run today.
