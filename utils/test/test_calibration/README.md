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

`variable` accepts `streamflow` (or `discharge`) today. Groundwater level, baseflow and stream
temperature arrive with the coupling work that produces them.

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

**The aquifer is not spun up with the soil column.** The calibration's one-year cold-start
spinup runs coupled, but MODFLOW restarts each sample from the initial head in its own input
(`STRT`), while SUMMA restarts from the spun-up state. Set `STRT` to a sensible water table
for the calibration period rather than relying on the spinup to settle it.

**This needs mizuRoute as well**, for the same reason every other calibration does: the
objective compares routed streamflow against gauge observations. A coupled calibration build
is therefore `-DUSE_MPI=ON -DUSE_MIZUROUTE=ON -DUSE_MODFLOW6=ON`.

There is no bundled end-to-end test of this path: it needs a domain that has a MODFLOW model,
a mizuRoute topology, streamflow observations, and a year of forcing before the calibration
period, and no domain in this repository has all four. The pieces it is built from are
covered separately - `utils/test/test_mflow` exercises the coupling itself, and
`test_calibration_bow.sh` below exercises the calibration machinery.

## `test_calibration_bow.sh`

A short, self-contained calibration on the Bow River at Banff (CAN_05BB001). Everything it needs
is already in the repository, under `utils/test/test_mizuroute/bow_real_data/`: lumped SUMMA
inputs and ten years of forcing, the mizuRoute topology and lumped-to-HRU remapping, and daily
streamflow observations.

```bash
./test_calibration_bow.sh [n_samples] [n_ranks]     # defaults: 4 samples, 2 ranks
```

It writes a TOML configuration pointing at that data, runs a DDS search over five parameters for
one simulated year after a year of spinup, and checks that finite, plausible objective values were
recorded for the trials.

The period and the sample count are deliberately tiny so the test finishes in minutes. It
exercises the machinery; it does **not** produce a calibrated parameter set. A real calibration
uses thousands of samples over several years — see `n_samples` in the generated config.

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
