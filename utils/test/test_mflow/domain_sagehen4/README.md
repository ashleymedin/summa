# domain_sagehen4 — four HRUs, lateral flow, coupled to MODFLOW 6

A synthetic 4-HRU case built from `domain_sagehen1` (one month of hourly forcing tiled
across four columns), run through the MODFLOW 6 coupler for 72 hourly steps, the same as
`ex-gwf-sagehen`. The forcing is the February 2017 wet period, 2017-02-07 to 2017-02-09,
written by `../tools/make_wet_forcing.py`. It does two jobs: it exercises
`groundwatr = modLatflow`, and it guards the HRU cascade ordering in `run_oneGRU`.

    hruId          1       2       3       4
    downHRUindex   0       1       1       0
    HRUarea       10%     40%     40%     10%    of the basin
    MODFLOW      lower   upper   upper   lower   elevation band
                 shared  private private shared

HRUs 2 and 3 drain into HRU 1. HRU 4 is an exact twin of HRU 1 — same cells, elevation and
area — that receives nothing, and serves as the control.

HRUs 2 and 3 each own a **private, disjoint** region of the grid; HRUs 1 and 4 share the
lower band because they must be exact twins for the check below to mean anything. Disjoint
and shared mappings both work, and an HRU may own as little as one cell.

## Why the ordering matters

The receiving HRU deliberately sits at index 1, i.e. **before** its contributors in the HRU
loop. That is the ordering that used to lose water: the outflow of HRUs 2 and 3 was added to
HRU 1's `mLayerColumnInflow` only after HRU 1 had already been solved, and the inflow array
is zeroed at the start of the next step, so the water was never applied to anything. Worse,
`mLayerColumnInflow` is written to output *after* the loop finishes accumulating, so it
reported the full inflow and hid the loss.

`run_oneGRU` now orders HRUs upslope-before-downslope (Kahn's algorithm over
`downHRUindex`) and stops on a cyclic network.

## Running

    ./run_sagehen4.sh              # with lateral flow    (groundwatr = modLatflow)
    ./run_sagehen4_noLatflow.sh    # without              (groundwatr = modflow)

Needs `srcextern/summa/bin/summa_modflow6.exe`, built with `-DUSE_MODFLOW6=ON`. The case
shares the `ex-gwf-sagehen` MODFLOW model and supplies its own coupler config via
`coupler_commands.sh -c`, which carries this case's `map_file`.

The two file managers point at decision files that differ on the `groundwatr` line **alone**,
so running the pair isolates what lateral flow contributes. They write `run1_latflow*` and
`run1_noLatflow*`, so neither overwrites the other.

## What to check

**This check is currently inert, and the case cannot detect a cascade regression.** Both
`run1_latflow` and `run1_noLatflow` report every `mLayerColumnOutflow` and
`mLayerColumnInflow` as exactly zero, so HRU 1 and its control twin HRU 4 are
bit-identical for want of any lateral water at all.

The reason is not the forcing. The case was moved from August 2019 -- the driest month in
the record, 3.8 mm -- to the February 2017 atmospheric river, 128 mm of precipitation in
the 72 hours run and 121 mm of it rain, and the outflow stayed at zero. The soil ends the
run wet but nowhere saturated:

    liquid fraction   0.3701 0.3700 0.3684 0.3606 0.3432 0.2723 0.2490 0.3392
    matric head (m)  -0.112 -0.113 -0.121 -0.159 -0.246 -0.732 -0.993 -0.267
                     fieldCapacity 0.20, theta_sat 0.55

`groundwatr.f90` zeroes the transmissivity of every layer above the saturated zone
(`if (ixSaturation>1) trSoil(1:ixSaturation-1) = 0`), so with no saturated layer there is
no lateral flow to cap or route. The drainable water is ample, well above field capacity,
so the outflow cap is not what closes this off.

What is odd, and is the thing to chase: with `bcLowrSoiH = presHead` the coupler sets
`lowerBoundHead = h_mf6 - (z_surface - soil_thickness)`, which for a water table 1.5 m below
the surface and a 4 m column is about +2.5 m -- the column base ought to be saturated. It is
not, the upward flux being conductivity-limited at roughly 1.8 mm/h. Until that is resolved
`groundwatr = modLatflow` moves no water in a coupled run, for any forcing.

For the record, the values this check reported before the lateral outflow cap of
`changes_fromV3Summa` entry 71:

    hru 1 (receiver    ): inflow 4.0458e-02   outflow 3.53643928e-03
    hru 4 (CONTROL twin): inflow 0.0000e+00   outflow 3.49913476e-03   (ratio 1.0107)

**HRU 1 and HRU 4 must not be equal** once lateral flow works again. If they are
bit-identical while HRU 1's inflow is non-zero, the cascade ordering has regressed and the
inflow is being discarded. Verified both ways at the time: with the ordering fix reverted,
the two were bit-identical at `3.49912619e-03`.

Comparing HRU 1 against HRU 2 or 3 does not work -- they have a different area and band, so
they differ for unrelated reasons. Only HRU 4 is a true control, which is why it shares
HRU 1's cells and elevation rather than owning a private region.

In `run1_noLatflow` every column outflow is zero, which is the expected contrast rather than
a second cascade check.
