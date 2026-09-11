# domain_sagehen4 — four HRUs, lateral flow, coupled to MODFLOW 6

A synthetic 4-HRU case built from `domain_sagehen1` (one month of hourly forcing tiled
across four columns), run through the MODFLOW 6 coupler for 72 hourly steps, the same as
`ex-gwf-sagehen`. It does two jobs: it exercises `groundwatr = modLatflow`, and it guards
the HRU cascade ordering in `run_oneGRU`.

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

An **exact-equality** check on `run1_latflow`, so the size of the difference does not matter:

    hru 1 (receiver    ): inflow 4.0458e-02   outflow 3.53643928e-03
    hru 4 (CONTROL twin): inflow 0.0000e+00   outflow 3.49913476e-03   (ratio 1.0107)

**HRU 1 and HRU 4 must not be equal.** If they are bit-identical while HRU 1's inflow is
non-zero, the cascade ordering has regressed and the inflow is being discarded. Verified both
ways: with the ordering fix reverted, the two are bit-identical at `3.49912619e-03`.

Comparing HRU 1 against HRU 2 or 3 does not work — they have a different area and band, so
they differ for unrelated reasons. Only HRU 4 is a true control, which is why it shares
HRU 1's cells and elevation rather than owning a private region.

In `run1_noLatflow` every column outflow is zero, which is the expected contrast rather than
a second cascade check.
