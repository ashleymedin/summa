# domain_cascade — HRU cascade ordering test, coupled to MODFLOW 6

A synthetic 4-HRU case that guards the lateral-flow cascade in `run_oneGRU`, run through the
MODFLOW 6 coupler with `groundwatr = modLatflow`.  72 hourly steps, the same as
`ex-gwf-sagehen`.

    hruId        1      2      3      4
    downHRUindex 0      1      1      0
    HRUarea     10%    40%    40%    10%     of the basin
    MODFLOW     lower  upper  upper  lower   elevation band

HRUs 2 and 3 drain into HRU 1.  HRU 4 is an exact twin of HRU 1 — same area, same cells,
same elevation — that receives nothing, and is the control.

The receiving HRU deliberately sits at index 1, i.e. **before** its contributors in the HRU
loop.  That is the ordering that used to lose water: the outflow of HRUs 2 and 3 was added
to HRU 1's `mLayerColumnInflow` only after HRU 1 had already been solved, and the inflow
array is zeroed at the start of the next step, so the water was never applied to anything.
`mLayerColumnInflow` is written to output after the loop has finished accumulating, so it
reported the full inflow and hid the loss.  `run_oneGRU` now orders HRUs
upslope-before-downslope (Kahn's algorithm over `downHRUindex`) and errors on a cyclic
network.

## Running

    ./domain_cascade/run.sh

Needs `srcextern/summa/bin/summa_modflow6.exe` (build with `-DUSE_MODFLOW6=ON`).  The case
passes its own coupler config with `coupler_commands.sh -c`, so it keeps its own
`map_file` while `mf6/` stays a plain copy of `ex-gwf-sagehen`.

## What to check

An **exact-equality** check, so the size of the difference does not matter:

    hru 1 (receiver    ): inflow 2.3509e-02   outflow 1.97601308e-03
    hru 4 (CONTROL twin): inflow 0.0000e+00   outflow 1.95963223e-03   (ratio 1.0084)

**HRU 1 and HRU 4 must not be equal.**  If they are bit-identical while HRU 1's inflow is
non-zero, the cascade ordering has regressed and the inflow is being discarded.  Comparing
against HRU 2 or 3 does not work — they have a different area and band, so they differ from
HRU 1 for unrelated reasons.  Only HRU 4 is a true control.

## Coupled-model stability, worth knowing

HRUs may cover part of the grid — this case does, two elevation bands — but **giving every
HRU its own private set of cells makes the run fail**.  That configuration was tried: four
disjoint bands, each HRU's elevation correctly set to the mean land surface of its own
cells, giving a perfectly sane `lowerBoundHead` of 2.5 m (water table 1.5 m below the
surface in a 4 m column).  SUMMA fails at the very first step, sub-stepping below 1 s.

It is not the head value.  The same 2.5 m head runs fine both with a whole-grid map and
with this two-region map.  What breaks is the feedback: with private cells, an HRU's
recharge lands only on the cells whose head it then reads back, so the explicit one-step
head/recharge loop is tight and undamped.  Sharing cells between HRUs averages over
contributors and damps it.

That is a property of the explicit coupling, not of this test case, and it will bite any
real application that maps one HRU per MODFLOW sub-area.  Mitigations would be relaxation
on the head feedback, a shorter coupling step, or iterating the exchange to convergence
within a step — none of which are implemented.
