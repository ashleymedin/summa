# gulkana_wolverine — a glacierized domain with GRACE storage observations

Two Alaskan glacier basins, Gulkana and Wolverine, as one SUMMA domain of 2 GRUs and 6 HRUs, with
the GRACE terrestrial-water-storage record for each. It is here to give the storage targets
something real to calibrate against: a glacier basin loses mass over the record, so basin storage
carries a signal a discharge series does not.

| | |
| --- | --- |
| Forcing | Oct 2008 – Dec 2009, 15 monthly files, double precision |
| Domain | 2 GRUs, 6 HRUs; GRU 1 = 26.6 km², GRU 2 = 24.7 km² |
| Observations | GRACE JPL mass anomaly per basin, mm equivalent water thickness, 2002-04 to 2025-11 |

The domain comes from `summa_test_cases/test_cases` (`settings/multiGruTestCases/gulkana_wolverine`
and `input_data/gulkana_wolverine`), taking the double-precision forcing only. The observations are the
JPL RL06.3 CRI mascon cell holding each GRU's centroid, written by
`utils/pre-processing/acquire_grace_tws.py` in the same form as the bundled streamflow observations,
so that one reader serves both:

```bash
cd utils/pre-processing
./acquire_grace_tws.py --attributes ../test/test_calibration/gulkana_wolverine/summa_inputs/attributes.nc \
                       --gru 1 ../test/test_calibration/gulkana_wolverine/observations/gulkana_grace_tws.nc
```

and `--gru 2` for Wolverine.

The forcing period sits inside the GRACE record, which is the point: the Bow domain bundled for the
streamflow test runs 1980–1990 and GRACE begins in 2002, so the two never overlap. This domain is
where a storage target can actually be exercised today.

## Scoring each basin against its own record

The two basins are nearly equal in area but their storage signals are not comparable:

| Basin | GRU | GRACE JPL anomaly range |
| --- | --- | --- |
| Gulkana | 1 | −997.8 to +519.9 mm |
| Wolverine | 2 | −4985.2 to +2375.0 mm |

A domain mean of the two would match neither record, so each target names its basin with `gru = 1`
or `gru = 2`. Over one Oct–Dec 2009 run, the Gulkana record scores 54 mm RMSE against GRU 1, 385 mm
against GRU 2, and 194 mm against the domain mean.

## Spinning it up

The calibration driver spins up for the year before `start_time`, and the forcing starts on
2008-10-01, so a calibration here starts no earlier than 2009-10-01 and covers at most three months, or
three monthly GRACE comparisons. That is enough to exercise the targets, not to calibrate anything.
Set `forcing_path` to the `forcing/` directory.

The GRACE anomalies are relative to their product's own baseline (2004–2009 for the JPL mascons),
which a one-year run cannot reproduce. Reference both series to the run period instead, by setting
the target's `baseline_start` and `baseline_end` to the calibration period: each series is then a
departure from its own mean over the same months, which is what makes them comparable.

No mizuRoute topology is bundled here. A storage-only calibration does not need one — the simulation
time axis is filled from the forcing rather than from the routing — but a streamflow target would.
