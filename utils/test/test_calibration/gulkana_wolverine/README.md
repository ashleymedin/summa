# gulkana_wolverine — a glacierized domain with GRACE storage observations

Two Alaskan glacier basins, Gulkana and Wolverine, as one SUMMA domain of 2 GRUs and 6 HRUs, with
the GRACE terrestrial-water-storage record for each. It is here to give the storage targets
something real to calibrate against: a glacier basin loses mass over the record, so basin storage
carries a signal a discharge series does not.

| | |
| --- | --- |
| Forcing | Oct 2008 – Dec 2009, 15 monthly files, double precision |
| Domain | 2 GRUs, 6 HRUs; GRU 1 = 26.6 km², GRU 2 = 24.7 km² |
| Observations | GRACE JPL mass anomaly per basin, mm equivalent water thickness, 2002-04 to 2025-03 |

The domain comes from `summa_test_cases/test_cases` (`settings/multiGruTestCases/gulkana_wolverine`
and `input_data/gulkana_wolverine`), taking the double-precision forcing only. The observations were
built from the SYMFLUENCE basin-average GRACE records with
`utils/pre-processing/grace_tws_to_netcdf.py`, which writes them in the same shape as the bundled
streamflow observations so that one reader serves both.

The forcing period sits inside the GRACE record, which is the point: the Bow domain bundled for the
streamflow test runs 1980–1990 and GRACE begins in 2002, so the two never overlap. This domain is
where a storage target can actually be exercised today.

## What still has to be built before a GRACE target runs on this

**A target has to be able to name which GRU it scores.** Simulated series are currently collected as
an area-weighted mean over the whole domain, which is right for a single-basin domain and wrong
here. The two basins are nearly equal in area but their storage signals are not comparable:

| Basin | GRACE JPL anomaly range |
| --- | --- |
| Gulkana | −105.8 to +45.0 mm |
| Wolverine | −506.9 to +212.6 mm |

A domain mean of those two would match neither record, so scoring it against either basin's GRACE
would be meaningless. Spatial selection is the missing feature, and the stream-temperature target
needs the same thing for a different reason — a reach outlet temperature is a property of one reach,
not an average over the basin.

Until that exists, this directory is the data, not a runnable calibration. A single-basin domain
would sidestep it, at the cost of the 2-GRU coverage.

## Spinning it up

15 months is enough for a short spinup and about a year of calibration, giving roughly twelve
monthly GRACE comparisons. That is thin for a calibration and fine for a test.

The GRACE anomalies are relative to their product's own baseline (2004–2009 for the JPL mascons),
which a one-year run cannot reproduce. Reference both series to the run period instead, by setting
the target's `baseline_start` and `baseline_end` to the calibration period: each series is then a
departure from its own mean over the same months, which is what makes them comparable.

No mizuRoute topology is bundled here. A storage-only calibration does not need one — the simulation
time axis is filled from the forcing rather than from the routing — but a streamflow target would.
