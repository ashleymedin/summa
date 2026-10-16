# gulkana_wolverine — two glacier basins as one domain of two GRUs

Gulkana and Wolverine glaciers, Alaska, as one SUMMA domain of 2 GRUs and 6 HRUs, with each basin's
glacier mass balance and GRACE storage. It is the minimal test of a multi-objective calibration over
several GRUs (`../test_calibration_gulkana_wolverine.sh`): each target names its basin by its GRU.

| | |
| --- | --- |
| Forcing | ERA5, Oct 2008 – Dec 2011, 39 monthly files, double precision |
| Domain | 2 GRUs, 6 HRUs; GRU 1 = Gulkana (HRUs 1–3, 26.6 km²), GRU 2 = Wolverine (HRUs 4–6, 24.7 km²) |
| Mass balance | USGS Benchmark Glacier glacier-wide seasonal and annual balances, mm w.e. |
| Storage | GRACE JPL mass anomaly per basin, mm equivalent water thickness, 2002-04 to 2025-11 |

The domain is written by `build_cases.py cases.yaml gulkana_wolverine` in `summa_glacier_workflow`, beside this
checkout, which also builds the Gulkana and Wolverine glacier cases: the two glacier domains side by
side, their HRUs, glacier grids and forcing concatenated, Wolverine's HRUs renumbered 4–6. The cold
state starts the soil and debris 1.5 °C above each HRU's mean air temperature over its first year,
October 2008 to September 2009. The text settings and `trialParams.nc` are this folder's own.

The mass balances were written by `acquire_glacier_mass_balance.py` (in `summa_glacier_workflow`) and the GRACE
records by `utils/pre-processing/acquire_grace_tws.py`, the JPL RL06.3 CRI mascon cell holding each
GRU's centroid:

```bash
cd utils/pre-processing
./acquire_grace_tws.py --attributes ../test/test_calibration/gulkana_wolverine/summa_inputs/attributes.nc \
                       --gru 1 ../test/test_calibration/gulkana_wolverine/observations/gulkana_grace_tws.nc
```

and `--gru 2` for Wolverine.

## Scoring each basin against its own record

The two basins are nearly equal in area but their storage signals are not comparable:

| Basin | GRU | GRACE JPL anomaly range |
| --- | --- | --- |
| Gulkana | 1 | −997.8 to +519.9 mm |
| Wolverine | 2 | −4985.2 to +2375.0 mm |

A domain mean of the two would match neither record, so each target names its basin with `gru = 1`
or `gru = 2`.

## Spinning it up

The calibration driver spins up for the year before `start_time`, and the forcing starts on
2008-10-01, so a calibration here starts no earlier than 2009-10-01. A seasonal balance is measured
from the opposite extreme a year before its observation, so the first balances inside the run are
2011's: May 2011 (winter) and September 2011 (summer).

A GRACE anomaly is relative to its product's own baseline (2004–2009 for the JPL mascons), which a short
run cannot reproduce. Reference both series to the run period instead, by setting the target's
`baseline_start` and `baseline_end` to the calibration period.

No mizuRoute topology is bundled here. A storage or mass-balance calibration does not need one, but a
streamflow target would.
