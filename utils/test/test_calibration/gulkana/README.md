# gulkana — a glacier basin with discharge, GRACE and mass-balance observations

Gulkana Glacier, Alaska, as one GRU of three HRUs (26.6 km², 17.3 km² of it glacier), with the three
records a glacier calibration can be scored against: discharge below the glacier, basin water storage
from GRACE, and the glacier's own mass balance.

| | |
| --- | --- |
| Forcing | ERA5, Jan 2008 – Dec 2012, 60 monthly files, UTC |
| Domain | 1 GRU, 3 HRUs, gridded glacier (`attributes.nc`, `coldState.nc`) |
| Discharge | USGS 15478040, Phelan Creek near Paxson, daily mean, m3 s-1 |
| Storage | GRACE JPL RL06.3 mascon cell holding the basin centroid, mm |
| Mass balance | USGS Benchmark Glacier glacier-wide seasonal (winter and summer) balance, mm w.e. |

The observations were written by the scripts in `utils/pre-processing/`:

```bash
cd utils/pre-processing
./acquire_streamflow.py usgs 15478040 2008-01-01 2012-12-31 \
    ../test/test_calibration/gulkana/observations/USGS_15478040_daily_flow.nc --model-utc-offset 0
./acquire_glacier_mass_balance.py Output_Gulkana_Glacier_Wide_solutions_calibrated.csv seasonal \
    ../test/test_calibration/gulkana/observations/gulkana_mass_balance_seasonal.nc
```

The mass-balance table is from the USGS data release at
https://www.sciencebase.gov/catalog/item/6441de03d34ee8d4ade7d2a5, downloaded by hand. The GRACE file
is the same cell as `gulkana_wolverine/observations/gulkana_grace_tws.nc`.

## Running it

The calibration spins up over the year before `start_time`, so start on 2009-01-01 and score
2009–2012. No mizuRoute topology is bundled: discharge is scored as `streamflow` over `gru = 1`, the
GRU's routed runoff times its area. The 2009 balances need an extreme from before the run and are left
out, so the winter and summer balances of 2010–2012 are scored. One evaluation takes about 100 s.

With the default parameters the discharge KGE is 0.46 and the seasonal balance scores 228 mm RMSE.
Score the mass balance or GRACE with discharge, not both. Against GRACE the RMSE is 1173 mm, nearly all
of it trend: the glacier loses 1.4–1.8 m w.e. a year over 64% of the basin, so the simulated basin
storage falls about 1000 mm a year, much faster than the GRACE cell does. Scored together, the
parameters that fit GRACE best fit the mass balance worst.
