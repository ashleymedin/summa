# wolverine — a glacier basin with a routed reach, scored on discharge, stream temperature and mass balance

Wolverine Glacier, Alaska, as one GRU of three land HRUs (24.7 km², 16.4 km² of it glacier) and one
stream HRU for the reach from the terminus to USGS 15236900, Wolverine Creek near Lawing. The gauge
records both discharge and water temperature, so a routed calibration can score all three records.

| | |
| --- | --- |
| Forcing | ERA5, Jan 2016 – Dec 2019, 48 monthly files, UTC |
| Domain | 1 GRU, 3 land HRUs + 1 stream HRU, gridded glacier; one mizuRoute reach (`topology.nc`) |
| Discharge | USGS 15236900, daily mean, m3 s-1, open-water season only |
| Water temperature | USGS 15236900, daily mean, degC, open-water season only: 144–209 days a year |
| Mass balance | USGS Benchmark Glacier glacier-wide seasonal (winter and summer) balance, mm w.e. |

The land HRUs, their parameters and forcing are the Wolverine glacier setup. The text settings are
the bundled Gulkana ones, with Wolverine's own `basinParamInfo.txt`. 
`utils/pre-processing/add_stream_hru.py` appends the stream HRU and writes the one-reach topology:

```bash
cd utils/pre-processing
S=~/Research/Symfluence/SYMFLUENCE_data/domain_Wolverine
./add_stream_hru.py $S/settings/SUMMA $S/forcing/SUMMA_input ../test/test_calibration/wolverine \
    --attributes attributes_glac.nc --cold-state coldstate_glac.nc \
    --first-month 2016-01 --last-month 2019-12 \
    --length 2000 --slope 0.03 --width 5 --elevation 390 --latitude 60.374 --longitude -148.912
./acquire_streamflow.py usgs 15236900 2016-01-01 2019-12-31 \
    ../test/test_calibration/wolverine/observations/USGS_15236900_daily_flow.nc --model-utc-offset 0
./acquire_stream_temperature.py 15236900 2016-01-01 2019-12-31 \
    ../test/test_calibration/wolverine/observations/USGS_15236900_daily_temperature.nc --model-utc-offset 0
./acquire_glacier_mass_balance.py $S/observations/storage/mass_balance/Output_Wolverine_Glacier_Wide_solutions_calibrated.csv \
    seasonal ../test/test_calibration/wolverine/observations/wolverine_mass_balance_seasonal.nc
```

The reach geometry is an estimate, since the domain has no river network: 2 km from the terminus to
the gauge (357 m), falling 60 m, and 5 m wide, the width mizuRoute's `wscale` gives this basin area
and within the 15–30 ft single-thread channel at the gauge. The stream HRU takes the forcing and soil
of HRU 1, the lowest land HRU.

`modelDecisions.txt` differs from Gulkana's in two decisions. `deepTherml = bedrockLyrs` puts 30 m
of thermal-only bedrock under the soil, and the groundwater reaching the channel leaves at its base;
with the default, it leaves at the base of the 4 m soil column, near freezing, and the water entering
the reach in July is 0.07 °C against 1.37 °C at the gauge. The cold state starts each land HRU's
soil 1.5 °C above its 2016–2019 mean annual air temperature (2.5, 1.2 and −0.7 °C), and the bedrock
starts from the deepest soil layer. Its base moves less than 0.2 °C by 2019, so that start sets the
groundwater temperature.
`hyporhTdyn = proxy` returns part of the reach flow at the temperature it had earlier, which damps the
daily swing.

## Running it

`test_calibration_wolverine.sh` spins up over 2016 and scores 2017–2019. The water temperature record
starts in April 2016, so this is the earliest run that scores it after a full year of spin-up. The
gauge is the reach outlet, so discharge and temperature both name `reach = 1`. The mass balance names
`gru = 1`. The 2017 balances need an extreme from before the run and are left out, so the winter and
summer balances of 2018–2019 are scored.

The glacier-fed creek runs at 0.0–2.6 °C over the scored years, so a temperature RMSE of a degree is
a large miss.
