#!/usr/bin/env python3
"""Write synthetic calibration targets from a reference trial of the calibration driver.

The reference trial writes hourly SUMMA/mizuRoute output (write_timeseries, write_Qreach) and MODFLOW's
OBS6 CSV.  Each target is taken as the driver would simulate it, then as daily means over UTC days,
stamped at the end of the day, the way the driver aligns an hourly series to daily observations:

    q_obs         discharge leaving the reach, m3 s-1
    t_obs         water temperature leaving the reach, degC
    h_<name>      head at OBS6 observation <name>, m
    storage_obs   basin__StorageChange over the domain, weighted by GRU area and accumulated, mm
    baseflow_obs  basin__AquiferBaseflow of one GRU, m s-1
    aquifer_obs   scalarAquiferStorage of one HRU, m

Usage:
    make_targets.py <timestep output .nc> <OBS6 .csv> <attributes.nc> <reach> <gru> <hru> <output .nc>
"""

import datetime as dt
import sys

import numpy as np
from netCDF4 import Dataset, num2date

EPOCH = dt.datetime(1950, 1, 1)


def daily_means(times, values):
    """Means over UTC days of the values stamped in (midnight, next midnight]; whole days only."""
    days = {}
    for t, v in zip(times, values):
        days.setdefault((t - dt.timedelta(microseconds=1)).date(), []).append(v)
    per_day = max(len(v) for v in days.values())
    return {d: float(np.mean(v)) for d, v in days.items() if len(v) == per_day}


def main(timestep_nc, obs_csv, attributes, reach, gru, hru, output):
    reach, gru, hru = int(reach), int(gru), int(hru)
    with Dataset(attributes) as a:
        hru2gru = np.asarray(a["hru2gruId"][:])
        gru_ids = np.asarray(a["gruId"][:])
        gru_area = np.array([np.asarray(a["HRUarea"][:])[hru2gru == g].sum() for g in gru_ids])

    with Dataset(timestep_nc) as d:
        tv = d["time"]
        times = [dt.datetime(*t.timetuple()[:6]) for t in num2date(tv[:], tv.units, getattr(tv, "calendar", "standard"))]
        iseg = int(np.where(np.asarray(d["seg"][:]) == reach)[0][0])
        igru = int(np.where(np.asarray(d["gruId"][:]) == gru)[0][0])
        ihru = int(np.where(np.asarray(d["hruId"][:]) == hru)[0][0])
        q = np.asarray(d["Q_reach"][:])[:, iseg, 0]
        temp = np.asarray(d["T_reach"][:])[:, iseg] - 273.15
        dS = np.asarray(d["basin__StorageChange"][:]) @ gru_area / gru_area.sum()
        baseflow = np.asarray(d["basin__AquiferBaseflow"][:])[:, igru]
        aquifer = np.asarray(d["scalarAquiferStorage"][:])[:, ihru, 0]   # the HRU's one domain; the rest is fill
        units = {v: d[v].units for v in ("basin__AquiferBaseflow", "scalarAquiferStorage")}

    # the driver integrates a rate step by step from the first, which adds nothing
    storage = np.concatenate([[0.0], np.cumsum(dS[1:] * 3600.0)])

    names = open(obs_csv).readline().strip().split(",")[1:]
    rows = np.loadtxt(obs_csv, delimiter=",", skiprows=1, ndmin=2)[-len(times):, 1:]
    if len(rows) != len(times):
        sys.exit(f"{obs_csv} holds {len(rows)} rows for {len(times)} time steps")

    series = {"q_obs": (q, "m3 s-1"), "t_obs": (temp, "degC"), "storage_obs": (storage, "mm"),
              "baseflow_obs": (baseflow, units["basin__AquiferBaseflow"]),
              "aquifer_obs": (aquifer, units["scalarAquiferStorage"])}
    for k, name in enumerate(names):
        series[f"h_{name.lower()}"] = (rows[:, k], "m")

    means = {v: daily_means(times, s) for v, (s, _) in series.items()}
    days = sorted(means["q_obs"])
    with Dataset(output, "w", format="NETCDF4") as o:
        o.createDimension("time", len(days))
        end = [int((dt.datetime.combine(d, dt.time()) + dt.timedelta(days=1) - EPOCH).total_seconds() // 60) for d in days]
        t = o.createVariable("time", "i8", ("time",))
        t.standard_name, t.units, t.calendar = "time", "minutes since 1950-01-01", "standard"
        t.time_zone, t.comment = "UTC+0", "end of the day each value covers"
        t[:] = end
        for vname, (_, u) in series.items():
            v = o.createVariable(vname, "f8", ("time",), fill_value=np.nan)
            v.units, v.cell_methods = u, "time: mean"
            v[:] = [means[vname].get(d, np.nan) for d in days]
        o.title = "Synthetic Sagehen calibration targets"
        o.source = f"reference trial {timestep_nc}, {obs_csv}; reach {reach}, gru {gru}, hru {hru}"
    print(f"{output}: {len(days)} days of {', '.join(series)}")


if __name__ == "__main__":
    if len(sys.argv) != 8:
        sys.exit(__doc__)
    main(*sys.argv[1:])
