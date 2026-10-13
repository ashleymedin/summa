#!/usr/bin/env python3
"""Correct SUMMA forcing toward on-site weather stations: air temperature, humidity, shortwave and pressure.

Reads the USGS Benchmark Glacier station records (<stations>/<elevation>/LVL2/*_15min_LVL2_*.csv, as
distributed with the glacier mass-balance data) and the monthly SUMMA forcing files, and writes corrected
copies of the forcing.  Precipitation, wind and longwave are copied unchanged.

  airtemp   K added, per HRU and calendar month: the mean station-minus-forcing difference, with the
            station moved to each HRU's elevation by that month's lapse rate between the two stations
  spechum   multiplied, per HRU: the ratio of mean station to mean forcing specific humidity, the
            station's computed from its relative humidity, temperature and pressure (not lapsed)
  SWRadAtm  multiplied, per HRU: the ratio of summed station to summed forcing incoming shortwave
  airpres   multiplied, per HRU: the ratio of mean station pressure, carried hypsometrically to the HRU
            elevation, to mean forcing pressure

A correction is applied only where the station records the variable; the others stay at 0 or 1.  Station
values are hourly means of the 15-minute samples in the hour ending at each forcing time.  The corrections
are written to corrections.csv beside the corrected files.

    ./correct_forcing_to_station.py $S/domain_Wolverine/forcing/SUMMA_input $S/domain_Wolverine/USGSweatherData \\
        $S/domain_Wolverine/settings/SUMMA/attributes.nc $S/domain_Wolverine/forcing/SUMMA_input/station_corrected
"""

import argparse
import csv
import glob
import os
from datetime import datetime, timedelta, timezone
from zoneinfo import ZoneInfo

import numpy as np
from netCDF4 import Dataset

EPOCH = datetime(1900, 1, 1, tzinfo=timezone.utc)   # forcing time is hours since 1900-01-01 UTC
LOCAL = ZoneInfo("America/Anchorage")
G, R_DRY = 9.80665, 287.05                           # m s-2, J kg-1 K-1
DEFAULT_LAPSE = -0.0065                              # K m-1
COLUMNS = {"site_temp": "T", "RelHum": "RH", "Barom": "P", "ShortwaveIn": "SW"}


def parse_time(text):
    """A station time stamp, either YYYY/MM/DD HH:MM or M/D/YYYY H:MM."""
    for fmt in ("%Y/%m/%d %H:%M", "%m/%d/%Y %H:%M"):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            pass
    raise ValueError(f"unrecognised station time {text!r}")


def read_station(folder):
    """Hourly means {variable: {hour: value}} of one station's 15-minute records; T degC, RH %, P kPa, SW W m-2."""
    sums, counts = {}, {}
    for f in sorted(glob.glob(os.path.join(folder, "LVL2", "*_15min_LVL2_*.csv"))):
        with open(f) as fh:
            rows = csv.reader(fh)
            head = next(rows)
            utc = head.index("UTC_time") if "UTC_time" in head else None
            loc = head.index("local_time")
            cols = {COLUMNS[c]: head.index(c) for c in COLUMNS if c in head}
            last = None
            for r in rows:
                if utc is not None:
                    t = parse_time(r[utc]).replace(tzinfo=timezone.utc)
                else:
                    # local-time-only records: the repeated hour at the end of daylight time takes the later fold
                    naive = parse_time(r[loc])
                    t = naive.replace(tzinfo=LOCAL, fold=0).astimezone(timezone.utc)
                    if last is not None and t <= last:
                        t = naive.replace(tzinfo=LOCAL, fold=1).astimezone(timezone.utc)
                    last = t
                    t = t + timedelta(minutes=15)  # these stamp the start of the sample, the others its end
                hour = int(np.ceil((t - EPOCH).total_seconds() / 3600.0))
                for v, i in cols.items():
                    try:
                        x = float(r[i])
                    except (ValueError, IndexError):
                        continue
                    if not np.isfinite(x):
                        continue
                    sums.setdefault(v, {})
                    counts.setdefault(v, {})
                    sums[v][hour] = sums[v].get(hour, 0.0) + x
                    counts[v][hour] = counts[v].get(hour, 0) + 1
    return {v: {h: s / counts[v][h] for h, s in d.items()} for v, d in sums.items()}


def specific_humidity(rh, t_c, p_kpa):
    """Specific humidity (kg kg-1) from relative humidity (%), temperature (degC) and pressure (kPa), Bolton (1980)."""
    e = rh / 100.0 * 0.6112 * np.exp(17.67 * t_c / (t_c + 243.5))
    return 0.622 * e / (p_kpa - 0.378 * e)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("forcing", help="directory of the monthly SUMMA forcing files")
    p.add_argument("stations", help="USGSweatherData directory, one subdirectory per station elevation (m)")
    p.add_argument("attributes", help="SUMMA attributes file, for the HRU elevations")
    p.add_argument("output", help="directory to write the corrected forcing to")
    a = p.parse_args()

    with Dataset(a.attributes) as d:
        hru_elev = dict(zip(np.array(d["hruId"][:]).astype(int), np.array(d["elevation"][:], float)))
    stations = {int(os.path.basename(s)): read_station(s)
                for s in sorted(glob.glob(os.path.join(a.stations, "*"))) if os.path.basename(s).isdigit()}
    if len(stations) < 2:
        raise SystemExit("two stations are needed for the temperature lapse rate")
    main_z = max(stations, key=lambda z: len(stations[z].get("T", {})))
    upper_z = [z for z in stations if z != main_z][0]
    st_main = stations[main_z]
    sw_z = max(stations, key=lambda z: len(stations[z].get("SW", {})))
    print(f"temperature from {main_z} m, lapse rate against {upper_z} m; shortwave from {sw_z} m "
          f"({len(stations[sw_z].get('SW', {}))} h)")

    # forcing, all files
    files = sorted(glob.glob(os.path.join(a.forcing, "*.nc")))
    hours, frc = [], {v: [] for v in ("airtemp", "spechum", "SWRadAtm", "airpres")}
    for f in files:
        with Dataset(f) as d:
            if "time" not in d.variables or "airtemp" not in d.variables:
                continue
            hours.append(np.array(d["time"][:]).astype(int))
            hru_ids = np.array(d["hruId"][:]).astype(int)
            for v in frc:
                frc[v].append(np.array(d[v][:], float))
    hours = np.concatenate(hours)
    frc = {v: np.concatenate(x) for v, x in frc.items()}
    month = np.array([(EPOCH + timedelta(hours=int(h))).month for h in hours])
    z = np.array([hru_elev[h] for h in hru_ids])

    def on_hours(series):
        return np.array([series.get(h, np.nan) for h in hours])

    # monthly lapse rate between the two stations
    t_main = on_hours(st_main.get("T", {}))
    t_up = on_hours(stations[upper_z].get("T", {}))
    lapse = np.full(13, DEFAULT_LAPSE)
    for m in range(1, 13):
        ok = (month == m) & np.isfinite(t_main) & np.isfinite(t_up)
        if ok.sum() > 24 * 30:
            lapse[m] = np.mean(t_up[ok] - t_main[ok]) / (upper_z - main_z)
    t_hru = t_main[:, None] + lapse[month][:, None] * (z[None, :] - main_z) + 273.15

    # temperature: monthly additive offset per HRU
    offset = np.zeros((13, len(hru_ids)))
    for m in range(1, 13):
        sel = month == m
        for j in range(len(hru_ids)):
            ok = sel & np.isfinite(t_hru[:, j])
            if ok.sum() > 24 * 30:
                offset[m, j] = np.mean(t_hru[ok, j] - frc["airtemp"][ok, j])

    ratio = {v: np.ones(len(hru_ids)) for v in ("spechum", "SWRadAtm", "airpres")}
    # shortwave: ratio of sums, no elevation adjustment
    sw = on_hours(stations[sw_z].get("SW", {}))
    sw = np.where(sw < 0, 0.0, sw)
    for j in range(len(hru_ids)):
        ok = np.isfinite(sw)
        if ok.sum() > 24 * 365:
            ratio["SWRadAtm"][j] = sw[ok].sum() / frc["SWRadAtm"][ok, j].sum()
    # humidity and pressure need the main station's barometer and hygrometer
    rh, pk = on_hours(st_main.get("RH", {})), on_hours(st_main.get("P", {}))
    q = specific_humidity(rh, t_main, pk)
    for j in range(len(hru_ids)):
        ok = np.isfinite(q)
        if ok.sum() > 24 * 365:
            ratio["spechum"][j] = q[ok].mean() / frc["spechum"][ok, j].mean()
        ok = np.isfinite(pk) & np.isfinite(t_hru[:, j])
        if ok.sum() > 24 * 365:
            lap = lapse[month][ok]
            p_hru = pk[ok] * 1000.0 * (t_hru[ok, j] / (t_main[ok] + 273.15)) ** (-G / (R_DRY * lap))
            ratio["airpres"][j] = p_hru.mean() / frc["airpres"][ok, j].mean()

    # write the corrected files
    os.makedirs(a.output, exist_ok=True)
    for f in files:
        with Dataset(f) as s, Dataset(os.path.join(a.output, os.path.basename(f)), "w", format="NETCDF4") as d:
            d.setncatts({k: s.getncattr(k) for k in s.ncattrs()})
            d.history = "corrected toward on-site weather stations by correct_forcing_to_station.py"
            for name, dim in s.dimensions.items():
                d.createDimension(name, None if dim.isunlimited() else len(dim))
            m_file = np.array([(EPOCH + timedelta(hours=int(h))).month for h in np.array(s["time"][:])])
            for name, v in s.variables.items():
                fill = v.getncattr("_FillValue") if "_FillValue" in v.ncattrs() else None
                o = d.createVariable(name, v.dtype, v.dimensions, zlib=True, complevel=4, fill_value=fill)
                o.setncatts({k: v.getncattr(k) for k in v.ncattrs() if k != "_FillValue"})
                x = v[:]
                if name == "airtemp":
                    x = x + offset[m_file]
                elif name in ratio:
                    x = x * ratio[name][None, :]
                o[:] = x

    with open(os.path.join(a.output, "corrections.csv"), "w") as fh:
        w = csv.writer(fh)
        w.writerow(["hruId", "elevation_m", "spechum_ratio", "SWRadAtm_ratio", "airpres_ratio"]
                   + [f"airtemp_offset_K_{m:02d}" for m in range(1, 13)])
        for j, h in enumerate(hru_ids):
            w.writerow([h, z[j], f"{ratio['spechum'][j]:.4f}", f"{ratio['SWRadAtm'][j]:.4f}",
                        f"{ratio['airpres'][j]:.4f}"] + [f"{offset[m, j]:.3f}" for m in range(1, 13)])
        w.writerow([])
        w.writerow(["lapse_rate_K_per_km"] + [f"{1000 * lapse[m]:.2f}" for m in range(1, 13)])

    print("lapse rate (K/km) by month:", " ".join(f"{1000 * lapse[m]:.1f}" for m in range(1, 13)))
    for j, h in enumerate(hru_ids):
        print(f"HRU {h} ({z[j]:.0f} m): q x{ratio['spechum'][j]:.3f}, SW x{ratio['SWRadAtm'][j]:.3f}, "
              f"P x{ratio['airpres'][j]:.3f}, T offset by month (K) "
              + " ".join(f"{offset[m, j]:+.1f}" for m in range(1, 13)))
    print(f"{len(files)} files written to {a.output}")


if __name__ == "__main__":
    main()
