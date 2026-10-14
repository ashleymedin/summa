#!/usr/bin/env python3
"""Correct SUMMA forcing toward on-site weather stations and, optionally, glacier stake winter balances.

Reads the USGS Benchmark Glacier station records (<stations>/<elevation>/LVL2/*_15min_LVL2_*.csv, as
distributed with the glacier mass-balance data) and the monthly SUMMA forcing files, and writes corrected
copies of the forcing.  Wind and longwave are copied unchanged.

  airtemp   K added, per HRU, calendar month and forcing hour type (wet or dry, by the forcing's own
            precipitation against --wet-threshold): the mean station-minus-forcing difference, with the
            station moved to each HRU's elevation by that month's lapse rate between the two stations
  spechum   multiplied, per HRU, calendar month and hour type: the ratio of mean station to mean forcing
            specific humidity, the station's computed from its relative humidity, temperature and pressure (not lapsed).
            Where the station has no humidity, --spechum-ratios supplies them (spechum_ratio_from_campaign.py)
  SWRadAtm  multiplied, per HRU: the ratio of summed station to summed forcing incoming shortwave
  airpres   multiplied, per HRU: the ratio of mean station pressure, carried hypsometrically to the HRU
            elevation, to mean forcing pressure
  pptrate   multiplied, per HRU, only with --stakes: measured winter balance over the corrected forcing's
            snowfall at the stakes within --stake-window of the HRU (or the nearest such group), over each
            stake's fall-to-spring interval.  Snowfall is split on wet-bulb temperature as SUMMA does, with
            --temp-crit-rain and --temp-range.  Only stakes above --stake-min-elev count, where winter
            ablation is small.  With --gauge-shield, every HRU the stakes leave out instead takes the ratio
            of the station's weighing-gauge catch, corrected for wind undercatch, to the forcing's
            precipitation at the gauge's elevation, summed over the full hours the gauge, wind and
            temperature all record.  The catch efficiency is Eq. 3 of Kochendorfer et al. (2017, HESS 21,
            3525) for gauge-height wind (capped at 7.2 m s-1) and air temperature; without --gauge-shield
            those HRUs keep their precipitation.

A correction is applied only where the data exist; the others stay at 0 or 1.  Station values are hourly
means of the 15-minute samples in the hour ending at each forcing time, precipitation the sum of all four.  The corrections are written to
corrections.csv beside the corrected files.

    ./correct_forcing_to_station.py $S/domain_Wolverine/forcing/SUMMA_input $S/domain_Wolverine/USGSweatherData \\
        $S/domain_Wolverine/settings/SUMMA/attributes.nc $S/domain_Wolverine/forcing/SUMMA_input/station_corrected \\
        --stakes $S/domain_Wolverine/observations/storage/mass_balance/Input_Wolverine_Glaciological_Data.csv
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
COLUMNS = {"site_temp": "T", "RelHum": "RH", "Barom": "P", "ShortwaveIn": "SW", "WindSpeed": "U",
           "Precip_Weighing_Incremental": "PR"}
CATCH = {"SA": (0.0348, 1.366, 0.779), "UN": (0.0785, 0.729, 0.407)}   # Kochendorfer et al. (2017) Table 2, gauge height
U_CATCH_MAX = 7.2                                                      # m s-1, where those fits end


def parse_time(text):
    """A station time stamp, either YYYY/MM/DD HH:MM or M/D/YYYY H:MM."""
    for fmt in ("%Y/%m/%d %H:%M", "%m/%d/%Y %H:%M"):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            pass
    raise ValueError(f"unrecognised station time {text!r}")


def read_station(folder):
    """Hourly means {variable: {hour: value}} of one station's 15-minute records; T degC, RH %, P kPa, SW W m-2,
    U m s-1; PR is the hour's total in mm, kept only where all four samples exist."""
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
    out = {v: {h: s / counts[v][h] for h, s in d.items()} for v, d in sums.items() if v != "PR"}
    if "PR" in sums:
        out["PR"] = {h: s for h, s in sums["PR"].items() if counts["PR"][h] == 4}
    return out


def catch_efficiency(u, t_c, shield):
    """Gauge catch efficiency (-) from gauge-height wind (m s-1) and air temperature (degC), Kochendorfer et al. (2017) Eq. 3."""
    a, b, c = CATCH[shield]
    return np.exp(-a * np.minimum(u, U_CATCH_MAX) * (1.0 - np.arctan(b * t_c) + c))


def specific_humidity(rh, t_c, p_kpa):
    """Specific humidity (kg kg-1) from relative humidity (%), temperature (degC) and pressure (kPa), Bolton (1980)."""
    e = rh / 100.0 * 0.6112 * np.exp(17.67 * t_c / (t_c + 243.5))
    return 0.622 * e / (p_kpa - 0.378 * e)


def wet_bulb(t_k, q, p_pa):
    """Wet-bulb temperature (K) from air temperature (K), specific humidity and pressure (Pa), Stull (2011)."""
    t = t_k - 273.15
    e = q * p_pa / (0.622 + 0.378 * q)
    rh = np.clip(100.0 * e / (611.2 * np.exp(17.67 * t / (t + 243.5))), 1.0, 100.0)
    return (t * np.arctan(0.151977 * np.sqrt(rh + 8.313659)) + np.arctan(t + rh) - np.arctan(rh - 1.676331)
            + 0.00391838 * rh ** 1.5 * np.arctan(0.023101 * rh) - 4.686035 + 273.15)


def along_elevation(z_hru, values, z):
    """values (time, hru) interpolated linearly in elevation to z, extrapolated from the end pair."""
    k = np.argsort(z_hru)
    zs, v = z_hru[k], values[:, k]
    i = int(np.clip(np.searchsorted(zs, z) - 1, 0, len(zs) - 2))
    w = (z - zs[i]) / (zs[i + 1] - zs[i])
    return v[:, i] + (v[:, i + 1] - v[:, i]) * w


def stake_ratios(path, hours, z, snow, min_elev, window):
    """Per HRU, summed stake winter balance over summed forcing snowfall (both mm) at the stakes near it."""
    rows = list(csv.DictReader(open(path)))
    fall = {(int(r["Year"]), r["site_name"]): r["fall_date"] for r in rows}
    stakes = []
    for r in rows:
        if r["bw"] in ("nan", "") or float(r["elevation"]) < min_elev:
            continue
        y = int(r["Year"])
        f0 = fall.get((y - 1, r["site_name"]), "nan")
        start = datetime.strptime(f0, "%Y/%m/%d") if f0 != "nan" else datetime(y - 1, 9, 30)
        end = datetime.strptime(r["spring_date"], "%Y/%m/%d")
        h0 = (start.replace(tzinfo=timezone.utc) - EPOCH).total_seconds() / 3600.0
        h1 = (end.replace(tzinfo=timezone.utc) - EPOCH).total_seconds() / 3600.0
        sel = (hours > h0) & (hours <= h1)
        if sel.sum() < 0.9 * (h1 - h0):
            continue
        ze = float(r["elevation"])
        stakes.append((ze, 1000.0 * float(r["bw"]), along_elevation(z, snow[sel], ze).sum()))
    ratio, used = np.ones(len(z)), [[] for _ in z]
    if not stakes:
        return ratio, used
    ze = np.array([s[0] for s in stakes])
    for j, zj in enumerate(z):
        if zj < min_elev:
            continue
        near = np.abs(ze - zj) <= window
        if not near.any():
            near = np.abs(ze - ze[np.argmin(np.abs(ze - zj))]) <= window
        used[j] = [s for s, n in zip(stakes, near) if n]
        ratio[j] = sum(s[1] for s in used[j]) / sum(s[2] for s in used[j])
    return ratio, used


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("forcing", help="directory of the monthly SUMMA forcing files")
    p.add_argument("stations", help="USGSweatherData directory, one subdirectory per station elevation (m)")
    p.add_argument("attributes", help="SUMMA attributes file, for the HRU elevations")
    p.add_argument("output", help="directory to write the corrected forcing to")
    p.add_argument("--stakes", help="Input_<glacier>_Glaciological_Data.csv; scale precipitation to the winter balances")
    p.add_argument("--stake-min-elev", type=float, default=1000.0, help="lowest stake used, m")
    p.add_argument("--stake-window", type=float, default=200.0, help="stakes within this of an HRU's elevation, m")
    p.add_argument("--temp-crit-rain", type=float, default=273.16, help="wet-bulb rain/snow threshold, K")
    p.add_argument("--temp-range", type=float, default=2.0, help="wet-bulb range over which the phase changes, K")
    p.add_argument("--wet-threshold", type=float, default=0.1,
                   help="forcing precipitation above which an hour is corrected as wet, mm/h")
    p.add_argument("--gauge-shield", choices=sorted(CATCH),
                   help="wind shield of the station's weighing gauge (SA single-Alter, UN none); scales the "
                        "precipitation of HRUs without stakes to the undercatch-corrected gauge")
    p.add_argument("--spechum-ratios", help="CSV of hruId, month, ratio for the months the station has no humidity for")
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
    hours, frc = [], {v: [] for v in ("airtemp", "spechum", "SWRadAtm", "airpres", "pptrate")}
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
    nh = len(hru_ids)

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

    # hour type per HRU: 1 where the forcing precipitates, 0 where it is dry
    wet = (frc["pptrate"] * 3600.0 > a.wet_threshold).astype(int)

    def by_type(fit, valid, neutral):
        """fit(ok, j) over the valid hours per hour type, month and HRU, as (2, 13, nh).
        A type with under 5 days of hours takes the month's all-hour fit; a month under 30 days keeps neutral."""
        out = np.full((2, 13, nh), neutral)
        for j in range(nh):
            for m in range(1, 13):
                sel = (month == m) & valid(j)
                if sel.sum() <= 24 * 30:
                    continue
                out[:, m, j] = fit(sel, j)
                for k in (0, 1):
                    ok = sel & (wet[:, j] == k)
                    if ok.sum() > 24 * 5:
                        out[k, m, j] = fit(ok, j)
        return out

    # temperature: additive offset per hour type, month and HRU
    offset = by_type(lambda ok, j: np.mean(t_hru[ok, j] - frc["airtemp"][ok, j]),
                     lambda j: np.isfinite(t_hru[:, j]), 0.0)

    ratio = {v: np.ones(nh) for v in ("SWRadAtm", "airpres", "pptrate")}
    # shortwave: ratio of sums, no elevation adjustment
    sw = on_hours(stations[sw_z].get("SW", {}))
    sw = np.where(sw < 0, 0.0, sw)
    for j in range(nh):
        ok = np.isfinite(sw)
        if ok.sum() > 24 * 365:
            ratio["SWRadAtm"][j] = sw[ok].sum() / frc["SWRadAtm"][ok, j].sum()
    # humidity (monthly) and pressure need the main station's barometer and hygrometer
    rh, pk = on_hours(st_main.get("RH", {})), on_hours(st_main.get("P", {}))
    q = specific_humidity(rh, t_main, pk)
    q_ratio = by_type(lambda ok, j: q[ok].mean() / frc["spechum"][ok, j].mean(), lambda j: np.isfinite(q), 1.0)
    if a.spechum_ratios:
        for r in csv.DictReader(open(a.spechum_ratios)):
            j, m = list(hru_ids).index(int(r["hruId"])), int(r["month"])
            if np.all(q_ratio[:, m, j] == 1.0):
                q_ratio[:, m, j] = float(r["ratio"])
    for j in range(nh):
        ok = np.isfinite(pk) & np.isfinite(t_hru[:, j])
        if ok.sum() > 24 * 365:
            lap = lapse[month][ok]
            p_hru = pk[ok] * 1000.0 * (t_hru[ok, j] / (t_main[ok] + 273.15)) ** (-G / (R_DRY * lap))
            ratio["airpres"][j] = p_hru.mean() / frc["airpres"][ok, j].mean()

    # precipitation: stake winter balance against the corrected forcing's snowfall
    used = [[] for _ in z]
    if a.stakes:
        cols = np.arange(nh)[None, :]
        t_c = frc["airtemp"] + offset[wet, month[:, None], cols]
        tw = wet_bulb(t_c, frc["spechum"] * q_ratio[wet, month[:, None], cols], frc["airpres"] * ratio["airpres"][None, :])
        snow_frac = np.clip((a.temp_crit_rain + a.temp_range / 2.0 - tw) / a.temp_range, 0.0, 1.0)
        snow = np.maximum(frc["pptrate"], 0.0) * 3600.0 * snow_frac
        ratio["pptrate"], used = stake_ratios(a.stakes, hours, z, snow, a.stake_min_elev, a.stake_window)

    # precipitation: undercatch-corrected gauge for the HRUs the stakes leave out
    gauge = ""
    if a.gauge_shield:
        gz = max(stations, key=lambda s: len(stations[s].get("PR", {})))
        pr, u, tg = (on_hours(stations[gz].get(v, {})) for v in ("PR", "U", "T"))
        ok = np.isfinite(pr) & np.isfinite(u) & np.isfinite(tg)
        if ok.sum() > 24 * 365:
            caught = pr[ok] / catch_efficiency(u[ok], tg[ok], a.gauge_shield)
            forced = along_elevation(z, np.maximum(frc["pptrate"][ok], 0.0) * 3600.0, gz)
            gauge = f"{gz} m gauge, {ok.sum()} h, catch {pr[ok].sum() / caught.sum():.2f}"
            for j in range(nh):
                if not used[j]:
                    ratio["pptrate"][j] = caught.sum() / forced.sum()

    # write the corrected files
    os.makedirs(a.output, exist_ok=True)
    for f in files:
        with Dataset(f) as s, Dataset(os.path.join(a.output, os.path.basename(f)), "w", format="NETCDF4") as d:
            d.setncatts({k: s.getncattr(k) for k in s.ncattrs()})
            d.history = "corrected toward on-site weather stations by correct_forcing_to_station.py"
            for name, dim in s.dimensions.items():
                d.createDimension(name, None if dim.isunlimited() else len(dim))
            m_file = np.array([(EPOCH + timedelta(hours=int(h))).month for h in np.array(s["time"][:])])[:, None]
            wet_file = (np.array(s["pptrate"][:], float) * 3600.0 > a.wet_threshold).astype(int)
            cols = np.arange(nh)[None, :]
            for name, v in s.variables.items():
                fill = v.getncattr("_FillValue") if "_FillValue" in v.ncattrs() else None
                o = d.createVariable(name, v.dtype, v.dimensions, zlib=True, complevel=4, fill_value=fill)
                o.setncatts({k: v.getncattr(k) for k in v.ncattrs() if k != "_FillValue"})
                x = v[:]
                if name == "airtemp":
                    x = x + offset[wet_file, m_file, cols]
                elif name == "spechum":
                    x = x * q_ratio[wet_file, m_file, cols]
                elif name in ratio:
                    x = x * ratio[name][None, :]
                o[:] = x

    with open(os.path.join(a.output, "corrections.csv"), "w") as fh:
        w = csv.writer(fh)
        types = ("dry", "wet")
        w.writerow(["hruId", "elevation_m", "SWRadAtm_ratio", "airpres_ratio", "pptrate_ratio", "stakes_used"]
                   + [f"airtemp_offset_{t}_K_{m:02d}" for t in types for m in range(1, 13)]
                   + [f"spechum_ratio_{t}_{m:02d}" for t in types for m in range(1, 13)])
        for j, h in enumerate(hru_ids):
            w.writerow([h, z[j], f"{ratio['SWRadAtm'][j]:.4f}", f"{ratio['airpres'][j]:.4f}",
                        f"{ratio['pptrate'][j]:.4f}", len(used[j])]
                       + [f"{offset[k, m, j]:.3f}" for k in (0, 1) for m in range(1, 13)]
                       + [f"{q_ratio[k, m, j]:.4f}" for k in (0, 1) for m in range(1, 13)])
        w.writerow([])
        w.writerow(["lapse_rate_K_per_km"] + [f"{1000 * lapse[m]:.2f}" for m in range(1, 13)])
        w.writerow(["wet_threshold_mm_per_h", f"{a.wet_threshold}"])
        if gauge:
            w.writerow(["pptrate_gauge", f"{a.gauge_shield} shield, {gauge}"])

    print("lapse rate (K/km) by month:", " ".join(f"{1000 * lapse[m]:.1f}" for m in range(1, 13)))
    if gauge:
        print(f"precipitation without stakes from the {a.gauge_shield}-shielded {gauge}")
    for j, h in enumerate(hru_ids):
        stk = (f", precip x{ratio['pptrate'][j]:.3f} from {len(used[j])} stake-years at "
               f"{min(s[0] for s in used[j]):.0f}-{max(s[0] for s in used[j]):.0f} m") if used[j] else (
               f", precip x{ratio['pptrate'][j]:.3f} from the gauge" if gauge else "")
        print(f"HRU {h} ({z[j]:.0f} m): SW x{ratio['SWRadAtm'][j]:.3f}, P x{ratio['airpres'][j]:.3f}{stk}")
        for k, t in enumerate(("dry", "wet")):
            print(f"    T offset by month, {t} (K) " + " ".join(f"{offset[k, m, j]:+.1f}" for m in range(1, 13)))
            print(f"    q ratio by month, {t}      " + " ".join(f"{q_ratio[k, m, j]:.2f}" for m in range(1, 13)))
    print(f"{len(files)} files written to {a.output}")


if __name__ == "__main__":
    main()
