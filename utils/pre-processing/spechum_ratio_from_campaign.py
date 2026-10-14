#!/usr/bin/env python3
"""Monthly specific-humidity ratios, station over forcing, from a station campaign outside the forcing years.

The campaign is compared against ARCO-ERA5 2 m humidity at the station's grid cell, and that comparison is
carried onto each HRU's forcing through the same grid cell over --bridge-years of the forcing record:

  ratio(hru, month) = [q_station / q_ERA5,2m]_campaign / [q_forcing,hru / q_ERA5,2m]_bridge

so the reanalysis-vs-station bias is assumed stable in time, and the bridge absorbs both the forcing's height
(lowest model level, not 2 m) and its remapping to HRUs.  Months with under --min-days of campaign hours are
left out.  Station q is from its relative humidity, temperature and own barometer (Bolton 1980, not lapsed);
station values are hourly means of the samples in the hour ending at each ERA5 time.

The station file is a Campbell logger table with UTC TIMESTAMP, AirTC_Avg (degC), RH (%) and BV_BP_Avg (hPa),
as in Wilson et al.'s Gulkana 2024 data release.  Writes hruId, month, ratio rows for
correct_forcing_to_station.py --spechum-ratios.  Needs xarray, gcsfs and zarr; the ERA5 pulls are cached.

    ./spechum_ratio_from_campaign.py $S/domain_Gulkana/forcing/SUMMA_input $S/domain_Gulkana/settings/SUMMA/attributes.nc \\
        $S/domain_Gulkana/Wilson_gulkana_2024_data/1_meteo/gulkanaB_AWS_2024.csv 63.285514 -145.410336 \\
        $S/domain_Gulkana/forcing/spechum_ratios_wilson2024.csv
"""

import argparse
import csv
import glob
import os
from datetime import datetime, timezone

import numpy as np
import xarray as xr

EPOCH = np.datetime64("1900-01-01T00:00")
ARCO = "gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"


def specific_humidity(e_kpa, p_kpa):
    """Specific humidity (kg kg-1) from vapour pressure and pressure (kPa)."""
    return 0.622 * e_kpa / (p_kpa - 0.378 * e_kpa)


def sat_vapour(t_c):
    """Saturation vapour pressure over water (kPa) at t_c (degC), Bolton (1980)."""
    return 0.6112 * np.exp(17.67 * t_c / (t_c + 243.5))


def read_campaign(path):
    """Hourly mean station q (kg kg-1) and RH (%) keyed by hours since 1900, from a Campbell logger table."""
    with open(path) as fh:
        rows = csv.reader(fh)
        head = next(rows)
        next(rows), next(rows)
        it, iT, iRH, iP = (head.index(c) for c in ("TIMESTAMP", "AirTC_Avg", "RH", "BV_BP_Avg"))
        sums = {}
        for r in rows:
            try:
                t_c, rh, p = float(r[iT]), float(r[iRH]), float(r[iP]) / 10.0
            except (ValueError, IndexError):
                continue
            if not (np.isfinite(t_c) and np.isfinite(rh) and np.isfinite(p)):
                continue
            t = datetime.strptime(r[it], "%m/%d/%Y %H:%M").replace(tzinfo=timezone.utc)
            hour = int(np.ceil((t - datetime(1900, 1, 1, tzinfo=timezone.utc)).total_seconds() / 3600.0))
            s = sums.setdefault(hour, [0.0, 0.0, 0])
            s[0] += specific_humidity(rh / 100.0 * sat_vapour(t_c), p)
            s[1] += rh
            s[2] += 1
    hours = np.array(sorted(sums))
    return hours, np.array([sums[h][0] / sums[h][2] for h in hours]), np.array([sums[h][1] / sums[h][2] for h in hours])


def era5_q(hours, lat, lon, cache):
    """ARCO-ERA5 2 m specific humidity at the grid cell nearest (lat, lon) on the given hours since 1900."""
    if os.path.exists(cache):
        c = np.load(cache)
        if np.array_equal(c["hours"], hours):
            return c["q"]
    import dask
    import gcsfs
    ds = xr.open_zarr(gcsfs.GCSFileSystem(token="anon").get_mapper(ARCO), consolidated=True, chunks={})
    cell = ds[["2m_dewpoint_temperature", "surface_pressure"]].sel(latitude=lat, longitude=lon % 360.0, method="nearest")
    times = EPOCH + hours.astype("timedelta64[h]")
    with dask.config.set(scheduler="threads", num_workers=32):
        x = cell.sel(time=times).compute()
    print(f"ERA5 cell {float(x.latitude):.2f} N {float(x.longitude) - 360:.2f} E, {len(hours)} h")
    q = specific_humidity(sat_vapour(x["2m_dewpoint_temperature"].values - 273.15), x["surface_pressure"].values / 1000.0)
    np.savez(cache, hours=hours, q=q)
    return q


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("forcing", help="directory of the monthly SUMMA forcing files")
    p.add_argument("attributes", help="SUMMA attributes file, for the HRU elevations")
    p.add_argument("station", help="campaign station file")
    p.add_argument("lat", type=float, help="station latitude, degrees N")
    p.add_argument("lon", type=float, help="station longitude, degrees E")
    p.add_argument("output", help="CSV of hruId, month, ratio")
    p.add_argument("--bridge-years", default="2020-2022", help="forcing years, first-last, compared to ERA5 2 m")
    p.add_argument("--bridge-stride", type=int, default=3, help="sample every this many bridge hours")
    p.add_argument("--min-days", type=float, default=10.0, help="fewest campaign days for a month's ratio")
    a = p.parse_args()

    st_hours, q_st, rh_st = read_campaign(a.station)
    st_month = (EPOCH + st_hours.astype("timedelta64[h]")).astype("datetime64[M]").astype(int) % 12 + 1
    months = [m for m in range(1, 13) if (st_month == m).sum() >= 24 * a.min_days]
    q_era = era5_q(st_hours, a.lat, a.lon, os.path.splitext(a.output)[0] + "_era5_campaign.npz")

    # forcing over the bridge months
    y0, y1 = (int(y) for y in a.bridge_years.split("-"))
    frc_hours, frc_q = [], []
    for f in sorted(glob.glob(os.path.join(a.forcing, "*.nc"))):
        with xr.open_dataset(f, engine="h5netcdf", decode_times=False) as d:
            if "spechum" not in d:
                continue
            h = d["time"].values.astype(int)
            mo = (EPOCH + h.astype("timedelta64[h]")).astype("datetime64[M]")
            yr, m = mo.astype("datetime64[Y]").astype(int) + 1970, mo.astype(int) % 12 + 1
            keep = (yr >= y0) & (yr <= y1) & np.isin(m, months)
            if keep.any():
                frc_hours.append(h[keep])
                frc_q.append(d["spechum"].values[keep])
                hru_ids = d["hruId"].values.astype(int)
    frc_hours, frc_q = np.concatenate(frc_hours), np.concatenate(frc_q)
    sub = np.arange(0, len(frc_hours), a.bridge_stride)
    frc_hours, frc_q = frc_hours[sub], frc_q[sub]
    frc_month = (EPOCH + frc_hours.astype("timedelta64[h]")).astype("datetime64[M]").astype(int) % 12 + 1
    q_bridge = era5_q(frc_hours, a.lat, a.lon, os.path.splitext(a.output)[0] + f"_era5_bridge_{y0}-{y1}.npz")

    with xr.open_dataset(a.attributes, engine="h5netcdf") as d:
        elev = dict(zip(d["hruId"].values.astype(int), d["elevation"].values))
    with open(a.output, "w") as fh:
        w = csv.writer(fh)
        w.writerow(["hruId", "month", "ratio", "station_over_era5", "forcing_over_era5", "campaign_hours", "station_rh"])
        for m in months:
            sel = (st_month == m) & np.isfinite(q_era)
            camp = q_st[sel].mean() / q_era[sel].mean()
            bsel = frc_month == m
            print(f"month {m:2d}: {sel.sum():4d} h, station RH {rh_st[sel].mean():.0f}%, station/ERA5 q {camp:.3f}")
            for j, h in enumerate(hru_ids):
                bridge = frc_q[bsel, j].mean() / q_bridge[bsel].mean()
                w.writerow([h, m, f"{camp / bridge:.4f}", f"{camp:.4f}", f"{bridge:.4f}", int(sel.sum()),
                            f"{rh_st[sel].mean():.1f}"])
                print(f"    HRU {h} ({elev[h]:.0f} m): forcing/ERA5 q {bridge:.3f} -> ratio {camp / bridge:.3f}")
    print(f"written {a.output}")


if __name__ == "__main__":
    main()
