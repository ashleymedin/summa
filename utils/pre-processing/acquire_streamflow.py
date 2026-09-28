#!/usr/bin/env python3
"""Fetch daily streamflow for a gauge and write it as a SUMMA calibration observation file.

Two sources, both public and needing no account:

  wsc   Water Survey of Canada, through the MSC GeoMet API (api.weather.gc.ca)
  usgs  U.S. Geological Survey, through the Water Data API (api.waterdata.usgs.gov)

The Bow at Banff, the gauge of the bundled test domain:

    ./acquire_streamflow.py wsc 05BB001 1980-01-01 1990-12-31 CAN_05BB001_daily_flow_observations.nc

A USGS gauge, given by site number:

    ./acquire_streamflow.py usgs 09380000 2000-01-01 2010-12-31 USGS_09380000_daily_flow_observations.nc

The file holds `q_obs`, the daily mean in m3 s-1; USGS values arrive in ft3 s-1 and are converted.
A day with no value is written as NaN rather than dropped, so the record stays evenly spaced.
"""

import argparse
import datetime as dt
import json
import sys
import urllib.parse
import urllib.request

import numpy as np

from observation_netcdf import write_observations

CFS_TO_CMS = 0.028316846592  # ft3 s-1 to m3 s-1

SOURCES = {
    "wsc": {
        "url": "https://api.weather.gc.ca/collections/hydrometric-daily-mean/items",
        "name": "Water Survey of Canada, MSC GeoMet hydrometric-daily-mean",
    },
    "usgs": {
        "url": "https://api.waterdata.usgs.gov/ogcapi/v0/collections/daily/items",
        "name": "U.S. Geological Survey, Water Data API daily values",
    },
}


def fetch_features(url, params):
    """Every feature of an OGC API query, following its `next` links."""
    features = []
    url = url + "?" + urllib.parse.urlencode(params)
    while url:
        with urllib.request.urlopen(url, timeout=120) as response:
            page = json.load(response)
        features.extend(page.get("features", []))
        url = next((link["href"] for link in page.get("links", []) if link.get("rel") == "next"), None)
    return features


def fetch_wsc(station, start, end):
    """Daily discharge (m3 s-1) by date for a WSC station."""
    features = fetch_features(SOURCES["wsc"]["url"], {
        "STATION_NUMBER": station,
        "datetime": f"{start}/{end}",
        "f": "json",
        "limit": 10000,
    })
    return {
        dt.date.fromisoformat(f["properties"]["DATE"]): f["properties"]["DISCHARGE"]
        for f in features
        if f["properties"]["DISCHARGE"] is not None
    }


def fetch_usgs(site, start, end):
    """Daily mean discharge (m3 s-1) by date for a USGS site."""
    features = fetch_features(SOURCES["usgs"]["url"], {
        "monitoring_location_id": f"USGS-{site}",
        "parameter_code": "00060",   # discharge
        "statistic_id": "00003",     # daily mean
        "time": f"{start}/{end}",
        "f": "json",
        "limit": 10000,
    })
    flows = {}
    for f in features:
        props = f["properties"]
        if props["value"] is None:
            continue
        if props["unit_of_measure"] != "ft^3/s":
            sys.exit(f"unexpected USGS unit {props['unit_of_measure']!r}")
        flows[dt.date.fromisoformat(props["time"][:10])] = float(props["value"]) * CFS_TO_CMS
    return flows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("source", choices=sorted(SOURCES))
    parser.add_argument("gauge", help="WSC station number or USGS site number")
    parser.add_argument("start", type=dt.date.fromisoformat, help="first day, YYYY-MM-DD")
    parser.add_argument("end", type=dt.date.fromisoformat, help="last day, YYYY-MM-DD")
    parser.add_argument("output", help="NetCDF file to write")
    args = parser.parse_args()

    fetch = fetch_wsc if args.source == "wsc" else fetch_usgs
    flows = fetch(args.gauge, args.start, args.end)
    if not flows:
        sys.exit(f"no discharge for {args.gauge} between {args.start} and {args.end}")

    days = [args.start + dt.timedelta(days=i) for i in range((args.end - args.start).days + 1)]
    starts = [dt.datetime.combine(day, dt.time()) for day in days]
    ends = [start + dt.timedelta(days=1) for start in starts]
    values = np.asarray([flows.get(day, np.nan) for day in days], dtype="f8")

    write_observations(
        args.output, starts, ends, values,
        varname="q_obs",
        units="m3 s-1",
        long_name="observed streamflow",
        attrs={
            "title": f"Daily streamflow at {args.source.upper()} {args.gauge}, for SUMMA calibration",
            "source": SOURCES[args.source]["name"],
            "history": f"created {dt.date.today()} by acquire_streamflow.py",
        },
    )
    print(f"{args.output}: {np.isfinite(values).sum()} of {len(values)} days")


if __name__ == "__main__":
    main()
