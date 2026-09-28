#!/usr/bin/env python3
"""Convert a GRACE terrestrial-water-storage record into the observation file SUMMA calibrates against.

A calibration target reads its observations from a NetCDF file holding a CF time coordinate and a
named variable carrying its own units - the same shape as the streamflow observations already in
this repository, `CAN_05BB001_daily_flow_observations.nc`.  Keeping every kind of observation in
that one shape is what lets the model read them all with one routine, instead of growing a reader
per product.

GRACE is distributed as monthly mass anomalies, and the basin-average records this was written for
arrive as a dated CSV with a column per processing centre.  This turns one of those columns into the
NetCDF form.

    ./grace_tws_to_netcdf.py Gulkana_HRUs_GRUs_grace_tws_anomaly.csv gulkana_grace_tws.nc \
        --column grace_jpl_anomaly

The anomalies are equivalent water thickness in millimetres, relative to the baseline of whichever
product they came from (2004-2009 for the JPL mascons).  That baseline matters: a calibration target
reading this file has to reference the simulated storage to the same period, which is what the
target's `baseline_start` and `baseline_end` are for.

Fetching the mascon grids and reducing them over a basin polygon is a separate, heavier job - see
`extract_grace_tws.py` in the SYMFLUENCE data tree, which produced the CSVs this reads.
"""

import argparse
import csv
import datetime as dt
import sys

import netCDF4 as nc
import numpy as np

# The reference the time coordinate is written against.  Any CF-style "<unit> since <date>" works;
# this one matches the bundled streamflow observations, so the two files read alike.
TIME_UNITS = "minutes since 1950-01-01"
TIME_REFERENCE = dt.datetime(1950, 1, 1)


def read_column(path, column):
    """Read one named column of a dated CSV, returning (datetimes, values).

    The date is the first column, named or not: these products are written with an unnamed index
    column.  A blank field is a month with no solution, and becomes NaN rather than a number.
    """
    times, values = [], []
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if column not in header:
            named = [name for name in header[1:] if name]
            sys.exit(f"column {column!r} is not in {path}\navailable: {', '.join(named)}")
        index = header.index(column)

        for row in reader:
            if not row or not row[0].strip():
                continue
            times.append(dt.datetime.strptime(row[0].strip()[:10], "%Y-%m-%d"))
            field = row[index].strip() if index < len(row) else ""
            values.append(float(field) if field else np.nan)

    if not times:
        sys.exit(f"{path} holds no dated rows")
    return times, np.asarray(values, dtype="f8")


def write_netcdf(path, times, values, varname, units, long_name, source_column, source_file):
    """Write the observations in the form a calibration target reads."""
    minutes = np.asarray(
        [(t - TIME_REFERENCE).total_seconds() / 60.0 for t in times], dtype="i8"
    )

    with nc.Dataset(path, "w", format="NETCDF4") as out:
        out.createDimension("time", len(minutes))

        time_var = out.createVariable("time", "i8", ("time",))
        time_var.units = TIME_UNITS
        time_var.long_name = "time"
        time_var.calendar = "standard"
        time_var[:] = minutes

        # NaN is what the metrics treat as missing, so it is also what the fill value has to be:
        # a month GRACE has no solution for must not read as a number.
        obs = out.createVariable(varname, "f8", ("time",), fill_value=np.nan)
        obs.units = units
        obs.long_name = long_name
        obs[:] = values

        out.title = "GRACE terrestrial water storage, for SUMMA calibration"
        out.source = f"{source_file}, column {source_column}"
        out.history = f"created {dt.datetime.now():%Y-%m-%d} by grace_tws_to_netcdf.py"
        out.comment = (
            "Monthly mass anomaly as equivalent water thickness, relative to the source product's "
            "own baseline. A calibration target reading this must reference the simulated storage "
            "to the same period."
        )

    finite = np.isfinite(values)
    print(f"wrote {path}")
    print(f"  variable:  {varname} [{units}]")
    print(f"  period:    {times[0]:%Y-%m} to {times[-1]:%Y-%m}, {len(times)} months")
    print(f"  finite:    {finite.sum()} of {len(values)}")
    if finite.any():
        print(f"  range:     {values[finite].min():.2f} to {values[finite].max():.2f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", help="dated CSV holding the GRACE record")
    parser.add_argument("netcdf", help="observation file to write")
    parser.add_argument("--column", default="grace_jpl_anomaly",
                        help="CSV column to take (default: %(default)s; JPL for glacierized basins)")
    parser.add_argument("--varname", default="tws_obs",
                        help="name of the variable to write (default: %(default)s)")
    parser.add_argument("--units", default="mm",
                        help="units of the anomaly (default: %(default)s, equivalent water thickness)")
    parser.add_argument("--long-name", default="GRACE terrestrial water storage anomaly",
                        help="long_name attribute of the variable")
    args = parser.parse_args()

    times, values = read_column(args.csv, args.column)
    write_netcdf(args.netcdf, times, values, args.varname, args.units, args.long_name,
                 args.column, args.csv)


if __name__ == "__main__":
    main()
