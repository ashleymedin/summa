#!/usr/bin/env python3
"""Write a USGS Benchmark Glacier mass balance as a SUMMA calibration observation file.

The source is the glacier-wide solutions table of the USGS data release "Glacier-Wide Mass Balance and
Compiled Data Inputs" (https://www.sciencebase.gov/catalog/item/6441de03d34ee8d4ade7d2a5), one file per
glacier, Output_<Glacier>_Glacier_Wide_solutions_calibrated.csv.  ScienceBase does not allow scripted
downloads, so the file is downloaded by hand and given here.

Each year carries the winter (Bw), summer (Bs) and annual (Ba) balance in m w.e., and the modelled
dates of the glacier's mass maximum (Bw_Date) and minimum (Ba_Date).  The balances are stratigraphic,
measured between those extremes, and each value spans:

  seasonal  winter, the previous year's Ba_Date to this year's Bw_Date, and summer, Bw_Date to Ba_Date,
            both in one series in date order
  annual    the previous year's Ba_Date to this year's Ba_Date

A winter or annual balance in the first year, or after a year with no Ba_Date, has no start and is left
out.  Values are written in mm w.e., each stamped at the end of the day its span ends on; the target
compares them against the simulation's own extremes near that date, so the hour does not matter.

    ./acquire_glacier_mass_balance.py Output_Gulkana_Glacier_Wide_solutions_calibrated.csv annual \\
        gulkana_mass_balance_annual.nc
"""

import argparse
import csv
import datetime as dt
import math

from observation_netcdf import write_observations


def read_solutions(path):
    """Rows of the glacier-wide solutions table, by year, with dates parsed and NaN kept as None."""
    def number(text):
        value = float(text)
        return None if math.isnan(value) else value

    def date(text):
        text = text.strip()
        if not text or text.lower() == "nan":
            return None
        return dt.datetime.strptime(text, "%Y/%m/%d")

    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            rows.append({
                "year": int(row["Year"]),
                "bw": number(row["Bw"]), "bs": number(row["Bs"]), "ba": number(row["Ba"]),
                "bw_date": date(row["Bw_Date"]), "ba_date": date(row["Ba_Date"]),
            })
    rows.sort(key=lambda r: r["year"])
    return rows


def spans(rows, balance):
    """(start, end, value in m w.e.) of each balance that has both ends and a value, in date order."""
    out = []
    previous = None
    for row in rows:
        # consecutive years only, so a gap in the record does not become a multi-year balance
        last_min = previous["ba_date"] if previous is not None and previous["year"] == row["year"] - 1 else None
        if balance == "seasonal":
            candidates = [(last_min, row["bw_date"], row["bw"]), (row["bw_date"], row["ba_date"], row["bs"])]
        else:
            candidates = [(last_min, row["ba_date"], row["ba"])]
        for start, end, value in candidates:
            if start is not None and end is not None and value is not None and end > start:
                out.append((start, end, value))
        previous = row
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("solutions", help="Output_<Glacier>_Glacier_Wide_solutions_calibrated.csv")
    parser.add_argument("balance", choices=("seasonal", "annual"))
    parser.add_argument("output", help="NetCDF file to write")
    args = parser.parse_args()

    balances = spans(read_solutions(args.solutions), args.balance)
    if not balances:
        raise SystemExit(f"no {args.balance} balances in {args.solutions}")

    one_day = dt.timedelta(days=1)
    starts = [start + one_day for start, _, _ in balances]
    ends = [end + one_day for _, end, _ in balances]
    values = [1000.0 * value for _, _, value in balances]

    write_observations(
        args.output, starts, ends, values,
        varname="mb_obs",
        units="mm",
        long_name=f"glacier-wide {args.balance} mass balance, water equivalent",
        cell_methods="time: sum",
        attrs={
            "title": f"Glacier-wide {args.balance} mass balance, for SUMMA calibration",
            "source": "USGS Benchmark Glacier Project, Glacier-Wide Mass Balance and Compiled Data Inputs, "
                      "https://www.sciencebase.gov/catalog/item/6441de03d34ee8d4ade7d2a5",
            "comment": "Stratigraphic balance between the modelled dates of mass minimum (Ba_Date) and "
                       "maximum (Bw_Date), geodetically calibrated, per unit glacier area.",
            "history": f"created {dt.date.today()} by acquire_glacier_mass_balance.py from "
                       f"{args.solutions.split('/')[-1]}",
        },
        time_zone="dates only",
    )
    print(f"{args.output}: {len(values)} {args.balance} balances, "
          f"{balances[0][1]:%Y} to {balances[-1][1]:%Y}")


if __name__ == "__main__":
    main()
