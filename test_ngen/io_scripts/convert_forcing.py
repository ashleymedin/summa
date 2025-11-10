#!/usr/bin/env python3
"""
Convert SUMMA forcing files to NGen forcing format or vice versa.
Usage:
  python3 convert_forcing.py input_file.nc output_file.nc typeOfConversion
"""
import argparse
import xarray as xr
import numpy as np
import re
import os
from ast import literal_eval

def _parse_id(fid):
    # try direct int, else extract trailing digits, else first digits, else return original string
    try:
        return int(fid)
    except Exception:
        s = str(fid).strip()
        m = re.search(r'(\d+)$', s) or re.search(r'(\d+)', s)
        return int(m.group(1)) if m else s

def main():
    p = argparse.ArgumentParser()
    p.add_argument("input_file", help="input_file.nc (SUMMA or NGen forcing file)")
    p.add_argument("output_file", help="output_file.nc (converted forcing file)")
    p.add_argument("typeOfConversion", help="type of conversion: summa_to_ngen or ngen_to_summa")
    args = p.parse_args()

    # open input file
    f_ds = xr.open_dataset(args.input_file)

    # do conversion
    if args.typeOfConversion == "ngen_to_summa":
        f_var = "ids"
        if f_var not in f_ds:
            raise RuntimeError(f"NGen forcing file missing '{f_var}' variable")
        f_ids_raw = np.array([str(x) for x in f_ds[f_var].values])
        f_ids = np.array([_parse_id(fid) for fid in f_ids_raw])



    elif args.typeOfConversion == "summa_to_ngen":
        f_var = "hruId"
        if f_var not in f_ds:
            raise RuntimeError(f"SUMMA forcing file missing '{f_var}' variable")
        f_ids = f_ds[f_var].values

            

    else:
        raise RuntimeError("typeOfConversion must be either 'summa_to_ngen' or 'ngen_to_summa'")






if __name__ == "__main__":
    main()