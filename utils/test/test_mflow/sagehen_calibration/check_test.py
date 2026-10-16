#!/usr/bin/env python3
"""Check run_test.sh's three runs: the truth reproduces its targets, the search ran, and each worker's
MODFLOW files hold its last trial's multipliers.

Usage:
    check_test.py <work dir> <PARAMS table> <population> <generations>
"""

import os
import re
import sys

import numpy as np
import netCDF4 as nc

work, params, pop, ngen = sys.argv[1], [l.split(None, 5) for l in sys.argv[2].splitlines()], int(sys.argv[3]), int(sys.argv[4])
mf6 = os.path.join(work, "domain", "mf6")

fail = False


def check(ok, what):
    global fail
    print(f"  {'ok  ' if ok else 'FAIL'} {what}")
    fail |= not ok


def read_array(path):
    return np.array(open(path).read().split(), float)


def read_column(path, block, column):
    """One column of the data rows of a block, 1-based."""
    text = open(path).read()
    body = re.search(rf"BEGIN\s+{block}\b.*?\n(.*?)\nEND\s+{block}", text, flags=re.S | re.I).group(1)
    rows = [l.split() for l in body.splitlines() if l.strip() and not l.lstrip().startswith("#")]
    return np.array([float(r[column - 1]) for r in rows])


# the true parameters reproduce their own targets
t = nc.Dataset(os.path.join(work, "truth", "output", "sagehen_calibration.nc"))
names = [b"".join(row).decode().strip() for row in np.asarray(t["target_name"][:])]
sense = [s.strip() for s in t["objective"].sense.split(",")]
f0 = np.asarray(t["objective"][:], float).ravel()
for name, s, v in zip(names, sense, f0):
    ok = (1.0 - v < 1e-6) if s == "maximize" else (abs(v) < 1e-6)
    check(ok, f"the true parameters reproduce {name}: {v:.3e}")

# every package the coupler drives moved water or heat in the truth trial: its budget is the log's last
budget = open(os.path.join(work, "truth", "run.log")).read().split("coupled water budget over the run")[-1]


def flow(label):
    return float(re.search(rf"{re.escape(label)}\s*:\s*(\S+)", budget).group(1))


check(flow("of which into UZF") > 0.0, f"UZF takes recharge: {flow('of which into UZF'):.1f} m3")
check(flow("base exchange into MODFLOW, returned from the column") != 0.0,
      f"the base exchange moves water: {flow('base exchange into MODFLOW, returned from the column'):.1f} m3")
check(flow("returned as baseflow") > 0.0, f"SFR returns baseflow: {flow('returned as baseflow'):.1f} m3")
check(flow("returned at land surface") > 0.0, f"DRN returns water at land surface: {flow('returned at land surface'):.1f} m3")
check(flow("taken as groundwater ET") > 0.0, f"EVT takes groundwater ET: {flow('taken as groundwater ET'):.1f} m3")
sent, got = flow("soil-base conduction sent (J)"), flow("GWE ESL loading received (J)")
check(abs(got - sent) <= 1e-6 * abs(sent), f"GWE's ESL receives the heat SUMMA sends: {got:.4e} J of {sent:.4e}")
check(abs(flow("mapping residual")) <= 1e-6 * flow("SUMMA drainage sent"), "the HRU-to-cell map conserves drainage")

# the search ran every trial and scored every target
d = nc.Dataset(os.path.join(work, "search", "output", "sagehen_calibration.nc"))
obj = np.ma.filled(d["objective"][:].astype(float), np.nan).T   # (target, sample)
failed = np.asarray(d["trial_failed"][:]) == 1
check(obj.shape == (len(names), pop * ngen), f"objective is {len(names)} targets x {pop * ngen} trials")
check(not failed.any(), f"{int((~failed).sum())} of {failed.size} trials ran")
check(bool(np.all(np.isfinite(obj[:, ~failed]))), "every objective of a trial that ran is finite")
check(bool((np.asarray(d["pareto_front"][:]) == 1).any()), "the front is not empty")

# a worker's MODFLOW files hold its last trial's multipliers on the model's own values
rank = np.asarray(d["worker_rank"][:])
tc = np.asarray(d["completion_time"][:]).reshape(rank.size, -1).astype(np.int64)   # y m d utc h min s ms
done = ((((((tc[:, 0] * 13 + tc[:, 1]) * 32 + tc[:, 2]) * 24 + tc[:, 4]) * 60 + tc[:, 5]) * 60 + tc[:, 6]) * 1000 + tc[:, 7])
zone3 = np.loadtxt(os.path.join(mf6, "syzone1.txt")).ravel() == 3
for name, truth, lo, hi, scale, spec in params:
    if spec == "summa":
        continue
    v = np.asarray(d[name][:], float)
    check(bool(np.all((v >= float(lo)) & (v <= float(hi))) and np.ptp(v) > 0.0), f"{name} varies within [{lo}, {hi}]")
    files = re.findall(r'"([^"]+)"', spec.split("|")[0])
    column = re.search(r"column = (\d+)", spec)
    for r in np.unique(rank):
        last = np.where(rank == r)[0][np.argmax(done[rank == r])]
        run_dir = os.path.join(work, "search", "output", f"modflow_rank{r:04d}")
        for f in files:
            if column:
                base = read_column(os.path.join(mf6, f), "PACKAGEDATA", int(column.group(1)))
                here = read_column(os.path.join(run_dir, f), "PACKAGEDATA", int(column.group(1)))
            else:
                base, here = read_array(os.path.join(mf6, f)), read_array(os.path.join(run_dir, f))
            want = np.full(base.shape, v[last])
            if "zone_file" in spec:
                want = np.where(zone3, v[last], 1.0)
            k = base != 0
            check(bool(np.allclose(here[k] / base[k], want[k], rtol=1e-6)), f"rank {r}'s {f} is its last trial's {name} {v[last]:.4f}")

sys.exit(1 if fail else 0)
