#!/usr/bin/env python3
"""Build ex-gwf-sagehen-gwe: ex-gwf-sagehen-sfr with a GWE model on the same grid.

GWE carries the aquifer temperature by advection and conduction (ADV, CND, EST) and along
the SFR reaches (SFE). Recharge brings SUMMA's drainage temperature in through RCH's
auxiliary TEMPERATURE, which the coupler overwrites every step and SSM names as the
recharge source. The coupler returns GWE's water-table temperature to SUMMA as
scalarAquiferTemp. GWE runs in degrees Celsius, in the GWF model's metres and seconds.

Usage:
    build_sagehen_gwe.py [<output model dir>]      (default ../ex-gwf-sagehen-gwe)
"""

import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SFR = os.path.join(HERE, "..", "ex-gwf-sagehen-sfr")
GWE = "SAGEHEN_GWE"

STRT = 6.0            # degC, initial aquifer and reach temperature
POROSITY = 0.25
RHO_SOLID = 2650.0    # kg m-3
CP_SOLID = 840.0      # J kg-1 degC-1
K_WATER = 0.58        # W m-1 degC-1
K_SOLID = 2.5         # W m-1 degC-1
K_BED = 2.0           # W m-1 degC-1, streambed
BED_THK = 0.5         # m, as build_sagehen_sfr.py's RBED_THK


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


def main(out_dir):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    shutil.copytree(SFR, out_dir, ignore=shutil.ignore_patterns(
        "*.lst", "*.hds", "*.cbc", "*.grb", "*.bin", ".DS_Store", "mfsim.stdout"))

    with open(os.path.join(SFR, "sagehen.sfr")) as f:
        nreach = int(re.search(r"NREACHES\s+(\d+)", f.read()).group(1))

    # recharge carries a temperature
    rch = os.path.join(out_dir, "sagehen.rcha")
    with open(rch) as f:
        text = f.read()
    text = text.replace("  READASARRAYS\n", "  READASARRAYS\n  AUXILIARY  TEMPERATURE\n")
    text = text.replace("    CONSTANT  1.0e-8\n", f"    CONSTANT  1.0e-8\n  TEMPERATURE\n    CONSTANT  {STRT}\n")
    write(rch, "# The coupler also overwrites TEMPERATURE (__INPUT__/SAGEHEN/RCHA/AUXVAR, degC) with SUMMA's\n"
               "# drainage temperature; GWE's SSM reads it as the recharge temperature.\n" + text)

    write(os.path.join(out_dir, "mfsim.nam"), f"""\
# MODFLOW 6 simulation name file - Sagehen with SFR and GWE, configured for the SUMMA coupler
#   * GWF model "SAGEHEN" and GWE model "{GWE}" on the same DIS grid
#   * TDIS in SECONDS, MODFLOW time step = SUMMA forcing data_step (3600 s here)
#   * RCH package "RCHA" with READASARRAYS and AUXILIARY TEMPERATURE
# Written by tools/build_sagehen_gwe.py.
BEGIN OPTIONS
END OPTIONS

BEGIN TIMING
  TDIS6  sagehen.tdis
END TIMING

BEGIN MODELS
  GWF6  sagehen.nam      SAGEHEN
  GWE6  sagehen_gwe.nam  {GWE}
END MODELS

BEGIN EXCHANGES
  GWF6-GWE6  sagehen.gwfgwe  SAGEHEN  {GWE}
END EXCHANGES

BEGIN SOLUTIONGROUP 1
  IMS6  sagehen.ims      SAGEHEN
  IMS6  sagehen_gwe.ims  {GWE}
END SOLUTIONGROUP 1
""")
    write(os.path.join(out_dir, "sagehen.gwfgwe"), "BEGIN OPTIONS\nEND OPTIONS\n")

    write(os.path.join(out_dir, "sagehen_gwe.nam"), """\
# GWE model name file: aquifer temperature on the GWF grid, degC.
BEGIN OPTIONS
  SAVE_FLOWS
END OPTIONS

BEGIN PACKAGES
  DIS6  sagehen_gwe.dis  dis
  IC6   sagehen_gwe.ic   ic
  ADV6  sagehen_gwe.adv  adv
  CND6  sagehen_gwe.cnd  cnd
  EST6  sagehen_gwe.est  est
  SSM6  sagehen_gwe.ssm  ssm
  SFE6  sagehen_gwe.sfe  SFE
  OC6   sagehen_gwe.oc   oc
END PACKAGES
""")

    with open(os.path.join(out_dir, "sagehen.dis")) as f:
        dis = f.read()
    write(os.path.join(out_dir, "sagehen_gwe.dis"),
          dis.replace("  LENGTH_UNITS  METERS\n", "  LENGTH_UNITS  METERS\n  NOGRB\n"))

    write(os.path.join(out_dir, "sagehen_gwe.ic"),
          f"BEGIN GRIDDATA\n  STRT\n    CONSTANT  {STRT}\nEND GRIDDATA\n")
    write(os.path.join(out_dir, "sagehen_gwe.adv"), "BEGIN OPTIONS\n  SCHEME  TVD\nEND OPTIONS\n")
    write(os.path.join(out_dir, "sagehen_gwe.cnd"), f"""\
BEGIN OPTIONS
  XT3D_OFF
END OPTIONS

BEGIN GRIDDATA
  KTW
    CONSTANT  {K_WATER}
  KTS
    CONSTANT  {K_SOLID}
END GRIDDATA
""")
    write(os.path.join(out_dir, "sagehen_gwe.est"), f"""\
BEGIN OPTIONS
  SAVE_FLOWS
  DENSITY_WATER        1000.0
  HEAT_CAPACITY_WATER  4184.0
END OPTIONS

BEGIN GRIDDATA
  POROSITY
    CONSTANT  {POROSITY}
  HEAT_CAPACITY_SOLID
    CONSTANT  {CP_SOLID}
  DENSITY_SOLID
    CONSTANT  {RHO_SOLID}
END GRIDDATA
""")
    # SFR is an advanced package, so SFE carries it; DRN only ever takes water out
    write(os.path.join(out_dir, "sagehen_gwe.ssm"), """\
BEGIN OPTIONS
  SAVE_FLOWS
END OPTIONS

BEGIN SOURCES
  RCHA  AUX  TEMPERATURE
END SOURCES
""")
    with open(os.path.join(out_dir, "sagehen_gwe.sfe"), "w") as f:
        f.write("# Streamflow energy along the SFR reaches, with conduction through the streambed.\n")
        f.write("BEGIN OPTIONS\n  FLOW_PACKAGE_NAME  SFR\n  SAVE_FLOWS\nEND OPTIONS\n\n")
        f.write("BEGIN PACKAGEDATA\n# rno strt ktf rbthcnd\n")
        for n in range(1, nreach + 1):
            f.write(f"  {n:4d} {STRT} {K_BED} {BED_THK}\n")
        f.write("END PACKAGEDATA\n")

    write(os.path.join(out_dir, "sagehen_gwe.ims"), """\
BEGIN OPTIONS
  PRINT_OPTION  SUMMARY
END OPTIONS

BEGIN NONLINEAR
  OUTER_DVCLOSE  1.0e-4
  OUTER_MAXIMUM  50
END NONLINEAR

BEGIN LINEAR
  INNER_MAXIMUM        300
  INNER_DVCLOSE        1.0e-5
  INNER_RCLOSE         1.0e-3
  LINEAR_ACCELERATION  BICGSTAB
END LINEAR
""")
    write(os.path.join(out_dir, "sagehen_gwe.oc"), """\
BEGIN OPTIONS
  BUDGET       FILEOUT  sagehen_gwe.cbc
  TEMPERATURE  FILEOUT  sagehen_gwe.ucn
END OPTIONS

BEGIN PERIOD 1
  SAVE  TEMPERATURE  FREQUENCY 24
  SAVE  BUDGET       FREQUENCY 24
  SAVE  TEMPERATURE  LAST
  SAVE  BUDGET       LAST
END PERIOD
""")


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "ex-gwf-sagehen-gwe")))
