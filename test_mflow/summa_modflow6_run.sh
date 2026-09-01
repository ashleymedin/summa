#!/bin/bash
# ---------------------------------------------------------------------------------------
# Minimal smoke test for the SUMMA <-> MODFLOW 6 coupler (summa_modflow6).
#
# Layout expected in $CASE_DIR:
#   mfsim.nam                     MODFLOW 6 simulation name file (single GWF model, DIS,
#                                 RCH with READASARRAYS, TDIS TIME_UNITS SECONDS with one
#                                 time step per SUMMA data step, length unit metres)
#   fileManager.txt              SUMMA file manager (its modelDecisions.txt must set
#                                 groundwatr = modflow  and  bcLowrSoiH = presHead)
#   summa_modflow6.config        &coupler namelist, e.g.
#                                   &coupler
#                                     mf6_model_name   = 'MYMODEL'
#                                     rch_package_name = 'RCHA'
#                                     soil_thickness   = 2.0
#                                     map_file         = ''        ! blank => nearest-cell map
#                                     feedback         = .true.
#                                   /
#
# Optional  map_file  (whitespace separated, '#' comments ignored), one line per SUMMA HRU
# in BMI flatten order:
#     iHRU   nPairs   cell_1 w_1   cell_2 w_2  ...
# where cell_k is the row-major horizontal MODFLOW cell index (irow-1)*ncol + icol and the
# weights w_k sum to 1.
#
# Usage:  ./summa_modflow6_run.sh  /path/to/CASE_DIR  [/path/to/summa_modflow6]
# ---------------------------------------------------------------------------------------
set -euo pipefail

CASE_DIR=${1:?"usage: $0 CASE_DIR [summa_modflow6_exe]"}
EXE=${2:-"$(cd "$(dirname "$0")/../bin" 2>/dev/null && pwd)/summa_modflow6"}

[ -x "$EXE" ] || { echo "coupler executable not found/executable: $EXE"; exit 1; }
[ -f "$CASE_DIR/mfsim.nam" ] || { echo "missing $CASE_DIR/mfsim.nam"; exit 1; }
[ -f "$CASE_DIR/fileManager.txt" ] || { echo "missing $CASE_DIR/fileManager.txt"; exit 1; }
[ -f "$CASE_DIR/summa_modflow6.config" ] || { echo "missing $CASE_DIR/summa_modflow6.config"; exit 1; }

# MODFLOW 6 is initialised from mfsim.nam in the working directory
cd "$CASE_DIR"
echo "running: $EXE fileManager.txt summa_modflow6.config   (in $CASE_DIR)"
"$EXE" fileManager.txt summa_modflow6.config

echo "done. check:"
echo "  - the SUMMA output NetCDF (scalarSoilDrainage, scalarAquiferStorage should be ~0)"
echo "  - the MODFLOW 6 listing budget (RCH inflow tracks SUMMA drainage)"
