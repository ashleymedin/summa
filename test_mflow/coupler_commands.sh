#!/bin/bash
# ---------------------------------------------------------------------------------------
# Expected in  MODFLOW_CASE  (the coupler runs here; MODFLOW reads mfsim.nam from cwd):
#   mfsim.nam                     MODFLOW 6 simulation name file (single GWF model, DIS,
#                                 RCH with READASARRAYS, TDIS TIME_UNITS SECONDS with one
#                                 time step per SUMMA data step, length unit metres)
#   (optional) the map_file named in summa_modflow6.config
#
# Expected one directory ABOVE MODFLOW_CASE:
#   summa_modflow6.config        &coupler namelist, e.g.
#                                   &coupler
#                                     mf6_model_name   = 'SAGEHEN'
#                                     rch_package_name = 'RCHA'
#                                     soil_thickness   = 2.0
#                                     map_file         = ''        ! blank => nearest-cell map
#                                     feedback         = .true.
#                                   /
#
# Expected in  SUMMA_DIR :
#   fileManager.txt              SUMMA file manager (its modelDecisions.txt must set
#                                 groundwatr = modflow  and  bcLowrSoiH = presHead)
#
# map_file format (whitespace separated, '#' comments ignored), one line per SUMMA HRU
# in BMI flatten order:
#     iHRU   nPairs   cell_1 w_1   cell_2 w_2  ...
# where cell_k is the row-major horizontal MODFLOW cell index (irow-1)*ncol + icol and the
# weights w_k sum to 1.
#
# Usage:  ./coupler_commands.sh  /path/to/MODFLOW_CASE  /path/to/SUMMA_FILEMANAGER  [/path/to/summa_modflow6.exe]
# ---------------------------------------------------------------------------------------
set -euo pipefail

MODFLOW_CASE=${1:?"usage: $0 MODFLOW_CASE SUMMA_FILEMANAGER [summa_modflow6_exe]"}
SUMMA_FILEMANAGER=${2:?"usage: $0 MODFLOW_CASE SUMMA_FILEMANAGER [summa_modflow6_exe]"}
if [ $# -ge 3 ]; then
  EXE=$3
else
  # default: the coupler built into srcextern/summa/bin (may not exist yet)
  BIN_DIR=$(cd "$(dirname "$0")/../bin" 2>/dev/null && pwd || true)
  EXE=${BIN_DIR:+$BIN_DIR/summa_modflow6.exe}
  EXE=${EXE:-"$(dirname "$0")/../bin/summa_modflow6.exe"}
fi

# Resolve to absolute paths so they survive the cd into MODFLOW_CASE.
# summa_modflow6.config lives one directory ABOVE mfsim.nam (i.e. above MODFLOW_CASE).
EXE=$(cd "$(dirname "$EXE")" 2>/dev/null && pwd)/$(basename "$EXE") || true
FILE_MANAGER=$(cd "$(dirname "$SUMMA_FILEMANAGER")" 2>/dev/null && pwd)/$(basename "$SUMMA_FILEMANAGER") || true
CONFIG=$(cd "$MODFLOW_CASE/.." 2>/dev/null && pwd)/summa_modflow6.config || true

[ -x "$EXE" ]                    || { echo "coupler executable not found/executable: $EXE"; exit 1; }
[ -f "$MODFLOW_CASE/mfsim.nam" ] || { echo "missing $MODFLOW_CASE/mfsim.nam"; exit 1; }
[ -f "$CONFIG" ]                 || { echo "missing $MODFLOW_CASE/../summa_modflow6.config"; exit 1; }
[ -f "$FILE_MANAGER" ]           || { echo "missing $SUMMA_FILEMANAGER"; exit 1; }

# MODFLOW 6 is initialized from mfsim.nam in the working directory
cd "$MODFLOW_CASE"
"$EXE" "$FILE_MANAGER" "$CONFIG"

echo "done. check:"
echo "  - the SUMMA output NetCDF scalarSoilDrainage should match the MODFLOW 6 RCH inflow"
echo "  - the SUMMA output NetCDF scalarAquiferStorage should match the MODFLOW 6 GWF storage change"
