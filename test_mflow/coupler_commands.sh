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
# map_file format: one "iHRU  cell  weight" triple per line (whitespace separated;
# blank lines and '#' comments ignored).  An HRU may span any number of lines.
# cell is the row-major horizontal MODFLOW index (irow-1)*ncol + icol; weights are
# normalised per HRU, so put e.g. 1.0 on every line to spread an HRU over its cells.
#
# Usage:  ./coupler_commands.sh [-c CONFIG] MODFLOW_CASE SUMMA_FILEMANAGER [summa_modflow6.exe]
#
#   -c, --config CONFIG   path to summa_modflow6.config.  Default: MODFLOW_CASE/../summa_modflow6.config,
#                         so several cases can share one MODFLOW model directory while each keeps its own
#                         coupler settings (model/package names, HRU->cell map_file, feedback).
# ---------------------------------------------------------------------------------------
set -euo pipefail

CONFIG_ARG=""
while [ $# -gt 0 ]; do
  case $1 in
    -c|--config) CONFIG_ARG=${2:?"$0: -c/--config needs a path"}; shift 2 ;;
    --)          shift; break ;;
    -*)          echo "$0: unknown option $1"; exit 1 ;;
    *)           break ;;
  esac
done

MODFLOW_CASE=${1:?"usage: $0 [-c CONFIG] MODFLOW_CASE SUMMA_FILEMANAGER [summa_modflow6_exe]"}
SUMMA_FILEMANAGER=${2:?"usage: $0 [-c CONFIG] MODFLOW_CASE SUMMA_FILEMANAGER [summa_modflow6_exe]"}
if [ $# -ge 3 ]; then
  EXE=$3
else
  # default: the coupler built into srcextern/summa/bin (may not exist yet)
  BIN_DIR=$(cd "$(dirname "$0")/../bin" 2>/dev/null && pwd || true)
  EXE=${BIN_DIR:+$BIN_DIR/summa_modflow6.exe}
  EXE=${EXE:-"$(dirname "$0")/../bin/summa_modflow6.exe"}
fi

# Resolve to absolute paths so they survive the cd into MODFLOW_CASE.
# Without -c, summa_modflow6.config is taken one directory ABOVE mfsim.nam (i.e. above MODFLOW_CASE).
EXE=$(cd "$(dirname "$EXE")" 2>/dev/null && pwd)/$(basename "$EXE") || true
FILE_MANAGER=$(cd "$(dirname "$SUMMA_FILEMANAGER")" 2>/dev/null && pwd)/$(basename "$SUMMA_FILEMANAGER") || true
if [ -n "$CONFIG_ARG" ]; then
  CONFIG=$(cd "$(dirname "$CONFIG_ARG")" 2>/dev/null && pwd)/$(basename "$CONFIG_ARG") || true
  [ -f "$CONFIG" ] || CONFIG=$CONFIG_ARG   # unresolvable: keep what was typed, for the error below
else
  CONFIG=$(cd "$MODFLOW_CASE/.." 2>/dev/null && pwd)/summa_modflow6.config || true
fi

[ -x "$EXE" ]                    || { echo "coupler executable not found/executable: $EXE"; exit 1; }
[ -f "$MODFLOW_CASE/mfsim.nam" ] || { echo "missing $MODFLOW_CASE/mfsim.nam"; exit 1; }
[ -f "$CONFIG" ]                 || { echo "missing coupler config: $CONFIG"; exit 1; }
[ -f "$FILE_MANAGER" ]           || { echo "missing $SUMMA_FILEMANAGER"; exit 1; }

# MODFLOW 6 is initialized from mfsim.nam in the working directory
cd "$MODFLOW_CASE"
"$EXE" "$FILE_MANAGER" "$CONFIG"

echo "done. check (SUMMA output NetCDF vs the MODFLOW 6 listing budget):"
echo "  - scalarSoilDrainage    <-> RCH (RCHA) inflow          [coupler imposes this on MODFLOW]"
echo "  - scalarAquiferBaseflow <-> bflow package (CHD/DRN/..) discharge over the same HRUs"
echo "  - scalarAquiferStorage   = Sy * (MODFLOW water table - soil-column base); a copy for"
echo "                             inspecting the water-table position, not an independent check"
