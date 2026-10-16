#!/bin/bash
# ---------------------------------------------------------------------------------------
# Coupled SUMMA / MODFLOW 6 runs on lumped Sagehen Creek, through the serial and the MPI coupler.
#
# The domain is sagehen_calibration's: 9 GRUs, each a land and a stream HRU, on a MODFLOW model with
# RCH, UZF and UZE, SFR, DRN, the base GHB/DRN exchange, EVT and a GWE model.  Ten days from a cold start:
#
#  1. the serial coupler with mizuRoute: its coupled budget must close and show water or heat through
#     every package, and the outlet reach must carry flow;
#  2. the serial coupler without mizuRoute, and the MPI coupler on 2 ranks: their SUMMA output must
#     agree bit for bit.
#
# Needs bin/summa_modflow6_sundials_mizuroute.exe and bin/summa_modflow6_mpi_sundials_mizuroute.exe,
# built with -DUSE_MPI=ON -DUSE_MIZUROUTE=ON -DUSE_MODFLOW6=ON, and a python3 with netCDF4.
#
# Usage:  ./run_test.sh
# ---------------------------------------------------------------------------------------
set -euo pipefail

TEST_DIR="$(cd "$(dirname "$0")" && pwd)"
SUMMA_ROOT="$(cd "${TEST_DIR}/../../../.." && pwd)"
DOMAIN="${TEST_DIR}/../sagehen_calibration/domain"
WORK="${TEST_DIR}/run"
PYTHON=${PYTHON:-python3}
SERIAL_EXE="${SUMMA_ROOT}/bin/summa_modflow6_sundials_mizuroute.exe"
MPI_EXE="${SUMMA_ROOT}/bin/summa_modflow6_mpi_sundials_mizuroute.exe"

for exe in "${SERIAL_EXE}" "${MPI_EXE}"; do
  if [ ! -x "${exe}" ]; then
    echo "ERROR: ${exe} not found"
    echo "       build it with -DUSE_MPI=ON -DUSE_MIZUROUTE=ON -DUSE_MODFLOW6=ON"
    exit 1
  fi
done

echo "======================================================================"
echo "SUMMA coupled run test (lumped Sagehen, serial and MPI couplers)"
echo "======================================================================"

rm -rf "${WORK}"
mkdir -p "${WORK}"
cp -r "${DOMAIN}" "${WORK}/domain"
sed -i.bak "s|^  map_file .*|  map_file           = '${WORK}/domain/hru2cell_map.txt'|" "${WORK}/domain/summa_modflow6.config"
printf 'scalarSoilDrainage | 1\nscalarAquiferStorage | 1\nscalarAquiferBaseflow | 1\n' \
  >> "${WORK}/domain/settings/SUMMA/outputControl.txt"

# run <name> <mizuRoute: yes|no> <command...>: each run in its own copy of the MODFLOW model
run() {
  local name=$1 routed=$2; shift 2
  local dir="${WORK}/${name}"
  mkdir -p "${dir}/output"
  cp -r "${WORK}/domain/mf6" "${dir}/mf6"
  cat > "${dir}/fileManager.txt" <<EOF
controlVersion       'SUMMA_FILE_MANAGER_V3.0.0'
simStartTime         '2017-10-01 01:00'
simEndTime           '2017-10-10 23:00'
tmZoneInfo           'utcTime'
outFilePrefix        '${name}'
settingsPath         '${WORK}/domain/settings/SUMMA/'
forcingPath          '${WORK}/domain/forcing/'
outputPath           '${dir}/output/'
initConditionFile    'coldState.nc'
attributeFile        'attributes.nc'
trialParamFile       'trialParams.nc'
forcingListFile      'forcingFileList.txt'
decisionsFile        'modelDecisions.txt'
outputControlFile    'outputControl.txt'
globalHruParamFile   'localParamInfo.txt'
globalGruParamFile   'basinParamInfo.txt'
vegTableFile         'TBL_VEGPARM.TBL'
soilTableFile        'TBL_SOILPARM.TBL'
generalTableFile     'TBL_GENPARM.TBL'
noahmpTableFile      'TBL_MPTABLE.TBL'
EOF
  local toml=()
  if [ "${routed}" = "yes" ]; then
    cat > "${dir}/mizuroute.toml" <<EOF
[simulation]

use_mizuroute = true

[mizuRoute]

namelist_path = "${WORK}/domain/settings/mizuroute/"
namelist_file = "mizuroute.param"
methods       = "3"
dt            = 3600
write_qbasin  = false
write_Qreach  = true

[hydrofabric]

hfabric_path = "${WORK}/domain/"
hfabric_file = "topology.nc"

dname_seg = "seg"
dname_hru = "hru"

varname_HRUid     = "hruId"
varname_segId     = "segId"
varname_hruSegId  = "hruToSegId"
varname_downSegId = "downSegId"

varname_area   = "area"
varname_slope  = "slope"
varname_length = "length"
EOF
    toml=("${dir}/mizuroute.toml")
  fi
  echo "Running ${name}..."
  local status
  set +e
  (cd "${dir}/mf6" && "$@" "${dir}/fileManager.txt" "${WORK}/domain/summa_modflow6.config" ${toml[@]+"${toml[@]}"}) \
    > "${dir}/run.log" 2>&1
  status=$?
  set -e
  if [ ${status} -ne 0 ]; then
    echo "FAIL: ${name} exited ${status}"
    tail -25 "${dir}/run.log"
    exit 1
  fi
}

run serial_routed yes "${SERIAL_EXE}"
run serial        no  "${SERIAL_EXE}"
run mpi           no  mpirun -np 2 "${MPI_EXE}"
echo

echo "Checking..."
"${PYTHON}" - "${WORK}" <<'PYEOF'
import glob, os, re, sys
import numpy as np
import netCDF4 as nc

work = sys.argv[1]
fail = False


def check(ok, what):
    global fail
    print(f"  {'ok  ' if ok else 'FAIL'} {what}")
    fail |= not ok


# every package the coupler drives moved water or heat
budget = open(os.path.join(work, "serial_routed", "run.log")).read().split("coupled water budget over the run")[-1]


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

# mizuRoute carries the outlet's flow
routed = glob.glob(os.path.join(work, "serial_routed", "output", "*_timestep.nc"))[0]
with nc.Dataset(routed) as d:
    q = np.asarray(d["Q_reach"][:])[:, int(np.argmax(np.asarray(d["seg"][:]) == 9)), 0]
check(bool(np.all(np.isfinite(q)) and q[-1] > 0.0), f"reach 9 carries flow: {q[-1]:.4f} m3 s-1 at the end")


# the MPI coupler reproduces the serial one: every variable of each output frequency, each rank's file
# matched to the serial file's HRUs and GRUs by id
for freq in ("timestep", "day"):
    serial = glob.glob(os.path.join(work, "serial", "output", f"serial_*_{freq}.nc"))
    ranks = sorted(glob.glob(os.path.join(work, "mpi", "output", f"mpi_*_{freq}.nc")))
    check(len(serial) == 1 and len(ranks) == 2, f"serial writes one {freq} file, MPI one per rank")
    if len(serial) != 1 or not ranks:
        continue
    differ = set()
    with nc.Dataset(serial[0]) as a:
        a.set_auto_mask(False)
        hru_a, gru_a = list(np.asarray(a["hruId"][:])), list(np.asarray(a["gruId"][:]))
        for f in ranks:
            with nc.Dataset(f) as b:
                b.set_auto_mask(False)
                rows = {"hru": [hru_a.index(i) for i in np.asarray(b["hruId"][:])],
                        "gru": [gru_a.index(i) for i in np.asarray(b["gruId"][:])]}
                for v in a.variables:
                    va = a[v]
                    if v not in b.variables:
                        differ.add(v)
                        continue
                    x, y = va[:], b[v][:]
                    for k, dim in enumerate(va.dimensions):
                        if dim in rows:
                            x = np.take(x, rows[dim], axis=k)
                    same = x.shape == y.shape and np.array_equal(x, y, equal_nan=x.dtype.kind == "f")
                    if not same:
                        differ.add(v)
    check(not differ, f"MPI matches serial bit for bit in the {freq} output"
          + (f": {', '.join(sorted(differ)[:6])} differ" if differ else ""))

sys.exit(1 if fail else 0)
PYEOF

echo
echo "PASS"
