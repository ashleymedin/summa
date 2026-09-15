#!/bin/bash
  
# Build modflow6 on HPC, from modflow6 directory put this one directory up and run this as ../build_mflow.cluster.bash
# Digital Resource Alliance of Canada settings
module load StdEnv/2023
module load gcc/12.3
module load openblas/0.3.24
module load openmpi/4.1.5
module load gdal/3.10.0
module load conda/2024.09
module load netcdf-fortran/4.6.1
#
# Purdue Anvil settings
#module load gcc/14.2.0
#module load openmpi/4.1.6
#module load gdal/3.10.0
#module load conda/2024.09
#module load openblas/0.3.17
#module load netcdf-fortran/4.5.3

# Only do this once
#conda env create -f environment.yml

# Environment variables may be set within this script (see examples below) or in the terminal environment before executing this script
# activate correct python environment, here is an example with conda environment named modflow6
: "${MODFLOW_CONDA_ENV:=modflow6}"
source ${HOME}/PythonEnv/${MODFLOW_CONDA_ENV}/bin/activate
# fallback: allow overriding python executable explicitly
: "${MODFLOW_PYTHON_EXECUTABLE:=$(which python 2>/dev/null || echo /usr/bin/python3)}"
# root of the active python environment (asked of the interpreter itself so a stale
# VIRTUAL_ENV/CONDA_PREFIX can't mislead it); used as a hint for find_package(Python)
: "${MODFLOW_PYTHON_ROOT:=$("${MODFLOW_PYTHON_EXECUTABLE}" -c 'import sys; print(sys.prefix)' 2>/dev/null || dirname "$(dirname "${MODFLOW_PYTHON_EXECUTABLE}")")}"
export PYTHONNOUSERSITE=1
python -m pip install --upgrade "pip<24.1" >/dev/null 2>&1 || true

# Build MODFLOW6
meson setup --prefix=`pwd` --libdir=bin builddir
meson install -C builddir

# Actors install of sundials
#export SUNDIALS_DIR="$CMAKE_PREFIX_PATH:$HOME/Summa-Actors/utils/dependencies/install/sundials/"
#
# Regular install of sundials
export SUNDIALS_DIR=$HOME/SummaSundials/sundials/instdir/

# May want to use this flag
#export FLAGS_OPT="-flto=1;-fuse-linker-plugin"

MF6_BIN="$(cd "$(pwd)/bin" && pwd)"
cmake -B srcextern/summa/build/cmake_build -S srcextern/summa/build \
    -DUSE_MODFLOW6=ON
    -DMODFLOW6_LIB_DIR="${MF6_BIN}" \
    -DUSE_SUNDIALS=ON \
    -DUSE_MPI=OFF \
    -DUSE_NEXTGEN=OFF \
    -DUSE_OPENWQ=OFF \
    -DSPECIFY_LAPACK_LINKS=ON \
    -DCMAKE_BUILD_TYPE=Release
    
cmake --build srcextern/summa/build/cmake_build --target all -j
