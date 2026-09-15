#!/bin/bash
  
# Build modflow6 on Mac, from modflow6 directory put this one directory up and run this as ../build_mflow.mac.bash

# Only do this once
#conda env create -f environment.yml

# Environment variables may be set within this script (see examples below) or in the terminal environment before executing this script
# activate correct python environment, here is an example with conda environment named modflow6
: "${MODFLOW_CONDA_ENV:=modflow6}"
# try common conda install locations; adjust if your conda is elsewhere
if [ -f "${HOME}/opt/anaconda3/etc/profile.d/conda.sh" ]; then
  . "${HOME}/opt/anaconda3/etc/profile.d/conda.sh"
elif [ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
  . "${HOME}/miniconda3/etc/profile.d/conda.sh"
elif command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)" || true
fi
# activate env if possible (non-fatal)
if command -v conda >/dev/null 2>&1; then
  conda activate "${MODFLOW_CONDA_ENV}" || true
fi
# fallback: allow overriding python executable explicitly
: "${DPython3_EXECUTABLE:=$(which python 2>/dev/null || echo /usr/bin/python3)}"
export DPython3_EXECUTABLE

# Build MODFLOW6
meson setup --prefix=`pwd` --libdir=bin builddir
meson install -C builddir

# Mac Example using MacPorts:
export FC=/opt/local/bin/gfortran                             # Fortran compiler family
#export FLAGS_OPT="-flto=1"                                   # -flto=1 is slow to compile, but might want to use
export LIBRARY_LINKS='-llapack'                               # list of library links
export SUNDIALS_DIR=../../SummaSundials/sundials/instdir/     # will not be used if -DUSE_SUNDIALS=OFF

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
