#!/bin/bash

# from {$maindir}/sundials/builddir, run
# cp ../../summa/build/cmake_external/build_cmakeSundials.bash build_cmake
# run script from the builddir directory with ./build_cmake
# run `make`, then `make install`
# Note, using -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_ENABLE_F2003=OFF -DEXAMPLES_ENABLE_CUDA=OFF -DEXAMPLES_ENABLE_CXX=OFF, if want to run examples should change

module load StdEnv/2023
module load gcc/12.3
module load openblas/0.3.24
module load openmpi/4.1.5
module load netcdf-fortran/4.6.1
module load cuda/12.2
module load magma/2.7.2

cmake ../../sundials-software/ -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_ENABLE_F2003=OFF -DEXAMPLES_ENABLE_CUDA=OFF -DEXAMPLES_ENABLE_CXX=OFF -DBUILD_FORTRAN_MODULE_INTERFACE=ON -DCMAKE_Fortran_COMPILER=gfortran -DENABLE_MAGMA=ON -DMAGMA_DIR=$EBROOTMAGMA -DENABLE_CUDA=ON -DCMAKE_INSTALL_PREFIX=../../sundials/instdir -DEXAMPLES_INSTALL_PATH=../../sundials/instdir/examples
