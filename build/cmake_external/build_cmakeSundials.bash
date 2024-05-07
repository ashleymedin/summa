#!/bin/bash

# from {$maindir}/sundials/builddir, run
# cp ../../summa/build/cmake_external/build_cmakeSundials.bash build_cmake
# run script from the builddir directory with ./build_cmake
# run `make`, then `make install`
# Note, using -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_ENABLE_F2003=OFF -DEXAMPLES_ENABLE_CUDA=OFF -DEXAMPLES_ENABLE_CXX=OFF, if want to run examples should change

module load StdEnv/2020
module load gcc/9.3.0
module load openblas/0.3.17
module load netcdf-fortran/4.5.2
module load cuda/11.0
module load magma/2.5.4

cmake ../../sundials-software/ -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_ENABLE_F2003=OFF -DEXAMPLES_ENABLE_CUDA=OFF -DEXAMPLES_ENABLE_CXX=OFF -DBUILD_FORTRAN_MODULE_INTERFACE=ON -DCMAKE_Fortran_COMPILER=gfortran -DENABLE_MAGMA=ON -DMAGMA_DIR=$EBROOTMAGMA -DENABLE_CUDA=ON -DCMAKE_INSTALL_PREFIX=../../sundials/instdir -DEXAMPLES_INSTALL_PATH=../../sundials/instdir/examples
