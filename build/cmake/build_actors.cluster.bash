#!/bin/bash
  
# build on Copernicus or Graham, from cmake directory run this as ./build_actors.cluster.bash
# for Summa
module load StdEnv/2020
module load gcc/9.3.0
module load openblas/0.3.17
module load netcdf-fortran/4.5.2
module load cuda/11.0
module load ginkgo/1.6.0
module load magma/2.5.4

# for Actors
module load caf/0.18.3

export FLAGS_OPT="-flto=1;-fuse-linker-plugin"
export SUNDIALS_PATH=/globalhome/kck540/HPC/Libraries/sundials/v7.0/instdir

cmake -B ../cmake_build -S . -DCMAKE_BUILD_TYPE=Sundials_Actors_Cluster
cmake --build ../cmake_build --target all
