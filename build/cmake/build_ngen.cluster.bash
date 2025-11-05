#!/bin/bash
  
# Build nextgen on HPC, from ngen directory put this one directory up and run this as ../build_ngen.cluster.bash
# required for SUMMA
module load StdEnv/2023
module load gcc/12.3
module load openblas/0.3.24
module load openmpi/4.1.5
module load netcdf-fortran/4.6.1

export FLAGS_OPT="-flto=1;-fuse-linker-plugin"                # -flto=1 is slow to compile, but might want to use
export SUNDIALS_DIR=../../../sundials/instdir/                # will not be used if -DUSE_SUNDIALS=OFF

# required NGEN
module load boost
module load udunits/2.2.28
module load sqlite

# Build SUMMA NGEN below
cmake -B extern/iso_c_fortran_bmi/cmake_build -S extern/iso_c_fortran_bmi
cmake --build extern/iso_c_fortran_bmi/cmake_build --target all

cmake -B extern/summa/cmake_build -S extern/summa -DUSE_NEXTGEN=ON -DUSE_SUNDIALS=OFF -DSPECIFY_LAPACK_LINKS=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build extern/summa/cmake_build --target all -j

cmake -S . -B cmake_build -DBoost_INCLUDE_DIR=/opt/local/libexec/boost/1.81/include -DPython_NumPy_INCLUDE_DIR=/opt/local/bin/python \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo            \
    -DNGEN_IS_MAIN_PROJECT=ON                    \
    -DNGEN_WITH_MPI:BOOL=OFF                     \
    -DNGEN_WITH_NETCDF:BOOL=ON                   \
    -DNGEN_WITH_SQLITE:BOOL=ON                   \
    -DNGEN_WITH_UDUNITS:BOOL=ON                  \
    -DNGEN_WITH_BMI_FORTRAN:BOOL=ON              \
    -DNGEN_WITH_BMI_C:BOOL=ON                    \
    -DNGEN_WITH_PYTHON:BOOL=ON                   \
    -DNGEN_WITH_ROUTING:BOOL=ON                  \
    -DNGEN_WITH_TESTS:BOOL=ON                    \
    -DNGEN_QUIET:BOOL=ON                         \
    -DNGEN_WITH_EXTERN_ALL:BOOL=ON
    
# Comments on above choices, and defaults
#    -DCMAKE_BUILD_TYPE=RelWithDebInfo:  to be able to run in gdb change to -DCMAKE_BUILD_TYPE=Debug
#    -DNGEN_IS_MAIN_PROJECT=ON        :  must be BOOL=ON for DNGEN_WITH_EXTERN_ALL:BOOL=ON
#    -DNGEN_WITH_MPI:BOOL=OFF         :  may want to turn this ON as well as uncommenting "make -j 8 -C cmake_build"
#    -DNGEN_WITH_NETCDF:BOOL=ON       :  must be BOOL=ON to build SUMMA NGEN
#    -DNGEN_WITH_SQLITE:BOOL=ON       :  must be BOOL=ON if planning to use GeoPackages (and not just geojsons)
#    -DNGEN_WITH_UDUNITS:BOOL=ON      :  must be BOOL=ON to build SUMMA NGEN
#    -DNGEN_WITH_BMI_FORTRAN:BOOL=ON  :  must be BOOL=ON to build SUMMA NGEN
#    -DNGEN_WITH_BMI_C:BOOL=ON        :  must be BOOL=ON for DNGEN_WITH_EXTERN_ALL:BOOL=ON
#    -DNGEN_WITH_PYTHON:BOOL=ON       :  must be BOOL=ON for DNGEN_WITH_EXTERN_ALL:BOOL=ON
#    -DNGEN_WITH_ROUTING:BOOL=ON      :  must have DNGEN_WITH_PYTHON:BOOL=ON for this to be ON
#    -DNGEN_WITH_TESTS:BOOL=ON        :  must have DNGEN_WITH_EXTERN_ALL:BOOL=ON for this to be ON
#    -DNGEN_QUIET:BOOL=ON             :  may want turn to this OFF, especially if debugging
#    -DNGEN_WITH_EXTERN_ALL:BOOL=ON   :  these submodules are not used with SUMMA, you may turn this off you don't want to use them

# make -j 8 -C cmake_build    # build w/ 8 parallel jobs, if uncomment then comment the next line and use DNGEN_WITH_MPI:BOOL=ON
make -C cmake_build

