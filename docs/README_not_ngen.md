# SummaSundials 

## About

This is configured with the [CMakeLists.txt](CMakeLists.txt) and other files in this directory.

#### Extra Outer Directory for Actors

Actors needs one


### Getting the Latest Changes

There are two steps to getting upstream submodule changes fully 
  1. fetching and locally checking out the changes from the remote
  2. committing the new checkout revision for the submodule

To fetch and check out the latest revision (for the [currently used branch](#viewing-the-current-branch)):

    git pull

# Usage

## Building Libraries

### Solver options in SUMMA build scripts
First, cd into the cmake directory in summa, without Actors:

    cd summa/build/cmake

or with Actors:

    cd summa/build/summa/build/cmake

If you want to use Sundials IDA or BE Kinsol, set -DCMAKE_BUILD_TYPE=Sundials* in the build script instead of BE*.  Then, before summa can be built, Sundials needs to be installed.

### Installing SUNDIALS
For reference, let the main summa directory be contained in the folder `top_dir`
    - e.g., summa executables would be located within `top_dir/summa/bin`
    
We maintain a Sundials version that is compatible with running on GPU and using the fortran interface. You need to install it from git: 
    within `top_dir`: 
    $ git clone https://github.com/ashleymedin/sundials sundials-software
    $ cd sundials-software
    $ git checkout -b gpu_fortran remotes/origin/gpu_fortran

Create new empty directories to prep for SUNDIALS installation
    within `top_dir`: 
    $ mkdir sundials
    $ cd sundials
    $ mkdir builddir instdir

Copy CMake build script from SUMMA files to properly configure SUNDIALS
    $ cd buildir
    $ cp ../../summa/build/cmake_external/build_cmakeSundials.bash .

Build SUNDIALS configured for SUMMA
    within `buildir`: 
    $ ./build_cmakeSundials.bash
    $ make
    $ make install
    
Note if you need to recompile after a system upgrade, delete the contents of sundials/instdir and sundials/buildir EXCEPT sundials/buildir/build_cmakeSundials.bash before building and installing.

Note that when there is an existing directory, it may sometimes be necessary to clear it and regenerate, especially if any changes were made to the CMakeLists.txt file.

### Building SUMMA
After there is build system directory, the shared library can be built using the `summa[actors]` CMake target. For example, the SummaSundials shared library file (i.e., the build config's `summa` target) can be built using:

    $ cmake --build ../cmake_build --target all

This will build a `cmake_build/libsumma.<ext>` file, where the extension depends on the local machine's operating system.    

There is an example of a bash script to build the summa libraries at /build/cmake/build[_actors].[system_type].bash. Sundials is turned on here. These need to be run in the cmake directory.

