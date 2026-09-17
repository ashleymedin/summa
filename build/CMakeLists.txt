cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project(summa VERSION 4.0.0 LANGUAGES Fortran C CXX)

if(DEFINED ENV{FC})
    set (CMAKE_Fortran_COMPILER $ENV{FC})  #for NextGen need, works for others
else()
    set (CMAKE_Fortran_COMPILER gfortran)
endif()

#=========================================================================================
# SUMMA can be compiled with debug flags by specifiying -DCMAKE_BUILD_TYPE=Debug
#
# There are multiple options for building SUMMA, including using the
# Sundials suite of solvers and the NextGen framework.
# The options are set by specifying -DOPTION=ON when running cmake.
# For example to compile SUMMA with Sundials, use
# cmake -B ../cmake_build -S ../. -DUSE_SUNDIALS=ON
# cmake --build ../cmake_build --target all -j
#
# To compile SUMMA with Actors, see the Summa-Actors repo:
# https://github.com/uofs-simlab/Summa-Actors
#=========================================================================================

# Add options for build type
set(CMAKE_CONFIGURATION_TYPES Release Debug)

# Options: Enable each with cmake -DOPTION=ON
option(USE_MIZUROUTE "Enable mizuRoute coupling" OFF)
option(USE_SUNDIALS "Use IDA solver from SUNDIALS suite" OFF)
option(USE_NEXTGEN "Use NextGen Framework" OFF)
option(USE_OPENWQ "Use OpenWQ Framework" OFF)
option(USE_MODFLOW6 "Build the summa_modflow6 SUMMA/MODFLOW 6 coupler executable" OFF)
# directory that contains the MODFLOW 6 shared library (libmf6); required with USE_MODFLOW6
set(MODFLOW6_LIB_DIR "" CACHE PATH "directory containing the MODFLOW 6 shared library (libmf6)")
option(USE_MPI "Use MPI parallelization" OFF)
option(USE_TOML "Build the TOML configuration reader (needed for -c)" OFF)

# The TOML reader is only wanted by things that actually read TOML: the mizuRoute
# configuration, and the calibration driver's run manifest (built with MPI). A plain run
# driven by a file manager never touches it, so it is off unless something needs it and
# the toml-f submodule is then not required either.
if (USE_MIZUROUTE OR USE_MPI)
    set(USE_TOML ON)
endif()
set(EXT_TARGETS) # list of external targets to link to

# Set top-level directory
if (USE_NEXTGEN)
    set(F_MASTER "${CMAKE_CURRENT_SOURCE_DIR}/summa")
else()
    set(F_MASTER "${CMAKE_CURRENT_SOURCE_DIR}/..")
endif()

set(PARENT_DIR "${F_MASTER}")
set(EXEC_DIR "${F_MASTER}/bin")

# Handle optional external libraries

# SUNDIALS
if (USE_SUNDIALS)
    message("ENABLING SUNDIALS")
    find_package(SUNDIALS REQUIRED)
    list(APPEND EXT_TARGETS SUNDIALS::fida_mod_shared
                            SUNDIALS::fkinsol_mod_shared)
    add_compile_definitions(SUNDIALS_ACTIVE)
endif()

# NextGen
if (USE_NEXTGEN)
    message("ENABLING NEXTGEN")
    add_compile_definitions(NGEN_ACTIVE NGEN_FORCING_ACTIVE
                            NGEN_OUTPUT_ACTIVE)
    add_subdirectory(../iso_c_fortran_bmi ${CMAKE_BINARY_DIR}/iso_c_bmi)
    list(INSERT CMAKE_MODULE_PATH 0 "${F_MASTER}/build/cmake/")
    list(INSERT CMAKE_PREFIX_PATH 0 "${F_MASTER}/build/cmake/")
    set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran)
    set(SUMMA_LIB_NAME_CMAKE summabmi)
    set(SUMMA_LIB_DESC_CMAKE "Summa-Sundials BMI Module Shared Library")
else()
    LIST(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_SOURCE_DIR}/cmake/")
    LIST(INSERT CMAKE_PREFIX_PATH 0 "${CMAKE_SOURCE_DIR}/cmake/")
endif()

# MODFLOW 6 coupler: builds the standard (non-NextGen) SUMMA shared library together with
# the SUMMA BMI module (summa_bmi.f90) and links the executable against libmf6.
if (USE_MODFLOW6)
    if (USE_NEXTGEN)
        message(FATAL_ERROR "USE_MODFLOW6 and USE_NEXTGEN are mutually exclusive")
    endif()
    message("ENABLING MODFLOW 6 COUPLER")
    add_compile_definitions(MODFLOW_ACTIVE)
    set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran)
    # The SUMMA BMI module (summa_bmi.f90, non-NextGen path) uses the CSDMS bmi-fortran
    # bindings, i.e. the module "bmif_2_0" from library target "bmif". Provide it either
    # from a sibling checkout of the NOAA-OWP iso_c_fortran_bmi repo (which vendors it, as
    # the NextGen build expects) or by fetching bmi-fortran directly.
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../iso_c_fortran_bmi/CMakeLists.txt")
        add_subdirectory(../iso_c_fortran_bmi ${CMAKE_BINARY_DIR}/iso_c_bmi)
    endif()
    if (NOT TARGET bmif)
        # Fetch the CSDMS bmi-fortran sources only; do NOT let its CMakeLists.txt
        # configure (as a subdirectory it names its targets "${CMAKE_PROJECT_NAME}_*"
        # which resolves to "summa_*", so no "bmif" target would ever appear).
        # SOURCE_SUBDIR points at a non-existent dir so MakeAvailable skips add_subdirectory.
        include(FetchContent)
        FetchContent_Declare(bmi_fortran
            GIT_REPOSITORY https://github.com/csdms/bmi-fortran.git
            GIT_TAG        v2.0
            SOURCE_SUBDIR  do-not-configure)
        FetchContent_MakeAvailable(bmi_fortran)
        # Build the single-source module ourselves as target "bmif" (module bmif_2_0).
        # CMAKE_Fortran_MODULE_DIRECTORY (set above) puts bmif_2_0.mod where summa looks.
        add_library(bmif STATIC ${bmi_fortran_SOURCE_DIR}/bmi.f90)
        set_target_properties(bmif PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif()
    if (NOT TARGET bmif)
        message(FATAL_ERROR "USE_MODFLOW6: could not provide the bmi-fortran 'bmif' target "
                            "(module bmif_2_0). Check out iso_c_fortran_bmi next to the summa "
                            "repo or make bmi-fortran available.")
    endif()
endif()

# OpenWQ Framework
if (USE_OPENWQ)
    message("ENABLING OPENWQ")
    add_compile_definitions(OPENWQ_ACTIVE)
endif()

# MPI
if (USE_MPI)
    message("ENABLING MPI")

    find_package(MPI REQUIRED COMPONENTS Fortran)

    add_library(PARALLEL_UTILS STATIC
        source/parallel/mpi_context.f90
        source/parallel/error_utils.f90
    )

    target_link_libraries(PARALLEL_UTILS PUBLIC
        SUMMA_NRTYPE
        MPI::MPI_Fortran
    )
endif()

# TOML parser
# NOTE: paths here are relative to F_MASTER, not to this file, because this same file is
#       also used as ngen/extern/summa/CMakeLists.txt, where the SUMMA source sits one
#       directory deeper. F_MASTER accounts for that; see above.
if (USE_TOML)
    message("ENABLING TOML CONFIGURATION READER")
    add_compile_definitions(TOML_ACTIVE)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    add_subdirectory(
        ${F_MASTER}/external/toml-f
        ${CMAKE_BINARY_DIR}/toml-f
    )
endif()

# mizuRoute
if (USE_MIZUROUTE)

    message("ENABLING MIZUROUTE")

    set(MIZU_DIR
        ${F_MASTER}/external/mizuRoute/route/build/src
    )

    set(SUMMA_MIZU_DIR
        ${F_MASTER}/build/source/mizuroute
    )

    include(${SUMMA_MIZU_DIR}/CMakeLists.txt)
endif()

# Executable names
#
# The name records which options the binary was built with, so that builds with
# different options can share bin/ without overwriting one another. Note that
# CMAKE_BUILD_TYPE is not part of the name: a Debug and a Release build of the
# same options still collide, so give them separate output directories if you
# need both at once.
set(EXEC_SUFFIX "")
if (USE_SUNDIALS)
    string(APPEND EXEC_SUFFIX "_sundials")
endif()
if (USE_OPENWQ)
    string(APPEND EXEC_SUFFIX "_openwq")
endif()
if (USE_MIZUROUTE)
    string(APPEND EXEC_SUFFIX "_mizuroute")
endif()

set(EXEC_NAME     "summa${EXEC_SUFFIX}.exe")
set(EXEC_NAME_MPI "summa${EXEC_SUFFIX}_mpi.exe")
set(EXEC_NAME_OPT "summa${EXEC_SUFFIX}_opt.exe")

message("Executable name: ${EXEC_NAME}")
if (USE_MPI)
    message("MPI executable name: ${EXEC_NAME_MPI}")
    message("calibration executable name: ${EXEC_NAME_OPT}")
endif()

get_filename_component(F_MASTER "${F_MASTER}" REALPATH)
get_filename_component(PARENT_DIR "${PARENT_DIR}" REALPATH)
get_filename_component(EXEC_DIR "${EXEC_DIR}" REALPATH)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXEC_DIR})

#=========================================================================================
# SET THE PATHS TO THE LIBRARIES AND INCLUDE FILES
#=========================================================================================

# NetCDF is found with a custom FindNetCDF.cmake file
find_package(NetCDF REQUIRED)
list(APPEND EXT_TARGETS NetCDF::NetCDF)

# Find LAPACK with user-specified links or automatic detection
option(SPECIFY_LAPACK_LINKS "Use specified LAPACK links" OFF)
if(SPECIFY_LAPACK_LINKS)
	message("Using LAPACK links from LIBRARY_LINKS environment variable")
	if(DEFINED ENV{LIBRARY_LINKS})
		set(LIBRARY_LINKS $ENV{LIBRARY_LINKS})
	else()
		message(FATAL_ERROR "LIBRARY_LINKS environment variable not found. Set this variable or try -DSPECIFY_LAPACK_LINKS=OFF in build file.")
	endif()
else()
	message("Automatically detect LAPACK links")
	set(LIBRARY_LINKS "")

	# OpenBLAS
	set(BLA_VENDOR OpenBLAS)
	find_package(OpenBLAS REQUIRED COMPONENTS serial)
	list(APPEND EXT_TARGETS OpenBLAS::OpenBLAS)
endif()

# Set compiler flags
set(FLAGS_OPT $ENV{FLAGS_OPT}) # get optional user-specified flags from environment variables
if(CMAKE_BUILD_TYPE MATCHES Debug)
	message("\nSetting SUMMA Debug Options")
	# Notes: - optimization -Og is selected here to maximize information displayed in debugging software (e.g., gdb)
	#        - however, it may be useful to vary the optimization level to test for optimization related issues
    set(FLAGS_NOAH -g -Og -fbacktrace -fbounds-check -ffree-form -ffree-line-length-none -fmax-errors=0 -fPIC -Wfatal-errors ${FLAGS_OPT})
    set(FLAGS_ALL -g -Og -fbacktrace -fbounds-check -ffree-line-length-none -fmax-errors=0 -fPIC -Wfatal-errors -cpp ${FLAGS_OPT})
    set(FLAGS_CXX -g -Og -fbounds-check -Wfatal-errors -std=c++17 ${FLAGS_OPT})
else()
	message("\nSetting SUMMA Release Options")
    set(FLAGS_NOAH -O3 -ffree-form -ffree-line-length-none -fmax-errors=0 -fPIC -Wfatal-errors ${FLAGS_OPT})
    set(FLAGS_ALL -O3 -ffree-line-length-none -fmax-errors=0 -fPIC -Wfatal-errors -cpp ${FLAGS_OPT})
    set(FLAGS_CXX -O3 -Wfatal-errors -std=c++17 ${FLAGS_OPT})
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # Suppress libgfortran end-of-run IEEE summaries (e.g., IEEE_DENORMAL).
    list(APPEND FLAGS_NOAH -ffpe-summary=none)
    list(APPEND FLAGS_ALL -ffpe-summary=none)
endif()

#=========================================================================================
# COMPILE PART 1: Define directory paths
#=========================================================================================

# Define directories that contains source code
set(DRIVER_DIR ${F_MASTER}/build/source/driver)
set(OBJFNC_DIR ${F_MASTER}/build/source/objfunc)
set(DSHARE_DIR ${F_MASTER}/build/source/dshare)
set(CALIBR_DIR ${F_MASTER}/build/source/calibration)
set(ENGINE_DIR ${F_MASTER}/build/source/engine)
set(HOOKUP_DIR ${F_MASTER}/build/source/hookup)
set(NETCDF_DIR ${F_MASTER}/build/source/netcdf)
set(NOAHMP_DIR ${F_MASTER}/build/source/noah-mp)
set(NUMREC_DIR ${F_MASTER}/build/source/numrec)

#=========================================================================================
# COMPILE PART 2: Assemble all of the SUMMA sub-routines
#=========================================================================================

# SUMMA Source Files are defined in the CMakeLists.txt file in the subdirectory
add_subdirectory(${F_MASTER}/build/source/)

# Noah-MP files use free-form syntax even though they use .F extension.
set_source_files_properties(${NOAHMP} PROPERTIES Fortran_FORMAT FREE)


#=========================================================================================
# COMPILE PART 3: Collect the subroutines into build groups depending on build type
#=========================================================================================

set(COMM_ALL ${NRPROC} ${HOOKUP} ${DATAMS} ${UTILMS})
set(SUMMA_ALL ${NETCDF} ${PRELIM} ${MODRUN} ${SOLVER} ${OBJFNC} ${CALIBR} ${DRIVER})

# the TOML configuration reader, only when one is being built
if (USE_TOML)
    set(SUMMA_ALL ${SUMMA_ALL} ${DRIVER_TOML})
endif()
set(MAIN_SUMMA ${DRIVER_DIR}/summa_driver.f90)

# Numerical type definitions used throughout SUMMA
add_library(SUMMA_NRTYPE STATIC ${NRTYPE})
target_compile_options(SUMMA_NRTYPE PRIVATE ${FLAGS_ALL})
set_target_properties(SUMMA_NRTYPE PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/summa_nrtype_mod)
target_include_directories(SUMMA_NRTYPE PUBLIC ${CMAKE_BINARY_DIR}/summa_nrtype_mod)

# Core SUMMA data type definitions
add_library(SUMMA_DATA_TYPES STATIC ${SUMMA_DATA_TYPES})
target_compile_options(SUMMA_DATA_TYPES PRIVATE ${FLAGS_ALL})
target_link_libraries(SUMMA_DATA_TYPES PUBLIC SUMMA_NRTYPE)
set_target_properties(SUMMA_DATA_TYPES PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/summa_data_types_mod)
target_include_directories(SUMMA_DATA_TYPES PUBLIC ${CMAKE_BINARY_DIR}/summa_data_types_mod)

if (USE_MPI)
    set(MAIN_SUMMA_MPI ${DRIVER_DIR}/summa_driver_mpi.f90)
    set(MAIN_SUMMA_OPT ${DRIVER_DIR}/summa_driver_opt.f90)

endif()

if (USE_NEXTGEN)
    set(SUMMA_ALL ${SUMMA_ALL} ${DRIVER_BMI})
endif()

if (USE_MODFLOW6)
    # include the SUMMA BMI module and the MODFLOW 6 coupling modules in the library, and
    # build the coupler as the main program
    set(SUMMA_ALL ${SUMMA_ALL} ${DRIVER_BMI} ${DRIVER_MODFLOW})
    set(MAIN_SUMMA ${DRIVER_DIR}/summa_modflow6.f90)
    set(EXEC_NAME "summa_modflow6${EXEC_SUFFIX}.exe")
    if (USE_MPI)
        set(MAIN_SUMMA_MODFLOW_MPI ${DRIVER_DIR}/summa_modflow6_mpi.f90)
        set(EXEC_NAME_MODFLOW_MPI "summa_modflow6_mpi${EXEC_SUFFIX}.exe")
        # The calibration driver is the same program either way: it evaluates each parameter
        # sample as a coupled SUMMA/MODFLOW 6 simulation when simulation.use_modflow is set in
        # the TOML configuration, and as a plain SUMMA run when it is not.  So this build's
        # calibration executable is simply renamed, the way EXEC_NAME is above, rather than
        # built a second time (two targets over one set of sources would also race for the
        # shared Fortran module directory).
        set(EXEC_NAME_OPT "summa_modflow6_opt${EXEC_SUFFIX}.exe")
    endif()
endif()

message("Executable name: ${EXEC_NAME}")
if (USE_MPI)
    message("MPI executable name: ${EXEC_NAME_MPI}")
    if (USE_MODFLOW6)
        message("MPI MODFLOW 6 coupler executable name: ${EXEC_NAME_MODFLOW_MPI}")
        message("MODFLOW 6 calibration executable name: ${EXEC_NAME_OPT}")
    endif()
endif()

if (USE_SUNDIALS)
    set(COMM_ALL ${COMM_ALL} ${DATAMS_SUNDIALS} ${UTILMS_SUNDIALS})
    set(SUMMA_ALL ${SUMMA_ALL} ${MODRUN_SUNDIALS} ${SOLVER_SUNDIALS})
endif()

if (USE_MIZUROUTE)

    # the TOML configuration reader is only built with mizuRoute, its only consumer

    # mizuRoute type and interface definitions
    add_library(MIZUROUTE_TYPES STATIC ${MIZU_TYPE})
    target_compile_options(MIZUROUTE_TYPES PRIVATE ${FLAGS_ALL})
    target_link_libraries(MIZUROUTE_TYPES PUBLIC SUMMA_NRTYPE)
    set_target_properties(MIZUROUTE_TYPES PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mizuroute_types_mod)
    target_include_directories(MIZUROUTE_TYPES PUBLIC ${CMAKE_BINARY_DIR}/mizuroute_types_mod)

    # Core mizuRoute implementation
    add_library(MIZUROUTE STATIC ${MIZU_BASE} ${MIZU_SHIM} ${MIZU_TOPO}
      ${MIZU_INIT} ${MIZU_RMAP} ${MIZU_ROUTE} ${MIZU_OUTPUT})
    target_compile_options(MIZUROUTE PRIVATE ${FLAGS_ALL})
    target_link_libraries(MIZUROUTE PUBLIC MIZUROUTE_TYPES NetCDF::NetCDF)
    set_target_properties(MIZUROUTE PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mizuroute_mod)

    # SUMMA-facing mizuRoute coupling layer
    add_library(MIZUROUTE_COUPLING STATIC ${MIZU_COUPLING})
    target_compile_options(MIZUROUTE_COUPLING PRIVATE ${FLAGS_ALL})
    target_link_libraries(MIZUROUTE_COUPLING PUBLIC MIZUROUTE)
    add_dependencies(MIZUROUTE_COUPLING SUMMA_TYPES toml-f-lib)
    set_target_properties(MIZUROUTE_COUPLING PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mizuroute_coupling_mod)
    target_include_directories(MIZUROUTE_COUPLING PRIVATE ${CMAKE_BINARY_DIR}/summa_types_mod ${CMAKE_BINARY_DIR}/mizuroute_mod ${CMAKE_BINARY_DIR}/toml-f/include)
    target_include_directories(MIZUROUTE_COUPLING PUBLIC ${CMAKE_BINARY_DIR}/mizuroute_coupling_mod)

endif()

# OpenWQ coupling
if (USE_OPENWQ)
    add_subdirectory(${F_MASTER}/build/source/openwq)
    if (DEFINED OPENWQ_ERROR)
        message(FATAL_ERROR "OpenWQ build failed")
    endif()
    set(SUMMA_ALL ${SUMMA_ALL} ${OPENWQ_COUPLER})
    list(APPEND EXT_TARGETS openWQ)
endif()

# Define version number
set(VERSIONFILE     ${DRIVER_DIR}/summaversion.inc)
execute_process(
    COMMAND bash -c "git tag | tail -n 1"
    OUTPUT_VARIABLE VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND date
    OUTPUT_VARIABLE BUILDTIM
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND bash -c "git describe --long --all --always | sed 's#heads/##'"
    OUTPUT_VARIABLE GITBRCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND git rev-parse HEAD
    OUTPUT_VARIABLE GITHASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

#=========================================================================================
# COMPILE PART 4: Do the compilation
#=========================================================================================

# update version information
file(WRITE  ${VERSIONFILE} "character(len=64), parameter     :: summaVersion = '${VERSION}'\n")
file(APPEND ${VERSIONFILE} "character(len=64), parameter     :: buildTime = '${BUILDTIM}'\n")
file(APPEND ${VERSIONFILE} "character(len=64), parameter     :: gitBranch = '${GITBRCH}'\n")
file(APPEND ${VERSIONFILE} "character(len=64), parameter     :: gitHash = '${GITHASH}'")

# Build SUMMA_NOAHMP Object
add_library(SUMMA_NOAHMP OBJECT ${NOAHMP} ${NRUTIL})
target_compile_options(SUMMA_NOAHMP PRIVATE ${FLAGS_NOAH})
set_target_properties(SUMMA_NOAHMP PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/summa_noahmp_mod)
target_include_directories(SUMMA_NOAHMP PUBLIC ${CMAKE_BINARY_DIR}/summa_noahmp_mod)
target_link_libraries(SUMMA_NOAHMP PUBLIC SUMMA_NRTYPE)
set_target_properties(SUMMA_NOAHMP PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/summa_noahmp_mod)
target_include_directories(SUMMA_NOAHMP PUBLIC ${CMAKE_BINARY_DIR}/summa_noahmp_mod)

# Build SUMMA_COMM Object
add_library(SUMMA_COMM OBJECT ${COMM_ALL})
target_compile_options(SUMMA_COMM PRIVATE ${FLAGS_ALL})
set_target_properties(SUMMA_COMM PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/summa_comm_mod)
target_include_directories(SUMMA_COMM PUBLIC ${CMAKE_BINARY_DIR}/summa_comm_mod)
target_link_libraries(SUMMA_COMM PUBLIC SUMMA_NOAHMP SUMMA_DATA_TYPES ${EXT_TARGETS} ${LIBRARY_LINKS})

# Main SUMMA model type used by the driver and coupling layers
add_library(SUMMA_TYPES STATIC ${SUMMA_TYPE})
target_compile_options(SUMMA_TYPES PRIVATE ${FLAGS_ALL})
target_link_libraries(SUMMA_TYPES PUBLIC SUMMA_COMM)
set_target_properties(SUMMA_TYPES PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/summa_types_mod)
target_include_directories(SUMMA_TYPES PUBLIC ${CMAKE_BINARY_DIR}/summa_types_mod)

if (USE_MIZUROUTE)
    target_compile_definitions(SUMMA_TYPES PRIVATE MIZUROUTE_ACTIVE)
    target_compile_definitions(SUMMA_COMM PRIVATE MIZUROUTE_ACTIVE)
    target_link_libraries(SUMMA_TYPES PUBLIC MIZUROUTE_TYPES)
endif()

# For NextGen, build SUMMA Shared Library and add the outside BMI libraries
if(USE_NEXTGEN)

  if(WIN32)
		add_library(summabmi ${SUMMA_ALL})
	else()
		add_library(summabmi SHARED ${SUMMA_ALL})
	endif()
	target_compile_options(summabmi PRIVATE ${FLAGS_ALL})
	target_link_libraries(summabmi PUBLIC SUMMA_TYPES ${EXT_TARGETS} iso_c_bmi SUMMA_NOAHMP SUMMA_COMM)

	set_target_properties(summabmi PROPERTIES VERSION ${PROJECT_VERSION})
	include(GNUInstallDirs)

	install(TARGETS summabmi
                LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
	configure_file(summabmi.pc.in summabmi.pc @ONLY)
	install(FILES ${CMAKE_BINARY_DIR}/summabmi.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

else()

  # SUMMA shared library
  add_library(summa SHARED ${SUMMA_ALL})
  target_compile_options(summa PRIVATE ${FLAGS_ALL})
  target_link_libraries(summa PUBLIC SUMMA_TYPES ${EXT_TARGETS} SUMMA_NOAHMP SUMMA_COMM ${LIBRARY_LINKS})
  if (USE_MODFLOW6)
    # summa_bmi.f90 (in the library) needs the CSDMS bmi-fortran module bmif_2_0, and
    # mf6_coupling.f90 (also in the library, so that summa_simulation.f90 can run a coupled
    # simulation per calibration sample) calls libmf6 directly.  Both are PUBLIC so that every
    # executable linking summa resolves them.
    if (NOT MODFLOW6_LIB_DIR)
      message(FATAL_ERROR "USE_MODFLOW6=ON requires -DMODFLOW6_LIB_DIR=<dir containing libmf6>")
    endif()
    find_library(MF6_LIB NAMES mf6 libmf6 PATHS "${MODFLOW6_LIB_DIR}" NO_DEFAULT_PATH)
    if (NOT MF6_LIB)
      message(FATAL_ERROR "could not find libmf6 in ${MODFLOW6_LIB_DIR}")
    endif()
    target_link_libraries(summa PUBLIC bmif "${MF6_LIB}")
    # so anything linking the library finds libmf6 at run time without LD_LIBRARY_PATH
    set_target_properties(summa PROPERTIES
      BUILD_RPATH "${MODFLOW6_LIB_DIR}"
      INSTALL_RPATH "${MODFLOW6_LIB_DIR}")
  endif()

  if (USE_TOML)
    target_link_libraries(summa PRIVATE toml-f-lib)
    target_include_directories(summa PRIVATE ${CMAKE_BINARY_DIR}/toml-f/include)
  endif()

  # Optional mizuRoute coupling
  if (USE_MIZUROUTE)
    target_compile_definitions(summa PRIVATE MIZUROUTE_ACTIVE)
    target_link_libraries(summa PRIVATE MIZUROUTE_COUPLING)
  endif()

  # Serial SUMMA executable
  add_executable(${EXEC_NAME} ${MAIN_SUMMA})
  target_compile_options(${EXEC_NAME} PRIVATE ${FLAGS_ALL})
  set_property(TARGET ${EXEC_NAME} PROPERTY LINKER_LANGUAGE Fortran)
  target_link_libraries(${EXEC_NAME} PRIVATE summa)
  if (USE_MODFLOW6)
    # bmif and libmf6 come through summa (PUBLIC, see above); the rpath does not, so set it here
    set_target_properties(${EXEC_NAME} PROPERTIES
      BUILD_RPATH "${MODFLOW6_LIB_DIR}"
      INSTALL_RPATH "${MODFLOW6_LIB_DIR}")
  endif()

  # MPI SUMMA executable
  if (USE_MPI)
    add_executable(${EXEC_NAME_MPI} ${MAIN_SUMMA_MPI})
    target_compile_options(${EXEC_NAME_MPI} PRIVATE ${FLAGS_ALL})
    set_property(TARGET ${EXEC_NAME_MPI} PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(${EXEC_NAME_MPI} PRIVATE summa PARALLEL_UTILS)
  endif()

  # Parallel optimization SUMMA executable.  With USE_MODFLOW6 this is the coupled calibration
  # driver (renamed above); each parameter sample is then evaluated as a coupled SUMMA/MODFLOW 6
  # simulation whenever simulation.use_modflow is set in the TOML configuration.
  if (USE_MPI)
    add_executable(${EXEC_NAME_OPT} ${CALIBR_MPI} ${MAIN_SUMMA_OPT})
    target_compile_options(${EXEC_NAME_OPT} PRIVATE ${FLAGS_ALL})
    set_property(TARGET ${EXEC_NAME_OPT} PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(${EXEC_NAME_OPT} PRIVATE summa PARALLEL_UTILS)
    if (USE_MODFLOW6)
      # so the executable can find libmf6 at run time without LD_LIBRARY_PATH
      set_target_properties(${EXEC_NAME_OPT} PROPERTIES
        BUILD_RPATH "${MODFLOW6_LIB_DIR}"
        INSTALL_RPATH "${MODFLOW6_LIB_DIR}")
    endif()
  endif()

  # MPI MODFLOW 6 coupler executable
  if (USE_MPI AND USE_MODFLOW6)
    add_executable(${EXEC_NAME_MODFLOW_MPI} ${MAIN_SUMMA_MODFLOW_MPI})
    target_compile_options(${EXEC_NAME_MODFLOW_MPI} PRIVATE ${FLAGS_ALL})
    set_property(TARGET ${EXEC_NAME_MODFLOW_MPI} PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(${EXEC_NAME_MODFLOW_MPI} PRIVATE summa PARALLEL_UTILS)
    # so the executable can find libmf6 at run time without LD_LIBRARY_PATH
    set_target_properties(${EXEC_NAME_MODFLOW_MPI} PROPERTIES
      BUILD_RPATH "${MODFLOW6_LIB_DIR}"
      INSTALL_RPATH "${MODFLOW6_LIB_DIR}")
  endif()

endif()
