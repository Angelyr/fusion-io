cmake_minimum_required(VERSION 3.24)
project(fusion_io VERSION 1.0.0 LANGUAGES C CXX Fortran)

enable_testing()
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(FindPkgConfig)
pkg_check_modules(HDF5 REQUIRED IMPORTED_TARGET hdf5)

option(FUSIONIO_ENABLE_PCMS OFF)
if (FUSIONIO_ENABLE_PCMS)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set(CMAKE_CXX_STANDARD 17)
    find_package(pcms REQUIRED)
    find_package(perfstubs REQUIRED)
    set_property(TARGET Omega_h::omega_h PROPERTY INTERFACE_COMPILE_OPTIONS "")
endif()

option(FUSIONIO_ENABLE_PYTHON "Enable Python" ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "-fdefault-real-8")
endif()

add_subdirectory(m3dc1_lib)
add_subdirectory(fusion_io)

option(FUSIONIO_ENABLE_TRACE "Enable Trace" ON)
if (FUSIONIO_ENABLE_TRACE)
    add_subdirectory(trace)
endif()

option(FUSIONIO_BUILD_TESTING "Build for CTest" OFF)
if (FUSIONIO_BUILD_TESTING)
    add_subdirectory(examples)
endif()