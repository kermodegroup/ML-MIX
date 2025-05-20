# -*- CMake -*- file for plugin examples.
# The is meant to be used as a template for plugins that are
# distributed independent from the LAMMPS package.
##########################################

cmake_minimum_required(VERSION 3.16)

# enforce out-of-source build
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds are not allowed. You must create and use a build directory. "
      "Please remove CMakeCache.txt and CMakeFiles first.")
endif()

project(plugins VERSION 1.0 LANGUAGES C)

get_filename_component(PLUGINS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)
set(PLUGINS_SOURCE_DIR     ${PLUGINS_DIR}/src)

set(PROJECT_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Modules")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_MODULE_PATH})

include("${PROJECT_MODULE_PATH}/FetchCPM.cmake")
include("${PROJECT_MODULE_PATH}/FindLAMMPSTools.cmake")

# when this file is included as subdirectory in the LAMMPS build, many settings are directly imported
if(LAMMPS_SOURCE_DIR)
    set(LAMMPS_HEADER_DIR ${LAMMPS_SOURCE_DIR})
else()
    # NOTE: the next line should be commented out when used outside of the LAMMPS package
    #get_filename_component(LAMMPS_SOURCE_DIR ${PROJECT_SOURCE_DIR}/../../src ABSOLUTE)
    set(LAMMPS_HEADER_DIR ${LAMMPS_SOURCE_DIR} CACHE PATH "Location of LAMMPS headers")
    if(NOT LAMMPS_HEADER_DIR)
        message(FATAL_ERROR "Must set LAMMPS_HEADER_DIR")
    endif()
    # by default, install into $HOME/.local (not /usr/local),
    # so that no root access (and sudo) is needed
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "Default install path" FORCE)
    endif()
    # ugly hacks for MSVC which by default always reports an old C++ standard in the __cplusplus macro
    # and prints lots of pointless warnings about "unsafe" functions
    if(MSVC)
        if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
            add_compile_options(/Zc:__cplusplus)
            add_compile_options(/wd4244)
            add_compile_options(/wd4267)
            add_compile_options(/wd4250)
            add_compile_options(/EHsc)
        endif()
        add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
    endif()
endif()

# find the build directory where liblammps.so is located
if(NOT LAMMPS_BUILD_DIR)
    # when built inside of the LAMMPS package LAMMPS_BUILD_DIR is the same as the build directory
    get_filename_component(LAMMPS_BUILD_DIR ${CMAKE_BUILD_DIRECTORY} ABSOLUTE)
else()
    # when built outside of the LAMMPS package LAMMPS_BUILD_DIR needs to be specified
    message(STATUS "LAMMPS build directory specified: ${LAMMPS_BUILD_DIR}")
endif()

# C++20 is required
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Need -restrict with Intel compilers
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -restrict")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
include(CheckIncludeFileCXX)
if(NOT LAMMPS_SOURCE_DIR)
  include(LAMMPSInterfaceCXX)
endif()

##########################
# building the plugins

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

add_library(mlmlplugin MODULE ${PLUGINS_SOURCE_DIR}/KOKKOS/mlmlpluginkokkos.cpp ${PLUGINS_SOURCE_DIR}/fix_mlml.cpp ${PLUGINS_SOURCE_DIR}/KOKKOS/fix_mlml_kokkos.cpp)
target_include_directories(mlmlplugin PRIVATE "${PLUGINS_SOURCE_DIR}" "${PLUGINS_SOURCE_DIR}/KOKKOS")
target_link_libraries(mlmlplugin PRIVATE LAMMPS::src LAMMPS::lammps)

if(Kokkos_FOUND)
    target_link_libraries(mlmlplugin PUBLIC Kokkos::src)
    target_compile_options(mlmlplugin PRIVATE ${KokkosCompileOptions})
endif()

set_target_properties(mlmlplugin PROPERTIES PREFIX "" SUFFIX ".so")

# MacOS seems to need this
if(CMAKE_SYSTEM_NAME STREQUAL Darwin)
    set_target_properties(mlmlplugin
        PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    # tell CMake to export all symbols to a .dll on Windows with special case for MinGW cross-compilers
    set_target_properties(mlmlplugin
        PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    if(CMAKE_CROSSCOMPILING)
        set_target_properties(mlmlplugin PROPERTIES LINK_FLAGS "-Wl,--export-all-symbols")
    endif()
else()
  set_target_properties(mlmlplugin PROPERTIES
      LINK_FLAGS "-rdynamic")
endif()

add_custom_target(plugins ALL ${CMAKE_COMMAND} -E echo "Building Plugins"
    DEPENDS mlmlplugin)
