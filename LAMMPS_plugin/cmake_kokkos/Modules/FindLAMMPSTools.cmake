# Global requirements
# -------------------

find_package(Git REQUIRED)


# Global variables
# ----------------

set(LAMMPS_URL "https://github.com/lammps/lammps.git")

if(NOT LAMMPS_ROOT)
    if(DEFINED ENV{LAMMPS_ROOT})
        set(LAMMPS_ROOT $ENV{LAMMPS_ROOT})
    elseif(CMAKE_PREFIX_PATH)
        find_path(LAMMPS_ROOT
            NAMES LAMMPS_Targets.cmake
            HINTS ${CMAKE_PREFIX_PATH}
            PATH_SUFFIXES LAMMPS
        )
        if("${LAMMPS_ROOT}" STREQUAL "LAMMPS_ROOT-NOTFOUND")
            message(FATAL_ERROR
                "Unable to find LAMMPS. Try setting the CMake "
                "variables LAMMPS_ROOT or CMAKE_PREFIX_PATH."
            )
        endif()
    endif()
endif()


# Utility functions
# -----------------

#     fetch_lammps
#
# Given a LAMMPS `tag` use CPM to retrieve the code for that release.
function(fetch_lammps tag)
    # We use lowercase lammps to avoid clashes with the main library
    CPMAddPackage(NAME lammps
        GIT_REPOSITORY  ${LAMMPS_URL}
        GIT_TAG         ${tag}
        GIT_SHALLOW     TRUE
        DOWNLOAD_ONLY   TRUE
    )
    if(lammps_ADDED)
        set(lammps_SOURCE_DIR ${lammps_SOURCE_DIR} PARENT_SCOPE)
    else()
        message(WARNING "Failed to download LAMMPS source")
    endif()
endfunction()

#     find_executable
#
# Looks for the IMPORTED_LOCATION of the first available IMPORTED_LOCATION_<config>
# for the given target and sets a variable `varname` on the parent scope.
function(find_executable target varname)
    get_target_property(var ${target} IMPORTED_LOCATION)
    if("${var}" STREQUAL "var-NOTFOUND")
        get_target_property(configs ${target} IMPORTED_CONFIGURATIONS)
        list(GET configs 0 config)
        get_target_property(var ${target} "IMPORTED_LOCATION_${config}")
    endif()
    set(${varname} ${var} PARENT_SCOPE)
endfunction()

#     find_lammps_cxx_compiler(path)
#
# Look for the  `nvcc_wrapper` shipped by lammps within the specified `path` and,
# if found, set its location to `LAMMPS_CXX_COMPILER` in the parent scope.
function(find_lammps_cxx_compiler path)
    get_filename_component(NVCC_WRAPPER "${path}/lib/kokkos/bin/nvcc_wrapper" ABSOLUTE)
    if(EXISTS ${NVCC_WRAPPER})
        set(LAMMPS_CXX_COMPILER ${NVCC_WRAPPER} PARENT_SCOPE)
    endif()
endfunction()

#     get_lammps_tag(version)
#
# Given a LAMMPS `version` as reported to CMake or Python, sets `LAMMPS_tag` in the
# parent scope to the latest git tag matching this version within the LAMMPS repo.
function(get_lammps_tag version)
    # Try to get the git tag or commit directly from LAMMPS' help
    execute_process(
        COMMAND ${LAMMPS_EXECUTABLE} -h
        RESULT_VARIABLE exit_code
        OUTPUT_VARIABLE LAMMPS_help
        ERROR_QUIET
    )

    if(exit_code EQUAL 0)
        string(REGEX MATCH "Git info [^ ]+ / ([^)]+)" _ "${LAMMPS_help}")

        if(
            (NOT ("${CMAKE_MATCH_1}" STREQUAL "")) AND
            (NOT ("${CMAKE_MATCH_1}" STREQUAL "(unknown)"))
        )
            set(LAMMPS_tag "${CMAKE_MATCH_1}" PARENT_SCOPE)
            return()
        endif()
    endif()

    # If we're unable to find it we search for the last tag that matches
    # the provided `version`.
    set(MONTHS _ Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec)

    string(REGEX MATCH "([0-9][0-9][0-9][0-9])([0-9][0-9])([0-9][0-9]).*" _ "${version}")
    set(year ${CMAKE_MATCH_1})
    list(GET MONTHS ${CMAKE_MATCH_2} month)
    math(EXPR day ${CMAKE_MATCH_3})
    set(pattern "*${day}${month}${year}*")

    execute_process(
        COMMAND ${GIT_EXECUTABLE} ls-remote --tags --refs ${LAMMPS_URL} ${pattern}
        RESULT_VARIABLE exit_code  # TODO: use this to check if it fails and inform the user
        OUTPUT_VARIABLE git_tags
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    string(REGEX MATCH ".*/(.+)$" _ "${git_tags}")

    set(LAMMPS_tag "${CMAKE_MATCH_1}" PARENT_SCOPE)
endfunction()

#     get_lammps_version(path)
#
# Tries loading the LAMMPSConfigVersion.cmake from within `path`,
# and sets LAMMPS_VERSION in the parent scope.
function(get_lammps_version path)
    if(EXISTS "${path}/lib64/cmake/LAMMPS/LAMMPSConfigVersion.cmake")
        include("${path}/lib64/cmake/LAMMPS/LAMMPSConfigVersion.cmake")
    endif()
    set(LAMMPS_VERSION "${PACKAGE_VERSION}" PARENT_SCOPE)
endfunction()

#     find_lammps()
#
# Macro equivalent to find_package(LAMMPS QUIET), but avoids looking for MPI,
# which is performed separately. It also looks for the NVCC wrapper that
# comes with LAMMPS and tries to find an appropriate git tag matching the
# LAMMPS version.
macro(find_lammps)
    include("${LAMMPS_ROOT}/lib64/cmake/LAMMPS/LAMMPS_Targets.cmake")
    find_executable(LAMMPS::lmp "LAMMPS_EXECUTABLE")
    find_lammps_cxx_compiler("${LAMMPS_EXECUTABLE}/..")
    get_lammps_version(${LAMMPS_ROOT})
    get_lammps_tag(${LAMMPS_VERSION})
endmacro()

#     append_paths(list)
#
# Given a list and any number of file names, appends to the list
# the absolute paths to the provided files.
function(append_paths list)
    list(POP_FRONT ARGV)
    foreach(file ${ARGV})
        get_filename_component(path ${file} ABSOLUTE)
        set(${list} ${${list}} ${path})
    endforeach()
    set(${list} ${${list}} PARENT_SCOPE)
endfunction()

#    copy_target_property
#
# Given a source target `src` and a destination target `dst`, reads the value of the
# `property` in `src` and if found, sets it in `dst`.
function(copy_target_property src dst property)
    get_target_property(var ${src} ${property})
    if(NOT ("${var}" STREQUAL "var-NOTFOUND"))
        set_target_properties(${dst} PROPERTIES ${property} "${var}")
    endif()
endfunction()

#    copy_target_property_fallback
#
# Given a source target `src` and a destination target `dst`, reads the value of the
# `property` in `src` and only if not found, sets it in `dst` using the corresponding
# value from the configuration `config` provided.
function(copy_target_property_fallback src dst property config)
    get_target_property(var ${src} ${property})
    if("${var}" STREQUAL "var-NOTFOUND")
        get_target_property(var ${src} "${property}_${config}")
        set_target_properties(${dst} PROPERTIES ${property} "${var}")
    endif()
endfunction()

#    copy_target_property
#
# Given a source target `src` and a destination target `dst`, finds all the
# IMPORTED_CONFIGURATIONS in `src` and for each of them tries to extract the value of the
# corresponding `property` in `src`, and if found, sets it in `dst`.
function(copy_target_property_configs src dst property)
    get_target_property(configs ${src} IMPORTED_CONFIGURATIONS)
    list(LENGTH configs nconfigs)

    if("${configs}" STREQUAL "configs-NOTFOUND")
        return()
    endif()

    if(${nconfigs})
        foreach(config ${configs})
            copy_target_property(${src} ${dst} "${property}_${config}")
        endforeach()
        # TODO: fix when multiple configurations are present
        # (we're assuming the first one is appropriate)
        list(GET configs 0 config)
        copy_target_property_fallback(${src} ${dst} ${property} ${config})
    endif()
endfunction()

#     set_python_module_path
#
# Tries finding the path where the python lammps module is installed and sets the
# variable `PYLAMMPS_PATH` to that path, otherwise `PYLAMMPS_PATH` will be empty.
function(set_python_module_path)
    find_package(Python QUIET COMPONENTS Interpreter)
    if(NOT (Python_FOUND AND Python_Interpreter_FOUND))
        message(FATAL_ERROR
            "Could not find Python interpreter, make sure it is installed and enabled"
        )
    endif()
    set(find_lammps_script "
from __future__ import print_function;
import os
try:
    import lammps
    print(os.path.dirname(lammps.__file__), end='')
except:
    print('', end='')"
    )
    set(find_liblammps_script "
from __future__ import print_function;
import os
try:
    import lammps
    lmp = lammps.lammps(cmdargs='-log none -screen none'.split())
    print(lmp.lib._name, end='')
except:
    print('', end='')"
    )
    execute_process(
        COMMAND ${Python_EXECUTABLE} -c "${find_lammps_script}"
        OUTPUT_VARIABLE PYLAMMPS_PATH
    )
    execute_process(
        COMMAND ${Python_EXECUTABLE} -c "${find_liblammps_script}"
        OUTPUT_VARIABLE PYLAMMPS_LIBRARY
    )
    if("${PYLAMMPS_PATH}" STREQUAL "" OR "${PYLAMMPS_LIBRARY}" STREQUAL "")
        unset(PYLAMMPS_LIBRARY)
        find_library(PYLAMMPS_LIBRARY
            NAMES lammps
            HINTS ${Python_SITELIB}
            PATH_SUFFIXES lammps
        )
        if("${PYLAMMPS_LIBRARY}" STREQUAL "PYLAMMPS_LIBRARY-NOTFOUND")
            message(FATAL_ERROR "Unable to locate LAMMPS python module")
        else()
            get_filename_component(PYLAMMPS_PATH "${PYLAMMPS_LIBRARY}" DIRECTORY)
        endif()
    endif()
    set(PYLAMMPS_PATH "${PYLAMMPS_PATH}" PARENT_SCOPE)
    set(PYLAMMPS_LIBRARY "${PYLAMMPS_LIBRARY}" PARENT_SCOPE)
endfunction()


# Setup LAMMPS
# ------------

# We use find_lammps() first instead of find_package(LAMMPS) to avoid finding
# MPI which requires CXX enabled, but we want to enable CXX after looking for
# the NVCC compiler wrapper that comes with LAMMPS.
find_lammps()

message(STATUS "Found LAMMPS at ${LAMMPS_ROOT} (version ${LAMMPS_VERSION})")

fetch_lammps(${LAMMPS_tag})

if(NOT CMAKE_BUILD_TYPE)
    if(${LAMMPS_VERSION} GREATER 20190618)
        set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Type of build" FORCE)
    else()
        set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
    endif()
endif()

set(CMAKE_CXX_EXTENSIONS OFF CACHE FILEPATH "Use compiler extensions")
if(LAMMPS_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER ${LAMMPS_CXX_COMPILER} CACHE FILEPATH "C++ compiler")
endif()

enable_language(CXX)

find_package(LAMMPS REQUIRED)

if(TARGET LAMMPS::mpi_stubs)  # LAMMPS was built without MPI support
    target_include_directories(LAMMPS::mpi_stubs SYSTEM INTERFACE
        "${lammps_SOURCE_DIR}/src/STUBS"
    )
else()
    find_package(MPI REQUIRED)
endif()

if(NOT LAMMPS_INSTALL_PREFIX)
    get_filename_component(LAMMPS_INSTALL_PREFIX "${LAMMPS_ROOT}/../../.." ABSOLUTE)
endif()

add_library(LAMMPS_src INTERFACE)
add_library(LAMMPS::src ALIAS LAMMPS_src)

target_include_directories(LAMMPS_src INTERFACE "${lammps_SOURCE_DIR}/src")

find_package(Kokkos QUIET)
if(Kokkos_FOUND)
    message(STATUS "Kokkos support has been enabled (version ${Kokkos_VERSION})")
    get_target_property(KokkosCore       Kokkos::kokkoscore INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(KokkosContainers Kokkos::kokkoscontainers INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(KokkosAlgorithms Kokkos::kokkosalgorithms INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(KokkosSIMD       Kokkos::kokkossimd INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(KokkosCompileOptions  Kokkos::kokkoscore INTERFACE_COMPILE_OPTIONS)

    add_library(Kokkos_src INTERFACE)
    add_library(Kokkos::src ALIAS Kokkos_src)
    target_include_directories(Kokkos_src INTERFACE "${lammps_SOURCE_DIR}/src/KOKKOS")
    target_include_directories(Kokkos_src INTERFACE "${KokkosCore}")
    target_include_directories(Kokkos_src INTERFACE "${KokkosContainers}")
    target_include_directories(Kokkos_src INTERFACE "${KokkosAlgorithms}")
    target_include_directories(Kokkos_src INTERFACE "${KokkosSIMD}")

    set_target_properties(Kokkos_src PROPERTIES INTERFACE_COMPILE_OPTIONS "${KokkosCompileOptions}")
else()
    message(STATUS
        "Kokkos support is not enabled. If you built LAMMPS with Kokkos,"
        "make sure that CMake if able to find it by providing the Kokkos_ROOT."
    )
endif()
