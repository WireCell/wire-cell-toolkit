# Wire-Cell Toolkit - modern CMake build (parallel to the waf/./wcb build).
#
# This is an additive, parallel build system.  It does not replace and must not
# require changes to the waf build (wscript, */wscript_build, waft/*).  Build it
# out-of-tree, e.g.:
#
#     cmake -S . -B build-cmake [-DWITH_<DEP>=...] [-DWCT_WITH_TESTS=ON]
#     cmake --build build-cmake -j
#     ctest --test-dir build-cmake
#     cmake --install build-cmake --prefix /path/to/install
#
# This file (epic wct-ike, task wct-ike.1) provides the project bootstrap:
# version/build-mode detection, the C++ standard, global compile flags, and the
# user-facing options.  Dependency discovery (wct-ike.2), BuildConfig.h
# (wct-ike.3), the wct_package() helper (wct-ike.4), the per-package
# subdirectories (wct-ike.5) and the rest are layered on top by later tasks and
# are intentionally not yet wired in here.

cmake_minimum_required(VERSION 3.16)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# ---------------------------------------------------------------------------
# Version and build mode (must run before project()).
# ---------------------------------------------------------------------------
# WCT_BUILD_MODE forces development/release; empty = auto-detect (see module).
set(WCT_BUILD_MODE "" CACHE STRING
    "Force build mode: 'development', 'release', or '' to auto-detect from git")
set_property(CACHE WCT_BUILD_MODE PROPERTY STRINGS "" "development" "release")

include(WCTVersion)

project(WireCellToolkit
        VERSION ${WCT_VERSION_TRIPLET}
        LANGUAGES CXX)

# ---------------------------------------------------------------------------
# User-facing options (scalar build settings; per-dependency WITH_<DEP>
# options are defined by the dependency module in wct-ike.2).
# ---------------------------------------------------------------------------
set(WCT_CXX_STANDARD "17" CACHE STRING "C++ standard to compile against (mirrors --cxxstd)")
set_property(CACHE WCT_CXX_STANDARD PROPERTY STRINGS "17" "20" "23")

option(WCT_SPDLOG_STATIC
       "spdlog is compiled as a library (defines SPDLOG_COMPILED_LIB); mirrors --with-spdlog-static" ON)

# Minimum compiled-in log level for SPDLOG_<LEVEL>() macros; baked into
# BuildConfig.h by wct-ike.3 as SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_<UPPER>.
set(WCT_SPDLOG_ACTIVE_LEVEL "trace" CACHE STRING
    "Compiled minimum spdlog level (mirrors --with-spdlog-active-level)")
set_property(CACHE WCT_SPDLOG_ACTIVE_LEVEL PROPERTY STRINGS
             trace debug info warn error critical off)

# Tests are off by default (mirrors waf's --tests opt-in); wired up in wct-ike.9.
option(WCT_WITH_TESTS "Build and register the test suite (mirrors --tests)" OFF)

# Comma-separated experiment names whose pgrapher/experiment configs to install
# ("all" for every experiment); mirrors --install-config.  Empty installs only
# the common/ui jsonnet.
set(WCT_INSTALL_CONFIG "" CACHE STRING
    "Experiment configs to install (comma list or 'all'); mirrors --install-config")

# ---------------------------------------------------------------------------
# C++ standard.
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD ${WCT_CXX_STANDARD})
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)          # -std=c++NN, not -std=gnu++NN (matches waf)

# Shared libraries need PIC; spdlog static use likewise (see wcb.py note).
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Helpful for editors/linters; replaces waft/clang_compilation_database.py.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Default to an optimized-with-symbols build, matching waf's --build-debug
# default of "-O2 -ggdb3", while still allowing -DCMAKE_BUILD_TYPE=... override.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
      "Build type (default RelWithDebInfo ~ waf -O2 -ggdb3)" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
               Debug Release RelWithDebInfo MinSizeRel)
endif()

# ---------------------------------------------------------------------------
# Install layout.
# ---------------------------------------------------------------------------
# WCT installs libraries under lib/, never lib64/ (issue #337), unless the user
# explicitly asked for something else on the command line.  Capture whether the
# user pre-set CMAKE_INSTALL_LIBDIR before GNUInstallDirs computes its default.
if(DEFINED CMAKE_INSTALL_LIBDIR)
  set(_wct_user_libdir TRUE)
else()
  set(_wct_user_libdir FALSE)
endif()
include(GNUInstallDirs)
if(NOT _wct_user_libdir AND CMAKE_INSTALL_LIBDIR STREQUAL "lib64")
  set(CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "Object code libraries (lib)" FORCE)
endif()

# ---------------------------------------------------------------------------
# Global compile flags (mirrors wscript build()/configure() and wcb.py).
# ---------------------------------------------------------------------------
# Always-on warnings.
add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wall> $<$<COMPILE_LANGUAGE:CXX>:-Wpedantic>)

# Eigen expectations carried by the waf build.
add_compile_definitions(EIGEN_HAS_CXX11 EIGEN_FFTW_DEFAULT=1)

if(WCT_SPDLOG_STATIC)
  add_compile_definitions(SPDLOG_COMPILED_LIB=1)
endif()

# Release builds are strict (mirrors the strict_flags in waft/wcb.py); a
# development checkout instead quells warnings selectively via #pragma.
if(NOT WCT_IS_DEVELOPMENT)
  add_compile_options(
    $<$<COMPILE_LANGUAGE:CXX>:-Werror>
    $<$<COMPILE_LANGUAGE:CXX>:-Werror=return-type>
    $<$<COMPILE_LANGUAGE:CXX>:-pedantic>
    $<$<COMPILE_LANGUAGE:CXX>:-Wno-unused-local-typedefs>)
endif()

# Directory holding generated headers (BuildConfig.h lands in
# ${WCT_GENERATED_INCLUDE_DIR}/WireCellUtil/ in wct-ike.3); exposed here so
# later tasks can add it to the include search path.
set(WCT_GENERATED_INCLUDE_DIR "${CMAKE_BINARY_DIR}/generated-include"
    CACHE INTERNAL "Root of build-time generated headers")

# ---------------------------------------------------------------------------
# External dependency discovery (wct-ike.2): creates WCT::<TOKEN> interface
# targets and WCT_HAVE_<TOKEN> flags consumed by later tasks.
# ---------------------------------------------------------------------------
include(WCTDependencies)

# ---------------------------------------------------------------------------
# Generated WireCellUtil/BuildConfig.h (wct-ike.3).
# ---------------------------------------------------------------------------
include(WCTBuildConfig)

# ---------------------------------------------------------------------------
# Per-package builds via wct_package() (wct-ike.4); the package list and
# optional-dependency gating are completed in wct-ike.5 / wct-ike.6.
# ---------------------------------------------------------------------------
include(WCTPackage)

# Package subdirectories.  CMake resolves inter-target links at generate time,
# so listing order is not significant; it is kept dependency-ordered for
# readability.  (wct-ike.5 fills in the remaining packages; wct-ike.6 adds the
# optional-dependency gating that conditions the optional ones.)
# Core, always-on packages (their externals are all mandatory).
set(WCT_PACKAGES util iface aux gen sigproc sig img pgraph sio quickhull)
# pyutil is Python-gated but is a build dependency of clus when Python is
# present, so it joins the core set here (the general optional-dependency gating
# is formalized in wct-ike.6).
if(WCT_HAVE_PYTHON)
  list(APPEND WCT_PACKAGES pyutil)
endif()
list(APPEND WCT_PACKAGES clus match apps test cfg docs)

foreach(_pkg IN LISTS WCT_PACKAGES)
  add_subdirectory(${_pkg})
endforeach()

# ---------------------------------------------------------------------------
# Optional, dependency-gated packages (root, tbb, hio, pytorch, triton, zio,
# cuda) are added conditionally in wct-ike.6.  Test integration (wct-ike.9) and
# install/export + pkg-config + CMake config (wct-ike.11) also still to land.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Configuration summary.
# ---------------------------------------------------------------------------
if(WCT_IS_DEVELOPMENT)
  set(_wct_mode "development (relaxed warnings)")
else()
  set(_wct_mode "release (strict warnings)")
endif()
message(STATUS "Wire-Cell Toolkit ${WCT_VERSION} [triplet ${WCT_VERSION_TRIPLET}]")
message(STATUS "  build mode      : ${_wct_mode}")
message(STATUS "  C++ standard    : c++${CMAKE_CXX_STANDARD}")
message(STATUS "  build type      : ${CMAKE_BUILD_TYPE}")
message(STATUS "  install prefix  : ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  libdir          : ${CMAKE_INSTALL_LIBDIR}")
message(STATUS "  spdlog level    : ${WCT_SPDLOG_ACTIVE_LEVEL} (static=${WCT_SPDLOG_STATIC})")
message(STATUS "  tests           : ${WCT_WITH_TESTS}")
