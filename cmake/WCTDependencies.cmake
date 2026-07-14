# Wire-Cell Toolkit - external dependency discovery (epic wct-ike, task wct-ike.2).
#
# Ports the dependency finding that waft/generic.py + waft/wcb.py
# package_descriptions do for the waf build, and follows the same find recipe
# (find_package vs pkg-config vs find_library) already encoded in
# waft/cmake.py:DEPMAP.
#
# For every external "use" token referenced by the package wscript_build files
# (SPDLOG BOOST FFTW EIGEN DYNAMO JSONCPP JSONNET ZLIB BZIP2 FFTWTHREADS TBB
# HDF5 ZMQ CZMQ ZYRE ZIO GRPC PROTOBUF TRITON PYTHON GLPK ROOTSYS LIBTORCH CUDA
# PTHREAD) this module, on success, creates an INTERFACE IMPORTED target
#
#     WCT::<TOKEN>
#
# carrying that dependency's include dirs and link libraries, and sets a cache
# boolean WCT_HAVE_<TOKEN> (and a mirror HAVE_<TOKEN> for BuildConfig.h, wct-ike.3).
# The wct_package() helper (wct-ike.4) maps a package's use= tokens to these
# targets; the optional-dependency gating (wct-ike.6) reads WCT_HAVE_<TOKEN>.
#
# Each token honors a cache option WITH_<TOKEN>:
#   ""              auto: mandatory -> required; optional -> quiet probe
#   no|off|false    skip (error if the token is mandatory)
#   yes|on|true     require
#   <path>          use <path> as an install-prefix hint, then require
# This mirrors the --with-NAME[/-include/-lib/-libs] semantics of generic.py.

include_guard(GLOBAL)

# Modern Boost ships BoostConfig.cmake; prefer it over the deprecated FindBoost.
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

find_package(PkgConfig QUIET)

# Mandatory tokens needed by the core library chain (mirrors generic.py
# mandatory=True); all others are optional and gate the packages that use them.
set(WCT_DEP_MANDATORY
    SPDLOG BOOST FFTW EIGEN DYNAMO JSONCPP JSONNET ZLIB BZIP2 PTHREAD)
set(WCT_DEP_OPTIONAL
    FFTWTHREADS GLPK TBB HDF5 ZMQ CZMQ ZYRE ZIO GRPC PROTOBUF TRITON
    PYTHON ROOTSYS LIBTORCH CUDA BACKTRACE)

# Opt-in dependencies: like waf's with_p() gating (waft/wcb.py), these are only
# probed when their WITH_<TOKEN> is set; an empty value means "do not use",
# never auto-detect.  (ROOT, libtorch and CUDA pull in heavy toolchains and are
# enabled deliberately, e.g. configit.sh passes --with-root.)
set(WCT_DEP_OPTIN ROOTSYS LIBTORCH CUDA KOKKOS)

# Declare the per-dependency options.  These mirror waf's family of
# --with-NAME[/-include/-lib/-libs] options:
#   WITH_<TOKEN>          '' (auto) | no | yes | <install-prefix path>
#   WITH_<TOKEN>_INCLUDE  include dir(s) (comma list), overrides <prefix>/include
#   WITH_<TOKEN>_LIB      library dir(s) (comma list), overrides <prefix>/lib
#   WITH_<TOKEN>_LIBS     exact library name(s) (comma list), overrides defaults
#                         e.g. -DWITH_JSONNET_LIBS=gojsonnet
# Run `cmake -LH <build>` (or ccmake) to see them all with their help strings.
foreach(_tok IN LISTS WCT_DEP_MANDATORY WCT_DEP_OPTIONAL)
  set(WITH_${_tok} "" CACHE STRING
      "Discovery for ${_tok}: ''=auto, no/yes, or an install-prefix path")
  set(WITH_${_tok}_INCLUDE "" CACHE STRING
      "Include dir(s) for ${_tok} (comma list; overrides <prefix>/include)")
  set(WITH_${_tok}_LIB "" CACHE STRING
      "Library dir(s) for ${_tok} (comma list; overrides <prefix>/lib)")
  set(WITH_${_tok}_LIBS "" CACHE STRING
      "Exact library name(s) for ${_tok} (comma list; overrides defaults)")
endforeach()

set(WCT_EXTERNAL_FOUND "" CACHE INTERNAL "External tokens discovered")
set(WCT_EXTERNAL_MISSING "" CACHE INTERNAL "Optional external tokens not found")

# Resolve WITH_<tok> into _wct_mode (SKIP|AUTO|REQUIRE) and apply any path hint
# to the search paths used by find_package / find_library / pkg-config.
macro(_wct_intent tok)
  list(FIND WCT_DEP_MANDATORY "${tok}" _wct_mand_idx)
  if(_wct_mand_idx GREATER -1)
    set(_wct_mandatory TRUE)
  else()
    set(_wct_mandatory FALSE)
  endif()

  set(_wct_val "${WITH_${tok}}")
  string(TOLOWER "${_wct_val}" _wct_vl)
  list(FIND WCT_DEP_OPTIN "${tok}" _wct_optin_idx)
  set(_wct_prefix "")
  if(_wct_vl STREQUAL "")
    if(_wct_mandatory)
      set(_wct_mode REQUIRE)
    elseif(_wct_optin_idx GREATER -1)
      set(_wct_mode SKIP)        # opt-in: not requested => not used
    else()
      set(_wct_mode AUTO)
    endif()
  elseif(_wct_vl MATCHES "^(no|off|false)$")
    if(_wct_mandatory)
      message(FATAL_ERROR "${tok} is mandatory; WITH_${tok}=${_wct_val} is not allowed")
    endif()
    set(_wct_mode SKIP)
  elseif(_wct_vl MATCHES "^(yes|on|true)$")
    set(_wct_mode REQUIRE)
  else()
    # A path: use it as an install-prefix hint for all search mechanisms.
    set(_wct_prefix "${_wct_val}")
    list(PREPEND CMAKE_PREFIX_PATH "${_wct_val}")
    if(EXISTS "${_wct_val}/lib/pkgconfig")
      set(ENV{PKG_CONFIG_PATH} "${_wct_val}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    endif()
    set(_wct_mode REQUIRE)
  endif()

  # Fine-grained overrides (mirror --with-NAME-include/-lib/-libs).  These set,
  # for the discovery block below:
  #   _wct_inc_hints  include dir(s) to search (find_path HINTS / CMAKE_INCLUDE_PATH)
  #   _wct_lib_hints  library dir(s) to search (find_library HINTS / CMAKE_LIBRARY_PATH)
  #   _wct_lib_names  explicit library name(s), or empty to use the built-in default
  set(_wct_inc_hints "")
  set(_wct_lib_hints "")
  set(_wct_lib_names "")
  if(NOT "${WITH_${tok}_INCLUDE}" STREQUAL "")
    string(REPLACE "," ";" _wct_inc_hints "${WITH_${tok}_INCLUDE}")
  elseif(_wct_prefix)
    set(_wct_inc_hints "${_wct_prefix}/include")
  endif()
  if(NOT "${WITH_${tok}_LIB}" STREQUAL "")
    string(REPLACE "," ";" _wct_lib_hints "${WITH_${tok}_LIB}")
  elseif(_wct_prefix)
    set(_wct_lib_hints "${_wct_prefix}/lib")
  endif()
  if(NOT "${WITH_${tok}_LIBS}" STREQUAL "")
    string(REPLACE "," ";" _wct_lib_names "${WITH_${tok}_LIBS}")
  endif()
  # Make find_package()/pkg-config honor explicit include/lib dirs too.
  if(_wct_inc_hints)
    list(APPEND CMAKE_INCLUDE_PATH ${_wct_inc_hints})
  endif()
  foreach(_lh IN LISTS _wct_lib_hints)
    list(APPEND CMAKE_LIBRARY_PATH "${_lh}")
    if(EXISTS "${_lh}/pkgconfig")
      set(ENV{PKG_CONFIG_PATH} "${_lh}/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    endif()
  endforeach()
endmacro()

# Create the WCT::<tok> interface target and record the token as found.
function(_wct_provide tok)
  cmake_parse_arguments(A "" "" "LINK;INCLUDE" ${ARGN})
  if(NOT TARGET WCT::${tok})
    add_library(WCT::${tok} INTERFACE IMPORTED GLOBAL)
  endif()
  if(A_LINK)
    set_property(TARGET WCT::${tok} PROPERTY INTERFACE_LINK_LIBRARIES "${A_LINK}")
  endif()
  if(A_INCLUDE)
    list(REMOVE_ITEM A_INCLUDE "")
    set_property(TARGET WCT::${tok} PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${A_INCLUDE}")
  endif()
  set(WCT_HAVE_${tok} TRUE CACHE INTERNAL "WCT has ${tok}")
  set(HAVE_${tok} 1 CACHE INTERNAL "BuildConfig mirror for ${tok}")
  list(APPEND WCT_EXTERNAL_FOUND ${tok})
  list(REMOVE_DUPLICATES WCT_EXTERNAL_FOUND)
  set(WCT_EXTERNAL_FOUND "${WCT_EXTERNAL_FOUND}" CACHE INTERNAL "")
endfunction()

function(_wct_missing tok)
  set(WCT_HAVE_${tok} FALSE CACHE INTERNAL "WCT has ${tok}")
  list(APPEND WCT_EXTERNAL_MISSING ${tok})
  list(REMOVE_DUPLICATES WCT_EXTERNAL_MISSING)
  set(WCT_EXTERNAL_MISSING "${WCT_EXTERNAL_MISSING}" CACHE INTERNAL "")
endfunction()

# ===========================================================================
# Mandatory dependencies
# ===========================================================================

# --- spdlog (config; pkg-config fallback) ---
_wct_intent(SPDLOG)
find_package(spdlog QUIET)
if(spdlog_FOUND)
  _wct_provide(SPDLOG LINK spdlog::spdlog)
elseif(PkgConfig_FOUND)
  pkg_check_modules(WCT_SPDLOG REQUIRED IMPORTED_TARGET spdlog)
  _wct_provide(SPDLOG LINK PkgConfig::WCT_SPDLOG)
else()
  message(FATAL_ERROR "spdlog not found (set WITH_SPDLOG=<prefix>)")
endif()

# --- Boost ---
_wct_intent(BOOST)
set(_wct_boost_comps filesystem graph thread program_options iostreams regex)
find_package(Boost CONFIG QUIET COMPONENTS ${_wct_boost_comps})
if(NOT Boost_FOUND)
  find_package(Boost REQUIRED COMPONENTS ${_wct_boost_comps})
endif()
set(_wct_boost_libs "")
foreach(_c IN LISTS _wct_boost_comps)
  list(APPEND _wct_boost_libs Boost::${_c})
endforeach()
if(TARGET Boost::headers)
  list(APPEND _wct_boost_libs Boost::headers)
endif()
_wct_provide(BOOST LINK "${_wct_boost_libs}")

# --- FFTW (single-precision, via pkg-config) ---
_wct_intent(FFTW)
pkg_check_modules(WCT_FFTW REQUIRED IMPORTED_TARGET fftw3f)
_wct_provide(FFTW LINK PkgConfig::WCT_FFTW)

# --- Eigen3 (header only; config, pkg-config fallback) ---
_wct_intent(EIGEN)
find_package(Eigen3 QUIET)
if(Eigen3_FOUND)
  _wct_provide(EIGEN LINK Eigen3::Eigen)
elseif(PkgConfig_FOUND)
  pkg_check_modules(WCT_EIGEN REQUIRED IMPORTED_TARGET eigen3)
  _wct_provide(EIGEN LINK PkgConfig::WCT_EIGEN)
else()
  message(FATAL_ERROR "Eigen3 not found (set WITH_EIGEN=<prefix>)")
endif()

# --- dl / dlfcn.h (the DYNAMO token) ---
_wct_intent(DYNAMO)
find_path(WCT_DLFCN_INCLUDE_DIR NAMES dlfcn.h)
if(NOT WCT_DLFCN_INCLUDE_DIR)
  message(FATAL_ERROR "dlfcn.h not found")
endif()
set(WCT_HAVE_DLFCN_H 1 CACHE INTERNAL "have dlfcn.h")
_wct_provide(DYNAMO LINK ${CMAKE_DL_LIBS} INCLUDE ${WCT_DLFCN_INCLUDE_DIR})

# --- JsonCpp (pkg-config) ---
_wct_intent(JSONCPP)
pkg_check_modules(WCT_JSONCPP REQUIRED IMPORTED_TARGET jsoncpp)
_wct_provide(JSONCPP LINK PkgConfig::WCT_JSONCPP)

# --- Jsonnet (find_library; jsonnet or gojsonnet) ---
# The C and Go implementations are ABI-compatible; select one with
# -DWITH_JSONNET_LIBS=gojsonnet (mirrors wcb's --with-jsonnet-libs=gojsonnet).
_wct_intent(JSONNET)
if(_wct_lib_names)
  set(_wct_jsonnet_names ${_wct_lib_names})
else()
  set(_wct_jsonnet_names jsonnet gojsonnet)
endif()
find_library(WCT_JSONNET_LIBRARY NAMES ${_wct_jsonnet_names} HINTS ${_wct_lib_hints})
find_path(WCT_JSONNET_INCLUDE_DIR NAMES libjsonnet.h HINTS ${_wct_inc_hints})
if(WCT_JSONNET_LIBRARY AND WCT_JSONNET_INCLUDE_DIR)
  _wct_provide(JSONNET LINK ${WCT_JSONNET_LIBRARY} INCLUDE ${WCT_JSONNET_INCLUDE_DIR})
else()
  message(FATAL_ERROR "jsonnet not found (set WITH_JSONNET=<prefix>)")
endif()

# --- zlib ---
_wct_intent(ZLIB)
find_package(ZLIB REQUIRED)
_wct_provide(ZLIB LINK ZLIB::ZLIB)

# --- bzip2 ---
_wct_intent(BZIP2)
find_package(BZip2 REQUIRED)
_wct_provide(BZIP2 LINK BZip2::BZip2)

# --- pthread ---
_wct_intent(PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
_wct_provide(PTHREAD LINK Threads::Threads)

# ===========================================================================
# Optional dependencies
# ===========================================================================

# Small helpers for the common optional patterns.
function(_wct_opt_pkgconfig tok pcname)
  _wct_intent(${tok})
  if(_wct_mode STREQUAL "SKIP" OR NOT PkgConfig_FOUND)
    _wct_missing(${tok})
    return()
  endif()
  set(_req "")
  if(_wct_mode STREQUAL "REQUIRE")
    set(_req REQUIRED)
  endif()
  pkg_check_modules(WCT_${tok} ${_req} QUIET IMPORTED_TARGET ${pcname})
  if(WCT_${tok}_FOUND)
    _wct_provide(${tok} LINK PkgConfig::WCT_${tok})
  else()
    _wct_missing(${tok})
  endif()
endfunction()

function(_wct_opt_findlib tok header)
  # remaining ARGN are the default candidate library names
  _wct_intent(${tok})
  if(_wct_mode STREQUAL "SKIP")
    _wct_missing(${tok})
    return()
  endif()
  # WITH_<tok>_LIBS overrides the library name(s); WITH_<tok>_LIB/_INCLUDE add
  # search hints (both computed by _wct_intent).
  if(_wct_lib_names)
    set(_names ${_wct_lib_names})
  else()
    set(_names ${ARGN})
  endif()
  set(_libs "")
  set(_ok TRUE)
  foreach(_nm IN LISTS _names)
    find_library(WCT_${tok}_${_nm}_LIB NAMES ${_nm} HINTS ${_wct_lib_hints})
    if(WCT_${tok}_${_nm}_LIB)
      list(APPEND _libs ${WCT_${tok}_${_nm}_LIB})
    else()
      set(_ok FALSE)
    endif()
  endforeach()
  set(_inc "")
  if(NOT header STREQUAL "")
    find_path(WCT_${tok}_INCLUDE_DIR NAMES ${header} HINTS ${_wct_inc_hints})
    if(WCT_${tok}_INCLUDE_DIR)
      set(_inc ${WCT_${tok}_INCLUDE_DIR})
    else()
      set(_ok FALSE)
    endif()
  endif()
  if(_ok)
    _wct_provide(${tok} LINK "${_libs}" INCLUDE "${_inc}")
  else()
    if(_wct_mode STREQUAL "REQUIRE")
      message(FATAL_ERROR "${tok} requested but not found")
    endif()
    _wct_missing(${tok})
  endif()
endfunction()

# --- FFTW threads (single lib) ---
_wct_opt_findlib(FFTWTHREADS "" fftw3f_threads)

# --- GLPK ---
_wct_opt_findlib(GLPK "glpk.h" glpk)

# --- TBB (config; pkg-config fallback) ---
_wct_intent(TBB)
if(_wct_mode STREQUAL "SKIP")
  _wct_missing(TBB)
else()
  find_package(TBB QUIET)
  if(TBB_FOUND)
    _wct_provide(TBB LINK TBB::tbb)
  elseif(PkgConfig_FOUND)
    pkg_check_modules(WCT_TBB QUIET IMPORTED_TARGET tbb)
    if(WCT_TBB_FOUND)
      _wct_provide(TBB LINK PkgConfig::WCT_TBB)
    else()
      _wct_missing(TBB)
    endif()
  else()
    _wct_missing(TBB)
  endif()
endif()

# --- HDF5 (FindHDF5 is unreliable; use the .pc as DEPMAP does) ---
_wct_opt_pkgconfig(HDF5 hdf5)

# --- ZeroMQ stack ---
_wct_opt_pkgconfig(ZMQ  libzmq)
_wct_opt_pkgconfig(CZMQ libczmq)
_wct_opt_pkgconfig(ZYRE libzyre)
_wct_opt_pkgconfig(ZIO  libzio)

# --- Protobuf (module/config; pkg-config fallback) ---
_wct_intent(PROTOBUF)
if(_wct_mode STREQUAL "SKIP")
  _wct_missing(PROTOBUF)
else()
  find_package(Protobuf QUIET)
  if(Protobuf_FOUND)
    _wct_provide(PROTOBUF LINK protobuf::libprotobuf)
  elseif(PkgConfig_FOUND)
    pkg_check_modules(WCT_PROTOBUF QUIET IMPORTED_TARGET protobuf)
    if(WCT_PROTOBUF_FOUND)
      _wct_provide(PROTOBUF LINK PkgConfig::WCT_PROTOBUF)
    else()
      _wct_missing(PROTOBUF)
    endif()
  else()
    _wct_missing(PROTOBUF)
  endif()
endif()

# --- gRPC ---
_wct_intent(GRPC)
if(_wct_mode STREQUAL "SKIP")
  _wct_missing(GRPC)
else()
  find_package(gRPC QUIET)
  if(gRPC_FOUND)
    _wct_provide(GRPC LINK gRPC::grpc++ gRPC::grpc gRPC::gpr)
  else()
    _wct_missing(GRPC)
  endif()
endif()

# --- Triton client (find_library set + header) ---
_wct_opt_findlib(TRITON "grpc_client.h"
  grpcclient tritoncommonerror tritoncommonmodelconfig tritoncommonlogging
  tritontableprinter tritonthreadpool tritonasyncworkqueue)

# --- Python3 (embedding) ---
_wct_intent(PYTHON)
if(_wct_mode STREQUAL "SKIP")
  _wct_missing(PYTHON)
else()
  find_package(Python3 QUIET COMPONENTS Development.Embed)
  if(Python3_FOUND)
    _wct_provide(PYTHON LINK Python3::Python)
  else()
    _wct_missing(PYTHON)
  endif()
endif()

# --- ROOT ---
_wct_intent(ROOTSYS)
if(_wct_mode STREQUAL "SKIP")
  _wct_missing(ROOTSYS)
else()
  find_package(ROOT QUIET)
  if(ROOT_FOUND)
    _wct_provide(ROOTSYS LINK ${ROOT_LIBRARIES} INCLUDE ${ROOT_INCLUDE_DIRS})
  else()
    _wct_missing(ROOTSYS)
  endif()
endif()

# --- libtorch ---
_wct_intent(LIBTORCH)
if(_wct_mode STREQUAL "SKIP")
  _wct_missing(LIBTORCH)
else()
  find_package(Torch QUIET)
  if(Torch_FOUND)
    _wct_provide(LIBTORCH LINK ${TORCH_LIBRARIES} INCLUDE ${TORCH_INCLUDE_DIRS})
  else()
    _wct_missing(LIBTORCH)
  endif()
endif()

# --- CUDA toolkit (enable_language(CUDA) is handled in wct-ike.8) ---
_wct_intent(CUDA)
if(_wct_mode STREQUAL "SKIP")
  _wct_missing(CUDA)
else()
  find_package(CUDAToolkit QUIET)
  if(CUDAToolkit_FOUND)
    _wct_provide(CUDA LINK CUDA::cudart)
  else()
    _wct_missing(CUDA)
  endif()
endif()

# --- Kokkos (opt-in; compiles .kokkos sources, see wct-ike.8) ---
_wct_intent(KOKKOS)
if(_wct_mode STREQUAL "SKIP")
  _wct_missing(KOKKOS)
else()
  find_package(Kokkos QUIET)
  if(Kokkos_FOUND)
    _wct_provide(KOKKOS LINK Kokkos::kokkos)
  else()
    _wct_missing(KOKKOS)
  endif()
endif()

# ===========================================================================
# Feature checks consumed by BuildConfig.h (wct-ike.3)
# ===========================================================================
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

# backtrace library (Exceptions.h): HAVE_BACKTRACE_LIB
_wct_intent(BACKTRACE)
if(NOT _wct_mode STREQUAL "SKIP")
  find_library(WCT_BACKTRACE_LIB NAMES backtrace)
  find_path(WCT_BACKTRACE_INCLUDE_DIR NAMES backtrace.h)
  if(WCT_BACKTRACE_LIB AND WCT_BACKTRACE_INCLUDE_DIR)
    cmake_push_check_state()
    set(CMAKE_REQUIRED_LIBRARIES ${WCT_BACKTRACE_LIB})
    set(CMAKE_REQUIRED_INCLUDES ${WCT_BACKTRACE_INCLUDE_DIR})
    check_cxx_source_compiles("
      #include <backtrace.h>
      int main(){ backtrace_create_state(nullptr,false,nullptr,nullptr); }
    " WCT_HAVE_BACKTRACE_LIB)
    cmake_pop_check_state()
  endif()
endif()
if(WCT_HAVE_BACKTRACE_LIB)
  _wct_provide(BACKTRACE LINK ${WCT_BACKTRACE_LIB} INCLUDE ${WCT_BACKTRACE_INCLUDE_DIR})
  set(HAVE_BACKTRACE_LIB 1 CACHE INTERNAL "")
else()
  _wct_missing(BACKTRACE)
endif()

# boost/core/span.hpp (until everyone is on boost >= 1.78): HAVE_BOOST_CORE_SPAN_HPP
cmake_push_check_state()
if(TARGET Boost::headers)
  get_target_property(_wct_boost_inc Boost::headers INTERFACE_INCLUDE_DIRECTORIES)
  if(_wct_boost_inc)
    set(CMAKE_REQUIRED_INCLUDES ${_wct_boost_inc})
  endif()
endif()
check_cxx_source_compiles("
  #include <boost/core/span.hpp>
  int main(){ return 0; }
" WCT_HAVE_BOOST_CORE_SPAN_HPP)
cmake_pop_check_state()
if(WCT_HAVE_BOOST_CORE_SPAN_HPP)
  set(HAVE_BOOST_CORE_SPAN_HPP 1 CACHE INTERNAL "")
endif()

# ===========================================================================
# Summary
# ===========================================================================
message(STATUS "WCT dependencies found  : ${WCT_EXTERNAL_FOUND}")
message(STATUS "WCT optional deps absent: ${WCT_EXTERNAL_MISSING}")
