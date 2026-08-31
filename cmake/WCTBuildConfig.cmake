# Wire-Cell Toolkit - generate WireCellUtil/BuildConfig.h (epic wct-ike, task wct-ike.3).
#
# Reproduces the header that the waf build writes via cfg.write_config_header()
# into build/WireCellUtil/BuildConfig.h, so the C++ code's #ifdef guards behave
# identically under the CMake build.  Include this after WCTDependencies
# (wct-ike.2), which sets the WCT_HAVE_<TOKEN> / WCT_HAVE_BACKTRACE_LIB /
# WCT_HAVE_BOOST_CORE_SPAN_HPP / WCT_HAVE_DLFCN_H flags it consumes.
#
# The generated header lands in ${WCT_GENERATED_INCLUDE_DIR}/WireCellUtil/ and
# is installed to <prefix>/include/WireCellUtil/, matching the waf layout.

include_guard(GLOBAL)

# --- version + development/release marker (mirrors wscript configure) ---
if(WCT_IS_DEVELOPMENT)
  set(WCT_WIRECELL_MODE_DEFINE "#define WIRECELL_DEVELOPMENT \"${WCT_VERSION}\"")
else()
  set(WCT_WIRECELL_MODE_DEFINE "#define WIRECELL_RELEASE \"${WCT_VERSION}\"")
endif()

# --- spdlog compiled-in level: SPDLOG_ACTIVE_LEVEL = SPDLOG_LEVEL_<UPPER> ---
string(TOUPPER "${WCT_SPDLOG_ACTIVE_LEVEL}" _wct_lvl)
set(WCT_SPDLOG_ACTIVE_LEVEL_MACRO "SPDLOG_LEVEL_${_wct_lvl}")

# --- per-dependency HAVE_* families ---
# For each external token, which of HAVE_<TOK> / HAVE_<TOK>_LIB / HAVE_<TOK>_INC
# the waf build emits when the dependency is found.  This matches the asymmetry
# in waft/wcb.py package_descriptions (e.g. EIGEN is header-only so only _INC;
# BZIP2/GLPK/JSONNET have no pcname so no bare HAVE_<TOK>; FFTWTHREADS is
# libs-only so HAVE_ + _LIB but no _INC).
set(_wct_have_map
    "SPDLOG:HAVE,LIB,INC"
    "ZLIB:HAVE,LIB,INC"
    "BZIP2:LIB,INC"
    "FFTW:HAVE,LIB,INC"
    "FFTWTHREADS:HAVE,LIB"
    "JSONCPP:HAVE,LIB,INC"
    "EIGEN:HAVE,INC"
    "GLPK:LIB,INC"
    "JSONNET:LIB,INC"
    "TBB:HAVE,LIB,INC"
    "HDF5:HAVE,LIB,INC"
    "ZMQ:HAVE,LIB,INC"
    "CZMQ:HAVE,LIB,INC"
    "ZYRE:HAVE,LIB,INC"
    "ZIO:HAVE,LIB,INC"
    "GRPC:HAVE,LIB,INC"
    "PROTOBUF:HAVE,LIB,INC"
    "TRITON:LIB,INC"
    "PYTHON:HAVE,LIB,INC"
    "ROOTSYS:HAVE"
    "LIBTORCH:HAVE"
    "CUDA:HAVE")

set(_lines "")
foreach(_entry IN LISTS _wct_have_map)
  string(REPLACE ":" ";" _parts "${_entry}")
  list(GET _parts 0 _tok)
  list(GET _parts 1 _flags)
  if(NOT WCT_HAVE_${_tok})
    continue()
  endif()
  string(REPLACE "," ";" _flaglist "${_flags}")
  foreach(_f IN LISTS _flaglist)
    if(_f STREQUAL "HAVE")
      string(APPEND _lines "#define HAVE_${_tok} 1\n")
    elseif(_f STREQUAL "LIB")
      string(APPEND _lines "#define HAVE_${_tok}_LIB 1\n")
    elseif(_f STREQUAL "INC")
      string(APPEND _lines "#define HAVE_${_tok}_INC 1\n")
    endif()
  endforeach()
endforeach()

# --- special checks (mirrors the explicit cfg.check* calls) ---
if(WCT_HAVE_BACKTRACE_LIB)
  string(APPEND _lines "#define HAVE_BACKTRACE_LIB 1\n")
endif()
if(WCT_HAVE_BOOST_CORE_SPAN_HPP)
  string(APPEND _lines "#define HAVE_BOOST_CORE_SPAN_HPP 1\n")
endif()
if(WCT_HAVE_ROOTSYS)
  string(APPEND _lines "#define HAVE_RTYPES_H 1\n")   # rootsys.py defines this
endif()
if(WCT_HAVE_DLFCN_H)
  string(APPEND _lines "#define HAVE_DLFCN_H 1\n")     # from the DYNAMO/dlfcn.h check
endif()

# Trim trailing newline so the template's blank line is tidy.
string(REGEX REPLACE "\n$" "" _lines "${_lines}")
set(WCT_BUILDCONFIG_HAVE_LINES "${_lines}")

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/BuildConfig.h.in"
  "${WCT_GENERATED_INCLUDE_DIR}/WireCellUtil/BuildConfig.h"
  @ONLY)

install(FILES "${WCT_GENERATED_INCLUDE_DIR}/WireCellUtil/BuildConfig.h"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/WireCellUtil")

message(STATUS "Generated WireCellUtil/BuildConfig.h (version ${WCT_VERSION})")
