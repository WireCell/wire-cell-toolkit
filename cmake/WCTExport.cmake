# Wire-Cell Toolkit - install/export + pkg-config + CMake package config
# (epic wct-ike, task wct-ike.11).
#
# Call once at the end of the top-level CMakeLists, after all packages are
# added (so the WireCellToolkitTargets export set and the WCT::<TOKEN> wrapper
# targets are complete).  Produces, installed under the conventional locations:
#
#   lib/cmake/WireCellToolkit/WireCellToolkitTargets.cmake        (install(EXPORT))
#   lib/cmake/WireCellToolkit/WireCellToolkitConfig.cmake         (find_package entry)
#   lib/cmake/WireCellToolkit/WireCellToolkitConfigVersion.cmake
#   lib/pkgconfig/wire-cell-toolkit.pc
#
# The native install(EXPORT) carries the WCT-internal target graph and library
# locations.  The Config recreates the external WCT::<TOKEN> wrappers the
# exported targets reference, reusing the exact link/include recipe captured on
# each wrapper at discovery time (get_target_property) -- no second copy of the
# dependency map.  This preserves the established downstream contract
# (WireCell::Util, WireCell::Aux, ..., aggregate WireCell::WireCell).

include_guard(GLOBAL)
include(CMakePackageConfigHelpers)

set(WCT_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/WireCellToolkit")

# How a consumer re-discovers each external token.  Sets, in the caller:
#   _cmd       commands (a ;-list of lines) to run before recreating the wrapper
#   _needpc    TRUE if the commands use pkg-config
#   _required  TRUE for mandatory deps (find_dependency / REQUIRED)
function(_wct_refind tok out_cmd out_needpc out_required)
  set(_c "")
  set(_pc FALSE)
  set(_req FALSE)
  if(tok STREQUAL "SPDLOG")
    set(_c "find_dependency(spdlog)")
    set(_req TRUE)
  elseif(tok STREQUAL "BOOST")
    set(_c "find_dependency(Boost COMPONENTS filesystem graph thread program_options iostreams regex)")
    set(_req TRUE)
  elseif(tok STREQUAL "EIGEN")
    set(_c "find_dependency(Eigen3)")
    set(_req TRUE)
  elseif(tok STREQUAL "ZLIB")
    set(_c "find_dependency(ZLIB)")
    set(_req TRUE)
  elseif(tok STREQUAL "BZIP2")
    set(_c "find_dependency(BZip2)")
    set(_req TRUE)
  elseif(tok STREQUAL "PTHREAD")
    set(_c "set(THREADS_PREFER_PTHREAD_FLAG ON)" "find_dependency(Threads)")
    set(_req TRUE)
  elseif(tok STREQUAL "FFTW")
    set(_c "pkg_check_modules(WCT_FFTW REQUIRED IMPORTED_TARGET fftw3f)")
    set(_pc TRUE)
    set(_req TRUE)
  elseif(tok STREQUAL "JSONCPP")
    set(_c "pkg_check_modules(WCT_JSONCPP REQUIRED IMPORTED_TARGET jsoncpp)")
    set(_pc TRUE)
    set(_req TRUE)
  elseif(tok STREQUAL "HDF5")
    set(_c "pkg_check_modules(WCT_HDF5 IMPORTED_TARGET hdf5)")
    set(_pc TRUE)
  elseif(tok STREQUAL "ZMQ")
    set(_c "pkg_check_modules(WCT_ZMQ IMPORTED_TARGET libzmq)")
    set(_pc TRUE)
  elseif(tok STREQUAL "CZMQ")
    set(_c "pkg_check_modules(WCT_CZMQ IMPORTED_TARGET libczmq)")
    set(_pc TRUE)
  elseif(tok STREQUAL "ZYRE")
    set(_c "pkg_check_modules(WCT_ZYRE IMPORTED_TARGET libzyre)")
    set(_pc TRUE)
  elseif(tok STREQUAL "ZIO")
    set(_c "pkg_check_modules(WCT_ZIO IMPORTED_TARGET libzio)")
    set(_pc TRUE)
  elseif(tok STREQUAL "TBB")
    set(_c "find_package(TBB QUIET)")
  elseif(tok STREQUAL "PROTOBUF")
    set(_c "find_package(Protobuf QUIET)")
  elseif(tok STREQUAL "GRPC")
    set(_c "find_package(gRPC QUIET)")
  elseif(tok STREQUAL "PYTHON")
    set(_c "find_package(Python3 QUIET COMPONENTS Development.Embed)")
  elseif(tok STREQUAL "ROOTSYS")
    set(_c "find_package(ROOT QUIET)")
  elseif(tok STREQUAL "LIBTORCH")
    set(_c "find_package(Torch QUIET)")
  elseif(tok STREQUAL "CUDA")
    set(_c "find_package(CUDAToolkit QUIET)")
  elseif(tok STREQUAL "KOKKOS")
    set(_c "find_package(Kokkos QUIET)")
  endif()
  # DYNAMO/JSONNET/FFTWTHREADS/GLPK/TRITON/BACKTRACE: no re-find; their wrapper
  # bakes the discovered absolute path(s)/CMAKE_DL_LIBS captured below.
  set(${out_cmd} "${_c}" PARENT_SCOPE)
  set(${out_needpc} ${_pc} PARENT_SCOPE)
  set(${out_required} ${_req} PARENT_SCOPE)
endfunction()

function(wct_install_exports)
  # ---- pkg-config .pc ----
  set(prefix "${CMAKE_INSTALL_PREFIX}")
  set(libdir "${CMAKE_INSTALL_FULL_LIBDIR}")
  set(VERSION "${WCT_VERSION}")
  set(LLIBS "-lWireCellAux -lWireCellIface -lWireCellUtil")
  # Requires = pkg-config names of the found deps that have one (mirrors the
  # REQUIRES set the waf build accumulates).
  set(_pcmap
      "SPDLOG=spdlog" "ZLIB=zlib" "FFTW=fftw3f" "JSONCPP=jsoncpp" "EIGEN=eigen3"
      "TBB=tbb" "HDF5=hdf5" "PROTOBUF=protobuf" "PYTHON=python3-embed"
      "ZMQ=libzmq" "CZMQ=libczmq" "ZYRE=libzyre" "ZIO=libzio" "GRPC=grpc++")
  set(_reqs "")
  foreach(_m IN LISTS _pcmap)
    string(REPLACE "=" ";" _kv "${_m}")
    list(GET _kv 0 _t)
    list(GET _kv 1 _pc)
    if(_t IN_LIST WCT_EXTERNAL_FOUND)
      list(APPEND _reqs "${_pc}")
    endif()
  endforeach()
  string(REPLACE ";" " " REQUIRES "${_reqs}")
  configure_file("${CMAKE_SOURCE_DIR}/wire-cell-toolkit.pc.in"
                 "${CMAKE_BINARY_DIR}/wire-cell-toolkit.pc" @ONLY)
  install(FILES "${CMAKE_BINARY_DIR}/wire-cell-toolkit.pc"
          DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

  # ---- native target export ----
  install(EXPORT WireCellToolkitTargets
          NAMESPACE WireCell::
          DESTINATION "${WCT_CMAKE_INSTALL_DIR}"
          FILE WireCellToolkitTargets.cmake)

  # ---- generate the external-wrapper recreation block for the Config ----
  set(_block "")
  set(_needpc FALSE)
  set(_seen "")
  foreach(_tok IN LISTS WCT_EXTERNAL_FOUND)
    if(NOT TARGET WCT::${_tok})
      continue()
    endif()
    get_target_property(_lnk WCT::${_tok} INTERFACE_LINK_LIBRARIES)
    get_target_property(_inc WCT::${_tok} INTERFACE_INCLUDE_DIRECTORIES)
    if(NOT _lnk)
      set(_lnk "")
    endif()
    if(NOT _inc)
      set(_inc "")
    endif()
    _wct_refind("${_tok}" _cmd _pc _req)
    if(_pc)
      set(_needpc TRUE)
    endif()
    foreach(_line IN LISTS _cmd)
      string(APPEND _block "${_line}\n")
    endforeach()
    string(APPEND _block "if(NOT TARGET WCT::${_tok})\n")
    string(APPEND _block "  add_library(WCT::${_tok} INTERFACE IMPORTED)\n")
    string(APPEND _block "  set_target_properties(WCT::${_tok} PROPERTIES\n")
    string(APPEND _block "    INTERFACE_LINK_LIBRARIES \"${_lnk}\"\n")
    string(APPEND _block "    INTERFACE_INCLUDE_DIRECTORIES \"${_inc}\")\n")
    string(APPEND _block "endif()\n")
  endforeach()
  if(_needpc)
    set(_block "find_package(PkgConfig REQUIRED)\n${_block}")
  endif()
  set(WCT_EXPORT_FIND_BLOCK "${_block}")

  configure_package_config_file(
    "${CMAKE_SOURCE_DIR}/cmake/WireCellToolkitConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/WireCellToolkitConfig.cmake"
    INSTALL_DESTINATION "${WCT_CMAKE_INSTALL_DIR}")
  write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/WireCellToolkitConfigVersion.cmake"
    VERSION "${WCT_VERSION_TRIPLET}"
    COMPATIBILITY SameMajorVersion)
  install(FILES
    "${CMAKE_BINARY_DIR}/WireCellToolkitConfig.cmake"
    "${CMAKE_BINARY_DIR}/WireCellToolkitConfigVersion.cmake"
    DESTINATION "${WCT_CMAKE_INSTALL_DIR}")
endfunction()
