# Wire-Cell Toolkit - wct_package() helper (epic wct-ike, task wct-ike.4).
#
# Reproduces the per-package conventions of waft/smplpkgs.py:smplpkg() so a
# package directory's CMakeLists.txt (wct-ike.5) is, like its wscript_build, a
# one-liner:
#
#     wct_package(WireCellUtil
#                 USE SPDLOG BOOST FFTW EIGEN DYNAMO JSONCPP JSONNET ZLIB)
#
# Conventions (matching smplpkg):
#   src/*.cxx                -> SHARED library target <NAME>  (libNAME.so)
#   inc/<NAME>/*.h           -> public headers, installed to include/<NAME>
#   apps/*.cxx               -> one executable each, installed to bin
#   USE/APP_USE/TEST_USE      tokens -> internal WCT lib targets (WireCell*/WCP*)
#                             or external WCT::<TOKEN> targets from WCTDependencies
#
# Tests (test/, doctest*, scripts) are handled by wct-ike.9, ROOT dictionaries
# (dict/) by wct-ike.7, and .cu/.kokkos sources by wct-ike.8; this helper calls
# into those when their hooks are present and otherwise leaves them alone.

include_guard(GLOBAL)

# Installed libraries find their siblings and external deps without LD_LIBRARY_PATH.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

# Resolve a list of "use" tokens to linkable targets.  Internal packages
# (WireCell*/WCP*) map to same-named targets (resolved at generate time, so
# subdirectory order does not matter); external tokens map to WCT::<UPPER>.
# Tokens for absent optional externals are dropped (the package gating in
# wct-ike.6 is what keeps a package from being built when a hard dep is absent).
function(wct_resolve_uses out_var)
  set(_res "")
  foreach(_tok IN LISTS ARGN)
    if(_tok MATCHES "^(WireCell|WCP)")
      list(APPEND _res ${_tok})
    else()
      string(TOUPPER "${_tok}" _up)
      if(TARGET WCT::${_up})
        list(APPEND _res WCT::${_up})
      elseif(DEFINED WCT_HAVE_${_up})
        message(STATUS "wct_package: optional '${_tok}' absent, not linking")
      else()
        message(WARNING "wct_package: unknown use token '${_tok}'")
      endif()
    endif()
  endforeach()
  set(${out_var} "${_res}" PARENT_SCOPE)
endfunction()

function(wct_package NAME)
  cmake_parse_arguments(P "" "HEADER_SUBDIR" "USE;APP_USE;TEST_USE" ${ARGN})

  set(_dir "${CMAKE_CURRENT_SOURCE_DIR}")
  if(P_HEADER_SUBDIR)
    set(_hsub "${P_HEADER_SUBDIR}")
  else()
    set(_hsub "${NAME}")
  endif()

  # Effective use-sets: apps/tests also see the library's own USE (as in smplpkg).
  set(_use      ${P_USE})
  set(_app_use  ${P_USE} ${P_APP_USE})
  set(_test_use ${P_USE} ${P_TEST_USE})

  wct_resolve_uses(_use_tgts ${_use})

  # ----- shared library from src/*.cxx (+ .cu/.kokkos when enabled) -----
  set(_have_lib FALSE)
  if(IS_DIRECTORY "${_dir}/src")
    file(GLOB _srcs CONFIGURE_DEPENDS "${_dir}/src/*.cxx")
    # CUDA sources compile via the CUDA language enabled at the top level
    # (wct-ike.8); only glob them when CUDA is in use.
    if(WCT_HAVE_CUDA)
      file(GLOB _cu CONFIGURE_DEPENDS "${_dir}/src/*.cu")
      list(APPEND _srcs ${_cu})
    endif()
    # Kokkos sources use a non-standard extension; compile them as C++ with the
    # Kokkos toolchain when Kokkos is in use.
    if(WCT_HAVE_KOKKOS)
      file(GLOB _kk CONFIGURE_DEPENDS "${_dir}/src/*.kokkos")
      if(_kk)
        set_source_files_properties(${_kk} PROPERTIES LANGUAGE CXX)
        list(APPEND _srcs ${_kk})
      endif()
    endif()
    if(_srcs)
      add_library(${NAME} SHARED ${_srcs})
      set(_have_lib TRUE)
      target_include_directories(${NAME} PUBLIC
        "$<BUILD_INTERFACE:${_dir}/inc>"
        "$<BUILD_INTERFACE:${WCT_GENERATED_INCLUDE_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
      if(_use_tgts)
        target_link_libraries(${NAME} PUBLIC ${_use_tgts})
      endif()
      # Exported name strips the WireCell/WCP prefix so downstream consumers use
      # WireCell::Util, WireCell::Quickhull, ... (the established contract from
      # the waf-generated config); the on-disk library stays libWireCellUtil.so.
      set(_exp "${NAME}")
      foreach(_pre WireCell WCP)
        if(_exp MATCHES "^${_pre}(.+)$")
          set(_exp "${CMAKE_MATCH_1}")
          break()
        endif()
      endforeach()
      set_target_properties(${NAME} PROPERTIES
        EXPORT_NAME "${_exp}"
        INSTALL_RPATH "$ORIGIN;${CMAKE_INSTALL_FULL_LIBDIR}")
      install(TARGETS ${NAME} EXPORT WireCellToolkitTargets
              LIBRARY  DESTINATION "${CMAKE_INSTALL_LIBDIR}"
              ARCHIVE  DESTINATION "${CMAKE_INSTALL_LIBDIR}"
              RUNTIME  DESTINATION "${CMAKE_INSTALL_BINDIR}")
      # ROOT dictionary (wct-ike.7) appends a generated source to this target.
      if(COMMAND wct_package_root_dict)
        wct_package_root_dict(${NAME} "${_dir}" "${_use}")
      endif()
    endif()
  endif()

  # ----- public headers -----
  if(IS_DIRECTORY "${_dir}/inc/${_hsub}")
    install(DIRECTORY "${_dir}/inc/${_hsub}"
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
            FILES_MATCHING PATTERN "*.h")
  endif()

  # ----- applications from apps/*.cxx -----
  if(IS_DIRECTORY "${_dir}/apps")
    file(GLOB _apps CONFIGURE_DEPENDS "${_dir}/apps/*.cxx")
    if(_apps)
      wct_resolve_uses(_appuse_tgts ${_app_use})
      foreach(_app IN LISTS _apps)
        get_filename_component(_aname "${_app}" NAME_WE)
        add_executable(${_aname} "${_app}")
        target_include_directories(${_aname} PRIVATE
          "${_dir}/apps"
          "$<BUILD_INTERFACE:${WCT_GENERATED_INCLUDE_DIR}>")
        set(_applink "")
        if(_have_lib)
          list(APPEND _applink ${NAME})
        endif()
        list(APPEND _applink ${_appuse_tgts})
        if(_applink)
          target_link_libraries(${_aname} PRIVATE ${_applink})
        endif()
        set_target_properties(${_aname} PROPERTIES
          INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR};${CMAKE_INSTALL_FULL_LIBDIR}")
        install(TARGETS ${_aname} RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
      endforeach()
    endif()
  endif()

  # ----- tests (wct-ike.9 provides the hook) -----
  if(WCT_WITH_TESTS AND COMMAND wct_register_tests)
    wct_register_tests("${NAME}" "${_dir}" "${_have_lib}" "${_test_use}")
  endif()
endfunction()
