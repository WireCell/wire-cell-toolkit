# Wire-Cell Toolkit - CTest integration (epic wct-ike, task wct-ike.9).
#
# Provides wct_register_tests(), the hook wct_package() (wct-ike.4) calls when
# WCT_WITH_TESTS is on.  Reproduces the test machinery of
# waft/smplpkgs.py:ValidationContext + waft/wcb_unit_test.py under CTest:
#
#   test/doctest*.cxx       -> aggregated into one wcdoctest-<pkg> runner
#   test/test_*.cxx         -> one program each, run as a test  (label: atomic)
#   test/atomic*.cxx        ->   "" (label: atomic)
#   test/check_*.cxx        -> built only, NOT run               (label: check)
#   test/history*.cxx       -> program + test                    (label: history)
#   test/report*.cxx        -> program + test                    (label: report)
#   test/{test,atomic,history,report}*.{py,sh,bats,jsonnet}
#                           -> run via python/bash/bats/wcsonnet (labelled)
#
# Test groups map to CTest LABELS: select with e.g. `ctest -L atomic`.  The waf
# --test-group {check,atomic,history,report} thus becomes label selection.
# Slow-test skipping (waf --test-duration) maps to per-test TIMEOUT/labels and
# can be driven with `ctest --timeout`.  A valgrind-style wrapper (waf
# --testcmd) maps to CMAKE_CROSSCOMPILING_EMULATOR-free `ctest -D
# ExperimentalMemCheck` or a TEST_LAUNCHER; we expose WCT_TEST_LAUNCHER for it.

include_guard(GLOBAL)
include(CTest)   # calls enable_testing()

find_program(WCT_BASH  NAMES bash)
find_program(WCT_BATS  NAMES bats PATHS "${CMAKE_SOURCE_DIR}/test/bats/bin")
# python3 for .py tests
if(NOT DEFINED Python3_EXECUTABLE)
  find_package(Python3 QUIET COMPONENTS Interpreter)
endif()

# Optional launcher prefix for every compiled test (e.g. "valgrind --error-exitcode=1");
# mirrors waf --testcmd.
set(WCT_TEST_LAUNCHER "" CACHE STRING "Command prefix to run each compiled test under (e.g. valgrind)")

# Common runtime environment for tests: find freshly built libs + apps and give
# bats its lib path (mirrors LD_LIBRARY_PATH injection + BATS_LIB_PATH in waf).
set(WCT_TEST_ENV
    "LD_LIBRARY_PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}:$ENV{LD_LIBRARY_PATH}"
    "BATS_LIB_PATH=${CMAKE_SOURCE_DIR}/test")

# Map a script extension to the command list that runs it.
function(_wct_script_command ext scriptfile out_cmd)
  if(ext STREQUAL ".py")
    set(${out_cmd} "${Python3_EXECUTABLE}" "${scriptfile}" PARENT_SCOPE)
  elseif(ext STREQUAL ".sh")
    set(${out_cmd} "${WCT_BASH}" "${scriptfile}" PARENT_SCOPE)
  elseif(ext STREQUAL ".bats")
    set(${out_cmd} "${WCT_BATS}" --jobs 1 "${scriptfile}" PARENT_SCOPE)
  elseif(ext STREQUAL ".jsonnet")
    # waf runs .jsonnet tests through the freshly built wcsonnet (see wcb.py).
    set(${out_cmd} "$<TARGET_FILE:wcsonnet>" "${scriptfile}" PARENT_SCOPE)
  else()
    set(${out_cmd} "" PARENT_SCOPE)
  endif()
endfunction()

function(wct_register_tests NAME pkgdir have_lib test_use)
  set(_testdir "${pkgdir}/test")
  if(NOT IS_DIRECTORY "${_testdir}")
    return()
  endif()

  # Libraries every test in this package links: the package lib + its test_use.
  wct_resolve_uses(_tlibs ${test_use})
  set(_link "")
  if(have_lib)
    list(APPEND _link ${NAME})
  endif()
  list(APPEND _link ${_tlibs})

  set(_inc "${pkgdir}/inc" "${_testdir}" "${WCT_GENERATED_INCLUDE_DIR}")

  # ---- aggregated doctest runner ----
  file(GLOB _dts CONFIGURE_DEPENDS "${_testdir}/doctest*.cxx")
  if(_dts)
    set(_main "${CMAKE_CURRENT_BINARY_DIR}/wcdoctest-${NAME}.cxx")
    configure_file("${CMAKE_SOURCE_DIR}/cmake/wcdoctest-main.cxx.in" "${_main}" @ONLY)
    set(_dtexe "wcdoctest-${NAME}")
    add_executable(${_dtexe} "${_main}" ${_dts})
    target_include_directories(${_dtexe} PRIVATE ${_inc})
    if(_link)
      target_link_libraries(${_dtexe} PRIVATE ${_link})
    endif()
    add_test(NAME "${NAME}/wcdoctest" COMMAND ${WCT_TEST_LAUNCHER} $<TARGET_FILE:${_dtexe}>)
    set_tests_properties("${NAME}/wcdoctest" PROPERTIES
      LABELS "atomic;doctest;${NAME}" ENVIRONMENT "${WCT_TEST_ENV}")
  endif()

  # ---- compiled programs by prefix/group ----
  # prefix -> (label, run?)
  foreach(_spec "test:atomic:1" "atomic:atomic:1" "check:check:0"
                "history:history:1" "report:report:1")
    string(REPLACE ":" ";" _s "${_spec}")
    list(GET _s 0 _prefix)
    list(GET _s 1 _label)
    list(GET _s 2 _run)
    file(GLOB _progs CONFIGURE_DEPENDS "${_testdir}/${_prefix}*.cxx")
    foreach(_src IN LISTS _progs)
      get_filename_component(_base "${_src}" NAME_WE)
      set(_tgt "${NAME}__${_base}")
      add_executable(${_tgt} "${_src}")
      target_include_directories(${_tgt} PRIVATE ${_inc})
      if(_link)
        target_link_libraries(${_tgt} PRIVATE ${_link})
      endif()
      if(_run)
        add_test(NAME "${NAME}/${_base}"
                 COMMAND ${WCT_TEST_LAUNCHER} $<TARGET_FILE:${_tgt}>
                 WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
        set_tests_properties("${NAME}/${_base}" PROPERTIES
          LABELS "${_label};${NAME}" ENVIRONMENT "${WCT_TEST_ENV}")
      endif()
    endforeach()
  endforeach()

  # ---- interpreted script tests (the check group runs no scripts) ----
  foreach(_prefix test atomic history report)
    foreach(_ext .py .sh .bats .jsonnet)
      file(GLOB _scripts CONFIGURE_DEPENDS "${_testdir}/${_prefix}*${_ext}")
      foreach(_script IN LISTS _scripts)
        get_filename_component(_base "${_script}" NAME)
        _wct_script_command("${_ext}" "${_script}" _cmd)
        if(_cmd)
          add_test(NAME "${NAME}/${_base}"
                   COMMAND ${_cmd}
                   WORKING_DIRECTORY "${pkgdir}")
          set_tests_properties("${NAME}/${_base}" PROPERTIES
            LABELS "${_prefix};script;${NAME}" ENVIRONMENT "${WCT_TEST_ENV}")
        endif()
      endforeach()
    endforeach()
  endforeach()
endfunction()
