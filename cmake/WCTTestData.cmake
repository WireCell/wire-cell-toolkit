# Wire-Cell Toolkit - test-data repository provisioning (epic wct-ike, task wct-ike.10).
#
# Mirrors waft/datarepo.py.  Provides:
#   - options WCT_TEST_DATA_URL / WCT_TEST_DATA_VERSIONS (--test-data-url/-versions)
#   - target `test-data`       : download+unpack input.tar and history-<ver>.tar
#                                into ${WCT_TEST_DATA_DIR} (build/tests)
#   - target `pack-test-data`  : re-archive them (./wcb packrepo)
#   - CTest fixture WCTTestData : a setup test so `ctest -L history` can provision
#                                on demand; data-dependent tests opt in by setting
#                                FIXTURES_REQUIRED WCTTestData.
# The history/report test groups and some atomic tests read build/tests/.

include_guard(GLOBAL)

set(WCT_TEST_DATA_URL "https://www.phy.bnl.gov/~bviren/tmp/wcttest/data_repo"
    CACHE STRING "Base URL of the WCT test-data repository (mirrors --test-data-url)")
set(WCT_TEST_DATA_VERSIONS "0.20.0,0.21.0,0.22.0,0.23.0,0.24.1"
    CACHE STRING "Release versions whose history archives to fetch (--test-data-versions)")
set(WCT_TEST_DATA_DIR "${CMAKE_BINARY_DIR}/tests"
    CACHE PATH "Directory under which the test-data repo is unpacked")

set(_dl "${CMAKE_SOURCE_DIR}/cmake/download-test-data.cmake")
set(_pk "${CMAKE_SOURCE_DIR}/cmake/pack-test-data.cmake")

add_custom_target(test-data
  COMMAND "${CMAKE_COMMAND}"
          -DURL=${WCT_TEST_DATA_URL}
          -DVERSIONS=${WCT_TEST_DATA_VERSIONS}
          -DDEST=${WCT_TEST_DATA_DIR}
          -P "${_dl}"
  COMMENT "Downloading WCT test-data repository into ${WCT_TEST_DATA_DIR}"
  VERBATIM USES_TERMINAL)

add_custom_target(pack-test-data
  COMMAND "${CMAKE_COMMAND}"
          -DSRC=${WCT_TEST_DATA_DIR}
          -DVERSIONS=${WCT_TEST_DATA_VERSIONS}
          -P "${_pk}"
  COMMENT "Packing WCT test-data archives from ${WCT_TEST_DATA_DIR}"
  VERBATIM)

# A CTest fixture-setup test that provisions the data; a data-dependent test
# requires it via set_tests_properties(<t> PROPERTIES FIXTURES_REQUIRED WCTTestData).
add_test(NAME WCTTestData-provision
         COMMAND "${CMAKE_COMMAND}"
                 -DURL=${WCT_TEST_DATA_URL}
                 -DVERSIONS=${WCT_TEST_DATA_VERSIONS}
                 -DDEST=${WCT_TEST_DATA_DIR}
                 -P "${_dl}")
set_tests_properties(WCTTestData-provision PROPERTIES FIXTURES_SETUP WCTTestData)
