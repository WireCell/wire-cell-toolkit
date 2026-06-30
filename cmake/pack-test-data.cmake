# Wire-Cell Toolkit - re-pack the test-data repo (epic wct-ike, task wct-ike.10).
#
# Run in script mode:
#   cmake -DSRC=<tests-dir> -DVERSIONS=<csv> -P pack-test-data.cmake
#
# Mirrors `./wcb packrepo`: tar up SRC/input -> input.tar and each
# SRC/history/<ver> -> history-<ver>.tar (only versions in VERSIONS, if given).

if(NOT SRC)
  message(FATAL_ERROR "pack-test-data: -DSRC=<tests-dir> is required")
endif()

if(EXISTS "${SRC}/input")
  message(STATUS "datarepo: packing input.tar")
  execute_process(COMMAND "${CMAKE_COMMAND}" -E tar cf "${SRC}/../input.tar" input
                  WORKING_DIRECTORY "${SRC}")
else()
  message(WARNING "datarepo: no input/ under ${SRC}")
endif()

if(VERSIONS)
  string(REPLACE "," ";" _vers "${VERSIONS}")
else()
  file(GLOB _vdirs RELATIVE "${SRC}/history" "${SRC}/history/*")
  set(_vers ${_vdirs})
endif()
foreach(_v IN LISTS _vers)
  if(_v AND IS_DIRECTORY "${SRC}/history/${_v}")
    message(STATUS "datarepo: packing history-${_v}.tar")
    execute_process(COMMAND "${CMAKE_COMMAND}" -E tar cf "${SRC}/../history-${_v}.tar" "history/${_v}"
                    WORKING_DIRECTORY "${SRC}")
  endif()
endforeach()
