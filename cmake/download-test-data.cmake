# Wire-Cell Toolkit - test-data downloader (epic wct-ike, task wct-ike.10).
#
# Run in script mode:
#   cmake -DURL=<repo> -DVERSIONS=<csv> -DDEST=<dir> -P download-test-data.cmake
#
# Mirrors waft/datarepo.py: fetches <URL>/input.tar and <URL>/history-<ver>.tar
# for each version and extracts them under DEST (so DEST/input/ and
# DEST/history/<ver>/ appear).  TLS verification is disabled to match the
# unverified SSL context datarepo uses (the BNL host cert does not parse
# cleanly), see waft/datarepo.py.

if(NOT DEST)
  message(FATAL_ERROR "download-test-data: -DDEST=<dir> is required")
endif()
file(MAKE_DIRECTORY "${DEST}")

function(_wct_get name)
  set(_tar "${DEST}/${name}.tar")
  message(STATUS "datarepo: downloading ${URL}/${name}.tar")
  file(DOWNLOAD "${URL}/${name}.tar" "${_tar}" STATUS _st TLS_VERIFY OFF SHOW_PROGRESS)
  list(GET _st 0 _rc)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "datarepo: download of ${name}.tar failed: ${_st}")
  endif()
  file(ARCHIVE_EXTRACT INPUT "${_tar}" DESTINATION "${DEST}")
  file(REMOVE "${_tar}")
endfunction()

# The input archive (unversioned).
_wct_get(input)

# The per-release history archives.
if(VERSIONS)
  string(REPLACE "," ";" _vers "${VERSIONS}")
  foreach(_v IN LISTS _vers)
    if(_v)
      _wct_get("history-${_v}")
    endif()
  endforeach()
endif()
