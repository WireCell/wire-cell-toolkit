# Determine the Wire-Cell Toolkit version and build mode.
#
# This mirrors determine_version()/is_development() in the top-level wscript so
# the CMake build bakes the same version string into BuildConfig.h and applies
# the same development-vs-release policy.  It is a pure helper: include() it
# before project() and it sets, in the including scope:
#
#   WCT_VERSION            full version string (e.g. 0.35.0-330-g61618538 or
#                          <branch>-<describe> on a feature branch)
#   WCT_VERSION_TRIPLET    numeric X.Y.Z extracted from WCT_VERSION, suitable
#                          for project(... VERSION ...)
#   WCT_IS_DEVELOPMENT     TRUE on a development (non-release) checkout
#
# Sources of the version, in priority order (same as wscript):
#   1. `git describe --tags` when building from a git clone
#   2. a version.txt file shipped in a source release archive
#
# WCT_BUILD_MODE (cache, default "") may force "development" or "release";
# otherwise the mode is detected from the git branch / version, exactly as the
# wscript does.

function(_wct_extract_triplet version out_var)
  # Pull a numeric X.Y.Z (or X.Y -> X.Y.0) out of a git-describe string.
  if(version MATCHES "([0-9]+)\\.([0-9]+)\\.([0-9]+)")
    set(${out_var} "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.${CMAKE_MATCH_3}" PARENT_SCOPE)
  elseif(version MATCHES "([0-9]+)\\.([0-9]+)")
    set(${out_var} "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.0" PARENT_SCOPE)
  else()
    set(${out_var} "0.0.0" PARENT_SCOPE)
  endif()
endfunction()

function(_wct_determine_version)
  set(_version "")
  set(_branch "")

  find_package(Git QUIET)
  if(GIT_FOUND)
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" describe --tags
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE _describe
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _git_rc
      ERROR_QUIET)
    if(_git_rc EQUAL 0 AND _describe)
      set(_version "${_describe}")
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE _branch
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    endif()
  endif()

  if(NOT _version)
    # Not a git clone: fall back to a shipped version.txt (release archive).
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/version.txt")
      file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/version.txt" _lines)
      list(GET _lines 0 _version)
      string(STRIP "${_version}" _version)
    else()
      message(FATAL_ERROR
        "Wire-Cell Toolkit must be built from a git clone or a source "
        "distribution containing version.txt")
    endif()
  endif()

  # On master or a numeric release branch the describe string is the version;
  # feature branches are tagged <branch>-<describe> (mirrors the wscript).
  string(SUBSTRING "${_branch}" 0 1 _branch0)
  if(_branch STREQUAL "master" OR _branch0 MATCHES "[0-9]")
    set(_full "${_version}")
  elseif(_branch)
    set(_full "${_branch}-${_version}")
  else()
    set(_full "${_version}")
  endif()

  # Development = not master and not a numeric release version (per wscript).
  string(SUBSTRING "${_full}" 0 1 _full0)
  if(_full STREQUAL "master" OR _full0 MATCHES "[0-9]")
    set(_dev FALSE)
  else()
    set(_dev TRUE)
  endif()

  # Allow an explicit override.
  if(WCT_BUILD_MODE STREQUAL "development")
    set(_dev TRUE)
  elseif(WCT_BUILD_MODE STREQUAL "release")
    set(_dev FALSE)
  elseif(NOT WCT_BUILD_MODE STREQUAL "")
    message(FATAL_ERROR "WCT_BUILD_MODE must be 'development', 'release' or empty (auto)")
  endif()

  _wct_extract_triplet("${_full}" _triplet)

  set(WCT_VERSION         "${_full}"    PARENT_SCOPE)
  set(WCT_VERSION_TRIPLET "${_triplet}" PARENT_SCOPE)
  set(WCT_IS_DEVELOPMENT  ${_dev}       PARENT_SCOPE)
endfunction()

_wct_determine_version()
