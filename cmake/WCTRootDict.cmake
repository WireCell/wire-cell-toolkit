# Wire-Cell Toolkit - ROOT dictionary generation (epic wct-ike, task wct-ike.7).
#
# Provides the wct_package_root_dict() hook that wct_package() (wct-ike.4) calls
# after it creates a package's library target.  It reproduces
# waft/rootsys.py:gen_rootcling_dict():
#
#   rootcling -f <NAME>Dict.cxx -rml lib<NAME>.so -rmf lib<NAME>.rootmap \
#             <include flags> dict/LinkDef.h
#
# then compiles <NAME>Dict.cxx into lib<NAME> and installs the generated
# lib<NAME>.rootmap and <NAME>Dict_rdict.pcm next to the library so ROOT can
# autoload it at runtime.  Only packages with a dict/LinkDef.h and ROOTSYS in
# their use-list get a dictionary (currently just WireCellRoot).

include_guard(GLOBAL)

if(WCT_HAVE_ROOTSYS)
  find_program(WCT_ROOTCLING NAMES rootcling
    HINTS "${ROOT_BINDIR}" "$ENV{ROOTSYS}/bin" "${ROOTSYS}/bin")
endif()

function(wct_package_root_dict NAME pkgdir use)
  if(NOT WCT_HAVE_ROOTSYS)
    return()
  endif()
  if(NOT "ROOTSYS" IN_LIST use)
    return()
  endif()
  set(_linkdef "${pkgdir}/dict/LinkDef.h")
  if(NOT EXISTS "${_linkdef}")
    return()
  endif()
  if(NOT WCT_ROOTCLING)
    message(FATAL_ERROR "ROOT dictionary for ${NAME} needs rootcling but it was not found")
  endif()

  set(_dictsrc "${CMAKE_CURRENT_BINARY_DIR}/${NAME}Dict.cxx")
  set(_rootmap "${CMAKE_CURRENT_BINARY_DIR}/lib${NAME}.rootmap")
  set(_pcm     "${CMAKE_CURRENT_BINARY_DIR}/${NAME}Dict_rdict.pcm")

  # Include flags: the package's own headers, the generated BuildConfig dir and
  # ROOT's headers (so rootcling can parse anything LinkDef.h pulls in).
  set(_incflags "-I${pkgdir}/inc" "-I${WCT_GENERATED_INCLUDE_DIR}")
  foreach(_i IN LISTS ROOT_INCLUDE_DIRS)
    list(APPEND _incflags "-I${_i}")
  endforeach()

  add_custom_command(
    OUTPUT "${_dictsrc}" "${_rootmap}" "${_pcm}"
    COMMAND "${WCT_ROOTCLING}" -f "${_dictsrc}"
            -rml "lib${NAME}.so" -rmf "${_rootmap}"
            ${_incflags} "${_linkdef}"
    DEPENDS "${_linkdef}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Generating ROOT dictionary for ${NAME}"
    VERBATIM)

  target_sources(${NAME} PRIVATE "${_dictsrc}")
  install(FILES "${_rootmap}" "${_pcm}" DESTINATION "${CMAKE_INSTALL_LIBDIR}")
endfunction()
