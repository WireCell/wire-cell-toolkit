# WireCellCfg (mirrors cfg/wscript_build): installs Jsonnet configuration into
# share/wirecell, preserving directory structure (relative_trick).  No library.
wct_package(WireCellCfg)

set(_share "${CMAKE_INSTALL_PREFIX}/share/wirecell")

# Always: top-level *.jsonnet plus the common/ and ui/ pgrapher trees.
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/wirecell.jsonnet
        DESTINATION "${_share}" OPTIONAL)
file(GLOB _top_jsonnet "${CMAKE_CURRENT_SOURCE_DIR}/*.jsonnet")
install(FILES ${_top_jsonnet} DESTINATION "${_share}")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/pgrapher/common"
        DESTINATION "${_share}/pgrapher" FILES_MATCHING PATTERN "*.jsonnet")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/pgrapher/ui"
        DESTINATION "${_share}/pgrapher" FILES_MATCHING PATTERN "*.jsonnet")

# Optionally: per-experiment configs selected by WCT_INSTALL_CONFIG.
if(WCT_INSTALL_CONFIG)
  string(REPLACE "," ";" _exps "${WCT_INSTALL_CONFIG}")
  foreach(_exp IN LISTS _exps)
    if(_exp STREQUAL "all")
      set(_exp "*")
    endif()
    file(GLOB _expdirs "${CMAKE_CURRENT_SOURCE_DIR}/pgrapher/experiment/${_exp}")
    foreach(_d IN LISTS _expdirs)
      get_filename_component(_dn "${_d}" NAME)
      install(DIRECTORY "${_d}"
              DESTINATION "${_share}/pgrapher/experiment"
              FILES_MATCHING PATTERN "*.jsonnet" PATTERN "*.fcl")
    endforeach()
  endforeach()
endif()
