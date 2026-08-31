# WireCellClus (mirrors clus/wscript_build): WCPPyUtil is added only when
# Python is available, exactly as the wscript_build keys off INCLUDES_PYTHON.
set(_clus_use WireCellAux WCPQuickhull)
if(WCT_HAVE_PYTHON)
  list(APPEND _clus_use WCPPyUtil)
endif()
wct_package(WireCellClus USE ${_clus_use})
