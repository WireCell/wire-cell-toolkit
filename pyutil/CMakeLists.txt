# WCPPyUtil (mirrors pyutil/wscript_build): a small library plus the Python
# modules that must be on PYTHONPATH at runtime.
# TEST_USE: the doctests (test/doctest_*.cxx) need WireCellUtil's doctest.h and logging.
wct_package(WCPPyUtil USE PYTHON TEST_USE WireCellUtil)

# SCN_Vertex.py is the top-level inference script the C++ SCN_Vertex wrapper
# calls; SCN/DeepVtx.py holds the sparseconvnet model class it imports.
# DV2_Vertex.py and its vendored package dv2max_vertexer/ are installed by the
# same rule (subdirectories included).
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python/"
        DESTINATION "${CMAKE_INSTALL_PREFIX}/python"
        FILES_MATCHING PATTERN "*.py")
