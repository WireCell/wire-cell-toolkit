# WCPPyUtil (mirrors pyutil/wscript_build): a small library plus the Python
# modules that must be on PYTHONPATH at runtime.
wct_package(WCPPyUtil USE PYTHON)

# SCN_Vertex.py is the top-level inference script the C++ SCN_Vertex wrapper
# calls; SCN/DeepVtx.py holds the sparseconvnet model class it imports.
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python/"
        DESTINATION "${CMAKE_INSTALL_PREFIX}/python"
        FILES_MATCHING PATTERN "*.py")
