# WireCellUtil - the base toolkit library (mirrors util/wscript_build).
# util is the foundation and must not "use" any other WireCell library.

wct_package(WireCellUtil
            USE SPDLOG BOOST FFTW EIGEN DYNAMO JSONCPP JSONNET ZLIB)

# util ships more than inc/WireCellUtil/*.h: also the vendored *.hpp headers and
# the custard/ and boost/ header trees (see util/wscript_build).  wct_package()
# installed the *.h; add the rest here.
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/inc/WireCellUtil"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        FILES_MATCHING
        PATTERN "*.hpp")
