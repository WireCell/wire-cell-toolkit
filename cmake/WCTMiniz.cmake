# Wire-Cell Toolkit - isolate the vendored miniz (zip/deflate) implementation.
#
# WCT vendors miniz (util/src/miniz.cxx, an amalgamated third-party zip
# library) for its .npz / custard I/O.  When WCT is built with LIBTORCH, this
# collides with libtorch, which *also* statically bundles its own copy of
# miniz.  Both expose the same global C symbols (mz_*, tinfl_*, tdefl_*).
# libWireCellUtil is loaded before libtorch, so torch's internal calls into
# what it believes is its own miniz get bound by the dynamic linker to WCT's
# copy instead.  The two miniz versions lay out mz_zip_archive differently, so
# torch::jit::load() dereferences a wild pointer and segfaults inside
# PyTorchStreamReader::init().
#
# Fix: compile miniz once into a static library with *hidden* symbol
# visibility, so the symbols resolve within each consuming WCT shared object at
# link time but never enter its dynamic symbol table.  Nothing WCT builds then
# exports miniz to the global scope, and torch binds to its own copy.  The
# static library is linked PRIVATE into every WCT library by wct_package();
# archive semantics mean only libraries that actually reference miniz
# (util, sio, clus, img, ...) pull in its objects.

include_guard(GLOBAL)

set(_wct_miniz_src "${CMAKE_SOURCE_DIR}/util/src/miniz.cxx")
if(EXISTS "${_wct_miniz_src}")
  add_library(wct_miniz STATIC "${_wct_miniz_src}")
  # miniz.cxx does #include "miniz.h" from the same directory.
  target_include_directories(wct_miniz PRIVATE "${CMAKE_SOURCE_DIR}/util/src")
  set_target_properties(wct_miniz PROPERTIES
    POSITION_INDEPENDENT_CODE ON      # linked into shared objects
    C_VISIBILITY_PRESET hidden
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON)
else()
  message(WARNING "WCTMiniz: ${_wct_miniz_src} not found; miniz isolation disabled")
endif()
