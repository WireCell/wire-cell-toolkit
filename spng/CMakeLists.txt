# WireCellSpng - libtorch-based signal processing (mirrors spng/wscript_build:
# use='WireCellAux SPDLOG LIBTORCH').  Gated on LIBTORCH by the top-level file.
wct_package(WireCellSpng USE WireCellAux SPDLOG LIBTORCH)

# The check_threads exploratory program (spng-fra.13) needs TBB to reproduce the
# "torch op on a TBB worker thread" scenario.  Link it in when available; the
# program itself is guarded by HAVE_TBB.
if(WCT_WITH_TESTS AND TARGET WireCellSpng__check_threads)
  find_package(TBB QUIET)
  if(TBB_FOUND)
    target_link_libraries(WireCellSpng__check_threads PRIVATE TBB::tbb)
    target_compile_definitions(WireCellSpng__check_threads PRIVATE HAVE_TBB=1)
  endif()
endif()
