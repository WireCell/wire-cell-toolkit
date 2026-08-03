# WireCellSpng - libtorch-based signal processing (mirrors spng/wscript_build:
# use='WireCellAux SPDLOG LIBTORCH').  Gated on LIBTORCH by the top-level file.
wct_package(WireCellSpng USE WireCellAux SPDLOG LIBTORCH)
