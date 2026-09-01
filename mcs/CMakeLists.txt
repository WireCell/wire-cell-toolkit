# WireCellMcs (mirrors mcs/wscript_build): a plain numeric library, NOT a
# plugin -- never list "WireCellMcs" in a jsonnet plugins: entry.  clus links
# it (see clus/CMakeLists.txt), so mcs must be listed before clus in
# WCT_PACKAGES.
wct_package(WireCellMcs USE WireCellUtil EIGEN)
