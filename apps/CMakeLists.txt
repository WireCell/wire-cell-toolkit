# WireCellApps (mirrors apps/wscript_build): builds the apps/ executables
# (wire-cell, wcsonnet, wcwires).  TBB is linked only when present -- the
# wscript_build appends it under 'if HAVE_TBB'; wct_resolve_uses drops an absent
# optional external automatically, so passing it unconditionally is equivalent.
wct_package(WireCellApps USE WireCellIface FFTWTHREADS TBB)
