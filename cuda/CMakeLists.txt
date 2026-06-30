# WireCellCuda (mirrors cuda/wscript_build).  The library is built from src/*.cu
# via the CUDA language support enabled in wct-ike.8; until that lands this
# package has no compilable sources under the plain wct_package() globbing.
wct_package(WireCellCuda USE CUDA)
