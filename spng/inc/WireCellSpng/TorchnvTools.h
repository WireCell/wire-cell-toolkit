#ifndef WIRECELLSPNG_TORCHNVTOOLS_H
#define WIRECELLSPNG_TORCHNVTOOLS_H

#ifdef HAVE_NVTX
    #include "nvToolsExt.h"
    #include <chrono>

    #define NVTX_RANGE_PUSH(name) nvtxRangePushA(name)
    #define NVTX_RANGE_POP() nvtxRangePop()
    #define NVTX_MARKER(name) nvtxMarkA(name)

    struct NvtxRange {
        NvtxRange(const char* name) {
            nvtxRangePushA(name);
        }
        ~NvtxRange() {
            nvtxRangePop();
        }
    };
    #define NVTX_SCOPED_RANGE(name) NvtxRange nvtxRangeInstance(name);
#else
    #define NVTX_RANGE_PUSH(name)
    #define NVTX_RANGE_POP()
    #define NVTX_MARKER(name)
    #define NVTX_SCOPED_RANGE(name)

#endif // HAVE_NVTX
#endif // WIRECELLSPNG_TORCHNVTOOLS_H