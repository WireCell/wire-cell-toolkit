#ifndef WIRECELLTBB_NVTXTOOLS_H
#define WIRECELLTBB_NVTXTOOLS_H

#ifdef HAVE_NVTX
    #include "nvToolsExt.h"

    #define NVTX_RANGE_PUSH(name) nvtxRangePushA(name)
    #define NVTX_RANGE_POP() nvtxRangePop()
#else
    #define NVTX_RANGE_PUSH(name)
    #define NVTX_RANGE_POP()
#endif // HAVE_NVTX
#endif // WIRECELLTBB_NVTXTOOLS_H
