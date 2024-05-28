// Multi-dimensional array support.
//
// For now, this simply brings in boost::multi_array in a way that ignores its
// internal use of deprecated functions.

#ifndef WIRECELL_MULTIARRAY
#define WIRECELL_MULTIARRAY

#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic warning "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic warning "-Warray-bounds"
#pragma GCC diagnostic ignored "-Warray-bounds"

#include <boost/multi_array.hpp>

#pragma GCC diagnostic pop

#endif
