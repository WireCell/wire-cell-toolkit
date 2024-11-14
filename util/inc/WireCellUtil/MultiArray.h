// Multi-dimensional array support.
//
// For now, this simply brings in boost::multi_array in a way that ignores its
// internal use of deprecated functions.

#ifndef WIRECELL_MULTIARRAY
#define WIRECELL_MULTIARRAY

#ifdef __clang__
#  if defined(__has_warning)
#    define HAS_WARNING(warning) __has_warning(warning)
#  else
#    define HAS_WARNING(warning) 1
#  endif
#else
#  define HAS_WARNING(warning) 1
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#if HAS_WARNING("-Wmaybe-uninitialized")
#pragma GCC diagnostic warning "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#pragma GCC diagnostic warning "-Warray-bounds"
#pragma GCC diagnostic ignored "-Warray-bounds"

#include <boost/multi_array.hpp>

#pragma GCC diagnostic pop

#endif
