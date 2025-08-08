// Main inclusion of Eigen

#ifndef WIRECELL_EIGEN
#define WIRECELL_EIGEN

#ifdef __clang__
#  if defined(__has_warning)
#    define HAS_WARNING(warning) __has_warning(warning)
#  else
#    define HAS_WARNING(warning) 1
#  endif
#else
#  define HAS_WARNING(warning) 1
#endif

#if HAS_WARNING("-Wmaybe-uninitialized")
#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#if HAS_WARNING("-Wignored-attributes")
#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wignored-attributes"
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
#if HAS_WARNING("-Wunused-but-set-variable")
#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues> 

#if HAS_WARNING("-Wunused-but-set-variable")
#pragma GCC diagnostic pop
#endif
#if HAS_WARNING("-Wignored-attributes")
#pragma GCC diagnostic pop
#endif
#if HAS_WARNING("-Wmaybe-uninitialized")
#pragma GCC diagnostic pop
#endif




#endif
