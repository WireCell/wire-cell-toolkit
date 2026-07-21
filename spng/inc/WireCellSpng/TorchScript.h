
#ifndef WIRECELL_SPNG_TORCHSCRIPT
#define WIRECELL_SPNG_TORCHSCRIPT

#include "WireCellSpng/Torch.h"

// Torch has a HUGE number of compiler warnings and apparently telling the
// compiler they are "system headers" aka could be found via "-isystem" instead
// of "-I" CLI flags, is a valid solution to quell them.  We can not easily tell
// waf to change to "-isystem" so rely on this pragma.

// #ifdef __clang__
// #  if defined(__has_warning)
// #    define HAS_WARNING(warning) __has_warning(warning)
// #  else
// #    define HAS_WARNING(warning) 1
// #  endif
// #else
// #  define HAS_WARNING(warning) 1
// #endif

// #if HAS_WARNING("-Wgnu-zero-variadic-macro-arguments")
// #pragma GCC diagnostic push
// #pragma GCC diagnostic warning "-Wgnu-zero-variadic-macro-arguments"
// #pragma GCC diagnostic ignored "-Wgnu-zero-variadic-macro-arguments"
// #endif

// #if HAS_WARNING("-Wc++20-extensions")
// #pragma GCC diagnostic push
// #pragma GCC diagnostic warning "-Wc++20-extensions"
// #pragma GCC diagnostic ignored "-Wc++20-extensions"
// #endif

#include <torch/script.h>

// #if HAS_WARNING("-Wc++20-extensions")
// #pragma GCC diagnostic pop
// #endif

// #if HAS_WARNING("-Wgnu-zero-variadic-macro-arguments")
// #pragma GCC diagnostic pop
// #endif

#endif
