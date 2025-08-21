/** SPNG library code must not #include torch headers directly and only #include this one.

    Note, any doctest code should see Testing.h.
*/

#ifndef WIRECELL_SPNG_TORCH
#define WIRECELL_SPNG_TORCH

// Torch has a HUGE number of compiler warnings and apparently telling the
// compiler they are "system headers" aka could be found via "-isystem" instead
// of "-I" CLI flags, is a valid solution to quell them.  We can not easily tell
// waf to change to "-isystem" so rely on this pragma.

#pragma GCC system_header
#include <torch/torch.h>  // should be only include site in entire spng!

#endif
