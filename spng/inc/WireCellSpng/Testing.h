/** Any SPNG unit test code using doctest MUST #include this header and SHOULD
 * NOT include the "Torch.h" header (its' included here).
 *
 * SPNG library could SHOULD NOT include this header but may include Torch.h.
 */

#ifndef WIRECELL_SPNG_TESTING
#define WIRECELL_SPNG_TESTING

#include "WireCellSpng/Torch.h"
#undef CHECK
#include "WireCellUtil/doctest.h"


#endif
