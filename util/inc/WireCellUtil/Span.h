// The "span" type.

#ifndef WIRECELL_SPAN
#define WIRECELL_SPAN

// fixme: Only need to keep this ifdef in place unil users upgrade to
// at least boost 1.78.  Can then remove this test from wscript.
#include "WireCellUtil/BuildConfig.h"
#ifdef HAVE_BOOST_CORE_SPAN_HPP
#include "boost/core/span.hpp"
#else
#include "WireCellUtil/boost/core/span.hpp"
#endif  

#endif
