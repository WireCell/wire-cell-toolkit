/** This API holds some general "mathy" stuff.
 */

#ifndef WIRECELL_UTIL_MATH
#define WIRECELL_UTIL_MATH

#include <boost/math/constants/constants.hpp>
#include <algorithm>

namespace WireCell {

    // Return the greatest common denominator between a and b.
    size_t GCD(size_t a, size_t b);

    // Return coprime of number nearest to target or zero if none
    // found smaller than number.  Search is bound by [0, number/2]
    // or [number/2, number] depending on which half target is in.
    size_t nearest_coprime(size_t number, size_t target);

    const double pi = boost::math::constants::pi<double>();

}
#endif
