/** Functions to help use ranges */

#ifndef WIRECELL_RANGE
#define WIRECELL_RANGE
#include "boost/range/irange.hpp"
namespace WireCell::Range {

    // convert from first/second to begin/end
    // https://stackoverflow.com/a/27140916
    template <typename I>
    struct iter_pair : std::pair<I, I>
    {
        using std::pair<I, I>::pair;
        I begin() { return this->first; }
        I end() { return this->second; }
    };
    // and use like:
    // for (auto e : make_range(my_pair)) { ... }
    template<typename I>
    iter_pair<I> make_range(std::pair<I, I> p) {
        return iter_pair<I>(p);
    }


    template<typename I>
    auto irange(std::pair<I, I> p) {
        return boost::irange(p.first, p.second);
    }

    using boost::irange;
    
}

#endif
