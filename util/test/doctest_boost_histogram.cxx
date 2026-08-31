#include "WireCellUtil/Logging.h"

#include "WireCellUtil/doctest.h"

#include <boost/histogram.hpp>
#include <boost/format.hpp>


#include <algorithm> 
#include <vector>
#include <string>

using spdlog::debug;

namespace bh = boost::histogram;

//int main() {  
TEST_CASE("exercise two dimensional boost histogram") {

    auto h = bh::make_histogram(bh::axis::regular<>(10, 0, 10),
                                bh::axis::regular<>(10, 0, 1));

    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0, 3.0};
    std::vector<double> y = {.10, .20, .30, .40, .50, 0.3};

    for (size_t ind=0; ind<x.size(); ++ind) {
        h(x[ind], y[ind]);
    }
    CHECK(h.at(0,0) == 0);
    CHECK(h.at(1,1) == 1);

    auto indexed = bh::indexed(h);

    // // iterate over bins
    // for (auto&& x : indexed) {
    //     std::cerr << boost::format("bin %i [ %.1f, %.1f ): %i\n")
    //         % x.index() % x.bin().lower() % x.bin().upper() % *x;
    // }

    {
        // Compare cell values explicitly: the accessor's own comparison operators
        // are ambiguous under Apple clang (boost::histogram detail::operators).
        auto it = std::max_element(indexed.begin(), indexed.end(), [](const auto& a, const auto& b) {
            return static_cast<double>(*a) < static_cast<double>(*b);
        });
        const auto& cell = *it;
        CHECK(cell.index(0) == 3);
        CHECK(cell.index(1) == 3);
        CHECK(*cell == 2);
    }
//    return 0;
}
