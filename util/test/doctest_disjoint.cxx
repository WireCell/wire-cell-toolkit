#include "WireCellUtil/DisjointRange.h"
// #define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Type.h"

#include <random>
#include <chrono>

using namespace WireCell;
using spdlog::debug;

template<typename C, typename R, typename I>
void dump_types(C& c, R& r, const I& i)
{
    debug("collection: {}", type(c));
    debug("range:      {}", type(r));
    debug("iterator:   {}", type(i));
}

template<typename V>
void dump_val(V& v)
{
    debug("value: ({}&) {}", type(v), v);
}

template<typename collection>
void disjoint_const_correctness(collection const& c)
{
    using inner_vector = typename collection::iterator::value_type;

    {
        collection col = c;
        disjoint_range<inner_vector> r(col);
        dump_types(col, r, r.begin());
        REQUIRE(r.size() == 5);
        for (auto& one : r) {
            dump_val(one);
            one = 42;
        }
        CHECK(42 == col[0][0]);
        CHECK(42 == r[0]);
        *r.begin() = 43;
        CHECK(43 == col[0][0]);
    }
    {
        collection const& col = c;
        disjoint_range<inner_vector const> const r(col);
        dump_types(col, r, r.begin());
        REQUIRE(r.size() == 5);
        for (auto& one : r) {   // implicitly const
            dump_val(one);
            // one = 42;
        }
        CHECK(0 == col[0][0]);
    }    
    {
        collection col = c;
        disjoint_range<inner_vector> const r(col);
        dump_types(col, r, r.begin());
        REQUIRE(r.size() == 5);
        for (auto& one : r) {   // implicitly const
            dump_val(one);
            // one = 42;
        }
    }
    // {
    //     collection const& col = c;
    //     disjoint_range<inner_vector const> r;
    //     for (const auto& one : col) {
    //         r.append(one.begin(), one.end());
    //     }
    //     std::cerr << "colection: " << type(col) << "\nrange: " << type(r) << "\niter: " << type(r.begin()) << "\n";
    //     REQUIRE(r.size() == 5);
    //     for (const auto& one : r) {
    //         std::cerr << "(" << type(one) << " const&) "<<one<<"\n";
    //         // one = 42;
    //     }
    // }
}

TEST_CASE("disjoint constness")
{
    using value_type = size_t;
    using collection = std::vector<std::vector<value_type>>;

    collection col = { {0,1,2}, {3}, {4} };
    SUBCASE("const test skipped") { WARN(false); }
    // disjoint_const_correctness(col);


}
