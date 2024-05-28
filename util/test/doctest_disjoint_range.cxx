#include "WireCellUtil/DisjointRange.h"
// #define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include <random>
#include <chrono>

using namespace WireCell;
using spdlog::debug;


TEST_CASE("boost range basics")
{
    using value_type = int;
    using inner_vector = std::vector<value_type>;
    using outer_vector = std::vector<inner_vector>;
    using range_type = boost::sub_range<inner_vector>;
    using const_range_type = boost::sub_range<const inner_vector>;

    outer_vector ov = { {0,1,2}, {3}, {4} };
    range_type r( ov[0] );
    REQUIRE ( r.begin() != r.end() );
    REQUIRE ( r.size() == ov[0].size() );
    const inner_vector& civ = ov[0];

    // const range_type const_incorrect_fails (civ);

    // const range of const container
    const const_range_type cr (civ);
    // const range of non-const container
    const const_range_type cr2 (ov[0]); 

    const outer_vector& cov = ov;
    const const_range_type cr3 (cov[0]);     
    CHECK(cr3.size() == 3);

    CHECK(3 == std::distance(cr3.begin(), cr3.end()));

    CHECK(0 == *cr3.begin());
}

TEST_CASE("disjoint range basics")
{
    using value_type = int;
    using inner_vector = std::vector<value_type>;
    using outer_vector = std::vector<inner_vector>;
    using disjoint_type = disjoint_range<inner_vector>;
    outer_vector hier = { {0,1,2}, {3}, {4} };

    auto dj = disjoint_type();
    bool empty = dj.begin() == dj.end();
    REQUIRE( empty );
    dj.append( hier[0] );
    REQUIRE( dj.begin() != dj.end() );    
    CHECK( dj.size() == 3);

    inner_vector& iv = hier[1];
    dj.append( iv );
    CHECK( dj.size() == 4);

    dj.append( hier[2] );
    CHECK( dj.size() == 5);

    debug("flat index distance: {} -> {}", dj.begin().flat_index(), dj.end().flat_index());
    CHECK( dj.size() == std::distance(dj.begin(), dj.end()));

    //SUBCASE("ignoring const disjoint range on coordinate array test") { WARN(false); };

    // const auto& djc = dj;
    // CHECK( 5 == djc.size());
    // CHECK( 5 == std::distance(djc.begin(), djc.end()));

    // {
    //     auto it0 = djc.begin();
    //     CHECK( 0 == *it0 );
    // }

    CHECK(dj[4] == 4);
    CHECK(dj[2] == 2);
}


static
void test_sequence(std::vector<std::vector<int>> hier,
                   std::vector<int> want)
{
    using value_type = int;
    using inner_vector = std::vector<value_type>;
    // using outer_vector = std::vector<inner_vector>;
    using disjoint_type = disjoint_range<inner_vector>;
    disjoint_type numbers(hier);
    // auto numbers = flatten(hier);
    debug("end is end");
    REQUIRE(numbers.end() == numbers.end());
    {
        REQUIRE(numbers.size() == want.size());
        auto beg = numbers.begin();
        auto end = numbers.end();
        REQUIRE(beg != end);
        REQUIRE(beg == numbers.begin());
        REQUIRE(end == numbers.end());
        REQUIRE(std::distance(beg, end) == want.size());
    }
    {
        auto n2 = numbers;
        REQUIRE(n2.size() == want.size());
        REQUIRE(std::distance(n2.begin(),n2.end()) == want.size());
    }
    {
        auto it = numbers.begin();
        REQUIRE(it == numbers.begin());
        REQUIRE(it.flat_index() == 0);
        debug("not at end is not at end");
        REQUIRE(it != numbers.end());
        it += want.size(); // now at end.
        REQUIRE(it.flat_index() == want.size());
        debug("jump to end is end, it={}, end={}", it.flat_index(), numbers.end().flat_index());
        REQUIRE(it == numbers.end());
    }
    {
        auto it = numbers.begin();
        for (size_t left=want.size(); left; --left) {
            REQUIRE(it != numbers.end());
            ++it;
        }
        REQUIRE(it == numbers.end());        
    }

    {
        // can also use major/minor indices.
        CHECK(numbers[{0,0}] == want[0]);
        CHECK(numbers.at({0,0}) == want.at(0));
    }


    size_t ind = 0;
    for (const auto& num : numbers) {
        CHECK(ind < 5);
        CHECK(num == want[ind]);
        ++ind;
    }
    disjoint_type dj(hier);
    // auto dj = flatten(hier);
    auto djit = dj.begin();     // 0
    auto end = dj.end();
    REQUIRE(djit == dj.begin());
    REQUIRE(djit.flat_index() == 0);
    ++djit;                     // 1
    REQUIRE(djit.flat_index() == 1);
    REQUIRE(djit != dj.begin());
    REQUIRE(*djit == want[1]);
    {
        auto& aa = dj.at(4);
        REQUIRE(aa == want[4]);
    }
    --djit;                     // 0
    REQUIRE(djit.flat_index() == 0);
    REQUIRE(*djit == want[0]);
    djit += 2;                  // 2
    REQUIRE(djit.flat_index() == 2);
    REQUIRE(*djit == want[2]);    
    djit -= 2;                  // 0
    REQUIRE(djit.flat_index() == 0);
    REQUIRE(*djit == want[0]);
    djit += want.size();        // 5
    REQUIRE(djit.flat_index() == 5);
    REQUIRE(djit == end);
    --djit;                     // backwards from end
    REQUIRE(djit.flat_index() == 4);
    REQUIRE(djit != end);
    debug("want.back={}", want.back());
    debug("last flat={} major={} minor={}", djit.flat_index(), djit.index().first, djit.index().second);
    debug("last djit={}", *djit);
    REQUIRE(*djit == want.back());
    djit -= want.size() - 1;
    REQUIRE(*djit == want.front());
    CHECK_THROWS_AS(--djit, std::out_of_range);
}

TEST_CASE("disjoint range iterator no empties") {
    test_sequence({ {0,1,2}, {3}, {4} }, {0,1,2,3,4});
}
TEST_CASE("disjoint range iterator with empties") {
    test_sequence({ {}, {0,1,2}, {}, {3}, {4}, {}}, {0,1,2,3,4});
}
TEST_CASE("disjoint range iterator mutate") {
    using value_type = int;
    using inner_vector = std::vector<value_type>;
    using disjoint_type = disjoint_range<inner_vector>;
    std::vector<std::vector<int>> hier = { {0,1,2}, {3}, {4} };
    disjoint_type numbers(hier);
    // auto djit = flatten(hier);
    disjoint_type::iterator djit = numbers.begin();
    CHECK(*djit == 0);
    *djit = 42;
    CHECK(*djit == 42);
    CHECK(hier[0][0] == 42);
}

TEST_CASE("disjoint range iterator const") {
    using vi = std::vector<int>;
    using vvi = std::vector<vi>;

    using dji_vvivc = disjoint_range<vi>;
    vvi hier = { {0,1,2}, {3}, {4} };
  
    auto vvivc_range = dji_vvivc(hier);
    auto vvivc = vvivc_range.begin();
    ++vvivc;
    int val = *vvivc;
    CHECK(val == 1);

    // should not compile!
    int& ref2 = *vvivc;
    CHECK(ref2 == 1);

    int const& ref = *vvivc;
    CHECK(ref == 1);

    // using dji_vvcvc = disjoint<vvi::const_iterator, vi::const_iterator>;
    // auto vvcvc = dji_vvcvc(hier.cbegin(), hier.cend());
    // ++vvcvc;

}

TEST_CASE("disjoint range iterator perf") {
    using value_type = size_t;
    using inner_type = std::vector<value_type>;
    using outer_type = std::vector<inner_type>;
    using disjoint_type = disjoint_range<inner_type>;

    const size_t onums =1000;
    const size_t inums =1000;
    const size_t nrands=10000;
    // const size_t nrands=1000;

    std::random_device rnd_device;
    std::mt19937 rnd_engine {rnd_device()};  // Generates random integers
    std::uniform_int_distribution<size_t> dist {0, onums * inums};

    size_t count = 0;
    outer_type outer;
    for (size_t onum=0; onum < onums; ++onum) {
        inner_type inner(inums);
        for (size_t inum=0; inum < inums; ++inum) {
            inner[inum] = count;
            ++count;
        }
        outer.push_back(inner);
    }
    
    disjoint_type dj(outer);
    size_t ntot = onums * inums;
    REQUIRE(ntot == count);
    REQUIRE(ntot == dj.size());

    {
        std::vector<size_t> indices = {0, 1, 2, inums-1, inums, 10*inums, ntot-1};
        for (size_t ind : indices) {
            size_t got = dj[ind];
            REQUIRE(ind == got);
        }
    }

    // random access flat dj many times
    std::vector<size_t> rands(nrands);
    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t dummy = 0; dummy<nrands; ++dummy) {
        const size_t rind = dist(rnd_engine);
        rands[dummy] = rind;
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt10 = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

    for (size_t rind : rands) {
        if (0 > rind || rind >= ntot) {
            REQUIRE (rind >= 0);
            REQUIRE (rind < ntot);
        }
        const size_t got = dj[rind];
        if (got == rind) continue;
        REQUIRE(got == rind);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto dt21 = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    std::sort(rands.begin(), rands.end());

    auto t3 = std::chrono::high_resolution_clock::now();
    auto dt32 = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
    for (size_t rind : rands) {
        if (0 > rind || rind >= ntot) {
            REQUIRE (rind >= 0);
            REQUIRE (rind < ntot);
        }
        const size_t got = dj[rind];
        if (got == rind) continue;
        REQUIRE(got == rind);
    }
    auto t4 = std::chrono::high_resolution_clock::now();
    auto dt43 = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();

    debug("{} us: make {} rands", dt10, nrands);
    debug("{} us: random index", dt21);
    debug("{} us: sorting", dt32);
    debug("{} us: sorted index", dt43);
    debug("{} MHz: sort and sum", (1.0*nrands)/(dt32+dt43));

    auto djit = dj.begin();
    for (size_t ind = 0; ind<ntot; ++ind) {
        if (ind == *djit) {
            ++djit;
            continue;
        }
        REQUIRE(ind == *djit);
    }
    REQUIRE( djit == dj.end() );
    
}

