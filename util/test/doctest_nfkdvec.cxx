#include "WireCellUtil/NFKDVec.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Logging.h"

// #define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "WireCellUtil/doctest.h"

#include <boost/container_hash/hash.hpp>

#include <random>
#include <vector>

using spdlog::debug;
using namespace WireCell;

using nfkdtree_static = NFKDVec::Tree<double, NFKDVec::IndexStatic>;
using nfkdtree_dynamic = NFKDVec::Tree<double>; // default

// these derived types do not depend on static/dynamic
using coordinates_type = nfkdtree_dynamic::coordinates_type;
using points_type = nfkdtree_dynamic::points_type;
using results_type = nfkdtree_dynamic::results_type;

using point_type = std::vector<double>;

static
bool same_coordinates(const coordinates_type& c1, const coordinates_type& c2)
{
    if (c1.size() != c2.size()) return false;

    const size_t num = c1.size();
    for (size_t ind=0; ind<num; ++ind) {
        if (c1.at(ind) != c2.at(ind)) return false;
    }
    return true;
}
static
bool same_points(const points_type& p1, const points_type& p2)
{
    if (p1.size() != p2.size()) return false;
    const size_t ndim = p1.size();
    for (size_t dim=0; dim<ndim; ++dim) {
        if (! same_coordinates(p1[dim], p2[dim])) return false;
    }
    return true;
}


static points_type make_points(size_t copies = 1)
{
    points_type ret(3), points = {
        {1.0,1.0,1.0},
        {2.0,1.0,3.0},
        {1.0,4.0,1.0}
    };
    
    while (copies) {
        for (size_t dim=0; dim<points.size(); ++dim) {
            for (size_t ind=0; ind<points[0].size(); ++ind) {
                ret[dim].push_back(copies*10 + points[dim][ind]);
            }
        }
        --copies;
    }
    return ret;
}
static points_type make_other_points()
{
    return {
        {2.0,1.0,1.0},
        {1.0,1.0,3.0},
        {1.0,4.0,5.0}
    };
}

static points_type random_points(size_t seed, size_t number = 3, size_t ndim = 3,
                                 double xmin=-1.0, double xmax=1.0)
{
    std::default_random_engine re{seed};
    std::uniform_real_distribution<double> dist(xmin, xmax);

    points_type points(ndim);
    for (size_t dim=0; dim<ndim; ++dim) {
        points[dim].resize(number);
        for (size_t ind=0; ind<number; ++ind) {
            points[dim][ind] = dist(re);
        }
    }
    return points;
}

template<typename nfkdtree_type>
void dump_kd_result(const nfkdtree_type& kd, const results_type& res)
{
    for (const auto& [ind,dist] : res) {
        auto pt = kd.point(ind);
        debug("query result: point=({},{},{}), dist={}, kd point calls={}",
              pt[0], pt[1], pt[2], dist, kd.point_calls());
    }
}

// by a fair roll of fishy dice
// $ for n in (seq 20)
//      printf "%d," (random) 
//  end
static std::vector<size_t> seeds = {
    31818,8536,32603,17366,502,20163,21205,19001,
    2242,30810,11876,4798,2771,4880,31270,29421,
    29210,14507,30280,30681
};

template<typename IndexType>
void appending()
{
    NFKDVec::Tree<double, IndexType> kd(3);
    REQUIRE(kd.npoints() == 0);
    for (size_t block=0; block<3; ++block) {
        points_type points = random_points(seeds[block],3);
        kd.append(points);
    }

    auto& kdpoints = kd.points();
    REQUIRE(kd.npoints() == 9);
    REQUIRE(kdpoints.size() == 3);
    REQUIRE(kdpoints[0].size() == 9);

    const point_type origin = {0.0,0.0,0.0};
    auto res = kd.radius(1000.0, origin);
    REQUIRE(res.size() == 9);

    for (const auto& [index, dist2] : res) {
        size_t major = kd.major_index(index);
        size_t minor = kd.minor_index(index);
        debug("index={} major={} minor={} dist2={}, pt=({},{},{})",
              index, major, minor, dist2,
              kdpoints[0][index], kdpoints[1][index], kdpoints[2][index]);
    }

    for (size_t block=3; block<6; ++block) {
        points_type points = random_points(seeds[block],3);
        kd.append(points);
    }

    REQUIRE(kd.npoints() == 18);
    REQUIRE(kdpoints.size() == 3);
    REQUIRE(kdpoints[0].size() == 18);

    std::unordered_set< size_t > unique_points;
    for (size_t ind=0; ind<18; ++ind) {
        std::size_t h = 0;
        boost::hash_combine(h, kdpoints[0][ind]);
        boost::hash_combine(h, kdpoints[1][ind]);
        boost::hash_combine(h, kdpoints[2][ind]);
        unique_points.insert(h);
    }
    REQUIRE(unique_points.size() == 18);
    

    auto res2 = kd.radius(1000.0, origin);
    REQUIRE(res2.size() == 18);

    for (size_t ind=0; ind<18; ++ind) {
        std::vector<double> pt = {kdpoints[0][ind],kdpoints[1][ind],kdpoints[2][ind]};

        auto res = kd.knn(1, pt);
        REQUIRE(res.size() == 1);
        CHECK(res[0].first == ind);
        CHECK(res[0].second == 0.0);

        auto res2 = kd.radius(0.00001, pt); // randoms could get lucky!
        bool found_it=false;
        for (const auto& [pind,dist] : res) {
            if (pind != ind) {
                continue;
            }
            CHECK(dist == 0);
            found_it = true;
        }
        CHECK(found_it);
    }
        

}

TEST_CASE("nfkdvec appending dynamic")
{
    appending<NFKDVec::IndexDynamic>();
}
TEST_CASE("nfkdvec appending static")
{
    appending<NFKDVec::IndexStatic>();
}

TEST_CASE("nfkdvec dynamic")
{
    points_type points = make_points();

    debug("nfkdvec: making k-d tree");
    nfkdtree_dynamic kd(points);

    REQUIRE(kd.npoints() == 3);
    debug("nfkdvec: {} point calls", kd.point_calls());

    CHECK(same_points(kd.points(), points));
    // CHECK_THROWS_AS(kd.points().at(3), std::out_of_range); // throws

    for (size_t count = 0; count < 3; ++count)
    {
        debug("nfkdvec: radius query #{}", count);
        point_type pt = {11.0, 12.0, 11.0};
        auto res = kd.radius(0.01, pt);
        CHECK(res.size() == 1);
        dump_kd_result(kd, res);
    }

    debug("nfkdvec: appending");
    kd.append(points);

    auto& kdpoints = kd.points();

    REQUIRE(kd.npoints() == 6);
    REQUIRE(kdpoints.size() == 3);
    REQUIRE(kdpoints[0].size() == 6);

    for (size_t kpt=0; kpt<6; ++kpt) {
        const size_t ipt = 0;   // points reused
        for (size_t dim=0; dim<3; ++dim) {
            const double kele = kdpoints.at(dim)[kpt];
            const double iele = points.at(dim)[ipt];
            debug("kpt={} dim={} {} =?= {}", kpt, dim, kele, iele);
            // CHECK(kele == iele);
        }
    }

    for (size_t N = 1; N <= 6; ++N) {
        debug("nfkdvec: knn={} query", N);
        point_type pt = {0,0,0};
        auto res = kd.knn(N, pt);
        CHECK(res.size() == N);
        dump_kd_result(kd, res);
    }

    debug("nfkdvec: appending more");
    kd.append(make_other_points());

    REQUIRE(kd.npoints() == 9);
    REQUIRE(kdpoints.size() == 3);
    REQUIRE(kdpoints[0].size() == 9);

    {
        point_type pt = {0.0,0.0,0.0};
        auto res = kd.knn(100, pt);
        REQUIRE(res.size() == 9);
    }
    {
        point_type pt = {0.0,0.0,0.0};
        auto res = kd.radius(1000.0, pt);
        REQUIRE(res.size() == 9);

        for (const auto& [index, dist2] : res) {
            size_t major = kd.major_index(index);
            size_t minor = kd.minor_index(index);
            debug("index={} major={} minor={} dist2={}, pt=({},{},{})",
                  index, major, minor, dist2,
                  kdpoints[0][index], kdpoints[1][index], kdpoints[2][index]);
        }
    }

}

TEST_CASE("nfkdvec static two")
{
    points_type points = make_points(2);

    debug("nfkdvec: making k-d tree");
    nfkdtree_static kd(points);

    REQUIRE(kd.npoints() == 6);
    debug("nfkdvec: {} point calls", kd.point_calls());

    CHECK(same_points(kd.points(), points));
    // CHECK_THROWS_AS(kd.points().at(6), std::out_of_range); // throws

    for (size_t count = 0; count < 3; ++count)
    {
        debug("nfkdvec: radius query #{}", count);
        point_type pt = {11.0, 12.0, 11.0};
        auto res = kd.radius(0.01, pt);
        CHECK(res.size() == 1);
        dump_kd_result(kd, res);
    }

    auto& kdpoints = kd.points();

    REQUIRE(kd.npoints() == 6);
    REQUIRE(kdpoints.size() == 3);
    REQUIRE(kdpoints[0].size() == 6);

    for (size_t kpt=0; kpt<6; ++kpt) {
        const size_t ipt = 0;   // points reused
        for (size_t dim=0; dim<3; ++dim) {
            const double kele = kdpoints.at(dim)[kpt];
            const double iele = points.at(dim)[ipt];
            debug("kpt={} dim={} {} =?= {}", kpt, dim, kele, iele);
            // CHECK(kele == iele);
        }
    }

    for (size_t N = 1; N <= 6; ++N) {
        debug("nfkdvec: knn={} query", N);
        point_type pt = {0,0,0};
        auto res = kd.knn(N, pt);
        CHECK(res.size() == N);
        dump_kd_result(kd, res);
    }
}

TEST_CASE("nfkdvec static")
{
    points_type points = make_points();
        
    debug("nfkdvec: making k-d tree");
    nfkdtree_static kd(points);

    CHECK(kd.npoints() == 3);
    debug("nfkdvec: {} point calls", kd.point_calls());

    CHECK(same_points(kd.points(), points));
    // CHECK_THROWS_AS(kd.points().at(3), std::out_of_range); // throws

    for (size_t count = 0; count < 3; ++count)
    {
        debug("nfkdvec: radius query #{}", count);
        point_type pt = {11.0, 12.0, 11.0};
        auto res = kd.radius(0.01, pt);
        CHECK(res.size() == 1);
        dump_kd_result(kd, res);
    }

    for (size_t N = 1; N <= 3; ++N) {
        debug("nfkdvec: knn={} query", N);
        point_type pt = {0,0,0};
        auto res = kd.knn(N, pt);
        CHECK(res.size() == N);
        dump_kd_result(kd, res);
    }

    {
        const auto& ckd = kd;
        point_type pt = {11.0, 12.0, 11.0};
        auto res = ckd.radius(0.01, pt);
        CHECK(res.size() == 1);
    }    

}
