#include "WireCellUtil/NFKDVec.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Logging.h"

// #define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "WireCellUtil/doctest.h"

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
            ret[dim].insert(ret[dim].end(), points[dim].begin(), points[dim].end());
        }
        --copies;
    }
    return ret;
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
        point_type pt = {1.0, 2.0, 1.0};
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
        point_type pt = {1.0, 2.0, 1.0};
        auto res = kd.radius(0.01, pt);
        CHECK(res.size() == 2);
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
        point_type pt = {1.0, 2.0, 1.0};
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
        point_type pt = {1.0, 2.0, 1.0};
        auto res = ckd.radius(0.01, pt);
        CHECK(res.size() == 1);
    }    

}
