#include "WireCellUtil/NFKD.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/ExecMon.h"
#include "WireCellUtil/String.h"

#include "WireCellUtil/PointCloudCoordinates.h"
#include "WireCellUtil/PointCloudDataset.h"

// #define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "WireCellUtil/doctest.h"

#include <vector>

using namespace WireCell;
using namespace WireCell::PointCloud;
using WireCell::String::format;


TEST_CASE("nfkd dataset atomic")
{
    Dataset ds({
            {"x", Array({1.0, 1.0, 1.0})},
            {"y", Array({2.0, 1.0, 3.0})},
            {"z", Array({1.0, 4.0, 1.0})},
            {"one", Array({1  ,2  ,3  })},
            {"two", Array({1.1,2.2,3.3})}});
    std::vector<std::string> names = {"x","y","z"};
    
    auto sel = ds.selection(names);
    using coords_type = coordinate_array<double>;
    using point_type = coords_type::point;
    coords_type coords(sel);
    std::vector<coords_type> vcoords = { coords };
    
    {
        const size_t npoints = coords.size();
        REQUIRE(npoints == 3);
        for (size_t ind=0; ind<npoints; ++ind) {
            const auto& pt = *(coords.begin() + ind);
            REQUIRE(pt.size() == 3);
            for (size_t dim=0; dim<3; ++dim) {
                SPDLOG_TRACE("nfkd dataset atomic:\tdirect[{},{}] = {}", ind, dim, pt.at(dim));
            }
        }
    }
    {
        using points_t = disjoint_range<coords_type>;
        points_t points(vcoords);
        SPDLOG_TRACE("nfkd dataset atomic: made disjoint points");
        const size_t npoints = points.size();
        REQUIRE(npoints == 3);
        for (size_t ind=0; ind<npoints; ++ind) {
            const point_type& pt = points.at(ind);
            REQUIRE(pt.size() == 3);
            for (size_t dim=0; dim<3; ++dim) {
                REQUIRE(pt.size() == 3);
                SPDLOG_TRACE("nfkd dataset atomic:\tdisjoint[{},{}] = {}", ind, dim, pt.at(dim));
            }
        }
    }
    {
        using points_t = disjoint_range<coords_type>;
        points_t points(vcoords);
        const size_t npoints = points.size();
        REQUIRE(npoints == 3);
        SPDLOG_TRACE("nfkd dataset atomic: made points");
        for (size_t ind=0; ind<npoints; ++ind) {
            SPDLOG_TRACE("nfkd dataset atomic: getting point {}", ind);
            const point_type& pt = points.at(ind); // ref
            SPDLOG_TRACE("nfkd dataset atomic: point {} size={}", ind, pt.size());
            REQUIRE(pt.size() == 3);
            for (size_t dim=0; dim<3; ++dim) {
                REQUIRE(pt.size() == 3);
                SPDLOG_TRACE("nfkd dataset atomic:\t[{},{}] = {}", ind, dim, pt.at(dim));
            }
        }
    }


    using kdtree_type = NFKD::Tree<coords_type, NFKD::IndexStatic>;
    kdtree_type kd(names.size(), vcoords);

    auto& points = kd.points();

    CHECK(points.size() == 3);
    SPDLOG_DEBUG("nfkd: {} point calls", kd.point_calls());

    {
        auto& point = points.at(0);
        CHECK(point.at(0) == 1.0);
        CHECK(point.at(1) == 2.0);
        CHECK(point.at(2) == 1.0);
    }
    {
        auto& point = points.at(1);        
        CHECK(point.at(0) == 1.0);
        CHECK(point.at(1) == 1.0);
        CHECK(point.at(2) == 4.0);
    }
    {
        auto& point = points.at(2);
        CHECK(point.at(0) == 1.0);
        CHECK(point.at(1) == 3.0);
        CHECK(point.at(2) == 1.0);
    }

    CHECK_THROWS_AS(points.at(3), std::out_of_range); // throws

    // Do three times and observe that calls to the point evaluation
    // method are repeated.
    for (size_t count = 0; count < 3; ++count)
    {
        SPDLOG_DEBUG("nfkd: radius query #{}", count);
        std::vector<double> pt = {1.0, 2.0, 1.0};
        auto res = kd.radius(0.01, pt);
        CHECK(res.size() == 1);
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_DEBUG
        for (const auto& [it,dist] : res) {
            SPDLOG_DEBUG("rad: #{} point=({},{},{}), dist={}, calls={}",
                         count, it->at(0), it->at(1), it->at(2), dist, kd.point_calls());
        }
#endif
    }
    for (size_t N = 1; N <= 3; ++N) {
        SPDLOG_DEBUG("nfkd: knn={} query", N);
        auto res = kd.knn<std::vector<double>>(N, {0,0,0});
        REQUIRE(res.size() == N);
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_DEBUG
        for (const auto& [it,dist] : res) {
            SPDLOG_DEBUG("knn: N={} point=({},{},{}), dist={}, calls={}",
                  N, it->at(0), it->at(1), it->at(2),
                  dist, kd.point_calls());
        }
#endif
    }

}


 

#include <random>
#include <algorithm>


// Smallest point cloud to consider
static const size_t start_size = 1024;

// How many times we double the point cloud size in tests.  Want to
// keep this somewhat small to allow the test to run quickly.  Feel
// free to increase it for a local build and test but keep it at ~3
// for commits.  see util/docs/nary-tree.org for some benchmarks.
static const size_t max_doubling = 3;

// A list of coordinate names and thus dimension of the k-d space.
static const std::vector<std::string> coord_names = {"x","y","z"};

Dataset make_dataset(size_t npts,
                     const std::vector<std::string>& names = coord_names)
{
    std::random_device rnd_device;
    // Specify the engine and distribution.
    std::mt19937 mersenne_engine {rnd_device()};  // Generates random integers
    std::uniform_real_distribution<double> dist {-1.0, 1.0};
    
    auto gen = [&dist, &mersenne_engine](){
                   return dist(mersenne_engine);
               };

    SPDLOG_DEBUG("make_dataset with {} points", npts);
    Dataset ds;

    for (const auto& name: names) {
        std::vector<double> arr(npts);
        SPDLOG_DEBUG("make_dataset make array {} with {} points", name, npts);
        generate(begin(arr), end(arr), gen);
        Array aa(arr);
        ds.add(name, std::move(aa));
        SPDLOG_DEBUG("make_dataset add array {} with {} points", name, npts);
    }
    SPDLOG_DEBUG("make_dataset with {} points returning", npts);
    return ds;
}


TEST_CASE("point cloud disjoint serialized construction")
{
    using coords_type = coordinate_array<double>;
    // using point_type = coords_type::point;

    std::vector<Dataset> ds_store;
    std::vector<Dataset::selection_t> se_store;
    std::vector<coords_type> co_store;

    const size_t nper = 74;
    const size_t ndses = 2;
    const size_t npts = nper*ndses;

    for (size_t count=0; count<ndses; ++count)  {
        SPDLOG_DEBUG("count={} make and emplace dataset", count);
        REQUIRE(count == ds_store.size());
        if (true) {
            ds_store.emplace_back(make_dataset(nper, coord_names));
        }
        else {
            Dataset ds_actual = make_dataset(nper, coord_names);
            ds_store.push_back(ds_actual);
        }
        Dataset& ds = ds_store.back();
        SPDLOG_DEBUG("count={} added ds with {}", count, ds.size_major());
        REQUIRE(nper == ds.size_major());
        REQUIRE(nper == ds_store[count].size_major());
    }
    for (size_t count=0; count<ndses; ++count)  {
        Dataset& ds = ds_store[count];
        SPDLOG_DEBUG("count={} add selection from ds with {}", count, ds.size_major());
        se_store.push_back(ds.selection(coord_names));
        auto& se = se_store.back();
        auto arr0 = se[0];
        SPDLOG_DEBUG("count={} add se with {} ({})", count, arr0->size_major(), se.size());
    }

    for (size_t count=0; count<ndses; ++count)  {
        auto& se = se_store[count];
        auto arr0 = se[0];
        SPDLOG_DEBUG("count={} add se with {} ({})", count, arr0->size_major(), se.size());
        co_store.emplace_back(se);
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_DEBUG
        auto& co = co_store.back();
        SPDLOG_DEBUG("count={} add co with {}", count, co.size());
#endif
    }

    for (size_t count=0; count<ndses; ++count) {
        SPDLOG_DEBUG("count={} check: ds={} se={}/{} co={}",
              count, 
              ds_store[count].size_major(),
              se_store[count].size(),
              se_store[count][0]->size_major(),
              co_store[count].size());

        REQUIRE(nper == ds_store[count].size_major());
        REQUIRE(nper == se_store[count][0]->size_major());
        REQUIRE(nper == co_store[count].size());
    }
    using kdtree_type = NFKD::Tree<coords_type, NFKD::IndexDynamic>;
    kdtree_type kd(coord_names.size());
    for (size_t count=0; count<ndses; ++count) {
        SPDLOG_DEBUG("count={} kd has {} adding: ds={} se={}/{} co={}",
              count, kd.kdtree_get_point_count(), 
              ds_store[count].size_major(),
              se_store[count].size(),
              se_store[count][0]->size_major(),
              co_store[count].size());

        kd.append(co_store[count]);
    }

    REQUIRE(kd.points().size() == npts);
    
}

TEST_CASE("point cloud disjoint staged construction")
{
    using coords_type = coordinate_array<double>;
    // using point_type = coords_type::point;

    const size_t nper = 74;
    const size_t ndses = 2;
    const size_t npts = nper*ndses;

    std::map<size_t, Dataset> ds_store;
    using selection_t = Dataset::selection_t;
    //std::vector<std::unique_ptr<selection_t>> se_store;
    // std::vector<selection_t> se_store;
    std::map<size_t, selection_t> se_store;
    std::vector<coords_type> co_store;

    for (size_t count=0; count<ndses; ++count)  {
        {
            SPDLOG_DEBUG("count={} make and emplace dataset", count);
            REQUIRE(count == ds_store.size());
            {
                Dataset ds_local = make_dataset(nper, coord_names);
                SPDLOG_DEBUG("count={} made dataset ({}) with {}",
                      count, (void*)&ds_local, ds_local.size_major());
                ds_store[count] = ds_local;
                SPDLOG_DEBUG("count={} placed dataset ({}) with {}",
                      count, (void*)&ds_store[count], ds_store[count].size_major());
            }
            Dataset& ds = ds_store[count];
            SPDLOG_DEBUG("count={} added ds with {}", count, ds.size_major());
            REQUIRE(nper == ds.size_major());
            REQUIRE(ds_store.size() == count+1);
            REQUIRE(nper == ds_store[count].size_major());
        }
        REQUIRE(ds_store.size() == count+1);
        REQUIRE(nper == ds_store[count].size_major());
        {
            Dataset& ds = ds_store[count];
            SPDLOG_DEBUG("count={} add selection from ds with {}", count, ds.size_major());
            // se_store.push_back(std::make_unique<selection_t>(ds.selection(coord_names)));
            se_store[count] = ds.selection(coord_names);
            REQUIRE(se_store.size() == count + 1);
            REQUIRE(se_store[count].size() > 0);
            const auto* septr = &se_store[count];
            auto arr0 = septr->at(0);
            SPDLOG_DEBUG("count={} add se with {} ({}) arr0@ {} se@ {}",
                  count, arr0->size_major(), septr->size(),
                  (void*)arr0.get(),
                  (void*)septr);
            for (size_t ind = 0; ind<=count; ++ind) {
                // const auto& selptr = se_store[ind].get();
                const auto& selptr = &se_store[ind];
                auto arr = selptr->at(0);
                Dataset& dsi = ds_store[ind];
                auto new_se = dsi.selection(coord_names);
                auto arr2 = new_se[0];
                SPDLOG_DEBUG("count={} ind={} add se with {} ({}) arr0@ {} arr2@ {} se@ {} ds@ {}",
                      count, ind, arr->size_major(), selptr->size(),
                      (void*)arr.get(), (void*)arr2.get(),
                      (void*)selptr,
                      (void*)(&dsi));
                REQUIRE(nper == arr->size_major());
            }
        }
        REQUIRE(ds_store.size() == count+1);
        REQUIRE(nper == ds_store[count].size_major());
        REQUIRE(se_store.size() == count + 1);
        REQUIRE(se_store[count].size() > 0);
        REQUIRE(nper == se_store[count].at(0)->size_major());
        {
            selection_t* septr = &se_store[count];
            auto arr0 = septr->at(0);
            SPDLOG_DEBUG("count={} add se with {} ({})", count, arr0->size_major(), septr->size());
            co_store.emplace_back(*septr);
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_DEBUG
            auto& co = co_store[count];
            SPDLOG_DEBUG("count={} add co with {}", count, co.size());
#endif
        }
        REQUIRE(nper == ds_store[count].size_major());
        REQUIRE(se_store[count].size() > 0);
        REQUIRE(nper == se_store[count].at(0)->size_major());
        REQUIRE(nper == co_store[count].size());
    }

    {
        REQUIRE(nper == se_store[ndses-1].at(0)->size_major());
        REQUIRE(nper == co_store[ndses-1].size());
    }

    for (size_t count=0; count<ndses; ++count) {
        SPDLOG_DEBUG("count={} check: ds={} se={}/{} co={}",
              count, 
              ds_store[count].size_major(),
              se_store[count].size(),
              se_store[count].at(0)->size_major(),
              co_store[count].size());

        REQUIRE(nper == ds_store[count].size_major());
        REQUIRE(se_store[count].size() > 0);
        REQUIRE(nper == se_store[count].at(0)->size_major());
        REQUIRE(nper == co_store[count].size());
    }
    using kdtree_type = NFKD::Tree<coords_type, NFKD::IndexDynamic>;
    kdtree_type kd(coord_names.size());
    for (size_t count=0; count<ndses; ++count) {
        SPDLOG_DEBUG("count={} kd has {} adding: ds={} se={}/{} co={}",
              count, kd.kdtree_get_point_count(), 
              ds_store[count].size_major(),
              se_store[count].size(),
              se_store[count].at(0)->size_major(),
              co_store[count].size());

        kd.append(co_store[count]);
    }

    REQUIRE(kd.points().size() == npts);
    
}

TEST_CASE("point cloud disjoint iteration")
{
    ExecMon em("point cloud disjoint iteration");

    // using point_type = std::vector<double>;

    for (size_t npts = start_size; npts <= (1<<max_doubling)*start_size; npts *= 2) {

        for (size_t nper = 16; nper <= npts; nper *= 8) {
            const size_t ndses = npts/nper;

            SPDLOG_DEBUG("\n{}", em(format("%d points, %d datasets: start test", npts, ndses)));

            std::vector<Dataset> ds_store;
            for (size_t count=0; count<ndses; ++count) {
                ds_store.emplace_back(make_dataset(nper, coord_names));
            }

            SPDLOG_DEBUG("\n{}", em(format("%d points, %d datasets: built datasets", npts, ndses)));
        }
    }
    SPDLOG_DEBUG("Time and memory usage:\n{}", em.summary());
}


template<typename IndexType>
void do_monolithic(const std::string& index_name)
{

    ExecMon em("nfkd dataset monolithic " + index_name);

    using coords_type = coordinate_array<double>;
    // using point_type = coords_type::point;
    using kdtree_type = NFKD::Tree<coords_type, IndexType>;

    std::vector<double> origin = {0,0,0};

    for (size_t npts = start_size; npts <= (1<<max_doubling)*start_size; npts *= 2) {

        em(format("%d points: start test", npts));

        auto ds = make_dataset(npts, coord_names);

        em(format("%d points: built dataset", npts));

        auto sel = ds.selection(coord_names);
        coords_type coords(sel);
        std::vector<coords_type> vcoords = {coords};
        kdtree_type kd(coord_names.size(), vcoords);
        CHECK(kd.points().size() == npts);

        em(format("%d points: built k-d tree", npts));

        for (size_t N = 1; N <= npts; N *= 2) {
            auto knn = kd.knn(N, origin);
            REQUIRE(knn.size() == N);

            em(format("%d points: knn query returns %d (expected %d)", npts, knn.size(), N));
        }

        for (double radius = 0.001; radius <= 1; radius *= 10) {
            auto rad = kd.radius(radius, origin);

            em(format("%d points: rad query returns %d (radius %f)", npts, rad.size(), radius));
        }
    }
    SPDLOG_DEBUG("Time and memory usage:\n{}", em.summary());
}


TEST_CASE("nfkd dataset performance monolithic static")
{
    do_monolithic<NFKD::IndexStatic>("static");
}
TEST_CASE("nfkd dataset performance monolithic dynamic")
{
    do_monolithic<NFKD::IndexDynamic>("dynamic");
}

TEST_CASE("nfkd dataset performance disjoint")
{
    ExecMon em("nfkd dataset disjoint");

    using coords_type = coordinate_array<double>;
    // using point_type = coords_type::point;
    using kdtree_type = NFKD::Tree<coords_type, NFKD::IndexDynamic>;
    const std::vector<double> origin = {0,0,0};

    for (size_t npts = start_size; npts <= (1<<max_doubling)*start_size; npts *= 2) {

        for (size_t nper = 16; nper <= npts; nper *= 8) {
            const size_t ndses = npts/nper;

            em(format("%d points, %d datasets: start test", npts, ndses));

            std::vector<Dataset> ds_store;
            for (size_t count=0; count<ndses; ++count) {
                ds_store.emplace_back(make_dataset(nper, coord_names));
            }

            em(format("%d points, %d datasets: built datasets", npts, ndses));

            std::vector<Dataset::selection_t> se_store;
            for (size_t count=0; count<ndses; ++count) {
                se_store.push_back(ds_store[count].selection(coord_names));
            }

            em(format("%d points, %d datasets: built selections", npts, ndses));

            std::vector<coords_type> co_store;
            for (size_t count=0; count<ndses; ++count) {
                Dataset::selection_t& sel = se_store[count];
                co_store.emplace_back(sel);
            }

            em(format("%d points, %d datasets: built coords", npts, ndses));

            kdtree_type kd(coord_names.size());
            for (size_t count=0; count<ndses; ++count) {
                auto& co = co_store.back();
                kd.append(co);
            }
            REQUIRE(kd.points().size() == npts);

            em(format("%d points, %d datasets: built k-d tree", npts, ndses));

            for (size_t N = 1; N <= npts; N *= 2) {
                auto knn = kd.knn(N, origin);
                CHECK(knn.size() == N);

                em(format("%d points, %d datasets: knn query returns %d (expected %d)",
                          npts, ndses, knn.size(), N));
            }

            for (double radius = 0.001; radius <= 1; radius *= 10) {
                auto rad = kd.radius(radius, origin);
                CHECK(rad.size() <= npts);

                em(format("%d points, %d datasets: rad query returns %d (radius %f)",
                          npts, ndses, rad.size(), radius));
            }
        }

    }
    spdlog::debug("Time and memory usage:\n{}", em.summary());
}
