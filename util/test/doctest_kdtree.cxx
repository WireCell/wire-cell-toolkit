#include "WireCellUtil/KDTree.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/TimeKeeper.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

#include <vector>
#include <string>
#include <iostream>


using namespace WireCell::PointCloud;
using namespace WireCell::KDTree;
using WireCell::TimeKeeper;
using namespace spdlog;

TEST_CASE("kdtree static")
{
    Dataset d({
        {"one", Array({1.0,1.0,3.0})},
        {"two", Array({1.1,2.2,3.3})},
    });

    auto one = d.get("one")->elements<double>();
    auto two = d.get("two")->elements<double>();

    auto kdq = query<double>(d, {"one", "two"});
    {
        auto a = kdq->knn(2, {2.0, 2.0});
        const size_t nfound = a.index.size();
        debug("{} found of 2 nearest neighbors", nfound);
        REQUIRE(nfound == 2);
        for (size_t ifound=0; ifound<nfound; ++ifound) {
            size_t ind = a.index[ifound];
            double dist = a.distance[ifound];
            debug("{}: [{}] p=({},{}) d={}", ifound, ind, one[ind], two[ind], dist);
        }
    }
    {
        // The extreme point is at linear r = 4.4598206.
        const double linear_radius = 4.46;
        // we use L2 metric!
        const double metric_radius = linear_radius*linear_radius; 

        auto a = kdq->radius(metric_radius, {0.0, 0.0});
        const size_t nfound = a.index.size();
        debug("{} found in radius", nfound);
        REQUIRE(nfound == 3);
        for (size_t ifound=0; ifound<nfound; ++ifound) {
            size_t ind = a.index[ifound];
            double dist = a.distance[ifound];
            debug("{}: [{}] p=({},{}) d={}", ifound, ind, one[ind], two[ind], dist);
        }
    }
    {
        // Catch just the first point
        const double linear_radius = 1.49;
        // we use L2 metric!
        const double metric_radius = linear_radius*linear_radius; 

        auto a = kdq->radius(metric_radius, {0.0, 0.0});
        const size_t nfound = a.index.size();
        debug("{} found in radius", nfound);
        REQUIRE(nfound == 1);
        for (size_t ifound=0; ifound<nfound; ++ifound) {
            size_t ind = a.index[ifound];
            double dist = a.distance[ifound];
            debug("{}: [{}] p=({},{}) d={}", ifound, ind, one[ind], two[ind], dist);
            REQUIRE(ind == 0);
        }
    }
}

TEST_CASE("kdtree dynamic")
{
    Dataset d({
        {"one", Array({1.0,1.0,3.0})},
        {"two", Array({1.1,2.2,3.3})},
    });

    // This must have registered update callback.
    auto kdq = query<double>(d, {"one", "two"}, true);
    REQUIRE(kdq);

    auto arrs = d.selection({"one", "two"});

    const Array& one = *arrs[0];
    const Array& two = *arrs[1];

    REQUIRE(one.num_elements() == 3);
    REQUIRE(two.num_elements() == 3);

    Dataset tail({
            {"one", Array({1.1})},
            {"two", Array({1.0})}});
    d.append(tail);

    REQUIRE(one.num_elements() == 4);
    REQUIRE(two.num_elements() == 4);

    const double linear_radius = 1.49;
    const double metric_radius = linear_radius*linear_radius; 

    auto ones = one.elements<double>();
    auto twos = two.elements<double>();

    auto a = kdq->radius(metric_radius, {0.0, 0.0});
    const size_t nfound = a.index.size();
    debug("{} found in radius", nfound);
    REQUIRE(nfound == 2);
    for (size_t ifound=0; ifound<nfound; ++ifound) {
        size_t ind = a.index[ifound];
        double dist = a.distance[ifound];
        debug("{}: [{}] p=({},{}) d={}", ifound, ind, ones[ind], twos[ind], dist);
    }

}

TEST_CASE("kdtree multi")
{
    Dataset d({
        {"one", Array({1.0,1.0,3.0})},
        {"two", Array({1.1,2.2,3.3})},
    });

    MultiQuery mq(d);

    // This must have registered update callback.
    Dataset::name_list_t onetwo = {"one", "two"};

    const bool dynamic = true;
    auto kdq = mq.get<double>(onetwo, dynamic);
    assert(kdq);
    REQUIRE(kdq->dynamic() == dynamic);
    REQUIRE(kdq->metric() == Metric::l2simple);

    auto kdq2 = mq.get<double>(onetwo, dynamic);
    REQUIRE(kdq2 == kdq);
    REQUIRE(kdq2);
    REQUIRE(kdq2->dynamic() == dynamic);
    REQUIRE(kdq2->metric() == Metric::l2simple);
}

#include <random>
void test_speed(size_t num=1000, size_t nlu = 1000, size_t kay=10,
                bool shared=true, bool dynamic=true);
void test_speed(size_t num, size_t nlu, size_t kay,
                bool shared, bool dynamic)
{
    std::stringstream label;
    label << "test KDTree: num=" << num
          << " nlu=" << nlu
          << " kay=" << kay
          << " shared=" << shared
          << " dynamic=" << dynamic;
    debug(label.str());

    TimeKeeper tk(label.str());

    const double xmax=1000;
    const double xmin=-xmax;
    std::random_device rd;
    std::default_random_engine re(rd());
    std::uniform_real_distribution<double> dist(xmin, xmax);

    std::vector<double> v1(num, 0), v2(num, 0), v3(num, 0);
    for (size_t ind=0; ind<num; ++ind) {
        v1[ind] = dist(re);
        v2[ind] = dist(re);
        v3[ind] = dist(re);
    }

    {
        Array junk(v1.data(), {num}, shared);
        auto aptr = &junk;
        REQUIRE(aptr);
        REQUIRE(aptr->size_major() == num);
        auto dat = aptr->bytes();
        REQUIRE(dat.size() == num * sizeof(double));
        REQUIRE(dat.data());
    }

    tk("randoms made");

    Dataset d({
        {"x", Array(v1.data(), {num}, shared)},
        {"y", Array(v2.data(), {num}, shared)},
        {"z", Array(v3.data(), {num}, shared)},
    });
    tk("dataset made");

    Dataset::name_list_t names = {"x", "y", "z"};
    for (const auto& key : names) {
        auto aptr = d.get(key);
        REQUIRE(aptr);
        REQUIRE(aptr->size_major() == num);
        auto dat = aptr->bytes();
        REQUIRE(dat.size() == num * sizeof(double));
        REQUIRE(dat.data());
    }

    auto sel = d.selection(names);
    REQUIRE(! sel.empty());
    for (auto aptr : sel) {
        REQUIRE(aptr);
        REQUIRE(aptr->size_major() == num);
        auto dat = aptr->bytes();
        REQUIRE(dat.size() == num * sizeof(double));
        REQUIRE(dat.data());
    }

    debug("making k-d tree query");
    auto kdq = query<double>(d, names, dynamic);
    debug("made k-d tree query");

    tk("kdtree made");

    for (size_t ind=0; ind<nlu; ++ind) {
        std::vector<double> qp = {dist(re), dist(re), dist(re)};
        auto res = kdq->knn(kay, qp);
        if (!ind) { tk("first query"); }            
    }
    tk("queries made");

    debug("summary:\n{}", tk.summary());
    
}

TEST_CASE("kdtree speed 100000")
{
    test_speed(10000, 1000, 10, false);
}    
