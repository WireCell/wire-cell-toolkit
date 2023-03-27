#include "WireCellUtil/KDTree.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/TimeKeeper.h"

#include <vector>
#include <string>
#include <iostream>


using namespace WireCell::PointCloud;
using namespace WireCell::KDTree;
using WireCell::TimeKeeper;

void test_static()
{
    Dataset::store_t s = {
        {"one", Array({1.0,1.0,3.0})},
        {"two", Array({1.1,2.2,3.3})},
    };
    Dataset d(s);

    auto one = s["one"].elements<double>();
    auto two = s["two"].elements<double>();

    auto kdq = query<double>(d, {"one", "two"});
    {
        auto a = kdq->knn(2, {2.0, 2.0});
        const size_t nfound = a.index.size();
        std::cerr << nfound << " found of 2 nearest neighbors\n";
        Assert(nfound == 2);
        for (size_t ifound=0; ifound<nfound; ++ifound) {
            size_t ind = a.index[ifound];
            double dist = a.distance[ifound];
            std::cerr << ifound << ": [" <<ind<< "] p=("<<one[ind]<<","<<two[ind]<<") d=" << dist << "\n";
        }
    }
    {
        // The extreme point is at linear r = 4.4598206.
        const double linear_radius = 4.46;
        // we use L2 metric!
        const double metric_radius = linear_radius*linear_radius; 

        auto a = kdq->radius(metric_radius, {0.0, 0.0});
        const size_t nfound = a.index.size();
        std::cerr << nfound << " found in radius\n";
        Assert(nfound == 3);
        for (size_t ifound=0; ifound<nfound; ++ifound) {
            size_t ind = a.index[ifound];
            double dist = a.distance[ifound];
            std::cerr << ifound << ": [" <<ind<< "] p=("<<one[ind]<<","<<two[ind]<<") d=" << dist << "\n";
        }
    }
    {
        // Catch just the first point
        const double linear_radius = 1.49;
        // we use L2 metric!
        const double metric_radius = linear_radius*linear_radius; 

        auto a = kdq->radius(metric_radius, {0.0, 0.0});
        const size_t nfound = a.index.size();
        std::cerr << nfound << " found in radius\n";
        Assert(nfound == 1);
        for (size_t ifound=0; ifound<nfound; ++ifound) {
            size_t ind = a.index[ifound];
            double dist = a.distance[ifound];
            std::cerr << ifound << ": [" <<ind<< "] p=("<<one[ind]<<","<<two[ind]<<") d=" << dist << "\n";
            Assert(ind == 0);
        }
    }
}

void test_dynamic()
{
    Dataset::store_t s = {
        {"one", Array({1.0,1.0,3.0})},
        {"two", Array({1.1,2.2,3.3})},
    };
    Dataset d(s);

    // This must have registered update callback.
    auto kdq = query<double>(d, {"one", "two"}, true);

    auto arrs = d.selection({"one", "two"});

    const Array& one = arrs[0];
    const Array& two = arrs[1];

    Assert(one.num_elements() == 3);
    Assert(two.num_elements() == 3);

    Dataset tail({
            {"one", Array({1.1})},
            {"two", Array({1.0})}});
    d.append(tail);

    Assert(one.num_elements() == 4);
    Assert(two.num_elements() == 4);

    const double linear_radius = 1.49;
    const double metric_radius = linear_radius*linear_radius; 

    auto ones = one.elements<double>();
    auto twos = two.elements<double>();

    auto a = kdq->radius(metric_radius, {0.0, 0.0});
    const size_t nfound = a.index.size();
    std::cerr << nfound << " found in radius\n";
    Assert(nfound == 2);
    for (size_t ifound=0; ifound<nfound; ++ifound) {
        size_t ind = a.index[ifound];
        double dist = a.distance[ifound];
        std::cerr << ifound << ": [" <<ind<< "] p=("<<ones[ind]<<","<<twos[ind]<<") d=" << dist << "\n";
    }

}

void test_multi()
{
    Dataset::store_t s = {
        {"one", Array({1.0,1.0,3.0})},
        {"two", Array({1.1,2.2,3.3})},
    };
    Dataset d(s);

    MultiQuery mq(d);

    // This must have registered update callback.
    name_list_t onetwo = {"one", "two"};

    const bool dynamic = true;
    auto kdq = mq.get<double>(onetwo, dynamic);
    assert(kdq);
    Assert(kdq->dynamic() == dynamic);
    Assert(kdq->metric() == Metric::l2simple);

    auto kdq2 = mq.get<double>(onetwo, dynamic);
    Assert(kdq2 == kdq);
    Assert(kdq2);
    Assert(kdq2->dynamic() == dynamic);
    Assert(kdq2->metric() == Metric::l2simple);
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

    TimeKeeper tk(label.str());

    const double xmax=1000;
    const double xmin=-xmax;
    std::random_device rd;
    std::default_random_engine re(rd());
    std::uniform_real_distribution<double> dist(xmin, xmax);

    std::vector<double> v1(num, 0), v2(num, 0), v3(num, 0);
    for (size_t ind=0; ind<num; ++ind) {
        v1[num] = dist(re);
        v2[num] = dist(re);
        v3[num] = dist(re);
    }

    tk("randoms made");

    Dataset::store_t s = {
        {"x", Array(v1.data(), {num}, shared)},
        {"y", Array(v2.data(), {num}, shared)},
        {"z", Array(v3.data(), {num}, shared)},
    };
    tk("arrays made");

    Dataset d(s);
    tk("dataset made");

    auto kdq = query<double>(d, {"x", "y", "z"}, dynamic);
    tk("kdtree made");
    
    for (size_t ind=0; ind<nlu; ++ind) {
        std::vector<double> qp = {dist(re), dist(re), dist(re)};
        auto res = kdq->knn(kay, qp);
        if (!ind) { tk("first query"); }            
    }
    tk("queries made");

    std::cerr << tk.summary() << "\n";
    
}

int main() {
    test_static();
    test_dynamic();
    test_multi();
    // test_speed(100);
    // test_speed(1000);
    // test_speed(10000);
    test_speed(100000);
    return 0;
}
