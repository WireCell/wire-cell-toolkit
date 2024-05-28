#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Exceptions.h"

#include <vector>
#include <iostream>
#include <memory>
#include <numeric>              // iota

using namespace WireCell;
using namespace WireCell::PointCloud;
using spdlog::debug;

Dataset make_one()
{
    Dataset ds;
    std::vector<int> v = {0, 1, 2};
    Array arr(v);
    REQUIRE(arr.is_type<int>());
    ds.add("arr",arr);
    return ds;
}


TEST_CASE("point cloud ownership")
{
    Dataset ds({
            {"x", Array({1.0, 1.0, 1.0})},
            {"y", Array({2.0, 1.0, 3.0})},
            {"z", Array({1.0, 4.0, 1.0})},
            {"one", Array({1  ,2  ,3  })},
            {"two", Array({1.1,2.2,3.3})}});
    Dataset::selection_t sel = ds.selection({"x","y","z"});

    Dataset ds2(std::move(ds));
    REQUIRE(ds.size() == 0);
    REQUIRE(ds2.size() == 5);
    REQUIRE(ds2.get("x")->size_major() == 3);


    Dataset::selection_t sel2 = ds2.selection({"x","y","z"});    

    CHECK(sel[0] == sel2[0]);

    {
        Dataset ds3 = make_one();
        REQUIRE(ds3.has("arr"));
        auto arr = ds3.get("arr")->elements<int>();
        REQUIRE(arr.size() == 3);
    }
}

template<typename ElementType=double>
void check_a_slice(const Array& a)
{
    auto bytes = a.bytes();

    debug("a slice size_major={} ndims={} shape=({}) nbytes={} ele_size={}({})",
          a.size_major(), a.shape().size(), a.shape()[0], bytes.size(), a.element_size(), sizeof(ElementType));

    REQUIRE(a.element_size() == sizeof(ElementType));

    CHECK(a.size_major() == 1);
    auto sa = a.elements<ElementType>();
    CHECK(sa.size() == 1);
    CHECK(sa[0] == 2.0);

    auto ma = a.indexed<ElementType, 1>();
    CHECK(ma.num_dimensions() == 1);
    CHECK(ma.num_elements() == 1);
}

template<typename ElementType=int>
void check_b_slice(const Array& b)
{
    REQUIRE(b.element_size() == sizeof(ElementType));

    CHECK(b.size_major() == 1);
    auto bshape = b.shape();
    CHECK(bshape.size() == 3);
    CHECK(bshape[0] == 1);
    CHECK(bshape[1] == 2);
    CHECK(bshape[2] == 4);

    // flattened
    auto ba = b.elements<ElementType>();
    CHECK(ba.size() == 1*2*4);
    CHECK(ba[0] == 8);
    CHECK(ba[1] == 9);
    CHECK(ba[2] == 10);
    CHECK(ba[3] == 11);
    CHECK(ba[4] == 12);
    CHECK(ba[5] == 13);
    CHECK(ba[6] == 14);
    CHECK(ba[7] == 15);
    
    // indexed
    auto ma = b.indexed<ElementType, 3>();
    CHECK(ma.num_dimensions() == 3);
    CHECK(ma.num_elements() == 1*2*4);

    auto mshape = ma.shape();
    CHECK(mshape[0] == 1);
    CHECK(mshape[1] == 2);
    CHECK(mshape[2] == 4);
    CHECK(ma[0][0][0] == 8);
    CHECK(ma[0][0][1] == 9);
    CHECK(ma[0][0][2] == 10);
    CHECK(ma[0][0][3] == 11);
    CHECK(ma[0][1][0] == 12);
    CHECK(ma[0][1][1] == 13);
    CHECK(ma[0][1][2] == 14);
    CHECK(ma[0][1][3] == 15);
}

TEST_CASE("point cloud slice")
{
    Array aa({1.0, 2.0, 3.0});
    auto aas = aa.slice(1,1);
    check_a_slice(aas);
    CHECK(aas.dtype() == aa.dtype());
    

    // major
    // axis
    // 0   [[[ 0  1  2  3] [ 4  5  6  7]]
    // 1    [[ 8  9 10 11] [12 13 14 15]]
    // 2    [[16 17 18 19] [20 21 22 23]]]
    //
    // Will slice out m.a. index 1 to get:
    //
    // [[[ 8  9 10 11] [12 13 14 15]]]
    std::vector<size_t> shape = {3,2,4};
    std::vector<int> counts(24);
    std::iota(counts.begin(), counts.end(), 0);
    Array bb(counts, shape, false);
    auto bbs = bb.slice(1,1);
    check_b_slice(bbs);
    CHECK(bbs.dtype() == bb.dtype());

    Dataset ds({
            {"a", aa},
            {"b", bb}});

    Dataset slc = ds.slice(1,1);
    CHECK(slc.size_major() == 1);
    CHECK(slc.has("a"));
    CHECK(slc.has("b"));

    auto a = slc.get("a");
    REQUIRE(a);
    check_a_slice(*a);

    auto b = slc.get("b");
    REQUIRE(b);
    check_b_slice(*b);
}

template<typename ElementType>
void assure_equal(const Array& a, const std::vector<ElementType>& v)
{
    debug("dtype={} num_elements={} size={}",
          a.dtype(), a.num_elements(), v.size());
    CHECK(a.is_type<ElementType>());
    CHECK(v.size() == a.num_elements()); // these two are 
    CHECK(v.size() == a.size_major());   // degenerate for
    CHECK(a.shape().size() == 1);        // a 1D array.

    auto ele = a.elements<ElementType>();

    for (size_t ind=0; ind<v.size(); ++ind) {
        debug("{}: {} =?= {}", ind, v[ind], ele[ind]);
        CHECK(v[ind] == ele[ind]);
    }
}

template<typename ElementType>
PointCloud::Array make_array(const ElementType* data, const Array::shape_t& shape, bool share)
{
    PointCloud::Array arr(data, shape, share);
    return arr;
}
template<typename ElementType>
auto make_arrays(const ElementType* data, const Array::shape_t& shape, bool share)
{
    std::map<std::string, PointCloud::Array> ret;
    PointCloud::Array arr(data, shape, share);
    ret["a"] = arr;
    return ret;
}

template<typename ElementType=int>
void test_array_noshare()
{
    std::vector<ElementType> v{1,2,3};

    {
        Array a = make_array<ElementType>(v.data(), {v.size()}, true);
        assure_equal(a, v);
    }
    {
        Array a = make_array<ElementType>(v.data(), {v.size()}, false);
        assure_equal(a, v);
    }
    {
        auto d = make_arrays<ElementType>(v.data(), {v.size()}, true);
        const Array& a = d["a"];
        assure_equal(a, v);
        Array aa = d["a"];
        assure_equal(aa, v);
    }
    {
        auto d = make_arrays<ElementType>(v.data(), {v.size()}, false);
        const Array& a = d["a"];
        assure_equal(a, v);
        Array aa = d["a"];
        assure_equal(aa, v);
    }


    // cross: (mutable, const) x (shared, copy)
    Array sa(v.data(), {v.size()}, true);
    Array ns(v.data(), {v.size()}, false);

    assure_equal(sa, v);
    assure_equal(ns, v);

    std::vector<ElementType> w(v.begin(), v.end());

    sa.assure_mutable();

    v[0] = v[1] = v[2] = 0;
    v.clear();

    assure_equal(sa, w);
    assure_equal(ns, w);
}

TEST_CASE("point cloud array 1d")
{

    std::vector<int> v{1,2,3};
    const std::vector<int>& w = v;

    // cross: (mutable, const) x (shared, copy)
    Array ms(v.data(), {v.size()}, true);
    Array mc(v.data(), {v.size()}, false);
    Array cs(w.data(), {w.size()}, true);
    Array cc(w.data(), {w.size()}, false);

    assure_equal(ms, v);
    assure_equal(mc, v);
    assure_equal(cs, v);
    assure_equal(cc, v);

    {
        auto flat_span = ms.elements<int>();
        CHECK(flat_span.size() == v.size());
        ms.elements<float>(); // same size okay (for now)
        bool caught = true;
        try {
            ms.elements<double>();
        }
        catch (const ValueError& err) {
            caught = true;
        }
        CHECK(caught);
    }

    v.push_back(4);
    ms.assign(v.data(), {v.size()}, true);
    mc.assign(v.data(), {v.size()}, false);
    cs.assign(w.data(), {w.size()}, true);
    cc.assign(w.data(), {w.size()}, false);

    CHECK(ms.size_major() == 4);
    CHECK(mc.size_major() == 4);
    CHECK(cs.size_major() == 4);
    CHECK(cc.size_major() == 4);

    assure_equal(ms, v);
    assure_equal(mc, v);
    assure_equal(cs, v);
    assure_equal(cc, v);

    v[0] = 42;

    CHECK(ms.element<int>(0) == 42);
    CHECK(mc.element<int>(0) == 1);
    CHECK(cs.element<int>(0) == 42);
    CHECK(cc.element<int>(0) == 1);

    std::vector<double> d{1., 2., 3.};
    const std::vector<double>& p=d;

    ms.assign(d.data(), {d.size()}, true);
    mc.assign(d.data(), {d.size()}, false);
    cs.assign(p.data(), {p.size()}, true);
    cc.assign(p.data(), {p.size()}, false);

    assure_equal(ms, d);
    assure_equal(mc, d);
    assure_equal(cs, d);
    assure_equal(cc, d);

    d.push_back(4);
    ms.assign(d.data(), {d.size()}, true);
    mc.assign(d.data(), {d.size()}, false);
    cs.assign(p.data(), {p.size()}, true);
    cc.assign(p.data(), {p.size()}, false);

    assure_equal(ms, d);
    assure_equal(mc, d);
    assure_equal(cs, d);
    assure_equal(cc, d);

    d[0] = 42;

    CHECK(ms.element<double>(0) == 42);
    CHECK(mc.element<double>(0) == 1);
    CHECK(cs.element<double>(0) == 42);
    CHECK(cc.element<double>(0) == 1);

    auto msarr = ms.indexed<double, 1>();
    CHECK(msarr[0] == 42);

    auto msele = ms.elements<double>();
    CHECK(msele[0] == 42);

    std::vector<double> d2{4.,5.,6.};
    {
        const size_t before = ms.size_major();
        ms.append(d2);
        const size_t after = ms.size_major();
        CHECK( before + d2.size() == after);
    }
    CHECK(ms.element<double>(3) == 4);
    {
        CHECK(ms.shape().size() == 1);
        CHECK(ms.shape()[0] == 7);
    }

    ms.clear();
    mc.clear();
    cs.clear();
    cc.clear();

    debug("ms.size_major after clear: {}", ms.size_major());
    CHECK(ms.size_major() == 0);
    CHECK(mc.size_major() == 0);
    CHECK(cs.size_major() == 0);
    CHECK(cc.size_major() == 0);

    debug("ms.num_elements after clear: {}", ms.num_elements());
    CHECK(ms.num_elements() == 0);
    CHECK(mc.num_elements() == 0);
    CHECK(cs.num_elements() == 0);
    CHECK(cc.num_elements() == 0);

}

TEST_CASE("point cloud array 2d")
{
    Array::shape_t shape = {2,5};
    std::vector<int> user(shape[0]*shape[1],0);
    Array arr(user, shape, true);
    CHECK(arr.shape() == shape);
    arr.append(user);
    shape[0] += shape[0];
    CHECK(arr.shape() == shape);    

    std::vector<double> fuser(10,0);
    bool caught = false;
    try {
        arr.append(fuser);      // should throw
        debug("failed to throw on type size error");
    }
    catch (const ValueError& err) {
        debug("error on type mismatch as expected");
        caught=true;
    }
    CHECK(caught);
    
    CHECK(arr.shape() == shape);    
    CHECK(arr.size_major() == 4);

    {
        auto ma = arr.indexed<int, 2>();
        CHECK(ma.num_dimensions() == 2);
        CHECK(ma.num_elements() == 20);
        CHECK(ma.shape()[0] == 4);
        CHECK(ma.shape()[1] == 5);
    }

    {                           // wrong ndims
        bool caught=false;
        try {
            auto ma = arr.indexed<int, 3>();
        }
        catch (const ValueError& err) {
            caught = true;
        }
        CHECK(caught);
    }

    {                           // wrong element size
        bool caught=false;
        try {
            auto ma = arr.indexed<double, 2>();
        }
        catch (const ValueError& err) {
            caught = true;
        }
        CHECK(caught);
    }


}

TEST_CASE("point cloud dataset")
{
    Dataset d({
        {"one", Array({1  ,2  ,3  })},
        {"two", Array({1.1,2.2,3.3})},
    });
    CHECK(d.size_major() == 3);
    CHECK(d.keys().size() == 2);
    CHECK(d.size() == 2);

    {
        auto d2 = d;
        Array::shape_t shape = {3,4};
        std::vector<int> user(shape[0]*shape[1],0);
        Array arr(user, shape, true);
        d2.add("twod", arr);
    }
    {
        auto d2 = d;
        Array::shape_t shape = {2,5};
        std::vector<int> user(shape[0]*shape[1],0);
        Array arr(user, shape, true);
        bool caught = false;
        try {
            d2.add("twod", arr);
        }
        catch (const ValueError& err) {
            caught = true;
        }
        CHECK(caught);
    }

    {
        auto sel = d.selection({"one", "two"});
        CHECK(sel.size() == 2);
        const auto arr = sel[0];
        CHECK(arr->size_major() == 3);
        CHECK(sel[1]->size_major() == 3);
    }

    {
        Dataset d1(d);
        CHECK(d.size_major() == 3);
        CHECK(d1.size_major() == 3);
        CHECK(d1.keys().size() == 2);
        
        Dataset d2 = d;
        CHECK(d.size_major() == 3);
        CHECK(d2.size_major() == 3);
        CHECK(d2.keys().size() == 2);
    
        Dataset d3 = std::move(d2);
        CHECK(d2.size_major() == 0);
        CHECK(d3.size_major() == 3);
        CHECK(d2.keys().size() == 0);
        CHECK(d3.keys().size() == 2);
    }
    CHECK(d.size_major() == 3);

    {
        Dataset d({
            {"one", Array({1  ,2  ,3  })},
            {"two", Array({1.1,2.2,3.3})},
        });
        CHECK(d.size_major() == 3);
        CHECK(d.keys().size() == 2);

        size_t beg=0, end=0;
        d.register_append([&](size_t b, size_t e) { beg=b; end=e; });

        Dataset tail;
        tail.add("one", Array({4  , 5}));
        CHECK(tail.size_major() == 2);

        bool caught = false;
        try {
            tail.add("two", Array({4.4}));
        }
        catch (const ValueError& err) {
            caught = true;
        }
        CHECK(caught);

        tail.add("two", Array({4.4, 5.4}));
        CHECK(tail.keys().size() == 2);
        CHECK(tail.size_major() == 2);

        d.append(tail);
        debug("HAS {}", d.size_major());

        CHECK(beg == 3);
        CHECK(end == 5);

        CHECK(d.keys().size() == 2);
        CHECK(d.size_major() == 5);

        {
            auto tail2 = d.zeros_like(7);
            CHECK(tail2.size_major() == 7);

            auto arr2 = tail2.get("one");
            auto one2 = arr2->indexed<int, 1>();
            one2[0] = 42;
            CHECK(arr2->element<int>(0) == 42);

            d.append(tail2);
            CHECK(d.size_major() == 12);
            auto arr22 = d.get("one");
            CHECK(arr22->size_major() == 12);
            CHECK(arr22->element<int>(5) == 42);
        }
    }        

}
