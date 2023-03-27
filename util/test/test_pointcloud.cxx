#include "WireCellUtil/PointCloud.h"
#include "WireCellUtil/Testing.h"

#include <vector>
#include <iostream>

using namespace WireCell;
using namespace WireCell::PointCloud;

template<typename ElementType>
void assure_equal(const Array& a, const std::vector<ElementType>& v)
{
    std::cerr << " dtype=" << a.dtype() << " num_elements=" << a.num_elements() << " size=" << v.size() << "\n";
    Assert(a.is_type<ElementType>());
    Assert(v.size() == a.num_elements()); // these two are 
    Assert(v.size() == a.size_major());   // degenerate for
    Assert(a.shape().size() == 1);        // a 1D array.

    auto ele = a.elements<ElementType>();

    for (size_t ind=0; ind<v.size(); ++ind) {
        std::cerr << ind << ": " << v[ind] << " =?= " << ele[ind] << "\n";
        Assert(v[ind] == ele[ind]);
    }
}

template<typename ElementType>
PointCloud::Array make_array(const ElementType* data, const Array::shape_t& shape, bool share)
{
    PointCloud::Array arr(data, shape, share);
    return arr;
}
template<typename ElementType>
PointCloud::Dataset::store_t make_arrays(const ElementType* data, const Array::shape_t& shape, bool share)
{
    PointCloud::Array arr(data, shape, share);
    PointCloud::Dataset::store_t ret;
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

void test_array()
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
        Assert(flat_span.size() == v.size());
        ms.elements<float>(); // same size okay (for now)
        bool caught = true;
        try {
            ms.elements<double>();
        }
        catch (const ValueError& err) {
            caught = true;
        }
        Assert(caught);
    }

    v.push_back(4);
    ms.assign(v.data(), {v.size()}, true);
    mc.assign(v.data(), {v.size()}, false);
    cs.assign(w.data(), {w.size()}, true);
    cc.assign(w.data(), {w.size()}, false);

    Assert(ms.size_major() == 4);
    Assert(mc.size_major() == 4);
    Assert(cs.size_major() == 4);
    Assert(cc.size_major() == 4);

    assure_equal(ms, v);
    assure_equal(mc, v);
    assure_equal(cs, v);
    assure_equal(cc, v);

    v[0] = 42;

    Assert(ms.element<int>(0) == 42);
    Assert(mc.element<int>(0) == 1);
    Assert(cs.element<int>(0) == 42);
    Assert(cc.element<int>(0) == 1);

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

    Assert(ms.element<double>(0) == 42);
    Assert(mc.element<double>(0) == 1);
    Assert(cs.element<double>(0) == 42);
    Assert(cc.element<double>(0) == 1);

    auto msarr = ms.indexed<double, 1>();
    Assert(msarr[0] == 42);

    auto msele = ms.elements<double>();
    Assert(msele[0] == 42);

    std::vector<double> d2{4.,5.,6.};
    ms.append(d2);
    Assert(ms.element<double>(3) == 4);
    {
        Assert(ms.shape().size() == 1);
        Assert(ms.shape()[0] == 7);
    }

    ms.clear();
    mc.clear();
    cs.clear();
    cc.clear();

    std::cerr << "ms.size_major after clear: " << ms.size_major() << "\n";
    Assert(ms.size_major() == 0);
    Assert(mc.size_major() == 0);
    Assert(cs.size_major() == 0);
    Assert(cc.size_major() == 0);

    std::cerr << "ms.num_elements after clear: " << ms.num_elements() << "\n";
    Assert(ms.num_elements() == 0);
    Assert(mc.num_elements() == 0);
    Assert(cs.num_elements() == 0);
    Assert(cc.num_elements() == 0);

}

void test_array2d()
{
    Array::shape_t shape = {2,5};
    std::vector<int> user(shape[0]*shape[1],0);
    Array arr(user, shape, true);
    Assert(arr.shape() == shape);
    arr.append(user);
    shape[0] += shape[0];
    Assert(arr.shape() == shape);    

    std::vector<double> fuser(10,0);
    bool caught = false;
    try {
        arr.append(fuser);      // should throw
        std::cerr << "failed to through on type size error\n";
    }
    catch (const ValueError& err) {
        std::cerr << "error on type mismatch as expected\n";
        caught=true;
    }
    Assert(caught);
    
    Assert(arr.shape() == shape);    
    Assert(arr.size_major() == 4);

    {
        auto ma = arr.indexed<int, 2>();
        Assert(ma.num_dimensions() == 2);
        Assert(ma.num_elements() == 20);
        Assert(ma.shape()[0] == 4);
        Assert(ma.shape()[1] == 5);
    }

    {                           // wrong ndims
        bool caught=false;
        try {
            auto ma = arr.indexed<int, 3>();
        }
        catch (const ValueError& err) {
            caught = true;
        }
        Assert(caught);
    }

    {                           // wrong element size
        bool caught=false;
        try {
            auto ma = arr.indexed<double, 2>();
        }
        catch (const ValueError& err) {
            caught = true;
        }
        Assert(caught);
    }


}

void test_dataset()
{
    Dataset::store_t s = {
        {"one", Array({1  ,2  ,3  })},
        {"two", Array({1.1,2.2,3.3})},
    };
    Dataset d(s);
    Assert(d.size_major() == 3);
    Assert(d.keys().size() == 2);
    Assert(d.store().size() == 2);

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
        Assert(caught);
    }

    {
        auto sel = d.selection({"one", "two"});
        Assert(sel.size() == 2);
        const Array& arr = sel[0];
        Assert(arr.size_major() == 3);
        Assert(sel[1].get().size_major() == 3);
    }

    {
        Dataset d1(d);
        Assert(d.size_major() == 3);
        Assert(d1.size_major() == 3);
        Assert(d1.keys().size() == 2);
        
        Dataset d2 = d;
        Assert(d.size_major() == 3);
        Assert(d2.size_major() == 3);
        Assert(d2.keys().size() == 2);
    
        Dataset d3 = std::move(d2);
        Assert(d2.size_major() == 0);
        Assert(d3.size_major() == 3);
        Assert(d2.keys().size() == 0);
        Assert(d3.keys().size() == 2);
    }
    Assert(d.size_major() == 3);

    {
        Dataset::store_t s = {
            {"one", Array({1  ,2  ,3  })},
            {"two", Array({1.1,2.2,3.3})},
        };
        Dataset d(s);
        Assert(d.size_major() == 3);
        Assert(d.keys().size() == 2);

        size_t beg=0, end=0;
        d.register_append([&](size_t b, size_t e) { beg=b; end=e; });

        Dataset tail;
        tail.add("one", Array({4  , 5}));
        Assert(tail.size_major() == 2);

        bool caught = false;
        try {
            tail.add("two", Array({4.4}));
        }
        catch (const ValueError& err) {
            caught = true;
        }
        Assert(caught);

        tail.add("two", Array({4.4, 5.4}));
        Assert(tail.keys().size() == 2);
        Assert(tail.size_major() == 2);

        d.append(tail);
        std::cerr << "HAS " << d.size_major() << "\n";

        Assert(beg == 3);
        Assert(end == 5);

        Assert(d.keys().size() == 2);
        Assert(d.size_major() == 5);

        {
            auto tail2 = d.zeros_like(7);
            Assert(tail2.size_major() == 7);

            auto& arr2 = tail2.get("one");
            auto one2 = arr2.indexed<int, 1>();
            one2[0] = 42;
            Assert(arr2.element<int>(0) == 42);

            d.append(tail2);
            Assert(d.size_major() == 12);
            auto& arr22 = d.get("one");
            Assert(arr22.size_major() == 12);
            Assert(arr22.element<int>(5) == 42);
        }
    }        

}

int main()
{
    test_array();
    test_array_noshare();
    test_array2d();
    test_dataset();

    return 0;
}
