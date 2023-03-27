#include "WireCellAux/TensorDM.h"
#include "WireCellAux/TensorTools.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellUtil/Testing.h"

#include <complex>
#include <vector>
#include <iostream>

using real_t = float;
using RV = std::vector<real_t>;
using complex_t = std::complex<real_t>;
using CV = std::vector<complex_t>;

// test fodder
const RV real_vector{0,1,2,3,4,5};
const RV real_vector_cw{0,3,1,4,2,5};
const CV complex_vector{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5}};
const WireCell::ITensor::shape_t shape{2,3};

using namespace WireCell;
using namespace WireCell::Aux;


void test_is_type()
{
    auto rt = std::make_shared<SimpleTensor>(shape, real_vector.data());
    Assert (Aux::is_type<real_t>(rt));
    Assert (!Aux::is_type<complex_t>(rt));
}

void test_is_row_major()
{
    // ST actually does not let us do anything but C-order/row-major
    auto rm = std::make_shared<SimpleTensor>(shape, real_vector.data());
    Assert(Aux::is_row_major(rm));
}

template<typename VectorType>
void assert_equal(const VectorType& v1, const VectorType& v2)
{
    Assert(v1.size() == v2.size());
    for (size_t ind=0; ind<v1.size(); ++ind) {
        Assert(v1[ind] == v2[ind]);
    }
}

// asvec 1) match type, 2) type mismatch
void test_asvec()
{
    auto rt = std::make_shared<SimpleTensor>(shape, real_vector.data());
    auto ct = std::make_shared<SimpleTensor>(shape, complex_vector.data());
    auto got_rt = Aux::asvec<real_t>(rt);
    auto got_ct = Aux::asvec<complex_t>(ct);
    assert_equal(real_vector, got_rt);
    assert_equal(complex_vector, got_ct);

    try {
        auto oops = Aux::asvec<complex_t>(rt);
    }
    catch (ValueError& err) {
    }
}

void test_asarray()
{
    // as array 2x2: (1d,2d) x (rw,cw)

    // make mutable copy to test that TT returns a copy
    RV my_vec(real_vector.begin(), real_vector.end());

    // test 2d
    auto rt = std::make_shared<SimpleTensor>(shape, my_vec.data());
    auto ra = Aux::asarray<real_t>(rt);
    auto shape = rt->shape();
    for (size_t irow = 0; irow < shape[0]; ++irow) {
        for (size_t icol = 0; icol < shape[1]; ++icol) {        
            Assert(ra(irow, icol) == my_vec[irow*shape[1] + icol]);
        }
    }

    // test 1d
    const WireCell::ITensor::shape_t shape1d{6,};
    auto rt1d = std::make_shared<SimpleTensor>(shape1d, my_vec.data());
    auto ra1d = Aux::asarray<real_t>(rt1d);
    for (size_t ind = 0; ind < shape[0]; ++ind) {        
        Assert(ra1d(ind) == my_vec[ind]);
    }

    // Assure the internal use of Eigen::Map leads to a copy on return
    my_vec[0] = 42;
    Assert(ra(0,0) == 0);
    Assert(ra1d(0) == 0);
}

template<typename ElementType>
void test_array_roundtrip()
{
    using namespace WireCell::PointCloud;
    using namespace WireCell::Aux;

    std::vector<ElementType> v{1,2,3};
    Array::shape_t shape = {3,};
    Array arr(v, shape, true);
    arr.metadata()["foo"] = "bar";

    Assert(arr.size_major() == 3);
    Assert(arr.is_type<ElementType>());

    Assert(arr.metadata()["foo"].asString() == "bar");
    auto itp = TensorDM::as_tensor(arr, "test");
    Assert(arr.bytes().data() != itp->data());
    Assert(arr.shape() == itp->shape());
    for (size_t ind=0; ind<v.size(); ++ind) {
        ElementType got = ((ElementType*)itp->data())[ind];

        std::cerr << ind
                  << ": orig = " << v[ind]
                  << " copy = " << got
                  << " type = " << arr.dtype() << "\n";

        Assert(v[ind] == got);
    }
    Assert(! itp->metadata()["foo"].isNull());
    Assert(itp->metadata()["foo"].isString());
    Assert(itp->metadata()["foo"].asString() == "bar");

    Array arr2 = TensorDM::as_array(itp, true);
    Assert(arr2.shape() == arr.shape());
    Assert(arr.size_major() == 3);
    Assert(arr.is_type<ElementType>());
    Assert(arr2.element<ElementType>(0) == arr.element<ElementType>(0));
    Assert(arr2.metadata()["foo"].asString() == "bar");
    Assert(arr2.metadata()["name"].isNull());

    for (size_t ind=0; ind<v.size(); ++ind) {
        Assert(v[ind] == arr2.element<ElementType>(ind));
    }

    arr2.assure_mutable();      // comment out and see asserts fail.
    itp = nullptr;              // trigger deletion

    Assert(arr2.shape() == arr.shape());
    Assert(arr.size_major() == 3);
    Assert(arr.is_type<ElementType>());
    Assert(arr2.element<ElementType>(0) == arr.element<ElementType>(0));
    Assert(arr2.metadata()["foo"].asString() == "bar");
    Assert(arr2.metadata()["name"].isNull());

    for (size_t ind=0; ind<v.size(); ++ind) {
        Assert(v[ind] == arr2.element<ElementType>(ind));
    }

}

void test_dataset_roundtrip()
{
    using namespace WireCell::PointCloud;
    using namespace WireCell::Aux;

    Array orig_one({1  ,2  ,3  });
    Array orig_two({1.1,2.2,3.3});

    Dataset::store_t s = {
        {"one", orig_one},
        {"two", orig_two},
    };

    Dataset d(s);
    {
        const auto& one = d.get("one");
        Assert(one.size_major() == 3);
        Assert(one.is_type<int>());
        const auto& two = d.get("two");
        Assert(two.size_major() == 3);
        Assert(two.is_type<double>());

        for (size_t ind = 0; ind<3; ++ind) {
            std::cerr << ind
                      << ": one: orig = " << orig_one.element<int>(ind)
                      << " copy = " << one.element<int>(ind) << "\n";
            Assert(one.element<int>(ind) == orig_one.element<int>(ind));
            std::cerr << ind
                      << ": two: orig = " << orig_two.element<double>(ind)
                      << " copy = " << two.element<double>(ind) << "\n";        
            Assert(two.element<double>(ind) == orig_two.element<double>(ind));
        }
    }

    const std::string datapath = "test";

    // We do not expect the round trip to step on our dataset metadata.
    const std::string not_an_ident_number = "this isn't an ident";
    d.metadata()["ident"] = not_an_ident_number;

    // to tensor

    auto itenvec = TensorDM::as_tensors(d, datapath);
    Assert(itenvec.size() > 0);
    auto itsp = TensorDM::as_tensorset(itenvec, 42);
    Assert(itsp->ident() == 42);
    auto itsmd = itenvec[0]->metadata();

    auto itv = itsp->tensors();
    {
        Assert(itv->size() == 2+1);
        for (size_t ind=0; ind<3; ++ind) {
            auto it = itv->at(ind);
            auto md = it->metadata();
            if (!ind) {         // the aggregator
                Assert(md["datatype"].asString() == "pcdataset");
            }
            else {              // an array
                Assert(md["datatype"].asString() == "pcarray");
            }
        }
    }
    {
        auto one = itv->at(1);
        auto two = itv->at(2);

        for (size_t ind=0; ind<3; ++ind) {
            Assert(orig_one.element<int>(ind) == ((int*)one->data())[ind]);
            Assert(orig_two.element<double>(ind) == ((double*)two->data())[ind]);
        }
    }
    

    // back to dataset

    bool share = false;
    Dataset d2 = TensorDM::as_dataset(itsp, datapath, share);

    Assert(d2.size_major() == 3);

    auto dmd = d2.metadata();
    Assert(dmd["ident"].asString() == not_an_ident_number);

    auto keys = d2.keys();
    Assert(keys.size() == 2);
    Assert(keys.find("one") != keys.end());
    Assert(keys.find("two") != keys.end());
    
    auto sel = d2.selection({"one", "two"});
    const Array& one = sel[0];
    const Array& two = sel[1];

    // const auto& one = d2.get("one");
    Assert(one.size_major() == 3);
    Assert(one.is_type<int>());
    Assert(one.dtype() == "i4");
    auto one_bytes = one.bytes();
    Assert(one_bytes.size() == 3*sizeof(int));
    // Assert(1 == ((int*)one_bytes.data())[0]);

    // const auto& two = d2.get("two");
    Assert(two.size_major() == 3);
    Assert(two.is_type<double>());

    Assert(one.metadata()["name"].isNull());

    for (size_t ind = 0; ind<3; ++ind) {
        std::cerr << ind
                  << ": one: orig = " << orig_one.element<int>(ind)
                  << " copy = " << one.element<int>(ind) << "\n";
        Assert(one.element<int>(ind) == orig_one.element<int>(ind));
        std::cerr << ind
                  << ": two: orig = " << orig_two.element<double>(ind)
                  << " copy = " << two.element<double>(ind) << "\n";        
        Assert(two.element<double>(ind) == orig_two.element<double>(ind));
    }

    return;
}

int main()
{
    test_is_type();
    test_is_row_major();
    test_asvec();
    test_asarray();

    test_array_roundtrip<int8_t>();
    test_array_roundtrip<int16_t>();
    test_array_roundtrip<int32_t>();
    test_array_roundtrip<int64_t>();
    test_array_roundtrip<uint8_t>();
    test_array_roundtrip<uint16_t>();
    test_array_roundtrip<uint32_t>();
    test_array_roundtrip<uint64_t>();
    test_array_roundtrip<float>();
    test_array_roundtrip<double>();
    test_array_roundtrip<std::complex<float>>();
    test_array_roundtrip<std::complex<double>>();

    test_dataset_roundtrip();

    return 0;
}
