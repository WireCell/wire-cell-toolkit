// TorchTensorSetService's conversions (wcfm doc 14): any-rank, any-of-four-dtype tensor sets go through a
// TorchScript module positionally, and a tuple result comes back as a tensor set with the same element types.
// The module is defined in-test (Module::define), so no model file is needed.

#include "WireCellUtil/doctest.h"

#include "WireCellPytorch/TorchTensorSetService.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"

#include <cstdint>
#include <vector>

using namespace WireCell;

namespace {
    template <typename T>
    ITensor::pointer mk(const ITensor::shape_t& shape, const std::vector<T>& v)
    {
        return std::make_shared<Aux::SimpleTensor>(shape, v.data());
    }
    template <typename T>
    std::vector<T> vals(const ITensor::pointer& t)
    {
        const T* p = reinterpret_cast<const T*>(t->data());
        return std::vector<T>(p, p + t->size() / sizeof(T));
    }
}  // namespace

TEST_CASE("tensor set through a scripted module: ranks, dtypes, tuple return, inputs untouched")
{
    torch::jit::script::Module m("m");
    m.define(R"(
def forward(self, a, b, c, d):
    # a: f4 [2,3], b: i8 [4], c: f8 [1], d: i4 [2,2]
    return (a.sum(1) * 2.0, b + 1, c * 3.0, d.long().sum().int().reshape(1), b.index_select(0, b[:2]))
)");
    const std::vector<float> a{1, 2, 3, 4, 5, 6};
    const std::vector<int64_t> b{1, 0, 3, 2};
    const std::vector<double> c{0.5};
    const std::vector<int32_t> d{1, 2, 3, 4};
    ITensor::vector tv{mk<float>({2, 3}, a), mk<int64_t>({4}, b), mk<double>({1}, c), mk<int32_t>({2, 2}, d)};
    Configuration md;
    md["tag"] = "x";
    auto in = std::make_shared<Aux::SimpleTensorSet>(7, md, std::make_shared<ITensor::vector>(tv));

    auto out = Pytorch::forward_tensorset(m, in, torch::Device(torch::kCPU));
    REQUIRE(out);
    CHECK(out->ident() == 7);
    CHECK(out->metadata()["tag"].asString() == "x");
    auto ot = *out->tensors();
    REQUIRE(ot.size() == 5);

    CHECK(ot[0]->dtype() == "f4");
    CHECK(ot[0]->shape() == ITensor::shape_t{2});
    CHECK((vals<float>(ot[0]) == std::vector<float>{12, 30}));
    CHECK(ot[1]->dtype() == "i8");
    CHECK((vals<int64_t>(ot[1]) == std::vector<int64_t>{2, 1, 4, 3}));
    CHECK(ot[2]->dtype() == "f8");
    CHECK(vals<double>(ot[2]) == std::vector<double>{1.5});
    CHECK(ot[3]->dtype() == "i4");
    CHECK(vals<int32_t>(ot[3]) == std::vector<int32_t>{10});
    CHECK((vals<int64_t>(ot[4]) == std::vector<int64_t>{0, 1}));

    // the inputs are wrapped, not copied, and must be unchanged
    CHECK(vals<float>(tv[0]) == a);
    CHECK(vals<int64_t>(tv[1]) == b);
}

TEST_CASE("a single tensor return, and an unsupported input type throws")
{
    torch::jit::script::Module m("m");
    m.define(R"(
def forward(self, a):
    return a * 2
)");
    const std::vector<float> a{1, 2};
    auto in = std::make_shared<Aux::SimpleTensorSet>(
        1, Configuration(), std::make_shared<ITensor::vector>(ITensor::vector{mk<float>({2}, a)}));
    auto out = Pytorch::forward_tensorset(m, in, torch::Device(torch::kCPU));
    REQUIRE(out->tensors()->size() == 1);
    CHECK((vals<float>(out->tensors()->at(0)) == std::vector<float>{2, 4}));

    const std::vector<uint8_t> u{1, 2};
    auto bad = std::make_shared<Aux::SimpleTensorSet>(
        1, Configuration(), std::make_shared<ITensor::vector>(ITensor::vector{mk<uint8_t>({2}, u)}));
    CHECK_THROWS(Pytorch::forward_tensorset(m, bad, torch::Device(torch::kCPU)));
}
