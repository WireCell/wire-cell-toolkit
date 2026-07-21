#include "WireCellSpng/Reduce.h"
#include "WireCellSpng/Testing.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellUtil/PluginManager.h"

using namespace WireCell;
using namespace WireCell::SPNG;


TEST_SUITE("spng reduce") {

    TEST_CASE("spng reduce sum zeros ones") {

        auto t1 = torch::zeros({2,3,4});
        auto t2 = torch::ones({2,3,4});

        SimpleTorchTensor::vector inv;
        inv.push_back(std::make_shared<SimpleTorchTensor>(t1));
        inv.push_back(std::make_shared<SimpleTorchTensor>(t2));

        Reduce r;
        auto cfg = r.default_configuration();
        cfg["operation"] = "sum";
        r.configure(cfg);

        ITorchTensor::pointer out;
        r.fanin_combine(inv, out);
        REQUIRE(out);
        auto out_tensor = out->tensor();
        // 0 + 1 = 1
        CHECK(torch::all( out_tensor == t2 ).item<bool>());
    }

    // --- Element-wise operations ---

    TEST_CASE("spng reduce sum 2d") {
        auto t1 = torch::tensor({{1.0, 2.0}, {3.0, 4.0}});
        auto t2 = torch::tensor({{10.0, 20.0}, {30.0, 40.0}});

        SimpleTorchTensor::vector inv;
        inv.push_back(std::make_shared<SimpleTorchTensor>(t1));
        inv.push_back(std::make_shared<SimpleTorchTensor>(t2));

        Reduce r;
        auto cfg = r.default_configuration();
        cfg["operation"] = "sum";
        r.configure(cfg);

        ITorchTensor::pointer out;
        r.fanin_combine(inv, out);
        REQUIRE(out);
        auto out_tensor = out->tensor();
        
        auto expected = torch::tensor({{11.0, 22.0}, {33.0, 44.0}});
        CHECK(torch::allclose(out_tensor, expected));
    }

    TEST_CASE("spng reduce min") {
        auto t1 = torch::tensor({{1.0, 2.0}, {3.0, 4.0}});
        auto t2 = torch::tensor({{10.0, 20.0}, {30.0, 40.0}});

        SimpleTorchTensor::vector inv;
        inv.push_back(std::make_shared<SimpleTorchTensor>(t1));
        inv.push_back(std::make_shared<SimpleTorchTensor>(t2));

        Reduce r;
        auto cfg = r.default_configuration();
        cfg["operation"] = "min";
        r.configure(cfg);

        ITorchTensor::pointer out;
        r.fanin_combine(inv, out);
        REQUIRE(out);
        auto out_tensor = out->tensor();
        
        // t1 is element-wise minimum
        CHECK(torch::allclose(out_tensor, t1));
    }

    TEST_CASE("spng reduce max") {
        auto t1 = torch::tensor({{1.0, 2.0}, {3.0, 4.0}});
        auto t2 = torch::tensor({{10.0, 20.0}, {30.0, 40.0}});

        SimpleTorchTensor::vector inv;
        inv.push_back(std::make_shared<SimpleTorchTensor>(t1));
        inv.push_back(std::make_shared<SimpleTorchTensor>(t2));

        Reduce r;
        auto cfg = r.default_configuration();
        cfg["operation"] = "max";
        r.configure(cfg);

        ITorchTensor::pointer out;
        r.fanin_combine(inv, out);
        REQUIRE(out);
        auto out_tensor = out->tensor();
        
        // t2 is element-wise maximum
        CHECK(torch::allclose(out_tensor, t2));
    }

    TEST_CASE("spng reduce mean") {
        auto t1 = torch::tensor({{1.0, 2.0}, {3.0, 4.0}});
        auto t2 = torch::tensor({{10.0, 20.0}, {30.0, 40.0}});

        SimpleTorchTensor::vector inv;
        inv.push_back(std::make_shared<SimpleTorchTensor>(t1));
        inv.push_back(std::make_shared<SimpleTorchTensor>(t2));

        Reduce r;
        auto cfg = r.default_configuration();
        cfg["operation"] = "mean";
        r.configure(cfg);

        ITorchTensor::pointer out;
        r.fanin_combine(inv, out);
        REQUIRE(out);
        auto out_tensor = out->tensor();
        
        auto expected = torch::tensor({{5.5, 11.0}, {16.5, 22.0}});
        CHECK(torch::allclose(out_tensor, expected));
    }

    TEST_CASE("spng reduce mul") {
        auto t1 = torch::tensor({{1.0, 2.0}, {3.0, 4.0}});
        auto t2 = torch::tensor({{10.0, 20.0}, {30.0, 40.0}});

        SimpleTorchTensor::vector inv;
        inv.push_back(std::make_shared<SimpleTorchTensor>(t1));
        inv.push_back(std::make_shared<SimpleTorchTensor>(t2));

        Reduce r;
        auto cfg = r.default_configuration();
        cfg["operation"] = "mul";
        r.configure(cfg);

        ITorchTensor::pointer out;
        r.fanin_combine(inv, out);
        REQUIRE(out);
        auto out_tensor = out->tensor();
        
        auto expected = torch::tensor({{10.0, 40.0}, {90.0, 160.0}});
        CHECK(torch::allclose(out_tensor, expected));
    }

    // --- Stacking/Concatenation operations ---

    TEST_CASE("spng reduce cat dim 1") {
        auto t1 = torch::tensor({{1.0, 2.0}, {3.0, 4.0}});
        auto t2 = torch::tensor({{10.0, 20.0}, {30.0, 40.0}});

        SimpleTorchTensor::vector inv;
        inv.push_back(std::make_shared<SimpleTorchTensor>(t1));
        inv.push_back(std::make_shared<SimpleTorchTensor>(t2));

        Reduce r;
        auto cfg = r.default_configuration();
        cfg["operation"] = "cat";
        cfg["dim"] = 1;
        r.configure(cfg);

        ITorchTensor::pointer out;
        r.fanin_combine(inv, out);
        REQUIRE(out);
        auto out_tensor = out->tensor();
        
        auto expected = torch::tensor({{1.0, 2.0, 10.0, 20.0}, {3.0, 4.0, 30.0, 40.0}});
        CHECK(torch::allclose(out_tensor, expected));
        CHECK(out_tensor.sizes().vec() == std::vector<int64_t>{2, 4});
    }

    TEST_CASE("spng reduce stack dim 0") {
        auto t1 = torch::tensor({{1.0, 2.0}, {3.0, 4.0}});
        auto t2 = torch::tensor({{10.0, 20.0}, {30.0, 40.0}});

        SimpleTorchTensor::vector inv;
        inv.push_back(std::make_shared<SimpleTorchTensor>(t1));
        inv.push_back(std::make_shared<SimpleTorchTensor>(t2));

        Reduce r;
        auto cfg = r.default_configuration();
        cfg["operation"] = "stack";
        cfg["dim"] = 0;
        r.configure(cfg);

        ITorchTensor::pointer out;
        r.fanin_combine(inv, out);
        REQUIRE(out);
        auto out_tensor = out->tensor();
        
        CHECK(out_tensor.sizes().vec() == std::vector<int64_t>{2, 2, 2});
        CHECK(torch::allclose(out_tensor[0], t1));
        CHECK(torch::allclose(out_tensor[1], t2));
    }

    TEST_CASE("spng reduce hstack") {
        auto t1 = torch::tensor({{1.0, 2.0}, {3.0, 4.0}});
        auto t2 = torch::tensor({{10.0, 20.0}, {30.0, 40.0}});

        SimpleTorchTensor::vector inv;
        inv.push_back(std::make_shared<SimpleTorchTensor>(t1));
        inv.push_back(std::make_shared<SimpleTorchTensor>(t2));

        Reduce r;
        auto cfg = r.default_configuration();
        cfg["operation"] = "hstack";
        r.configure(cfg);

        ITorchTensor::pointer out;
        r.fanin_combine(inv, out);
        REQUIRE(out);
        auto out_tensor = out->tensor();
        
        // hstack for 2D tensors is cat along dim 1
        auto expected = torch::tensor({{1.0, 2.0, 10.0, 20.0}, {3.0, 4.0, 30.0, 40.0}});
        CHECK(torch::allclose(out_tensor, expected));
        CHECK(out_tensor.sizes().vec() == std::vector<int64_t>{2, 4});
    }

    TEST_CASE("spng reduce vstack") {
        auto t1 = torch::tensor({{1.0, 2.0}, {3.0, 4.0}});
        auto t2 = torch::tensor({{10.0, 20.0}, {30.0, 40.0}});

        SimpleTorchTensor::vector inv;
        inv.push_back(std::make_shared<SimpleTorchTensor>(t1));
        inv.push_back(std::make_shared<SimpleTorchTensor>(t2));

        Reduce r;
        auto cfg = r.default_configuration();
        cfg["operation"] = "vstack";
        r.configure(cfg);

        ITorchTensor::pointer out;
        r.fanin_combine(inv, out);
        REQUIRE(out);
        auto out_tensor = out->tensor();
        
        // vstack for 2D tensors is cat along dim 0
        auto expected = torch::tensor({{1.0, 2.0}, {3.0, 4.0}, {10.0, 20.0}, {30.0, 40.0}});
        CHECK(torch::allclose(out_tensor, expected));
        CHECK(out_tensor.sizes().vec() == std::vector<int64_t>{4, 2});
    }

    TEST_CASE("spng reduce dstack") {
        auto t1 = torch::tensor({{1.0, 2.0}, {3.0, 4.0}});
        auto t2 = torch::tensor({{10.0, 20.0}, {30.0, 40.0}});

        SimpleTorchTensor::vector inv;
        inv.push_back(std::make_shared<SimpleTorchTensor>(t1));
        inv.push_back(std::make_shared<SimpleTorchTensor>(t2));

        Reduce r;
        auto cfg = r.default_configuration();
        cfg["operation"] = "dstack";
        r.configure(cfg);

        ITorchTensor::pointer out;
        r.fanin_combine(inv, out);
        REQUIRE(out);
        auto out_tensor = out->tensor();
        
        // dstack results in a 3D tensor (2, 2, 2)
        CHECK(out_tensor.sizes().vec() == std::vector<int64_t>{2, 2, 2});
        
        // Check slices along dim 2
        CHECK(torch::allclose(out_tensor.select(2, 0), t1));
        CHECK(torch::allclose(out_tensor.select(2, 1), t2));
    }
        
}
