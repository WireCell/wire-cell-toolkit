
#include "WireCellSpng/Reduce.h"
#include "WireCellSpng/Testing.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellUtil/PluginManager.h"

using namespace WireCell;
using namespace WireCell::SPNG;


TEST_SUITE("spng reduce") {

    TEST_CASE("spng reduce simple add") {

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
        CHECK(torch::all( out_tensor == t2 ).item<bool>());
    }
        
}
