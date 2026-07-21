#include "WireCellSpng/Testing.h"
#include "WireCellSpng/TensorIndex.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/TensorSelector.h"
#include "WireCellSpng/TensorRenaming.h"
#include "WireCellSpng/Stopwatch.h"
#include "WireCellUtil/Exceptions.h"

using namespace WireCell;
using namespace WireCell::SPNG;

static TensorIndex make_ti()
{
    torch::Tensor dummy = torch::ones({2,2});
    TensorIndex ti;
    {
        Configuration md;
        md["datatype"] = "dummy"; // doesn't need to be known, just defined
        md["datapath"] = "/path/to/dummy";
        auto iten = std::make_shared<SimpleTorchTensor>(dummy, md);
        ti.add(iten);           // doesn't throw
    }
    {
        Configuration md;
        md["datatype"] = "kid"; // doesn't need to be known, just defined
        md["datapath"] = "/path/to/kid";
        md["parent"] = "/path/to/dummy";
        auto iten = std::make_shared<SimpleTorchTensor>(dummy, md);
        ti.add(iten);           // doesn't throw
    }
    {
        Configuration md;
        md["datatype"] = "kid"; // doesn't need to be known, just defined
        md["datapath"] = "/path/to/kid2";
        md["parent"] = "/path/to/dummy";
        auto iten = std::make_shared<SimpleTorchTensor>(dummy, md);
        ti.add(iten);           // doesn't throw
    }
    return ti;
}

TEST_CASE("spng tensor selector speed")
{
    Configuration cfg;
    cfg["tensor_selection"][0]["accept"] = "/path/to/dummy";
    TensorSelector ts;
    ts.configure(cfg);

    TensorIndex ti = make_ti();

    Stopwatch timer;

    const int ntries=1000;
    for (int count=0; count < ntries; ++count) { 
        auto ti2 = ts.apply(ti);
    }
    std::cerr << ntries << " selections in " << timer.lap() << " us\n";
}

TEST_CASE("spng tensor renaming speed")
{
    Configuration cfg;
    cfg["tensor_renaming"][0]["match"] = "(.*)/dummy";
    cfg["tensor_renaming"][0]["replace"] = "$1/smarty";
    TensorRenaming tr;
    tr.configure(cfg);

    TensorIndex ti = make_ti();

    Stopwatch timer;

    const int ntries=1000;
    for (int count=0; count < ntries; ++count) { 
        auto ti2 = tr.apply(ti);
    }
    std::cerr << ntries << " renamings in " << timer.lap() << " us\n";        
}

TEST_CASE("spng tensor index speed")
{
    const int ntries=1000;
    Stopwatch timer;

    for (int count=0; count < ntries; ++count) { 
        TensorIndex ti = make_ti();
    }

    std::cerr << ntries << " indexing in " << timer.lap() << " us\n";        
}

TEST_CASE("spng tensor index") {

    TensorIndex ti;             // start empty
    REQUIRE(ti.ident() < 0);    // we gave no ident, so the ident is not legit
    REQUIRE(ti.nparents() == 0);
    REQUIRE(ti.tensors().size() == 0);


    torch::Tensor dummy = torch::ones({2,2});

    {
        // no TDM datatype and datapath metadata attributes cause exception
        auto iten = std::make_shared<SimpleTorchTensor>(dummy);
        CHECK_THROWS_AS(ti.add(iten), ValueError);
    }
    {
        Configuration md;
        md["datatype"] = "dummy"; // doesn't need to be known, just defined
        md["datapath"] = "/path/to/dummy";
        auto iten = std::make_shared<SimpleTorchTensor>(dummy, md);
        ti.add(iten);           // doesn't throw
        REQUIRE(ti.nparents() == 1);
    }
    {
        Configuration md;
        md["datatype"] = "kid"; // doesn't need to be known, just defined
        md["datapath"] = "/path/to/kid";
        md["parent"] = "/path/to/does not exist";
        auto iten = std::make_shared<SimpleTorchTensor>(dummy, md);
        CHECK_THROWS_AS(ti.add(iten), ValueError); // parent not known
    }
    {
        Configuration md;
        md["datatype"] = "kid"; // doesn't need to be known, just defined
        md["datapath"] = "/path/to/kid";
        md["parent"] = "/path/to/dummy";
        auto iten = std::make_shared<SimpleTorchTensor>(dummy, md);
        ti.add(iten);           // doesn't throw
        REQUIRE(ti.nparents() == 1);
        REQUIRE(ti.parents().size() == 1);
    }
    {
        Configuration md;
        md["datatype"] = "kid"; // doesn't need to be known, just defined
        md["datapath"] = "/path/to/kid2";
        md["parent"] = "/path/to/dummy";
        auto iten = std::make_shared<SimpleTorchTensor>(dummy, md);
        ti.add(iten);           // doesn't throw
        REQUIRE(ti.nparents() == 1);
        REQUIRE(ti.parents().size() == 1);

        auto parent = ti.parent(iten);
        REQUIRE(parent);
        REQUIRE(parent->metadata()["datapath"] == "/path/to/dummy");
    }
    {
        auto tens = ti.tensors();
        REQUIRE(tens.size() == 3);
        REQUIRE(ti.of_type("dummy").size() == 1);
        REQUIRE(ti.of_type("kid").size() == 2);
        REQUIRE(ti.at_path("/path/to/dummy"));
        REQUIRE(ti.at_path("/path/to/kid"));
        REQUIRE(ti.at_path("/path/to/kid2"));
    }
    
    ///
    /// Test selector and renaming.
    ///
    /// Don't pull you hair out testing regex in compiled code.  Use util's
    /// check-regex CLI.
    ///

    {                           // kids live with parent
        Configuration cfg;
        cfg["tensor_selection"][0]["accept"] = "/path/to/dummy";
        TensorSelector ts;
        ts.configure(cfg);
        auto ti2 = ts.apply(ti);
        REQUIRE(ti2.at_path("/path/to/dummy"));
        REQUIRE(ti2.at_path("/path/to/kid"));
        REQUIRE(ti2.at_path("/path/to/kid2"));
    }
    {                           // kids die with parent
        Configuration cfg;
        cfg["tensor_selection"][0]["reject"] = "/path/to/dummy";
        TensorSelector ts;
        ts.configure(cfg);
        auto ti2 = ts.apply(ti);
        REQUIRE(nullptr == ti2.at_path("/path/to/dummy"));
        REQUIRE(nullptr == ti2.at_path("/path/to/kid"));
        REQUIRE(nullptr == ti2.at_path("/path/to/kid2"));
    }
    {                           // default accept
        Configuration cfg;
        TensorSelector ts;
        ts.configure(cfg);
        auto ti2 = ts.apply(ti, true);
        REQUIRE(ti2.at_path("/path/to/dummy"));
    }
    {                           // default reject
        Configuration cfg;
        TensorSelector ts;
        ts.configure(cfg);
        auto ti2 = ts.apply(ti, false);
        REQUIRE(nullptr == ti2.at_path("/path/to/dummy"));
    }

    {
        Configuration cfg;
        cfg["tensor_renaming"][0]["match"] = "(.*)/dummy";
        cfg["tensor_renaming"][0]["replace"] = "$1/smarty";
        TensorRenaming tr;
        tr.configure(cfg);
        auto ti2 = tr.apply(ti);

        REQUIRE(nullptr == ti2.at_path("/path/to/dummy"));
        REQUIRE(ti2.at_path("/path/to/smarty"));
        REQUIRE(ti2.at_path("/path/to/kid"));
        REQUIRE(ti2.at_path("/path/to/kid2"));
        REQUIRE(ti2.at_path("/path/to/kid")->metadata()["parent"] == "/path/to/smarty");
        REQUIRE(ti2.at_path("/path/to/kid2")->metadata()["parent"] == "/path/to/smarty");
    }

}
