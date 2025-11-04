#include "WireCellSpng/TorchFans.h"
#include "WireCellSpng/Testing.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/INamed.h"
#include "WireCellUtil/PluginManager.h"

using namespace WireCell;
using namespace WireCell::SPNG;

TEST_CASE("spng fanbase get")
{
    auto& pm = WireCell::PluginManager::instance();
    pm.add("WireCellSpng");

    auto icfg = Factory::lookup<IConfigurable>("SPNGFaninTensorSets");
    auto cfg = icfg->default_configuration();
    std::cout << "Default for FaninTensorSets cfg:\n" << cfg << "\n";
    cfg["verbosity"] = 2;
    cfg["device"] = "gpu";
    std::cout << "Configured FaninTensorSets cfg:\n" << cfg << "\n";
    icfg->configure(cfg);
    
    auto inamed = Factory::find_tn<INamed>("SPNGFaninTensorSets");
    std::cout << inamed->get_name() << "\n";
    inamed->set_name("SPNGFaninTensorSets");
    std::cout << inamed->get_name() << "\n";

    auto ifan = std::dynamic_pointer_cast<FaninTensorSets>(inamed);
    REQUIRE(ifan != nullptr);
    std::cout << "device: " << ifan->device() << "\n";    
    CHECK(ifan->verbosity() == 2);

}

TEST_CASE("spng fanintensorsets config")
{
    FaninTensorSets fan;
    auto cfg = fan.default_configuration();
    std::cout << "Default for FaninTensorSets cfg:\n" << cfg << "\n";
    cfg["verbosity"] = 2;
    cfg["device"] = "gpu";
    fan.configure(cfg);
    CHECK(fan.verbosity() == 2);
    fan.logit("test");
    std::cout << "device: " << fan.device() << "\n";
}
