#include "WireCellSpng/TorchFans.h"
#include "WireCellSpng/Testing.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/INamed.h"
#include "WireCellUtil/PluginManager.h"

using namespace WireCell;
using namespace WireCell::SPNG;

TEST_CASE("spng fanbase config")
{
    auto& pm = WireCell::PluginManager::instance();
    pm.add("WireCellSpng");

    auto icfg = Factory::lookup<IConfigurable>("SPNGFaninTensorSets");
    auto cfg = icfg->default_configuration();
    std::cout << "Default for FaninTensorSets cfg:\n" << cfg << "\n";
    icfg->configure(cfg);
    
    auto inamed = Factory::find_tn<INamed>("SPNGFaninTensorSets");
    std::cout << inamed->get_name() << "\n";
    inamed->set_name("SPNGFaninTensorSets");
    std::cout << inamed->get_name() << "\n";
}
