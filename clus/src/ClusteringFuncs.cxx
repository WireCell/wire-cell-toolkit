#include "WireCellClus/ClusteringFuncs.h"

#include "WireCellUtil/NamedFactory.h"

using namespace WireCell;


WireCell::Clus::IClusteringMethod::pointer WireCell::Clus::Facade::getClusteringMethod(const WireCell::Configuration& config)
{
    auto tn = config["name"].asString(); // this requires config to go from snake_case to CamelCase
    auto icfg = Factory::lookup_tn<IConfigurable>(tn);
    auto cfg = icfg->default_configuration();
    auto cfg2 = config; // fixme: make update() api const correct for 2nd arg....
    icfg->configure(update(cfg, cfg2));
    return Factory::find_tn<IClusteringMethod>(tn);
}
