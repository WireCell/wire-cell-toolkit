#include "WireCellApps/SchemaDumper.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Type.h"
#include "WireCellUtil/Logging.h"

#include "WireCellIface/INode.h"
#include "WireCellIface/ISourceNode.h"
#include "WireCellIface/ISinkNode.h"
#include "WireCellIface/IFunctionNode.h"
#include "WireCellIface/IQueuedoutNode.h"
#include "WireCellIface/IJoinNode.h"
#include "WireCellIface/ISplitNode.h"
#include "WireCellIface/IFaninNode.h"
#include "WireCellIface/IFanoutNode.h"
#include "WireCellIface/IHydraNode.h"

#include "WireCellIface/IDepo.h"
#include "WireCellIface/IDepoSource.h"
#include "WireCellIface/IDepoSink.h"
#include "WireCellIface/IDepoFilter.h"
#include "WireCellIface/IDepoFramer.h"
#include "WireCellIface/IDepoSet.h"
#include "WireCellIface/IDepoSetSource.h"
#include "WireCellIface/IDepoSetSink.h"
#include "WireCellIface/IDepoSetFilter.h"

#include "WireCellIface/IFrame.h"
#include "WireCellIface/IFrameSource.h"
#include "WireCellIface/IFrameSink.h"
#include "WireCellIface/IFrameFilter.h"

#include "WireCellIface/ICluster.h"
#include "WireCellIface/IClusterSource.h"
#include "WireCellIface/IClusterSink.h"

#include "WireCellIface/ITensorSet.h"
#include "WireCellIface/ITensorSetSource.h"
#include "WireCellIface/ITensorSetSink.h"

#include "WireCellIface/IBlobSet.h"
#include "WireCellIface/IBlobSetSource.h"
#include "WireCellIface/IBlobSetSink.h"

#include <set>
#include <map>

WIRECELL_FACTORY(SchemaDumper, WireCellApps::SchemaDumper, WireCell::IApplication, WireCell::IConfigurable)

using spdlog::info;
using spdlog::warn;

using namespace std;
using namespace WireCell;
using namespace WireCellApps;

SchemaDumper::SchemaDumper()
  : m_cfg(default_configuration())
{
}

SchemaDumper::~SchemaDumper() {}

void SchemaDumper::configure(const Configuration& config) { m_cfg = config; }

WireCell::Configuration SchemaDumper::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = "/dev/stdout";
    return cfg;
}

// Helper function to collect factory information for a given interface
template<typename IType>
void collect_factories(const std::string& interface_name,
                      std::map<std::string, Json::Value>& factories)
{
    auto known = Factory::known_types<IType>();
    for (const auto& classname : known) {
        // Initialize factory entry if it doesn't exist
        if (factories.find(classname) == factories.end()) {
            factories[classname] = Json::objectValue;
            factories[classname]["classname"] = classname;
            factories[classname]["interfaces"] = Json::arrayValue;
        }
        // Add this interface to the list
        factories[classname]["interfaces"].append(interface_name);

        // Try to get type information for the concrete class
        // We do this only once (first time we encounter this class)
        if (!factories[classname].isMember("concrete_type")) {
            try {
                auto instance = Factory::lookup<IType>(classname, "", true, true);
                if (instance) {
                    factories[classname]["concrete_type"] = type(*instance);
                }
            }
            catch (...) {
                // If we can't instantiate, that's okay - just skip the concrete type
            }
        }
    }
}

void SchemaDumper::execute()
{
    std::map<std::string, Json::Value> factories;

    // Collect factories from all major interface types
    collect_factories<IApplication>("IApplication", factories);
    collect_factories<IConfigurable>("IConfigurable", factories);

    // Node interfaces
    collect_factories<INode>("INode", factories);
    collect_factories<ISourceNodeBase>("ISourceNodeBase", factories);
    collect_factories<ISinkNodeBase>("ISinkNodeBase", factories);
    collect_factories<IFunctionNodeBase>("IFunctionNodeBase", factories);
    collect_factories<IQueuedoutNodeBase>("IQueuedoutNodeBase", factories);
    collect_factories<IJoinNodeBase>("IJoinNodeBase", factories);
    collect_factories<ISplitNodeBase>("ISplitNodeBase", factories);
    collect_factories<IFaninNodeBase>("IFaninNodeBase", factories);
    collect_factories<IFanoutNodeBase>("IFanoutNodeBase", factories);
    collect_factories<IHydraNodeBase>("IHydraNodeBase", factories);

    // Depo interfaces
    collect_factories<IDepoSource>("IDepoSource", factories);
    collect_factories<IDepoSink>("IDepoSink", factories);
    collect_factories<IDepoFilter>("IDepoFilter", factories);
    collect_factories<IDepoFramer>("IDepoFramer", factories);
    collect_factories<IDepoSetSource>("IDepoSetSource", factories);
    collect_factories<IDepoSetSink>("IDepoSetSink", factories);
    collect_factories<IDepoSetFilter>("IDepoSetFilter", factories);

    // Frame interfaces
    collect_factories<IFrameSource>("IFrameSource", factories);
    collect_factories<IFrameSink>("IFrameSink", factories);
    collect_factories<IFrameFilter>("IFrameFilter", factories);

    // Cluster interfaces
    collect_factories<IClusterSource>("IClusterSource", factories);
    collect_factories<IClusterSink>("IClusterSink", factories);

    // TensorSet interfaces
    collect_factories<ITensorSetSource>("ITensorSetSource", factories);
    collect_factories<ITensorSetSink>("ITensorSetSink", factories);

    // BlobSet interfaces
    collect_factories<IBlobSetSource>("IBlobSetSource", factories);
    collect_factories<IBlobSetSink>("IBlobSetSink", factories);

    // Build the final JSON structure
    Configuration output;
    output["factories"] = Json::objectValue;

    for (const auto& entry : factories) {
        output["factories"][entry.first] = entry.second;
    }

    // Add metadata
    output["metadata"] = Json::objectValue;
    output["metadata"]["generator"] = "WireCell::SchemaDumper";
    output["metadata"]["num_factories"] = (int)factories.size();

    // Dump to file
    Persist::dump(get<string>(m_cfg, "filename"), output);

    info("SchemaDumper: dumped {} factories to {}",
         factories.size(), get<string>(m_cfg, "filename"));
}
