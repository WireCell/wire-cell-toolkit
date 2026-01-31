#include "WireCellApps/SchemaDumper.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Type.h"
#include "WireCellUtil/Logging.h"

#include "WireCellIface/INode.h"

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

// Helper to convert INode::NodeCategory enum to string
static std::string node_category_to_string(INode::NodeCategory cat)
{
    switch (cat) {
        case INode::unknown: return "unknown";
        case INode::sourceNode: return "sourceNode";
        case INode::sinkNode: return "sinkNode";
        case INode::functionNode: return "functionNode";
        case INode::queuedoutNode: return "queuedoutNode";
        case INode::joinNode: return "joinNode";
        case INode::splitNode: return "splitNode";
        case INode::faninNode: return "faninNode";
        case INode::fanoutNode: return "fanoutNode";
        case INode::multioutNode: return "multioutNode";
        case INode::hydraNode: return "hydraNode";
        default: return "unknown";
    }
}

void SchemaDumper::execute()
{
    std::map<std::string, Json::Value> factories;

    // Walk all registered interface types using the global registry
    auto all_interfaces = Factory::all_interfaces();

    info("SchemaDumper: found {} registered interface types", all_interfaces.size());

    for (const auto& iface_info : all_interfaces) {
        const std::string& interface_name = iface_info.interface_name;

        // Get all known types for this interface
        std::vector<std::string> known_types;
        try {
            known_types = iface_info.get_known_types();
        }
        catch (const std::exception& e) {
            warn("SchemaDumper: failed to get known types for interface {}: {}",
                 interface_name, e.what());
            continue;
        }

        info("SchemaDumper: interface {} has {} known types",
             interface_name, known_types.size());

        for (const auto& classname : known_types) {
            // Initialize factory entry if it doesn't exist
            if (factories.find(classname) == factories.end()) {
                factories[classname] = Json::objectValue;
                factories[classname]["classname"] = classname;
                factories[classname]["interfaces"] = Json::arrayValue;
            }

            // Add this interface to the list
            factories[classname]["interfaces"].append(interface_name);

            // Try to instantiate to get concrete type and INode info
            // We do this only once (first time we encounter this class)
            if (!factories[classname].isMember("concrete_type")) {
                try {
                    auto instance = iface_info.instantiate(classname, true);
                    if (instance) {
                        factories[classname]["concrete_type"] = type(*instance);

                        // Check if this is an INode and extract node-specific information
                        auto node = std::dynamic_pointer_cast<INode>(instance);
                        if (node) {
                            Json::Value node_info = Json::objectValue;

                            // Get node category
                            auto cat = node->category();
                            node_info["category"] = node_category_to_string(cat);

                            // Get input types
                            auto input_types = node->input_types();
                            if (!input_types.empty()) {
                                node_info["input_types"] = Json::arrayValue;
                                for (const auto& itype : input_types) {
                                    node_info["input_types"].append(demangle(itype));
                                }
                            }

                            // Get output types
                            auto output_types = node->output_types();
                            if (!output_types.empty()) {
                                node_info["output_types"] = Json::arrayValue;
                                for (const auto& otype : output_types) {
                                    node_info["output_types"].append(demangle(otype));
                                }
                            }

                            // Get signature
                            std::string sig = node->signature();
                            if (!sig.empty()) {
                                node_info["signature"] = demangle(sig);
                            }

                            // Get concurrency
                            node_info["concurrency"] = node->concurrency();

                            factories[classname]["node"] = node_info;
                        }
                    }
                }
                catch (const std::exception& e) {
                    warn("SchemaDumper: failed to instantiate {}: {}", classname, e.what());
                    // Continue anyway - we still have the interface information
                }
                catch (...) {
                    warn("SchemaDumper: failed to instantiate {}: unknown error", classname);
                    // Continue anyway
                }
            }
        }
    }

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
    output["metadata"]["num_interfaces"] = (int)all_interfaces.size();

    // Dump to file
    Persist::dump(get<string>(m_cfg, "filename"), output);

    info("SchemaDumper: dumped {} factories across {} interfaces to {}",
         factories.size(), all_interfaces.size(), get<string>(m_cfg, "filename"));
}
