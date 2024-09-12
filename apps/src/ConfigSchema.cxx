#include "WireCellApps/ConfigSchema.h"

#include "WireCellIface/INode.h"

#include "WireCellUtil/String.h"
#include "WireCellUtil/Type.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/ConfigManager.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Exceptions.h"


WIRECELL_FACTORY(ConfigSchema,
                 WireCellApps::ConfigSchema,
                 WireCell::IApplication,
                 WireCell::IConfigurable)

using spdlog::info;
using spdlog::warn;

using namespace WireCell;
using namespace WireCellApps;

ConfigSchema::ConfigSchema()
  : m_cfg(default_configuration())
{
}

ConfigSchema::~ConfigSchema() {}

void ConfigSchema::configure(const Configuration& config) { m_cfg = config; }

WireCell::Configuration ConfigSchema::default_configuration() const
{
    // yo dawg, I heard you liked dumping so I made a dumper that dumps the dumper.
    Configuration cfg;
    cfg["filename"] = "/dev/stdout";
    cfg["components"] = Json::arrayValue;
    return cfg;
}

// the "name" is the local type name.  We don't really have that available so
// instead if the cfg is a value of an object, the name will be the attribute
// key name.  O.w., empty.  This is something to fix by hand after the dump.
Configuration dump_schema(const Configuration& cfg, const std::string& name="");
Configuration dump_object(const Configuration& cfg, const std::string& name="");
Configuration dump_number(const Configuration& cfg, const std::string& name="");
Configuration dump_string(const Configuration& cfg, const std::string& name="");
Configuration dump_array(const Configuration& cfg, const std::string& name="");
Configuration dump_bool(const Configuration& cfg, const std::string& name="");
Configuration dump_null(const Configuration& cfg, const std::string& name="");

static
void maybe_name(Configuration& schema, const std::string& name)
{
    if (name.empty()) return;
    schema["name"] = name;
}

Configuration dump_object(const Configuration& cfg, const std::string& name)
{
    Configuration schema;
    schema["schema"] = "record";
    maybe_name(schema, name);
    schema["fields"] = Json::arrayValue;

    for (const auto& key : cfg.getMemberNames()) {
        auto val = cfg[key];
        Configuration field;
        field["name"] = key;
        field["item"] = dump_schema(val, key);
        if (! val.empty()) {
            field["default"] = val;
        }
        schema["fields"].append(field);
    }

    return schema;
}

Configuration dump_number(const Configuration& cfg, const std::string& name)
{
    Configuration schema;
    schema["schema"] = "number";
    maybe_name(schema, name);
    // omit: path, deps.

    // How can we possibly guess these these distinctions?
    if(cfg.isInt()) {
        schema["dtype"] = "i4";
    }
    else {
        schema["dtype"] = "f8";
    }

    return schema;
}

Configuration dump_string(const Configuration& cfg, const std::string& name)
{
    Configuration schema;
    schema["schema"] = "string";
    maybe_name(schema, name);
    // omit: path, deps, pattern/format.

    return schema;
}

Configuration dump_array(const Configuration& cfg, const std::string& name)
{
    Configuration schema;
    schema["schema"] = "sequence";
    maybe_name(schema, name);

    if (cfg.empty()) {
        schema["item"] = "null";
    }
    else {
        schema["item"] = dump_schema(cfg[0]);
    }
    return schema;
}

Configuration dump_bool(const Configuration& cfg, const std::string& name)
{
    Configuration schema;
    schema["schema"] = "boolean";
    maybe_name(schema, name);

    return schema;
}

Configuration dump_null(const Configuration& cfg, const std::string& name)
{
    Configuration schema;
    schema["schema"] = "null";
    maybe_name(schema, name);

    return schema;
}


Configuration dump_schema(const Configuration& cfg, const std::string& name)
{
    if (cfg.isObject()) {
        return dump_object(cfg);
    }
    if (cfg.isNumeric()) {
        return dump_number(cfg);
    }
    if (cfg.isString()) {
        return dump_string(cfg);
    }
    if (cfg.isArray()) {
        return dump_array(cfg);
    }
    if (cfg.isBool()) {
        return dump_bool(cfg);
    }
    if (cfg.isNull()) {
        return dump_null(cfg);
    }
    raise<ValueError>("unknown type: %s", cfg);
    // quell compiler warning about no return value.
    return Configuration();
}

void ConfigSchema::execute()
{
    // ConfigManager cm;
    int nfailed = 0;

    std::vector<std::string> comps;
    for (auto jone : m_cfg["components"]) {
        comps.push_back(jone.asString());
    }
    
    if (comps.empty()) {
        comps = Factory::known_types<IConfigurable>();
    }
    
    std::unordered_map<INode::NodeCategory, std::string> categories = {
        { INode::unknown,"unknown" },
        { INode::sourceNode, "source" },
        { INode::sinkNode, "sink" },
        { INode::functionNode, "function" }, 
        { INode::queuedoutNode, "queuedout" },
        { INode::joinNode, "join" },
        { INode::splitNode, "split" },
        { INode::faninNode, "fanin" },
        { INode::fanoutNode, "fanout" },
        { INode::multioutNode, "multiout" },
        { INode::hydraNode, "hydra" },
    };

    Configuration all_schema;

    for (auto const& [wctype, tinfo] : NamedFactoryRegistryBase::known_type_info()) {
        Configuration schema;
        schema["cppname"] = tinfo.cppname;

        schema["interfaces"] = Json::arrayValue;
        for (auto const& intname : tinfo.intnames) {
            schema["interfaces"].append(intname);
        }

        auto inode = Factory::lookup_tn<INode>(wctype, true, true);
        if (inode) {
            Configuration node;

            node["category"] = categories[inode->category()];
            node["signature"]= demangle( inode->signature() );
            node["concurrency"] = inode->concurrency();

            {
                Configuration ports = Json::arrayValue;
                for (const auto& port : inode->input_types()) {
                    ports.append( demangle(port) );
                }
                node["iports"] = ports;
            }
            {
                Configuration ports = Json::arrayValue;
                for (const auto& port : inode->output_types()) {
                    ports.append( demangle(port) );
                }
                node["oports"] = ports;
            }
            schema["node"] = node;
        }
        all_schema[wctype] = schema;
    }

    for (auto c : comps) {
        auto [wctype, instname] = String::parse_pair(convert<std::string>(c));

        Configuration cfg;
        try {
            auto cfgobj = Factory::lookup<IConfigurable>(wctype, instname);
            cfg = cfgobj->default_configuration();
        }
        catch (FactoryException& fe) {
            warn("failed lookup component: \"{}\":\"{}\"", wctype, instname);
            ++nfailed;
            continue;
        }

        all_schema[wctype]["config"] = dump_schema(cfg);
    }

    Persist::dump(get<std::string>(m_cfg, "filename"), all_schema);
}
