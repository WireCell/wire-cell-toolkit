#include "WireCellApps/ConfigDumper.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/ConfigManager.h"
#include "WireCellUtil/Logging.h"

WIRECELL_FACTORY(ConfigDumper, WireCellApps::ConfigDumper, WireCell::IApplication, WireCell::IConfigurable)

using spdlog::info;
using spdlog::warn;

using namespace std;
using namespace WireCell;
using namespace WireCellApps;

ConfigDumper::ConfigDumper()
  : m_cfg(default_configuration())
{
}

ConfigDumper::~ConfigDumper() {}

void ConfigDumper::configure(const Configuration& config) { m_cfg = config; }

WireCell::Configuration ConfigDumper::default_configuration() const
{
    // yo dawg, I heard you liked dumping so I made a dumper that dumps the dumper.
    Configuration cfg;
    cfg["filename"] = "/dev/stdout";
    cfg["components"] = Json::arrayValue;
    return cfg;
}

void ConfigDumper::execute()
{
    ConfigManager cm;
    int nfailed = 0;

    std::vector<std::string> comps;
    for (auto jone : m_cfg["components"]) {
        comps.push_back(jone.asString());
    }
    if (comps.empty()) {
        comps = Factory::known_types<IConfigurable>();
    }

    for (auto c : comps) {
        string type, name;
        tie(type, name) = String::parse_pair(convert<string>(c));

        Configuration cfg;
        try {
            auto cfgobj = Factory::lookup<IConfigurable>(type, name);
            cfg = cfgobj->default_configuration();
        }
        catch (FactoryException& fe) {
            warn("failed lookup component: \"{}\":\"{}\"", type, name);
            ++nfailed;
            continue;
        }
        cm.add(cfg, type, name);
    }

    Persist::dump(get<string>(m_cfg, "filename"), cm.all());
    (void)nfailed; // unused, but useful for debugging
}
