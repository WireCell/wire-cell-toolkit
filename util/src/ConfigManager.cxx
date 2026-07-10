#include "WireCellUtil/ConfigManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"

#include <unordered_map>

using namespace std;
using namespace WireCell;

// A named logger so these messages flow through WCT's configured sinks
// instead of relying on the process-global default logger.
static Log::logptr_t logger() {
    static Log::logptr_t l = Log::logger("cfgmgr");
    return l;
}


ConfigManager::ConfigManager()
  : m_top(Json::arrayValue)
{
}
ConfigManager::~ConfigManager() {}

void ConfigManager::extend(Configuration more)
{
    // O(N+M) instead of O(N*M): build an index map of existing entries once,
    // then probe per new entry.  Bulk loads of large configs (thousands of
    // entries) used to scale quadratically through add()->index() linear scan.
    std::unordered_map<std::string, int> by_tn;
    by_tn.reserve(m_top.size() + more.size());
    for (int i = 0; i < static_cast<int>(m_top.size()); ++i) {
        const auto& c = m_top[i];
        by_tn[get<string>(c, "type") + "\t" + get<string>(c, "name")] = i;
    }

    // We own `more` (by-value parameter): move each entry into m_top instead
    // of deep-copying.  Read type/name into locals first because we may move
    // `one` afterwards.  Callers that std::move()-in their argument also skip
    // the caller-side deep-copy at the call boundary.
    for (auto& one : more) {
        const std::string type = get<string>(one, "type");
        const std::string name = get<string>(one, "name");
        const std::string key = type + "\t" + name;
        auto it = by_tn.find(key);
        if (it == by_tn.end()) {
            const int ind = m_top.size();
            m_top[ind] = std::move(one);
            by_tn[key] = ind;
        }
        else {
            logger()->warn("ConfigManager::extend() overwriting type=\"{}\" name=\"{}\"",
                           type, name);
            m_top[it->second] = std::move(one);
        }
    }
}

int ConfigManager::index(const std::string& type, const std::string& name) const
{
    int ind = -1;
    for (const auto& c : m_top) {
        ++ind;
        if (get<string>(c, "type") != type) {
            continue;
        }
        if (get<string>(c, "name") != name) {
            continue;
        }
        return ind;
    }
    return -1;
}

int ConfigManager::add(Configuration& cfg)
{
    const std::string type = get<string>(cfg, "type");
    const std::string name = get<string>(cfg, "name");

    int ind = this->index(type, name);
    if (ind < 0) {
        ind = m_top.size();
    }
    else {
        logger()->warn("ConfigManager:add() overwriting existing type=\"{}\" name=\"{}\"", type, name);
    }
    m_top[ind] = cfg;
    return ind;
}

int ConfigManager::add(Configuration& payload, const std::string& type, const std::string& name)
{
    Configuration cfg;
    cfg["data"] = payload;
    cfg["type"] = type;
    cfg["name"] = name;
    return add(cfg);
}

Configuration ConfigManager::at(int ind) const
{
    if (ind < 0 || ind >= size()) {
        return Configuration();
    }
    return m_top[ind];
}

std::vector<ConfigManager::ClassInstance> ConfigManager::configurables() const
{
    std::vector<ConfigManager::ClassInstance> ret;
    for (auto c : m_top) {
        ret.push_back(make_pair(get<string>(c, "type"), get<string>(c, "name")));
    }
    return ret;
}

Configuration ConfigManager::pop(int ind)
{
    if (ind < 0 || ind >= size()) {
        return Configuration();
    }
    Configuration ret;
    // In-place removal via jsoncpp; avoids cloning every non-popped entry into
    // a fresh `reduced` array and then re-assigning that array back to m_top.
    m_top.removeIndex(static_cast<Json::ArrayIndex>(ind), &ret);
    return ret;
}

void ConfigManager::clear_data()
{
    // For each top-level entry, drop its "data" sub-tree.  finalize() only
    // looks at "type" and "name" so the bulk payload can be released after
    // configure() has run.
    for (auto& c : m_top) {
        if (c.isObject()) {
            c.removeMember("data");
        }
    }
}
