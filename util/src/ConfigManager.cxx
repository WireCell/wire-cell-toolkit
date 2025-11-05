#include "WireCellUtil/ConfigManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"

using namespace std;
using namespace WireCell;


ConfigManager::ConfigManager()
  : m_top(Json::arrayValue)
{
}
ConfigManager::~ConfigManager() {}

void ConfigManager::extend(Configuration more)
{
    // m_top = append(m_top, more);

    for (auto one : more) {
        add(one);
    }
}

int ConfigManager::index(const std::string& type, const std::string& name) const
{
    int ind = -1;
    for (auto c : m_top) {
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
        spdlog::warn("ConfigManager:add() overwriting existing type=\"{}\" name=\"{}\"", type, name);
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
    Configuration reduced(Json::arrayValue);
    int siz = size();
    for (int i = 0; i < siz; ++i) {
        if (i == ind) {
            ret = m_top[i];
        }
        else {
            reduced.append(m_top[i]);
        }
    }
    m_top = reduced;
    return ret;
}
