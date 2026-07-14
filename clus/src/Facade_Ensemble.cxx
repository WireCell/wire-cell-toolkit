#include "WireCellClus/Facade_Grouping.h"

#include "WireCellClus/Facade_Ensemble.h"

using namespace WireCell::Clus::Facade; 

bool Ensemble::has(const std::string& name) const
{
    std::vector<std::string> ret;
    for (const auto* child : children()) {
        if (name == child->get_name()) {
            return true;
        }
    }
    return false;
}

std::vector<std::string> Ensemble::names() const
{
    std::vector<std::string> ret;
    for (const auto* child : children()) {
        ret.push_back(child->get_name());
    }
    return ret;
}

std::set<std::string> Ensemble::unique_names() const
{
    std::set<std::string> ret;
    for (const auto* child : children()) {
        ret.insert(child->get_name());
    }
    return ret;
}


std::vector<const Grouping*> Ensemble::with_name(const std::string& name) const
{
    std::vector<const Grouping*> ret;
    for (const auto* child : children()) {
        if (child->get_name() == name) {
            ret.push_back(child);
        }
    }
    return ret;
}
std::vector<Grouping*> Ensemble::with_name(const std::string& name) 
{
    std::vector<Grouping*> ret;
    for (auto* child : children()) {
        if (child->get_name() == name) {
            ret.push_back(child);
        }
    }
    return ret;
}

Grouping& Ensemble::make_grouping(const std::string& name)
{
    auto* pnode = m_node->insert();
    Grouping* grouping = pnode->value.facade<Grouping>();
    grouping->set_name(name);
    return *grouping;
}

Grouping& Ensemble::add_grouping_node(const std::string& name, points_t::node_ptr&& gnode)
{
    auto* pnode = m_node->insert(std::move(gnode));
    Grouping* grouping = pnode->value.facade<Grouping>();
    grouping->set_name(name);
    return *grouping;
}

std::map<std::string, Grouping*> Ensemble::groupings_by_name()
{
    std::map<std::string, Grouping*> ret;
    for (auto* child : children()) {
        auto name = child->get_name();
        if (ret.find(name) == ret.end()) {
            ret[name] = child;
        }
    }
    return ret;
}

