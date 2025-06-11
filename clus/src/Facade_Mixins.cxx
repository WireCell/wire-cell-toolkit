#include "WireCellClus/Facade_Mixins.h"

using namespace WireCell::Clus::Facade;

Mixins::Graphs::graph_type& Mixins::Graphs::make_graph(const std::string& name, size_t nvertices)
{
    m_graph_store[name] = graph_type(nvertices);
    return m_graph_store[name];
}

Mixins::Graphs::graph_type& Mixins::Graphs::give_graph(const std::string& name, Mixins::Graphs::graph_type&& gr)
{
    auto it = m_graph_store.find(name);
    if (it != m_graph_store.end()) {
        m_graph_store.erase(it);
    }
    auto it2 = m_graph_store.emplace(name, std::move(gr));
    return it2.first->second;
}

Mixins::Graphs::graph_type& Mixins::Graphs::get_graph(const std::string& name)
{
    auto it = m_graph_store.find(name);
    if (it == m_graph_store.end()) {
        return make_graph(name);
    }
    return it->second;
}

const Mixins::Graphs::graph_type& Mixins::Graphs::get_graph(const std::string& name) const
{
    auto it = m_graph_store.find(name);
    if (it == m_graph_store.end()) {
        raise<ValueError>("no graph with name " + name);
    }
    return it->second;
}

Mixins::Graphs::graph_type Mixins::Graphs::take_graph(const std::string& name)
{
    // not actually in C++17 for GCC at least?
    // auto entry = m_graph_store.extract(name);
    // if (entry) {
    //     return std::move(entry.value());
    // }
    auto it = m_graph_store.find(name);
    if (it == m_graph_store.end()) {
        return graph_type{};
    }
    auto g = std::move(it->second);
    m_graph_store.erase(it);
    return g;
}


