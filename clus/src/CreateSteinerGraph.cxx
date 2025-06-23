#include "CreateSteinerGraph.h"
#include "SteinerGrapher.h"

#include "WireCellClus/Graphs.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(CreateSteinerGraph, WireCell::Clus::Steiner::CreateSteinerGraph,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)


using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

Steiner::CreateSteinerGraph::CreateSteinerGraph()
    : Aux::Logger("CreateSteinerGraph", "clus")
{
}
Steiner::CreateSteinerGraph::~CreateSteinerGraph()
{
}


void Steiner::CreateSteinerGraph::configure(const WireCell::Configuration& cfg)
{
    m_grouping_name = get(cfg, "grouping", m_grouping_name);
    m_graph_name = get(cfg, "graph", m_graph_name);
    m_replace = get(cfg, "replace", m_replace);

    NeedDV::configure(cfg);
    NeedPCTS::configure(cfg);

    m_grapher_config.dv = m_dv;
    m_grapher_config.pcts = m_pcts;

    /// Do we even need samplers?
    // m_grapher_config.samplers.clear()
    // if (cfg.isMember("samplers") && cfg["samplers"].isArray()) {
    //     // Process array of samplers
    //     for (const auto& sampler_cfg : cfg["samplers"]) {
    //         int apa = sampler_cfg["apa"].asInt();
    //         int face = sampler_cfg["face"].asInt();
    //         std::string sampler_name = sampler_cfg["name"].asString();
            
    //         if (sampler_name.empty()) {
    //             raise<ValueError>("RetileCluster requires an IBlobSampler name for APA %d face %d", apa, face);
    //         }
    //         auto sampler_ptr = Factory::find_tn<IBlobSampler>(sampler_name);
    //         m_grapher_config.samplers[apa][face] = sampler_ptr;
    //     }
    // }
}

Configuration Steiner::CreateSteinerGraph::default_configuration() const
{
    Configuration cfg;
    // Build the Steiner graph for clusters in this grouping.
    cfg["grouping"] = m_grouping_name;
    // Name of the resulting graph on the cluster
    cfg["graph"] = m_graph_name;
    // If true, replace any pre-existing graph with that name, else do
    // nothing if one already exists.
    cfg["replace"] = m_replace;

    return cfg;
}


void Steiner::CreateSteinerGraph::visit(Ensemble& ensemble) const
{
    auto& grouping = *ensemble.with_name(m_grouping_name).at(0);
    for (auto* cluster : grouping.children()) {
        
        bool already = cluster->has_graph(m_graph_name);
        if (already || m_replace) {
            Steiner::Grapher sg(*cluster, m_grapher_config, log);
            auto gr = sg.create_steiner_graph();
            // Do we do any tests, eg on num_vertices()?
            cluster->give_graph(m_graph_name, std::move(gr));
        }
    }
}
