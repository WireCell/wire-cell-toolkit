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
    m_grapher_config.sampler = Factory::find_tn<IBlobSampler>(cfg["blob_sampler"].asString());
    m_grapher_config.dv = Factory::find_tn<IDetectorVolumes>(cfg["detector_volumes"].asString());
    m_grapher_config.pcts = Factory::find_tn<WireCell::Clus::IPCTransformSet>(cfg["pc_transforms"].asString());
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
