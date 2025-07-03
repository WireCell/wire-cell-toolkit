#include "CreateSteinerGraph.h"
#include "SteinerGrapher.h"

#include "WireCellClus/Graphs.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellClus/ClusteringFuncs.h"


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
    const std::string retiler_tn = get<std::string>(cfg, "retiler", "RetileCluster");
    m_grapher_config.retile = Factory::find_tn<IPCTreeMutate>(retiler_tn);
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
        Steiner::Grapher sg(*cluster, m_grapher_config, log);
        
        bool already = cluster->has_graph(m_graph_name);
        if (already || m_replace) {

            if (cluster->get_flag(Flags::beam_flash)){
                // separate the clusters into separated pieces ...

                auto cc =cluster->get_pcarray("isolated", "perblob");
                // convert span to vector
                std::vector<int> cc_vec(cc.begin(), cc.end());
                std::cout << "Xin: " << cluster->ident() << " " << cc_vec.size() << std::endl;
                std::cout << "Xin: " << cluster->get_flash().time()/units::us << " " << cluster->nchildren() << " " << cluster->npoints() <<  " " << std::endl;
                
                // if (cc_vec.size() < 2) continue;
                // auto scope = cluster->get_default_scope();
                // auto scope_transform = cluster->get_scope_transform(scope);
                // // // origi_cluster still have the original main cluster ... 
                // // std::cout << "Start: " << orig_cluster->kd_blobs().size() << " " << orig_cluster->nchildren() << std::endl;
                auto splits = grouping.separate(cluster, cc_vec);

                // Apply the scope filter settings to all new clusters
                for (auto& [id, new_cluster] : splits) {
                    std::cout << "Xin: " << new_cluster->ident() << " " << new_cluster->kd_blobs().size() << " " << new_cluster->nchildren() << " " << id << std::endl;
                    //     new_cluster->set_scope_filter(scope, true);
                    //     new_cluster->set_default_scope(scope);
                    //     new_cluster->set_scope_transform(scope,scope_transform);
                }

            }
                      
        }

        

            //auto gr = sg.create_steiner_graph();
            // Do we do any tests, eg on num_vertices()?
            //cluster->give_graph(m_graph_name, std::move(gr));
    }
}
