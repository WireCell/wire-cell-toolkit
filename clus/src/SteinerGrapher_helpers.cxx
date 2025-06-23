#include "WireCellUtil/Exceptions.h"
#include "SteinerGrapher.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

Steiner::Grapher::Grapher(Cluster& cluster, const Steiner::Grapher::Config& cfg, Log::logptr_t log)
    : log(log), m_cluster(cluster), m_config(cfg)
{

}



void Steiner::Grapher::put_point_cloud(PointCloud::Dataset&& pc, const std::string& name)
{
    m_cluster.local_pcs().emplace(name, pc);
}
void Steiner::Grapher::put_point_cloud(const PointCloud::Dataset& pc, const std::string& name)
{
    m_cluster.local_pcs().emplace(name, pc);
}


PointCloud::Dataset& Steiner::Grapher::get_point_cloud(const std::string& name)
{
    if (m_cluster.has_pc(name)) {
        return m_cluster.get_pc(name);
    }

    // Fixme? configure the scope?  for now, the default.
    const auto& sv = m_cluster.sv();

    put_point_cloud(sv.flat_coords(), name);

    // Note, if more than the x,y,z coordinates are needed we would replace
    // flat_coords() with something like:
    // sv.flat_pc("3d", {"x","y","z","charge"})

    // Return the in-place reference
    return m_cluster.get_pc(name);
}

Steiner::Grapher::graph_type& Steiner::Grapher::get_graph(const std::string& flavor)
{
    // If graph of given flavor does not exist, the Cluster knows how to make
    // three "reserved" flavors, "basic", "ctpc" and "relaxed".
    return m_cluster.find_graph(flavor, m_config.dv, m_config.pcts); // throws if no flavor
}
