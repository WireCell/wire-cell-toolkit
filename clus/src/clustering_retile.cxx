// This defines the component ClusteringRetile, a "clustering function" aka
// "ensemble visitor" which delegates to an IPCTreeTransform to produce a new
// set of clusters from a subset of an input set.  See retile_cluster.cxx.


#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/RayHelpers.h"

#include "WireCellIface/IBlob.h"
#include "WireCellIface/IBlobSampler.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IPCTreeMutate.h"

#include "WireCellAux/PlaneTools.h"

#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/ClusteringFuncsMixins.h"


#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"

#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SamplingHelpers.h"

#include "WireCellUtil/PointTree.h"

#include "WireCellAux/SimpleSlice.h"
#include "WireCellClus/GroupingHelper.h"



#include <vector>

class ClusteringRetile;
WIRECELL_FACTORY(ClusteringRetile, ClusteringRetile,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;



class ClusteringRetile : public IConfigurable, public Clus::IEnsembleVisitor, private Clus::NeedScope {
public:
    ClusteringRetile() {};
    virtual ~ClusteringRetile() {};

    void configure(const WireCell::Configuration& cfg) {
        NeedScope::configure(cfg);
        m_retiler = Factory::find_tn<IPCTreeMutate>(get<std::string>(cfg, "retiler", "RetileCluster"));
    }
    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["retiler"] = "RetileCluster";
        return cfg;
    }

    void visit(Ensemble& ensemble) const;

private:
    IPCTreeMutate::pointer m_retiler;

    // std::map<Cluster*, std::tuple<Cluster*, int, Cluster*>> process_groupings(
    //     Grouping& original, 
    //     Grouping& shadow, 
    //     const std::string& aname = "isolated", 
    //     const std::string& pname = "perblob") const;
        

};

void ClusteringRetile::visit(Ensemble& ensemble) const
{
    // fixme: make grouping names configurable

    auto& orig_grouping = *ensemble.with_name("live").at(0);
    auto& shad_grouping = ensemble.make_grouping("shadow");
    shad_grouping.from(orig_grouping);
    auto* shad_root = shad_grouping.node();

    for (auto* orig_cluster : orig_grouping.children()) {

        if (!orig_cluster->get_scope_filter(m_scope)) {
            // move on if the cluster is not in the scope filter ...
            continue;
        }

        if (orig_cluster->get_default_scope().hash() != m_scope.hash()) {
            orig_cluster->set_default_scope(m_scope);
        }

        auto shad_node = m_retiler->mutate(*orig_cluster->node());
        if (! shad_node) {
            continue;
        }
        shad_root->insert(std::move(shad_node));
    }
}


// std::map<WCF::Cluster*, std::tuple<WCF::Cluster*, int, WCF::Cluster*>> 
// ClusteringRetile::process_groupings(
//     WCF::Grouping& original,
//     WCF::Grouping& shadow,
//     const std::string& aname,
//     const std::string& pname) const
// {
//     return process_groupings_helper(original, shadow, aname, pname);
// }

