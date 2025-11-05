#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/IEnsembleVisitor.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PointTree.h"

#include <vector>

class ClusteringPointed;
WIRECELL_FACTORY(ClusteringPointed, ClusteringPointed,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;

class ClusteringPointed : public IConfigurable, public Clus::IEnsembleVisitor {
public:
    ClusteringPointed() {};
    virtual ~ClusteringPointed() {};

    void configure(const WireCell::Configuration& cfg) {
        auto jgroupings = cfg["groupings"];
        if (! jgroupings.isArray()) {
            return;
        }
        m_groupings.clear();
        for (const auto& one : jgroupings) {
            m_groupings.push_back(one.asString());
        }
    }
    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["groupings"][0] = "live";
        return cfg;
    }

    void visit(Ensemble& ensemble) const {

        for (const auto& name : m_groupings) {
            auto got = ensemble.with_name(name);
            if (got.empty()) { continue; }
            auto* grouping = got[0];

            std::vector<Cluster*> doomed_clusters;

            for (auto* cluster : grouping->children()) {

                std::vector<Blob*> doomed_blobs;

                for (auto* blob : cluster->children()) {
                    if (! blob->npoints()) {
                        doomed_blobs.push_back(blob);
                    }
                }
                
                for (auto* dead : doomed_blobs) {
                    cluster->destroy_child(dead);
                }

                if (! cluster->nchildren()) {
                    doomed_clusters.push_back(cluster);
                }
            }

            for (auto* dead : doomed_clusters) {
                grouping->destroy_child(dead);
            }

        }
    }

private:

    std::vector<std::string> m_groupings = {"live"};

};
