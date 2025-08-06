#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"

class TaggerCheckSTM;
WIRECELL_FACTORY(TaggerCheckSTM, TaggerCheckSTM,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

/**
 * Clustering function that checks the main cluster from clustering_recovering_bundle
 * for Short Track Muon (STM) characteristics and sets the STM flag when conditions are met.
 * This function works on clusters that have already been processed by clustering_recovering_bundle.
 */
class TaggerCheckSTM : public IConfigurable, public Clus::IEnsembleVisitor, private Clus::NeedDV {
public:
    TaggerCheckSTM() {}
    virtual ~TaggerCheckSTM() {}

    virtual void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        m_grouping_name = get<std::string>(config, "grouping", "live");
    }
    
    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping_name;
        cfg["detector_volumes"] = "DetectorVolumes";
        return cfg;
    }

    virtual void visit(Ensemble& ensemble) const {
        
        // Get the specified grouping (default: "live")
        auto groupings = ensemble.with_name(m_grouping_name);
        if (groupings.empty()) {
            return;
        }
        
        auto& grouping = *groupings.at(0);
        
        // Find clusters that have the main_cluster flag (set by clustering_recovering_bundle)
        std::vector<Cluster*> main_clusters;

        for (auto* cluster : grouping.children()) {
            if (cluster->get_flag(Flags::main_cluster)) {
                main_clusters.push_back(cluster);
            }
        }

        std::cout << "TaggerCheckSTM: Found " << main_clusters.size() 
                  << " main clusters to check for STM conditions." << std::endl;

        // Process each main cluster
        size_t stm_count = 0;
        for (auto* cluster : main_clusters) {
            if (check_stm_conditions(*cluster)) {
                cluster->set_flag(Flags::STM);
                stm_count++;
            }
        }
        
        (void)stm_count;
    }

private:
    std::string m_grouping_name{"live"};

   
    /**
     * Check if a cluster meets the conditions for STM (Short Track Muon) tagging.
     * This is where you'll implement your specific STM detection algorithm.
     * 
     * @param cluster The main cluster to analyze
     * @return true if cluster should be flagged as STM
     */
    bool check_stm_conditions(const Cluster& cluster) const {
        // get all the angles ...

        // Get all the wire plane IDs from the grouping
        const auto& wpids = cluster.grouping()->wpids();
        // Key: pair<APA, face>, Value: drift_dir, angle_u, angle_v, angle_w
        std::map<WirePlaneId , std::tuple<geo_point_t, double, double, double>> wpid_params;
        std::map<WirePlaneId, std::pair<geo_point_t, double> > wpid_U_dir;
        std::map<WirePlaneId, std::pair<geo_point_t, double> > wpid_V_dir;
        std::map<WirePlaneId, std::pair<geo_point_t, double> > wpid_W_dir;
        std::set<int> apas;
        compute_wireplane_params(wpids, m_dv, wpid_params, wpid_U_dir, wpid_V_dir, wpid_W_dir, apas);

        std::cout << "TaggerCheckSTM: Checking cluster with " << wpids.size() 
                  << " wire plane IDs and " << apas.size() << " APAs." << std::endl;

        

        return false;
    }
    
 
};