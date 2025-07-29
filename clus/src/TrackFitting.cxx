#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"

class TrackFitting;
WIRECELL_FACTORY(TrackFitting, TrackFitting,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

/**
 * Track fitting function that processes clusters marked with the do_tracking flag.
 * Supports both single and multiple track fitting modes based on configuration.
 */
class TrackFitting : public IConfigurable, public Clus::IEnsembleVisitor {
public:
    TrackFitting() {}
    virtual ~TrackFitting() {}

    virtual void configure(const WireCell::Configuration& config) {
        m_grouping_name = get<std::string>(config, "grouping", "live");
        m_fitting_type = get<std::string>(config, "fitting_type", "single");
        
        // Validate fitting_type parameter
        if (m_fitting_type != "single" && m_fitting_type != "multiple") {
            std::cerr << "TrackFitting: Invalid fitting_type '" << m_fitting_type 
                      << "'. Must be 'single' or 'multiple'. Using 'single' as default." << std::endl;
            m_fitting_type = "single";
        }
    }
        
    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping_name;
        cfg["fitting_type"] = m_fitting_type;
        return cfg;
    }

    virtual void visit(Ensemble& ensemble) const {
        
        // Get the specified grouping (default: "live")
        auto groupings = ensemble.with_name(m_grouping_name);
        if (groupings.empty()) {
            return;
        }
        
        auto& grouping = *groupings.at(0);
        
        // Find clusters that have the do_tracking flag
        std::vector<Cluster*> tracking_clusters;
        
        for (auto* cluster : grouping.children()) {
            if (cluster->get_flag(Flags::do_tracking)) {
                tracking_clusters.push_back(cluster);
            }
        }

        std::cout << "TrackFitting: Found " << tracking_clusters.size()
                  << " clusters marked for tracking (fitting_type: " << m_fitting_type << ")." << std::endl;

        // Process each cluster marked for tracking
        for (auto* cluster : tracking_clusters) {
            perform_track_fitting(*cluster);
        }
    }

private:
    std::string m_grouping_name{"live"};
    std::string m_fitting_type{"single"};
    
    /**
     * Perform track fitting on a cluster marked with do_tracking flag.
     * Implementation depends on fitting_type configuration.
     * 
     * @param cluster The cluster to perform track fitting on
     */
    void perform_track_fitting(Cluster& cluster) const {
        if (m_fitting_type == "single") {
            // TODO: Implement single track fitting algorithm
            std::cout << "TrackFitting: Performing single track fitting on cluster " 
                      << cluster.ident() << std::endl;
            
        } else if (m_fitting_type == "multiple") {
            // TODO: Implement multiple track fitting algorithm
            std::cout << "TrackFitting: Performing multiple track fitting on cluster " 
                      << cluster.ident() << std::endl;
        }
        
        // TODO: Set appropriate flags or store results based on fitting outcome
        // Example: cluster.set_flag(Flags::fitted_track);
    }
};