#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"

class ClusteringTaggerFlagTransfer;
WIRECELL_FACTORY(ClusteringTaggerFlagTransfer, ClusteringTaggerFlagTransfer,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

/**
 * Lightweight visitor that transfers tagger information from point clouds to cluster flags
 * This should run FIRST in the clustering pipeline, right after clusters are created
 */
class ClusteringTaggerFlagTransfer : public IConfigurable, public Clus::IEnsembleVisitor {
public:
    ClusteringTaggerFlagTransfer() {}
    virtual ~ClusteringTaggerFlagTransfer() {}

    virtual void configure(const WireCell::Configuration& config) {
        // No configuration needed - this is a pure transfer operation
    }
    
    virtual Configuration default_configuration() const {
        return Configuration{};
    }

    virtual void visit(Ensemble& ensemble) const {
        using spdlog::debug;
        
        // Process all groupings (live, dead, etc.)
        for (auto* grouping : ensemble.children()) {
            auto clusters = grouping->children();
            
            for (auto* cluster : clusters) {
                transfer_tagger_flags(*cluster);
            }
            
            debug("ClusteringTaggerFlagTransfer: Processed {} clusters in grouping '{}'", 
                  clusters.size(), grouping->get_name());
        }
    }

private:
    void transfer_tagger_flags(Cluster& cluster) const {
        const auto& lpc = cluster.value().local_pcs();
        auto it = lpc.find("tagger_info");
        if (it == lpc.end()) {
            return; // No tagger info available
        }
        
        const auto& tagger_pc = it->second;
        
        // Debug: List all available keys in tagger_info
        // std::cout << "Xin: Cluster " << cluster.ident() << " tagger_info keys: ";
        // for (const auto& key : tagger_pc.keys()) {
        //     std::cout << key << " ";
        // }
        // std::cout << std::endl;
        
        // Helper lambda to check if a flag should be set
        auto should_set_flag = [&](const std::string& flag_name) -> bool {
            auto arr = tagger_pc.get(flag_name);
            // std::cout << "Xin: " << cluster.ident() << " checking " << flag_name << " - arr ptr: " << arr;
            if (arr) {
                // std::cout << ", size: " << arr->size_major();
                if (arr->size_major() > 0) {
                    int value = arr->element<int>(0);
                    // std::cout << ", value: " << value;
                    // std::cout << std::endl;
                    return value > 0;
                }
            }
            // std::cout << " - NO DATA" << std::endl;
            return false;
        };
        
      
        
        // Set flags based on stored tagger information
        if (should_set_flag("has_beam_flash")) {
            cluster.set_flag(Flags::beam_flash);

            // Try individual flags first, fall back to array method
            // Only set other flags if beam flash is true (following prototype logic)
            if (should_set_flag("has_tgm") ) {
                cluster.set_flag(Flags::tgm);
            }
            if (should_set_flag("has_low_energy") ) {
                cluster.set_flag(Flags::low_energy);
            }
            if (should_set_flag("has_light_mismatch") ) {
                cluster.set_flag(Flags::light_mismatch);
            }
            if (should_set_flag("has_fully_contained") ) {
                cluster.set_flag(Flags::fully_contained);
            }
            if (should_set_flag("has_short_track_muon") ) {
                cluster.set_flag(Flags::short_track_muon);
            }
            if (should_set_flag("has_full_detector_dead") ) {
                cluster.set_flag(Flags::full_detector_dead);
            }
        }
        // If not beam flash coincident, no additional flags are set

        // if (cluster.get_flag(Flags::beam_flash))
        // std::cout << "Xin: " << cluster.ident() << " has beam_flash: "
        //           << cluster.get_flag(Flags::beam_flash) << " " << should_set_flag("has_beam_flash") << ", tgm: "
        //           << cluster.get_flag(Flags::tgm) << " " << should_set_flag("has_tgm") << ", low_energy: "
        //           << cluster.get_flag(Flags::low_energy) << " " << should_set_flag("has_low_energy") << ", light_mismatch: "
        //           << cluster.get_flag(Flags::light_mismatch) << " " << should_set_flag("has_light_mismatch") << ", fully_contained: "
        //           << cluster.get_flag(Flags::fully_contained) << " " << should_set_flag("has_fully_contained") << ", short_track_muon: "
        //           << cluster.get_flag(Flags::short_track_muon) << " " << should_set_flag("has_short_track_muon") << ", full_detector_dead: "
        //           << cluster.get_flag(Flags::full_detector_dead) << " " << should_set_flag("has_full_detector_dead") << std::endl;
    }
};

// Now in ClusteringFuncs.h, we can have simple flag checking utilities:

namespace WireCell::Clus::Facade {
    namespace TaggerUtils {
        
        // Simple flag checking using get_flag() accessor
        inline bool is_beam_flash(const Cluster& cluster) {
            return cluster.get_flag(Flags::beam_flash);
        }

        inline bool is_tgm(const Cluster& cluster) {
            return cluster.get_flag(Flags::tgm);
        }

        inline bool is_low_energy(const Cluster& cluster) {
            return cluster.get_flag(Flags::low_energy);
        }

        inline bool is_light_mismatch(const Cluster& cluster) {
            return cluster.get_flag(Flags::light_mismatch);
        }

        inline bool is_fully_contained(const Cluster& cluster) {
            return cluster.get_flag(Flags::fully_contained);
        }

        inline bool is_short_track_muon(const Cluster& cluster) {
            return cluster.get_flag(Flags::short_track_muon);
        }

        inline bool is_full_detector_dead(const Cluster& cluster) {
            return cluster.get_flag(Flags::full_detector_dead);
        }
        
        /**
         * Get event type from tagger metadata (still need point cloud for this)
         */
        inline int get_event_type(const Cluster& cluster) {
            const auto& lpc = cluster.value().local_pcs();
            auto it = lpc.find("tagger_info");
            if (it == lpc.end()) return -1;
            
            auto arr = it->second.get("event_type");
            return (arr && arr->size_major() > 0) ? arr->element<int>(0) : -1;
        }
        
        /**
         * Get cluster length from tagger metadata
         */
        inline double get_cluster_length(const Cluster& cluster) {
            const auto& lpc = cluster.value().local_pcs();
            auto it = lpc.find("tagger_info");
            if (it == lpc.end()) return -1.0;
            
            auto arr = it->second.get("cluster_length");
            return (arr && arr->size_major() > 0) ? arr->element<double>(0) : -1.0;
        }
        
        /**
         * Check if cluster has any tagger flags set
         */
        inline bool has_any_tagger_flags(const Cluster& cluster) {
            return is_beam_flash(cluster) || is_tgm(cluster) || is_low_energy(cluster) ||
                   is_light_mismatch(cluster) || is_fully_contained(cluster) ||
                   is_short_track_muon(cluster) || is_full_detector_dead(cluster);
        }
    }
}