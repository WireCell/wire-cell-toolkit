#include "WireCellClus/TrackFitting.h"
#include "WireCellUtil/Logging.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

TrackFitting::TrackFitting(FittingType fitting_type) 
    : m_fitting_type(fitting_type) 
{
}

bool TrackFitting::fit_track(Cluster& cluster) const {
    // if (!validate_cluster_for_fitting(cluster)) {
        // return false;
    // }

    bool success = false;
    switch (m_fitting_type) {
        case FittingType::Single:
            // success = perform_single_track_fitting(cluster);
            break;
        case FittingType::Multiple:
            // success = perform_multiple_track_fitting(cluster);
            break;
    }

    // update_cluster_after_fitting(cluster, success);
    return success;
}

size_t TrackFitting::fit_tracks(std::vector<Cluster*>& clusters) const {
    size_t successful_fits = 0;
    
    for (auto* cluster : clusters) {
        if (cluster && fit_track(*cluster)) {
            successful_fits++;
        }
    }
    
    return successful_fits;
}


