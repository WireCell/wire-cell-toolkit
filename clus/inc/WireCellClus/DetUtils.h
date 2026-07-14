/** This API is for various utility functions related to detector information */

#ifndef WIRECELL_CLUS_DETUTILS
#define WIRECELL_CLUS_DETUTILS

#include "WireCellClus/DynamicPointCloud.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellUtil/Point.h"
#include <map>
#include <vector>
#include <tuple>

namespace WireCell::Clus {

    /// Get all APA ident numbers known to the detector volumes.
    std::set<int> apa_idents(IDetectorVolumes::pointer dv);

    /// Get map from "layer" level WPID to a drift "dirx" vector and three wire angles.
    using wpid_faceparams_map = std::map<WirePlaneId, std::tuple<Vector, double, double, double>>;
    wpid_faceparams_map face_parameters(IDetectorVolumes::pointer dv);
    
    /// Create a new and empty DynamicPointCloud.  The DV is needed for the wire angles.
    std::shared_ptr<Facade::DynamicPointCloud> make_dynamicpointcloud(IDetectorVolumes::pointer dv);


}

#endif
