/// This API provides some mixin classes for "clustering classes" to handle
/// common behavior.

#ifndef WIRECELLCLUS_CLUSTERINGFUNCSMIXINS
#define WIRECELLCLUS_CLUSTERINGFUNCSMIXINS

#include "WireCellClus/IPCTransform.h"

#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellUtil/Configuration.h"

#include "WireCellUtil/PointTree.h"

namespace WireCell::Clus {

    // A mixin to get an IDetectorVolumes
    // 
    // Configuration:
    //
    // "detector_volumes" : string type/name of an IDetectorVolumes
    struct NeedDV {
        IDetectorVolumes::pointer m_dv;
        void configure(const WireCell::Configuration &cfg);
    };

    // A mixin to get an IPCTransformSet
    // 
    // Configuration:
    //
    // "pc_transforms" : string type/name of an IPCTransformSet
    struct NeedPCTS {
        Clus::IPCTransformSet::pointer m_pcts;
        void configure(const WireCell::Configuration &cfg);
    };

    // A mixin for things that need to be configured for a scope (pc name and coord names)
    // 
    // Configuration:
    //
    // "pc_name" : string name for a "node-local" PointCloud::Dataset
    // "coords" : array of string giving names of arrays in the PC to use as coordinates
    struct NeedScope {
        PointCloud::Tree::Scope m_scope;
        NeedScope(const std::string &pcname = "3d", 
                  const std::vector<std::string> &coords = {"x", "y", "z"},
                  size_t depth = 0);
        void configure(const WireCell::Configuration &cfg);
    };

}
#endif
