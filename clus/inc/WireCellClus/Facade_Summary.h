#ifndef WIRECELLCLUS_FACADE_SUMMARY
#define WIRECELLCLUS_FACADE_SUMMARY

#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"

#include "WireCellUtil/Configuration.h"

// fixme: this namespace here and elsewhere in clus should be in
// WireCell::Clus::Facade or something.
namespace WireCell::PointCloud::Facade {

    // Summarize facades as JSON.  These recur down the type hierarchy and into the tree via 
    Configuration json_summary(const Grouping& grp);
    Configuration json_summary(const Cluster& cls);
    Configuration json_summary(const Blob& blb);


}

#endif
