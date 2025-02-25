#ifndef WIRECELLUTIL_POINTSUMMARY
#define WIRECELLUTIL_POINTSUMMARY

#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/PointCloudArray.h"

#include "WireCellUtil/Configuration.h"

namespace WireCell::PointCloud {

    /// Summarize objects as JSON.  These are NOT full serializations but can be
    /// useful for debugging and validating.  Summary of points recurs to
    /// summary of datasets, datasets to arrays.  
    /// 
    /// If recur is true the summary recurs does DFS walk of the tree.
    Configuration json_summary(const Tree::Points& value, bool recur=true);
    Configuration json_summary(const Dataset& ds);
    Configuration json_summary(const Array& arr);
}

#endif
