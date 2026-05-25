#ifndef WIRECELL_MATCH_UTIL
#define WIRECELL_MATCH_UTIL

#include "WireCellMatch/TimingTPCBundle.h"
#include "WireCellIface/ITensorSet.h"
#include "WireCellUtil/PointTree.h"

#include <map>
#include <string>

namespace WireCell::Match {

    void dump_bee_3d(const PointCloud::Tree::Points::node_t& root,
                     const std::string& fn);

    void dump_bee_flash(const ITensorSet::pointer& ts,
                        const std::string& fn);

    void dump_bee_bundle(const FlashBundlesMap& f2bundle,
                         const std::map<WireCell::Clus::Facade::Cluster*, int>& cluster_idx_map,
                         const std::string& fn);

} // namespace WireCell::Match

#endif
