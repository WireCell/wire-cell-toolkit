
#include "WireCellImg/GeomClusteringUtil.h"

#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/RayClustering.h"

using namespace WireCell;
using namespace WireCell::Img;

// return the time between two blob sets relative to the slice span of
// the first.
static int rel_time_diff(const IBlobSet::pointer& one, const IBlobSet::pointer& two)
{
    const auto here = one->slice();
    return std::round((two->slice()->start() - here->start()) / here->span());
}

static bool adjacent(const std::set<double>& times, const double& time1, const double& time2)
{
    double tmin = time1;
    double tmax = time2;
    if (time1 > time2) std::swap(tmin, tmax);
    auto iter1 = times.find(tmin);
    if (iter1 == times.end()) THROW(ValueError() << errmsg{String::format("tmin \"%d\" not in \"times\"", tmin)});
    auto iter2 = times.find(tmax);
    if (iter2 == times.end()) THROW(ValueError() << errmsg{String::format("tmax \"%d\" not in \"times\"", tmin)});
    return (std::distance(iter1, iter2) == 1);
}

void WireCell::Img::geom_clustering(cluster_indexed_graph_t& grind, IBlobSet::vector::iterator beg,
                                    IBlobSet::vector::iterator end, std::string policy)
{
    auto next = beg + 1;

    if (next == end) {
        return;
    }

    std::unordered_set<std::string> known_policies = {"simple", "uboone", "uboone_local", "dead_clus"};
    if (known_policies.find(policy) == known_policies.end()) {
        THROW(ValueError() << errmsg{String::format("policy \"%s\" not implemented!", policy)});
    }

    // uboone policy
    int max_rel_diff = 2;
    std::map<int, RayGrid::grid_index_t> map_gap_tol = {{1, 2}, {2, 1}};

    if (policy == "uboone_local") {
        max_rel_diff = 2;
        map_gap_tol = {{1, 2}, {2, 2}};
    }

    if (policy == "simple") {
        max_rel_diff = 1;
        map_gap_tol = {{1, 0}};
    }

    // overlap + tolerance
    struct TolerantVisitor {
        RayGrid::grid_index_t tolerance{0};
        bool verbose{false};
        RayGrid::blobvec_t operator()(const RayGrid::blobref_t& blob, const RayGrid::blobproj_t& proj,
                                      RayGrid::layer_index_t layer)
        {
            return RayGrid::overlap(blob, proj, layer, tolerance, verbose);
        }
    };

    IBlob::vector iblobs2 = (*beg)->blobs();
    RayGrid::blobs_t blobs2 = (*beg)->shapes();
    const auto beg2 = blobs2.begin();
    for (auto test = next; test != end and rel_time_diff(*beg, *test) <= max_rel_diff; ++test) {
        int rel_diff = rel_time_diff(*beg, *test);
        if (map_gap_tol.find(rel_diff) == map_gap_tol.end()) continue;

        // handle each face separately faces
        IBlob::vector iblobs1 = (*test)->blobs();

        RayGrid::blobs_t blobs1 = (*test)->shapes();

        const auto beg1 = blobs1.begin();

        auto assoc = [&](RayGrid::blobref_t& a, RayGrid::blobref_t& b) {
            int an = a - beg1;
            int bn = b - beg2;
            grind.edge(iblobs1[an], iblobs2[bn]);
        };
        bool verbose = false;
        TolerantVisitor tv{map_gap_tol[rel_diff], verbose};
        RayGrid::associate(blobs1, blobs2, assoc, tv);
    }
}

// void WireCell::Img::geom_clustering(cluster_graph_t& cg, std::string policy, gc_filter_t filter)
// {
//     std::unordered_set<std::string> known_policies = {"simple", "uboone", "uboone_local"};
//     if (known_policies.find(policy) == known_policies.end()) {
//         THROW(ValueError() << errmsg{String::format("policy \"%s\" not implemented!", policy)});
//     }

//     // uboone policy
//     int max_rel_diff = 2;
//     std::map<int, RayGrid::grid_index_t> map_gap_tol = {{1, 2}, {2, 1}};

//     if (policy == "uboone_local") {
//         max_rel_diff = 2;
//         map_gap_tol = {{1, 2}, {2, 2}};
//     }

//     if (policy == "simple") {
//         max_rel_diff = 1;
//         map_gap_tol = {{1, 0}};
//     }

//     // overlap + tolerance
//     struct TolerantVisitor {
//         RayGrid::grid_index_t tolerance{0};
//         bool verbose{false};
//         RayGrid::blobvec_t operator()(const RayGrid::blobref_t& blob, const RayGrid::blobproj_t& proj,
//                                       RayGrid::layer_index_t layer)
//         {
//             return RayGrid::overlap(blob, proj, layer, tolerance, verbose);
//         }
//     };

//     /// collect all the slices
//     /// TODO: need sorting?
//     /// assumes all slices are unique
//     std::vector<cluster_vertex_t> slices;
//     std::unordered_map<cluster_vertex_t, std::vector<cluster_vertex_t> > s2b;
//     std::set<double> slice_times;
//     for (const auto& bvtx : GraphTools::mir(boost::vertices(cg))) {
//         if (cg[bvtx].code() == 'b') {
//             if (!filter(bvtx)) continue;
//             for (const auto& svtx : neighbors_oftype<cluster_node_t::slice_t>(cg, bvtx)) {
//                 slices.push_back(svtx);
//                 s2b[svtx].push_back(bvtx);
//                 const auto islice = get<cluster_node_t::slice_t>(cg[svtx].ptr);
//                 slice_times.insert(islice->start());
//             }
//         }
//     }

//     /// make SimpleBlobSet for slice pairs and make edges
//     for (auto siter1 = slices.begin(); siter1 != slices.end(); ++siter1) {
//         const auto islice1 = get<cluster_node_t::slice_t>(cg[*siter1].ptr);
//         auto& bdescs1 = s2b[*siter1];
//         auto gen = [&](const std::vector<cluster_vertex_t>& descs) {
//             RayGrid::blobs_t ret;
//             for (auto desc : descs) {
//                 const auto iblob = get<cluster_node_t::blob_t>(cg[desc].ptr);
//                 ret.push_back(iblob->shape());
//             }
//             return ret;
//         };
//         RayGrid::blobs_t blobs1 = gen(bdescs1);
//         const auto beg1 = blobs1.begin();
//         for (auto siter2 = siter1 + 1; siter2 != slices.end(); ++siter2) {
//             const auto islice2 = get<cluster_node_t::slice_t>(cg[*siter2].ptr);
//             int rel_diff = std::round(fabs((islice1->start() - islice2->start()) / islice1->span()));
//             if (map_gap_tol.find(rel_diff) == map_gap_tol.end()) continue;
//             if (policy == "uboone_local" && !adjacent(slice_times, islice1->start(), islice2->start())) continue;

//             // auto allbdescs1 = neighbors_oftype<cluster_node_t::blob_t>(cg, *siter1);
//             // auto allbdescs2 = neighbors_oftype<cluster_node_t::blob_t>(cg, *siter2);
//             // auto vfilter = [&](const std::vector<cluster_vertex_t>& descs) {
//             //     std::vector<cluster_vertex_t> ret;
//             //     for (auto desc : descs) {
//             //         if (!filter(desc)) continue;
//             //         ret.push_back(desc);
//             //     }
//             //     return ret;
//             // };
//             // auto bdescs1 = vfilter(allbdescs1);
//             // auto bdescs2 = vfilter(allbdescs2);

//             auto& bdescs2 = s2b[*siter2];
//             RayGrid::blobs_t blobs2 = gen(bdescs2);
//             const auto beg2 = blobs2.begin();

//             auto assoc = [&](RayGrid::blobref_t& a, RayGrid::blobref_t& b) {
//                 int an = a - beg1;
//                 int bn = b - beg2;
//                 boost::add_edge(bdescs1[an], bdescs2[bn], cg);
//             };
//             bool verbose = false;
//             TolerantVisitor tv{map_gap_tol[rel_diff], verbose};
//             RayGrid::associate(blobs1, blobs2, assoc, tv);
//         }
//     }
// }

void WireCell::Img::grouped_geom_clustering(cluster_graph_t& cg, std::string policy,
                                            const std::unordered_map<cluster_vertex_t, int> groups)
{
    std::unordered_set<std::string> known_policies = {"simple", "uboone", "uboone_local", "dead_clus"};
    if (known_policies.find(policy) == known_policies.end()) {
        THROW(ValueError() << errmsg{String::format("policy \"%s\" not implemented!", policy)});
    }

    // uboone policy
    // int max_rel_diff = 2;
    std::map<int, RayGrid::grid_index_t> map_gap_tol = {{1, 2}, {2, 1}};

    if (policy == "uboone_local") {
        // max_rel_diff = 2;
        map_gap_tol = {{1, 2}, {2, 2}};
    }

    if (policy == "simple") {
        // max_rel_diff = 1;
        map_gap_tol = {{1, 0}};
    }

    if (policy == "dead_clus") {
        // max_rel_diff = 1;
        map_gap_tol = {{0, 1}, {1, 1}};
    }

    // overlap + tolerance
    struct TolerantVisitor {
        RayGrid::grid_index_t tolerance{0};
        bool verbose{false};
        RayGrid::blobvec_t operator()(const RayGrid::blobref_t& blob, const RayGrid::blobproj_t& proj,
                                      RayGrid::layer_index_t layer)
        {
            return RayGrid::overlap(blob, proj, layer, tolerance, verbose);
        }
    };

    /// collect all the slices
    /// TODO: need sorting?
    /// assumes all slices are unique
    std::vector<cluster_vertex_t> slices;
    std::set<double> slice_times;
    for (const auto& svtx : GraphTools::mir(boost::vertices(cg))) {
        if (cg[svtx].code() == 's') {
            slices.push_back(svtx);
            const auto islice = get<cluster_node_t::slice_t>(cg[svtx].ptr);
            slice_times.insert(islice->start());
        }
    }

    /// make SimpleBlobSet for slice pairs and make edges
    for (auto siter1 = slices.begin(); siter1 != slices.end(); ++siter1) {
        auto& islice1 = get<cluster_node_t::slice_t>(cg[*siter1].ptr);
        auto bdescs1 = neighbors_oftype<cluster_node_t::blob_t>(cg, *siter1);
        auto gen = [&](const std::vector<cluster_vertex_t>& descs) {
            RayGrid::blobs_t ret;
            for (auto desc : descs) {
                const auto iblob = get<cluster_node_t::blob_t>(cg[desc].ptr);
                ret.push_back(iblob->shape());
            }
            return ret;
        };
        RayGrid::blobs_t blobs1 = gen(bdescs1);
        const auto beg1 = blobs1.begin();

        for (auto siter2 = siter1 + 1; siter2 != slices.end(); ++siter2) {
            auto& islice2 = get<cluster_node_t::slice_t>(cg[*siter2].ptr);
            int rel_diff = std::round(fabs((islice1->start() - islice2->start()) / islice1->span()));
            if (map_gap_tol.find(rel_diff) == map_gap_tol.end()) continue;
            if (policy == "uboone_local" && !adjacent(slice_times, islice1->start(), islice2->start())) continue;
            auto bdescs2 = neighbors_oftype<cluster_node_t::blob_t>(cg, *siter2);

            RayGrid::blobs_t blobs2 = gen(bdescs2);
            const auto beg2 = blobs2.begin();

            auto assoc = [&](RayGrid::blobref_t& a, RayGrid::blobref_t& b) {
                int an = a - beg1;
                int bn = b - beg2;
                if (groups.size() != 0) {
                    const auto& iter1 = groups.find(bdescs1[an]);
                    if (iter1 == groups.end()) return;
                    const auto& iter2 = groups.find(bdescs2[bn]);
                    if (iter2 == groups.end()) return;
                    if (iter1->second != iter2->second) return;
                }
                boost::add_edge(bdescs1[an], bdescs2[bn], cg);
            };
            bool verbose = false;
            TolerantVisitor tv{map_gap_tol[rel_diff], verbose};
            RayGrid::associate(blobs1, blobs2, assoc, tv);
        }
    }
}
