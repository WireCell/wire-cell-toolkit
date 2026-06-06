#include "WireCellImg/Projection2D.h"
#include "WireCellUtil/Exceptions.h"

#include <boost/graph/filtered_graph.hpp>
#include "WireCellUtil/Array.h"
#include "WireCellUtil/Stream.h"
#include "WireCellUtil/String.h"

#include <fstream>
#include <limits>
#include <unordered_set>
#include <utility>
#include <algorithm>

using namespace WireCell;
using namespace WireCell::Img;
using namespace WireCell::Img::Projection2D;

// Here, "node" implies a cluster graph vertex payload object.
using channel_t = cluster_node_t::channel_t;
using wire_t = cluster_node_t::wire_t;
using blob_t = cluster_node_t::blob_t;
using slice_t = cluster_node_t::slice_t;
using meas_t = cluster_node_t::meas_t;

// template <typename It> boost::iterator_range<It> mir(std::pair<It, It> const& p) {
//     return boost::make_iterator_range(p.first, p.second);
// }
// template <typename It> boost::iterator_range<It> mir(It b, It e) {
//     return boost::make_iterator_range(b, e);
// }

// std::vector<cluster_vertex_t> WireCell::Img::Projection2D::neighbors(const WireCell::cluster_graph_t& cg, const
// cluster_vertex_t& vd)
// {
//     std::vector<cluster_vertex_t> ret;
//     for (auto edge : boost::make_iterator_range(boost::out_edges(vd, cg))) {
//         cluster_vertex_t neigh = boost::target(edge, cg);
//         ret.push_back(neigh);
//     }
//     return ret;
// }

// template <typename Type>
// std::vector<cluster_vertex_t> WireCell::Img::Projection2D::neighbors_oftype(const WireCell::cluster_graph_t& cg,
// const cluster_vertex_t& vd)
// {
//     std::vector<cluster_vertex_t> ret;
//     for (const auto& vp : neighbors(cg, vd)) {
//         if (std::holds_alternative<Type>(cg[vp].ptr)) {
//             ret.push_back(vp);
//         }
//     }
//     return ret;
// }

std::unordered_map<int, std::set<cluster_vertex_t> > WireCell::Img::Projection2D::get_geom_clusters(
    const WireCell::cluster_graph_t& cg)
{
    std::unordered_map<int, std::set<cluster_vertex_t> > groups;

    // Blob-only subgraph as a lazy view: keep 'b' vertices (a b-b edge survives in a
    // vertex-filtered graph only when both of its endpoints do, so this is exactly the
    // blob + b-b-edge subgraph).  Was a materialized copy (add_vertex/add_edge into a
    // fresh cg_blob) + connected_components; the filtered_graph yields identical
    // components with no graph allocation, and the descriptors are already the
    // original cg vertices so no new2old remap is needed.
    using BFiltered =
        boost::filtered_graph<cluster_graph_t, boost::keep_all, std::function<bool(cluster_vertex_t)> >;
    BFiltered bcg(cg, {}, [&](auto vtx) { return cg[vtx].code() == 'b'; });

    std::unordered_map<cluster_vertex_t, int> desc2id;
    boost::connected_components(bcg, boost::make_assoc_property_map(desc2id));
    for (auto& [desc, id] : desc2id) {  // invert: component id -> original blob descriptors
        groups[id].insert(desc);
    }

    return groups;
}

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(std::pair<T1, T2> const& pair) const
    {
        std::size_t h1 = std::hash<T1>()(pair.first);
        std::size_t h2 = std::hash<T2>()(pair.second);
        // boost::hash_combine style mixing to avoid h1^h2 symmetry collisions
        h1 ^= h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2);
        return h1;
    }
};

LayerProjection2DMap WireCell::Img::Projection2D::get_projection(const WireCell::cluster_graph_t& cg,
                                                                 const std::set<cluster_vertex_t>& group,
                                                                 const size_t nchan, const size_t nslice,
                                                                 double uncer_cut, double dead_default_charge)
{
    using triplet_t = Eigen::Triplet<scaler_t>;
    using triplet_vec_t = std::vector<triplet_t>;
    std::unordered_map<WirePlaneLayer_t, triplet_vec_t> lcoeff;
    std::unordered_set<std::pair<int, int>, pair_hash> filled;
    //    std::unordered_set<int> filled_slices;

    // record the number of slices for each layer ...
    std::unordered_map<WirePlaneLayer_t, std::unordered_set<int> > layer_nslices;

    // layer_projection_map_t ret;
    LayerProjection2DMap ret;

    // blob_est_min_charge = min(u,v,w)
    // cluster_est_min_charge = sum(blob_est_min_charge )
    double estimated_minimum_charge = 0;
    double estimated_total_charge = 0;
    //    int saved_flag = 0;
    // int saved_flag_1 = 0;
    int number_blobs = 0;
    //    int number_slices = 0;

    // assumes one blob linked to one slice
    // use b-w-c to find all channels linked to the blob
    std::unordered_map<slice_t, std::vector<cluster_vertex_t> > map_s2vb;
    std::unordered_map<cluster_vertex_t, std::vector<cluster_vertex_t> > map_b2c;
    for (const auto& bdesc : group) {
        const auto& node = cg[bdesc];
        if (node.code() == 'b') {
            number_blobs++;
            for (auto edge : boost::make_iterator_range(boost::out_edges(bdesc, cg))) {
                cluster_vertex_t neigh = boost::target(edge, cg);
                /// ASSUMPTION: only 1 b-s for each b
                if (cg[neigh].code() == 's') {
                    auto& slice = std::get<slice_t>(cg[neigh].ptr);
                    map_s2vb[slice].push_back(bdesc);
                }
                /// ASSUMPTION: only 1 w-c for each w
                if (cg[neigh].code() == 'w') {
                    for (const auto& wedge : boost::make_iterator_range(boost::out_edges(neigh, cg))) {
                        cluster_vertex_t cdesc = boost::target(wedge, cg);
                        if (cg[cdesc].code() == 'c') {
                            map_b2c[bdesc].push_back(cdesc);
                            break;
                        }
                    }
                }
            }
        }
    }

    for (const auto& [slice, bdescs] : map_s2vb) {
        int start = slice->start() / slice->span();
        auto activity = slice->activity();

        for (const auto& bdesc : bdescs) {
            std::unordered_map<WirePlaneLayer_t, double> layer_charge;
            // initialization ...
            layer_charge[kUlayer] = 0;
            layer_charge[kVlayer] = 0;
            layer_charge[kWlayer] = 0;

            for (const auto& chan_desc : map_b2c[bdesc]) {
                const auto& chan = std::get<channel_t>(cg[chan_desc].ptr);
                WirePlaneLayer_t layer = chan->planeid().layer();
                int cident = chan->index();
                auto charge = activity[chan].value();
                auto unc = activity[chan].uncertainty();
                // TODO: make this configurable and robust
                if (unc > uncer_cut) {
                    charge = dead_default_charge;
                }
                else {
                    // TODO: double check this
                    layer_charge[layer] += charge;
                    layer_nslices[layer].insert(start);
                }
                // if filled, skip
                if (filled.find({cident, start}) != filled.end()) {
                    continue;
                }
                filled.insert({cident, start});
                // TODO: validate this
                lcoeff[layer].push_back({cident, start, charge});
            }  // loop over channel

            double sum_charge = 0;
            int sum_n = 0;
            double min_charge = std::numeric_limits<double>::max();
            for (auto it = layer_charge.begin(); it != layer_charge.end(); it++) {
                if (it->second != 0) {
                    sum_charge += it->second;
                    sum_n++;
                    if (it->second < min_charge) min_charge = it->second;
                }
            }
            // protection: no layer had non-zero charge
            if (sum_n == 0) min_charge = 0;
            if (sum_n > 0) estimated_total_charge = sum_charge / sum_n;
            estimated_minimum_charge += min_charge;  // min_iter->second;

        }  // loop over blobs
    }      // loop over slices

    for (auto it = layer_nslices.begin(); it != layer_nslices.end(); it++) {
        ret.m_number_layer_slices[it->first] = it->second.size();
    }

    ret.m_layer_proj.insert({kUlayer, Projection2D(nchan, nslice)});
    ret.m_layer_proj.insert({kVlayer, Projection2D(nchan, nslice)});
    ret.m_layer_proj.insert({kWlayer, Projection2D(nchan, nslice)});

    for (auto lc : lcoeff) {
        auto& l = lc.first;
        auto& c = lc.second;
        // if (ret.m_layer_proj.find(l) == ret.m_layer_proj.end()) {
        //     ret.m_layer_proj.insert({l, Projection2D(nchan, nslice)});
        // }
        ret.m_layer_proj[l].m_proj.setFromTriplets(c.begin(), c.end());
    }

    ret.m_estimated_minimum_charge = estimated_minimum_charge;
    ret.m_estimated_total_charge = estimated_total_charge;
    ret.m_number_blobs = number_blobs;
    ret.m_number_slices = map_s2vb.size();  // filled_slices.size();
    return ret;
}

std::string WireCell::Img::Projection2D::dump(const Projection2D& proj2d, bool verbose)
{
    std::stringstream ss;
    ss << "Projection2D:";

    auto ctq = proj2d.m_proj;
    size_t counter = 0;
    for (int k = 0; k < ctq.outerSize(); ++k) {
        for (sparse_mat_t::InnerIterator it(ctq, k); it; ++it) {
            if (verbose) {
                ss << " {" << it.row() << ", " << it.col() << "}->" << it.value();
            }
            ++counter;
        }
    }
    ss << " #: " << counter;

    return ss.str();
}

bool WireCell::Img::Projection2D::write(const Projection2D& proj2d, const std::string& fname)
{
    using namespace WireCell::Stream;
    boost::iostreams::filtering_ostream fout;
    fout.clear();
    output_filters(fout, fname);
    Eigen::MatrixXf dense_m(proj2d.m_proj);
    Array::array_xxf dense_a = dense_m.array();
    const std::string aname = String::format("proj2d_%d.npy", 0);
    WireCell::Stream::write(fout, aname, dense_a);
    fout.flush();
    fout.pop();
    return 0;
}

// some local utility functions
namespace {

    // return a sparse matrix mask for elements evaluate true for function f
    sparse_mat_t mask(const sparse_mat_t& sm, std::function<bool(scaler_t)> f)
    {
        sparse_mat_t smm = sm;
        for (int k = 0; k < smm.outerSize(); ++k) {
            for (sparse_mat_t::InnerIterator it(smm, k); it; ++it) {
                if (f(it.value())) {
                    it.valueRef() = 1;
                }
                else {
                    it.valueRef() = 0;
                }
            }
        }
        return smm;
    }

}  // namespace

// 1: tar is part of ref: REF_COVERS_TAR
// -1: ref is part of tar: TAR_COVERAS_REF
// 2: tar is equal to ref: REF_EQ_TAR
// -2: tar is equal to ref and both are empty: BOTH EMPTY ...
// 0 ref and tar do not belong to each other:   OTHER
WireCell::Img::Projection2D::Coverage WireCell::Img::Projection2D::judge_coverage(const Projection2D& ref,
                                                                                  const Projection2D& tar,
                                                                                  double uncer_cut)
{
    // Single union-pass over the two projections' stored entries, replacing the
    // previous mask(ref)+mask(tar)+(ref_mask-tar_mask) (three sparse allocations).
    // Bit-identical: mask() sets every STORED entry to 1 (live: value > -uncer_cut)
    // or 0, keeping the same sparsity, so ref_mask-tar_mask is in {-1,0,1} over the
    // UNION of stored positions, and the four booleans reduce exactly to:
    //   ref_non_zero  = exists a position live in ref
    //   tar_non_zero  = exists a position live in tar
    //   ref_m_tar_pos = exists a position live in ref but NOT live in tar  (diff +1)
    //   ref_m_tar_neg = exists a position live in tar but NOT live in ref  (diff -1)
    // where "live in M" == stored in M.m_proj with value > -uncer_cut, and
    // "not live" includes both not-stored and stored-but-<=-uncer_cut. Both
    // projections share dims (nchan x nslice) and are column-major, so a per-column
    // merge of the (row-sorted) InnerIterators visits each union position once.
    bool ref_non_zero = false, tar_non_zero = false;
    bool ref_m_tar_pos = false, ref_m_tar_neg = false;
    const int ncols = ref.m_proj.outerSize();
    for (int k = 0; k < ncols; ++k) {
        sparse_mat_t::InnerIterator rit(ref.m_proj, k);
        sparse_mat_t::InnerIterator tit(tar.m_proj, k);
        while (rit || tit) {
            const int rrow = rit ? (int) rit.row() : std::numeric_limits<int>::max();
            const int trow = tit ? (int) tit.row() : std::numeric_limits<int>::max();
            bool r_live = false, t_live = false;
            if (rrow == trow) {
                r_live = rit.value() > -uncer_cut;
                t_live = tit.value() > -uncer_cut;
                ++rit;
                ++tit;
            }
            else if (rrow < trow) {
                r_live = rit.value() > -uncer_cut;
                ++rit;
            }
            else {
                t_live = tit.value() > -uncer_cut;
                ++tit;
            }
            ref_non_zero  = ref_non_zero  || r_live;
            tar_non_zero  = tar_non_zero  || t_live;
            ref_m_tar_pos = ref_m_tar_pos || (r_live && !t_live);
            ref_m_tar_neg = ref_m_tar_neg || (t_live && !r_live);
            // all flags saturated: both non-empty with exclusive live cells -> OTHER,
            // and no later cell can change that.
            if (ref_non_zero && tar_non_zero && ref_m_tar_pos && ref_m_tar_neg) {
                return OTHER;
            }
        }
    }

    if ((!ref_non_zero) && (!tar_non_zero)) {
        return BOTH_EMPTY;
    }
    else if (ref_non_zero && (!tar_non_zero)) {
        return REF_COVERS_TAR;
    }
    else if ((!ref_non_zero) && tar_non_zero) {
        return TAR_COVERS_REF;
    }
    else if ((!ref_m_tar_neg) && (!ref_m_tar_pos)) {
        return REF_EQ_TAR;
    }
    else if (ref_m_tar_neg && (!ref_m_tar_pos)) {
        return TAR_COVERS_REF;
    }
    else if ((!ref_m_tar_neg) && ref_m_tar_pos) {
        return REF_COVERS_TAR;
    }
    else {
        return OTHER;
    }

    // ref * tar
    // sparse_mat_t ref_t_tar = ref.m_proj.cwiseProduct(tar.m_proj);

    // non overlapping
    // bool all_zero = !loop_exist(ref_t_tar, [](scaler_t x){return (x!=0 && x > -1e8);});
    //  if (all_zero) return OTHER;

    // partial overlapping
    // if (ref_m_tar_neg && ref_m_tar_pos) {
    //    return OTHER;
    // }

    // no non-overlapping pixels
    // TODO: all dead is also non-zero?

    // if (!ref_m_tar_neg && !ref_m_tar_pos) {
    //      if (ref_non_zero) {
    //          return REF_EQ_TAR;
    //      }
    //      return BOTH_EMPTY;
    //  }

    // ref > tar
    //  if (ref_m_tar_pos) {
    //      return REF_COVERS_TAR;
    //  }
}

// 1: tar is part of ref
// -1: ref is part of tar
// 2: tar is equal to ref
// 0 ref and tar do not belong to each other
WireCell::Img::Projection2D::Coverage WireCell::Img::Projection2D::judge_coverage_alt(const Projection2D& ref,
                                                                                      const Projection2D& tar,
                                                                                      std::vector<double>& cut_values,
                                                                                      double uncer_cut)
{
    // ref * tar
    //    sparse_mat_t ref_t_tar = ref.m_proj.cwiseProduct(tar.m_proj);
    // sparse_mat_t inter_mask = mask(ref_t_tar, [](scaler_t x){return x>0;});
    /// sparse_mat_t inter_proj = ref.m_proj.cwiseProduct(inter_mask);

    sparse_mat_t ref_mask = mask(ref.m_proj, [](scaler_t x) { return x > 0; });
    sparse_mat_t tar_mask = mask(tar.m_proj, [](scaler_t x) { return x > 0; });
    sparse_mat_t inter_mask = ref_mask.cwiseProduct(tar_mask);
    sparse_mat_t inter_proj = ref.m_proj.cwiseProduct(inter_mask);

    // Single-pass stats for ref: count live, count dead, sum charge
    int num_ref = 0, num_dead_ref = 0;
    scaler_t charge_ref = 0;
    for (int k = 0; k < ref.m_proj.outerSize(); ++k) {
        for (sparse_mat_t::InnerIterator it(ref.m_proj, k); it; ++it) {
            if (it.value() > 0) { ++num_ref; charge_ref += it.value(); }
            else if (it.value() < (-1) * uncer_cut) { ++num_dead_ref; }
        }
    }
    // Single-pass stats for tar: count live, count dead, sum charge
    int num_tar = 0, num_dead_tar = 0;
    scaler_t charge_tar = 0;
    for (int k = 0; k < tar.m_proj.outerSize(); ++k) {
        for (sparse_mat_t::InnerIterator it(tar.m_proj, k); it; ++it) {
            if (it.value() > 0) { ++num_tar; charge_tar += it.value(); }
            else if (it.value() < (-1) * uncer_cut) { ++num_dead_tar; }
        }
    }
    // Single-pass stats for inter_proj: count and charge
    int num_inter = 0;
    scaler_t charge_inter = 0;
    for (int k = 0; k < inter_proj.outerSize(); ++k) {
        for (sparse_mat_t::InnerIterator it(inter_proj, k); it; ++it) {
            if (it.value() > 0) { ++num_inter; charge_inter += it.value(); }
        }
    }

    if (num_ref != 0 && num_tar != 0 && num_inter == 0) return OTHER;

    //    if (num_ref == 0 && num_tar == 0) {
    //    return BOTH_EMPTY;
    // }
    // else if (num_ref == num_tar && num_dead_ref == num_dead_tar) {
    //    return REF_EQ_TAR;
    // }
    if (num_inter == num_ref) {
        return TAR_COVERS_REF;
    }
    else if (num_inter == num_tar) {
        return REF_COVERS_TAR;
    }
    else {
        Coverage value;

        float small_counts;
        float small_charge;
        float dead_counts;
        if (num_ref < num_tar) {
            value = TAR_COVERS_REF;
            small_counts = num_ref;
            small_charge = charge_ref;
            dead_counts = num_dead_ref;
        }
        else {
            value = REF_COVERS_TAR;
            small_counts = num_tar;
            small_charge = charge_tar;
            dead_counts = num_dead_tar;
        }

        float common_counts = num_inter;
        float common_charge = charge_inter;

        if (small_charge == 0 || small_counts == 0) return OTHER;

        if ((1 - common_charge / small_charge) <
                std::min(cut_values[0] * (small_counts + dead_counts) / small_counts, cut_values[1]) &&
            (1 - common_counts / small_counts) <
                std::min(cut_values[2] * (small_counts + dead_counts) / small_counts, cut_values[3])) {
            return value;
        }
        else {
            return OTHER;
        }
    }
}
