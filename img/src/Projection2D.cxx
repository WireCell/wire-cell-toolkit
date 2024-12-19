#include "WireCellImg/Projection2D.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/Array.h"
#include "WireCellUtil/Stream.h"
#include "WireCellUtil/String.h"

#include <fstream>
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
    cluster_graph_t cg_blob;

    size_t nblobs = 0;
    std::unordered_map<cluster_vertex_t, cluster_vertex_t> old2new;
    std::unordered_map<cluster_vertex_t, cluster_vertex_t> new2old;
    for (const auto& vtx : GraphTools::mir(boost::vertices(cg))) {
        const auto& node = cg[vtx];
        if (node.code() == 'b') {
            ++nblobs;
            auto newvtx = boost::add_vertex(node, cg_blob);
            old2new[vtx] = newvtx;
            new2old[newvtx] = vtx;
        }
    }
    // debug
    // std::cout << "nblobs: " << nblobs << std::endl;

    if (!nblobs) {
        return groups;
    }

    for (auto edge : GraphTools::mir(boost::edges(cg))) {
        auto old_tail = boost::source(edge, cg);
        auto old_head = boost::target(edge, cg);

        auto old_tit = old2new.find(old_tail);
        if (old_tit == old2new.end()) {
            continue;
        }
        auto old_hit = old2new.find(old_head);
        if (old_hit == old2new.end()) {
            continue;
        }
        const auto& hnode = cg_blob[old_hit->second];
        const auto& tnode = cg_blob[old_tit->second];
        if (hnode.code() == 'b' && tnode.code() == 'b') {
            boost::add_edge(old_tit->second, old_hit->second, cg_blob);
        }
    }

    // for debugging, return as one group
    // for (const auto& desc : GraphTools::mir(boost::vertices(cg_blob))) {
    //     groups[0].insert(new2old[desc]);
    // }
    // return groups;

    std::unordered_map<cluster_vertex_t, int> desc2id;
    boost::connected_components(cg_blob, boost::make_assoc_property_map(desc2id));
    for (auto& [desc, id] : desc2id) {  // invert
        groups[id].insert(new2old[desc]);
    }

    // debug
    // for (auto id2desc : groups) {
    //     auto & descvec = id2desc.second;
    //     if (descvec.size()>10) {
    //         std::cout << id2desc.first << ": ";
    //         for (auto &desc : descvec) {
    //             std::cout << desc << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    // }

    return groups;
}

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(std::pair<T1, T2> const& pair) const
    {
        std::size_t h1 = std::hash<T1>()(pair.first);
        std::size_t h2 = std::hash<T2>()(pair.second);
        return h1 ^ h2;
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
            double min_charge = 1e12;
            for (auto it = layer_charge.begin(); it != layer_charge.end(); it++) {
                if (it->second != 0) {
                    sum_charge += it->second;
                    sum_n++;
                    if (it->second < min_charge) min_charge = it->second;
                }
            }
            // protection ...
            if (min_charge == 1e12) min_charge = 0;
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

    bool loop_exist(const sparse_mat_t& sm, std::function<bool(scaler_t)> f)
    {
        for (int k = 0; k < sm.outerSize(); ++k) {
            for (sparse_mat_t::InnerIterator it(sm, k); it; ++it) {
                if (f(it.value())) {
                    return true;
                }
            }
        }
        return false;
    }

    int loop_count(const sparse_mat_t& sm, std::function<bool(scaler_t)> f)
    {
        int ret = 0;
        for (int k = 0; k < sm.outerSize(); ++k) {
            for (sparse_mat_t::InnerIterator it(sm, k); it; ++it) {
                if (f(it.value())) {
                    ret += 1;
                }
            }
        }
        return ret;
    }

    scaler_t loop_sum(const sparse_mat_t& sm, std::function<scaler_t(scaler_t)> f)
    {
        scaler_t ret = 0;
        for (int k = 0; k < sm.outerSize(); ++k) {
            for (sparse_mat_t::InnerIterator it(sm, k); it; ++it) {
                ret += f(it.value());
            }
        }
        return ret;
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
    sparse_mat_t ref_mask = mask(ref.m_proj, [&](scaler_t x) { return x > -uncer_cut; });
    sparse_mat_t tar_mask = mask(tar.m_proj, [&](scaler_t x) { return x > -uncer_cut; });

    // non zero component in both
    // bool ref_non_zero = loop_exist(ref.m_proj, [&](scaler_t x){return (x!=0 && x > -uncer_cut);}); // number of live
    // elements bool tar_non_zero = loop_exist(tar.m_proj, [&](scaler_t x){return (x!=0 && x > -uncer_cut);}); // number
    // of live elements ...

    bool ref_non_zero = loop_exist(ref_mask, [](scaler_t x) { return (x != 0); });  // number of live elements
    bool tar_non_zero = loop_exist(tar_mask, [](scaler_t x) { return (x != 0); });  // number of live elements ...

    if ((!ref_non_zero) && (!tar_non_zero)) {
        return BOTH_EMPTY;
    }
    else if (ref_non_zero && (!tar_non_zero)) {
        return REF_COVERS_TAR;
    }
    else if ((!ref_non_zero) && tar_non_zero) {
        return TAR_COVERS_REF;
    }
    else {
        // ref - tar
        //    sparse_mat_t ref_m_tar = ref.m_proj - tar.m_proj;
        sparse_mat_t ref_m_tar = ref_mask - tar_mask;

        // bool ref_m_tar_neg = loop_exist(ref_m_tar, [&](scaler_t x){return (x<0 && fabs(x)<fabs(uncer_cut));}); //
        // number of live elements belong to tar, not belong to ref bool ref_m_tar_pos = loop_exist(ref_m_tar,
        // [&](scaler_t x){return (x>0 && fabs(x)<fabs(uncer_cut));}); // number of live elements belong to ref, not
        // belong to tar

        bool ref_m_tar_neg = loop_exist(ref_m_tar, [](scaler_t x) {
            return x < -0.01;
        });  // number of live elements belong to tar, not belong to ref
        bool ref_m_tar_pos = loop_exist(ref_m_tar, [](scaler_t x) {
            return x > 0.01;
        });  // number of live elements belong to ref, not belong to tar

        if ((!ref_m_tar_neg) && (!ref_m_tar_pos)) {
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
    }

    return OTHER;

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

    int num_ref = loop_count(ref.m_proj, [](scaler_t x) { return x > 0; });
    int num_tar = loop_count(tar.m_proj, [](scaler_t x) { return x > 0; });
    int num_dead_ref = loop_count(ref.m_proj, [&](scaler_t x) {
        return x < (-1) * uncer_cut;
    });  // -1e12 for dead pixels, TODO: make sure this is robust
    int num_dead_tar = loop_count(tar.m_proj, [&](scaler_t x) { return x < (-1) * uncer_cut; });
    int num_inter = loop_count(inter_proj, [](scaler_t x) { return x > 0; });

    scaler_t charge_ref = loop_sum(ref.m_proj, [](scaler_t x) { return (x > 0) * x; });
    scaler_t charge_tar = loop_sum(tar.m_proj, [](scaler_t x) { return (x > 0) * x; });
    scaler_t charge_inter = loop_sum(inter_proj, [](scaler_t x) { return (x > 0) * x; });

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
