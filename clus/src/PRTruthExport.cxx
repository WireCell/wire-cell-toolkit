// PRTruthExport -- the pattern-recognition chain's point assignment, exported
// for MC-truth labels of the reconstructed objects (wcp-porting-img
// pr_ml/docs/23 sec 11, WP0.2 "R" exporter; label definitions
// pr_ml/docs/15 sec 2.4).  Doc 19 sec 3.7 called this exporter
// PRTruthMatcher; as for the Steiner graph (SteinerTruthLabeler), the C++
// side only exports what exists nowhere else and the matching is offline
// (truth/python/pr_truth.py).
//
// What exists nowhere else is the PR's own assignment of cluster points to
// segments: clustering_points_segments() builds each segment's
// "associate_points" cloud together with the cluster point indices it came
// from (Segment::global_indices).  The segments' fits, vertices, showers and
// the neutrino vertex are already written by PrDisplayDump (pr_display,
// calib-pr-evt<N>.json) with the join keys used here.
//
// For every TrackFitting the neutrino tagger left on the grouping (the named
// slots "nu0", "nu1", ... when present, else the unnamed one), every segment
// of its PR graph (PR::ordered_edges) adds rows to two cluster-node local
// clouds of the segment's cluster:
//
//   pr_segments    one row per segment:
//     seg_id       cluster_id*1000 + segment id (graph index when unset):
//                  the id of PrDisplayDump's segments[] / showers[].id and of
//                  the Bee particle-flow tree
//     seg_gidx     graph index
//     nu_index     candidate index (0 for the unnamed slot)
//     n_fit        number of fit points
//     n_cloud      rows of the associate_points cloud (0 without one)
//     n_index      stored cluster point indices (0 without them)
//     n_exact      cloud rows equal to the cluster point at their index
//   pr_seg_points  one row per associate_points row:
//     seg_id, nu_index
//     p_index      default-scope k-d index of the cluster point (the same
//                  index SteinerTruthLabeler's p_index uses)
//     src          0: the stored index, the point equals the cloud row
//                  1: stored index stale or absent; nearest cluster point
//     d_join       distance cloud row -> that point (0 for src 0)
//     x_raw, y_raw, z_raw   raw "3d" coordinates of the point (join check
//                  against the dump's blob clouds)
//
// No zero-row cloud is written; a cluster without segments gets nothing.
//
// A NEW component: absent from every existing pipeline => no output changes.
// It reads the PR graph and never touches it or any existing cloud.  Place it
// after the last stage that changes the graph (after tagger_check_neutrino).
// Deterministic: candidates in slot order, segments in graph order, clusters
// in tree order.

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRSegment.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Logging.h"

#include <cmath>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

class PRTruthExport;
WIRECELL_FACTORY(PRTruthExport, PRTruthExport,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

static Log::logptr_t logger() {
    static Log::logptr_t l = Log::logger("clus.PRTruthExport");
    return l;
}

namespace {
    // PrDisplayDump.cxx pf_node_id (the Bee producer's encoding), duplicated.
    int seg_display_id(const PR::SegmentPtr& seg)
    {
        int sid = seg->id();
        if (sid < 0) sid = static_cast<int>(seg->get_graph_index());
        const auto* cl = seg->cluster();
        return cl ? cl->get_cluster_id() * 1000 + sid : sid;
    }

    struct Rows {
        std::vector<int> seg_id, seg_gidx, nu_index, n_fit, n_cloud, n_index, n_exact;
        std::vector<int> p_seg_id, p_nu_index, p_index, p_src;
        std::vector<double> p_d_join, p_x_raw, p_y_raw, p_z_raw;
    };
}

class PRTruthExport : public IConfigurable, public Clus::IEnsembleVisitor {
public:
    PRTruthExport() {}
    virtual ~PRTruthExport() {}

    void configure(const WireCell::Configuration& config) {
        m_grouping = get<std::string>(config, "grouping", m_grouping);
        m_assoc_cloud = get<std::string>(config, "assoc_cloud", m_assoc_cloud);
        m_out_segments = get<std::string>(config, "out_segments", m_out_segments);
        m_out_points = get<std::string>(config, "out_points", m_out_points);
    }

    WireCell::Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping;
        cfg["assoc_cloud"] = m_assoc_cloud;
        cfg["out_segments"] = m_out_segments;
        cfg["out_points"] = m_out_points;
        return cfg;
    }

    void visit(Ensemble& ensemble) const {
        auto vec = ensemble.with_name(m_grouping);
        if (vec.empty()) {
            logger()->warn("PRTruthExport: no '{}' grouping found, skipping", m_grouping);
            return;
        }
        Grouping& grouping = *vec.at(0);

        // the candidates' fitters: named slots nu0, nu1, ... else the unnamed one
        std::vector<std::shared_ptr<TrackFitting>> tfs;
        for (int i = 0;; ++i) {
            auto tf = grouping.get_track_fitting("nu" + std::to_string(i));
            if (!tf) break;
            tfs.push_back(tf);
        }
        if (tfs.empty()) {
            if (auto tf = grouping.get_track_fitting()) tfs.push_back(tf);
        }
        if (tfs.empty()) {
            logger()->info("PRTruthExport: no TrackFitting on grouping '{}' (no neutrino candidate), nothing exported",
                           m_grouping);
            return;
        }

        std::unordered_map<const Cluster*, Rows> rows;   // lookup only; output in tree order
        std::unordered_map<const Cluster*, bool> in_grouping;
        for (const Cluster* cl : grouping.children()) in_grouping[cl] = true;

        size_t nseg = 0, nseg_nocluster = 0, nseg_nocloud = 0, nseg_noindex = 0, npts = 0, npts_src0 = 0;
        for (size_t inu = 0; inu < tfs.size(); ++inu) {
            auto graph = tfs[inu]->get_graph();
            if (!graph) continue;
            for (auto ed : PR::ordered_edges(*graph)) {
                auto seg = (*graph)[ed].segment;
                if (!seg) continue;
                Cluster* cl = seg->cluster();
                if (!cl || !in_grouping.count(cl)) { ++nseg_nocluster; continue; }
                ++nseg;
                Rows& r = rows[cl];
                const int sid = seg_display_id(seg);
                const auto dpc = seg->dpcloud(m_assoc_cloud);
                const auto& gind = seg->global_indices(m_assoc_cloud);
                const size_t ncloud = dpc ? dpc->npoints() : 0;
                if (!ncloud) ++nseg_nocloud;
                if (gind.empty()) ++nseg_noindex;

                const auto& sv = cl->sv3d();
                const auto& kd = cl->kd3d();
                const size_t ncl = cl->npoints();
                int nexact = 0;
                for (size_t ip = 0; ip < ncloud; ++ip) {
                    const auto& pts = dpc->points();
                    const geo_point_t q(pts.x[ip], pts.y[ip], pts.z[ip]);
                    size_t ind = 0;
                    int src = 1;
                    double d = 0;
                    if (gind.size() == ncloud && gind[ip] < ncl) {
                        const auto p = cl->point3d(gind[ip]);
                        if (p.x() == q.x() && p.y() == q.y() && p.z() == q.z()) {
                            ind = gind[ip];
                            src = 0;
                        }
                    }
                    if (src) {
                        const auto nn = cl->kd_knn(1, q);
                        ind = nn.at(0).first;
                        d = std::sqrt(nn.at(0).second);
                    }
                    else {
                        ++nexact;
                    }
                    const auto* node = sv.nodes().at(kd.major_index(ind));
                    const size_t row = kd.minor_index(ind);
                    double xr = 0, yr = 0, zr = 0;
                    const auto& lpcs = node->value.local_pcs();
                    auto it3 = lpcs.find("3d");
                    if (it3 != lpcs.end()) {
                        xr = it3->second.get("x")->elements<double>()[row];
                        yr = it3->second.get("y")->elements<double>()[row];
                        zr = it3->second.get("z")->elements<double>()[row];
                    }
                    r.p_seg_id.push_back(sid);
                    r.p_nu_index.push_back(inu);
                    r.p_index.push_back(ind);
                    r.p_src.push_back(src);
                    r.p_d_join.push_back(d);
                    r.p_x_raw.push_back(xr);
                    r.p_y_raw.push_back(yr);
                    r.p_z_raw.push_back(zr);
                }
                npts += ncloud;
                npts_src0 += nexact;
                r.seg_id.push_back(sid);
                r.seg_gidx.push_back(seg->get_graph_index());
                r.nu_index.push_back(inu);
                r.n_fit.push_back(seg->fits().size());
                r.n_cloud.push_back(ncloud);
                r.n_index.push_back(gind.size());
                r.n_exact.push_back(nexact);
            }
        }

        size_t ncl_out = 0;
        for (Cluster* cl : grouping.children()) {   // tree order: deterministic
            auto it = rows.find(cl);
            if (it == rows.end()) continue;
            Rows& r = it->second;
            auto& lpcs = cl->local_pcs();
            PointCloud::Dataset segs;
            segs.add("seg_id", PointCloud::Array(r.seg_id));
            segs.add("seg_gidx", PointCloud::Array(r.seg_gidx));
            segs.add("nu_index", PointCloud::Array(r.nu_index));
            segs.add("n_fit", PointCloud::Array(r.n_fit));
            segs.add("n_cloud", PointCloud::Array(r.n_cloud));
            segs.add("n_index", PointCloud::Array(r.n_index));
            segs.add("n_exact", PointCloud::Array(r.n_exact));
            lpcs.erase(m_out_segments);
            lpcs.emplace(m_out_segments, std::move(segs));
            lpcs.erase(m_out_points);
            if (!r.p_index.empty()) {
                PointCloud::Dataset pts;
                pts.add("seg_id", PointCloud::Array(r.p_seg_id));
                pts.add("nu_index", PointCloud::Array(r.p_nu_index));
                pts.add("p_index", PointCloud::Array(r.p_index));
                pts.add("src", PointCloud::Array(r.p_src));
                pts.add("d_join", PointCloud::Array(r.p_d_join));
                pts.add("x_raw", PointCloud::Array(r.p_x_raw));
                pts.add("y_raw", PointCloud::Array(r.p_y_raw));
                pts.add("z_raw", PointCloud::Array(r.p_z_raw));
                lpcs.emplace(m_out_points, std::move(pts));
            }
            ++ncl_out;
        }
        logger()->info("PRTruthExport: {} candidate(s), {} segment(s) in {} cluster(s) ({} without an '{}' cloud, "
                       "{} without stored indices, {} outside the grouping); {} assigned points, {} joined by the "
                       "stored index",
                       tfs.size(), nseg, ncl_out, nseg_nocloud, m_assoc_cloud, nseg_noindex, nseg_nocluster,
                       npts, npts_src0);
    }

private:
    std::string m_grouping{"live"};
    std::string m_assoc_cloud{"associate_points"};
    std::string m_out_segments{"pr_segments"};
    std::string m_out_points{"pr_seg_points"};
};

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
