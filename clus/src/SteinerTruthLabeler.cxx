// SteinerTruthLabeler -- MC-truth labels on the Steiner graph (wcp-porting-img
// pr_ml/docs/23, WP0.2 "G" exporter; label definitions pr_ml/docs/15 sec 2.2
// and 12 sec 4.3).
//
// For every cluster of the grouping that carries a "steiner_graph" and its
// "steiner_pc" (CreateSteinerGraph), adds two cluster-node local point clouds:
//
//   steiner_truth  one row per Steiner vertex, same order as steiner_pc
//                  (vertex i of steiner_graph == row i of steiner_pc):
//     tid          majority track key of the cluster points within
//                  join_radius of the vertex (count vote, ties to the smaller
//                  key); the nearest point's key when none is within the radius
//     tid_frac     share of the joined points carrying tid (0 without a join)
//     npts_join    number of cluster points within join_radius
//     d_near       distance to the nearest cluster point (default-scope coords)
//     tid_near     track key of that nearest point
//     orig_tid, pdg, cls, iid, pf_node, d_traj, pid_frac
//                  copied from the "truth" row of the nearest point carrying tid
//     p_index      default-scope k-d index of the nearest point (node -> P map)
//     x_raw, y_raw, z_raw   raw "3d" x/y/z of the nearest point
//     dx_cor, dy_cor, dz_cor  default-scope minus raw coordinates of that point
//                  (the frame shift a polyline in the raw frame needs)
//   steiner_edges  one row per edge of steiner_graph, boost::edges order:
//     src, dst     vertex indices (rows of steiner_pc / steiner_truth)
//     weight       the routing price stored on the edge
//     length       Euclidean length in the default-scope coordinates
//
// The per-point "truth" cloud is the one TruthPointLabeler (WireCellTruth)
// hangs on every live blob node in the clustering job; its names and
// sentinels are those of truth/inc/WireCellTruth/Schema.h (tid -1 ghost,
// -2 shared).  This file only reads them by name: clus does not link the
// truth plugin.  A blob without a "truth" cloud yields tid_none (-3).
//
// The retiled cloud the Steiner graph is built on carries no truth (the
// retiled blobs are new), hence the geometric join against the original
// cluster's points in the same default scope (doc 23 sec 2).  The edge
// classes (on-trajectory / junction / false bridge / ghost) and the flow
// sign need the depo polylines and the parent relations and are computed
// offline from these two clouds (truth/python/steiner_truth.py).
//
// A NEW component: absent from every existing pipeline => no output changes.
// It never touches the graph or any existing cloud.  Deterministic: clusters
// in tree order, vertices and edges in graph order.

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Logging.h"

#include <boost/graph/adjacency_list.hpp>

#include <cmath>
#include <map>

class SteinerTruthLabeler;
WIRECELL_FACTORY(SteinerTruthLabeler, SteinerTruthLabeler,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

static Log::logptr_t logger() {
    static Log::logptr_t l = Log::logger("clus.SteinerTruthLabeler");
    return l;
}

namespace {
    constexpr int tid_none = -3;   // no "truth" cloud on the blob

    // The "truth" row of one cluster point: the blob's truth arrays and the
    // point's index in them.  null arrays when the blob has no truth cloud.
    struct TruthRow {
        const PointCloud::Dataset* ds{nullptr};
        size_t row{0};
        double x{0}, y{0}, z{0};   // raw "3d" coordinates
    };

    template <typename T>
    T value(const PointCloud::Dataset* ds, const std::string& key, size_t row, T def)
    {
        if (!ds) return def;
        auto arr = ds->get(key);
        if (!arr) return def;
        return arr->elements<T>()[row];   // elements() checks the element size
    }
}

class SteinerTruthLabeler : public IConfigurable, public Clus::IEnsembleVisitor {
public:
    SteinerTruthLabeler() {}
    virtual ~SteinerTruthLabeler() {}

    void configure(const WireCell::Configuration& config) {
        m_grouping = get<std::string>(config, "grouping", m_grouping);
        m_graph = get<std::string>(config, "graph", m_graph);
        m_steiner_pc = get<std::string>(config, "steiner_pc", m_steiner_pc);
        m_truth_pc = get<std::string>(config, "truth_pc", m_truth_pc);
        m_out_nodes = get<std::string>(config, "out_nodes", m_out_nodes);
        m_out_edges = get<std::string>(config, "out_edges", m_out_edges);
        m_join_radius = get<double>(config, "join_radius", m_join_radius);
    }

    WireCell::Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping;
        cfg["graph"] = m_graph;
        cfg["steiner_pc"] = m_steiner_pc;
        cfg["truth_pc"] = m_truth_pc;
        cfg["out_nodes"] = m_out_nodes;
        cfg["out_edges"] = m_out_edges;
        cfg["join_radius"] = m_join_radius;
        return cfg;
    }

    void visit(Ensemble& ensemble) const {
        auto vec = ensemble.with_name(m_grouping);
        if (vec.empty()) {
            logger()->warn("SteinerTruthLabeler: no '{}' grouping found, skipping", m_grouping);
            return;
        }
        Grouping& grouping = *vec.at(0);
        size_t ncl = 0, nvtx = 0, nedge = 0, njoined = 0, nnotruth = 0, nghost = 0, nshared = 0;
        for (Cluster* cluster : grouping.children()) {   // tree order: deterministic
            if (!cluster->has_graph(m_graph) || !cluster->has_pc(m_steiner_pc)) continue;
            const auto& graph = cluster->get_graph(m_graph);
            const auto& spc = cluster->get_pc(m_steiner_pc);
            const auto& coords = cluster->get_default_scope().coords;
            const auto xa = spc.get(coords.at(0));
            const auto ya = spc.get(coords.at(1));
            const auto za = spc.get(coords.at(2));
            if (!xa || !ya || !za) continue;
            const auto& xs = xa->elements<double>();
            const auto& ys = ya->elements<double>();
            const auto& zs = za->elements<double>();
            const size_t nv = boost::num_vertices(graph);
            if (nv != xs.size()) {
                raise<ValueError>("SteinerTruthLabeler: cluster %d: %d graph vertices but %d %s rows",
                                  cluster->ident(), nv, xs.size(), m_steiner_pc);
            }

            // default-scope k-d index -> (blob node, row of its clouds)
            const auto& sv = cluster->sv3d();
            const auto& kd = cluster->kd3d();
            auto truth_row = [&](size_t ind) {
                TruthRow tr;
                const auto* node = sv.nodes().at(kd.major_index(ind));
                tr.row = kd.minor_index(ind);
                const auto& lpcs = node->value.local_pcs();
                auto it = lpcs.find(m_truth_pc);
                if (it != lpcs.end()) tr.ds = &it->second;
                auto it3 = lpcs.find("3d");
                if (it3 != lpcs.end()) {
                    tr.x = value<double>(&it3->second, "x", tr.row, 0.0);
                    tr.y = value<double>(&it3->second, "y", tr.row, 0.0);
                    tr.z = value<double>(&it3->second, "z", tr.row, 0.0);
                }
                return tr;
            };
            auto tid_of = [&](const TruthRow& tr) { return tr.ds ? value<int>(tr.ds, "tid", tr.row, tid_none) : tid_none; };

            std::vector<int> tid(nv), npts(nv), tid_near(nv), orig_tid(nv), pdg(nv), cls(nv), iid(nv), pf_node(nv), p_index(nv);
            std::vector<double> tid_frac(nv), d_near(nv), d_traj(nv), pid_frac(nv);
            std::vector<double> x_raw(nv), y_raw(nv), z_raw(nv), dx_cor(nv), dy_cor(nv), dz_cor(nv);
            for (size_t iv = 0; iv < nv; ++iv) {
                const geo_point_t p(xs[iv], ys[iv], zs[iv]);
                const auto nn = cluster->kd_knn(1, p);
                const size_t inear = nn.at(0).first;
                const TruthRow near = truth_row(inear);
                const auto pnear = cluster->point3d(inear);

                // count vote over the points within the join radius
                const auto res = cluster->kd_radius(m_join_radius, p);
                std::map<int, int> votes;   // int-keyed: ordered, deterministic
                for (const auto& [ind, d2] : res) ++votes[tid_of(truth_row(ind))];
                int best = tid_of(near), nbest = 0;
                for (const auto& [k, n] : votes) {
                    if (n > nbest) { best = k; nbest = n; }   // ties keep the smaller key
                }
                // the nearest point carrying the chosen key supplies its attributes
                TruthRow rep = near;
                if (!res.empty()) {
                    double dbest = -1;
                    for (const auto& [ind, d2] : res) {
                        const TruthRow tr = truth_row(ind);
                        if (tid_of(tr) != best) continue;
                        if (dbest < 0 || d2 < dbest || (d2 == dbest && ind < inear)) { dbest = d2; rep = tr; }
                    }
                    ++njoined;
                }
                tid[iv] = best;
                npts[iv] = res.size();
                tid_frac[iv] = res.empty() ? 0.0 : double(nbest) / double(res.size());
                d_near[iv] = std::sqrt(nn.at(0).second);
                tid_near[iv] = tid_of(near);
                orig_tid[iv] = value<int>(rep.ds, "orig_tid", rep.row, tid_none);
                pdg[iv] = value<int>(rep.ds, "pdg", rep.row, 0);
                cls[iv] = value<int>(rep.ds, "cls", rep.row, -1);
                iid[iv] = value<int>(rep.ds, "iid", rep.row, -1);
                pf_node[iv] = value<int>(rep.ds, "pf_node", rep.row, -1);
                d_traj[iv] = value<double>(rep.ds, "d_traj", rep.row, -1.0);
                pid_frac[iv] = value<double>(rep.ds, "pid_frac", rep.row, 0.0);
                p_index[iv] = inear;
                x_raw[iv] = near.x;
                y_raw[iv] = near.y;
                z_raw[iv] = near.z;
                dx_cor[iv] = pnear.x() - near.x;
                dy_cor[iv] = pnear.y() - near.y;
                dz_cor[iv] = pnear.z() - near.z;
                if (best == tid_none) ++nnotruth;
                else if (best == -1) ++nghost;
                else if (best == -2) ++nshared;
            }

            std::vector<int> src, dst;
            std::vector<double> weight, length;
            const auto wmap = boost::get(boost::edge_weight, graph);
            for (auto [eit, eend] = boost::edges(graph); eit != eend; ++eit) {
                const size_t a = boost::source(*eit, graph), b = boost::target(*eit, graph);
                src.push_back(a);
                dst.push_back(b);
                weight.push_back(wmap[*eit]);
                length.push_back(std::sqrt((xs[a] - xs[b]) * (xs[a] - xs[b]) + (ys[a] - ys[b]) * (ys[a] - ys[b]) +
                                           (zs[a] - zs[b]) * (zs[a] - zs[b])));
            }

            PointCloud::Dataset nodes;
            nodes.add("tid", PointCloud::Array(tid));
            nodes.add("tid_frac", PointCloud::Array(tid_frac));
            nodes.add("npts_join", PointCloud::Array(npts));
            nodes.add("d_near", PointCloud::Array(d_near));
            nodes.add("tid_near", PointCloud::Array(tid_near));
            nodes.add("orig_tid", PointCloud::Array(orig_tid));
            nodes.add("pdg", PointCloud::Array(pdg));
            nodes.add("cls", PointCloud::Array(cls));
            nodes.add("iid", PointCloud::Array(iid));
            nodes.add("pf_node", PointCloud::Array(pf_node));
            nodes.add("d_traj", PointCloud::Array(d_traj));
            nodes.add("pid_frac", PointCloud::Array(pid_frac));
            nodes.add("p_index", PointCloud::Array(p_index));
            nodes.add("x_raw", PointCloud::Array(x_raw));
            nodes.add("y_raw", PointCloud::Array(y_raw));
            nodes.add("z_raw", PointCloud::Array(z_raw));
            nodes.add("dx_cor", PointCloud::Array(dx_cor));
            nodes.add("dy_cor", PointCloud::Array(dy_cor));
            nodes.add("dz_cor", PointCloud::Array(dz_cor));

            PointCloud::Dataset edges;
            edges.add("src", PointCloud::Array(src));
            edges.add("dst", PointCloud::Array(dst));
            edges.add("weight", PointCloud::Array(weight));
            edges.add("length", PointCloud::Array(length));

            auto& lpcs = cluster->local_pcs();
            lpcs.erase(m_out_nodes);
            lpcs.emplace(m_out_nodes, std::move(nodes));
            // a graph without edges gets no edge cloud (no zero-row local
            // cloud; cf. the empty steiner_pc guard of CreateSteinerGraph)
            lpcs.erase(m_out_edges);
            if (!src.empty()) lpcs.emplace(m_out_edges, std::move(edges));

            ++ncl;
            nvtx += nv;
            nedge += src.size();
        }
        logger()->info("SteinerTruthLabeler: {} cluster(s), {} vertices ({} joined within {:.2f} cm; {} ghost, {} shared, "
                       "{} without truth), {} edges",
                       ncl, nvtx, njoined, m_join_radius / units::cm, nghost, nshared, nnotruth, nedge);
    }

private:
    std::string m_grouping{"live"};
    std::string m_graph{"steiner_graph"};
    std::string m_steiner_pc{"steiner_pc"};
    std::string m_truth_pc{"truth"};
    std::string m_out_nodes{"steiner_truth"};
    std::string m_out_edges{"steiner_edges"};
    double m_join_radius{0.5 * units::cm};   // doc 15 sec 2.2
};

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
