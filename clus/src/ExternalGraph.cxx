// pr_ml doc 79 (P3): read a segment graph built outside the toolkit and build it as a PR::Graph for one cluster.
// See WireCellClus/ExternalGraph.h for the file format.  Called only from TaggerCheckNeutrino's pr_fit_only = "external"
// mode (default off), so no existing output changes.  Deterministic: vertices, segments and clumps in file order.

#include "WireCellClus/ExternalGraph.h"
#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRVertex.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/String.h"

#include <algorithm>
#include <map>

using namespace WireCell;
using namespace WireCell::Clus;

namespace {
    WireCell::Point point_of(const Json::Value& j, const std::string& what)
    {
        if (!j.isArray() || j.size() != 3) {
            raise<ValueError>("external graph: %s is not a 3-vector", what);
        }
        return WireCell::Point(j[0].asDouble() * units::mm, j[1].asDouble() * units::mm, j[2].asDouble() * units::mm);
    }

    std::vector<WireCell::Point> points_of(const Json::Value& j, const std::string& what)
    {
        std::vector<WireCell::Point> out;
        if (!j.isArray()) raise<ValueError>("external graph: %s has no point list", what);
        out.reserve(j.size());
        for (Json::ArrayIndex i = 0; i < j.size(); ++i) out.push_back(point_of(j[i], what));
        return out;
    }
}

PR::ExternalGraph PR::read_external_graph(const std::string& path)
{
    Json::Value top = Persist::load(path);
    if (top.get("version", 0).asInt() != 1) raise<ValueError>("external graph %s: version must be 1", path);
    if (top.get("units", "").asString() != "mm") raise<ValueError>("external graph %s: units must be \"mm\"", path);
    ExternalGraph eg;
    std::map<int, size_t> vidx;
    for (const auto& jv : top["vertices"]) {
        ExternalGraph::Vtx v;
        v.id = jv["id"].asInt();
        v.p = point_of(jv["p"], "vertex " + std::to_string(v.id));
        if (vidx.count(v.id)) raise<ValueError>("external graph %s: vertex id %d repeated", path, v.id);
        vidx[v.id] = eg.vertices.size();
        eg.vertices.push_back(v);
    }
    const double tol = 0.01 * units::mm;
    for (const auto& js : top["segments"]) {
        ExternalGraph::Seg s;
        s.id = js["id"].asInt();
        const std::string what = "segment " + std::to_string(s.id);
        if (!js["v"].isArray() || js["v"].size() != 2) raise<ValueError>("external graph %s: %s needs two vertices", path, what);
        s.v0 = js["v"][0].asInt(); s.v1 = js["v"][1].asInt();
        if (!vidx.count(s.v0) || !vidx.count(s.v1)) raise<ValueError>("external graph %s: %s names an unknown vertex", path, what);
        s.pts = points_of(js["pts"], what);
        if (s.pts.size() < 2) raise<ValueError>("external graph %s: %s has fewer than 2 points", path, what);
        if ((s.pts.front() - eg.vertices[vidx[s.v0]].p).magnitude() > tol ||
            (s.pts.back() - eg.vertices[vidx[s.v1]].p).magnitude() > tol) {
            raise<ValueError>("external graph %s: %s does not end on its vertices", path, what);
        }
        s.type = js.get("type", "").asString();
        eg.segments.push_back(s);
    }
    for (const auto& jc : top["clumps"]) {
        ExternalGraph::Clump c;
        c.id = jc["id"].asInt();
        c.pts = points_of(jc["pts"], "clump " + std::to_string(c.id));
        if (c.pts.empty()) raise<ValueError>("external graph %s: clump %d has no point", path, c.id);
        c.type = jc.get("type", "").asString();
        eg.clumps.push_back(c);
    }
    return eg;
}

namespace {
    double median_distance(const Facade::Cluster& cluster, const std::vector<WireCell::Point>& pts)
    {
        std::vector<double> d;
        d.reserve(pts.size());
        for (const auto& p : pts) d.push_back((cluster.get_closest_wcpoint(p).second - p).magnitude());
        std::sort(d.begin(), d.end());
        const size_t n = d.size();
        return n % 2 ? d[n / 2] : 0.5 * (d[n / 2 - 1] + d[n / 2]);
    }
}

PR::ExternalGraphStats PR::build_external_graph(PatternAlgorithms& pa, Graph& graph, Facade::Cluster& cluster, TrackFitting& tf,
                                                IDetectorVolumes::pointer dv, const ExternalGraph& eg, double max_med_dist,
                                                const std::vector<const Facade::Cluster*>& others)
{
    // 0: used; 1: farther than max_med_dist; 2: another cluster is nearer
    auto verdict = [&](const std::vector<WireCell::Point>& pts, double& md) {
        md = median_distance(cluster, pts);
        if (md > max_med_dist) return 1;
        for (const auto* o : others) {
            if (!o || o == &cluster || o->npoints() == 0) continue;
            if (median_distance(*o, pts) < md) return 2;
        }
        return 0;
    };
    ExternalGraphStats st;
    std::map<int, const ExternalGraph::Vtx*> vin;
    for (const auto& v : eg.vertices) vin[v.id] = &v;
    std::map<int, VertexPtr> vmade;
    auto vertex = [&](int id) {
        auto it = vmade.find(id);
        if (it != vmade.end()) return it->second;
        VertexPtr v = make_vertex(graph);
        v->wcpt().point = vin.at(id)->p;
        v->cluster(&cluster);
        vmade[id] = v;
        return v;
    };
    for (const auto& s : eg.segments) {
        ++st.seg_read;
        double md = 0;
        const int v = verdict(s.pts, md);
        if (v == 1) { ++st.seg_skipped; continue; }
        if (v == 2) { ++st.seg_other_nearer; continue; }
        auto seg = pa.create_segment_for_cluster(cluster, dv, s.pts);
        if (!seg) { ++st.seg_refused; continue; }
        add_segment(graph, seg, vertex(s.v0), vertex(s.v1));
        tf.add_segment(seg);
        ++st.seg_used;
        st.seg_med_dist.push_back(md);
    }
    std::vector<std::vector<WireCell::Point>> clumps;
    for (const auto& c : eg.clumps) {
        ++st.clump_read;
        double md = 0;
        const int v = verdict(c.pts, md);
        if (v == 2) ++st.clump_other_nearer;
        if (v != 0) continue;
        clumps.push_back(c.pts);
        ++st.clump_used;
    }
    tf.set_charge_clumps(&cluster, std::move(clumps));
    return st;
}
