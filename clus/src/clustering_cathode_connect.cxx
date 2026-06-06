// SBND cathode-crossing connector (default-OFF, retireable).
//
// After Q/L matching, a cosmic that crosses the central cathode is reconstructed
// as two halves (one per TPC).  The generic merge passes (regular/parallel_prolong/
// ...) only connect such a pair when its closest-point distance is <= 3 cm; a data
// cathode-crosser's distance = physical cathode gap (~1.5 cm) + a drift-x calibration
// residual, which pushes the worst cases just past 3 cm so they are left as two
// clusters (see clus/docs/cathode-crossing-clustering.md).
//
// This pass merges exactly that population using a narrow, cathode-specific cut set
// that cannot fire within a single TPC (it requires the closest points to sit in
// DIFFERENT TPCs).  ALWAYS required:
//   1. the two local track directions are collinear      (vhough_transform)
//   2. both closest points end at the cathode plane      (|x - cathode_x| < cathode_x_cut)
//   3. the two closest points are in separate TPCs       (wpid.apa() differ)
//   4. the closest points are at the same DRIFT depth    (|x1 - x2| < drift_cut, tight)
// The 3D closest-point distance is then handled in two regimes:
//   - close (dis < dis_cut): the p1->p2 connection vector is dominated by the drift-x
//     offset (the calibration artifact, ~along the drift axis), not the track, so the
//     generic passes' connection-alignment test rejects these crossers; here we accept
//     on the local (Hough) track-direction collinearity alone -- the hole the pass fills.
//   - far  (dis_cut <= dis < max_dis): a shallow-angle crosser travels within the cathode
//     plane, so the halves are offset transversely (large dis) and the p1->p2 connection
//     vector becomes a reliable direction.  Two refinements vs close: (a) the local Hough
//     direction can be unreliable when one half is a dense blob at the cathode, so the
//     cluster PCA principal axis is added as an ALTERNATIVE direction estimate (Hough OR
//     PCA may satisfy collinearity -- additional, it does not replace the Hough test);
//     (b) require the connection vector to align with the track (Hough or PCA) within
//     conn_far_cut, which rejects parallel-offset cosmics (connection ~perpendicular,
//     >=50 deg) while passing a real crosser (a few deg to ~20 deg).
// plus the all-APA flash-time gate (same matched flash group), so only coincident
// clusters are paired.  It runs after the generic merge passes, so it can only ADD
// merges they missed; it is the same logic QLMatching uses to flag cross-TPC pairs
// (flag_xtpc_consistent / cull_cross_tpc), here used to connect rather than flag.
//
// The misalignment it compensates is a calibration artifact: when the pos_offset /
// space-charge transverse calibration tightens, the residual shrinks, the generic
// 3 cm path catches these crossers, and this pass can be turned off and retired.

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include <unordered_map>
#include <cmath>

class ClusteringCathodeConnect;
WIRECELL_FACTORY(ClusteringCathodeConnect, ClusteringCathodeConnect,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

static void clustering_cathode_connect(
    Grouping& live_grouping,
    const Tree::Scope& scope,
    double drift_cut,        // |x1 - x2| below which the two points are at the same drift depth
    double dis_cut,          // 3D distance below which we accept on track collinearity alone
    double max_dis,          // 3D distance ceiling (connection-aligned far regime)
    double angle_cut,        // degrees; collinearity of the two half-directions
    double conn_far_cut,     // degrees; far-regime alignment of the p1->p2 connection with the track
    double cathode_x,        // cathode position in the (T0-corrected) clustering frame
    double cathode_x_cut,    // |x - cathode_x| below which a closest point "ends at the cathode"
    double hough_radius,
    double min_length,
    double flash_t0_window);

class ClusteringCathodeConnect : public IConfigurable, public Clus::IEnsembleVisitor, private NeedScope {
public:
    ClusteringCathodeConnect() {}
    virtual ~ClusteringCathodeConnect() {}

    void configure(const WireCell::Configuration& config) {
        NeedScope::configure(config);

        drift_cut_       = get(config, "drift_cut", 5*units::cm);
        dis_cut_         = get(config, "dis_cut", 5*units::cm);
        max_dis_         = get(config, "max_dis", 25*units::cm);
        angle_cut_       = get(config, "angle_cut", 10.0);          // degrees
        conn_far_cut_    = get(config, "conn_far_cut", 30.0);       // degrees
        cathode_x_       = get(config, "cathode_x", 0.0);
        cathode_x_cut_   = get(config, "cathode_x_cut", 3.5*units::cm);
        hough_radius_    = get(config, "hough_radius", 20*units::cm);
        min_length_      = get(config, "min_length", 10*units::cm);
        flash_t0_window_ = get(config, "flash_t0_window", 80*units::ns);
    }
    virtual Configuration default_configuration() const {
        Configuration cfg;
        return cfg;
    }

    void visit(Ensemble& ensemble) const {
        auto& live = *ensemble.with_name("live").at(0);
        clustering_cathode_connect(live, m_scope, drift_cut_, dis_cut_, max_dis_, angle_cut_,
                                   conn_far_cut_, cathode_x_, cathode_x_cut_, hough_radius_,
                                   min_length_, flash_t0_window_);
    }

private:
    double drift_cut_{5*units::cm};
    double dis_cut_{5*units::cm};
    double max_dis_{25*units::cm};
    double angle_cut_{10.0};
    double conn_far_cut_{30.0};
    double cathode_x_{0.0};
    double cathode_x_cut_{3.5*units::cm};
    double hough_radius_{20*units::cm};
    double min_length_{10*units::cm};
    double flash_t0_window_{80*units::ns};
};


// Fold an angle (deg, 0..180) about 180 -> collinearity (0 = parallel or anti-parallel).
static inline double collinear_deg(double angle_rad)
{
    double a = angle_rad / 3.1415926 * 180.0;
    return std::min(a, 180.0 - a);
}

// Return true if (cluster1, cluster2) is a cathode-crossing pair to connect.
static bool is_cathode_crossing_pair(
    const Cluster& cluster1,
    const Cluster& cluster2,
    double length_1,
    double length_2,
    double drift_cut,
    double dis_cut,
    double max_dis,
    double angle_cut,
    double conn_far_cut,
    double cathode_x,
    double cathode_x_cut,
    double hough_radius)
{
    geo_point_t p1;
    geo_point_t p2;
    // Find_Closest_Points ignores its length_cut argument (it returns the global
    // closest 3D point pair); we apply our own cuts below.
    double dis = WireCell::Clus::Facade::Find_Closest_Points(cluster1, cluster2,
                                                             length_1, length_2,
                                                             max_dis, p1, p2);

    // (3) the two closest points must sit in DIFFERENT TPCs (this is what makes the
    //     pass incapable of acting within a single TPC).
    auto wpid_p1 = cluster1.wpid(p1);
    auto wpid_p2 = cluster2.wpid(p2);
    if (wpid_p1.apa() == wpid_p2.apa()) return false;

    // (2) both closest points must end at the cathode plane.
    if (std::fabs(p1.x() - cathode_x) >= cathode_x_cut) return false;
    if (std::fabs(p2.x() - cathode_x) >= cathode_x_cut) return false;

    // (4) the two points must be at the same DRIFT depth (always tight): the only
    //     thing between the two halves is the ~1.5 cm cathode gap plus a drift-x
    //     calibration residual (observed up to ~4.1 cm in data).
    if (std::fabs(p1.x() - p2.x()) >= drift_cut) return false;

    // outer 3D distance ceiling.
    if (dis >= max_dis) return false;

    // (1) local track directions, from the Hough transform at the closest points.
    geo_point_t dir1 = cluster1.vhough_transform(p1, hough_radius);
    geo_point_t dir2 = cluster2.vhough_transform(p2, hough_radius);
    double tt_hough = collinear_deg(dir1.angle(dir2));

    if (dis < dis_cut) {
        // CLOSE regime: the p1->p2 connection vector is dominated by the drift-x offset
        //   (the calibration artifact, ~along the drift axis), not the track, so the
        //   generic passes' connection-alignment test rejects these crossers -- the hole
        //   this pass fills.  Accept on the local half-track collinearity alone.
        return tt_hough < angle_cut;
    }

    // FAR regime: the two halves are offset within the cathode plane (large transverse
    //   separation from a shallow-angle crosser travelling along the cathode).  Two
    //   things change vs CLOSE: (a) the local Hough direction can be unreliable when one
    //   half is a dense blob near the cathode, so we ADD the cluster PCA principal axis
    //   as an alternative direction estimate (Hough OR PCA may satisfy collinearity);
    //   (b) the now-long p1->p2 connection vector is meaningful, so we require it to
    //   align with the track (PCA or Hough) -- this rejects parallel-offset cosmics,
    //   whose connection is ~perpendicular (>=50 deg), while passing a real crosser.
    geo_point_t pca1, pca2;
    bool have_pca = false;
    {
        const auto& ax1 = cluster1.get_pca().axis;
        const auto& ax2 = cluster2.get_pca().axis;
        if (!ax1.empty() && !ax2.empty()) {
            pca1.set(ax1.at(0).x(), ax1.at(0).y(), ax1.at(0).z());
            pca2.set(ax2.at(0).x(), ax2.at(0).y(), ax2.at(0).z());
            have_pca = true;
        }
    }
    double tt_pca = have_pca ? collinear_deg(pca1.angle(pca2)) : 999.0;
    // track-track collinear: Hough OR (additional) PCA.
    if (tt_hough >= angle_cut && tt_pca >= angle_cut) return false;

    // connection vector aligned with the track: Hough OR (additional) PCA.
    geo_point_t conn(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
    double cc_hough = collinear_deg(conn.angle(dir1));
    double cc_pca = have_pca ? std::min(collinear_deg(conn.angle(pca1)),
                                        collinear_deg(conn.angle(pca2))) : 999.0;
    if (cc_hough >= conn_far_cut && cc_pca >= conn_far_cut) return false;

    return true;
}


static void clustering_cathode_connect(
    Grouping& live_grouping,
    const Tree::Scope& scope,
    double drift_cut,
    double dis_cut,
    double max_dis,
    double angle_cut,
    double conn_far_cut,
    double cathode_x,
    double cathode_x_cut,
    double hough_radius,
    double min_length,
    double flash_t0_window)
{
    // prepare graph ... (same skeleton as the other merge passes)
    typedef cluster_connectivity_graph_t Graph;
    Graph g;
    std::unordered_map<int, int> ilive2desc;
    std::unordered_map<const Cluster*, int> map_cluster_index;
    auto live_clusters = live_grouping.children();

    // Cross-TPC connection only makes sense for clusters coincident in flash time;
    // unmatched clusters get unique singleton groups (never linked).
    std::map<const Cluster*, int> flash_t0_group = assign_flash_t0_groups(live_clusters, flash_t0_window);

    // Vertex index in children() order (merge_clusters dereferences against children()).
    for (size_t ilive = 0; ilive < live_clusters.size(); ++ilive) {
        const auto& live = live_clusters.at(ilive);
        if (live->get_default_scope().hash() != scope.hash()) {
            live->set_default_scope(scope);
        }
        map_cluster_index[live] = ilive;
        ilive2desc[ilive] = boost::add_vertex(ilive, g);
    }
    sort_clusters(live_clusters);  // deterministic edge-building order

    for (size_t i = 0; i != live_clusters.size(); i++) {
        auto cluster_1 = live_clusters.at(i);
        if (!cluster_1->get_scope_filter(scope)) continue;
        if (cluster_1->get_length() < min_length) continue;
        for (size_t j = i + 1; j < live_clusters.size(); j++) {
            auto cluster_2 = live_clusters.at(j);
            if (!cluster_2->get_scope_filter(scope)) continue;
            if (cluster_2->get_length() < min_length) continue;
            if (flash_t0_group.at(cluster_1) != flash_t0_group.at(cluster_2)) continue;
            if (is_cathode_crossing_pair(*cluster_1, *cluster_2,
                                         cluster_1->get_length(), cluster_2->get_length(),
                                         drift_cut, dis_cut, max_dis, angle_cut, conn_far_cut,
                                         cathode_x, cathode_x_cut, hough_radius)) {
                boost::add_edge(ilive2desc[map_cluster_index[cluster_1]],
                                ilive2desc[map_cluster_index[cluster_2]], g);
            }
        }
    }

    merge_clusters(g, live_grouping);
}


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
