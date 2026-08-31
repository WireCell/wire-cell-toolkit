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
//     BOTH-LONG sub-case (default OFF, short_dir_len > 0): when both halves are long
//     enough for their cluster PCA to be reliable (>= short_dir_len) but the local Hough
//     at the dense cathode tips is noisy (e.g. 183096: tt_hough=24 deg, tt_pca=1.1 deg),
//     accept on the FAR-regime test applied to this close pair -- PCA-axis collinearity
//     AND a connection vector that continues the track (cc_pca < conn_far_cut); the
//     connection term rejects two distinct parallel cosmics passing close at the cathode.
//     SHORT-STUB sub-case (default OFF, short_dir_len > 0): when one half is too short
//     (< short_dir_len) for its own direction to be trusted and the other is a genuinely
//     long anchor (>= short_dir_len), the collinearity test is unreliable, so instead
//     require the anchor->stub connection to continue the ANCHOR's direction (Hough or
//     PCA) within conn_short_cut -- a gap-splintered stub attaches to its anchor while a
//     transversely-offset short cluster (connection ~perpendicular) is rejected.
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
// clusters are paired.  The gate is configurable (use_flash_t0, default on) so the
// pass can act on detectors without flash matching (e.g. PDHD).  It runs after the
// generic merge passes, so it can only ADD merges they missed; it is the same logic QLMatching uses to flag cross-TPC pairs
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
#include <cstdlib>   // std::getenv  (CATHODE_CONNECT_DEBUG instrumentation, removable)
#include <cstdio>    // std::fprintf

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
    double min_length,       // the longer half of a pair must reach this (the anchor track)
    double min_length_short, // the shorter member only needs this (admits a bridge fragment)
    double short_dir_len,    // below this a half's own direction is untrusted (prolong instead)
    double conn_short_cut,   // degrees; short-stub anchor->stub connection vs anchor direction
    double tip_touch_cut,    // below this 3D gap the cc connection term is dropped (both-long PCA branch)
    double tip_touch_angle_cut, // degrees; relaxed local-Hough collinearity in the tip-touch branch
    bool use_flash_t0,
    double flash_t0_window,
    double crosser_conn_relax, // degrees; loose cc_pca bound for drift-gated 6cm-cathode crossers (0=off)
    double crosser_pca_angle,  // degrees; raised tt_pca bound for drift-gated bent crossers (0=off)
    double cathode_band_dis);  // distance; near-cathode closest-approach retry for hard-gated crossers (0=off)

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
        // default min_length_short == min_length: both members must reach min_length,
        // i.e. byte-identical to the pre-asymmetric behaviour unless explicitly set.
        min_length_short_= get(config, "min_length_short", min_length_);
        // default short_dir_len == 0: short-stub prolongation branch OFF, i.e.
        // byte-identical to the pre-prolongation behaviour unless explicitly set.
        short_dir_len_   = get(config, "short_dir_len", 0.0);
        conn_short_cut_  = get(config, "conn_short_cut", 30.0);     // degrees
        // default tip_touch_cut == 0: the both-long PCA branch keeps its cc_pca
        // connection-alignment gate (byte-identical to before).  When > 0, a pair
        // whose closest 3D distance is below this is treated as "tips touching" and
        // accepted on PCA collinearity alone (the connection vector over a ~1 cm gap
        // is dominated by sub-cm transverse jitter and is uninformative there).
        tip_touch_cut_   = get(config, "tip_touch_cut", 0.0);
        // default tip_touch_angle_cut == angle_cut: the tip-touch local-Hough fallback
        // reduces to the primary close test (byte-identical, no new accepts).  Set
        // larger to also merge a tip-touching crosser whose curved half inflates its
        // GLOBAL PCA above angle_cut while its LOCAL Hough arm stays collinear.
        tip_touch_angle_cut_ = get(config, "tip_touch_angle_cut", angle_cut_);
        // default true: keep the flash-time coincidence gate (the original
        // behaviour).  Set false on detectors without flash matching, where no
        // cluster carries a matched flash and the gate would veto every pair.
        use_flash_t0_    = get(config, "use_flash_t0", true);
        flash_t0_window_ = get(config, "flash_t0_window", 80*units::ns);
        // default crosser_conn_relax == 0: the both-long PCA branch keeps its cc_pca
        // gate at conn_far_cut (byte-identical).  When > conn_far_cut, a drift-gated
        // cathode crosser (already past the drift_cut + cathode_x_cut + same-flash
        // gates, i.e. a same-depth opposite-TPC pair whose tips end at the cathode)
        // is accepted on PCA collinearity with the connection-vector bound RELAXED to
        // this looser value.  Rationale: for a real crosser the two truncated cathode
        // tips are displaced transversely by space charge, so the p1->p2 connection
        // vector is a poor (SCE-noisy) track-direction estimate (observed 34-69 deg on
        // confirmed PDVD crossers) while the cluster PCA stays tight (<10 deg); a full
        // drop would readmit perpendicular-connection parallel cosmics (~90 deg), so we
        // relax rather than drop.
        crosser_conn_relax_ = get(config, "crosser_conn_relax", 0.0);  // degrees
        // default crosser_pca_angle == 0: the both-long PCA branch keeps its tt_pca
        // collinearity bound at angle_cut (byte-identical).  When > angle_cut, a
        // drift-gated crosser is accepted with the tt_pca bound raised to this value:
        // a genuine crosser can bend (delta ray / SCE curvature) so its two halves'
        // GLOBAL PCA axes differ by 10-15 deg while remaining one track; the QL-pin
        // truth shows real crossers reach ttP~18 deg (p90) whereas coincidences sit
        // at ttP p50~24 / p90~51, so raising the bound to ~15 recovers bent crossers
        // without admitting the higher-ttP coincidences.
        crosser_pca_angle_ = get(config, "crosser_pca_angle", 0.0);    // degrees
        // default cathode_band_dis == 0: use only the global closest 3D point pair
        // (byte-identical).  When > 0, if that global pair hard-gates (a tip far from
        // the cathode or a large drift separation -- as happens for a long inclined
        // crosser whose global closest approach falls mid-track), retry with the
        // closest approach RESTRICTED to points within this distance of the cathode
        // plane, then apply the SAME gates to the near-cathode pair.  Additive: only
        // ever converts a hard-gate reject into a candidate that must still pass the
        // collinearity/connection gates; the legacy pair and all accepts are unchanged.
        cathode_band_dis_ = get(config, "cathode_band_dis", 0.0);
    }
    virtual Configuration default_configuration() const {
        Configuration cfg;
        return cfg;
    }

    void visit(Ensemble& ensemble) const {
        auto& live = *ensemble.with_name("live").at(0);
        clustering_cathode_connect(live, m_scope, drift_cut_, dis_cut_, max_dis_, angle_cut_,
                                   conn_far_cut_, cathode_x_, cathode_x_cut_, hough_radius_,
                                   min_length_, min_length_short_, short_dir_len_, conn_short_cut_,
                                   tip_touch_cut_, tip_touch_angle_cut_, use_flash_t0_, flash_t0_window_,
                                   crosser_conn_relax_, crosser_pca_angle_, cathode_band_dis_);
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
    double min_length_short_{10*units::cm};
    double short_dir_len_{0.0};
    double conn_short_cut_{30.0};
    double tip_touch_cut_{0.0};
    double tip_touch_angle_cut_{10.0};
    bool use_flash_t0_{true};
    double flash_t0_window_{80*units::ns};
    double crosser_conn_relax_{0.0};
    double crosser_pca_angle_{0.0};
    double cathode_band_dis_{0.0};
};


// Fold an angle (deg, 0..180) about 180 -> collinearity (0 = parallel or anti-parallel).
static inline double collinear_deg(double angle_rad)
{
    double a = angle_rad / 3.1415926 * 180.0;
    return std::min(a, 180.0 - a);
}

// Closest 3D approach between the two clusters RESTRICTED to points within `band` of
// the cathode plane.  Find_Closest_Points returns the GLOBAL closest pair, which for a
// long inclined crosser can fall mid-track (far from the cathode); this finds the pair
// at the two cathode tips instead.  Brute force over the (small) near-cathode subsets;
// only invoked when the global pair hard-gated and the knob is on.  Returns 1e9 if
// either cluster has < 3 points in the band.
static double cathode_band_closest(const Cluster& c1, const Cluster& c2,
                                   double cathode_x, double band,
                                   geo_point_t& q1, geo_point_t& q2)
{
    std::vector<geo_point_t> b1, b2;
    const int n1 = c1.npoints();
    for (int i = 0; i < n1; ++i) {
        geo_point_t p = c1.point3d(i);
        if (std::fabs(p.x() - cathode_x) < band) b1.push_back(p);
    }
    const int n2 = c2.npoints();
    for (int i = 0; i < n2; ++i) {
        geo_point_t p = c2.point3d(i);
        if (std::fabs(p.x() - cathode_x) < band) b2.push_back(p);
    }
    if (b1.size() < 3 || b2.size() < 3) return 1e9;
    double best = 1e9;
    for (const auto& a : b1) {
        for (const auto& b : b2) {
            const double dx = a.x() - b.x(), dy = a.y() - b.y(), dz = a.z() - b.z();
            const double d = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (d < best) { best = d; q1 = a; q2 = b; }
        }
    }
    return best;
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
    double hough_radius,
    double short_dir_len,
    double conn_short_cut,
    double tip_touch_cut,
    double tip_touch_angle_cut,
    double crosser_conn_relax,
    double crosser_pca_angle,
    double cathode_band_dis)
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
    // --- CATHODE_CONNECT_DEBUG instrumentation (removable; env-gated) ---
    static const bool cc_dbg = std::getenv("CATHODE_CONNECT_DEBUG") != nullptr;
    if (cc_dbg && dis < 10*units::cm && std::min(length_1, length_2) >= 10*units::cm) {
        std::fprintf(stderr,
            "[cc] c%d<->c%d dis=%.2f p1=(%.2f,%.1f,%.1f) p2=(%.2f,%.1f,%.1f) apa=%d/%d driftsep=%.2f len=%.1f/%.1f t0=%.3f/%.3fus dt0=%.3fus\n",
            (int)cluster1.ident(), (int)cluster2.ident(), dis/units::cm,
            p1.x()/units::cm, p1.y()/units::cm, p1.z()/units::cm,
            p2.x()/units::cm, p2.y()/units::cm, p2.z()/units::cm,
            wpid_p1.apa(), wpid_p2.apa(), std::fabs(p1.x()-p2.x())/units::cm,
            length_1/units::cm, length_2/units::cm,
            cluster1.get_cluster_t0()/units::us, cluster2.get_cluster_t0()/units::us,
            std::fabs(cluster1.get_cluster_t0()-cluster2.get_cluster_t0())/units::us);
    }
    if (wpid_p1.apa() == wpid_p2.apa()) return false;

    // --- CATHODE-BAND closest-approach retry (default OFF: cathode_band_dis == 0).
    //     If the GLOBAL closest pair hard-gates (a tip far from the cathode, or a large
    //     drift separation -- a long inclined crosser whose global closest approach
    //     falls mid-track), retry with the closest approach restricted to the cathode
    //     band and adopt it, so the gates below re-test the near-cathode tips.  Additive:
    //     can only rescue a would-be hard-gate reject; an already-passing pair keeps its
    //     global closest.  Both band points are cross-TPC by construction (each is its
    //     own cluster's near-cathode point), so the apa relation above still holds. ---
    if (cathode_band_dis > 0) {
        const bool global_hardgated =
            std::fabs(p1.x() - cathode_x) >= cathode_x_cut ||
            std::fabs(p2.x() - cathode_x) >= cathode_x_cut ||
            std::fabs(p1.x() - p2.x())    >= drift_cut;
        if (global_hardgated) {
            geo_point_t q1, q2;
            double bd = cathode_band_closest(cluster1, cluster2, cathode_x,
                                             cathode_band_dis, q1, q2);
            if (bd < max_dis) { p1 = q1; p2 = q2; dis = bd; }
        }
    }

    // --- CATHODE_CONNECT_DEBUG per-return tracer (removable; env-gated) ---
    double _ttH = -1, _ttP = -1, _ccb = -1;
    auto DBG = [&](const char* why, bool acc) -> bool {
        if (cc_dbg && dis < max_dis) {
            std::fprintf(stderr,
                "[ccx] c%d<->c%d %-16s dis=%.2f dX=%.2f tip=%.2f/%.2f ttH=%.1f ttP=%.1f cc=%.1f len=%.0f/%.0f -> %s\n",
                (int)cluster1.ident(), (int)cluster2.ident(), why, dis/units::cm,
                std::fabs(p1.x()-p2.x())/units::cm,
                std::fabs(p1.x()-cathode_x)/units::cm, std::fabs(p2.x()-cathode_x)/units::cm,
                _ttH, _ttP, _ccb, length_1/units::cm, length_2/units::cm,
                acc ? "ACCEPT" : "reject");
        }
        return acc;
    };

    // --- CC_FEATURE_DUMP: full labeled feature vector for EVERY cross-APA candidate
    //     pair (removable; env-gated).  Computes all directions unconditionally (local
    //     copies, does not touch the real logic below) + p1/p2 world coords so an
    //     offline join can label each pair by QL-pin truth (point-membership).  Used
    //     to data-drive a geometric separator of real crossers vs coincidences. ---
    static const bool cc_feat = std::getenv("CC_FEATURE_DUMP") != nullptr;
    if (cc_feat && dis < max_dis) {
        geo_point_t fd1 = cluster1.vhough_transform(p1, hough_radius);
        geo_point_t fd2 = cluster2.vhough_transform(p2, hough_radius);
        double f_ttH = collinear_deg(fd1.angle(fd2));
        double f_ttP = 999, f_ccH = 999, f_ccP = 999;
        const auto& fax1 = cluster1.get_pca().axis;
        const auto& fax2 = cluster2.get_pca().axis;
        geo_point_t fconn(p1.x()-p2.x(), p1.y()-p2.y(), p1.z()-p2.z());
        f_ccH = collinear_deg(fconn.angle(fd1));
        if (!fax1.empty() && !fax2.empty()) {
            geo_point_t fp1(fax1.at(0).x(), fax1.at(0).y(), fax1.at(0).z());
            geo_point_t fp2(fax2.at(0).x(), fax2.at(0).y(), fax2.at(0).z());
            f_ttP = collinear_deg(fp1.angle(fp2));
            f_ccP = std::min(collinear_deg(fconn.angle(fp1)), collinear_deg(fconn.angle(fp2)));
        }
        std::fprintf(stderr,
            "[feat] dis=%.2f dX=%.2f tip1=%.2f tip2=%.2f ttH=%.1f ttP=%.1f ccH=%.1f ccP=%.1f "
            "len1=%.1f len2=%.1f p1=%.1f,%.1f,%.1f p2=%.1f,%.1f,%.1f apa=%d/%d\n",
            dis/units::cm, std::fabs(p1.x()-p2.x())/units::cm,
            std::fabs(p1.x()-cathode_x)/units::cm, std::fabs(p2.x()-cathode_x)/units::cm,
            f_ttH, f_ttP, f_ccH, f_ccP, length_1/units::cm, length_2/units::cm,
            p1.x()/units::cm, p1.y()/units::cm, p1.z()/units::cm,
            p2.x()/units::cm, p2.y()/units::cm, p2.z()/units::cm,
            wpid_p1.apa(), wpid_p2.apa());
    }

    // (2) both closest points must end at the cathode plane.
    if (std::fabs(p1.x() - cathode_x) >= cathode_x_cut) return DBG("tipx1", false);
    if (std::fabs(p2.x() - cathode_x) >= cathode_x_cut) return DBG("tipx2", false);

    // (4) the two points must be at the same DRIFT depth (always tight): the only
    //     thing between the two halves is the ~1.5 cm cathode gap plus a drift-x
    //     calibration residual (observed up to ~4.1 cm in data).
    if (std::fabs(p1.x() - p2.x()) >= drift_cut) return DBG("drift", false);

    // outer 3D distance ceiling.
    if (dis >= max_dis) return DBG("maxdis", false);

    // (1) local track directions, from the Hough transform at the closest points.
    geo_point_t dir1 = cluster1.vhough_transform(p1, hough_radius);
    geo_point_t dir2 = cluster2.vhough_transform(p2, hough_radius);
    double tt_hough = collinear_deg(dir1.angle(dir2));
    _ttH = tt_hough;

    // cluster PCA principal axes (an alternative, global direction estimate used by both
    // the short-stub branch and the far regime; empty/invalid axis is guarded).
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

    if (dis < dis_cut) {
        // CLOSE regime: the p1->p2 connection vector is dominated by the drift-x offset
        //   (the calibration artifact, ~along the drift axis), not the track, so the
        //   generic passes' connection-alignment test rejects these crossers -- the hole
        //   this pass fills.  Accept on the local half-track collinearity alone.
        if (tt_hough < angle_cut) return DBG("close_primary", true);

        // BOTH-LONG PCA fallback (default OFF: short_dir_len == 0).  When BOTH halves
        //   are long enough for their cluster PCA to be reliable (>= short_dir_len) but
        //   the local Hough at the dense cathode tips is noisy, the Hough collinearity
        //   above is untrustworthy (e.g. 183096: tt_hough=24 deg while tt_pca=1.1 deg).
        //   Accept on the FAR-regime test applied to this close pair: the PCA principal
        //   axes are collinear (tt_pca < angle_cut) AND the p1->p2 connection vector
        //   continues the track (cc_pca < conn_far_cut).  The connection-alignment term
        //   is essential -- two distinct parallel cosmics that merely pass close at the
        //   cathode are also PCA-collinear, but their connection is ~perpendicular, so
        //   cc_pca rejects them (a PCA-collinearity-only test would wrongly merge them).
        if (short_dir_len > 0 && have_pca &&
            std::min(length_1, length_2) >= short_dir_len) {
            double tt_pca = collinear_deg(pca1.angle(pca2));
            geo_point_t conn(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
            double cc_pca = std::min(collinear_deg(conn.angle(pca1)),
                                     collinear_deg(conn.angle(pca2)));
            _ttP = tt_pca; _ccb = cc_pca;
            // TIP-TOUCH relaxation (default OFF: tip_touch_cut == 0).  When the two
            //   halves nearly touch (dis < tip_touch_cut, e.g. a genuine crosser whose
            //   cathode tips meet within ~1 cm), the p1->p2 connection vector spans only
            //   that ~1 cm gap and is dominated by sub-cm transverse jitter -- cc_pca is
            //   ~perpendicular (>=75 deg) even for a real crosser, so it is uninformative
            //   and wrongly vetoes the merge.  There we drop the cc_pca term and accept on
            //   PCA collinearity alone; the small closest distance itself is the evidence
            //   the tips meet (a transversely-offset parallel cosmic has a larger dis and
            //   still goes through the cc_pca gate below).  Gate pairs on flash-time
            //   coincidence (use_flash_t0) so only same-flash tips are relaxed this way.
            const bool tip_touch = (tip_touch_cut > 0 && dis < tip_touch_cut);
            // 6cm-CATHODE crosser relaxation (default OFF: crosser_conn_relax == 0).
            //   Reaching this branch already means a drift-gated crosser: opposite TPC,
            //   both tips within cathode_x_cut of the cathode, same DRIFT depth
            //   (|dX| < drift_cut), same flash.  The two half-tracks stop at the two
            //   FACES of the 6cm-thick cathode, so the p1->p2 connection spans ~6cm of
            //   pure cathode geometry (drift gap + SCE transverse shift), not track --
            //   an SCE-noisy direction estimate.  Relax the cc_pca bound from
            //   conn_far_cut to the looser crosser_conn_relax (still rejecting the
            //   ~perpendicular connection of two distinct parallel cosmics).
            const double conn_bound =
                (crosser_conn_relax > conn_far_cut) ? crosser_conn_relax : conn_far_cut;
            // crosser_pca_angle raises the tt_pca collinearity bound for bent crossers
            // (default 0 => angle_cut => byte-identical).
            const double pca_ang =
                (crosser_pca_angle > angle_cut) ? crosser_pca_angle : angle_cut;
            // global-PCA acceptance (original; cc_pca dropped when tips touch, else
            // tested against conn_bound == conn_far_cut unless the crosser relax is on).
            const bool acc_pca = (tt_pca < pca_ang && (tip_touch || cc_pca < conn_bound));
            // LOCAL-Hough acceptance (default OFF: tip_touch_angle_cut == angle_cut, so
            //   this reduces to the primary close test that already ran).  When tips
            //   touch, a curved half can inflate its GLOBAL PCA above angle_cut while its
            //   LOCAL cathode arm still continues the partner (e.g. 29107 evt991 cl26:
            //   global tt_pca=15 deg, but Hough tt_hough=10 deg).  Accept on the local
            //   Hough with the relaxed bound.  Hough is charge-weighted so it correctly
            //   REJECTS oblique/perpendicular touchers (tt_hough 49/90 deg) that an
            //   unweighted local SVD would wrongly pass; the strong spatial + same-flash
            //   gates make a non-crossing coincidence at the cathode unlikely.
            const bool acc_hough = (tip_touch && tt_hough < tip_touch_angle_cut);
            if (cc_dbg && dis < 10*units::cm) {
                std::fprintf(stderr,
                    "[cc]   c%d<->c%d CLOSE both-long: tt_hough=%.1f tt_pca=%.1f cc_pca=%.1f tip_touch=%d(cut=%.1f,ang=%.1f) -> %s\n",
                    (int)cluster1.ident(), (int)cluster2.ident(), tt_hough, tt_pca, cc_pca,
                    (int)tip_touch, tip_touch_cut/units::cm, tip_touch_angle_cut,
                    (acc_pca || acc_hough) ? "ACCEPT" : "reject");
            }
            if (acc_pca || acc_hough) return DBG("close_bothlong", true);
        }

        // SHORT-STUB prolongation (default OFF: short_dir_len == 0).  When one half is
        //   too short for its own direction to be trusted (a gap-splintered cathode
        //   stub), the collinearity test above is unreliable.  Instead use only the LONG
        //   (anchor) half's reliable direction: the anchor->stub connection vector must
        //   continue the anchor track (within conn_short_cut).  The anchor's Hough (at its
        //   closest point) OR cluster PCA axis may satisfy it.  Require EXACTLY one short
        //   member (< short_dir_len) and the other a genuinely long anchor (>= short_dir_len,
        //   not merely >= min_length) -- a mid-length "anchor" has an unreliable direction
        //   too, so trusting its extrapolation would test conn against junk.  A coincidental
        //   short cluster offset transversely from the anchor's extrapolation has a
        //   ~perpendicular connection and is rejected.
        if (short_dir_len > 0 &&
            std::min(length_1, length_2) < short_dir_len &&
            std::max(length_1, length_2) >= short_dir_len) {
            const bool one_is_anchor = (length_1 >= length_2);
            const geo_point_t& anchor_pt    = one_is_anchor ? p1 : p2;
            const geo_point_t& stub_pt      = one_is_anchor ? p2 : p1;
            const geo_point_t& anchor_hough = one_is_anchor ? dir1 : dir2;
            const geo_point_t& anchor_pca   = one_is_anchor ? pca1 : pca2;
            geo_point_t conn(anchor_pt.x() - stub_pt.x(),
                             anchor_pt.y() - stub_pt.y(),
                             anchor_pt.z() - stub_pt.z());
            double cc_hough = collinear_deg(conn.angle(anchor_hough));
            double cc_pca = have_pca ? collinear_deg(conn.angle(anchor_pca)) : 999.0;
            _ccb = std::min(cc_hough, cc_pca);
            return DBG("close_shortstub", std::min(cc_hough, cc_pca) < conn_short_cut);
        }
        return DBG("close_fallthrough", false);
    }

    // FAR regime: the two halves are offset within the cathode plane (large transverse
    //   separation from a shallow-angle crosser travelling along the cathode).  Two
    //   things change vs CLOSE: (a) the local Hough direction can be unreliable when one
    //   half is a dense blob near the cathode, so we ADD the cluster PCA principal axis
    //   as an alternative direction estimate (Hough OR PCA may satisfy collinearity);
    //   (b) the now-long p1->p2 connection vector is meaningful, so we require it to
    //   align with the track (PCA or Hough) -- this rejects parallel-offset cosmics,
    //   whose connection is ~perpendicular (>=50 deg), while passing a real crosser.
    double tt_pca = have_pca ? collinear_deg(pca1.angle(pca2)) : 999.0;
    _ttP = tt_pca;
    // track-track collinear: Hough OR (additional) PCA.
    if (tt_hough >= angle_cut && tt_pca >= angle_cut) return DBG("far_notcollinear", false);

    // connection vector aligned with the track: Hough OR (additional) PCA.
    geo_point_t conn(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
    double cc_hough = collinear_deg(conn.angle(dir1));
    double cc_pca = have_pca ? std::min(collinear_deg(conn.angle(pca1)),
                                        collinear_deg(conn.angle(pca2))) : 999.0;
    _ccb = std::min(cc_hough, cc_pca);
    if (cc_hough >= conn_far_cut && cc_pca >= conn_far_cut) return DBG("far_conn", false);

    return DBG("far_accept", true);
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
    double min_length_short,
    double short_dir_len,
    double conn_short_cut,
    double tip_touch_cut,
    double tip_touch_angle_cut,
    bool use_flash_t0,
    double flash_t0_window,
    double crosser_conn_relax,
    double crosser_pca_angle,
    double cathode_band_dis)
{
    // prepare graph ... (same skeleton as the other merge passes)
    typedef cluster_connectivity_graph_t Graph;
    Graph g;
    std::unordered_map<int, int> ilive2desc;
    std::unordered_map<const Cluster*, int> map_cluster_index;
    auto live_clusters = live_grouping.children();

    // Cross-TPC connection only makes sense for clusters coincident in flash time;
    // unmatched clusters get unique singleton groups (never linked).  When flash
    // matching is unavailable (use_flash_t0 false), skip the gate entirely.
    std::map<const Cluster*, int> flash_t0_group;
    if (use_flash_t0) {
        flash_t0_group = assign_flash_t0_groups(live_clusters, flash_t0_window);
    }

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
        if (cluster_1->get_length() < min_length_short) continue;
        for (size_t j = i + 1; j < live_clusters.size(); j++) {
            auto cluster_2 = live_clusters.at(j);
            if (!cluster_2->get_scope_filter(scope)) continue;
            if (cluster_2->get_length() < min_length_short) continue;
            // Asymmetric length gate: the longer member must be a real anchor track
            // (>= min_length); the shorter only needs min_length_short.  This admits a
            // short bridge fragment that attaches to a long half (a gap-splintered
            // crosser) while forbidding short<->short pairs.  With min_length_short ==
            // min_length (the default) this is the original "both >= min_length" gate.
            if (std::max(cluster_1->get_length(), cluster_2->get_length()) < min_length) continue;
            if (use_flash_t0 && flash_t0_group.at(cluster_1) != flash_t0_group.at(cluster_2)) {
                static const bool cc_dbg2 = std::getenv("CATHODE_CONNECT_DEBUG") != nullptr;
                if (cc_dbg2 && std::min(cluster_1->get_length(), cluster_2->get_length()) >= 15*units::cm) {
                    std::fprintf(stderr, "[ccflash] c%d<->c%d FLASH-GATE reject len=%.0f/%.0f t0=%.3f/%.3fus\n",
                        (int)cluster_1->ident(), (int)cluster_2->ident(),
                        cluster_1->get_length()/units::cm, cluster_2->get_length()/units::cm,
                        cluster_1->get_cluster_t0()/units::us, cluster_2->get_cluster_t0()/units::us);
                }
                continue;
            }
            if (is_cathode_crossing_pair(*cluster_1, *cluster_2,
                                         cluster_1->get_length(), cluster_2->get_length(),
                                         drift_cut, dis_cut, max_dis, angle_cut, conn_far_cut,
                                         cathode_x, cathode_x_cut, hough_radius,
                                         short_dir_len, conn_short_cut, tip_touch_cut,
                                         tip_touch_angle_cut, crosser_conn_relax, crosser_pca_angle,
                                         cathode_band_dis)) {
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
