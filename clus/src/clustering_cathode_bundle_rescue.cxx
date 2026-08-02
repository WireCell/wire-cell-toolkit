// SBND cathode BUNDLE rescue (default-OFF, isolated, removable).
//
// NOT a WCP port -- a new algorithm.  It patches, on the TPC side, a known SBND
// light-reconstruction defect: an early flash opens a readout window (~8 us) that
// can ABSORB a later flash.  A beam interaction whose track crosses the central
// cathode scintillates in both drift volumes and should yield a beam-coincident
// flash on each side; when one side's flash is absorbed, that side's charge half
// gets Q/L-matched to a DIFFERENT flash (2-13 us away, outside the beam window)
// with a perfectly good light pattern.  The two halves then sit in different
// bundles, cathode_connect's flash gate correctly refuses to join them, and the
// beam-window gate hands the PR chain only half the track (doc pr/12 sec 6:
// 7 of 12 "fit stops at the cathode" candidates; sbnd_xin/docs/pr/14).
//
// What it does, per event, after cathode_connect and BEFORE examine_bundles:
//   1. find pairs (K_beam, K_far): K_beam matched in the beam window with a tip
//      at the cathode; K_far in a DIFFERENT flash bundle within
//      [-rescue_t0_early, +rescue_t0_late] of K_beam's t0, outside the beam
//      window; the pair passes the cathode-crossing geometric test (duplicated
//      from clustering_cathode_connect.cxx -- see note below) with K_far's x
//      evaluated under the destination-T0 hypothesis.
//   2. decide WHICH flash the joined crosser belongs to from the bundle
//      structure: a/b = K_beam's length vs the longest OTHER member of its
//      bundle; c/d = same for K_far.  The bundle that is dominated by its
//      crosser half keeps the crosser; ties go to the longer half ("small
//      merges to large").  Both directions occur in the hand-scanned truth set
//      (4 into the beam bundle, 4 out of it -- doc pr/14 sec 3).
//   3. merge the pair into ONE cluster (merge_clusters, exactly as
//      cathode_connect) and re-stamp the merged cluster's flash identity
//      (cluster_t0 / flash / matched_flash_gid / lm_flag) from the destination
//      bundle, then re-materialize the corrected coordinates under the new T0.
//      Merging BEFORE examine_bundles is what makes the crosser one
//      flash-collapse member, so the PR job's unmerge_bundle keeps it whole
//      (same invariant as cathode_connect; doc 53).
//   4. (rescue_unmatched, default OFF; sbnd_xin/docs/pr/17) a second pass
//      adopts NON-MATCHED clusters: when the far half is one whole cluster
//      (e.g. after the pr/15 vertex veto), no wrong flash may claim it and
//      Q/L matching leaves it flashless -- invisible downstream.  If it
//      geometrically continues a beam-window cluster across the cathode it
//      is merged into the beam bundle (forced direction: it has no flash of
//      its own).  Its geometry is evaluated in the RAW scope (its x_t0cor
//      carries the sentinel T0).
//
// Isolation: this file is self-contained (all helpers file-static); the
// geometric pair test is DUPLICATED from clustering_cathode_connect.cxx
// (M10 fork-by-duplication -- the production connector stays byte-identical).
// When the light reconstruction fixes the absorbing-window defect upstream,
// retire this pass by flipping its single config flag (cathode_rescue_on in
// cfg/pgrapher/experiment/sbnd/clus.jsonnet).
//
// Inert by default: with beam_window_low == beam_window_high (C++ default 0/0)
// no cluster is "in beam" and visit() is a no-op even if mistakenly inserted.

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/Logger.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>   // std::getenv  (CATHODE_RESCUE_DEBUG instrumentation, removable)
#include <cstdio>    // std::fprintf

class ClusteringCathodeBundleRescue;
WIRECELL_FACTORY(ClusteringCathodeBundleRescue, ClusteringCathodeBundleRescue,
                 WireCell::INamed, WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

namespace {

// All geometry/window parameters of the pass, bundled to keep signatures sane.
struct CbrParams {
    // beam window (half-open [low, high), same convention as CreateSteinerGraph).
    double beam_window_low{0.0};
    double beam_window_high{0.0};
    // asymmetric bundle search window around the beam cluster's t0:
    // accept t0_far - t0_beam in [-rescue_t0_early, +rescue_t0_late].
    double rescue_t0_early{8*units::us};
    double rescue_t0_late{13*units::us};
    bool require_far_out_of_beam{true};
    // geometric pair test -- duplicated knob set of clustering_cathode_connect.cxx.
    double drift_cut{5*units::cm};
    double dis_cut{5*units::cm};
    double max_dis{25*units::cm};
    double angle_cut{10.0};        // degrees
    double conn_far_cut{30.0};     // degrees
    double cathode_x{0.0};
    double cathode_x_cut{3.5*units::cm};
    double hough_radius{20*units::cm};
    double min_length{10*units::cm};
    double min_length_short{10*units::cm};
    double short_dir_len{0.0};
    double conn_short_cut{30.0};   // degrees
    double tip_touch_cut{0.0};
    double tip_touch_angle_cut{10.0};
    double cathode_band_dis{0.0};
    // unmatched-cluster adoption pass (default OFF; sbnd_xin/docs/pr/17).
    // A cluster the Q/L matching left with NO flash at all (56463 veto-ON
    // mode: the true flash was absorbed and no wrong flash claims the whole
    // cluster) is invisible downstream.  If it geometrically continues a
    // beam-window cluster across the cathode, adopt it into the beam bundle.
    bool rescue_unmatched{false};
    double unmatched_min_length{30*units::cm};
    int unmatched_min_npts{200};
};

// Fold an angle (deg, 0..180) about 180 -> collinearity (0 = parallel or anti-parallel).
// (duplicated: clustering_cathode_connect.cxx collinear_deg)
inline double collinear_deg(double angle_rad)
{
    double a = angle_rad / 3.1415926 * 180.0;
    return std::min(a, 180.0 - a);
}

// Closest 3D approach RESTRICTED to points within `band` of the cathode plane.
// (duplicated: clustering_cathode_connect.cxx cathode_band_closest)
double cathode_band_closest(const Cluster& c1, const Cluster& c2,
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

// Cathode-crossing pair test, duplicated from clustering_cathode_connect.cxx
// is_cathode_crossing_pair (see that file's long header comment for the full cut
// rationale) with ONE addition: `xshift2` translates cluster2's coordinates in x
// before every cut.  cluster2 (the far half) carries a WRONG matched T0, so its
// stored x_t0cor is off by v_drift * dt0 (up to ~2 cm at 12 us); xshift2 puts it
// in the destination-T0 frame so the tip/drift/distance cuts are unbiased.  The
// relative pair geometry is the SAME in either half's frame (the two TPCs drift
// in opposite x), so one shift serves both possible move directions.  Residual
// caveat: the closest-pair SELECTION (Find_Closest_Points / the band retry) runs
// on unshifted coordinates; a <= 2 cm rigid x-shift does not change which tips
// face each other at the cathode.  Local Hough / PCA directions are translation
// invariant and are evaluated at the ORIGINAL points.
bool is_cathode_crossing_pair(
    const Cluster& cluster1,
    const Cluster& cluster2,
    double length_1,
    double length_2,
    double xshift2,
    const CbrParams& P)
{
    geo_point_t p1;
    geo_point_t p2;
    // Find_Closest_Points ignores its length_cut argument (it returns the global
    // closest 3D point pair); we apply our own cuts below.
    double dis = WireCell::Clus::Facade::Find_Closest_Points(cluster1, cluster2,
                                                             length_1, length_2,
                                                             P.max_dis, p1, p2);

    // the two closest points must sit in DIFFERENT TPCs (this is what makes the
    // pass incapable of acting within a single TPC).
    auto wpid_p1 = cluster1.wpid(p1);
    auto wpid_p2 = cluster2.wpid(p2);
    if (wpid_p1.apa() == wpid_p2.apa()) return false;

    // CATHODE-BAND closest-approach retry (default OFF: cathode_band_dis == 0);
    // see clustering_cathode_connect.cxx for rationale.  Runs on unshifted
    // coordinates; the shift is applied to whichever pair is adopted.
    if (P.cathode_band_dis > 0) {
        const bool global_hardgated =
            std::fabs(p1.x() - P.cathode_x) >= P.cathode_x_cut ||
            std::fabs(p2.x() + xshift2 - P.cathode_x) >= P.cathode_x_cut ||
            std::fabs(p1.x() - p2.x() - xshift2) >= P.drift_cut;
        if (global_hardgated) {
            geo_point_t q1, q2;
            double bd = cathode_band_closest(cluster1, cluster2, P.cathode_x,
                                             P.cathode_band_dis, q1, q2);
            if (bd < P.max_dis) { p1 = q1; p2 = q2; }
        }
    }

    // Destination-T0 hypothesis: shift cluster2's tip in x; keep the original
    // point for the (translation-invariant) local Hough direction below.
    const geo_point_t p2_orig = p2;
    geo_point_t p2s(p2.x() + xshift2, p2.y(), p2.z());
    dis = std::sqrt((p1.x() - p2s.x())*(p1.x() - p2s.x()) +
                    (p1.y() - p2s.y())*(p1.y() - p2s.y()) +
                    (p1.z() - p2s.z())*(p1.z() - p2s.z()));

    // --- CATHODE_RESCUE_DEBUG per-return tracer (removable; env-gated) ---
    static const bool cbr_dbg = std::getenv("CATHODE_RESCUE_DEBUG") != nullptr;
    double _ttH = -1, _ttP = -1, _ccb = -1;
    auto DBG = [&](const char* why, bool acc) -> bool {
        if (cbr_dbg && dis < P.max_dis) {
            std::fprintf(stderr,
                "[cbrx] c%d<->c%d %-16s dis=%.2f dX=%.2f tip=%.2f/%.2f xsh=%.2f ttH=%.1f ttP=%.1f cc=%.1f len=%.0f/%.0f -> %s\n",
                (int)cluster1.ident(), (int)cluster2.ident(), why, dis/units::cm,
                std::fabs(p1.x()-p2s.x())/units::cm,
                std::fabs(p1.x()-P.cathode_x)/units::cm, std::fabs(p2s.x()-P.cathode_x)/units::cm,
                xshift2/units::cm, _ttH, _ttP, _ccb,
                length_1/units::cm, length_2/units::cm,
                acc ? "ACCEPT" : "reject");
        }
        return acc;
    };

    // both closest points must end at the cathode plane.
    if (std::fabs(p1.x() - P.cathode_x) >= P.cathode_x_cut) return DBG("tipx1", false);
    if (std::fabs(p2s.x() - P.cathode_x) >= P.cathode_x_cut) return DBG("tipx2", false);

    // same DRIFT depth (always tight).
    if (std::fabs(p1.x() - p2s.x()) >= P.drift_cut) return DBG("drift", false);

    // outer 3D distance ceiling.
    if (dis >= P.max_dis) return DBG("maxdis", false);

    // local track directions, from the Hough transform at the closest points.
    geo_point_t dir1 = cluster1.vhough_transform(p1, P.hough_radius);
    geo_point_t dir2 = cluster2.vhough_transform(p2_orig, P.hough_radius);
    double tt_hough = collinear_deg(dir1.angle(dir2));
    _ttH = tt_hough;

    // cluster PCA principal axes (guarded).
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

    if (dis < P.dis_cut) {
        // CLOSE regime: accept on local half-track collinearity alone.
        if (tt_hough < P.angle_cut) return DBG("close_primary", true);

        // BOTH-LONG PCA fallback (OFF unless short_dir_len > 0); tip-touch
        // relaxation as in the connector.  The crosser_conn_relax /
        // crosser_pca_angle knobs of the connector are NOT duplicated here
        // (PDVD-specific; SBND leaves them off).
        if (P.short_dir_len > 0 && have_pca &&
            std::min(length_1, length_2) >= P.short_dir_len) {
            double tt_pca = collinear_deg(pca1.angle(pca2));
            geo_point_t conn(p1.x() - p2s.x(), p1.y() - p2s.y(), p1.z() - p2s.z());
            double cc_pca = std::min(collinear_deg(conn.angle(pca1)),
                                     collinear_deg(conn.angle(pca2)));
            _ttP = tt_pca; _ccb = cc_pca;
            const bool tip_touch = (P.tip_touch_cut > 0 && dis < P.tip_touch_cut);
            const bool acc_pca = (tt_pca < P.angle_cut && (tip_touch || cc_pca < P.conn_far_cut));
            const bool acc_hough = (tip_touch && tt_hough < P.tip_touch_angle_cut);
            if (acc_pca || acc_hough) return DBG("close_bothlong", true);
        }

        // SHORT-STUB prolongation (OFF unless short_dir_len > 0).
        if (P.short_dir_len > 0 &&
            std::min(length_1, length_2) < P.short_dir_len &&
            std::max(length_1, length_2) >= P.short_dir_len) {
            const bool one_is_anchor = (length_1 >= length_2);
            const geo_point_t& anchor_pt    = one_is_anchor ? p1 : p2s;
            const geo_point_t& stub_pt      = one_is_anchor ? p2s : p1;
            const geo_point_t& anchor_hough = one_is_anchor ? dir1 : dir2;
            const geo_point_t& anchor_pca   = one_is_anchor ? pca1 : pca2;
            geo_point_t conn(anchor_pt.x() - stub_pt.x(),
                             anchor_pt.y() - stub_pt.y(),
                             anchor_pt.z() - stub_pt.z());
            double cc_hough = collinear_deg(conn.angle(anchor_hough));
            double cc_pca = have_pca ? collinear_deg(conn.angle(anchor_pca)) : 999.0;
            _ccb = std::min(cc_hough, cc_pca);
            return DBG("close_shortstub", std::min(cc_hough, cc_pca) < P.conn_short_cut);
        }
        return DBG("close_fallthrough", false);
    }

    // FAR regime.
    double tt_pca = have_pca ? collinear_deg(pca1.angle(pca2)) : 999.0;
    _ttP = tt_pca;
    if (tt_hough >= P.angle_cut && tt_pca >= P.angle_cut) return DBG("far_notcollinear", false);

    geo_point_t conn(p1.x() - p2s.x(), p1.y() - p2s.y(), p1.z() - p2s.z());
    double cc_hough = collinear_deg(conn.angle(dir1));
    double cc_pca = have_pca ? std::min(collinear_deg(conn.angle(pca1)),
                                        collinear_deg(conn.angle(pca2))) : 999.0;
    _ccb = std::min(cc_hough, cc_pca);
    if (cc_hough >= P.conn_far_cut && cc_pca >= P.conn_far_cut) return DBG("far_conn", false);

    return DBG("far_accept", true);
}

} // anonymous namespace

class ClusteringCathodeBundleRescue : public IConfigurable,
                                      public Clus::IEnsembleVisitor,
                                      public Aux::Logger,
                                      private NeedDV, private NeedPCTS, private NeedScope {
public:
    ClusteringCathodeBundleRescue()
        : Aux::Logger("CathodeBundleRescue", "clus")
    {}
    virtual ~ClusteringCathodeBundleRescue() {}

    void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedPCTS::configure(config);
        NeedScope::configure(config);

        // Window knobs.  beam_window defaults to an EMPTY window: no beam
        // cluster, no candidate pair, visit() is a no-op unless configured.
        p_.beam_window_low        = get(config, "beam_window_low", p_.beam_window_low);
        p_.beam_window_high       = get(config, "beam_window_high", p_.beam_window_high);
        p_.rescue_t0_early        = get(config, "rescue_t0_early", p_.rescue_t0_early);
        p_.rescue_t0_late         = get(config, "rescue_t0_late", p_.rescue_t0_late);
        p_.require_far_out_of_beam= get(config, "require_far_out_of_beam", p_.require_far_out_of_beam);

        // Geometric pair test knobs (same names and defaults as the connector).
        p_.drift_cut       = get(config, "drift_cut", p_.drift_cut);
        p_.dis_cut         = get(config, "dis_cut", p_.dis_cut);
        p_.max_dis         = get(config, "max_dis", p_.max_dis);
        p_.angle_cut       = get(config, "angle_cut", p_.angle_cut);
        p_.conn_far_cut    = get(config, "conn_far_cut", p_.conn_far_cut);
        p_.cathode_x       = get(config, "cathode_x", p_.cathode_x);
        p_.cathode_x_cut   = get(config, "cathode_x_cut", p_.cathode_x_cut);
        p_.hough_radius    = get(config, "hough_radius", p_.hough_radius);
        p_.min_length      = get(config, "min_length", p_.min_length);
        p_.min_length_short= get(config, "min_length_short", p_.min_length);
        p_.short_dir_len   = get(config, "short_dir_len", p_.short_dir_len);
        p_.conn_short_cut  = get(config, "conn_short_cut", p_.conn_short_cut);
        p_.tip_touch_cut   = get(config, "tip_touch_cut", p_.tip_touch_cut);
        p_.tip_touch_angle_cut = get(config, "tip_touch_angle_cut", p_.angle_cut);
        p_.cathode_band_dis= get(config, "cathode_band_dis", p_.cathode_band_dis);

        // Unmatched-cluster adoption pass (C++ default false => the pass is
        // absent and behavior is byte-identical to pre-knob builds).
        p_.rescue_unmatched     = get(config, "rescue_unmatched", p_.rescue_unmatched);
        p_.unmatched_min_length = get(config, "unmatched_min_length", p_.unmatched_min_length);
        p_.unmatched_min_npts   = get(config, "unmatched_min_npts", p_.unmatched_min_npts);
    }

    virtual Configuration default_configuration() const {
        Configuration cfg;
        return cfg;
    }

    void visit(Ensemble& ensemble) const {
        auto groupings = ensemble.with_name("live");
        if (groupings.empty()) return;
        rescue(*groupings.at(0));
    }

private:
    CbrParams p_;

    struct Member {           // one flash-matched cluster, snapshot per round
        Cluster* cluster;
        int ident;
        int gid;
        double t0;
        double length;
    };

    bool in_beam(double t0) const {
        return t0 >= p_.beam_window_low && t0 < p_.beam_window_high;
    }

    // x-shift that moves the far half into the destination-T0 frame:
    // x_hyp = x_t0cor(old t0) - dirx * (t0_dest - t0_old) * v_drift.
    // dirx/speed are per (apa, face) of the far tip; trigger offsets cancel in
    // the difference.  Falls back to 0 (no shift) if the speed map lacks the
    // face -- the cuts then simply run on the stored coordinates, which the
    // pr/12 measurements show still pass for the known population.
    double far_xshift(const Grouping& grouping, const Cluster& far, double dt0_dest) const {
        // any point of the far cluster gives its apa/face (single-TPC cluster).
        if (far.npoints() == 0) return 0.0;
        auto wpid = far.wpid(far.point3d(0));
        const auto speeds = grouping.get_drift_speed();
        auto ait = speeds.find(wpid.apa());
        if (ait == speeds.end()) return 0.0;
        auto fit = ait->second.find(wpid.face());
        if (fit == ait->second.end()) return 0.0;
        const double dirx = m_dv->face_dirx(WirePlaneId(kAllLayers, wpid.face(), wpid.apa()));
        return -dirx * dt0_dest * fit->second;
    }

    void rescue(Grouping& live_grouping) const {
        if (!(p_.beam_window_high > p_.beam_window_low)) return;

        // One accepted move per round, then FULL re-enumeration: merge_clusters
        // destroys both inputs and appends a fresh cluster, so every pointer,
        // index and bundle tally is stale after it.  Progress is guaranteed
        // (each round removes one cluster and the merged pair shares one gid,
        // so it can never re-qualify against itself).
        int nrounds = 0;
        while (true) {
            auto live_clusters = live_grouping.children();
            for (auto* live : live_clusters) {
                if (live->get_default_scope().hash() != m_scope.hash()) {
                    live->set_default_scope(m_scope);
                }
            }

            // Flash-matched membership snapshot (by value; ident order).
            std::vector<Member> members;
            for (auto* c : live_clusters) {
                if (c->get_scalar<int>("flash", -1) < 0) continue;
                const int gid = c->get_scalar<int>("matched_flash_gid", -1);
                if (gid < 0) continue;
                members.push_back({c, (int)c->ident(), gid,
                                   c->get_cluster_t0(), c->get_length()});
            }
            std::sort(members.begin(), members.end(),
                      [](const Member& a, const Member& b) { return a.ident < b.ident; });

            // First qualifying pair in (ident_beam, ident_far) order.
            const Member* best_beam = nullptr;
            const Member* best_far = nullptr;
            for (const auto& kb : members) {
                if (best_beam) break;
                if (!in_beam(kb.t0)) continue;
                if (!kb.cluster->get_scope_filter(m_scope)) continue;
                if (kb.length < p_.min_length_short) continue;
                for (const auto& kf : members) {
                    if (kf.gid == kb.gid) continue;
                    const double dt0 = kf.t0 - kb.t0;
                    if (dt0 < -p_.rescue_t0_early || dt0 > p_.rescue_t0_late) continue;
                    if (p_.require_far_out_of_beam && in_beam(kf.t0)) continue;
                    if (!kf.cluster->get_scope_filter(m_scope)) continue;
                    if (kf.length < p_.min_length_short) continue;
                    if (std::max(kb.length, kf.length) < p_.min_length) continue;
                    const double xshift = far_xshift(live_grouping, *kf.cluster, kb.t0 - kf.t0);
                    if (!is_cathode_crossing_pair(*kb.cluster, *kf.cluster,
                                                  kb.length, kf.length, xshift, p_)) continue;
                    best_beam = &kb;
                    best_far = &kf;
                    break;
                }
            }
            if (!best_beam) break;

            // Bundle-structure inputs for the direction rule (a/b/c/d).
            double b_len = 0, d_len = 0;
            for (const auto& m : members) {
                if (m.cluster == best_beam->cluster || m.cluster == best_far->cluster) continue;
                if (m.gid == best_beam->gid) b_len = std::max(b_len, m.length);
                if (m.gid == best_far->gid) d_len = std::max(d_len, m.length);
            }
            const double a_len = best_beam->length;
            const double c_len = best_far->length;
            const bool beam_dom = a_len >= b_len;
            const bool far_dom = c_len >= d_len;
            // Which flash does the joined crosser adopt?  The bundle dominated
            // by its crosser half keeps it; ties go to the longer half.
            bool dest_is_beam;
            const char* rule;
            if (beam_dom && !far_dom)      { dest_is_beam = true;  rule = "beam-dominant"; }
            else if (far_dom && !beam_dom) { dest_is_beam = false; rule = "far-dominant"; }
            else                           { dest_is_beam = (a_len >= c_len); rule = "longer-half"; }

            // Snapshot destination identity and the SOURCE bundle id by value
            // (merge_clusters destroys both member clusters).
            const Member& dest = dest_is_beam ? *best_beam : *best_far;
            const Member& src = dest_is_beam ? *best_far : *best_beam;
            const double dest_t0 = dest.t0;
            const int dest_flash = dest.cluster->get_scalar<int>("flash", -1);
            const int dest_gid = dest.gid;
            const int dest_lm = dest.cluster->get_scalar<int>("lm_flag", -1);
            const int src_gid = src.gid;

            log->info("rescue round {}: c{} (gid {}, t0 {:.3f} us, {:.1f} cm) + "
                      "c{} (gid {}, t0 {:.3f} us, {:.1f} cm) -> gid {} ({}; "
                      "a={:.1f} b={:.1f} c={:.1f} d={:.1f} cm)",
                      nrounds, best_beam->ident, best_beam->gid,
                      best_beam->t0/units::us, best_beam->length/units::cm,
                      best_far->ident, best_far->gid,
                      best_far->t0/units::us, best_far->length/units::cm,
                      dest_gid, rule,
                      a_len/units::cm, b_len/units::cm, c_len/units::cm, d_len/units::cm);

            // Merge exactly as cathode_connect does (2-vertex graph, default
            // merge_clusters args carry the assoc_cluster_id/main provenance).
            typedef cluster_connectivity_graph_t Graph;
            Graph g;
            std::map<int, int> desc;
            int idx_beam = -1, idx_far = -1;
            for (size_t ilive = 0; ilive < live_clusters.size(); ++ilive) {
                if (live_clusters[ilive] == best_beam->cluster) idx_beam = (int)ilive;
                if (live_clusters[ilive] == best_far->cluster) idx_far = (int)ilive;
            }
            for (size_t ilive = 0; ilive < live_clusters.size(); ++ilive) {
                desc[(int)ilive] = boost::add_vertex((int)ilive, g);
            }
            boost::add_edge(desc[idx_beam], desc[idx_far], g);
            auto fresh = merge_clusters(g, live_grouping);
            if (fresh.size() != 1) {
                log->warn("rescue round {}: expected 1 merged cluster, got {} -- stopping",
                          nrounds, fresh.size());
                break;
            }
            Cluster* merged = fresh.at(0);

            // Re-stamp the flash identity AFTER the merge (merge_clusters derives
            // it from the longest member, which may carry the source flash), then
            // re-materialize the corrected coordinates under the destination T0
            // (T0 first, then correct -- the improvecluster_2.cxx idiom; no
            // geometry access on `merged` before add_corrected_points).
            merged->set_cluster_t0(dest_t0);
            merged->set_scalar<int>("flash", dest_flash);
            merged->set_scalar<int>("matched_flash_gid", dest_gid);
            merged->set_scalar<int>("lm_flag", dest_lm);
            // The joined crosser is a bundle main by construction (multi-main
            // flash groups are legal downstream).
            merged->set_flag(Flags::main_cluster, 1);
            merged->set_flag(Flags::associated_cluster, 0);
            const auto& default_scope = merged->get_default_scope();
            const auto& raw_scope = merged->get_raw_scope();
            if (default_scope.hash() != raw_scope.hash()) {
                auto correction_name = merged->get_scope_transform(default_scope);
                merged->add_corrected_points(m_pcts, correction_name);
            }

            // Source-bundle main repair: the source bundle lost its crosser
            // half; if members remain but none is a main, promote the longest
            // so the bundle stays visible to unmerge/TGM/STM/Neutrino.
            Cluster* promote = nullptr;
            bool src_has_main = false;
            for (auto* c : live_grouping.children()) {
                if (c == merged) continue;
                if (c->get_scalar<int>("matched_flash_gid", -1) != src_gid) continue;
                if (c->get_flag(Flags::main_cluster)) { src_has_main = true; break; }
                if (!promote || c->get_length() > promote->get_length()) promote = c;
            }
            if (!src_has_main && promote) {
                promote->set_flag(Flags::main_cluster, 1);
                promote->set_flag(Flags::associated_cluster, 0);
                log->debug("rescue round {}: promoted c{} ({:.1f} cm) to main of source gid {}",
                           nrounds, (int)promote->ident(),
                           promote->get_length()/units::cm, src_gid);
            }

            ++nrounds;
        }

        // ---- Pass 2 (rescue_unmatched; sbnd_xin/docs/pr/17): adopt NON-MATCHED
        // clusters.  Same absorbing-window defect, different failure shape: when
        // the whole far half is ONE cluster (e.g. after the pr/15 vertex veto
        // rejoined it), no wrong flash may claim it at all and Q/L matching
        // leaves it with no flash.  Such a cluster is invisible downstream (no
        // bundle, no beam-window evaluation) AND its corrected coordinates were
        // materialized with the sentinel T0 (56463: x_t0cor off by 1.5e6 m), so
        // it must be evaluated in its RAW scope under the destination-T0
        // hypothesis (raw x is the t0=0 frame; the shift for beam t0 is sub-mm).
        // Direction is FORCED into the beam bundle -- an unmatched cluster has
        // no flash to offer -- so this pass can only ADD charge to the beam
        // bundle; a wrongly adopted cosmic crosser becomes a joined cathode
        // crosser the TGM/STM taggers evaluate downstream (the self-correcting
        // pattern of doc pr/14 sec 7.4, evts 288952/352365).  Runs AFTER the
        // bundle pass, on the post-move state.
        int nadopt = 0;
        if (p_.rescue_unmatched) {
            while (true) {
                auto live_clusters = live_grouping.children();
                // Matched clusters in the analysis scope; unmatched candidates
                // in their raw scope (see above).  Restored after the pass.
                for (auto* live : live_clusters) {
                    const bool matched = live->get_scalar<int>("flash", -1) >= 0 &&
                                         live->get_scalar<int>("matched_flash_gid", -1) >= 0;
                    const auto& want = matched ? m_scope : live->get_raw_scope();
                    if (live->get_default_scope().hash() != want.hash()) {
                        live->set_default_scope(want);
                    }
                }

                std::vector<Member> members;    // flash-matched, ident order
                std::vector<Member> orphans;    // unmatched candidates, ident order
                for (auto* c : live_clusters) {
                    const int gid = c->get_scalar<int>("matched_flash_gid", -1);
                    if (c->get_scalar<int>("flash", -1) >= 0 && gid >= 0) {
                        members.push_back({c, (int)c->ident(), gid,
                                           c->get_cluster_t0(), c->get_length()});
                        continue;
                    }
                    // Size floors keep specks and debris out of the beam bundle.
                    if (c->npoints() < p_.unmatched_min_npts) continue;
                    const double len = c->get_length();
                    if (len < p_.unmatched_min_length) continue;
                    // No scope-filter gate here: the in/out-of-volume partition
                    // needs a T0 the unmatched cluster does not have.  The
                    // tip-at-the-cathode geometry below is the containment gate.
                    orphans.push_back({c, (int)c->ident(), -1, 0.0, len});
                }
                std::sort(members.begin(), members.end(),
                          [](const Member& a, const Member& b) { return a.ident < b.ident; });
                std::sort(orphans.begin(), orphans.end(),
                          [](const Member& a, const Member& b) { return a.ident < b.ident; });

                const Member* best_beam = nullptr;
                const Member* best_orphan = nullptr;
                for (const auto& kb : members) {
                    if (best_beam) break;
                    if (!in_beam(kb.t0)) continue;
                    if (!kb.cluster->get_scope_filter(m_scope)) continue;
                    if (kb.length < p_.min_length_short) continue;
                    for (const auto& ko : orphans) {
                        if (std::max(kb.length, ko.length) < p_.min_length) continue;
                        // Raw x is the t0=0 frame: dt0_dest = t0_beam - 0.
                        const double xshift = far_xshift(live_grouping, *ko.cluster, kb.t0);
                        if (!is_cathode_crossing_pair(*kb.cluster, *ko.cluster,
                                                      kb.length, ko.length, xshift, p_)) continue;
                        best_beam = &kb;
                        best_orphan = &ko;
                        break;
                    }
                }
                if (!best_beam) break;

                // Snapshot the (beam) destination identity by value.
                const double dest_t0 = best_beam->t0;
                const int dest_flash = best_beam->cluster->get_scalar<int>("flash", -1);
                const int dest_gid = best_beam->gid;
                const int dest_lm = best_beam->cluster->get_scalar<int>("lm_flag", -1);

                log->info("unmatched rescue round {}: c{} (gid {}, t0 {:.3f} us, {:.1f} cm) + "
                          "c{} (unmatched, {:.1f} cm, {} pts) -> gid {}",
                          nadopt, best_beam->ident, best_beam->gid,
                          best_beam->t0/units::us, best_beam->length/units::cm,
                          best_orphan->ident, best_orphan->length/units::cm,
                          best_orphan->cluster->npoints(), dest_gid);

                typedef cluster_connectivity_graph_t Graph;
                Graph g;
                std::map<int, int> desc;
                int idx_beam = -1, idx_orphan = -1;
                for (size_t ilive = 0; ilive < live_clusters.size(); ++ilive) {
                    if (live_clusters[ilive] == best_beam->cluster) idx_beam = (int)ilive;
                    if (live_clusters[ilive] == best_orphan->cluster) idx_orphan = (int)ilive;
                }
                for (size_t ilive = 0; ilive < live_clusters.size(); ++ilive) {
                    desc[(int)ilive] = boost::add_vertex((int)ilive, g);
                }
                boost::add_edge(desc[idx_beam], desc[idx_orphan], g);
                auto fresh = merge_clusters(g, live_grouping);
                if (fresh.size() != 1) {
                    log->warn("unmatched rescue round {}: expected 1 merged cluster, got {} -- stopping",
                              nadopt, fresh.size());
                    break;
                }
                Cluster* merged = fresh.at(0);

                // Beam identity, then corrected coordinates under the beam T0
                // (this also replaces the orphan points' sentinel-T0 x_t0cor).
                // The merge inputs had MIXED scopes (orphan raw), so pin the
                // merged cluster's default scope explicitly before correcting.
                merged->set_cluster_t0(dest_t0);
                merged->set_scalar<int>("flash", dest_flash);
                merged->set_scalar<int>("matched_flash_gid", dest_gid);
                merged->set_scalar<int>("lm_flag", dest_lm);
                merged->set_flag(Flags::main_cluster, 1);
                merged->set_flag(Flags::associated_cluster, 0);
                merged->set_default_scope(m_scope);
                if (m_scope.hash() != merged->get_raw_scope().hash()) {
                    auto correction_name = merged->get_scope_transform(m_scope);
                    merged->add_corrected_points(m_pcts, correction_name);
                }
                // No source-bundle repair: the orphan belonged to no bundle.

                ++nadopt;
            }

            // Leave every cluster in the analysis scope, as the bundle pass does
            // (unaccepted orphans were parked in their raw scope above).
            for (auto* live : live_grouping.children()) {
                if (live->get_default_scope().hash() != m_scope.hash()) {
                    live->set_default_scope(m_scope);
                }
            }
        }

        if (nrounds) {
            log->info("cathode bundle rescue: {} move(s) applied", nrounds);
        }
        if (nadopt) {
            log->info("unmatched rescue: {} adoption(s) applied", nadopt);
        }
    }
};


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
