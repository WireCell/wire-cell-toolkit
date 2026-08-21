#include "WireCellUtil/Array.h"  // Eigen, for the local-PCA direction fallback; must precede
                                 // other headers that transitively pull in SIMD intrinsics
                                 // (matches the include order in Facade_Cluster.cxx)

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Util.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/Units.h"
#include "WireCellAux/Logger.h"

#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <unordered_set>
#include <vector>

class ClusteringProtectBundle;
WIRECELL_FACTORY(ClusteringProtectBundle, ClusteringProtectBundle,
                 WireCell::INamed, WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

/**
 * PR-stage overclustering protection: the toolkit counterpart of the
 * prototype's second graph-examination round,
 *
 *     prototype (ProtectOverClustering.cxx lines 6-160, called from
 *     wire-cell-prod-nue.cxx line 1322 after Q/L matching):
 *       for every beam-window matched bundle:
 *         for every member cluster:
 *           Examine_graph(ct_point_cloud)      // rebuild graph, ct-test the
 *                                              // inter-fragment bridges
 *           -> split into surviving connected components;
 *         the main cluster's largest component keeps the main cluster id
 *         (lines 57-121), every other component becomes a new cluster in the
 *         same parent bundle (lines 104-137), which NeutrinoID then fits
 *         separately (wire-cell-prod-nue.cxx lines 1345, 1360).
 *
 * Without this stage the PR graph's uncapped MST bridges charge fragments
 * that Clustering_neutrino merged into the nu cluster (photons near the
 * vertex), and the track fitter interpolates fit points across genuinely
 * charge-free space (sbnd_xin doc pr/22 section 8, SBND evt 386948).
 *
 * Placement (doc pr/23 ordering decision): after the cosmic taggers
 * (tagger_check_tgm/stm/fc) and before tagger_check_neutrino, with a second
 * steiner pass in between.  This is the prototype-faithful order: uboone's
 * cosmic verdicts were computed on the UNSPLIT clusters (the TGM flag prod-stm
 * reads at wire-cell-prod-stm.cxx lines 806-810 comes from the Q/L-stage
 * event_type) and Protect_Over_Clustering ran only in the nue executable
 * immediately before NeutrinoID (wire-cell-prod-nue.cxx line 1322), with the
 * point clouds recomputed right after (lines 1325-1330).  Running the split
 * before the taggers instead defeats TGM (one boundary-touching end split off
 * => both-ends test fails) -- measured at 5 cosmic->nu promotions + 3
 * demotions per 572 data events (doc pr/23 section 6).  The visitor therefore
 * (a) skips bundles whose in-window main already carries a cosmic conviction
 * (skip_convicted, the TaggerCheckNeutrino nu_skip_cosmic condition), and
 * (b) invalidates the retained cluster's pre-split steiner products so the
 * follow-up steiner pass (replace=false) rebuilds exactly the touched
 * clusters.
 *
 * The split/flag/ident idiom follows ClusteringUnmergeBundle mode
 * "component": connected_blobs() on the configured graph flavor, longest
 * component by get_length() keeps the retained cluster (the prototype keeps
 * the component with the most mcells, ProtectOverClustering.cxx lines 60-70;
 * length is the established toolkit idiom -- divergence recorded in doc
 * pr/23), fragments get flag_associated_cluster and a collision-free ident.
 *
 * SBND divergence, knob-gated OFF by default: the relaxed graph does not join
 * the two halves of a cathode crosser (the ~0.9 cm inactive band at x=0 shows
 * up as a ~4-5 cm charge gap with a ~1 cm transverse offset, sbnd_xin doc
 * pr/20), a geometry uboone did not have -- splitting there would undo the
 * cathode_connect/B0 work.  When cathode_rejoin_xcut > 0, component pairs
 * whose closest points both lie within cathode_rejoin_xcut of cathode_x,
 * within cathode_rejoin_dis in 3D and cathode_rejoin_dyz transversely, are
 * re-united before the split.
 *
 * cathode_rejoin_dyz is frame-aligned (a bound on transverse offset in the
 * detector's y-z plane) and is right for a crosser nearly parallel to the
 * drift direction, but for an OBLIQUE crosser the transverse offset across
 * the ~2 cm cathode x-gap is real track travel (dyz ~ |dx|*tan(theta)), so
 * the frame-aligned bound rejects it by construction (sbnd_xin doc pr/25,
 * SBND evt 489327: dyz 6.50 cm on two halves collinear to 1.5 deg at a 75
 * deg angle to drift).  When cathode_rejoin_perp > 0, a dyz failure falls
 * back to a direction-agreement test: local-PCA collinearity
 * (cathode_rejoin_angle) plus a perpendicular-offset bound
 * (cathode_rejoin_perp) between each tip and the other component's local
 * direction line.  dis and xcut are unchanged by this fallback -- a 572-event
 * population census (pr25_rejoin_census.py) found exactly 6 pairs failing
 * ONLY dyz, separating cleanly at 1.5-7.1 deg / perp 1.50-2.68 cm for 5
 * genuine crossers vs 89.2 deg / perp 6.11 cm for the one same-side junk
 * stub (which also fails the point-count floor independently).
 */
class ClusteringProtectBundle : public IConfigurable,
                                public Clus::IEnsembleVisitor,
                                public Aux::Logger,
                                private NeedDV,
                                private NeedPCTS {
public:
    ClusteringProtectBundle()
        : Aux::Logger("ClusteringProtectBundle", "clus")
    {}
    virtual ~ClusteringProtectBundle() {}

    virtual void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedPCTS::configure(config);
        m_grouping_name = get<std::string>(config, "grouping", m_grouping_name);
        m_graph_name = get<std::string>(config, "graph_name", m_graph_name);
        m_pcarray_name = get<std::string>(config, "pcarray_name", m_pcarray_name);
        m_require_in_scope = get<bool>(config, "require_in_scope", m_require_in_scope);
        // Restrict to the beam-coincident bundle(s), the prototype's scope:
        // to_be_checked there is built from beam-window flashes only
        // (wire-cell-prod-nue.cxx lines 1313-1320).  Same gate and same key
        // names as CreateSteinerGraph.  Default false = every scope-passing
        // cluster.
        m_beam_window_only = get<bool>(config, "beam_window_only", m_beam_window_only);
        m_beam_window_low = get<double>(config, "beam_window_low", m_beam_window_low);
        m_beam_window_high = get<double>(config, "beam_window_high", m_beam_window_high);
        // Only meaningful with beam_window_only: a convicted in-window main
        // (TaggerCheckTGM/STM flags or the Q/L LM tagger's lm_flag > 0 -- the
        // exact TaggerCheckNeutrino nu_skip_cosmic condition) does not open
        // its bundle for splitting.  Default true: at the post-tagger pipeline
        // position, splitting a convicted cosmic bundle cannot change any
        // verdict but churns the archives of cosmic-tagged events; skipping
        // keeps them byte-identical to the stage-off chain.  An unconvicted
        // in-window main sharing the gid (the doc pr/16 spared-mate case)
        // still opens the bundle.
        m_skip_convicted = get<bool>(config, "skip_convicted", m_skip_convicted);
        m_open_convicted_bundles = get<bool>(config, "open_convicted_bundles", m_open_convicted_bundles);
        // Cathode re-join (SBND).  cathode_rejoin_xcut <= 0 (default) disables
        // the pass entirely => prototype-faithful behavior.  Internal units.
        m_cathode_x = get<double>(config, "cathode_x", m_cathode_x);
        m_cathode_rejoin_xcut = get<double>(config, "cathode_rejoin_xcut", m_cathode_rejoin_xcut);
        m_cathode_rejoin_dyz = get<double>(config, "cathode_rejoin_dyz", m_cathode_rejoin_dyz);
        m_cathode_rejoin_dis = get<double>(config, "cathode_rejoin_dis", m_cathode_rejoin_dis);
        // Direction-agreement fallback for a dyz-only failure (SBND doc
        // pr/25).  cathode_rejoin_perp <= 0 (default) disables the fallback
        // entirely => byte-identical to the dyz-only behavior above.
        // Internal units except cathode_rejoin_angle (degrees).
        m_cathode_rejoin_perp = get<double>(config, "cathode_rejoin_perp", m_cathode_rejoin_perp);
        m_cathode_rejoin_angle = get<double>(config, "cathode_rejoin_angle", m_cathode_rejoin_angle);
        m_cathode_rejoin_dir_radius = get<double>(config, "cathode_rejoin_dir_radius", m_cathode_rejoin_dir_radius);
        m_cathode_rejoin_dir_npts = get<int>(config, "cathode_rejoin_dir_npts", m_cathode_rejoin_dir_npts);
    }

    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["detector_volumes"] = "DetectorVolumes";
        cfg["pc_transforms"] = "PCTransformSet";
        cfg["grouping"] = m_grouping_name;
        cfg["graph_name"] = m_graph_name;
        cfg["pcarray_name"] = m_pcarray_name;
        cfg["require_in_scope"] = m_require_in_scope;
        cfg["beam_window_only"] = m_beam_window_only;
        cfg["beam_window_low"] = m_beam_window_low;
        cfg["beam_window_high"] = m_beam_window_high;
        cfg["skip_convicted"] = m_skip_convicted;
        cfg["open_convicted_bundles"] = m_open_convicted_bundles;   // doc pr/94 round 3; false = a convicted main does not open its bundle (legacy)
        cfg["cathode_x"] = m_cathode_x;
        cfg["cathode_rejoin_xcut"] = m_cathode_rejoin_xcut;
        cfg["cathode_rejoin_dyz"] = m_cathode_rejoin_dyz;
        cfg["cathode_rejoin_dis"] = m_cathode_rejoin_dis;
        cfg["cathode_rejoin_perp"] = m_cathode_rejoin_perp;
        cfg["cathode_rejoin_angle"] = m_cathode_rejoin_angle;
        cfg["cathode_rejoin_dir_radius"] = m_cathode_rejoin_dir_radius;
        cfg["cathode_rejoin_dir_npts"] = m_cathode_rejoin_dir_npts;
        return cfg;
    }

    virtual void visit(Ensemble& ensemble) const {
        auto groupings = ensemble.with_name(m_grouping_name);
        if (groupings.empty()) {
            log->warn("no '{}' grouping found, skipping", m_grouping_name);
            return;
        }
        auto& grouping = *groupings.at(0);

        // Snapshot: separate() adds children to the grouping.
        const std::vector<Cluster*> clusters = grouping.children();

        std::set<int> taken;
        for (const auto* c : clusters) taken.insert(c->ident());

        // Beam gate (CreateSteinerGraph idiom): gids of in-window mains; a
        // bundle member is the main itself or a companion sharing its gid.
        std::unordered_set<int> beam_gids;
        const bool gate = m_beam_window_only && m_beam_window_low < m_beam_window_high;
        size_t n_convicted = 0;
        if (gate) {
            for (const auto* cluster : clusters) {
                if (!cluster->get_flag(Flags::main_cluster)) continue;
                const double t0 = cluster->get_cluster_t0();
                if (t0 < m_beam_window_low || t0 >= m_beam_window_high) continue;
                const int gid = cluster->get_scalar<int>("matched_flash_gid", -1);
                if (gid < 0) continue;
                // skip_convicted: honor the upstream cosmic verdicts
                // (TaggerCheckNeutrino nu_skip_cosmic condition, lines
                // 339-350: TGM/STM flags, lm_flag 1=low-energy
                // 2=light-mismatch, absent scalar never skips).
                if (m_skip_convicted
                    && (cluster->get_flag(Flags::TGM) || cluster->get_flag(Flags::STM)
                        || cluster->get_scalar<int>("lm_flag", -1) > 0)) {
                    // doc pr/94 round 3 -- open_convicted_bundles.  A convicted
                    // main suppresses the second graph examination for its WHOLE
                    // bundle, so a co-bundled secondary activity -- the one the
                    // demoted-main fallback / per-bundle mode will go on to call
                    // the neutrino -- never gets the split that every cluster in
                    // an unconvicted bundle gets (SBND 18255/395148: the 198.9 cm
                    // secondary keeps a graph bridge that fits 17 cm of
                    // trajectory through empty space).  When on, the convicted
                    // main still OPENS its bundle; it is still never SPLIT
                    // itself -- the per-member guard below is untouched, so the
                    // cosmic-tagged tree is left exactly as it is.
                    if (m_open_convicted_bundles) {
                        log->debug("OC94OPEN main ident={} nblobs={} gid={} t0={:.2f}us "
                                   "convicted TGM={} STM={} lm={} -- bundle OPENED for its "
                                   "unconvicted members (open_convicted_bundles)",
                                   cluster->ident(), cluster->nchildren(), gid, t0/units::us,
                                   cluster->get_flag(Flags::TGM), cluster->get_flag(Flags::STM),
                                   cluster->get_scalar<int>("lm_flag", -1));
                        beam_gids.insert(gid);
                        continue;   // counted as opened, not as skipped
                    }
                    ++n_convicted;
                    // doc pr/53 round 6: name the gate per cluster (log-only)
                    log->debug("OC53SKIP main ident={} nblobs={} gid={} t0={:.2f}us "
                               "convicted TGM={} STM={} lm={} -- bundle not opened",
                               cluster->ident(), cluster->nchildren(), gid, t0/units::us,
                               cluster->get_flag(Flags::TGM), cluster->get_flag(Flags::STM),
                               cluster->get_scalar<int>("lm_flag", -1));
                    continue;   // this main does not open its bundle
                }
                beam_gids.insert(gid);
            }
        }

        size_t n_split = 0, n_parts = 0, n_rejoin = 0;
        for (Cluster* cluster : clusters) {
            if (m_require_in_scope
                && !cluster->get_scope_filter(cluster->get_default_scope())) continue;
            if (gate) {
                const int gid = cluster->get_scalar<int>("matched_flash_gid", -1);
                if (gid < 0 || beam_gids.count(gid) == 0) continue;
            }
            // A convicted cluster is never split, even as a member of a
            // bundle an unconvicted mate opened (SBND evt 52195: gid 6 holds
            // the TGM-tagged 513 cm crosser AND an untagged shard main --
            // the shard opens the bundle, but splitting the convicted
            // crosser would only churn the cosmic-tagged tree).
            if (m_skip_convicted
                && (cluster->get_flag(Flags::TGM) || cluster->get_flag(Flags::STM)
                    || cluster->get_scalar<int>("lm_flag", -1) > 0)) {
                ++n_convicted;
                // doc pr/53 round 6: name the gate per cluster (log-only)
                log->debug("OC53SKIP member ident={} nblobs={} main={} "
                           "convicted TGM={} STM={} lm={} -- never split",
                           cluster->ident(), cluster->nchildren(),
                           cluster->get_flag(Flags::main_cluster),
                           cluster->get_flag(Flags::TGM), cluster->get_flag(Flags::STM),
                           cluster->get_scalar<int>("lm_flag", -1));
                continue;
            }
            const size_t before = grouping.nchildren();
            if (process_cluster(grouping, cluster, taken, n_rejoin)) {
                ++n_split;
                n_parts += grouping.nchildren() - before;
            }
        }
        log->debug("split {} bundle cluster(s) into {} extra cluster(s) "
                   "({} cathode re-join(s), {} convicted main(s) skipped, graph '{}')",
                   n_split, n_parts, n_rejoin, n_convicted, m_graph_name);
    }

private:
    std::string m_grouping_name{"live"};
    std::string m_graph_name{"relaxed"};
    std::string m_pcarray_name{"perblob"};
    bool m_require_in_scope{true};
    bool m_beam_window_only{false};
    bool m_skip_convicted{true};
    // doc pr/94 round 3.  When true, a cosmic-convicted main still contributes
    // its gid to beam_gids -- i.e. it opens its bundle so the bundle's OTHER,
    // unconvicted members are graph-examined and split -- while the per-member
    // skip_convicted guard keeps the convicted cluster itself from ever being
    // split.  Default false = the pre-pr/94 behaviour, byte-identical.
    bool m_open_convicted_bundles{false};
    double m_beam_window_low{0};
    double m_beam_window_high{0};
    double m_cathode_x{0};
    double m_cathode_rejoin_xcut{0};       // <= 0: pass disabled (prototype behavior)
    double m_cathode_rejoin_dyz{4 * units::cm};
    double m_cathode_rejoin_dis{8 * units::cm};
    double m_cathode_rejoin_perp{0};       // <= 0: direction fallback disabled (default)
    double m_cathode_rejoin_angle{20.0};   // degrees
    double m_cathode_rejoin_dir_radius{15 * units::cm};
    int m_cathode_rejoin_dir_npts{20};

    // Fold an angle (radians, 0..pi) about pi -> collinearity in degrees
    // (0 = parallel or anti-parallel).  Same convention as
    // clustering_cathode_connect.cxx / clustering_cathode_bundle_rescue.cxx
    // (file-local by the repo's fork-by-duplication convention, not shared).
    static double collinear_deg(double angle_rad) {
        const double a = angle_rad / 3.1415926 * 180.0;
        return std::min(a, 180.0 - a);
    }

    /// Local principal-axis direction of `cloud` within
    /// m_cathode_rejoin_dir_radius of `tip`, oriented away from the local
    /// centroid (i.e. pointing outward from the tip into the component).
    /// Returns false if fewer than m_cathode_rejoin_dir_npts points fall in
    /// the radius (a floor against a direction estimated from a handful of
    /// points -- the doc pr/25 census's rejected junk stub is 12 points
    /// total and fails this independently of the angle test).
    bool local_direction(const Simple3DPointCloud& cloud, const geo_point_t& tip,
                         geo_point_t& dir) const {
        const auto near = cloud.get_closest_wcpoints_radius(tip, m_cathode_rejoin_dir_radius);
        if ((int) near.size() < m_cathode_rejoin_dir_npts) return false;

        geo_point_t center;
        for (const auto& [ind, p] : near) center += p;
        center /= (double) near.size();

        // Same 3x3 covariance + SelfAdjointEigenSolver pattern as
        // Facade_Cluster::get_pca(), restricted to the local neighborhood.
        Eigen::MatrixXd cov(3, 3);
        for (int i = 0; i != 3; ++i) {
            for (int j = i; j != 3; ++j) {
                double s = 0;
                for (const auto& [ind, p] : near) s += (p[i] - center[i]) * (p[j] - center[j]);
                cov(i, j) = s;
            }
        }
        cov(1, 0) = cov(0, 1);
        cov(2, 0) = cov(0, 2);
        cov(2, 1) = cov(1, 2);

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov);
        const auto vecs = es.eigenvectors();
        // Eigen returns ascending eigenvalues; the principal axis is the last column.
        geo_point_t v(vecs(0, 2), vecs(1, 2), vecs(2, 2));
        const double norm = v.magnitude();
        if (norm <= 0) return false;
        v = v / norm;
        // Orient outward from the tip (away from the local centroid), so the
        // two tips' directions can be compared for anti-collinearity without
        // a sign ambiguity flipping the perpendicular-offset test below.
        if ((tip - center).dot(v) < 0) v = v * -1;
        dir = v;
        return true;
    }

    /// Direction-agreement fallback for a dyz-only cathode_rejoin failure:
    /// local-PCA collinearity at the two tips plus a perpendicular-offset
    /// bound, in place of the frame-aligned dyz cut.  cloud1/cloud2 MUST be
    /// the same per-component point clouds cathode_rejoin() already built --
    /// do not query Cluster::vhough_transform() here, since at re-join time
    /// both components are still blobs of the same not-yet-split cluster and
    /// a cluster-level direction query would mix in the OTHER component's
    /// points.
    bool direction_agrees(const Simple3DPointCloud& cloud1, const Simple3DPointCloud& cloud2,
                          const geo_point_t& p1, const geo_point_t& p2) const {
        geo_point_t v1, v2;
        if (!local_direction(cloud1, p1, v1)) return false;
        if (!local_direction(cloud2, p2, v2)) return false;
        if (collinear_deg(v1.angle(v2)) >= m_cathode_rejoin_angle) return false;

        const geo_point_t d12 = p2 - p1;
        const double perp1 = (d12 - v1 * d12.dot(v1)).magnitude();
        const geo_point_t d21 = p1 - p2;
        const double perp2 = (d21 - v2 * d21.dot(v2)).magnitude();
        return std::max(perp1, perp2) < m_cathode_rejoin_perp;
    }

    /// Union components that meet across the cathode band.  `b2g` is edited
    /// in place; returns the number of pair unions applied.
    size_t cathode_rejoin(const Cluster* cluster, std::vector<int>& b2g) const {
        std::set<int> ids;
        for (int id : b2g) if (id >= 0) ids.insert(id);
        if (ids.size() < 2) return 0;

        // One point cloud per component, points mapped through their blob
        // (same construction as connect_graph_relaxed.cxx lines 71-82; the
        // point coordinates are in the default scope, i.e. the T0-corrected
        // frame that cathode_x is defined in).
        std::map<int, Simple3DPointCloud> clouds;
        const auto& points = cluster->points();
        const size_t npts = points[0].size();
        for (size_t i = 0; i < npts; ++i) {
            const int comp = b2g[cluster->kd3d().major_index(i)];
            if (comp < 0) continue;
            clouds[comp].add({points[0][i], points[1][i], points[2][i]});
        }

        // Union-find over component ids.
        std::map<int, int> parent;
        for (int id : ids) parent[id] = id;
        std::function<int(int)> find = [&](int a) {
            while (parent[a] != a) { parent[a] = parent[parent[a]]; a = parent[a]; }
            return a;
        };

        size_t n_union = 0;
        for (auto j = clouds.begin(); j != clouds.end(); ++j) {
            for (auto k = std::next(j); k != clouds.end(); ++k) {
                const auto [i1, i2, dis] = j->second.get_closest_points(k->second);
                if (dis >= m_cathode_rejoin_dis) continue;
                const auto p1 = j->second.point(i1);
                const auto p2 = k->second.point(i2);
                if (std::abs(p1.x() - m_cathode_x) >= m_cathode_rejoin_xcut) continue;
                if (std::abs(p2.x() - m_cathode_x) >= m_cathode_rejoin_xcut) continue;
                const double dyz = std::hypot(p1.y() - p2.y(), p1.z() - p2.z());
                bool via_direction = false;
                if (dyz >= m_cathode_rejoin_dyz) {
                    // Oblique crosser: dyz is real track travel across the
                    // cathode x-gap, not misalignment (doc pr/25, SBND evt
                    // 489327).  Disabled by default (m_cathode_rejoin_perp
                    // <= 0) => byte-identical to the check above.
                    if (m_cathode_rejoin_perp <= 0) continue;
                    if (!direction_agrees(j->second, k->second, p1, p2)) continue;
                    via_direction = true;
                }
                const int ra = find(j->first), rb = find(k->first);
                if (ra == rb) continue;
                parent[std::max(ra, rb)] = std::min(ra, rb);   // keep lowest id
                ++n_union;
                SPDLOG_LOGGER_DEBUG(log,
                    "cluster {}: cathode re-join comp {}+{} (gap {:.2f} cm, dyz {:.2f} cm, "
                    "x {:.2f}/{:.2f} cm{})",
                    cluster->ident(), j->first, k->first, dis / units::cm,
                    dyz / units::cm, p1.x() / units::cm, p2.x() / units::cm,
                    via_direction ? ", via direction fallback" : "");
            }
        }
        if (n_union) {
            for (int& id : b2g) if (id >= 0) id = find(id);
        }
        return n_union;
    }

    /// Split one bundle cluster at its surviving component boundaries.
    /// Returns true if a split happened.
    bool process_cluster(Grouping& grouping, Cluster*& cluster,
                         std::set<int>& taken, size_t& n_rejoin) const {
        const size_t nb = cluster->nchildren();
        if (nb < 2) return false;

        // prototype: Examine_graph(ct_point_cloud), ProtectOverClustering.cxx:56
        auto b2g = cluster->connected_blobs(m_dv, m_pcts, m_graph_name);
        if (b2g.size() != nb) return false;

        if (m_cathode_rejoin_xcut > 0) n_rejoin += cathode_rejoin(cluster, b2g);

        std::set<int> ids;
        for (int id : b2g) if (id >= 0) ids.insert(id);
        if (ids.size() < 2) return false;

        // prototype keeps the component with the most mcells
        // (ProtectOverClustering.cxx:60-70); here longest by get_length(),
        // ascending iteration so ties keep the lowest id.
        double max_length = -1;
        int longest = -1;
        for (int id : ids) {
            const double len = get_length(cluster, b2g, id);
            if (len > max_length) { max_length = len; longest = id; }
        }

        std::map<int, int> id2grp;
        int g = 0;
        for (int id : ids) if (id != longest) id2grp[id] = g++;

        std::vector<int> groups(nb, -1);
        for (size_t i = 0; i < nb; ++i) {
            if (b2g[i] < 0 || b2g[i] == longest) continue;   // unassigned stays put
            groups[i] = id2grp.at(b2g[i]);
        }

        // Snapshot the per-blob PC before the blobs move away (same landmine
        // as ClusteringUnmergeBundle: separate(remove=false) leaves the
        // retained cluster holding a full-length array over a shorter blob
        // list).
        PointCloud::Dataset perblob;
        bool have_perblob = false;
        {
            auto& lpcs = cluster->value().local_pcs();
            auto it = lpcs.find(m_pcarray_name);
            if (it != lpcs.end()) {
                if (it->second.size_major() == nb) {
                    perblob = it->second;      // copy
                    have_perblob = true;
                }
                else {
                    log->warn("cluster {}: '{}' has {} rows for {} blobs, dropping it",
                              cluster->ident(), m_pcarray_name,
                              it->second.size_major(), nb);
                }
            }
        }

        const int main_ident = cluster->ident();
        const bool was_main = cluster->get_flag(Flags::main_cluster);
        auto splits = grouping.separate(cluster, groups);   // remove=false: largest stays put

        auto carve = [&](Cluster* part, int want) {
            auto& lpcs = part->value().local_pcs();
            lpcs.erase(m_pcarray_name);
            if (!have_perblob) return;
            std::vector<size_t> rows;
            rows.reserve(nb);
            for (size_t i = 0; i < nb; ++i) if (groups[i] == want) rows.push_back(i);
            if (rows.size() != part->nchildren()) {
                log->warn("cluster {}: carved {} rows for {} blobs, dropping '{}'",
                          part->ident(), rows.size(), part->nchildren(), m_pcarray_name);
                return;
            }
            lpcs[m_pcarray_name] = perblob.subset(rows);
        };

        carve(cluster, -1);
        // Invalidate the retained cluster's pre-split products.  At the
        // production position (after the cosmic taggers, doc pr/23 ordering)
        // the steiner stage and the taggers have already run on this facade:
        // the node-local steiner_pc, the facade-owned steiner_graph, AND every
        // other facade-owned graph (STM's track fit builds 'basic_pid', the
        // taggers build ctpc flavors, connected_blobs above built
        // m_graph_name) were all made on the pre-split blob set, so their
        // vertex indices exceed the post-split point count
        // (blob_with_point(3906) on 2738 points -- evt 386948 abort,
        // TrackFitting::form_point_association via graph_algorithms
        // 'basic_pid').  Named flavors rebuild lazily
        // (Cluster::graph_algorithms), steiner_graph/steiner_pc by the
        // follow-up steiner pass; fragments come out of separate() with fresh
        // facades carrying none of these (the ClusteringUnmergeBundle lines
        // 78-84 observation).
        if (cluster->value().local_pcs().erase("steiner_pc")) {
            log->debug("cluster {}: dropped the pre-split steiner_pc "
                       "(follow-up steiner pass rebuilds it)", main_ident);
        }
        {
            std::vector<std::string> gnames;
            for (const auto& [gname, gr] : cluster->graph_store()) gnames.push_back(gname);
            for (const auto& gname : gnames) {
                cluster->take_graph(gname);                 // erases the store entry
                cluster->remove_graph_algorithms(gname);
            }
            // A GraphAlgorithms cache can exist without a stored graph
            // (reference-graph constructions) -- drop those too.
            for (const auto& gname : cluster->get_cached_graph_algorithms()) {
                cluster->remove_graph_algorithms(gname);
            }
            if (!gnames.empty()) {
                log->debug("cluster {}: dropped {} pre-split graph(s)",
                           main_ident, gnames.size());
            }
        }

        // prototype: the largest component keeps the main cluster id and role
        // (ProtectOverClustering.cxx:101-110); the retained cluster already
        // carries the flags and cluster_scalar (matched_flash_gid, cluster_t0)
        // -- fragments get them via separate()->from(), so a fragment of the
        // bundle main must have flag_main_cluster cleared or every fragment
        // would look like a bundle main to STM/TGM/FC.
        size_t nb_parts = 0;
        for (const auto& [gid, part] : splits) nb_parts += part->nchildren();
        log->debug("cluster {}{}: {} blobs -> retained {} + {} fragment(s) holding {}",
                   main_ident, was_main ? " (main)" : "", nb, cluster->nchildren(),
                   splits.size(), nb_parts);
        if (cluster->nchildren() + nb_parts != nb) {
            log->error("cluster {}: BLOB LOSS in separate(): {} + {} != {}",
                       main_ident, cluster->nchildren(), nb_parts, nb);
        }
        int sub_id = 1;
        for (auto& [gid, part] : splits) {
            part->set_flag(Flags::main_cluster, 0);
            part->set_flag(Flags::associated_cluster);
            // prototype: fresh ids via get_next_cluster_id
            // (ProtectOverClustering.cxx:113-124); here the collision-avoiding
            // ClusteringUnmergeBundle scheme.
            part->set_ident(alloc_ident(taken, main_ident * 100 + sub_id));
            ++sub_id;
            carve(part, gid);
        }
        return true;
    }

    /// `want` if free, else the lowest free integer above every taken ident.
    static int alloc_ident(std::set<int>& taken, int want) {
        if (!taken.count(want)) { taken.insert(want); return want; }
        int cand = taken.empty() ? 0 : *taken.rbegin() + 1;
        while (taken.count(cand)) ++cand;
        taken.insert(cand);
        return cand;
    }
};
