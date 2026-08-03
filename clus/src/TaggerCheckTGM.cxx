/// TaggerCheckTGM: through-going-muon (TGM) tagger.
///
/// Port of the WCP prototype WCPPID::ToyFiducial::check_tgm
/// (prototype pid/src/Cosmic_tagger.h:1331; algorithm walkthrough in
/// clus/docs/tgm/check_tgm.html).  A cluster is a TGM when its two ends both
/// exit the fiducial volume (case A), or when an apparently-inside end is
/// explained by a prolonged-signal artefact or a dead region (case B).
///
/// Differences from the prototype, by design:
///  - The FV inside/outside tests run against a dedicated IFiducial
///    ("fiducial" config; e.g. a BoxFiducial spanning BOTH SBND TPCs so a
///    cathode-crossing track does not look like an exiter at x=0), with an
///    optional tolerance vector applying the boundary margins.  The
///    dead-region / signal-processing checks still use the grouping's
///    FiducialUtils (per-(apa,face) dead-channel logic).
///  - offset_x = 0 everywhere: this runs after switch_scope, so point
///    coordinates are already T0-corrected by the matched flash time.
///  - The prototype's main_flash->get_type()==2 (beam flash) protection is
///    replaced by a beam-window test on cluster_t0.  With the
///    check_neutrino_candidate knob OFF (default, v1 behavior) the protected
///    branches never tag an in-beam-window bundle.  With the knob ON the
///    ported check_neutrino_candidate() (prototype Cosmic_tagger.h:1677)
///    arbitrates as in the prototype; see that method below for the
///    two-TPC generalizations.
///  - Multi-main: every flagged main cluster is checked (uBooNE = one).

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellClus/FiducialUtils.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Units.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numeric>
#include <set>

class TaggerCheckTGM;
WIRECELL_FACTORY(TaggerCheckTGM, TaggerCheckTGM,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

static auto t_log = WireCell::Log::logger("clus.NeutrinoPattern");

class TaggerCheckTGM : public IConfigurable, public Clus::IEnsembleVisitor,
                       private Clus::NeedDV, private Clus::NeedPCTS, private Clus::NeedFiducial {
public:
    TaggerCheckTGM() {}
    virtual ~TaggerCheckTGM() {}

    virtual void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedPCTS::configure(config);
        NeedFiducial::configure(config);
        m_grouping_name = get<std::string>(config, "grouping", m_grouping_name);
        m_beam_window_low = get(config, "beam_window_low", m_beam_window_low);
        m_beam_window_high = get(config, "beam_window_high", m_beam_window_high);
        // beam_window_only (default false = evaluate every main cluster):
        // restrict the tagger to the beam-coincident bundle, i.e. to mains
        // whose matched flash time (cluster_t0) lies in the SAME
        // [beam_window_low, beam_window_high) window the in-beam protection
        // below already uses.  Off, or with an empty window (low >= high), the
        // cluster selection is byte-identical to the pre-knob one.  Verdicts on
        // the surviving mains are untouched: this only removes clusters from
        // the loop, and each main's check_tgm() reads no other cluster.
        m_beam_window_only = get(config, "beam_window_only", m_beam_window_only);
        if (m_beam_window_only && !(m_beam_window_low < m_beam_window_high)) {
            SPDLOG_LOGGER_WARN(t_log, "configure: TaggerCheckTGM: beam_window_only set but the window is empty ([{}, {}) us) -- gate disabled, every main evaluated",
                               m_beam_window_low/units::us, m_beam_window_high/units::us);
        }
        m_length_limit_frac = get(config, "length_limit_frac", m_length_limit_frac);
        m_enable_case_b = get(config, "enable_case_b", m_enable_case_b);
        // require_in_scope (default false = historical behavior): also require
        // the cluster to pass the default-scope filter set by switch_scope, i.e.
        // to have at least one blob whose T0-corrected points land inside the
        // active volume.  switch_scope SEPARATES the failing blobs into their own
        // cluster which stays in the grouping and inherits flag_main_cluster, so
        // without this the tagger also evaluates non-physical out-of-volume
        // shards -- which sit outside the FV by construction and therefore
        // satisfy the CASE-A through-going test almost automatically.  The Bee
        // writer (filter:1) and clustering_examine_bundles already honor this
        // same flag; the taggers were the only consumers ignoring it.
        m_require_in_scope = get(config, "require_in_scope", m_require_in_scope);
        // evaluate_demoted_mains (default false = historical behavior): also
        // evaluate a cluster carrying Flags::demoted_main -- a split part that
        // was ITSELF a matched Q/L bundle main before the flash-time merge
        // demoted it (ClusteringUnmergeBundle restore_demoted_mains, doc pr/20
        // Part I).  Such a cluster is not a fragment; it is exactly the
        // population these cuts were tuned on, and today no cosmic tagger ever
        // looks at it.  The scope filter and beam-window gate below apply to it
        // unchanged.  Default false => no cluster is added => byte-identical.
        m_evaluate_demoted_mains = get<bool>(config, "evaluate_demoted_mains", m_evaluate_demoted_mains);

        // check_neutrino_candidate (default false = v1 behavior): enable the
        // ported prototype neutrino-candidate veto in the beam-protected
        // branches.  OFF keeps the conservative never-tag-in-beam behavior
        // byte-identical.
        m_check_neutrino_candidate = get(config, "check_neutrino_candidate", m_check_neutrino_candidate);
        // require_chord_charge (default false = historical behavior): also
        // require the cluster to carry charge ALONG the chord between the two
        // extreme points before that pair may tag.  See chord_has_charge().
        m_require_chord_charge = get(config, "require_chord_charge", m_require_chord_charge);
        m_chord_support_radius = get(config, "chord_support_radius", m_chord_support_radius);
        m_chord_max_gap = get(config, "chord_max_gap", m_chord_max_gap);
        // chord_charge_mode (default "chord" = the straight-chord sampling
        // test, preserved for reproducibility of earlier runs): "path"
        // replaces it with the piecewise charge-path test, robust to track
        // curvature -- see path_components().  Only consulted when
        // require_chord_charge is on.
        m_chord_charge_mode = get<std::string>(config, "chord_charge_mode", m_chord_charge_mode);
        if (m_chord_charge_mode != "chord" && m_chord_charge_mode != "path") {
            SPDLOG_LOGGER_WARN(t_log, "TaggerCheckTGM: unknown chord_charge_mode '{}', using 'chord'",
                               m_chord_charge_mode);
            m_chord_charge_mode = "chord";
        }
        // component_extremes (default false = historical behavior): find the
        // extreme points PER connected component and union them, instead of
        // taking 8 global extremes over a flash-merged cluster.  See
        // component_extreme_wcps().  Intended to be used together with
        // require_chord_charge.
        m_component_extremes = get(config, "component_extremes", m_component_extremes);
        m_component_min_length = get(config, "component_min_length", m_component_min_length);
        m_component_graph = get<std::string>(config, "component_graph", m_component_graph);
        // component_rescue (default false = historical behavior): a component
        // SHORTER than component_min_length still donates its extreme points
        // when it is path-connected (same path_components() component, i.e.
        // reachable through the cluster's own charge with no jump longer than
        // chord_max_gap) to a component that passed the length cut.  The
        // min-length cut exists so a detached merge-grafted speck cannot fake
        // a TGM end (evt285185 cluster 20, 608 cm away); but a genuine track
        // END that fragments into a sub-10 cm piece behind small gaps is
        // dropped by the same cut and the track loses its pair (SBND
        // evt286681 cluster 7: anode->top corner clipper whose last 2.5 cm
        // before the top wall sits behind 4-6 cm gaps).  Path connectivity
        // separates the two cases with the same 30 cm-step rule the chord
        // guard uses.  Only consulted when component_extremes is on.
        m_component_rescue = get(config, "component_rescue", m_component_rescue);
        // rescue_chord_check (default false = rescue-as-introduced behavior):
        // a CASE-A/CASE-B pair with an end donated by a RESCUED component must
        // ALSO pass the straight-chord support test (chord_has_charge), even
        // in path mode.  Rescue's target -- a genuine track end fragmented
        // behind small gaps -- lies ON the straight line to the track's other
        // end, so it passes.  But path mode alone lets a rescued speck pair
        // across TWO merged tracks: SBND evt288727 cluster 6, two touching
        // cosmics (one enters downstream, one enters the bottom via a
        // drift-parallel fragmented track), pair joined by an L-shaped detour
        // with jumps < 30 cm while 88% of the straight chord is > 6 cm from
        // any charge.  Only consulted when component_rescue is on.
        m_rescue_chord_check = get(config, "rescue_chord_check", m_rescue_chord_check);
        // main_component_pairs (default false = historical behavior): a
        // CASE-A/CASE-B pair may tag only when at least one end lies in the
        // cluster's MAIN charge component -- the path_components() component
        // (30 cm-step linkage, so a cathode crosser stays one component)
        // holding the most points.  On a flash-merged Cluster a merged-in
        // fragment that is itself through-going otherwise tags the whole
        // bundle on its own within-component pair, which the chord guard by
        // design does not reject: SBND evt289343 cluster 9 (in beam window),
        // a 26 cm bottom->downstream corner-clipping cosmic fragment ~450 cm
        // from the 52 cm main track tagged the bundle TGM, and the
        // check_neutrino_candidate veto could not protect because it walks
        // the path between the PAIR's own endpoints.  The TGM verdict should
        // follow the bundle's main cluster; the largest path component is
        // its proxy (per-point provenance of the pre-merge main cluster does
        // not survive the merge).
        m_main_component_pairs = get(config, "main_component_pairs", m_main_component_pairs);
        // main_component_mode (default "path" = doc-36 behavior): how the
        // "main" is identified when main_component_pairs is on.
        //  "path": proxy -- the largest 30 cm path component.  Wrong when a
        //          merged-in cosmic is larger than the bundle's main cluster.
        //  "real": exact -- the per-blob "real_cluster_main" perblob array
        //          written by the QL flash merge (marks the REPRESENTATIVE
        //          member: the flash/flags donor, i.e. the one main the
        //          merged cluster reports), persisted through the pctree
        //          tarball when the QL job saves with save_real_cluster_id.
        //          Falls back to the "path" proxy when the array is absent
        //          (old tarballs, unmerged-only events).
        m_main_component_mode = get<std::string>(config, "main_component_mode", m_main_component_mode);
        if (m_main_component_mode != "path" && m_main_component_mode != "real") {
            SPDLOG_LOGGER_WARN(t_log, "TaggerCheckTGM: unknown main_component_mode '{}', using 'path'",
                               m_main_component_mode);
            m_main_component_mode = "path";
        }
        // exempt_demoted_main_pairs (default false = historical behavior):
        // skip main_pair_rejects entirely for a cluster carrying
        // Flags::demoted_main.  A demoted main's OWN real_cluster_main /
        // path-component provenance is all-zero by construction after
        // ClusteringUnmergeBundle::carve (it was A bundle main, but never
        // the merge's chosen representative, and nothing re-stamps the array
        // post-split) -- so with main_component_pairs on, every demoted-main
        // pair was vetoed unconditionally, independent of geometry (sbnd_xin
        // doc pr/25, SBND evt 320029/18255-1: cluster 30, a 37 cm
        // corner-clipping CASE-A shape with both ends on a boundary, rejected
        // before the boundary geometry ran).  Only meaningful together with
        // evaluate_demoted_mains (which is what lets such a cluster reach
        // this tagger at all).  Default off => byte-identical.
        m_exempt_demoted_main_pairs = get<bool>(config, "exempt_demoted_main_pairs", m_exempt_demoted_main_pairs);
        auto tol = config["fv_tolerance"];
        if (!tol.isNull() && tol.isArray()) {
            m_fv_tolerance.clear();
            for (const auto& t : tol) m_fv_tolerance.push_back(t.asDouble());
        }
        // interior_fv_tolerance (default empty = use fv_tolerance,
        // byte-identical): separate tolerance vector for the CASE-A
        // interior-support tests only (flag_check midpoints and
        // flag_check_again waypoint chords).  The ENDPOINT outside/inside
        // tests and CASE-B keep fv_tolerance, so "contained" keeps one
        // meaning for exit qualification while a widened endpoint inset
        // (doc 32 tgm_fv_zmax_margin) no longer starves a wall-hugging
        // corner clipper of its midpoint support.  See inside_fv_interior().
        auto itol = config["interior_fv_tolerance"];
        if (!itol.isNull() && itol.isArray()) {
            m_interior_fv_tolerance.clear();
            for (const auto& t : itol) m_interior_fv_tolerance.push_back(t.asDouble());
        }
    }

    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping_name;
        cfg["detector_volumes"] = "DetectorVolumes";
        cfg["pc_transforms"] = "PCTransformSet";
        cfg["fiducial"] = "DetectorVolumes";
        // Boundary margins for the FV tests, FiducialUtils tolerance-vec
        // convention: [x_lo, x_hi, y_lo, y_hi, z_lo, z_hi], negative = inset.
        cfg["fv_tolerance"] = Json::Value(Json::arrayValue);
        cfg["interior_fv_tolerance"] = Json::Value(Json::arrayValue);
        cfg["beam_window_low"] = m_beam_window_low;   // window on cluster_t0; low >= high
        cfg["beam_window_high"] = m_beam_window_high; // disables the beam protection
        cfg["beam_window_only"] = m_beam_window_only; // evaluate only in-window mains
        cfg["length_limit_frac"] = m_length_limit_frac;
        cfg["enable_case_b"] = m_enable_case_b;
        cfg["require_in_scope"] = m_require_in_scope;
        cfg["evaluate_demoted_mains"] = m_evaluate_demoted_mains;
        cfg["check_neutrino_candidate"] = m_check_neutrino_candidate;
        cfg["require_chord_charge"] = m_require_chord_charge;
        cfg["chord_support_radius"] = m_chord_support_radius;
        cfg["chord_max_gap"] = m_chord_max_gap;
        cfg["chord_charge_mode"] = m_chord_charge_mode;
        cfg["component_extremes"] = m_component_extremes;
        cfg["component_min_length"] = m_component_min_length;
        cfg["component_graph"] = m_component_graph;
        cfg["component_rescue"] = m_component_rescue;
        cfg["rescue_chord_check"] = m_rescue_chord_check;
        cfg["main_component_pairs"] = m_main_component_pairs;
        cfg["main_component_mode"] = m_main_component_mode;
        cfg["exempt_demoted_main_pairs"] = m_exempt_demoted_main_pairs;
        return cfg;
    }

    virtual void visit(Ensemble& ensemble) const {
        auto groupings = ensemble.with_name(m_grouping_name);
        if (groupings.empty()) return;
        auto& grouping = *groupings.at(0);

        const bool beam_gate = m_beam_window_only && m_beam_window_low < m_beam_window_high;

        std::vector<Cluster*> main_clusters;
        int n_out_of_scope = 0;
        int n_out_of_window = 0;
        int n_demoted = 0;
        for (auto* cluster : grouping.children()) {
            const bool demoted = m_evaluate_demoted_mains
                && !cluster->get_flag(Flags::main_cluster)
                && cluster->get_flag(Flags::demoted_main);
            if (!cluster->get_flag(Flags::main_cluster) && !demoted) continue;
            if (m_require_in_scope && !cluster->get_scope_filter(cluster->get_default_scope())) {
                ++n_out_of_scope;
                continue;
            }
            if (beam_gate) {
                const double t0 = cluster->get_cluster_t0();
                if (t0 < m_beam_window_low || t0 >= m_beam_window_high) {
                    ++n_out_of_window;
                    continue;
                }
            }
            main_clusters.push_back(cluster);
            if (demoted) ++n_demoted;
        }
        if (n_demoted) {
            SPDLOG_LOGGER_INFO(t_log, "visit: TaggerCheckTGM: evaluate_demoted_mains: {} demoted main(s) added",
                               n_demoted);
        }
        if (n_out_of_scope) {
            SPDLOG_LOGGER_INFO(t_log, "visit: TaggerCheckTGM: skipped {} out-of-scope main cluster(s)",
                               n_out_of_scope);
        }
        if (beam_gate) {
            SPDLOG_LOGGER_INFO(t_log, "visit: TaggerCheckTGM: beam_window_only [{:.3f}, {:.3f}) us: {} main(s) evaluated, {} out of window",
                               m_beam_window_low/units::us, m_beam_window_high/units::us,
                               main_clusters.size(), n_out_of_window);
        }
        if (main_clusters.empty()) return;

        for (auto* main_cluster : main_clusters) {
            if (main_cluster->get_flag(Flags::TGM)) continue;
            bool is_tgm = false;
            try {
                is_tgm = check_tgm(*main_cluster);
            }
            catch (const std::exception& err) {
                SPDLOG_LOGGER_WARN(t_log, "visit: TaggerCheckTGM: cluster {} check failed: {}",
                                   main_cluster->ident(), err.what());
            }
            if (is_tgm) main_cluster->set_flag(Flags::TGM);
            SPDLOG_LOGGER_INFO(t_log, "visit: TaggerCheckTGM: cluster {} → TGM={}",
                               main_cluster->ident(), is_tgm);
        }
    }

private:
    std::string m_grouping_name{"live"};
    double m_beam_window_low{0};
    double m_beam_window_high{0};
    bool m_beam_window_only{false};
    double m_length_limit_frac{0.45};
    bool m_require_in_scope{false};
    bool m_evaluate_demoted_mains{false};
    bool m_enable_case_b{true};
    bool m_check_neutrino_candidate{false};
    bool m_require_chord_charge{false};
    double m_chord_support_radius{6 * units::cm};
    double m_chord_max_gap{30 * units::cm};
    std::string m_chord_charge_mode{"chord"};
    bool m_component_extremes{false};
    double m_component_min_length{10 * units::cm};
    std::string m_component_graph{"relaxed"};
    bool m_component_rescue{false};
    bool m_rescue_chord_check{false};
    bool m_main_component_pairs{false};
    std::string m_main_component_mode{"path"};
    bool m_exempt_demoted_main_pairs{false};
    std::vector<double> m_fv_tolerance;
    std::vector<double> m_interior_fv_tolerance;

    // A through-going muon deposits charge ALONG its whole path.  The CASE-A
    // and CASE-B pair tests below only ask whether the straight line between
    // the two extreme points crosses the fiducial volume -- never whether the
    // cluster has any charge on it.  That was safe in the prototype, where
    // check_tgm sees one connected PR3DCluster (bundle->get_main_cluster(),
    // 2dtoy/src/ToyFiducial.cxx:905-928).  Here the flash-time merge in
    // clustering_examine_bundles (use_flash_t0) collapses a whole flash group
    // into ONE Cluster, so a detached fragment hundreds of cm away is inside
    // the object under test and its chord to the real track passes the FV test
    // trivially.  SBND evt 285185 cluster 20: a 14-point / 1.1 cm speck on the
    // far top field cage, 608 cm chord, 96% of it empty, tagged TGM.
    //
    // Sample the chord ~every cm and reject the pair if any contiguous stretch
    // longer than m_chord_max_gap is unsupported, where "supported" means a
    // cluster point lies within m_chord_support_radius.  A cathode-crossing
    // track stays supported -- the CPA gap is a few cm, below the radius --
    // which a connectivity-based rule would NOT achieve: no connected
    // component ever spans the SBND cathode (the relaxed graph stops at
    // |x| ~ 1 cm), so "main component only" would drop every crosser.
    // Measured on the 10-event MCP2025C reco1 sample, support radius 6 cm:
    // genuine through-goers (incl. every cathode crosser) 0.0 cm longest
    // unsupported run, merge artefacts 93-583 cm.  Clean for radii >= 4 cm.
    bool chord_has_charge(const Cluster& cluster,
                          const geo_point_t& a, const geo_point_t& b) const {
        const double len = (b - a).magnitude();
        if (len <= 0) return true;
        const int nsamp = std::max(20, (int)(len / (1 * units::cm)));
        const double step = len / nsamp;
        const double r2 = m_chord_support_radius * m_chord_support_radius;
        double run = 0;  // running unsupported length
        for (int is = 0; is <= nsamp; ++is) {
            const geo_point_t p = a + (b - a) * ((double)is / nsamp);
            const auto res = cluster.kd_knn(1, p);
            // kd_knn returns (index, SQUARED distance), cf. TaggerCheckSTM.cxx:1497.
            if (!res.empty() && res[0].second <= r2) { run = 0; continue; }
            run += step;
            if (run > m_chord_max_gap) return false;
        }
        return true;
    }

    // Piecewise-path variant of the chord guard (chord_charge_mode: "path").
    //
    // The straight-chord test above rejects a pair whenever the STRAIGHT
    // segment between the extremes leaves the charge -- which a genuinely
    // curved track does: SBND evt285185 cluster 16, a continuous 480 cm
    // top->anode crosser (largest real point gap 3.1 cm), bows up to 10 cm
    // away from its own end-to-end chord, giving a 126 cm "unsupported" run
    // at the 6 cm support radius and losing a real TGM.  Path mode asks the
    // intended question directly: can pe1 reach pe2 through the cluster's
    // own points without a jump longer than m_chord_max_gap?  Robust to
    // curvature by construction.  The flash-merge artefacts stay rejected --
    // their grafted fragments sit >= 93 cm from the rest of the bundle's
    // charge on the 10-event MCP2025C sample -- and a cathode crosser still
    // passes: the ~3 cm CPA gap is far below the 30 cm linking distance.
    // (This is NOT the pre-built "relaxed"-graph connectivity, whose ~1 cm
    // linking stops at the cathode; the linking distance here is
    // m_chord_max_gap.)  m_chord_support_radius is not used in this mode.
    //
    // Implementation: single-linkage components over a voxel downsample of
    // the cluster points (cell = max_gap/6; the representative of a voxel is
    // its lowest-index point), linking representatives closer than max_gap.
    // Representatives of adjacent occupied voxels are at most 2 cell
    // diagonals ~= 17 cm apart at the 30 cm default, well under the linking
    // distance, so a contiguous track can never be split; a real gap breaks
    // the path at an effective tolerance of max_gap -+ 2 cell diagonals
    // ~= [13, 47] cm -- both edges far from the ~3 cm (genuine) and
    // >= 93 cm (artefact) populations this knob separates.  Deterministic:
    // point-index and voxel-index iteration only.
    //
    // comp[i] = component id of point i, -1 for excluded points.
    std::vector<int> path_components(const Cluster& cluster) const {
        const int npts = cluster.npoints();
        std::vector<int> comp(npts, -1);
        if (npts <= 0) return comp;
        const double cell = m_chord_max_gap / 6.0;
        std::map<std::array<long, 3>, int> vox;  // voxel key -> voxel index
        std::vector<geo_point_t> rep;            // per voxel: first (lowest-index) point
        std::vector<std::array<long, 3>> keys;   // per voxel: its key
        std::vector<int> pt_vox(npts, -1);       // per point: voxel index
        for (int i = 0; i < npts; ++i) {
            if (cluster.is_point_excluded(i)) continue;
            const geo_point_t p = cluster.point3d(i);
            const std::array<long, 3> key{(long)std::floor(p.x() / cell),
                                          (long)std::floor(p.y() / cell),
                                          (long)std::floor(p.z() / cell)};
            auto [it, fresh] = vox.try_emplace(key, (int)rep.size());
            if (fresh) {
                rep.push_back(p);
                keys.push_back(key);
            }
            pt_vox[i] = it->second;
        }
        const int nvox = (int)rep.size();
        if (!nvox) return comp;
        std::vector<int> parent(nvox);
        std::iota(parent.begin(), parent.end(), 0);
        auto find = [&parent](int a) {
            while (parent[a] != a) {
                parent[a] = parent[parent[a]];
                a = parent[a];
            }
            return a;
        };
        const long nspan = 6;  // = max_gap / cell by construction
        const double gap2 = m_chord_max_gap * m_chord_max_gap;
        for (int a = 0; a < nvox; ++a) {
            for (int b = a + 1; b < nvox; ++b) {
                if (std::abs(keys[a][0] - keys[b][0]) > nspan) continue;
                if (std::abs(keys[a][1] - keys[b][1]) > nspan) continue;
                if (std::abs(keys[a][2] - keys[b][2]) > nspan) continue;
                const int ra = find(a), rb = find(b);
                if (ra == rb) continue;
                const geo_vector_t d = rep[a] - rep[b];
                if (d.dot(d) > gap2) continue;
                parent[std::max(ra, rb)] = std::min(ra, rb);
            }
        }
        for (int i = 0; i < npts; ++i) {
            if (pt_vox[i] >= 0) comp[i] = find(pt_vox[i]);
        }
        return comp;
    }

    // Are a and b (extreme points, i.e. actual cluster points) in the same
    // charge-path component?  Fails OPEN -- the guard may only suppress.
    bool path_connected(const Cluster& cluster, const std::vector<int>& comp,
                        const geo_point_t& a, const geo_point_t& b) const {
        if (comp.empty()) return true;
        const auto ra = cluster.kd_knn(1, a);
        const auto rb = cluster.kd_knn(1, b);
        if (ra.empty() || rb.empty()) return true;
        const int ca = comp[ra[0].first];
        const int cb = comp[rb[0].first];
        if (ca < 0 || cb < 0) return true;  // nearest point excluded: don't veto
        return ca == cb;
    }

    // Component-aware extreme points (knob off => never called).
    //
    // get_extreme_wcps() (Facade_Cluster.cxx:3042) scans the WHOLE cluster for
    // 8 GLOBAL extremes.  On a flash-merged Cluster each slot (max z, min x,
    // ...) is claimed by whichever merge component reaches furthest, so a
    // second component's own wall-exit never becomes a candidate and the
    // legitimate WITHIN-component through-going pair is never even tested:
    //   SBND evt284657 cluster 25 -- comp 11 runs top wall (y=199.9) to
    //   downstream wall (z=500.5), but comp 10 also reaches z=500.5 and took
    //   the global max-z slot, leaving comp 11 with only its top end;
    //   cluster 26 -- comp 20 runs top wall to the anode (x=-201.3), but comp
    //   16 also reaches x=-201.3 and took the global min-x slot.
    // Both are real through-going muons that lost their tag.
    //
    // Fix: run the same scan PER connected component and union the results.
    // Components shorter than m_component_min_length contribute nothing (a
    // few-point speck would otherwise donate two "ends" of its own).  Cross-
    // component pairs are still formed -- they are exactly what the
    // chord-charge test rejects -- so a cathode crosser, whose two halves ARE
    // separate components here, keeps its tag: its chord across the ~3 cm CPA
    // gap stays charge-supported.  Use the two knobs together.
    //
    // The scan and the 5 cm proximity grouping are duplicated from
    // get_extreme_wcps() rather than shared (CLAUDE.md M10).
    bool component_extremes(const Cluster& cluster, const std::vector<size_t>& idxs,
                            std::array<geo_point_t, 8>& ex) const {
        if (idxs.size() < 2) return false;
        std::vector<geo_point_t> pts;
        pts.reserve(idxs.size());
        geo_point_t center(0, 0, 0);
        for (size_t i : idxs) {
            const geo_point_t p = cluster.point3d(i);
            pts.push_back(p);
            center = center + p;
        }
        center = center * (1.0 / (double)pts.size());
        geo_vector_t main_axis = calc_pca_dir(center, pts);
        if (main_axis.y() < 0) main_axis = main_axis * (-1);

        double val[8];
        for (int i = 0; i < 8; ++i) ex[i] = pts[0];
        val[0] = val[1] = pts[0].dot(main_axis);
        val[2] = val[3] = pts[0].y();
        val[4] = val[5] = pts[0].z();
        val[6] = val[7] = pts[0].x();
        for (const auto& p : pts) {
            const double mp = p.dot(main_axis);
            if (mp > val[0]) { ex[0] = p; val[0] = mp; }
            if (mp < val[1]) { ex[1] = p; val[1] = mp; }
            if (p.y() > val[2]) { ex[2] = p; val[2] = p.y(); }
            if (p.y() < val[3]) { ex[3] = p; val[3] = p.y(); }
            if (p.z() > val[4]) { ex[4] = p; val[4] = p.z(); }
            if (p.z() < val[5]) { ex[5] = p; val[5] = p.z(); }
            if (p.x() > val[6]) { ex[6] = p; val[6] = p.x(); }
            if (p.x() < val[7]) { ex[7] = p; val[7] = p.x(); }
        }
        return (ex[0] - ex[1]).magnitude() >= m_component_min_length;
    }

    // rescued_out (optional): parallel to the returned groups, one flag per
    // point, 1 = the point was donated by a rescue-pass component.  Consumed
    // by the rescue_chord_check pair guard in check_tgm().
    std::vector<std::vector<geo_point_t>> component_extreme_wcps(
        const Cluster& cluster,
        std::vector<std::vector<char>>* rescued_out = nullptr) const {
        std::vector<std::vector<geo_point_t>> out;
        std::vector<std::vector<char>> res;
        const int npts = cluster.npoints();
        if (npts <= 0) return out;
        std::vector<int> b2c;
        try {
            b2c = cluster.connected_blobs(m_dv, m_pcts, m_component_graph);
        }
        catch (const std::exception& err) {
            SPDLOG_LOGGER_WARN(t_log, "component_extreme_wcps: cluster {} connected_blobs failed: {}",
                               cluster.ident(), err.what());
            return out;
        }
        if (b2c.empty()) return out;

        // std::map keyed by component id => deterministic iteration (no pointer keys).
        std::map<int, std::vector<size_t>> comp_pts;
        for (int i = 0; i < npts; ++i) {
            if (cluster.is_point_excluded(i)) continue;
            const int bind = cluster.kd3d().major_index(i);
            if (bind < 0 || bind >= (int)b2c.size()) continue;
            comp_pts[b2c[bind]].push_back((size_t)i);
        }
        // Largest component first so ITS main-axis pair seeds groups 0 and 1:
        // CASE B special-cases (i==0 && k==1) as "the" main-axis pair.  Ties
        // break on component id (stable_sort over the ordered map).
        std::vector<const std::vector<size_t>*> order;
        for (const auto& [cid, v] : comp_pts) order.push_back(&v);
        std::stable_sort(order.begin(), order.end(),
                         [](const std::vector<size_t>* a, const std::vector<size_t>* b) {
                             return a->size() > b->size();
                         });

        const double min_separation = 5.0 * units::cm;
        int n_used = 0;
        std::vector<const std::vector<size_t>*> kept;     // rescue anchors (knob on)
        std::vector<const std::vector<size_t>*> skipped;  // rescue candidates (knob on)
        for (const auto* v : order) {
            std::array<geo_point_t, 8> ex;
            if (!component_extremes(cluster, *v, ex)) {
                // >= 2 points so a rescue re-scan can fill ex.
                if (m_component_rescue && v->size() >= 2) skipped.push_back(v);
                continue;
            }
            ++n_used;
            if (m_component_rescue) kept.push_back(v);
            int first = 0;
            if (out.empty()) {
                out.push_back({ex[0]});
                out.push_back({ex[1]});
                res.push_back({0});
                res.push_back({0});
                first = 2;
            }
            for (int i = first; i < 8; ++i) {
                bool distinct = true;
                for (size_t g = 0; g < out.size(); ++g) {
                    if ((ex[i] - out[g][0]).magnitude() < min_separation) {
                        out[g].push_back(ex[i]);
                        res[g].push_back(0);
                        distinct = false;
                        break;
                    }
                }
                if (distinct) {
                    out.push_back({ex[i]});
                    res.push_back({0});
                }
            }
        }
        // Rescue pass (knob off => never taken): a sub-min-length component
        // donates its extremes after all when it shares a path_components()
        // component (30 cm-step charge path, same rule as the path-mode chord
        // guard) with a component that passed the length cut.  This cannot
        // re-admit the evt285185 speck: that fragment is path-DISCONNECTED
        // from the track by construction (>= 93 cm from any other charge).
        // All points of a "relaxed"-graph component (~1 cm linking) share one
        // path component (30 cm linking), so any single point identifies it.
        if (m_component_rescue && n_used > 0 && !skipped.empty()) {
            const std::vector<int> pcomp = path_components(cluster);
            std::set<int> anchors;
            for (const auto* v : kept) {
                const int pid = pcomp[(*v)[0]];
                if (pid >= 0) anchors.insert(pid);
            }
            int n_rescued = 0;
            for (const auto* v : skipped) {
                const int pid = pcomp[(*v)[0]];
                if (pid < 0 || !anchors.count(pid)) continue;
                std::array<geo_point_t, 8> ex;
                component_extremes(cluster, *v, ex);  // fills ex; length verdict ignored
                ++n_rescued;
                for (int i = 0; i < 8; ++i) {
                    bool distinct = true;
                    for (size_t g = 0; g < out.size(); ++g) {
                        if ((ex[i] - out[g][0]).magnitude() < min_separation) {
                            out[g].push_back(ex[i]);
                            res[g].push_back(1);
                            distinct = false;
                            break;
                        }
                    }
                    if (distinct) {
                        out.push_back({ex[i]});
                        res.push_back({1});
                    }
                }
            }
            if (n_rescued) {
                SPDLOG_LOGGER_DEBUG(t_log,
                    "component_rescue: cluster {} rescued {} of {} short component(s) path-connected to a kept one",
                    cluster.ident(), n_rescued, skipped.size());
            }
        }
        SPDLOG_LOGGER_DEBUG(t_log, "component_extreme_wcps: cluster {} {} component(s), {} above {:.1f} cm -> {} extreme group(s)",
                            cluster.ident(), comp_pts.size(), n_used,
                            m_component_min_length / units::cm, out.size());
        if (rescued_out) *rescued_out = std::move(res);
        return out;
    }

    // FiducialUtils::inside_fiducial_volume logic against our own IFiducial.
    bool inside_fv(const Point& p) const { return inside_fv_tol(p, m_fv_tolerance); }
    // Interior-support variant: the CASE-A "does the track pass through the
    // detector interior" tests (flag_check chord midpoints and the
    // flag_check_again waypoint chords) may use their OWN tolerance vector.
    // A widened endpoint inset (tgm_fv_zmax_margin 3->5 cm, doc 32) also
    // shrinks the interior these tests sample, so a genuine corner clipper
    // RUNNING ALONG the downstream wall inside the widened band loses its
    // midpoint support and the waypoint re-check vetoes the tag (SBND
    // evt287517 cluster 16, evt289805 cluster 9: top->downstream-z corner
    // clippers, whole chord in the z = 493-501 cm band).  Empty (default)
    // => fall back to fv_tolerance, byte-identical legacy behavior.
    bool inside_fv_interior(const Point& p) const {
        return inside_fv_tol(p, m_interior_fv_tolerance.empty() ? m_fv_tolerance
                                                                : m_interior_fv_tolerance);
    }
    bool inside_fv_tol(const Point& p, const std::vector<double>& tv) const {
        if (tv.empty()) return m_fiducial->contained(p);
        double txl, txh, tyl, tyh, tzl, tzh;
        if (tv.size() >= 6) { txl = tv[0]; txh = tv[1]; tyl = tv[2]; tyh = tv[3]; tzl = tv[4]; tzh = tv[5]; }
        else if (tv.size() >= 3) { txl = txh = tv[0]; tyl = tyh = tv[1]; tzl = tzh = tv[2]; }
        else { txl = txh = tyl = tyh = tzl = tzh = tv[0]; }
        if (!m_fiducial->contained(Point(p.x() - txl, p.y(), p.z()))) return false;
        if (!m_fiducial->contained(Point(p.x() + txh, p.y(), p.z()))) return false;
        if (!m_fiducial->contained(Point(p.x(), p.y() - tyl, p.z()))) return false;
        if (!m_fiducial->contained(Point(p.x(), p.y() + tyh, p.z()))) return false;
        if (!m_fiducial->contained(Point(p.x(), p.y(), p.z() - tzl))) return false;
        if (!m_fiducial->contained(Point(p.x(), p.y(), p.z() + tzh))) return false;
        return true;
    }

    bool check_tgm(Cluster& cluster) const {
        auto fiducial_utils = cluster.grouping()->get_fiducialutils();
        if (!fiducial_utils) {
            SPDLOG_LOGGER_WARN(t_log, "check_tgm: no FiducialUtils on grouping; run fiducialutils first");
            return false;
        }

        // Component-aware extremes when the knob is on; fall back to the
        // global scan if it yields fewer than two groups (tiny cluster), so
        // coverage is never LOST relative to the knob-off path.
        std::vector<std::vector<char>> rescued_grps;
        auto out_vec_wcps = m_component_extremes ? component_extreme_wcps(cluster, &rescued_grps)
                                                 : std::vector<std::vector<geo_point_t>>();
        if (out_vec_wcps.size() < 2) {
            out_vec_wcps = cluster.get_extreme_wcps();
            rescued_grps.clear();
        }
        if (out_vec_wcps.size() < 2 || out_vec_wcps[0].empty() || out_vec_wcps[1].empty()) return false;
        for (const auto& grp : out_vec_wcps) {
            if (grp.empty()) return false;
        }

        // Beam protection: the prototype gates its tags on
        // main_flash->get_type()==2 (beam flash) and check_neutrino_candidate.
        // Here: beam-coincident bundle == cluster_t0 in the window; the
        // protected branches never tag (check_neutrino_candidate unported).
        const bool in_beam_window = (m_beam_window_low < m_beam_window_high)
            && cluster.get_cluster_t0() >= m_beam_window_low
            && cluster.get_cluster_t0() < m_beam_window_high;

        const double length_limit = (out_vec_wcps[0][0] - out_vec_wcps[1][0]).magnitude();

        // Per-(apa,face) U/V/W wire directions for the prolonged-signal test.
        std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> wpid_params;
        std::map<WirePlaneId, std::pair<geo_point_t, double>> wpid_U_dir, wpid_V_dir, wpid_W_dir;
        std::set<int> apas;
        Facade::compute_wireplane_params(cluster.grouping()->wpids(), m_dv,
                                         wpid_params, wpid_U_dir, wpid_V_dir, wpid_W_dir, apas);
        // The |dir.x| construction makes the drift-angle test sign-agnostic,
        // so a single abs drift axis serves both SBND drift directions.
        const geo_vector_t drift_dir_abs(1, 0, 0);
        if (wpid_U_dir.empty()) return false;
        auto uvw_dirs_at = [&](const Point& p) {
            WirePlaneId wpid = m_dv->contained_by(p);
            auto pick = [&](const std::map<WirePlaneId, std::pair<geo_point_t, double>>& m) {
                if (wpid.apa() >= 0 && wpid.face() >= 0) {
                    for (const auto& [k, v] : m) {
                        if (k.apa() == wpid.apa() && k.face() == wpid.face()) return v.first;
                    }
                }
                return m.begin()->second.first;
            };
            return std::array<geo_vector_t, 3>{pick(wpid_U_dir), pick(wpid_V_dir), pick(wpid_W_dir)};
        };
        // Prototype prolonged-signal angle: project dir into the (y,z) plane,
        // measure against the wire direction, rebuild the (|dx|, transverse)
        // vector and take its angle to the drift axis, in degrees.
        auto drift_angle_deg = [&](const geo_vector_t& dir, const geo_vector_t& wire_dir) {
            geo_vector_t dir_1(0, dir.y(), dir.z());
            const double angle = dir_1.angle(wire_dir);
            geo_vector_t tempV(std::fabs(dir.x()),
                               std::sqrt(dir.y()*dir.y() + dir.z()*dir.z()) * std::sin(angle), 0);
            return tempV.angle(drift_dir_abs) / 3.1415926 * 180.;
        };

        const geo_point_t dir_main = cluster.get_pca().axis.at(0);

        // Chord-charge guard state + dispatch (knob off => never taken).  The
        // path-mode components are computed once per cluster; the lambda
        // keeps the two call sites terse and the chord-mode log line
        // byte-identical to the pre-path-mode code.
        const std::vector<int> path_comp =
            ((m_require_chord_charge && m_chord_charge_mode == "path") || m_main_component_pairs)
                ? path_components(cluster) : std::vector<int>();
        auto chord_guard_rejects = [&](const geo_point_t& a, const geo_point_t& b,
                                       const char* case_tag, size_t gi, size_t gk) {
            if (!m_require_chord_charge) return false;
            if (m_chord_charge_mode == "path") {
                if (path_connected(cluster, path_comp, a, b)) return false;
                SPDLOG_LOGGER_DEBUG(t_log, "check_tgm: cluster {} {} pair ({},{}) rejected: no {:.1f} cm-step charge path between the ends ({:.1f} cm chord)",
                                    cluster.ident(), case_tag, gi, gk,
                                    m_chord_max_gap / units::cm,
                                    (a - b).magnitude() / units::cm);
                return true;
            }
            if (chord_has_charge(cluster, a, b)) return false;
            SPDLOG_LOGGER_DEBUG(t_log, "check_tgm: cluster {} {} pair ({},{}) rejected: chord {:.1f} cm has an unsupported run > {:.1f} cm",
                                cluster.ident(), case_tag, gi, gk,
                                (a - b).magnitude() / units::cm,
                                m_chord_max_gap / units::cm);
            return true;
        };
        // rescue_chord_check pair guard (knob off => never taken): a pair
        // whose end point was donated by a RESCUED component must also pass
        // the STRAIGHT-chord support test.  Path mode alone lets a rescued
        // wall-touching speck pair with the other cosmic of a two-track
        // merge through an L-shaped charge detour (evt288727 cluster 6);
        // a genuine fragmented track end lies on its own pair's chord.
        auto end_rescued = [&](size_t g, int idx) -> bool {
            if (g >= rescued_grps.size()) return false;
            const auto& v = rescued_grps[g];
            if (idx < 0 || idx >= (int)v.size()) return false;
            return v[idx] != 0;
        };
        auto rescued_pair_rejects = [&](const geo_point_t& a, const geo_point_t& b,
                                        size_t gi, int ai, size_t gk, int bi,
                                        const char* case_tag) {
            if (!m_rescue_chord_check || rescued_grps.empty()) return false;
            if (!end_rescued(gi, ai) && !end_rescued(gk, bi)) return false;
            if (chord_has_charge(cluster, a, b)) return false;
            SPDLOG_LOGGER_DEBUG(t_log, "check_tgm: cluster {} {} pair ({},{}) rejected: rescued end, straight chord {:.1f} cm has an unsupported run > {:.1f} cm",
                                cluster.ident(), case_tag, gi, gk,
                                (a - b).magnitude() / units::cm,
                                m_chord_max_gap / units::cm);
            return true;
        };
        // main_component_pairs guard (knob off => never taken): the MAIN
        // charge component is the path_components() label holding the most
        // points; ascending-id iteration with a strict > keeps the lowest id
        // on a tie (deterministic).  A pair with NEITHER end in it is a
        // fragment-only pair -- on a flash-merged cluster that is a merged-in
        // cosmic tagging the bundle, not the main track (evt289343 cluster 9).
        // A genuine main-track pair has both ends in it (a cathode crosser is
        // ONE path component -- the 30 cm linkage bridges the CPA gap), and a
        // main<->fragment cross pair keeps one end, so only the fragment's
        // own pairs are vetoed.  Fails OPEN like path_connected().
        int main_comp = -1;
        if (m_main_component_pairs && !path_comp.empty()) {
            std::map<int, int> comp_npts;
            for (int c : path_comp) {
                if (c >= 0) ++comp_npts[c];
            }
            int best_npts = 0;
            for (const auto& [c, n] : comp_npts) {
                if (n > best_npts) { best_npts = n; main_comp = c; }
            }
        }
        // Exact mode ("real"): per-blob provenance from the QL flash merge --
        // the contributing member's "main_cluster" flag, persisted through the
        // pctree tarball by save_real_cluster_id.  The array is parallel to
        // the cluster's blob (children) order, which is also the scoped-view
        // blob order that kd().major_index() indexes (same invariant the Bee
        // writer relies on).  Absent array (old tarballs) => fall back to the
        // largest-path-component proxy.
        std::vector<int> real_main;
        if (m_main_component_pairs && m_main_component_mode == "real"
            && cluster.has_pcarray<int>("real_cluster_main", "perblob")) {
            auto sp = cluster.get_pcarray<int>("real_cluster_main", "perblob");
            real_main.assign(sp.begin(), sp.end());
        }
        // 1 = end is in the main cluster/component, 0 = definitely not,
        // -1 = cannot tell (fails OPEN: only a definite 0/0 pair is vetoed).
        auto end_in_main = [&](const geo_point_t& p) -> int {
            const auto r = cluster.kd_knn(1, p);
            if (r.empty()) return -1;
            if (!real_main.empty()) {
                const size_t bidx = cluster.kd().major_index(r[0].first);
                if (bidx >= real_main.size()) return -1;
                return real_main[bidx] != 0 ? 1 : 0;
            }
            if (main_comp < 0) return -1;
            return path_comp[r[0].first] == main_comp ? 1 : 0;
        };
        auto main_pair_rejects = [&](const geo_point_t& a, const geo_point_t& b,
                                     const char* case_tag, size_t gi, size_t gk) {
            if (!m_main_component_pairs) return false;
            if (m_exempt_demoted_main_pairs && cluster.get_flag(Flags::demoted_main)) return false;
            if (real_main.empty() && main_comp < 0) return false;
            if (end_in_main(a) != 0 || end_in_main(b) != 0) return false;
            SPDLOG_LOGGER_DEBUG(t_log, "check_tgm: cluster {} {} pair ({},{}) rejected: neither end in the {} ({:.1f} cm chord)",
                                cluster.ident(), case_tag, gi, gk,
                                real_main.empty() ? "main charge component"
                                                  : "pre-merge main cluster",
                                (a - b).magnitude() / units::cm);
            return true;
        };

        for (size_t i = 0; i != out_vec_wcps.size(); i++) {
            bool flag_p1_inside = true;
            int p1_index = -1;
            for (size_t j = 0; j != out_vec_wcps[i].size(); j++) {
                flag_p1_inside = flag_p1_inside && inside_fv(out_vec_wcps[i][j]);
                if (!flag_p1_inside) { p1_index = j; break; }
            }

            for (size_t k = i + 1; k != out_vec_wcps.size(); k++) {
                bool flag_p2_inside = true;
                int p2_index = -1;
                for (size_t j = 0; j != out_vec_wcps[k].size(); j++) {
                    flag_p2_inside = flag_p2_inside && inside_fv(out_vec_wcps[k][j]);
                    if (!flag_p2_inside) { p2_index = j; break; }
                }

                if (!flag_p1_inside && !flag_p2_inside) {
                    // CASE A: both ends exit the FV.
                    const geo_point_t& pe1 = out_vec_wcps[i][p1_index];
                    const geo_point_t& pe2 = out_vec_wcps[k][p2_index];
                    bool flag_check = false;
                    for (int kk = 0; kk != 3; kk++) {
                        const geo_point_t p3 = pe1 + (pe2 - pe1) * ((kk + 1) / 4.);
                        flag_check = flag_check || inside_fv_interior(p3);
                    }
                    if (std::getenv("WCT_TGM_DEBUG")) {
                        SPDLOG_LOGGER_INFO(t_log, "check_tgm dbg: cluster {} pair ({},{}) ngrp {} pe1 ({:.1f},{:.1f},{:.1f}) pe2 ({:.1f},{:.1f},{:.1f}) mid_inside {} len {:.1f}/{:.1f} cm",
                            cluster.ident(), i, k, out_vec_wcps.size(),
                            pe1.x()/units::cm, pe1.y()/units::cm, pe1.z()/units::cm,
                            pe2.x()/units::cm, pe2.y()/units::cm, pe2.z()/units::cm,
                            flag_check, (pe1-pe2).magnitude()/units::cm, length_limit/units::cm);
                    }
                    // Chord-charge guard (knob off => never taken).  Placed
                    // before every CASE-A tagging branch so it covers both the
                    // flag_check path and the out_vec_wcps.size()==2 /
                    // flag_check_again early exits below.  It can only SUPPRESS
                    // a tag; the loop continues to the remaining pairs, so a
                    // sub-track that is itself through-going still tags on its
                    // own charge-supported pair.
                    if (chord_guard_rejects(pe1, pe2, "CASE-A", i, k)) continue;
                    if (rescued_pair_rejects(pe1, pe2, i, p1_index, k, p2_index, "CASE-A")) continue;
                    if (main_pair_rejects(pe1, pe2, "CASE-A", i, k)) continue;
                    if (flag_check) {
                        if (in_beam_window) {
                            if (m_check_neutrino_candidate) {
                                // prototype (Cosmic_tagger.h lines 1421-1437): beam flash →
                                // tag only when not a neutrino candidate and long enough.
                                const double temp_length = (pe1 - pe2).magnitude();
                                if (!check_neutrino_candidate(cluster, pe1, pe2)
                                    && temp_length > m_length_limit_frac * length_limit) return true;
                                continue;
                            }
                            // knob off: conservative NO TAG (v1 behavior).
                            continue;
                        }
                        return true;
                    }
                    else {
                        if (out_vec_wcps.size() == 2) return true;
                        // Re-check chords through the other extreme groups.
                        bool flag_check_again = false;
                        for (size_t kkk = 0; kkk != out_vec_wcps.size(); kkk++) {
                            if (kkk == i || kkk == k) continue;
                            for (int kk = 0; kk != 4; kk++) {
                                const geo_point_t p3 = pe1 + (out_vec_wcps[kkk][0] - pe1) * ((kk + 1) / 4.);
                                flag_check_again = flag_check_again || inside_fv_interior(p3);
                            }
                            for (int kk = 0; kk != 3; kk++) {
                                // NB: mixed endpoints as in the prototype.
                                const geo_point_t p3 = out_vec_wcps[kkk][0]
                                    + (pe2 - out_vec_wcps[i][0]) * ((kk + 1) / 4.);
                                flag_check_again = flag_check_again || inside_fv_interior(p3);
                            }
                        }
                        if (!flag_check_again) {
                            const double temp_length = (pe1 - pe2).magnitude();
                            if (m_check_neutrino_candidate) {
                                // prototype (Cosmic_tagger.h lines 1466-1477): NO beam gate
                                // here -- the neutrino-candidate veto arbitrates for every
                                // cluster, in or out of the beam window.
                                if (!check_neutrino_candidate(cluster, pe1, pe2)
                                    && temp_length > m_length_limit_frac * length_limit) return true;
                            }
                            else {
                                if (in_beam_window) continue;  // v1: needs check_neutrino_candidate
                                if (temp_length > m_length_limit_frac * length_limit) return true;
                            }
                        }
                    }
                }
                else if (m_enable_case_b) {
                    // CASE B: at least one end looks inside — test whether it is
                    // really a prolonged-signal artefact or exits into a dead region.
                    const geo_point_t p1 = out_vec_wcps[i][0];
                    const geo_point_t p2 = out_vec_wcps[k][0];
                    geo_vector_t dir_test = p1 - p2;

                    const double perp_deg =
                        std::fabs((3.1415926/2. - dir_test.angle(dir_main)) / 3.1415926 * 180.);
                    if (!(perp_deg > 75 || (i == 0 && k == 1))) continue;

                    // Chord-charge guard, CASE B (knob off => never taken).
                    // Before the hough / dead-volume work so it also saves that cost.
                    if (chord_guard_rejects(p1, p2, "CASE-B", i, k)) continue;
                    if (rescued_pair_rejects(p1, p2, i, 0, k, 0, "CASE-B")) continue;
                    if (main_pair_rejects(p1, p2, "CASE-B", i, k)) continue;

                    bool skip_pair = false;
                    bool flag_p1_inside_p = flag_p1_inside;
                    if (flag_p1_inside_p) {
                        geo_vector_t dir = cluster.vhough_transform(p1, 30*units::cm);
                        dir = dir * (-1);
                        if (dir.angle(dir_test) > 3.1415926*2./3.) { skip_pair = true; }
                        else {
                            const auto uvw = uvw_dirs_at(p1);
                            if (drift_angle_deg(dir, uvw[0]) < 10 || drift_angle_deg(dir, uvw[1]) < 10
                                || drift_angle_deg(dir, uvw[2]) < 5) {
                                flag_p1_inside_p = flag_p1_inside_p
                                    && fiducial_utils->check_signal_processing(cluster, p1, dir, 1*units::cm);
                            }
                            if (std::fabs((3.1415926/2. - dir.angle(dir_main)) / 3.1415926 * 180.) > 60) {
                                flag_p1_inside_p = flag_p1_inside_p
                                    && fiducial_utils->check_dead_volume(cluster, p1, dir, 1*units::cm);
                            }
                        }
                    }
                    if (skip_pair) continue;

                    bool flag_p2_inside_p = flag_p2_inside;
                    if (flag_p2_inside_p) {
                        geo_vector_t dir = cluster.vhough_transform(p2, 30*units::cm);
                        dir = dir * (-1);
                        if (dir.angle(dir_test) < 3.1415926/3.) { skip_pair = true; }
                        else {
                            const auto uvw = uvw_dirs_at(p2);
                            if (drift_angle_deg(dir, uvw[0]) < 10 || drift_angle_deg(dir, uvw[1]) < 10
                                || drift_angle_deg(dir, uvw[2]) < 5) {
                                flag_p2_inside_p = flag_p2_inside_p
                                    && fiducial_utils->check_signal_processing(cluster, p2, dir, 1*units::cm);
                            }
                            if (std::fabs((3.1415926/2. - dir.angle(dir_main)) / 3.1415926 * 180.) > 60) {
                                flag_p2_inside_p = flag_p2_inside_p
                                    && fiducial_utils->check_dead_volume(cluster, p2, dir, 1*units::cm);
                            }
                        }
                    }
                    if (skip_pair) continue;

                    if (!flag_p1_inside_p && !flag_p2_inside_p) {
                        if (in_beam_window) {
                            if (m_check_neutrino_candidate) {
                                // prototype (Cosmic_tagger.h lines 1564-1581): beam flash →
                                // veto arbitration on the first extreme points of the pair.
                                const double temp_length = (p1 - p2).magnitude();
                                if (!check_neutrino_candidate(cluster, p1, p2)
                                    && temp_length > m_length_limit_frac * length_limit) return true;
                            }
                            continue;  // knob off: v1 conservative NO TAG
                        }
                        return true;
                    }
                }
            }
        }
        return false;
    }

    // Port of the prototype's Dijkstra path-topology neutrino veto,
    // WCPPID::ToyFiducial::check_neutrino_candidate
    // (prototype pid/src/Cosmic_tagger.h lines 1677-1985).  True when the
    // in-cluster shortest path between the two extreme points shows a
    // neutrino-like topology:
    //  (a) "gap" veto: >7 consecutive path samples inside the FV close
    //      (<25 cm) to an endpoint with no 2-view ctpc support and >7 bad
    //      points, not explained by a dead region;
    //  (b) "kink" veto: a sustained (>=3-sample, or 1 with tighter length
    //      cuts) direction break whose apex is inside the FV and not in a
    //      dead region.
    //
    // Built on the existing toolkit ports: graph_algorithms("ctpc") is the
    // prototype's Create_graph + dijkstra_shortest_paths/cal_shortest_path;
    // Grouping::get_closest_points / is_good_point are the ToyCTPointCloud
    // queries (same defaults 0.6 cm / ch_range 1 / allowed_bad 1);
    // FiducialUtils::inside_dead_region is the prototype's dead-region test;
    // calc_pca_dir is the prototype's calc_PCA_dir.
    //
    // Two-TPC generalizations vs the single-TPC prototype:
    //  - the path and all FV tests run in the T0-corrected default scope
    //    (offset_x = 0 as everywhere in this component); every per-point
    //    ctpc/dead query first resolves the point's (apa, face) via
    //    DetectorVolumes and backward-transforms into that volume's raw
    //    coordinates (the FiducialUtils::check_signal_processing pattern);
    //  - path points outside every sensitive volume (e.g. the cathode gap
    //    of a two-TPC crosser path) are treated like dead regions: no hits
    //    can exist there, so they reset the gap counters;
    //  - the prolonged-signal angle test uses the per-(apa,face) wire
    //    directions of the endpoints' volumes (prototype: hard-coded uBooNE
    //    U/V/W); when the endpoints are in different volumes either frame
    //    may disable the 2-view requirement;
    //  - drift_dir = (1,0,0) is kept from the prototype: every use is of the
    //    form |90° − angle(drift, v)| or |dir.x|, both invariant under the
    //    drift-sign flip between the two drift volumes.
    bool check_neutrino_candidate(Cluster& cluster,
                                  const geo_point_t& wcp1, const geo_point_t& wcp2,
                                  bool flag_2view_check = true) const {
        auto* grouping = cluster.grouping();
        auto fiducial_utils = grouping->get_fiducialutils();
        if (!fiducial_utils) return false;

        // prototype: Create_graph(ct_point_cloud); dijkstra_shortest_paths(wcp1);
        // cal_shortest_path(wcp2).
        const size_t src_idx = cluster.get_closest_point_index(wcp1);
        const size_t dst_idx = cluster.get_closest_point_index(wcp2);
        const auto& path_indices =
            cluster.graph_algorithms("ctpc", m_dv, m_pcts).shortest_path(src_idx, dst_idx);
        const auto path_pts = cluster.indices_to_points(path_indices);
        if (path_pts.empty()) return false;

        // Resample the path (prototype lines 1687-1722): path_wcps_vec keeps
        // points >0.5 cm apart; path_wcps_vec1 caps spacing at 1 cm by
        // interpolation (NB the prototype's back() is re-read after each push,
        // giving its characteristic shrinking-step interpolation -- kept).
        const double low_dis_limit = 0.5 * units::cm;
        std::vector<geo_point_t> path_wcps_vec;
        std::vector<geo_point_t> path_wcps_vec1;
        for (const auto& pt : path_pts) {
            if (path_wcps_vec.empty()) {
                path_wcps_vec.push_back(pt);
                path_wcps_vec1.push_back(pt);
                continue;
            }
            double dis = (pt - path_wcps_vec.back()).magnitude();
            if (dis > low_dis_limit) path_wcps_vec.push_back(pt);
            dis = (pt - path_wcps_vec1.back()).magnitude();
            if (dis <= 2 * low_dis_limit) {
                path_wcps_vec1.push_back(pt);
            }
            else {
                const int nseg = dis / 2. / low_dis_limit + 1;
                for (int i = 0; i != nseg; i++) {
                    const geo_point_t temp_p = path_wcps_vec1.back()
                        + (pt - path_wcps_vec1.back()) * ((i + 1.) / nseg);
                    path_wcps_vec1.push_back(temp_p);
                }
            }
        }

        // Raw-coordinate transform for the per-point ctpc/dead queries
        // (FiducialUtils::check_signal_processing pattern).
        const auto transform =
            m_pcts->pc_transform(cluster.get_scope_transform(cluster.get_default_scope()));
        const double cluster_t0 = cluster.get_cluster_t0();
        auto to_raw = [&](const geo_point_t& p, const WirePlaneId& wpid) {
            return transform->backward(p, cluster_t0, wpid.face(), wpid.apa());
        };

        // Check whether the path is good (prototype lines 1725-1830).
        {
            int num_nth = 0;
            double min_dis = 1e9;

            // U/V/W induction-view check: drop the 2-view requirement when the
            // chord is quasi-parallel to a wire plane or quasi-isochronous.
            if (flag_2view_check) {
                std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> wpid_params;
                std::map<WirePlaneId, std::pair<geo_point_t, double>> wpid_U_dir, wpid_V_dir, wpid_W_dir;
                std::set<int> apas;
                Facade::compute_wireplane_params(grouping->wpids(), m_dv,
                                                 wpid_params, wpid_U_dir, wpid_V_dir, wpid_W_dir, apas);
                if (!wpid_U_dir.empty()) {
                    const geo_vector_t drift_dir_abs(1, 0, 0);
                    const geo_vector_t dir = wcp2 - wcp1;
                    auto angle_deg = [&](const geo_vector_t& wire_dir) {
                        // prototype lines 1743-1755 (same construction as
                        // check_tgm's prolonged-signal test).
                        geo_vector_t dir_1(0, dir.y(), dir.z());
                        const double angle = dir_1.angle(wire_dir);
                        geo_vector_t tempV(std::fabs(dir.x()),
                                           std::sqrt(dir.y()*dir.y() + dir.z()*dir.z()) * std::sin(angle), 0);
                        return tempV.angle(drift_dir_abs) / 3.1415926 * 180.;
                    };
                    auto pick = [&](const std::map<WirePlaneId, std::pair<geo_point_t, double>>& m,
                                    const WirePlaneId& wpid) {
                        if (wpid.apa() >= 0 && wpid.face() >= 0) {
                            for (const auto& [k, v] : m) {
                                if (k.apa() == wpid.apa() && k.face() == wpid.face()) return v.first;
                            }
                        }
                        return m.begin()->second.first;
                    };
                    auto quasi_parallel = [&](const WirePlaneId& wpid) {
                        const double angle1_1 = angle_deg(pick(wpid_U_dir, wpid));
                        const double angle2_1 = angle_deg(pick(wpid_V_dir, wpid));
                        const double angle3_1 = angle_deg(pick(wpid_W_dir, wpid));
                        const double angle4 =
                            std::fabs(3.1415926/2. - drift_dir_abs.angle(dir)) / 3.1415926 * 180.;
                        return angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5 || angle4 < 5.;
                    };
                    const WirePlaneId wpid1 = m_dv->contained_by(wcp1);
                    const WirePlaneId wpid2 = m_dv->contained_by(wcp2);
                    bool quasi = quasi_parallel(wpid1);
                    if (wpid2.apa() != wpid1.apa() || wpid2.face() != wpid1.face()) {
                        quasi = quasi || quasi_parallel(wpid2);
                    }
                    if (quasi) flag_2view_check = false;
                }
            }

            int num_bad = 0;
            for (size_t i = 0; i != path_wcps_vec1.size(); i++) {
                const auto& p = path_wcps_vec1[i];
                bool flag_reset = false;
                const WirePlaneId wpid = m_dv->contained_by(p);
                if (wpid.apa() < 0 || wpid.face() < 0) {
                    // Outside every sensitive volume (cathode gap / exited):
                    // no hits are possible -- treat like a dead region.
                    flag_reset = true;
                }
                else {
                    const int apa = wpid.apa();
                    const int face = wpid.face();
                    const auto p_raw = to_raw(p, wpid);
                    const size_t nu = grouping->get_closest_points(p_raw, low_dis_limit * 2, apa, face, 0).size();
                    const size_t nv = grouping->get_closest_points(p_raw, low_dis_limit * 2, apa, face, 1).size();
                    const size_t nw = grouping->get_closest_points(p_raw, low_dis_limit * 2, apa, face, 2).size();

                    if (flag_2view_check) {
                        if ((nu > 0 && nv > 0) ||  // require two planes to be good ...
                            (nu > 0 && nw > 0) ||
                            (nv > 0 && nw > 0)) {
                            flag_reset = true;
                        }
                        else if (fiducial_utils->inside_dead_region(p_raw, apa, face)) {
                            flag_reset = true;
                        }
                    }
                    else {
                        if (nu > 0 ||  // require one plane to be good ...
                            nv > 0 ||
                            nw > 0) {
                            flag_reset = true;
                        }
                        else if (fiducial_utils->inside_dead_region(p_raw, apa, face)) {
                            flag_reset = true;
                        }
                    }

                    if (!grouping->is_good_point(p_raw, apa, face)) num_bad++;
                }

                if (flag_reset) {
                    num_nth = 0;
                    min_dis = 1e9;
                    num_bad = 0;
                }
                else {
                    if (inside_fv(p)) {
                        const double dis1 = (p - wcp1).magnitude();
                        const double dis2 = (p - wcp2).magnitude();
                        if (dis1 < min_dis) min_dis = dis1;
                        if (dis2 < min_dis) min_dis = dis2;
                        num_nth++;
                    }
                    if (num_nth > 7 && min_dis < 25*units::cm && num_bad > 7) return true;  // too big a gap
                }
            }
        }

        // Kink search (prototype lines 1834-1983).
        int count = 0;
        double max_angle = 0;
        geo_point_t max_point(0, 0, 0);
        const geo_vector_t drift_dir(1, 0, 0);
        // Apex must be in the FV and not in a dead region for the veto to fire.
        auto apex_in_dead = [&](const geo_point_t& pt) {
            const WirePlaneId wpid = m_dv->contained_by(pt);
            if (wpid.apa() < 0 || wpid.face() < 0) return true;  // outside sensitive volume
            return fiducial_utils->inside_dead_region(to_raw(pt, wpid), wpid.apa(), wpid.face());
        };
        for (size_t i = 5; i + 5 < path_wcps_vec.size(); i++) {
            geo_vector_t dir1 = path_wcps_vec[i] - path_wcps_vec[i - 5];
            geo_vector_t dir2 = path_wcps_vec[i] - path_wcps_vec[i + 5];

            geo_vector_t dir3, dir4, dir5, dir6;
            {
                std::vector<geo_point_t> pts;
                double temp_x = 0, temp_y = 0, temp_z = 0, temp_count = 0;
                for (size_t j = 1; j != 15; j++) {
                    if (i >= j) {
                        const auto& pt = path_wcps_vec[i - j];
                        if (j <= 12 && j > 2) {
                            temp_x += pt.x();
                            temp_y += pt.y();
                            temp_z += pt.z();
                            temp_count++;
                        }
                        pts.push_back(pt);
                    }
                }
                // The prototype recomputes dir3/dir5 on every j iteration;
                // only the final all-points value survives -- computed once.
                dir3 = calc_pca_dir(path_wcps_vec[i], pts);
                dir5 = geo_vector_t(temp_x / temp_count - path_wcps_vec[i].x(),
                                    temp_y / temp_count - path_wcps_vec[i].y(),
                                    temp_z / temp_count - path_wcps_vec[i].z());
                if (dir3.angle(dir1) > 3.1415926/2.) dir3 = dir3 * (-1);
            }
            {
                std::vector<geo_point_t> pts;
                double temp_x = 0, temp_y = 0, temp_z = 0, temp_count = 0;
                for (size_t j = 1; j != 15; j++) {
                    if (i + j < path_wcps_vec.size()) {
                        const auto& pt = path_wcps_vec[i + j];
                        if (j <= 12 && j > 2) {
                            temp_x += pt.x();
                            temp_y += pt.y();
                            temp_z += pt.z();
                            temp_count++;
                        }
                        pts.push_back(pt);
                    }
                }
                dir4 = calc_pca_dir(path_wcps_vec[i], pts);
                dir6 = geo_vector_t(temp_x / temp_count - path_wcps_vec[i].x(),
                                    temp_y / temp_count - path_wcps_vec[i].y(),
                                    temp_z / temp_count - path_wcps_vec[i].z());
                if (dir4.angle(dir2) > 3.1415926/2.) dir4 = dir4 * (-1);
            }

            int cut1 = 0;
            if ((3.1415926 - dir1.angle(dir2)) / 3.1415926 * 180. > 25) cut1++;
            if ((3.1415926 - dir3.angle(dir4)) / 3.1415926 * 180. > 25) cut1++;
            if ((3.1415926 - dir5.angle(dir6)) / 3.1415926 * 180. > 25) cut1++;
            int cut2 = 0;
            if (std::fabs(3.1415926/2. - drift_dir.angle(dir1 - dir2)) / 3.1415926 * 180. > 5) cut2++;
            if (std::fabs(3.1415926/2. - drift_dir.angle(dir3 - dir4)) / 3.1415926 * 180. > 5) cut2++;
            if (std::fabs(3.1415926/2. - drift_dir.angle(dir5 - dir6)) / 3.1415926 * 180. > 5) cut2++;

            if (cut1 >= 3 && cut2 >= 2) {
                if ((3.1415926 - dir3.angle(dir4)) / 3.1415926 * 180. > max_angle) {
                    max_angle = (3.1415926 - dir3.angle(dir4)) / 3.1415926 * 180.;
                    max_point = path_wcps_vec[i];
                }

                count++;
                if (count >= 3) {
                    const geo_vector_t temp1 = path_wcps_vec[i] - wcp1;
                    const geo_vector_t temp2 = path_wcps_vec[i] - wcp2;
                    const double open_deg = (3.1415926 - temp1.angle(temp2)) / 3.1415926 * 180.;
                    const double sym_deg =
                        std::fabs(3.1415926/2. - drift_dir.angle(temp1 + temp2)) / 3.1415926 * 180.;
                    // prototype lines 1935-1939, operator precedence preserved:
                    // G1 || (G2 && >10cm legs) || (G3 && >15cm legs), each
                    // Gn = (open > thresh && sym > 5.5) || open > 60.
                    if (((open_deg > 35 && sym_deg > 5.5) || open_deg > 60) ||
                        (((open_deg > 32 && sym_deg > 5.5) || open_deg > 60)
                         && temp1.magnitude() > 10*units::cm && temp2.magnitude() > 10*units::cm) ||
                        (((open_deg > 25 && sym_deg > 5.5) || open_deg > 60)
                         && temp1.magnitude() > 15*units::cm && temp2.magnitude() > 15*units::cm)) {
                        if ((!inside_fv(max_point)) ||           // must be in fiducial
                            (apex_in_dead(max_point) && open_deg < 45)) {  // not in dead volume
                        }
                        else {
                            return true;
                        }
                    }
                }
                else if (count >= 1) {
                    const geo_vector_t temp1 = path_wcps_vec[i] - wcp1;
                    const geo_vector_t temp2 = path_wcps_vec[i] - wcp2;
                    const double open_deg = (3.1415926 - temp1.angle(temp2)) / 3.1415926 * 180.;
                    const double sym_deg =
                        std::fabs(3.1415926/2. - drift_dir.angle(temp1 + temp2)) / 3.1415926 * 180.;
                    // prototype lines 1957-1961.
                    if (((open_deg > 35 && sym_deg > 5.5) || open_deg > 60)
                        && temp1.magnitude() > 5*units::cm && temp2.magnitude() > 5*units::cm) {
                        if ((!inside_fv(max_point)) ||           // must be in fiducial
                            (apex_in_dead(max_point) && open_deg < 45)) {  // not in dead volume
                        }
                        else {
                            return true;
                        }
                    }
                }
            }
            else {
                count = 0;
                max_angle = 0;
                max_point = geo_point_t(0, 0, 0);
            }
        }

        return false;
    }
};
