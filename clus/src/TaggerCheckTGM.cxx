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
        // check_neutrino_candidate (default false = v1 behavior): enable the
        // ported prototype neutrino-candidate veto in the beam-protected
        // branches.  OFF keeps the conservative never-tag-in-beam behavior
        // byte-identical.
        m_check_neutrino_candidate = get(config, "check_neutrino_candidate", m_check_neutrino_candidate);
        auto tol = config["fv_tolerance"];
        if (!tol.isNull() && tol.isArray()) {
            m_fv_tolerance.clear();
            for (const auto& t : tol) m_fv_tolerance.push_back(t.asDouble());
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
        cfg["beam_window_low"] = m_beam_window_low;   // window on cluster_t0; low >= high
        cfg["beam_window_high"] = m_beam_window_high; // disables the beam protection
        cfg["length_limit_frac"] = m_length_limit_frac;
        cfg["enable_case_b"] = m_enable_case_b;
        cfg["require_in_scope"] = m_require_in_scope;
        cfg["check_neutrino_candidate"] = m_check_neutrino_candidate;
        return cfg;
    }

    virtual void visit(Ensemble& ensemble) const {
        auto groupings = ensemble.with_name(m_grouping_name);
        if (groupings.empty()) return;
        auto& grouping = *groupings.at(0);

        std::vector<Cluster*> main_clusters;
        int n_out_of_scope = 0;
        for (auto* cluster : grouping.children()) {
            if (!cluster->get_flag(Flags::main_cluster)) continue;
            if (m_require_in_scope && !cluster->get_scope_filter(cluster->get_default_scope())) {
                ++n_out_of_scope;
                continue;
            }
            main_clusters.push_back(cluster);
        }
        if (n_out_of_scope) {
            SPDLOG_LOGGER_INFO(t_log, "visit: TaggerCheckTGM: skipped {} out-of-scope main cluster(s)",
                               n_out_of_scope);
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
    double m_length_limit_frac{0.45};
    bool m_require_in_scope{false};
    bool m_enable_case_b{true};
    bool m_check_neutrino_candidate{false};
    std::vector<double> m_fv_tolerance;

    // FiducialUtils::inside_fiducial_volume logic against our own IFiducial.
    bool inside_fv(const Point& p) const {
        if (m_fv_tolerance.empty()) return m_fiducial->contained(p);
        double txl, txh, tyl, tyh, tzl, tzh;
        const auto& tv = m_fv_tolerance;
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

        const auto out_vec_wcps = cluster.get_extreme_wcps();
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
                        flag_check = flag_check || inside_fv(p3);
                    }
                    if (std::getenv("WCT_TGM_DEBUG")) {
                        SPDLOG_LOGGER_INFO(t_log, "check_tgm dbg: cluster {} pair ({},{}) ngrp {} pe1 ({:.1f},{:.1f},{:.1f}) pe2 ({:.1f},{:.1f},{:.1f}) mid_inside {} len {:.1f}/{:.1f} cm",
                            cluster.ident(), i, k, out_vec_wcps.size(),
                            pe1.x()/units::cm, pe1.y()/units::cm, pe1.z()/units::cm,
                            pe2.x()/units::cm, pe2.y()/units::cm, pe2.z()/units::cm,
                            flag_check, (pe1-pe2).magnitude()/units::cm, length_limit/units::cm);
                    }
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
                                flag_check_again = flag_check_again || inside_fv(p3);
                            }
                            for (int kk = 0; kk != 3; kk++) {
                                // NB: mixed endpoints as in the prototype.
                                const geo_point_t p3 = out_vec_wcps[kkk][0]
                                    + (pe2 - out_vec_wcps[i][0]) * ((kk + 1) / 4.);
                                flag_check_again = flag_check_again || inside_fv(p3);
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
