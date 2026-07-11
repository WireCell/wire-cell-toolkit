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
///    replaced by a beam-window test on cluster_t0.  The protected branches
///    require the unported check_neutrino_candidate(); v1 is conservative
///    and never tags an in-beam-window bundle through those branches.
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
                       private Clus::NeedDV, private Clus::NeedFiducial {
public:
    TaggerCheckTGM() {}
    virtual ~TaggerCheckTGM() {}

    virtual void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedFiducial::configure(config);
        m_grouping_name = get<std::string>(config, "grouping", m_grouping_name);
        m_beam_window_low = get(config, "beam_window_low", m_beam_window_low);
        m_beam_window_high = get(config, "beam_window_high", m_beam_window_high);
        m_length_limit_frac = get(config, "length_limit_frac", m_length_limit_frac);
        m_enable_case_b = get(config, "enable_case_b", m_enable_case_b);
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
        cfg["fiducial"] = "DetectorVolumes";
        // Boundary margins for the FV tests, FiducialUtils tolerance-vec
        // convention: [x_lo, x_hi, y_lo, y_hi, z_lo, z_hi], negative = inset.
        cfg["fv_tolerance"] = Json::Value(Json::arrayValue);
        cfg["beam_window_low"] = m_beam_window_low;   // window on cluster_t0; low >= high
        cfg["beam_window_high"] = m_beam_window_high; // disables the beam protection
        cfg["length_limit_frac"] = m_length_limit_frac;
        cfg["enable_case_b"] = m_enable_case_b;
        return cfg;
    }

    virtual void visit(Ensemble& ensemble) const {
        auto groupings = ensemble.with_name(m_grouping_name);
        if (groupings.empty()) return;
        auto& grouping = *groupings.at(0);

        std::vector<Cluster*> main_clusters;
        for (auto* cluster : grouping.children()) {
            if (cluster->get_flag(Flags::main_cluster)) main_clusters.push_back(cluster);
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
    bool m_enable_case_b{true};
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
                            // prototype: !check_neutrino_candidate && len > frac*limit
                            // → conservative NO TAG until that port exists.
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
                            if (in_beam_window) continue;  // needs check_neutrino_candidate
                            const double temp_length = (pe1 - pe2).magnitude();
                            if (temp_length > m_length_limit_frac * length_limit) return true;
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
                        if (in_beam_window) continue;  // needs check_neutrino_candidate
                        return true;
                    }
                }
            }
        }
        return false;
    }
};
