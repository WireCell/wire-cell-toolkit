// TaggerCheckTGM -- Through-Going-Muon (TGM) tagger.
//
// Faithful port of the WCP (wire-cell-prototype) ToyFiducial::check_tgm
// (pid/src/Cosmic_tagger.h, function check_tgm) into the WCT clustering
// IEnsembleVisitor framework.  Follows the IO / framework conventions of
// TaggerCheckSTM.cxx.
//
// SBND differences from the WCP MicroBooNE original:
//   - Runs AFTER the two-APA Q/L matching + stitching stage, on the
//     T0-corrected, stitched all-APA clusters (default scope "3d" with
//     x_t0cor/y_cor/z_cor coords).  cluster_t0 already folded into x_t0cor.
//   - Each endpoint is SCE-corrected back to its true position with
//     SCECorrection::forward(p, t0=0, face, apa) BEFORE the fiducial test
//     (t0=0 because x_t0cor already carries the drift correction; forward then
//     only adds the SCE spatial displacement).  When no ISCEField is wired the
//     transform is a no-op and the test runs on the x_t0cor point, which is the
//     current SBND configuration.
//   - With SCE-corrected (true) points the fiducial boundary is a simple BOX,
//     which is exactly what FiducialUtils::inside_fiducial_volume() tests
//     (m_sd.fiducial->contained(p)); the box bounds come from
//     cfg/pgrapher/experiment/sbnd/clus.jsonnet (dvm.overall) via DetectorVolumes.
//   - ALL live clusters are checked and flagged (not just the main cluster).
//
// NOT ported: WCP check_neutrino_candidate() (Dijkstra path-topology neutrino
// veto) has no WCT equivalent.  Here all flashes are treated as normal (WCP
// flash type != 2): the type-2 neutrino-veto branches are skipped and the
// neutrino-candidate guards are treated as "not a neutrino" (i.e. they pass),
// with the WCP length guard (chord > 0.45*length_limit) retained where it
// gated those branches.  This is the conservative TGM-tagging choice for
// post-Q/L cosmic clusters in a box fiducial volume.

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellClus/FiducialUtils.h"
#include "WireCellClus/IPCTransform.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Units.h"

#include <cmath>
#include <vector>

class TaggerCheckTGM;
WIRECELL_FACTORY(TaggerCheckTGM, TaggerCheckTGM,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

static auto s_log = WireCell::Log::logger("clus.TGM");

class TaggerCheckTGM : public IConfigurable,
                       public Clus::IEnsembleVisitor,
                       private Clus::NeedDV,
                       private Clus::NeedPCTS {
public:
    TaggerCheckTGM() {}
    virtual ~TaggerCheckTGM() {}

    virtual void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedPCTS::configure(config);
        m_grouping_name = get<std::string>(config, "grouping", "live");
        m_debug = get<bool>(config, "debug", false);
        m_debug_charge_array = get<std::string>(config, "debug_charge_array", "tgm_charge");
        m_debug_charge_pcname = get<std::string>(config, "debug_charge_pcname", "tgm_debug");
        m_debug_endpoint_charge = get<double>(config, "debug_endpoint_charge", 10000.0);
        m_debug_body_charge = get<double>(config, "debug_body_charge", 100.0);
    }

    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping_name;
        cfg["detector_volumes"] = "DetectorVolumes";
        cfg["pc_transforms"] = "PCTransformSet";
        // Debug mode: write a per-point charge array on each tagged cluster so
        // MultiAlgBlobClustering can dump a Bee point set highlighting the
        // tagged-track endpoints.  Default OFF.
        cfg["debug"] = m_debug;
        cfg["debug_charge_array"] = m_debug_charge_array;
        cfg["debug_charge_pcname"] = m_debug_charge_pcname;
        cfg["debug_endpoint_charge"] = m_debug_endpoint_charge;
        cfg["debug_body_charge"] = m_debug_body_charge;
        return cfg;
    }

    virtual void visit(Ensemble& ensemble) const {
        auto groupings = ensemble.with_name(m_grouping_name);
        if (groupings.empty()) return;
        auto& grouping = *groupings.at(0);

        auto fiducial_utils = grouping.get_fiducialutils();
        if (!fiducial_utils) {
            SPDLOG_LOGGER_WARN(s_log, "visit: no FiducialUtils on grouping '{}', skipping", m_grouping_name);
            return;
        }

        // SCE transform (no-op when no ISCEField is wired).
        IPCTransform::pointer sce = m_pcts ? m_pcts->pc_transform("SCECorrection") : nullptr;

        int ntagged = 0;
        for (auto* cluster : grouping.children()) {
            if (check_tgm(*cluster, *fiducial_utils, sce)) {
                cluster->set_flag(Flags::TGM);
                ++ntagged;
                if (m_debug) {
                    mark_debug_points(*cluster);
                }
            }
        }
        SPDLOG_LOGGER_DEBUG(s_log, "visit: tagged {}/{} live clusters as TGM",
                            ntagged, (int)grouping.children().size());
    }

private:
    std::string m_grouping_name{"live"};
    bool m_debug{false};
    std::string m_debug_charge_array{"tgm_charge"};
    std::string m_debug_charge_pcname{"tgm_debug"};
    double m_debug_endpoint_charge{10000.0};
    double m_debug_body_charge{100.0};

    // SCE-correct a (default-scope, x_t0cor) point back to true position and box-test.
    bool inside_fv(const geo_point_t& p, const FiducialUtils& fid, IPCTransform::pointer sce) const {
        geo_point_t p_cor = p;
        if (sce) {
            WirePlaneId wpid = m_dv->contained_by(p);
            if (wpid.apa() >= 0 && wpid.face() >= 0) {
                // t0=0: x_t0cor already carries the drift correction; forward()
                // here only adds the SCE spatial displacement.
                p_cor = sce->forward(p, 0.0, wpid.face(), wpid.apa());
            }
        }
        return fid.inside_fiducial_volume(p_cor);
    }

    // Port of WCP ToyFiducial::check_tgm.  Returns true if the cluster is a
    // through-going muon (both ends exit the fiducial volume).
    bool check_tgm(const Cluster& cluster, const FiducialUtils& fid, IPCTransform::pointer sce) const {

        std::vector<std::vector<geo_point_t>> out_vec_wcps = cluster.get_extreme_wcps();
        if (out_vec_wcps.size() < 2) return false;
        for (const auto& g : out_vec_wcps) {
            if (g.empty()) return false;
        }

        const double pi = 3.1415926;
        const geo_point_t drift_dir(1, 0, 0);
        // Hard coded for U and V planes (SBND induction at +/-60 deg, W vertical).
        const geo_point_t U_dir(0, std::cos(60. / 180. * pi),  std::sin(60. / 180. * pi));
        const geo_point_t V_dir(0, std::cos(60. / 180. * pi), -std::sin(60. / 180. * pi));
        const geo_point_t W_dir(0, 1, 0);

        auto& g0 = out_vec_wcps.at(0).at(0);
        auto& g1 = out_vec_wcps.at(1).at(0);
        const double length_limit =
            std::sqrt(std::pow(g0.x() - g1.x(), 2) + std::pow(g0.y() - g1.y(), 2) + std::pow(g0.z() - g1.z(), 2));

        for (size_t i = 0; i != out_vec_wcps.size(); i++) {
            bool flag_p1_inside = true;
            int p1_index = -1;
            for (size_t j = 0; j != out_vec_wcps.at(i).size(); j++) {
                const geo_point_t& p1 = out_vec_wcps.at(i).at(j);
                flag_p1_inside = flag_p1_inside && inside_fv(p1, fid, sce);
                if (!flag_p1_inside) { p1_index = (int)j; break; }
            }

            for (size_t k = i + 1; k != out_vec_wcps.size(); k++) {
                bool flag_p2_inside = true;
                int p2_index = -1;
                for (size_t j = 0; j != out_vec_wcps.at(k).size(); j++) {
                    const geo_point_t& p2 = out_vec_wcps.at(k).at(j);
                    flag_p2_inside = flag_p2_inside && inside_fv(p2, fid, sce);
                    if (!flag_p2_inside) { p2_index = (int)j; break; }
                }

                if ((!flag_p1_inside) && (!flag_p2_inside)) {
                    // CASE A: both endpoint groups appear to exit the FV.
                    const geo_point_t& e1 = out_vec_wcps.at(i).at(p1_index);
                    const geo_point_t& e2 = out_vec_wcps.at(k).at(p2_index);

                    // Sample 3 points between the two exit points; require at
                    // least one to lie inside the FV (track passes through it).
                    bool flag_check = false;
                    for (int kk = 0; kk != 3; kk++) {
                        geo_point_t p3(e1.x() + (e2.x() - e1.x()) / 4. * (kk + 1),
                                       e1.y() + (e2.y() - e1.y()) / 4. * (kk + 1),
                                       e1.z() + (e2.z() - e1.z()) / 4. * (kk + 1));
                        flag_check = flag_check || inside_fv(p3, fid, sce);
                    }

                    if (flag_check) {
                        // Normal flash (no type-2 neutrino veto in WCT) -> TGM.
                        return true;
                    }
                    else {
                        if (out_vec_wcps.size() == 2) {
                            return true;
                        }
                        else {
                            // Re-test passing-through using the OTHER endpoint
                            // groups as intermediate waypoints.
                            bool flag_check_again = false;
                            for (int kkk = 0; kkk != (int)out_vec_wcps.size(); kkk++) {
                                if (kkk == (int)i || kkk == (int)k) continue;
                                const geo_point_t& em = out_vec_wcps.at(kkk).at(0);
                                for (int kk = 0; kk != 4; kk++) {
                                    geo_point_t p3(e1.x() + (em.x() - e1.x()) / 4. * (kk + 1),
                                                   e1.y() + (em.y() - e1.y()) / 4. * (kk + 1),
                                                   e1.z() + (em.z() - e1.z()) / 4. * (kk + 1));
                                    flag_check_again = flag_check_again || inside_fv(p3, fid, sce);
                                }
                                for (int kk = 0; kk != 3; kk++) {
                                    geo_point_t p3(em.x() + (e2.x() - g0.x()) / 4. * (kk + 1),
                                                   em.y() + (e2.y() - g0.y()) / 4. * (kk + 1),
                                                   em.z() + (e2.z() - g0.z()) / 4. * (kk + 1));
                                    flag_check_again = flag_check_again || inside_fv(p3, fid, sce);
                                }
                            }
                            if (!flag_check_again) {
                                double temp_length =
                                    std::sqrt(std::pow(e1.x() - e2.x(), 2) + std::pow(e1.y() - e2.y(), 2) +
                                              std::pow(e1.z() - e2.z(), 2));
                                // (neutrino-candidate guard treated as pass)
                                if (temp_length > 0.45 * length_limit) return true;
                            }
                        }
                    }
                }
                else {
                    // CASE B: at least one endpoint group appears to stay inside.
                    geo_point_t dir_main = cluster.get_pca().axis.at(0);
                    geo_point_t dir_test(g0_at(out_vec_wcps, i).x() - g0_at(out_vec_wcps, k).x(),
                                         g0_at(out_vec_wcps, i).y() - g0_at(out_vec_wcps, k).y(),
                                         g0_at(out_vec_wcps, i).z() - g0_at(out_vec_wcps, k).z());

                    // Only consider near-axial endpoint pairs.
                    if (std::fabs((pi / 2. - dir_test.angle(dir_main)) / pi * 180.) > 75 || (i == 0 && k == 1)) {

                        bool flag_p1_inside_p = flag_p1_inside;
                        if (flag_p1_inside_p) {
                            geo_point_t p1 = g0_at(out_vec_wcps, i);
                            geo_point_t dir = cluster.vhough_transform(p1, 30 * units::cm);
                            dir = dir * (-1.0);
                            if (dir.angle(dir_test) > pi * 2. / 3.) continue;

                            if (prolonged_track(dir, drift_dir, U_dir, V_dir, W_dir, pi)) {
                                flag_p1_inside_p = flag_p1_inside_p &&
                                                   fid.check_signal_processing(cluster, p1, dir, 1 * units::cm);
                            }
                            if (std::fabs((pi / 2. - dir.angle(dir_main)) / pi * 180.) > 60)
                                flag_p1_inside_p = flag_p1_inside_p &&
                                                   fid.check_dead_volume(cluster, p1, dir, 1 * units::cm);
                        }

                        bool flag_p2_inside_p = flag_p2_inside;
                        if (flag_p2_inside_p) {
                            geo_point_t p2 = g0_at(out_vec_wcps, k);
                            geo_point_t dir = cluster.vhough_transform(p2, 30 * units::cm);
                            dir = dir * (-1.0);
                            if (dir.angle(dir_test) < pi / 3.) continue;

                            if (prolonged_track(dir, drift_dir, U_dir, V_dir, W_dir, pi)) {
                                flag_p2_inside_p = flag_p2_inside_p &&
                                                   fid.check_signal_processing(cluster, p2, dir, 1 * units::cm);
                            }
                            if (std::fabs((pi / 2. - dir.angle(dir_main)) / pi * 180.) > 60)
                                flag_p2_inside_p = flag_p2_inside_p &&
                                                   fid.check_dead_volume(cluster, p2, dir, 1 * units::cm);
                        }

                        if ((!flag_p1_inside_p) && (!flag_p2_inside_p)) {
                            // Normal flash (no type-2 neutrino veto) -> TGM.
                            return true;
                        }
                    }
                }
            }
        }

        return false;
    }

    // Convenience: first point of endpoint group g.
    static const geo_point_t& g0_at(const std::vector<std::vector<geo_point_t>>& v, size_t g) {
        return v.at(g).at(0);
    }

    // WCP prolonged-track / isochronous test: project the Hough direction onto
    // each wire plane and check whether the track is close to isochronous in
    // U (<10 deg), V (<10 deg) or W (<5 deg).
    static bool prolonged_track(const geo_point_t& dir, const geo_point_t& drift_dir,
                                const geo_point_t& U_dir, const geo_point_t& V_dir,
                                const geo_point_t& W_dir, double pi) {
        const double yz = std::sqrt(dir.y() * dir.y() + dir.z() * dir.z());
        geo_point_t dir_1(0, dir.y(), dir.z());

        double angle1 = dir_1.angle(U_dir);
        geo_point_t tempV1(std::fabs(dir.x()), yz * std::sin(angle1), 0);
        double angle1_1 = tempV1.angle(drift_dir) / pi * 180.;

        double angle2 = dir_1.angle(V_dir);
        geo_point_t tempV2(std::fabs(dir.x()), yz * std::sin(angle2), 0);
        double angle2_1 = tempV2.angle(drift_dir) / pi * 180.;

        double angle3 = dir_1.angle(W_dir);
        geo_point_t tempV3(std::fabs(dir.x()), yz * std::sin(angle3), 0);
        double angle3_1 = tempV3.angle(drift_dir) / pi * 180.;

        return (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5);
    }

    // Debug: write a per-point charge array on a tagged cluster -- endpoints
    // (the two main-axis extreme points) get the endpoint charge, all other
    // points the body charge.  Consumed by MultiAlgBlobClustering's bee dump
    // (a bee_points_sets entry with charge_array == m_debug_charge_array).
    void mark_debug_points(Cluster& cluster) const {
        const int n = cluster.npoints();
        if (n <= 0) return;
        std::vector<double> q(n, m_debug_body_charge);

        std::vector<std::vector<geo_point_t>> ev = cluster.get_extreme_wcps();
        if (ev.size() >= 2 && !ev.at(0).empty() && !ev.at(1).empty()) {
            for (size_t g : {size_t(0), size_t(1)}) {
                size_t idx = cluster.get_closest_point_index(ev.at(g).at(0));
                if ((int)idx < n) q[idx] = m_debug_endpoint_charge;
            }
        }
        cluster.put_pcarray(q, m_debug_charge_array, m_debug_charge_pcname);
    }
};
