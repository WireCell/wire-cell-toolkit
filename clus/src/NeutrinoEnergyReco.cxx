#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRShowerFunctions.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Logging.h"

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;
using namespace WireCell;


double PatternAlgorithms::cal_corr_factor(WireCell::Point& pt, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    double corr_factor = 1.0;
    // So far this is an empty class that needs to be filled with actual logic ...

    // Example 1: Find APA and face using detector volumes
    // The WirePlaneId contains apa and face information
    WirePlaneId wpid = dv->contained_by(pt);
    int apa = wpid.apa();
    int face = wpid.face();
    int plane = wpid.index();  // 0=U, 1=V, 2=W

    // Example 2: Find the grouping from track_fitter
    // The track_fitter contains a reference to the grouping
    auto grouping = track_fitter.grouping();

    (void)apa;
    (void)face;
    (void)plane;
    (void)grouping;

    return corr_factor;
}


namespace {

// ChargeMap and WireMap are defined in NeutrinoPatternBase.h (WireCell::Clus::PR namespace).
using WireCell::Clus::PR::ChargeMap;
using WireCell::Clus::PR::WireMap;

// Core charge-to-energy conversion given pre-collected 2D charge maps and point clouds.
// CorrFn: callable with signature double(WireCell::Point&).
// Both cal_kine_charge overloads and calculate_shower_kinematics use this.
template<typename CorrFn>
static double kine_charge_from_maps(
    std::shared_ptr<Facade::DynamicPointCloud> pcloud1,
    std::shared_ptr<Facade::DynamicPointCloud> pcloud2,
    double fudge_factor,
    double recom_factor,
    const ChargeMap& charge_2d_u,
    const ChargeMap& charge_2d_v,
    const ChargeMap& charge_2d_w,
    const WireMap& map_apa_ch_plane_wires,
    Facade::Grouping* grouping,
    CorrFn&& corr_fn,
    double dis_cut)
{
    const ChargeMap* maps[3] = {&charge_2d_u, &charge_2d_v, &charge_2d_w};
    double sums[3] = {0, 0, 0};
    int n_hits_total = 0, n_hits_within_cut = 0;

    for (int plane_id = 0; plane_id < 3; ++plane_id) {
        int plane_sample = 0;  // print first 3 hits per plane for diagnosis
        for (const auto& [coord_key, charge_data] : *maps[plane_id]) {
            int time_slice = coord_key.time;
            int channel    = coord_key.channel;
            int apa        = coord_key.apa;

            auto wire_it = map_apa_ch_plane_wires.find({apa, channel});
            if (wire_it == map_apa_ch_plane_wires.end()) continue;

            int face = -1, local_wire = -1;
            for (const auto& [f, plane, wire] : wire_it->second) {
                if (plane == plane_id) { face = f; local_wire = wire; break; }
            }
            if (face < 0 || local_wire < 0) continue;

            // Use local_wire (plane-local index from map_apa_ch_plane_wires), NOT the global
            // channel number. V channels start at ~2400 and W at ~4800, so passing channel
            // as the wire index gives enormous wrong p2d.second for V/W planes.
            auto p2d = grouping->convert_time_wire_2Dpoint(time_slice, local_wire, apa, face, plane_id);

            double dis = 1e9;
            size_t point_index = 0;
            const Facade::Cluster* closest_cluster = nullptr;

            // Use direct (drift, wire_perp) query — p2d already provides coordinates in the
            // wire-perpendicular space matching the KD2D tree storage; applying the angle
            // projection again in get_closest_2d_point_info would corrupt the coordinates.
            if (pcloud1) {
                auto res    = pcloud1->get_closest_2d_point_info_direct(p2d.first, p2d.second, plane_id, face, apa);
                dis             = std::get<0>(res);
                closest_cluster = std::get<1>(res);
                point_index     = std::get<2>(res);
            }

            ++n_hits_total;
            if (dis < dis_cut && closest_cluster) ++n_hits_within_cut;

            if (plane_sample < 3) {
                ++plane_sample;
                auto pcloud1_pts = pcloud1 ? pcloud1->get_points().size() : 0;
                // Also sample first point in pcloud1 for coordinate comparison
                double pc_x = 0, pc_y = 0, pc_z = 0;
                if (pcloud1 && pcloud1_pts > 0) {
                    const auto& pt0 = pcloud1->get_points()[0];
                    pc_x = pt0.x; pc_y = pt0.y; pc_z = pt0.z;
                }
                SPDLOG_LOGGER_DEBUG(s_log,
                    "kine_charge_from_maps:   plane={} ts={} ch={} wire={} apa={} face={}"
                    " p2d=({:.3f},{:.3f}) dis={:.4f}cm dis_cut={:.4f}cm cluster={}"
                    " pc0_3d=({:.3f},{:.3f},{:.3f})cm",
                    plane_id, time_slice, channel, local_wire, apa, face,
                    p2d.first / units::cm, p2d.second / units::cm,
                    dis / units::cm, dis_cut / units::cm,
                    closest_cluster ? "ok" : "null",
                    pc_x / units::cm, pc_y / units::cm, pc_z / units::cm);
            }

            // Accumulate charge from pc at the already-looked-up point_index/closest_cluster/dis.
            // Returns true if the charge was added.
            auto try_add = [&](std::shared_ptr<Facade::DynamicPointCloud> pc) -> bool {
                if (!pc || dis >= dis_cut || !closest_cluster) return false;
                const auto& pts = pc->get_points();
                if (point_index >= pts.size()) return false;
                WireCell::Point tp(pts[point_index].x, pts[point_index].y, pts[point_index].z);
                sums[plane_id] += charge_data.charge * corr_fn(tp);
                return true;
            };

            if (!try_add(pcloud1) && pcloud2) {
                auto res    = pcloud2->get_closest_2d_point_info_direct(p2d.first, p2d.second, plane_id, face, apa);
                dis             = std::get<0>(res);
                closest_cluster = std::get<1>(res);
                point_index     = std::get<2>(res);
                try_add(pcloud2);
            }
        }
    }

    SPDLOG_LOGGER_DEBUG(s_log,
        "kine_charge_from_maps:   hits total={} within_cut={} dis_cut={:.4f}cm sums=[{:.1f},{:.1f},{:.1f}]",
        n_hits_total, n_hits_within_cut, dis_cut / units::cm,
        sums[0], sums[1], sums[2]);

    // Find min / med / max plane indices by charge.
    int min_idx = 0, max_idx = 0, med_idx = 0;
    double min_q = 1e9, max_q = -1e9;
    for (int i = 0; i < 3; ++i) {
        if (sums[i] < min_q) { min_q = sums[i]; min_idx = i; }
        if (sums[i] > max_q) { max_q = sums[i]; max_idx = i; }
    }
    if (min_idx != max_idx) {
        for (int i = 0; i < 3; ++i) {
            if (i != min_idx && i != max_idx) { med_idx = i; break; }
        }
    } else {
        min_idx = 0; med_idx = 1; max_idx = 2;
    }

    const double weight[3] = {0.25, 0.25, 1.0};
    const double weight_sum = weight[0] + weight[1] + weight[2];

    double max_asy = 0;
    if (sums[med_idx] + sums[max_idx] > 0)
        max_asy = std::abs(sums[med_idx] - sums[max_idx]) / (sums[med_idx] + sums[max_idx]);

    double overall = (weight[0]*sums[0] + weight[1]*sums[1] + weight[2]*sums[2]) / weight_sum;
    if (max_asy > 0.04)
        overall = (weight[med_idx]*sums[med_idx] + weight[min_idx]*sums[min_idx])
                  / (weight[med_idx] + weight[min_idx]);

    return overall / recom_factor / fudge_factor * 23.6 / 1e6 * units::MeV;
}

} // anonymous namespace


double PatternAlgorithms::cal_kine_charge(ShowerPtr shower,
    const ChargeMap& charge_2d_u, const ChargeMap& charge_2d_v, const ChargeMap& charge_2d_w,
    const WireMap& map_apa_ch_plane_wires,
    TrackFitting& track_fitter, IDetectorVolumes::pointer dv)
{
    if (!shower) return 0.0;
    auto grouping = track_fitter.grouping();
    if (!grouping) return 0.0;

    double fudge_factor = 0.95, recom_factor = 0.7;
    if (shower->get_flag_shower()) {
        recom_factor = 0.5;
        fudge_factor = 0.8;
    } else if (std::abs(shower->get_particle_type()) == 2212) {
        recom_factor = 0.35;
    }

    auto pcloud1 = shower->get_pcloud("associate_points");
    auto pcloud2 = shower->get_pcloud("fit");
    if (!pcloud1 && !pcloud2) return 0.0;
    if (!pcloud1) pcloud1 = pcloud2;
    if (!pcloud2) pcloud2 = pcloud1;

    return kine_charge_from_maps(
        pcloud1, pcloud2, fudge_factor, recom_factor,
        charge_2d_u, charge_2d_v, charge_2d_w, map_apa_ch_plane_wires,
        grouping,
        [&](WireCell::Point& pt) { return cal_corr_factor(pt, track_fitter, dv); },
        0.6 * units::cm);
}


void PatternAlgorithms::collect_charge_maps(TrackFitting& track_fitter)
{
    m_charge_2d_u.clear();
    m_charge_2d_v.clear();
    m_charge_2d_w.clear();
    m_map_apa_ch_plane_wires.clear();
    track_fitter.collect_2D_charge(m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires);
}


double PatternAlgorithms::cal_kine_charge(ShowerPtr shower, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    (void)graph;
    if (!shower) return 0.0;
    if (!track_fitter.grouping()) return 0.0;

    // Use cached maps if available, otherwise collect (standalone call outside shower_clustering_with_nv).
    if (m_charge_2d_u.empty()) collect_charge_maps(track_fitter);
    return cal_kine_charge(shower, m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires, track_fitter, dv);
}


double PatternAlgorithms::cal_kine_charge(SegmentPtr segment, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    (void)graph;
    if (!segment) return 0.0;

    auto grouping = track_fitter.grouping();
    if (!grouping) return 0.0;

    double fudge_factor = 0.95, recom_factor = 0.7;
    if (segment->flags_any(PR::SegmentFlags::kShowerTopology)) {
        recom_factor = 0.5;
        fudge_factor = 0.8;
    } else if (segment->has_particle_info() && std::abs(segment->particle_info()->pdg()) == 2212) {
        recom_factor = 0.35;
    }

    ChargeMap charge_2d_u, charge_2d_v, charge_2d_w;
    WireMap   map_apa_ch_plane_wires;
    track_fitter.collect_2D_charge(charge_2d_u, charge_2d_v, charge_2d_w, map_apa_ch_plane_wires);

    auto pcloud1 = segment->dpcloud("associate_points");
    auto pcloud2 = segment->dpcloud("fit");
    if (!pcloud1 && !pcloud2) return 0;
    if (!pcloud1) pcloud1 = pcloud2;
    if (!pcloud2) pcloud2 = pcloud1;

    return kine_charge_from_maps(
        pcloud1, pcloud2, fudge_factor, recom_factor,
        charge_2d_u, charge_2d_v, charge_2d_w, map_apa_ch_plane_wires,
        grouping,
        [&](WireCell::Point& pt) { return cal_corr_factor(pt, track_fitter, dv); },
        0.6 * units::cm);
}


void PatternAlgorithms::calculate_shower_kinematics(std::set<ShowerPtr>& showers, std::set<VertexPtr>& vertices_in_long_muon, std::set<SegmentPtr>& segments_in_long_muon, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
    (void)vertices_in_long_muon;
    (void)graph;

    auto grouping = track_fitter.grouping();
    if (!grouping) {
        SPDLOG_LOGGER_DEBUG(s_log, "calculate_shower_kinematics: grouping is null, returning early");
        return;
    }

    // Use pre-collected maps if available (populated once by shower_clustering_with_nv
    // via collect_charge_maps()), otherwise collect here for standalone calls.
    if (m_charge_2d_u.empty()) collect_charge_maps(track_fitter);

    SPDLOG_LOGGER_DEBUG(s_log,
        "calculate_shower_kinematics: {} shower(s), charge maps U={} V={} W={} hits",
        showers.size(), m_charge_2d_u.size(), m_charge_2d_v.size(), m_charge_2d_w.size());

    auto corr_fn = [&](WireCell::Point& pt) { return cal_corr_factor(pt, track_fitter, dv); };
    const double dis_cut = 0.6 * units::cm;

    for (auto& shower : showers) {
        if (!shower || shower->get_flag_kinematics()) continue;

        if (std::abs(shower->get_particle_type()) != 13) {
            shower->calculate_kinematics(particle_data, recomb_model);
        } else {
            shower->calculate_kinematics_long_muon(segments_in_long_muon, particle_data, recomb_model);
        }

        double fudge_factor = 0.95, recom_factor = 0.7;
        if (shower->get_flag_shower()) {
            recom_factor = 0.5;
            fudge_factor = 0.8;
        } else if (std::abs(shower->get_particle_type()) == 2212) {
            recom_factor = 0.35;
        }

        auto pcloud1 = shower->get_pcloud("associate_points");
        auto pcloud2 = shower->get_pcloud("fit");
        if (!pcloud1 && !pcloud2) {
            SPDLOG_LOGGER_DEBUG(s_log,
                "calculate_shower_kinematics:   shower pdg={} nseg={} — no pclouds, kine_charge=0",
                shower->get_particle_type(), shower->get_num_segments());
            shower->set_flag_kinematics(true);
            continue;
        }
        if (!pcloud1) pcloud1 = pcloud2;
        if (!pcloud2) pcloud2 = pcloud1;

        SPDLOG_LOGGER_DEBUG(s_log,
            "calculate_shower_kinematics:   shower pdg={} nseg={} pcloud1_pts={} pcloud2_pts={}",
            shower->get_particle_type(), shower->get_num_segments(),
            pcloud1->get_points().size(), pcloud2->get_points().size());

        double kine_charge = kine_charge_from_maps(
            pcloud1, pcloud2, fudge_factor, recom_factor,
            m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires,
            grouping, corr_fn, dis_cut);

        SPDLOG_LOGGER_DEBUG(s_log,
            "calculate_shower_kinematics:   shower pdg={} nseg={} kine_charge={:.1f}MeV",
            shower->get_particle_type(), shower->get_num_segments(), kine_charge / units::MeV);

        shower->set_kine_charge(kine_charge);
        shower->set_flag_kinematics(true);
    }
}
