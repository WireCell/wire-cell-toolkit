#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRShowerFunctions.h"
#include "WireCellUtil/Units.h"

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

using ChargeMap = std::map<TrackFitting::CoordReadout, TrackFitting::ChargeMeasurement>;
using WireMap   = std::map<std::pair<int,int>, std::vector<std::tuple<int,int,int>>>;

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

    for (int plane_id = 0; plane_id < 3; ++plane_id) {
        for (const auto& [coord_key, charge_data] : *maps[plane_id]) {
            int time_slice = coord_key.time;
            int channel    = coord_key.channel;
            int apa        = coord_key.apa;

            auto wire_it = map_apa_ch_plane_wires.find({apa, channel});
            if (wire_it == map_apa_ch_plane_wires.end()) continue;

            int face = -1;
            for (const auto& [f, plane, wire] : wire_it->second) {
                if (plane == plane_id) { face = f; break; }
            }
            if (face < 0) continue;

            auto p2d = grouping->convert_time_ch_2Dpoint(time_slice, channel, apa, face, plane_id);
            WireCell::Point test_p2d(p2d.first, p2d.second, 0);

            double dis = 1e9;
            size_t point_index = 0;
            const Facade::Cluster* closest_cluster = nullptr;

            if (pcloud1) {
                auto res    = pcloud1->get_closest_2d_point_info(test_p2d, plane_id, face, apa);
                dis             = std::get<0>(res);
                closest_cluster = std::get<1>(res);
                point_index     = std::get<2>(res);
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
                auto res    = pcloud2->get_closest_2d_point_info(test_p2d, plane_id, face, apa);
                dis             = std::get<0>(res);
                closest_cluster = std::get<1>(res);
                point_index     = std::get<2>(res);
                try_add(pcloud2);
            }
        }
    }

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


double PatternAlgorithms::cal_kine_charge(ShowerPtr shower, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    (void)graph;
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

    ChargeMap charge_2d_u, charge_2d_v, charge_2d_w;
    WireMap   map_apa_ch_plane_wires;
    track_fitter.collect_2D_charge(charge_2d_u, charge_2d_v, charge_2d_w, map_apa_ch_plane_wires);

    auto pcloud1 = shower->get_pcloud("associate_points");
    auto pcloud2 = shower->get_pcloud("fit");
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
    if (!grouping) return;

    // Collect 2D charge maps once for all showers.  The 2D charge data is
    // event-global (not shower-specific), so collecting inside the per-shower
    // loop would redundantly repeat this O(N_hits) operation N_shower times.
    ChargeMap charge_2d_u, charge_2d_v, charge_2d_w;
    WireMap   map_apa_ch_plane_wires;
    track_fitter.collect_2D_charge(charge_2d_u, charge_2d_v, charge_2d_w, map_apa_ch_plane_wires);

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
        if (!pcloud1 && !pcloud2) { shower->set_flag_kinematics(true); continue; }
        if (!pcloud1) pcloud1 = pcloud2;
        if (!pcloud2) pcloud2 = pcloud1;

        double kine_charge = kine_charge_from_maps(
            pcloud1, pcloud2, fudge_factor, recom_factor,
            charge_2d_u, charge_2d_v, charge_2d_w, map_apa_ch_plane_wires,
            grouping, corr_fn, dis_cut);

        shower->set_kine_charge(kine_charge);
        shower->set_flag_kinematics(true);
    }
}
